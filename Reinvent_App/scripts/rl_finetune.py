"""REINFORCE-style RL fine-tune of the REINVENT 4 prior on a target-aware
reward stack.

This is the canonical REINVENT RL algorithm (Olivecrona et al. 2017):

    1. Two copies of the LSTM: ``prior`` (frozen) and ``agent`` (trainable).
    2. Each step:
         a) Agent samples N SMILES (no_grad). Get sequences + agent NLLs.
         b) Compute reward R(smiles) ∈ [0, 1] for each sample.
         c) Re-evaluate the agent's NLL FOR THE SAMPLED SEQUENCES with
            gradients enabled (this is the policy gradient path).
         d) Compute the prior's NLL for the same sequences (frozen).
         e) Augmented likelihood:  L_aug = -prior_NLL + sigma * R
         f) Loss:                  ((L_aug - (-agent_NLL))**2).mean()
            i.e. the agent's log-prob should track (prior_log_prob + sigma*R).
         g) Backprop + Adam step.
    3. Periodically: log mean reward, validity %, mean QED, save checkpoint.

Reward stack (free, runs locally, no external service required):

    R = w_qed * QED  +  w_pains * (1 if not PAINS) +  w_tani * max_tani_to_actives

with weights normalised to sum to 1. Each component lands in [0, 1].

If you want ADMET-in-the-loop or Vina-in-the-loop the trivial extension is to
hook them into _compute_reward as additional terms, gated on a CLI flag.
They are NOT enabled by default because each adds 10+ seconds per molecule
× batch_size × steps = hours of wall-clock — not the right default for
"researcher tries it on Tuesday afternoon".

Usage:
    python scripts/rl_finetune.py EGFR \
        --steps 200 --batch-size 64 --sigma 120 --lr 1e-4 \
        --weight-qed 1 --weight-tani 1 --weight-pains 1 \
        --device cpu --max-len 256

Inputs:
    prior/reinvent.prior                            (base prior; starting weights for BOTH copies)
    prior/training_data/<TARGET>_actives.smi        (optional Tanimoto reference set)

Outputs:
    prior/finetuned/<TARGET>_rl_<UTC_timestamp>.prior      (RL-fine-tuned agent)
    prior/finetuned/<TARGET>_rl_<UTC_timestamp>_log.json   (per-step metrics)

The output file is intentionally NOT named ``<TARGET>_finetuned.prior``
because that path is the supervised-TL artefact slot — overwriting it
would silently swap models out from under the De Novo pipeline. After
manual evaluation the user can rename / promote.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import types
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

# Force UTF-8 console output (Turkish Windows cp1254 trips on Unicode arrows)
try:
    sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
    sys.stderr.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
except Exception:
    pass

import torch

_PKG_ROOT = Path(__file__).resolve().parent.parent
_PRIOR_PATH = _PKG_ROOT / "prior" / "reinvent.prior"
_TRAINING_DIR = _PKG_ROOT / "prior" / "training_data"
_FINETUNED_DIR = _PKG_ROOT / "prior" / "finetuned"

if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))


def _install_pickle_aliases() -> None:
    """Same shim as the sampler — let torch.load resolve REINVENT 4 pickles
    against our local reinvent_core."""
    if "reinvent.models.reinvent.models.vocabulary" in sys.modules:
        return
    from reinvent_core import vocabulary, meta_data

    def _stub(name: str) -> types.ModuleType:
        if name in sys.modules:
            return sys.modules[name]
        m = types.ModuleType(name)
        m.__path__ = []
        sys.modules[name] = m
        return m

    _stub("reinvent")
    _stub("reinvent.models")
    sys.modules["reinvent.models.meta_data"] = meta_data
    _stub("reinvent.models.reinvent")
    _stub("reinvent.models.reinvent.models")
    sys.modules["reinvent.models.reinvent.models.vocabulary"] = vocabulary


# ---------------------------------------------------------------------------
# Reward components
# ---------------------------------------------------------------------------
#
# Each component returns a float in [0, 1]. Invalid SMILES short-circuit
# to 0.0 before any of these run (handled by the caller). Components are
# pure functions of (mol, optional ref_fps); no global state.


def _qed_score(mol) -> float:
    from rdkit.Chem import QED  # local import keeps cold-start light

    try:
        return float(QED.qed(mol))
    except Exception:  # noqa: BLE001
        return 0.0


_PAINS_CATALOG = None  # lazily-initialised PAINS_A+B+C filter


def _pains_pass(mol) -> float:
    """Return 1.0 when the molecule passes PAINS, 0.0 when it hits an alert."""
    global _PAINS_CATALOG
    if _PAINS_CATALOG is None:
        from rdkit.Chem import FilterCatalog

        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        _PAINS_CATALOG = FilterCatalog.FilterCatalog(params)
    try:
        # HasMatch returns True when ANY catalog filter matches → fail.
        return 0.0 if _PAINS_CATALOG.HasMatch(mol) else 1.0
    except Exception:  # noqa: BLE001
        return 1.0


def _max_tanimoto(mol, ref_fps) -> float:
    """Max Tanimoto similarity of `mol` to any reference fingerprint."""
    if not ref_fps:
        return 0.5  # neutral when no references provided
    try:
        from rdkit.Chem import AllChem
        from rdkit import DataStructs

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        best = 0.0
        for ref in ref_fps:
            t = DataStructs.TanimotoSimilarity(fp, ref)
            if t > best:
                best = t
        return float(best)
    except Exception:  # noqa: BLE001
        return 0.0


def _build_reference_fingerprints(smiles_list: list[str]):
    """Compute Morgan ECFP4 fingerprints for the Tanimoto reference set."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    fps = []
    for s in smiles_list:
        s = (s or "").strip()
        if not s:
            continue
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048))
    return fps


# ---------------------------------------------------------------------------
# Reward composition
# ---------------------------------------------------------------------------
def _compose_reward(
    smiles_batch: list[str],
    weights: dict[str, float],
    ref_fps,
) -> tuple[torch.Tensor, dict[str, float]]:
    """Per-SMILES reward + a small per-step stats dict for the log.

    Composition: weighted ARITHMETIC mean over components (NOT geometric).
    Geometric mean would zero out the whole reward when any one component
    is zero (a single PAINS alert collapses the entire signal). Arithmetic
    is more forgiving and lets the policy gradient still move the model
    in a useful direction when one component is hard to hit early in
    training.
    """
    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")

    norm = sum(max(0.0, w) for w in weights.values()) or 1.0
    w_qed = max(0.0, weights.get("qed", 0.0)) / norm
    w_pains = max(0.0, weights.get("pains", 0.0)) / norm
    w_tani = max(0.0, weights.get("tani", 0.0)) / norm

    rewards = []
    valid = 0
    qed_sum = 0.0
    pains_pass = 0
    tani_sum = 0.0
    for smi in smiles_batch:
        # An empty / non-string sample never survives the LSTM decoder but
        # guard defensively — the sampler can emit truncated junk if the
        # max sequence length is reached without an end token.
        if not smi or not isinstance(smi, str):
            rewards.append(0.0)
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rewards.append(0.0)
            continue
        valid += 1
        r_q = _qed_score(mol)
        r_p = _pains_pass(mol)
        r_t = _max_tanimoto(mol, ref_fps)
        qed_sum += r_q
        if r_p > 0.5:
            pains_pass += 1
        tani_sum += r_t
        rewards.append(w_qed * r_q + w_pains * r_p + w_tani * r_t)
    bs = max(1, len(smiles_batch))
    stats = {
        "validity": valid / bs,
        "mean_qed": qed_sum / max(1, valid),
        "pains_pass_rate": pains_pass / max(1, valid),
        "mean_tanimoto": tani_sum / max(1, valid),
    }
    return torch.tensor(rewards, dtype=torch.float32), stats


# ---------------------------------------------------------------------------
# RL training loop
# ---------------------------------------------------------------------------
def _sample_with_grad(model, batch_size: int):
    """Sample sequences from ``model`` AND return their NLLs with gradients.

    REINVENT's vanilla ``Model.sample`` is wrapped in @torch.no_grad which
    is fine for inference but blocks the policy gradient. We re-implement
    the same auto-regressive sampling here without the no_grad wrapper so
    the agent NLLs carry the backward graph the optimizer needs.
    """
    from reinvent_core import vocabulary as mv

    nll_loss = torch.nn.NLLLoss(reduction="none")
    sequences = [
        torch.full(
            (batch_size, 1),
            model.vocabulary[mv.START_TOKEN],
            dtype=torch.long,
            device=model.device,
        )
    ]
    input_vector = torch.full(
        (batch_size,), model.vocabulary[mv.START_TOKEN],
        dtype=torch.long, device=model.device,
    )
    hidden_state = None
    nlls = torch.zeros(batch_size, device=model.device)

    for _ in range(model.max_sequence_length - 1):
        logits, hidden_state = model.network(input_vector.unsqueeze(1), hidden_state)
        logits = logits.squeeze(1)
        log_probs = logits.log_softmax(dim=1)
        probabilities = logits.softmax(dim=1)
        # Multinomial sample is non-differentiable — that's OK; the gradient
        # path flows through the log_prob of the SAMPLED token (REINFORCE).
        input_vector = torch.multinomial(probabilities.detach(), num_samples=1).view(-1)
        sequences.append(input_vector.view(-1, 1))
        nlls = nlls + nll_loss(log_probs, input_vector)
        # Early termination is unsafe here because every position needs to
        # contribute to the gradient — we keep going to max length and let
        # the end-of-seq token padding-style behaviour handle the tail.
        if input_vector.sum() == 0 and not model.network.training:
            break

    concat = torch.cat(sequences, dim=1)
    smiles = [
        model.tokenizer.untokenize(model.vocabulary.decode(seq))
        for seq in concat.detach().cpu().numpy()
    ]
    return concat, smiles, nlls


def _save_checkpoint(model, path: Path) -> None:
    """Mirror finetune_target.py's persistence shape so the sampler can
    load this checkpoint with no special handling."""
    ckpt = {
        "model_type": "Reinvent",
        "version": 1,
        "vocabulary": model.vocabulary,
        "tokenizer": model.tokenizer,
        "meta_data": model.meta_data,
        "max_sequence_length": model.max_sequence_length,
        "network": model.network.state_dict(),
        "network_params": getattr(model.network, "_network_params", {})
        or {
            "layer_size": getattr(model.network, "_layer_size", 512),
            "num_layers": getattr(model.network, "_num_layers", 3),
            "cell_type": getattr(model.network, "_cell_type", "lstm"),
            "embedding_layer_size": getattr(model.network, "_embedding_layer_size", 256),
            "dropout": getattr(model.network, "_dropout", 0.0),
            "layer_normalization": getattr(model.network, "_layer_normalization", False),
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(ckpt, path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="REINFORCE RL fine-tune of REINVENT 4 prior on a target-aware reward.",
    )
    parser.add_argument("target", help="Target name (case-insensitive). Used for output filename + Tanimoto reference lookup.")
    parser.add_argument("--steps", type=int, default=200, help="Number of policy-gradient updates.")
    parser.add_argument("--batch-size", type=int, default=64, help="SMILES per RL step.")
    parser.add_argument("--sigma", type=float, default=120.0, help="REINFORCE augmentation strength.")
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--weight-qed", type=float, default=1.0)
    parser.add_argument("--weight-pains", type=float, default=1.0)
    parser.add_argument("--weight-tani", type=float, default=1.0)
    parser.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    parser.add_argument("--max-len", type=int, default=256)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--checkpoint-every", type=int, default=50,
        help="Save intermediate checkpoint every N steps (final is always saved).",
    )
    parser.add_argument(
        "--log-path", type=str, default=None,
        help="Optional explicit path for the training-log JSON. Defaults next to the checkpoint.",
    )
    parser.add_argument(
        "--out-path", type=str, default=None,
        help="Optional explicit path for the final checkpoint. Defaults to prior/finetuned/<TARGET>_rl_<utc>.prior.",
    )
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    target = args.target.strip().upper()
    device = torch.device(args.device)

    # 1. Load Tanimoto reference set if available — silently neutral
    # (returns 0.5 for every molecule) when missing, so a target with no
    # actives can still be fine-tuned on QED + PAINS alone.
    smi_path = _TRAINING_DIR / f"{target}_actives.smi"
    ref_smiles: list[str] = []
    if smi_path.is_file():
        ref_smiles = [s.strip() for s in smi_path.read_text(encoding="utf-8").splitlines() if s.strip()]
    print(f"[rl] target={target} ref_smiles={len(ref_smiles)} steps={args.steps} batch={args.batch_size} sigma={args.sigma} lr={args.lr} device={device}")
    ref_fps = _build_reference_fingerprints(ref_smiles) if ref_smiles else []

    # 2. Load TWO copies of the prior — agent gets trained, prior is frozen.
    _install_pickle_aliases()
    from reinvent_core.model import Model

    print(f"[rl] loading prior from {_PRIOR_PATH}")
    ckpt = torch.load(_PRIOR_PATH, map_location=device, weights_only=False)
    prior = Model.create_from_dict(ckpt, mode="inference", device=device)
    agent = Model.create_from_dict(ckpt, mode="training", device=device)
    # Switch agent to training mode but keep prior eval (no dropout etc)
    agent.network.train()
    prior.network.eval()
    for p in prior.network.parameters():
        p.requires_grad_(False)

    # Optimizer only on the AGENT's parameters.
    optimizer = torch.optim.Adam(agent.network.parameters(), lr=args.lr)

    weights = {"qed": args.weight_qed, "pains": args.weight_pains, "tani": args.weight_tani}

    # 3. Output paths.
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    if args.out_path:
        out_path = Path(args.out_path)
    else:
        out_path = _FINETUNED_DIR / f"{target}_rl_{stamp}.prior"
    log_path = Path(args.log_path) if args.log_path else out_path.with_suffix(".json").with_name(out_path.stem + "_log.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[rl] checkpoint => {out_path}")
    print(f"[rl] log        => {log_path}")

    history: list[dict] = []
    t0 = time.time()
    best_mean_reward = -math.inf

    for step in range(1, args.steps + 1):
        # 3a. Sample with gradients enabled (REINFORCE).
        seqs, smiles, agent_nlls = _sample_with_grad(agent, args.batch_size)
        # 3b. Compute reward.
        rewards_t, comp_stats = _compose_reward(smiles, weights, ref_fps)
        rewards_t = rewards_t.to(device)
        # 3c. Prior NLLs — frozen, no grad.
        with torch.no_grad():
            prior_nlls = prior.likelihood(seqs)
        # 3d. REINFORCE augmented loss.
        # The standard form is L_aug = -prior_NLL + sigma * R.
        # The agent's log_prob = -agent_NLL should match that target.
        # MSE between (-agent_NLL) and (-prior_NLL + sigma * R).
        augmented = -prior_nlls + args.sigma * rewards_t
        agent_logp = -agent_nlls
        loss = ((augmented - agent_logp) ** 2).mean()

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(agent.network.parameters(), max_norm=1.0)
        optimizer.step()

        mean_reward = float(rewards_t.mean().item())
        rec = {
            "step": step,
            "loss": float(loss.item()),
            "mean_reward": mean_reward,
            "validity": comp_stats["validity"],
            "mean_qed": comp_stats["mean_qed"],
            "pains_pass_rate": comp_stats["pains_pass_rate"],
            "mean_tanimoto": comp_stats["mean_tanimoto"],
            "elapsed_s": round(time.time() - t0, 2),
        }
        history.append(rec)

        # Single-line progress so a background runner can tail-grep without
        # parsing JSON every line. Flushed so the orchestrator sees it now.
        print(
            f"[rl] step={step}/{args.steps} loss={rec['loss']:.4f} "
            f"reward={mean_reward:.4f} valid={comp_stats['validity']:.2f} "
            f"qed={comp_stats['mean_qed']:.3f} pains={comp_stats['pains_pass_rate']:.2f} "
            f"tani={comp_stats['mean_tanimoto']:.3f} t={rec['elapsed_s']}s",
            flush=True,
        )

        # 3e. Intermediate checkpoint + log dump.
        if args.checkpoint_every and step % args.checkpoint_every == 0:
            _save_checkpoint(agent, out_path)
            log_path.write_text(json.dumps({
                "target": target,
                "args": vars(args),
                "history": history,
            }, indent=2), encoding="utf-8")

        if mean_reward > best_mean_reward:
            best_mean_reward = mean_reward

    # 4. Final checkpoint + log.
    _save_checkpoint(agent, out_path)
    log_path.write_text(json.dumps({
        "target": target,
        "args": vars(args),
        "history": history,
        "best_mean_reward": best_mean_reward,
    }, indent=2), encoding="utf-8")
    print(f"[rl] DONE — saved {out_path}, best mean reward {best_mean_reward:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
