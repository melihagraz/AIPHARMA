"""Forensic inspection of Reinvent_App's .prior files.

Verifies empirically whether the prior weights are TRAINED or RANDOM:

1. Load the .prior pickle (REINVENT 4 format) using the sampler's alias setup.
2. Walk the network parameters; compute mean, std, min, max, |w|<0.001 ratio.
3. Untrained / random weights typically:
     - Initialized from a uniform/normal distribution
     - std ≈ 0.05–0.1 (Xavier/Kaiming init range)
     - No bias toward near-zero values
   Trained weights typically:
     - std much smaller for most layers (0.01–0.05)
     - Long tails (some weights become very large or near-zero)
     - High concentration near 0 (sparsity emerges through learning)

4. Sample 50 SMILES via the existing sampler. A trained prior should give
   ~80–99 % valid SMILES; a random prior gives ~0–10 %.

Run:
    .venv\\Scripts\\python.exe scripts\\inspect_prior.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Reinvent_App"))

import torch  # noqa: E402

from app.sampler import (  # noqa: E402
    _install_pickle_aliases,
    resolve_prior_path,
    get_default_sampler,
)


def inspect_state(name: str, prior_path: Path) -> None:
    print(f"\n{'='*70}")
    print(f"INSPECTING: {name}")
    print(f"PATH:       {prior_path}")
    print(f"SIZE:       {prior_path.stat().st_size:,} bytes "
          f"({prior_path.stat().st_size / 1024 / 1024:.1f} MB)")
    print('='*70)

    _install_pickle_aliases()
    try:
        obj = torch.load(prior_path, map_location="cpu", weights_only=False)
    except Exception as exc:
        print(f"  LOAD FAILED: {exc}")
        return

    print(f"  Object type: {type(obj).__name__}")
    # REINVENT 4 prior is a dict whose 'network' entry is a state_dict
    # (OrderedDict[str, Tensor]) — not an nn.Module. Iterate as a state_dict.
    network = None
    metadata = None
    network_params = None
    vocabulary = None
    if isinstance(obj, dict):
        print(f"  Dict keys: {list(obj.keys())}")
        if obj.get("version"):
            print(f"  REINVENT version: {obj.get('version')!r}")
        if obj.get("model_type"):
            print(f"  Model type:       {obj.get('model_type')!r}")
        metadata = obj.get("metadata")
        network_params = obj.get("network_params")
        vocabulary = obj.get("vocabulary")
        network = obj.get("network")
        if metadata:
            print(f"  Metadata: {metadata!r}")
        if network_params:
            print(f"  Network params: {network_params!r}")
        if vocabulary is not None:
            try:
                vocab_len = len(vocabulary)
            except TypeError:
                vocab_len = "?"
            print(f"  Vocabulary size:  {vocab_len}")
    elif hasattr(obj, "named_parameters"):
        network = {k: v for k, v in obj.named_parameters()}

    if network is None:
        print("  Could not locate network weights. Aborting.")
        return

    # network here is an OrderedDict[str, Tensor] (state_dict)
    items = list(network.items())
    total = sum(int(v.numel()) for _, v in items if hasattr(v, "numel"))
    print(f"  Total parameters: {total:,}")
    print(f"  Tensor count:    {len(items)}")
    params = items  # rename for downstream code

    print("\n  Per-tensor weight statistics:")
    print(f"  {'tensor':<40} {'shape':<22} {'mean':>10} {'std':>10} {'|w|<.001 %':>12}")
    print(f"  {'-'*40} {'-'*22} {'-'*10} {'-'*10} {'-'*12}")
    untrained_signals = 0
    trained_signals = 0
    for pname, p in params:
        if not hasattr(p, "numel") or p.numel() == 0:
            continue
        w = p.detach().float() if hasattr(p, "detach") else torch.as_tensor(p, dtype=torch.float32)
        mean = w.mean().item()
        std = w.std().item()
        near_zero_pct = (w.abs() < 0.001).float().mean().item() * 100
        shape = "x".join(str(s) for s in tuple(w.shape))
        print(f"  {pname[:40]:<40} {shape:<22} {mean:>10.4f} {std:>10.4f} {near_zero_pct:>11.1f}%")
        # Cheap heuristic — Xavier init for LSTM weights gives std ~ 1/sqrt(hidden)
        # so for hidden=512 we'd see std ~ 0.044 across the board.
        # Trained LSTMs typically show LARGER std (1-1.5) in input-to-hidden weights
        # because gradient descent stretches the distribution.
        if "weight" in pname and "norm" not in pname.lower() and std > 0:
            if 0.03 <= std <= 0.06 and near_zero_pct < 5:
                untrained_signals += 1
            elif std > 0.1 or near_zero_pct > 8:
                trained_signals += 1

    print(f"\n  Heuristic verdict:")
    print(f"    Layers with untrained-init signature: {untrained_signals}")
    print(f"    Layers with trained signature:        {trained_signals}")
    if trained_signals > untrained_signals:
        print(f"    ==> Looks TRAINED.")
    elif untrained_signals > trained_signals * 2:
        print(f"    ==> Looks UNTRAINED / random.")
    else:
        print(f"    ==> Mixed signature; can't tell from weights alone.")


def sampling_test(prior_name: str = "reinvent", n: int = 50) -> None:
    print(f"\n{'='*70}")
    print(f"BEHAVIORAL TEST: sample {n} SMILES from prior '{prior_name}'")
    print('='*70)

    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception as exc:
        print(f"  RDKit not importable: {exc}")
        return

    try:
        sampler = get_default_sampler(prior_name=prior_name)
    except Exception as exc:
        print(f"  Sampler init failed: {exc}")
        return

    try:
        smiles_list = sampler.sample(batch_size=n, seed=42)
    except Exception as exc:
        print(f"  Sampling failed: {exc}")
        return

    valid = 0
    examples = []
    for s in smiles_list:
        smi = s if isinstance(s, str) else s.get("smiles", "")
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is not None and mol.GetNumAtoms() > 0:
            valid += 1
            if len(examples) < 5:
                examples.append(smi)

    pct = 100.0 * valid / max(1, len(smiles_list))
    print(f"  Generated:        {len(smiles_list)}")
    print(f"  Valid (RDKit):    {valid}  ({pct:.1f}%)")
    print(f"  Verdict:")
    if pct >= 80:
        print(f"    [OK] TRAINED — random priors give ~0–10%, you got {pct:.1f}%.")
    elif pct >= 40:
        print(f"    [!] PARTIALLY trained or unusual setup.")
    else:
        print(f"    [X] UNTRAINED — only {pct:.1f}% valid is consistent with random weights.")
    if examples:
        print(f"\n  Sample valid SMILES (first 5):")
        for s in examples:
            print(f"    {s}")


def main() -> int:
    base = resolve_prior_path(None)  # default: prior/reinvent.prior
    inspect_state("BASE REINVENT prior", base)

    ft_dir = ROOT / "Reinvent_App" / "prior" / "finetuned"
    for ft in sorted(ft_dir.glob("*.prior")):
        inspect_state(f"FINE-TUNED: {ft.name}", ft)

    sampling_test(prior_name="reinvent", n=50)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
