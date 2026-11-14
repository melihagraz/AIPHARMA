"""RL convergence study — does the REINFORCE loop actually learn?

This test answers the "is the reward signal real" question by running the
RL script for N steps and asserting that:
    * the reward distribution moves UP over training (last quartile mean
      strictly greater than first quartile mean)
    * validity stays above a floor (training shouldn't trash basic SMILES
      grammar — if it does, we've over-rewarded a degenerate mode)
    * the trained checkpoint loads back through the sampler and produces
      valid molecules (the file isn't corrupt)

It's marked ``@pytest.mark.slow`` because a meaningful run is ~100 steps
× batch=32 on CPU — roughly 5-10 minutes wall clock. The default ``pytest``
invocation skips it; opt in with ``pytest -m slow`` or set ``RL_CONV_STEPS``
in the env to override the step count for a faster (less reliable) check.

What this test does NOT prove:
    * That the prior generalises to held-out targets
    * That the molecules satisfy ANY real bioactivity criterion (Tanimoto-
      to-actives is a structural proxy, not a binding assay)
    * That hyperparameter choices (sigma, lr, weights) are optimal
The point is purely "does the training math move the reward needle". A
more rigorous evaluation lives in /api/reinvent/evaluate-prior (which
samples a large batch and reports the population-level distribution).
"""

from __future__ import annotations

import importlib.util
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

import pytest


pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def rl_script_path() -> Path:
    """Resolve the RL script path once per session."""
    return (
        Path(__file__).resolve().parents[2]
        / "Reinvent_App"
        / "scripts"
        / "rl_finetune.py"
    )


@pytest.fixture(scope="module")
def reinvent_prior_path() -> Path:
    """Skip the whole convergence module when the prior file is missing —
    it's a 23 MB artefact that may not ship with the repo on all
    machines, and there's nothing meaningful to test without it."""
    p = (
        Path(__file__).resolve().parents[2]
        / "Reinvent_App"
        / "prior"
        / "reinvent.prior"
    )
    if not p.is_file():
        pytest.skip(f"Reinvent prior file missing at {p} — can't run RL convergence study.")
    return p


def _run_rl(
    script: Path,
    out_dir: Path,
    target: str,
    steps: int,
    batch_size: int,
    timeout: int,
) -> tuple[Path, list[dict]]:
    """Spawn the RL script as a subprocess and return (checkpoint_path, history).

    Why subprocess instead of importing + calling main():
      The script's top-level code mutates sys.modules (pickle aliases)
      and sys.path. Running it inline would leak that state into the
      rest of the test session. A subprocess gives clean isolation —
      same pattern the production backend uses to run RL jobs.
    """
    out_path = out_dir / "agent.prior"
    log_path = out_dir / "log.json"
    cmd = [
        sys.executable, str(script),
        target,
        "--steps", str(steps),
        "--batch-size", str(batch_size),
        # Keep weights all-on so the test covers the full reward stack.
        "--weight-qed", "1",
        "--weight-pains", "1",
        "--weight-tani", "1",
        "--device", "cpu",
        "--checkpoint-every", str(max(1, steps // 4)),
        "--out-path", str(out_path),
        "--log-path", str(log_path),
        # Smaller max_len for snappier sampling — we don't need 256-char
        # SMILES to validate the training math.
        "--max-len", "120",
    ]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=str(script.parent.parent),
    )
    assert proc.returncode == 0, (
        f"RL script exited rc={proc.returncode}\n"
        f"--- stdout (last 2000) ---\n{proc.stdout[-2000:]}\n"
        f"--- stderr (last 2000) ---\n{proc.stderr[-2000:]}"
    )
    assert out_path.is_file(), f"Checkpoint missing at {out_path}"
    assert log_path.is_file(), f"Log missing at {log_path}"
    log = json.loads(log_path.read_text(encoding="utf-8"))
    history = log.get("history") or []
    assert history, "Empty training history — script wrote no per-step records."
    return out_path, history


def _quartile_means(values: list[float]) -> tuple[float, float]:
    """Return (mean of first quartile, mean of last quartile).

    Used to check 'reward went up across training' without assuming
    a monotonic curve — the trajectory is noisy and a strict 'every
    step is better than the previous one' assertion would flake.
    """
    n = len(values)
    if n < 4:
        return statistics.fmean(values), statistics.fmean(values)
    q = max(1, n // 4)
    first = values[:q]
    last = values[-q:]
    return statistics.fmean(first), statistics.fmean(last)


def test_reward_distribution_moves_upward(
    rl_script_path, reinvent_prior_path, tmp_path
):
    """Last-quartile mean reward must be strictly greater than the first-
    quartile mean. With 60 steps × batch 32 the signal-to-noise is high
    enough to make this a reliable assertion on CPU."""
    steps = int(os.environ.get("RL_CONV_STEPS", "60"))
    batch = int(os.environ.get("RL_CONV_BATCH", "32"))
    # 60 steps × ~6 s/step on CPU = ~6 min; pad the timeout for slower
    # CI runners. Override locally via RL_CONV_TIMEOUT for faster boxes.
    timeout = int(os.environ.get("RL_CONV_TIMEOUT", "1500"))
    t0 = time.time()
    out_path, history = _run_rl(
        rl_script_path, tmp_path,
        target="JAK2",  # base prior, no actives needed → neutral 0.5 Tanimoto.
        steps=steps, batch_size=batch, timeout=timeout,
    )
    elapsed = time.time() - t0
    rewards = [float(r["mean_reward"]) for r in history if "mean_reward" in r]
    assert len(rewards) >= steps - 1, f"Expected ≥{steps - 1} reward records, got {len(rewards)}"

    q1_mean, q4_mean = _quartile_means(rewards)
    print(
        f"[rl-conv] steps={steps} batch={batch} elapsed={elapsed:.0f}s "
        f"q1={q1_mean:.4f} q4={q4_mean:.4f} delta={q4_mean - q1_mean:+.4f}",
        flush=True,
    )
    # The minimum meaningful improvement is ~1% — any smaller and we
    # can't distinguish a real signal from sampling noise. If you
    # tighten this floor make sure the step count rises too.
    assert q4_mean > q1_mean + 0.01, (
        f"Reward distribution did not move upward: q1={q1_mean:.4f} q4={q4_mean:.4f}. "
        f"Training may be misconfigured (sigma too small? lr too large?), or the "
        f"reward components may be saturated already. Inspect {out_path.with_suffix('.json')}."
    )


def test_validity_floor_preserved(rl_script_path, reinvent_prior_path, tmp_path):
    """Validity rate must stay above 0.5 throughout training.

    A common RL failure mode is "policy collapses to a degenerate token
    that maximises reward but isn't valid SMILES" — guarding against
    that is the whole reason validity is also in the reward. We assert
    the floor on every step rather than just the final step so an
    intermediate collapse fails fast."""
    steps = int(os.environ.get("RL_CONV_STEPS", "60"))
    batch = int(os.environ.get("RL_CONV_BATCH", "32"))
    timeout = int(os.environ.get("RL_CONV_TIMEOUT", "1500"))
    _, history = _run_rl(
        rl_script_path, tmp_path,
        target="JAK2",
        steps=steps, batch_size=batch, timeout=timeout,
    )
    validity = [float(r["validity"]) for r in history if "validity" in r]
    assert validity, "No validity records in training history."
    bad_steps = [(i, v) for i, v in enumerate(validity, start=1) if v < 0.5]
    assert not bad_steps, (
        f"Validity dropped below 0.5 at {len(bad_steps)} steps. "
        f"First few: {bad_steps[:5]}. The agent may have collapsed onto a "
        f"degenerate token sequence."
    )


def test_trained_checkpoint_round_trips_through_sampler(
    rl_script_path, reinvent_prior_path, tmp_path
):
    """The whole point of training is producing a usable checkpoint.
    After ~30 steps (less than the upward-move test) we should be
    able to load the saved .prior back through the sampler and draw
    valid SMILES from it. This catches serialisation regressions
    (a checkpoint that 'looks fine' on disk but can't be re-loaded
    is worse than a clearly-broken one)."""
    steps = max(20, int(os.environ.get("RL_CONV_STEPS_LOAD", "20")))
    batch = int(os.environ.get("RL_CONV_BATCH", "16"))
    timeout = int(os.environ.get("RL_CONV_TIMEOUT", "600"))
    out_path, _ = _run_rl(
        rl_script_path, tmp_path,
        target="JAK2",
        steps=steps, batch_size=batch, timeout=timeout,
    )

    # Load the script as a module (same trick as test_rl_script.py) so
    # we can call its install_pickle_aliases without polluting the
    # test session's sys.modules.
    spec = importlib.util.spec_from_file_location("rl_under_test_load", rl_script_path)
    rl_mod = importlib.util.module_from_spec(spec)
    sys.modules["rl_under_test_load"] = rl_mod
    try:
        spec.loader.exec_module(rl_mod)
        rl_mod._install_pickle_aliases()

        import torch
        from reinvent_core.model import Model
        from rdkit import Chem
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")

        ckpt = torch.load(out_path, map_location="cpu", weights_only=False)
        model = Model.create_from_dict(ckpt, mode="inference", device=torch.device("cpu"))
        _, smiles, _ = model.sample(batch_size=20)
        n_valid = sum(1 for s in smiles if Chem.MolFromSmiles(s) is not None)
        # Most samples should be valid; if a freshly-trained checkpoint
        # produces <50 % valid, something is wrong with the save/load
        # path or the training cratered the model.
        assert n_valid >= 10, (
            f"Reloaded checkpoint produced only {n_valid}/20 valid SMILES — "
            f"checkpoint may be corrupt or training catastrophically failed."
        )
    finally:
        sys.modules.pop("rl_under_test_load", None)
