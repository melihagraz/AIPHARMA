#!/usr/bin/env python
"""Validate REINVENT generation stability across multiple random seeds.

This is a generation-level validation, not evidence of biological activity.
It reports validity, uniqueness, PAH drug-likeness pass rate, Bemis-Murcko
scaffold diversity, and cross-seed overlap in a machine-readable JSON file.

Example:
    python scripts/validate_generation_seeds.py --seeds 1 2 3 4 5 6 7 8 9 10
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Reinvent_App"))
sys.path.insert(0, str(REPO / "scripts"))

from app.sampler import ReinventSampler  # noqa: E402
from run_pah_pipeline import FILT, drug_like  # noqa: E402


def canonicalize(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, canonical=True) if mol is not None else None


def scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    core = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(core, canonical=True) if core.GetNumAtoms() else "ACYCLIC"


def summarize_seed(seed: int, requested: int) -> tuple[dict, set[str]]:
    sampled = ReinventSampler().sample(batch_size=requested, seed=seed)
    raw = [item.smiles for item in sampled if item.smiles]
    canonical = [value for smi in raw if (value := canonicalize(smi))]
    unique = set(canonical)
    passing = {smi for smi in unique if (props := drug_like(smi)) and props["pass"]}
    scaffolds = {scaffold(smi) for smi in passing}
    scaffolds.discard("")
    return {
        "seed": seed,
        "requested": requested,
        "returned": len(raw),
        "valid": len(canonical),
        "unique": len(unique),
        "drug_like": len(passing),
        "scaffolds": len(scaffolds),
        "validity_rate": round(len(canonical) / requested, 4),
        "uniqueness_rate": round(len(unique) / len(canonical), 4) if canonical else 0.0,
        "drug_like_rate": round(len(passing) / len(unique), 4) if unique else 0.0,
        "scaffold_to_drug_like_ratio": round(len(scaffolds) / len(passing), 4) if passing else 0.0,
    }, passing


def mean_sd(rows: list[dict], key: str) -> dict:
    values = [float(row[key]) for row in rows]
    return {
        "mean": round(statistics.mean(values), 4),
        "sd": round(statistics.stdev(values), 4) if len(values) > 1 else 0.0,
        "min": round(min(values), 4),
        "max": round(max(values), 4),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(1, 11)))
    parser.add_argument("--generate", type=int, default=100, help="Molecules requested per seed")
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO / "reports" / "generation_10seed_validation.json",
    )
    args = parser.parse_args()
    if len(set(args.seeds)) != len(args.seeds):
        parser.error("Seed values must be unique")
    if args.generate < 1:
        parser.error("--generate must be positive")

    rows: list[dict] = []
    passing_by_seed: dict[int, set[str]] = {}
    for seed in args.seeds:
        row, passing = summarize_seed(seed, args.generate)
        rows.append(row)
        passing_by_seed[seed] = passing
        print(
            f"seed={seed:>3} valid={row['valid']:>4}/{args.generate} "
            f"unique={row['unique']:>4} drug_like={row['drug_like']:>4} "
            f"scaffolds={row['scaffolds']:>4}",
            flush=True,
        )

    all_passing = set().union(*passing_by_seed.values()) if rows else set()
    all_scaffolds = {scaffold(smi) for smi in all_passing}
    all_scaffolds.discard("")
    recurring = {
        smi
        for smi in all_passing
        if sum(smi in values for values in passing_by_seed.values()) > 1
    }
    metrics = [
        "validity_rate",
        "uniqueness_rate",
        "drug_like_rate",
        "scaffold_to_drug_like_ratio",
    ]
    result = {
        "validation_type": "multi_seed_generation",
        "target_claim": "generation robustness only; not biological activity",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "model": "REINVENT prior",
        "seeds": args.seeds,
        "requested_per_seed": args.generate,
        "filter_thresholds": FILT,
        "per_seed": rows,
        "aggregate": {
            key: mean_sd(rows, key) for key in metrics
        } | {
            "total_requested": args.generate * len(rows),
            "unique_drug_like_across_seeds": len(all_passing),
            "unique_scaffolds_across_seeds": len(all_scaffolds),
            "drug_like_repeated_across_seeds": len(recurring),
            "cross_seed_repeat_rate": round(len(recurring) / len(all_passing), 4)
            if all_passing else 0.0,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
