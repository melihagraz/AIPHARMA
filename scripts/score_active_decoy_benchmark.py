#!/usr/bin/env python
"""Score an ALK5 active/decoy docking benchmark from a CSV file.

Required columns: ``label`` (1 active, 0 decoy) and ``score``. Docking scores
are assumed to be lower-is-better unless ``--higher-is-better`` is supplied.
Optional ``smiles`` values are used to report active/decoy scaffold leakage.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Scoring.Scoring import CalcBEDROC
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parent.parent


def scaffold(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    core = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(core, canonical=True) if core.GetNumAtoms() else "ACYCLIC"


def ef_at_fraction(labels: list[int], ranks: list[float], fraction: float) -> float:
    n = max(1, math.ceil(len(labels) * fraction))
    order = sorted(range(len(labels)), key=lambda idx: ranks[idx], reverse=True)
    found = sum(labels[idx] for idx in order[:n])
    total = sum(labels)
    return (found / n) / (total / len(labels)) if total else 0.0


def metrics(labels: list[int], ranks: list[float]) -> dict:
    ordered = sorted(zip(ranks, labels), reverse=True)
    return {
        "roc_auc": float(roc_auc_score(labels, ranks)),
        "pr_auc": float(average_precision_score(labels, ranks)),
        "ef_1pct": float(ef_at_fraction(labels, ranks, 0.01)),
        "bedroc_alpha_20": float(CalcBEDROC([[score, label] for score, label in ordered], 1, 20.0)),
    }


def bootstrap(labels: list[int], ranks: list[float], repeats: int, seed: int) -> dict:
    rng = random.Random(seed)
    samples: dict[str, list[float]] = {key: [] for key in metrics(labels, ranks)}
    for _ in range(repeats):
        idx = [rng.randrange(len(labels)) for _ in labels]
        ys = [labels[i] for i in idx]
        if len(set(ys)) < 2:
            continue
        values = metrics(ys, [ranks[i] for i in idx])
        for key, value in values.items():
            samples[key].append(value)
    out = {}
    for key, values in samples.items():
        values.sort()
        if not values:
            out[key] = [None, None]
            continue
        lo = values[int(0.025 * (len(values) - 1))]
        hi = values[int(0.975 * (len(values) - 1))]
        out[key] = [round(lo, 4), round(hi, 4)]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--higher-is-better", action="store_true")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--output", type=Path, default=REPO / "reports" / "alk5_active_decoy_benchmark.json"
    )
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not {"label", "score"}.issubset(rows[0]):
        parser.error("CSV must contain label and score columns")
    total_input = len(rows)
    failed_rows = [row for row in rows if not str(row.get("score", "")).strip()]
    rows = [row for row in rows if str(row.get("score", "")).strip()]
    if not rows:
        parser.error("CSV contains no successfully scored rows")
    labels = [int(row["label"]) for row in rows]
    if set(labels) != {0, 1}:
        parser.error("label must contain both 0 (decoy) and 1 (active)")
    raw_scores = [float(row["score"]) for row in rows]
    ranks = raw_scores if args.higher_is_better else [-value for value in raw_scores]
    values = metrics(labels, ranks)

    active_scaffolds: set[str] = set()
    decoy_scaffolds: set[str] = set()
    invalid_smiles = 0
    for row, label in zip(rows, labels):
        if not row.get("smiles"):
            continue
        value = scaffold(row["smiles"])
        if value is None:
            invalid_smiles += 1
        elif label:
            active_scaffolds.add(value)
        else:
            decoy_scaffolds.add(value)
    overlap = active_scaffolds & decoy_scaffolds
    result = {
        "validation_type": "active_decoy_virtual_screening",
        "input": str(args.input.resolve()),
        "score_direction": "higher_is_better" if args.higher_is_better else "lower_is_better",
        "n": len(rows),
        "input_rows": total_input,
        "excluded_missing_score": len(failed_rows),
        "excluded_by_label": {
            "active": sum(int(row.get("label", 0)) == 1 for row in failed_rows),
            "decoy": sum(int(row.get("label", 0)) == 0 for row in failed_rows),
        },
        "actives": sum(labels),
        "decoys": len(labels) - sum(labels),
        "metrics": {key: round(value, 4) for key, value in values.items()},
        "bootstrap_95ci": bootstrap(labels, ranks, args.bootstrap, args.seed),
        "scaffold_leakage": {
            "active_scaffolds": len(active_scaffolds),
            "decoy_scaffolds": len(decoy_scaffolds),
            "overlapping_scaffolds": len(overlap),
            "overlap_fraction_of_active_scaffolds": round(
                len(overlap) / len(active_scaffolds), 4
            ) if active_scaffolds else None,
            "invalid_smiles": invalid_smiles,
        },
        "interpretation_guardrail": (
            "Enrichment validates ranking on this retrospective set only; it does not establish "
            "binding, inhibition, safety, or efficacy."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
