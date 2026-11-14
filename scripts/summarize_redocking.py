#!/usr/bin/env python
"""Aggregate redocking JSON files without overstating mixed-target results."""
from __future__ import annotations

import argparse
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def wilson(successes: int, total: int, z: float = 1.96) -> list[float]:
    if total == 0:
        return [0.0, 0.0]
    p = successes / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    margin = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / den
    return [round(max(0.0, center - margin), 4), round(min(1.0, center + margin), 4)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reports-dir", type=Path, default=REPO / "reports")
    parser.add_argument("--target", default="ALK5")
    parser.add_argument("--target-pdb", nargs="+", default=["1PY5", "1VJY"])
    parser.add_argument("--threshold", type=float, default=2.0)
    parser.add_argument(
        "--output", type=Path, default=REPO / "reports" / "redocking_validation_summary.json"
    )
    args = parser.parse_args()

    rows = []
    for path in sorted(args.reports_dir.glob("redock_*.json")):
        if path.name == args.output.name:
            continue
        row = json.loads(path.read_text(encoding="utf-8"))
        rmsd = row.get("rmsd")
        if rmsd is None:
            continue
        row["source"] = path.name
        row["passes_threshold"] = float(rmsd) <= args.threshold
        row["target_specific"] = str(row.get("pdb", "")).upper() in {
            value.upper() for value in args.target_pdb
        }
        rows.append(row)

    def summary(items: list[dict]) -> dict:
        rmsds = [float(row["rmsd"]) for row in items]
        passed = sum(bool(row["passes_threshold"]) for row in items)
        return {
            "n": len(items),
            "passed": passed,
            "success_rate": round(passed / len(items), 4) if items else 0.0,
            "success_rate_wilson_95ci": wilson(passed, len(items)),
            "median_rmsd_angstrom": round(statistics.median(rmsds), 3) if rmsds else None,
            "mean_rmsd_angstrom": round(statistics.mean(rmsds), 3) if rmsds else None,
        }

    target_rows = [row for row in rows if row["target_specific"]]
    result = {
        "validation_type": "redocking",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "success_threshold_angstrom": args.threshold,
        "target": args.target,
        "target_pdb_ids": [value.upper() for value in args.target_pdb],
        "target_specific": summary(target_rows),
        "mixed_target_protocol_check": summary(rows),
        "interpretation_guardrail": (
            "Mixed-target success estimates validate the docking protocol broadly and must not be "
            f"reported as {args.target}-specific performance. A single target structure is insufficient "
            "for a general target-specific claim."
        ),
        "runs": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"target_specific": result["target_specific"],
                      "mixed_target": result["mixed_target_protocol_check"]}, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
