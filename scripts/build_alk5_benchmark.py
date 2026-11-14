#!/usr/bin/env python
"""Build a property-matched human ALK5 active/control benchmark from ChEMBL JSON.

The input files are ChEMBL activity API responses. Potent actives and weak
controls are canonicalized, deduplicated by molecule, restricted to exact
measurements, and greedily matched on basic physicochemical properties.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

REPO = Path(__file__).resolve().parent.parent


def describe(smiles: str) -> dict | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    canonical = Chem.MolToSmiles(mol, canonical=True)
    core = MurckoScaffold.GetScaffoldForMol(mol)
    return {
        "smiles": canonical,
        "mw": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rot": Lipinski.NumRotatableBonds(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "scaffold": Chem.MolToSmiles(core, canonical=True) if core.GetNumAtoms() else "ACYCLIC",
    }


def load(path: Path, active: bool) -> list[dict]:
    body = json.loads(path.read_text(encoding="utf-8"))
    best: dict[str, dict] = {}
    for row in body.get("activities", []):
        if row.get("standard_relation") != "=":
            continue
        pchembl = row.get("pchembl_value")
        smi = row.get("canonical_smiles")
        if pchembl is None or not smi:
            continue
        desc = describe(smi)
        if desc is None:
            continue
        desc.update({
            "molecule_chembl_id": row.get("molecule_chembl_id"),
            "pchembl_value": float(pchembl),
            "standard_type": row.get("standard_type"),
            "label": int(active),
        })
        key = desc["smiles"]
        previous = best.get(key)
        if previous is None or (
            active and desc["pchembl_value"] > previous["pchembl_value"]
        ) or (
            not active and desc["pchembl_value"] < previous["pchembl_value"]
        ):
            best[key] = desc
    return list(best.values())


def distance(a: dict, b: dict) -> float:
    scales = {"mw": 50.0, "logp": 1.0, "hbd": 1.0, "hba": 2.0, "rot": 2.0, "tpsa": 30.0}
    return math.sqrt(sum(((a[key] - b[key]) / scale) ** 2 for key, scale in scales.items()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actives-json", type=Path, required=True)
    parser.add_argument("--controls-json", type=Path, required=True)
    parser.add_argument("--max-pairs", type=int, default=100)
    parser.add_argument("--allow-scaffold-overlap", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=REPO / "reports" / "alk5_active_control_benchmark.csv"
    )
    args = parser.parse_args()

    actives = load(args.actives_json, True)
    controls = load(args.controls_json, False)
    active_scaffolds = {row["scaffold"] for row in actives}
    if not args.allow_scaffold_overlap:
        controls = [row for row in controls if row["scaffold"] not in active_scaffolds]

    available = list(actives)
    pairs: list[tuple[dict, dict]] = []
    for control in sorted(controls, key=lambda row: (row["mw"], row["logp"])):
        if not available or len(pairs) >= args.max_pairs:
            break
        active = min(available, key=lambda row: distance(row, control))
        available.remove(active)
        pairs.append((active, control))

    fields = [
        "pair_id", "label", "molecule_chembl_id", "smiles", "pchembl_value",
        "standard_type", "mw", "logp", "hbd", "hba", "rot", "tpsa", "scaffold", "score",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for pair_id, (active, control) in enumerate(pairs, 1):
            for row in (active, control):
                writer.writerow({key: row.get(key, "") for key in fields} | {"pair_id": pair_id})

    active_rows = [pair[0] for pair in pairs]
    control_rows = [pair[1] for pair in pairs]
    diagnostics = {
        "target": "human ALK5/TGFBR1",
        "target_chembl_id": "CHEMBL4439",
        "pairs": len(pairs),
        "active_threshold": "pChEMBL >= 7, exact relation",
        "control_threshold": "pChEMBL <= 5, exact relation",
        "exact_scaffold_overlap_removed": not args.allow_scaffold_overlap,
        "mean_property_delta": {
            key: round(sum(abs(a[key] - c[key]) for a, c in pairs) / len(pairs), 3)
            if pairs else None
            for key in ("mw", "logp", "hbd", "hba", "rot", "tpsa")
        },
        "active_scaffolds": len({row["scaffold"] for row in active_rows}),
        "control_scaffolds": len({row["scaffold"] for row in control_rows}),
    }
    meta = args.output.with_suffix(".metadata.json")
    meta.write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
    print(json.dumps(diagnostics, indent=2))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
