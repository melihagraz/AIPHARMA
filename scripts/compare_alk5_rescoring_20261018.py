"""Exploratory, fixed-pose scoring comparison. Never selects a production winner.

Run in the existing RDKit environment with isolated vina==1.2.7 on PYTHONPATH.
Inputs are immutable artifacts from the previous ensemble experiment.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from vina import Vina
import vina


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(labels, values):
    y, s = np.asarray(labels), np.asarray(values)
    auc = float(roc_auc_score(y, s))
    rng = np.random.default_rng(20261018)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    boot = []
    for _ in range(2000):
        ix = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        boot.append(roc_auc_score(y[ix], s[ix]))
    null = []
    for indices in itertools.combinations(range(len(y)), int(y.sum())):
        yp = np.zeros(len(y)); yp[list(indices)] = 1
        null.append(roc_auc_score(yp, s))
    return {"auc": auc, "stratified_molecule_bootstrap_95": np.quantile(boot, [.025, .975]).tolist(),
            "exact_label_permutation_p_one_sided": float(np.mean(np.asarray(null) >= auc - 1e-12)),
            "permutations": len(null)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--audit", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    source, audit, out = args.source, args.audit, args.output
    out.mkdir(parents=True, exist_ok=True)
    if (out / "metrics.json").exists():
        raise RuntimeError("Completed comparison is immutable; choose a new output directory")
    items = json.loads((source / "locked_items.json").read_text())
    assays = {r["assay_id"]: r for r in json.loads((audit / "assay_classification.json").read_text())}
    protocol = {"scope": "Post-hoc exploratory fixed-pose rescoring, not an independent test or Vinardo redocking benchmark",
                "vina_version": vina.__version__, "seed": 20261018,
                "pose_selection": "MODEL 1 for each of the two previously fixed Vina seeds; no new pose search or minimization",
                "methods": ["historical_Vina", "historical_Vina_per_heavy_atom", "Vina_1.2.7_fixed_pose", "Vinardo_fixed_pose"],
                "direction": "negative score and negative score/heavy_atoms: higher is ranked first",
                "aggregation": "Mean of two seed scores within each receptor, never pooled across receptors",
                "control_limit": "3 active + 3 weak molecules, not 36 independent controls; shared documents reduce independence further",
                "labels": "Previously exposed; endpoint classification based on assay descriptions, not original full protocols",
                "promotion_rule": "No production ranking method promoted from this small exposed comparison",
                "input_hashes": {"locked_items": digest(source / "locked_items.json"), "historical_scores": digest(source / "evaluation/scores.csv"), "assay_classification": digest(audit / "assay_classification.json")}}
    (out / "protocol.json").write_text(json.dumps(protocol, indent=2))
    records, inputs, control_audit = [], {}, []
    for item in items:
        mol = Chem.MolFromSmiles(item["smiles"])
        item["heavy_atoms"] = mol.GetNumHeavyAtoms()
        item["mw"] = Descriptors.MolWt(mol)
        item["logp"] = Crippen.MolLogP(mol)
        if item["kind"] != "candidate":
            descriptions = [assays[a] for a in item["assay_ids"]]
            assert all(a["endpoint"] == "enzyme_catalytic" and a["mutation"] == "no_explicit_mutant" for a in descriptions)
            control_audit.append({**item, "assay_review": descriptions})
    (out / "control_audit.json").write_text(json.dumps(control_audit, indent=2))
    for pdb in ["1PY5", "1VJY", "2WOT"]:
        box = json.loads((source / pdb / "box.json").read_text())
        receptor = source / pdb / "receptor.pdbqt"
        inputs[str(receptor)] = digest(receptor)
        for sf in ["vina", "vinardo"]:
            engine = Vina(sf_name=sf, cpu=1, seed=20261018, verbosity=0)
            engine.set_receptor(str(receptor))
            engine.compute_vina_maps(center=box["center"], box_size=box["size"])
            for item in items:
                for seed in [20261018, 20261019]:
                    pose = source / "jobs" / f"{pdb}_{item['id']}_{seed}" / "pose.pdbqt"
                    text = pose.read_text()
                    inputs[str(pose)] = digest(pose)
                    assert "MODEL 1" in text and "ENDMDL" in text
                    first = text.split("MODEL 1", 1)[1].split("ENDMDL", 1)[0].strip() + "\n"
                    dest = out / "poses" / f"{pdb}_{item['id']}_{seed}.pdbqt"
                    dest.parent.mkdir(exist_ok=True)
                    if not dest.exists(): dest.write_text(first)
                    engine.set_ligand_from_string(first)
                    energies = engine.score()
                    records.append({"pdb": pdb, "id": item["id"], "kind": item["kind"], "seed": seed,
                                    "scoring": sf, "score": float(energies[0]), "energy_terms": energies.tolist()})
            print(f"scored {pdb} {sf}", flush=True)
    (out / "raw_scores.json").write_text(json.dumps(records, indent=2))
    historical = list(csv.DictReader((source / "evaluation/scores.csv").open()))
    rows, metrics = [], {}
    for pdb in ["1PY5", "1VJY", "2WOT"]:
        for item in items:
            old = next(r for r in historical if r["pdb"] == pdb and r["id"] == item["id"])
            row = {k: item[k] for k in ["id", "kind", "heavy_atoms", "mw", "logp"]}
            row.update({"pdb": pdb, "historical_Vina": -float(old["mean_score"]),
                        "historical_Vina_per_heavy_atom": -float(old["mean_score"]) / item["heavy_atoms"]})
            for sf, name in [("vina", "Vina_1.2.7_fixed_pose"), ("vinardo", "Vinardo_fixed_pose")]:
                scores = [r["score"] for r in records if r["pdb"] == pdb and r["id"] == item["id"] and r["scoring"] == sf]
                assert len(scores) == 2
                row[name] = -float(np.mean(scores))
            rows.append(row)
        controls = [r for r in rows if r["pdb"] == pdb and r["kind"] != "candidate"]
        y = [int(r["kind"] == "active_control") for r in controls]
        metrics[pdb] = {name: evaluate(y, [r[name] for r in controls]) for name in protocol["methods"]}
        metrics[pdb]["size_baselines"] = {name: evaluate(y, [r[name] for r in controls]) for name in ["heavy_atoms", "mw", "logp"]}
        metrics[pdb]["correlations_controls_only"] = {name: float(spearmanr([r["historical_Vina"] for r in controls], [r[name] for r in controls]).statistic) for name in ["heavy_atoms", "mw", "logp"]}
    with (out / "comparison.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2))
    assert len(records) == 108 and len(rows) == 27 and len(control_audit) == 6
    assert all(digest(Path(path)) == expected for path, expected in inputs.items())
    (out / "verification.json").write_text(json.dumps({"passed": True, "fixed_pose_scores": len(records), "aggregated_rows": len(rows), "controls": len(control_audit), "input_files_unchanged": len(inputs), "input_hashes": inputs, "script_sha256": digest(Path(__file__))}, indent=2))
    print(json.dumps({p: {m: round(v["auc"], 3) for m, v in values.items() if "auc" in v} for p, values in metrics.items()}), flush=True)


if __name__ == "__main__":
    main()
