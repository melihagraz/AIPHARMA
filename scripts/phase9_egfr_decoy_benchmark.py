#!/usr/bin/env python3
"""Prospective, frozen EGFR DUD-E active/decoy Vina pilot.

Runs on the *newly crystallographically prepared 1M17 receptor*, not on the
older dockstring EGFR bundled receptor. Requires a completed Phase 9 redocking
report with pose RMSD <2 A; otherwise this program fails closed. DUD-E decoys
are presumed decoys (NOT experimentally tested EGFR inactives).
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
from urllib.request import urlopen, Request

URLS = {
    "active": "https://dude.docking.org/targets/egfr/actives_final.ism",
    "decoy": "https://dude.docking.org/targets/egfr/decoys_final.ism",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def acquire(path: Path, url: str, offline=False):
    if path.exists() and path.stat().st_size > 50:
        return
    if offline:
        raise FileNotFoundError(f"Missing DUD-E source: {path}; download from {url}")
    req = Request(url, headers={"User-Agent": "AIPHARMA-phase9-benchmark/1.0"})
    with urlopen(req, timeout=35) as response:
        content = response.read(25_000_000)
    if len(content) < 50 or b"<html" in content[:1000].lower():
        raise ValueError(f"Download is not an ISM molecular file: {url}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def parse_ism(path: Path):
    from rdkit import Chem
    good = {}
    rejected = 0
    for line in path.read_text(errors="replace").splitlines():
        tok = line.split()
        if not tok or tok[0].startswith("#"):
            continue
        smi = tok[0]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rejected += 1
            continue
        canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        good.setdefault(canon, {"smiles": canon, "source_id": tok[1] if len(tok) > 1 else ""})
    return list(good.values()), rejected


def auc_pairwise(pos, neg):
    # Compares binding-ranking score = -Vina affinity, where higher is better.
    if not pos or not neg:
        return None
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


def stratified_ci(pos, neg, *, reps=1000, seed=20261109):
    if not pos or not neg:
        return None
    rng = random.Random(seed)
    vals = []
    for _ in range(reps):
        ps = rng.choices(pos, k=len(pos))
        ns = rng.choices(neg, k=len(neg))
        vals.append(auc_pairwise(ps, ns))
    vals.sort()
    return [round(vals[int(0.025 * reps)], 4), round(vals[min(reps - 1, int(0.975 * reps))], 4)]


def build_sample(actives, decoys, *, n, seed):
    active_ids = {x["smiles"] for x in actives}
    decoys = [x for x in decoys if x["smiles"] not in active_ids]
    if len(actives) < n or len(decoys) < n:
        raise ValueError(f"Insufficient unique non-overlapping inputs for n={n}")
    rng = random.Random(seed)
    a = rng.sample(sorted(actives, key=lambda x: x["smiles"]), n)
    d = rng.sample(sorted(decoys, key=lambda x: x["smiles"]), n)
    rows = [{**x, "class": "active"} for x in a] + [{**x, "class": "decoy"} for x in d]
    rng.shuffle(rows)
    return rows


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--n-per-class", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261109)
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--offline", action="store_true")
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--out", type=Path, default=Path("reports/phase9/egfr_benchmark"))
    a = p.parse_args(argv)
    if not (4 <= a.n_per_class <= 1000):
        p.error("n-per-class must be 4..1000")
    root = a.root.resolve()
    redock = root / "reports/phase9/redocking/1M17/redock_report.json"
    if not redock.exists():
        p.error("Run scripts/phase9_redock.py --target EGFR first")
    rdr = json.loads(redock.read_text())
    if not rdr.get("passed") or not rdr.get("pose_reproduced_within_2A"):
        p.error("EGFR crystal redocking RMSD must be measured and <2 A before benchmarking")
    receptor = redock.parent / "prepared_receptor.pdbqt"
    config = redock.parent / "vina_grid.txt"
    hashes = rdr.get("files_sha256", {})
    if not receptor.exists() or sha(receptor) != hashes.get(receptor.name):
        p.error("EGFR receptor missing or modified since validated redocking")
    if not config.exists() or sha(config) != hashes.get(config.name):
        p.error("EGFR grid missing or modified since validated redocking")
    source_dir = a.out / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    input_paths = {c: source_dir / f"{c}.ism" for c in URLS}
    for klass, file in input_paths.items():
        acquire(file, URLS[klass], offline=a.offline)
    actives, failed_a = parse_ism(input_paths["active"])
    decoys, failed_d = parse_ism(input_paths["decoy"])
    panel = build_sample(actives, decoys, n=a.n_per_class, seed=a.seed)
    a.out.mkdir(parents=True, exist_ok=True)
    # Selection is frozen and saved BEFORE scoring. Never redraw failed molecules.
    with (a.out / "frozen_panel.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["smiles", "source_id", "class"])
        writer.writeheader(); writer.writerows(panel)
    from Docking_App.app.dock import dock_one
    box = tuple(rdr["center"] + rdr["size"])
    rows = []
    for i, item in enumerate(panel):
        outcome = dock_one(item["smiles"], "EGFR", exhaustiveness=a.exhaustiveness,
                           num_modes=3, cpu=1, timeout=360,
                           pocket_box=box, custom_receptor_pdbqt=receptor)
        rows.append({**item, "score_kcal_mol": outcome.affinity,
                     "error": outcome.error, "num_modes": outcome.num_modes})
        print(f"{i + 1}/{len(panel)} {item['class']} {outcome.affinity} {outcome.error or ''}", flush=True)
    with (a.out / "actual_vina_scores.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["smiles", "source_id", "class", "score_kcal_mol", "error", "num_modes"])
        writer.writeheader(); writer.writerows(rows)
    vals = {c: [-float(r["score_kcal_mol"]) for r in rows if r["class"] == c and r["score_kcal_mol"] is not None]
            for c in ("active", "decoy")}
    valid = {k: len(v) for k, v in vals.items()}
    enough = all(valid[k] / a.n_per_class >= .90 for k in valid)
    report = {"created_utc": datetime.now(timezone.utc).isoformat(), "target": "EGFR",
              "dataset": "DUD-E EGFR", "source_urls": URLS,
              "source_hashes": {c: sha(x) for c, x in input_paths.items()},
              "receptor_sha256": sha(receptor), "grid_sha256": sha(config),
              "source_rdkit_valid_counts": {"active": len(actives), "decoy": len(decoys)},
              "rejected_invalid_smiles": {"active": failed_a, "decoy": failed_d},
              "fixed_seed": a.seed, "requested_per_class": a.n_per_class,
              "valid_docking_counts": valid, "coverage_sufficient_for_primary_auc": enough,
              "roc_auc_complete_case": round(auc_pairwise(vals["active"], vals["decoy"]), 4) if enough else None,
              "stratified_bootstrap95_conditional": stratified_ci(vals["active"], vals["decoy"], seed=a.seed) if enough else None,
              "limitations": ["Decoys are computer-generated presumed negatives, not confirmed inactive by assay",
                              "Small pilot and possible source/chemotype bias; not prospective experiment",
                              "No activity estimates from Vina kcal/mol", "Score-missingness may be class-dependent",
                              "If sample n is small, AUC uncertainty and stability are limited"],
              "actual_docking_executed": True}
    (a.out / "analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"coverage": valid, "roc_auc": report["roc_auc_complete_case"],
                      "95_ci": report["stratified_bootstrap95_conditional"]}, indent=2))
    return 0 if enough else 3


if __name__ == "__main__":
    raise SystemExit(main())
