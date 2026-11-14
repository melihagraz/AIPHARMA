#!/usr/bin/env python
"""End-to-end PAH discovery run — reproducible, single script.

Pipeline:  PDB -> fpocket pocket -> REINVENT de novo -> RDKit drug-likeness
           (PAH profile) -> optional ADMET-AI Ames gate -> AutoDock Vina
           docking -> ranked candidates.

Designed to run inside the WSL conda env ``aipharma`` (fpocket + rdkit +
openbabel + autodock-vina + torch + admet_ai all present there). REINVENT is
imported from the repo's Reinvent_App; the PAH filter thresholds come from
backend/app/data/pah_filter_profile.py so code and manuscript stay in lock-step.

Usage (inside WSL):
    mamba run -n aipharma python scripts/run_pah_pipeline.py \
        --pdb 1PY5 --generate 150 --dock-top 8 --exhaustiveness 8
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REINVENT_DIR = REPO / "Reinvent_App"
RUN_ROOT = Path.home() / ".aipharma_run"

# PAH drug-likeness thresholds (base_filters subset — the RDKit-computable ones).
FILT = {"mw_max": 500.0, "logp_max": 5.0, "hbd_max": 5, "hba_max": 10,
        "rotatable_bonds_max": 10, "tpsa_max": 140.0, "qed_min": 0.50}


def log(msg: str) -> None:
    print(f"[pipeline] {msg}", flush=True)


def fetch_pdb(pdb: str, workdir: Path) -> Path:
    dest = workdir / f"{pdb}.pdb"
    if not dest.is_file() or dest.stat().st_size == 0:
        url = f"https://files.rcsb.org/download/{pdb.upper()}.pdb"
        log(f"fetching {pdb} from RCSB")
        urllib.request.urlretrieve(url, dest)  # noqa: S310
    return dest


# --------------------------------------------------------------------------- pocket
def detect_pocket(pdb_path: Path, workdir: Path) -> dict:
    log("fpocket cavity detection")
    subprocess.run(["fpocket", "-f", str(pdb_path)], cwd=workdir,
                   check=True, capture_output=True)
    stem = pdb_path.stem
    out = workdir / f"{stem}_out"
    info = (out / f"{stem}_info.txt").read_text()
    # best pocket by druggability
    best_num, best_drug = None, -1.0
    cur = None
    drug = {}
    for line in info.splitlines():
        s = line.strip()
        if s.startswith("Pocket ") and s.endswith(":"):
            cur = int(s.split()[1])
        elif cur is not None and "Druggability Score" in s:
            drug[cur] = float(s.split(":")[1])
    best_num = max(drug, key=drug.get)
    best_drug = drug[best_num]
    # center + box from that pocket's alpha spheres
    vert = (out / "pockets" / f"pocket{best_num}_vert.pqr").read_text()
    xs, ys, zs = [], [], []
    for line in vert.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            p = line.split()
            try:
                xs.append(float(p[-5])); ys.append(float(p[-4])); zs.append(float(p[-3]))
            except (ValueError, IndexError):
                continue
    def box(v):
        lo, hi = min(v), max(v)
        c = (lo + hi) / 2
        s = max(16.0, min(30.0, (hi - lo) + 8.0))
        return round(c, 3), round(s, 3)
    cx, sx = box(xs); cy, sy = box(ys); cz, sz = box(zs)
    log(f"  pocket {best_num}: druggability {best_drug:.3f}, "
        f"center ({cx},{cy},{cz}), box ({sx},{sy},{sz})")
    return {"pocket": best_num, "druggability": round(best_drug, 3),
            "center": [cx, cy, cz], "size": [sx, sy, sz]}


# --------------------------------------------------------------------------- receptor
# Catalytic metal ions: resname -> (AutoDock type, formal charge). Kept in the
# receptor (with the right charge) so metalloprotein pockets dock correctly.
METALS = {
    "ZN": ("Zn", 2.0), "MG": ("Mg", 2.0), "MN": ("Mn", 2.0), "FE": ("Fe", 2.0),
    "CA": ("Ca", 2.0), "NI": ("Ni", 2.0), "CU": ("Cu", 2.0), "CO": ("Co", 2.0),
    "NA": ("Na", 1.0), "K": ("K", 1.0),
}


def _fix_metal_charges(pdbqt: Path) -> int:
    types = {t: q for (t, q) in METALS.values()}
    n, out = 0, []
    for line in pdbqt.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line.split():
            adtype = line.split()[-1]
            if adtype in types:
                line, k = re.subn(
                    r"([+-]?\d+\.\d{3})(\s+" + re.escape(adtype) + r")\s*$",
                    f"{types[adtype]:+.3f}" + r"\2", line)
                n += k
        out.append(line)
    pdbqt.write_text("\n".join(out) + "\n")
    return n


def prep_receptor(pdb_path: Path, workdir: Path, pocket: dict) -> Path:
    import math
    log("preparing rigid receptor PDBQT (protein + active-site metals)")
    center = tuple(pocket["center"])
    clean = workdir / "receptor_clean.pdb"
    metals = []
    with pdb_path.open() as fh, clean.open("w") as out:
        for line in fh:
            if line.startswith(("ATOM", "TER")):
                out.write(line)
            elif line.startswith("HETATM") and line[17:20].strip().upper() in METALS:
                try:
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue
                if math.dist(xyz, center) <= 12.0:
                    out.write(line)
                    metals.append(line[17:20].strip().upper())
        out.write("END\n")
    pdbqt = workdir / "receptor.pdbqt"
    subprocess.run(["obabel", str(clean), "-O", str(pdbqt), "-xr", "-p", "7.4"],
                   check=True, capture_output=True)
    n_fixed = _fix_metal_charges(pdbqt)
    if metals:
        log(f"  kept {len(metals)} active-site metal ion(s): {', '.join(sorted(set(metals)))} "
            f"(formal charges assigned to {n_fixed})")
    return pdbqt


# --------------------------------------------------------------------------- generate
def generate(n: int, seed: int) -> list[tuple[str, float]]:
    log(f"REINVENT de novo generation (n={n}, seed={seed})")
    sys.path.insert(0, str(REINVENT_DIR))
    from app.sampler import ReinventSampler
    mols = ReinventSampler().sample(batch_size=n, seed=seed)
    return [(m.smiles, m.nll) for m in mols if m.smiles]


# --------------------------------------------------------------------------- filter
def drug_like(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, QED, Crippen, rdMolDescriptors, Lipinski
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    d = {
        "mw": Descriptors.MolWt(mol), "logp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "rot": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "hbd": Lipinski.NumHDonors(mol), "hba": Lipinski.NumHAcceptors(mol),
    }
    try:
        d["qed"] = QED.qed(mol)
    except Exception:
        d["qed"] = 0.0
    ok = (d["mw"] <= FILT["mw_max"] and d["logp"] <= FILT["logp_max"]
          and d["tpsa"] <= FILT["tpsa_max"] and d["rot"] <= FILT["rotatable_bonds_max"]
          and d["hbd"] <= FILT["hbd_max"] and d["hba"] <= FILT["hba_max"]
          and d["qed"] >= FILT["qed_min"])
    d["pass"] = ok
    return d


# --------------------------------------------------------------------------- dock
def dock(smiles: str, receptor: Path, pocket: dict, workdir: Path,
         idx: int, exhaustiveness: int) -> float | None:
    lig = workdir / f"lig_{idx}.pdbqt"
    out = workdir / f"dock_{idx}.pdbqt"
    try:
        subprocess.run(["obabel", f"-:{smiles}", "-O", str(lig), "--gen3d", "-p", "7.4"],
                       check=True, capture_output=True, timeout=120)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    if not lig.is_file() or lig.stat().st_size == 0:
        return None
    cx, cy, cz = pocket["center"]; sx, sy, sz = pocket["size"]
    cmd = ["vina", "--receptor", str(receptor), "--ligand", str(lig),
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
           "--exhaustiveness", str(exhaustiveness), "--cpu", "4", "--out", str(out)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return None
    m = re.search(r"VINA RESULT:\s*(-?\d+\.\d+)", out.read_text())
    return float(m.group(1)) if m else None


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="1PY5")
    ap.add_argument("--target", default="ALK5")
    ap.add_argument("--generate", type=int, default=150)
    ap.add_argument("--dock-top", type=int, default=8)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    workdir = RUN_ROOT / args.pdb
    workdir.mkdir(parents=True, exist_ok=True)
    pdb_path = fetch_pdb(args.pdb, workdir)

    pocket = detect_pocket(pdb_path, workdir)
    receptor = prep_receptor(pdb_path, workdir, pocket)

    raw = generate(args.generate, args.seed)
    log(f"  generated {len(raw)} non-empty SMILES")

    cands = []
    seen = set()
    for smi, nll in raw:
        if smi in seen:
            continue
        seen.add(smi)
        d = drug_like(smi)
        if d and d["pass"]:
            d.update(smiles=smi, nll=round(nll, 1))
            cands.append(d)
    cands.sort(key=lambda c: -c["qed"])
    log(f"  {len(cands)} unique drug-like candidates (PAH profile)")

    top = cands[: args.dock_top]
    log(f"docking top {len(top)} with AutoDock Vina (exhaustiveness={args.exhaustiveness})")
    for i, c in enumerate(top):
        score = dock(c["smiles"], receptor, pocket, workdir, i, args.exhaustiveness)
        c["dock_score"] = score
        log(f"  [{i + 1}/{len(top)}] QED {c['qed']:.2f}  dock "
            f"{score if score is not None else 'FAIL':>6}  {c['smiles'][:48]}")

    docked = [c for c in top if c.get("dock_score") is not None]
    docked.sort(key=lambda c: c["dock_score"])  # most negative first

    result = {
        "target": args.target, "pdb": args.pdb, "pocket": pocket,
        "generated": len(raw), "drug_like": len(cands),
        "docked": len(docked),
        "ranked": [
            {"rank": i + 1, "smiles": c["smiles"], "dock_score": c["dock_score"],
             "qed": round(c["qed"], 2), "mw": round(c["mw"], 1),
             "logp": round(c["logp"], 2), "tpsa": round(c["tpsa"], 1),
             "rot": c["rot"], "nll": c["nll"]}
            for i, c in enumerate(docked)
        ],
    }
    out_json = REPO / "reports" / f"pah_run_{args.pdb}.json"
    out_json.parent.mkdir(exist_ok=True)
    out_json.write_text(json.dumps(result, indent=2))
    log(f"wrote {out_json}")
    print("\n=== RANKED CANDIDATES (" + args.target + " / " + args.pdb + ") ===")
    print(f"{'rank':>4} {'dock':>7} {'QED':>5} {'MW':>6} {'cLogP':>6} {'TPSA':>6}  SMILES")
    for r in result["ranked"]:
        print(f"{r['rank']:>4} {r['dock_score']:>7} {r['qed']:>5} {r['mw']:>6} "
              f"{r['logp']:>6} {r['tpsa']:>6}  {r['smiles']}")
    print(json.dumps({"summary": {k: result[k] for k in
          ('target', 'pdb', 'generated', 'drug_like', 'docked')}}))


if __name__ == "__main__":
    main()
