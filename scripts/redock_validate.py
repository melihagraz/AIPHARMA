#!/usr/bin/env python
"""Redocking (self-docking) validation.

Takes a PDB, removes its co-crystallized ligand, re-docks that same ligand
into the binding site, and measures the RMSD between the docked pose and the
original crystallographic pose. A low RMSD (< 2.0 A, the accepted community
threshold) demonstrates that the docking protocol reproduces the experimental
binding mode -- i.e. the method is reproducible.

Runs inside the WSL conda env ``aipharma`` (obabel + autodock-vina + obrms).

Usage:
    mamba run -n aipharma python scripts/redock_validate.py --pdb 1PY5 --exhaustiveness 16
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
RUN = Path.home() / ".aipharma_redock"
NON_LIGAND = {
    "HOH", "WAT", "DOD", "NA", "CL", "MG", "ZN", "CA", "K", "FE", "CU", "MN", "NI",
    "SO4", "PO4", "GOL", "EDO", "PEG", "DMS", "TRS", "ACE", "FMT", "BME", "MES",
    "EPE", "BTB", "CIT", "TRI", "MPD", "NH2", "ACT", "EOH", "MOH", "IOD", "BR",
}


# Catalytic metal ions: resname -> (AutoDock atom type, formal charge). AutoDock
# receptor prep drops or mis-charges these, so metalloproteins (Zn/Mg/Mn active
# sites) dock wrong unless the ion is kept and given its formal charge.
METALS = {
    "ZN": ("Zn", 2.0), "MG": ("Mg", 2.0), "MN": ("Mn", 2.0), "FE": ("Fe", 2.0),
    "CA": ("Ca", 2.0), "NI": ("Ni", 2.0), "CU": ("Cu", 2.0), "CO": ("Co", 2.0),
    "NA": ("Na", 1.0), "K": ("K", 1.0),
}


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def fix_metal_charges(pdbqt_path: Path) -> int:
    """obabel keeps metals with the right AutoDock type but a 0.000 charge.
    Overwrite that with the formal charge (Zn/Mg/... +2, Na/K +1). Returns the
    number of metal atoms fixed."""
    types = {t: q for (t, q) in METALS.values()}
    n = 0
    out = []
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line.split():
            adtype = line.split()[-1]
            if adtype in types:
                newc = f"{types[adtype]:+.3f}"
                line, k = re.subn(
                    r"([+-]?\d+\.\d{3})(\s+" + re.escape(adtype) + r")\s*$",
                    newc + r"\2", line)
                n += k
        out.append(line)
    pdbqt_path.write_text("\n".join(out) + "\n")
    return n


def fetch(pdb, wd):
    dest = wd / f"{pdb}.pdb"
    if not dest.is_file():
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb.upper()}.pdb", dest)  # noqa: S310
    return dest


def split_structure(pdb_path, wd):
    """Write receptor (protein + active-site metal ions) and the largest
    non-buffer HETATM ligand. Metal ions within 12 A of the ligand are kept in
    the receptor so metalloprotein pockets dock correctly."""
    import math
    groups = {}
    protein = []
    metals = []  # (resname, (x,y,z), raw_line)
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM", "TER")):
            protein.append(line)
        elif line.startswith("HETATM"):
            res = line[17:20].strip().upper()
            if res in METALS:
                try:
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    metals.append((res, xyz, line))
                except ValueError:
                    pass
                continue
            if res in NON_LIGAND:
                continue
            chain, seq = line[21:22], line[22:26]
            groups.setdefault((res, chain, seq), []).append(line)
    if not groups:
        raise SystemExit("No co-crystal ligand found (apo structure) -- use fpocket pocket detection instead.")
    (res, chain, seq), lig_lines = max(groups.items(), key=lambda kv: len(kv[1]))

    lc = tuple(sum(float(l[c:c + 8]) for l in lig_lines) / len(lig_lines) for c in (30, 38, 46))
    near = [m for m in metals if math.dist(m[1], lc) <= 12.0]

    rec = wd / "receptor.pdb"
    rec.write_text("\n".join(protein + [m[2] for m in near]) + "\nEND\n")
    lig = wd / "crystal_lig.pdb"
    lig.write_text("\n".join(lig_lines) + "\nEND\n")
    return rec, lig, res, len(lig_lines), [m[0] for m in near]


def prep_ligand_meeko(lig_pdb_path, resname, out_pdbqt):
    """Proper Vina-ready ligand PDBQT: CCD bond orders + explicit H + Meeko
    (Gasteiger charges, torsion tree). obabel's naive prep loses the correct
    protonation/charge pattern, which is what makes the docked orientation
    wrong even when the location is right."""
    sys.path.insert(0, str(REPO / "Pocket_App"))
    from app.pharmacophore import _build_ligand_mol
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = _build_ligand_mol(lig_pdb_path.read_text(), resname)  # 3D coords + bond orders
    mol = Chem.AddHs(mol, addCoords=True)
    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise RuntimeError(f"Meeko ligand prep failed: {err}")
    Path(out_pdbqt).write_text(pdbqt)


def rmsd_inplace(ref_pdb, docked_pdb):
    """Symmetry-corrected, in-place (no superposition) heavy-atom RMSD — the
    correct metric for redocking pose reproduction."""
    from spyrmsd import io as sio, rmsd as srmsd
    ref = sio.loadmol(str(ref_pdb))
    mol = sio.loadmol(str(docked_pdb))
    return float(srmsd.symmrmsd(
        ref.coordinates, mol.coordinates,
        ref.atomicnums, mol.atomicnums,
        ref.adjacency_matrix, mol.adjacency_matrix, minimize=False))


def centroid_box(lig_path, margin=8.0):
    """Snug per-axis box centered on the ligand (redocking convention).

    A box much larger than the ligand lets the search wander off-site; here
    each axis is sized to the ligand's extent plus a small margin and clamped
    to a sensible [16, 26] A window.
    """
    xs, ys, zs = [], [], []
    for line in lig_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
    cx, cy, cz = (sum(xs) / len(xs)), (sum(ys) / len(ys)), (sum(zs) / len(zs))

    def ax(v):
        return round(max(16.0, min(26.0, (max(v) - min(v)) + 2 * margin)), 1)

    return (round(cx, 3), round(cy, 3), round(cz, 3)), (ax(xs), ax(ys), ax(zs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", default="1PY5")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    args = ap.parse_args()

    wd = RUN / args.pdb
    wd.mkdir(parents=True, exist_ok=True)
    pdb_path = fetch(args.pdb, wd)
    rec_pdb, lig_pdb, resname, n_atoms, metals = split_structure(pdb_path, wd)
    metal_note = f"; active-site metals: {', '.join(metals)}" if metals else ""
    print(f"[redock] {args.pdb}: ligand {resname} ({n_atoms} atoms) extracted{metal_note}")

    # receptor + ligand prep
    rec_qt = wd / "receptor.pdbqt"
    sh(["obabel", str(rec_pdb), "-O", str(rec_qt), "-xr", "-p", "7.4"])
    n_fixed = fix_metal_charges(rec_qt)
    if n_fixed:
        print(f"[redock] assigned formal charges to {n_fixed} metal ion(s) in the receptor")
    lig_qt = wd / "ligand.pdbqt"
    prep_ligand_meeko(lig_pdb, resname, lig_qt)

    (cx, cy, cz), (sx, sy, sz) = centroid_box(lig_pdb)
    print(f"[redock] box center ({cx},{cy},{cz}) size ({sx},{sy},{sz}) A; docking (exhaustiveness={args.exhaustiveness})")

    out_qt = wd / "docked.pdbqt"
    sh(["vina", "--receptor", str(rec_qt), "--ligand", str(lig_qt),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(args.exhaustiveness), "--cpu", "4", "--out", str(out_qt)])
    score = None
    m = re.search(r"VINA RESULT:\s*(-?\d+\.\d+)", out_qt.read_text())
    if m:
        score = float(m.group(1))

    # top pose -> heavy-atom pdb, then RMSD vs crystal (heavy atoms)
    pose_pdb = wd / "docked_pose.pdb"
    sh(["obabel", str(out_qt), "-O", str(pose_pdb), "-f", "1", "-l", "1", "-d"])
    ref_pdb = wd / "crystal_lig_noH.pdb"
    sh(["obabel", str(lig_pdb), "-O", str(ref_pdb), "-d"])
    rmsd = round(rmsd_inplace(ref_pdb, pose_pdb), 3)

    ok = rmsd is not None and rmsd < 2.0
    result = {"pdb": args.pdb, "ligand": resname, "ligand_atoms": n_atoms,
              "dock_score": score, "rmsd": rmsd, "reproduced": ok,
              "metals": metals,
              "box_center": [cx, cy, cz], "box_size": [sx, sy, sz],
              "exhaustiveness": args.exhaustiveness}
    outp = Path(__file__).resolve().parent.parent / "reports" / f"redock_{args.pdb}.json"
    outp.parent.mkdir(exist_ok=True)
    outp.write_text(json.dumps(result, indent=2))
    # copy poses next to the JSON so the UI (Windows side) can load them
    import shutil
    shutil.copy(pose_pdb, outp.parent / f"redock_{args.pdb}_docked.pdb")
    shutil.copy(ref_pdb, outp.parent / f"redock_{args.pdb}_crystal.pdb")

    print("\n=== REDOCKING VALIDATION ===")
    print(f"  PDB              : {args.pdb}")
    print(f"  Ligand           : {resname} ({n_atoms} atoms)")
    print(f"  Dock score       : {score} kcal/mol")
    print(f"  Pose RMSD to xtal : {rmsd} A")
    print(f"  Reproduced (<2 A): {'YES' if ok else 'NO'}")
    print(f"  -> {outp}")


if __name__ == "__main__":
    main()
