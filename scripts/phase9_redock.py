#!/usr/bin/env python3
"""Phase 9: real crystallographic redocking of EGFR and ALK5.

Strictly experimental: requires actual RCSB PDB, Meeko, AutoDock Vina,
Open Babel, spyrmsd, RDKit. Does not report success without real runs.
Run from the complete AIPHARMA repository with this script in scripts/.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.request import urlopen, Request

TARGETS = {
    "EGFR": {"pdb": "1M17", "ligand": "AQ4", "source": "https://www.rcsb.org/structure/1M17"},
    "ALK5": {"pdb": "1VJY", "ligand": "460", "source": "https://www.rcsb.org/structure/1VJY"},
}


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ligand_centroid_box(lines, pad=8.0, min_size=20.0, max_size=30.0):
    coords = []
    for line in lines:
        if line.startswith("HETATM"):
            coords.append(tuple(float(line[i:i + 8]) for i in (30, 38, 46)))
    if len(coords) < 5:
        raise ValueError("Too few co-crystallized ligand atoms")
    xyz = list(zip(*coords))
    center = [round(sum(axis) / len(axis), 3) for axis in xyz]
    size = [round(max(min_size, min(max_size, max(axis) - min(axis) + 2 * pad)), 3) for axis in xyz]
    return center, size


def split_pdb(source: Path, ligand_code: str, protein_out: Path, ligand_out: Path) -> dict:
    lines = source.read_text().splitlines()
    # PDB asymmetric unit can contain several chain/ligand copies. Choose the
    # ligand with the requested CCD ID in chain A; never select a random HETATM.
    residue_groups = {}
    protein = []
    first_model_only = False
    for raw in lines:
        if raw.startswith("MODEL") and first_model_only:
            break
        if raw.startswith("ENDMDL"):
            first_model_only = True
            continue
        if raw.startswith("ATOM  ") and raw[21:22] == "A" and raw[16:17] in (" ", "A"):
            protein.append(raw)
        if raw.startswith("HETATM") and raw[17:20].strip() == ligand_code:
            if raw[16:17] not in (" ", "A"):
                continue
            key = (raw[21:22], raw[22:26], raw[26:27])
            residue_groups.setdefault(key, []).append(raw)
    if not protein:
        raise ValueError("No receptor protein chain A ATOM records")
    if not residue_groups:
        raise ValueError(f"Expected crystallographic ligand {ligand_code} not found")
    # No silent ligand switching: ambiguity is recorded, preferred chain A.
    keys = sorted(residue_groups, key=lambda k: (k[0] != "A", -len(residue_groups[k]), k))
    key = keys[0]
    ligand = residue_groups[key]
    protein_out.write_text("\n".join(protein + ["TER", "END", ""]) )
    ligand_out.write_text("\n".join(ligand + ["END", ""]))
    center, size = ligand_centroid_box(ligand)
    return {"chosen_ligand_residue": key, "ligand_copies": len(residue_groups),
            "ligand_heavy_atom_count": len(ligand), "protein_atom_count": len(protein),
            "center": center, "size": size}


def command(argv, cwd=None, timeout=240):
    proc = subprocess.run([str(x) for x in argv], cwd=cwd, text=True,
                          capture_output=True, timeout=timeout, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed (exit {proc.returncode}): {' '.join(map(str, argv))}\n"
                           f"STDERR: {proc.stderr[-1300:]}\nSTDOUT: {proc.stdout[-700:]}")
    return proc.stdout[-2500:]


def fetch_pdb(pdb: str, target_file: Path, offline=False):
    if target_file.is_file() and target_file.stat().st_size > 10000:
        return
    if offline:
        raise FileNotFoundError(f"Missing cached PDB: {target_file}")
    url = f"https://files.rcsb.org/download/{pdb}.pdb"
    req = Request(url, headers={"User-Agent": "AIPHARMA-reproducibility-phase9/1.0"})
    with urlopen(req, timeout=35) as response:
        content = response.read(10_000_000)
    if b"ATOM" not in content or b"HEADER" not in content:
        raise ValueError("RCSB response does not look like a PDB file")
    target_file.write_bytes(content)


def run_target(target: str, out_dir: Path, prior_root: Path, offline=False, exhaustiveness=16):
    from importlib.util import spec_from_file_location, module_from_spec
    redock_path = prior_root / "scripts" / "redock_validate.py"
    spec = spec_from_file_location("aipharma_legacy_redock_helpers", redock_path)
    helper = module_from_spec(spec)
    spec.loader.exec_module(helper)
    specinfo = TARGETS[target]
    pdb, ligand_code = specinfo["pdb"], specinfo["ligand"]
    wd = out_dir / pdb
    wd.mkdir(parents=True, exist_ok=True)
    rep = {"target": target, "pdb_id": pdb, "ccd_ligand": ligand_code,
           "date_utc": datetime.now(timezone.utc).isoformat(),
           "crystal_source": specinfo["source"], "passed": False,
           "experiment": "redocking; not binding affinity or experimental efficacy validation",
           "commands": [], "blocker": None}
    report_path = wd / "redock_report.json"
    try:
        needed = ["vina", "obabel", "mk_prepare_receptor.py"]
        missing = [exe for exe in needed if shutil.which(exe) is None]
        if missing:
            raise RuntimeError(f"Missing executables: {', '.join(missing)}")
        for module in ["spyrmsd", "meeko", "rdkit"]:
            __import__(module)
        pdbfile = wd / f"{pdb}.pdb"
        fetch_pdb(pdb, pdbfile, offline=offline)
        recpdb, ligpdb = wd / "protein_chainA.pdb", wd / "co_crystal_ligand.pdb"
        geom = split_pdb(pdbfile, ligand_code, recpdb, ligpdb)
        rep.update(geom)
        rep["hashes"] = {name: digest(p) for name, p in
                          [("original_pdb", pdbfile), ("protein_pdb", recpdb), ("ligand_pdb", ligpdb)]}
        # The default Meeko strict mode fails rather than quietly deleting any
        # incomplete kinase residues. This is preferable to an unreviewed receptor.
        receptor = wd / "prepared_receptor.pdbqt"
        prep_argv = ["mk_prepare_receptor.py", "--read_pdb", recpdb,
                     "--write_pdbqt", receptor]
        # The EGFR crystal structure 1M17 contains alternate A/B locations at
        # residues A:751 and A:831. split_pdb() keeps altloc A explicitly; pass
        # these matching choices to Meeko rather than deleting residues.
        if target == "EGFR":
            prep_argv += ["--wanted_altloc", "A:751=A,A:831=A"]
            rep["selected_altlocs"] = {"A:751": "A", "A:831": "A"}
        rep["commands"].append(" ".join(map(str, prep_argv)))
        rep["receptor_prep_output"] = command(prep_argv, timeout=240)
        if not receptor.is_file() or not any(x.startswith("ATOM") for x in receptor.read_text().splitlines()):
            raise ValueError("Receptor PDBQT missing or empty")
        # PDB coordinates alone cannot establish the aromatic/heteroatom bond
        # orders; require a live CCD chemical template, do not accept a guessed
        # fallback from the helper in an experiment presented as redocking.
        sys.path.insert(0, str(prior_root / "Pocket_App"))
        from app.pharmacophore import _fetch_ccd_smiles
        from rdkit import Chem
        ccd_smiles = _fetch_ccd_smiles(ligand_code)
        if not ccd_smiles:
            raise RuntimeError(f"RCSB CCD bond-order template unavailable for {ligand_code}")
        ccd_mol = Chem.MolFromSmiles(ccd_smiles)
        if ccd_mol is None or ccd_mol.GetNumHeavyAtoms() != geom["ligand_heavy_atom_count"]:
            raise RuntimeError("Co-crystallized ligand atom count differs from CCD template; manual review required")
        rep["ccd_template_smiles"] = ccd_smiles
        ligqt = wd / "prepared_ligand.pdbqt"
        helper.prep_ligand_meeko(ligpdb, ligand_code, ligqt)
        if not ligqt.is_file():
            raise RuntimeError("Meeko ligand conversion did not produce PDBQT")
        center, size = geom["center"], geom["size"]
        config = wd / "vina_grid.txt"
        config.write_text("\n".join([f"center_{ax} = {v:.3f}" for ax, v in zip("xyz", center)] +
                                  [f"size_{ax} = {v:.3f}" for ax, v in zip("xyz", size)]) + "\n")
        out_pose = wd / "redocked.pdbqt"
        argv = ["vina", "--receptor", receptor, "--ligand", ligqt,
                "--config", config, "--exhaustiveness", str(exhaustiveness),
                "--cpu", "1", "--seed", "20261109", "--num_modes", "9", "--out", out_pose]
        rep["commands"].append(" ".join(map(str, argv)))
        rep["vina_stdout"] = command(argv, timeout=1200)
        if not out_pose.is_file() or "VINA RESULT:" not in out_pose.read_text():
            raise RuntimeError("Docking produced no scored poses")
        top = wd / "top_pose.pdb"
        ref = wd / "crystal_noH.pdb"
        command(["obabel", out_pose, "-O", top, "-f", "1", "-l", "1", "-d"])
        command(["obabel", ligpdb, "-O", ref, "-d"])
        rmsd = round(float(helper.rmsd_inplace(ref, top)), 4)
        if not math.isfinite(rmsd):
            raise ValueError("Invalid redocking RMSD")
        score = re.search(r"VINA RESULT:\s*(-?\d+(?:\.\d+)?)", out_pose.read_text())
        rep["pose_rmsd_angstrom"] = rmsd
        rep["docking_score_kcal_mol"] = float(score.group(1)) if score else None
        rep["pose_reproduced_within_2A"] = rmsd < 2.0
        rep["passed"] = True  # executed and RMSD computed, not successful pose recovery
        rep["files_sha256"] = {p.name: digest(p) for p in (receptor, ligqt, config, out_pose, top, ref)}
    except Exception as e:
        rep["blocker"] = f"{type(e).__name__}: {str(e)[:1800]}"
    finally:
        report_path.write_text(json.dumps(rep, indent=2) + "\n")
    return rep


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--target", choices=["EGFR", "ALK5", "BOTH"], default="BOTH")
    p.add_argument("--out", type=Path, default=Path("reports/phase9/redocking"))
    p.add_argument("--root", type=Path, default=Path("."))
    p.add_argument("--offline", action="store_true", help="Use already-downloaded PDB inputs")
    p.add_argument("--exhaustiveness", type=int, default=16)
    a = p.parse_args(argv)
    if not (1 <= a.exhaustiveness <= 128):
        p.error("exhaustiveness must be 1 to 128")
    targets = ["EGFR", "ALK5"] if a.target == "BOTH" else [a.target]
    results = [run_target(t, a.out, a.root.resolve(), a.offline, a.exhaustiveness) for t in targets]
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "index.json").write_text(json.dumps(results, indent=2) + "\n")
    for r in results:
        print(f"{r['target']}: executed={r['passed']}, RMSD={r.get('pose_rmsd_angstrom')}, "
              f"pose<2A={r.get('pose_reproduced_within_2A')}, blocker={r.get('blocker')}")
    return 0 if all(r["passed"] for r in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
