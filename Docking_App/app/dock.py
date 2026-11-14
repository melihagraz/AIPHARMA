"""AutoDock Vina docking driver.

End-to-end pipeline:
    SMILES -> RDKit 3D conformer -> Meeko PDBQT -> Vina subprocess -> parse affinity

Receptor PDBQT files and grid-box configs are pre-bundled under
``Docking_App/receptors/<TARGET>_target.pdbqt`` and ``<TARGET>_conf.txt``
(sourced from dockstring; Apache 2.0).
"""

from __future__ import annotations

import logging
import hashlib
import os
import re
import shutil
import sys
import subprocess
import tempfile
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import List, Optional, Tuple

# Silence Meeko's deprecation warning (non-fatal, library-level)
warnings.filterwarnings("ignore", category=DeprecationWarning, module="meeko")

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Filesystem layout
# ---------------------------------------------------------------------------
_PKG_ROOT = Path(__file__).resolve().parent.parent
_VINA = _PKG_ROOT / "bin" / "vina.exe"  # Windows-shipped; on Linux replace with "vina"
_RECEPTORS_DIR = _PKG_ROOT / "receptors"


def _vina_binary() -> str:
    """Return the platform-appropriate Vina executable path.

    Windows: bin/vina.exe
    Linux/macOS: bin/vina (binary downloaded from CCSB releases page)
    """
    win = _PKG_ROOT / "bin" / "vina.exe"
    nix = _PKG_ROOT / "bin" / "vina"
    if os.name == "nt":
        if win.is_file():
            return str(win)
        found = shutil.which("vina.exe") or shutil.which("vina")
    else:
        if nix.is_file() and os.access(nix, os.X_OK):
            return str(nix)
        found = shutil.which("vina")
    if found:
        return found
    raise FileNotFoundError(
        f"No executable AutoDock Vina for {sys.platform}; install platform-appropriate "
        "Vina or provide Docking_App/bin/vina (Linux/macOS) or vina.exe (Windows)"
    )


def list_targets() -> List[str]:
    """Return targets that have BOTH a *_target.pdbqt and *_conf.txt available."""
    if not _RECEPTORS_DIR.is_dir():
        return []
    targets = set()
    for f in _RECEPTORS_DIR.iterdir():
        m = re.match(r"^(.+)_target\.pdbqt$", f.name)
        if m:
            t = m.group(1)
            conf = _RECEPTORS_DIR / f"{t}_conf.txt"
            if conf.exists():
                targets.add(t)
    return sorted(targets)


def _resolve_target(target: str) -> tuple[Path, Path]:
    """Map a target name (e.g. 'EGFR' or 'egfr') to (pdbqt_path, conf_path).

    Tries exact match first, then case-insensitive. Raises ValueError if missing.
    """
    direct = _RECEPTORS_DIR / f"{target}_target.pdbqt"
    direct_conf = _RECEPTORS_DIR / f"{target}_conf.txt"
    if direct.exists() and direct_conf.exists():
        return direct, direct_conf
    # case-insensitive
    lookup = target.upper()
    for t in list_targets():
        if t.upper() == lookup:
            return _RECEPTORS_DIR / f"{t}_target.pdbqt", _RECEPTORS_DIR / f"{t}_conf.txt"
    raise ValueError(
        f"Receptor for target {target!r} not found. Available: {list_targets()}"
    )


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------
def smiles_to_pdbqt(smiles: str) -> Optional[str]:
    """Convert SMILES → PDBQT string. Returns None on parse / embed failure."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    RDLogger.DisableLog("rdApp.*")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Stable ligand conformer initialization: the Vina seed below is fixed,
    # but it cannot make docking reproducible if the INPUT conformer changes.
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    seed = (int.from_bytes(hashlib.sha256(canonical.encode()).digest()[:4], "big") % 2147483646) + 1
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Deterministic fallback, not RDKit's random default.
        if AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True) != 0:
            return None
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:  # noqa: BLE001
        pass  # A failed minimization is disclosed only via ligand validation.

    # Meeko >=0.5 returns MoleculeSetup objects from prepare(). The old
    # prep.setup attribute is deprecated and may disappear in future releases.
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if len(setups) != 1:
        # Never silently choose an arbitrary reactive/multiple setup.
        return None
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok or not pdbqt or not any(line.startswith(("ATOM", "HETATM")) for line in pdbqt.splitlines()):
        logger.warning("Meeko ligand preparation failed: %s", str(err)[:180])
        return None
    return pdbqt


# ---------------------------------------------------------------------------
# Docking
# ---------------------------------------------------------------------------
@dataclass
class DockResult:
    smiles: str
    target: str
    affinity: Optional[float]  # kcal/mol; None on failure
    error: Optional[str] = None
    num_modes: int = 0
    # Top-pose 3D coordinates as a vanilla PDB block (the AutoDock-specific
    # charge + atom-type columns are stripped). Allows the calling backend
    # to persist the docked pose alongside the affinity so the Results
    # page can render the actual binding geometry instead of an RDKit
    # free-state conformer. ~3-15 KB per ligand; for a 30-molecule batch
    # we add ~150 KB to the AuditLog payload, well within budget.
    pose_pdb: Optional[str] = None
    # All Vina poses concatenated into a single multi-MODEL PDB block —
    # roughly num_modes × top-pose size, so 5 modes ≈ 15-75 KB. Lets the
    # frontend cycle through alternate binding poses without an extra
    # round-trip; None when capture failed.
    pose_pdb_all: Optional[str] = None


def _pdbqt_block_to_pdb(pdbqt_text: str, *, only_first_model: bool = True) -> str:
    """Convert a Vina out.pdbqt (multi-MODEL) block to vanilla PDB text.

    PDB columns 1-66 carry the standard ATOM/HETATM record (name, resname,
    chain, xyz, occupancy, tempFactor). AutoDock's PDBQT extension lives
    in columns 67+ (partial charge + atom type) and confuses most PDB
    parsers including NGL.js. Truncating each line at column 66 yields
    a clean PDB block while preserving every coordinate.

    Vina writes multiple MODEL blocks (one per pose) ranked top-down by
    affinity. ``only_first_model=True`` keeps just MODEL 1, which is the
    pose whose affinity the Results table already shows. Pass False to
    return all poses in a single multi-model PDB if the frontend ever
    wants a pose-cycler widget.
    """
    out_lines: list[str] = []
    saw_first_model = False
    for raw in pdbqt_text.splitlines():
        if raw.startswith("MODEL"):
            if only_first_model and saw_first_model:
                break
            saw_first_model = True
            out_lines.append(raw)
            continue
        if raw.startswith("ENDMDL"):
            out_lines.append(raw)
            if only_first_model:
                break
            continue
        if raw.startswith(("ATOM", "HETATM")):
            out_lines.append(raw[:66].rstrip())
        elif raw.startswith(("REMARK", "TER", "CONECT")):
            out_lines.append(raw)
    if not any(l.startswith(("ATOM", "HETATM")) for l in out_lines):
        return ""
    if not out_lines or not out_lines[-1].startswith("END"):
        out_lines.append("END")
    return "\n".join(out_lines)


# Vina prints rmsd columns as either "0" or "1.23" — the integer-zero on the
# top mode (which has no rmsd vs itself) means we cannot require decimal points.
_AFFINITY_LINE_RE = re.compile(r"^\s*1\s+(-?\d+\.\d+)\s+\S+\s+\S+\s*$")


def _parse_best_affinity(stdout: str) -> Optional[float]:
    """Extract the top mode's affinity (kcal/mol) from Vina's stdout table."""
    for line in stdout.splitlines():
        m = _AFFINITY_LINE_RE.match(line)
        if m:
            return float(m.group(1))
    return None


def dock_one(
    smiles: str,
    target: str,
    exhaustiveness: int = 8,
    num_modes: int = 5,
    cpu: int = 4,
    timeout: float = 180.0,
    pocket_box: Optional[Tuple[float, float, float, float, float, float]] = None,
    custom_receptor_pdbqt: Optional[Path] = None,
) -> DockResult:
    """Dock a single ligand against a pre-bundled target. Returns DockResult.

    ``pocket_box`` is an optional ``(center_x, center_y, center_z, size_x,
    size_y, size_z)`` tuple. When supplied, the corresponding Vina CLI flags
    override whatever centre/size the bundled ``conf.txt`` declares — so the
    docking happens inside the *user-detected* pocket rather than the default
    ATP-site grid. When ``None`` we fall back to the bundled grid (legacy).

    ``custom_receptor_pdbqt`` overrides the receptor file lookup entirely:
    when supplied, Vina is given this PDBQT (typically produced upstream by
    ``receptor_prep.prepare_custom_receptor`` from a user-uploaded PDB) and
    the bundled receptor for ``target`` is ignored. ``target`` is still used
    as a logging / cache label.
    """
    if custom_receptor_pdbqt is not None:
        receptor_path = custom_receptor_pdbqt
        # No bundled conf for custom receptors; pocket_box is *mandatory* in
        # this mode because Vina needs a search box from somewhere.
        if pocket_box is None:
            return DockResult(
                smiles=smiles, target=target, affinity=None,
                error="custom_receptor requires pocket_box (centre + size)",
            )
        conf_path = None
    else:
        try:
            receptor_path, conf_path = _resolve_target(target)
        except ValueError as e:
            return DockResult(smiles=smiles, target=target, affinity=None, error=str(e))

    pdbqt = smiles_to_pdbqt(smiles)
    if pdbqt is None:
        return DockResult(smiles=smiles, target=target, affinity=None, error="ligand_prep_failed")

    with tempfile.TemporaryDirectory(prefix="dock_") as td:
        td = Path(td)
        lig_path = td / "ligand.pdbqt"
        out_path = td / "out.pdbqt"
        lig_path.write_text(pdbqt)
        cmd = [
            _vina_binary(),
            "--receptor", str(receptor_path),
            "--ligand", str(lig_path),
            "--out", str(out_path),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--cpu", str(cpu),
            "--seed", "42",  # reproducibility
        ]
        # Bundled receptors come with a conf.txt that declares the search
        # grid (centre + size). Custom receptors skip --config because they
        # have no bundled conf — the caller must supply pocket_box.
        if conf_path is not None:
            cmd.extend(["--config", str(conf_path)])
        # Pocket-box override: append CLI flags AFTER --config so they win.
        # Vina's argument parser is "last-flag-wins" for center_x/size_x/etc.,
        # which means a user-supplied pocket from Pocket Analysis can target
        # an allosteric or mutant-specific site that the bundled grid doesn't
        # cover.
        if pocket_box is not None:
            cx, cy, cz, sx, sy, sz = pocket_box
            cmd.extend([
                "--center_x", f"{cx:.3f}",
                "--center_y", f"{cy:.3f}",
                "--center_z", f"{cz:.3f}",
                "--size_x",   f"{sx:.3f}",
                "--size_y",   f"{sy:.3f}",
                "--size_z",   f"{sz:.3f}",
            ])
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return DockResult(smiles=smiles, target=target, affinity=None, error="vina_timeout")
        if r.returncode != 0:
            return DockResult(
                smiles=smiles,
                target=target,
                affinity=None,
                error=f"vina_failed_rc{r.returncode}: {r.stderr[:300]}",
            )
        affinity = _parse_best_affinity(r.stdout)
        if affinity is None:
            return DockResult(smiles=smiles, target=target, affinity=None,
                              error="vina_output_missing_valid_affinity")
        # Capture docked poses so the API caller can persist them
        # alongside the affinity. We grab BOTH the top-1 pose (cheap,
        # always rendered) and the full multi-model block (lets the
        # frontend cycle through alternative binding modes). The
        # out_path is still inside the temp dir so we MUST read it
        # before this `with` block exits and the dir is wiped. Any
        # failure to read / convert is swallowed — docking itself
        # succeeded, the missing pose just means the Results page
        # falls back to the RDKit free-state conformer.
        pose_pdb: Optional[str] = None
        pose_pdb_all: Optional[str] = None
        try:
            if out_path.exists() and out_path.stat().st_size > 0:
                raw_pdbqt = out_path.read_text(encoding="utf-8", errors="ignore")
                converted = _pdbqt_block_to_pdb(raw_pdbqt, only_first_model=True)
                if converted:
                    pose_pdb = converted
                # All-pose block keeps every MODEL in the same PDB so the
                # NGL viewer can switch between them without an extra
                # round-trip. Skipped (None) if the converter returned
                # an empty block.
                full = _pdbqt_block_to_pdb(raw_pdbqt, only_first_model=False)
                if full and full != pose_pdb:
                    pose_pdb_all = full
        except Exception as exc:  # noqa: BLE001
            logger.warning("dock_one pose capture failed: %s", exc)
        return DockResult(
            smiles=smiles, target=target, affinity=affinity,
            num_modes=num_modes, pose_pdb=pose_pdb, pose_pdb_all=pose_pdb_all,
        )


# Module-level lock to serialise heavy concurrent prep work; Vina subprocess
# itself is fine to parallelise but RDKit's embedding can be racy in some builds.
_prep_lock = Lock()


def dock_batch(
    smiles_list: List[str],
    target: str,
    exhaustiveness: int = 8,
    num_modes: int = 5,
    cpu_per_dock: int = 2,
    workers: int = 4,
    timeout_per_ligand: float = 180.0,
    pocket_box: Optional[Tuple[float, float, float, float, float, float]] = None,
    custom_receptor_pdbqt: Optional[Path] = None,
) -> List[DockResult]:
    """Dock a batch of ligands in parallel. Returns results in input order.

    ``pocket_box`` is forwarded to each ``dock_one`` call — see its docstring
    for the override semantics. ``custom_receptor_pdbqt`` forwards a
    user-prepared receptor PDBQT (typically produced by
    ``receptor_prep.prepare_custom_receptor``) so every worker docks against
    the SAME custom receptor without re-resolving the bundled file each time.
    """
    if not smiles_list:
        return []
    # Only validate the bundled target when we'll actually be using it.
    if custom_receptor_pdbqt is None:
        _resolve_target(target)

    results: List[Optional[DockResult]] = [None] * len(smiles_list)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(
                dock_one, smi, target,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                cpu=cpu_per_dock,
                timeout=timeout_per_ligand,
                pocket_box=pocket_box,
                custom_receptor_pdbqt=custom_receptor_pdbqt,
            ): i
            for i, smi in enumerate(smiles_list)
        }
        for fut in as_completed(futures):
            i = futures[fut]
            try:
                results[i] = fut.result()
            except Exception as e:  # noqa: BLE001
                results[i] = DockResult(
                    smiles=smiles_list[i], target=target, affinity=None, error=str(e)
                )
    return [r for r in results if r is not None]
