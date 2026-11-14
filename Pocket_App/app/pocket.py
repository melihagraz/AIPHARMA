"""Pocket detection — pure-Python SBDD pocket extraction.

We avoid Java-based tools (p2rank, fpocket-windows, etc.) and instead use
two pragmatic strategies that cover real-world drug discovery:

1. **Bundled receptor lookup** — the 8 bundled targets (EGFR/JAK2/...) all
   ship with a Vina ``*_conf.txt`` defining the validated pocket center +
   bounding box. These come from the dockstring repository and were each
   curated against published co-crystal structures. For any of these
   targets we just read the conf file. Druggability score = 1.0.

2. **Co-crystal ligand center** — for any user-provided PDB with a bound
   small-molecule ligand, the bound ligand's center of geometry IS the
   pocket center, and the bounding box can be sized from the ligand's
   coordinates with a margin. This is how Schrödinger's Glide grid box
   is typically defined. Druggability score = 0.95 (high — a real ligand
   binds here, but pocket may be too tight without padding).

3. **Fallback** — if neither path works (apo PDB, no ligand) we return a
   centered guess based on protein center of mass, marked as low
   druggability so the caller knows to upload a co-crystal instead.

Skipped HETATM residues (treated as non-ligand by default):
    HOH, WAT, NA, CL, MG, ZN, CA, K, SO4, PO4, GOL, EDO, PEG, ACE, NH2,
    common buffer / cryoprotectant codes.
"""

from __future__ import annotations

import io
import logging
import re
import urllib.request
from pathlib import Path
from typing import Iterable, Optional

from .schemas import PocketResult

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

POCKET_APP_DIR = Path(__file__).resolve().parents[1]
RECEPTORS_DIR = POCKET_APP_DIR.parent / "Docking_App" / "receptors"
STRUCTURES_DIR = POCKET_APP_DIR / "structures"
STRUCTURES_DIR.mkdir(exist_ok=True)

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Common HETATMs that are NOT ligands of interest.
NON_LIGAND_RESNAMES = {
    # Water + ions
    "HOH", "WAT", "DOD",
    "NA", "CL", "MG", "ZN", "CA", "K", "FE", "CU", "MN", "NI",
    # Buffers / cryoprotectants
    "SO4", "PO4", "GOL", "EDO", "PEG", "DMS", "TRS", "ACE", "FMT",
    "BME", "MES", "EPE", "BTB", "CIT", "TRI", "MPD",
    # Modifiers / caps
    "NH2", "ACT", "EOH", "MOH",
}


# ---------------------------------------------------------------------------
# Bundled receptor lookup
# ---------------------------------------------------------------------------


def list_bundled_targets() -> list[str]:
    """Names of bundled targets that have a ``*_conf.txt`` file."""
    if not RECEPTORS_DIR.is_dir():
        return []
    out = []
    for f in RECEPTORS_DIR.glob("*_conf.txt"):
        name = f.stem.removesuffix("_conf")
        if name:
            out.append(name)
    return sorted(out)


def _parse_conf(path: Path) -> Optional[dict[str, float]]:
    """Parse a Vina ``conf.txt`` into a dict of floats."""
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    out: dict[str, float] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = re.match(r"^\s*(\w+)\s*=\s*([-+]?\d*\.?\d+)\s*$", line)
        if m:
            try:
                out[m.group(1).lower()] = float(m.group(2))
            except ValueError:
                continue
    return out


def pocket_from_bundled_target(target: str) -> Optional[PocketResult]:
    """Read the Vina conf for a bundled target and return its pocket."""
    if not target:
        return None
    target_norm = target.strip()
    # case-insensitive match against bundled list
    available = list_bundled_targets()
    match = None
    for t in available:
        if t.lower() == target_norm.lower():
            match = t
            break
    if match is None:
        return None
    conf = _parse_conf(RECEPTORS_DIR / f"{match}_conf.txt")
    if not conf:
        return None
    try:
        return PocketResult(
            center_x=conf["center_x"],
            center_y=conf["center_y"],
            center_z=conf["center_z"],
            size_x=conf.get("size_x", 24.0),
            size_y=conf.get("size_y", 24.0),
            size_z=conf.get("size_z", 24.0),
            druggability_score=1.0,
            source="bundled_receptor",
            notes=f"Bundled validated pocket for {match} (dockstring source).",
        )
    except KeyError:
        return None


# ---------------------------------------------------------------------------
# PDB fetching + parsing
# ---------------------------------------------------------------------------


def _fetch_pdb(pdb_id: str) -> Path:
    """Download a PDB from RCSB if not cached locally. Returns the local path."""
    pid = pdb_id.strip().lower()
    if not re.fullmatch(r"[a-z0-9]{4}", pid):
        raise ValueError(f"Invalid PDB ID: {pdb_id!r}")
    local = STRUCTURES_DIR / f"{pid}.pdb"
    if local.is_file() and local.stat().st_size > 0:
        return local
    url = RCSB_URL.format(pdb_id=pid.upper())
    logger.info("Fetching PDB %s from %s", pid, url)
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:  # noqa: S310
            # Remote archives are not trusted; do not fetch unbounded data.
            data = resp.read(20 * 1024 * 1024 + 1)
        if len(data) > 20 * 1024 * 1024:
            raise RuntimeError(f"PDB {pid} exceeds the 20 MiB download limit")
        if not any(record in data for record in (b"ATOM  ", b"HETATM")):
            raise RuntimeError(f"PDB {pid} contains no ATOM/HETATM records")
        local.write_bytes(data)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Failed to fetch PDB {pid}: {exc}") from exc
    return local


def _iter_hetatm_records(pdb_text: str) -> Iterable[tuple[str, list[tuple[float, float, float]]]]:
    """Yield (resname, [coords...]) for each HETATM residue.

    Groups atoms by (resname, chain, resseq).
    """
    groups: dict[tuple[str, str, str], list[tuple[float, float, float]]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        # PDB fixed-column format
        try:
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        groups.setdefault((resname, chain, resseq), []).append((x, y, z))
    for (resname, _chain, _resseq), coords in groups.items():
        yield resname, coords


def _iter_atom_records(pdb_text: str) -> Iterable[tuple[str, str, tuple[float, float, float]]]:
    """Yield (resname, atom_name, coords) for each protein ATOM record."""
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        yield resname, atom_name, (x, y, z)


def _bbox_and_center(
    coords: list[tuple[float, float, float]],
    margin: float = 8.0,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Compute (center, box_size) from a list of XYZ coords with margin."""
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    cx, cy, cz = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2
    sx = (max(xs) - min(xs)) + 2 * margin
    sy = (max(ys) - min(ys)) + 2 * margin
    sz = (max(zs) - min(zs)) + 2 * margin
    # Vina expects 1–30 Å typical; clamp to [16, 30] so we don't undersize tiny ligands
    sx = max(16.0, min(30.0, sx))
    sy = max(16.0, min(30.0, sy))
    sz = max(16.0, min(30.0, sz))
    return (cx, cy, cz), (sx, sy, sz)


def _pocket_residues_near(
    pdb_text: str, center: tuple[float, float, float], radius: float = 6.0
) -> list[str]:
    """Return up to ~25 protein residues whose any atom is within ``radius`` Å of ``center``."""
    cx, cy, cz = center
    seen: set[tuple[str, str]] = set()
    out: list[str] = []
    r2 = radius * radius
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        dx, dy, dz = x - cx, y - cy, z - cz
        if dx * dx + dy * dy + dz * dz <= r2:
            key = (chain, resseq)
            if key not in seen:
                seen.add(key)
                out.append(f"{resname}{resseq}")
                if len(out) >= 25:
                    break
    return out


# ---------------------------------------------------------------------------
# Co-crystal ligand → pocket
# ---------------------------------------------------------------------------


def pocket_from_pdb(
    pdb_id: str,
    ligand_resname: Optional[str] = None,
    prepared_pdb_path: Optional[str] = None,
) -> PocketResult:
    """Detect a pocket from a PDB file using the bound ligand's coordinates.

    When ``prepared_pdb_path`` is provided, that file is parsed directly and
    the RCSB fetch is skipped. The caller is expected to have run the file
    through a preparation pipeline (PDBFixer cleanup, protonation, water
    handling) before passing it in.
    """
    if prepared_pdb_path:
        prep_path = Path(prepared_pdb_path)
        if not prep_path.is_file():
            raise ValueError(
                f"prepared_pdb_path does not exist or is not a file: {prepared_pdb_path}"
            )
        pdb_path = prep_path
    else:
        pdb_path = _fetch_pdb(pdb_id)
    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")

    chosen_resname: Optional[str] = None
    chosen_coords: list[tuple[float, float, float]] = []

    candidates = list(_iter_hetatm_records(pdb_text))
    if ligand_resname:
        target_res = ligand_resname.strip().upper()
        for resname, coords in candidates:
            if resname.upper() == target_res:
                chosen_resname = resname
                chosen_coords = coords
                break
        if chosen_resname is None:
            available = sorted({r for r, _ in candidates})
            raise ValueError(
                f"Ligand {ligand_resname!r} not found in PDB {pdb_id}. "
                f"Available HETATM resnames: {available}"
            )
    else:
        # Pick the largest non-trivial HETATM (most atoms = the bound drug).
        best_size = 0
        for resname, coords in candidates:
            if resname.upper() in NON_LIGAND_RESNAMES:
                continue
            if len(coords) > best_size and len(coords) >= 6:
                best_size = len(coords)
                chosen_resname = resname
                chosen_coords = coords

    if not chosen_coords:
        # Apo / no recognisable ligand → centre-of-mass fallback
        all_atoms = [c for _r, _n, c in _iter_atom_records(pdb_text)]
        if not all_atoms:
            raise ValueError(f"PDB {pdb_id} has no usable ATOM/HETATM records.")
        center, size = _bbox_and_center(all_atoms, margin=0.0)
        return PocketResult(
            center_x=center[0],
            center_y=center[1],
            center_z=center[2],
            size_x=24.0,
            size_y=24.0,
            size_z=24.0,
            druggability_score=0.30,
            source="fallback",
            ligand_resname=None,
            pocket_residues=[],
            notes=(
                f"No bound ligand detected in {pdb_id}. Using protein "
                "centre-of-mass as a guess — accuracy is low. Upload a "
                "co-crystal PDB with a bound inhibitor for real pocket detection."
            ),
        )

    center, size = _bbox_and_center(chosen_coords, margin=8.0)
    residues = _pocket_residues_near(pdb_text, center, radius=6.0)

    return PocketResult(
        center_x=round(center[0], 3),
        center_y=round(center[1], 3),
        center_z=round(center[2], 3),
        size_x=round(size[0], 3),
        size_y=round(size[1], 3),
        size_z=round(size[2], 3),
        druggability_score=0.95,
        source="cocrystal_ligand",
        ligand_resname=chosen_resname,
        pocket_residues=residues,
        notes=(
            f"Pocket inferred from bound ligand {chosen_resname!r} "
            f"({len(chosen_coords)} atoms) in PDB {pdb_id}. "
            f"{len(residues)} residues lining pocket."
        ),
    )
