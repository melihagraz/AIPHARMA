"""Real cavity detection via fpocket (Voronoi tessellation / alpha spheres).

This is the scientifically-rigorous pocket-detection path used for the
manuscript pipeline. Unlike the co-crystal-ligand shortcut in ``pocket.py``
(which assumes the binding site is already known and assigns a fixed
druggability confidence), fpocket *detects* cavities de novo and reports a
computed Druggability Score per pocket (Le Guilloux et al., 2009;
Schmidtke & Barril, 2010).

fpocket is a Linux C tool with no reliable native-Windows build, so on
Windows we drive the copy installed inside WSL (Ubuntu). The service itself
still runs as a normal Windows Python process; only the fpocket call is
delegated to WSL. On Linux/macOS hosts where ``fpocket`` is on ``PATH``,
the same code runs it directly with no WSL indirection.

Pocket geometry (center + Vina bounding box) is derived from the alpha-sphere
vertices of the top pocket; pocket-lining residues reuse the 6.0 Å cutoff
shared with the rest of the SBDD stack.
"""

from __future__ import annotations

import logging
import platform
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Optional

from .pocket import (
    _bbox_and_center,
    _fetch_pdb,
    _pocket_residues_near,
)
from .schemas import PocketResult

logger = logging.getLogger(__name__)

# Where fpocket lives inside WSL when installed via the bundled Miniforge /
# bioconda path. ``command -v fpocket`` is tried first so a system-wide or
# conda-env install is picked up automatically.
_WSL_FPOCKET_FALLBACK = "$HOME/miniforge3/bin/fpocket"

# Bash that runs fpocket in an isolated work dir and streams the info file +
# every pocket's alpha-sphere vertices back on stdout. Read from stdin
# (``bash -s``) with the PDB path passed as ``$1`` so no shell-quoting or
# CRLF issues can corrupt it.
_RUN_SCRIPT = """
set -e
SRC="$1"
BASE="$(basename "$SRC")"
STEM="${BASE%.*}"
WORK="$HOME/.aipharma_fpocket/${STEM}"
rm -rf "$WORK"; mkdir -p "$WORK"
cp "$SRC" "$WORK/${STEM}.pdb"
cd "$WORK"
FPOCKET="$(command -v fpocket || echo "%FALLBACK%")"
"$FPOCKET" -f "${STEM}.pdb" >/dev/null 2>&1
echo "===INFO==="
cat "${STEM}_out/${STEM}_info.txt"
for f in "${STEM}_out"/pockets/pocket*_vert.pqr; do
  [ -f "$f" ] || continue
  echo "===VERT:$(basename "$f")==="
  cat "$f"
done
echo "===END==="
""".replace("%FALLBACK%", _WSL_FPOCKET_FALLBACK)


def _is_windows() -> bool:
    return platform.system() == "Windows"


def _wsl_prefix() -> list[str]:
    """Command prefix to run a POSIX command, via WSL on Windows."""
    return ["wsl"] if _is_windows() else []


@lru_cache(maxsize=1)
def fpocket_available() -> bool:
    """True if fpocket can be invoked (directly, or through WSL on Windows)."""
    try:
        if _is_windows():
            probe = (
                'command -v fpocket >/dev/null 2>&1 '
                f'|| test -x {_WSL_FPOCKET_FALLBACK}'
            )
            res = subprocess.run(
                ["wsl", "bash", "-lc", probe],
                capture_output=True,
                timeout=30,
            )
            return res.returncode == 0
        return shutil.which("fpocket") is not None
    except (OSError, subprocess.SubprocessError):
        return False


def _to_run_path(pdb_path: Path) -> str:
    """Path as seen by the fpocket runtime (translated to /mnt/... under WSL)."""
    if not _is_windows():
        return str(pdb_path)
    # Backslashes get eaten crossing the Windows->WSL argv boundary; wslpath
    # accepts a forward-slash Windows path just fine ("C:/a/b" -> /mnt/c/a/b).
    win = str(pdb_path).replace("\\", "/")
    res = subprocess.run(
        ["wsl", "wslpath", "-a", win],
        capture_output=True,
        text=True,
        timeout=30,
    )
    out = res.stdout.strip()
    if res.returncode != 0 or not out:
        raise RuntimeError(f"wslpath translation failed for {pdb_path}: {res.stderr}")
    return out


def _run_fpocket(pdb_path: Path) -> str:
    """Run fpocket on ``pdb_path`` and return the raw marker-delimited stdout."""
    run_path = _to_run_path(pdb_path)
    cmd = [*_wsl_prefix(), "bash", "-s", run_path]
    # Send the script as bytes with LF endings. Text mode would re-translate
    # '\n' to '\r\n' on Windows, and bash rejects the stray '\r'.
    script = _RUN_SCRIPT.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")
    try:
        res = subprocess.run(
            cmd,
            input=script,
            capture_output=True,
            timeout=300,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("fpocket timed out (>300s).") from exc
    stdout = res.stdout.decode("utf-8", "replace")
    if res.returncode != 0 or "===INFO===" not in stdout:
        stderr = res.stderr.decode("utf-8", "replace")
        raise RuntimeError(
            "fpocket run failed "
            f"(exit {res.returncode}): {stderr.strip()[:500] or 'no stderr'}"
        )
    return stdout


def _parse_info(info_text: str) -> dict[int, dict[str, float]]:
    """Parse ``*_info.txt`` into {pocket_number: {descriptor: value}}."""
    pockets: dict[int, dict[str, float]] = {}
    current: Optional[int] = None
    for line in info_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("Pocket ") and stripped.endswith(":"):
            try:
                current = int(stripped.split()[1])
                pockets[current] = {}
            except (ValueError, IndexError):
                current = None
            continue
        if current is None or ":" not in stripped:
            continue
        key, _, val = stripped.partition(":")
        key = key.strip().lower()
        try:
            value = float(val.strip())
        except ValueError:
            continue
        if "druggability" in key:
            pockets[current]["druggability"] = value
        elif key.startswith("volume :") or key == "volume":
            pockets[current]["volume"] = value
        elif "number of alpha spheres" in key:
            pockets[current]["n_alpha"] = value
    return pockets


def _parse_vertices(stdout: str) -> dict[int, list[tuple[float, float, float]]]:
    """Parse the ``===VERT:pocketN_vert.pqr===`` blocks into per-pocket coords.

    PQR is whitespace-delimited (not fixed-column like PDB), so x/y/z are read
    as the three numeric fields preceding the trailing charge + radius.
    """
    verts: dict[int, list[tuple[float, float, float]]] = {}
    current: Optional[int] = None
    for line in stdout.splitlines():
        if line.startswith("===VERT:"):
            token = line[len("===VERT:"):].rstrip("=").strip()  # pocketN_vert.pqr
            digits = "".join(c for c in token.split("_")[0] if c.isdigit())
            current = int(digits) if digits else None
            if current is not None:
                verts.setdefault(current, [])
            continue
        if line.startswith("==="):
            current = None
            continue
        if current is None or not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            x, y, z = float(parts[-5]), float(parts[-4]), float(parts[-3])
        except (ValueError, IndexError):
            continue
        verts[current].append((x, y, z))
    return verts


def pocket_from_fpocket(
    pdb_id: Optional[str] = None,
    prepared_pdb_path: Optional[str] = None,
    min_druggability: float = 0.60,
) -> PocketResult:
    """Detect the top druggable pocket with fpocket.

    The pocket with the highest Druggability Score is selected. Per the
    manuscript convention, a score ≥ ``min_druggability`` marks a druggable
    candidate; if no pocket clears the threshold the top-scoring pocket is
    still returned (with ``meets_druggability_threshold=False``).
    """
    if prepared_pdb_path:
        pdb_path = Path(prepared_pdb_path)
        if not pdb_path.is_file():
            raise ValueError(f"prepared_pdb_path does not exist: {prepared_pdb_path}")
    elif pdb_id:
        pdb_path = _fetch_pdb(pdb_id)
    else:
        raise ValueError("pocket_from_fpocket needs pdb_id or prepared_pdb_path.")

    stdout = _run_fpocket(pdb_path)
    info_part = stdout.split("===INFO===", 1)[1]
    info_text = info_part.split("===VERT:", 1)[0]
    pockets = _parse_info(info_text)
    verts = _parse_vertices(stdout)

    if not pockets:
        raise RuntimeError(f"fpocket found no pockets in {pdb_path.name}.")

    best = max(pockets, key=lambda p: pockets[p].get("druggability", -1.0))
    descr = pockets[best]
    drug = descr.get("druggability", 0.0)

    coords = verts.get(best, [])
    if not coords:
        raise RuntimeError(
            f"fpocket pocket {best} in {pdb_path.name} has no alpha-sphere vertices."
        )
    center, size = _bbox_and_center(coords, margin=4.0)

    pdb_text = pdb_path.read_text(encoding="utf-8", errors="replace")
    residues = _pocket_residues_near(pdb_text, center, radius=6.0)

    meets = drug >= min_druggability
    label = pdb_id or Path(prepared_pdb_path or "").stem
    return PocketResult(
        center_x=round(center[0], 3),
        center_y=round(center[1], 3),
        center_z=round(center[2], 3),
        size_x=round(size[0], 3),
        size_y=round(size[1], 3),
        size_z=round(size[2], 3),
        druggability_score=round(drug, 3),
        source="fpocket",
        ligand_resname=None,
        pocket_residues=residues,
        volume=round(descr["volume"], 2) if "volume" in descr else None,
        n_alpha_spheres=int(descr["n_alpha"]) if "n_alpha" in descr else None,
        meets_druggability_threshold=meets,
        notes=(
            f"fpocket detected {len(pockets)} pockets in {label}; selected "
            f"pocket {best} (druggability {drug:.3f}"
            f"{'' if meets else f', below {min_druggability:.2f} threshold'}). "
            f"{len(residues)} residues lining pocket."
        ),
    )
