"""Redocking validation driver — runs the self-docking pipeline and returns
the crystal + docked poses (as PDB strings) plus the RMSD, for the UI to
overlay in the 3D viewer.

The heavy toolchain (Meeko + AutoDock Vina + spyRMSD) lives in the WSL conda
env ``aipharma``; on Windows this shells out to it exactly like
``fpocket_detect``. ``scripts/redock_validate.py`` writes the JSON summary and
the two pose PDBs into ``reports/``, which this module then reads back.
"""

from __future__ import annotations

import json
import logging
import platform
import re
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[2]
REPORTS = REPO / "reports"
_PDB_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")


def _wsl_repo() -> str:
    win = str(REPO).replace("\\", "/")
    res = subprocess.run(["wsl", "wslpath", "-a", win],
                         capture_output=True, text=True, timeout=30)
    out = res.stdout.strip()
    if res.returncode != 0 or not out:
        raise RuntimeError(f"wslpath failed for {REPO}: {res.stderr}")
    return out


def run_redock(pdb_id: str, exhaustiveness: int = 12) -> dict:
    """Self-dock the co-crystal ligand of ``pdb_id`` and return poses + RMSD."""
    pid = (pdb_id or "").strip().upper()
    if not _PDB_RE.match(pid):
        raise ValueError(f"Invalid PDB id: {pdb_id!r}")
    exh = max(1, min(32, int(exhaustiveness)))

    if platform.system() == "Windows":
        repo = _wsl_repo()
        inner = (
            f'cd "{repo}" && ~/miniforge3/bin/mamba run -n aipharma '
            f"python scripts/redock_validate.py --pdb {pid} --exhaustiveness {exh}"
        )
        cmd = ["wsl", "bash", "-lc", inner]
    else:
        cmd = ["python", "scripts/redock_validate.py", "--pdb", pid,
               "--exhaustiveness", str(exh)]

    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=str(REPO))
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("Redocking timed out (>600s).") from exc
    if res.returncode != 0:
        raise RuntimeError(f"Redocking failed: {res.stderr.strip()[:400] or 'no stderr'}")

    summary_path = REPORTS / f"redock_{pid}.json"
    crystal_path = REPORTS / f"redock_{pid}_crystal.pdb"
    docked_path = REPORTS / f"redock_{pid}_docked.pdb"
    if not summary_path.is_file():
        raise RuntimeError(f"Redocking produced no summary for {pid}.")

    summary = json.loads(summary_path.read_text())
    return {
        **summary,
        "crystal_pdb": crystal_path.read_text() if crystal_path.is_file() else None,
        "docked_pdb": docked_path.read_text() if docked_path.is_file() else None,
    }
