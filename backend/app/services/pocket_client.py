"""HTTP client for the standalone Pocket_App service (Phase 3 / SBDD).

Two endpoints over the production pipeline:

    POST /detect          → druggable pocket center + bounding box
    POST /pharmacophore   → 3D pharmacophore feature list

Both calls are best-effort: if Pocket_App is unreachable, the upstream
caller (production_pipeline / agent_service) should fall back to the
existing LBDD-only flow rather than fail the whole job.
"""

from __future__ import annotations

import logging
from typing import Any, Optional

import httpx
from fastapi import HTTPException

from app.core.config import get_settings
from app.core.internal_service import service_headers

logger = logging.getLogger(__name__)


def _base_url() -> str:
    cfg = get_settings()
    base = getattr(cfg, "pocket_app_url", None) or "http://127.0.0.1:8050"
    return base.rstrip("/")


# ---------------------------------------------------------------------------
# Pocket detection
# ---------------------------------------------------------------------------


async def detect_pocket(
    target: Optional[str] = None,
    pdb_id: Optional[str] = None,
    ligand_resname: Optional[str] = None,
    timeout: float = 30.0,
) -> Optional[dict[str, Any]]:
    """Return the pocket dict ``{center_x, center_y, center_z, size_x, ..., source, druggability_score, pocket_residues, notes}``.

    Returns ``None`` if both target and pdb_id are missing or the service
    is unreachable; raises :class:`HTTPException` on a 4xx/5xx with detail.
    """
    if not target and not pdb_id:
        return None
    url = f"{_base_url()}/detect"
    payload: dict[str, Any] = {}
    if target:
        payload["target"] = str(target)
    if pdb_id:
        payload["pdb_id"] = str(pdb_id)
    if ligand_resname:
        payload["ligand_resname"] = str(ligand_resname)
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.post(url, json=payload, headers=service_headers())
    except httpx.ConnectError as exc:
        logger.warning("Pocket_App unreachable at %s: %s", url, exc)
        return None
    if r.status_code >= 400:
        text = r.text[:500]
        raise HTTPException(status_code=r.status_code, detail=f"Pocket_App: {text}")
    body = r.json()
    return body.get("pocket")


# ---------------------------------------------------------------------------
# Pharmacophore extraction
# ---------------------------------------------------------------------------


async def extract_pharmacophore(
    target: Optional[str] = None,
    pdb_id: Optional[str] = None,
    ligand_resname: Optional[str] = None,
    timeout: float = 30.0,
) -> Optional[dict[str, Any]]:
    """Return the full pharmacophore response ``{features, ligand_smiles, ...}``.

    Returns ``None`` when the service is unreachable or no source is given;
    raises :class:`HTTPException` on a non-2xx error from the service.
    """
    if not target and not pdb_id:
        return None
    url = f"{_base_url()}/pharmacophore"
    payload: dict[str, Any] = {}
    if target:
        payload["target"] = str(target)
    if pdb_id:
        payload["pdb_id"] = str(pdb_id)
    if ligand_resname:
        payload["ligand_resname"] = str(ligand_resname)
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.post(url, json=payload, headers=service_headers())
    except httpx.ConnectError as exc:
        logger.warning("Pocket_App unreachable at %s: %s", url, exc)
        return None
    if r.status_code >= 400:
        text = r.text[:500]
        raise HTTPException(status_code=r.status_code, detail=f"Pocket_App: {text}")
    return r.json()


async def redock_validate(
    pdb_id: str,
    exhaustiveness: int = 12,
    timeout: float = 700.0,
) -> Optional[dict[str, Any]]:
    """Run self-docking validation on ``pdb_id`` via Pocket_App.

    Returns ``{rmsd, reproduced, ligand, dock_score, crystal_pdb, docked_pdb, ...}``
    or ``None`` if the service is unreachable. Slow (Vina): the default timeout
    is generous. Raises :class:`HTTPException` on a service-side 4xx/5xx.
    """
    if not pdb_id:
        return None
    url = f"{_base_url()}/redock"
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.post(url, json={"pdb_id": str(pdb_id), "exhaustiveness": exhaustiveness}, headers=service_headers())
    except httpx.ConnectError as exc:
        logger.warning("Pocket_App unreachable at %s: %s", url, exc)
        return None
    if r.status_code >= 400:
        raise HTTPException(status_code=r.status_code, detail=f"Pocket_App: {r.text[:500]}")
    return r.json()


async def health() -> dict[str, Any]:
    url = f"{_base_url()}/health"
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(5.0)) as client:
            r = await client.get(url, headers=service_headers())
        return r.json()
    except Exception:  # noqa: BLE001
        return {"status": "unreachable", "rdkit_pharm3d": False, "bio_pdb": False}
