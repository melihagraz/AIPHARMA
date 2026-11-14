"""Thin proxy routes for the Pocket_App service (Phase 3 / SBDD).

Exposes synchronous endpoints under ``/api/pocket/*`` that the New Job UI
can call directly. Pocket detection and pharmacophore extraction are fast
(< 2 s for bundled targets, < 10 s for PDB fetch + parse) so we run them
synchronously rather than queuing a job — the user gets the result inline.

Routes:
    GET  /api/pocket/health             — service liveness + bundled target list
    POST /api/pocket/detect             — locate druggable pocket
    POST /api/pocket/pharmacophore      — extract 3D pharmacophore features
    POST /api/pocket/best-pdb-for-drug  — receptor name → ranked drug-ready PDBs
                                          (Smart PDB Picker; ~7-12 s wall clock)
"""

from __future__ import annotations

import logging
from typing import Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from app.services import pocket_client

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/pocket", tags=["pocket"])


# ---------------------------------------------------------------------------
# Schemas (lightweight; the full structured response lives in Pocket_App)
# ---------------------------------------------------------------------------


class PocketDetectIn(BaseModel):
    target: Optional[str] = None
    pdb_id: Optional[str] = None
    ligand_resname: Optional[str] = None


class PharmacophoreIn(BaseModel):
    target: Optional[str] = None
    pdb_id: Optional[str] = None
    ligand_resname: Optional[str] = None


class RedockIn(BaseModel):
    pdb_id: str
    exhaustiveness: int = 12


class BestPdbIn(BaseModel):
    """Body for /best-pdb-for-drug.

    target_name is the free-form receptor name the chemist typed in Step 1
    (e.g. "EGFR", "ADRB2"). mode is the optional drug class — when
    "agonist"/"antagonist" is supplied we pre-filter for the matching
    conformational state, which matters for GPCRs but is irrelevant for
    kinases.
    """
    target_name: str
    mode: Optional[str] = None  # "agonist" | "antagonist" | None
    max_candidates_to_check: int = 8
    max_results: int = 3


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------


@router.get("/health")
async def health() -> dict:
    return await pocket_client.health()


@router.post("/detect")
async def detect(body: PocketDetectIn) -> dict:
    if not body.target and not body.pdb_id:
        raise HTTPException(status_code=400, detail="Provide 'target' or 'pdb_id'.")
    try:
        pocket = await pocket_client.detect_pocket(
            target=body.target,
            pdb_id=body.pdb_id,
            ligand_resname=body.ligand_resname,
        )
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        logger.exception("Pocket detection failed")
        raise HTTPException(status_code=500, detail=f"detect failed: {exc}") from None
    if pocket is None:
        raise HTTPException(
            status_code=503,
            detail=(
                "Pocket_App unreachable on port 8050. Start with: "
                "cd Pocket_App && uvicorn app.main:app --port 8050"
            ),
        )
    return {"target": body.target, "pdb_id": body.pdb_id, "pocket": pocket}


@router.post("/pharmacophore")
async def pharmacophore(body: PharmacophoreIn) -> dict:
    if not body.target and not body.pdb_id:
        raise HTTPException(status_code=400, detail="Provide 'target' or 'pdb_id'.")
    try:
        resp = await pocket_client.extract_pharmacophore(
            target=body.target,
            pdb_id=body.pdb_id,
            ligand_resname=body.ligand_resname,
        )
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        logger.exception("Pharmacophore extraction failed")
        raise HTTPException(status_code=500, detail=f"pharmacophore failed: {exc}") from None
    if resp is None:
        raise HTTPException(
            status_code=503,
            detail=(
                "Pocket_App unreachable on port 8050. Start with: "
                "cd Pocket_App && uvicorn app.main:app --port 8050"
            ),
        )
    return resp


@router.post("/redock")
async def redock(body: RedockIn) -> dict:
    """Self-docking validation: re-dock the co-crystal ligand and return the
    crystal + docked poses (PDB strings) plus the pose RMSD. Slow (Vina)."""
    try:
        resp = await pocket_client.redock_validate(body.pdb_id, body.exhaustiveness)
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        logger.exception("Redocking failed")
        raise HTTPException(status_code=500, detail=f"redock failed: {exc}") from None
    if resp is None:
        raise HTTPException(status_code=503, detail="Pocket_App unreachable on port 8050.")
    return resp


@router.post("/best-pdb-for-drug")
async def best_pdb_for_drug(body: BestPdbIn) -> dict:
    """End-to-end ranker that picks the most drug-discovery-ready PDB
    for a free-form receptor name. Calls rcsb_search to enumerate
    candidates, then runs Pocket_App's drug-pocket detection on each so
    only structures with a real co-crystal ligand (not waters, ions,
    buffers, or cofactors) survive the filter.

    Wall clock: ~0.3 s per candidate × max_candidates_to_check ≈ 3 s
    minimum, more typically 7-12 s when RCSB pagination misses on the
    first page and the UniProt fallback kicks in.
    """
    from app.services.best_pdb_picker import pick_best_pdb

    name = (body.target_name or "").strip()
    if len(name) < 2:
        raise HTTPException(
            status_code=400,
            detail="target_name is required (≥2 chars).",
        )
    mode = (body.mode or "").strip().lower() or None
    if mode and mode not in {"agonist", "antagonist"}:
        raise HTTPException(
            status_code=400,
            detail=f"mode must be 'agonist' or 'antagonist' (got {body.mode!r}).",
        )
    try:
        result = await pick_best_pdb(
            target_name=name,
            mode=mode,  # type: ignore[arg-type]
            max_candidates_to_check=max(1, min(int(body.max_candidates_to_check), 30)),
            max_results=max(1, min(int(body.max_results), 10)),
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        logger.exception("best_pdb_for_drug failed")
        raise HTTPException(status_code=500, detail=f"picker failed: {exc}") from None
    return result
