"""FastAPI app — pocket detection + pharmacophore extraction.

Run:
    cd Pocket_App
    uvicorn app.main:app --host 127.0.0.1 --port 8050 --reload
"""

from __future__ import annotations

import asyncio
import logging

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from .security import check_internal_access, validate_prepared_pdb_path
from fastapi.middleware.cors import CORSMiddleware

from .pocket import (
    list_bundled_targets,
    pocket_from_bundled_target,
    pocket_from_pdb,
)
from .fpocket_detect import fpocket_available, pocket_from_fpocket
from .redock import run_redock
from .pharmacophore import (
    pharmacophore_from_pdb,
    pharmacophore_from_target,
)
from .schemas import (
    HealthResponse,
    PharmacophoreRequest,
    PharmacophoreResponse,
    PocketDetectRequest,
    PocketDetectResponse,
)

logger = logging.getLogger(__name__)

app = FastAPI(
    title="AI Pharma — Pocket & Pharmacophore Service",
    description=(
        "Phase 3 SBDD component. Provides druggable-pocket detection and 3D "
        "pharmacophore extraction from bound ligands. Pocket detection uses "
        "fpocket (Voronoi/alpha-sphere cavity detection with computed "
        "druggability) when available, with a pure-Python RDKit + Biopython "
        "co-crystal-ligand fallback. Pharmacophore extraction is RDKit-based."
    ),
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.middleware("http")
async def internal_service_auth(request: Request, call_next):
    denied = check_internal_access(request)
    if denied is not None:
        return denied
    return await call_next(request)


@app.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    rdkit_pharm3d_ok = False
    try:
        from rdkit.Chem import ChemicalFeatures  # noqa: F401
        from rdkit import RDConfig  # noqa: F401
        rdkit_pharm3d_ok = True
    except Exception:  # noqa: BLE001
        pass
    bio_pdb_ok = False
    try:
        import Bio.PDB  # noqa: F401
        bio_pdb_ok = True
    except Exception:  # noqa: BLE001
        pass
    return HealthResponse(
        status="ok",
        rdkit_pharm3d=rdkit_pharm3d_ok,
        bio_pdb=bio_pdb_ok,
        fpocket=fpocket_available(),
        bundled_targets=list_bundled_targets(),
    )


@app.get("/targets")
def targets() -> dict:
    return {"targets": list_bundled_targets()}


@app.post("/detect", response_model=PocketDetectResponse)
async def detect(req: PocketDetectRequest) -> PocketDetectResponse:
    """Locate a druggable pocket. See ``PocketDetectRequest`` for resolution order."""
    if not req.target and not req.pdb_id and not req.prepared_pdb_path:
        raise HTTPException(
            status_code=400, detail="Provide 'target', 'pdb_id', or approved 'prepared_pdb_path'."
        )
    try:
        local_path = validate_prepared_pdb_path(req.prepared_pdb_path)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
    if local_path is not None:
        req = req.model_copy(update={"prepared_pdb_path": local_path})
    loop = asyncio.get_event_loop()

    # 1. Bundled target
    if req.target:
        pocket = await loop.run_in_executor(None, pocket_from_bundled_target, req.target)
        if pocket is not None:
            return PocketDetectResponse(target=req.target, pocket=pocket)
        logger.info(
            "Target %r not bundled; falling back to PDB resolution", req.target
        )

    # 2. PDB ID + (optional) ligand_resname + (optional) prepared_pdb_path
    if req.pdb_id or req.prepared_pdb_path:
        # 2a. Real cavity detection with fpocket (preferred when available).
        want_fpocket = req.use_fpocket if req.use_fpocket is not None else fpocket_available()
        if want_fpocket:
            try:
                pocket = await loop.run_in_executor(
                    None,
                    lambda: pocket_from_fpocket(
                        req.pdb_id, req.prepared_pdb_path, req.min_druggability
                    ),
                )
                return PocketDetectResponse(
                    target=req.target, pdb_id=req.pdb_id, pocket=pocket
                )
            except (RuntimeError, ValueError) as exc:
                if req.use_fpocket is True:
                    # Caller explicitly demanded fpocket — surface the failure.
                    raise HTTPException(status_code=502, detail=str(exc)) from None
                logger.warning(
                    "fpocket detection failed (%s); falling back to co-crystal ligand.",
                    exc,
                )

        # 2b. Co-crystal-ligand shortcut (fallback, or when use_fpocket=False).
        try:
            pocket = await loop.run_in_executor(
                None,
                pocket_from_pdb,
                req.pdb_id,
                req.ligand_resname,
                req.prepared_pdb_path,
            )
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from None
        except RuntimeError as exc:
            raise HTTPException(status_code=502, detail=str(exc)) from None
        return PocketDetectResponse(
            target=req.target, pdb_id=req.pdb_id, pocket=pocket
        )

    raise HTTPException(
        status_code=400,
        detail=(
            f"Target {req.target!r} is not bundled and no pdb_id provided. "
            "Either pick a bundled target or supply a PDB ID."
        ),
    )


@app.post("/redock")
async def redock(req: dict) -> dict:
    """Self-docking validation: remove the co-crystal ligand, re-dock it, and
    return crystal + docked poses (PDB strings) with the pose RMSD. Slow
    (Vina, ~15-60s) — the caller should show a progress state."""
    pdb_id = (req or {}).get("pdb_id")
    if not isinstance(pdb_id, str) or not pdb_id:
        raise HTTPException(status_code=400, detail="Provide 'pdb_id'.")
    exh = (req or {}).get("exhaustiveness", 12)
    loop = asyncio.get_event_loop()
    try:
        return await loop.run_in_executor(None, lambda: run_redock(pdb_id, exh))
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
    except RuntimeError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from None


@app.post("/pharmacophore", response_model=PharmacophoreResponse)
async def pharmacophore(req: PharmacophoreRequest) -> PharmacophoreResponse:
    """Extract a 3D pharmacophore from a co-crystal ligand.

    Either provide a bundled ``target`` (we'll use a known co-crystal PDB)
    or a ``pdb_id`` + optional ``ligand_resname`` directly.
    """
    if not req.target and not req.pdb_id and not req.prepared_pdb_path:
        raise HTTPException(
            status_code=400, detail="Provide 'target', 'pdb_id', or approved 'prepared_pdb_path'."
        )
    try:
        local_path = validate_prepared_pdb_path(req.prepared_pdb_path)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
    if local_path is not None:
        req = req.model_copy(update={"prepared_pdb_path": local_path})
    loop = asyncio.get_event_loop()

    if req.target:
        try:
            features, pdb_id, ligand_resname, smiles = await loop.run_in_executor(
                None, pharmacophore_from_target, req.target
            )
        except ValueError as exc:
            # Target not in bundled cocrystal map — fall through if pdb_id given
            if not req.pdb_id:
                raise HTTPException(status_code=400, detail=str(exc)) from None
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(status_code=500, detail=f"Pharmacophore extraction failed: {exc}") from None
        else:
            return PharmacophoreResponse(
                target=req.target,
                pdb_id=pdb_id,
                ligand_resname=ligand_resname,
                ligand_smiles=smiles,
                feature_count=len(features),
                features=features,
                source="cocrystal",
                notes=(
                    f"Extracted {len(features)} 3D features from {ligand_resname!r} "
                    f"in co-crystal PDB {pdb_id} for target {req.target}."
                ),
            )

    # Direct PDB path
    try:
        features, ligand_resname, smiles = await loop.run_in_executor(
            None,
            pharmacophore_from_pdb,
            req.pdb_id,
            req.ligand_resname,
            req.prepared_pdb_path,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from None
    except RuntimeError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from None
    except Exception as exc:  # noqa: BLE001
        logger.exception("Pharmacophore extraction failed")
        raise HTTPException(status_code=500, detail=f"Pharmacophore extraction failed: {exc}") from None

    return PharmacophoreResponse(
        target=req.target,
        pdb_id=req.pdb_id,
        ligand_resname=ligand_resname,
        ligand_smiles=smiles,
        feature_count=len(features),
        features=features,
        source="cocrystal",
        notes=(
            f"Extracted {len(features)} 3D features from {ligand_resname!r} in PDB {req.pdb_id}."
        ),
    )
