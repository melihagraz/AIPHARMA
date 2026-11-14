"""FastAPI app exposing AutoDock Vina docking over HTTP.

Run:
    cd Docking_App
    uvicorn app.main:app --host 127.0.0.1 --port 8040 --reload
"""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from importlib.util import find_spec

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
import os
import hmac
from fastapi.middleware.cors import CORSMiddleware

from .dock import dock_batch, list_targets, _vina_binary
from .schemas import (
    DockRequest,
    DockResponse,
    DockResultOut,
    HealthResponse,
    TargetsResponse,
)

logger = logging.getLogger(__name__)

app = FastAPI(
    title="AI Pharma — Docking Service",
    description=(
        "AutoDock Vina docking against pre-bundled receptor PDBQT files. "
        "Targets: EGFR, JAK2, BACE1, KDR, MET, AKT1, LCK, KIT (extendable). "
        "Vina is Apache 2.0; receptors from dockstring (Apache 2.0). See NOTICE."
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
    """Deny direct unauthenticated scientific microservice calls in production.

    Service endpoints must ALSO be restricted to internal networks. This
    secret is for service-to-service requests, not a browser API token.
    """
    environment = os.getenv("APP_ENV", "development").lower().strip()
    if environment not in {"development", "test", "production"}:
        return JSONResponse(status_code=503, content={"detail": "Invalid APP_ENV configuration"})
    if environment == "production":
        secret = os.getenv("INTERNAL_SERVICE_KEY", "")
        if len(secret) < 32:
            return JSONResponse(status_code=503, content={"detail":"Internal service key not configured"})
        supplied = request.headers.get("X-Internal-Service-Key", "")
        if not hmac.compare_digest(supplied, secret):
            return JSONResponse(status_code=401, content={"detail":"Internal service authentication required"})
    return await call_next(request)


@app.get("/health", response_model=HealthResponse)
def health() -> HealthResponse:
    try:
        binary_ok = Path(_vina_binary()).is_file()
    except FileNotFoundError:
        binary_ok = False
    meeko_ok = find_spec("meeko") is not None
    targets = list_targets()
    return HealthResponse(
        status="ok" if binary_ok and meeko_ok else "dependencies_missing",
        vina_binary=binary_ok,
        meeko_available=meeko_ok,
        ready_for_docking=binary_ok and meeko_ok and bool(targets),
        available_targets=targets,
    )


@app.get("/targets", response_model=TargetsResponse)
def targets() -> TargetsResponse:
    return TargetsResponse(targets=list_targets())


@app.post("/dock", response_model=DockResponse)
async def dock(req: DockRequest) -> DockResponse:
    """Dock a batch of SMILES against a target. Runs Vina in a thread pool
    so the FastAPI event loop is not blocked while the C++ binary works."""
    # Assemble the optional pocket-box override — Pydantic gave us six
    # individual fields, but dock_batch expects either a 6-tuple or None.
    box_fields = (req.center_x, req.center_y, req.center_z, req.size_x, req.size_y, req.size_z)
    pocket_box = box_fields if all(v is not None for v in box_fields) else None

    # Custom receptor: prep PDBQT on the fly via meeko, cache by content
    # hash. If prep fails OR the caller didn't supply pocket_box, we surface
    # a clear 400 rather than silently falling back to the bundled receptor
    # (which would dock against the WRONG protein and quietly mislead the
    # user). pocket_box is required because custom receptors don't carry a
    # bundled conf.txt search grid.
    custom_receptor_pdbqt = None
    custom_receptor_hash = None
    if req.custom_receptor_pdb:
        if pocket_box is None:
            raise HTTPException(
                status_code=400,
                detail=(
                    "custom_receptor_pdb requires pocket_box "
                    "(center_x/y/z + size_x/y/z) — no bundled grid to fall back on."
                ),
            )
        from .receptor_prep import prepare_custom_receptor
        try:
            custom_receptor_pdbqt, custom_receptor_hash = prepare_custom_receptor(
                req.custom_receptor_pdb, pdb_id=req.target,
            )
        except ValueError as exc:
            raise HTTPException(
                status_code=400,
                detail=f"custom receptor prep failed: {exc}",
            ) from exc

    try:
        # dock_batch is sync + CPU-bound; offload to default executor
        loop = asyncio.get_event_loop()
        results = await loop.run_in_executor(
            None,
            lambda: dock_batch(
                smiles_list=req.smiles,
                target=req.target,
                exhaustiveness=req.exhaustiveness,
                num_modes=req.num_modes,
                cpu_per_dock=req.cpu_per_dock,
                workers=req.workers,
                pocket_box=pocket_box,
                custom_receptor_pdbqt=custom_receptor_pdbqt,
            ),
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from None
    except Exception as e:  # noqa: BLE001
        logger.exception("dock_batch failed")
        raise HTTPException(status_code=500, detail=f"docking failed: {e}") from None

    out = [
        DockResultOut(
            smiles=r.smiles,
            target=r.target,
            affinity=r.affinity,
            error=r.error,
            num_modes=r.num_modes,
            pose_pdb=r.pose_pdb,
            pose_pdb_all=r.pose_pdb_all,
        )
        for r in results
    ]
    success_count = sum(1 for r in out if r.affinity is not None)
    return DockResponse(target=req.target, count=len(out), success_count=success_count, results=out)
