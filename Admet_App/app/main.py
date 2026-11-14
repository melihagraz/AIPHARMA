"""FastAPI app exposing ADMET-AI predictions over HTTP.

Run:
    cd Admet_App
    uvicorn app.main:app --host 127.0.0.1 --port 8030 --reload
"""

from __future__ import annotations

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
import os
import hmac
from fastapi.middleware.cors import CORSMiddleware

from .predictor import get_default_predictor
from .schemas import (
    HealthResponse,
    MoleculeAdmet,
    PredictRequest,
    PredictResponse,
)

app = FastAPI(
    title="AI Pharma — ADMET Predictor",
    description=(
        "FastAPI wrapper around ADMET-AI (Stanford, Swanson et al. 2024). "
        "Predicts 52 ADMET endpoints + DrugBank percentiles per molecule. "
        "ADMET-AI is MIT licensed (see NOTICE)."
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
    p = get_default_predictor()
    return HealthResponse(status="ok", model_loaded=p.is_loaded)


@app.post("/predict", response_model=PredictResponse)
def predict(req: PredictRequest) -> PredictResponse:
    p = get_default_predictor()
    try:
        results = p.predict(req.smiles)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from None
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"prediction failed: {e}") from None

    out = [
        MoleculeAdmet(smiles=r.smiles, raw=r.raw, percentile=r.percentile) for r in results
    ]
    return PredictResponse(count=len(out), predictions=out)
