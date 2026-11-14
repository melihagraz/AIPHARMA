"""FastAPI app exposing the REINVENT 4 prior as an HTTP service.

Run:
    cd Reinvent_App
    uvicorn app.main:app --host 127.0.0.1 --port 8020 --reload
"""

from __future__ import annotations

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
import os
import hmac
from fastapi.middleware.cors import CORSMiddleware

from .sampler import get_default_sampler, list_available_priors, resolve_prior_path
from .schemas import (
    HealthResponse,
    MoleculeOut,
    PriorsResponse,
    SampleRequest,
    SampleResponse,
)

app = FastAPI(
    title="AI Pharma — Reinvent Generator",
    description=(
        "FastAPI wrapper around the REINVENT 4 SMILES prior. "
        "The model architecture and prior weights are licensed under "
        "Apache 2.0 by AstraZeneca (see NOTICE)."
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
    s = get_default_sampler()
    if s.is_loaded:
        return HealthResponse(
            status="ok",
            prior_loaded=True,
            vocab_size=s.vocab_size,
            max_sequence_length=s.max_sequence_length,
        )
    return HealthResponse(status="ok", prior_loaded=False)


@app.get("/priors", response_model=PriorsResponse)
def priors() -> PriorsResponse:
    """List available priors (base + any fine-tuned target models)."""
    return PriorsResponse(available=sorted(list_available_priors().keys()))


@app.post("/sample", response_model=SampleResponse)
def sample(req: SampleRequest) -> SampleResponse:
    s = get_default_sampler()
    # Resolve which prior is actually used so the response is honest about
    # whether we hit a fine-tuned target or fell back to the base.
    available = list_available_priors()
    requested = (req.target_prior or "").strip()
    if not requested or requested.lower() in {"base", "default", "generic"}:
        prior_used = "base"
    else:
        match = next((k for k in available if k.lower() == requested.lower()), None)
        prior_used = match if match else "base"
    try:
        results = s.sample(
            batch_size=req.batch_size,
            seed=req.seed,
            target_prior=req.target_prior,
        )
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from None
    except Exception as e:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"sampling failed: {e}") from None

    # Pre-compute reference fingerprints once if a target-relevance set was provided.
    reference_fps = None
    if req.reference_smiles:
        from .scoring import build_reference_fingerprints

        reference_fps = build_reference_fingerprints(req.reference_smiles)
        if not reference_fps:
            # Provided list parsed to zero valid fingerprints — keep going but flag.
            reference_fps = None

    candidates: list[MoleculeOut] = []
    valid_count = 0

    if req.validate_with_rdkit:
        # Lazy import: keep RDKit out of the import path until needed.
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED
        from rdkit import RDLogger

        from .scoring import max_tanimoto_to_references

        RDLogger.DisableLog("rdApp.*")  # silence parser warnings on bad SMILES

        for r in results:
            mol = Chem.MolFromSmiles(r.smiles)
            valid = mol is not None
            if valid:
                valid_count += 1
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                qed = QED.qed(mol)
                tpsa = Descriptors.TPSA(mol)
                tani = (
                    max_tanimoto_to_references(r.smiles, reference_fps)
                    if reference_fps
                    else None
                )
                candidates.append(
                    MoleculeOut(
                        smiles=r.smiles,
                        nll=r.nll,
                        valid=True,
                        mw=round(mw, 3),
                        logp=round(logp, 3),
                        qed=round(qed, 3),
                        tpsa=round(tpsa, 3),
                        max_tanimoto=round(tani, 4) if tani is not None else None,
                    )
                )
            else:
                candidates.append(MoleculeOut(smiles=r.smiles, nll=r.nll, valid=False))
    else:
        from .scoring import max_tanimoto_to_references

        for r in results:
            tani = (
                max_tanimoto_to_references(r.smiles, reference_fps)
                if reference_fps
                else None
            )
            candidates.append(
                MoleculeOut(
                    smiles=r.smiles,
                    nll=r.nll,
                    max_tanimoto=round(tani, 4) if tani is not None else None,
                )
            )

    return SampleResponse(
        requested=req.batch_size,
        returned=len(candidates),
        valid_count=valid_count if req.validate_with_rdkit else None,
        candidates=candidates,
        prior_vocab_size=s.vocab_size,
        prior_max_seq_len=s.max_sequence_length,
        prior_used=prior_used,
    )
