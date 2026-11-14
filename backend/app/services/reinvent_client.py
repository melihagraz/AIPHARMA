"""HTTP client for the standalone Reinvent_App service.

Calls ``POST /sample`` on the Reinvent_App FastAPI server and returns
generated SMILES with NLL + RDKit basic descriptors.
"""

from __future__ import annotations

import logging
import math
from typing import Any, Optional

import httpx
from fastapi import HTTPException

from app.core.config import get_settings
from app.core.internal_service import service_headers

logger = logging.getLogger(__name__)


async def sample_smiles(
    batch_size: int = 400,
    seed: Optional[int] = None,
    validate_with_rdkit: bool = True,
    reference_smiles: Optional[list[str]] = None,
    target_prior: Optional[str] = None,
    timeout: float = 300.0,
) -> tuple[list[dict[str, Any]], str]:
    """Call Reinvent_App and return (candidate list, prior_actually_used).

    Each candidate: {smiles, nll, valid, mw, logp, qed, tpsa, max_tanimoto?}

    When ``reference_smiles`` is provided, each generated molecule gets a
    ``max_tanimoto`` (ECFP4) field measuring chemical similarity to the
    closest reference active.

    When ``target_prior`` is provided (e.g. "EGFR"), Reinvent_App will use
    the matching fine-tuned prior if one exists, else fall back to base.
    The second return value reports which one was actually used.

    Raises HTTPException on connection / non-2xx errors.
    """
    cfg = get_settings()
    base = cfg.reinvent_app_url.rstrip("/")
    url = f"{base}/sample"
    payload: dict[str, Any] = {
        "batch_size": int(batch_size),
        "seed": seed,
        "validate_with_rdkit": bool(validate_with_rdkit),
    }
    if reference_smiles:
        payload["reference_smiles"] = list(reference_smiles)
    if target_prior:
        payload["target_prior"] = str(target_prior)
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.post(url, headers=service_headers(), json=payload)
    except httpx.ConnectError as exc:
        raise HTTPException(
            status_code=503,
            detail=(
                f"Reinvent_App unreachable at {url}. "
                "Start the service: cd Reinvent_App && uvicorn app.main:app --port 8020"
            ),
        ) from exc
    if r.status_code >= 400:
        text = r.text[:1000]
        raise HTTPException(status_code=r.status_code, detail=f"Reinvent_App error: {text}")
    try:
        body = r.json()
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=502, detail=f"Reinvent_App non-JSON response: {exc}") from exc
    if not isinstance(body, dict) or not isinstance(body.get("candidates"), list):
        raise HTTPException(status_code=502, detail="REINVENT service returned invalid candidate schema")
    candidates = body["candidates"]
    if body.get("returned") != len(candidates) or body.get("requested") != batch_size:
        raise HTTPException(status_code=502, detail="REINVENT candidate counts are inconsistent")
    for index, item in enumerate(candidates):
        if not isinstance(item, dict) or not isinstance(item.get("smiles"), str) or not item["smiles"].strip():
            raise HTTPException(status_code=502, detail=f"REINVENT candidate {index} lacks valid SMILES field")
        if not isinstance(item.get("valid"), bool):
            raise HTTPException(status_code=502, detail=f"REINVENT candidate {index} lacks validity status")
        nll = item.get("nll")
        if item["valid"] and (isinstance(nll, bool) or not isinstance(nll, (int, float)) or not math.isfinite(nll)):
            raise HTTPException(status_code=502, detail=f"REINVENT candidate {index} lacks valid model likelihood")
    prior_used = body.get("prior_used")
    if not isinstance(prior_used, str) or not prior_used.strip():
        raise HTTPException(status_code=502, detail="REINVENT service did not identify the loaded prior")
    return candidates, prior_used


async def list_available_priors(timeout: float = 5.0) -> list[str]:
    """Probe Reinvent_App for the set of currently available priors."""
    cfg = get_settings()
    base = cfg.reinvent_app_url.rstrip("/")
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.get(f"{base}/priors", headers=service_headers())
            return (r.json() or {}).get("available") or []
    except Exception:  # noqa: BLE001
        return []


async def health() -> dict[str, Any]:
    cfg = get_settings()
    base = cfg.reinvent_app_url.rstrip("/")
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(5.0)) as client:
            r = await client.get(f"{base}/health", headers=service_headers())
        return r.json()
    except Exception:  # noqa: BLE001
        return {"status": "unreachable", "prior_loaded": False}
