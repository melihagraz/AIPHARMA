"""HTTP client for the standalone Admet_App service.

Calls ``POST /predict`` and returns per-molecule ADMET predictions
(52 raw endpoints + 52 DrugBank-approved percentiles).
"""

from __future__ import annotations

import logging
import math
from typing import Any

import httpx
from fastapi import HTTPException

from app.core.config import get_settings
from app.core.internal_service import service_headers

logger = logging.getLogger(__name__)


async def predict_admet(smiles: list[str], timeout: float = 600.0) -> list[dict[str, Any]]:
    """Call Admet_App and return predictions aligned with the input order.

    Each item: {smiles, raw: {...52 keys...}, percentile: {...52 keys...}}
    """
    if not smiles:
        return []
    cfg = get_settings()
    base = cfg.admet_app_url.rstrip("/")
    url = f"{base}/predict"
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.post(url, headers=service_headers(), json={"smiles": list(smiles)})
    except httpx.ConnectError as exc:
        raise HTTPException(
            status_code=503,
            detail=(
                f"Admet_App unreachable at {url}. "
                "Start the service: cd Admet_App && uvicorn app.main:app --port 8030"
            ),
        ) from exc
    if r.status_code >= 400:
        text = r.text[:1000]
        raise HTTPException(status_code=r.status_code, detail=f"Admet_App error: {text}")
    try:
        body = r.json()
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=502, detail=f"Admet_App non-JSON response: {exc}") from exc
    if not isinstance(body, dict) or not isinstance(body.get("predictions"), list):
        raise HTTPException(status_code=502, detail="ADMET service returned an invalid prediction schema")
    predictions = body["predictions"]
    if len(predictions) != len(smiles):
        raise HTTPException(status_code=502, detail="ADMET prediction count differs from requested molecules")
    if "count" in body and body["count"] != len(smiles):
        raise HTTPException(status_code=502, detail="ADMET service count field is inconsistent")
    for i, (requested, item) in enumerate(zip(smiles, predictions)):
        if not isinstance(item, dict) or item.get("smiles") != requested:
            raise HTTPException(status_code=502, detail=f"ADMET molecule order/identity mismatch at index {i}")
        raw = item.get("raw")
        if not isinstance(raw, dict) or not raw:
            raise HTTPException(status_code=502, detail=f"ADMET raw predictions are missing for index {i}")
        if all(val is None for val in raw.values()):
            raise HTTPException(status_code=502, detail=f"ADMET raw endpoints all unavailable at index {i}")
        for name, val in raw.items():
            if val is not None and (not isinstance(val, (float, int)) or not math.isfinite(val)):
                raise HTTPException(status_code=502, detail=f"Invalid ADMET value for {name!r} at index {i}")
    return predictions


async def health() -> dict[str, Any]:
    cfg = get_settings()
    base = cfg.admet_app_url.rstrip("/")
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(5.0)) as client:
            r = await client.get(f"{base}/health", headers=service_headers())
        return r.json()
    except Exception:  # noqa: BLE001
        return {"status": "unreachable", "model_loaded": False}
