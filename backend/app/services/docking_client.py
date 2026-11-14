"""HTTP client for the standalone Docking_App service (AutoDock Vina)."""

from __future__ import annotations

import logging
import math
from typing import Any

import httpx
from fastapi import HTTPException

from app.core.config import get_settings
from app.core.internal_service import service_headers

logger = logging.getLogger(__name__)


async def dock_batch(
    smiles: list[str],
    target: str,
    exhaustiveness: int = 8,
    num_modes: int = 5,
    cpu_per_dock: int = 2,
    workers: int = 4,
    timeout: float = 1800.0,
) -> list[dict[str, Any]]:
    """Call Docking_App and return per-molecule docking results.

    Each item: {smiles, target, affinity (kcal/mol), error?, num_modes}
    """
    if not smiles:
        return []
    cfg = get_settings()
    base = cfg.docking_app_url.rstrip("/")
    url = f"{base}/dock"
    payload = {
        "smiles": list(smiles),
        "target": target,
        "exhaustiveness": int(exhaustiveness),
        "num_modes": int(num_modes),
        "cpu_per_dock": int(cpu_per_dock),
        "workers": int(workers),
    }
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.post(url, headers=service_headers(), json=payload)
    except httpx.ConnectError as exc:
        raise HTTPException(
            status_code=503,
            detail=(
                f"Docking_App unreachable at {url}. "
                "Start the service: cd Docking_App && uvicorn app.main:app --port 8040"
            ),
        ) from exc
    if r.status_code >= 400:
        text = r.text[:1000]
        raise HTTPException(status_code=r.status_code, detail=f"Docking_App error: {text}")
    try:
        body = r.json()
    except Exception as exc:
        raise HTTPException(status_code=502, detail="Docking service returned non-JSON results") from exc
    if not isinstance(body, dict) or not isinstance(body.get("results"), list):
        raise HTTPException(status_code=502, detail="Docking service returned invalid results schema")
    results = body["results"]
    if len(results) != len(smiles) or any(
        not isinstance(item, dict) or item.get("smiles") != expected
        for item, expected in zip(results, smiles)
    ):
        raise HTTPException(status_code=502, detail="Docking results missing or out of input SMILES order")
    if "count" in body and body["count"] != len(smiles):
        raise HTTPException(status_code=502, detail="Docking count metadata inconsistent")
    if "target" in body and body["target"].casefold() != target.casefold():
        raise HTTPException(status_code=502, detail="Docking receptor identity mismatch")
    successful = 0
    for i, item in enumerate(results):
        if item.get("target", target).casefold() != target.casefold():
            raise HTTPException(status_code=502, detail=f"Docking row target mismatch at index {i}")
        affinity = item.get("affinity")
        if affinity is None:
            if not item.get("error"):
                raise HTTPException(status_code=502, detail=f"Docking row {i} has neither score nor error")
        elif isinstance(affinity, bool) or not isinstance(affinity, (float, int)) or not math.isfinite(affinity):
            raise HTTPException(status_code=502, detail=f"Docking row {i} has invalid affinity")
        else:
            successful += 1
    if "success_count" in body and body["success_count"] != successful:
        raise HTTPException(status_code=502, detail="Docking success_count metadata inconsistent")
    return results


async def list_available_targets(timeout: float = 5.0) -> list[str]:
    """Probe Docking_App for the set of pre-bundled targets."""
    cfg = get_settings()
    base = cfg.docking_app_url.rstrip("/")
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            r = await client.get(f"{base}/targets", headers=service_headers())
            return (r.json() or {}).get("targets") or []
    except Exception:  # noqa: BLE001
        return []


async def health() -> dict[str, Any]:
    cfg = get_settings()
    base = cfg.docking_app_url.rstrip("/")
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(5.0)) as client:
            r = await client.get(f"{base}/health", headers=service_headers())
        return r.json()
    except Exception:  # noqa: BLE001
        return {"status": "unreachable", "vina_binary": False, "available_targets": []}
