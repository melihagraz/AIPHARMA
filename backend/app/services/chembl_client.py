"""ChEMBL REST client for fetching known active SMILES per target.

Used by the production discovery pipeline to assemble a target-relevance
reference set: when the user types "EGFR", we pull a small batch of
high-potency known inhibitors from ChEMBL and pass them to the Reinvent_App
sampler as ``reference_smiles``. Generated molecules are then ranked partly
by Tanimoto similarity to that reference set, making the system genuinely
target-aware (vs. the generic prior alone).

Resolution flow (all read-only HTTPS, no auth):
    1. /target/search?q=NAME   → top match → target_chembl_id (e.g. CHEMBL203)
    2. /molecule?molecule_chembl_id__in=...  is too coarse; instead we use
       /activity?target_chembl_id=...&pchembl_value__gte=N  to filter for
       potent activities, then pull the SMILES per molecule_chembl_id.

Cached in-memory by target name to avoid hammering ChEMBL across requests.
"""

from __future__ import annotations

import asyncio
import logging
from threading import Lock
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
DEFAULT_PCHEMBL_MIN = 7.0  # IC50 ≤ 100 nM ish — "potent" threshold
DEFAULT_LIMIT = 50

# Module-level cache: { (target_name_lower, pchembl_min, limit): [SMILES, ...] }
_cache: dict[tuple[str, float, int], list[str]] = {}
_cache_lock = Lock()


async def _resolve_target_id(client: httpx.AsyncClient, target_name: str) -> Optional[str]:
    """Find the best-matching ChEMBL target_chembl_id for a free-text name.

    Strategy: prefer human SINGLE PROTEIN targets (cleanest activity binding),
    then human PROTEIN COMPLEX targets, before considering other organisms.
    This avoids silently selecting the mouse TGFBR1 record for ALK5/TGFBR1.
    """
    url = f"{CHEMBL_BASE}/target/search.json"
    r = await client.get(url, params={"q": target_name, "limit": 25})
    r.raise_for_status()
    body = r.json()
    targets = body.get("targets") or []
    if not targets:
        return None
    # Organism is scientifically material: search results for TGFBR1 currently
    # list mouse before human even though the query is otherwise identical.
    for preferred_organism in ("Homo sapiens", None):
        for preferred_type in ("SINGLE PROTEIN", "PROTEIN COMPLEX", None):
            for t in targets:
                organism_ok = (
                    preferred_organism is None
                    or t.get("organism") == preferred_organism
                )
                type_ok = preferred_type is None or t.get("target_type") == preferred_type
                if organism_ok and type_ok:
                    tid = t.get("target_chembl_id")
                    if tid:
                        return tid
    return None


async def _activities_for_target(
    client: httpx.AsyncClient,
    target_chembl_id: str,
    pchembl_min: float,
    limit: int,
) -> list[str]:
    """Fetch potent activities and return the deduplicated SMILES set."""
    url = f"{CHEMBL_BASE}/activity.json"
    params = {
        "target_chembl_id": target_chembl_id,
        "pchembl_value__gte": pchembl_min,
        "limit": limit * 4,  # over-fetch to compensate for missing SMILES
    }
    r = await client.get(url, params=params)
    r.raise_for_status()
    body = r.json()
    activities = body.get("activities") or []
    seen: set[str] = set()
    out: list[str] = []
    for act in activities:
        smi = (act.get("canonical_smiles") or "").strip()
        if not smi or smi in seen:
            continue
        seen.add(smi)
        out.append(smi)
        if len(out) >= limit:
            break
    return out


async def fetch_known_actives(
    target_name: str,
    pchembl_min: float = DEFAULT_PCHEMBL_MIN,
    limit: int = DEFAULT_LIMIT,
    timeout: float = 30.0,
) -> list[str]:
    """Return up to ``limit`` deduplicated SMILES of known potent actives.

    Returns an empty list if the target name does not resolve in ChEMBL or
    no activities meet the pchembl threshold. Never raises on network /
    HTTP errors — it is a best-effort enrichment, not a hard dependency.
    """
    name = (target_name or "").strip()
    if not name:
        return []
    cache_key = (name.lower(), pchembl_min, limit)
    with _cache_lock:
        if cache_key in _cache:
            return list(_cache[cache_key])

    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(timeout)) as client:
            tid = await _resolve_target_id(client, name)
            if not tid:
                logger.info("ChEMBL: no target match for %r", name)
                with _cache_lock:
                    _cache[cache_key] = []
                return []
            smiles = await _activities_for_target(client, tid, pchembl_min, limit)
            logger.info(
                "ChEMBL: target=%r → %s, %d known actives @ pChEMBL>=%.1f",
                name, tid, len(smiles), pchembl_min,
            )
            with _cache_lock:
                _cache[cache_key] = list(smiles)
            return smiles
    except (httpx.HTTPError, asyncio.TimeoutError) as exc:
        logger.warning("ChEMBL lookup failed for %r: %s", name, exc)
        return []


def clear_cache() -> None:
    with _cache_lock:
        _cache.clear()
