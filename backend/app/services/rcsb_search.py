"""State-aware PDB selection via RCSB Search API.

Given a receptor/target name and an activity mode (``agonist`` | ``antagonist``),
queries RCSB's full-text search and ranks the results so that the binding
pocket's *conformational state* matches the intended drug class:

    * Agonist  → active-state structures (agonist-bound, G-protein complex, ...)
    * Antagonist → inactive-state structures (antagonist-bound, inverse agonist, ...)

Picking the wrong state biases the De Novo model in a misleading direction —
training on an active-state pocket while designing an antagonist (or vice
versa) teaches the model the wrong binding geometries.

The ranker combines RCSB's text-relevance score with:
    * a state-consistency bonus (does the title contain the expected keywords?)
    * a resolution penalty (lower resolution is worse)
"""

from __future__ import annotations

import asyncio
import logging
import re
from dataclasses import dataclass
from typing import Literal, Optional

import httpx

logger = logging.getLogger(__name__)

_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

ActivityMode = Literal["agonist", "antagonist"]

# Keywords used for the RCSB full-text query (search side).
_AGONIST_KEYWORDS = [
    "agonist-bound",
    "agonist bound",
    "active state",
    "active-state",
    "G-protein complex",
    "G protein complex",
    "signaling complex",
    "fully active",
]
_ANTAGONIST_KEYWORDS = [
    "antagonist-bound",
    "antagonist bound",
    "inactive state",
    "inactive-state",
    "apo inactive",
    "inverse agonist",
    "inverse-agonist",
    "antagonist complex",
]

# Regex patterns for state inference (title-based, word-boundary aware).
# IMPORTANT: antagonist patterns are checked BEFORE agonist ones so that the
# word "antagonist" is never misread as the word "agonist" (the suffix
# distinguishes them). Plain word-boundary matching also handles arbitrary
# word order ("BOUND ANTAGONIST" vs "ANTAGONIST-BOUND") which substring
# matching cannot.
_ANTAGONIST_RE = re.compile(
    r"\b(?:antagonist|inverse[\s-]?agonist|inactive(?:[\s-]?state)?|apo\s+inactive)\b",
    re.IGNORECASE,
)
_AGONIST_RE = re.compile(
    r"\b(?:agonist|active(?:[\s-]?state)?|g[\s-]?protein\s+(?:complex|coupled)|signaling\s+complex|fully\s+active)\b",
    re.IGNORECASE,
)


@dataclass
class PDBCandidate:
    """One ranked PDB suggestion."""

    pdb_id: str
    title: str
    resolution: Optional[float]
    method: Optional[str]
    inferred_state: str  # "active" | "inactive" | "unclear"
    state_consistent: bool
    name_match: float  # 0..1, fraction of target-name tokens found in title
    rcsb_score: float
    combined_score: float


# Tokens that look like name fragments but are too generic to count as a match.
_NAME_STOPWORDS = {
    "the", "and", "of", "with", "a", "an", "by", "for",
    "receptor", "protein", "complex", "human", "structure", "crystal", "cryoem", "cryo",
}


def _tokenize_target_name(name: str) -> list[str]:
    """Lowercase + split on non-word chars; drop stopwords and 1-char tokens."""
    raw = re.split(r"[\s/\-,_()]+", name.lower())
    return [w for w in raw if len(w) >= 2 and w not in _NAME_STOPWORDS]


def _name_match_score(target_tokens: list[str], title: str) -> float:
    """Fraction of target tokens that appear as substrings in the title."""
    if not target_tokens or not title:
        return 0.0
    t = title.lower()
    hits = sum(1 for tok in target_tokens if tok in t)
    return hits / len(target_tokens)


def _infer_state(title: str) -> str:
    """Pick a state label from the title via word-boundary regex matching.

    Antagonist patterns are checked first because the word "antagonist"
    contains the substring "agonist" — without the antagonist-first ordering
    a naive scanner would mislabel an antagonist-bound structure as "active".
    """
    if not title:
        return "unclear"
    if _ANTAGONIST_RE.search(title):
        return "inactive"
    if _AGONIST_RE.search(title):
        return "active"
    return "unclear"


def _build_query(target_name: str, mode: ActivityMode) -> dict:
    """Build an RCSB search query: target_name AND (any of state-keywords)."""
    keywords = _AGONIST_KEYWORDS if mode == "agonist" else _ANTAGONIST_KEYWORDS
    # First 4 keywords give a reasonable OR-of-synonyms; more would dilute relevance.
    state_nodes = [
        {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": kw},
        }
        for kw in keywords[:4]
    ]
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": target_name},
                },
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": state_nodes,
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": 20},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }


async def _fetch_entry_metadata(client: httpx.AsyncClient, pdb_id: str) -> Optional[dict]:
    url = _ENTRY_URL.format(pdb_id=pdb_id)
    try:
        response = await client.get(url, timeout=10.0)
    except Exception as exc:  # noqa: BLE001
        logger.warning("RCSB entry fetch failed for %s: %s", pdb_id, exc)
        return None
    if response.status_code != 200:
        return None
    try:
        return response.json()
    except Exception:  # noqa: BLE001
        return None


def _build_candidate(
    entry: dict,
    pdb_id: str,
    rcsb_score: float,
    mode: ActivityMode,
    target_tokens: list[str],
) -> PDBCandidate:
    title = str((entry.get("struct") or {}).get("title") or "")
    info = entry.get("rcsb_entry_info") or {}

    resolution: Optional[float] = None
    raw_res = info.get("resolution_combined")
    if isinstance(raw_res, list) and raw_res:
        try:
            resolution = float(raw_res[0])
        except (TypeError, ValueError):
            resolution = None
    elif isinstance(raw_res, (int, float)):
        resolution = float(raw_res)

    method: Optional[str] = None
    exptl = entry.get("exptl") or []
    if isinstance(exptl, list) and exptl:
        method = exptl[0].get("method")

    inferred = _infer_state(title)
    if mode == "agonist":
        state_consistent = inferred == "active"
    else:
        state_consistent = inferred == "inactive"

    # Name match: how much of the target_name (significant tokens) appears in
    # the title? Heavily weighted — wrong-receptor matches are the second
    # biggest scientific risk after wrong-state. A 0.0 match (no overlap)
    # gets -3.0; 1.0 match gets +0.0. We avoid a positive bonus here because
    # the RCSB relevance score already rewards name matches; this term only
    # *demotes* off-target entries that RCSB ranked too high.
    name_match = _name_match_score(target_tokens, title)

    combined = rcsb_score
    if state_consistent:
        combined += 1.0
    elif inferred == "unclear":
        combined -= 0.2
    else:
        combined -= 1.5  # opposite state — strongly penalize
    if resolution is not None:
        combined -= max(0.0, (resolution - 2.0) * 0.2)
    # Apply name penalty: -3.0 when none of the target tokens appear in the title,
    # scaling linearly back to 0.0 when all of them do.
    combined -= (1.0 - name_match) * 3.0

    return PDBCandidate(
        pdb_id=pdb_id,
        title=title,
        resolution=resolution,
        method=method,
        inferred_state=inferred,
        state_consistent=state_consistent,
        name_match=round(name_match, 2),
        rcsb_score=rcsb_score,
        combined_score=round(combined, 3),
    )


async def search_pdb_for_target(
    target_name: str,
    mode: ActivityMode,
    max_candidates: int = 5,
) -> list[PDBCandidate]:
    """Return ranked PDB candidates for ``target_name`` matching the activity mode.

    Best first. Empty list = no usable matches.
    """
    if not target_name or len(target_name.strip()) < 2:
        return []
    if mode not in ("agonist", "antagonist"):
        raise ValueError(f"mode must be 'agonist' or 'antagonist', got {mode!r}")

    target_tokens = _tokenize_target_name(target_name)
    query = _build_query(target_name.strip(), mode)
    async with httpx.AsyncClient(timeout=httpx.Timeout(15.0)) as client:
        try:
            response = await client.post(_SEARCH_URL, json=query)
        except httpx.RequestError as exc:
            from app.core.error_log import log_error

            log_error(
                "rcsb_search_failed",
                f"RCSB search network error for target={target_name!r} mode={mode}: {exc.__class__.__name__}",
                target_name=target_name,
                mode=mode,
                error=str(exc),
            )
            logger.warning("RCSB search request failed: %s", exc)
            return []
        # RCSB returns 204 No Content for "valid request, zero results" —
        # treat that the same as a 200 with an empty result_set so the
        # UniProt fallback below has a chance to fire. Any other non-2xx
        # is a real upstream problem.
        if response.status_code not in (200, 204):
            from app.core.error_log import log_error

            log_error(
                "rcsb_search_failed",
                f"RCSB search HTTP {response.status_code} for target={target_name!r} mode={mode}",
                target_name=target_name,
                mode=mode,
                http_status=response.status_code,
            )
            logger.info("RCSB search HTTP %s: %s", response.status_code, response.text[:200])
            return []
        if response.status_code == 204:
            data = {}
        else:
            try:
                data = response.json()
            except Exception:  # noqa: BLE001
                data = {}

        result_set = data.get("result_set") or []
        used_uniprot_fallback = False
        if not result_set:
            # RCSB full-text search misses many short / generic gene names
            # ("MET" collides with methionine, "JAK" with miscellaneous JAK*
            # constructs, "AKT" with non-protein hits). Fall back to UniProt
            # which indexes by exact gene symbol and carries an authoritative
            # PDB cross-reference list for every reviewed entry.
            from .uniprot_client import find_pdbs_for_gene

            # Pull a generous slice (40) because UniProt orders PDB
            # cross-refs roughly by deposition year — the oldest entries
            # are usually apo / fragment structures with no real co-crystal
            # ligand, which the drug-pocket validator downstream rejects.
            # Drug-co-crystal entries tend to surface later in the list, so
            # asking for 40 means the picker has a realistic shot at finding
            # one that survives the validation gate.
            uniprot_pdbs, uniprot_accession = await find_pdbs_for_gene(
                target_name.strip(), max_pdbs=40,
            )
            if not uniprot_pdbs:
                return []
            logger.info(
                "RCSB returned 0 hits for %r — using UniProt %s with %d PDB IDs",
                target_name, uniprot_accession, len(uniprot_pdbs),
            )
            used_uniprot_fallback = True
            # UniProt has no per-entry relevance score; treat them as
            # equal-priority and let the local re-ranker (state-consistency
            # + resolution penalty) sort them out. Score = 1.0 placeholder.
            # Cap at 20 to keep the metadata-fetch fan-out reasonable; the
            # picker downstream verifies each one with Pocket_App detect.
            pid_to_score = {pid: 1.0 for pid in uniprot_pdbs[:20]}
            entries = await asyncio.gather(*[
                _fetch_entry_metadata(client, pid) for pid in pid_to_score
            ])
        else:
            # Take top 10 by RCSB relevance, then fetch metadata in parallel
            # and re-rank locally with state-consistency + resolution penalty.
            top = result_set[:10]
            pid_to_score = {item["identifier"]: float(item.get("score", 0.0)) for item in top}

            entries = await asyncio.gather(*[
                _fetch_entry_metadata(client, pid) for pid in pid_to_score
            ])

    candidates: list[PDBCandidate] = []
    for pid, entry in zip(pid_to_score.keys(), entries):
        if entry is None:
            continue
        candidates.append(_build_candidate(entry, pid, pid_to_score[pid], mode, target_tokens))

    # Drop candidates where NONE of the target tokens appear in the title —
    # these are off-target hits (a different receptor with overlapping
    # state keywords) and surfacing them to the user is more harmful than
    # returning an empty list, because the alternative is silently feeding
    # an unrelated structure into the pipeline.
    #
    # UniProt fallback exception: when we resolved the target through
    # UniProt's gene-symbol index, every PDB ID in the list is already
    # guaranteed to belong to that exact receptor — but the RCSB *title*
    # may not contain the gene symbol at all (e.g. MET kinase entries
    # often read "hepatocyte growth factor receptor" with no "MET"
    # substring). Dropping by title token here would erase the entire
    # fallback set; skip the filter when fallback fired.
    if target_tokens and not used_uniprot_fallback:
        candidates = [c for c in candidates if c.name_match > 0.0]

    # Second-chance UniProt fallback — RCSB matched something but every hit
    # got dropped by the title-token filter (classic "MET" case: 200 with
    # methionine-bound apo structures, none of which actually carries the
    # gene symbol in the title). UniProt's gene-symbol index gives us the
    # canonical receptor entry instead.
    if not candidates and not used_uniprot_fallback:
        from .uniprot_client import find_pdbs_for_gene
        # _fetch_entry_metadata is defined at module scope above, no import
        # needed — the earlier ``from .rcsb_search import …`` line tripped
        # Python's UnboundLocalError because the symbol was already a local
        # binding in this function's frame.

        uniprot_pdbs, uniprot_accession = await find_pdbs_for_gene(
            target_name.strip(), max_pdbs=40,
        )
        if uniprot_pdbs:
            logger.info(
                "RCSB hits all dropped by title filter for %r — second-chance UniProt %s with %d PDB IDs",
                target_name, uniprot_accession, len(uniprot_pdbs),
            )
            pid_to_score_fb = {pid: 1.0 for pid in uniprot_pdbs[:20]}
            async with httpx.AsyncClient(timeout=httpx.Timeout(15.0)) as client:
                entries_fb = await asyncio.gather(*[
                    _fetch_entry_metadata(client, pid) for pid in pid_to_score_fb
                ])
            for pid, entry in zip(pid_to_score_fb.keys(), entries_fb):
                if entry is None:
                    continue
                # Build with empty target_tokens so the local re-ranker
                # doesn't punish entries again for the missing gene symbol.
                candidates.append(
                    _build_candidate(entry, pid, pid_to_score_fb[pid], mode, [])
                )

    candidates.sort(key=lambda c: c.combined_score, reverse=True)
    return candidates[:max_candidates]
