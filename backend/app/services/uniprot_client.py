"""UniProt gene-name → PDB cross-reference lookup.

Used as a fallback for ``rcsb_search`` when RCSB's full-text search returns
zero hits for a short or ambiguous target token (e.g. ``"MET"`` — also the
amino-acid methionine code, which buries the receptor under noise).
UniProt's gene-symbol index is much more reliable for canonical receptor
names, and every reviewed (Swiss-Prot) entry carries an explicit list of
PDB cross-references so we can hand the IDs back to RCSB for metadata
enrichment without paying another text-search round trip.

Resolution chain:

    1. ``gene_exact:<TARGET> AND organism_id:9606 AND reviewed:true``
       → human, Swiss-Prot reviewed entry.
    2. If empty, drop the ``gene_exact`` constraint and try ``gene:<TARGET>``
       (loose match, picks up the receptor when the gene field uses a
       slightly different casing).
    3. If still empty, return ``[]``.
"""

from __future__ import annotations

import logging
from typing import Optional

import httpx

logger = logging.getLogger(__name__)

_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_TIMEOUT = 12.0


async def find_pdbs_for_gene(target_name: str, max_pdbs: int = 25) -> tuple[list[str], Optional[str]]:
    """Return ``(pdb_ids, uniprot_accession)`` for the canonical reviewed
    human entry of ``target_name``. Empty list when nothing matches.

    The accession (e.g. ``"P08581"``) is returned too so the caller can
    surface it in the response — useful when the user wants to verify
    that ``"MET"`` resolved to the MET kinase rather than methionine
    aminopeptidase or anything similar.
    """
    target = target_name.strip()
    if len(target) < 2:
        return [], None

    queries = [
        f"gene_exact:{target} AND organism_id:9606 AND reviewed:true",
        f"gene:{target} AND organism_id:9606 AND reviewed:true",
    ]

    async with httpx.AsyncClient(timeout=httpx.Timeout(_TIMEOUT)) as client:
        for query in queries:
            params = {
                "query": query,
                "format": "json",
                "fields": "accession,protein_name,xref_pdb",
                "size": "1",
            }
            try:
                response = await client.get(_UNIPROT_SEARCH, params=params)
            except httpx.RequestError as exc:
                logger.info("UniProt search network error for %r: %s", target, exc)
                continue
            if response.status_code != 200:
                continue
            try:
                data = response.json()
            except Exception:  # noqa: BLE001
                continue
            results = data.get("results") or []
            if not results:
                continue
            entry = results[0]
            accession = entry.get("primaryAccession")
            xrefs = entry.get("uniProtKBCrossReferences") or []
            pdb_ids: list[str] = []
            for xref in xrefs:
                if xref.get("database") == "PDB":
                    pid = xref.get("id")
                    if isinstance(pid, str) and pid:
                        pdb_ids.append(pid.upper())
                        if len(pdb_ids) >= max_pdbs:
                            break
            if pdb_ids:
                logger.info(
                    "UniProt fallback resolved %s → %s (%d PDB cross-refs)",
                    target, accession, len(pdb_ids),
                )
                return pdb_ids, accession

    return [], None
