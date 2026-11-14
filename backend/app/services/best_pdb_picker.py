"""Receptor name → best PDB for drug discovery.

A target's RCSB hits can include dozens of crystal structures, but only a
fraction are useful for designing a drug:

    * The crystal must contain an actual **drug-like co-crystal ligand**,
      not just a buffer (sulfate, glycerol) or cofactor (NAG, NAD).
    * Pocket detection must succeed — Pocket_App's ``source: cocrystal_ligand``,
      not ``fallback``.
    * Resolution should be drug-design-grade (X-ray ≤ 2.5 Å or Cryo-EM ≤ 3.5 Å
      preferred; we don't hard-reject, just down-rank).
    * State must match the intended drug class (agonist → active, antagonist
      → inactive). Only when the caller specifies ``mode``.

This picker stitches together ``rcsb_search`` and ``pocket_client.detect`` to
do an end-to-end check on each candidate before ranking. Slower than a pure
text search (~0.3 s per candidate), but the answer is real-world useful:
"this PDB has a real drug pocket you can dock against" rather than "this
PDB shares words with your query".

Filtering rules:
    * Drop candidates whose detected ligand is in ``_NON_DRUG_HETERO``
      (waters, ions, common buffers/cryoprotectants).
    * Drop candidates with < ``_MIN_LIGAND_HEAVY_ATOMS`` atoms in the pocket
      (rules out small fragments that aren't drug-like).
    * Drop ``source == "fallback"`` — Pocket_App couldn't even find a ligand.

Ranking score (higher is better):
    + druggability * 1.0       (Pocket_App's own confidence)
    + name_match * 0.5         (title vs target_name overlap)
    + state_consistent * 0.4   (when mode is given)
    + min(0, 2.5 - resolution) * 0.2   (resolution bonus)
"""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass
from typing import Literal, Optional

from app.services import pocket_client
from app.services.rcsb_search import search_pdb_for_target

logger = logging.getLogger(__name__)

#: Heterogens that are not drug ligands (waters, ions, buffer molecules,
#: cryoprotectants, cofactors). Mirrors protein_prep's non-ligand set plus
#: a few commonly mis-identified small molecules.
_NON_DRUG_HETERO = {
    # waters / heavy water
    "HOH", "WAT", "DOD",
    # monatomic ions
    "NA", "K", "MG", "CA", "ZN", "FE", "MN", "CU", "CO", "NI",
    "CL", "BR", "I", "F",
    # buffer / cryoprotectant / salt
    "SO4", "PO4", "NO3", "ACT", "EDO", "GOL", "PEG", "MPD",
    "DMS", "TRS", "BME", "BCT", "FMT", "OXL",
    # crystallographic additives / sugars
    "NAG", "MAN", "BMA", "FUC", "GAL",  # common N-glycans (not drug ligands)
    # generic placeholders
    "UNL", "UNK",
    # phosphate / nucleotide cofactors (debatable for drug discovery —
    # included because they almost never mark a drug pocket)
    "STP", "GTP", "GDP", "ATP", "ADP", "AMP", "CMP",
    "FMN", "FAD", "COA",
}

#: Below this heavy-atom count the "ligand" is more likely a fragment or
#: a buffer chunk than a drug.
_MIN_LIGAND_HEAVY_ATOMS = 15


@dataclass
class DrugPdbCandidate:
    pdb_id: str
    title: str
    resolution: Optional[float]
    method: Optional[str]
    ligand_resname: Optional[str]
    pocket_residue_count: int
    druggability_score: float
    pocket_source: str
    state_inferred: str
    state_consistent: bool
    name_match: float
    rcsb_score: float
    drug_pocket_score: float
    rejection_reason: Optional[str]  # populated for filtered-out candidates
    notes: list[str]


async def _verify_candidate(
    pdb_id: str,
    pdb_metadata: dict,
    mode: Optional[str],
) -> DrugPdbCandidate:
    """Run pocket detection on a single PDB and assemble a DrugPdbCandidate."""
    notes: list[str] = []
    rejection: Optional[str] = None
    ligand_resname: Optional[str] = None
    pocket_source: str = "unknown"
    druggability: float = 0.0
    residue_count: int = 0
    ligand_atom_count: int = 0

    try:
        pocket = await pocket_client.detect_pocket(pdb_id=pdb_id)
    except Exception as exc:  # noqa: BLE001
        logger.info("Pocket detect failed for %s during best-pdb scan: %s", pdb_id, exc)
        pocket = None

    if pocket is None:
        rejection = "Pocket_App unreachable or detection failed."
    else:
        pocket_source = str(pocket.get("source") or "unknown")
        ligand_resname = pocket.get("ligand_resname")
        druggability = float(pocket.get("druggability_score") or 0.0)
        residue_count = len(pocket.get("pocket_residues") or [])
        ligand_atom_count = 0
        raw_notes = pocket.get("notes") or ""
        # Pocket_App embeds the atom count in its notes string like "ligand 'AQ4' (29 atoms)"
        import re as _re
        m = _re.search(r"\((\d+)\s*atoms?\)", raw_notes)
        if m:
            try:
                ligand_atom_count = int(m.group(1))
            except ValueError:
                pass

        # Filtering rules
        if pocket_source == "fallback":
            rejection = "Pocket_App fell back to protein centre of mass — no co-crystal ligand."
        elif ligand_resname and ligand_resname.upper() in _NON_DRUG_HETERO:
            rejection = f"Bound HETATM {ligand_resname!r} is a buffer/cofactor, not a drug."
        elif ligand_atom_count and ligand_atom_count < _MIN_LIGAND_HEAVY_ATOMS:
            rejection = (
                f"Bound ligand {ligand_resname!r} has only {ligand_atom_count} atoms — "
                f"likely a fragment, not a drug."
            )

    # Score the candidate (even if rejected; useful for reporting)
    drug_score = druggability * 1.0
    if pdb_metadata.get("name_match"):
        drug_score += float(pdb_metadata["name_match"]) * 0.5
    if pdb_metadata.get("state_consistent"):
        drug_score += 0.4
    res = pdb_metadata.get("resolution")
    if isinstance(res, (int, float)) and res > 0:
        # Up to +0.5 bonus for resolution below 2.5 Å; never a penalty.
        drug_score += min(0.5, max(0.0, (2.5 - res)) * 0.2)
    drug_score = round(drug_score, 3)

    return DrugPdbCandidate(
        pdb_id=pdb_id,
        title=pdb_metadata.get("title") or "",
        resolution=pdb_metadata.get("resolution"),
        method=pdb_metadata.get("method"),
        ligand_resname=ligand_resname,
        pocket_residue_count=residue_count,
        druggability_score=druggability,
        pocket_source=pocket_source,
        state_inferred=pdb_metadata.get("inferred_state", "unclear"),
        state_consistent=bool(pdb_metadata.get("state_consistent")),
        name_match=float(pdb_metadata.get("name_match") or 0.0),
        rcsb_score=float(pdb_metadata.get("rcsb_score") or 0.0),
        drug_pocket_score=drug_score,
        rejection_reason=rejection,
        notes=notes,
    )


async def pick_best_pdb(
    target_name: str,
    mode: Optional[Literal["agonist", "antagonist"]] = None,
    max_candidates_to_check: int = 8,
    max_results: int = 3,
) -> dict:
    """End-to-end: receptor name → ranked drug-ready PDBs.

    When ``mode`` is given, candidates are pre-filtered for the right
    conformational state. When ``mode`` is None, drug-pocket presence is the
    only filter — useful for kinases and other non-GPCR targets where the
    agonist/antagonist distinction is moot.
    """
    if not target_name or len(target_name.strip()) < 2:
        raise ValueError("target_name is required and must be ≥2 chars.")

    # 1) Get RCSB candidates. When no mode given, query both states and merge.
    # We oversample (×2.5) so that — when the search engine had to fall back
    # to UniProt and grabbed mostly apo / fragment structures from the head
    # of the cross-ref list — the drug-pocket validator downstream has a
    # realistic chance of finding at least one entry with a real co-crystal
    # ligand. Verification is still capped by max_candidates_to_check below.
    inner_cap = max(max_candidates_to_check * 3, 12)
    if mode:
        suggestions = await search_pdb_for_target(target_name, mode, max_candidates=inner_cap)
    else:
        agonist_task = search_pdb_for_target(target_name, "agonist", max_candidates=inner_cap)
        antagonist_task = search_pdb_for_target(target_name, "antagonist", max_candidates=inner_cap)
        agonist_hits, antagonist_hits = await asyncio.gather(agonist_task, antagonist_task)
        merged: dict[str, object] = {}
        for hit in (agonist_hits + antagonist_hits):
            if hit.pdb_id not in merged:
                merged[hit.pdb_id] = hit
        suggestions = list(merged.values())[:inner_cap]

    if not suggestions:
        return {
            "target_name": target_name,
            "mode": mode,
            "checked": 0,
            "candidates": [],
            "rejected": [],
            "warning": (
                f"No RCSB matches for target {target_name!r}. Try a more specific name "
                f"(e.g. 'dopamine D2 receptor' instead of 'D2'), or supply a PDB ID directly."
            ),
        }

    # 2) Verify each candidate with Pocket_App detect, in parallel
    verify_tasks = [
        _verify_candidate(
            hit.pdb_id,
            {
                "title": hit.title,
                "resolution": hit.resolution,
                "method": hit.method,
                "inferred_state": hit.inferred_state,
                "state_consistent": hit.state_consistent,
                "name_match": hit.name_match,
                "rcsb_score": hit.rcsb_score,
            },
            mode,
        )
        for hit in suggestions
    ]
    verified = await asyncio.gather(*verify_tasks)

    accepted = [v for v in verified if v.rejection_reason is None]
    rejected = [v for v in verified if v.rejection_reason is not None]

    accepted.sort(key=lambda c: c.drug_pocket_score, reverse=True)

    return {
        "target_name": target_name,
        "mode": mode,
        "checked": len(verified),
        "candidates": [_candidate_to_dict(c) for c in accepted[:max_results]],
        "rejected": [_candidate_to_dict(c) for c in rejected],
        "top_pick": _candidate_to_dict(accepted[0]) if accepted else None,
        "warning": (
            None if accepted
            else f"None of the {len(verified)} candidates checked had a real drug pocket — "
                 f"the target may not have any drug co-crystal in RCSB. "
                 f"See ``rejected`` for the reasons."
        ),
    }


def _candidate_to_dict(c: DrugPdbCandidate) -> dict:
    return {
        "pdb_id": c.pdb_id,
        "title": c.title,
        "resolution": c.resolution,
        "method": c.method,
        "ligand_resname": c.ligand_resname,
        "pocket_residue_count": c.pocket_residue_count,
        "druggability_score": c.druggability_score,
        "pocket_source": c.pocket_source,
        "state_inferred": c.state_inferred,
        "state_consistent": c.state_consistent,
        "name_match": c.name_match,
        "drug_pocket_score": c.drug_pocket_score,
        "rejection_reason": c.rejection_reason,
    }
