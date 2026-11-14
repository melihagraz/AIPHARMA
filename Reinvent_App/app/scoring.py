"""Target-relevance scoring via Tanimoto similarity to known actives.

Given a generated molecule and a reference set of known target actives,
compute the maximum Tanimoto similarity (Morgan/ECFP4 fingerprints) as a
proxy for "this molecule is in the chemical neighbourhood of known binders".

This is a real, defensible mechanism used widely in AI drug discovery
papers; it is not RL-fine-tuning, but it makes the system *target-aware*
when the prior alone is not.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

# RDKit imports are heavy; we keep them at module load time because RDKit
# is a required dependency for any production motor work.
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")  # silence parse warnings on noisy SMILES

# ECFP4 ≡ Morgan radius 2; 2048-bit folded fingerprint is the standard
# benchmark choice in cheminformatics literature.
_FP_RADIUS = 2
_FP_BITS = 2048


def _fingerprint(smiles: str):
    """Return Morgan/ECFP4 fingerprint or None if SMILES does not parse."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, _FP_RADIUS, nBits=_FP_BITS)


def build_reference_fingerprints(reference_smiles: Sequence[str]) -> list:
    """Pre-compute fingerprints for the reference active set.

    Invalid SMILES are dropped silently. Caller can compare returned length
    to input length to detect parse failures.
    """
    out = []
    for smi in reference_smiles:
        fp = _fingerprint(smi)
        if fp is not None:
            out.append(fp)
    return out


def max_tanimoto_to_references(query_smiles: str, reference_fps: List) -> Optional[float]:
    """Return max Tanimoto similarity between ``query_smiles`` and the reference set.

    Returns None when the query fails to parse or the reference set is empty.
    """
    if not reference_fps:
        return None
    qfp = _fingerprint(query_smiles)
    if qfp is None:
        return None
    sims = DataStructs.BulkTanimotoSimilarity(qfp, reference_fps)
    return float(max(sims)) if sims else None


def score_batch_against_references(
    smiles_list: Sequence[str], reference_fps: List
) -> list[Optional[float]]:
    """Vectorized version: return max-Tanimoto for each input SMILES."""
    if not reference_fps:
        return [None] * len(smiles_list)
    out: list[Optional[float]] = []
    for smi in smiles_list:
        qfp = _fingerprint(smi)
        if qfp is None:
            out.append(None)
            continue
        sims = DataStructs.BulkTanimotoSimilarity(qfp, reference_fps)
        out.append(float(max(sims)) if sims else None)
    return out
