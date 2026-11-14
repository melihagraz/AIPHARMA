"""Pydantic schemas for the REINVENT sampler API."""

from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class SampleRequest(BaseModel):
    batch_size: int = Field(32, ge=1, le=4096, description="Number of SMILES to draw")
    seed: Optional[int] = Field(None, description="Optional torch seed for reproducibility")
    validate_with_rdkit: bool = Field(
        True, description="If true, attach RDKit validity + basic descriptors per molecule"
    )
    reference_smiles: Optional[List[str]] = Field(
        None,
        description=(
            "Optional list of known active SMILES (target-relevance reference set). "
            "When provided, each generated molecule gets a max_tanimoto field — "
            "highest Tanimoto similarity (ECFP4) vs the reference set. "
            "Higher = closer to known binders for the target."
        ),
    )
    target_prior: Optional[str] = Field(
        None,
        description=(
            "Optional fine-tuned prior name (e.g. 'EGFR'). Maps to "
            "prior/finetuned/<name>_finetuned.prior. None / 'base' → base "
            "ChEMBL 25 prior. Unknown names silently fall back to base."
        ),
    )


class MoleculeOut(BaseModel):
    smiles: str
    nll: float = Field(..., description="Negative log-likelihood under the prior; lower = more likely")
    valid: Optional[bool] = None
    mw: Optional[float] = None
    logp: Optional[float] = None
    qed: Optional[float] = None
    tpsa: Optional[float] = None
    max_tanimoto: Optional[float] = Field(
        None,
        description="Max ECFP4 Tanimoto similarity to the reference active set (0..1); only set if reference_smiles was provided",
    )


class SampleResponse(BaseModel):
    requested: int
    returned: int
    valid_count: Optional[int] = Field(
        None, description="Number of molecules parseable by RDKit (only set when validate=True)"
    )
    candidates: List[MoleculeOut]
    engine: str = "reinvent4-prior (AI Pharma adapted)"
    prior_vocab_size: int
    prior_max_seq_len: int
    prior_used: str = Field(
        default="base",
        description="Which prior was actually used: 'base' or '<TARGET>' (e.g. 'EGFR')",
    )


class PriorsResponse(BaseModel):
    available: List[str] = Field(
        ..., description="Available prior names. 'base' is always present; fine-tuned targets follow."
    )


class HealthResponse(BaseModel):
    status: str
    prior_loaded: bool
    vocab_size: Optional[int] = None
    max_sequence_length: Optional[int] = None
