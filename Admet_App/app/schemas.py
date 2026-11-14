"""Pydantic schemas for the ADMET prediction API."""

from __future__ import annotations

from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class PredictRequest(BaseModel):
    smiles: List[str] = Field(
        ..., min_length=1, max_length=1000, description="List of SMILES to predict (1..1000)"
    )


class MoleculeAdmet(BaseModel):
    smiles: str
    raw: Dict[str, Optional[float]] = Field(
        ..., description="Raw model outputs: 52 endpoints incl. BBB, hERG, CYPs, Caco2, half-life, tox classes"
    )
    percentile: Dict[str, Optional[float]] = Field(
        ...,
        description="DrugBank-approved percentile per endpoint (0–100). e.g. hERG=65 means safer than 65% of approved drugs.",
    )


class PredictResponse(BaseModel):
    count: int
    predictions: List[MoleculeAdmet]
    engine: str = "admet-ai (Stanford, MIT)"


class HealthResponse(BaseModel):
    status: str
    model_loaded: bool
