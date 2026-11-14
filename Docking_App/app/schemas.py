"""Pydantic schemas for the Docking API."""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field, model_validator
import math


class DockRequest(BaseModel):
    target: str = Field(..., description="Pre-bundled target name, e.g. 'EGFR' (case-insensitive)")
    smiles: List[str] = Field(..., min_length=1, max_length=500, description="Ligand SMILES list")
    exhaustiveness: int = Field(8, ge=1, le=32, description="Vina exhaustiveness (search effort)")
    num_modes: int = Field(5, ge=1, le=20, description="How many top docking poses to keep")
    cpu_per_dock: int = Field(2, ge=1, le=16, description="CPU threads per individual docking")
    workers: int = Field(4, ge=1, le=16, description="Parallel docking jobs")
    # Optional pocket-box override — when supplied (typically by Pocket Analysis
    # upstream), Vina is told to search inside this box instead of the bundled
    # default grid. All six fields must be present together; any missing one
    # falls back to the receptor's pre-bundled conf.txt grid (legacy behaviour).
    center_x: Optional[float] = Field(None, description="Pocket centre X (Å); overrides bundled grid")
    center_y: Optional[float] = Field(None, description="Pocket centre Y (Å)")
    center_z: Optional[float] = Field(None, description="Pocket centre Z (Å)")
    size_x: Optional[float] = Field(None, description="Box X size (Å); overrides bundled grid")
    size_y: Optional[float] = Field(None, description="Box Y size (Å)")
    size_z: Optional[float] = Field(None, description="Box Z size (Å)")

    @model_validator(mode="after")
    def validate_box(self):
        keys = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
        values = [getattr(self, k) for k in keys]
        if any(v is not None for v in values):
            if not all(v is not None and math.isfinite(v) for v in values):
                raise ValueError("all six docking box fields must be supplied and finite")
            if any(getattr(self, k) <= 0 or getattr(self, k) > 60 for k in keys[3:]):
                raise ValueError("docking grid dimensions must be in (0, 60] Å")
        return self
    # Custom receptor override — when the caller has its own PDB (typically
    # user-uploaded via /api/pocket/upload-pdb, e.g. for an allosteric or
    # mutant pocket that's not in the bundled 8) it supplies the raw PDB
    # text here. The server runs meeko to prep a PDBQT on the fly, caches
    # by content hash, and docks against THAT receptor instead of the
    # bundled file resolved by ``target``. ``target`` is still required —
    # it's used as the cache label and as a fallback if prep fails.
    custom_receptor_pdb: Optional[str] = Field(
        None,
        max_length=20_000_000,
        description=(
            "Raw PDB text. When supplied, server preps a custom PDBQT via "
            "meeko and uses it instead of the bundled receptor."
        ),
    )


class DockResultOut(BaseModel):
    smiles: str
    target: str
    affinity: Optional[float] = Field(
        None, description="Best (most negative) binding affinity in kcal/mol; None on failure"
    )
    error: Optional[str] = None
    num_modes: int = 0
    # Top-pose PDB block (AutoDock charge/atom-type columns stripped).
    # Lets the calling backend persist the docked pose alongside the
    # affinity, so the Results page can render the real binding geometry
    # in NGL.js instead of an RDKit free-state conformer.
    pose_pdb: Optional[str] = Field(
        None,
        description=(
            "Top docked pose as a vanilla PDB block (truncated at col 66 "
            "to strip PDBQT charge + atom-type extensions)."
        ),
    )
    # All Vina poses concatenated as a multi-MODEL PDB block. Lets the
    # frontend cycle through alternate binding modes without an extra
    # round-trip. None when capture failed or only one pose was produced.
    pose_pdb_all: Optional[str] = Field(
        None,
        description="All Vina poses as a multi-MODEL PDB block (top mode first).",
    )


class DockResponse(BaseModel):
    target: str
    count: int
    success_count: int
    results: List[DockResultOut]
    engine: str = "AutoDock Vina (runtime version not reported by API)"


class HealthResponse(BaseModel):
    status: str
    vina_binary: bool
    meeko_available: bool = False
    ready_for_docking: bool = False
    available_targets: List[str]


class TargetsResponse(BaseModel):
    targets: List[str]
