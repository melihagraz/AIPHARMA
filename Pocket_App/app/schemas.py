"""Pydantic schemas for Pocket_App."""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import BaseModel, Field

# ---------------------------------------------------------------------------
# Pocket detection
# ---------------------------------------------------------------------------


class PocketDetectRequest(BaseModel):
    """Locate a druggable pocket.

    Resolution order:
      1. ``target`` — one of the bundled receptors (EGFR, JAK2, ...). Pocket
         center is read from the existing Vina conf file (validated, fast).
      2. ``pdb_id`` + ``ligand_resname`` — fetch the PDB from RCSB and use
         the bound ligand's center of geometry as the pocket center.
      3. ``pdb_id`` only — same as above, auto-detect a non-trivial HETATM
         (skips waters, ions, common buffers).
    """

    target: Optional[str] = Field(
        default=None,
        description="Bundled target name (e.g. 'EGFR'). Fastest path.",
    )
    pdb_id: Optional[str] = Field(
        default=None,
        description="4-char PDB ID (e.g. '1M17'). Fetched from RCSB.",
    )
    ligand_resname: Optional[str] = Field(
        default=None,
        description="3-char HETATM residue name (e.g. 'IRE' for gefitinib). "
        "If omitted, the first non-trivial HETATM is used.",
    )
    prepared_pdb_path: Optional[str] = Field(
        default=None,
        description="Absolute path to a pre-prepared PDB (PDBFixer-cleaned: "
        "missing residues repaired, protonated, waters cleaned). When given, "
        "the file is read directly and RCSB fetch is skipped.",
    )
    use_fpocket: Optional[bool] = Field(
        default=None,
        description="Cavity detection engine for the pdb_id path. True forces "
        "fpocket (real Voronoi/alpha-sphere detection + computed druggability); "
        "False forces the co-crystal-ligand shortcut. None (default) = use "
        "fpocket when available, else fall back to the ligand shortcut.",
    )
    min_druggability: float = Field(
        default=0.60,
        description="fpocket druggability threshold marking a pocket as a "
        "druggable candidate (manuscript convention). The top-scoring pocket "
        "is returned regardless; this only sets meets_druggability_threshold.",
    )


class PocketResult(BaseModel):
    center_x: float
    center_y: float
    center_z: float
    size_x: float = Field(default=24.0, description="Bounding box X (Å)")
    size_y: float = Field(default=24.0, description="Bounding box Y (Å)")
    size_z: float = Field(default=24.0, description="Bounding box Z (Å)")
    druggability_score: float = Field(
        default=0.0,
        description="0..1 druggability. For source='fpocket' this is fpocket's "
        "computed Druggability Score (Voronoi/alpha-sphere geometry). For "
        "'bundled_receptor'/'cocrystal_ligand' it is a fixed confidence "
        "(1.0/0.95); 'fallback' is a low-confidence guess.",
    )
    source: Literal["fpocket", "bundled_receptor", "cocrystal_ligand", "fallback"]
    ligand_resname: Optional[str] = None
    pocket_residues: list[str] = Field(
        default_factory=list,
        description="Up to ~25 residues lining the pocket (e.g. ['LYS745','MET793']).",
    )
    volume: Optional[float] = Field(
        default=None, description="Pocket volume in Å³ (fpocket only)."
    )
    n_alpha_spheres: Optional[int] = Field(
        default=None, description="Number of Voronoi alpha spheres (fpocket only)."
    )
    meets_druggability_threshold: Optional[bool] = Field(
        default=None,
        description="Whether druggability_score ≥ the requested threshold "
        "(fpocket only). False = highest-scoring pocket but below threshold.",
    )
    notes: Optional[str] = None


class PocketDetectResponse(BaseModel):
    target: Optional[str] = None
    pdb_id: Optional[str] = None
    pocket: PocketResult


# ---------------------------------------------------------------------------
# Pharmacophore extraction
# ---------------------------------------------------------------------------


class PharmacophoreFeature(BaseModel):
    """One 3D pharmacophore feature."""

    type: Literal[
        "Donor", "Acceptor", "Aromatic", "Hydrophobe", "PosIonizable", "NegIonizable"
    ] = Field(description="RDKit BaseFeatures category.")
    x: float
    y: float
    z: float
    tolerance: float = Field(
        default=1.5,
        description="Acceptable spatial deviation in Å when matching (default 1.5).",
    )
    label: str = Field(default="", description="Human-readable feature name.")


class PharmacophoreRequest(BaseModel):
    """Extract a 3D pharmacophore.

    Strategy:
      1. ``target`` (bundled) → parse bundled co-crystal PDB → extract bound ligand.
      2. ``pdb_id`` + optional ``ligand_resname`` → fetch + extract.
      3. ``prepared_pdb_path`` (when given) → read directly, skip RCSB fetch.
    """

    target: Optional[str] = None
    pdb_id: Optional[str] = None
    ligand_resname: Optional[str] = None
    prepared_pdb_path: Optional[str] = None


class PharmacophoreResponse(BaseModel):
    target: Optional[str] = None
    pdb_id: Optional[str] = None
    ligand_smiles: Optional[str] = None
    ligand_resname: Optional[str] = None
    feature_count: int
    features: list[PharmacophoreFeature]
    source: Literal["cocrystal", "fallback"]
    notes: Optional[str] = None


# ---------------------------------------------------------------------------
# Health
# ---------------------------------------------------------------------------


class HealthResponse(BaseModel):
    status: str = "ok"
    rdkit_pharm3d: bool
    bio_pdb: bool
    fpocket: bool = Field(
        default=False,
        description="Whether fpocket cavity detection is reachable "
        "(directly or via WSL on Windows).",
    )
    bundled_targets: list[str]
