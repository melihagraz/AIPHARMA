from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class DiscoveryRequest(BaseModel):
    model_config = ConfigDict(extra="ignore")
    target_name: str = Field(..., min_length=2)
    pdb_id: str | None = None
    sequence: str | None = None
    disease_module: str = "general"
    desired_route: str = "oral"
    smiles_seed: str | None = None
    # Integer seed for the real REINVENT sampling service (if supported by the model).
    seed: int | None = Field(default=None, ge=0, le=2**31-1)
    generate_count: int = Field(default=400, ge=50, le=3000)
    filter_overrides: dict[str, dict[str, float | bool]] | None = None
    #: Farmakoför + cep kısıtları (ana platformdan); meta ve skor tuzaklaması için
    pharmacophore_context: dict[str, Any] | None = None
    #: Per-request engine override: "demo" (sentetik motor) | "production" (REINVENT 4 + ADMET-AI).
    #: None ise sunucudaki universal_discovery_mode kullanılır.
    engine_mode: str | None = Field(default=None, description="demo | production | None=use server default")
    #: Target-relevance: explicit list of known active SMILES. When provided,
    #: skips ChEMBL lookup. Used by production pipeline for Tanimoto-based scoring.
    known_actives_smiles: list[str] | None = Field(default=None)
    #: When true (default), production pipeline fetches known actives from ChEMBL
    #: by ``target_name`` if no explicit ``known_actives_smiles`` was given.
    use_chembl: bool = Field(default=True)
    #: ChEMBL filter: only activities with pChEMBL ≥ this value (7.0 ≈ IC50 100 nM).
    chembl_pchembl_min: float = Field(default=7.0, ge=0.0, le=15.0)
    #: ChEMBL: max known actives to pull as the reference set.
    chembl_limit: int = Field(default=50, ge=1, le=500)
    #: Run AutoDock Vina docking on top pre-ranked candidates. Default true in
    #: production mode. Disable for fast preview runs.
    enable_docking: bool = Field(default=True)
    #: How many top pre-ranked candidates to send through Vina (docking is slow).
    dock_top_k: int = Field(default=30, ge=1, le=200)
    #: Optional explicit Vina receptor name (e.g. "EGFR"). If unset, falls back
    #: to ``target_name``. Only used if a matching pre-bundled PDBQT exists.
    dock_target: str | None = None
    #: Bypass the target_name validation gate. When True, the production
    #: pipeline will run even if ``target_name`` has zero recognised signal
    #: (no ChEMBL actives, no fine-tuned prior, no bundled Vina receptor).
    #: The resulting batch is target-agnostic — only Lipinski + ADMET filters
    #: apply — and ``meta.validation.degraded`` is set so the UI can warn.
    force_generic: bool = Field(default=False)


class Candidate(BaseModel):
    rank: int
    smiles: str
    module: str
    score: float
    # None when Vina docking was unavailable for this run (target not bundled,
    # Docking_App down, individual ligand prep failed). UI must treat as "—".
    dock_score: float | None = None
    alerts: list[str]
    component_scores: dict[str, float]
    descriptors: dict[str, float | None]


class DiscoveryResponse(BaseModel):
    target: str
    module: str
    total_generated: int
    after_filter_count: int
    top_n: int
    candidates: list[Candidate]
    meta: dict[str, Any]
