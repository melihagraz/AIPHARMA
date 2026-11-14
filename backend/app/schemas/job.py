from datetime import datetime
from enum import Enum
from pydantic import BaseModel, Field


class FunctionalMode(str, Enum):
    agonist = "agonist"
    antagonist = "antagonist"
    partial_agonist = "partial_agonist"


class JobStatus(str, Enum):
    queued = "queued"
    running = "running"
    completed = "completed"
    failed = "failed"
    cancelled = "cancelled"


class JobCreateRequest(BaseModel):
    tenant_id: str = Field(default="default", min_length=1)
    project_id: str = Field(default="default", min_length=1)
    pdb_id: str | None = None
    uniprot_id: str | None = None
    fasta_sequence: str | None = None
    receptor_name: str | None = None
    mode: FunctionalMode = FunctionalMode.agonist
    candidate_count: int = Field(default=100, ge=1, le=1000)


class JobResponse(BaseModel):
    id: str
    tenant_id: str
    project_id: str
    status: JobStatus
    mode: FunctionalMode
    receptor_name: str | None = None
    created_at: datetime
    updated_at: datetime


class JobDetailResponse(JobResponse):
    pdb_id: str | None = None
    uniprot_id: str | None = None
    fasta_sequence: str | None = None
    candidate_count: int


class JobStatusUpdate(BaseModel):
    status: JobStatus
