from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from app.db.session import get_db
from app.core.security import Principal, require_session
from app.schemas.job import JobCreateRequest, JobDetailResponse, JobResponse
from app.services.job_service import JobService

router = APIRouter(prefix="/jobs", tags=["jobs"])


@router.post("", response_model=JobResponse)
async def create_job(payload: JobCreateRequest, db: Session = Depends(get_db), principal: Principal = Depends(require_session)) -> JobResponse:
    if not any([payload.pdb_id, payload.uniprot_id, payload.fasta_sequence, payload.receptor_name]):
        raise HTTPException(status_code=400, detail="Provide at least one receptor input field.")

    # Scope new jobs to authenticated user; ignore user-supplied tenant_id.
    payload = payload.model_copy(update={"tenant_id": principal.user_id})
    service = JobService(db)
    job = service.create_job(payload)
    return JobResponse(
        id=job.id,
        tenant_id=job.tenant_id,
        project_id=job.project_id,
        status=job.status,
        mode=job.mode,
        receptor_name=job.receptor_name,
        created_at=job.created_at,
        updated_at=job.updated_at,
    )


@router.get("", response_model=list[JobResponse])
async def list_jobs(
    tenant_id: str | None = Query(default=None),
    project_id: str | None = Query(default=None),
    db: Session = Depends(get_db),
    principal: Principal = Depends(require_session),
) -> list[JobResponse]:
    service = JobService(db)
    jobs = service.list_jobs(tenant_id=principal.user_id, project_id=project_id)
    return [
        JobResponse(
            id=j.id,
            tenant_id=j.tenant_id,
            project_id=j.project_id,
            status=j.status,
            mode=j.mode,
            receptor_name=j.receptor_name,
            created_at=j.created_at,
            updated_at=j.updated_at,
        )
        for j in jobs
    ]


@router.get("/{job_id}", response_model=JobDetailResponse)
async def get_job(job_id: UUID, db: Session = Depends(get_db), principal: Principal = Depends(require_session)) -> JobDetailResponse:
    """UUID-only so ``/jobs/start-denovo`` is never captured as a job id (POST would otherwise see 405)."""
    service = JobService(db)
    job = service.get_job(str(job_id))
    if not job or job.tenant_id != principal.user_id:
        raise HTTPException(status_code=404, detail="Job not found.")

    return JobDetailResponse(
        id=job.id,
        tenant_id=job.tenant_id,
        project_id=job.project_id,
        status=job.status,
        mode=job.mode,
        receptor_name=job.receptor_name,
        created_at=job.created_at,
        updated_at=job.updated_at,
        pdb_id=job.pdb_id,
        uniprot_id=job.uniprot_id,
        fasta_sequence=job.fasta_sequence,
        candidate_count=job.candidate_count,
    )


@router.delete("/{job_id}")
async def delete_job(job_id: UUID, db: Session = Depends(get_db), principal: Principal = Depends(require_session)) -> dict:
    """UUID-only so paths like ``/jobs/start-denovo`` are never captured here (405 on POST)."""
    sid = str(job_id)
    service = JobService(db)
    existing = service.get_job(sid)
    if not existing or existing.tenant_id != principal.user_id:
        raise HTTPException(status_code=404, detail="Job not found.")
    if not service.delete_job(sid):
        raise HTTPException(status_code=404, detail="Job not found.")
    return {"ok": True, "job_id": sid}
