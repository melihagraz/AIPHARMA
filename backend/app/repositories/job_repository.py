from sqlalchemy.orm import Session

from app.models.audit_log import AuditLog
from app.models.job import Job
from app.schemas.job import JobCreateRequest


class JobRepository:
    def __init__(self, db: Session):
        self.db = db

    def create(self, job_id: str, payload: JobCreateRequest) -> Job:
        entity = Job(
            id=job_id,
            tenant_id=payload.tenant_id,
            project_id=payload.project_id,
            status="queued",
            mode=payload.mode.value,
            receptor_name=payload.receptor_name,
            pdb_id=payload.pdb_id,
            uniprot_id=payload.uniprot_id,
            fasta_sequence=payload.fasta_sequence,
            candidate_count=payload.candidate_count,
        )
        self.db.add(entity)
        self.db.commit()
        self.db.refresh(entity)
        return entity

    def list(self, tenant_id: str | None = None, project_id: str | None = None) -> list[Job]:
        query = self.db.query(Job)
        if tenant_id:
            query = query.filter(Job.tenant_id == tenant_id)
        if project_id:
            query = query.filter(Job.project_id == project_id)
        return query.order_by(Job.created_at.desc()).all()

    def get(self, job_id: str) -> Job | None:
        return self.db.query(Job).filter(Job.id == job_id).first()

    def delete(self, job_id: str) -> bool:
        """Remove job and any audit rows keyed by this job id (outputs, etc.)."""
        job = self.get(job_id)
        if not job:
            return False
        self.db.query(AuditLog).filter(AuditLog.resource_id == job_id).delete(synchronize_session=False)
        self.db.delete(job)
        self.db.commit()
        return True
