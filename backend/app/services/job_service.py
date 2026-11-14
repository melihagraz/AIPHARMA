import logging
from datetime import datetime
from uuid import uuid4
from sqlalchemy.orm import Session

from app.repositories.job_repository import JobRepository
from app.core.config import get_settings
from app.schemas.job import JobCreateRequest


class JobService:
    def __init__(self, db: Session):
        self.repo = JobRepository(db)
        self.settings = get_settings()

    def create_job(self, payload: JobCreateRequest):
        job_id = str(uuid4())
        created = self.repo.create(job_id=job_id, payload=payload)
        if self.settings.enable_async_dispatch:
            try:
                # Optional worker dependency: do not break the API import in
                # an offline environment with dispatch disabled.
                from app.tasks.job_tasks import run_pipeline_async
                run_pipeline_async.delay(job_id)
            except Exception as exc:
                created.status = "failed"
                self.repo.db.commit()
                logging.getLogger(__name__).error("Async dispatch failed for %s: %s", job_id, exc)
        else:
            # No worker is attached. Leave job queued, not falsely completed.
            logging.getLogger(__name__).info("Job %s queued without async worker", job_id)
        return created

    def list_jobs(self, tenant_id: str | None = None, project_id: str | None = None):
        jobs = self.repo.list(tenant_id=tenant_id, project_id=project_id)
        changed = False
        for job in jobs:
            changed = self._advance_job_state(job) or changed
        if changed:
            self.repo.db.commit()
        return jobs

    def get_job(self, job_id: str):
        job = self.repo.get(job_id)
        if not job:
            return None
        if self._advance_job_state(job):
            self.repo.db.commit()
        return job

    def delete_job(self, job_id: str) -> bool:
        return self.repo.delete(job_id)

    def _advance_job_state(self, job) -> bool:
        """No synthetic completion: only a worker may change real job status."""
        return False
