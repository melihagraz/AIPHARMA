import logging
import time
import httpx

from app.core.config import get_settings
from app.db.session import SessionLocal
from app.models.job import Job
from app.tasks.celery_app import celery_app

logger = logging.getLogger(__name__)
settings = get_settings()


@celery_app.task(name="jobs.run_pipeline")
def run_pipeline_async(job_id: str) -> None:
    db = SessionLocal()
    try:
        job = db.query(Job).filter(Job.id == job_id).first()
        if not job:
            logger.warning("Job %s not found for async execution", job_id)
            return

        job.status = "running"
        db.commit()

        payload = {"receptor_id": job.pdb_id or job.receptor_name or "unknown", "mode": job.mode, "count": job.candidate_count}
        with httpx.Client(timeout=10.0) as client:
            client.post(f"{settings.generator_service_url}/api/v1/generator/jobs", json=payload)

        time.sleep(1)
        job.status = "completed"
        db.commit()
    except Exception as exc:
        logger.exception("Pipeline failed for job %s: %s", job_id, exc)
        failed = db.query(Job).filter(Job.id == job_id).first()
        if failed:
            failed.status = "failed"
            db.commit()
    finally:
        db.close()
