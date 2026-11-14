import json
from uuid import uuid4

from fastapi import APIRouter, Depends, Header, HTTPException, Query, Request
from sqlalchemy.orm import Session

from app.core.config import get_settings
from app.core.security import issue_session_token, require_admin
from app.core.login_throttle import check_login_throttle, record_login_failure
from app.db.session import get_db
from app.models.audit_log import AuditLog
from app.models.job import Job
from app.models.user_account import UserAccount
from app.schemas.admin import (
    AdminAuditLogResponse,
    AdminLoginRequest,
    AdminLoginResponse,
    AdminUserCreateRequest,
    AdminUserResponse,
    AdminUserUpdateRequest,
)

router = APIRouter(prefix="/admin", tags=["admin"])
settings = get_settings()


def write_audit(db: Session, actor_id: str, action: str, resource_type: str, resource_id: str, payload: dict):
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id=actor_id,
            action=action,
            resource_type=resource_type,
            resource_id=resource_id,
            payload=json.dumps(payload, ensure_ascii=True),
        )
    )
    db.commit()


def to_user_response(user: UserAccount) -> AdminUserResponse:
    return AdminUserResponse(
        id=user.id,
        email=user.email,
        full_name=user.full_name,
        is_active=user.is_active,
        is_admin=user.is_admin,
        pocket_quota=user.pocket_quota,
        denovo_quota=user.denovo_quota,
        admet_quota=user.admet_quota,
        docking_quota=user.docking_quota,
        created_at=user.created_at,
        updated_at=user.updated_at,
    )


@router.post("/login", response_model=AdminLoginResponse)
async def admin_login(payload: AdminLoginRequest, request: Request) -> AdminLoginResponse:
    import hmac
    check_login_throttle(request, payload.username)
    if not settings.admin_password or not settings.session_signing_key or (
        payload.username != settings.admin_username or
        not hmac.compare_digest(payload.password, settings.admin_password)
    ):
        record_login_failure(request, payload.username)
        raise HTTPException(status_code=401, detail="Invalid admin credentials.")
    return AdminLoginResponse(access_token=issue_session_token("admin-local", "admin"))


@router.get("/logs", response_model=list[AdminAuditLogResponse], dependencies=[Depends(require_admin)])
async def list_audit_logs(limit: int = Query(default=100, ge=1, le=500), db: Session = Depends(get_db)):
    logs = db.query(AuditLog).order_by(AuditLog.created_at.desc()).limit(limit).all()
    return [
        AdminAuditLogResponse(
            id=log.id,
            actor_id=log.actor_id,
            action=log.action,
            resource_type=log.resource_type,
            resource_id=log.resource_id,
            payload=log.payload,
            created_at=log.created_at,
        )
        for log in logs
    ]


@router.get("/users", response_model=list[AdminUserResponse], dependencies=[Depends(require_admin)])
async def list_users(db: Session = Depends(get_db)):
    users = db.query(UserAccount).order_by(UserAccount.created_at.desc()).all()
    return [to_user_response(user) for user in users]


@router.post("/users", response_model=AdminUserResponse, dependencies=[Depends(require_admin)])
async def create_user(payload: AdminUserCreateRequest, db: Session = Depends(get_db)):
    existing = db.query(UserAccount).filter(UserAccount.email == payload.email).first()
    if existing:
        raise HTTPException(status_code=409, detail="User already exists.")
    user = UserAccount(
        id=str(uuid4()),
        email=payload.email,
        full_name=payload.full_name,
        is_admin=payload.is_admin,
        pocket_quota=payload.pocket_quota,
        denovo_quota=payload.denovo_quota,
        admet_quota=payload.admet_quota,
        docking_quota=payload.docking_quota,
    )
    db.add(user)
    db.commit()
    db.refresh(user)
    write_audit(db, "admin", "user_created", "user", user.id, {"email": user.email})
    return to_user_response(user)


@router.patch("/users/{user_id}", response_model=AdminUserResponse, dependencies=[Depends(require_admin)])
async def update_user(user_id: str, payload: AdminUserUpdateRequest, db: Session = Depends(get_db)):
    user = db.query(UserAccount).filter(UserAccount.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found.")
    updates = payload.model_dump(exclude_none=True)
    for field, value in updates.items():
        setattr(user, field, value)
    db.commit()
    db.refresh(user)
    write_audit(db, "admin", "user_updated", "user", user.id, updates)
    return to_user_response(user)


@router.post("/jobs/{job_id}/cancel", dependencies=[Depends(require_admin)])
async def cancel_job(job_id: str, db: Session = Depends(get_db)):
    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")
    if job.status in {"completed", "failed", "cancelled"}:
        raise HTTPException(status_code=409, detail=f"Job already {job.status}.")
    job.status = "cancelled"
    db.commit()
    write_audit(db, "admin", "job_cancelled", "job", job.id, {"previous_status": "running_or_queued"})
    return {"ok": True, "job_id": job.id, "status": job.status}
