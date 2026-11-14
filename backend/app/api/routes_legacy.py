import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional
from uuid import UUID, uuid4

import httpx
from fastapi import APIRouter, Depends, File, Form, HTTPException, Query, Request, UploadFile
from pydantic import ValidationError
from sqlalchemy import desc
from sqlalchemy.orm import Session

from app.core.config import get_settings
from app.core.security import hash_password, verify_password, issue_session_token, require_session, Principal
from app.core.login_throttle import check_login_throttle, record_login_failure
from app.db.session import get_db
from app.models.audit_log import AuditLog
from app.models.auth_user import AuthUser
from app.models.job import Job
from app.models.user_account import UserAccount
from app.repositories.job_repository import JobRepository
from app.schemas.job import JobCreateRequest
from app.services.job_service import JobService
from app.universal_discovery.models import DiscoveryRequest, DiscoveryResponse

router = APIRouter(prefix="/api", tags=["legacy"])
settings = get_settings()
logger = logging.getLogger(__name__)


def _normalize_pocket_token(raw) -> str:
    """Align client pocket ids with DB (e.g. JSON float 1.0 -> '1')."""
    if raw is None:
        return ""
    s = str(raw).strip()
    if not s:
        return ""
    try:
        f = float(s)
        if f == int(f):
            return str(int(f))
    except (TypeError, ValueError):
        pass
    return s


async def _discover_upstream(denovo_request: dict) -> tuple[dict, str]:
    """Run discovery; returns (payload_json, engine_label).

    Per-request override: if ``denovo_request["engine_mode"]`` is set to
    "demo" or "production", that wins over the server default.
    """
    cfg = get_settings()
    # Per-request override wins over server default
    requested_mode = (denovo_request.get("engine_mode") or "").strip().lower()
    if requested_mode == "production":
        mode = "production"
    elif requested_mode == "demo":
        mode = "embedded"  # the synthetic legacy motor lives in embedded/subprocess paths
    else:
        mode = (cfg.universal_discovery_mode or "production").strip().lower()
    # Legacy embedded, subprocess and HTTP integrations are demo-only. An
    # unannotated request must never silently switch to fabricated chemistry.
    if mode != "production" and requested_mode != "demo":
        raise HTTPException(status_code=503, detail=(
            "Non-production discovery is disabled by default. Explicitly request "
            "engine_mode='demo' for synthetic examples, or configure production services."
        ))
    base = (cfg.universal_ai_service_url or "").strip().rstrip("/")

    # PRODUCTION: real REINVENT 4 generator + ADMET-AI predictor + Vina pipeline.
    # The engine label here must match what ``production_pipeline.py`` writes
    # into ``meta.discovery_engine`` so the UI gets a single consistent string.
    # The label is "_vina"-suffixed because the production path runs Vina by
    # default whenever the target maps to a bundled receptor (see
    # production_pipeline.py:553). For unbundled targets the inner pipeline
    # still returns "reinvent4_admet_vina" because the Vina *attempt* is part
    # of the engine — whether it succeeded for that run is recorded separately
    # in ``meta.docking_status``.
    if mode == "production":
        from app.services.production_pipeline import (
            TargetValidationError,
            UnsupportedFilterError,
            run_production_discovery,
        )

        try:
            return await run_production_discovery(denovo_request), "reinvent4_admet_vina"
        except UnsupportedFilterError as exc:
            raise HTTPException(status_code=422, detail=str(exc)) from exc
        except TargetValidationError as exc:
            # Translate the pipeline's silent-failure gate into a 422 the
            # frontend can render with a friendly error. The detail object
            # carries the per-axis findings (ChEMBL hits, fine-tuned prior
            # availability, bundled receptor availability) so the UI can
            # tell the user exactly what's wrong with the target_name —
            # rather than "the pipeline failed" or worse, silently
            # producing 15 generic drug-like molecules pretending to be
            # target-specific.
            raise HTTPException(
                status_code=422,
                detail={
                    "error": "target_validation_failed",
                    "message": (
                        f"Target {exc.target!r} has no recognised signal: "
                        f"ChEMBL returned {exc.findings.get('chembl_hits', 0)} potent actives, "
                        f"no fine-tuned REINVENT prior is available, and no Vina docking receptor "
                        f"is bundled for this name. Either correct the target name, or rerun with "
                        f"force_generic=true to acknowledge that the pipeline will produce "
                        f"target-agnostic drug-like molecules only."
                    ),
                    "target": exc.target,
                    **exc.findings,
                },
            ) from exc

    async def _http_discover() -> dict:
        if not base:
            raise HTTPException(
                status_code=500,
                detail="UNIVERSAL_AI_SERVICE_URL / universal_ai_service_url is not set for HTTP discovery.",
            )
        url = f"{base}/api/discover"
        try:
            async with httpx.AsyncClient(timeout=httpx.Timeout(300.0)) as client:
                r = await client.post(url, json=denovo_request)
        except httpx.ConnectError as exc:
            raise HTTPException(
                status_code=503,
                detail=f"Discovery motor unreachable at {url}. The upstream service must be started. Synthetic demo is available only via engine_mode=demo.",
            ) from exc
        if r.status_code >= 400:
            detail_txt = r.text[:2000]
            try:
                body = r.json()
                if isinstance(body, dict) and body.get("detail") is not None:
                    detail_txt = str(body.get("detail"))
            except Exception:
                pass
            raise HTTPException(status_code=r.status_code, detail=f"Discovery HTTP error: {detail_txt}") from None
        try:
            return r.json()
        except Exception as exc:
            raise HTTPException(status_code=502, detail=f"Discovery response is not JSON: {exc}") from exc

    if mode == "http":
        return await _http_discover(), "http"

    if mode in ("yeni_motor_subprocess", "yeni_motor_script", "subprocess_yeni_motor"):
        from app.services.yeni_motor_discover_subprocess import run_discovery_via_yeni_motor_script

        return await run_discovery_via_yeni_motor_script(denovo_request), "universal_discovery_engine"

    from app.universal_discovery import run_universal_discovery

    try:
        return run_universal_discovery(denovo_request), "embedded"
    except ValidationError:
        raise
    except FileNotFoundError as exc:
        logger.warning("Embedded discovery config missing: %s", exc)
        if mode == "http_fallback" and base:
            return await _http_discover(), "http"
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    except Exception as exc:
        logger.exception("Embedded universal discovery failed")
        if mode == "http_fallback" and base:
            logger.warning("Retrying discovery via HTTP at %s", base)
            return await _http_discover(), "http"
        raise HTTPException(status_code=500, detail=f"Discovery pipeline failed: {exc}") from exc


async def motor_discover_impl(body: DiscoveryRequest) -> DiscoveryResponse:
    """POST /api/motor/discover — API_KULLANIM.txt bölüm 2; yeni_motor /api/discover ile aynı gövde (DiscoveryRequest)."""
    try:
        data, engine_label = await _discover_upstream(body.model_dump(exclude_none=False))
    except HTTPException:
        raise
    meta = data.get("meta")
    if not isinstance(meta, dict):
        meta = {}
    # Prefer the pipeline's own engine name when it set one (production_pipeline
    # writes "reinvent4_admet_vina" on success and "reinvent4_admet" on the
    # empty-result fallback path). The outer ``engine_label`` is the coarse
    # transport-level tag (production / http / embedded / universal_discovery_engine)
    # and is only used as a fallback when the inner pipeline didn't self-identify.
    data["meta"] = {**meta}
    if not data["meta"].get("discovery_engine"):
        data["meta"]["discovery_engine"] = engine_label
    is_synthetic = engine_label != "reinvent4_admet_vina"
    data["meta"].update({
        "engine_type": "synthetic_demo" if is_synthetic else "computational_prediction",
        "is_synthetic": is_synthetic,
        "research_validated": False,
    })
    data["meta"].setdefault("run_id", str(uuid4()))
    data["meta"].setdefault("seed", body.smiles_seed if is_synthetic else body.seed)
    data["meta"].setdefault("model_versions", {})
    data["meta"].setdefault("reference_hash", None)
    data["meta"].setdefault("source_data", "synthetic_demo" if is_synthetic else "unspecified")
    data["meta"].setdefault("degraded", bool(is_synthetic or
                                    (data["meta"].get("validation") or {}).get("degraded", False)))
    try:
        return DiscoveryResponse.model_validate(data)
    except ValidationError as exc:
        raise HTTPException(status_code=502, detail=exc.errors()) from exc


def _progress_for_job(job) -> int:
    """Best-available 0-100 progress percentage for a Job row.

    Pipeline writes ``Job.progress_pct`` (0..100) as it advances through
    the REINVENT → ADMET → Vina chain. When that value is present we
    return it verbatim; otherwise we fall back to coarse status-only
    estimates so older jobs (no granular tracking) still render
    something in the progress bar.
    """
    if job.status == "completed":
        return 100
    if job.status in {"failed", "cancelled"}:
        return 0
    pct = getattr(job, "progress_pct", None)
    if isinstance(pct, (int, float)) and 0 < pct <= 100:
        return int(pct)
    if job.status == "queued":
        return 5
    if job.status == "running":
        return 35
    return 0


def _progress_detail_for_job(job) -> dict:
    """Granular progress dict: {pct, step, message}.

    Used by /api/jobs/queue so the UI progress tracker can render a
    "REINVENT sampling 42%" style line instead of just a generic bar.
    All three fields default to None when the pipeline never updated
    them (e.g. legacy jobs from before progress tracking landed).
    """
    return {
        "pct": _progress_for_job(job),
        "step": getattr(job, "current_step", None) or None,
        "message": getattr(job, "progress_message", None) or None,
    }


def _message_for_job(job, section: str = "analysis") -> str:
    label = (
        "Pocket analysis"
        if section == "pocket_analysis"
        else (
            "Pharmacophore analysis"
            if section in {"farmakofor", "pharmacophore"}
            else ("De Novo analysis" if section in {"denovo", "de_novo"} else "Analysis")
        )
    )
    if job.status == "queued":
        return f"{label} queued..."
    if job.status == "running":
        return f"{label} in progress..."
    if job.status == "completed":
        return f"{label} completed."
    return job.status


def _job_outputs(db: Session, job_ids: list[str]) -> dict[str, dict]:
    if not job_ids:
        return {}
    rows = (
        db.query(AuditLog)
        .filter(
            AuditLog.resource_type == "job_output",
            AuditLog.resource_id.in_(job_ids),
        )
        .all()
    )
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        try:
            payload = json.loads(row.payload or "{}")
        except Exception:
            payload = {}
        payload["__action"] = row.action
        # Preserve the actor_id so the queue/hierarchy responses can surface
        # provenance ("did the agent produce this row, or did the user click
        # the button?"). The legacy UI writer uses actor_id="system"; the
        # Phase 2 agent uses actor_id="agent". The frontend renders a small
        # badge on each row based on this field — no DB schema change,
        # just propagating existing audit_log metadata up to the API surface.
        payload["__actor_id"] = row.actor_id or "system"
        grouped.setdefault(row.resource_id, []).append(payload)
    out: dict[str, dict] = {}
    for job_id, payloads in grouped.items():
        preferred = next((p for p in payloads if p.get("__action") == "pharmacophore_output"), None)
        if not preferred:
            preferred = next((p for p in payloads if p.get("__action") == "denovo_output"), None)
        if not preferred:
            preferred = next((p for p in payloads if p.get("__action") == "pocket_output"), None)
        if not preferred:
            preferred = payloads[-1]
        out[job_id] = preferred
    return out


def _source_label(actor_id: str | None) -> str:
    """Map audit_log.actor_id → a stable "source" string for the frontend.

    We normalise to two values the UI can switch on:
      - "agent" → chat agent (Phase 2 LLM) generated this run
      - "ui"    → button-driven UI workflow generated this run
    """
    if (actor_id or "").lower() == "agent":
        return "agent"
    return "ui"


def _pocket_payload_for_job(db: Session, job_id: str) -> dict:
    if not job_id or str(job_id).strip() in {"", "0"}:
        return {}
    row = (
        db.query(AuditLog)
        .filter(
            AuditLog.resource_type == "job_output",
            AuditLog.resource_id == str(job_id),
            AuditLog.action == "pocket_output",
        )
        .order_by(desc(AuditLog.created_at))
        .first()
    )
    if not row:
        return {}
    try:
        return json.loads(row.payload or "{}")
    except Exception:
        return {}


def _pharmacophore_payload_for_job(db: Session, job_id: str) -> dict | None:
    row = (
        db.query(AuditLog)
        .filter(
            AuditLog.resource_type == "job_output",
            AuditLog.resource_id == str(job_id),
            AuditLog.action == "pharmacophore_output",
        )
        .order_by(desc(AuditLog.created_at))
        .first()
    )
    if not row:
        return None
    try:
        return json.loads(row.payload or "{}")
    except Exception:
        return None


@router.get("/jobs/queue")
async def legacy_jobs_queue(db: Session = Depends(get_db)):
    service = JobService(db)
    jobs = service.list_jobs()
    outputs = _job_outputs(db, [job.id for job in jobs])
    return [
        (lambda section, actor: {
            "id": job.id,
            "job_id": job.id,
            "status": job.status,
            "message": _message_for_job(job, section),
            "created_at": job.created_at.strftime("%d.%m.%Y %H:%M") if isinstance(job.created_at, datetime) else str(job.created_at),
            "updated_at": job.updated_at.strftime("%d.%m.%Y %H:%M")
            if isinstance(job.updated_at, datetime)
            else (str(job.updated_at) if job.updated_at else None),
            "pdb_id": job.pdb_id or job.receptor_name or "---",
            "analysis_section": section,
            "progress": _progress_for_job(job),
            "progress_detail": _progress_detail_for_job(job),
            "results_json": outputs.get(job.id, {}).get("results_json") or outputs.get(job.id, {}).get("pockets", []),
            # "agent" if the chat-driven agent created this row, else "ui".
            # The frontend ResultsModal renders a 🤖 / 👤 badge per row and
            # exposes a 3-way filter toggle (All / Agent / UI) above the
            # tables. No DB schema change — we just propagate the existing
            # audit_log.actor_id up to the JSON surface.
            "source": _source_label(actor),
        })(
            outputs.get(job.id, {}).get("analysis_section", "pocket_analysis"),
            outputs.get(job.id, {}).get("__actor_id"),
        )
        for job in jobs
    ]


@router.get("/jobs/hierarchy")
async def legacy_jobs_hierarchy(db: Session = Depends(get_db)):
    service = JobService(db)
    # Sort newest-first so the UI dropdown shows the most recent analysis on top.
    # SQLite stores DATETIME columns as TEXT and SQLAlchemy hands them back as
    # ``datetime`` objects whose tzinfo depends on the writer: some pipeline
    # paths use timezone-aware ``datetime.now(timezone.utc)`` (e.g. RL +
    # selectivity persistence) while older code writes naive
    # ``datetime.utcnow()``. Mixing the two in one sorted() call raises
    # ``can't compare offset-naive and offset-aware datetimes`` and the
    # entire /jobs/hierarchy endpoint 500's — which then breaks the
    # ``fetchData`` chain in NewJobForm so the ADMET dropdown stays empty
    # ("ADMET sonuçları bulunamadı" even though they exist). Coerce every
    # value to naive-UTC for the sort key so the two writers can coexist
    # in the same table without breaking the read path.
    def _sort_key(j):
        d = getattr(j, "created_at", None)
        if d is None:
            return datetime.min
        # tz-aware → convert to UTC, then drop tzinfo to match the naive rows.
        if getattr(d, "tzinfo", None) is not None:
            try:
                from datetime import timezone as _tz
                return d.astimezone(_tz.utc).replace(tzinfo=None)
            except Exception:  # noqa: BLE001
                return datetime.min
        return d
    jobs = sorted(service.list_jobs(), key=_sort_key, reverse=True)
    outputs = _job_outputs(db, [job.id for job in jobs])

    def _build_entry(job, section: str, actor: str | None) -> dict:
        """Project a job into the hierarchy-dropdown shape.

        For De Novo / Docking rows we surface the same total_generated /
        admet_passed / candidates_saved trio that /jobs/admet-source-list
        exposes — Steps 4/5/6 dropdowns then share a single visual idiom
        ("EGFR — 30 ADMET ✓ / 50") instead of each tab having a different
        view of the same job. Pocket / pharmacophore rows leave the pass
        fields as None because the concept doesn't apply.
        """
        out = outputs.get(job.id, {}) or {}
        total_generated = None
        admet_passed = None
        candidates_saved = None
        if section in ("denovo", "denovo_design", "docking"):
            results_json = out.get("results_json") or {}
            saved_rows = results_json.get("data") or []
            if isinstance(saved_rows, list):
                candidates_saved = len(saved_rows) or None
            tg = out.get("total_generated")
            ap = out.get("after_filter_count")
            if isinstance(tg, (int, float)):
                total_generated = int(tg)
            if isinstance(ap, (int, float)):
                admet_passed = int(ap)
        return {
            "id": job.id,
            "job_id": job.id,
            "pdb_id": job.pdb_id or job.receptor_name or "---",
            "status": job.status,
            "message": _message_for_job(job, section),
            "date": job.created_at.strftime("%d.%m.%Y %H:%M"),
            "pockets": out.get("pockets", []),
            # Step 6 candidate-analysis dropdown filters on this field
            # (h.analysis_section === "denovo_design" || ... === "docking").
            # Previously omitted, so the dropdown was always empty.
            "analysis_section": section,
            # Same provenance hint as /jobs/queue — lets the Step 6 dropdown
            # show a small 🤖/👤 marker so the user can tell at a glance
            # whether a candidate row came from chat or from a UI run.
            "source": _source_label(actor),
            "short_id": job.id[:8],
            # Pass-count trio for Step 6 dropdown parity with Steps 4/5.
            "total_generated": total_generated,
            "admet_passed": admet_passed,
            "candidates_saved": candidates_saved,
        }

    return [
        _build_entry(
            job,
            outputs.get(job.id, {}).get("analysis_section", "pocket_analysis"),
            outputs.get(job.id, {}).get("__actor_id"),
        )
        for job in jobs
    ]


@router.post("/jobs/log-pharmacophore-analysis")
async def legacy_log_pharmacophore_analysis(body: dict, db: Session = Depends(get_db)):
    """Persist a real Pocket_App pharmacophore result as a completed Job + AuditLog.

    Same pattern as ``log-pocket-analysis``: the frontend's "Extract from
    co-crystal" button hits ``/api/pocket/pharmacophore`` synchronously,
    shows the features inline, then calls this endpoint so the result lands
    in the Job history. Without this hook, only the legacy
    ``/jobs/start-pharmacophore`` placeholder (residue_count // 3 = acceptor,
    etc.) would show up in the Analysis Results table — which is exactly the
    "all entries have 1/1/1/1/1 counts" symptom the user spotted.

    Replaces the synthetic Acceptor/Donor/Hydrophobic counts with the actual
    per-family tallies from RDKit's BaseFeatures.fdef (after the within-
    family dedup we added to Pocket_App/pharmacophore.py).
    """
    pdb_id = str(body.get("pdb_id") or "").strip().upper()
    target_name = str(body.get("target_name") or "").strip()
    pharma_result = body.get("pharmacophore_result") or {}
    if not pdb_id and not target_name:
        raise HTTPException(status_code=400, detail="pdb_id or target_name required")

    features = pharma_result.get("features") if isinstance(pharma_result, dict) else None
    if not isinstance(features, list):
        features = []

    # Tally features by family so the existing Analysis Results table renders
    # the same shape it already expects (acceptor / donor / hydrophobic /
    # positive / negative counts).
    by_family: dict[str, int] = {}
    for f in features:
        if not isinstance(f, dict):
            continue
        ftype = str(f.get("type") or "").strip().lower()
        if ftype:
            by_family[ftype] = by_family.get(ftype, 0) + 1

    receptor_value = pdb_id or target_name or "UNKNOWN"
    repo = JobRepository(db)
    job = repo.create(
        job_id=str(uuid4()),
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=100,
        ),
    )
    # Pharmacophore extraction already happened (sync) — mark completed.
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned is not None:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()

    # Compute the per-family tally once and reuse in BOTH top-level fields
    # (for direct backend consumers) and the nested ``results_json`` (which
    # is what /api/jobs/queue exposes to the UI's Pharmacophore table — the
    # legacy ``start-pharmacophore`` endpoint also uses this nested shape,
    # so we mirror it for table-rendering parity).
    counts = {
        "acceptor": int(by_family.get("acceptor", 0)),
        "donor": int(by_family.get("donor", 0)),
        "hydrophobic": int(
            by_family.get("hydrophobe", 0)
            + by_family.get("hydrophobic", 0)
            + by_family.get("lumpedhydrophobe", 0)
        ),
        "aromatic": int(by_family.get("aromatic", 0)),
        "positive": int(by_family.get("posionizable", 0) + by_family.get("positive", 0)),
        "negative": int(by_family.get("negionizable", 0) + by_family.get("negative", 0)),
    }
    audit_payload = {
        "analysis_section": "farmakofor",
        "job_id": job.id,
        "pdb_id": pdb_id or "---",
        "target_name": target_name or None,
        "ligand_resname": pharma_result.get("ligand_resname"),
        "ligand_smiles": pharma_result.get("ligand_smiles"),
        "source": pharma_result.get("source") or "cocrystal",
        # Top-level fields for any direct backend consumer
        **counts,
        "feature_count": len(features),
        "features": features,  # full coordinate list for downstream consumers
        # Nested results_json — this is what /api/jobs/queue serves to the
        # frontend Pharmacophore table. Without this, the row renders as
        # 0/0/0/0/0 even though the data is in the payload.
        "results_json": {
            **counts,
            "feature_count": len(features),
            "source": pharma_result.get("source") or "cocrystal",
            "ligand_resname": pharma_result.get("ligand_resname"),
        },
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="pharmacophore_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(audit_payload, ensure_ascii=True),
        )
    )
    db.commit()
    return {"job_id": job.id, "status": "completed", "pdb_id": pdb_id, "feature_count": len(features)}


@router.post("/jobs/log-pocket-analysis")
async def legacy_log_pocket_analysis(body: dict, db: Session = Depends(get_db)):
    """Persist a real Pocket_App result as a completed Job + AuditLog entry.

    The frontend's sync path (handleAnalyze → /api/pocket/detect) shows the
    pocket result immediately for UX, then calls this endpoint to record
    the analysis in the Job history. That way the COMPLETED ANALYSES
    dropdown in the UI is populated with *real* pocket data, not the
    hardcoded placeholders the legacy ``_ensure_pocket_output`` used to
    inject when Jobs auto-completed.
    """
    pdb_id = str(body.get("pdb_id") or "").strip().upper()
    target_name = str(body.get("target_name") or "").strip()
    pocket_result = body.get("pocket_result") or {}
    if not pdb_id and not target_name:
        raise HTTPException(status_code=400, detail="pdb_id or target_name required")

    # Pocket_App returns { "target": ..., "pocket": {...}, "pdb_id": ... } —
    # accept either the wrapped or unwrapped form for resilience.
    pocket_data = pocket_result.get("pocket") if isinstance(pocket_result, dict) else None
    if not isinstance(pocket_data, dict):
        pocket_data = pocket_result if isinstance(pocket_result, dict) else {}

    receptor_value = pdb_id or target_name or "UNKNOWN"
    repo = JobRepository(db)
    job = repo.create(
        job_id=str(uuid4()),
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=100,
        ),
    )
    # Real analysis already happened (sync) — mark completed immediately.
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned is not None:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()

    audit_payload = {
        "analysis_section": "pocket_analysis",
        "job_id": job.id,
        "pdb_id": pdb_id or "---",
        "target_name": target_name or None,
        "best_pocket": 1,
        "pockets": [
            {
                "pocket_number": 1,
                "druggability_score": pocket_data.get("druggability_score", 0.0),
                "source": pocket_data.get("source"),
                "ligand_resname": pocket_data.get("ligand_resname"),
                "ligand_count": 1 if pocket_data.get("ligand_resname") else 0,
                "residues": pocket_data.get("pocket_residues") or [],
                "center": [
                    pocket_data.get("center_x"),
                    pocket_data.get("center_y"),
                    pocket_data.get("center_z"),
                ],
                "box": [
                    pocket_data.get("size_x"),
                    pocket_data.get("size_y"),
                    pocket_data.get("size_z"),
                ],
                "notes": pocket_data.get("notes"),
            }
        ],
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="pocket_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(audit_payload, ensure_ascii=True),
        )
    )
    db.commit()
    return {"job_id": job.id, "status": "completed", "pdb_id": pdb_id}


@router.get("/jobs/denovo-source-list")
async def legacy_denovo_source_list(db: Session = Depends(get_db)):
    service = JobService(db)
    jobs = service.list_jobs()
    output_rows = (
        db.query(AuditLog)
        .filter(AuditLog.resource_type == "job_output", AuditLog.action == "pharmacophore_output")
        .all()
    )
    # Capture actor_id along with payload so each pharmacophore source entry
    # can declare its provenance ("UI button" vs "agent"). The Step 2 → Step
    # 3 chain dropdown can then show a 🤖 prefix for agent-produced
    # pharmacophore runs the same way Step 4 / Step 5 dropdowns do.
    output_map: dict[str, tuple[dict, str]] = {}
    for row in output_rows:
        try:
            output_map[row.resource_id] = (json.loads(row.payload or "{}"), row.actor_id or "system")
        except Exception:
            output_map[row.resource_id] = ({}, row.actor_id or "system")

    items: list[dict] = []
    for job in jobs:
        if job.status != "completed" or job.id not in output_map:
            continue
        pharma, actor_id = output_map[job.id]
        src = _source_label(actor_id)
        parent_job_id = str(pharma.get("parent_job_id") or "").strip()
        parent_pocket_payload = _pocket_payload_for_job(db, parent_job_id) if parent_job_id and parent_job_id != "0" else {}
        pockets = list(parent_pocket_payload.get("pockets") or [])
        designated = str(pharma.get("pocket_number") or "0").strip()
        residues = pharma.get("selected_residues") or []
        if not pockets:
            pockets = [
                {
                    "pocket_number": designated or "1",
                    "druggability_score": 1.0,
                    "ligand_count": 0,
                    "residues": residues if isinstance(residues, list) else [],
                    "lining_residues": residues if isinstance(residues, list) else [],
                }
            ]

        nice_target = _short_target_label(job.pdb_id or job.receptor_name, src)
        items.append(
            {
                "id": job.id,
                "job_id": job.id,
                "pdb_id": job.pdb_id or job.receptor_name or "---",
                "parent_job_id": parent_job_id if parent_job_id != "0" else None,
                "designated_pocket_number": designated,
                "selected_residues": residues if isinstance(residues, list) else [],
                "pockets": pockets,
                "display_name": (
                    f"[FARMAKOFOR OUTPUT] {nice_target} — "
                    f"Cep {designated or '?'}"
                ),
                "counts": {
                    "acceptor": pharma.get("acceptor", 0),
                    "donor": pharma.get("donor", 0),
                    "positive": pharma.get("positive", 0),
                    "negative": pharma.get("negative", 0),
                    "hydrophobic": pharma.get("hydrophobic", 0),
                    "aromatic": 0,
                },
                "date": job.created_at.strftime("%d.%m.%Y %H:%M"),
                "source": src,
                "short_id": job.id[:8],
            }
        )
    return items


def _short_target_label(pdb_id: str | None, source: str) -> str:
    """Render a friendly target label for dropdowns.

    UI-button runs use the raw PDB id (e.g. "4ZAU"). Agent runs use a long
    receptor_value like ``AGENT-DENOVO-EGFR-2070e038`` — the suffix is a
    DB tracking id, not researcher-meaningful. Strip the prefix so agent
    rows read e.g. "EGFR (De Novo)" rather than the wall-of-text UUID.
    """
    raw = (pdb_id or "").strip() or "Unknown"
    if source != "agent" or not raw.startswith("AGENT-"):
        return raw
    # AGENT-{KIND}-{TARGET}-{8charid} → take the {TARGET} segment
    parts = raw.split("-")
    if len(parts) >= 4:
        kind = parts[1].title()  # DENOVO → Denovo, ANALOG → Analog
        target = parts[2]
        return f"{target} ({kind})"
    return raw


@router.get("/jobs/admet-source-list")
async def legacy_admet_source_list(db: Session = Depends(get_db)):
    service = JobService(db)
    jobs = service.list_jobs()
    output_rows = (
        db.query(AuditLog)
        .filter(AuditLog.resource_type == "job_output", AuditLog.action == "denovo_output")
        .all()
    )
    # output_map now carries (payload, actor_id) so we can stamp each item
    # with a provenance hint the Step 4 ADMET dropdown can use to render
    # a 🤖 prefix for agent-driven runs.
    output_map: dict[str, tuple[dict, str]] = {}
    for row in output_rows:
        try:
            output_map[row.resource_id] = (json.loads(row.payload or "{}"), row.actor_id or "system")
        except Exception:
            output_map[row.resource_id] = ({}, row.actor_id or "system")

    items: list[dict] = []
    for job in jobs:
        if job.status != "completed" or job.id not in output_map:
            continue
        payload, actor_id = output_map[job.id]
        src = _source_label(actor_id)
        nice_target = _short_target_label(job.pdb_id or job.receptor_name, src)
        # Pass-count surfacing for Steps 4/5/6 dropdowns: derive the trio
        # from the De Novo output payload. `total_generated` = upstream
        # REINVENT/discovery total; `after_filter_count` = ADMET pass
        # count; `results_json.data` length = final candidates persisted.
        # All three may be missing on legacy demo / embedded-engine runs,
        # in which case the dropdown silently falls back to "no pass info"
        # rather than rendering garbage zeros.
        results_json = payload.get("results_json") or {}
        saved_rows = results_json.get("data") or []
        candidates_saved = len(saved_rows) if isinstance(saved_rows, list) else 0
        total_generated_raw = payload.get("total_generated")
        admet_passed_raw = payload.get("after_filter_count")
        items.append(
            {
                "job_id": job.id,
                "display_name": f"{nice_target} - De Novo Output",
                "date": job.created_at.timestamp(),
                "pocket_id": "1",
                "pdb_id_base": job.pdb_id or "6LU7",
                "residues": ["ALA205", "GLN203", "ILE207"],
                "druggability_score": 0.81,
                "candidate_count": candidates_saved,
                # The trio Steps 4/5 dropdowns format as
                # "EGFR — N ADMET ✓ / M" + the inline summary card.
                "total_generated": int(total_generated_raw) if isinstance(total_generated_raw, (int, float)) else None,
                "admet_passed": int(admet_passed_raw) if isinstance(admet_passed_raw, (int, float)) else None,
                "candidates_saved": candidates_saved or None,
                # Provenance hint for the frontend dropdown (3-char emoji
                # prefix vs default cleanliness). Same field name as the
                # /jobs/queue + /jobs/hierarchy responses for consistency.
                "source": src,
                "short_id": job.id[:8],
            }
        )
    return {"status": "success", "data": items}


# ---------------------------------------------------------------------------
# Iterative Mode endpoints — standalone ADMET + standalone Vina docking.
# These let the user re-analyze a completed De Novo job with different
# filter thresholds, or bring an external SMILES list (CSV / paste) for
# ADMET prediction and/or docking against a chosen receptor. Closes the
# "Step 4 / Step 5 buttons return Method Not Allowed" gap that previously
# made the UI 6-step workflow misleading.
#
# Workflow:
#   /api/jobs/log-* — record real Pocket_App results (sync path)
#   /api/admet/run  — re-ADMET a De Novo job OR external SMILES list
#   /api/docking/run — Vina dock a De Novo / ADMET job OR external SMILES,
#                      with optional target receptor override (cross-target
#                      / reverse-virtual-screening per Egemen's suggestion)
# ---------------------------------------------------------------------------


def _smiles_list_from_job_output(db: Session, job_id: str) -> list[str]:
    """Pull the SMILES list from a previous De Novo or ADMET output AuditLog.

    Tries denovo_output first (the standard De Novo result shape), then
    admet_output (for re-docking an ADMET-filtered set). Returns empty list
    if no usable output is found.
    """
    smiles, _ = _job_output_smiles_and_target(db, job_id)
    return smiles


def _job_output_smiles_and_target(db: Session, job_id: str) -> tuple[list[str], Optional[str]]:
    """Pull SMILES list + target_name from a job's denovo/admet output.

    The target_name (e.g. "EGFR") is what the bundled receptor list keys on
    for Vina docking — NOT the PDB structure id (e.g. "4ZAU") which the
    Job row's ``pdb_id`` column happens to store when both fields are
    supplied at De Novo time. Pulling the target from the denovo_output
    payload's ``target`` field is the source of truth.

    Returns (smiles_list, target_name_or_None). Empty / None on failure.
    """
    row = (
        db.query(AuditLog)
        .filter(
            AuditLog.resource_type == "job_output",
            AuditLog.resource_id == str(job_id),
            # docking_output added so Step 6 candidate-analysis (and any
            # cross-target re-dock) can pull SMILES off a finished Vina run.
            # All three writer endpoints (/api/jobs/start-denovo, /api/admet/run,
            # /api/docking/run) share the {results_json.data: [{smiles, ...}]}
            # shape, so the parser below handles all three uniformly.
            AuditLog.action.in_(["denovo_output", "admet_output", "docking_output"]),
        )
        .order_by(desc(AuditLog.created_at))
        .first()
    )
    if not row:
        return [], None
    try:
        payload = json.loads(row.payload or "{}")
    except Exception:
        return [], None
    # De Novo shape: results_json.data = [{smiles, ...}, ...]
    results = payload.get("results_json") or {}
    data = results.get("data") or []
    if not isinstance(data, list):
        data = []
    smiles_out: list[str] = []
    for item in data:
        if isinstance(item, dict):
            smi = item.get("smiles")
            if isinstance(smi, str) and smi.strip():
                smiles_out.append(smi.strip())
    # Target lives at multiple levels depending on writer; check both
    target = (
        payload.get("target")
        or results.get("target")
        or payload.get("target_name")
    )
    if isinstance(target, str):
        target = target.strip() or None
    else:
        target = None
    return smiles_out, target


@router.post("/admet/run")
async def iterative_admet_run(body: dict, db: Session = Depends(get_db)):
    """Standalone ADMET prediction — re-analyze a De Novo job or external SMILES.

    Body:
        {
          "job_id":  "<existing De Novo job id>"   # mode A: pull SMILES from job
          # OR
          "smiles":  ["CC(=O)Oc1...", ...]          # mode B: external list
          "user_id": "...",                         # optional, audit only
        }

    Returns:
        { "job_id": "<new ADMET job>", "status": "completed",
          "candidate_count": N }

    The new job's AuditLog (action="admet_output") payload mirrors the
    De Novo output shape so the UI Step 4 table renders the per-candidate
    ADMET descriptors without UI changes.
    """
    from app.services.admet_client import predict_admet

    source_job_id = str(body.get("job_id") or "").strip()
    raw_smiles = body.get("smiles") or []
    if not isinstance(raw_smiles, list):
        raw_smiles = []

    smiles_list: list[str] = []
    if source_job_id:
        smiles_list = _smiles_list_from_job_output(db, source_job_id)
        if not smiles_list:
            raise HTTPException(
                status_code=404,
                detail=f"No SMILES found in job_output for job_id={source_job_id!r}",
            )
    elif raw_smiles:
        # External list — sanitize entries
        smiles_list = [str(s).strip() for s in raw_smiles if isinstance(s, str) and str(s).strip()]
    else:
        raise HTTPException(
            status_code=400,
            detail="Provide either 'job_id' (existing De Novo) or 'smiles' (list of SMILES).",
        )

    if not smiles_list:
        raise HTTPException(status_code=400, detail="Empty SMILES list after sanitization.")

    # Call Admet_App
    try:
        admet_preds = await predict_admet(smiles_list)
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=502, detail=f"ADMET service failed: {exc}") from exc

    # Build per-candidate rows (UI Step 4 table reads these via /jobs/queue
    # → results_json.data). Keep the field names aligned with denovo_output
    # so the existing render path Just Works.
    #
    # ResultsModal.jsx Step 4 ADMET table reads item.mw, item.logp, item.tpsa,
    # item.herg directly off each row — NOT from item.descriptors. Compute
    # MW/LogP/TPSA via RDKit (same approach as /docking/run); derive a
    # categorical herg label from the ADMET-AI hERG raw probability so the
    # red/green badge actually has a value to render. Same fix shape that
    # solved the empty-column problem on the Docking results table.
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors, Crippen  # type: ignore
        _rdkit_ok = True
    except Exception:  # noqa: BLE001
        _rdkit_ok = False

    def _props_for(smi: str | None) -> tuple[float | None, float | None, float | None]:
        if not _rdkit_ok or not smi:
            return None, None, None
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None, None, None
            return (
                round(Descriptors.MolWt(mol), 1),
                round(Crippen.MolLogP(mol), 2),
                round(Descriptors.TPSA(mol), 1),
            )
        except Exception:  # noqa: BLE001
            return None, None, None

    candidates: list[dict] = []
    for idx, pred in enumerate(admet_preds, start=1):
        smi = pred.get("smiles") if isinstance(pred, dict) else None
        raw = (pred.get("raw") if isinstance(pred, dict) else None) or {}
        pct = (pred.get("percentile") if isinstance(pred, dict) else None) or {}
        descriptors = {
            "hERG": raw.get("hERG"),
            "DILI": raw.get("DILI"),
            "AMES": raw.get("AMES"),
            "BBB": raw.get("BBB_Martins"),
            "Bioavailability": raw.get("Bioavailability_Ma"),
            "Caco2": raw.get("Caco2_Wang"),
            "PPBR": raw.get("PPBR_AZ"),
            "CYP3A4_inhib": raw.get("CYP3A4_Veith"),
            "CYP2D6_inhib": raw.get("CYP2D6_Veith"),
            "CYP2C9_inhib": raw.get("CYP2C9_Veith"),
            "Half_Life_h": raw.get("Half_Life_Obach"),
        }
        mw, logp, tpsa = _props_for(smi)
        # ADMET-AI hERG raw score is a probability of blockade [0,1].
        # >0.5 is the standard cutoff in TDC for "cardiotoxic". Render as
        # "VAR" / "YOK" so the existing red/green badge styling lights up.
        herg_raw = raw.get("hERG")
        if isinstance(herg_raw, (int, float)):
            herg_label = "VAR" if float(herg_raw) > 0.5 else "YOK"
        else:
            herg_label = None
        candidates.append(
            {
                "rank": idx,
                "smiles": smi,
                "module": "general",
                "score": None,
                "dock_score": None,
                # Top-level fields the Step 4 table renderer reads.
                "mw": mw,
                "logp": logp,
                "tpsa": tpsa,
                "herg": herg_label,
                "alerts": [],
                "component_scores": {},
                "descriptors": {
                    **{k: v for k, v in descriptors.items() if v is not None},
                    **({"mw": mw} if mw is not None else {}),
                    **({"logp": logp} if logp is not None else {}),
                    **({"tpsa": tpsa} if tpsa is not None else {}),
                },
                "admet_raw": raw,
                "admet_percentile": pct,
            }
        )

    # Create the new job
    source_label = source_job_id[:8] if source_job_id else "EXTERNAL"
    receptor_value = f"ADMET-RERUN-{source_label}"
    repo = JobRepository(db)
    job = repo.create(
        job_id=str(uuid4()),
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=max(1, len(candidates)),
        ),
    )
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned is not None:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()

    admet_result = {
        "analysis_section": "admet",
        "source_job_id": source_job_id or None,
        "source_kind": "denovo_job" if source_job_id else "external_smiles",
        "input_count": len(smiles_list),
        "output_count": len(candidates),
        "results_json": {
            "section": "admet",
            "total_input": len(smiles_list),
            "total_output": len(candidates),
            "data": candidates,
        },
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="admet_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(admet_result, ensure_ascii=True),
        )
    )
    db.commit()
    return {
        "job_id": job.id,
        "status": "completed",
        "candidate_count": len(candidates),
        "source": "denovo_job" if source_job_id else "external_smiles",
    }


@router.post("/docking/run")
async def iterative_docking_run(body: dict, db: Session = Depends(get_db)):
    """Standalone Vina docking — dock candidates from a previous job or
    external SMILES list against a chosen bundled receptor.

    Body:
        {
          "job_id":  "<existing De Novo / ADMET job id>"   # mode A
          # OR
          "smiles":  ["...", ...]                           # mode B (external)
          "target":  "EGFR"                                 # which bundled
                                                            # receptor to dock
                                                            # against. Defaults
                                                            # to the source
                                                            # job's pdb_id.
          "pocket_id": "1"                                  # (informational)
        }

    Returns:
        { "job_id": "<new docking job>", "status": "completed",
          "docked_count": N, "target": "EGFR", "succeeded": M }

    Target override enables **cross-target / reverse virtual screening**:
    dock an EGFR-designed candidate against CYP2C9 / JAK2 / etc. to probe
    off-target toxicity or selectivity. (Egemen's #2 recommendation —
    reverse virtual screening for toxicity targets.)
    """
    from app.services.docking_client import dock_batch, list_available_targets

    source_job_id = str(body.get("job_id") or "").strip()
    raw_smiles = body.get("smiles") or []
    if not isinstance(raw_smiles, list):
        raw_smiles = []
    target_override = str(body.get("target") or "").strip()

    smiles_list: list[str] = []
    inferred_target: str | None = None
    if source_job_id:
        # Pull both SMILES list and the target_name from the source job's
        # denovo_output payload. The payload's "target" field stores the
        # bundled receptor key (e.g. "EGFR") which is what Vina dispatches on
        # — NOT the Job row's pdb_id (e.g. "4ZAU") which is just the PDB
        # structure id and won't match any bundled PDBQT.
        smiles_list, inferred_target = _job_output_smiles_and_target(db, source_job_id)
        if not smiles_list:
            raise HTTPException(
                status_code=404,
                detail=f"No SMILES found in job_output for job_id={source_job_id!r}",
            )
        # Fallback to Job row's stored receptor_name / pdb_id when the
        # payload didn't carry a target — older denovo runs may not have
        # populated the field. The pdb_id-only path will still fail bundled-
        # receptor validation below, which is the correct behavior (asks user
        # to pick a target explicitly).
        if not inferred_target:
            source_job = db.query(Job).filter(Job.id == source_job_id).first()
            if source_job is not None:
                inferred_target = source_job.receptor_name or source_job.pdb_id
    elif raw_smiles:
        smiles_list = [str(s).strip() for s in raw_smiles if isinstance(s, str) and str(s).strip()]
    else:
        raise HTTPException(
            status_code=400,
            detail="Provide either 'job_id' (existing job) or 'smiles' (list of SMILES).",
        )

    if not smiles_list:
        raise HTTPException(status_code=400, detail="Empty SMILES list after sanitization.")

    # Resolve docking target — explicit override wins, then inferred, then 400.
    docking_target = (target_override or inferred_target or "").strip()
    if not docking_target:
        raise HTTPException(
            status_code=400,
            detail="No docking target — provide 'target' field or pass a job_id with a known pdb_id.",
        )

    # Validate against bundled receptors (Vina needs a PDBQT file)
    try:
        available = await list_available_targets()
    except Exception:  # noqa: BLE001
        available = []
    case_match = next((t for t in available if t.lower() == docking_target.lower()), None)
    if case_match is None:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Target {docking_target!r} is not in bundled receptors. "
                f"Available: {sorted(available)}. "
                "Add the receptor PDBQT to Docking_App/receptors/ to enable docking."
            ),
        )
    docking_target = case_match

    # Run Vina (top 30 by default to keep wall-clock reasonable; user can
    # request larger by submitting a smaller SMILES list)
    DOCK_CAP = 30
    docking_input = smiles_list[:DOCK_CAP]
    try:
        dock_results = await dock_batch(docking_input, docking_target)
    except HTTPException:
        raise
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=502, detail=f"Docking service failed: {exc}") from exc

    # Build candidate rows the UI Step 5 table can render
    by_smi: dict[str, float | None] = {}
    # Parallel pose-PDB map (smiles -> PDB block) from Docking_App's Phase B
    # response shape. Stored per-row below so the NGL viewer can fetch the
    # real Vina pose via /api/jobs/{id}/poses/{rank}.pdb instead of falling
    # back to the RDKit free-state conformer.
    pose_by_smi: dict[str, str] = {}
    # Multi-pose map (all Vina MODELs in one PDB) — used by NGL's pose
    # cycler at /api/jobs/{id}/poses/{rank}.pdb?mode=all
    pose_all_by_smi: dict[str, str] = {}
    for r in dock_results or []:
        if isinstance(r, dict):
            smi = r.get("smiles")
            aff = r.get("affinity")
            if isinstance(smi, str) and isinstance(aff, (int, float)):
                by_smi[smi] = float(aff)
            if isinstance(smi, str):
                pose = r.get("pose_pdb")
                if isinstance(pose, str) and pose.strip():
                    pose_by_smi[smi] = pose
                pose_all = r.get("pose_pdb_all")
                if isinstance(pose_all, str) and pose_all.strip():
                    pose_all_by_smi[smi] = pose_all

    # ResultsModal.jsx (lines 380-384) reads item.mw, item.logp, item.affinity
    # directly off each row — not from descriptors{}. Compute them via RDKit
    # so the Docking results table actually shows numbers instead of blanks.
    # (The user just saw all three columns empty despite "5/5 succeeded".)
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors, Crippen  # type: ignore
        _rdkit_ok = True
    except Exception:  # noqa: BLE001
        _rdkit_ok = False

    def _props_for(smi: str) -> tuple[float | None, float | None]:
        if not _rdkit_ok:
            return None, None
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None, None
            return round(Descriptors.MolWt(mol), 1), round(Crippen.MolLogP(mol), 2)
        except Exception:  # noqa: BLE001
            return None, None

    candidates: list[dict] = []
    for idx, smi in enumerate(docking_input, start=1):
        aff = by_smi.get(smi)
        mw, logp = _props_for(smi)
        aff_rounded = round(aff, 3) if aff is not None else None
        row: dict = {
            "rank": idx,
            "smiles": smi,
            "module": "general",
            "score": None,
            "dock_score": aff_rounded,
            # Mirror dock_score into `affinity` — ResultsModal.jsx Step 5
            # column reads item.affinity directly.
            "affinity": aff_rounded,
            "mw": mw,
            "logp": logp,
            "alerts": [],
            "component_scores": {},
            "descriptors": {
                "vina_affinity_kcal_mol": aff_rounded,
                "mw": mw,
                "logp": logp,
            },
        }
        pose = pose_by_smi.get(smi)
        if pose:
            row["pose_pdb"] = pose
        pose_all = pose_all_by_smi.get(smi)
        if pose_all:
            row["pose_pdb_all"] = pose_all
        candidates.append(row)
    candidates.sort(key=lambda c: (c["dock_score"] is None, c["dock_score"] or 0))
    # Re-stamp rank after the sort so the persisted rank matches the row
    # order (NGL viewer asks for /poses/{rank}.pdb where rank is 1-based
    # and resolves to the same row the user clicked in the table).
    for new_rank, c in enumerate(candidates, start=1):
        c["rank"] = new_rank

    # Create the new job + AuditLog
    source_label = source_job_id[:8] if source_job_id else "EXTERNAL"
    receptor_value = f"DOCK-{docking_target}-{source_label}"
    repo = JobRepository(db)
    job = repo.create(
        job_id=str(uuid4()),
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=max(1, len(candidates)),
        ),
    )
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned is not None:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()

    succeeded = sum(1 for c in candidates if c["dock_score"] is not None)
    docking_result = {
        "analysis_section": "docking",
        "source_job_id": source_job_id or None,
        "source_kind": "denovo_job" if source_job_id else "external_smiles",
        "docking_target": docking_target,
        "cross_target": bool(target_override and inferred_target
                             and target_override.lower() != (inferred_target or "").lower()),
        "input_count": len(smiles_list),
        "docked_count": len(docking_input),
        "succeeded": succeeded,
        "results_json": {
            "section": "docking",
            "target": docking_target,
            "input_count": len(smiles_list),
            "docked_count": len(docking_input),
            "succeeded": succeeded,
            "data": candidates,
        },
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="docking_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(docking_result, ensure_ascii=True),
        )
    )
    db.commit()
    return {
        "job_id": job.id,
        "status": "completed",
        "target": docking_target,
        "docked_count": len(docking_input),
        "succeeded": succeeded,
        "cross_target": docking_result["cross_target"],
    }


@router.post("/jobs/start-candidate-analysis")
async def legacy_start_candidate_analysis(
    analysis_section: str = Form(default="candidate_analysis"),
    parent_job_id: str = Form(default="0"),
    smiles_input: str = Form(default=""),
    db: Session = Depends(get_db),
):
    """Step 6 — Single-molecule deep-dive.

    Modes:
      A) ``smiles_input`` non-empty → analyze that SMILES literally.
      B) ``parent_job_id`` set + ``smiles_input`` empty → pull the top
         candidate from the parent job (De Novo or Docking output, the only
         two analysis sections the frontend list filter exposes for Step 6).

    What "deep-dive" means here:
      - RDKit descriptors: MW, logP, TPSA, HBD/HBA, rotatable bonds, formal
        charge, QED, SA score
      - Lipinski Rule of Five pass/fail breakdown
      - PAINS hard-flag + BRENK soft-flag substructure alerts
      - ADMET-AI prediction (all 52 endpoints) for the one molecule

    The output is stored under action ``candidate_analysis_output`` so the
    Analysis Results page's "Candidate Molecule Analysis" tab can render it
    without changes to the queue/hierarchy plumbing. Returns
    ``{ job_id, status, smiles, lipinski_pass, alert_count }`` — the
    frontend already calls ``alert(t('signal_sent_success', { jobId }))``
    on success.
    """
    from app.services.admet_client import predict_admet
    from app.services.production_pipeline import (
        _extra_rdkit_descriptors,
        _structural_alerts,
    )

    # ---- 1. Resolve the SMILES to analyze -----------------------------------
    smi_raw = (smiles_input or "").strip()
    parent_id = (parent_job_id or "").strip()
    parent_pdb: str | None = None

    if smi_raw:
        # User pasted SMILES — split on whitespace/comma/newline and take
        # the first non-empty token. (The textarea accepts free text.)
        for tok in smi_raw.replace(",", "\n").split():
            tok = tok.strip()
            if tok:
                smi_raw = tok
                break
    elif parent_id and parent_id != "0":
        pulled, _t = _job_output_smiles_and_target(db, parent_id)
        if not pulled:
            raise HTTPException(
                status_code=404,
                detail=f"No SMILES found in parent job_output for parent_job_id={parent_id!r}",
            )
        smi_raw = pulled[0]  # top-ranked candidate (audit is rank-ordered)
        parent = db.query(Job).filter(Job.id == parent_id).first()
        if parent is not None:
            parent_pdb = parent.pdb_id or parent.receptor_name
    else:
        raise HTTPException(
            status_code=400,
            detail="Provide either 'smiles_input' (manual SMILES) or 'parent_job_id' (existing De Novo / docking job).",
        )

    # ---- 2. RDKit validation + descriptors ----------------------------------
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors, Crippen, Lipinski, QED  # type: ignore
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"RDKit not available: {exc}") from exc

    mol = Chem.MolFromSmiles(smi_raw)
    if mol is None:
        raise HTTPException(
            status_code=400,
            detail=f"Could not parse SMILES: {smi_raw!r}",
        )

    canonical = Chem.MolToSmiles(mol)
    try:
        mw = float(Descriptors.MolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))
        qed = float(QED.qed(mol))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Descriptor calc failed: {exc}") from exc

    # Lipinski Rule of Five — explicit per-criterion breakdown so the UI can
    # explain *why* a molecule fails (vs. a single bool).
    lipinski_checks = {
        "mw_le_500": mw <= 500.0,
        "logp_le_5": logp <= 5.0,
        "hbd_le_5": hbd <= 5,
        "hba_le_10": hba <= 10,
    }
    lipinski_pass = all(lipinski_checks.values())
    # Veber rules (oral bioavailability)
    veber_checks = {
        "rotb_le_10": rotb <= 10,
        "tpsa_le_140": tpsa <= 140.0,
    }
    veber_pass = all(veber_checks.values())

    extras = _extra_rdkit_descriptors(canonical)  # sa_score + formal_charge
    hard_alerts, soft_alerts = _structural_alerts(canonical)

    # ---- 3. ADMET prediction (single molecule) ------------------------------
    admet_raw: dict = {}
    admet_pct: dict = {}
    try:
        preds = await predict_admet([canonical])
        if preds and isinstance(preds[0], dict):
            admet_raw = preds[0].get("raw") or {}
            admet_pct = preds[0].get("percentile") or {}
    except HTTPException:
        # ADMET service errors are non-fatal here — the descriptor / alert
        # block is still useful on its own. Surface the failure in the row.
        admet_raw = {}
    except Exception:  # noqa: BLE001
        admet_raw = {}

    # ---- 4. Persist as a completed Job + audit row --------------------------
    parent_label = parent_id[:8] if parent_id and parent_id != "0" else "MANUAL"
    receptor_value = f"CANDIDATE-{parent_label}"
    repo = JobRepository(db)
    job = repo.create(
        job_id=str(uuid4()),
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=1,
        ),
    )
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned is not None:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()

    row = {
        "rank": 1,
        "smiles": canonical,
        "input_smiles": smi_raw if smi_raw != canonical else None,
        "module": "general",
        "score": None,
        "dock_score": None,
        "mw": round(mw, 2),
        "logp": round(logp, 3),
        "tpsa": round(tpsa, 2),
        "qed": round(qed, 4),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "sa_score": (round(extras["sa_score"], 2) if "sa_score" in extras else None),
        "formal_charge_abs": (int(extras["formal_charge_abs"]) if "formal_charge_abs" in extras else None),
        "lipinski_pass": lipinski_pass,
        "lipinski_checks": lipinski_checks,
        "veber_pass": veber_pass,
        "veber_checks": veber_checks,
        "alerts": hard_alerts,
        "soft_alerts": soft_alerts,
        "admet_raw": admet_raw,
        "admet_percentile": admet_pct,
        "descriptors": {
            "mw": round(mw, 2),
            "logp": round(logp, 3),
            "qed": round(qed, 4),
            "tpsa": round(tpsa, 2),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": rotb,
        },
        "parent_job_id": parent_id if parent_id != "0" else None,
        "parent_pdb": parent_pdb,
    }
    result = {
        "analysis_section": "candidate_analysis",
        "parent_job_id": parent_id if parent_id != "0" else None,
        "smiles": canonical,
        "lipinski_pass": lipinski_pass,
        "veber_pass": veber_pass,
        "alert_count": len(hard_alerts),
        "results_json": {
            "section": "candidate_analysis",
            "total_input": 1,
            "total_output": 1,
            "data": [row],
        },
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="candidate_analysis_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(result, ensure_ascii=True),
        )
    )
    db.commit()
    return {
        "job_id": job.id,
        "status": "completed",
        "smiles": canonical,
        "lipinski_pass": lipinski_pass,
        "veber_pass": veber_pass,
        "alert_count": len(hard_alerts),
    }


@router.post("/jobs/start-analysis")
async def legacy_start_analysis(
    analysis_section: str = Form(default="pocket_analysis"),
    pdb_id: str | None = Form(default=None),
    entry_type: str = Form(default="PDB ID"),
    file: UploadFile | None = File(default=None),
    db: Session = Depends(get_db),
):
    receptor_value = (pdb_id or "").strip().upper() or (file.filename if file else "UPLOADED")
    payload = JobCreateRequest(
        tenant_id="default",
        project_id="legacy-ui",
        pdb_id=receptor_value,
        receptor_name=receptor_value,
        mode="agonist",
        candidate_count=100,
    )
    service = JobService(db)
    job = service.create_job(payload)
    return {
        "job_id": job.id,
        "status": "processing",
        "pdb_id": receptor_value,
        "membership": "gold",
        "analysis_section": analysis_section,
        "entry_type": entry_type,
    }


@router.post("/jobs/start-pharmacophore")
async def legacy_start_pharmacophore(
    parent_job_id: str = Form(default="0"),
    pocket_number: str = Form(default="0"),
    selected_residues: str = Form(default="[]"),
    db: Session = Depends(get_db),
):
    receptor_value = "PHARMA-TARGET"
    if parent_job_id and parent_job_id != "0":
        base_job = db.query(Job).filter(Job.id == parent_job_id).first()
        if base_job:
            receptor_value = base_job.pdb_id or base_job.receptor_name or receptor_value
    payload = JobCreateRequest(
        tenant_id="default",
        project_id="legacy-ui",
        pdb_id=receptor_value,
        receptor_name=receptor_value,
        mode="agonist",
        candidate_count=100,
    )
    # Farmakoför: gerçek zamanlı pipeline gerekmez; create_job Celery/8s simülasyonuna
    # sokup işi "completed" olmadan kilitlemesin (De Novo liste + sonuçlar ekranı için anında hazır).
    repo = JobRepository(db)
    job = repo.create(job_id=str(uuid4()), payload=payload)
    try:
        residue_list = json.loads(selected_residues or "[]")
    except Exception:
        residue_list = []

    # Manual mode means the user picked pocket-lining residues by name (e.g.
    # ["ALA12", "LYS45", ...]) without a bound co-crystal ligand. A genuine
    # pharmacophore is the 3D arrangement of chemical features a *ligand*
    # must present — we don't have one here, so we can't compute it.
    # Instead we derive residue-side-chain hints from the canonical 3-letter
    # codes: acidic side chains imply H-bond ACCEPTOR opportunities at the
    # pocket wall, basic ones imply DONOR opportunities, etc. This is an
    # honest residue-type tally — NOT a real ligand-derived pharmacophore —
    # and we label it as such in the saved payload so downstream consumers
    # don't mistake it for the proper Pocket_App output.
    ACIDIC = {"ASP", "GLU"}
    BASIC = {"LYS", "ARG", "HIS"}
    AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
    HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}

    def _resname(s: object) -> str:
        if not isinstance(s, str):
            return ""
        # Strip residue numbers ("ALA12" → "ALA") and uppercase.
        out = "".join(ch for ch in s if ch.isalpha()).upper()
        return out[:3]

    res3 = [_resname(r) for r in residue_list]
    acceptor_count = sum(1 for r in res3 if r in ACIDIC)
    donor_count = sum(1 for r in res3 if r in BASIC)
    aromatic_count = sum(1 for r in res3 if r in AROMATIC)
    hydrophobic_count = sum(1 for r in res3 if r in HYDROPHOBIC)

    pharma_result = {
        "analysis_section": "farmakofor",
        "parent_job_id": parent_job_id,
        "pocket_number": pocket_number,
        # Residue-type tally — labelled clearly as a hint, not a pharmacophore.
        "acceptor": acceptor_count,
        "donor": donor_count,
        "hydrophobic": hydrophobic_count,
        "aromatic": aromatic_count,
        "positive": donor_count,  # basic side chains are positively charged at pH 7
        "negative": acceptor_count,  # acidic side chains are negatively charged at pH 7
        "selected_residues": residue_list,
        "source": "manual_residue_hint",
        "note": (
            "Residue-type hint from selected pocket residues. NOT a ligand-"
            "derived pharmacophore — use 'Extract from co-crystal (PDB)' for "
            "the genuine 3D pharmacophore from Pocket_App."
        ),
        "results_json": {
            "acceptor": acceptor_count,
            "donor": donor_count,
            "hydrophobic": hydrophobic_count,
            "aromatic": aromatic_count,
            "positive": donor_count,
            "negative": acceptor_count,
            "residue_count": len(residue_list),
            "source": "manual_residue_hint",
        },
    }
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="pharmacophore_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(pharma_result, ensure_ascii=True),
        )
    )
    db.commit()
    return {"job_id": job.id, "status": "completed", "analysis_section": "farmakofor"}


# POST /api/jobs/start-denovo and POST {api_prefix}/jobs/start-denovo are registered on the FastAPI app
# in main.py *before* routers so no /jobs/{job_id} route can shadow them (405).
async def legacy_start_denovo(payload: dict, db: Session):
    disc = payload.get("discovery") if isinstance(payload.get("discovery"), dict) else payload
    target_name = str(disc.get("target_name") or "").strip()
    if len(target_name) < 2:
        raise HTTPException(status_code=400, detail="target_name must be at least 2 characters")

    # Production mode (REINVENT 4 + ADMET-AI + Vina + ChEMBL Tanimoto) does NOT
    # need a pharmacophore/pocket prerequisite — target enrichment is fetched
    # automatically from ChEMBL by ``target_name``.
    requested_engine = str(disc.get("engine_mode") or "").strip().lower()
    if requested_engine not in {"", "production", "demo"}:
        raise HTTPException(status_code=422, detail="engine_mode must be production or demo")
    is_production = requested_engine != "demo"

    pharma_job_id = str(disc.get("pharmacophore_job_id") or "").strip()
    if not is_production and not pharma_job_id:
        raise HTTPException(status_code=400, detail="pharmacophore_job_id required — select a completed pharmacophore analysis")

    pocket_number_raw = disc.get("pocket_number")
    pocket_num = _normalize_pocket_token(pocket_number_raw)
    if not is_production and not pocket_num:
        raise HTTPException(status_code=400, detail="pocket_number required — select the binding pocket for De Novo")

    # Skip the pharma DB lookups + pocket validation entirely for production mode.
    pharma = None
    pharma_row_job = None
    parent_job_id = ""
    pockets = []
    designated = ""
    residues = []
    chosen_pocket_detail = None

    if not is_production:
        pharma = _pharmacophore_payload_for_job(db, pharma_job_id)
        if not pharma:
            raise HTTPException(status_code=404, detail="Pharmacophore output not found for the selected job")

        pharma_row_job = db.query(Job).filter(Job.id == pharma_job_id).first()
        if not pharma_row_job or pharma_row_job.status != "completed":
            raise HTTPException(status_code=404, detail="Pharmacophore job not found or not completed")

        parent_job_id = str(pharma.get("parent_job_id") or "").strip()
        parent_pocket_payload = _pocket_payload_for_job(db, parent_job_id) if parent_job_id and parent_job_id != "0" else {}
        pockets = list(parent_pocket_payload.get("pockets") or [])
        designated = str(pharma.get("pocket_number") or "0").strip()
        residues = pharma.get("selected_residues") or []
        if not pockets:
            pockets = [
                {
                    "pocket_number": designated or "1",
                    "druggability_score": 1.0,
                    "residues": residues if isinstance(residues, list) else [],
                    "lining_residues": residues if isinstance(residues, list) else [],
                }
            ]

        valid_pockets = {
            _normalize_pocket_token(p.get("pocket_number"))
            for p in pockets
            if p.get("pocket_number") is not None and _normalize_pocket_token(p.get("pocket_number"))
        }
        if pocket_num not in valid_pockets:
            raise HTTPException(
                status_code=400,
                detail=f"Pocket {pocket_num} is not valid for this analysis. Choose one of: {sorted(valid_pockets)}",
            )

        chosen_pocket_detail = next(
            (p for p in pockets if _normalize_pocket_token(p.get("pocket_number")) == pocket_num),
            None,
        )

    try:
        gen_count = int(disc.get("generate_count", 400))
    except (TypeError, ValueError):
        gen_count = 400
    gen_count = max(50, min(gen_count, 3000))

    pdb_raw = disc.get("pdb_id")
    pdb_val = pdb_raw.strip().upper() if isinstance(pdb_raw, str) and pdb_raw.strip() else pdb_raw
    if isinstance(pdb_val, str) and not pdb_val:
        pdb_val = None
    if not pdb_val and pharma_row_job is not None:
        pdb_val = (pharma_row_job.pdb_id or pharma_row_job.receptor_name or "").strip().upper() or None

    seq_raw = disc.get("sequence")
    seq_val = seq_raw.strip() if isinstance(seq_raw, str) else seq_raw

    seed_raw = disc.get("smiles_seed")
    seed_val = seed_raw.strip() if isinstance(seed_raw, str) else seed_raw

    if is_production:
        # No pharmacophore context in production mode — pipeline uses ChEMBL
        # actives for target relevance instead.
        pharmacophore_context = None
    else:
        pharmacophore_context = {
            "pharmacophore_job_id": pharma_job_id,
            "pocket_number": pocket_num,
            "parent_job_id": parent_job_id if parent_job_id != "0" else None,
            "designated_pocket_at_pharma": designated,
            "selected_residues": residues if isinstance(residues, list) else [],
            "pharmacophore_features": {
                "acceptor": pharma.get("acceptor"),
                "donor": pharma.get("donor"),
                "hydrophobic": pharma.get("hydrophobic"),
                "positive": pharma.get("positive"),
                "negative": pharma.get("negative"),
            },
            "pocket_druggability_score": (chosen_pocket_detail or {}).get("druggability_score"),
            "pocket_residues": (chosen_pocket_detail or {}).get("residues")
            or (chosen_pocket_detail or {}).get("lining_residues")
            or [],
        }

    denovo_request: dict = {
        "target_name": target_name,
        "pdb_id": pdb_val or None,
        "sequence": seq_val or None,
        "smiles_seed": seed_val or None,
        "disease_module": str(disc.get("disease_module") or "general"),
        "desired_route": str(disc.get("desired_route") or "oral"),
        "generate_count": gen_count,
        "pharmacophore_context": pharmacophore_context,
    }
    if is_production:
        # Carry the engine override forward so _discover_upstream picks the
        # production_pipeline branch (REINVENT 4 + ADMET-AI + ChEMBL + Vina).
        denovo_request["engine_mode"] = "production"
    # Forward force_generic so production_pipeline can bypass the target
    # validation gate when the user explicitly accepts a target-agnostic
    # run (frontend prompts for confirmation after a 422 from the gate).
    if disc.get("force_generic"):
        denovo_request["force_generic"] = True

    fos = disc.get("filter_overrides")
    if isinstance(fos, dict) and fos:
        denovo_request["filter_overrides"] = fos

    try:
        disc_req = DiscoveryRequest.model_validate(denovo_request)
    except ValidationError as exc:
        raise HTTPException(status_code=422, detail=exc.errors()) from exc
    cfg_denovo = get_settings()
    scripts_dir = Path(__file__).resolve().parents[3] / "yeni_motor" / "scripts"
    script_path = scripts_dir / "run_universal_discovery_stdin.cjs"
    # Production mode never goes through the synthetic Node subprocess; it must
    # reach motor_discover_impl → _discover_upstream → run_production_discovery.
    use_client_script = (not is_production) and bool(cfg_denovo.denovo_use_yeni_motor_script and script_path.is_file())
    if cfg_denovo.denovo_use_yeni_motor_script and not script_path.is_file():
        logger.warning(
            "denovo_use_yeni_motor_script enabled but Node runner missing: %s — using motor_discover_impl",
            script_path,
        )

    try:
        if use_client_script:
            from app.services.yeni_motor_discover_subprocess import run_discovery_via_yeni_motor_script

            upstream = await run_discovery_via_yeni_motor_script(disc_req.model_dump(exclude_none=False))
            try:
                DiscoveryResponse.model_validate(upstream)
            except ValidationError as exc:
                raise HTTPException(status_code=502, detail=exc.errors()) from exc
        else:
            disc_resp = await motor_discover_impl(disc_req)
            upstream = disc_resp.model_dump()
    except HTTPException:
        raise
    meta = upstream.get("meta") or {}
    if not isinstance(meta, dict):
        meta = {}
    discovery_engine = str(
        meta.get("discovery_engine") or ("universal_discovery_engine" if use_client_script else "embedded")
    )
    is_synthetic = not is_production
    # Keep provenance attached to persisted outputs, not just the ephemeral API response.
    meta.update({"is_synthetic": is_synthetic,
                 "engine_type": "synthetic_demo" if is_synthetic else "computational_prediction",
                 "research_validated": False})
    meta.setdefault("run_id", str(uuid4()))
    meta.setdefault("reference_hash", None)
    meta.setdefault("model_versions", {})
    meta.setdefault("degraded", bool(is_synthetic or (meta.get("validation") or {}).get("degraded", False)))
    candidates = upstream.get("candidates", []) or []
    # Persist full candidate rows (same shape as yeni_motor /api/discover) for UI: QED, MW, logP in descriptors.
    normalized = [
        {
            "rank": item.get("rank", idx + 1),
            "smiles": item.get("smiles"),
            "module": item.get("module"),
            "score": item.get("score"),
            "dock_score": item.get("dock_score"),
            "alerts": item.get("alerts", []),
            "component_scores": item.get("component_scores") or {},
            "descriptors": item.get("descriptors") or {},
        }
        for idx, item in enumerate(candidates)
    ]

    display_receptor = (denovo_request.get("pdb_id") or target_name)[:255]
    job_payload = JobCreateRequest(
        tenant_id="default",
        project_id="legacy-ui",
        pdb_id=display_receptor,
        receptor_name=display_receptor,
        mode="agonist",
        candidate_count=max(1, len(normalized)),
    )
    service = JobService(db)
    job = service.create_job(job_payload)

    denovo_result = {
        "analysis_section": "denovo",
        "pharmacophore_job_id": pharma_job_id,
        "pocket_number": pocket_num,
        "pharmacophore_context": pharmacophore_context,
        "target": upstream.get("target"),
        "module": upstream.get("module"),
        "total_generated": upstream.get("total_generated"),
        "after_filter_count": upstream.get("after_filter_count"),
        "top_n": upstream.get("top_n"),
        "results_json": {
            "target": upstream.get("target"),
            "module": upstream.get("module"),
            "total_generated": upstream.get("total_generated"),
            "after_filter_count": upstream.get("after_filter_count"),
            "top_n": upstream.get("top_n"),
            "data": normalized,
            "meta": meta,
        },
        "meta": meta,
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="denovo_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(denovo_result, ensure_ascii=True),
        )
    )
    db.commit()
    # Forward the production pipeline's validation block so the UI banner
    # can warn the researcher when a run is degraded (e.g. force_generic=true
    # path). Missing on demo/embedded engines — degraded defaults to False.
    validation_block = (upstream.get("meta") or {}).get("validation") or {}
    return {
        "job_id": job.id,
        "status": "processing",
        "analysis_section": "denovo",
        "discovery_engine": discovery_engine,
        "is_synthetic": is_synthetic,
        "discovery_summary": {
            "target": upstream.get("target"),
            "module": upstream.get("module"),
            "total_generated": upstream.get("total_generated"),
            "after_filter_count": upstream.get("after_filter_count"),
            "top_n": upstream.get("top_n"),
            "candidates_saved": len(normalized),
            "validation": validation_block,
        },
    }


@router.delete("/jobs/{job_id}")
async def legacy_delete_job(job_id: UUID, db: Session = Depends(get_db)):
    """Remove a job and its output audit rows so it disappears from the queue.

    `job_id` must be a UUID so this route never shadows ``POST /jobs/start-denovo`` (405).
    """
    sid = str(job_id)
    service = JobService(db)
    if not service.delete_job(sid):
        raise HTTPException(status_code=404, detail="Job not found")
    return {"ok": True, "job_id": sid}


@router.get("/user/quota-status")
async def legacy_quota_status(user_id: str, feature: str):
    return {
        "user_id": user_id,
        "feature": feature,
        "limit": "Sınırsız",
        "used": 0,
        "remaining": 9999,
        "can_start": True,
        "status": "ok",
    }


@router.post("/auth/register")
async def legacy_auth_register(payload: dict, local_kw: str = Query(default="tr"), db: Session = Depends(get_db)):
    _ = local_kw
    email = (payload.get("email") or "").strip().lower()
    full_name = (payload.get("full_name") or "").strip()
    password = payload.get("password") or ""
    if not email or not password:
        raise HTTPException(status_code=400, detail="email and password required")
    existing = db.query(AuthUser).filter(AuthUser.email == email).first()
    if existing:
        raise HTTPException(status_code=409, detail="User already exists")
    try:
        password_hash = hash_password(password)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    user = AuthUser(id=str(uuid4()), email=email, full_name=full_name, password=password_hash)
    db.add(user)
    # Keep admin/user management screen in sync.
    managed = db.query(UserAccount).filter(UserAccount.email == email).first()
    if not managed:
        db.add(
            UserAccount(
                id=str(uuid4()),
                email=email,
                full_name=full_name,
                is_admin=False,
                pocket_quota=10,
                denovo_quota=10,
                admet_quota=10,
                docking_quota=10,
            )
        )
    db.commit()
    return {"ok": True, "email": email}


@router.post("/auth/login")
async def legacy_auth_login(payload: dict, request: Request, db: Session = Depends(get_db)):
    email = (payload.get("username") or payload.get("email") or "").strip().lower()
    password = payload.get("password") or ""
    check_login_throttle(request, email)
    # Allow configured platform admin to log in from legacy auth screen.
    import hmac
    if (settings.admin_password and settings.session_signing_key and
        email == settings.admin_username.lower() and
        hmac.compare_digest(password, settings.admin_password)):
        managed = db.query(UserAccount).filter(UserAccount.email == "admin@local").first()
        if not managed:
            db.add(
                UserAccount(
                    id=str(uuid4()),
                    email="admin@local",
                    full_name="Platform Admin",
                    is_admin=True,
                    pocket_quota=9999,
                    denovo_quota=9999,
                    admet_quota=9999,
                    docking_quota=9999,
                )
            )
            db.commit()
        return {
            "access_token": issue_session_token("admin-local", "admin"),
            "token_type": "bearer",
            "user": {"id": "admin-local", "email": "admin@local", "full_name": "Platform Admin"},
        }
    if email == settings.admin_username.lower() and settings.admin_password:
        record_login_failure(request, email)
    user = db.query(AuthUser).filter(AuthUser.email == email).first()
    valid = bool(user and verify_password(user.password, password))
    if user and not valid and settings.app_env == "development" and user.password and not user.password.startswith("$argon2id$"):
        # Migration on verified login for old local databases. Production rejects
        # legacy plaintext rows, which require an administrator-assisted reset.
        valid = hmac.compare_digest(user.password, password)
        if valid:
            user.password = hash_password(password)
            db.commit()
    if not valid:
        record_login_failure(request, email)
        raise HTTPException(status_code=401, detail="Invalid credentials")
    return {
        "access_token": issue_session_token(user.id),
        "token_type": "bearer",
        "user": {"id": user.id, "email": user.email, "full_name": user.full_name},
    }


@router.get("/auth/profile")
async def legacy_auth_profile(principal: Principal = Depends(require_session), db: Session = Depends(get_db)):
    # Never trust an email query parameter as authentication.
    if principal.role == "admin":
        return {"id": principal.user_id, "email": "admin@local", "full_name": "Platform Admin"}
    user = db.query(AuthUser).filter(AuthUser.id == principal.user_id).first()
    if not user:
        raise HTTPException(status_code=401, detail="Session user no longer exists")
    return {"id": user.id, "email": user.email, "full_name": user.full_name}


# ---------------------------------------------------------------------------
# Util / export endpoints for the Results page redesign
# ---------------------------------------------------------------------------
# Both endpoints intentionally live in routes_legacy because the Results
# modal already calls /api/jobs/* on this same router, so we keep the
# CORS / mount surface identical. Heavy lifting (RDKit MolToImage, CSV
# assembly) is small (<50 LoC) so a dedicated module would be overkill.


# Per-process LRU keyed on the rendered PNG payload. RDKit MolToImage is
# ~10-30 ms per molecule on a 2024 laptop; 5k unique SMILES * 1 KB PNG =
# ~5 MB resident — well within budget for the gain of skipping a re-render
# every time the user scrolls a 200-row table.
from functools import lru_cache  # noqa: E402


@lru_cache(maxsize=8192)
def _render_mol_png(smiles: str, width: int, height: int) -> bytes:
    """SMILES -> 2D PNG bytes. Returns a transparent 1x1 PNG on failure.

    Failure modes (un-parseable SMILES, broken templates, etc.) return the
    blank sentinel rather than raising — the Results table renders 200+
    rows in one shot and a single bad SMILES should not 500 the whole page.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("unparseable SMILES")
        # Skip expensive 3D embedding — 2D layout is enough for a thumbnail
        # and ~50x faster. Compute2DCoords gives a deterministic depiction.
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:  # noqa: BLE001
            pass
        img = Draw.MolToImage(mol, size=(int(width), int(height)))
        import io
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:  # noqa: BLE001
        # 1x1 transparent PNG — keeps the <img> tag rendering without an
        # alt-text shrug, and the upstream cache entry prevents repeated
        # failed renders of the same bad SMILES from hitting RDKit again.
        return (
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
            b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89\x00\x00\x00\rIDATx\x9cc\xf8"
            b"\x0f\x00\x00\x01\x01\x00\x01\\\xcdF\xd7\x00\x00\x00\x00IEND\xaeB`\x82"
        )


@lru_cache(maxsize=2048)
def _smiles_to_pdb_3d(smiles: str) -> str | None:
    """Generate a 3D PDB conformer from SMILES via RDKit ETKDGv3 + UFF.

    Used by the NGL viewer in the Results page to render the ligand
    stick. Returns ``None`` on parse / embed failure (the frontend then
    skips the ligand and just shows the receptor cartoon).

    NOTE: this is NOT the Vina-docked pose — it's just a low-energy
    free-state 3D conformer of the SMILES. Showing the real docked pose
    requires plumbing Docking_App's out.pdbqt back through the audit log,
    which is a separate (bigger) change. For now the viewer surfaces the
    structure-aware conformer next to the binding pocket so the chemist
    can at least eyeball steric fit / functional group positioning.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            # Fallback to random-coords seeding for stubborn ring systems
            if AllChem.EmbedMolecule(mol, useRandomCoords=True) != 0:
                return None
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:  # noqa: BLE001
            pass
        return Chem.MolToPDBBlock(mol)
    except Exception:  # noqa: BLE001
        return None


@router.get("/util/mol-pdb")
async def util_mol_pdb(smiles: str = Query(..., min_length=1, max_length=2048)):
    """SMILES -> 3D PDB text. Used by the NGL viewer.

    Returns 404 when the SMILES can't be parsed / embedded so the NGL
    viewer can fall back to receptor-only rendering instead of trying
    to load an empty model.
    """
    from fastapi.responses import PlainTextResponse

    pdb = _smiles_to_pdb_3d(smiles)
    if not pdb:
        raise HTTPException(status_code=404, detail="Could not generate 3D structure")
    return PlainTextResponse(
        content=pdb,
        media_type="chemical/x-pdb",
        headers={"Cache-Control": "public, max-age=86400"},
    )


@router.get("/util/mol-image")
async def util_mol_image(
    smiles: str = Query(..., min_length=1, max_length=2048),
    w: int = Query(default=120, ge=24, le=600),
    h: int = Query(default=80, ge=24, le=600),
):
    """Render a SMILES string as a 2D PNG thumbnail.

    Cheap RDKit-backed sketch endpoint the Results table uses to surface
    structures inline. ``w`` / ``h`` are capped at 600 px so a malicious
    client cannot waste server time rendering arbitrarily large images.
    Cached per (smiles, w, h) for the process lifetime — the common case
    is the same 50-200 candidate SMILES being requested as the table
    re-renders.
    """
    from fastapi.responses import Response

    png = _render_mol_png(smiles, w, h)
    # Encourage browser-side caching so even a hard page refresh skips the
    # network round-trip for already-rendered structures.
    return Response(
        content=png,
        media_type="image/png",
        headers={"Cache-Control": "public, max-age=86400, immutable"},
    )


def _candidates_from_output(payload: dict) -> tuple[str, list[dict]]:
    """Pull the (section, candidate rows) tuple out of an AuditLog payload.

    Centralised so the CSV exporter and any future "results detail"
    endpoint share one definition of "where the molecules live".
    """
    section = (payload.get("analysis_section") or "").strip().lower()
    results_json = payload.get("results_json") or {}
    rows = results_json.get("data") or []
    if not isinstance(rows, list):
        rows = []
    return section, rows


# ---------------------------------------------------------------------------
# REINVENT — fine-tuned prior discovery
# ---------------------------------------------------------------------------
# Status of REINVENT fine-tuning in this codebase (as of this commit):
#
#   ✓ Base prior (Reinvent_App/prior/reinvent.prior) loaded by sampler.
#   ✓ Supervised fine-tuning script exists at
#     Reinvent_App/scripts/finetune_target.py and ships an
#     EGFR_finetuned.prior under prior/finetuned/.
#   ✓ End-to-end target_prior plumbing already works in this backend:
#     production_pipeline.run() passes target_prior=target into
#     reinvent_client.sample_smiles(), which Reinvent_App resolves to
#     the fine-tuned file when a matching <TARGET>_finetuned.prior
#     exists, falling back to base otherwise.
#   ✗ No UI surface exposing "you are using a fine-tuned prior" yet.
#   ✗ No RL (reward-shaped) fine-tune loop — current script is
#     supervised next-token modeling on ChEMBL actives only.
#   ✗ No backend endpoint to kick off a new fine-tune from the UI.
#
# The endpoint below closes the smallest gap (UI discovery of priors).
# RL fine-tuning is genuinely multi-week work — it needs reward
# shaping (ADMET + Vina + Tanimoto), GPU training infra, eval harness,
# and a model registry. That's out of scope for a single endpoint
# add and should land as a dedicated "Reinvent_App/scripts/rl_finetune.py"
# wrapper around REINVENT 4's built-in RL module.


# ---------------------------------------------------------------------------
# REINVENT RL fine-tune — async job runner
# ---------------------------------------------------------------------------
# The actual RL training loop lives in Reinvent_App/scripts/rl_finetune.py.
# This backend layer adds:
#   - POST /api/reinvent/rl-finetune       launch as subprocess, returns job_id
#   - GET  /api/reinvent/rl-finetune/{id}  status + parsed progress (last log step)
#   - GET  /api/reinvent/rl-finetune/{id}/log  raw stdout tail
#   - GET  /api/reinvent/rl-finetune       list known fine-tune jobs
# Subprocess strategy (asyncio.create_subprocess_exec) was picked over an
# in-process worker because the RL loop uses PyTorch which doesn't play
# well with FastAPI's event loop on CPU — a misbehaving step would block
# every other request. Spawning a separate Python process gives us hard
# isolation and the ability to kill a runaway training without bouncing
# the API.

import asyncio  # noqa: E402

_RL_JOBS_DIR = (Path(__file__).resolve().parents[2] / "storage" / "rl_jobs")
_RL_RUNTIME: dict[str, dict] = {}  # job_id -> {"proc": asyncio.subprocess.Process, "started": iso}


def _rl_job_paths(job_id: str) -> tuple[Path, Path, Path, Path]:
    """Only canonical UUID job IDs may resolve private RL artifacts.

    Reject traversal / arbitrary paths even when called from a JSON body,
    where URL routing does not prevent embedded directory separators.
    """
    from uuid import UUID
    try:
        if not isinstance(job_id, str) or str(UUID(job_id)) != job_id:
            raise ValueError("RL job identifier must be a canonical UUID")
    except (ValueError, TypeError, AttributeError) as exc:
        raise ValueError("Invalid RL job identifier") from exc
    base = _RL_JOBS_DIR / job_id
    return (
        base / "agent.prior",      # final RL agent checkpoint
        base / "train_log.json",   # structured per-step metrics
        base / "stdout.log",       # raw stdout/stderr
        base,
    )


async def _rl_subprocess_runner(
    job_id: str,
    cmd: list[str],
    workdir: Path,
    stdout_path: Path,
) -> None:
    """Spawn the RL training subprocess and stream stdout to a log file.

    Persists a final summary to the AuditLog when the process exits so the
    Results modal can surface the trained checkpoint + log.
    """
    from app.db.session import SessionLocal

    workdir.mkdir(parents=True, exist_ok=True)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
        cwd=str(workdir.parent.parent),  # project root so reinvent_core resolves
    )
    _RL_RUNTIME[job_id] = {
        "proc": proc,
        "pid": proc.pid,
        "started": datetime.utcnow().isoformat(),
    }
    # Stream stdout to the log file line by line so the status endpoint
    # can tail the file as the run progresses. We could use `tee` but
    # this gives us cross-platform behaviour and an explicit flush point.
    try:
        with stdout_path.open("w", encoding="utf-8") as f:
            while True:
                line = await proc.stdout.readline()
                if not line:
                    break
                try:
                    f.write(line.decode("utf-8", errors="replace"))
                    f.flush()
                except Exception:  # noqa: BLE001
                    pass
        rc = await proc.wait()
    except Exception as exc:  # noqa: BLE001
        rc = -1
        logger.exception("RL subprocess failed: %s", exc)
    finally:
        _RL_RUNTIME.pop(job_id, None)

    # Update the job row + write a docking_output-style audit row so the
    # Results modal can pick it up. Action label intentionally distinct
    # ("rl_finetune_output") so the hierarchy endpoint can section it
    # separately from De Novo / ADMET / Docking jobs.
    db = SessionLocal()
    try:
        job = db.query(Job).filter(Job.id == job_id).first()
        if job is not None:
            job.status = "completed" if rc == 0 else "failed"
            job.progress_pct = 100
            job.current_step = "completed" if rc == 0 else "failed"
            job.updated_at = datetime.utcnow()
        db.add(
            AuditLog(
                id=str(uuid4()),
                actor_id="system",
                action="rl_finetune_output",
                resource_type="job_output",
                resource_id=job_id,
                payload=json.dumps({
                    "analysis_section": "rl_finetune",
                    "return_code": rc,
                    "stdout_path": str(stdout_path),
                    "checkpoint_path": str(_rl_job_paths(job_id)[0]),
                    "log_path": str(_rl_job_paths(job_id)[1]),
                }, ensure_ascii=True),
            )
        )
        db.commit()
    except Exception:  # noqa: BLE001
        db.rollback()
    finally:
        db.close()


@router.get("/reinvent/rl-finetune")
async def legacy_reinvent_rl_list(db: Session = Depends(get_db)):
    """List all RL fine-tune jobs (running + completed + cancelled).

    Resolved by joining the ``jobs`` table on the ``RL-`` receptor_name
    convention. For each job we surface the headline metrics from the
    persisted training log (or stdout fallback for jobs that ended
    before their first structured flush). Returned newest-first so the
    Lab page can render a "recent runs" table.
    """
    # Pull the candidate jobs in one query, then enrich from on-disk
    # artefacts. We deliberately read the artefact paths here (instead
    # of always loading the full JSON history) so the list stays cheap
    # for the common case of a Lab page with dozens of past runs.
    candidates = (
        db.query(Job)
        .filter(Job.receptor_name.like("RL-%"))
        .order_by(Job.created_at.desc())
        .all()
    )
    out: list[dict] = []
    for job in candidates:
        ck_path, log_path, stdout_path, _ = _rl_job_paths(job.id)
        # Parse the structured log when present — same helper as the
        # single-job status endpoint, returns {} when missing.
        progress = _parse_progress_from_log_file(log_path)
        # Derive target / steps from the audit row when present, else
        # fall back to the receptor_name convention "RL-<TARGET>-<short>".
        target = None
        if isinstance(progress, dict):
            target = progress.get("target")
        if not target and job.receptor_name and job.receptor_name.startswith("RL-"):
            parts = job.receptor_name.split("-", 2)
            if len(parts) >= 2:
                target = parts[1]
        out.append({
            "job_id": job.id,
            "target": target,
            "status": job.status,
            "created_at": job.created_at.strftime("%d.%m.%Y %H:%M") if job.created_at else None,
            "checkpoint_exists": ck_path.is_file(),
            "checkpoint_size_bytes": ck_path.stat().st_size if ck_path.is_file() else 0,
            "log_path_exists": log_path.is_file(),
            "stdout_path_exists": stdout_path.is_file(),
            "is_running": job.id in _RL_RUNTIME and _RL_RUNTIME[job.id]["proc"].returncode is None,
            "progress": {
                "current_step": progress.get("current_step"),
                "total_steps": progress.get("total_steps"),
                "best_mean_reward": progress.get("best_mean_reward"),
                "mean_reward": progress.get("mean_reward"),
                "validity": progress.get("validity"),
                "mean_qed": progress.get("mean_qed"),
                "pains_pass_rate": progress.get("pains_pass_rate"),
                "mean_tanimoto": progress.get("mean_tanimoto"),
                "elapsed_s": progress.get("elapsed_s"),
            },
        })
    return out


@router.get("/reinvent/rl-finetune/{job_id}/checkpoint")
async def legacy_reinvent_rl_download_checkpoint(job_id: str, db: Session = Depends(get_db)):
    """Stream the trained agent.prior checkpoint as a download.

    Lets the user export the model for offline analysis / sharing without
    having to navigate to the storage directory on the server.
    """
    from fastapi.responses import FileResponse

    job = db.query(Job).filter(Job.id == job_id).first()
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    ck_path, _, _, _ = _rl_job_paths(job_id)
    if not ck_path.is_file():
        raise HTTPException(status_code=404, detail="No checkpoint persisted for this job.")
    short = str(job_id)[:8]
    fname = f"rl_agent_{short}.prior"
    return FileResponse(ck_path, media_type="application/octet-stream", filename=fname)


@router.post("/reinvent/rl-finetune")
async def legacy_reinvent_rl_start(body: dict, db: Session = Depends(get_db)):
    """Kick off a REINFORCE RL fine-tune as a background subprocess.

    Body (all optional except ``target``):
        {
          "target": "EGFR",
          "steps": 200,
          "batch_size": 64,
          "sigma": 120.0,
          "lr": 0.0001,
          "weight_qed": 1.0,
          "weight_pains": 1.0,
          "weight_tani": 1.0,
          "device": "cpu" | "cuda",
          "checkpoint_every": 50
        }

    Returns ``{"job_id": "<uuid>", "status": "running", "log_url": "..."}``.

    Cost reality check: each RL step at batch=64 is ~2 min on CPU, ~5 s
    on GPU. A 200-step run is ~7 h CPU / ~17 min GPU. The endpoint does
    NOT block — it spawns the subprocess via asyncio and returns
    immediately so the caller can poll the status endpoint.
    """
    target = str(body.get("target") or "").strip().upper()
    if not target:
        raise HTTPException(status_code=400, detail="`target` is required.")
    steps = max(1, min(int(body.get("steps", 200)), 10000))
    batch_size = max(1, min(int(body.get("batch_size", 64)), 512))
    sigma = float(body.get("sigma", 120.0))
    lr = float(body.get("lr", 1e-4))
    w_qed = float(body.get("weight_qed", 1.0))
    w_pains = float(body.get("weight_pains", 1.0))
    w_tani = float(body.get("weight_tani", 1.0))
    device = str(body.get("device") or "cpu").lower()
    if device not in {"cpu", "cuda"}:
        device = "cpu"
    checkpoint_every = max(1, min(int(body.get("checkpoint_every", 50)), 1000))

    # 1) Create the Job row up front so the status endpoint has something
    # to query even before the subprocess emits its first log line.
    job_id = str(uuid4())
    receptor_value = f"RL-{target}-{job_id[:8]}"
    repo = JobRepository(db)
    repo.create(
        job_id=job_id,
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=steps,
        ),
    )
    pinned = db.query(Job).filter(Job.id == job_id).first()
    if pinned is not None:
        pinned.status = "running"
        pinned.progress_pct = 0
        pinned.current_step = "warming_up"
        pinned.progress_message = f"steps={steps} batch={batch_size} device={device}"
        pinned.updated_at = datetime.utcnow()
    db.commit()

    # 2) Resolve paths + assemble the subprocess argv.
    project_root = Path(__file__).resolve().parents[3]  # backend/app/api/.. -> project root
    script = project_root / "Reinvent_App" / "scripts" / "rl_finetune.py"
    if not script.is_file():
        raise HTTPException(
            status_code=500,
            detail=f"RL script missing at {script}",
        )
    venv_python = project_root / ".venv" / "Scripts" / "python.exe"
    # Cross-platform Python resolution: Linux uses .venv/bin/python
    if not venv_python.is_file():
        alt = project_root / ".venv" / "bin" / "python"
        if alt.is_file():
            venv_python = alt
    if not venv_python.is_file():
        venv_python = Path(sys.executable)
    out_path, log_path, stdout_path, work_base = _rl_job_paths(job_id)
    work_base.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(venv_python), str(script),
        target,
        "--steps", str(steps),
        "--batch-size", str(batch_size),
        "--sigma", str(sigma),
        "--lr", str(lr),
        "--weight-qed", str(w_qed),
        "--weight-pains", str(w_pains),
        "--weight-tani", str(w_tani),
        "--device", device,
        "--checkpoint-every", str(checkpoint_every),
        "--out-path", str(out_path),
        "--log-path", str(log_path),
    ]

    # 3) Spawn — fire-and-forget. asyncio.create_task keeps the subprocess
    # supervisor coroutine alive after this handler returns.
    asyncio.create_task(_rl_subprocess_runner(job_id, cmd, work_base, stdout_path))

    return {
        "job_id": job_id,
        "status": "running",
        "target": target,
        "steps": steps,
        "batch_size": batch_size,
        "device": device,
        "log_url": f"/api/reinvent/rl-finetune/{job_id}/log",
        "status_url": f"/api/reinvent/rl-finetune/{job_id}",
    }


def _parse_progress_from_log_file(log_path: Path) -> dict:
    """Read the per-step JSON history written by rl_finetune.py.

    Returns the most recent step's metrics + overall summary. Empty dict
    when the file doesn't exist yet (subprocess hasn't reached its first
    checkpoint flush) — caller falls back to stdout tail in that case.
    """
    if not log_path.is_file():
        return {}
    try:
        data = json.loads(log_path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {}
    hist = data.get("history") or []
    if not hist:
        return {"target": data.get("target")}
    last = hist[-1]
    return {
        "target": data.get("target"),
        "args": data.get("args"),
        "best_mean_reward": data.get("best_mean_reward"),
        "current_step": last.get("step"),
        "total_steps": (data.get("args") or {}).get("steps"),
        "loss": last.get("loss"),
        "mean_reward": last.get("mean_reward"),
        "validity": last.get("validity"),
        "mean_qed": last.get("mean_qed"),
        "pains_pass_rate": last.get("pains_pass_rate"),
        "mean_tanimoto": last.get("mean_tanimoto"),
        "elapsed_s": last.get("elapsed_s"),
        "history_len": len(hist),
    }


def _parse_progress_from_stdout(stdout_path: Path) -> dict:
    """Fallback when the structured log hasn't been written yet.

    Parses the last ``[rl] step=X/Y ...`` line for current progress.
    Cheap regex; rl_finetune.py format is stable.
    """
    if not stdout_path.is_file():
        return {}
    try:
        # Read last ~8 KB; the per-step line is short so this is plenty.
        with stdout_path.open("rb") as f:
            try:
                f.seek(-8192, os.SEEK_END)
            except OSError:
                f.seek(0)
            tail = f.read().decode("utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        return {}
    import re
    last_step = None
    last_total = None
    last_reward = None
    for m in re.finditer(
        r"\[rl\]\s+step=(\d+)/(\d+)\s+loss=([\d.eE+-]+)\s+reward=([\d.eE+-]+)",
        tail,
    ):
        last_step = int(m.group(1))
        last_total = int(m.group(2))
        last_reward = float(m.group(4))
    if last_step is None:
        return {}
    return {
        "current_step": last_step,
        "total_steps": last_total,
        "mean_reward": last_reward,
    }


@router.get("/reinvent/rl-finetune/{job_id}")
async def legacy_reinvent_rl_status(job_id: str, db: Session = Depends(get_db)):
    """Status + progress for an RL fine-tune job.

    Reads the structured log file first (preferred — has every metric)
    and falls back to parsing the stdout tail when the log hasn't been
    written yet (first ~50 steps before the first checkpoint flush).
    """
    job = db.query(Job).filter(Job.id == job_id).first()
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    out_path, log_path, stdout_path, _ = _rl_job_paths(job_id)
    progress = _parse_progress_from_log_file(log_path) or _parse_progress_from_stdout(stdout_path)
    is_running = job_id in _RL_RUNTIME and _RL_RUNTIME[job_id]["proc"].returncode is None
    return {
        "job_id": job_id,
        "status": "running" if is_running else (job.status or "unknown"),
        "checkpoint_exists": out_path.is_file(),
        "checkpoint_size_bytes": out_path.stat().st_size if out_path.is_file() else 0,
        "log_url": f"/api/reinvent/rl-finetune/{job_id}/log",
        "progress": progress,
    }


@router.get("/reinvent/rl-finetune/{job_id}/log")
async def legacy_reinvent_rl_log(job_id: str, tail: int = 200):
    """Tail of the RL subprocess stdout. Caller can poll for live updates."""
    _, _, stdout_path, _ = _rl_job_paths(job_id)
    if not stdout_path.is_file():
        return {"job_id": job_id, "lines": []}
    try:
        text = stdout_path.read_text(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        return {"job_id": job_id, "lines": []}
    lines = text.splitlines()
    tail = max(1, min(int(tail), 5000))
    return {"job_id": job_id, "lines": lines[-tail:], "total_lines": len(lines)}


@router.delete("/reinvent/rl-finetune/{job_id}")
async def legacy_reinvent_rl_cancel(job_id: str, db: Session = Depends(get_db)):
    """Terminate a running RL fine-tune job. The checkpoint produced so
    far (if any) is kept on disk — only the subprocess is killed."""
    rec = _RL_RUNTIME.get(job_id)
    if rec is None:
        raise HTTPException(status_code=404, detail="No running RL job with that id.")
    proc = rec["proc"]
    try:
        proc.terminate()
    except Exception:  # noqa: BLE001
        pass
    # Best-effort SIGKILL if termiante didn't take effect within a few sec.
    try:
        await asyncio.wait_for(proc.wait(), timeout=5.0)
    except asyncio.TimeoutError:
        try:
            proc.kill()
        except Exception:  # noqa: BLE001
            pass
    job = db.query(Job).filter(Job.id == job_id).first()
    if job is not None:
        job.status = "cancelled"
        job.updated_at = datetime.utcnow()
        db.commit()
    return {"ok": True, "job_id": job_id}


def _pains_score_smiles(smiles_list: list[str]) -> tuple[int, int]:
    """Returns (n_pains_pass, n_evaluated). Skips invalid SMILES.

    Same PAINS_A+B+C catalog as the RL reward — keeps the eval consistent
    with what the RL training optimised for.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import FilterCatalog
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:  # noqa: BLE001
        return 0, 0

    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    catalog = FilterCatalog.FilterCatalog(params)
    n_pass = 0
    n_eval = 0
    for smi in smiles_list:
        if not isinstance(smi, str) or not smi.strip():
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        n_eval += 1
        try:
            if not catalog.HasMatch(mol):
                n_pass += 1
        except Exception:  # noqa: BLE001
            n_pass += 1  # benefit of the doubt on catalog parse errors
    return n_pass, n_eval


def _aggregate_eval_metrics(candidates: list[dict]) -> dict:
    """Aggregate per-molecule fields returned by Reinvent_App /sample into
    population-level metrics suitable for prior A/B comparison."""
    n_total = len(candidates)
    valid = [c for c in candidates if c.get("valid") is not False]
    n_valid = len(valid)
    # Reinvent_App returns mw/logp/qed only when validate_with_rdkit=True
    # and only for valid molecules; mean over valid only.
    def _mean(field: str):
        xs = [float(c[field]) for c in valid if isinstance(c.get(field), (int, float))]
        return round(sum(xs) / len(xs), 4) if xs else None
    smiles = [c.get("smiles") for c in valid if isinstance(c.get("smiles"), str)]
    unique = len({s for s in smiles if s})
    pains_pass, pains_eval = _pains_score_smiles(smiles)
    # Max Tanimoto is per-molecule from Reinvent_App; compute mean for
    # population-level "how close to actives" signal.
    tani_xs = [float(c["max_tanimoto"]) for c in valid if isinstance(c.get("max_tanimoto"), (int, float))]
    return {
        "n_requested": n_total,
        "n_valid": n_valid,
        "n_unique_valid": unique,
        "validity_rate": round(n_valid / n_total, 4) if n_total else None,
        "uniqueness_rate": round(unique / n_valid, 4) if n_valid else None,
        "mean_qed": _mean("qed"),
        "mean_mw": _mean("mw"),
        "mean_logp": _mean("logp"),
        "mean_tpsa": _mean("tpsa"),
        "pains_pass_rate": round(pains_pass / pains_eval, 4) if pains_eval else None,
        "mean_max_tanimoto": round(sum(tani_xs) / len(tani_xs), 4) if tani_xs else None,
    }


@router.post("/reinvent/evaluate-prior")
async def legacy_reinvent_evaluate_prior(body: dict):
    """Sample N SMILES from a prior and report population-level metrics.

    Body:
        {
          "prior": "EGFR" | "base",
          "n_samples": 200,
          "reference_smiles": ["..."] (optional, enables Tanimoto),
          "compare_to_base": true (optional, samples base in the same call)
        }

    Returns metrics for the chosen prior (and base, when comparing), plus
    per-metric delta. Use this to verify a fine-tune actually moved the
    distribution toward the reward stack instead of just memorising the
    training set.

    NOTE: RL checkpoints written under storage/rl_jobs/<id>/agent.prior are
    not discoverable by Reinvent_App until promoted via
    /api/reinvent/promote-checkpoint. This endpoint only evaluates priors
    that Reinvent_App's sampler can resolve by name.
    """
    from app.services.reinvent_client import sample_smiles

    prior_in = str(body.get("prior") or "base").strip()
    n_samples = max(20, min(int(body.get("n_samples") or 200), 2000))
    refs = body.get("reference_smiles") or []
    if not isinstance(refs, list):
        refs = []
    compare = bool(body.get("compare_to_base"))

    target_arg = None if prior_in.lower() == "base" else prior_in
    candidates, prior_used = await sample_smiles(
        batch_size=n_samples,
        validate_with_rdkit=True,
        reference_smiles=refs if refs else None,
        target_prior=target_arg,
    )
    metrics = _aggregate_eval_metrics(candidates)
    out: dict = {
        "prior": prior_in,
        "prior_used": prior_used,
        "metrics": metrics,
        "sample_smiles": [c.get("smiles") for c in candidates[:10] if c.get("smiles")],
    }

    if compare and prior_used.lower() != "base":
        base_candidates, _ = await sample_smiles(
            batch_size=n_samples,
            validate_with_rdkit=True,
            reference_smiles=refs if refs else None,
            target_prior=None,
        )
        base_metrics = _aggregate_eval_metrics(base_candidates)
        # Per-metric delta (candidate − base). Positive delta means the
        # fine-tune moved the distribution in the direction that metric
        # measures. For QED / validity / Tanimoto bigger=better → positive
        # delta is good. For MW / LogP delta sign is preference-dependent.
        delta: dict = {}
        for k in ("validity_rate", "uniqueness_rate", "mean_qed", "mean_mw",
                  "mean_logp", "mean_tpsa", "pains_pass_rate", "mean_max_tanimoto"):
            a = metrics.get(k)
            b = base_metrics.get(k)
            if isinstance(a, (int, float)) and isinstance(b, (int, float)):
                delta[k] = round(a - b, 4)
        out["base_metrics"] = base_metrics
        out["delta"] = delta
        out["base_sample_smiles"] = [c.get("smiles") for c in base_candidates[:10] if c.get("smiles")]
    return out


@router.post("/reinvent/promote-checkpoint")
async def legacy_reinvent_promote_checkpoint(body: dict, db: Session = Depends(get_db)):
    """Copy an RL-trained checkpoint into the Reinvent_App finetuned slot
    so the sampler can resolve it by name.

    Body:
        {"job_id": "<rl_job_id>", "target_label": "EGFR_rl"}

    The target_label is what /api/reinvent/priors will surface and what
    callers pass as target_prior to use this model. We deliberately
    require a distinct label (e.g. "EGFR_rl") so promoting an RL model
    does NOT overwrite a previously-promoted supervised TL prior named
    plain "EGFR" — the chemist can A/B them and pick.
    """
    import shutil

    job_id = str(body.get("job_id") or "").strip()
    label = str(body.get("target_label") or "").strip()
    if not job_id or not label:
        raise HTTPException(status_code=400, detail="job_id + target_label required.")
    # Sanitise the label — only alnum + underscore + hyphen. This becomes
    # the filename prefix the sampler globs against.
    import re
    safe_label = re.sub(r"[^A-Za-z0-9_-]+", "_", label).strip("_")
    if not safe_label:
        raise HTTPException(status_code=400, detail="target_label is empty after sanitisation.")
    try:
        src, _, _, _ = _rl_job_paths(job_id)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid RL job identifier") from None
    job = db.query(Job).filter(Job.id == job_id).first()
    if job is None or not (job.receptor_name or "").startswith("RL-"):
        raise HTTPException(status_code=404, detail="RL job not found")
    if job.status != "completed":
        raise HTTPException(status_code=409, detail="RL training must complete before promotion")
    if not src.is_file() or src.is_symlink():
        raise HTTPException(status_code=404, detail="Checkpoint not found or not a regular trusted file")
    project_root = Path(__file__).resolve().parents[3]
    dst = project_root / "Reinvent_App" / "prior" / "finetuned" / f"{safe_label}_finetuned.prior"
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise HTTPException(
            status_code=409,
            detail=f"Destination {dst.name} already exists. Pick a different target_label.",
        )
    try:
        shutil.copy2(src, dst)
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=500, detail=f"Copy failed: {exc}") from exc
    # Audit the promotion so the operator can track it.
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="rl_checkpoint_promoted",
            resource_type="job_output",
            resource_id=job_id,
            payload=json.dumps({
                "analysis_section": "rl_finetune",
                "src": str(src),
                "dst": str(dst),
                "target_label": safe_label,
            }, ensure_ascii=True),
        )
    )
    db.commit()
    return {
        "ok": True,
        "job_id": job_id,
        "src": str(src),
        "dst": str(dst),
        "available_as": safe_label,
    }


# ---------------------------------------------------------------------------
# Extended receptor catalog — on-demand RCSB → meeko → Vina docking.
# ---------------------------------------------------------------------------
# The 8 entries under Docking_App/receptors/ are *bundled* (PDBQT + Vina
# conf.txt grid box shipped). The catalog under app.data.receptor_catalog
# adds ~75 more receptors the chemist can dock against on-demand:
#
#   1) User picks a receptor in the UI from the extended catalog.
#   2) Frontend POSTs /api/receptors/extended/prepare?name=EGFR.
#      Backend fetches the canonical PDB from RCSB and caches it under
#      storage/raw_pdbs/<PDB>.pdb. Idempotent — already-cached PDBs
#      short-circuit without a network round-trip.
#   3) The chemist runs Step 1 (Pocket Analysis) on the cached PDB to
#      get pocket centre + size — we do NOT pre-curate pockets for the
#      extended catalog because each project has its own preferred site.
#   4) Docking proceeds via the existing custom_receptor_pdb path
#      (production_pipeline reads the cached PDB and passes its text
#      content to Docking_App, which runs meeko per-call). No new
#      pipeline plumbing required.
#
# Why no pre-prep into PDBQT? Meeko produces deterministic output for
# the same PDB content, and the Docking_App already caches by content
# hash (Docking_App/receptors/custom_cache/<PDB>_<hash>.pdbqt). Doing
# the prep here would just duplicate that cache without saving wall
# time on the first dock.


@router.get("/receptors/extended-catalog")
async def legacy_receptors_extended_catalog():
    """List every receptor in the on-demand catalog + family histogram.

    Used by the frontend dropdown to render the extended-target picker
    grouped by family. Bundled receptors (Docking_App/receptors/) are
    served by ``/api/docking/targets`` — these two lists are
    complementary, not overlapping.
    """
    from app.data.receptor_catalog import all_entries, families_with_counts

    return {
        "entries": all_entries(),
        "families": families_with_counts(),
        "count": len(all_entries()),
    }


@router.get("/receptors/extended-catalog/{name}")
async def legacy_receptors_extended_lookup(name: str):
    """Resolve a free-form name / gene / alias to a catalog entry.

    Case-insensitive. Returns 404 with the original query echoed back so
    the UI can render "no extended target named 'XYZ'" without guessing
    which lookup failed.
    """
    from app.data.receptor_catalog import lookup, to_dict

    entry = lookup(name)
    if entry is None:
        raise HTTPException(
            status_code=404,
            detail=f"No extended catalog entry for {name!r}. Try the gene symbol or a common alias.",
        )
    return to_dict(entry)


@router.post("/receptors/extended/prepare")
async def legacy_receptors_extended_prepare(body: dict):
    """Fetch + cache the canonical PDB for an extended-catalog receptor.

    Body shape::

        {"name": "HER2"}    # resolves through catalog.lookup
        {"pdb":  "3PP0"}    # bypass catalog, just download this PDB

    Returns ``{"pdb": "...", "cached_path": "...", "size_bytes": ...}``.
    Idempotent — if the file already exists under storage/raw_pdbs/
    we skip the network call entirely and just report size.

    On RCSB failure we surface a 502 with the response body (truncated)
    so the UI can show the underlying error instead of a generic 500.
    Common failure modes: invalid PDB id (404 from RCSB), obsolete
    entry (returns the replacement id in headers, which we don't
    follow).
    """
    from app.data.receptor_catalog import lookup

    raw_pdb_id = (body.get("pdb") or "").strip().upper()
    name = (body.get("name") or "").strip()

    pdb_id: str | None = None
    catalog_match = None
    if raw_pdb_id:
        # Direct PDB override — just validate the 4-char shape.
        if not _is_valid_pdb_id(raw_pdb_id):
            raise HTTPException(status_code=400, detail=f"Invalid PDB id {raw_pdb_id!r}")
        pdb_id = raw_pdb_id
    elif name:
        catalog_match = lookup(name)
        if catalog_match is None:
            raise HTTPException(
                status_code=404,
                detail=f"No extended catalog entry for {name!r}",
            )
        pdb_id = catalog_match.pdb
    else:
        raise HTTPException(status_code=400, detail="Provide 'name' or 'pdb' in the body.")

    cache_dir = _raw_pdb_cache_dir()
    project_root = Path(__file__).resolve().parents[2]
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{pdb_id}.pdb"

    if cached.is_file() and cached.stat().st_size > 0:
        # Already cached — return existing file metadata without hitting
        # RCSB. The chemist can call this endpoint as a cheap pre-flight
        # before kicking off docking.
        return {
            "pdb": pdb_id,
            "cached": True,
            "size_bytes": cached.stat().st_size,
            "path": _safe_relative(cached, project_root),
            "catalog_entry": _catalog_entry_dict(catalog_match),
        }

    # Network fetch from RCSB's static files mirror. We use the legacy
    # /download/{id}.pdb endpoint because it returns vanilla PDB text;
    # the newer .cif endpoint would force a format conversion downstream.
    rcsb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(rcsb_url)
    except httpx.HTTPError as exc:
        raise HTTPException(
            status_code=502,
            detail=f"RCSB fetch failed for {pdb_id}: {exc!s}",
        ) from exc
    if r.status_code != 200:
        raise HTTPException(
            status_code=502,
            detail=f"RCSB returned {r.status_code} for {pdb_id}: {r.text[:300]}",
        )
    pdb_text = r.text
    if not pdb_text.startswith(("HEADER", "ATOM", "REMARK")):
        raise HTTPException(
            status_code=502,
            detail=f"RCSB returned non-PDB content for {pdb_id} (first 80 chars: {pdb_text[:80]!r})",
        )
    # write_bytes (not write_text) — Path.write_text on Windows would do
    # newline translation (\n → \r\n), bloating the file by ~1 byte per
    # ATOM record and breaking the byte-exact contract Docking_App's
    # meeko prep expects.
    cached.write_bytes(pdb_text.encode("utf-8"))
    return {
        "pdb": pdb_id,
        "cached": False,
        "size_bytes": cached.stat().st_size,
        "path": _safe_relative(cached, project_root),
        "catalog_entry": _catalog_entry_dict(catalog_match),
    }


def _is_valid_pdb_id(raw: str) -> bool:
    """RCSB PDB id format: first char 1-9, three more alphanumeric."""
    import re
    return bool(re.fullmatch(r"[1-9][A-Za-z0-9]{3}", raw or ""))


def _raw_pdb_cache_dir() -> Path:
    """Resolve the on-demand PDB cache dir.

    ``RAW_PDB_CACHE_DIR`` env var wins so tests can redirect to tmp_path
    without monkeypatching the route handler. Default lives under
    ``backend/storage/raw_pdbs/`` (gitignored — see .gitignore).
    """
    env_override = os.environ.get("RAW_PDB_CACHE_DIR")
    if env_override:
        return Path(env_override)
    project_root = Path(__file__).resolve().parents[2]
    return project_root / "storage" / "raw_pdbs"


def _safe_relative(path: Path, base: Path) -> str:
    """``path.relative_to(base)`` that falls back to the absolute path
    when ``path`` isn't under ``base``. Used for the cached PDB
    response so a test-overridden cache dir doesn't crash with
    ValueError when the path is in /tmp instead of the project root.
    """
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def _catalog_entry_dict(entry) -> dict | None:
    """Inline catalog entry → dict to keep the prepare response self-contained.

    Returns None when the user passed a direct ``pdb`` override (no
    catalog match). The frontend can use this to decide whether to
    show "manual receptor" vs "{name} (gene)".
    """
    if entry is None:
        return None
    from app.data.receptor_catalog import to_dict
    return to_dict(entry)


@router.get("/reinvent/priors")
async def legacy_reinvent_priors():
    """List base + fine-tuned priors discoverable by Reinvent_App.

    The De Novo pipeline auto-selects a fine-tuned prior when its name
    matches the requested target. This endpoint lets the UI surface
    "EGFR — fine-tuned ✓" vs "MET — base prior" badges so the chemist
    knows which model is generating their candidates.
    """
    from app.services.reinvent_client import list_available_priors

    try:
        priors = await list_available_priors()
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(
            status_code=503,
            detail=f"Reinvent_App unreachable: {exc}",
        ) from exc
    # Normalise the response: base is always implicitly available, so
    # callers always see at least one entry. Fine-tuned priors are
    # surfaced with their target label exactly as Reinvent_App reports.
    fine_tuned = [p for p in (priors or []) if p.lower() not in {"base", "default", "generic"}]
    return {
        "base_available": True,
        "fine_tuned": sorted(fine_tuned),
        "all": sorted(set(priors or []) | {"base"}),
    }


# ---------------------------------------------------------------------------
# Phase 3 — Selectivity / off-target cross-dock panel
# ---------------------------------------------------------------------------
# Idea: take a small set of candidate SMILES, dock them against EVERY
# bundled receptor (currently 8: AKT1, BACE1, EGFR, JAK2, KDR, KIT, LCK,
# MET) and emit a kcal/mol matrix the UI can render as a heatmap. The
# difference between "intended-target affinity" and "best off-target
# affinity" is a first-pass selectivity signal — a candidate that binds
# EGFR at -10 kcal/mol but also JAK2 at -9.5 is probably going to be
# promiscuous in vivo, while one that binds EGFR at -10 and everything
# else at -6 is a much cleaner lead.
#
# Cost: ~10-15 s per (molecule, target) on the bundled Vina path. Capping
# top_k at 5 (default 3) keeps a single panel call under 10 minutes for
# the 8-target bundle. For batches above this the caller should switch
# to an async job pattern, which is a separate (future) endpoint.


# In-process registry of running selectivity jobs so the status endpoint
# can read partial progress between Docking_App calls. Mirrors the RL
# subprocess registry pattern. {job_id -> {"current_target": str,
# "targets": [str], "matrix": dict, "started": iso, "rows_partial": list}}
_SEL_RUNTIME: dict[str, dict] = {}


async def _run_cross_dock(
    job_id: str,
    smiles_top: list[str],
    targets: list[str],
    primary_target: Optional[str],
    source_job_id: Optional[str],
    top_k: int,
) -> None:
    """Background worker: dock smiles_top against every target sequentially.

    Updates _SEL_RUNTIME[job_id] after each target so the status endpoint
    can return partial progress (n_targets_done / total) without waiting
    for the full 5-10 min wall-clock. Persists the final result as a
    selectivity_output AuditLog row so the Results modal can re-fetch.
    """
    from app.db.session import SessionLocal
    from app.services.docking_client import dock_batch

    state = _SEL_RUNTIME[job_id]
    matrix: dict[str, dict[str, Optional[float]]] = {smi: {} for smi in smiles_top}
    failed_targets: list[str] = []
    for target in targets:
        state["current_target"] = target
        state["targets_done"] = len([t for t in targets if t in state.get("completed_targets", [])])
        try:
            results = await dock_batch(smiles_top, target)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Cross-dock failed for target %s: %s", target, exc)
            failed_targets.append(target)
            for smi in smiles_top:
                matrix[smi][target] = None
            state.setdefault("completed_targets", []).append(target)
            continue
        by_smi: dict[str, Optional[float]] = {smi: None for smi in smiles_top}
        for r in results or []:
            if not isinstance(r, dict):
                continue
            smi = r.get("smiles")
            aff = r.get("affinity")
            if isinstance(smi, str) and isinstance(aff, (int, float)):
                by_smi[smi] = float(aff)
        for smi in smiles_top:
            matrix[smi][target] = by_smi.get(smi)
        state.setdefault("completed_targets", []).append(target)
        # Update progress matrix so the status endpoint can stream partial
        # heatmap data instead of waiting until the final flush.
        state["matrix_partial"] = {smi: dict(scores) for smi, scores in matrix.items()}

    # Per-row summary metrics — same shape as the sync endpoint.
    rows = []
    for smi in smiles_top:
        scores = matrix.get(smi, {})
        finite = [v for v in scores.values() if isinstance(v, (int, float))]
        primary_aff = scores.get(primary_target) if primary_target else None
        off_target_finite = [
            v for k, v in scores.items()
            if isinstance(v, (int, float)) and k != primary_target
        ]
        best_off = min(off_target_finite) if off_target_finite else None
        delta = None
        if isinstance(primary_aff, (int, float)) and isinstance(best_off, (int, float)):
            delta = round(best_off - primary_aff, 3)
        rows.append({
            "smiles": smi,
            "scores": scores,
            "primary_affinity": primary_aff,
            "best_off_target": best_off,
            "selectivity_delta": delta,
            "tested_targets": len(finite),
        })

    # Persist to the AuditLog + mark the Job row completed.
    db = SessionLocal()
    try:
        job = db.query(Job).filter(Job.id == job_id).first()
        if job is not None:
            job.status = "completed"
            job.progress_pct = 100
            job.current_step = "completed"
            job.updated_at = datetime.utcnow()
        audit_payload = {
            "analysis_section": "selectivity",
            "source_job_id": source_job_id,
            "primary_target": primary_target,
            "targets": targets,
            "failed_targets": failed_targets,
            "top_k": top_k,
            "results_json": {
                "section": "selectivity",
                "primary_target": primary_target,
                "targets": targets,
                "failed_targets": failed_targets,
                "matrix": matrix,
                "data": rows,
            },
        }
        db.add(
            AuditLog(
                id=str(uuid4()),
                actor_id="system",
                action="selectivity_output",
                resource_type="job_output",
                resource_id=job_id,
                payload=json.dumps(audit_payload, ensure_ascii=True),
            )
        )
        db.commit()
    except Exception:  # noqa: BLE001
        db.rollback()
    finally:
        db.close()
        # Keep state.matrix_partial around for late polls
        state["final"] = True


@router.post("/selectivity/cross-dock-async")
async def legacy_selectivity_cross_dock_async(body: dict, db: Session = Depends(get_db)):
    """Async variant of /selectivity/cross-dock — returns immediately with a
    job_id; caller polls /selectivity/cross-dock/{job_id} for progress.

    Body shape is identical to the sync endpoint. Use this when wall-clock
    matters — the sync version blocks the caller for the full 5-10 min run
    which is bad UX for users hitting it from the browser.
    """
    from app.services.docking_client import list_available_targets

    source_job_id = str(body.get("job_id") or "").strip()
    raw_smiles = body.get("smiles") or []
    primary_target_in = str(body.get("primary_target") or "").strip()
    top_k = max(1, min(int(body.get("top_k") or 3), 5))

    smiles_list: list[str] = []
    inferred_target: Optional[str] = None
    if source_job_id:
        smiles_list, inferred_target = _job_output_smiles_and_target(db, source_job_id)
    elif isinstance(raw_smiles, list):
        smiles_list = [s.strip() for s in raw_smiles if isinstance(s, str) and s.strip()]
    if not smiles_list:
        raise HTTPException(status_code=400, detail="No SMILES to cross-dock.")
    primary_target = (primary_target_in or inferred_target or "").upper() or None
    smiles_top = smiles_list[:top_k]

    targets = await list_available_targets()
    if not targets:
        raise HTTPException(status_code=503, detail="Docking_App returned no targets.")

    # Job row up front so polling has something to read.
    job_id = str(uuid4())
    receptor_value = f"SEL-{primary_target or 'MIXED'}-{(source_job_id[:8] if source_job_id else 'EXT')}"
    repo = JobRepository(db)
    repo.create(
        job_id=job_id,
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=len(smiles_top),
        ),
    )
    pinned = db.query(Job).filter(Job.id == job_id).first()
    if pinned is not None:
        pinned.status = "running"
        pinned.progress_pct = 0
        pinned.current_step = "warming_up"
        pinned.progress_message = f"top_k={top_k} targets={len(targets)}"
        pinned.updated_at = datetime.utcnow()
    db.commit()

    _SEL_RUNTIME[job_id] = {
        "smiles_top": list(smiles_top),
        "targets": list(targets),
        "primary_target": primary_target,
        "started": datetime.utcnow().isoformat(),
        "completed_targets": [],
        "matrix_partial": {smi: {} for smi in smiles_top},
        "final": False,
    }
    asyncio.create_task(_run_cross_dock(job_id, smiles_top, targets, primary_target, source_job_id or None, top_k))
    return {
        "job_id": job_id,
        "status": "running",
        "primary_target": primary_target,
        "targets": targets,
        "top_k": top_k,
        "smiles_count": len(smiles_top),
        "status_url": f"/api/selectivity/cross-dock/{job_id}",
    }


@router.get("/selectivity/cross-dock/{job_id}")
async def legacy_selectivity_cross_dock_status(job_id: str, db: Session = Depends(get_db)):
    """Poll status of an async cross-dock job.

    Returns the partial matrix while running so the UI can render rows
    progressively (greying out columns that haven't completed yet) instead
    of waiting for the full 5-10 min run.
    """
    job = db.query(Job).filter(Job.id == job_id).first()
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    state = _SEL_RUNTIME.get(job_id, {})
    completed = state.get("completed_targets", [])
    targets = state.get("targets", [])
    # Done = either runtime says final OR the job row is no longer "running".
    final = bool(state.get("final")) or (job.status not in ("running", "queued"))
    out = {
        "job_id": job_id,
        "status": "running" if not final else (job.status or "completed"),
        "primary_target": state.get("primary_target"),
        "targets": targets,
        "completed_targets": completed,
        "current_target": state.get("current_target"),
        "progress_pct": int(100 * len(completed) / max(1, len(targets))),
        "matrix_partial": state.get("matrix_partial", {}),
    }
    # When final and the audit row was persisted, surface the full row
    # summary so the UI can render ΔSel without re-deriving from matrix.
    if final:
        row = (
            db.query(AuditLog)
            .filter(
                AuditLog.resource_type == "job_output",
                AuditLog.resource_id == str(job_id),
                AuditLog.action == "selectivity_output",
            )
            .order_by(desc(AuditLog.created_at))
            .first()
        )
        if row:
            try:
                payload = json.loads(row.payload or "{}")
                res = payload.get("results_json") or {}
                out["data"] = res.get("data") or []
                out["matrix"] = res.get("matrix") or {}
                out["failed_targets"] = res.get("failed_targets") or []
            except Exception:  # noqa: BLE001
                pass
    return out


@router.post("/selectivity/cross-dock")
async def legacy_selectivity_cross_dock(body: dict, db: Session = Depends(get_db)):
    """Dock a small set of candidate SMILES against every bundled receptor.

    Body (one of):
        {"job_id": "<denovo or docking job>", "top_k": 3}    # mode A: pull SMILES from job
        {"smiles": ["CCO", ...], "top_k": 3}                  # mode B: external list
        {"smiles": [...], "primary_target": "EGFR"}           # primary highlight on the heatmap

    Returns ``{"job_id": <new selectivity job>, "matrix": {smi: {target: kcal_mol|null}}}``.

    The matrix is also persisted as a ``selectivity_output`` AuditLog so
    the Results modal can re-fetch it after a refresh without redoing the
    8x docking (which can take 5-10 minutes).
    """
    from app.services.docking_client import dock_batch, list_available_targets

    # Resolve the SMILES list — either pull from an existing De Novo /
    # ADMET / Docking job's audit log or take an external list verbatim.
    source_job_id = str(body.get("job_id") or "").strip()
    raw_smiles = body.get("smiles") or []
    primary_target_in = str(body.get("primary_target") or "").strip()
    top_k = max(1, min(int(body.get("top_k") or 3), 5))

    smiles_list: list[str] = []
    inferred_target: Optional[str] = None
    if source_job_id:
        smiles_list, inferred_target = _job_output_smiles_and_target(db, source_job_id)
    elif isinstance(raw_smiles, list):
        smiles_list = [s.strip() for s in raw_smiles if isinstance(s, str) and s.strip()]
    if not smiles_list:
        raise HTTPException(
            status_code=400,
            detail="No SMILES to cross-dock (need job_id or non-empty smiles list).",
        )
    primary_target = (primary_target_in or inferred_target or "").upper() or None
    smiles_top = smiles_list[:top_k]

    # Discover the bundled receptor list at call time so the panel
    # automatically widens when new receptors are added to Docking_App.
    targets = await list_available_targets()
    if not targets:
        raise HTTPException(
            status_code=503,
            detail="Docking_App returned no bundled targets — is the service running?",
        )

    # Run docking one target at a time so we get progressive results in
    # the audit log even if the user cancels halfway. We use the Docking_App's
    # batched /dock so the per-target call already parallelises ligands.
    matrix: dict[str, dict[str, Optional[float]]] = {smi: {} for smi in smiles_top}
    failed_targets: list[str] = []
    for target in targets:
        try:
            results = await dock_batch(smiles_top, target)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Cross-dock failed for target %s: %s", target, exc)
            failed_targets.append(target)
            for smi in smiles_top:
                matrix[smi][target] = None
            continue
        # Build a per-target {smiles -> affinity} map, then stamp into matrix.
        by_smi: dict[str, Optional[float]] = {smi: None for smi in smiles_top}
        for r in results or []:
            if not isinstance(r, dict):
                continue
            smi = r.get("smiles")
            aff = r.get("affinity")
            if isinstance(smi, str) and isinstance(aff, (int, float)):
                by_smi[smi] = float(aff)
        for smi in smiles_top:
            matrix[smi][target] = by_smi.get(smi)

    # Per-row summary fields so the UI can rank rows by selectivity
    # without re-computing on every render.
    rows = []
    for smi in smiles_top:
        scores = matrix.get(smi, {})
        finite = [v for v in scores.values() if isinstance(v, (int, float))]
        primary_aff = scores.get(primary_target) if primary_target else None
        # Best off-target affinity = strongest binding to ANY target other
        # than the primary. "Strongest" = most negative kcal/mol. None when
        # we have no off-target hits at all (e.g. all-failed).
        off_target_finite = [
            v for k, v in scores.items()
            if isinstance(v, (int, float)) and k != primary_target
        ]
        best_off = min(off_target_finite) if off_target_finite else None
        # Selectivity delta — positive means "binds primary stronger than
        # every off-target". Bigger delta = cleaner ligand.
        delta = None
        if isinstance(primary_aff, (int, float)) and isinstance(best_off, (int, float)):
            delta = round(best_off - primary_aff, 3)
        rows.append({
            "smiles": smi,
            "scores": scores,
            "primary_affinity": primary_aff,
            "best_off_target": best_off,
            "selectivity_delta": delta,
            "tested_targets": len(finite),
        })

    # Persist as a new Job + AuditLog so the Results modal can read it.
    source_label = source_job_id[:8] if source_job_id else "EXTERNAL"
    receptor_value = f"SEL-{primary_target or 'MIXED'}-{source_label}"
    repo = JobRepository(db)
    job = repo.create(
        job_id=str(uuid4()),
        payload=JobCreateRequest(
            tenant_id="default",
            project_id="legacy-ui",
            pdb_id=receptor_value,
            receptor_name=receptor_value,
            mode="agonist",
            candidate_count=max(1, len(smiles_top)),
        ),
    )
    pinned = db.query(Job).filter(Job.id == job.id).first()
    if pinned is not None:
        pinned.status = "completed"
        pinned.updated_at = datetime.utcnow()
    audit_payload = {
        "analysis_section": "selectivity",
        "source_job_id": source_job_id or None,
        "primary_target": primary_target,
        "targets": targets,
        "failed_targets": failed_targets,
        "top_k": top_k,
        "results_json": {
            "section": "selectivity",
            "primary_target": primary_target,
            "targets": targets,
            "failed_targets": failed_targets,
            "matrix": matrix,
            "data": rows,
        },
    }
    db.add(
        AuditLog(
            id=str(uuid4()),
            actor_id="system",
            action="selectivity_output",
            resource_type="job_output",
            resource_id=job.id,
            payload=json.dumps(audit_payload, ensure_ascii=True),
        )
    )
    db.commit()
    return {
        "job_id": job.id,
        "status": "completed",
        "primary_target": primary_target,
        "targets": targets,
        "failed_targets": failed_targets,
        "matrix": matrix,
        "data": rows,
    }


@router.get("/jobs/{job_id}/poses/{rank}.pdb")
async def legacy_job_pose_pdb(
    job_id: str,
    rank: int,
    mode: str = Query(default="first", description="'first' = top pose, 'all' = every Vina MODEL"),
    db: Session = Depends(get_db),
):
    """Serve the persisted Vina docking pose (PDB block) for one candidate.

    Phase B captured the top pose (``pose_pdb`` field on each candidate).
    Phase B+1 adds ``pose_pdb_all`` — a multi-MODEL block with every Vina
    pose for that candidate. ``?mode=all`` serves the multi-MODEL block
    when available so the NGL viewer can cycle through alternate binding
    modes; missing fall-back is the single top pose.

    Returns 404 when:
      - the job has no docking_output / denovo_output entry
      - the rank exceeds the persisted row count
      - the row has no pose (e.g. legacy data written before Phase B, or a
        Vina failure that recorded affinity=None / no pose).

    Resolution order on the audit log: docking_output first (the iterative
    Step 5 writer), then denovo_output (the De Novo pipeline persists poses
    when target_aware docking ran). admet_output never has poses.
    """
    from fastapi.responses import PlainTextResponse

    row = (
        db.query(AuditLog)
        .filter(
            AuditLog.resource_type == "job_output",
            AuditLog.resource_id == str(job_id),
            AuditLog.action.in_(["docking_output", "denovo_output"]),
        )
        .order_by(desc(AuditLog.created_at))
        .first()
    )
    if not row:
        raise HTTPException(status_code=404, detail="No docking output for this job.")
    try:
        payload = json.loads(row.payload or "{}")
    except Exception:  # noqa: BLE001
        payload = {}
    results_json = payload.get("results_json") or {}
    candidates = results_json.get("data") or []
    if not isinstance(candidates, list):
        raise HTTPException(status_code=404, detail="Output payload has no candidate list.")

    # Caller indexes rank 1..N; the persisted row may use 0-based or
    # already-1-based rank field. Prefer the row's own "rank" field when
    # present (single source of truth) and fall back to list position.
    target_row = None
    for idx, item in enumerate(candidates, start=1):
        if not isinstance(item, dict):
            continue
        row_rank = item.get("rank")
        if isinstance(row_rank, int) and row_rank == rank:
            target_row = item
            break
        if row_rank is None and idx == rank:
            target_row = item
            break
    if target_row is None and 1 <= rank <= len(candidates) and isinstance(candidates[rank - 1], dict):
        # Last-resort fallback: list position. Helps with legacy payloads
        # whose rank field was never stamped.
        target_row = candidates[rank - 1]
    if target_row is None:
        raise HTTPException(status_code=404, detail=f"Rank {rank} out of range.")

    # Multi-pose mode preference: serve pose_pdb_all when requested AND
    # available, else fall back to the top pose. Frontend can detect by
    # counting MODEL records to enable/disable the pose-cycler controls.
    use_all = (mode or "first").lower() in {"all", "multi", "every"}
    pose = None
    if use_all:
        pose = target_row.get("pose_pdb_all")
    if not isinstance(pose, str) or not pose.strip():
        pose = target_row.get("pose_pdb")
    if not isinstance(pose, str) or not pose.strip():
        raise HTTPException(
            status_code=404,
            detail=(
                "No persisted Vina pose for this row. Either the docking "
                "failed for this molecule or the job pre-dates Phase B."
            ),
        )
    return PlainTextResponse(
        content=pose,
        media_type="chemical/x-pdb",
        headers={"Cache-Control": "public, max-age=86400"},
    )


@router.get("/jobs/{job_id}/export.csv")
async def legacy_job_export_csv(job_id: str, db: Session = Depends(get_db)):
    """Stream the candidate list of a denovo / admet / docking job as CSV.

    Researchers asked to round-trip candidate rows into Spotfire / ChemDraw
    / Excel — copy-pasting from the Results table loses the SMILES quoting
    and the floating-point precision. This endpoint pulls the latest
    job_output AuditLog row and emits a single sheet:

        rank,smiles,score,dock_score,qed,mw,logp,...

    Column set is the union of all keys across the rows so we don't drop
    descriptors that appear on some molecules but not others. Output
    filename is stamped with the short id so a chemist downloading three
    job CSVs in a row can tell them apart on disk.
    """
    import csv
    import io

    from fastapi.responses import StreamingResponse

    row = (
        db.query(AuditLog)
        .filter(
            AuditLog.resource_type == "job_output",
            AuditLog.resource_id == str(job_id),
            AuditLog.action.in_(
                ["denovo_output", "admet_output", "docking_output"]
            ),
        )
        .order_by(desc(AuditLog.created_at))
        .first()
    )
    if not row:
        raise HTTPException(
            status_code=404,
            detail="No exportable output for this job (denovo/admet/docking).",
        )
    try:
        payload = json.loads(row.payload or "{}")
    except Exception:  # noqa: BLE001
        payload = {}
    section, candidates = _candidates_from_output(payload)
    if not candidates:
        raise HTTPException(
            status_code=404,
            detail="Job output has no candidate rows to export.",
        )

    # Flatten nested dicts (e.g. row["descriptors"]) so each key becomes a
    # top-level column. Standard chemist-friendly order (rank, smiles first,
    # then numeric descriptors) — anything we don't recognise lands at the
    # tail in stable sorted order.
    preferred = [
        "rank", "smiles", "module", "score", "dock_score", "qed",
        "mw", "logp", "tpsa", "affinity",
    ]
    flat_rows: list[dict[str, object]] = []
    union_keys: set[str] = set()
    for item in candidates:
        if not isinstance(item, dict):
            continue
        flat: dict[str, object] = {}
        for k, v in item.items():
            if isinstance(v, dict):
                # one level of flattening: descriptors.qed -> qed
                for sk, sv in v.items():
                    flat.setdefault(sk, sv)
            elif isinstance(v, list):
                flat[k] = "|".join(str(x) for x in v)
            else:
                flat[k] = v
        flat_rows.append(flat)
        union_keys.update(flat.keys())
    if not flat_rows:
        raise HTTPException(status_code=404, detail="No rows to export.")

    cols = [c for c in preferred if c in union_keys] + sorted(
        k for k in union_keys if k not in preferred
    )

    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=cols, extrasaction="ignore")
    writer.writeheader()
    for r in flat_rows:
        writer.writerow({c: r.get(c, "") for c in cols})
    body = buf.getvalue().encode("utf-8-sig")  # BOM so Excel opens UTF-8 right

    short = str(job_id)[:8]
    filename = f"{section or 'job'}_{short}.csv"
    return StreamingResponse(
        iter([body]),
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
