"""Shared pytest fixtures for the backend test suite.

Design goals:
    * Hermetic by default — every test gets a fresh in-memory SQLite, no
      writes to backend/enterprise.db.
    * No external service required — Docking_App / Reinvent_App / Admet_App
      HTTP calls are intercepted with httpx.MockTransport and given canned responses
      that mirror the real schemas (DockResponse, SampleResponse, etc.).
    * Background tasks short-circuited — asyncio.create_subprocess_exec
      and asyncio.create_task patches let RL fine-tune / selectivity
      async endpoints return without actually spawning a process or
      hammering Docking_App.
    * Cheap to compose — seed_* helpers exist as fixtures and the test
      writes a single line ``job = seed_completed_denovo_job(jobs_db, ...)``
      to get a fully-shaped job row + AuditLog payload, instead of every
      test re-implementing the same INSERT statements.

Fixture layering (highest first; later layers depend on earlier):
    test_engine     — SQLAlchemy engine bound to in-memory SQLite + schema
    test_db_session — single Session for the test (rolled back in teardown)
    client          — FastAPI TestClient with get_db dependency overridden
    mocked_*        — httpx routers for each upstream service
"""

from __future__ import annotations

import json
import os
import sys
import uuid
from collections.abc import Iterator
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pytest

# ---------------------------------------------------------------------------
# Path setup — let tests `import app.*` from backend/app without forcing
# the user to set PYTHONPATH manually. Mirrors how uvicorn loads the app
# via --app-dir backend.
# ---------------------------------------------------------------------------
_BACKEND_ROOT = Path(__file__).resolve().parents[1]
if str(_BACKEND_ROOT) not in sys.path:
    sys.path.insert(0, str(_BACKEND_ROOT))

# Force a known database URL BEFORE app.* gets imported so the engine
# is built against the in-memory DB, not whatever .env / settings.py
# resolves to in dev. The test_engine fixture below also overrides
# get_db, but this belt-and-braces approach prevents the import-time
# engine from grabbing the dev enterprise.db when settings are read.
os.environ.setdefault("DATABASE_URL", "sqlite:///:memory:")


@pytest.fixture(scope="session")
def app():
    """The FastAPI app instance — imported once per session.

    Importing inside the fixture (not at module top) keeps pytest
    collection fast when filtering down to a subset of tests, and means
    a misconfigured app doesn't crash test discovery — it raises during
    the first test that requests this fixture.
    """
    from app.main import app as fastapi_app  # type: ignore[no-redef]

    return fastapi_app


@pytest.fixture()
def test_engine(monkeypatch):
    """Per-test in-memory SQLite engine with the full schema created.

    Each test gets a fresh engine (function scope) so audit log rows
    from one test don't leak into another. We patch the GLOBAL
    ``app.db.session.SessionLocal`` to bind to this engine — that's
    necessary because some background workers (e.g.
    routes_legacy._run_cross_dock, _rl_subprocess_runner) instantiate
    SessionLocal directly instead of going through get_db, so the
    FastAPI dependency override in the client fixture isn't enough
    on its own.
    """
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.pool import StaticPool

    from app import models  # noqa: F401  trigger metadata registration
    from app.db import session as session_module
    from app.db.base import Base

    # StaticPool so every Session sees the same in-memory DB. Without it
    # SQLite spins up a per-connection DB and tables vanish between the
    # create_all and the first query.
    engine = create_engine(
        "sqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    Base.metadata.create_all(bind=engine)
    # Redirect the global SessionLocal so background workers that bypass
    # get_db (selectivity worker, RL subprocess supervisor) write into
    # this test's DB rather than the dev enterprise.db.
    test_session_factory = sessionmaker(bind=engine, autocommit=False, autoflush=False)
    monkeypatch.setattr(session_module, "SessionLocal", test_session_factory)
    monkeypatch.setattr(session_module, "engine", engine)
    yield engine
    engine.dispose()


@pytest.fixture()
def test_db(test_engine):
    """Yield a SQLAlchemy Session bound to the test engine.

    Tests that want to insert rows directly (seed fixtures, etc.) use
    this. The session is closed on teardown — no rollback semantics are
    enforced here because the engine itself is per-test, so leaks are
    bounded to the test that produced them.
    """
    from sqlalchemy.orm import sessionmaker

    Session = sessionmaker(bind=test_engine, autocommit=False, autoflush=False)
    db = Session()
    try:
        yield db
    finally:
        db.close()


@pytest.fixture()
async def async_client(app, test_engine):
    """Async httpx client with ASGITransport — the only way to test
    endpoints that fire-and-forget via asyncio.create_task.

    TestClient creates a fresh event loop per request; any background
    task scheduled inside the request handler gets orphaned when that
    loop closes. With ASGITransport on a long-lived AsyncClient, every
    request shares the test's loop, so create_task'd workers actually
    drain. Use this for /selectivity/cross-dock-async,
    /reinvent/rl-finetune POST, and any other endpoint whose effects
    you need to observe across multiple HTTP calls.
    """
    import httpx
    from sqlalchemy.orm import sessionmaker

    from app.db.session import get_db

    Session = sessionmaker(bind=test_engine, autocommit=False, autoflush=False)

    def _override():
        s = Session()
        try:
            yield s
        finally:
            s.close()

    app.dependency_overrides[get_db] = _override
    transport = httpx.ASGITransport(app=app)
    try:
        async with httpx.AsyncClient(transport=transport, base_url="http://test") as ac:
            yield ac
    finally:
        app.dependency_overrides.pop(get_db, None)


@pytest.fixture()
def client(app, test_engine):
    """FastAPI TestClient with the get_db dependency overridden.

    Every request inside the client gets a fresh Session bound to the
    same in-memory engine, so HTTP-level tests see the data that fixtures
    have already INSERT'ed via test_db.
    """
    from fastapi.testclient import TestClient
    from sqlalchemy.orm import sessionmaker

    from app.db.session import get_db

    Session = sessionmaker(bind=test_engine, autocommit=False, autoflush=False)

    def _override():
        s = Session()
        try:
            yield s
        finally:
            s.close()

    app.dependency_overrides[get_db] = _override
    try:
        with TestClient(app) as tc:
            yield tc
    finally:
        app.dependency_overrides.pop(get_db, None)


# ---------------------------------------------------------------------------
# Seed helpers — composable fixtures that insert canonical Job + AuditLog
# rows for the various analysis sections (denovo, docking, selectivity,
# rl_finetune). Each returns the job_id string so tests can chain.
# ---------------------------------------------------------------------------


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(sep=" ", timespec="seconds")


@pytest.fixture()
def seed_denovo_job(test_db):
    """Insert a completed De Novo job with N candidate rows.

    Usage::

        def test_x(client, seed_denovo_job):
            job_id = seed_denovo_job(target="EGFR", n_candidates=10)
            r = client.get(f"/api/jobs/{job_id}/export.csv")
            assert r.status_code == 200
    """
    from app.models.audit_log import AuditLog
    from app.models.job import Job

    def _seed(target: str = "EGFR", n_candidates: int = 5, total_generated: int = 10) -> str:
        job_id = str(uuid.uuid4())
        now = datetime.now(timezone.utc)
        job = Job(
            id=job_id,
            tenant_id="default",
            project_id="legacy-ui",
            status="completed",
            mode="agonist",
            receptor_name=target,
            pdb_id=target,
            uniprot_id=None,
            fasta_sequence=None,
            candidate_count=n_candidates,
            progress_pct=100,
            current_step="completed",
            progress_message="done",
            pipeline_resolution="embedded",
            error_message=None,
            created_at=now,
            updated_at=now,
        )
        test_db.add(job)
        candidates = [
            {
                "rank": i + 1,
                "smiles": f"C{'C' * i}O",
                "score": round(0.9 - i * 0.05, 3),
                "dock_score": -7.5 + i * 0.3,
                "qed": round(0.7 - i * 0.02, 3),
                "mw": 100 + i * 30,
                "logp": 1.5 + i * 0.2,
                "descriptors": {"qed": round(0.7 - i * 0.02, 3), "mw": 100 + i * 30, "logp": 1.5 + i * 0.2},
                "alerts": [],
            }
            for i in range(n_candidates)
        ]
        payload = {
            "analysis_section": "denovo",
            "target": target,
            "total_generated": total_generated,
            "after_filter_count": n_candidates,
            "results_json": {
                "target": target,
                "total_generated": total_generated,
                "after_filter_count": n_candidates,
                "data": candidates,
            },
        }
        test_db.add(
            AuditLog(
                id=str(uuid.uuid4()),
                actor_id="system",
                action="denovo_output",
                resource_type="job_output",
                resource_id=job_id,
                payload=json.dumps(payload, ensure_ascii=True),
                created_at=now,
            )
        )
        test_db.commit()
        return job_id

    return _seed


@pytest.fixture()
def seed_docking_job(test_db):
    """Insert a completed docking job with optional pose_pdb / pose_pdb_all per row."""
    from app.models.audit_log import AuditLog
    from app.models.job import Job

    # Minimal vanilla PDB MODEL block — enough for the poses endpoint test
    # to read non-empty content back out. Real Vina blocks are larger but
    # shape-only tests don't need that.
    _POSE_SINGLE = (
        "MODEL        1\n"
        "HETATM    1  C1  LIG A   1      12.345  56.789   1.234  1.00  0.00           C\n"
        "ENDMDL\nEND"
    )
    _POSE_MULTI = (
        "MODEL        1\n"
        "HETATM    1  C1  LIG A   1      12.345  56.789   1.234  1.00  0.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "HETATM    1  C1  LIG A   1      13.000  57.000   1.500  1.00  0.00           C\n"
        "ENDMDL\nEND"
    )

    def _seed(target: str = "EGFR", with_pose: bool = True, with_multi_pose: bool = False) -> str:
        job_id = str(uuid.uuid4())
        now = datetime.now(timezone.utc)
        job = Job(
            id=job_id,
            tenant_id="default",
            project_id="legacy-ui",
            status="completed",
            mode="agonist",
            receptor_name=f"DOCK-{target}",
            pdb_id=f"DOCK-{target}",
            uniprot_id=None,
            fasta_sequence=None,
            candidate_count=3,
            progress_pct=100,
            current_step="completed",
            progress_message="done",
            pipeline_resolution="embedded",
            error_message=None,
            created_at=now,
            updated_at=now,
        )
        test_db.add(job)
        rows: list[dict[str, Any]] = []
        for i in range(3):
            row: dict[str, Any] = {
                "rank": i + 1,
                "smiles": f"C{'C' * i}O",
                "module": "general",
                "score": None,
                "dock_score": -8.0 + i,
                "affinity": -8.0 + i,
                "mw": 100 + i * 30,
                "logp": 1.5 + i * 0.2,
                "alerts": [],
            }
            # Sprinkle pose on rank 1 + 3, skip rank 2 to test "missing pose" path.
            if with_pose and i != 1:
                row["pose_pdb"] = _POSE_SINGLE
            if with_multi_pose and i == 0:
                row["pose_pdb_all"] = _POSE_MULTI
            rows.append(row)
        payload = {
            "analysis_section": "docking",
            "target": target,
            "results_json": {
                "section": "docking",
                "target": target,
                "input_count": 3,
                "docked_count": 3,
                "succeeded": 3,
                "data": rows,
            },
        }
        test_db.add(
            AuditLog(
                id=str(uuid.uuid4()),
                actor_id="system",
                action="docking_output",
                resource_type="job_output",
                resource_id=job_id,
                payload=json.dumps(payload, ensure_ascii=True),
                created_at=now,
            )
        )
        test_db.commit()
        return job_id

    return _seed


@pytest.fixture()
def seed_rl_job(test_db):
    """Insert an RL fine-tune Job row (without spawning a real subprocess)."""
    from app.models.job import Job

    def _seed(target: str = "EGFR", status: str = "completed") -> str:
        job_id = str(uuid.uuid4())
        now = datetime.now(timezone.utc)
        test_db.add(
            Job(
                id=job_id,
                tenant_id="default",
                project_id="legacy-ui",
                status=status,
                mode="agonist",
                receptor_name=f"RL-{target}-{job_id[:8]}",
                pdb_id=f"RL-{target}-{job_id[:8]}",
                uniprot_id=None,
                fasta_sequence=None,
                candidate_count=50,
                progress_pct=100 if status == "completed" else 0,
                current_step=status,
                progress_message=None,
                pipeline_resolution="embedded",
                error_message=None,
                created_at=now,
                updated_at=now,
            )
        )
        test_db.commit()
        return job_id

    return _seed


# ---------------------------------------------------------------------------
# External-service mocks — built-in httpx.MockTransport for AsyncClient
# calls. Each test that needs one of these as_dep'ing a fixture, no
# need for autouse patches that hide what's being mocked.
# ---------------------------------------------------------------------------


@pytest.fixture()
def mocked_reinvent_app(monkeypatch):
    """httpx.MockTransport context that fakes Reinvent_App /priors + /sample.

    The /sample stub is request-aware: when target_prior matches a name
    in the (mocked) prior catalog we echo it back as prior_used so the
    eval endpoint's "compare_to_base only fires when prior_used !=
    'base'" branch behaves like it would against the real service.

    Yields a small HTTPX test router so a test can override individual routes
    if it needs a custom response (e.g. mid-test failure injection).
    """
    import json as _json

    from tests.httpx_stub import install
    from httpx import Response

    _AVAILABLE = {"EGFR", "JAK2"}
    _CANDIDATES = [
        {"smiles": "CCO", "nll": 12.0, "valid": True, "mw": 46.07, "logp": -0.31, "qed": 0.42, "tpsa": 20.2},
        {"smiles": "CCN", "nll": 13.5, "valid": True, "mw": 45.08, "logp": -0.13, "qed": 0.39, "tpsa": 26.0},
        {"smiles": "c1ccccc1", "nll": 11.2, "valid": True, "mw": 78.11, "logp": 1.69, "qed": 0.51, "tpsa": 0.0},
    ]

    def _sample_side_effect(request):
        try:
            body = _json.loads(request.content or b"{}")
        except Exception:  # noqa: BLE001
            body = {}
        requested = (body.get("target_prior") or "").strip()
        if requested and requested.upper() in _AVAILABLE:
            prior_used = requested.upper()
        else:
            prior_used = "base"
        return Response(200, json={
            "requested": int(body.get("batch_size", len(_CANDIDATES))),
            "returned": len(_CANDIDATES),
            "valid_count": len(_CANDIDATES),
            "candidates": _CANDIDATES,
            "prior_vocab_size": 34,
            "prior_max_seq_len": 256,
            "prior_used": prior_used,
        })

    mock = install(monkeypatch, "http://127.0.0.1:8020")
    mock.get("/priors").mock(
        return_value=Response(200, json={"available": sorted(_AVAILABLE)})
    )
    mock.post("/sample").mock(side_effect=_sample_side_effect)
    yield mock


@pytest.fixture()
def mocked_docking_app(monkeypatch):
    """HTTPX transport that fakes Docking_App /targets + /dock + /health."""
    from tests.httpx_stub import install
    from httpx import Response

    mock = install(monkeypatch, "http://127.0.0.1:8040")
    mock.get("/targets").mock(
        return_value=Response(200, json={"targets": ["EGFR", "JAK2", "AKT1"]})
    )
    mock.get("/health").mock(
        return_value=Response(200, json={"status": "ok", "vina_binary": True, "available_targets": ["EGFR"]})
    )
    # Service responses must match the requested target, number and order of
    # SMILES; otherwise the production fail-closed data-integrity guard
    # correctly rejects them. This exercise used to be skipped when respx
    # was unavailable and had a mismatched canned response.
    def _dock_side_effect(request):
        body = json.loads(request.content)
        target = body["target"]
        results = []
        for smiles in body["smiles"]:
            if smiles == "CCO":
                results.append({"smiles": smiles, "target": target,
                                "affinity": -7.5, "num_modes": 5,
                                "pose_pdb": "MODEL 1\nHETATM 1 C 1 0.0 0.0 0.0 1.00\nENDMDL\nEND"})
            else:
                results.append({"smiles": smiles, "target": target,
                                "affinity": None, "error": "ligand_prep_failed", "num_modes": 0})
        return Response(200, json={
            "target": target,
            "count": len(results),
            "success_count": sum(r["affinity"] is not None for r in results),
            "engine": "httpx-test-stub",
            "results": results,
        })

    mock.post("/dock").mock(side_effect=_dock_side_effect)
    yield mock


@pytest.fixture()
def isolated_rl_storage(monkeypatch, tmp_path):
    """Re-point storage/rl_jobs at tmp_path so promote / checkpoint tests
    don't litter the real backend/storage directory."""
    from app.api import routes_legacy

    monkeypatch.setattr(routes_legacy, "_RL_JOBS_DIR", tmp_path / "rl_jobs")
    return tmp_path / "rl_jobs"
