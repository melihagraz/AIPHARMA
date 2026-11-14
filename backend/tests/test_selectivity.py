"""Tests for /api/selectivity/cross-dock (sync + async).

Cross-dock makes real HTTP calls into Docking_App so every test here
goes through the mocked_docking_app fixture. We exercise:
    * 400 on empty SMILES
    * 503 when Docking_App reports no targets
    * sync mode returns the selectivity_delta row summaries
    * async mode returns a job_id + status URL, lets a polling caller
      see the matrix grow column-by-column as targets complete
"""

from __future__ import annotations

import pytest
from httpx import Response


def test_cross_dock_400_when_no_smiles(client):
    """Empty SMILES + no job_id is a client error, not a 5xx."""
    r = client.post("/api/selectivity/cross-dock", json={"top_k": 1})
    assert r.status_code == 400
    assert "no smiles" in r.json()["detail"].lower()


def test_cross_dock_400_for_async_variant_too(client):
    """Same validation on the async path so callers see the same error
    surface whether they use sync or async."""
    r = client.post("/api/selectivity/cross-dock-async", json={"smiles": [], "top_k": 1})
    assert r.status_code == 400


def test_cross_dock_503_when_no_targets(client, mocked_docking_app):
    """A Docking_App with no bundled receptors should surface as 503."""
    mocked_docking_app.get("/targets").mock(return_value=Response(200, json={"targets": []}))
    r = client.post("/api/selectivity/cross-dock", json={"smiles": ["CCO"], "top_k": 1})
    assert r.status_code == 503


def test_cross_dock_sync_returns_matrix(client, mocked_docking_app):
    """Sync endpoint runs the full cross-dock inline. mocked_docking_app
    returns a canned -7.5 affinity for CCO so every target ends up
    with the same value — and selectivity_delta is therefore 0.0."""
    r = client.post("/api/selectivity/cross-dock", json={
        "smiles": ["CCO"], "top_k": 1, "primary_target": "EGFR",
    })
    assert r.status_code == 200
    data = r.json()
    assert data["status"] == "completed"
    assert "EGFR" in data["targets"]
    # One row per SMILES.
    assert len(data["data"]) == 1
    row = data["data"][0]
    assert row["smiles"] == "CCO"
    # Every target hit the same -7.5 mock so primary == best off-target;
    # selectivity_delta = 0.0 — a sanity check that the math runs at all.
    assert row["primary_affinity"] == pytest.approx(-7.5)
    assert row["selectivity_delta"] == pytest.approx(0.0)


def test_cross_dock_async_returns_job_id_immediately(client, mocked_docking_app):
    """Async path returns 200 + job_id without waiting for docking to
    finish. Mock the Docking_App /dock to return immediately so the
    background task drains quickly within the test session."""
    r = client.post("/api/selectivity/cross-dock-async", json={
        "smiles": ["CCO"], "top_k": 1, "primary_target": "EGFR",
    })
    assert r.status_code == 200
    data = r.json()
    assert "job_id" in data
    assert data["status"] == "running"
    assert data["status_url"].startswith("/api/selectivity/cross-dock/")


async def test_cross_dock_async_status_polling(async_client, mocked_docking_app, monkeypatch):
    """After firing the async job we should be able to poll the status
    endpoint and see progress (completed_targets, current_target).

    Implementation note on the test scaffolding:
      The POST handler kicks off the worker via asyncio.create_task,
      which on FastAPI's TestClient orphans (the loop closes when the
      response returns). We patch create_task to a no-op for this test
      so the POST sets up the runtime state cleanly, then we drive the
      worker ourselves with a single direct await. Result: the test
      observes exactly one full pass through every target.
    """
    import asyncio as _asyncio

    from app.api import routes_legacy

    captured: dict = {}

    def _capture_no_run(coro):
        # Pull out the coroutine so we can await it explicitly later,
        # and silence the "coroutine was never awaited" warning by
        # closing the duplicate one the test won't drive.
        captured.setdefault("coro", coro)
        # Return a benign sentinel that satisfies asyncio.create_task's
        # callers (they only use it via _RL_RUNTIME book-keeping which
        # selectivity does not touch).
        async def _noop():
            return None
        return _asyncio.ensure_future(_noop())

    monkeypatch.setattr(routes_legacy.asyncio, "create_task", _capture_no_run)

    start = (await async_client.post(
        "/api/selectivity/cross-dock-async",
        json={"smiles": ["CCO"], "top_k": 1, "primary_target": "EGFR"},
    )).json()
    job_id = start["job_id"]
    assert job_id in routes_legacy._SEL_RUNTIME
    # Drive the captured worker coroutine to completion.
    await captured["coro"]
    s = (await async_client.get(f"/api/selectivity/cross-dock/{job_id}")).json()
    assert s["status"] == "completed"
    assert s["progress_pct"] == 100
    assert "data" in s and len(s["data"]) == 1


def test_cross_dock_async_status_404_for_unknown_job(client):
    r = client.get("/api/selectivity/cross-dock/00000000-0000-0000-0000-000000000000")
    assert r.status_code == 404
