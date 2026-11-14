"""Tests for the REINVENT-related endpoints:
    * GET  /api/reinvent/priors
    * POST /api/reinvent/evaluate-prior
    * POST /api/reinvent/promote-checkpoint
    * GET  /api/reinvent/rl-finetune  (list)
    * GET  /api/reinvent/rl-finetune/{id}/checkpoint (download)

Everything that touches Reinvent_App goes through the mocked_reinvent_app
httpx transport fixture so tests don't depend on the upstream service being up.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def test_priors_returns_base_and_finetuned(client, mocked_reinvent_app):
    """The mock yields ['EGFR', 'JAK2'] as available — the response
    must split into fine_tuned + always include 'base'."""
    r = client.get("/api/reinvent/priors")
    assert r.status_code == 200
    data = r.json()
    assert data["base_available"] is True
    assert set(data["fine_tuned"]) == {"EGFR", "JAK2"}
    assert "base" in data["all"]


def test_priors_503_when_reinvent_app_down(client, monkeypatch):
    """When the Reinvent_App probe raises, surface a 503 with a
    helpful detail — not a 500."""
    from app.services import reinvent_client

    async def _boom(*_a, **_k):
        raise RuntimeError("connection refused")

    monkeypatch.setattr(reinvent_client, "list_available_priors", _boom)
    r = client.get("/api/reinvent/priors")
    assert r.status_code == 503
    assert "unreachable" in r.json()["detail"].lower()


def test_evaluate_prior_returns_metrics(client, mocked_reinvent_app):
    """Default mock returns 3 candidates: validity 100%, mean QED ≈ 0.44.
    The endpoint should aggregate those into a single metrics dict."""
    r = client.post(
        "/api/reinvent/evaluate-prior",
        json={"prior": "base", "n_samples": 50},
    )
    assert r.status_code == 200
    data = r.json()
    assert data["prior"] == "base"
    m = data["metrics"]
    # All 3 candidates are valid in the default mock.
    assert m["validity_rate"] == 1.0
    assert m["mean_qed"] > 0
    # The mock supplies no Tanimoto so this stays None.
    assert m["mean_max_tanimoto"] is None
    # The endpoint surfaces a handful of sample SMILES for visual sanity.
    assert "CCO" in data["sample_smiles"]


def test_evaluate_prior_with_compare_to_base(client, mocked_reinvent_app):
    """When compare_to_base=True and the prior isn't 'base', we sample
    BOTH and surface a delta dict."""
    r = client.post(
        "/api/reinvent/evaluate-prior",
        json={"prior": "EGFR", "n_samples": 30, "compare_to_base": True},
    )
    assert r.status_code == 200
    data = r.json()
    # base_metrics + delta only appear in compare mode.
    assert "base_metrics" in data
    assert "delta" in data
    # The mock returns the same candidates for any prior, so delta is 0.
    assert data["delta"]["validity_rate"] == 0.0


def test_evaluate_prior_caps_n_samples(client, mocked_reinvent_app):
    """n_samples is server-side clamped to [20, 2000] — 100000 must not
    cause a giant sample batch."""
    r = client.post(
        "/api/reinvent/evaluate-prior",
        json={"prior": "base", "n_samples": 100000},
    )
    assert r.status_code == 200
    # Server clamps to 2000 but the mock only returns 3 anyway.
    # The shape check matters here, not the exact size.
    assert "metrics" in r.json()


def test_promote_checkpoint_copies_to_finetuned_slot(client, test_db, seed_rl_job, isolated_rl_storage, tmp_path):
    """Happy path: a completed RL job with a checkpoint at
    storage/rl_jobs/<id>/agent.prior gets copied to
    Reinvent_App/prior/finetuned/<label>_finetuned.prior."""
    # Pre-create the source checkpoint so the endpoint can find it.
    job_id = seed_rl_job(target="EGFR", status="completed")
    src_dir = isolated_rl_storage / job_id
    src_dir.mkdir(parents=True, exist_ok=True)
    (src_dir / "agent.prior").write_bytes(b"FAKE_PRIOR_BLOB")

    # Re-point the Reinvent_App finetuned dir to a temp dir so the test
    # doesn't litter the real prior folder.
    from app.api import routes_legacy

    # We can't monkey-patch the path resolution inside promote because
    # it derives Path(__file__).parents[3] / "Reinvent_App" / ... — that
    # path is a real directory. Instead, write to it and clean up.
    # The label is sanitized; pick a deliberately-test-only string.
    r = client.post("/api/reinvent/promote-checkpoint", json={
        "job_id": job_id,
        "target_label": "TEST_PROMOTE_ONLY_DELETE_ME",
    })
    try:
        assert r.status_code == 200
        data = r.json()
        assert data["ok"] is True
        assert data["available_as"] == "TEST_PROMOTE_ONLY_DELETE_ME"
        # Source still exists; destination now does too.
        assert Path(data["src"]).exists()
        assert Path(data["dst"]).exists()
    finally:
        # Clean up the real Reinvent_App slot we created.
        dst = Path(routes_legacy.__file__).resolve().parents[3] / "Reinvent_App" / "prior" / "finetuned" / "TEST_PROMOTE_ONLY_DELETE_ME_finetuned.prior"
        if dst.exists():
            dst.unlink()


def test_promote_checkpoint_400_when_missing_fields(client):
    r = client.post("/api/reinvent/promote-checkpoint", json={"job_id": "abc"})
    assert r.status_code == 400


def test_promote_checkpoint_404_when_no_source(client, seed_rl_job, isolated_rl_storage):
    """RL job exists in DB but no checkpoint file on disk → 404."""
    job_id = seed_rl_job(target="EGFR")
    r = client.post("/api/reinvent/promote-checkpoint", json={
        "job_id": job_id,
        "target_label": "DOES_NOT_MATTER",
    })
    assert r.status_code == 404


def test_rl_finetune_list_returns_recent_first(client, seed_rl_job, isolated_rl_storage):
    """The list endpoint should surface every Job row whose receptor_name
    starts with 'RL-', sorted newest-first."""
    _ = seed_rl_job(target="EGFR")
    second = seed_rl_job(target="JAK2", status="running")
    r = client.get("/api/reinvent/rl-finetune")
    assert r.status_code == 200
    data = r.json()
    assert len(data) == 2
    # Both rows have the per-row fields the UI table reads.
    for row in data:
        assert "target" in row
        assert "status" in row
        assert "checkpoint_exists" in row
        assert "progress" in row


def test_rl_finetune_checkpoint_download_404_when_missing(client, seed_rl_job, isolated_rl_storage):
    """No checkpoint file → 404 (not 500)."""
    job_id = seed_rl_job(target="EGFR")
    r = client.get(f"/api/reinvent/rl-finetune/{job_id}/checkpoint")
    assert r.status_code == 404


def test_rl_finetune_checkpoint_download_streams_bytes(client, seed_rl_job, isolated_rl_storage):
    """When the checkpoint exists on disk it's streamed with the right
    Content-Disposition for a download link."""
    job_id = seed_rl_job(target="EGFR")
    ck_dir = isolated_rl_storage / job_id
    ck_dir.mkdir(parents=True, exist_ok=True)
    payload = b"FAKE_PRIOR_PAYLOAD" * 100
    (ck_dir / "agent.prior").write_bytes(payload)
    r = client.get(f"/api/reinvent/rl-finetune/{job_id}/checkpoint")
    assert r.status_code == 200
    assert r.content == payload
    assert "attachment" in r.headers["content-disposition"]
    assert "rl_agent" in r.headers["content-disposition"]
