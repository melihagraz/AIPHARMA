"""Phase 5 fail-closed tests. These are security/data-integrity tests,
not tests of biological accuracy or live ADMET/Vina inference.
"""
from __future__ import annotations

import importlib
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


@pytest.mark.parametrize("module", [
    "Reinvent_App.app.main",
    "Admet_App.app.main",
    "Docking_App.app.main",
    "Pocket_App.app.main",
])
def test_microservices_reject_unrecognized_environment(module):
    """Isolate each service to avoid REINVENT pickle import aliases leaking.

    Test all deployment-typo variants against its actual HTTP middleware.
    """
    import subprocess
    import sys
    code = r'''
import importlib, os
from fastapi.testclient import TestClient
app = importlib.import_module('{module}').app
for bad in ('prod', 'staging', 'production-with-typo'):
    os.environ['APP_ENV'] = bad
    with TestClient(app) as client:
        r = client.get('/health')
    assert r.status_code == 503, (bad, r.status_code, r.text)
    assert r.json()['detail'] == 'Invalid APP_ENV configuration'
'''
    completed = subprocess.run(
        [sys.executable, "-c", code.format(module=module)],
        cwd=str(Path(__file__).resolve().parents[2]),
        capture_output=True, text=True, timeout=30,
    )
    assert completed.returncode == 0, completed.stderr[-2000:] + completed.stdout[-2000:]


@pytest.mark.parametrize("malicious", [
    "../../tmp", "/etc/passwd", "..", "fake", "",
    "11111111-1111-1111-1111-111111111111/../../",
])
def test_rl_storage_rejects_non_uuid_paths(malicious):
    from app.api.routes_legacy import _rl_job_paths
    with pytest.raises(ValueError, match="Invalid RL job identifier"):
        _rl_job_paths(malicious)


def test_rl_storage_accepts_canonical_uuid():
    from app.api.routes_legacy import _rl_job_paths
    canonical = "11111111-1111-4111-8111-111111111111"
    paths = _rl_job_paths(canonical)
    assert len(paths) == 4
    assert paths[0].name == "agent.prior"
    assert paths[3].name == canonical


def test_promote_rejects_path_traversal_before_io(client):
    resp = client.post("/api/reinvent/promote-checkpoint", json={
        "job_id": "../../arbitrary", "target_label": "BAD"
    })
    assert resp.status_code == 400
    assert "Invalid RL job identifier" in resp.json()["detail"]


def test_promote_rejects_incomplete_training(client, seed_rl_job, isolated_rl_storage):
    job_id = seed_rl_job(status="running")
    target = isolated_rl_storage / job_id
    target.mkdir(parents=True, exist_ok=True)
    (target / "agent.prior").write_bytes(b"checkpoint")
    resp = client.post("/api/reinvent/promote-checkpoint", json={
        "job_id": job_id, "target_label": "NO_PROMOTE_RUNNING"
    })
    assert resp.status_code == 409


def test_promote_rejects_symlink_checkpoint(client, seed_rl_job, isolated_rl_storage, tmp_path):
    job_id = seed_rl_job(status="completed")
    target = isolated_rl_storage / job_id
    target.mkdir(parents=True, exist_ok=True)
    real_file = tmp_path / "untrusted.prior"
    real_file.write_bytes(b"not trusted")
    (target / "agent.prior").symlink_to(real_file)
    resp = client.post("/api/reinvent/promote-checkpoint", json={
        "job_id": job_id, "target_label": "NO_PROMOTE_SYMLINK"
    })
    assert resp.status_code == 404


@pytest.mark.parametrize("module", [
    "Reinvent_App.app.main",
    "Admet_App.app.main",
    "Docking_App.app.main",
    "Pocket_App.app.main",
])
def test_microservice_production_key_blocks_unauthorized_requests(module):
    """Actual FastAPI middleware: valid internal key is required in production."""
    import subprocess
    import sys
    code = r'''
import importlib, os
from fastapi.testclient import TestClient
app = importlib.import_module('{module}').app
os.environ['APP_ENV'] = 'production'
os.environ['INTERNAL_SERVICE_KEY'] = 'test-key-' * 6
with TestClient(app) as client:
    missing = client.get('/openapi.json')
    incorrect = client.get('/openapi.json', headers={{'X-Internal-Service-Key':'invalid'}})
    correct = client.get('/openapi.json', headers={{'X-Internal-Service-Key':'test-key-' * 6}})
assert missing.status_code == 401, (missing.status_code, missing.text)
assert incorrect.status_code == 401, (incorrect.status_code, incorrect.text)
assert correct.status_code == 200, (correct.status_code, correct.text)
'''
    completed = subprocess.run([sys.executable, '-c', code.format(module=module)],
                               cwd=str(Path(__file__).resolve().parents[2]),
                               capture_output=True,text=True,timeout=30)
    assert completed.returncode == 0, completed.stderr[-2000:] + completed.stdout[-2000:]
