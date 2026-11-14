"""Tests for /api/jobs/{id}/poses/{rank}.pdb (Phase B persistence).

The NGL viewer's "Vina pose vs RDKit fallback" badge depends on this
endpoint behaving cleanly in 4 cases:
    1. Rank has a real pose_pdb  →  200 + PDB block
    2. Rank exists but no pose   →  404 + clear message
    3. Rank out of range          →  404
    4. Unknown job                →  404
plus the ?mode=all extension that serves the multi-MODEL block when
available.
"""

from __future__ import annotations


def test_pose_rank_1_returns_persisted_pdb(client, seed_docking_job):
    """seed_docking_job stamps pose on ranks 1+3, skips rank 2."""
    job_id = seed_docking_job(with_pose=True)
    r = client.get(f"/api/jobs/{job_id}/poses/1.pdb")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("chemical/x-pdb")
    body = r.text
    # Sanity: it's a valid-looking PDB block.
    assert "HETATM" in body
    assert "END" in body


def test_pose_rank_without_pdb_returns_404(client, seed_docking_job):
    """Rank 2 is seeded without pose_pdb to exercise this path."""
    job_id = seed_docking_job(with_pose=True)
    r = client.get(f"/api/jobs/{job_id}/poses/2.pdb")
    assert r.status_code == 404
    assert "No persisted Vina pose" in r.json()["detail"]


def test_pose_rank_out_of_range_returns_404(client, seed_docking_job):
    """Asking for rank=99 when only 3 candidates exist must 404 with a
    "out of range" message — not 500 from an IndexError."""
    job_id = seed_docking_job(with_pose=True)
    r = client.get(f"/api/jobs/{job_id}/poses/99.pdb")
    assert r.status_code == 404
    assert "out of range" in r.json()["detail"].lower()


def test_pose_unknown_job_returns_404(client):
    r = client.get("/api/jobs/00000000-0000-0000-0000-000000000000/poses/1.pdb")
    assert r.status_code == 404
    assert "no docking output" in r.json()["detail"].lower()


def test_pose_mode_all_returns_multi_model(client, seed_docking_job):
    """When the row has pose_pdb_all (rank 1 in our seed), ?mode=all
    serves the multi-MODEL block. The frontend counts MODEL records to
    enable the pose cycler."""
    job_id = seed_docking_job(with_pose=True, with_multi_pose=True)
    r = client.get(f"/api/jobs/{job_id}/poses/1.pdb", params={"mode": "all"})
    assert r.status_code == 200
    body = r.text
    # Two MODELs means the cycler will activate; this is the contract.
    assert body.count("MODEL ") >= 2


def test_pose_mode_all_falls_back_to_top_pose(client, seed_docking_job):
    """When pose_pdb_all isn't persisted (rank 2 setup), ?mode=all
    silently falls back to the top pose so legacy data still renders."""
    # Rank 1 has top-pose but NO multi-pose in this seed config.
    job_id = seed_docking_job(with_pose=True, with_multi_pose=False)
    r = client.get(f"/api/jobs/{job_id}/poses/1.pdb", params={"mode": "all"})
    assert r.status_code == 200
    body = r.text
    # Single MODEL — no cycler activation, but still rendered.
    assert body.count("MODEL ") == 1
