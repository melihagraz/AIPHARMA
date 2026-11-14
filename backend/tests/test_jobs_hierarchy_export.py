"""Tests for /api/jobs/hierarchy + /api/jobs/{id}/export.csv.

These are the "what shows up in the UI dropdowns" endpoints. The contract
the frontend depends on:
    * hierarchy returns total_generated / admet_passed / candidates_saved
      for denovo/docking rows (drives the Step 4/5/6 dropdown enrichment)
    * export.csv is BOM-prefixed UTF-8 (Excel reads it correctly), has
      preferred-column ordering (rank, smiles first), and 404s cleanly
      when there's no docking/denovo/admet output to export
"""

from __future__ import annotations

import csv
import io


def test_hierarchy_empty_when_no_jobs(client):
    r = client.get("/api/jobs/hierarchy")
    assert r.status_code == 200
    assert r.json() == []


def test_hierarchy_returns_denovo_pass_counts(client, seed_denovo_job):
    """Step 6 dropdown depends on these three fields being non-None for
    denovo rows. Regression for the 'all candidates show as empty in
    the picker' bug from the Step 6 chain-closing work."""
    job_id = seed_denovo_job(target="EGFR", n_candidates=7, total_generated=20)
    r = client.get("/api/jobs/hierarchy")
    assert r.status_code == 200
    data = r.json()
    assert len(data) == 1
    row = data[0]
    assert row["job_id"] == job_id
    assert row["analysis_section"] == "denovo"
    assert row["total_generated"] == 20
    assert row["admet_passed"] == 7  # after_filter_count from payload
    assert row["candidates_saved"] == 7
    assert row["short_id"] == job_id[:8]


def test_hierarchy_returns_docking_rows(client, seed_docking_job):
    """Docking jobs surface in hierarchy with analysis_section='docking'."""
    job_id = seed_docking_job(target="EGFR")
    r = client.get("/api/jobs/hierarchy")
    assert r.status_code == 200
    rows = r.json()
    docking = [row for row in rows if row["job_id"] == job_id]
    assert len(docking) == 1
    assert docking[0]["analysis_section"] == "docking"
    assert docking[0]["candidates_saved"] == 3


def test_export_csv_returns_bom_prefixed_utf8(client, seed_denovo_job):
    """Excel needs the UTF-8 BOM at byte 0 to render UTF-8 correctly."""
    job_id = seed_denovo_job(target="EGFR", n_candidates=4)
    r = client.get(f"/api/jobs/{job_id}/export.csv")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/csv")
    assert r.content.startswith(b"\xef\xbb\xbf")
    assert "attachment" in r.headers.get("content-disposition", "")


def test_export_csv_preferred_column_order(client, seed_denovo_job):
    """rank + smiles must come first so a chemist scanning the file at a
    glance sees the identifier before the descriptor block."""
    job_id = seed_denovo_job(target="EGFR", n_candidates=3)
    r = client.get(f"/api/jobs/{job_id}/export.csv")
    # Strip BOM, parse, assert header order.
    text = r.content.decode("utf-8-sig")
    reader = csv.reader(io.StringIO(text))
    header = next(reader)
    assert header[0] == "rank"
    assert header[1] == "smiles"
    # Numeric descriptors land in the preferred slots before the alpha tail.
    for k in ("score", "qed", "mw", "logp"):
        assert k in header


def test_export_csv_404_for_unknown_job(client):
    """No audit row → 404 + a clear detail message (not 500)."""
    r = client.get("/api/jobs/00000000-0000-0000-0000-000000000000/export.csv")
    assert r.status_code == 404
    assert "exportable output" in r.json()["detail"].lower()


def test_export_csv_has_one_row_per_candidate(client, seed_denovo_job):
    """Sanity: the file row count matches the seeded candidate count."""
    job_id = seed_denovo_job(target="JAK2", n_candidates=12)
    r = client.get(f"/api/jobs/{job_id}/export.csv")
    text = r.content.decode("utf-8-sig")
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)
    # 1 header + 12 data rows = 13.
    assert len(rows) == 13


def test_export_csv_flattens_nested_descriptors(client, seed_denovo_job):
    """Descriptors live in row['descriptors'] in the payload. The CSV
    flattens one level so descriptors.qed lands as plain 'qed' column —
    not as 'descriptors.qed' or stringified JSON."""
    job_id = seed_denovo_job(target="EGFR", n_candidates=2)
    r = client.get(f"/api/jobs/{job_id}/export.csv")
    text = r.content.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    rows = list(reader)
    assert len(rows) == 2
    # Each row carries the flat scalar form, not a serialised dict.
    assert "{" not in rows[0]["qed"]
    assert float(rows[0]["qed"]) > 0
