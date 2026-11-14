"""Tests for the extended receptor catalog + on-demand RCSB prepare flow.
The catalog itself is a pure-data module (no network) — tests assert the
contract the frontend depends on (75 entries, family histogram, alias
lookup). The prepare endpoint hits RCSB so we mock it via httpx.MockTransport, plus
check the idempotent cache short-circuit and the error paths the UI
surfaces inline.
"""
from __future__ import annotations
import pytest
from tests.httpx_stub import install
from httpx import Response
def test_catalog_list_returns_75_entries_grouped_by_family(client):
    """Contract test for the dropdown render. If the count or family
    histogram changes the UI optgroup labels need updating in lockstep."""
    r = client.get("/api/receptors/extended-catalog")
    assert r.status_code == 200
    data = r.json()
    assert data["count"] == 75
    # Every entry has the four fields the UI option label uses.
    for e in data["entries"]:
        assert {"name", "pdb", "gene", "family", "aliases"} <= set(e.keys())
    # Family count totals must match entry count — guards against a
    # family being added in FAMILIES without any catalog entries.
    fam_total = sum(f["count"] for f in data["families"])
    assert fam_total == data["count"]
def test_catalog_lookup_resolves_canonical_name(client):
    r = client.get("/api/receptors/extended-catalog/HER2")
    assert r.status_code == 200
    data = r.json()
    assert data["name"] == "HER2"
    assert data["gene"] == "ERBB2"
    assert data["pdb"] == "3PP0"
def test_catalog_lookup_resolves_aliases_case_insensitive(client):
    """Alias resolution is the whole point of having the catalog —
    the user types whatever name they remember (NEU, erbb2, neu, HER2)
    and they all land on the same entry."""
    for q in ["HER2", "her2", "ERBB2", "erbb2", "NEU", "neu"]:
        r = client.get(f"/api/receptors/extended-catalog/{q}")
        assert r.status_code == 200, q
        assert r.json()["name"] == "HER2", q
def test_catalog_lookup_404_for_unknown(client):
    r = client.get("/api/receptors/extended-catalog/NotARealReceptor")
    assert r.status_code == 404
    assert "NotARealReceptor" in r.json()["detail"]
def test_prepare_400_when_body_empty(client):
    r = client.post("/api/receptors/extended/prepare", json={})
    assert r.status_code == 400
    assert "name" in r.json()["detail"] or "pdb" in r.json()["detail"]
def test_prepare_400_for_invalid_pdb_id_shape(client):
    """RCSB ids are 4 chars starting with 1-9; reject anything else
    before hitting the network."""
    r = client.post("/api/receptors/extended/prepare", json={"pdb": "XYZ"})
    assert r.status_code == 400
    r2 = client.post("/api/receptors/extended/prepare", json={"pdb": "0M17"})
    assert r2.status_code == 400
def test_prepare_404_for_unknown_name(client):
    r = client.post("/api/receptors/extended/prepare", json={"name": "Wibble"})
    assert r.status_code == 404
import pytest
@pytest.fixture()
def isolated_raw_pdb_cache(monkeypatch, tmp_path):
    """Redirect the raw_pdbs cache to a tmp dir so prepare tests don't
    litter backend/storage/raw_pdbs/. The route consults RAW_PDB_CACHE_DIR
    on every call so a single setenv is enough."""
    monkeypatch.setenv("RAW_PDB_CACHE_DIR", str(tmp_path / "raw_pdbs"))
    return tmp_path / "raw_pdbs"
async def test_prepare_downloads_then_short_circuits_when_cached(async_client, isolated_raw_pdb_cache, monkeypatch):
    """First call fetches RCSB and writes the file; second call sees
    the cache and skips the network entirely. Catch the absence of the
    second HTTP call by counting handler invocations.
    Uses async_client because the prepare endpoint creates an
    httpx.AsyncClient internally — httpx.MockTransport intercepts that cleanly on the
    test's own event loop, but the sync TestClient creates a fresh loop
    per request which trips a "client has been closed" race.
    """
    fake_pdb = (
        "HEADER    TEST                                    01-JAN-00   1ABC              \n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    call_count = {"n": 0}
    def _rcsb_handler(request):
        call_count["n"] += 1
        return Response(200, text=fake_pdb)
    mock = install(monkeypatch, "https://files.rcsb.org")
    mock.get("/download/1M17.pdb").mock(side_effect=_rcsb_handler)
    r1 = await async_client.post("/api/receptors/extended/prepare", json={"name": "EGFR"})
    assert r1.status_code == 200, r1.text
    d1 = r1.json()
    assert d1["pdb"] == "1M17"
    assert d1["cached"] is False
    assert d1["size_bytes"] == len(fake_pdb.encode("utf-8"))
    assert call_count["n"] == 1
    r2 = await async_client.post("/api/receptors/extended/prepare", json={"name": "EGFR"})
    assert r2.status_code == 200
    d2 = r2.json()
    # Second call must short-circuit before hitting RCSB.
    assert d2["cached"] is True
    assert call_count["n"] == 1, "RCSB was hit twice — cache did not short-circuit"
def test_prepare_502_when_rcsb_returns_non_pdb_content(client, isolated_raw_pdb_cache, monkeypatch):
    """RCSB sometimes returns an HTML 200 (e.g. CDN cache miss page).
    Guard against that landing as a corrupt PDB on disk."""
    mock = install(monkeypatch, "https://files.rcsb.org")
    mock.get("/download/1M17.pdb").mock(
        return_value=Response(200, text="<html>not a PDB</html>")
    )
    r = client.post("/api/receptors/extended/prepare", json={"name": "EGFR"})
    assert r.status_code == 502
    assert "non-PDB content" in r.json()["detail"]
def test_prepare_502_when_rcsb_404(client, isolated_raw_pdb_cache, monkeypatch):
    mock = install(monkeypatch, "https://files.rcsb.org")
    mock.get("/download/1M17.pdb").mock(
        return_value=Response(404, text="not found")
    )
    r = client.post("/api/receptors/extended/prepare", json={"name": "EGFR"})
    assert r.status_code == 502
    assert "404" in r.json()["detail"]
async def test_prepare_with_pdb_override_skips_catalog(async_client, isolated_raw_pdb_cache, monkeypatch):
    """Direct PDB id bypasses the catalog — useful when the user wants
    a different structure than the curated default. catalog_entry is
    None in the response so the UI can render "manual" instead of
    "{name} ({gene})"."""
    fake_pdb = (
        "HEADER    TEST\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
    )
    mock = install(monkeypatch, "https://files.rcsb.org")
    mock.get("/download/2RH1.pdb").mock(return_value=Response(200, text=fake_pdb))
    r = await async_client.post("/api/receptors/extended/prepare", json={"pdb": "2RH1"})
    assert r.status_code == 200
    d = r.json()
    assert d["pdb"] == "2RH1"
    assert d["catalog_entry"] is None
