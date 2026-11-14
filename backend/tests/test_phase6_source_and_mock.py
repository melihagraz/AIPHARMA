"""Phase 6: strict network stubbing and document-level fixed-score sensitivity."""
import importlib.util
from pathlib import Path

import httpx
import pytest

from tests.httpx_stub import install

_ROOT = Path(__file__).resolve().parents[2]
_SPEC = importlib.util.spec_from_file_location("phase6_source", _ROOT / "scripts/alk5_document_sensitivity_phase6.py")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


def test_real_archive_source_impact_and_provenance():
    import json
    data = json.loads((_ROOT / "reports/alk5_endpoint_archive_20261016/dataset.json").read_text())
    preds = json.loads((_ROOT / "reports/alk5_endpoint_archive_20261016/predictions.json").read_text())
    result = _MODULE.evaluate_source_deletions(data, preds)
    assert result["dataset_n"] == 356
    assert result["negative_n"] == 11
    assert result["negative_containing_documents"] == 7
    lm = result["splits"]["document_and_scaffold"]["logistic"]
    assert lm["original_roc_auc"] == pytest.approx(0.9175230566534915)
    assert lm["deletion_auc_min"] == pytest.approx(0.8746355685131195)
    assert any(x["excluded_document"] == "CHEMBL3352089" and x["removed_compounds"] == 7
               for x in lm["deletions"])


def test_source_deletion_respects_multi_document_provenance():
    # A multi-source positive compound is excluded if either of its
    # cited documents is deleted; compound identity is never duplicated.
    ds = [{"smiles": "A", "label": 1, "documents": ["n1", "p"]},
          {"smiles": "B", "label": 1, "documents": ["p"]},
          {"smiles": "C", "label": 1, "documents": ["p"]},
          {"smiles": "D", "label": 0, "documents": ["n1"]},
          {"smiles": "E", "label": 0, "documents": ["n2"]}]
    pred = [{"smiles": r["smiles"], "label": r["label"], "score": s,
             "split": "source", "method": "toy"}
            for r,s in zip(ds, [0.8, 0.7, 0.6, 0.2, 0.1])]
    report = _MODULE.evaluate_source_deletions(ds, pred)
    d = report["splits"]["source"]["toy"]["deletions"]
    removed = {r["excluded_document"]:r["removed_compounds"] for r in d}
    assert removed == {"n1": 2, "n2": 1}


def test_source_deletion_rejects_mismatched_labels():
    ds = [{"smiles":"A", "label":1,"documents":["a"]},
          {"smiles":"B", "label":0,"documents":["b"]}]
    pred = [{"smiles":"A", "label":0,"score":0.9,"split":"s","method":"m"},
            {"smiles":"B", "label":0,"score":0.1,"split":"s","method":"m"}]
    with pytest.raises(ValueError,match="Label mismatch"):
        _MODULE.evaluate_source_deletions(ds,pred)


@pytest.mark.asyncio
async def test_http_mock_rejects_undocumented_external_host(monkeypatch):
    router=install(monkeypatch,"https://files.rcsb.org")
    router.get("/download/1M17.pdb").mock(return_value=httpx.Response(200,text="HEADER test"))
    async with httpx.AsyncClient() as client:
        with pytest.raises(AssertionError,match="Unexpected network host"):
            await client.get("https://example.com/secrets")


@pytest.mark.asyncio
async def test_http_mock_never_returns_unregistered_request(monkeypatch):
    router=install(monkeypatch,"http://127.0.0.1:8040")
    router.get("/health").mock(return_value=httpx.Response(200,json={"ok":True}))
    async with httpx.AsyncClient() as client:
        with pytest.raises(AssertionError, match="Unexpected external HTTP request"):
            await client.post("http://127.0.0.1:8040/dock", json={"smiles": ["CCO"]})


@pytest.mark.asyncio
@pytest.mark.parametrize("path", ["/api/v1/jobs-export", "/api/v1/admin-debug", "/api/agent/healthcheck", "/api/v1/jobs/start-denovo"])
async def test_production_guard_requires_exact_route_segments(monkeypatch, path):
    from fastapi import Request
    from fastapi.responses import JSONResponse
    from app import main
    monkeypatch.setattr(main.settings, "app_env", "production")
    scope = {"type": "http", "path": path, "raw_path": path.encode(),
             "headers": [], "query_string": b"", "method": "GET", "scheme": "http", "server": ("test",80)}
    called = []
    async def nxt(_request):
        called.append(True)
        return JSONResponse({"allowed": True})
    r = await main.block_unreviewed_legacy_api_in_production(Request(scope), nxt)
    assert r.status_code == 403
    assert not called


@pytest.mark.asyncio
@pytest.mark.parametrize("path", ["/api/v1/jobs", "/api/v1/jobs/11111111-1111-1111-1111-111111111111", "/api/v1/admin/users", "/api/auth/login", "/api/agent/health"])
async def test_production_guard_preserves_intentionally_allowed_routes(monkeypatch, path):
    from fastapi import Request
    from fastapi.responses import JSONResponse
    from app import main
    monkeypatch.setattr(main.settings, "app_env", "production")
    scope = {"type": "http", "path": path, "raw_path": path.encode(),
             "headers": [], "query_string": b"", "method": "GET", "scheme": "http", "server": ("test",80)}
    called = []
    async def nxt(_request):
        called.append(True)
        return JSONResponse({"allowed": True})
    r = await main.block_unreviewed_legacy_api_in_production(Request(scope), nxt)
    assert r.status_code == 200
    assert called
