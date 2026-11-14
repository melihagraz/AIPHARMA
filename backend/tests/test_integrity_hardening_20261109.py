"""Offline integrity/security tests. No model endpoints or docking are simulated as science."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import httpx
import pandas as pd
import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Admet_App.app.predictor import AdmetPredictor
from Docking_App.app.dock import _parse_best_affinity
from Docking_App.app.schemas import DockRequest
from Docking_App.app import dock as docking_driver
from Pocket_App.app.main import app as pocket_app
from Pocket_App.app.security import validate_prepared_pdb_path
from app.services import pocket_client


@pytest.fixture
def production_pocket(monkeypatch):
    monkeypatch.setenv("APP_ENV", "production")
    monkeypatch.setenv("INTERNAL_SERVICE_KEY", "testing-secret-very-long-and-not-for-production")
    with TestClient(pocket_app) as client:
        yield client


def test_pocket_requires_internal_key_in_production(production_pocket):
    client = production_pocket
    assert client.get("/targets").status_code == 401
    assert client.get("/targets", headers={"X-Internal-Service-Key": "wrong"}).status_code == 401
    result = client.get(
        "/targets", headers={"X-Internal-Service-Key": "testing-secret-very-long-and-not-for-production"}
    )
    assert result.status_code == 200
    assert "EGFR" in result.json()["targets"]


def test_pocket_rejects_arbitrary_path_before_access(production_pocket, tmp_path):
    private = tmp_path / "private.pdb"
    private.write_text("ATOM random secret file\n")
    payload = {"pdb_id": "1M17", "prepared_pdb_path": str(private), "use_fpocket": False}
    r = production_pocket.post(
        "/detect", json=payload,
        headers={"X-Internal-Service-Key": "testing-secret-very-long-and-not-for-production"},
    )
    assert r.status_code == 400
    assert "approved PDB" in r.text


def test_pocket_approved_pdb_path_and_symlink_escape(tmp_path, monkeypatch):
    allowed = ROOT / "Pocket_App/structures/1py5.pdb"
    assert Path(validate_prepared_pdb_path(str(allowed))) == allowed.resolve()
    monkeypatch.setenv("POCKET_ALLOWED_PDB_DIR", str(tmp_path))
    target = tmp_path / "inside.pdb"
    target.write_text("ATOM test\n")
    assert Path(validate_prepared_pdb_path(str(target))) == target.resolve()
    outside = tmp_path.parent / "outside.pdb"
    outside.write_text("ATOM cannot access\n")
    link = tmp_path / "escape.pdb"
    link.symlink_to(outside)
    with pytest.raises(ValueError, match="approved PDB"):
        validate_prepared_pdb_path(str(link))
    outside.unlink()


def test_pocket_path_rejects_wrong_extension(tmp_path, monkeypatch):
    monkeypatch.setenv("POCKET_ALLOWED_PDB_DIR", str(tmp_path))
    fake = tmp_path / "secret.txt"
    fake.write_text("secret")
    with pytest.raises(ValueError, match="approved PDB"):
        validate_prepared_pdb_path(str(fake))


@pytest.mark.asyncio
async def test_pocket_http_client_forwards_internal_key(monkeypatch):
    seen = []
    def response(request):
        seen.append(request.headers.get("X-Internal-Service-Key"))
        return httpx.Response(200, json={"pocket": {"center_x": 0}})
    original = httpx.AsyncClient
    transport = httpx.MockTransport(response)
    monkeypatch.setattr(pocket_client, "service_headers", lambda: {"X-Internal-Service-Key": "expected"})
    monkeypatch.setattr(httpx, "AsyncClient", lambda *a, **k: original(*a, transport=transport, **k))
    assert (await pocket_client.detect_pocket(target="EGFR"))["center_x"] == 0
    assert seen == ["expected"]


def _predict_with_frame(frame, smiles):
    model = AdmetPredictor()
    model._model = SimpleNamespace(predict=lambda smiles: frame)
    return model.predict(smiles)


def test_admet_reorders_by_smiles_instead_of_mislabeling():
    frame = pd.DataFrame({"hERG": [0.8, 0.1], "AMES": [0.3, 0.2]}, index=pd.Index(["O", "CCO"], name="smiles"))
    result = _predict_with_frame(frame, ["CCO", "O"])
    assert [x.smiles for x in result] == ["CCO", "O"]
    assert [x.raw["hERG"] for x in result] == [0.1, 0.8]


def test_admet_duplicates_are_preserved():
    frame = pd.DataFrame({"smiles": ["O", "CCO", "CCO"], "hERG": [0.9, 0.1, 0.2]})
    result = _predict_with_frame(frame, ["CCO", "CCO", "O"])
    assert [x.raw["hERG"] for x in result] == [0.1, 0.2, 0.9]


def test_admet_rejects_missing_and_unknown_identities():
    bad = pd.DataFrame({"smiles": ["O", "N"], "hERG": [0.2, 0.3]})
    with pytest.raises(ValueError, match="SMILES identities"):
        _predict_with_frame(bad, ["O", "CCO"])
    with pytest.raises(ValueError, match="row count"):
        _predict_with_frame(bad, ["O"])


def test_admet_nonfinite_values_remain_null():
    frame = pd.DataFrame({"smiles": ["CCO"], "hERG": [float("inf")], "DILI": [float("nan")]})
    result = _predict_with_frame(frame, ["CCO"])
    assert result[0].raw == {"hERG": None, "DILI": None}


@pytest.mark.parametrize("fields", [
    {"center_x": 1.0},
    {"center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": 24, "size_y": 24},
    {"center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": -1, "size_y": 24, "size_z": 24},
    {"center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": 1000, "size_y": 24, "size_z": 24},
])
def test_docking_rejects_partial_or_unsafe_box(fields):
    with pytest.raises(ValidationError):
        DockRequest(target="EGFR", smiles=["CCO"], **fields)


def test_vina_no_affinity_is_explicit_failure(monkeypatch, tmp_path):
    receptor = tmp_path / "receptor.pdbqt"
    config = tmp_path / "conf.txt"
    receptor.write_text("RECEPTOR\n")
    config.write_text("center_x=0\n")
    monkeypatch.setattr(docking_driver, "_resolve_target", lambda target: (receptor, config))
    monkeypatch.setattr(docking_driver, "_vina_binary", lambda: "/usr/bin/true")
    monkeypatch.setattr(docking_driver, "smiles_to_pdbqt", lambda smiles: "MOCK PDBQT")
    monkeypatch.setattr(docking_driver.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess([], 0, "no modes", ""))
    result = docking_driver.dock_one("CCO", "EGFR")
    assert result.affinity is None
    assert result.error == "vina_output_missing_valid_affinity"


def test_vina_affinity_parser_accepts_actual_table():
    stdout = "mode | affinity | dist from best mode\n   1    -7.346     0.000     0.000\n"
    assert _parse_best_affinity(stdout) == -7.346

from app.services.production_pipeline import (
    _check_filters, validate_supported_filters, UnsupportedFilterError,
)
from app.services import admet_client, docking_client
from fastapi import HTTPException


def test_requested_model_filter_never_passes_missing_admet_prediction():
    failed, _ = _check_filters(
        {"mw": 250, "qed": 0.5}, {},
        {"safety_filters": {"herg_risk_max": 0.6},
         "pk_filters": {"oral_bioavailability_min": 0.3}},
        smiles="CCO",
    )
    assert "hERG_missing" in failed
    assert "Bioavailability_Ma_missing" in failed


def test_missing_requested_rdkit_descriptor_cannot_pass():
    failed, _ = _check_filters(
        {"mw": None, "qed": 0.55}, {"hERG": 0.1},
        {"base_filters": {"mw_max": 500}}, smiles="CCO"
    )
    assert "mw_missing" in failed


@pytest.mark.parametrize("bad_filters", [
    {"safety_filters": {"unimplemented_toxicity_max": 0.1}},
    {"base_filters": {"hbd_max": float("nan")}},
    {"docking": {"score_threshold": -7.5}},
    {"pk_filters": {"half_life_min_hours": 10}},
    {"unknown_group": {"foo": 1}},
])
def test_unknown_or_unsupported_scientific_filters_rejected(bad_filters):
    with pytest.raises(UnsupportedFilterError):
        validate_supported_filters(bad_filters)


def test_supported_ui_filters_accepted():
    assert validate_supported_filters({
        "base_filters": {"qed_min": .4, "mw_min": 150},
        "safety_filters": {"herg_risk_max": .6},
        "pk_filters": {"caco2_min": .4},
        "docking": {"enabled": True, "score_threshold": 0.0},
    }) is None


@pytest.mark.asyncio
async def test_admet_rejects_null_only_response(monkeypatch):
    original = httpx.AsyncClient
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json={
        "count": 1, "predictions": [{"smiles": "CCO", "raw": {"hERG": None}, "percentile": {}}]
    }))
    monkeypatch.setattr(httpx, "AsyncClient", lambda *a, **k: original(*a, transport=transport, **k))
    with pytest.raises(HTTPException, match="all unavailable"):
        await admet_client.predict_admet(["CCO"])


@pytest.mark.asyncio
async def test_docking_rejects_rogue_receptor_identity(monkeypatch):
    original = httpx.AsyncClient
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json={
        "target": "JAK2", "count": 1, "success_count": 1,
        "results": [{"smiles": "CCO", "target": "JAK2", "affinity": -7.0}]
    }))
    monkeypatch.setattr(httpx, "AsyncClient", lambda *a, **k: original(*a, transport=transport, **k))
    with pytest.raises(HTTPException, match="receptor identity"):
        await docking_client.dock_batch(["CCO"], target="EGFR")


@pytest.mark.asyncio
async def test_docking_rejects_missing_affinity_with_no_error(monkeypatch):
    original = httpx.AsyncClient
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json={
        "target": "EGFR", "count": 1, "success_count": 0,
        "results": [{"smiles": "CCO", "target": "EGFR", "affinity": None}]
    }))
    monkeypatch.setattr(httpx, "AsyncClient", lambda *a, **k: original(*a, transport=transport, **k))
    with pytest.raises(HTTPException, match="neither score nor error"):
        await docking_client.dock_batch(["CCO"], target="EGFR")


def test_parent_cluster_bootstrap_is_deterministic_and_archival():
    from scripts.alk5_cluster_bootstrap_20261109 import analyze
    first = analyze(n_boot=125, seed=91)
    second = analyze(n_boot=125, seed=91)
    assert first == second
    assert (first["n_designs"], first["n_distinct_parent_compounds"],
            first["n_comparisons"]) == (30, 16, 180)
    assert first["endpoints"]["DILI"]["improved_direction_count"] == 13
    assert first["endpoints"]["Caco2_Wang"]["improved_direction_count"] == 23
    assert first["not_biological_or_clinical_validation"] if "not_biological_or_clinical_validation" in first else bool(first["interpretation_limitations"])
