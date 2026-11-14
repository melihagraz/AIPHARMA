#!/usr/bin/env python3
"""REAL local REINVENT -> ADMET-AI -> AutoDock Vina integration experiment.

No fake downstream scores, cached stand-ins, or model approximations. On a
missing runtime, writes an explicit blocked report rather than fabricating
scores. Default EGFR uses a bundled receptor; ALK5 has no bundled receptor.

Invoke from repository root after installing each service's real dependencies:
    python scripts/run_true_three_stage_e2e_20261109.py --target EGFR --count 4

This is *software integration*, not target activity / safety validation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata, util
import json
import os
from pathlib import Path
import platform
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package_version(name: str):
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def _local_service_headers() -> dict[str, str]:
    """Supply the same internal service credential a production caller needs.

    This script uses in-process FastAPI TestClient, not public network listeners.
    Local development/test requires no credential; production fails closed.
    """
    environment = os.getenv("APP_ENV", "development").strip().lower()
    if environment == "production":
        secret = os.getenv("INTERNAL_SERVICE_KEY", "")
        if len(secret) < 32:
            raise RuntimeError("APP_ENV=production needs INTERNAL_SERVICE_KEY >=32 chars")
        return {"X-Internal-Service-Key": secret}
    if environment not in {"development", "test"}:
        raise RuntimeError("APP_ENV must be development, test, or production")
    return {}


def record(report: dict[str, Any], out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str) + "\n")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", default="EGFR", help="Bundled receptor; ALK5 requires a separately validated receptor")
    p.add_argument("--count", type=int, default=4, help="Small fixed-seed generator batch (2..32)")
    p.add_argument("--seed", type=int, default=20261109)
    p.add_argument("--out", type=Path, default=ROOT / "reports/live_e2e_20261109/run.json")
    args = p.parse_args()
    if not (2 <= args.count <= 32):
        p.error("--count must be between 2 and 32")

    from rdkit import Chem
    from fastapi.testclient import TestClient
    from Docking_App.app.dock import _vina_binary, _resolve_target, list_targets

    report: dict[str, Any] = {
        "name": "AIPHARMA actual-model service integration test",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "not_biological_or_clinical_validation": True,
        "no_simulated_or_replayed_admet_or_docking_results": True,
        "full_three_stage_success": False,
        "target": args.target,
        "seed": args.seed,
        "requested_count": args.count,
        "environment": {"python": sys.version.split()[0], "platform": platform.platform(),
                        "torch": package_version("torch"), "rdkit": package_version("rdkit"),
                        "admet-ai": package_version("admet-ai"), "meeko": package_version("meeko"),
                        "vina": package_version("vina")},
        "stages": {}, "blockers": [],
    }
    # Independently probe every required executable, even if the first missing
    # dependency prevents the full run. Never conflate available with executed.
    available_admet = util.find_spec("admet_ai") is not None
    available_meeko = util.find_spec("meeko") is not None
    try:
        native_vina = _vina_binary()
    except FileNotFoundError:
        native_vina = None
    try:
        receptor_path, conf_path = _resolve_target(args.target)
    except ValueError:
        receptor_path = conf_path = None
    report["prerequisites"] = {
        "admet_ai_installed": available_admet,
        "meeko_installed": available_meeko,
        "native_vina_available": native_vina is not None,
        "bundled_receptor_and_grid_available": receptor_path is not None,
        "bundled_targets": list_targets(),
    }
    if not available_meeko or native_vina is None or receptor_path is None:
        report["stages"]["vina"] = {"status": "not_run_missing_prerequisites"}
        if not available_meeko:
            report["blockers"].append("Install Meeko for actual ligand PDBQT preparation")
        if native_vina is None:
            report["blockers"].append("Install a native AutoDock Vina executable (Windows .exe is unusable on Linux)")
        if receptor_path is None:
            report["blockers"].append(f"No validated receptor/grid for {args.target}; bundled targets: {list_targets()}")
    status = 2
    try:
        # Real local generator; independent of optional ADMET and docking installs.
        from Reinvent_App.app.main import app as reinvent_app
        prior = ROOT / "Reinvent_App/prior/reinvent.prior"
        report["model_prior_sha256"] = sha256(prior)
        payload = {"batch_size": args.count, "seed": args.seed,
                   "target_prior": "base", "validate_with_rdkit": True}
        t0 = time.perf_counter()
        with TestClient(reinvent_app) as client:
            response = client.post("/sample", json=payload, headers=_local_service_headers())
        if response.status_code != 200:
            raise RuntimeError(f"real REINVENT /sample failed with HTTP {response.status_code}: {response.text[:150]}")
        candidate_data = response.json()
        candidates = candidate_data["candidates"]
        smiles = [item["smiles"] for item in candidates if item.get("valid")]
        if not smiles or any(Chem.MolFromSmiles(s) is None for s in smiles):
            raise RuntimeError("no usable RDKit-valid molecules from actual REINVENT generator")
        report["stages"]["reinvent"] = {
            "status": "real_inference_passed", "runtime_s": round(time.perf_counter() - t0, 2),
            "returned": len(candidates), "valid": len(smiles), "smiles": smiles,
            "prior_used": candidate_data.get("prior_used")}

        # No cached/synthetic ADMET substitute is permitted.
        if not available_admet:
            report["stages"]["admet"] = {"status": "blocked_missing_admet_ai"}
            report["blockers"].append("Install admet-ai and its actual pretrained prediction models")
            return status

        from Admet_App.app.main import app as admet_app
        t0 = time.perf_counter()
        with TestClient(admet_app) as client:
            response = client.post("/predict", json={"smiles": smiles}, headers=_local_service_headers())
        if response.status_code != 200:
            raise RuntimeError(f"real ADMET /predict failed HTTP {response.status_code}: {response.text[:300]}")
        admet = response.json()
        preds = admet.get("predictions") or []
        if len(preds) != len(smiles) or [x.get("smiles") for x in preds] != smiles:
            raise RuntimeError("ADMET predictions are absent/out of molecule order")
        if any(not any(v is not None for v in x.get("raw", {}).values()) for x in preds):
            raise RuntimeError("ADMET produced one or more entirely empty endpoint rows")
        report["stages"]["admet"] = {
            "status": "real_model_inference_passed", "runtime_s": round(time.perf_counter() - t0, 2),
            "count": len(preds), "predictions": preds, "engine_claimed": admet.get("engine")}

        # Verify ligand prep AND native docking executable; .exe not usable on Linux.
        if util.find_spec("meeko") is None:
            report["stages"]["vina"] = {"status": "blocked_missing_meeko"}
            report["blockers"].append("Install Meeko in the Docking_App environment")
            return status
        try:
            vina_binary = _vina_binary()
        except FileNotFoundError as exc:
            report["stages"]["vina"] = {"status": "blocked_missing_native_vina"}
            report["blockers"].append(str(exc))
            return status
        try:
            receptor, conf = _resolve_target(args.target)
        except ValueError:
            report["stages"]["vina"] = {"status": "blocked_missing_bundled_receptor"}
            report["blockers"].append(f"No validated receptor/grid for {args.target}; bundled targets: {list_targets()}")
            return status

        from Docking_App.app.main import app as docking_app
        t0 = time.perf_counter()
        # Deliberately small subset: real CPU docking can be expensive.
        dock_smiles = smiles[:min(3, len(smiles))]
        with TestClient(docking_app) as client:
            response = client.post("/dock", json={"target": args.target, "smiles": dock_smiles,
                                                  "workers": 1, "cpu_per_dock": 1,
                                                  "exhaustiveness": 4, "num_modes": 3},
                                   headers=_local_service_headers())
        if response.status_code != 200:
            raise RuntimeError(f"real Vina /dock failed HTTP {response.status_code}: {response.text[:300]}")
        dock_result = response.json()
        rows = dock_result.get("results") or []
        if [x.get("smiles") for x in rows] != dock_smiles:
            raise RuntimeError("Vina result identity mismatch")
        successes = [x for x in rows if x.get("affinity") is not None]
        report["stages"]["vina"] = {
            "status": "real_docking_passed" if successes else "real_docking_failed_all",
            "runtime_s": round(time.perf_counter()-t0, 2), "requested": len(dock_smiles),
            "succeeded": len(successes), "results": rows,
            "binary_path": vina_binary, "receptor_sha256": sha256(receptor),
            "grid_sha256": sha256(conf)}
        if len(successes) != len(dock_smiles):
            report["blockers"].append("Vina did not produce valid affinity for all submitted molecules")
            return status
        report["full_three_stage_success"] = True
        status = 0
    except Exception as exc:
        report["stages"]["error"] = {"status": "failed", "kind": type(exc).__name__, "detail": str(exc)[:500]}
        report["blockers"].append(f"Execution failed: {type(exc).__name__}: {str(exc)[:160]}")
    finally:
        record(report, args.out)
        print(json.dumps({"full_three_stage_success": report["full_three_stage_success"],
                          "stages": {k: v["status"] for k, v in report["stages"].items()},
                          "blockers": report["blockers"], "report_path": str(args.out)},
                         ensure_ascii=False, indent=2))
    return status


if __name__ == "__main__":
    raise SystemExit(main())
