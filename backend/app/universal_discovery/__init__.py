"""In-process De Novo discovery (formerly separate yeni_motor HTTP service)."""

import os
from pathlib import Path

from pydantic import ValidationError

from .config import AppConfig, load_config
from .models import DiscoveryRequest
from .pipeline import run_discovery

_UNIVERSAL_CFG: AppConfig | None = None


def _ilac_master_repo_root() -> Path:
    """Parent of `backend/` (this file: backend/app/universal_discovery/__init__.py)."""
    return Path(__file__).resolve().parents[3]


def universal_discovery_config_path() -> Path:
    """Resolve TOML: env override → yeni_motor copy → embedded fallback."""
    override = (os.environ.get("UNIVERSAL_DISCOVERY_TOML") or "").strip()
    if override:
        p = Path(override).expanduser()
        if p.is_file():
            return p.resolve()
    external = _ilac_master_repo_root() / "yeni_motor" / "universal_ai_drug_discovery.toml"
    if external.is_file():
        return external.resolve()
    return (Path(__file__).resolve().parent / "universal_ai_drug_discovery.toml").resolve()


def get_universal_discovery_config() -> AppConfig:
    global _UNIVERSAL_CFG
    if _UNIVERSAL_CFG is None:
        path = universal_discovery_config_path()
        if not path.is_file():
            raise FileNotFoundError(f"Discovery config missing: {path}")
        _UNIVERSAL_CFG = load_config(path)
    return _UNIVERSAL_CFG


def run_universal_discovery(denovo_request: dict) -> dict:
    """Validate payload, run pipeline, return dict (same shape as former /api/discover JSON)."""
    req = DiscoveryRequest.model_validate(denovo_request)
    cfg = get_universal_discovery_config()
    return run_discovery(req, cfg).model_dump()


__all__ = [
    "DiscoveryRequest",
    "ValidationError",
    "get_universal_discovery_config",
    "run_universal_discovery",
    "run_discovery",
    "universal_discovery_config_path",
]
