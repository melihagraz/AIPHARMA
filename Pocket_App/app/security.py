"""Guard internal pocket endpoints and prevent arbitrary filesystem reads.

Only server-managed PDB directories are accepted for prepared PDB paths.
The browser never receives the internal service secret.
"""
from __future__ import annotations

import hmac
import os
from pathlib import Path

from fastapi import Request
from fastapi.responses import JSONResponse

STRUCTURES_DIR = Path(__file__).resolve().parents[1] / "structures"
MAX_PDB_BYTES = 20 * 1024 * 1024


def check_internal_access(request: Request):
    environment = os.getenv("APP_ENV", "development").strip().lower()
    if environment not in {"development", "test", "production"}:
        return JSONResponse(status_code=503, content={"detail": "Invalid APP_ENV configuration"})
    if environment != "production":
        return None
    secret = os.getenv("INTERNAL_SERVICE_KEY", "")
    if len(secret) < 32:
        return JSONResponse(status_code=503, content={"detail": "Internal service key not configured"})
    provided = request.headers.get("X-Internal-Service-Key", "")
    if not hmac.compare_digest(provided, secret):
        return JSONResponse(status_code=401, content={"detail": "Internal service authentication required"})
    return None


def validate_prepared_pdb_path(raw_path: str | None) -> str | None:
    """Reject paths outside configured PDB storage, symlink escapes or oversized files.

    An additional private storage directory may be set via POCKET_ALLOWED_PDB_DIR.
    Arbitrary user-provided paths are not treated as trusted local resources.
    """
    if raw_path is None:
        return None
    raw = Path(raw_path)
    allowed_roots = [STRUCTURES_DIR.resolve()]
    configured = os.getenv("POCKET_ALLOWED_PDB_DIR", "").strip()
    if configured:
        allowed_roots.append(Path(configured).resolve())
    try:
        candidate = raw.resolve(strict=True)
        if not candidate.is_file() or candidate.suffix.lower() not in {".pdb", ".ent"}:
            raise ValueError("not a regular .pdb/.ent file")
        if not any(candidate.is_relative_to(base) for base in allowed_roots):
            raise ValueError("path is outside approved PDB storage")
        if candidate.stat().st_size > MAX_PDB_BYTES:
            raise ValueError("PDB exceeds 20 MiB limit")
    except (OSError, RuntimeError, ValueError) as exc:
        raise ValueError(f"prepared_pdb_path is not an approved PDB file: {exc}") from None
    return str(candidate)
