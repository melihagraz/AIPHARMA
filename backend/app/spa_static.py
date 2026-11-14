"""Tek süreçte UI: Vite build ``dist/`` + FastAPI ``/api`` aynı kökten sunulur."""

from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles


def mount_spa(app: FastAPI) -> None:
    """``dist/index.html`` varsa statik SPA ve ``/`` fallback ekler. ``SERVE_SPA=0`` ile kapatılır."""
    repo_root = Path(__file__).resolve().parents[2]
    dist = repo_root / "dist"
    index = dist / "index.html"
    if not index.is_file():
        return

    flag = os.environ.get("SERVE_SPA", "").strip().lower()
    if flag in ("0", "false", "no", "off"):
        return

    assets = dist / "assets"
    if assets.is_dir():
        app.mount("/assets", StaticFiles(directory=str(assets)), name="spa_assets")

    dist_resolved = dist.resolve()

    @app.get("/", include_in_schema=False)
    async def spa_index() -> FileResponse:
        return FileResponse(index)

    @app.get("/{full_path:path}", include_in_schema=False)
    async def spa_files(full_path: str) -> FileResponse:
        reserved = frozenset({"api", "docs", "redoc", "openapi.json", "assets"})
        seg0 = full_path.split("/", 1)[0] if full_path else ""
        if seg0 in reserved:
            raise HTTPException(status_code=404, detail="Not Found")
        target = (dist / full_path).resolve()
        try:
            target.relative_to(dist_resolved)
        except ValueError:
            raise HTTPException(status_code=404, detail="Not Found") from None
        if target.is_file():
            return FileResponse(target)
        return FileResponse(index)
