"""De Novo: ``yeni_motor/scripts/run_universal_discovery_stdin.cjs`` → ``universal_discovery_engine.cjs`` + config JSON."""

from __future__ import annotations

import asyncio
import json
import logging
import os
import shutil
from pathlib import Path

from fastapi import HTTPException

logger = logging.getLogger(__name__)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _node_bin() -> str:
    return (os.environ.get("NODE") or "").strip() or shutil.which("node") or "node"


async def run_discovery_via_yeni_motor_script(denovo_request: dict) -> dict:
    """Eski isim korunur; artık Node universal motor çalışır (Python ``run_discover_client`` yok)."""
    repo = _repo_root()
    script = repo / "yeni_motor" / "scripts" / "run_universal_discovery_stdin.cjs"
    engine = repo / "yeni_motor" / "scripts" / "universal_discovery_engine.cjs"
    cfg = repo / "yeni_motor" / "scripts" / "universal_discovery_config.json"
    if not script.is_file():
        raise HTTPException(
            status_code=500,
            detail=f"Node discovery stdin runner not found: {script}",
        )
    if not engine.is_file() or not cfg.is_file():
        raise HTTPException(
            status_code=500,
            detail=f"universal_discovery_engine.cjs or universal_discovery_config.json missing under yeni_motor/scripts",
        )
    yeni_root = repo / "yeni_motor"
    payload = json.dumps(denovo_request, ensure_ascii=False).encode("utf-8")
    node = _node_bin()
    proc = await asyncio.create_subprocess_exec(
        node,
        str(script),
        stdin=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=str(yeni_root),
    )
    try:
        out, err = await asyncio.wait_for(proc.communicate(input=payload), timeout=900.0)
    except asyncio.TimeoutError:
        proc.kill()
        raise HTTPException(status_code=504, detail="Node discovery subprocess timed out (900s).") from None

    if proc.returncode != 0:
        msg = (err or b"").decode("utf-8", errors="replace").strip() or f"exit {proc.returncode}"
        logger.error("run_universal_discovery_stdin.cjs failed: %s", msg[:4000])
        raise HTTPException(status_code=502, detail=f"Node universal discovery failed: {msg[:2000]}")

    text = (out or b"").decode("utf-8", errors="replace").strip()
    if not text:
        raise HTTPException(status_code=502, detail="Node discovery returned empty stdout.")
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        logger.error("Invalid JSON from Node discovery (first 500 chars): %s", text[:500])
        raise HTTPException(
            status_code=502,
            detail=f"Node discovery stdout is not JSON: {exc}",
        ) from exc
