"""Chat agent route — SSE streaming over the production pipeline.

POST /api/agent/chat
    body: ChatRequest (message + history + optional model override)
    response: text/event-stream — newline-delimited ``data: {json}\\n\\n`` frames
              terminated by ``event: end\\n\\n``

GET /api/agent/health
    quick liveness probe + reports whether the API key is configured.

The route is intentionally thin: it only adapts ``stream_chat`` (an async
generator of Pydantic events) to the Server-Sent Events transport. All
orchestration logic lives in ``app.services.agent_service``.
"""

from __future__ import annotations

import json
import logging
from typing import AsyncIterator

from fastapi import APIRouter, Depends
from app.core.security import require_admin
from fastapi.responses import StreamingResponse

from app.core.config import get_settings
from app.schemas.agent import ChatRequest
from app.services.agent_service import TOOLS, stream_chat

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/agent", tags=["agent"])


@router.get("/health")
async def agent_health() -> dict[str, object]:
    cfg = get_settings()
    return {
        "status": "ok",
        "api_key_configured": bool(cfg.anthropic_api_key),
        "model": cfg.agent_model,
        "max_iterations": cfg.agent_max_iterations,
        "tool_count": len(TOOLS),
    }


@router.get("/tools")
async def list_tools() -> dict[str, object]:
    """Expose the tool catalogue to the UI (for debugging / sidebar display)."""
    return {
        "count": len(TOOLS),
        "tools": [{"name": t["name"], "description": t["description"]} for t in TOOLS],
    }


def _sse(event: str, data: object) -> bytes:
    """Encode one Server-Sent Events frame (UTF-8 bytes, \\n\\n terminator)."""
    payload = json.dumps(data, ensure_ascii=False, default=str)
    # ``event:`` is optional but helpful for client debugging.
    return f"event: {event}\ndata: {payload}\n\n".encode("utf-8")


async def _sse_stream(req: ChatRequest) -> AsyncIterator[bytes]:
    yield _sse("start", {"ok": True})
    try:
        async for ev in stream_chat(req):
            # Pydantic v2 ``model_dump`` keeps the ``type`` discriminator field.
            payload = ev.model_dump()
            event_name = payload.get("type", "message")
            yield _sse(event_name, payload)
    except Exception as exc:  # noqa: BLE001
        logger.exception("agent SSE stream failed")
        yield _sse("error", {"type": "error", "message": f"stream failed: {exc}"})
    finally:
        yield _sse("end", {"ok": True})


@router.post("/chat")
async def chat(req: ChatRequest, _principal=Depends(require_admin)) -> StreamingResponse:
    return StreamingResponse(
        _sse_stream(req),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache, no-transform",
            "Connection": "keep-alive",
            # Disable buffering on nginx / hosted reverse proxies.
            "X-Accel-Buffering": "no",
        },
    )
