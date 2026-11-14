"""Pydantic schemas for the LLM-driven chat agent (Phase 2).

The agent is a chat-first orchestration layer over the existing production
discovery pipeline. The user types a free-text request ("EGFR'a karşı 200
molekül üret, top 5'i dock et"), Claude routes that to the right tool calls
(generate_molecules → predict_admet → dock_molecules → ...), then summarises.

Wire format (frontend ↔ /api/agent/chat):
    Request:  ChatRequest with message + history
    Response: SSE stream of ``ChatStreamEvent`` JSON, terminated by `event: end`
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Conversation history (mirrors Anthropic's messages-API role/content shape,
# but kept JSON-simple for transport — no SDK type leak into the wire format).
# ---------------------------------------------------------------------------

class ChatMessageIn(BaseModel):
    """One past message replayed by the client to reconstruct context."""

    role: Literal["user", "assistant"] = "user"
    content: str = Field(..., description="Plain text. Tool-call traces are not replayed.")


class ChatRequest(BaseModel):
    message: str = Field(..., min_length=1, description="Latest user message.")
    history: list[ChatMessageIn] = Field(
        default_factory=list,
        description="Prior turns (user + assistant text). Keep ≤ ~20 for cost.",
    )
    model: str | None = Field(
        default=None,
        description="Override the configured agent_model. Use sparingly.",
    )


# ---------------------------------------------------------------------------
# SSE stream event shapes — emitted by /api/agent/chat
# ---------------------------------------------------------------------------

class ChatTextDelta(BaseModel):
    """Incremental assistant text. Concatenate in order."""

    type: Literal["text"] = "text"
    text: str


class ChatThinkingDelta(BaseModel):
    """Adaptive-thinking summary chunk (display:summarized). May be empty."""

    type: Literal["thinking"] = "thinking"
    text: str


class ChatToolUse(BaseModel):
    """Claude is about to call a tool."""

    type: Literal["tool_use"] = "tool_use"
    name: str
    tool_use_id: str
    input: dict[str, Any]


class ChatToolResult(BaseModel):
    """Tool finished executing. ``content`` is the raw JSON string returned."""

    type: Literal["tool_result"] = "tool_result"
    tool_use_id: str
    name: str
    is_error: bool = False
    content: str  # JSON-encoded payload (truncated/summarised)


class ChatError(BaseModel):
    type: Literal["error"] = "error"
    message: str


class ChatDone(BaseModel):
    type: Literal["done"] = "done"
    stop_reason: str | None = None
    iterations: int = 0
    input_tokens: int | None = None
    output_tokens: int | None = None


ChatStreamEvent = (
    ChatTextDelta | ChatThinkingDelta | ChatToolUse | ChatToolResult | ChatError | ChatDone
)
