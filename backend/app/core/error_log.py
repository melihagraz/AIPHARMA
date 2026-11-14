"""Structured error log persisted as JSON-Lines at ``storage/logs/errors.jsonl``.

Why a JSONL file instead of just relying on the stdlib ``logging`` module:
    * Admin needs to read past errors over HTTP without tailing stdout.
    * Each entry stays parseable forever — easy to filter by ``category``,
      build dashboards, or feed back into the model.
    * Append-only with automatic rotation (one rolled-over file kept) so
      the file never grows unbounded.

Categories used so far:
    * ``rcsb_validation_failed`` — both RCSB and PDBe rejected/unreachable
    * ``pdb_fetch_failed``       — could not download a PDB from any mirror
    * ``protein_prep_failed``    — PDBFixer/OpenMM pipeline crashed
    * ``prolif_failed``          — ProLIF interaction fingerprint crashed
    * ``rcsb_search_failed``     — state-aware PDB suggestion failed
"""

from __future__ import annotations

import json
import logging
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

logger = logging.getLogger(__name__)

_LOG_DIR = Path(__file__).resolve().parents[2] / "storage" / "logs"
_LOG_FILE = _LOG_DIR / "errors.jsonl"
_ROTATED_FILE = _LOG_DIR / "errors.1.jsonl"
_MAX_BYTES = 5 * 1024 * 1024  # 5 MB before rotation

_LOCK = threading.Lock()


def log_error(category: str, message: str, **fields: Any) -> None:
    """Append a structured error event. Safe to call from any thread.

    Failures inside this function never propagate — error logging must not
    itself raise. We log to stdlib ``logging`` as a last resort.
    """
    entry = {
        "ts": time.time(),
        "ts_iso": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "category": category,
        "message": message,
    }
    if fields:
        entry["fields"] = {k: v for k, v in fields.items() if v is not None}

    with _LOCK:
        try:
            _LOG_DIR.mkdir(parents=True, exist_ok=True)
            if _LOG_FILE.exists() and _LOG_FILE.stat().st_size > _MAX_BYTES:
                if _ROTATED_FILE.exists():
                    _ROTATED_FILE.unlink()
                _LOG_FILE.rename(_ROTATED_FILE)
            with _LOG_FILE.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps(entry, default=str, ensure_ascii=False) + "\n")
        except Exception:  # noqa: BLE001
            logger.exception("Failed to write error log entry for category=%s", category)


def read_recent(limit: int = 50, category: Optional[str] = None) -> list[dict]:
    """Return the most-recent ``limit`` entries (newest first).

    Reads both the current log file and the most recent rotated file so
    rotation doesn't blind the admin view. Malformed lines are skipped.
    """
    entries: list[dict] = []
    for path in (_ROTATED_FILE, _LOG_FILE):
        if not path.exists():
            continue
        try:
            with path.open("r", encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        entry = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if category and entry.get("category") != category:
                        continue
                    entries.append(entry)
        except Exception:  # noqa: BLE001
            logger.exception("Failed to read error log file %s", path)
    return entries[-limit:][::-1]


def count_by_category() -> dict[str, int]:
    """Aggregate by category across both log files."""
    counts: dict[str, int] = {}
    for path in (_ROTATED_FILE, _LOG_FILE):
        if not path.exists():
            continue
        try:
            with path.open("r", encoding="utf-8") as fh:
                for line in fh:
                    try:
                        entry = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    cat = entry.get("category", "unknown")
                    counts[cat] = counts.get(cat, 0) + 1
        except Exception:  # noqa: BLE001
            continue
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))
