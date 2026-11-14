"""Process-local authentication throttling for research deployments.

Not a replacement for Redis/load-balancer-level rate limiting. Keys use the
actual peer address, never untrusted X-Forwarded-For. Separate buckets per peer
and peer/account restrict password guessing and credential stuffing.
"""
from __future__ import annotations

from collections import deque
from threading import RLock
import time

from fastapi import HTTPException, Request

_WINDOW_SECONDS = 300
_LIMIT = 8
_BUCKETS: dict[str, deque[float]] = {}
_LOCK = RLock()


def _keys(request: Request, account: str) -> tuple[str, str]:
    ip = request.client.host if request.client else "unknown-peer"
    return f"ip:{ip}", f"account:{ip}:{account.strip().lower()}"


def _prune(bucket: deque[float], now: float) -> None:
    while bucket and now - bucket[0] >= _WINDOW_SECONDS:
        bucket.popleft()


def check_login_throttle(request: Request, account: str) -> None:
    now = time.monotonic()
    with _LOCK:
        for key in _keys(request, account):
            bucket = _BUCKETS.get(key)
            if bucket:
                _prune(bucket, now)
                if len(bucket) >= _LIMIT:
                    raise HTTPException(status_code=429, detail="Too many login attempts. Retry later.")


def record_login_failure(request: Request, account: str) -> None:
    now = time.monotonic()
    with _LOCK:
        for key in _keys(request, account):
            bucket = _BUCKETS.setdefault(key, deque())
            _prune(bucket, now)
            bucket.append(now)
        if len(_BUCKETS) > 10000:
            # Bound memory under username-spray attempts; evict oldest
            # inactive buckets, retaining live rate-limited buckets.
            for key in list(_BUCKETS):
                _prune(_BUCKETS[key], now)
                if not _BUCKETS[key]:
                    _BUCKETS.pop(key, None)
