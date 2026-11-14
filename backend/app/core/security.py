"""Small, dependency-light signed sessions, with Argon2id credential hashing.

Tokens use HMAC-SHA256, carry an explicit expiration and role, and are NOT
OAuth/JWT. Use HTTPS, revocation and a managed identity provider for public
multi-user deployments. Keys must never be committed to source control.
"""
from __future__ import annotations

import base64
from dataclasses import dataclass
import hashlib
import hmac
import json
import secrets
import time

from argon2 import PasswordHasher
from argon2.exceptions import InvalidHashError, VerifyMismatchError, VerificationError
from fastapi import Header, HTTPException

from app.core.config import get_settings

_hasher = PasswordHasher(time_cost=3, memory_cost=65536, parallelism=2)


@dataclass(frozen=True)
class Principal:
    user_id: str
    role: str


def hash_password(plain: str) -> str:
    if not isinstance(plain, str) or len(plain) < 12 or len(plain) > 1024:
        raise ValueError("Password must contain 12 to 1024 characters")
    return _hasher.hash(plain)


def verify_password(stored: str, plain: str) -> bool:
    if not stored or not stored.startswith("$argon2id$"):
        return False
    try:
        return bool(_hasher.verify(stored, plain))
    except (InvalidHashError, VerifyMismatchError, VerificationError):
        return False


def _key() -> bytes:
    key = get_settings().session_signing_key
    if len(key) < 32:
        raise HTTPException(status_code=503, detail="Authentication is not configured; set SESSION_SIGNING_KEY")
    return key.encode("utf-8")


def issue_session_token(user_id: str, role: str = "user") -> str:
    if role not in {"user", "admin"} or not user_id:
        raise ValueError("Invalid session identity")
    config = get_settings()
    payload = {"sub": user_id, "role": role, "exp": int(time.time()) + config.session_ttl_seconds,
               "nonce": secrets.token_hex(8)}
    raw = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode()
    data = base64.urlsafe_b64encode(raw).rstrip(b"=").decode()
    mac = hmac.new(_key(), data.encode(), hashlib.sha256).hexdigest()
    return f"{data}.{mac}"


def decode_session_token(token: str) -> Principal:
    try:
        data, mac = token.rsplit(".", 1)
        expected = hmac.new(_key(), data.encode(), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(expected, mac):
            raise ValueError("Invalid signature")
        payload = json.loads(base64.urlsafe_b64decode(data + "=" * (-len(data) % 4)))
        if int(payload["exp"]) <= int(time.time()) or payload["role"] not in {"user", "admin"}:
            raise ValueError("Expired or invalid role")
        sub = payload["sub"]
        if not isinstance(sub, str) or not sub:
            raise ValueError("Invalid user")
        return Principal(user_id=sub, role=payload["role"])
    except HTTPException:
        raise
    except (ValueError, TypeError, KeyError, UnicodeError, json.JSONDecodeError) as exc:
        raise HTTPException(status_code=401, detail="Invalid or expired session token") from exc


def require_session(authorization: str | None = Header(default=None)) -> Principal:
    if not authorization or not authorization.startswith("Bearer "):
        raise HTTPException(status_code=401, detail="Bearer token required")
    return decode_session_token(authorization[7:].strip())


def require_admin(authorization: str | None = Header(default=None)) -> Principal:
    principal = require_session(authorization)
    if principal.role != "admin":
        raise HTTPException(status_code=403, detail="Administrator permissions required")
    return principal
