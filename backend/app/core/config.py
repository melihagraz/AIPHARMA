from functools import lru_cache
from pathlib import Path

from pydantic import AnyHttpUrl, field_validator, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

# Load repo-root .env when API is started from backend/ (uvicorn cwd).
_BACKEND_DIR = Path(__file__).resolve().parents[2]
_REPO_ROOT = _BACKEND_DIR.parent


def _env_file_tuple() -> tuple[str, ...]:
    paths = [_BACKEND_DIR / ".env", _REPO_ROOT / ".env"]
    found = [str(p) for p in paths if p.is_file()]
    return tuple(found) if found else (".env",)


class Settings(BaseSettings):
    app_name: str = "AI Drug Discovery Platform API"
    app_env: str = "development"
    api_prefix: str = "/api/v1"
    debug: bool = False
    cors_origins: list[AnyHttpUrl] | list[str] = [
        "http://localhost:3000",
        "http://localhost:5173",
        "http://localhost:5174",
        "http://localhost:5175",
        "http://127.0.0.1:5173",
        "http://127.0.0.1:5174",
        "http://127.0.0.1:5175",
    ]

    database_url: str = "sqlite:///./enterprise.db"
    redis_url: str = "redis://localhost:6379/0"
    enable_async_dispatch: bool = False

    docking_service_url: str = "http://localhost:8090"
    generator_service_url: str = "http://localhost:8091"
    universal_ai_service_url: str = "http://127.0.0.1:8010"
    #: Reinvent_App (REINVENT 4 SMILES generator) — ``cd Reinvent_App && uvicorn app.main:app --port 8020``
    reinvent_app_url: str = "http://127.0.0.1:8020"
    #: Admet_App (ADMET-AI predictor) — ``cd Admet_App && uvicorn app.main:app --port 8030``
    admet_app_url: str = "http://127.0.0.1:8030"
    #: Docking_App (AutoDock Vina) — ``cd Docking_App && uvicorn app.main:app --port 8040``
    docking_app_url: str = "http://127.0.0.1:8040"
    #: Pocket_App (Phase 3 SBDD: pocket detection + 3D pharmacophore extraction)
    pocket_app_url: str = "http://127.0.0.1:8050"
    #: embedded = ilac-master ``app.universal_discovery``; http = yeni_motor HTTP ``/api/discover``;
    #: http_fallback = embedded sonra HTTP; yeni_motor_subprocess = Node ``run_universal_discovery_stdin.cjs`` (universal_discovery_engine.cjs);
    #: production = Reinvent_App + Admet_App pipeline (gerçek motor)
    universal_discovery_mode: str = "production"
    #: De Novo: ``yeni_motor/scripts/run_universal_discovery_stdin.cjs`` varsa Node universal motor (``DENOVO_USE_YENI_MOTOR_SCRIPT=0`` ile kapat)
    denovo_use_yeni_motor_script: bool = True
    admin_username: str = "admin"
    admin_password: str = ""  # disabled until explicitly configured
    # Static admin bearer tokens are no longer accepted. Only signed, expiring sessions.
    admin_token: str = ""  # deprecated legacy setting (ignored)
    session_signing_key: str = ""  # must be explicitly set to permit login
    internal_service_key: str = ""  # shared secret for internal scientific microservices in production
    session_ttl_seconds: int = 3600

    # ------------------------------------------------------------------
    # Agent layer (Phase 2) — LLM-driven chat orchestration
    # ------------------------------------------------------------------
    #: Anthropic API key for the chat agent. Set via ``.env`` or env var
    #: ``ANTHROPIC_API_KEY``. Empty string disables the agent route at runtime.
    anthropic_api_key: str = ""
    #: Default model. Skill default — adaptive thinking + effort:"high".
    agent_model: str = "claude-opus-4-7"
    #: Hard ceiling per assistant turn (model-side). The agent loop also has
    #: its own max_iterations safety stop in agent_service.
    agent_max_tokens: int = 8192
    #: Max number of tool-call → tool-result roundtrips inside one user turn.
    agent_max_iterations: int = 12

    model_config = SettingsConfigDict(env_file=_env_file_tuple(), env_file_encoding="utf-8", extra="ignore")

    @model_validator(mode="after")
    def validate_release_secrets(self):
        if self.app_env == "production" and self.debug:
            raise ValueError("DEBUG must be disabled in production")
        if self.app_env == "production" and len(self.session_signing_key) < 32:
            raise ValueError("SESSION_SIGNING_KEY must be at least 32 characters in production")
        if self.app_env == "production" and len(self.internal_service_key) < 32:
            raise ValueError("INTERNAL_SERVICE_KEY must be at least 32 characters in production")
        return self

    @field_validator("app_env")
    @classmethod
    def normalize_env(cls, value: str) -> str:
        """Unknown environment labels must not silently disable production guards."""
        normalized = value.lower().strip()
        if normalized not in {"development", "test", "production"}:
            raise ValueError("APP_ENV must be development, test, or production")
        return normalized

    @field_validator("database_url", mode="before")
    @classmethod
    def resolve_sqlite_path(cls, value):
        """Keep the development SQLite database stable across launch modes.

        The API is started with ``backend/`` as its working directory by the
        orchestrator and Vite auto-starter, but some developers run uvicorn
        from the repository root with ``--app-dir backend``.  SQLite normally
        resolves relative URLs against the process working directory, which
        made those launch modes select different files.  Resolve a relative
        SQLite URL against ``backend/`` once so every entry point uses the
        same database and the repository remains movable between machines.
        """
        if not isinstance(value, str):
            return value
        prefix = "sqlite:///./"
        if not value.startswith(prefix):
            return value
        relative_path = value[len(prefix):]
        absolute_path = (_BACKEND_DIR / relative_path).resolve()
        return f"sqlite:///{absolute_path.as_posix()}"

    @field_validator("cors_origins", mode="before")
    @classmethod
    def parse_cors_origins(cls, value):
        if value is None:
            return value
        if isinstance(value, str):
            parts = [p.strip() for p in value.split(",") if p.strip()]
            return parts if parts else value
        return value

    @field_validator("anthropic_api_key", mode="before")
    @classmethod
    def _fallback_anthropic_key(cls, value):
        """Defend against an empty ``ANTHROPIC_API_KEY`` injected by the shell.

        Some launchers (Claude Code CLI, certain CI runners) export an empty
        ``ANTHROPIC_API_KEY`` into the parent shell. Pydantic-settings normally
        prefers ``os.environ`` over ``env_file``, so an empty shell var would
        silently override the real key in ``.env``. If the value reaching this
        validator is empty/whitespace, fall back to a direct ``.env`` read so
        the configured key always wins.
        """
        if isinstance(value, str) and value.strip():
            return value
        for env_file in _env_file_tuple():
            try:
                with open(env_file, encoding="utf-8") as fh:
                    for raw in fh:
                        line = raw.strip().lstrip("﻿")
                        if line.startswith("ANTHROPIC_API_KEY="):
                            v = line.split("=", 1)[1].strip().strip('"').strip("'")
                            if v:
                                return v
            except OSError:
                continue
        return value or ""


@lru_cache
def get_settings() -> Settings:
    return Settings()
