"""Headers for backend-to-internal model services."""
from app.core.config import get_settings


def service_headers() -> dict[str,str]:
    cfg=get_settings()
    if cfg.app_env != "production":
        return {}
    if len(cfg.internal_service_key) < 32:
        raise RuntimeError("INTERNAL_SERVICE_KEY is required for production service calls")
    return {"X-Internal-Service-Key":cfg.internal_service_key}
