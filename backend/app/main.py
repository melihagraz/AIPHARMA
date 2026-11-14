from contextlib import asynccontextmanager

from fastapi import Depends, FastAPI, Request
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.orm import Session

from app.api.routes_health import router as health_router
from app.api.routes_jobs import router as jobs_router
from app.api.routes_admin import router as admin_router
from app.api.routes_legacy import router as legacy_router
from app.api.routes_agent import router as agent_router
from app.api.routes_pocket import router as pocket_router
from app.api.routes_alk5_design import router as alk5_design_router
from app.core.config import get_settings
from app.spa_static import mount_spa
from app.core.logging import configure_logging
from app.db.base import Base
from app.db.session import engine, get_db
from app import models  # noqa: F401
from app.universal_discovery.models import DiscoveryRequest, DiscoveryResponse

settings = get_settings()


@asynccontextmanager
async def lifespan(_: FastAPI):
    configure_logging()
    Base.metadata.create_all(bind=engine)
    yield


app = FastAPI(title=settings.app_name, debug=settings.debug, lifespan=lifespan)
@app.middleware("http")
async def block_unreviewed_legacy_api_in_production(request: Request, call_next):
    """Fail closed for legacy scientific routes that lack full tenant security.

    This intentionally limits public deployment functionality. Legacy job,
    download and scientific routes must be explicitly audited before enabling.
    Local controlled research mode (development) retains these features.
    """
    path = request.url.path
    if settings.app_env == "production" and path.startswith("/api/"):
        # Match exact path segments, not a loose startswith check. A future
        # /api/v1/jobs-export or /api/agent/healthcheck must not inherit
        # access merely because its name begins with an allowed prefix.
        scoped = (settings.api_prefix + "/jobs", settings.api_prefix + "/admin")
        safe = (
            path.startswith("/api/auth/")
            or path == "/api/agent/health"
            or any(path == base or path.startswith(base + "/") for base in scoped)
        )
        if not safe or path == settings.api_prefix + "/jobs/start-denovo":
            return JSONResponse(status_code=403, content={
                "detail": "Legacy research API is disabled for production until access controls are audited"
            })
    return await call_next(request)


app.add_middleware(
    CORSMiddleware,
    allow_origins=[str(origin) for origin in settings.cors_origins],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.post(f"{settings.api_prefix}/jobs/start-denovo")
@app.post("/api/jobs/start-denovo")
async def denovo_start_embedded(payload: dict, db: Session = Depends(get_db)):
    """Registered on the app before include_router so path /jobs/start-denovo is never matched as /jobs/{{job_id}} (405)."""
    from app.api.routes_legacy import legacy_start_denovo

    return await legacy_start_denovo(payload, db)


@app.post("/api/motor/discover", response_model=DiscoveryResponse)
async def motor_discover_embedded(body: DiscoveryRequest) -> DiscoveryResponse:
    """API_KULLANIM.txt — ilac-master köprüsü; De Novo da aynı ``motor_discover_impl`` yolunu kullanır."""
    from app.api.routes_legacy import motor_discover_impl

    return await motor_discover_impl(body)


app.include_router(health_router, prefix=settings.api_prefix)
app.include_router(jobs_router, prefix=settings.api_prefix)
app.include_router(admin_router, prefix=settings.api_prefix)
app.include_router(legacy_router)
# Agent layer (Phase 2) — chat-driven LLM orchestration over the pipeline.
# The router declares its own ``/api/agent`` prefix so it's mounted as-is.
app.include_router(agent_router)
# SBDD layer (Phase 3) — synchronous pocket detection + pharmacophore.
app.include_router(pocket_router)
app.include_router(alk5_design_router)

mount_spa(app)
