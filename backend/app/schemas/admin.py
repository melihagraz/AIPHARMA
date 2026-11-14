from datetime import datetime
from pydantic import BaseModel, Field


class AdminLoginRequest(BaseModel):
    username: str
    password: str


class AdminLoginResponse(BaseModel):
    access_token: str
    token_type: str = "bearer"


class AdminUserCreateRequest(BaseModel):
    email: str
    full_name: str = ""
    is_admin: bool = False
    pocket_quota: int = Field(default=10, ge=0)
    denovo_quota: int = Field(default=10, ge=0)
    admet_quota: int = Field(default=10, ge=0)
    docking_quota: int = Field(default=10, ge=0)


class AdminUserUpdateRequest(BaseModel):
    full_name: str | None = None
    is_active: bool | None = None
    is_admin: bool | None = None
    pocket_quota: int | None = Field(default=None, ge=0)
    denovo_quota: int | None = Field(default=None, ge=0)
    admet_quota: int | None = Field(default=None, ge=0)
    docking_quota: int | None = Field(default=None, ge=0)


class AdminUserResponse(BaseModel):
    id: str
    email: str
    full_name: str
    is_active: bool
    is_admin: bool
    pocket_quota: int
    denovo_quota: int
    admet_quota: int
    docking_quota: int
    created_at: datetime
    updated_at: datetime


class AdminAuditLogResponse(BaseModel):
    id: str
    actor_id: str
    action: str
    resource_type: str
    resource_id: str
    payload: str
    created_at: datetime
