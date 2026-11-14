from datetime import datetime
from sqlalchemy import DateTime, Integer, String, Text
from sqlalchemy.orm import Mapped, mapped_column

from app.db.base import Base


class Job(Base):
    __tablename__ = "jobs"

    id: Mapped[str] = mapped_column(String(36), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(String(64), index=True, default="default")
    project_id: Mapped[str] = mapped_column(String(64), index=True, default="default")
    status: Mapped[str] = mapped_column(String(32), index=True)
    mode: Mapped[str] = mapped_column(String(32), index=True)
    receptor_name: Mapped[str | None] = mapped_column(String(255), nullable=True)
    pdb_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    uniprot_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    fasta_sequence: Mapped[str | None] = mapped_column(Text, nullable=True)
    candidate_count: Mapped[int] = mapped_column(Integer, default=100)
    # Progress telemetry — populated by the production pipeline's progress
    # reporter while a long-running De Novo job is in flight. Lets the UI
    # show a real progress bar instead of a 12-minute spinner.
    progress_pct: Mapped[int] = mapped_column(Integer, default=0)
    current_step: Mapped[str | None] = mapped_column(String(64), nullable=True)
    progress_message: Mapped[str | None] = mapped_column(Text, nullable=True)
    pipeline_resolution: Mapped[str | None] = mapped_column(String(32), nullable=True)
    error_message: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
