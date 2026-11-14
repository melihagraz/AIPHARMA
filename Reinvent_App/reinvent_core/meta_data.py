# ---------------------------------------------------------------------------
# This file is adapted from REINVENT 4 (https://github.com/MolecularAI/REINVENT4)
# Copyright (c) AstraZeneca AB
# Licensed under the Apache License, Version 2.0 (see LICENSE & NOTICE in
# the parent directory of this package).
#
# Source path in REINVENT 4 repo: reinvent/models/meta_data.py
# Modifications: stripped to the dataclass needed for pickle compatibility;
# hash/update helpers (update_model_data, check_valid_hash) removed because
# they are not required for inference and brought in the xxhash dependency.
# ---------------------------------------------------------------------------

from __future__ import annotations

__all__ = ["ModelMetaData"]

from dataclasses import dataclass, field, asdict
from uuid import UUID
from typing import List, Tuple


@dataclass
class ModelMetaData:
    """Meta data for model files to support data provenance and lineage."""

    hash_id: str | None
    hash_id_format: str
    model_id: UUID | str
    origina_data_source: str  # original typo preserved for pickle compat
    creation_date: float
    date_format: str = "UNIX epoch"
    model_id_format: str = ""
    updates: List[Tuple] = field(default_factory=list)
    comments: List[str] = field(default_factory=list)

    def as_dict(self):
        return asdict(self)
