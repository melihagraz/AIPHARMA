# ---------------------------------------------------------------------------
# This file is adapted from REINVENT 4 (https://github.com/MolecularAI/REINVENT4)
# Copyright (c) AstraZeneca AB
# Licensed under the Apache License, Version 2.0 (see LICENSE & NOTICE in
# the parent directory of this package).
#
# Source path in REINVENT 4 repo: reinvent/models/model_mode_enum.py
# Modifications: imports rewritten to local package paths; otherwise unchanged.
# ---------------------------------------------------------------------------

from dataclasses import dataclass


@dataclass(frozen=True)
class ModelModeEnum:
    INFERENCE = "inference"
    TRAINING = "training"
