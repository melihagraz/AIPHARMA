# ---------------------------------------------------------------------------
# This file is adapted from REINVENT 4 (https://github.com/MolecularAI/REINVENT4)
# Copyright (c) AstraZeneca AB
# Licensed under the Apache License, Version 2.0 (see LICENSE & NOTICE in
# the parent directory of this package).
#
# Source path in REINVENT 4 repo: reinvent/models/reinvent/utils.py
# Modifications: imports rewritten to local package paths; otherwise unchanged.
# ---------------------------------------------------------------------------

"""Helper routines

A set of common auxiliary functionality.
"""

from typing import List

import torch


def collate_fn(encoded_seqs: List[torch.Tensor]) -> torch.Tensor:
    """Converts a list of encoded sequences into a padded tensor

    :param encoded_seqs: encodes sequences to be padded with zeroes
    :return: padded tensor
    """

    max_length = max([seq.size(0) for seq in encoded_seqs])
    collated_arr = torch.zeros(
        len(encoded_seqs), max_length, dtype=torch.long
    )  # padded with zeroes

    for i, seq in enumerate(encoded_seqs):
        collated_arr[i, : seq.size(0)] = seq

    return collated_arr
