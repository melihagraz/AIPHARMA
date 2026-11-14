# ---------------------------------------------------------------------------
# This file is adapted from REINVENT 4 (https://github.com/MolecularAI/REINVENT4)
# Copyright (c) AstraZeneca AB
# Licensed under the Apache License, Version 2.0 (see LICENSE & NOTICE in
# the parent directory of this package).
#
# Source path in REINVENT 4 repo: reinvent/models/reinvent/models/dataset.py
# Modifications: imports rewritten to local package paths; otherwise unchanged.
# ---------------------------------------------------------------------------

"""SMILES dataset for PyTorch DataLoader"""

import torch
import torch.utils.data as tud

from . import conversions
from .utils import collate_fn


class Dataset(tud.Dataset):
    """Custom PyTorch Dataset that takes a file containing \n separated SMILES"""

    def __new__(cls, *args, **kwargs):
        if "randomize" in kwargs and kwargs["randomize"]:
            cls.__getitem__ = cls._getitem_with_randomization
        else:
            cls.__getitem__ = cls._getitem

        return super().__new__(cls)

    def __init__(self, smiles_list, vocabulary, tokenizer, randomize=False):
        self._vocabulary = vocabulary
        self._tokenizer = tokenizer
        self._smiles_list = list(smiles_list)

    def _getitem(self, i):
        smiles = self._smiles_list[i]
        tokens = self._tokenizer.tokenize(smiles)
        encoded = self._vocabulary.encode(tokens)

        return torch.tensor(encoded, dtype=torch.long)

    def _getitem_with_randomization(self, i):
        smiles = conversions.randomize_smiles(self._smiles_list[i])
        tokens = self._tokenizer.tokenize(smiles)
        encoded = self._vocabulary.encode(tokens)

        return torch.tensor(encoded, dtype=torch.long)

    def __len__(self):
        return len(self._smiles_list)


Dataset.collate_fn = collate_fn
