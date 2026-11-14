"""REINVENT 4 prior loader + sampler.

This module is original AI Pharma code; it imports the adapted REINVENT 4
components from ``reinvent_core`` (see ../NOTICE for attribution) and exposes
a small, well-typed interface for SMILES generation.

The prior checkpoint was originally pickled with class paths under the
``reinvent.*`` namespace. We register import-time aliases in ``sys.modules``
so torch.load can rehydrate the SMILESTokenizer / Vocabulary / ModelMetaData
instances against our local ``reinvent_core`` package without taking a hard
dependency on the upstream REINVENT 4 distribution.
"""

from __future__ import annotations

import sys
import types
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import List, Optional

import torch

# Resolve filesystem layout: this file lives at Reinvent_App/app/sampler.py;
# the package root and prior live one directory up.
_PKG_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_PRIOR = _PKG_ROOT / "prior" / "reinvent.prior"
_FINETUNED_DIR = _PKG_ROOT / "prior" / "finetuned"


def list_available_priors() -> dict[str, Path]:
    """Return {prior_name -> path} mapping for the base prior + any fine-tuned ones.

    Naming convention:
        prior/reinvent.prior                          → "base"
        prior/finetuned/<TARGET>_finetuned.prior      → "<TARGET>" (e.g. "EGFR")
    """
    out: dict[str, Path] = {}
    if _DEFAULT_PRIOR.is_file():
        out["base"] = _DEFAULT_PRIOR
    if _FINETUNED_DIR.is_dir():
        for f in _FINETUNED_DIR.iterdir():
            name = f.name
            if name.endswith("_finetuned.prior"):
                target = name[: -len("_finetuned.prior")]
                if target:
                    out[target] = f
    return out


def resolve_prior_path(prior_name: Optional[str]) -> Path:
    """Map a prior name (case-insensitive) to a checkpoint path.

    None / "" / "base" → base prior. Otherwise look for a fine-tuned model
    matching the target (e.g. "EGFR" → prior/finetuned/EGFR_finetuned.prior).
    Falls back to base prior if the requested target is not yet fine-tuned.
    """
    available = list_available_priors()
    if not prior_name or prior_name.lower() in {"base", "default", "generic"}:
        return available.get("base", _DEFAULT_PRIOR)
    # case-insensitive match
    for key, path in available.items():
        if key.lower() == prior_name.lower():
            return path
    # fallback
    return available.get("base", _DEFAULT_PRIOR)

# Make `reinvent_core` importable when this module is loaded as part of an
# external app (e.g. the AI Pharma backend) that may not have Reinvent_App on
# its sys.path yet.
if str(_PKG_ROOT) not in sys.path:
    sys.path.insert(0, str(_PKG_ROOT))


def _install_pickle_aliases() -> None:
    """Register module aliases so the prior pickle resolves to ``reinvent_core``.

    Idempotent: safe to call more than once.
    """
    if "reinvent.models.reinvent.models.vocabulary" in sys.modules:
        return  # already installed

    from reinvent_core import vocabulary, meta_data

    def _stub(name: str) -> types.ModuleType:
        if name in sys.modules:
            return sys.modules[name]
        m = types.ModuleType(name)
        m.__path__ = []  # mark as package so submodule lookup works
        sys.modules[name] = m
        return m

    _stub("reinvent")
    _stub("reinvent.models")
    sys.modules["reinvent.models.meta_data"] = meta_data
    _stub("reinvent.models.reinvent")
    _stub("reinvent.models.reinvent.models")
    sys.modules["reinvent.models.reinvent.models.vocabulary"] = vocabulary


@dataclass
class SamplerConfig:
    prior_path: Path = _DEFAULT_PRIOR
    device: str = "cpu"
    seed: Optional[int] = None  # None → fresh randomness each call


@dataclass
class SampledMolecule:
    smiles: str
    nll: float  # negative log likelihood under the prior; lower = more likely


class ReinventSampler:
    """Thread-safe wrapper around the REINVENT 4 prior — caches one model
    per prior path so multiple targets can be served by the same process."""

    def __init__(self, config: Optional[SamplerConfig] = None):
        self.config = config or SamplerConfig()
        # path → loaded Model. We allow multiple priors to coexist so a single
        # process can serve EGFR, JAK2, ... fine-tuned variants without restart.
        self._models: dict[Path, "object"] = {}
        self._lock = Lock()

    def _ensure_loaded(self, prior_path: Optional[Path] = None) -> "object":
        path = (prior_path or self.config.prior_path).resolve()
        if path in self._models:
            return self._models[path]
        with self._lock:
            if path in self._models:
                return self._models[path]
            _install_pickle_aliases()
            from reinvent_core.model import Model

            ckpt = torch.load(path, map_location=self.config.device, weights_only=False)
            model = Model.create_from_dict(
                ckpt, mode="inference", device=torch.device(self.config.device)
            )
            self._models[path] = model
            return model

    @property
    def _model(self):
        """Backwards-compatible alias for the default-prior model (only loaded
        once is_loaded has been triggered for the default path)."""
        return self._models.get(self.config.prior_path.resolve())

    @property
    def is_loaded(self) -> bool:
        return len(self._models) > 0

    @property
    def vocab_size(self) -> int:
        m = self._ensure_loaded()
        return len(m.vocabulary)

    @property
    def max_sequence_length(self) -> int:
        m = self._ensure_loaded()
        return m.max_sequence_length

    def sample(
        self,
        batch_size: int = 32,
        seed: Optional[int] = None,
        target_prior: Optional[str] = None,
    ) -> List[SampledMolecule]:
        """Sample ``batch_size`` molecules from a prior.

        :param batch_size: number of SMILES to draw (1..1024 reasonable)
        :param seed: optional torch seed for reproducibility
        :param target_prior: optional fine-tuned target name (e.g. "EGFR").
            None / "base" → use the base prior. Unknown names also fall back
            to base, so callers can safely pass any target_name.
        """
        if batch_size < 1 or batch_size > 4096:
            raise ValueError("batch_size must be between 1 and 4096")

        prior_path = resolve_prior_path(target_prior)
        model = self._ensure_loaded(prior_path)
        if seed is None:
            seed = self.config.seed
        if seed is not None:
            torch.manual_seed(seed)

        # The model holds parameters under torch grad mode; inference only.
        with torch.no_grad():
            _seqs, smiles_list, nlls = model.sample(batch_size=batch_size)

        out: List[SampledMolecule] = []
        for smi, nll in zip(smiles_list, nlls.tolist()):
            out.append(SampledMolecule(smiles=smi, nll=float(nll)))
        return out


# Module-level default sampler — most callers only need one.
_default_sampler: Optional[ReinventSampler] = None


def get_default_sampler() -> ReinventSampler:
    global _default_sampler
    if _default_sampler is None:
        _default_sampler = ReinventSampler()
    return _default_sampler
