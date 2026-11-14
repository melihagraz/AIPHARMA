"""ADMET predictor wrapper around ADMET-AI (MIT licensed, pip dependency).

Loads the ensemble Chemprop models once and exposes a small thread-safe
``predict(smiles_list)`` interface. Output structure splits the 104 raw
columns into two groups (raw values + drugbank percentiles) so callers can
present them sensibly.
"""

from __future__ import annotations

from dataclasses import dataclass
from threading import Lock
from typing import Any, Dict, List, Optional
from collections import defaultdict, deque
import math


PERCENTILE_SUFFIX = "_drugbank_approved_percentile"


@dataclass
class AdmetPrediction:
    smiles: str
    raw: Dict[str, float]            # 52 raw model outputs
    percentile: Dict[str, float]     # 52 DrugBank-approved percentiles


class AdmetPredictor:
    """Thread-safe lazy-loaded wrapper around ``admet_ai.ADMETModel``.

    First call to ``predict`` triggers model download (one-time, cached by
    ADMET-AI under the user cache dir). Subsequent calls reuse the loaded
    ensemble.
    """

    def __init__(self, num_workers: int = 0):
        self._model = None
        self._lock = Lock()
        self._num_workers = num_workers

    def _ensure_loaded(self) -> None:
        if self._model is not None:
            return
        with self._lock:
            if self._model is not None:
                return
            from admet_ai import ADMETModel  # heavy import — defer

            self._model = ADMETModel(num_workers=self._num_workers)

    @property
    def is_loaded(self) -> bool:
        return self._model is not None

    def predict(self, smiles: List[str]) -> List[AdmetPrediction]:
        if not smiles:
            return []
        if any(not isinstance(s, str) or not s.strip() for s in smiles):
            raise ValueError("All SMILES must be non-empty strings")

        self._ensure_loaded()
        df = self._model.predict(smiles=smiles)

        # Verify row identity instead of assuming predictions preserve input order.
        # An incorrect mapping makes all molecule/ADMET conclusions invalid.
        df = _align_predictions_to_inputs(df, smiles)

        # Split columns into raw vs percentile
        all_cols = list(df.columns)
        percentile_cols = [c for c in all_cols if c.endswith(PERCENTILE_SUFFIX)]
        raw_cols = [c for c in all_cols if not c.endswith(PERCENTILE_SUFFIX) and c != "smiles"]

        if not raw_cols:
            raise ValueError("ADMET model returned no raw endpoints")
        out: List[AdmetPrediction] = []
        for i, smi in enumerate(smiles):
            row = df.iloc[i]
            raw = {c: _coerce_float(row[c]) for c in raw_cols}
            pct = {
                c[: -len(PERCENTILE_SUFFIX)]: _coerce_float(row[c]) for c in percentile_cols
            }
            out.append(AdmetPrediction(smiles=smi, raw=raw, percentile=pct))
        return out


def _align_predictions_to_inputs(df, smiles: List[str]):
    """Require one prediction per input, including duplicates and reordered rows.

    ADMET-AI commonly returns SMILES in the index. A RangeIndex without a
    SMILES column is supported only as an explicitly position-ordered result.
    Unexpected indices and partially missing rows fail closed.
    """
    if len(df) != len(smiles):
        raise ValueError(f"ADMET row count mismatch: requested {len(smiles)}, got {len(df)}")
    if "smiles" in df.columns:
        identities = [str(v) for v in df["smiles"].tolist()]
    elif df.index.name == "smiles" or all(isinstance(v, str) for v in df.index):
        identities = [str(v) for v in df.index.tolist()]
    elif list(df.index) == list(range(len(smiles))):
        # No row identifiers: the library's output order is the only possible mapping.
        return df.reset_index(drop=True)
    else:
        raise ValueError("ADMET model returned rows without verifiable molecule order")
    slots = defaultdict(deque)
    for i, smi in enumerate(identities):
        slots[smi].append(i)
    try:
        positions = [slots[smi].popleft() for smi in smiles]
    except (IndexError, KeyError):
        raise ValueError("ADMET model returned missing/unexpected SMILES identities") from None
    if any(slots.values()):
        raise ValueError("ADMET model returned extra/duplicate prediction identities")
    return df.iloc[positions].reset_index(drop=True)


def _coerce_float(v: Any) -> Optional[float]:
    """Convert numpy scalars / pandas NA to plain float (or None)."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    # Missing and non-finite results are explicitly null, never "zero".
    if not math.isfinite(f):
        return None
    return f


# Module-level default
_default_predictor: Optional[AdmetPredictor] = None


def get_default_predictor() -> AdmetPredictor:
    global _default_predictor
    if _default_predictor is None:
        _default_predictor = AdmetPredictor()
    return _default_predictor
