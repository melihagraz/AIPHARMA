"""Regression tests for Phase 4: no silent acceptance of unusable scientific values.

No external models or HTTP mocks; this tests software contracts, NOT
REINVENT/ADMET/Vina predictive performance.
"""
from __future__ import annotations

import math
import pytest
from pydantic import ValidationError

from app.core.config import Settings
from app.services.production_pipeline import _check_filters, _normalize_nll


def test_unknown_app_environment_is_not_implicitly_development():
    with pytest.raises(ValidationError, match="APP_ENV must be"):
        Settings(_env_file=None, app_env="prod")


def test_production_debug_is_rejected_even_with_valid_secrets():
    with pytest.raises(ValidationError, match="DEBUG must be disabled"):
        Settings(_env_file=None, app_env="production", debug=True,
                 session_signing_key="s"*48, internal_service_key="i"*48)

@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_reinvent_nll_rejected(invalid):
    with pytest.raises(ValueError, match="likelihood must be finite"):
        _normalize_nll(invalid)

@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_admet_prediction_never_satisfies_safety_threshold(invalid):
    failures, _ = _check_filters({"mw": 300.0}, {"hERG": invalid},
                                 {"safety_filters": {"herg_risk_max": .5}}, smiles="CCO")
    assert "hERG_missing" in failures

@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -float('inf')])
def test_nonfinite_rdkit_descriptor_never_satisfies_base_threshold(invalid):
    failures, _ = _check_filters({"mw": invalid}, {},
                                 {"base_filters": {"mw_max": 450.0}}, smiles="CCO")
    assert "mw_missing" in failures


def test_missing_smiles_with_hbd_filter_fails_closed():
    failures, _ = _check_filters({"mw": 330}, {},
                                 {"base_filters": {"hbd_max": 3}}, smiles=None)
    assert "hbd_missing" in failures


def test_missing_smiles_with_sa_filter_fails_closed():
    failures, _ = _check_filters({"mw": 330}, {},
                                 {"base_filters": {"sa_score_max": 4.0}})
    assert "sa_score_missing" in failures


def test_valid_finite_base_descriptor_kept():
    failures, _ = _check_filters({"mw": 350}, {},
                                 {"base_filters": {"mw_max": 450.0}})
    assert failures == []
