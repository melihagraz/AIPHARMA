"""Offline regression guards for security and evidence labeling.

These tests do not need REINVENT, docking, Celery, external datasets or network.
"""
import os
from pathlib import Path
import subprocess
import sys

import pytest
from fastapi import HTTPException

from app.core.config import Settings, get_settings
from app.core.security import (hash_password, verify_password,
                               issue_session_token, decode_session_token, require_session)
from app.services.production_pipeline import UnsupportedFilterError, validate_supported_filters
from app.universal_discovery.pipeline import _make_smiles


@pytest.fixture()
def session_key(monkeypatch):
    monkeypatch.setenv('SESSION_SIGNING_KEY', 'offline-test-secret-ALWAYS-NOT-FOR-PRODUCTION-123456')
    get_settings.cache_clear()
    yield
    get_settings.cache_clear()


def test_admin_defaults_not_open():
    cfg = Settings(_env_file=None)
    assert cfg.admin_password == ''
    assert cfg.admin_token == ''
    assert cfg.universal_discovery_mode == 'production'


def test_production_rejects_missing_session_key():
    with pytest.raises(ValueError, match='SESSION_SIGNING_KEY'):
        Settings(_env_file=None, app_env='production', session_signing_key='')


def test_argon2id_passwords_and_rejections():
    pwd = 'long-and-private-passphrase'
    hashed = hash_password(pwd)
    assert hashed.startswith('$argon2id$')
    assert pwd not in hashed
    assert verify_password(hashed, pwd)
    assert not verify_password(hashed, pwd+'!')
    assert not verify_password('a-legacy-plaintext-password', 'a-legacy-plaintext-password')
    with pytest.raises(ValueError):
        hash_password('short')


def test_signed_session_and_reject_tampering(session_key):
    token = issue_session_token('alice')
    assert decode_session_token(token).user_id == 'alice'
    assert decode_session_token(token).role == 'user'
    assert require_session('Bearer '+token).user_id == 'alice'
    with pytest.raises(HTTPException) as e:
        decode_session_token(token+'Z')
    assert e.value.status_code == 401
    with pytest.raises(HTTPException) as e:
        require_session(None)
    assert e.value.status_code == 401


def test_unsupported_filters_cannot_silently_pass():
    with pytest.raises(UnsupportedFilterError, match='bsep_risk_max'):
        validate_supported_filters({'safety_filters': {'bsep_risk_max': 0.4}})
    with pytest.raises(UnsupportedFilterError, match='half_life_range'):
        validate_supported_filters({'pk_filters': {'half_life_range': [1,4]}})
    assert validate_supported_filters({'safety_filters': {'herg_risk_max': 0.4}}) is None


def test_demo_seed_determinism_across_processes():
    a = _make_smiles('demo', 4, 'ALK5')
    code = "from app.universal_discovery.pipeline import _make_smiles; print(_make_smiles('demo',4,'ALK5'))"
    env = os.environ.copy()
    env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1])
    for hash_seed in ('1', '2'):
        env['PYTHONHASHSEED'] = hash_seed
        b = subprocess.check_output([sys.executable, '-c', code], env=env, text=True).strip()
        assert a == b


def test_missing_admet_descriptors_remain_null():
    from app.services.production_pipeline import _format_candidate
    row = _format_candidate(1, 'CCO', 25.0, {'qed':0.6, 'mw':46.07, 'logp':-0.1, 'tpsa':20.2},
                            {}, [], 'general', max_tanimoto=None, dock_affinity=None)
    assert row['descriptors']['hERG'] is None
    assert row['descriptors']['DILI'] is None
    assert row['dock_score'] is None
    assert row['descriptors']['qed'] == .6


@pytest.mark.asyncio
async def test_explicit_demo_response_is_marked_synthetic(monkeypatch):
    from app.api import routes_legacy
    from app.universal_discovery.models import DiscoveryRequest
    async def synthetic(_req):
        return ({'target':'ALK5','module':'general','total_generated':0,
                 'after_filter_count':0,'top_n':0,'candidates':[],
                 'meta':{}},'embedded')
    monkeypatch.setattr(routes_legacy,'_discover_upstream',synthetic)
    body=DiscoveryRequest(target_name='ALK5',engine_mode='demo')
    result=await routes_legacy.motor_discover_impl(body)
    assert result.meta['is_synthetic'] is True
    assert result.meta['research_validated'] is False
    assert result.meta['engine_type'] == 'synthetic_demo'
    assert result.meta['run_id']
    assert result.meta['degraded'] is True
