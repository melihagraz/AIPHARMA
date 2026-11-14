"""Strict upstream-response, rate-limit and platform-binary regression tests.

MockTransport checks network contracts only; it does NOT execute ML inference.
"""
import httpx
import pytest
from fastapi import HTTPException
from starlette.requests import Request

from app.services import admet_client, reinvent_client, docking_client
from app.core import login_throttle


def _fake_httpx(monkeypatch, json_body):
    original = httpx.AsyncClient
    transport = httpx.MockTransport(lambda request: httpx.Response(200, json=json_body))
    monkeypatch.setattr(httpx, 'AsyncClient', lambda *args, **kwargs: original(*args, transport=transport, **kwargs))


@pytest.mark.asyncio
async def test_admet_rejects_missing_rows(monkeypatch):
    _fake_httpx(monkeypatch, {'count': 1, 'predictions': [{'smiles':'CCO', 'raw':{'hERG':0.3},'percentile':{}}]})
    with pytest.raises(HTTPException, match='count differs') as err:
        await admet_client.predict_admet(['CCO','O'])
    assert err.value.status_code == 502


@pytest.mark.asyncio
async def test_admet_rejects_wrong_smiles_and_invalid_prediction(monkeypatch):
    _fake_httpx(monkeypatch, {'count':1,'predictions':[{'smiles':'N','raw':{'hERG':0.3},'percentile':{}}]})
    with pytest.raises(HTTPException, match='identity mismatch'):
        await admet_client.predict_admet(['CCO'])


@pytest.mark.asyncio
async def test_admet_good_response(monkeypatch):
    _fake_httpx(monkeypatch, {'count':1,'predictions':[{'smiles':'CCO','raw':{'hERG':0.3,'AMES':None},'percentile':{}}]})
    assert (await admet_client.predict_admet(['CCO']))[0]['raw']['hERG'] == 0.3


@pytest.mark.asyncio
async def test_reinvent_rejects_invalid_sample_count(monkeypatch):
    _fake_httpx(monkeypatch, {'requested':5,'returned':0,'candidates':[],'prior_used':'base'})
    with pytest.raises(HTTPException, match='counts are inconsistent'):
        await reinvent_client.sample_smiles(batch_size=4)


@pytest.mark.asyncio
async def test_reinvent_accepts_well_formed_sample(monkeypatch):
    _fake_httpx(monkeypatch, {'requested':1,'returned':1,'candidates':[{'smiles':'CCO','valid':True,'nll':21.0}],'prior_used':'base'})
    candidates,used = await reinvent_client.sample_smiles(batch_size=1)
    assert used=='base' and candidates[0]['smiles']=='CCO'


@pytest.mark.asyncio
async def test_docking_rejects_misaligned_result(monkeypatch):
    _fake_httpx(monkeypatch, {'results':[{'smiles':'N','target':'EGFR','affinity':-6.2}]})
    with pytest.raises(HTTPException, match='out of input SMILES order'):
        await docking_client.dock_batch(['CCO'],target='EGFR')


@pytest.mark.asyncio
async def test_docking_accepts_failed_pose_as_explicit_result(monkeypatch):
    _fake_httpx(monkeypatch, {'results':[{'smiles':'CCO','target':'EGFR','affinity':None,'error':'embedding failed'}]})
    assert (await docking_client.dock_batch(['CCO'],target='EGFR'))[0]['error']=='embedding failed'


def test_login_attempt_throttle_blocks_brute_force_and_ignores_forwarded_for():
    login_throttle._BUCKETS.clear()
    scope = {'type':'http','client':('192.0.2.44',50000), 'method':'POST','path':'/api/auth/login',
             'headers': [(b'x-forwarded-for',b'127.0.0.1')]}
    req = Request(scope)
    for i in range(8):
        login_throttle.check_login_throttle(req,'alice@example.org')
        login_throttle.record_login_failure(req,'alice@example.org')
    with pytest.raises(HTTPException) as err:
        login_throttle.check_login_throttle(req,'alice@example.org')
    assert err.value.status_code == 429
    assert '192.0.2.44' in ''.join(login_throttle._BUCKETS)
    assert '127.0.0.1' not in ''.join(login_throttle._BUCKETS)
    login_throttle._BUCKETS.clear()


def test_linux_does_not_pick_windows_vina(monkeypatch,tmp_path):
    import os
    from pathlib import Path
    import sys
    sys.path.insert(0,str(__import__("pathlib").Path(__file__).resolve().parents[2]))
    from Docking_App.app import dock
    (tmp_path / 'bin').mkdir()
    (tmp_path / 'bin' / 'vina.exe').write_bytes(b'MZ')
    monkeypatch.setattr(dock, '_PKG_ROOT',tmp_path)
    monkeypatch.setattr(dock.shutil,'which',lambda binary:None)
    if os.name != 'nt':
        with pytest.raises(FileNotFoundError,match='platform-appropriate'):
            dock._vina_binary()


def test_production_requires_internal_service_secret():
    from app.core.config import Settings
    from pydantic import ValidationError
    with pytest.raises(ValidationError,match='INTERNAL_SERVICE_KEY'):
        Settings(_env_file=None, app_env='production', session_signing_key='x'*40, internal_service_key='')


def test_production_microservice_endpoints_fail_closed(monkeypatch):
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from fastapi.testclient import TestClient
    from Reinvent_App.app.main import app as reinv
    from Admet_App.app.main import app as admet
    from Docking_App.app.main import app as dock
    monkeypatch.setenv('APP_ENV','production')
    monkeypatch.setenv('INTERNAL_SERVICE_KEY','offline-test-internal-key-length-more-than-32')
    for app,route in [(reinv,'/priors'), (admet,'/health'), (dock,'/targets')]:
        with TestClient(app) as client:
            assert client.get(route).status_code==401
            assert client.get(route,headers={'X-Internal-Service-Key':'bad'}).status_code==401
            ok=client.get(route,headers={'X-Internal-Service-Key':'offline-test-internal-key-length-more-than-32'})
            assert ok.status_code==200,ok.text
