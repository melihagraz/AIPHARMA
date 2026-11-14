"""Real main API + real ADMET HTTP service; isolated DB; leave running for UI QA.

Does not copy .env files or use cloud LLM services. No mocked prediction route.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from uuid import uuid4

ap = argparse.ArgumentParser()
ap.add_argument('root', type=Path)
ap.add_argument('--api-port', type=int, default=18765)
ap.add_argument('--admet-port', type=int, default=18766)
args = ap.parse_args()
root, source = args.root, args.root/'source'
out = root/'evidence'/uuid4().hex[:8]
out.mkdir(parents=True)
base = f'http://127.0.0.1:{args.api_port}'
admet = f'http://127.0.0.1:{args.admet_port}'
db = out/'fresh.sqlite3'
assert not db.exists()
assert not (source/'.env').exists() and not (source/'backend/.env').exists()
protocol = dict(scope='Fresh main-backend venv and DB; production main app and real ADMET HTTP service; built React UI', scientific_environment='Existing aipharma conda environment and cached ADMET weights; NOT a clean scientific-model installation', excluded=['LLM chat', 'REINVENT sampling', 'pocket detection', 'docking', 'cross-machine or public deployment'], repeated_design='exact full JSON equality', repeated_admet='same design identities; maximum absolute endpoint difference <=1e-5, fixed before execution', fresh_database=True, secrets_copied=False, api_base=base, admet_base=admet, python=platform.python_version())
(out/'protocol.json').write_text(json.dumps(protocol, indent=2))
processes, handles = [], []
def start(name, command, cwd, env):
    log = (out/f'{name}.log').open('w')
    handles.append(log)
    p = subprocess.Popen(command, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
    processes.append(p)
    return p
env = os.environ.copy()
env.update(DATABASE_URL=f'sqlite:///{db}', ANTHROPIC_API_KEY='', ADMIN_USERNAME='delivery-admin', ADMIN_PASSWORD='local-delivery-only', ADMIN_TOKEN='delivery-local-token', ADMET_APP_URL=admet, ENABLE_ASYNC_DISPATCH='false', UNIVERSAL_DISCOVERY_MODE='production', OMP_NUM_THREADS='4', MKL_NUM_THREADS='4')
# The scientific runtime already contains ADMET model dependencies; only
# the HTTP framework lives in an isolated target directory.
scientific_env = env.copy()
scientific_env['PYTHONPATH'] = '/home/cobre/tools/aipharma-design-api-test'
evidence = []
def request(path, body=None, name=None, expected=200):
    data = json.dumps(body).encode() if body is not None else None
    started = time.perf_counter()
    req = Request(base+path, data=data, headers={'Content-Type': 'application/json'} if data else {})
    try:
        with urlopen(req, timeout=240) as response:
            status, raw = response.status, response.read()
    except HTTPError as exc:
        status, raw = exc.code, exc.read()
    evidence.append(dict(path=path.split('?')[0], status=status, seconds=time.perf_counter()-started, bytes=len(raw)))
    assert status == expected, (path, status, raw[:300])
    if name:
        (out/name).write_bytes(raw)
    return raw
def wait_health(url, process):
    deadline = time.monotonic()+90
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError('Server exited; inspect isolated logs')
        try:
            with urlopen(url, timeout=2) as response:
                if response.status == 200:
                    return
        except (URLError, TimeoutError):
            pass
        time.sleep(.5)
    raise RuntimeError('Server health timeout')
def shutdown(*_):
    raise KeyboardInterrupt
signal.signal(signal.SIGTERM, shutdown)
try:
    backend = start('main_api', [str(root/'venv/bin/python'), '-m', 'uvicorn', 'app.main:app', '--host', '127.0.0.1', '--port', str(args.api_port), '--app-dir', str(source/'backend')], source, env)
    predictor = start('admet_service', ['/home/cobre/miniforge3/envs/aipharma/bin/python', '-m', 'uvicorn', 'app.main:app', '--host', '127.0.0.1', '--port', str(args.admet_port), '--app-dir', str(source/'Admet_App')], source/'Admet_App', scientific_env)
    wait_health(base+'/api/v1/health', backend)
    wait_health(admet+'/health', predictor)
    request('/api/v1/health', name='main_health.json')
    assert json.loads(request('/api/v1/jobs', name='empty_jobs.json')) == []
    assert b'<html' in request('/', name='served_index.html').lower()
    request('/openapi.json', name='openapi.json')
    first = json.loads(request('/api/motor/alk5/analogues', dict(count=30), 'designs_http_1.json'))
    second = json.loads(request('/api/motor/alk5/analogues', dict(count=30), 'designs_http_2.json'))
    assert first == second and first['returned'] == 30
    for n in [0, 101, True]:
        request('/api/motor/alk5/analogues', dict(count=n), expected=422)
    png = request('/api/util/mol-image?'+urlencode(dict(smiles=first['candidates'][0]['smiles'], w=160, h=112)), name='served_molecule.png')
    assert png.startswith(b'\x89PNG') and len(png) > 200
    profiles = []
    for repeat in [1, 2]:
        profile = json.loads(request('/api/motor/alk5/analogues/admet', dict(count=30), f'profiles_http_{repeat}.json'))
        assert len(profile['pairs']) == 30 and profile['meta']['missing_prediction_structures'] == 0
        assert all(p['smiles'] == c['smiles'] for p, c in zip(profile['pairs'], first['candidates']))
        profiles.append(profile)
    deltas = [abs(a['endpoints'][e][k]-b['endpoints'][e][k]) for a, b in zip(profiles[0]['pairs'], profiles[1]['pairs']) for e in profiles[0]['endpoints'] for k in ['parent', 'candidate', 'delta']]
    assert max(deltas) <= 1e-5
    # Verify existing sign-in without saving a bearer token in evidence.
    login = json.loads(request('/api/auth/login', dict(username='delivery-admin', password='local-delivery-only')))
    assert login.get('access_token') and login.get('user')
    result = dict(status='passed', real_main_application=True, mocked_predictions=False, design_count=30, repeated_designs_equal=True, endpoint_pairs=180, max_repeated_admet_absolute_difference=max(deltas), main_backend_clean_venv=True, model_environment_reused=True, fresh_database=True, sign_in_http=True, frontend_served=True, molecule_png_served=True, requests=evidence)
    (out/'http_verification.json').write_text(json.dumps(result, indent=2))
    (out/'backend_pip_freeze.txt').write_bytes(subprocess.check_output([str(root/'venv/bin/python'), '-m', 'pip', 'freeze']))
    (root/'active_evidence.json').write_text(json.dumps(dict(path=str(out), api_base=base, pids=[p.pid for p in processes]), indent=2))
    print(json.dumps(dict(status='HTTP checks passed; servers ready for browser QA', evidence=str(out), api_base=base, max_repeated_admet_absolute_difference=max(deltas))), flush=True)
    while all(p.poll() is None for p in processes):
        time.sleep(1)
except KeyboardInterrupt:
    pass
finally:
    for process in processes:
        if process.poll() is None:
            process.terminate()
    for process in processes:
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
    for handle in handles:
        handle.close()
