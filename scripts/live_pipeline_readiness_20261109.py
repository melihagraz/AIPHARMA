#!/usr/bin/env python3
"""Actually exercise available REINVENT weights; inspect ADMET/Vina prerequisites.

Do not represent status 'blocked' as an end-to-end scientific validation.
Use local TestClient for REINVENT so model weights are loaded and /sample is
executed without starting network listeners or deploying public services.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import platform
from pathlib import Path
import shutil
import sys
import time
from importlib.metadata import version,PackageNotFoundError

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
OUT=ROOT/'reports/validation_20261109_extended/live_service_readiness.json'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    from fastapi.testclient import TestClient
    from Reinvent_App.app.main import app as reinvent_app
    from Docking_App.app.dock import _vina_binary, list_targets
    from rdkit import Chem
    from rdkit.Chem import QED
    prior=ROOT/'Reinvent_App/prior/reinvent.prior'
    r={'type':'direct_local_model_execution_and_prerequisite_audit',
       'environment':{'python':sys.version.split()[0], 'platform':platform.platform(),
                      'rdkit':version('rdkit'),'torch':version('torch')},
       'stages':{},'full_live_reinvent_admet_vina_e2e_verified':False,
       'methodological_limitations':['No external ADMET inference or Vina docking was executed',
           'Molecular validity and QED do not establish ALK5 activity, safety or synthesizability']}
    start=time.perf_counter()
    request={'batch_size':16,'seed':20261109,'target_prior':'base','validate_with_rdkit':True}
    with TestClient(reinvent_app) as client:
        a=client.post('/sample',json=request)
        b=client.post('/sample',json=request)
    assert a.status_code==200,(a.status_code,a.text[:500]);assert b.status_code==200
    data=a.json();second=b.json()
    smi=[c['smiles'] for c in data['candidates']]
    assert smi == [c['smiles'] for c in second['candidates']]
    valid=[Chem.MolFromSmiles(x) for x in smi]
    valid_count=sum(x is not None for x in valid)
    assert valid_count==data['valid_count']
    r['stages']['reinvent']={
        'status':'passed_real_model_inference',
        'model_prior_sha256':sha(prior),
        'prior_used':data['prior_used'],
        'requested':16,'returned':len(smi),'valid_rdkit':valid_count,
        'unique_smiles':len(set(smi)),'replication_same_seed_smiles_identical':True,
        'mean_qed_valid':sum(QED.qed(m) for m in valid if m is not None)/valid_count if valid_count else None,
        'runtime_two_batches_seconds':round(time.perf_counter()-start,3),
        'not_an_alk5_activity_test':True,
    }
    try:
        present=importlib.util.find_spec('admet_ai') is not None
    except ImportError:
        present=False
    r['stages']['admet_ai']={'status':'available_for_execution' if present else 'blocked_missing_python_package',
                              'admet_ai_installed':present,'real_inference_performed':False}
    binary=None; error=None
    try:
        binary=_vina_binary()
    except FileNotFoundError as exc:
        error=str(exc)
    r['stages']['vina']={'status':'available_for_execution' if binary else 'blocked_missing_platform_executable',
                         'vina_executable':binary, 'error':error,
                         'meeko_installed':importlib.util.find_spec('meeko') is not None,
                         'available_receptor_targets':list_targets(),
                         'ALK5_bundled_receptor': 'ALK5' in list_targets(),
                         'real_docking_performed':False}
    r['full_live_reinvent_admet_vina_e2e_verified']=False
    OUT.parent.mkdir(exist_ok=True,parents=True)
    OUT.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v['status'] for k,v in r['stages'].items()},indent=2))
    print('Wrote',OUT)

if __name__=='__main__': main()
