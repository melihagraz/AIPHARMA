#!/usr/bin/env python3
"""Offline, independently rerunnable ALK5 case-study verification.

Reads archived prediction outputs. Does not invoke ADMET-AI, REINVENT or Vina.
Outputs are retrospective arithmetic checks, not new molecule predictions.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import statistics
import sys
from importlib.metadata import version, PackageNotFoundError


def sha256(file: Path) -> str:
    return hashlib.sha256(file.read_bytes()).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=root/'reports/software_validation_20261109')
    args=parser.parse_args()
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(root/'backend'))
    from app.services.alk5_analogue_design import design_analogues, ASSET
    d1=design_analogues(30)
    d2=design_analogues(30)
    assert [x['smiles'] for x in d1['candidates']]==[x['smiles'] for x in d2['candidates']]
    assert d1['returned']==30
    assert len(set(x['normalized_parent'] for x in d1['candidates']))==30
    assert all(x['activity_status']=='unmeasured' for x in d1['candidates'])
    report_dir=root/'reports/alk5_design_admet_20261019'
    csvpath=report_dir/'paired_outputs.csv'
    rows=list(csv.DictReader(csvpath.open(newline='',encoding='utf-8')))
    groups={}
    for row in rows:
        endpoint=row['endpoint']
        a,b=float(row['parent']),float(row['candidate'])
        delta=float(row['delta'])
        assert abs((b-a)-delta)<1e-9, f'invalid paired delta {row}'
        groups.setdefault(endpoint,[]).append(delta)
    assert len(rows)==180 and len(groups)==6
    summaries={}
    old=json.loads((report_dir/'descriptive_summary.json').read_text())
    for name,vals in sorted(groups.items()):
        item=dict(lower=sum(v < -1e-10 for v in vals),higher=sum(v > 1e-10 for v in vals),unchanged=sum(abs(v)<=1e-10 for v in vals),
                  median_delta=statistics.median(vals),min_delta=min(vals),max_delta=max(vals),n=len(vals))
        assert item['n']==30, name
        for key in ('lower','higher','unchanged'):
            assert item[key]==old[name][key],(name,key)
        assert abs(item['median_delta']-old[name]['median_delta'])<1e-9,name
        summaries[name]=item
    packages={}
    for m in ('rdkit','fastapi','pydantic','argon2-cffi','pytest'):
        try:packages[m]=version(m)
        except PackageNotFoundError:packages[m]=None
    result={
        'verification_type':'offline_source_recalculation',
        'scientific_boundary':'does not validate biological ALK5 inhibition, safety, synthesis or clinical efficacy',
        'design':{'returned':d1['returned'],**d1['metrics'],
                  'identity_prefix_repeated':True,
                  'generated_smiles_sha256':hashlib.sha256('\n'.join(x['smiles'] for x in d1['candidates']).encode()).hexdigest(),
                  'reference_asset_sha256':sha256(ASSET)},
        'paired_admet':{'paired_observations':len(rows),'endpoint_count':len(groups),'summary':summaries,
                        'source_sha256':sha256(csvpath)},
        'environment':{'python':sys.version.split()[0],'platform':platform.platform(),'packages':packages},
        'code_sha256':{str((root/'backend/app/services/alk5_analogue_design.py').relative_to(root)):sha256(root/'backend/app/services/alk5_analogue_design.py'),
                       str((root/'backend/app/services/alk5_design_admet.py').relative_to(root)):sha256(root/'backend/app/services/alk5_design_admet.py')},
        'run_status':'passed',
    }
    (out/'revalidation_results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    (out/'revalidated_designs.json').write_text(json.dumps(d1,indent=2,ensure_ascii=False)+'\n')
    print('ALK5 revalidation PASS:',result['design']['returned'],'designs,',len(rows),'paired endpoint rows; files in',out)

if __name__=='__main__':main()
