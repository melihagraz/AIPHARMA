#!/usr/bin/env python3
"""AIPHARMA Phase 16: on-copy, no-docking, fault-injection audit comparison.

Direct comparison of five *local control policies*, NOT a performance test of
Icolos or DDA. All perturbations affect temporary copies of existing real
benchmark records, never the original records. Requires the Phase 13/14 EAVE
modules already installed in the AIPHARMA project.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aipharma_eave.audit import Audit, CASES  # noqa: E402
from aipharma_eave.phase14 import gate, trusted_reference, digest_comparison  # noqa: E402

POLICIES = ('presence_only', 'metric_only', 'snapshot_only', 'eave_without_snapshot', 'eave_full')
SCENARIOS = (
    ('clean', False),
    ('unrelated_note', False),
    ('source_modified', True),
    ('wrong_receptor', True),
    ('duplicate_score', True),
    ('falsified_auc', True),
    ('changed_score_same_auc', True),
    ('missing_source', True),
    ('unsupported_potency_claim', True),
)

def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, data: dict):
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')


def csv_load(path: Path):
    with path.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        return reader.fieldnames,list(reader)


def csv_save(path: Path, fields, rows):
    with path.open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def auc_by_pairwise(active: list[float], decoy: list[float]) -> float:
    if not active or not decoy:raise ValueError('empty class')
    return sum((a < d) + .5*(a == d) for a in active for d in decoy)/(len(active)*len(decoy))


def stage_case(root: Path, case: str, destination: Path):
    """Only copy pinned input/output files into scratch; no data distributed."""
    baseline, error=trusted_reference(root,case)
    if error: raise ValueError(error)
    if not baseline['artifact_sha256']:raise ValueError('empty reference')
    for rel, old in baseline['artifact_sha256'].items():
        src=root/rel
        if not src.is_file():raise FileNotFoundError(str(src))
        if digest(src)!=old:raise ValueError(f'Input differs from trusted baseline before experiments: {rel}')
        dest=destination/rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest)
    # Also copy case redocking/selection records if not pinned by original baseline.
    cfg=CASES[case]
    for rel in (cfg['panel'],cfg['descriptors'],cfg['selection'],
                cfg['output']+'/actual_vina_scores.csv',
                cfg['output']+'/execution_manifest.json',
                cfg['output']+'/analysis.json',
                cfg['redock']+'/redock_report.json'):
        src=root/rel;dest=destination/rel
        if not dest.exists() and src.is_file():
            dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dest)
    bfile=destination/'phase14_eave_baselines'/f'{case}.json'
    bfile.parent.mkdir(parents=True,exist_ok=True)
    write_json(bfile,baseline)
    # The staged clean case must be accepted by the real engine, not a mock.
    rc,report=gate(destination,case,bootstrap=0,write=False)
    if rc!=0:raise RuntimeError(f'Clean staged records not verified: {case}: {report["gate_verdict"]}')
    return baseline


def mutate(staged: Path, case: str, scenario: str):
    cfg=CASES[case]
    if scenario=='clean':return
    if scenario=='unrelated_note':
        p=staged/'test_notes/independent_readme.txt'
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text('Unrelated new text; not part of frozen scientific evidence.\n')
    elif scenario=='source_modified':
        p=staged/cfg['sources']/'active.ism'
        with p.open('a',encoding='utf8') as f:f.write('\n# INJECTED TEST MUTATION (NOT SCIENTIFIC DATA)\n')
    elif scenario=='wrong_receptor':
        p=staged/cfg['output']/'execution_manifest.json'
        j=json.loads(p.read_text());j['receptor_pdb']='INTENTIONALLY_WRONG_PDB';write_json(p,j)
    elif scenario in ('duplicate_score','changed_score_same_auc'):
        p=staged/cfg['output']/'actual_vina_scores.csv'
        fields,rows=csv_load(p)
        if scenario=='duplicate_score':
            rows.append(dict(rows[0]))
        else:
            # Select an UN-TIED score; the previous first-row perturbation
            # broke an EGFR active--decoy tie (-8.656/-8.656) and raised AUC
            # from 0.6144 to 0.6146, contradicting the scenario name.
            # No opposite-class ordering or ties may change in this scenario.
            eps = 1e-8
            values = [float(r['score_kcal_mol']) for r in rows]
            active = [v for r,v in zip(rows,values) if r['class']=='active']
            decoy = [v for r,v in zip(rows,values) if r['class']=='decoy']
            original_auc = auc_by_pairwise(active,decoy)
            changed = False
            for idx, val in enumerate(values):
                if not math.isfinite(val):
                    continue
                # Require a gap > 100 perturbation units from all scores in
                # the OPPOSITE class; within-class ties cannot alter AUC.
                if any(abs(val-other) <= 100*eps for j,other in enumerate(values)
                       if rows[j]['class'] != rows[idx]['class']):
                    continue
                rows[idx]['score_kcal_mol'] = format(val-eps,'.10f')
                modified = [float(r['score_kcal_mol']) for r in rows]
                new_active = [v for r,v in zip(rows,modified) if r['class']=='active']
                new_decoy = [v for r,v in zip(rows,modified) if r['class']=='decoy']
                if auc_by_pairwise(new_active,new_decoy) == original_auc:
                    changed = True
                    break
                rows[idx]['score_kcal_mol'] = format(val,'.10f')
            if not changed:
                raise ValueError('No suitable rank-neutral score perturbation found')
        csv_save(p,fields,rows)
    elif scenario=='falsified_auc':
        p=staged/cfg['output']/'analysis.json'
        j=json.loads(p.read_text());key='roc_auc_3HMM' if case=='alk5-3hmm' else 'roc_auc'
        j[key]=0.1111 if float(j[key])>0.3 else 0.9999
        write_json(p,j)
    elif scenario=='missing_source':(staged/cfg['sources']/'active.ism').unlink()
    elif scenario=='unsupported_potency_claim':return # An explicit new claim request, not a file modification.
    else:raise ValueError(scenario)


def evaluate(root:Path,case:str,scenario:str,ref:dict):
    cfg=CASES[case]
    def presence_only():
        required=[cfg['sources']+'/active.ism',cfg['sources']+'/decoy.ism',
                  cfg['panel'],cfg['output']+'/execution_manifest.json',
                  cfg['output']+'/analysis.json',cfg['output']+'/actual_vina_scores.csv',
                  cfg['redock']+'/redock_report.json']
        missing=[x for x in required if not (root/x).is_file()]
        return bool(missing),'missing='+','.join(missing)
    def metric_only():
        """Deliberately narrow: compare only reported AUC against score-derived AUC."""
        p=root/cfg['output']/'actual_vina_scores.csv'
        data=json.loads((root/cfg['output']/'analysis.json').read_text())
        _, rows=csv_load(p)
        act=[];dec=[]
        for row in rows:
            x=float(row['score_kcal_mol'])
            if not math.isfinite(x):raise ValueError('nonfinite')
            if row['class']=='active':act.append(x)
            elif row['class']=='decoy':dec.append(x)
            else:raise ValueError('unexpected class')
        measured=auc_by_pairwise(act,dec)
        metric='roc_auc_3HMM' if case=='alk5-3hmm' else 'roc_auc'
        reported=float(data[metric])
        return abs(measured-reported)>0.00011,f'computed={measured:.4f}, reported={reported:.4f}'
    def snapshot_only():
        changed=digest_comparison({p:digest(root/p) for p in ref['artifact_sha256'] if (root/p).is_file()},ref['artifact_sha256'])
        return bool(changed),f'changed_files={len(changed)}'
    def eave_without_snapshot():
        a=Audit(root,case,0).run()
        claim= a.get('claims',{}).get('experimental_potency',{}) if scenario=='unsupported_potency_claim' else None
        detected=(a['audit_status']!='verified' or (claim is not None and claim.get('status')!='supported'))
        return detected, f"audit={a['audit_status']}; claim={claim.get('status') if claim is not None else 'not requested'}"
    def eave_full():
        claim='experimental_potency' if scenario=='unsupported_potency_claim' else None
        rc,out=gate(root,case,bootstrap=0,claim=claim,write=False)
        return rc!=0, f"verdict={out['gate_verdict']}; changes={len(out['baseline_changed_artifacts'])}"
    impl={'presence_only':presence_only,'metric_only':metric_only,
          'snapshot_only':snapshot_only,'eave_without_snapshot':eave_without_snapshot,
          'eave_full':eave_full}
    result={}
    for name in POLICIES:
        start=time.perf_counter()
        try:
            detected, detail=impl[name]()
            error=None
        except (OSError,KeyError,ValueError,TypeError,ZeroDivisionError) as exc:
            detected=True
            detail='invalid/unreadable evidence'
            error=f'{type(exc).__name__}: {str(exc)[:130]}'
        ms=(time.perf_counter()-start)*1000
        result[name]={'flagged':bool(detected),'runtime_ms':round(ms,3),'detail':detail,'runtime_error':error}
    return result


def summarise(trials:list[dict]):
    rows=[]
    for name in POLICIES:
        tp=fn=fp=tn=0
        runtimes=[]
        for t in trials:
            truth=t['is_fault']; flag=t['policies'][name]['flagged']
            runtimes.append(t['policies'][name]['runtime_ms'])
            if truth and flag:tp+=1
            elif truth and not flag:fn+=1
            elif not truth and flag:fp+=1
            else:tn+=1
        rows.append({'policy':name,'tp':tp,'fn':fn,'fp':fp,'tn':tn,
                     'sensitivity_on_injected_faults':round(tp/(tp+fn),4) if tp+fn else None,
                     'false_positive_rate_on_controls':round(fp/(fp+tn),4) if fp+tn else None,
                     'median_audit_ms':round(statistics.median(runtimes),3)})
    return rows


def run_suite(project_root:Path, output:Path, cases:list[str]):
    project_root=project_root.resolve(); output=output.resolve()
    for c in cases:
        if c not in CASES:raise ValueError(f'Unsupported case {c}')
    # Explicitly refuse writing results into a source evidence directory.
    for c in cases:
        cfg=CASES[c]
        for rel in (cfg['output'],cfg['sources'],cfg['panel']):
            p=(project_root/rel).resolve()
            if output==p or p in output.parents:raise ValueError('Output must not be under frozen experimental files')
    output.mkdir(parents=True,exist_ok=True)
    records=[]; baseline_summaries=[]
    for case in cases:
        with tempfile.TemporaryDirectory(prefix='aipharma_phase16_clean_') as tmp:
            baseline=stage_case(project_root,case,Path(tmp))
            baseline_summaries.append({'case':case,'original_verified':True,
               'auc':Audit(Path(tmp),case,0).run()['metrics']['auc_recomputed'],
               'baseline_sha256':digest(project_root/'phase14_eave_baselines'/f'{case}.json')})
        for scenario,is_fault in SCENARIOS:
            with tempfile.TemporaryDirectory(prefix='aipharma_phase16_fault_') as td:
                stage_case(project_root,case,Path(td))
                mutate(Path(td),case,scenario)
                policies=evaluate(Path(td),case,scenario,baseline)
                records.append({'case':case,'scenario':scenario,'is_fault':is_fault,'policies':policies})
            print(f'{case:12s} {scenario:27s} '+ ' '.join(f'{p}={int(policies[p]["flagged"])}' for p in POLICIES),flush=True)
    summaries=summarise(records)
    report={'schema':'aipharma-phase16-controlled-audit-comparison-v1',
            'generated_utc':datetime.now(timezone.utc).isoformat(),
            'cases':cases,'n_scenarios_per_case':len(SCENARIOS),
            'n_fault_scenarios_per_case':sum(f for _,f in SCENARIOS),
            'n_negative_controls_per_case':sum(not f for _,f in SCENARIOS),
            'data_basis':'Mutations of temporary COPIES of existing local benchmark artifacts; no docking performed',
            'comparison_scope':'Five local control policies; Icolos and DDA NOT executed or measured',
            'baseline_trust_note':'Reference manifests are local unsigned files; hash comparisons are not cryptographic attestation',
            'timing_note':'Single-process wall-time for local policy calls, not system-level throughput nor cross-tool timing benchmark',
            'labeling_note':'Injected fault/clean labels defined by this protocol; limited faults and negative controls',
            'baseline_summary':baseline_summaries,'policy_summary':summaries,'trials':records}
    write_json(output/'benchmark_results.json',report)
    with (output/'policy_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    with (output/'trial_matrix.csv').open('w',newline='') as f:
        fields=['case','scenario','ground_truth_fault',*POLICIES]
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for rec in records:
            writer.writerow({'case':rec['case'],'scenario':rec['scenario'],'ground_truth_fault':int(rec['is_fault']),
                            **{p:int(rec['policies'][p]['flagged']) for p in POLICIES}})
    print('\nSUMMARY')
    for s in summaries:
        print(f"{s['policy']:24s} TP={s['tp']} FN={s['fn']} FP={s['fp']} TN={s['tn']} median_ms={s['median_audit_ms']:.2f}")
    print('Report:',output/'benchmark_results.json')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd(),help='AIPHARMA repo with completed Phase 14 baselines')
    p.add_argument('--output',type=Path,default=None,help='Report directory (default reports/eave/phase16_comparison)')
    p.add_argument('--case',choices=('all',*sorted(CASES)),default='alk5-1vjy')
    args=p.parse_args(argv)
    root=args.root.resolve()
    output=args.output or root/'reports/eave/phase16_comparison'
    cases=list(CASES) if args.case=='all' else [args.case]
    try:run_suite(root,output,cases)
    except (ValueError,RuntimeError,OSError) as e:
        print('ERROR:',e,file=sys.stderr)
        return 2
    return 0

if __name__=='__main__':raise SystemExit(main())
