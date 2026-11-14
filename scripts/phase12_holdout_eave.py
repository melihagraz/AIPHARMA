#!/usr/bin/env python3
"""AIPHARMA Phase 12: read-only-to-originals, fault-injection benchmark.

Uses the existing Phase 16 staging and audit policies. New fault mechanisms and
separate semantically benign representation-drift diagnostics. Does not dock,
retrain, or modify frozen experiment inputs. This is NOT a cross-platform test.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import statistics
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = 'aipharma-eave-phase12-scenario-isolation-v1'
SCENARIOS = (
    ('clean', 'clean', 'Original unchanged staged record'),
    ('touch_mtime', 'clean', 'Alter file modification time only, not its bytes'),
    ('unrelated_aux_file', 'clean', 'Create a non-evidence auxiliary log outside the pinned artifacts'),
    ('swap_within_class_scores', 'semantic_fault', 'Swap scores between two identities in the same class; AUC unchanged'),
    ('swap_within_class_source_ids', 'semantic_fault', 'Swap source identifiers in the same class; scores/AUC unchanged'),
    ('drop_row_and_recompute_auc', 'semantic_fault', 'Delete scored molecule and update claimed AUC coherently'),
    ('corrupt_source_id', 'semantic_fault', 'Replace a scored molecule source identifier with a non-existent value'),
    ('zero_modes_with_finite_score', 'semantic_fault', 'Claim zero docking modes for a row with finite score'),
    ('equivalent_score_format', 'representation_drift', 'Rewrite one numeric score with equivalent extra zero'),
    ('reverse_csv_row_order', 'representation_drift', 'Permute row order without changing field values'),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load_phase16(root: Path):
    path = root / 'scripts' / 'phase16_compare_audit.py'
    if not path.is_file():
        raise FileNotFoundError(f'Existing Phase 16 runner missing: {path}')
    sys.path.insert(0, str(root))
    spec = importlib.util.spec_from_file_location('aipharma_original_phase16_for_phase12', path)
    if spec is None or spec.loader is None:
        raise RuntimeError('Could not load the existing Phase 16 runner')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_score_rows(phase16, staged: Path, case: str):
    cfg = phase16.CASES[case]
    path = staged / cfg['output'] / 'actual_vina_scores.csv'
    fields, rows = phase16.csv_load(path)
    if not fields or not rows:
        raise ValueError(f'Missing score fields or rows in {path}')
    required = {'score_kcal_mol', 'source_id', 'class', 'num_modes'}
    if not required.issubset(fields):
        raise ValueError(f'Missing score columns: {sorted(required - set(fields))}')
    return path, fields, rows


def save_score_rows(phase16, path: Path, fields, rows):
    phase16.csv_save(path, fields, rows)


def calculated_auc(phase16, rows):
    a = [float(r['score_kcal_mol']) for r in rows if r['class'] == 'active']
    d = [float(r['score_kcal_mol']) for r in rows if r['class'] == 'decoy']
    return phase16.auc_by_pairwise(a, d)


def first_distinct_pair(rows, key, *, same_class=True):
    for i, a in enumerate(rows):
        for j in range(i+1, len(rows)):
            b = rows[j]
            if same_class and a['class'] != b['class']:
                continue
            if key == 'source_id' and a.get('pair_id') and a.get('pair_id') == b.get('pair_id'):
                continue
            if a[key] != b[key]:
                return i, j
    raise ValueError(f'No two suitable rows have distinct {key}')


def mutate(phase16, staged: Path, case: str, scenario: str) -> dict:
    if scenario == 'clean':
        return {'change': 'none'}
    path, fields, rows = read_score_rows(phase16, staged, case)
    if scenario == 'touch_mtime':
        stat = path.stat()
        os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        return {'change':'timestamp only', 'content_sha256_unchanged':True}
    if scenario == 'unrelated_aux_file':
        p = staged / '_phase12_auxiliary_unpinned' / 'run_note.log'
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('Auxiliary text is outside the frozen reference.\n', encoding='utf8')
        return {'change':'unrelated non-evidence file'}
    if scenario == 'swap_within_class_scores':
        i,j = first_distinct_pair(rows, 'score_kcal_mol')
        old_auc = calculated_auc(phase16, rows)
        rows[i]['score_kcal_mol'],rows[j]['score_kcal_mol'] = rows[j]['score_kcal_mol'],rows[i]['score_kcal_mol']
        if calculated_auc(phase16, rows) != old_auc:
            raise AssertionError('Same-class score swap unexpectedly changed AUC')
        save_score_rows(phase16,path,fields,rows)
        return {'change':'score assigned to wrong molecule identity', 'auc_unchanged':True}
    if scenario == 'swap_within_class_source_ids':
        i,j = first_distinct_pair(rows, 'source_id')
        old_auc = calculated_auc(phase16, rows)
        rows[i]['source_id'],rows[j]['source_id'] = rows[j]['source_id'],rows[i]['source_id']
        if calculated_auc(phase16, rows) != old_auc:
            raise AssertionError('Source ID swap unexpectedly changed AUC')
        save_score_rows(phase16,path,fields,rows)
        return {'change':'wrong mapping from score to source compound', 'auc_unchanged':True}
    if scenario == 'drop_row_and_recompute_auc':
        old_n = len(rows)
        idx = next(i for i,r in enumerate(rows) if r['class']=='active')
        rows.pop(idx)
        new_auc = calculated_auc(phase16, rows)
        save_score_rows(phase16,path,fields,rows)
        analysis_path = staged/phase16.CASES[case]['output']/'analysis.json'
        analysis = json.loads(analysis_path.read_text(encoding='utf8'))
        key = 'roc_auc_3HMM' if case == 'alk5-3hmm' else 'roc_auc'
        if key not in analysis:
            raise KeyError(f'Missing published AUC key {key}')
        analysis[key] = new_auc
        phase16.write_json(analysis_path,analysis)
        return {'change':'one active scored row dropped, reported AUC recalculated',
                'original_n':old_n,'new_n':len(rows),'new_reported_auc':new_auc}
    if scenario == 'corrupt_source_id':
        original = rows[0]['source_id']
        rows[0]['source_id'] = 'PHASE12_INVALID_SOURCE_ID__NOT_IN_PANEL'
        save_score_rows(phase16,path,fields,rows)
        return {'change':'unsupported molecule-to-source mapping','previous_id':original}
    if scenario == 'zero_modes_with_finite_score':
        idx = next(i for i,r in enumerate(rows) if r['score_kcal_mol'].strip() and r['num_modes'].strip() != '0')
        rows[idx]['num_modes']='0'
        save_score_rows(phase16,path,fields,rows)
        return {'change':'scored pose row now declares zero docking modes'}
    if scenario == 'equivalent_score_format':
        idx = next(i for i,r in enumerate(rows)
                   if re.fullmatch(r'[+-]?\d+(?:\.\d+)?', r['score_kcal_mol'].strip()))
        old = rows[idx]['score_kcal_mol']
        updated = old+'0' if '.' in old else old+'.0'
        assert float(updated)==float(old), (old,updated)
        rows[idx]['score_kcal_mol'] = updated
        save_score_rows(phase16,path,fields,rows)
        return {'change':'same numeric value, different CSV bytes','old':old,'new':updated}
    if scenario == 'reverse_csv_row_order':
        save_score_rows(phase16,path,fields,list(reversed(rows)))
        return {'change':'same scored records, reverse CSV order'}
    raise KeyError(f'Unexpected scenario: {scenario}')


def changed_artifacts(stage: Path, ref: dict):
    changes=[]
    for rel,original in ref['artifact_sha256'].items():
        p=stage/rel
        now=sha256(p) if p.is_file() else None
        if now != original:
            changes.append(rel)
    return sorted(changes)


def integrity_table(records, policies):
    out=[]
    for policy in policies:
        for label in ('clean','semantic_fault','representation_drift'):
            subset=[r for r in records if r['category']==label]
            tp=sum(bool(r['policies'][policy]['flagged']) for r in subset)
            err=sum(bool(r['policies'][policy].get('runtime_error')) for r in subset)
            out.append({'policy':policy,'category':label,'n':len(subset),
                        'flagged':tp,'not_flagged':len(subset)-tp,'runtime_errors':err,
                        'flagged_fraction':round(tp/len(subset),4) if subset else None})
    return out


def emit_csv(path, records, policies):
    with path.open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f)
        w.writerow(['case','scenario','category','n_changed_pinned_files','mutation_detail',*policies])
        for rec in records:
            w.writerow([rec['case'],rec['scenario'],rec['category'],
                        len(rec['pinned_changes']),json.dumps(rec['mutation'],sort_keys=True),
                        *[int(rec['policies'][p]['flagged']) for p in policies]])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--output',type=Path)
    p.add_argument('--cases',nargs='+',default=['egfr-1m17','alk5-1vjy','alk5-3hmm'])
    p.add_argument('--mode',choices=['preflight','run'],default='run')
    args=p.parse_args()
    root=args.root.resolve()
    out=(args.output or root/'reports/eave/phase12_holdout').resolve()
    phase16=load_phase16(root)
    cases=args.cases
    for c in cases:
        if c not in phase16.CASES:
            p.error(f'Unknown existing EAVE case: {c}')
    if out == root or root not in out.parents:
        p.error('Output must be within the project report directory')
    if out == root/'reports' or out == root/'reports/eave':
        p.error('Output must be a dedicated Phase 12 subdirectory')
    for c in cases:
        cfg=phase16.CASES[c]
        for rel in (cfg['output'],cfg['sources'],cfg['panel'],cfg['descriptors'],cfg['selection'],cfg['redock']):
            evidence=(root/rel).resolve()
            if out == evidence or evidence in out.parents or out in evidence.parents:
                p.error(f'Output would conflict with evidence: {out} vs {evidence}')
    print('PHASE 12: staged copies only. No original evidence files modified.')
    print('Comparison: five LOCAL Phase 16 policies; no external platform comparison.')
    print('New scenario families: integrity, no-byte-change controls, representation drift; no biological validation.')
    baselines={}
    for c in cases:
        with tempfile.TemporaryDirectory(prefix='aipharma_phase12_preflight_') as tmp:
            ref=phase16.stage_case(root,c,Path(tmp))
            baselines[c]=ref
            print(f'PREFLIGHT OK: {c} pinned_artifacts={len(ref["artifact_sha256"])}')
    if args.mode=='preflight':
        print('No benchmark scenarios executed (preflight only).')
        return 0
    records=[]
    for case in cases:
        for scenario,category,desc in SCENARIOS:
            with tempfile.TemporaryDirectory(prefix='aipharma_phase12_case_') as tmp:
                stage=Path(tmp)
                ref=phase16.stage_case(root,case,stage)
                info=mutate(phase16,stage,case,scenario)
                changed=changed_artifacts(stage,ref)
                if category=='clean' and changed:
                    raise AssertionError(f'Clean control {scenario} changed pinned evidence: {changed}')
                if category != 'clean' and not changed:
                    raise AssertionError(f'Mutation {scenario} did not change any pinned file')
                policies=phase16.evaluate(stage,case,scenario,ref)
                records.append({'case':case,'scenario':scenario,'category':category,
                                'description':desc,'mutation':info,'pinned_changes':changed,'policies':policies})
                print(f'{case:12s} {scenario:31s} {category:22s} ' +
                      ' '.join(f'{k}:{int(policies[k]["flagged"])}' for k in phase16.POLICIES),flush=True)
    # Check that the on-disk original evidence still matches the same trusted hashes.
    for case,ref in baselines.items():
        for rel, expected in ref['artifact_sha256'].items():
            src=root/rel
            if not src.is_file() or sha256(src)!=expected:
                raise RuntimeError(f'ORIGINAL ARTIFACT CHANGED after benchmark: {rel}')
    out.mkdir(parents=True,exist_ok=True)
    policies=list(phase16.POLICIES)
    summary=integrity_table(records,policies)
    report={'schema':SCHEMA,'created_utc':datetime.now(timezone.utc).isoformat(),
            'cases':cases,'scenario_manifest':[{'name':s,'category':c,'description':d} for s,c,d in SCENARIOS],
            'categories':['clean','semantic_fault','representation_drift'],
            'n_trials':len(records),'n_unique_workflows':len(cases),
            'data_scope':'Three previously used EAVE workflows; new local fault mechanisms on staged copies',
            'comparison_scope':'Five local policies, NOT other software products; no docking',
            'integrity_truth':'Pinned artifact byte changes; representation drifts are reported separately from semantic faults',
            'claim_truth':'Not tested. Phase 16 claim checks are not pooled with integrity outcomes.',
            'time_note':'Single-process local call times; not comparable end-to-end pipeline runtimes',
            'dependency_sha256':{'phase16_compare_audit.py':sha256(root/'scripts/phase16_compare_audit.py'),
                                 'eave_audit.py':sha256(root/'scripts/eave_audit.py')},
            'summary':summary,'trials':records,
            'limitations':[
               'The three workflow cases previously appeared in Phase 16; new fault mechanisms alone do not establish generalization to new targets.',
               'The tests were designed with knowledge of EAVE code, not independently preregistered or assessed by third-party auditors.',
               'Snapshot hashes detect changes in pinned bytes; semantically equivalent reserialization also changes hashes.',
               'A fully signed external baseline and competing platforms have not been assessed.',
               'Scenarios are correlated within workflow; naive molecule-level confidence intervals would be invalid.',
               'Legacy comparator operational exceptions are recorded and counted as alerts; inspect runtime_error before claims.',
            ]}
    (out/'phase12_results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf8')
    emit_csv(out/'phase12_trial_matrix.csv',records,policies)
    with (out/'phase12_summary.csv').open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=list(summary[0]))
        writer.writeheader();writer.writerows(summary)
    lines=['# Phase 12 — scenario-isolated EAVE integrity check','',
           '**Experimental scope:** 3 previously used workflow cases; new synthetic mutations on temporary copies.','',
           '**No independent platform benchmark, biological validity, or prospective generalization claim.**','',
           '| Policy | Clean flags | Semantic faults detected | Representation drifts flagged |',
           '|---|---:|---:|---:|']
    for policy in policies:
        lookup={r['category']:r for r in summary if r['policy']==policy}
        lines.append(f"| {policy} | {lookup['clean']['flagged']}/{lookup['clean']['n']} | "
                     f"{lookup['semantic_fault']['flagged']}/{lookup['semantic_fault']['n']} | "
                     f"{lookup['representation_drift']['flagged']}/{lookup['representation_drift']['n']} |")
    lines.extend(['','**Interpretation:** A hash alert on unchanged numerical meaning is a byte-integrity alert, not proof of scientific corruption.',
                  'Claim-level decisions and file-integrity sensitivity must be reported separately.',
                  'For errors thrown by comparator policies see each `runtime_error` in phase12_results.json.',''])
    (out/'phase12_summary.md').write_text('\n'.join(lines),encoding='utf8')
    print('\n'+'\n'.join(lines))
    print('REPORT:',out/'phase12_summary.md')
    return 0

if __name__=='__main__':
    try:
        sys.exit(main())
    except (RuntimeError,ValueError,FileNotFoundError,AssertionError,KeyError,OSError) as ex:
        print(f'PHASE 12 STOPPED WITHOUT CLAIMING A RESULT: {type(ex).__name__}: {ex}',file=sys.stderr)
        sys.exit(2)
