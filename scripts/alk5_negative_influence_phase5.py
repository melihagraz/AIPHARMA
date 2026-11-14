#!/usr/bin/env python3
"""Audit how 11 weak/inactive records influence archived ALK5 ROC-AUC.

Not a retraining exercise, independent validation, or potency experiment.
Uses the fixed predictions already reported in October 2026. It quantifies
class/source imbalance and individual weak-compound leverage on pooled AUC.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import argparse

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'reports/alk5_endpoint_archive_20261016/dataset.json'
PREDS=ROOT/'reports/alk5_endpoint_archive_20261016/predictions.json'
DEFAULT=ROOT/'reports/phase5/alk5_negative_influence.json'


def file_hash(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ranking_share_of_one_negative(positives:np.ndarray, negative:float)->float:
    """Fraction of correctly ordered positive-negative pairs, ties half."""
    return float((np.count_nonzero(positives>negative)+.5*np.count_nonzero(positives==negative))/len(positives))


def analyze(data_path=DATA, pred_path=PREDS):
    rows=json.loads(Path(data_path).read_text())
    predictions=json.loads(Path(pred_path).read_text())
    smiles=[r['smiles'] for r in rows]
    if len(smiles)!=len(set(smiles)):
        raise ValueError('Duplicate input SMILES; score alignment would be ambiguous')
    labels=np.asarray([int(r['label']) for r in rows])
    if not set(labels).issubset({0,1}) or np.count_nonzero(labels==0)<2:
        raise ValueError('Expected positive and multiple weak/inactive records')
    negatives=np.flatnonzero(labels==0)
    positives=np.flatnonzero(labels==1)
    document_counts=Counter(doc for row in rows for doc in row['documents'])
    weak_documents=Counter(doc for row in rows if row['label']==0 for doc in row['documents'])
    pred_dict=defaultdict(dict)
    for row in predictions:
        key=(row['split'],row['method'])
        if row['smiles'] in pred_dict[key]:
            raise ValueError(f'Duplicate predictions in {key}: {row["smiles"]}')
        pred_dict[key][row['smiles']]=float(row['score'])
    results={}
    for (split,method), scored in sorted(pred_dict.items()):
        if set(scored)!=set(smiles):
            raise ValueError(f'Incomplete/mismatched score record for {split}/{method}')
        scores=np.asarray([scored[s] for s in smiles],dtype=float)
        if not np.isfinite(scores).all():
            raise ValueError('Non-finite archived scores')
        roc=float(roc_auc_score(labels,scores))
        ap=float(average_precision_score(labels,scores))
        weak=[]
        for idx in negatives:
            s=ranking_share_of_one_negative(scores[positives],scores[idx])
            keep=np.ones(len(rows),dtype=bool);keep[idx]=False
            auc_minus=float(roc_auc_score(labels[keep],scores[keep]))
            weak.append({'smiles':smiles[idx], 'documents': rows[idx]['documents'],
                         'scaffold':rows[idx]['scaffold'],
                         'positive_pair_ranking_share':s,
                         'roc_auc_after_removal':auc_minus,
                         'roc_auc_change_after_removal':float(auc_minus-roc)})
        expected=float(np.mean([r['positive_pair_ranking_share'] for r in weak]))
        if not np.isclose(roc,expected,atol=1e-12):
            raise AssertionError('Pairwise decomposition does not reconcile AUC')
        loo=np.asarray([r['roc_auc_after_removal'] for r in weak])
        results[f'{split}/{method}']={
            'roc_auc':roc,'average_precision':ap,
            'active_only_average_precision_baseline':float(len(positives)/len(rows)),
            'average_precision_above_active_only_baseline':float(ap-len(positives)/len(rows)),
            'negative_leave_one_out_auc_range':[float(loo.min()),float(loo.max())],
            'weak_compound_contributions':weak,
            'interpretation':'One weak/inactive record contributes 1/11 of the pooled AUC pairwise average; this is not a new out-of-sample estimate.'}
    return {
        'name':'Phase 5 archived ALK5 weak-negative leverage audit',
        'evidence_category':'offline reanalysis of fixed out-of-fold archived predictions',
        'independent_external_test':False,
        'retrained_models':False,
        'input_sha256':{'dataset.json':file_hash(Path(data_path)),'predictions.json':file_hash(Path(pred_path))},
        'counts':{'total':len(rows),'active':len(positives),'weak_or_inactive':len(negatives),
                  'active_fraction':float(len(positives)/len(rows)),
                  'distinct_scaffolds':len(set(r['scaffold'] for r in rows)),
                  'distinct_documents':len(document_counts),
                  'documents_containing_weak':len(weak_documents),
                  'weak_molecules_from_most_common_document':max(weak_documents.values()),
                  'most_common_weak_document':weak_documents.most_common(1)[0][0]},
        'cautions':['Only 11 weak/inactive molecules; each is highly influential.',
                    'Some documents contain most weak examples; scaffold and source confounding may remain.',
                    'AP baseline equals archive active prevalence, not a deployable threshold-dependent baseline.',
                    'Leave-one-negative-out is a retrospective influence analysis; neither retraining nor external validation.',
                    'Endpoint classification is annotation-derived; not confirmed by new wet-lab assays.'],
        'results':results}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=DEFAULT)
    args=parser.parse_args()
    report=analyze()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print('Archived AKL5 weak leverage:',json.dumps({
      'counts':report['counts'],
      'summary':{k:{name: v[name] for name in ('roc_auc','average_precision','average_precision_above_active_only_baseline','negative_leave_one_out_auc_range')} for k,v in report['results'].items()}},indent=2))

if __name__=='__main__':
    main()
