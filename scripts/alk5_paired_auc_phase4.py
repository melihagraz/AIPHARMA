#!/usr/bin/env python3
"""Paired scaffold-cluster bootstrap for fixed archived ALK5 OOF predictions.

This is *not* retraining, a test-set analysis, or model-selection uncertainty.
The purpose is to quantify how much the DIFFERENCE in AUC changes when
resampling the same scaffolds for both retrospective score vectors.
"""
from __future__ import annotations

import json
from pathlib import Path
from collections import defaultdict
import numpy as np
from sklearn.metrics import roc_auc_score

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'reports/alk5_endpoint_archive_20261016/dataset.json'
PREDS=ROOT/'reports/alk5_endpoint_archive_20261016/predictions.json'
OUTPUT=ROOT/'reports/phase4/alk5_paired_auc_bootstrap.json'


def analyze(seed:int=20261109, n_boot:int=1500):
    dataset=json.loads(DATA.read_text())
    preds=json.loads(PREDS.read_text())
    if len(dataset)!=356: raise ValueError('Unexpected reference sample count')
    if len({d['smiles'] for d in dataset})!=len(dataset):
        raise ValueError('Ambiguous duplicate SMILES in archived dataset')
    y=np.array([int(d['label']) for d in dataset], dtype=int)
    if int(y.sum()) != 345 or len(y)-int(y.sum()) != 11:
        raise ValueError('Archive labels changed')
    groups=defaultdict(list)
    for idx, row in enumerate(dataset):
        groups[row['scaffold']].append(idx)
    groups=[np.asarray(idx,dtype=int) for idx in groups.values()]
    lookup=defaultdict(dict)
    for row in preds:
        lookup[(row['split'], row['method'])][row['smiles']]=row
    rng=np.random.default_rng(seed)
    index_resamples=[]
    for _ in range(n_boot):
        selected=np.concatenate([groups[k] for k in rng.integers(0,len(groups),size=len(groups))])
        if len(set(y[selected]))==2:
            index_resamples.append(selected)
    results={}
    for split in ['document_and_scaffold','scaffold']:
        values=[]
        for method in ['logistic','similarity']:
            p=lookup[(split,method)]
            if len(p)!=len(dataset): raise ValueError('Missing prediction vector')
            score=np.asarray([float(p[row['smiles']]['score']) for row in dataset])
            if not np.isfinite(score).all():raise ValueError('Nonfinite score')
            values.append(score)
        a,b=values
        point=float(roc_auc_score(y,a)-roc_auc_score(y,b))
        deltas=np.array([roc_auc_score(y[ind],a[ind])-roc_auc_score(y[ind],b[ind]) for ind in index_resamples])
        lo,hi=np.quantile(deltas,[.025,.975])
        results[split]={
            'logistic_auc':float(roc_auc_score(y,a)),
            'similarity_auc':float(roc_auc_score(y,b)),
            'paired_auc_difference':point,
            'paired_scaffold_bootstrap_95pct_interval':[float(lo),float(hi)],
            'valid_bootstrap_replicates':len(deltas),
            'positive_difference_fraction':float(np.mean(deltas>0)),
        }
    return {'analysis':'paired_scaffold_cluster_bootstrap_of_fixed_archived_oof_predictions',
            'not_retrained_not_external_validation':True,
            'methods_compared':['logistic','similarity'],
            'sample_count':356,'positive_count':345,'weak_inactive_count':11,
            'scaffold_count':len(groups),'seed':seed,'requested_bootstrap_replicates':n_boot,
            'interpretation':'Conditional, exploratory within-archive score comparison. Confidence-like percentile intervals do not incorporate model selection, patient/source sampling, assay labels or domain shift.',
            'splits':results}

if __name__=='__main__':
    r=analyze()
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r['splits'],indent=2))
