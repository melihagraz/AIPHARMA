#!/usr/bin/env python3
"""Retrospective, fixed-OOF ALK5 enzyme prediction audit.

Recalculates ROC-AUC, average precision, class counts, and a group bootstrap
from previously produced OOF predictions. Does NOT refit, rerun cross-validation,
validate therapeutic efficacy, or establish prospective generalization.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import random
from collections import Counter,defaultdict

import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'reports/alk5_endpoint_archive_20261016'
OUTPUT=ROOT/'reports/validation_20261109_extended/alk5_fixed_oof_analysis.json'


def main():
    data=json.loads((INPUT/'dataset.json').read_text())
    pred=json.loads((INPUT/'predictions.json').read_text())
    assert len(data)==356 and len(pred)==1424
    by_smiles={x['smiles']:x for x in data}
    assert len(by_smiles)==356
    combos={}
    for p in pred:
        key=(p['split'],p['method'])
        combos.setdefault(key,{})[p['smiles']]=p
    counts=Counter(x['label'] for x in data)
    summary={'type':'retrospective_fixed_out_of_fold_prediction_recalculation',
             'warning':'No retraining or prospective validation. Group bootstrap resamples fixed predictions, not model-selection uncertainty.',
             'molecules':len(data),'positive':counts[1], 'weak_inactive':counts[0],
             'majority_class_accuracy':counts[1]/len(data),
             'random_top10_expected_active':10*counts[1]/len(data),
             'scaffolds':len(set(x['scaffold'] for x in data)),
             'weak_scaffolds':len(set(x['scaffold'] for x in data if x['label']==0)),
             'methods':{},
             'inputs_sha256':{f:hashlib.sha256((INPUT/f).read_bytes()).hexdigest() for f in ['dataset.json','predictions.json','metrics.json','protocol.json']}}
    # Pre-group by scaffold so the uncertainty analysis does not assume molecules
    # from the same scaffold are independent.
    grp=defaultdict(list)
    for i,d in enumerate(data):
        grp[d['scaffold']].append(i)
    groups=list(grp.values())
    for (split,method),lookup in sorted(combos.items()):
        assert len(lookup)==356 and set(lookup)==set(by_smiles)
        y=np.asarray([d['label'] for d in data],dtype=int)
        score=np.asarray([float(lookup[d['smiles']]['score']) for d in data],dtype=float)
        assert np.isfinite(score).all()
        assert all(lookup[d['smiles']]['label']==d['label'] for d in data)
        rng=random.Random(20261109)
        boot_auc=[];boot_ap=[]
        for i in range(1500):
            selected=[idx for g in rng.choices(groups,k=len(groups)) for idx in g]
            labels=y[selected]
            if len(set(labels))!=2:continue
            boot_auc.append(float(roc_auc_score(labels,score[selected])))
            boot_ap.append(float(average_precision_score(labels,score[selected])))
        key=f'{split}/{method}'
        summary['methods'][key]={
            'roc_auc':float(roc_auc_score(y,score)),
            'average_precision':float(average_precision_score(y,score)),
            'average_precision_baseline':counts[1]/len(data),
            'group_bootstrap_fixed_prediction_auc_95ci':np.quantile(boot_auc,[.025,.975]).tolist() if boot_auc else None,
            'group_bootstrap_fixed_prediction_ap_95ci':np.quantile(boot_ap,[.025,.975]).tolist() if boot_ap else None,
            'bootstrap_valid_replicates':len(boot_auc),
            'negative_samples':counts[0],
        }
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(json.dumps(summary,indent=2)+'\n')
    print('Fixed OOF ALK5 audit: 356 observations, 11 weak/inactive; outputs:',OUTPUT)
    for k,v in summary['methods'].items():print(k,'AUC',round(v['roc_auc'],4),'AP',round(v['average_precision'],4),'CI',list(map(lambda a:round(a,4),v['group_bootstrap_fixed_prediction_auc_95ci'])))

if __name__=='__main__': main()
