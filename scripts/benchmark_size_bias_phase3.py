#!/usr/bin/env python3
"""Phase 3: exploratory molecule-size sensitivity audit on existing ALK5 DUD-E panel.

No docking, remote downloads or experimental-activity claims. Depends on phase2
script in scripts/benchmark_alk5_phase2.py and its installed RDKit/sklearn.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def cv_predict(X, y, groups):
    out = np.full(len(y), np.nan)
    folds = GroupKFold(n_splits=5)
    for tr, te in folds.split(X, y, groups):
        if len(set(y[tr])) != 2 or len(set(y[te])) != 2:
            raise RuntimeError('One-class fold: cannot compute valid CV predictions')
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=1.0, class_weight='balanced'))
        model.fit(X[tr], y[tr])
        out[te] = model.predict_proba(X[te])[:, 1]
    return out


def auc(y, s):
    return float(roc_auc_score(y, s))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default='.')
    args = ap.parse_args()
    repo = Path(args.repo).resolve()
    sys.path.insert(0, str(repo/'scripts'))
    from benchmark_alk5_phase2 import features, get_rows, make_groups
    a = get_rows(repo/'reports/phase11/alk5_matched50_docking/actual_vina_scores.csv')
    b = get_rows(repo/'reports/phase12/alk5_3hmm_matched50_docking/actual_vina_scores.csv')
    if set(a) != set(b):
        raise RuntimeError('Different pair memberships: cannot compare')
    keys = sorted(a)
    if len(keys) != 100 or len(set(k[0] for k in keys)) != 50:
        raise RuntimeError('Unexpected panel size')
    for k in keys:
        if a[k]['canonical_smiles'] != b[k]['canonical_smiles'] or a[k]['source_id'] != b[k]['source_id']:
            raise RuntimeError('Molecule identity mismatch')
    y = np.array([int(k[1]=='active') for k in keys])
    descriptors = np.array([features(a[k]['mol']) for k in keys], dtype=float)
    heavy = np.array([a[k]['mol'].GetNumHeavyAtoms() for k in keys], dtype=float)
    scores = {
        'Vina_1VJY': np.array([-a[k]['score'] for k in keys]),
        'Vina_3HMM': np.array([-b[k]['score'] for k in keys]),
    }
    groups = make_groups(keys, a)
    if len(set(groups)) < 5:
        raise RuntimeError('Too few leakage groups')
    comparisons = {
        'heavy_atoms_only_OOF': heavy[:, None],
        'RDKit_6_properties_OOF': descriptors,
        'Vina_1VJY_only_OOF': scores['Vina_1VJY'][:, None],
        'Vina_3HMM_only_OOF': scores['Vina_3HMM'][:, None],
        'Vina_1VJY_plus_heavy_OOF': np.column_stack([scores['Vina_1VJY'], heavy]),
        'Vina_3HMM_plus_heavy_OOF': np.column_stack([scores['Vina_3HMM'], heavy]),
        'RDKit_plus_Vina_1VJY_OOF': np.column_stack([descriptors, scores['Vina_1VJY']]),
        'RDKit_plus_Vina_3HMM_OOF': np.column_stack([descriptors, scores['Vina_3HMM']]),
    }
    out = {
        'scope':'Exploratory DUD-E TGFR1 50 active/50 presumed decoy; no independent source or measured biological activity',
        'n_molecules':len(keys), 'n_pairs':50, 'n_group_components':len(set(groups)),
        'model_rule':'All classifiers use the same 5-fold pair/scaffold/source grouped CV, standardization inside training folds',
        'n_heavy_atoms_by_class': {
            cls: {'mean':float(np.mean(heavy[y==label])), 'median':float(np.median(heavy[y==label]))}
            for cls,label in [('active',1),('decoy',0)]
        },
        'raw_heavy_atom_auc':auc(y, heavy),
        'vina_raw_auc': {name:auc(y,s) for name,s in scores.items()},
        'vina_score_heavy_atom_spearman': {
            name: {'rho':float(spearmanr(s,heavy).statistic),'p':float(spearmanr(s,heavy).pvalue)}
            for name,s in scores.items()
        },
        'grouped_cv_auc':{},
        'paired_size_sensitivity':{},
    }
    for name,X in comparisons.items():
        pred = cv_predict(X,y,groups)
        out['grouped_cv_auc'][name] = auc(y, pred)
    index = {k:i for i,k in enumerate(keys)}
    for window in (0,1,2,3):
        selected = []
        for pid in sorted(set(k[0] for k in keys)):
            ia=index[(pid,'active')]
            idc=index[(pid,'decoy')]
            if abs(heavy[ia]-heavy[idc]) <= window:
                selected.append((ia,idc))
        entry = {'pairs':len(selected)}
        for name,s in scores.items():
            entry[name+'_active_better'] = int(sum(bool(s[ia]>s[idc]) for ia,idc in selected))
        out['paired_size_sensitivity'][f'abs_delta_heavy_le_{window}']=entry
    save = repo/'reports/validation_phase1/alk5_size_bias_phase3.json'
    save.parent.mkdir(parents=True,exist_ok=True)
    save.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print('RESULT FILE:',save)
    print('SIZE-ONLY AUC:',round(out['raw_heavy_atom_auc'],4))
    print('MEAN HEAVY ATOMS:',out['n_heavy_atoms_by_class'])
    for k,entry in out['vina_score_heavy_atom_spearman'].items():
        print('SIZE CORRELATION',k,'rho=',round(entry['rho'],3))
    print('\nGROUPED 5-FOLD OOF AUC (small exploratory panel)')
    for k,v in out['grouped_cv_auc'].items():
        print(f'{k:34} {v:.4f}')
    print('\nWITHIN-PAIR SIMILAR-SIZE SENSITIVITY')
    for k,entry in out['paired_size_sensitivity'].items():
        print(k,entry)
    print('\nCAUTION: DUD-E decoys are not experimentally proven inactive; no external validation.')

if __name__=='__main__':
    main()
