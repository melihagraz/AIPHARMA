#!/usr/bin/env python3
"""AIPHARMA phase 2: paired ALK5 docking vs leakage-audited RDKit descriptor baseline.

Uses existing archived 50-active/50-decoy DUD-E TGFR1 files.  No docking,
model fitting to experimental potency, network access, or report-file mutation.

DO NOT interpret this synthetic-decoy benchmark as independent source or
experimental activity validation, or as a benchmark against other software.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

FEATURE_NAMES = ["MolWt", "MolLogP", "TPSA", "HBA", "HBD", "NumRotatableBonds"]


def get_rows(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    needed = {'pair_id', 'class', 'smiles', 'source_id', 'score_kcal_mol', 'error'}
    if not rows or not needed.issubset(rows[0]):
        raise ValueError(f"Missing columns in {path.name}; expected {sorted(needed)}")
    out = {}
    for r in rows:
        cls = r['class'].strip().lower()
        if cls not in {'active', 'decoy'}:
            raise ValueError(f"Unexpected class {cls!r}")
        k = (r['pair_id'].strip(), cls)
        if k in out:
            raise ValueError(f"Duplicate pair/class {k}")
        if r['error'].strip() or not r['score_kcal_mol'].strip():
            raise ValueError(f"Missing or failed docking {k}: {r['error']}")
        score = float(r['score_kcal_mol'])
        if not math.isfinite(score):
            raise ValueError(f"Non-finite docking score {k}")
        mol = Chem.MolFromSmiles(r['smiles'])
        if mol is None:
            raise ValueError(f"Invalid SMILES at {k}")
        r['canonical_smiles'] = Chem.MolToSmiles(mol, isomericSmiles=True)
        r['score'] = score
        r['mol'] = mol
        out[k] = r
    return out


def features(mol):
    return [Descriptors.MolWt(mol), Crippen.MolLogP(mol),
            rdMolDescriptors.CalcTPSA(mol), Lipinski.NumHAcceptors(mol),
            Lipinski.NumHDonors(mol), Lipinski.NumRotatableBonds(mol)]


def make_groups(keys, rows):
    n = len(keys)
    parent = list(range(n))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    def union(a, b):
        x, y = root(a), root(b)
        if x != y:
            parent[y] = x
    first_pair, first_scaffold, first_identity, first_source = {}, {}, {}, {}
    for i, k in enumerate(keys):
        r = rows[k]
        for table, token in ((first_pair, k[0]),
                             (first_identity, r['canonical_smiles']),
                             (first_source, r['source_id'].strip())):
            if token:
                if token in table:
                    union(i, table[token])
                else:
                    table[token] = i
        # Exact Murcko scaffold shared across the two classes or across pairs.
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=r['mol'], includeChirality=False)
        # Do not treat all acyclic structures as one connected scaffold.
        if scaffold:
            if scaffold in first_scaffold:
                union(i, first_scaffold[scaffold])
            else:
                first_scaffold[scaffold] = i
    roots = {root(i) for i in range(n)}
    idmap = {r: j for j, r in enumerate(sorted(roots))}
    return np.array([idmap[root(i)] for i in range(n)])


def ci_and_deltas(y, preds, pair_labels, nboot=2000, seed=20261112):
    pair_ids = sorted(set(pair_labels))
    row_indices = {p: np.where(pair_labels == p)[0] for p in pair_ids}
    rng = np.random.default_rng(seed)
    draws = {name: [] for name in preds}
    diffs = {name: [] for name in preds if name != 'RDKit_descriptor_logistic_OOF'}
    for _ in range(nboot):
        sampled = rng.choice(pair_ids, size=len(pair_ids), replace=True)
        ix = np.concatenate([row_indices[p] for p in sampled])
        vals = {name: roc_auc_score(y[ix], pr[ix]) for name, pr in preds.items()}
        for name, v in vals.items():
            draws[name].append(v)
        for name in diffs:
            diffs[name].append(vals[name] - vals['RDKit_descriptor_logistic_OOF'])
    return ({k: list(map(float, np.quantile(v, [0.025, 0.975]))) for k, v in draws.items()},
            {k: list(map(float, np.quantile(v, [0.025, 0.975]))) for k, v in diffs.items()})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, default=Path('.'))
    ap.add_argument('--out', type=Path, default=None)
    args = ap.parse_args()
    repo = args.repo.expanduser().resolve()
    file1 = repo / 'reports/phase11/alk5_matched50_docking/actual_vina_scores.csv'
    file2 = repo / 'reports/phase12/alk5_3hmm_matched50_docking/actual_vina_scores.csv'
    a, b = get_rows(file1), get_rows(file2)
    if set(a) != set(b):
        raise ValueError('Panels differ in pair_id/class membership; do not compare!')
    keys = sorted(a)
    classes = Counter(k[1] for k in keys)
    pair_counts = Counter(k[0] for k in keys)
    if classes != {'active': 50, 'decoy': 50} or len(pair_counts) != 50 or any(v != 2 for v in pair_counts.values()):
        raise ValueError(f'Unexpected panel: classes={dict(classes)}, pairs={len(pair_counts)}')
    for k in keys:
        if a[k]['canonical_smiles'] != b[k]['canonical_smiles'] or a[k]['source_id'] != b[k]['source_id']:
            raise ValueError(f'Panel identity mismatch: {k}')
    # Fail if the exact same canonical chemical structure carries conflicting labels.
    label_by_identity = {}
    for k in keys:
        identity = a[k]['canonical_smiles']
        if identity in label_by_identity and label_by_identity[identity] != k[1]:
            raise ValueError(f'Conflicting active/decoy labels for the same structure: {identity}')
        label_by_identity[identity] = k[1]
    y = np.array([int(k[1] == 'active') for k in keys])
    X = np.asarray([features(a[k]['mol']) for k in keys], dtype=float)
    groups = make_groups(keys, a)
    group_counts = Counter(groups)
    ng = len(group_counts)
    if ng < 5:
        raise ValueError(f'Only {ng} independent pair/source/scaffold groups; no 5-fold group-safe baseline')
    fold = np.full(len(keys), -1)
    pred = np.full(len(keys), np.nan)
    cv = GroupKFold(n_splits=5)
    for j, (train, test) in enumerate(cv.split(X, y, groups=groups)):
        if len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2:
            raise ValueError(f'Fold {j} has only one class; baseline would be misleading')
        pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, class_weight='balanced', C=1.0))
        pipe.fit(X[train], y[train])
        pred[test] = pipe.predict_proba(X[test])[:, 1]
        fold[test] = j
    assert np.isfinite(pred).all()
    heavy = np.array([a[k]['mol'].GetNumHeavyAtoms() for k in keys])
    methods = {
        'Vina_1VJY': np.array([-a[k]['score'] for k in keys]),
        'Vina_3HMM': np.array([-b[k]['score'] for k in keys]),
        'Vina_1VJY_per_heavy_atom': np.array([-a[k]['score'] for k in keys]) / heavy,
        'Vina_3HMM_per_heavy_atom': np.array([-b[k]['score'] for k in keys]) / heavy,
        'RDKit_descriptor_logistic_OOF': pred,
    }
    pairs = np.array([k[0] for k in keys])
    auc_intervals, delta_intervals = ci_and_deltas(y, methods, pairs)
    result = {
        'interpretation_scope': 'Same DUD-E ALK5 matched active/presumed-decoy panel; NOT independent source or biological validation',
        'n_molecules': len(keys), 'n_pairs': len(pair_counts),
        'n_grouped_components_pair_scaffold_source': ng,
        'largest_component_size_molecules': max(group_counts.values()),
        'cv': '5-fold GroupKFold with pairs, exact normalized identity, source ID, and nonempty exact Murcko scaffold kept together',
        'descriptor_features': FEATURE_NAMES,
        'auc_uncertainty': 'Conditional pair-cluster bootstrap (2000 repeats); OOF predictions fixed; no retraining in bootstrap',
        'methods': {},
    }
    for name, scores in methods.items():
        wins = sum(float(scores[np.where((pairs == p) & (y == 1))[0][0]]) >
                   float(scores[np.where((pairs == p) & (y == 0))[0][0]]) for p in sorted(set(pairs)))
        ties = sum(float(scores[np.where((pairs == p) & (y == 1))[0][0]]) ==
                   float(scores[np.where((pairs == p) & (y == 0))[0][0]]) for p in sorted(set(pairs)))
        entry = {
            'roc_auc': float(roc_auc_score(y, scores)),
            'average_precision': float(average_precision_score(y, scores)),
            'roc_auc_95ci_pair_bootstrap': auc_intervals[name],
            'active_better_within_pair': wins,
            'ties_within_pair': ties,
            'n_pairs': 50,
        }
        if name in delta_intervals:
            entry['auc_minus_RDKit_descriptor_logistic_OOF'] = float(roc_auc_score(y, scores) - roc_auc_score(y, pred))
            entry['auc_delta_95ci_pair_bootstrap'] = delta_intervals[name]
        result['methods'][name] = entry
    out = args.out or repo / 'reports/validation_phase1/alk5_vina_vs_descriptors.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print('RESULT FILE:', out)
    print('PAIRED 50+50 ALK5 PANEL | scaffold/source/pair grouped descriptor OOF')
    print('Number of grouped components:', ng, '| largest component size:', max(group_counts.values()))
    print(f"{'method':37} {'AUC':>6} {'95% CI':>17} {'within-pair wins':>16}")
    for name, entry in result['methods'].items():
        print(f"{name:37} {entry['roc_auc']:.3f} {str([round(t,3) for t in entry['roc_auc_95ci_pair_bootstrap']]):>17} {entry['active_better_within_pair']:>10}/50")
    if max(group_counts.values()) > len(keys) * 0.3:
        print('CAUTION: A group contains >30% of molecules; grouped-CV baseline may be unstable.')
    print('CAUTION: DUD-E decoys are presumed negatives, not experimental inactivity measurements.')
    print('This is neither an independent source-validation result nor a claim of software superiority.')

if __name__ == '__main__':
    main()
