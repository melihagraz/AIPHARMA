#!/usr/bin/env python3
"""AIPHARMA Phase 9: paired stratified bootstrap of existing Phase 8 predictions.

No ADMET inference, no model training, no new datasets. Analyses are retrospective
TDC same-data-family reference comparisons, NOT independent external validation.
Paired resampling retains identical test molecules for both scores; class-stratified
resampling holds observed label prevalence fixed. CIs treat rows as independent and
thus do not account for chemical scaffold/source clustering or model-training reuse.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

ENDPOINTS = ('DILI', 'hERG', 'AMES')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def metric_values(y, a, b):
    vals = {}
    for tag, x in (('ADMET_AI', a), ('RDKit_baseline', b)):
        vals[tag] = {
            'roc_auc': float(roc_auc_score(y, x)),
            'average_precision': float(average_precision_score(y, x)),
            'brier': float(brier_score_loss(y, x)),
        }
    vals['difference_ADMET_minus_baseline'] = {
        'roc_auc': vals['ADMET_AI']['roc_auc'] - vals['RDKit_baseline']['roc_auc'],
        'average_precision': vals['ADMET_AI']['average_precision'] - vals['RDKit_baseline']['average_precision'],
        'brier_improvement_baseline_minus_ADMET': vals['RDKit_baseline']['brier'] - vals['ADMET_AI']['brier'],
    }
    return vals


def ci(a):
    z = np.asarray(a, dtype=float)
    return [float(np.quantile(z, .025)), float(np.quantile(z, .975))]


def analyze(folder: Path, n_boot: int, seed: int):
    p = folder / 'comparison_predictions.csv'
    e = folder / 'evaluation.json'
    if not p.exists() or not e.exists():
        raise FileNotFoundError(f'{folder}: phase8 comparison_predictions.csv and evaluation.json required')
    recorded = json.loads(e.read_text(encoding='utf-8'))
    got_hash = sha256(p)
    if got_hash != recorded.get('prediction_sha256'):
        raise RuntimeError(f'{folder.name}: prediction SHA mismatch. File has changed; aborting.')
    with p.open(newline='', encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))
    required = {'canonical_smiles', 'label', 'admet_ai_score', 'rdkit_8_descriptor_logistic_score'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f'{folder.name}: prediction CSV lacks required columns')
    y = np.array([int(r['label']) for r in rows], dtype=int)
    a = np.array([float(r['admet_ai_score']) for r in rows], dtype=float)
    b = np.array([float(r['rdkit_8_descriptor_logistic_score']) for r in rows], dtype=float)
    ids = [r['canonical_smiles'] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError(f'{folder.name}: duplicate canonical SMILES in test comparison')
    if not set(y.tolist()) <= {0, 1} or min(int(y.sum()), int((y == 0).sum())) < 5:
        raise ValueError(f'{folder.name}: insufficient valid binary class labels')
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        raise ValueError(f'{folder.name}: non-finite scores')
    if not (np.all((a >= 0) & (a <= 1)) and np.all((b >= 0) & (b <= 1))):
        raise ValueError(f'{folder.name}: Brier score requires probabilities in [0,1]')
    if len(y) != recorded.get('test_n'):
        raise RuntimeError(f'{folder.name}: test count differs from archived evaluation')
    point = metric_values(y, a, b)
    for name, phase8 in [('ADMET_AI', 'ADMET_AI'), ('RDKit_baseline', 'RDKit_8_feature_logistic')]:
        for k in ('roc_auc', 'average_precision', 'brier'):
            if not math.isclose(point[name][k], recorded['models'][phase8][k], abs_tol=1e-9, rel_tol=1e-9):
                raise RuntimeError(f'{folder.name}: recomputed {name} {k} differs from Phase 8')
    negative_idx, positive_idx = np.flatnonzero(y == 0), np.flatnonzero(y == 1)
    rng = np.random.default_rng(seed)
    sample_metrics = {t: {k: [] for k in values} for t, values in point.items()}
    for _ in range(n_boot):
        take = np.concatenate((rng.choice(negative_idx, size=len(negative_idx), replace=True),
                               rng.choice(positive_idx, size=len(positive_idx), replace=True)))
        x = metric_values(y[take], a[take], b[take])
        for t, metrics in x.items():
            for k, v in metrics.items():
                sample_metrics[t][k].append(v)
    intervals = {t: {k: ci(v) for k, v in vals.items()} for t, vals in sample_metrics.items()}
    return {
        'n': len(y), 'positive': int(y.sum()), 'negative': int((y == 0).sum()),
        'prediction_sha256': got_hash,
        'point_estimates': point,
        'bootstrap95_percentile_intervals': intervals,
        'bootstrap_replicates': n_boot,
        'resampling': 'same indices for both models; within-class row resampling',
        'limitations': [
            'ADMET-AI used TDC during training/development: NOT an independent external validation.',
            'CIs are retrospective conditional intervals, not independent/model-development-adjusted CIs.',
            'Rows are resampled independently; chemical scaffold/assay/source correlation is NOT accounted for.',
            'RDKit comparator is a deliberately basic 8-feature logistic model, not a state-of-the-art challenger.',
            'Brier score is uncalibrated; baseline used class_weight=balanced; do not conclude clinically validated calibration.',
            'Point estimates and bootstrap intervals measure test-set sampling uncertainty only.',
        ]
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, default=Path('.'))
    ap.add_argument('--bootstrap', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20261112)
    ap.add_argument('--only', nargs='*', choices=ENDPOINTS)
    args = ap.parse_args()
    if args.bootstrap < 200:
        ap.error('--bootstrap must be at least 200')
    endpoints = tuple(dict.fromkeys(args.only)) if args.only else ENDPOINTS
    root = args.repo.resolve() / 'reports' / 'validation_phase1' / 'phase8_admet'
    report = {
        'kind': 'AIPHARMA_Phase9_paired_stratified_bootstrap_TDC_reference_not_external',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'seed': args.seed, 'bootstrap_replicates': args.bootstrap,
        'endpoints': {},
    }
    for i, endpoint in enumerate(endpoints):
        report['endpoints'][endpoint] = analyze(root / endpoint, args.bootstrap, args.seed + i)
        r = report['endpoints'][endpoint]
        d = r['point_estimates']['difference_ADMET_minus_baseline']
        interval = r['bootstrap95_percentile_intervals']['difference_ADMET_minus_baseline']['roc_auc']
        print(f"{endpoint:5s} n={r['n']:4d} | ADMET AUC={r['point_estimates']['ADMET_AI']['roc_auc']:.4f} "
              f"RDKit AUC={r['point_estimates']['RDKit_baseline']['roc_auc']:.4f} "
              f"ΔAUC={d['roc_auc']:+.4f} (paired 95% CI [{interval[0]:+.4f}, {interval[1]:+.4f}])", flush=True)
    out = root / 'phase9_bootstrap_summary.json'
    out.write_text(json.dumps(report, indent=2, allow_nan=False, ensure_ascii=False) + '\n', encoding='utf-8')
    print('SAVED:', out)
    print('WARNING: TDC reference reuse; NOT independent validation or AIPHARMA model superiority.')

if __name__ == '__main__':
    main()
