#!/usr/bin/env python3
"""AIPHARMA Phase 8: auditable TDC ADMET reference benchmarking.

IMPORTANT: ADMET-AI 2 is trained on TDC data. TDC testing does NOT establish
independent external validation of ADMET-AI, even with official scaffold splits.

Modes:
  audit    Download official TDC ADMET_Group partitions, save SHA256 & overlap audit.
  evaluate Score only previously audited test inputs with installed ADMET-AI and a
           fixed RDKit descriptor baseline (descriptive comparison, NOT external).
  all      Run audit then evaluate.

Use --local-data-dir only for OFFLINE TESTING or pre-existing trusted TDC CSVs.
Expected <dir>/<ENDPOINT>/train_val.csv and test.csv, columns Drug,Y.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ENDPOINTS = ('DILI', 'hERG', 'AMES')
REQUIRED_COLUMNS = {'Drug', 'Y'}
FEATURE_NAMES = ('MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotB', 'Heavy', 'FractionCSP3')
PROVENANCE_WARNING = (
    'ADMET-AI v2 was trained on Therapeutics Data Commons (TDC) datasets. '
    'The TDC benchmark test partition cannot be assumed external to its training/model development. '
    'These metrics are retrospective, same-data-family reference evaluations, NOT independent validation, '
    'NOT evidence of clinical or experimental compound safety, and NOT a fair head-to-head superiority claim.'
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')


def versions() -> dict:
    result = {'python': platform.python_version(), 'rdkit': rdBase.rdkitVersion,
              'numpy': np.__version__, 'pandas': pd.__version__}
    for p in ('scikit-learn', 'PyTDC', 'admet-ai'):
        try:
            result[p] = importlib.metadata.version(p)
        except importlib.metadata.PackageNotFoundError:
            result[p] = 'not_installed'
    return result


def canonical(s: str) -> str | None:
    if not isinstance(s, str) or not s.strip():
        return None
    m = Chem.MolFromSmiles(s)
    if m is None:
        return None
    return Chem.MolToSmiles(m, isomericSmiles=True)


def prep_frame(frame: pd.DataFrame, label: str):
    if not REQUIRED_COLUMNS <= set(frame.columns):
        raise ValueError(f'{label}: expected Drug,Y; got {list(frame.columns)}')
    df = frame[['Drug', 'Y']].copy()
    df['Drug'] = df['Drug'].astype(str)
    df['Y'] = pd.to_numeric(df['Y'], errors='coerce')
    start = len(df)
    df = df[df['Y'].isin([0, 1])].copy()
    invalid_labels = start - len(df)
    # Enforce valid structures; preserve original SMILES alongside canonical.
    df['canonical_smiles'] = [canonical(x) for x in df.Drug]
    invalid_smiles = int(df.canonical_smiles.isna().sum())
    df = df[df.canonical_smiles.notna()].copy()
    df['Y'] = df['Y'].astype(int)
    conflicting = set(df.groupby('canonical_smiles')['Y'].nunique().loc[lambda x: x > 1].index)
    n_conflict_rows = int(df.canonical_smiles.isin(conflicting).sum())
    # Conservative: exclude both conflicting classes for repeated molecule identity.
    df = df[~df.canonical_smiles.isin(conflicting)].drop_duplicates('canonical_smiles', keep='first').copy()
    df.sort_values(['canonical_smiles'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    stats = {'n_raw': start, 'n_bad_or_missing_label': int(invalid_labels),
             'n_invalid_smiles': int(invalid_smiles), 'n_rows_with_conflicting_identity_labels': n_conflict_rows,
             'n_deduplicated_valid': len(df), 'positive': int(df.Y.sum()),
             'negative': int((df.Y == 0).sum())}
    return df, stats


def get_tdc_group(repo: Path):
    try:
        from tdc.benchmark_group import admet_group
    except ImportError as exc:
        raise RuntimeError('PyTDC not installed. For data acquisition use a separate '
                           'environment if possible: python -m pip install PyTDC. '
                           'Then rerun --mode audit. The --mode evaluate command needs '
                           'the AIPHARMA venv with ADMET-AI.') from exc
    cache = repo / 'reports' / 'validation_phase1' / 'phase8_admet' / '_tdc_cache'
    cache.mkdir(parents=True, exist_ok=True)
    return admet_group(path=str(cache))


def get_partitions(group, endpoint):
    choices = [endpoint, endpoint.lower(), endpoint.upper()]
    errors = []
    for name in choices:
        try:
            result = group.get(name)
            if 'train_val' not in result or 'test' not in result:
                raise RuntimeError(f'Unexpected benchmark keys: {list(result)}')
            return result['train_val'], result['test'], name
        except (ValueError, KeyError, RuntimeError) as exc:
            errors.append(f'{name}: {exc}')
    raise RuntimeError(f'Cannot retrieve official TDC ADMET benchmark for {endpoint}: {errors}')


def audit(args, root: Path):
    group = None if args.local_data_dir else get_tdc_group(args.repo)
    total = {'kind': 'TDC_ADMET_reference_data_audit_NOT_external_validation',
             'status': 'AUDITED_NO_MODEL_SCORING', 'generated_utc': datetime.now(timezone.utc).isoformat(),
             'dataset_source': str(args.local_data_dir) if args.local_data_dir else 'official_PyTDC_ADMET_Group',
             'versions': versions(), 'warnings': [PROVENANCE_WARNING,
                 'Cross-dataset model-training identities and model-selection exposure remain unresolved.',
                 'Exact canonical-SMILES overlap is not a salt/tautomer/stereoisomer or scaffold independence audit.'],
             'endpoints': {}}
    for endpoint in args.endpoints:
        folder = root / endpoint
        folder.mkdir(parents=True, exist_ok=True)
        if args.local_data_dir:
            src = args.local_data_dir / endpoint
            train = pd.read_csv(src / 'train_val.csv')
            test = pd.read_csv(src / 'test.csv')
            used_name = f'local:{endpoint}'
        else:
            train, test, used_name = get_partitions(group, endpoint)
        # Preserve EXACT input data as a separate signed-by-digest snapshot.
        train_path, test_path = folder / 'train_val_original.csv', folder / 'test_original.csv'
        train.to_csv(train_path, index=False)
        test.to_csv(test_path, index=False)
        a, sa = prep_frame(train, f'{endpoint}/train_val')
        b, sb = prep_frame(test, f'{endpoint}/test')
        overlap = set(a.canonical_smiles) & set(b.canonical_smiles)
        b_clean = b[~b.canonical_smiles.isin(overlap)].copy()
        a.to_csv(folder / 'train_val_curated.csv', index=False)
        b.to_csv(folder / 'test_curated_before_overlap_filter.csv', index=False)
        b_clean.to_csv(folder / 'test_curated.csv', index=False)
        summary = {'tdc_benchmark_name': used_name,
                   'train_val': sa, 'test': sb, 'exact_canonical_train_test_overlap_n': len(overlap),
                   'test_after_exact_overlap_screen': {'n': len(b_clean), 'positive': int(b_clean.Y.sum()),
                      'negative': int((b_clean.Y == 0).sum())},
                   'original_sha256': {'train_val': sha256_file(train_path), 'test': sha256_file(test_path)},
                   'curated_sha256': {'train_val': sha256_file(folder/'train_val_curated.csv'),
                                      'test': sha256_file(folder/'test_curated.csv')},
                   'caveat': 'Official TDC test is NOT independently external to ADMET-AI training.'}
        write_json(folder / 'audit.json', summary)
        total['endpoints'][endpoint] = summary
        print(f'{endpoint:5s} train={len(a):5d} test={len(b):5d} '
              f'overlap={len(overlap):4d} eligible={len(b_clean):5d} '
              f'positive={int(b_clean.Y.sum()):5d} negative={int((b_clean.Y==0).sum()):5d}')
    write_json(root / 'audit_summary.json', total)
    print('AUDIT OUTPUT:', root / 'audit_summary.json')
    print('NOT EXTERNAL VALIDATION. TDC dataset overlaps model training provenance.')


def descriptors(smi: str) -> list[float]:
    m = Chem.MolFromSmiles(smi)
    return [Descriptors.MolWt(m), Descriptors.MolLogP(m), rdMolDescriptors.CalcTPSA(m),
            Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m), Lipinski.NumRotatableBonds(m),
            m.GetNumHeavyAtoms(), rdMolDescriptors.CalcFractionCSP3(m)]


def calculate_metrics(y, scores):
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores, dtype=float)
    good = np.isfinite(scores)
    y, scores = y[good], scores[good]
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    summary = {'n': int(len(y)), 'positive': n_pos, 'negative': n_neg,
               'prevalence': float(n_pos / len(y)) if len(y) else None,
               'average_precision_no_skill_baseline': float(n_pos / len(y)) if len(y) else None}
    if not len(y) or min(n_pos, n_neg) < 5:
        summary['status'] = 'INSUFFICIENT_CLASS_SUPPORT'
        return summary
    summary['status'] = 'DESCRIPTIVE_TDC_REFERENCE_NOT_EXTERNAL'
    summary['roc_auc'] = float(roc_auc_score(y, scores))
    summary['average_precision'] = float(average_precision_score(y, scores))
    if np.min(scores) >= 0 and np.max(scores) <= 1:
        summary['brier'] = float(brier_score_loss(y, scores))
    return summary


def evaluate(args, root: Path):
    # No hidden modifications of inputs: audit digests must still match.
    summary = {'kind': 'ADMET_AI_vs_fixed_RDKit_baseline_TDC_reference_NOT_external_validation',
               'generated_utc': datetime.now(timezone.utc).isoformat(),
               'model_name': 'Installed AIPHARMA AdmetPredictor (ADMET-AI v2)',
               'versions': versions(), 'warnings': [PROVENANCE_WARNING], 'endpoints': {}}
    sys.path.insert(0, str(args.repo.resolve()))
    try:
        from Admet_App.app.predictor import AdmetPredictor
    except Exception as exc:
        raise RuntimeError('Cannot import Admet_App.app.predictor from --repo. '
                           'Run --mode evaluate using the project .venv from ~/AIPHARMA.') from exc
    predictor = AdmetPredictor(num_workers=0)
    for endpoint in args.endpoints:
        folder = root / endpoint
        info_path = folder / 'audit.json'
        if not info_path.exists():
            raise RuntimeError(f'Run --mode audit first; missing {info_path}')
        info = json.loads(info_path.read_text())
        train_path, test_path = folder/'train_val_curated.csv', folder/'test_curated.csv'
        for name, path in [('train_val', train_path), ('test', test_path)]:
            actual, expected = sha256_file(path), info['curated_sha256'][name]
            if actual != expected:
                raise RuntimeError(f'DATA CHANGED: {path} sha256={actual}, expected={expected}')
        train, test = pd.read_csv(train_path), pd.read_csv(test_path)
        if len(test) == 0:
            summary['endpoints'][endpoint] = {'status': 'EMPTY_TEST_AFTER_AUDIT'}
            continue
        if train.Y.nunique() < 2:
            summary['endpoints'][endpoint] = {'status': 'TRAINING_LABELS_SINGLE_CLASS'}
            continue
        model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0,
                                        class_weight='balanced', random_state=20261112))
        X_train = np.array([descriptors(s) for s in train.Drug], dtype=float)
        X_test = np.array([descriptors(s) for s in test.Drug], dtype=float)
        model.fit(X_train, train.Y.to_numpy(dtype=int))
        baseline_scores = model.predict_proba(X_test)[:, 1]
        # Predict in fixed batches without altering model settings or thresholds.
        admet_scores = []
        print(f'{endpoint}: predicting {len(test)} structures; batch={args.batch}...', flush=True)
        for start in range(0, len(test), args.batch):
            batch = test.Drug.iloc[start:start+args.batch].tolist()
            records = predictor.predict(batch)
            if len(records) != len(batch):
                raise RuntimeError(f'{endpoint}: ADMET-AI returned {len(records)}/{len(batch)} predictions')
            for pred in records:
                val = pred.raw.get(endpoint)
                if val is None or not math.isfinite(float(val)):
                    raise RuntimeError(f'ADMET-AI missing numeric endpoint {endpoint}; keys={list(pred.raw)[:30]}')
                admet_scores.append(float(val))
        if not len(admet_scores) == len(test):
            raise RuntimeError(f'Prediction count mismatch: {len(admet_scores)} != {len(test)}')
        outcome = pd.DataFrame({
            'smiles': test.Drug, 'canonical_smiles': test.canonical_smiles,
            'label': test.Y.astype(int), 'admet_ai_score': admet_scores,
            'rdkit_8_descriptor_logistic_score': baseline_scores})
        output = folder / 'comparison_predictions.csv'
        outcome.to_csv(output, index=False)
        results = {'test_n': len(test), 'test_sha256': sha256_file(test_path),
                   'prediction_sha256': sha256_file(output),
                   'models': {'ADMET_AI': calculate_metrics(test.Y, admet_scores),
                              'RDKit_8_feature_logistic': calculate_metrics(test.Y, baseline_scores)},
                   'warning': PROVENANCE_WARNING}
        summary['endpoints'][endpoint] = results
        write_json(folder/'evaluation.json', results)
        m1, m2 = results['models']['ADMET_AI'], results['models']['RDKit_8_feature_logistic']
        print(endpoint, 'ADMET-AI ROC-AUC=', m1.get('roc_auc'),
              'RDKit ROC-AUC=', m2.get('roc_auc'), 'AP baseline=', m1.get('prevalence'))
    write_json(root/'evaluation_summary.json', summary)
    print('EVALUATION OUTPUT:', root/'evaluation_summary.json')
    print('CAUTION: This is TDC reference reuse. No independent validation was established.')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--repo', type=Path, default=Path('.'))
    ap.add_argument('--mode', choices=('audit', 'evaluate', 'all'), default='audit')
    ap.add_argument('--only', nargs='*', choices=ENDPOINTS,
                    help='Select subset of endpoints (e.g. --only DILI hERG)')
    ap.add_argument('--batch', type=int, default=32)
    ap.add_argument('--local-data-dir', type=Path, default=None,
                    help='Offline fixture only. Folder with DILI/hERG/AMES/train_val.csv,test.csv')
    args = ap.parse_args()
    args.repo = args.repo.resolve()
    if args.local_data_dir:
        args.local_data_dir = args.local_data_dir.resolve()
    args.endpoints = tuple(dict.fromkeys(args.only)) if args.only else ENDPOINTS
    if args.batch < 1:
        ap.error('--batch must be >0')
    root = args.repo / 'reports' / 'validation_phase1' / 'phase8_admet'
    root.mkdir(parents=True, exist_ok=True)
    if args.mode in ('audit', 'all'):
        audit(args, root)
    if args.mode in ('evaluate', 'all'):
        evaluate(args, root)


if __name__ == '__main__':
    main()
