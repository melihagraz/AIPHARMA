"""Leave-one-connected-source/scaffold-component-out audit; no new-data claim."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import rdkit
import sklearn
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def component_groups(rows, scaffolds):
    parent = list(range(len(rows)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}
    for i, r in enumerate(rows):
        assert r['documents'], 'Missing source'
        keys = [('document', d) for d in r['documents']]
        if scaffolds[i]:
            keys.append(('scaffold', scaffolds[i]))
        for key in keys:
            if key in seen:
                parent[root(i)] = root(seen[key])
            seen[key] = i
    return np.array([root(i) for i in range(len(rows))])


def run(input_path, output):
    output.mkdir(parents=True, exist_ok=True)
    data = input_path.read_bytes()
    rows = json.loads(data)
    protocol = {
        'input_sha256': hashlib.sha256(data).hexdigest(),
        'purpose': 'Retrospective source robustness stress test of previously inspected data, NOT fresh external validation',
        'model': 'Unchanged Morgan radius=2, bits=2048; balanced logistic C=1 liblinear, max_iter=2000',
        'splits': 'Leave each connected document-family OR exact Murcko scaffold component out; every row tested once',
        'baselines': ['maximum training-active Tanimoto', 'standardized MW/logP/heavy atoms balanced logistic'],
        'primary_comparison': 'Equal-component mean AUC and paired mean AUC differences across test components containing both labels',
        'uncertainty': '2000 paired bootstrap samples of mixed-class test components; fixed fitted predictions, no refit uncertainty',
        'seed': 20261019,
        'sparse_data_guard': 'Fewer than 10 mixed-class independent components or fewer than 5 weak molecules in any mixed test component prohibits validated deployment; this is a conservative operational guard, not a statistical theorem',
        'other_gates': 'Every mixed component model AUC >0.5; paired mean difference CI lower bound >0 versus both baselines',
        'single_class_tests': 'Retain scores, report counts; AUC undefined, excluded from equal-component mean',
        'no_tuning': True,
        'fresh_external_data': False,
    }
    (output/'protocol.json').write_text(json.dumps(protocol, indent=2))
    assert len({r['smiles'] for r in rows}) == len(rows)
    mols = [Chem.MolFromSmiles(r['smiles']) for r in rows]
    assert all(m is not None for m in mols)
    scaffolds = [Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) for m in mols]
    groups = component_groups(rows, scaffolds)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(m) for m in mols]
    X = np.array([np.asarray(f) for f in fps])
    props = np.array([[r['mw'], r['logp'], r['heavy_atoms']] for r in rows])
    y = np.array([r['label'] for r in rows])
    scores = {k: np.full(len(rows), np.nan) for k in ['model', 'similarity', 'properties']}
    lr = lambda: LogisticRegression(C=1, class_weight='balanced', solver='liblinear', max_iter=2000, random_state=20261019)
    folds = []
    for group in sorted(set(groups)):
        tr, te = np.flatnonzero(groups != group), np.flatnonzero(groups == group)
        sources = lambda ix: {d for i in ix for d in rows[i]['documents']}
        assert not (sources(tr) & sources(te))
        assert not (({scaffolds[i] for i in tr} & {scaffolds[i] for i in te}) - {''})
        fold = dict(component=int(group), test_n=len(te), active=int(y[te].sum()), weak=int((y[te] == 0).sum()), sources=sorted(sources(te)), source_overlap=0, scaffold_overlap=0)
        if len(set(y[tr])) < 2:
            fold['error'] = 'single-class training set'
            folds.append(fold)
            continue
        scores['model'][te] = lr().fit(X[tr], y[tr]).decision_function(X[te])
        scores['properties'][te] = make_pipeline(StandardScaler(), lr()).fit(props[tr], y[tr]).decision_function(props[te])
        act = [fps[i] for i in tr if y[i] == 1]
        scores['similarity'][te] = [max(DataStructs.BulkTanimotoSimilarity(fps[i], act)) for i in te]
        fold['auc'] = {k: float(roc_auc_score(y[te], s[te])) if len(set(y[te])) == 2 else None for k, s in scores.items()}
        folds.append(fold)
    mixed = [f for f in folds if f.get('auc', {}).get('model') is not None]
    result = dict(n=len(rows), active=int(y.sum()), weak=int((y == 0).sum()), components=len(folds), mixed_components=len(mixed), folds=folds, fresh_external_validation=False)
    means = {k: float(np.mean([f['auc'][k] for f in mixed])) for k in scores} if mixed else {}
    result['equal_component_mean_auc'] = means
    rng = np.random.default_rng(20261019)
    intervals = {}
    if mixed:
        for key in ['similarity', 'properties']:
            differences = np.array([f['auc']['model'] - f['auc'][key] for f in mixed])
            boot = differences[rng.integers(0, len(mixed), size=(2000, len(mixed)))].mean(axis=1)
            intervals[key] = dict(mean_difference=float(differences.mean()), conditional_ci95=np.quantile(boot, [.025, .975]).tolist())
    result['paired_component_mean_differences'] = intervals
    reasons = []
    if len(mixed) < 10:
        reasons.append('fewer than 10 mixed-class independent test components')
    if any(f['weak'] < 5 for f in mixed):
        reasons.append('one or more mixed test components have fewer than 5 weak examples')
    if any(f['auc']['model'] <= .5 for f in mixed):
        reasons.append('model fails to exceed chance in one or more mixed test components')
    if any(v['conditional_ci95'][0] <= 0 for v in intervals.values()) or not intervals:
        reasons.append('no demonstrated equal-component advantage against both baselines')
    if not all(np.isfinite(s).all() for s in scores.values()):
        reasons.append('incomplete out-of-component predictions')
    result['validated_deployment_allowed'] = not bool(reasons)
    result['blocking_reasons'] = reasons
    result['versions'] = dict(python=platform.python_version(), rdkit=rdkit.__version__, sklearn=sklearn.__version__, numpy=np.__version__)
    (output/'metrics.json').write_text(json.dumps(result, indent=2))
    (output/'predictions.json').write_text(json.dumps([dict(smiles=r['smiles'], label=int(y[i]), component=int(groups[i]), **{k: float(s[i]) if np.isfinite(s[i]) else None for k, s in scores.items()}) for i, r in enumerate(rows)], indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != 'folds'}, indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    run(args.input, args.output)
