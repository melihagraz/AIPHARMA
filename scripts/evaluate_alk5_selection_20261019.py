"""Fixed, exploratory source-family AND scaffold separated selection audit."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', type=Path)
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    protocol = dict(seed=20261019, folds=5,
        purpose='Exploratory public-data comparison; not blind or biological validation',
        exclusions='Relative-stereo new records and monomer 667163 (measurement ceiling)',
        labels='Existing catalytic enzyme labels; interval 5–10 nM used only as active class',
        grouping='Connected components of source-family OR exact recomputed Murcko scaffold',
        model='Morgan radius 2 / 2048 bits; logistic C=1 balanced liblinear; no tuning',
        baselines=['maximum training-active Tanimoto', 'standardized MW/logP/heavy atoms logistic'],
        promotion='All test folds contain both classes; zero source/scaffold overlap; all fold model AUC >0.5; paired component-bootstrap AUC difference lower bound >0 for BOTH baselines',
        caveat='Bootstrap conditions on fixed out-of-fold predictions, not model refitting. No calibrated activity probability. No proof of safety or efficacy.')
    (args.output/'protocol.json').write_text(json.dumps(protocol, indent=2), encoding='utf-8')
    raw = args.input.read_bytes()
    original = json.loads(raw)
    rows = [r for r in original if r.get('pair_eligible', True) and str(r.get('monomer')) != '667163']
    assert len({r['smiles'] for r in rows}) == len(rows), 'Duplicate parents'
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols = [Chem.MolFromSmiles(r['smiles']) for r in rows]
    assert all(m is not None for m in mols)
    scaffolds = [Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) for m in mols]
    fps = [generator.GetFingerprint(m) for m in mols]
    X = np.array([np.asarray(f) for f in fps])
    props = np.array([[r['mw'], r['logp'], r['heavy_atoms']] for r in rows])
    y = np.array([r['label'] for r in rows])
    parent = list(range(len(rows)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}
    for i, r in enumerate(rows):
        assert r['documents'], 'Missing source'
        keys = ['D:'+d for d in r['documents']]
        if scaffolds[i]:
            keys.append('S:'+scaffolds[i])
        for key in keys:
            if key in seen:
                parent[root(i)] = root(seen[key])
            seen[key] = i
    groups = np.array([root(i) for i in range(len(rows))])
    result = dict(input_sha256=hashlib.sha256(raw).hexdigest(), n=len(rows),
        active=int(y.sum()), weak=int((y == 0).sum()), components=len(set(groups)),
        class_components=[len(set(groups[y == c])) for c in [0, 1]], promoted=False)
    (args.output/'dataset.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
    scores = {k: np.full(len(rows), np.nan) for k in ['model', 'similarity', 'properties']}
    folds = []
    if min(result['class_components']) >= 5:
        splits = StratifiedGroupKFold(5, shuffle=True, random_state=20261019).split(X, y, groups)
        for k, (tr, te) in enumerate(splits):
            sources = lambda ix: {d for i in ix for d in rows[i]['documents']}
            overlap = sources(tr) & sources(te)
            scaffold_overlap = (set(scaffolds[i] for i in tr) & set(scaffolds[i] for i in te)) - {''}
            assert not overlap and not scaffold_overlap
            fold = dict(fold=k, active=int(y[te].sum()), weak=int((y[te] == 0).sum()), source_overlap=0, scaffold_overlap=0)
            if len(set(y[tr])) < 2:
                fold['error'] = 'single training class'
                folds.append(fold)
                break
            lr = lambda: LogisticRegression(C=1, class_weight='balanced', solver='liblinear', max_iter=2000, random_state=20261019)
            scores['model'][te] = lr().fit(X[tr], y[tr]).decision_function(X[te])
            scores['properties'][te] = make_pipeline(StandardScaler(), lr()).fit(props[tr], y[tr]).decision_function(props[te])
            act = [fps[i] for i in tr if y[i] == 1]
            scores['similarity'][te] = [max(DataStructs.BulkTanimotoSimilarity(fps[i], act)) for i in te]
            fold['auc'] = {key: float(roc_auc_score(y[te], v[te])) if len(set(y[te])) == 2 else None for key, v in scores.items()}
            folds.append(fold)
    result['folds'] = folds
    if all(np.isfinite(v).all() for v in scores.values()):
        result['pooled_auc'] = {key: float(roc_auc_score(y, v)) for key, v in scores.items()}
        rng = np.random.default_rng(20261019)
        unique = sorted(set(groups))
        differences = {key: [] for key in ['similarity', 'properties']}
        for _ in range(2000):
            ix = np.concatenate([np.flatnonzero(groups == g) for g in rng.choice(unique, len(unique), replace=True)])
            if len(set(y[ix])) < 2:
                continue
            for key in differences:
                differences[key].append(roc_auc_score(y[ix], scores['model'][ix]) - roc_auc_score(y[ix], scores[key][ix]))
        result['conditional_auc_difference_ci95'] = {k: np.quantile(v, [.025, .975]).tolist() for k, v in differences.items()}
        result['promoted'] = bool(all(f.get('auc', {}).get('model') is not None and f['auc']['model'] > .5 for f in folds) and all(v[0] > 0 for v in result['conditional_auc_difference_ci95'].values()))
        (args.output/'predictions.json').write_text(json.dumps([dict(smiles=r['smiles'], label=int(y[i]), component=int(groups[i]), **{k: float(v[i]) for k, v in scores.items()}) for i, r in enumerate(rows)], indent=2))
    result['decision'] = 'Eligible for further independent evaluation, not validated discovery' if result['promoted'] else 'Do not deploy as validated activity ranking'
    (args.output/'metrics.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
