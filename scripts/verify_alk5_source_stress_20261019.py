"""Independent pairwise-AUC arithmetic and coverage verification (stdlib only)."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
out = root/'reports/alk5_source_stress_20261019'
source = root/'reports/alk5_selection_20261019/dataset.json'
data = json.loads(source.read_text())
metrics = json.loads((out/'metrics.json').read_text())
protocol = json.loads((out/'protocol.json').read_text())
pred = json.loads((out/'predictions.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest() == protocol['input_sha256']
assert len(pred) == len(data) == metrics['n']
by_smiles = {r['smiles']: r for r in data}
assert len({r['smiles'] for r in pred}) == len(data)
assert set(by_smiles) == {r['smiles'] for r in pred}
seen_sources, seen_scaffolds = {}, {}
for r in pred:
    original = by_smiles[r['smiles']]
    assert r['label'] == original['label']
    for key in original['documents']:
        assert seen_sources.setdefault(key, r['component']) == r['component']
    scaffold = original.get('scaffold')
    if scaffold and scaffold != 'ACYCLIC':
        assert seen_scaffolds.setdefault(scaffold, r['component']) == r['component']
    assert all(r[k] is not None for k in ['model', 'similarity', 'properties'])
methods = ['model', 'similarity', 'properties']
mixed = []
for fold in metrics['folds']:
    rows = [r for r in pred if r['component'] == fold['component']]
    assert len(rows) == fold['test_n']
    assert sum(r['label'] == 0 for r in rows) == fold['weak']
    assert sum(r['label'] == 1 for r in rows) == fold['active']
    positive = [r for r in rows if r['label'] == 1]
    negative = [r for r in rows if r['label'] == 0]
    for key in methods:
        if not positive or not negative:
            assert fold['auc'][key] is None
        else:
            auc = sum((p[key] > n[key]) + .5*(p[key] == n[key]) for p in positive for n in negative)/(len(positive)*len(negative))
            assert abs(auc-fold['auc'][key]) < 1e-12
    if positive and negative:
        mixed.append(fold)
for key in methods:
    mean = sum(f['auc'][key] for f in mixed)/len(mixed)
    assert abs(mean-metrics['equal_component_mean_auc'][key]) < 1e-12
assert metrics['fresh_external_validation'] is False
assert metrics['validated_deployment_allowed'] is False
result = dict(status='passed', verified_rows=len(pred), components=len(metrics['folds']), mixed_components=len(mixed), auc_method='Independent positive-negative pair enumeration', source_scaffold_group_consistency=True)
(out/'verification.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
