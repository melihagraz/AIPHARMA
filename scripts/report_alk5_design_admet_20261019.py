"""Independent paired arithmetic verification and descriptive plots; no ranking."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from statistics import median

ap = argparse.ArgumentParser()
ap.add_argument('report', type=Path)
args = ap.parse_args()
p = args.report
designs = json.loads((p/'input_designs.json').read_text())
protocol = json.loads((p/'protocol.json').read_text())
result = json.loads((p/'comparison.json').read_text())
assert hashlib.sha256((p/'input_designs.json').read_bytes()).hexdigest() == protocol['designs_sha256']
with (p/'all_admet_predictions.csv').open() as f:
    reader = csv.DictReader(f)
    structure_column = reader.fieldnames[0]
    raw = {r[structure_column]: r for r in reader}
assert len(raw) == protocol['unique_structures'] == 46
assert len(result['pairs']) == len(designs['candidates']) == 30
rows, summaries, deviations = [], {}, []
for endpoint in result['endpoints']:
    counts = Counter(lower=0, higher=0, unchanged=0, missing=0)
    deltas = []
    for c, pair in zip(designs['candidates'], result['pairs']):
        assert c['candidate_id'] == pair['candidate_id'] and c['smiles'] == pair['smiles']
        parent = float(raw[c['provenance']['parent_smiles']][endpoint])
        candidate = float(raw[c['smiles']][endpoint])
        delta = candidate-parent
        e = pair['endpoints'][endpoint]
        # Pandas wrote float32 predictions to CSV with decimal rounding;
        # JSON retains their full Python-float conversion. Audit the
        # difference explicitly instead of claiming bitwise CSV equality.
        deviations += [abs(e['parent']-parent), abs(e['candidate']-candidate), abs(e['delta']-delta)]
        assert max(deviations) <= 1e-6
        counts['unchanged' if abs(delta) <= 1e-6 else 'lower' if delta < 0 else 'higher'] += 1
        deltas.append(delta)
        rows.append(dict(candidate_id=c['candidate_id'], parent_id=c['provenance']['parent_id'], endpoint=endpoint, parent=parent, candidate=candidate, delta=delta))
    assert dict(counts) == result['summary'][endpoint]
    summaries[endpoint] = dict(**dict(counts), median_delta=median(deltas), min_delta=min(deltas), max_delta=max(deltas))
with (p/'paired_outputs.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(p/'descriptive_summary.json').write_text(json.dumps(summaries, indent=2))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 3, figsize=(12, 7))
for ax, endpoint in zip(axes.flat, result['endpoints']):
    values = [r['delta'] for r in rows if r['endpoint'] == endpoint]
    ax.scatter(range(1, len(values)+1), values, s=22, color='#2563eb', alpha=.8)
    ax.axhline(0, color='#475569', linewidth=.8)
    ax.set_title(endpoint)
    ax.set_xlabel('Design index (not activity rank)')
    ax.set_ylabel('Design - parent model output')
fig.suptitle('Paired ADMET predictions: 30 designs / 16 parents\nSeparate endpoint scales; model changes do not establish safety', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(p/'paired_admet_changes.png', dpi=160)
fig.savefig(p/'paired_admet_changes.svg')
plt.close(fig)
verification = dict(status='passed', raw_structures=46, paired_designs=30, endpoint_pairs_verified=len(rows), independent_arithmetic='CSV parent/design subtraction and direction-count recalculation', csv_rounding_absolute_tolerance=1e-6, max_csv_json_deviation=max(deviations), rounding_note='Original float32 CSV formatting differs from full-precision JSON; direction counts also independently agree', all_pairs_preserved=True, validated_safety=False)
(p/'verification.json').write_text(json.dumps(verification, indent=2))
print(json.dumps(dict(summary=summaries, verification=verification), indent=2))
