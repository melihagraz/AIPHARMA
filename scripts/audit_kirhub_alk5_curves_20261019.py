"""Describe reserved curves without fitting or using them for model selection."""
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

p = Path(sys.argv[1])
rows = json.loads((p / 'reserved_alk5_dose_series.json').read_text())
groups = defaultdict(list)
for r in rows:
    groups[r['compound_id']].append(r)
out = []
for identity, points in sorted(groups.items()):
    points.sort(key=lambda r: r['dose_uM'])
    assert len(points) == 10
    assert len({r['dose_uM'] for r in points}) == 10
    crossings = [dict(lower_dose_uM=a['dose_uM'], upper_dose_uM=b['dose_uM'])
                 for a, b in zip(points, points[1:])
                 if a['response_raw'] < 50 <= b['response_raw']]
    out.append(dict(compound_id=identity, compound=points[0]['compound'],
        points=len(points), minimum_dose_uM=points[0]['dose_uM'],
        maximum_dose_uM=points[-1]['dose_uM'],
        maximum_observed_inhibition_percent=max(r['response_raw'] for r in points),
        inhibition_at_maximum_dose_percent=points[-1]['response_raw'],
        upward_50_percent_crossing_brackets=crossings,
        decreases_exceeding_10_percentage_points=sum(a['response_raw']-b['response_raw'] > 10
            for a,b in zip(points, points[1:])),
        interpretation='Observed dose-response values only; no fitted IC50 or potency bound asserted',
        source_cells=[r['source_cell'] for r in points]))
single = json.loads((p/'reserved_alk5_single_dose.json').read_text())
assert len(rows) == 140 and len(groups) == 14 and len(single) == 92
assert len({r['source_cell'] for r in rows}) == 140
assert not json.loads((p/'reserved_missing_values.json').read_text())['series']
source_hash = hashlib.sha256((p/'41587_2026_3090_MOESM4_ESM.xlsx').read_bytes()).hexdigest()
assert source_hash == json.loads((p/'extraction_summary.json').read_text())['source_sha256']
result=dict(verification_passed=True, source_sha256=source_hash,
    compound_curves=out, curves_with_observed_upward_50_percent_crossing=sum(bool(r['upward_50_percent_crossing_brackets']) for r in out),
    curves_never_reaching_50_percent=sum(r['maximum_observed_inhibition_percent'] < 50 for r in out),
    source_reserved=True, models_fitted=False,
    single_dose_source_range=f"A9:L{max(int(r['source_cell'][1:]) for r in single)}")
(p/'curve_audit.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
