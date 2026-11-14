#!/usr/bin/env python3
"""Verify archived ALK5 ADMET numerical provenance without re-running ADMET-AI.

This DOES NOT generate fresh model predictions or validate molecular safety.
"""
from __future__ import annotations
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / 'reports/alk5_design_admet_20261019'
OUT = ROOT / 'reports/phase7/archived_admet_provenance.json'


def main():
    comparison = json.loads((DIR / 'comparison.json').read_text())
    manifest = json.loads((DIR / 'MANIFEST_SHA256.json').read_text())
    inputs = json.loads((DIR / 'input_smiles.json').read_text())
    endpoint_names = comparison['endpoints']
    pairs = comparison['pairs']
    failed = []
    for rel, expected_sha in manifest.items():
        p = DIR / rel
        actual = hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        if actual != expected_sha:
            failed.append(rel)
    counts = Counter()
    numerical_failed = []
    for i, pair in enumerate(pairs):
        for endpoint in endpoint_names:
            row = pair['endpoints'].get(endpoint)
            if not row:
                numerical_failed.append(f'{i}/{endpoint}/missing')
                continue
            diff = float(row['candidate']) - float(row['parent'])
            if abs(diff - float(row['delta'])) > 1e-6:
                numerical_failed.append(f'{i}/{endpoint}/delta')
            counts[endpoint] += 1
    meta = comparison['meta']
    version = meta.get('versions', {}).get('admet-ai')
    accepted = (
        not failed and not numerical_failed and len(pairs) == 30 and
        len(endpoint_names) == 6 and len(inputs) == 46 and
        len(set(inputs)) == 46 and all(counts[x] == 30 for x in endpoint_names) and
        version == '2.0.1'
    )
    result = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'archive': str(DIR.relative_to(ROOT)),
        'archive_validated': accepted,
        'manifest_file_count': len(manifest),
        'manifest_mismatches': failed,
        'design_pairs': len(pairs),
        'distinct_input_molecules': len(set(inputs)),
        'endpoint_count': len(endpoint_names),
        'numerical_pairs_checked': sum(counts.values()),
        'numerical_errors': numerical_failed,
        'archived_admet_ai_version': version,
        'archived_torch_version': meta.get('versions',{}).get('torch'),
        'archived_rdkit_version': meta.get('versions',{}).get('rdkit'),
        'this_run_inference_performed': False,
        'safety_or_alk5_activity_validation': False,
        'limitations': 'Hashes and arithmetic only; no fresh ADMET model inference or laboratory validation',
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return 0 if accepted else 1

if __name__ == '__main__':
    raise SystemExit(main())
