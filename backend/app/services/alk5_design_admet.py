"""Pair model outputs with measured-parent structures; never infer safety."""
from collections import Counter
import hashlib
import json
import math

ENDPOINTS = ('DILI', 'hERG', 'AMES', 'Bioavailability_Ma', 'Caco2_Wang', 'Solubility_AqSolDB')
CLASSIFICATION = {'DILI', 'hERG', 'AMES', 'Bioavailability_Ma'}


def profile_inputs(design):
    return list(dict.fromkeys(s for c in design['candidates'] for s in (c['smiles'], c['provenance']['parent_smiles'])))


def compare_profiles(design, predictions):
    lookup = {}
    for row in predictions:
        if not isinstance(row, dict) or not isinstance(row.get('smiles'), str) or not isinstance(row.get('raw') or {}, dict):
            raise ValueError('Malformed prediction row')
        if row['smiles'] in lookup:
            raise ValueError('Duplicate prediction SMILES; ambiguous alignment')
        lookup[row['smiles']] = row.get('raw') or {}
    def value(smiles, endpoint):
        v = lookup.get(smiles, {}).get(endpoint)
        if v is None or isinstance(v, bool):
            return None
        try:
            v = float(v)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(v) or (endpoint in CLASSIFICATION and not 0 <= v <= 1):
            return None
        return v
    pairs = []
    counts = {e: Counter(lower=0, higher=0, unchanged=0, missing=0) for e in ENDPOINTS}
    for c in design['candidates']:
        endpoints = {}
        for endpoint in ENDPOINTS:
            parent = value(c['provenance']['parent_smiles'], endpoint)
            candidate = value(c['smiles'], endpoint)
            delta = candidate-parent if candidate is not None and parent is not None else None
            direction = 'missing' if delta is None else 'unchanged' if abs(delta) <= 1e-6 else 'higher' if delta > 0 else 'lower'
            counts[endpoint][direction] += 1
            endpoints[endpoint] = dict(parent=parent, candidate=candidate, delta=delta, direction=direction)
        pairs.append(dict(candidate_id=c['candidate_id'], smiles=c['smiles'], parent_id=c['provenance']['parent_id'], parent_smiles=c['provenance']['parent_smiles'], endpoints=endpoints))
    return dict(pairs=pairs, endpoints=list(ENDPOINTS), summary={k: dict(v) for k, v in counts.items()},
        meta=dict(model='ADMET-AI', design_set_sha256=hashlib.sha256(json.dumps([(c['candidate_id'], c['smiles']) for c in design['candidates']], separators=(',', ':')).encode()).hexdigest(), reference_asset_sha256=design['meta']['reference_asset_sha256'], requested=design['requested'], unique_input_structures=len(profile_inputs(design)), missing_prediction_structures=sum(s not in lookup for s in profile_inputs(design)), interpretation='Differences in uncalibrated model outputs; not measured improvement, clinical risk or safety validation', decision='no_pass_fail_gate_no_activity_ranking', numeric_equality_tolerance=1e-6, activity_validated=False, safety_validated=False, endpoints_not_statistically_independent=True, parent_prediction_is_not_experimental_admet=True))
