"""Locked 30-design/parent paired ADMET screen, preserving all outputs."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time

ap = argparse.ArgumentParser()
ap.add_argument('root', type=Path)
ap.add_argument('designs', type=Path)
ap.add_argument('output', type=Path)
args = ap.parse_args()
sys.path.insert(0, str(args.root/'backend'))
from app.services.alk5_design_admet import ENDPOINTS, compare_profiles, profile_inputs

args.output.mkdir(parents=True, exist_ok=True)
if (args.output/'comparison.json').exists():
    raise RuntimeError('Completed comparison exists; do not overwrite')
data = args.designs.read_bytes()
design = json.loads(data)
inputs = profile_inputs(design)
protocol = dict(designs_sha256=hashlib.sha256(data).hexdigest(), unique_structures=len(inputs), designs=len(design['candidates']), endpoints=list(ENDPOINTS), selection='All locked designs and their parents, no post-prediction selection', comparison='Candidate minus own parent model output', gate='none; no 0.5 cutoff or safety decision', caveat='Predictions for both designs and measured-active parents; measured parent activity is not measured ADMET. No model retraining, docking or measured activity evaluation.')
(args.output/'protocol.json').write_text(json.dumps(protocol, indent=2))
(args.output/'input_designs.json').write_bytes(data)
(args.output/'input_smiles.json').write_text(json.dumps(inputs, indent=2))
import torch
torch.set_num_threads(4)
from admet_ai import ADMETModel
started = time.perf_counter()
model = ADMETModel(num_workers=0)
df = model.predict(smiles=inputs)
df.to_csv(args.output/'all_admet_predictions.csv', float_format='%.17g')
# Match by returned structure, never silently by row position.
assert df.index.is_unique
assert set(inputs) <= set(df.index), 'Prediction index differs from input structures'
predictions = [dict(smiles=s, raw={k: float(df.loc[s, k]) for k in df.columns if not k.endswith('_drugbank_approved_percentile')}) for s in inputs]
comparison = compare_profiles(design, predictions)
comparison['meta']['runtime_seconds'] = time.perf_counter()-started
comparison['meta']['versions'] = {p: importlib.metadata.version(p) for p in ['admet-ai', 'torch', 'rdkit']}
(args.output/'comparison.json').write_text(json.dumps(comparison, indent=2, allow_nan=False))
print(json.dumps(dict(summary=comparison['summary'], meta=comparison['meta']), indent=2), flush=True)
