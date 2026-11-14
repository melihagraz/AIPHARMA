"""Run and archive the explicit design-only branch without ADMET or docking."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ap = argparse.ArgumentParser()
ap.add_argument('root', type=Path)
ap.add_argument('output', type=Path)
args = ap.parse_args()
sys.path.insert(0, str(args.root/'backend'))
from app.services.alk5_analogue_design import ASSET, design_analogues

args.output.mkdir(parents=True, exist_ok=True)
protocol = dict(method='reference_anchored_brics_v1', reference_asset_sha256=hashlib.sha256(ASSET.read_bytes()).hexdigest(), count=30,
    design='Single BRICS peripheral substitution using donor fragments from measured catalytic ALK5 actives; exact parent Murcko scaffold retained',
    boundaries='Parent Morgan radius2/2048 Tanimoto 0.5–0.9, MW<=550, logP in [-2,6], heavy atoms<=parent+8; these are design bounds, not validated activity or safety criteria',
    archive_exclusion='2514 normalized archived ChEMBL parents; no global novelty claim',
    selection='Deterministic enumeration, rotated donor ordering across parents; up to8 per parent, round robin output, not activity ranking',
    statistical_model='none', admet='not run', docking='not run', synthesis='not demonstrated',
    success='Valid unique complete structures; all retain parent scaffold and traceable measurement/fragment provenance; deterministic output',
    previous_analogue_test='Earlier F/Cl atom-edit series remains preserved; no claim this design branch passed its failed ADMET thresholds')
(args.output/'protocol.json').write_text(json.dumps(protocol, indent=2))
result = design_analogues(30)
(args.output/'designs.json').write_text(json.dumps(result, indent=2))
print(json.dumps({k: v for k, v in result.items() if k != 'candidates'}, indent=2))
