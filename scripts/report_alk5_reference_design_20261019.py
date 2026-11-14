"""Verify measurement provenance and render parent/design structure examples."""
import argparse
import hashlib
import json
from pathlib import Path
from statistics import mean
from rdkit import Chem
from rdkit.Chem import Draw

ap = argparse.ArgumentParser()
ap.add_argument('asset', type=Path)
ap.add_argument('report', type=Path)
ap.add_argument('raw', type=Path)
args = ap.parse_args()
asset = json.loads(args.asset.read_text())
result = json.loads((args.report/'designs.json').read_text())
raw = {r['activity_id']: r for name in ['all_activities.json', 'weak_boundary_all.json'] for r in json.loads((args.raw/name).read_text())['activities']}
for ref in asset['references']:
    r = raw[ref['representative_activity_id']]
    assert r['standard_type'] == 'IC50' and r['standard_units'] == 'nM' and r['standard_relation'] == '='
    assert float(r['pchembl_value']) >= 7
    assert r['molecule_chembl_id'] == ref['id']
    assert ref['representative_activity_id'] in ref['activity_ids']
assert hashlib.sha256(args.asset.read_bytes()).hexdigest() == result['meta']['reference_asset_sha256']
rows = result['candidates']
measurements = dict(reference_count=len(asset['references']), measured_active_representatives_verified=True, design_parent_count=sum(r['design_parent'] for r in asset['references']), archived_parent_exclusion_count=len(asset['excluded_normalized_parents']))
summary = dict(requested=result['requested'], returned=len(rows), represented_parent_actives=len({r['provenance']['parent_id'] for r in rows}), unique_scaffolds=result['metrics']['unique_scaffolds'], replacement_fragments_used=result['metrics']['replacement_fragments_used'], parent_similarity_min=min(r['parent_similarity'] for r in rows) if rows else None, parent_similarity_max=max(r['parent_similarity'] for r in rows) if rows else None, mean_pairwise_tanimoto=result['metrics']['mean_pairwise_tanimoto'], mean_mw_change=mean(r['property_changes']['mw'] for r in rows) if rows else None, mean_logp_change=mean(r['property_changes']['logp'] for r in rows) if rows else None, **measurements, activity_validated=False, safety_validated=False, synthesis_validated=False)
(args.report/'summary.json').write_text(json.dumps(summary, indent=2))
mols, labels = [], []
for r in rows[:6]:
    mols += [Chem.MolFromSmiles(r['provenance']['parent_smiles']), Chem.MolFromSmiles(r['smiles'])]
    labels += [r['provenance']['parent_id']+' (measured parent)', r['candidate_id']+' (unmeasured design)']
img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(460, 260), legends=labels)
img.save(str(args.report/'parent_design_examples.png'))
print(json.dumps(summary, indent=2))
