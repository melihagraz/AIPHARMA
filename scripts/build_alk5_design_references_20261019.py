"""Create traceable, stereo-resolved catalytic-active design inputs before generation."""
import argparse
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import BRICS, Crippen, Descriptors, QED
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

ap = argparse.ArgumentParser()
ap.add_argument('dataset', type=Path)
ap.add_argument('raw_directory', type=Path)
ap.add_argument('output', type=Path)
args = ap.parse_args()
enum = rdMolStandardize.TautomerEnumerator()
uncharger = rdMolStandardize.Uncharger()
def normalized(m):
    m = enum.Canonicalize(uncharger.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))))
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)
raw = {r['activity_id']: r for name in ['all_activities.json', 'weak_boundary_all.json'] for r in json.loads((args.raw_directory/name).read_text())['activities']}
known = set()
for r in raw.values():
    if r.get('canonical_smiles'):
        m = Chem.MolFromSmiles(r['canonical_smiles'])
        if m is not None:
            known.add(normalized(m))
eligible = []
for r in json.loads(args.dataset.read_text()):
    if r['label'] != 1 or r.get('origin') != 'prior_exact_eligible':
        continue
    for activity_id in sorted(r['activity_ids']):
        activity = raw.get(activity_id)
        if not activity or activity.get('standard_relation') != '=' or not activity.get('canonical_smiles'):
            continue
        if activity.get('standard_type') != 'IC50' or activity.get('standard_units') != 'nM' or float(activity.get('pchembl_value') or 0) < 7:
            continue
        m = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(Chem.MolFromSmiles(activity['canonical_smiles'])))
        if normalized(m) != r['smiles']:
            continue
        if any(str(s.specified) == 'Unspecified' for s in Chem.FindPotentialStereo(m)):
            continue
        core = MurckoScaffold.GetScaffoldForMol(m)
        if core.GetNumHeavyAtoms() < 10 or Descriptors.MolWt(m) > 500 or Crippen.MolLogP(m) > 5 or QED.qed(m) < .3:
            continue
        eligible.append(dict(id=activity['molecule_chembl_id'], smiles=Chem.MolToSmiles(m), normalized_parent=r['smiles'], scaffold=Chem.MolToSmiles(core), documents=r['documents'], activity_ids=r['activity_ids'], representative_activity_id=activity_id, parent_median_pchembl=r['median_pchembl'], endpoint=r['stratum'], support=len(r['activity_ids'])))
        break
eligible.sort(key=lambda r: (-r['support'], r['id'], r['smiles']))
refs, seen = [], set()
def has_peripheral_cut(r):
    m = Chem.MolFromSmiles(r['smiles'])
    core = MurckoScaffold.GetScaffoldForMol(m)
    for bond in BRICS.FindBRICSBonds(m):
        parts = Chem.GetMolFrags(BRICS.BreakBRICSBonds(m, bonds=[bond]), asMols=True)
        if len(parts) != 2:
            continue
        for part, side in [parts, parts[::-1]]:
            if part.HasSubstructMatch(core) and not any(a.IsInRing() for a in side.GetAtoms()):
                return True
    return False
for r in eligible:
    if r['scaffold'] in seen or not has_peripheral_cut(r):
        continue
    refs.append(r)
    seen.add(r['scaffold'])
    if len(refs) == 16:
        break
assert refs
parent_ids = {r['id'] for r in refs}
# Donors retain their own measured provenance; do not limit the fragment
# vocabulary to the small scaffold-diverse set of editable parent seeds.
donors = [r for r in eligible if r['id'] not in parent_ids][:128]
refs = [dict(r, design_parent=True) for r in refs] + [dict(r, design_parent=False) for r in donors]
asset = dict(provenance=dict(dataset_sha256=hashlib.sha256(args.dataset.read_bytes()).hexdigest(), source='alk5_structure_mapping_20261018/completed_panel_nodes.json and archived ChEMBL raw activity rows', endpoint='catalytic enzyme IC50, no explicitly annotated mutant; not confirmed wild type', selection='Up to16 distinct editable parent scaffolds plus up to128 measured donors; exact active, raw stereo resolved, MW<=500, logP<=5, QED>=0.3; descending measurement support then ChEMBL ID; parent requires a peripheral BRICS cut preserving scaffold; no candidate activity or ADMET outcome used', parent_activity_not_transferred_to_design=True), references=refs, excluded_normalized_parents=sorted(known))
args.output.write_text(json.dumps(asset, indent=2))
print(json.dumps(dict(references=len(refs), design_parents=sum(r['design_parent'] for r in refs), excluded_archived_parents=len(known))))
