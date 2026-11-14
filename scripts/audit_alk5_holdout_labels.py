"""Verify every revealed label against archived raw experimental records."""
from pathlib import Path
import json,statistics
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
p=Path.home()/'aipharma_holdout_20261016'
if not (p/'docking/summary.json').exists():raise RuntimeError('Do not open labels before scoring completes')
labels=json.loads((p/'locked_bounded/label_key.json').read_text());allrows={r['activity_id']:r for fn in ['all_activities.json','weak_boundary_all.json'] for r in json.loads((p/'raw'/fn).read_text())['activities']};assays={a['assay_chembl_id']:a for a in json.loads((p/'raw/boundary_assays_all.json').read_text())['assays']};te=rdMolStandardize.TautomerEnumerator();uc=rdMolStandardize.Uncharger()
def norm(s):
 m=te.Canonicalize(uc.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(Chem.MolFromSmiles(s)))));Chem.RemoveStereochemistry(m);return Chem.MolToSmiles(m)
verified=[]
for case,row in labels.items():
 records=[allrows[a] for a in row['activity_ids']]
 for r in records:
  a=assays[r['assay_chembl_id']];assert a['target_chembl_id']=='CHEMBL4439' and a['assay_type']=='B' and a['confidence_score']==9
  assert r['standard_type']=='IC50' and r['standard_units']=='nM';assert not r.get('data_validity_comment') and r.get('potential_duplicate') not in [True,1,'1'];assert norm(r['canonical_smiles'])==row['parent_identity']
 if row['label']==1:
  exact=[float(r['pchembl_value']) for r in records if r['standard_relation']=='=' and r.get('pchembl_value') is not None];assert statistics.median(exact)>=7;assert row['median_pchembl']==statistics.median(exact)
 else:
  assert row['median_pchembl'] is None and all(r['standard_relation'] in ['>','>='] and float(r['standard_value'])>=10000 for r in records)
 verified.append(case)
report={'verified_case_count':len(verified),'all_labels_supported_by_archived_measurements':True,'human_target':'CHEMBL4439','high_confidence_biochemical_assay_only':True,'censored_weak_values_not_assigned_exact_potency':True,'verified_case_ids':sorted(verified)}
(p/'evaluation/ground_truth_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
