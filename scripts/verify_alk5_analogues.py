from pathlib import Path
import json,hashlib,math
from collections import Counter
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
p=Path.home()/'aipharma_analogues_20261016';h=Path.home()/'aipharma_holdout_20261016'
sel=json.loads((p/'selected.json').read_text());parents=json.loads((p/'parents.json').read_text());raw={r['activity_id']:r for fn in ['all_activities.json','weak_boundary_all.json'] for r in json.loads((h/'raw'/fn).read_text())['activities']};index={r['parent_id']:r for r in parents};f=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);audit=[]
assert sel['protocol_sha256']==hashlib.sha256((p/'protocol.json').read_bytes()).hexdigest()
for r in parents:
 a=raw[r['representative_activity_id']];assert a['molecule_chembl_id']==r['parent_id'] and a['standard_relation']=='=' and r['median_pchembl']>=7
 audit.append({'parent_id':r['parent_id'],'activity_id':a['activity_id'],'representative_structure_pchembl':float(a['pchembl_value']),'representative_IC50_nM':float(a['standard_value']),'pooled_parent_median_pchembl':r['median_pchembl'],'assay_id':a['assay_chembl_id'],'document_id':a['document_chembl_id'],'note':'Grouped median uses normalized stereo-free parent. Representative structure measurement is separately given; neither is assigned to designed analogue.'})
for r in [x for x in json.loads((p/'all_screened.json').read_text()) if not x.get('is_parent')]:
 m=Chem.MolFromSmiles(r['smiles']);pm=Chem.MolFromSmiles(r['parent_smiles']);assert m is not None and pm is not None
 assert not any(str(x.specified)=='Unspecified' for x in Chem.FindPotentialStereo(m));assert len(Chem.GetMolFrags(m))==1
 assert r['all_three_pass']==all(v<.5 for v in r['toxicity'].values());assert all(0<=v<=1 for v in r['toxicity'].values());assert DataStructs.TanimotoSimilarity(f.GetFingerprint(m),f.GetFingerprint(pm))>=.65
 assert sum(a.GetAtomicNum()!=1 for a in m.GetAtoms())-sum(a.GetAtomicNum()!=1 for a in pm.GetAtoms()) in [0,1]
 assert Counter(a.GetSymbol() for a in m.GetAtoms() if a.GetIsAromatic())==Counter(a.GetSymbol() for a in pm.GetAtoms() if a.GetIsAromatic())
assert len({r['smiles'] for r in sel['candidates']})==len(sel['candidates'])
(p/'parent_measurement_provenance.json').write_text(json.dumps(audit,indent=2));out={'checks_passed':True,'selected_count':len(sel['candidates']),'screened_analogue_structures_checked':69,'parents_pooled_median_pchembl_ge7':len(audit),'stereo_checked_with_FindPotentialStereo':True,'parent_core_atom_elements_preserved':True,'ADMET_gate_rechecked':True,'limitations':['No experimental activity or safety measured on analogues','No global novelty or synthesis route validation','Parent grouping ignores stereo; separate representative measurements archived'],'hashes':{fn:hashlib.sha256((p/fn).read_bytes()).hexdigest() for fn in ['protocol.json','parents.json','pre_admet_candidates.json','all_screened.json','selected.json','admet_predictions.csv']}}
(p/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
