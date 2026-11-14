from pathlib import Path
import json,hashlib,pandas as pd
from rdkit import Chem
p=Path.home()/'aipharma_generation_20261016';r=p/'results';data=json.loads((r/'locked_top20.json').read_text());sel=data['candidates'];assert len(sel)==20
assert len({x['smiles'] for x in sel})==20
assert all(all(v<0.5 for v in x['toxicity'].values()) for x in sel)
assert not any(label=='?' for x in sel for _,label in Chem.FindMolChiralCenters(Chem.MolFromSmiles(x['smiles']),includeUnassigned=True))
assert hashlib.sha256((r/'locked_top20.json').read_bytes()).hexdigest()==(r/'locked_top20.sha256').read_text().strip()
old=json.loads((p/'pah_admet_1PY5.json').read_text())['rows'];df=pd.read_csv(r/'baseline12_admet_recomputed.csv',index_col=0);diff={}
for oldkey,col in [('ames','AMES'),('herg','hERG'),('dili','DILI')]:diff[oldkey]=max(abs(float(df.loc[x['smiles'],col])-x[oldkey]) for x in old)
report={'unique_selected':20,'all_three_thresholds_verified':True,'legacy_chiral_center_check_passed':True,'all_potential_stereogenic_units_resolved':False,'audit_note':'Subsequent FindPotentialStereo audit found unresolved units in four locked candidates; see stereochemistry_and_conformer_audit.json','lock_hash_verified':True,'previous12_prediction_max_absolute_change':diff,'interpretation':'New set is selected on model predictions; no claim of experimental safety or retrained generator improvement'}
(r/'selection_quality_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
