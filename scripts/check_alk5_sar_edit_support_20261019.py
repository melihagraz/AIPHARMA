"""Check whether the observed series actually perturbs the newly edited site."""
import json,sys
from pathlib import Path
from rdkit import Chem
root=Path(sys.argv[1])
rows=json.loads((root/'local_sar.json').read_text())
queries={
 'pyridyl_methyl':'[n;r6]1[c;r6]([CH3])[c;r6][c;r6][c;r6][c;r6]1',
 'pyridyl_hydroxymethyl':'[n;r6]1[c;r6]([CH2][OX2H])[c;r6][c;r6][c;r6][c;r6]1',
 'aryl_ethyl_carbamate':'[c][OX2][CX3](=[OX1])[NX3][CH2][CH3]'}
queries={k:Chem.MolFromSmarts(v) for k,v in queries.items()}
assert all(queries.values())
results=[]
for r in rows:
    pool={x['analogue_smiles'] for x in r['same_assay_parent_pairs']}
    pool.update(x['smiles'] for x in r['neighbors'] if x['is_measured_parent'])
    molecules=[Chem.MolFromSmiles(s) for s in sorted(pool)]
    results.append(dict(candidate_id=r['candidate_id'],candidate_features={k:Chem.MolFromSmiles(r['smiles']).HasSubstructMatch(q) for k,q in queries.items()},
        distinct_measured_parent_series=len(pool),measured_feature_counts={k:sum(m.HasSubstructMatch(q) for m in molecules) for k,q in queries.items()},
        warning='Feature counts are site-support diagnostics, not causal matched-pair analysis; absence of the edited feature provides no direct measured support.'))
(root/'edit_support.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
