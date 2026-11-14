"""Actual structures and model deltas; no implied activity ranking."""
import json,sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
p=Path(sys.argv[1])
base=next(r for r in json.loads((p/'input_designs.json').read_text())['candidates'] if r['candidate_id']=='ALK5_design_020')
old=json.loads((p/'refinement/results.json').read_text())[0]
new=next(r for r in json.loads((p/'refinement_sar_guided/results.json').read_text()) if r['candidate_id']=='ALK5_refine_003')
structures=[base['provenance']['parent_smiles'],base['smiles'],old['smiles'],new['smiles']]
legends=['Measured reference CHEMBL1835522\nIC50 2.26 nM in CHEMBL1837692\nReference activity does not transfer to designs',
    'Initial design 020: UNMEASURED\nEthyl carbamate substitution\nNo direct measured support for new group',
    'Round 2: UNMEASURED\nPyridyl CH3 -> CH2OH\nEdited site invariant in measured series',
    'Round 3: UNMEASURED\nCarbamate ethyl -> methyl; pyridyl CH3 retained\nSame-model property improvement only']
mols=[Chem.MolFromSmiles(s) for s in structures]; assert all(mols)
im=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(650,440),legends=legends)
im.save(p/'design_lineage.png')
(p/'design_lineage.svg').write_text(Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(650,440),legends=legends,useSVG=True))
