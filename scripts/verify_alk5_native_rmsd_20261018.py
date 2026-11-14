"""Independent RDKit in-place native RMSD audit; never substitutes docking poses."""
from pathlib import Path
import json
import sys
from rdkit import Chem
from rdkit.Chem import rdMolAlign,rdqueries
from run_alk5_five_block_validation_20261018 import coordinates

out=Path(sys.argv[1]) if len(sys.argv)>1 else Path.home()/'aipharma_five_block_validation_20261018';records=[]
previous=out/'native_rmsd_independent_audit.json'
if previous.exists() and not (out/'native_rmsd_independent_audit_initial.json').exists():
    (out/'native_rmsd_independent_audit_initial.json').write_bytes(previous.read_bytes())
for folder in sorted((out/'jobs').glob('*_native_*')):
    result=json.loads((folder/'result.json').read_text())
    if result.get('native_rmsd_A') is None:continue
    template=(out/'ligands'/result['id']/'ligand.pdbqt').read_text()
    mol,xyz=coordinates(template,(folder/'poses.pdbqt').read_text())
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,point in enumerate(xyz):conf.SetAtomPosition(i,tuple(point))
    mol.AddConformer(conf)
    ref=Chem.SDMolSupplier(str(out/'receptors'/result['pdb']/'native_ref.sdf'),removeHs=True)[0]
    row={k:result[k] for k in ['id','pdb','function','seed','native_rmsd_A','native_pose_pass']}
    try:
        # Match the primary element/connectivity graph convention, independent of
        # OpenBabel versus SMILES aromatic bond orders and protonation flags.
        query=Chem.RWMol(mol);anybond=Chem.MolFromSmarts('[*]~[*]').GetBondWithIdx(0)
        for atom in mol.GetAtoms():query.ReplaceAtom(atom.GetIdx(),rdqueries.AtomNumEqualsQueryAtom(atom.GetAtomicNum()))
        for bond in mol.GetBonds():query.ReplaceBond(bond.GetIdx(),anybond)
        matches=ref.GetSubstructMatches(query.GetMol(),uniquify=False,maxMatches=100000,useChirality=False)
        assert matches and len(matches)<100000
        maps=[list(enumerate(match)) for match in matches]
        value=float(rdMolAlign.CalcRMS(mol,ref,map=maps,symmetrizeConjugatedTerminalGroups=False));delta=abs(value-result['native_rmsd_A'])
        row.update({'rdkit_in_place_RMSD_A':value,'absolute_difference_A':delta,'same_pass_fail_decision':(value<2)==result['native_pose_pass'],
                    'agreement_within_0_01A':delta<.01,'explicit_graph_isomorphisms':len(matches)})
    except Exception as exc:row['audit_error']=str(exc)
    records.append(row)
summary={'native_records':len(records),'graph_convention':'Atomic element and heavy-atom connectivity; explicit RDKit graph isomorphisms, no bond-order/aromaticity requirement, same primary convention',
         'agreements_within_0_01A':sum(r.get('agreement_within_0_01A',False) for r in records),
         'same_pass_fail_decisions':sum(r.get('same_pass_fail_decision',False) for r in records),
         'superposition_used':False,'records':records}
(out/'native_rmsd_independent_audit.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='records'},indent=2),flush=True)
