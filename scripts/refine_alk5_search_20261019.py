"""Local second-round edits; development only, no reserved-set feedback."""
import json, sys
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Crippen, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
root=Path(sys.argv[1]); out=root/'refinement_sar_guided'
out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
protocol=json.loads((root/'protocol.json').read_text())
save('protocol.json',dict(inherited=protocol,edits=['aryl methoxy demethylation','terminal N-demethylation','N-ethyl shortening'],max_rounds=1,selection='Development engineering hits; preserve pyridyl methyl after SAR audit found it invariant across 31 measured compounds. Keep original measured parent comparator. New edits still unmeasured.',stop='All eligible single edits evaluated; no ADMET threshold changes',policy_change='Chemical edit policy updated after local SAR audit; not an independent validation round'))
design=json.loads((root/'input_designs.json').read_text())
lookup={c['candidate_id']:c for c in design['candidates']}
hits=[r for f in sorted(root.glob('development_*.json')) for r in json.loads(f.read_text()) if r['engineering_objectives_met']]
seen={c['normalized_parent'] for c in design['candidates']}
asset=json.loads((root.parent/'backend/app/data/alk5_design_references.json').read_text())
seen.update(asset['excluded_normalized_parents'])
enum=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger()
fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def ident(m):
    m=enum.Canonicalize(un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))))
    Chem.RemoveStereochemistry(m)
    return Chem.MolToSmiles(m)
proposals=[]; rejected=[]
for hit in hits:
    c=lookup[hit['candidate_id']]; m=Chem.MolFromSmiles(c['smiles']); parent=Chem.MolFromSmiles(c['provenance']['parent_smiles'])
    for atom in m.GetAtoms():
        if atom.GetAtomicNum()!=6 or atom.GetDegree()!=1 or atom.GetTotalNumHs()!=3: continue
        neighbor=atom.GetNeighbors()[0]
        rw=Chem.RWMol(m)
        if neighbor.GetAtomicNum()==7 and not neighbor.GetIsAromatic():
            rw.RemoveAtom(atom.GetIdx()); method='N-demethylation'
        elif neighbor.GetAtomicNum()==8 and any(n.GetIsAromatic() for n in neighbor.GetNeighbors()):
            rw.RemoveAtom(atom.GetIdx()); method='aryl methoxy to hydroxy'
        elif neighbor.GetAtomicNum()==6 and neighbor.GetDegree()==2 and any(n.GetAtomicNum()==7 for n in neighbor.GetNeighbors()):
            rw.RemoveAtom(atom.GetIdx()); method='N-ethyl to N-methyl'
        else: continue
        try:
            prod=rw.GetMol(); Chem.SanitizeMol(prod)
            key=ident(prod)
            sim=DataStructs.TanimotoSimilarity(fp.GetFingerprint(parent),fp.GetFingerprint(prod))
            valid=(key not in seen and sim>=.6 and Descriptors.MolWt(prod)<=550 and Crippen.MolLogP(prod)<=6 and Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(prod))==c['provenance']['parent_scaffold'] and not any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(prod)))
            if not valid:
                rejected.append(dict(seed=c['candidate_id'],smiles=Chem.MolToSmiles(prod),reason='duplicate/archive/scaffold/stereo/property/similarity constraint')); continue
            seen.add(key)
            proposals.append(dict(candidate_id=f'ALK5_refine_{len(proposals)+1:03d}',seed=c['candidate_id'],smiles=Chem.MolToSmiles(prod),parent_smiles=c['provenance']['parent_smiles'],parent_id=c['provenance']['parent_id'],parent_similarity=sim,edit=method))
        except Exception as exc: rejected.append(dict(seed=c['candidate_id'],reason=str(exc)))
save('proposals.json',proposals); save('rejected.json',rejected)
if not proposals:
    save('summary.json',dict(evaluated=0,reason='No eligible edits')); sys.exit()
import torch
torch.set_num_threads(4)
from admet_ai import ADMETModel
inputs=sorted({s for c in proposals for s in (c['smiles'],c['parent_smiles'])})
df=ADMETModel(num_workers=0).predict(smiles=inputs)
assert df.index.is_unique and set(inputs)<=set(df.index)
df.to_csv(out/'all_predictions.csv',float_format='%.17g')
keys=['DILI','hERG','AMES','Solubility_AqSolDB','Bioavailability_Ma','Caco2_Wang']
rows=[]
for c in proposals:
    delta={k:float(df.loc[c['smiles'],k]-df.loc[c['parent_smiles'],k]) for k in keys}
    passed=(delta['Solubility_AqSolDB']>=.3 and all(delta[k]<=.02 for k in keys[:3]) and any(delta[k]<=-.05 for k in keys[:3]) and delta['Bioavailability_Ma']>=-.05 and delta['Caco2_Wang']>=-.2)
    rows.append(dict(**c,delta=delta,engineering_objectives_met=passed,activity_validated=False,safety_validated=False))
save('results.json',rows)
summary=dict(evaluated=len(rows),engineering_hits=[r['candidate_id'] for r in rows if r['engineering_objectives_met']],biological_validation=False)
save('summary.json',summary); print(json.dumps(summary),flush=True)
