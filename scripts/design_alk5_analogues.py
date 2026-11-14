"""Local analogue proposals from experimentally measured ALK5 actives.
No claim of novelty, synthesis feasibility, activity or safety.
"""
from pathlib import Path
import json,hashlib,sys,math
from collections import Counter
import numpy as np
from rdkit import Chem,DataStructs,RDLogger,RDConfig
from rdkit.Chem import Descriptors,Crippen,QED,Lipinski,rdMolDescriptors,rdFingerprintGenerator,FilterCatalog,Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
sys.path.insert(0,str(Path(RDConfig.RDContribDir)/'SA_Score'));import sascorer
RDLogger.DisableLog('rdApp.*')
h=Path.home()/'aipharma_holdout_20261016';o=Path.home()/'aipharma_analogues_20261016';o.mkdir(exist_ok=True)
if (o/'selected.json').exists():raise RuntimeError('Do not overwrite completed selection')
protocol={'design':'single aromatic C-H to F OR single terminal aromatic F/Cl swap; core retained','parent_selection':'up to12 distinct scaffolds, >=7 median pChEMBL, no unresolved stereo, MW<=450 logP<=4, PAINS-free; descending measurement support then median potency then QED','candidate_selection':'up to6 variants per parent before ADMET; predicted AMES/hERG/DILI each<0.5, then up to3 distinct parent scaffolds x3 candidates ordered by parent support and minimum max predicted toxicity','candidate_filters':'MW<=500 logP<=5 HBD<=5 HBA<=10 rot<=10 TPSA<=140 SA<=6 QED>=0.3; no new unresolved stereo; same normalized Murcko scaffold; parent similarity>=0.65','exclude':'all structures in archived ChEMBL rows and both prior generated pools after normalization; no claim of patent novelty','seed':20261016,'warning':'Local atom edits are design hypotheses, NOT validated synthetic routes or activity predictions. Models not experimentally calibrated for these candidates.'}
(o/'protocol.json').write_text(json.dumps(protocol,indent=2))
e=rdMolStandardize.TautomerEnumerator();u=rdMolStandardize.Uncharger();cache={}
def norm(s):
 if s not in cache:
  m=Chem.MolFromSmiles(s)
  if m is None:return None
  m=e.Canonicalize(u.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))));Chem.RemoveStereochemistry(m);cache[s]=Chem.MolToSmiles(m)
 return cache[s]
f=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
params=FilterCatalog.FilterCatalogParams();params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS);pains=FilterCatalog.FilterCatalog(params)
def props(m):return {'mw':Descriptors.MolWt(m),'logp':Crippen.MolLogP(m),'hbd':Lipinski.NumHDonors(m),'hba':Lipinski.NumHAcceptors(m),'rot':rdMolDescriptors.CalcNumRotatableBonds(m),'tpsa':rdMolDescriptors.CalcTPSA(m),'qed':QED.qed(m),'sa':sascorer.calculateScore(m)}
def unspecified(m):return any(str(x.specified)=='Unspecified' for x in Chem.FindPotentialStereo(m))
def scaffold(m):return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
raw=[r for name in ['all_activities.json','weak_boundary_all.json'] for r in json.loads((h/'raw'/name).read_text())['activities']];lookup={r['activity_id']:r for r in raw};known={norm(r['canonical_smiles']) for r in raw if r.get('canonical_smiles')}
known.update(norm(r['smiles']) for r in json.loads((h/'generation_records.json').read_text()) if r.get('smiles'));known.update(norm(r['smiles']) for r in json.loads((h/'pah_run_1PY5.json').read_text())['ranked'])
rows=json.loads((h/'qsar_exploratory/dataset.json').read_text());eligible=[]
for r in rows:
 if r['label']!=1:continue
 # Dataset parents removed stereochemistry: use representative raw structure and check its stereo.
 rr=[lookup[i] for i in r['activity_ids'] if i in lookup and lookup[i].get('standard_relation')=='='];rep=min(rr,key=lambda z:(abs(float(z['pchembl_value'])-r['median_pchembl']),int(z['activity_id'])))
 m=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(Chem.MolFromSmiles(rep['canonical_smiles'])))
 Chem.SanitizeMol(m);pp=props(m)
 if unspecified(m) or pp['mw']>450 or pp['logp']>4 or pains.HasMatch(m):continue
 eligible.append({**r,'smiles':Chem.MolToSmiles(m),'parent_id':rep['molecule_chembl_id'],'representative_activity_id':rep['activity_id'],'properties':pp,'support':len(r['activity_ids'])})
eligible.sort(key=lambda r:(-r['support'],-r['median_pchembl'],-r['properties']['qed'],r['smiles']));parents=[];scafs=set()
for r in eligible:
 if r['scaffold'] in scafs:continue
 parents.append(r);scafs.add(r['scaffold'])
 if len(parents)==12:break
(o/'parents.json').write_text(json.dumps(parents,indent=2));seen=set();candidates=[];counts=Counter()
for parent in parents:
 m=Chem.MolFromSmiles(parent['smiles']);pf=f.GetFingerprint(m);local=[]
 for atom in m.GetAtoms():
  edits=[]
  if atom.GetAtomicNum()==6 and atom.GetIsAromatic() and atom.GetTotalNumHs()==1:edits.append('aromatic_H_to_F')
  if atom.GetAtomicNum() in [9,17] and atom.GetDegree()==1 and atom.GetNeighbors()[0].GetIsAromatic():edits.append('aromatic_F_to_Cl' if atom.GetAtomicNum()==9 else 'aromatic_Cl_to_F')
  for edit in edits:
   counts['attempted']+=1;rw=Chem.RWMol(m)
   if edit=='aromatic_H_to_F':
    idx=rw.AddAtom(Chem.Atom(9));rw.AddBond(atom.GetIdx(),idx,Chem.BondType.SINGLE);rw.GetAtomWithIdx(atom.GetIdx()).SetNumExplicitHs(0)
   else:rw.GetAtomWithIdx(atom.GetIdx()).SetAtomicNum(17 if atom.GetAtomicNum()==9 else 9)
   mm=rw.GetMol()
   try:Chem.SanitizeMol(mm)
   except Exception:counts['invalid']+=1;continue
   smi=Chem.MolToSmiles(mm);identity=norm(smi)
   if identity in known or identity in seen:counts['known_or_duplicate']+=1;continue
   pp=props(mm);sim=DataStructs.TanimotoSimilarity(pf,f.GetFingerprint(mm))
   ok=pp['mw']<=500 and pp['logp']<=5 and pp['hbd']<=5 and pp['hba']<=10 and pp['rot']<=10 and pp['tpsa']<=140 and pp['sa']<=6 and pp['qed']>=.3 and sim>=.65 and not unspecified(mm) and not pains.HasMatch(mm) and norm(scaffold(mm))==norm(scaffold(m))
   if not ok:counts['quality_reject']+=1;continue
   seen.add(identity);local.append({'smiles':smi,'parent_id':parent['parent_id'],'parent_smiles':parent['smiles'],'parent_pchembl':parent['median_pchembl'],'parent_support':parent['support'],'parent_activity_ids':parent['activity_ids'],'parent_documents':parent['documents'],'edit':edit,'edited_parent_atom_index':atom.GetIdx(),'parent_similarity':sim,'properties':pp,'scaffold':parent['scaffold']})
 local.sort(key=lambda r:(-r['parent_similarity'],-r['properties']['qed'],r['smiles']));candidates.extend(local[:6])
(o/'pre_admet_candidates.json').write_text(json.dumps(candidates,indent=2));print(json.dumps({'parents':len(parents),'admet_inputs':len(candidates),'counts':dict(counts)}),flush=True)
import torch;torch.set_num_threads(4)
from admet_ai import ADMETModel
model=ADMETModel(num_workers=0);inputs=candidates+[{'smiles':r['smiles'],'is_parent':True,'parent_id':r['parent_id']} for r in parents];unique=list(dict.fromkeys(r['smiles'] for r in inputs));df=model.predict(smiles=unique);df.to_csv(o/'admet_predictions.csv')
for r in inputs:
 r['toxicity']={k:float(df.loc[r['smiles'],k]) for k in ['AMES','hERG','DILI']};assert all(math.isfinite(v) and 0<=v<=1 for v in r['toxicity'].values());r['all_three_pass']=all(v<.5 for v in r['toxicity'].values())
(o/'all_screened.json').write_text(json.dumps(inputs,indent=2));passing=[r for r in candidates if r['all_three_pass']];families=[]
for par in parents:
 local=[r for r in passing if r['parent_id']==par['parent_id']];local.sort(key=lambda r:(max(r['toxicity'].values()),-r['parent_similarity'],-r['properties']['qed'],r['smiles']))
 if local:families.append((par,local))
selected=[]
for par,local in families[:3]:
 for r in local[:3]:selected.append({**r,'candidate_id':f'analogue_{len(selected)+1:02d}'})
(o/'selected.json').write_text(json.dumps({'protocol_sha256':hashlib.sha256((o/'protocol.json').read_bytes()).hexdigest(),'candidates':selected},indent=2));summary={'parent_count':len(parents),'screened_analogues':len(candidates),'passing_analogues':len(passing),'selected':len(selected),'selected_parent_families':len({r['parent_id'] for r in selected}),'counts':dict(counts),'novelty':'Excluded from archived structures only; no global literature/patent novelty search','synthesis':'No route or reaction yield established','activity':'Unmeasured for proposed analogues; parent potency not assigned to analogue','safety':'Predicted gates only'};(o/'summary.json').write_text(json.dumps(summary,indent=2))
if selected:
 img=Draw.MolsToGridImage([Chem.MolFromSmiles(r['smiles']) for r in selected],molsPerRow=3,subImgSize=(420,300),legends=[r['candidate_id']+' / '+r['parent_id'] for r in selected]);img.save(str(o/'selected_structures.png'))
print(json.dumps(summary,indent=2),flush=True)
