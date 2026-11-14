from pathlib import Path
import json,csv,hashlib
import numpy as np
from rdkit import Chem
from spyrmsd import rmsd
from sklearn.metrics import roc_auc_score
p=Path.home()/'aipharma_ensemble_20261018';o=p/'evaluation';o.mkdir(exist_ok=True)
if not (p/'summary.json').exists():raise RuntimeError('Wait until all scheduled jobs finish')
items=json.loads((p/'locked_items.json').read_text());protocol=json.loads((p/'protocol.json').read_text());native=json.loads((p/'native_pose_audit.json').read_text());rows=[json.loads(f.read_text()) for f in sorted((p/'jobs').glob('*/result.json'))];assert len(rows)==57

def pose_coords(itemid,pdb,seed):
 text=(p/'jobs'/f'{pdb}_{itemid}_{seed}'/'pose.pdbqt').read_text().split('ENDMDL',1)[0];template=(p/'ligands'/itemid/'ligand.pdbqt').read_text();smi=next(l[len('REMARK SMILES '):] for l in template.splitlines() if l.startswith('REMARK SMILES ') and not l.startswith('REMARK SMILES IDX'));m=Chem.MolFromSmiles(smi);mapping=[]
 for l in template.splitlines():
  if l.startswith('REMARK SMILES IDX '):mapping.extend(map(int,l[len('REMARK SMILES IDX '):].split()))
 pos={int(l[6:11]):[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in text.splitlines() if l.startswith(('ATOM','HETATM'))};xyz=np.full((m.GetNumAtoms(),3),np.nan)
 for i,j in zip(mapping[::2],mapping[1::2]):xyz[i-1]=pos[j]
 assert np.isfinite(xyz).all();return m,xyz
aggregate=[];receptors={}
for pdb in protocol['targets']:
 prof=[]
 for item in items:
  rr=[r for r in rows if r['pdb']==pdb and r['id']==item['id']];assert len(rr)==2;good=[r for r in rr if r.get('score') is not None and r['status']=='complete'];row={'pdb':pdb,'id':item['id'],'kind':item['kind'],'scored_seeds':len(good),'scores':[r['score'] for r in good]}
  if good:row.update(mean_score=float(np.mean(row['scores'])),seed_score_range=float(np.ptp(row['scores'])))
  if len(good)==2:
   m,x=pose_coords(item['id'],pdb,20261018);_,y=pose_coords(item['id'],pdb,20261019);nums=np.array([a.GetAtomicNum() for a in m.GetAtoms()]);adj=Chem.GetAdjacencyMatrix(m);row['seed_pose_symmetry_RMSD_A']=float(rmsd.symmrmsd(x,y,nums,nums,adj,adj,minimize=False))
  prof.append(row);aggregate.append(row)
 controls=[r for r in prof if r['kind'] in ['active_control','weak_control'] and r.get('mean_score') is not None];labels=[int(r['kind']=='active_control') for r in controls];auc=float(roc_auc_score(labels,[-r['mean_score'] for r in controls])) if len(set(labels))==2 else None
 receptors[pdb]={'native_control':next(r for r in native if r['pdb']==pdb),'control_count':len(controls),'control_ROC_AUC_descriptive_only':auc,'active_scores':[r['mean_score'] for r in controls if r['kind']=='active_control'],'weak_scores':[r['mean_score'] for r in controls if r['kind']=='weak_control']}
for pdb in protocol['targets']:
 part=sorted([r for r in aggregate if r['pdb']==pdb and 'mean_score' in r],key=lambda r:(r['mean_score'],r['id']))
 for rank,r in enumerate(part,1):r['rank_among_nine']=rank
out={'protocol':protocol,'receptors':receptors,'per_structure':aggregate,'completed_scoring_jobs':sum(r.get('score') is not None for r in rows),'scheduled_jobs':57,'MD':{'status':'not_run','reason':'OpenMM available but ligand parameterization dependencies unavailable; two conda-forge requests reset by peer. No stand-in dynamics claimed.'},'audit':'Native pose report corrected atom mapping using SMILES IDX and MODEL1; no redocking or score retuning. Initial audit exceptions remain in original logs/summary. Per-job corrected records are authoritative.','limits':['Only3 active and3 weak controls; no inferential claims from descriptive AUC','Receptors include previously studied controls, not external blind validation','Score agreement and repeated poses do not establish affinity, biological efficacy or safety','No best candidate chosen for MD from new scores; planned MD trio remained fixed']}
(o/'metrics.json').write_text(json.dumps(out,indent=2))
fields=['pdb','id','kind','scored_seeds','mean_score','seed_score_range','seed_pose_symmetry_RMSD_A','rank_among_nine']
with (o/'scores.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(aggregate)
checks={'scheduled_jobs':57,'result_files':len(rows),'all_three_native_controls_pass':all(r['pose_control_pass'] for r in native),'selection_locked_sha256':hashlib.sha256((p/'locked_items.json').read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((p/'protocol.json').read_bytes()).hexdigest(),'no_MD_claim':True};(o/'verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps({'receptors':receptors,'candidates':[r for r in aggregate if r['kind']=='candidate'],'scoring_completed':out['completed_scoring_jobs']},indent=2),flush=True)
