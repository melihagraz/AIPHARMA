from pathlib import Path
import json,csv,hashlib
import numpy as np
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
RDLogger.DisableLog('rdApp.*')
p=Path.home()/'aipharma_holdout_20261016/qsar_exploratory';rows=json.loads((p/'dataset.json').read_text());source=Path.home()/'aipharma_generation_20261016/results/locked_top20.json';cc=json.loads(source.read_text());print('candidate container',type(cc).__name__,flush=True)
if isinstance(cc,dict):
 print('keys',list(cc),flush=True)
 cc=cc.get('selected',cc.get('candidates',cc.get('structures',cc.get('ranked',[]))))
f=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);fps=[f.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in rows];X=np.array([np.asarray(x) for x in fps]);y=np.array([r['label'] for r in rows]);model=LogisticRegression(C=1,class_weight='balanced',max_iter=2000,solver='liblinear',random_state=20261016).fit(X,y)
parent=list(range(len(rows)))
def root(i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
seen={}
for i,r in enumerate(rows):
 for key in ['S:'+r['scaffold']]+['D:'+d for d in r['documents']]:
  if key in seen:parent[root(i)]=root(seen[key])
  seen[key]=i
groups=[root(i) for i in range(len(rows))];ref=[]
for tr,te in StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=20261016).split(X,y,groups):
 train=[fps[j] for j in tr]
 for j in te:ref.append(max(DataStructs.BulkTanimotoSimilarity(fps[j],train)))
e=rdMolStandardize.TautomerEnumerator();u=rdMolStandardize.Uncharger();out=[]
for j,r in enumerate(cc):
 m=u.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(Chem.MolFromSmiles(r['smiles']))));m=e.Canonicalize(m);Chem.RemoveStereochemistry(m);ff=f.GetFingerprint(m);sim=DataStructs.BulkTanimotoSimilarity(ff,fps);ix=int(np.argmax(sim));xx=np.asarray(ff).reshape(1,-1);stereo=Chem.FindPotentialStereo(Chem.MolFromSmiles(r['smiles']));out.append({'candidate':r.get('candidate',r.get('id',r.get('candidate_id',f'candidate_{j+1:02d}'))),'smiles':r['smiles'],'logistic_decision_score_not_probability':float(model.decision_function(xx)[0]),'nearest_training_similarity':float(sim[ix]),'max_active_similarity':max(v for v,yy in zip(sim,y) if yy==1),'max_weak_similarity':max(v for v,yy in zip(sim,y) if yy==0),'similarity_percentile_vs_source_disjoint_test':float(np.mean(np.asarray(ref)<=sim[ix])*100),'nearest_training_label':int(y[ix]),'nearest_training_smiles':rows[ix]['smiles'],'unresolved_stereo':sum(str(s.specified)=='Unspecified' for s in stereo)})
assert len(out)==20
(p/'candidate_applicability.json').write_text(json.dumps({'interpretation':'Exploratory scores only; not calibrated probabilities, IC50 predictions or validated ranking. Similarity percentile describes relative coverage and is not an error guarantee. No threshold selected post hoc.','candidate_file_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'reference_similarity_quantiles':dict(zip(['min','q05','q25','median','q75','max'],np.quantile(ref,[0,.05,.25,.5,.75,1]).tolist())),'candidates':out},indent=2))
with (p/'candidate_applicability.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print(json.dumps({'n':len(out),'candidate_similarity_range':[min(r['nearest_training_similarity'] for r in out),max(r['nearest_training_similarity'] for r in out)],'percentile_range':[min(r['similarity_percentile_vs_source_disjoint_test'] for r in out),max(r['similarity_percentile_vs_source_disjoint_test'] for r in out)],'reference_median':float(np.median(ref))},indent=2))
