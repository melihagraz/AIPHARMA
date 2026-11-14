"""Fixed exploratory ALK5 model comparison; old revealed holdout excluded entirely."""
from pathlib import Path
import json,hashlib,statistics,csv
from collections import defaultdict,Counter
import numpy as np
from rdkit import Chem,DataStructs,RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score
RDLogger.DisableLog('rdApp.*')
p=Path.home()/'aipharma_holdout_20261016';o=p/'qsar_exploratory';o.mkdir(exist_ok=True)
protocol={'purpose':'Exploratory fixed-model comparison, not prospective validation','seed':20261016,'active':'median exact pChEMBL >=7','weak':'median exact pChEMBL <=5 or only IC50 right bounds >=10000 nM','model':'Morgan radius2 2048bits; logistic regression C1 balanced, no tuning','baseline':'max similarity to training-fold actives','splits':['5-fold scaffold grouped','5-fold connected document+scaffold grouped when feasible'],'old_holdout':'all 36 parents and scaffolds excluded','caveat':'Historical public data already inspected; no new blind-trial claim; censoring and publication confounding remain'}
(o/'protocol.json').write_text(json.dumps(protocol,indent=2))
e=rdMolStandardize.TautomerEnumerator();u=rdMolStandardize.Uncharger();cache={}
def norm(s):
 if s not in cache:
  m=Chem.MolFromSmiles(s)
  if m is None:return None
  m=u.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m)));m=e.Canonicalize(m);Chem.RemoveStereochemistry(m)
  cache[s]=(Chem.MolToSmiles(m),Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) or 'ACYCLIC')
 return cache[s]
old=json.loads((p/'locked_bounded/label_key.json').read_text());excluded={r['scaffold'] for r in old.values()}
a={r['assay_chembl_id']:r for r in json.loads((p/'raw/boundary_assays_all.json').read_text())['assays']}
raw={r['activity_id']:r for fn in ['all_activities.json','weak_boundary_all.json'] for r in json.loads((p/'raw'/fn).read_text())['activities']};g=defaultdict(list)
for r in raw.values():
 aa=a[r['assay_chembl_id']]
 if aa.get('assay_type')!='B' or aa.get('confidence_score')!=9 or aa.get('target_chembl_id')!='CHEMBL4439':continue
 if r.get('standard_type')!='IC50' or r.get('standard_units')!='nM' or r.get('data_validity_comment') or r.get('potential_duplicate') in [True,1,'1']:continue
 exact=r.get('standard_relation')=='=' and r.get('pchembl_value') is not None
 bound=r.get('standard_relation') in ['>','>='] and float(r.get('standard_value') or 0)>=10000
 if not (exact or bound) or not r.get('canonical_smiles'):continue
 n=norm(r['canonical_smiles'])
 if n and n[1] not in excluded:g[n[0]].append((r,exact,n[1]))
rows=[]
for smi,rr in sorted(g.items()):
 vals=[float(r['pchembl_value']) for r,x,s in rr if x];bounded=any(not x for r,x,s in rr)
 if vals and ((min(vals)<=5 and max(vals)>=7) or (bounded and max(vals)>=7)):continue
 med=statistics.median(vals) if vals else None
 label=1 if med is not None and med>=7 else 0 if (med is not None and med<=5) or (med is None and bounded) else None
 if label is None:continue
 rows.append({'smiles':smi,'scaffold':rr[0][2],'label':label,'median_pchembl':med,'evidence':'exact' if vals else 'censored','documents':sorted({r['document_chembl_id'] for r,x,s in rr if r.get('document_chembl_id')}),'activity_ids':[r['activity_id'] for r,x,s in rr]})
(o/'dataset.json').write_text(json.dumps(rows,indent=2));fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);fps=[fp.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in rows];X=np.array([np.asarray(x) for x in fps]);y=np.array([r['label'] for r in rows]);N=len(rows)
# Union documents AND scaffolds to prohibit either crossing the source-aware split.
parent=list(range(N))
def root(i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
seen={}
for i,r in enumerate(rows):
 for key in ['S:'+r['scaffold']]+['D:'+d for d in r['documents']]:
  if key in seen:parent[root(i)]=root(seen[key])
  seen[key]=i
reports={};preds=[]
for name,groups in [('scaffold',[r['scaffold'] for r in rows]),('document_and_scaffold',[str(root(i)) for i in range(N)])]:
 counts=[len({groups[i] for i in range(N) if y[i]==c}) for c in [0,1]]
 if min(counts)<5:reports[name]={'status':'infeasible','class_group_counts':counts};continue
 splits=list(StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=20261016).split(X,y,groups));lr=np.full(N,np.nan);sim=np.full(N,np.nan);folds=[]
 for k,(tr,te) in enumerate(splits):
  if len(set(y[tr]))<2:break
  model=LogisticRegression(C=1,class_weight='balanced',max_iter=2000,solver='liblinear',random_state=20261016).fit(X[tr],y[tr]);lr[te]=model.decision_function(X[te]);act=[fps[j] for j in tr if y[j]==1]
  for j in te:sim[j]=max(DataStructs.BulkTanimotoSimilarity(fps[j],act))
  folds.append({'fold':k,'train':len(tr),'test':len(te),'test_active':int(y[te].sum()),'test_weak':int((y[te]==0).sum()),'scaffold_overlap':len({rows[j]['scaffold'] for j in tr}&{rows[j]['scaffold'] for j in te}),'document_overlap':len({d for j in tr for d in rows[j]['documents']}&{d for j in te for d in rows[j]['documents']})})
 if np.isnan(lr).any():reports[name]={'status':'infeasible_training_class','folds':folds};continue
 report={'status':'complete','folds':folds,'class_group_counts':counts,'methods':{}}
 for method,score in [('logistic',lr),('similarity',sim)]:
  order=np.argsort(-score,kind='stable');report['methods'][method]={'roc_auc':float(roc_auc_score(y,score)),'average_precision':float(average_precision_score(y,score)),'top10_active':int(y[order[:10]].sum()),'fold_AUC':[float(roc_auc_score(y[te],score[te])) if len(set(y[te]))==2 else None for tr,te in splits]}
  for i,r in enumerate(rows):preds.append({'split':name,'method':method,'smiles':r['smiles'],'label':r['label'],'score':float(score[i])})
 # Paired group bootstrap for AUC difference; prevents row-wise pseudo-replication.
 rng=np.random.default_rng(20261016);ug=sorted(set(groups));indices={gg:np.array([i for i in range(N) if groups[i]==gg]) for gg in ug};diff=[]
 for _ in range(1000):
  ix=np.concatenate([indices[gg] for gg in rng.choice(ug,len(ug),replace=True)])
  if len(set(y[ix]))==2:diff.append(roc_auc_score(y[ix],lr[ix])-roc_auc_score(y[ix],sim[ix]))
 report['auc_logistic_minus_similarity_group_bootstrap_CI']=np.quantile(diff,[.025,.975]).tolist();reports[name]=report
summary={'protocol':protocol,'n':N,'classes':dict(Counter(map(int,y))),'scaffolds':len(set(r['scaffold'] for r in rows)),'documents':len({d for r in rows for d in r['documents']}),'weak_evidence':dict(Counter(r['evidence'] for r in rows if r['label']==0)),'evaluations':reports,'warning':'OOF pooled scores may vary by fold; consult per-fold AUC. Data-source/censoring confounding remains. Not biological validation.'}
(o/'metrics.json').write_text(json.dumps(summary,indent=2));(o/'predictions.json').write_text(json.dumps(preds));print(json.dumps(summary,indent=2),flush=True)
