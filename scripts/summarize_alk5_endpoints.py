from pathlib import Path
import json
from collections import Counter
import numpy as np
from sklearn.metrics import roc_auc_score
p=Path.home()/'aipharma_holdout_20261016/endpoint_audit';summ=json.loads((p/'comparison_summary.json').read_text());out={}
for name,data in summ.items():
 rows=json.loads((p/name/'dataset.json').read_text());pred=json.loads((p/name/'predictions.json').read_text());N=len(rows);par=list(range(N))
 def root(i):
  while par[i]!=i:par[i]=par[par[i]];i=par[i]
  return i
 seen={}
 for i,r in enumerate(rows):
  for key in ['S:'+r['scaffold']]+['D:'+d for d in r['documents']]:
   if key in seen:par[root(i)]=root(seen[key])
   seen[key]=i
 groups=[root(i) for i in range(N)];y=np.array([r['label'] for r in rows]);smi=[r['smiles'] for r in rows];method={}
 for split,res in data['evaluations'].items():
  if res['status']!='complete':continue
  gg=groups if split=='document_and_scaffold' else [r['scaffold'] for r in rows];ug=sorted(set(gg),key=str);ixs={g:np.array([i for i in range(N) if gg[i]==g]) for g in ug};score={}
  for met in res['methods']:
   lookup={r['smiles']:r['score'] for r in pred if r['split']==split and r['method']==met};score[met]=np.array([lookup[s] for s in smi])
  rng=np.random.default_rng(20261016);boots={met:[] for met in score}
  for _ in range(2000):
   ix=np.concatenate([ixs[g] for g in rng.choice(ug,len(ug),replace=True)])
   if len(set(y[ix]))<2:continue
   for met,ss in score.items():boots[met].append(float(roc_auc_score(y[ix],ss[ix])))
  method[split]={'fixed_OOF_group_bootstrap_AUC_95CI':{met:np.quantile(b,[.025,.975]).tolist() for met,b in boots.items()},'valid_bootstrap_replicates':len(next(iter(boots.values()))),'macro_fold_AUC':{met:float(np.mean([x for x in mm['fold_AUC'] if x is not None])) for met,mm in res['methods'].items()},'random_top10_expected_active':10*float(y.mean()),'censored_weak_count':sum(r['label']==0 and r['evidence']=='censored' for r in rows),'weak_source_scaffold_components':len({groups[i] for i in range(N) if y[i]==0})}
 out[name]=method
(p/'uncertainty_supplement.json').write_text(json.dumps({'caveat':'Fixed OOF prediction bootstrap; does not retrain models or account for all selection/classification uncertainty. Small negative groups; exploratory CIs only.','strata':out},indent=2));print(json.dumps(out,indent=2))
