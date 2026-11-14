from pathlib import Path
import json,hashlib,platform
import numpy as np,sklearn,rdkit
from sklearn.metrics import roc_auc_score
p=Path.home()/'aipharma_holdout_20261016/qsar_exploratory';m=json.loads((p/'metrics.json').read_text());r=json.loads((p/'dataset.json').read_text());pp=json.loads((p/'predictions.json').read_text());old=json.loads((p.parent/'locked_bounded/label_key.json').read_text())
assert len(r)==542 and len({x['smiles'] for x in r})==542
assert not {x['scaffold'] for x in r}&{x['scaffold'] for x in old.values()}
for split,e in m['evaluations'].items():
 assert sum(f['test'] for f in e['folds'])==len(r)
 assert all(f['scaffold_overlap']==0 for f in e['folds'])
 if split=='document_and_scaffold':assert all(f['document_overlap']==0 for f in e['folds'])
 for method,mm in e['methods'].items():
  z=[x for x in pp if x['split']==split and x['method']==method]
  assert len(z)==len(r) and len({x['smiles'] for x in z})==len(r)
  assert abs(roc_auc_score([x['label'] for x in z],[x['score'] for x in z])-mm['roc_auc'])<1e-12
c=json.loads((p/'candidate_applicability.json').read_text())['candidates'];assert len(c)==20
out={'checks_passed':True,'dataset_unique_parents':542,'old_holdout_scaffold_overlap':0,'OOF_predictions_per_method':542,'reported_AUC_independently_recomputed':True,'source_split_document_and_scaffold_overlap':0,'candidate_count':20,'candidates_below_reference_5th_percentile':sum(x['similarity_percentile_vs_source_disjoint_test']<5 for x in c),'python':platform.python_version(),'sklearn':sklearn.__version__,'rdkit':rdkit.__version__,'hashes':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ['dataset.json','protocol.json','metrics.json','predictions.json','candidate_applicability.json']}}
(p/'verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
