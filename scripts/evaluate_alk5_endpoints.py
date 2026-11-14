from pathlib import Path
import json,hashlib,math
from collections import Counter
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score
p=Path.home()/'aipharma_holdout_20261016';base=p/'endpoint_audit';datasets=json.loads((base/'datasets_all.json').read_text());protocol=json.loads((base/'protocol.json').read_text());fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
# Reuse unchanged fixed-model evaluation implementation from the previous run.
source=(p/'compare_alk5_qsar.py').read_text();core=source[source.index('# Union documents AND scaffolds'):]
combined={}
for stratum in ['enzyme_catalytic__no_explicit_mutant','cell_signaling__no_explicit_mutant']:
 rows=[r for r in datasets[stratum] if not r['old_holdout_scaffold']];o=base/stratum;o.mkdir(exist_ok=True);(o/'dataset.json').write_text(json.dumps(rows,indent=2));(o/'protocol.json').write_text(json.dumps(protocol,indent=2));N=len(rows);fps=[fp.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in rows];X=np.array([np.asarray(x) for x in fps]);y=np.array([r['label'] for r in rows]);print('STRATUM '+stratum,flush=True)
 if min(Counter(y).values())<5 or len(set(y))<2:combined[stratum]={'status':'insufficient_classes','n':N};continue
 exec(compile(core,'unchanged_fixed_evaluator','exec'))
 combined[stratum]=summary
(base/'comparison_summary.json').write_text(json.dumps(combined,indent=2))
# Verify input/output membership, source separation and every reported AUC.
old=json.loads((p/'locked_bounded/label_key.json').read_text());oldsc={r['scaffold'] for r in old.values()};checked=0
for stratum in combined:
 out=base/stratum
 if not (out/'predictions.json').exists():continue
 rows=json.loads((out/'dataset.json').read_text());pred=json.loads((out/'predictions.json').read_text());met=combined[stratum];assert len({r['smiles'] for r in rows})==len(rows);assert not {r['scaffold'] for r in rows}&oldsc
 for r in rows:
  for aid in r['assay_ids']:
   c=next(a for a in json.loads((base/'assay_classification.json').read_text()) if a['assay_id']==aid);assert c['endpoint']+'__'+c['mutation']==stratum
 for split,res in met['evaluations'].items():
  if res['status']!='complete':continue
  assert sum(z['test'] for z in res['folds'])==len(rows);assert all(z['scaffold_overlap']==0 for z in res['folds'])
  if split=='document_and_scaffold':assert all(z['document_overlap']==0 for z in res['folds'])
  for method,m in res['methods'].items():
   z=[r for r in pred if r['split']==split and r['method']==method];assert len(z)==len(rows) and len({r['smiles'] for r in z})==len(rows);assert abs(roc_auc_score([r['label'] for r in z],[r['score'] for r in z])-m['roc_auc'])<1e-12;checked+=1
(base/'verification.json').write_text(json.dumps({'checks_passed':True,'recomputed_AUC_count':checked,'old_holdout_scaffold_overlap':0,'endpoint_and_mutation_cross_pooling':False,'classification_sha256':hashlib.sha256((base/'assay_classification.json').read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((base/'protocol.json').read_bytes()).hexdigest(),'evaluation_core_sha256':hashlib.sha256(core.encode()).hexdigest(),'limitations':'Assay descriptions not complete original protocols. Per-fold class scarcity and model-fitting uncertainty remain.'},indent=2));print('VERIFIED',checked,flush=True)
