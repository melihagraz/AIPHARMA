"""Frozen nearest-reference potency baseline on source/scaffold-connected holdouts."""
import json,sys,hashlib,math,statistics,random
from pathlib import Path
from collections import defaultdict
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
p=Path(sys.argv[1]); out=p/'continuous_transfer'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(method='Nearest Morgan radius2 2048 measured reference; no fitted weights or tuning',
    split='Connected components joining document ID, normalized compound identity or exact Murcko scaffold',
    data='All exact catalytic IC50 values, including intermediate potency; repeated assay/structure median',
    baseline='Training-set median pIC50', metrics=['component-balanced MAE','pooled MAE','RMSE'],
    engineering_acceptance=dict(max_component_mean_MAE=.7,min_component_mean_MAE_gain=.1,bootstrap_gain_lower_bound_above_zero=True),
    bootstrap='10000 component resamples, seed 20261019; fixed predictions, not refitted models',
    independence='Source/scaffold-separated retrospective evaluation; dataset already inspected, not prospective or blind',
    caveat='Measured IC50 varies with assay conditions; exact scaffold separation does not remove all chemical similarity. Acceptance limits are chosen engineering goals, not validated biological error margins.',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
rows=json.loads((p/'local_sar/measured_records.json').read_text()); n=len(rows)
roots=list(range(n))
def find(i):
    while roots[i]!=i: roots[i]=roots[roots[i]]; i=roots[i]
    return i
def join(i,j): roots[find(j)]=find(i)
seen={}
for i,r in enumerate(rows):
    keys=[('identity',r['normalized_parent']),('scaffold',r['scaffold'])]+[('document',m['document_id']) for m in r['measurements'] if m['document_id']]
    for k in keys:
        if k in seen: join(i,seen[k])
        else: seen[k]=i
groups=defaultdict(list)
for i in range(n): groups[find(i)].append(i)
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
fps=[gen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in rows]
predictions=[]; components=[]
for number,ids in enumerate(sorted(groups.values(),key=lambda ids:min(ids))):
    held=set(ids); train=[i for i in range(n) if i not in held]
    if not train: continue
    baseline=statistics.median(rows[i]['median_pic50'] for i in train)
    errs=[]; baseerrs=[]
    for i in ids:
        sims=DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in train]); best=max(range(len(train)),key=lambda j:sims[j]); j=train[best]
        actual=rows[i]['median_pic50']; prediction=rows[j]['median_pic50']
        error=abs(prediction-actual); be=abs(baseline-actual); errs.append(error); baseerrs.append(be)
        predictions.append(dict(component=number,assay_id=rows[i]['assay_id'],normalized_parent=rows[i]['normalized_parent'],actual_pic50=actual,nearest_prediction=prediction,train_median_prediction=baseline,nearest_similarity=sims[best],nearest_normalized_parent=rows[j]['normalized_parent']))
    components.append(dict(component=number,count=len(ids),nearest_MAE=statistics.mean(errs),median_MAE=statistics.mean(baseerrs)))
save('predictions.json',predictions); save('components.json',components)
gains=[c['median_MAE']-c['nearest_MAE'] for c in components]
rng=random.Random(20261019)
boot=sorted(statistics.mean(rng.choices(gains,k=len(gains))) for _ in range(10000)) if gains else []
mae=statistics.mean(c['nearest_MAE'] for c in components) if components else None
gain=statistics.mean(gains) if gains else None
ci=[boot[250],boot[9749]] if boot else None
passed=bool(mae is not None and mae<=.7 and gain>=.1 and ci[0]>0)
summary=dict(records=n,components=len(components),evaluated=len(predictions),component_mean_nearest_MAE=mae,
    component_mean_median_MAE=statistics.mean(c['median_MAE'] for c in components) if components else None,
    component_mean_gain=gain,fixed_prediction_component_bootstrap_gain_95=ci,
    pooled_nearest_MAE=statistics.mean(abs(r['nearest_prediction']-r['actual_pic50']) for r in predictions) if predictions else None,
    pooled_nearest_RMSE=math.sqrt(statistics.mean((r['nearest_prediction']-r['actual_pic50'])**2 for r in predictions)) if predictions else None,
    engineering_acceptance_met=passed,candidate_activity_validated=False)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
