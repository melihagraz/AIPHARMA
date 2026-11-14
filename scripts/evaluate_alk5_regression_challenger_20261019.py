"""One frozen regression challenger; all source/scaffold holdouts retained."""
import hashlib,json,math,random,statistics,sys,time
from collections import defaultdict
from pathlib import Path
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor

p=Path(sys.argv[1]); out=p/'regression_challenger'; out.mkdir(exist_ok=False)
def save(name,data): (out/name).write_text(json.dumps(data,indent=2,allow_nan=False))
parameters=dict(n_estimators=300,min_samples_leaf=3,max_features=.5,random_state=20261019,n_jobs=4)
save('protocol.json',dict(model='RandomForestRegressor',parameters=parameters,
    input='Morgan radius2 2048 bits',endpoint='Exact catalytic IC50, continuous pIC50, intermediates retained',
    splits='Reuse fixed 24 document/identity/scaffold components; no training identity appears in held component',
    training_weight='Equal total weight per normalized identity within each training split; repeated assay records do not multiply compound weight',
    comparison='Same split training median baseline and nearest-reference baseline',
    gates=dict(max_component_mean_MAE=.7,min_component_mean_gain=.1,component_bootstrap_gain_lower_bound_positive=True),
    candidate_prediction='Only if every gate passes; no candidate-based hyperparameter tuning',
    status='Exploratory new model after previous baseline failure; reused validation set, not independent confirmation',
    model_selection='Single configuration; preserve failure, no parameter search on same benchmark',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
rows=json.loads((p/'local_sar/measured_records.json').read_text())
old=json.loads((p/'continuous_transfer/predictions.json').read_text())
assignment={(r['assay_id'],r['normalized_parent']):r['component'] for r in old}
groups=np.array([assignment[r['assay_id'],r['normalized_parent']] for r in rows])
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
X=np.zeros((len(rows),2048),dtype=np.uint8)
for i,r in enumerate(rows): DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])),X[i])
y=np.array([r['median_pic50'] for r in rows]); saved=[]; components=[]
start=time.perf_counter()
for group in sorted(set(groups)):
    train=np.flatnonzero(groups!=group); held=np.flatnonzero(groups==group)
    assert set(rows[i]['normalized_parent'] for i in train).isdisjoint(rows[i]['normalized_parent'] for i in held)
    counts=defaultdict(int)
    for i in train: counts[rows[i]['normalized_parent']]+=1
    weights=np.array([1/counts[rows[i]['normalized_parent']] for i in train])
    model=RandomForestRegressor(**parameters).fit(X[train],y[train],sample_weight=weights)
    pred=model.predict(X[held]); baseline=float(np.median(y[train]))
    for i,z in zip(held,pred): saved.append(dict(component=int(group),assay_id=rows[i]['assay_id'],normalized_parent=rows[i]['normalized_parent'],actual_pic50=float(y[i]),prediction=float(z),baseline=baseline))
    c=dict(component=int(group),count=len(held),model_MAE=float(np.mean(np.abs(pred-y[held]))),baseline_MAE=float(np.mean(np.abs(baseline-y[held]))))
    components.append(c)
    save('progress.json',dict(completed_components=len(components),total_components=len(set(groups)),elapsed_seconds=time.perf_counter()-start))
    print(json.dumps(c),flush=True)
gain=[c['baseline_MAE']-c['model_MAE'] for c in components]
rng=random.Random(20261019)
boot=sorted(statistics.mean(rng.choices(gain,k=len(gain))) for _ in range(10000))
mae=statistics.mean(c['model_MAE'] for c in components); mean_gain=statistics.mean(gain); ci=[boot[250],boot[9749]]
passed=mae<=.7 and mean_gain>=.1 and ci[0]>0
save('predictions.json',saved); save('components.json',components)
summary=dict(records=len(saved),components=len(components),component_mean_MAE=mae,component_mean_baseline_MAE=statistics.mean(c['baseline_MAE'] for c in components),component_mean_gain=mean_gain,fixed_prediction_component_bootstrap_gain_95=ci,
    pooled_MAE=statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in saved),
    pooled_RMSE=math.sqrt(statistics.mean((r['prediction']-r['actual_pic50'])**2 for r in saved)),
    engineering_acceptance_met=passed,candidate_predictions_permitted_by_protocol=passed,
    independent_validation=False,activity_validated=False,runtime_seconds=time.perf_counter()-start)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
