"""Evaluate fixed assay-offset GP on the exact same measured-reference task."""
import hashlib,json,random,statistics,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from threadpoolctl import threadpool_limits
threadpool_limits(4)
root=Path(sys.argv[1]); report=Path(sys.argv[2]); out=report/'reference_gp'; out.mkdir(exist_ok=False)
sys.path.insert(0,str(root/'backend'))
from app.services.alk5_delta_model import molecular_features
from app.services.alk5_reference_gp import AssayReferenceGP
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(model='Tanimoto GP plus assay random offset, conditioned on one measured medoid anchor',
    molecular_variance=1.,assay_variance=1.,measurement_noise_variance=.04,
    optimizer='None; single fixed configuration',
    task='Same 501 anchor-assisted queries and source-excluded training sets as reference_delta',
    anchors='Same label-free structural medoids; their measured labels enter conditioning only',
    abstention='Anchor similarity<0.65 or model predictive SD>0.7 pIC50',
    acceptance='Coverage>=0.5, component-balanced MAE<=0.7, gain over unchanged reference>=0.1 and fixed-prediction bootstrap95 lower gain>0; interval coverage also reported',
    caveat='Development after previous failures on already inspected data; not independent confirmation. Gaussian uncertainty may be miscalibrated; fixed noise is not measured experimental variance.',
    candidate_use='No activity ranking unless all gates pass',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model_sha256=hashlib.sha256((root/'backend/app/services/alk5_reference_gp.py').read_bytes()).hexdigest()))
rows=json.loads((report/'local_sar/measured_records.json').read_text())
old=json.loads((report/'reference_delta/predictions.json').read_text())
splits=json.loads((report/'reference_delta/splits.json').read_text())
X=molecular_features([r['smiles'] for r in rows])[:,:2048]; y=np.array([r['median_pic50'] for r in rows]); assays=np.array([r['assay_id'] for r in rows])
pred=[]
for split in splits:
    train=split['training_indices']; model=AssayReferenceGP().fit(X[train],y[train],assays[train])
    for r in [r for r in old if r['component']==split['component']]:
        i=r['index']; anchor=r['anchor_index']
        assert i not in train and anchor not in train and assays[i] not in set(assays[train])
        estimate=model.predict_with_reference(X[i],X[anchor],float(y[anchor]),assays[i])
        pred.append(dict(component=r['component'],index=i,anchor_index=anchor,assay_id=r['assay_id'],normalized_parent=r['normalized_parent'],actual_pic50=float(y[i]),anchor_pic50=float(y[anchor]),**estimate))
    print(json.dumps(dict(component=split['component'],evaluated=len(pred))),flush=True)
assert {(r['index'],r['anchor_index']) for r in pred}=={(r['index'],r['anchor_index']) for r in old}
def metrics(rs):
    groups=defaultdict(list)
    for r in rs: groups[r['component']].append(r)
    cs=[dict(component=g,count=len(v),MAE=statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in v),baseline_MAE=statistics.mean(abs(r['anchor_pic50']-r['actual_pic50']) for r in v)) for g,v in sorted(groups.items())]
    if not cs: return dict(queries=0,components=0,component_mean_MAE=None,component_mean_gain=None,bootstrap_gain_95=None)
    gains=[c['baseline_MAE']-c['MAE'] for c in cs]; rng=random.Random(20261019)
    bootstrap=sorted(statistics.mean(rng.choices(gains,k=len(gains))) for _ in range(10000))
    return dict(queries=len(rs),components=len(cs),component_mean_MAE=statistics.mean(c['MAE'] for c in cs),component_mean_baseline_MAE=statistics.mean(c['baseline_MAE'] for c in cs),component_mean_gain=statistics.mean(gains),bootstrap_gain_95=[bootstrap[250],bootstrap[9749]],interval_95_coverage=statistics.mean(r['interval_95'][0]<=r['actual_pic50']<=r['interval_95'][1] for r in rs),component_details=cs)
allm=metrics(pred); accepted=metrics([r for r in pred if not r['abstain']]); coverage=accepted['queries']/len(pred)
passed=bool(coverage>=.5 and accepted['component_mean_MAE'] is not None and accepted['component_mean_MAE']<=.7 and accepted['component_mean_gain']>=.1 and accepted['bootstrap_gain_95'][0]>0)
summary=dict(all_queries=allm,non_abstained=accepted,coverage=coverage,engineering_acceptance_met=passed,independent_validation=False,activity_validated=False,safety_validated=False)
save('predictions.json',pred); save('summary.json',summary)
print(json.dumps({k:v for k,v in summary.items() if k not in ['all_queries','non_abstained']},indent=2),flush=True)
print(json.dumps(dict(all_MAE=allm['component_mean_MAE'],all_baseline=allm['component_mean_baseline_MAE'],non_abstained_MAE=accepted['component_mean_MAE'],non_abstained_baseline=accepted.get('component_mean_baseline_MAE'),non_abstained_gain=accepted['component_mean_gain']),indent=2),flush=True)
