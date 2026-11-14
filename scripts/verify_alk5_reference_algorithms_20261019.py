"""Split/arithmetic audit and independent joint-Gaussian posterior check."""
import json,statistics,sys
from pathlib import Path
import numpy as np
root=Path(sys.argv[1]); p=Path(sys.argv[2]); sys.path.insert(0,str(root/'backend'))
from app.services.alk5_delta_model import ReferenceDeltaModel
from app.services.alk5_reference_gp import AssayReferenceGP
def read(f): return json.loads((p/f).read_text())
rows=read('local_sar/measured_records.json'); splits=read('reference_delta/splits.json'); pairs=read('reference_delta/pair_indices.json')
assign={(r['assay_id'],r['normalized_parent']):r['component'] for r in read('continuous_transfer/predictions.json')}
for r in pairs:
    i,j=r['first'],r['second']
    assert rows[i]['assay_id']==rows[j]['assay_id']
    assert assign[rows[i]['assay_id'],rows[i]['normalized_parent']]==r['component']==assign[rows[j]['assay_id'],rows[j]['normalized_parent']]
bygroup={r['component']:set(r['training_indices']) for r in splits}
for stage in ['reference_delta','reference_gp']:
    records=read(stage+'/predictions.json'); s=read(stage+'/summary.json')
    assert len(records)==501
    for r in records:
        train=bygroup[r['component']]
        assert r['index'] not in train and r['anchor_index'] not in train
        assert all(rows[i]['assay_id']!=r['assay_id'] for i in train)
        assert r['normalized_parent']==rows[r['index']]['normalized_parent']
    for label,subset in [('all_queries',records),('non_abstained',[r for r in records if not r['abstain']])]:
        groups={g:[r for r in subset if r['component']==g] for g in sorted({r['component'] for r in subset})}
        mae=statistics.mean(statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in rs) for rs in groups.values())
        baseline=statistics.mean(statistics.mean(abs(r['anchor_pic50']-r['actual_pic50']) for r in rs) for rs in groups.values())
        assert abs(mae-s[label]['component_mean_MAE'])<1e-12
        assert abs(baseline-s[label]['component_mean_baseline_MAE'])<1e-12
    coverage=sum(not r['abstain'] for r in records)/len(records)
    assert abs(coverage-s['coverage'])<1e-12
    m=s['non_abstained']; passed=coverage>=.5 and m['component_mean_MAE']<=.7 and m['component_mean_gain']>=.1 and m['bootstrap_gain_95'][0]>0
    assert passed==s['engineering_acceptance_met']
keys=lambda stage:{(r['index'],r['anchor_index']) for r in read(stage+'/predictions.json')}
assert keys('reference_delta')==keys('reference_gp')

# Independent posterior calculation by a single joint solve, rather than the
# implementation's two-step conditional update.
bits=np.array([[1,1,0,0],[1,0,1,0],[0,0,1,1]],dtype=float)
y=np.array([6.,6.7,8.1]); assays=np.array(['a','a','b'])
anchor=np.array([1,1,1,0],dtype=float); query=np.array([1,1,0,1],dtype=float); actual_anchor=7.3
gp=AssayReferenceGP().fit(bits,y,assays)
pred=gp.predict_with_reference(query,anchor,actual_anchor,'new')
joint=np.vstack([bits,anchor]); labels=list(assays)+['new']
def kernel(a,b):
    intersection=float(np.minimum(a,b).sum()); union=float(np.maximum(a,b).sum())
    return intersection/union
cov=np.array([[kernel(a,b)+(labels[i]==labels[j]) for j,b in enumerate(joint)] for i,a in enumerate(joint)])
cov+=np.eye(4)*.04
kq=np.array([kernel(query,b)+(label=='new') for b,label in zip(joint,labels)])
mean=float(np.median(y)); expected=mean+float(kq@np.linalg.solve(cov,np.r_[y,actual_anchor]-mean))
sd=float(np.sqrt(2.04-kq@np.linalg.solve(cov,kq)))
assert abs(expected-pred['prediction'])<1e-10 and abs(sd-pred['standard_deviation'])<1e-10
try: gp.predict_with_reference(query,anchor,actual_anchor,'a')
except ValueError: pass
else: raise AssertionError('Known training assay was not rejected')
try: ReferenceDeltaModel().fit(np.eye(3),np.array([1.,2.,3.]),[(0,2,0)],{0,1})
except ValueError: pass
else: raise AssertionError('Forbidden training endpoint was not rejected')
gp_rows=read('reference_gp/predictions.json')
result=dict(status='passed',checks=['all fitted pair endpoints share eligible assay and component','501 matched targets in both models','held-out assay and reference excluded from fitted training sets','component-balanced error, baseline and rejection decisions recomputed','GP posterior mean and variance agree with independent joint covariance solve','training-assay conditioning rejected','forbidden pair endpoint rejected'],posterior_mean_absolute_error=abs(expected-pred['prediction']),posterior_sd_absolute_error=abs(sd-pred['standard_deviation']),mean_model_interval_width_pic50=statistics.mean(r['interval_95'][1]-r['interval_95'][0] for r in gp_rows),mean_non_abstained_interval_width_pic50=statistics.mean(r['interval_95'][1]-r['interval_95'][0] for r in gp_rows if not r['abstain']),independent_biological_validation=False)
(p/'reference_algorithm_verification.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
