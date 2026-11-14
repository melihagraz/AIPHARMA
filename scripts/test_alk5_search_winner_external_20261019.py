"""Frozen exploratory winner, evaluated on a separate already-inspected source."""
import json,sys,hashlib,statistics
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
threadpool_limits(2)
root=Path(sys.argv[1]); p=Path(sys.argv[2]); out=p/'search_winner_external_completed'; out.mkdir(exist_ok=False)
sys.path.insert(0,str(root/'backend'))
from app.services.alk5_delta_model import molecular_features
from app.services.alk5_reference_gp import AssayReferenceGP,tanimoto
search=json.loads((p/'nested_method_search/summary.json').read_text())
winner=search['exploratory_configuration_leaderboard'][0]['method']
config=next(c for c in json.loads((p/'nested_method_search/protocol.json').read_text())['options'] if c['name']==winner)
assert config['kind']=='gp'
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(frozen_model=config,selection='Lowest exploratory table MAE; nested selection did not pass acceptance, so this remains diagnostic',
    source='US20250179094A1 / BindingDB entry12979, already inspected and active-only; not blind',
    anchor='One label-free structural medoid among 51 normalized source structures',
    train='Original 541 exact catalytic assay/structure groups, with no normalized identity or exact scaffold overlap with source',
    evaluation='50 remaining structures vs unchanged measured anchor; no threshold or parameter adaptation',
    interpretation='One external source cannot establish generalization across sources; active-only cannot establish active-weak discrimination',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
rows=json.loads((p/'local_sar/measured_records.json').read_text())
external=json.loads((p/'external_evaluation_deduplicated/predictions.json').read_text())
assert not ({r['normalized_parent'] for r in rows}&{r['normalized_parent'] for r in external})
assert not ({r['scaffold'] for r in rows}&{r['scaffold'] for r in external})
X=molecular_features([r['smiles'] for r in rows])[:,:2048]
Q=molecular_features([r['smiles'] for r in external])[:,:2048]
similarity=tanimoto(Q,Q)
anchor=sorted(range(len(external)),key=lambda j:(-float(similarity[j].sum()),external[j]['normalized_parent']))[0]
save('anchor.json',external[anchor])
model=AssayReferenceGP(); model.assay_variance=config['assay_variance']; model.noise_variance=config['noise_variance']
model.fit(X,np.array([r['median_pic50'] for r in rows]),np.array([r['assay_id'] for r in rows]))
pred=[]
for i,r in enumerate(external):
    if i==anchor: continue
    estimate=model.predict_with_reference(Q[i],Q[anchor],external[anchor]['pic50'],'external_entry12979')
    pred.append(dict(normalized_parent=r['normalized_parent'],compound_number=r['compound_number'],actual_pic50=r['pic50'],anchor_pic50=external[anchor]['pic50'],**estimate))
save('predictions.json',pred)
def metrics(rs):
    if not rs: return dict(count=0)
    rho=float(spearmanr([r['actual_pic50'] for r in rs],[r['prediction'] for r in rs]).statistic) if len(rs)>2 else None
    return dict(count=len(rs),MAE=statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in rs),baseline_MAE=statistics.mean(abs(r['anchor_pic50']-r['actual_pic50']) for r in rs),Spearman=rho if rho is not None and np.isfinite(rho) else None,mean_interval_width_pic50=statistics.mean(r['interval_95'][1]-r['interval_95'][0] for r in rs))
summary=dict(model=winner,all=metrics(pred),non_abstained=metrics([r for r in pred if not r['abstain']]),coverage=sum(not r['abstain'] for r in pred)/len(pred),independent_blinded_validation=False,candidate_activity_validated=False,prior_search_acceptance_passed=False)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
