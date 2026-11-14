"""Nested source-component selection among 14 frozen reference-assisted options."""
import hashlib,json,random,statistics,sys,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from sklearn.linear_model import Ridge
from threadpoolctl import threadpool_limits
threadpool_limits(2)
root=Path(sys.argv[1]); p=Path(sys.argv[2]); out=p/'nested_method_search'; out.mkdir(exist_ok=False)
sys.path.insert(0,str(root/'backend'))
from app.services.alk5_delta_model import molecular_features
from app.services.alk5_reference_gp import AssayReferenceGP
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
configs=[dict(name='unchanged_reference',kind='baseline')]
configs += [dict(name=f'delta_ridge_alpha_{a}',kind='ridge',alpha=a) for a in [1,10,100,1000]]
configs += [dict(name=f'reference_gp_assay_{a}_noise_{n}',kind='gp',assay_variance=a,noise_variance=n) for a in [.25,1.,4.] for n in [.04,.25,1.]]
assert len(configs)==14
save('protocol.json',dict(options=configs,task='Same 501 measured-reference-assisted queries, not zero-shot prediction',
    outer='24 fixed document/identity/scaffold components, held out one at a time',
    inner='Remaining components assigned deterministically modulo3; choose minimum component-balanced MAE, ties simplest option first',
    selection='Outer query labels never used for choice; all 14 outer outcomes retained as exploratory comparisons, not independent evidence for their posthoc winner',
    test='Evaluate selected option on held outer component and compare unchanged measured reference on exactly the same queries',
    anchors='Fixed label-free structural medoids, supplied measurements only at prediction; neither query nor reference is fitted in its validation split',
    score='All queries retained; baseline fallback participates explicitly; this search does not establish confidence or candidate applicability',
    acceptance='Nested component-mean MAE<=0.7, gain>=0.1 and fixed-prediction component-bootstrap95 gain lower bound>0',
    limits='Repeatedly inspected retrospective development data; not blinded independent confirmation; no candidate activity/safety claims',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
rows=json.loads((p/'local_sar/measured_records.json').read_text())
queries=json.loads((p/'reference_delta/predictions.json').read_text())
old=json.loads((p/'continuous_transfer/predictions.json').read_text())
assignment={(r['assay_id'],r['normalized_parent']):r['component'] for r in old}
groups=np.array([assignment[r['assay_id'],r['normalized_parent']] for r in rows])
X=molecular_features([r['smiles'] for r in rows]); y=np.array([r['median_pic50'] for r in rows]); assays=np.array([r['assay_id'] for r in rows])
pairs=json.loads((p/'reference_delta/pair_indices.json').read_text())
allgroups=sorted(set(groups.tolist())); started=time.perf_counter()
def evaluate(config,train_groups,test_queries):
    train=np.flatnonzero(np.isin(groups,list(train_groups)))
    assert all(r['component'] not in train_groups and r['index'] not in train and r['anchor_index'] not in train for r in test_queries)
    if config['kind']=='baseline': return [r['anchor_pic50'] for r in test_queries]
    if config['kind']=='ridge':
        ps=[r for r in pairs if r['component'] in train_groups]
        counts=Counter(r['component'] for r in ps)
        dx=np.array([X[r['second']]-X[r['first']] for r in ps])
        dy=np.array([y[r['second']]-y[r['first']] for r in ps])
        weights=[len(ps)/(len(counts)*counts[r['component']]) for r in ps]
        model=Ridge(alpha=config['alpha'],fit_intercept=False,solver='lsqr',tol=1e-5).fit(dx,dy,sample_weight=weights)
        differences=np.array([X[r['index']]-X[r['anchor_index']] for r in test_queries])
        if not len(differences): return []
        return (np.array([r['anchor_pic50'] for r in test_queries])+model.predict(differences)).tolist()
    model=AssayReferenceGP(); model.assay_variance=config['assay_variance']; model.noise_variance=config['noise_variance']
    model.fit(X[train,:2048],y[train],assays[train])
    return [model.predict_with_reference(X[r['index'],:2048],X[r['anchor_index'],:2048],r['anchor_pic50'],r['assay_id'])['prediction'] for r in test_queries]
def group_errors(qs,preds):
    gs=defaultdict(list)
    for r,pred in zip(qs,preds): gs[r['component']].append(abs(pred-r['actual_pic50']))
    return {g:statistics.mean(v) for g,v in gs.items()}
outer_results=[]; selections=[]; leaderboard=defaultdict(list)
for held in allgroups:
    test=[r for r in queries if r['component']==held]
    if not test: continue
    available=set(allgroups)-{held}; inner_scores=[]
    for config in configs:
        errs={}
        for fold in range(3):
            validation={g for g in available if g%3==fold}; training=available-validation
            val=[r for r in queries if r['component'] in validation]
            errs.update(group_errors(val,evaluate(config,training,val)))
        inner_scores.append(statistics.mean(errs.values()))
    best=min(range(len(configs)),key=lambda j:(inner_scores[j],j)); chosen=configs[best]
    outer_predictions={}
    for config in configs:
        pred=evaluate(config,available,test); outer_predictions[config['name']]=pred
        leaderboard[config['name']].append(dict(component=held,MAE=statistics.mean(abs(z-r['actual_pic50']) for z,r in zip(pred,test))))
    for r,pred in zip(test,outer_predictions[chosen['name']]):
        outer_results.append(dict(component=held,index=r['index'],anchor_index=r['anchor_index'],actual_pic50=r['actual_pic50'],anchor_pic50=r['anchor_pic50'],prediction=pred,selected_method=chosen['name']))
    selection=dict(held_component=held,selected=chosen['name'],inner_scores={c['name']:v for c,v in zip(configs,inner_scores)},outer_all_method_predictions=outer_predictions,outer_query_indices=[r['index'] for r in test])
    selections.append(selection); save(f'component_{held:02d}.json',selection)
    print(json.dumps(dict(component=held,selected=chosen['name'],evaluated=len(outer_results),elapsed_seconds=time.perf_counter()-started)),flush=True)
save('nested_predictions.json',outer_results)
cs=[]
for g in sorted({r['component'] for r in outer_results}):
    rs=[r for r in outer_results if r['component']==g]
    cs.append(dict(component=g,count=len(rs),MAE=statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in rs),baseline_MAE=statistics.mean(abs(r['anchor_pic50']-r['actual_pic50']) for r in rs)))
gains=[r['baseline_MAE']-r['MAE'] for r in cs]; rng=random.Random(20261019)
boot=sorted(statistics.mean(rng.choices(gains,k=len(gains))) for _ in range(10000))
mae=statistics.mean(r['MAE'] for r in cs); gain=statistics.mean(gains)
ranked=sorted([dict(method=k,component_mean_MAE=statistics.mean(r['MAE'] for r in v),component_results=v) for k,v in leaderboard.items()],key=lambda r:r['component_mean_MAE'])
summary=dict(queries=len(outer_results),outer_evaluated_components=len(cs),options=len(configs),selection_counts=dict(Counter(r['selected'] for r in selections)),nested_component_MAE=mae,baseline_component_MAE=statistics.mean(r['baseline_MAE'] for r in cs),nested_gain=gain,conditional_component_bootstrap_gain_95=[boot[250],boot[9749]],engineering_acceptance_met=bool(mae<=.7 and gain>=.1 and boot[250]>0),exploratory_configuration_leaderboard=ranked,independent_validation=False,candidate_activity_validated=False)
save('summary.json',summary); print(json.dumps({k:v for k,v in summary.items() if k!='exploratory_configuration_leaderboard'},indent=2),flush=True)
