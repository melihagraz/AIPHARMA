"""Source-excluded evaluation with one measured structural medoid anchor per assay."""
import hashlib,json,random,statistics,sys,time
from collections import defaultdict
from pathlib import Path
import numpy as np

root=Path(sys.argv[1]); report=Path(sys.argv[2]); out=report/'reference_delta'
out.mkdir(exist_ok=False)
sys.path.insert(0,str(root/'backend'))
from app.services.alk5_delta_model import ReferenceDeltaModel,molecular_features
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(task='Predict activity CHANGE given an experimentally measured same-assay anchor; distinct from zero-shot absolute potency prediction',
    model='Ridge alpha10 no intercept; within-assay differences of Morgan2/2048 plus seven scaled physicochemical descriptors; equal total source-component weights',
    split='24 fixed document/identity/scaffold components; all fitted pair endpoints outside held component',
    anchor='One structural medoid per held assay, ties normalized identity; chosen without any potency labels; anchor potency supplied only at prediction',
    baseline='Measured anchor potency unchanged, evaluated on exactly the same query rows',
    pairs='Same assay, Morgan similarity >=0.5; <=200 pairs per assay selected by SHA256 identities without outcome access',
    abstention='Anchor similarity <0.65 or edit cosine >=0.3 supported by fewer than 3 training source components',
    acceptance=dict(coverage_at_least=.5,component_mean_MAE_at_most=.7,component_mean_gain_at_least=.1,component_bootstrap_gain_lower_bound_positive=True),
    constraints='No threshold or model search after results; report both all-query and non-abstained error; no discarded failures',
    independence='Retrospective development on previously inspected data; anchors intentionally share assay/source with queries but never train the difference model. Not blinded or prospective validation.',
    candidate_use='No candidate activity ranking unless acceptance passes; even a pass does not establish candidate biology',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    model_sha256=hashlib.sha256((root/'backend/app/services/alk5_delta_model.py').read_bytes()).hexdigest()))
records=json.loads((report/'local_sar/measured_records.json').read_text())
old=json.loads((report/'continuous_transfer/predictions.json').read_text())
assignment={(r['assay_id'],r['normalized_parent']):r['component'] for r in old}
groups=np.array([assignment[r['assay_id'],r['normalized_parent']] for r in records])
X=molecular_features([r['smiles'] for r in records]); y=np.array([r['median_pic50'] for r in records])
assays=defaultdict(list)
for i,r in enumerate(records): assays[r['assay_id']].append(i)
anchors={}; pairs=[]
for assay,ids in sorted(assays.items()):
    if len(ids)<2: continue
    bits=X[ids,:2048]; intersection=bits@bits.T; union=bits.sum(axis=1)[:,None]+bits.sum(axis=1)[None,:]-intersection
    sim=intersection/np.maximum(union,1)
    best=sorted(range(len(ids)),key=lambda j:(-float(sim[j].sum()),records[ids[j]]['normalized_parent']))[0]
    anchors[assay]=ids[best]
    eligible=[]
    for ii,i in enumerate(ids):
        for jj in range(ii+1,len(ids)):
            if sim[ii,jj]<.5: continue
            j=ids[jj]; assert groups[i]==groups[j]
            key=hashlib.sha256((records[i]['normalized_parent']+'|'+records[j]['normalized_parent']).encode()).hexdigest()
            eligible.append((key,i,j,int(groups[i])))
    pairs.extend((i,j,g) for key,i,j,g in sorted(eligible)[:200])
save('anchors.json',[dict(assay_id=a,index=i,normalized_parent=records[i]['normalized_parent'],pic50=float(y[i])) for a,i in anchors.items()])
save('pair_indices.json',[dict(first=i,second=j,component=g) for i,j,g in pairs])
results=[]; splits=[]; started=time.perf_counter()
for group in sorted(set(groups)):
    train=np.flatnonzero(groups!=group).tolist(); fit_pairs=[r for r in pairs if r[2]!=group]
    model=ReferenceDeltaModel().fit(X,y,fit_pairs,train)
    split=dict(component=int(group),training_indices=train,training_pairs=len(fit_pairs)); splits.append(split)
    for i in np.flatnonzero(groups==group):
        r=records[i]; anchor=anchors.get(r['assay_id'])
        if anchor is None or anchor==i: continue
        diagnostic=model.predict_change(X[anchor],X[i])
        results.append(dict(component=int(group),assay_id=r['assay_id'],index=int(i),normalized_parent=r['normalized_parent'],anchor_index=anchor,anchor_identity=records[anchor]['normalized_parent'],actual_pic50=float(y[i]),anchor_pic50=float(y[anchor]),prediction=float(y[anchor]+diagnostic['delta_pic50']),**diagnostic))
    print(json.dumps(dict(component=int(group),evaluated=len(results),elapsed_seconds=time.perf_counter()-started)),flush=True)
save('splits.json',splits); save('predictions.json',results)
def metrics(rows):
    gs=defaultdict(list)
    for r in rows: gs[r['component']].append(r)
    components=[dict(component=g,count=len(rs),MAE=statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in rs),baseline_MAE=statistics.mean(abs(r['anchor_pic50']-r['actual_pic50']) for r in rs)) for g,rs in sorted(gs.items())]
    if not components: return dict(queries=0,components=0,component_mean_MAE=None,component_mean_gain=None,bootstrap_gain_95=None)
    gains=[r['baseline_MAE']-r['MAE'] for r in components]; rng=random.Random(20261019)
    bootstrap=sorted(statistics.mean(rng.choices(gains,k=len(gains))) for _ in range(10000))
    return dict(queries=len(rows),components=len(components),component_mean_MAE=statistics.mean(r['MAE'] for r in components),component_mean_baseline_MAE=statistics.mean(r['baseline_MAE'] for r in components),component_mean_gain=statistics.mean(gains),bootstrap_gain_95=[bootstrap[250],bootstrap[9749]],component_details=components)
all_metrics=metrics(results); accepted=metrics([r for r in results if not r['abstain']]); coverage=accepted['queries']/len(results) if results else 0
passed=bool(coverage>=.5 and accepted['component_mean_MAE'] is not None and accepted['component_mean_MAE']<=.7 and accepted['component_mean_gain']>=.1 and accepted['bootstrap_gain_95'][0]>0)
summary=dict(training_pairs_total=len(pairs),all_queries=all_metrics,non_abstained=accepted,coverage=coverage,engineering_acceptance_met=passed,activity_validated=False,safety_validated=False,independent_validation=False)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
