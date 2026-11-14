"""Fixed leave-one-identity-out prediction with multiple same-assay analogues."""
import hashlib,json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from threadpoolctl import threadpool_limits
threadpool_limits(2)
p=Path(sys.argv[1]); out=p/'within_series'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(task='Multiple measured same-assay same-Murcko analogues available; leave each identity out. Different from earlier single-anchor unseen-series task.',
    eligibility='Exact catalytic IC50 only; assay plus exact normalized Murcko scaffold series, >=8 unique identities. No potency-range or outcome-based inclusion.',
    input='Morgan radius2 2048 bit Tanimoto',model='Kernel ridge, fixed lambda 0.25; centered on median of training analogues only; no tuning',
    comparison='Nearest measured analogue unchanged; ties resolved by normalized structure lexical order. Also report structural-medoid unchanged.',
    subset='Similarity >=0.65 to nearest measured training analogue, fixed before evaluation; retain all-query results too',
    dependence='Union series sharing publication or chemical identity before component-balanced errors and component bootstrap',
    acceptance='Supported queries >=30, >=5 dependence components, component mean MAE<=0.7, gain over nearest unchanged>=0.1 and conditional bootstrap lower bound>0',
    interpretation='Retrospective development, previously inspected data; no external confirmation and no candidate activity/safety evidence',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
rows=json.loads((p/'curated_complete_identity/records.json').read_text()); series=defaultdict(list)
for r in rows:
    if r['median_pic50'] is not None: series[r['assay_id'],r['scaffold']].append(r)
series={k:sorted(v,key=lambda r:r['normalized_parent']) for k,v in series.items() if len(v)>=8}
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
predictions=[]; summaries=[]
for sid,(key,rs) in enumerate(sorted(series.items())):
    assert len({r['normalized_parent'] for r in rs})==len(rs)
    fps=[gen.GetFingerprint(Chem.MolFromSmiles(r['normalized_parent'])) for r in rs]
    K=np.array([DataStructs.BulkTanimotoSimilarity(f,fps) for f in fps]); y=np.array([r['median_pic50'] for r in rs])
    for i,r in enumerate(rs):
        ix=np.array([j for j in range(len(rs)) if j!=i]); assert r['normalized_parent'] not in {rs[j]['normalized_parent'] for j in ix}
        nearest=int(ix[np.argmax(K[i,ix])]); medoid=int(ix[np.argmax(K[np.ix_(ix,ix)].mean(axis=1))])
        center=float(np.median(y[ix])); alpha=np.linalg.solve(K[np.ix_(ix,ix)]+.25*np.eye(len(ix)),y[ix]-center)
        prediction=float(center+K[i,ix]@alpha)
        predictions.append(dict(series=sid,assay_id=key[0],scaffold=key[1],normalized_parent=r['normalized_parent'],documents=r['documents'],
            actual_pic50=float(y[i]),prediction=prediction,nearest_unchanged=float(y[nearest]),medoid_unchanged=float(y[medoid]),nearest_similarity=float(K[i,nearest]),
            nearest_identity=rs[nearest]['normalized_parent'],training_analogues=len(ix),supported=bool(K[i,nearest]>=.65)))
    summaries.append(dict(series=sid,assay_id=key[0],scaffold=key[1],records=len(rs),documents=sorted({d for r in rs for d in r['documents']})))
# Whole-source/identity components for uncertainty; series folds do not imply independent sources.
parents=list(range(len(summaries)))
def find(i):
    while parents[i]!=i: parents[i]=parents[parents[i]]; i=parents[i]
    return i
seen={}
for r in predictions:
    for key in [('identity',r['normalized_parent'])]+[('document',d) for d in r['documents']]:
        if key in seen: parents[find(r['series'])]=find(seen[key])
        else: seen[key]=r['series']
for r in predictions: r['component']=find(r['series'])
rng=np.random.default_rng(20261019)
def metrics(rs):
    by=defaultdict(list)
    for r in rs: by[r['component']].append(r)
    cs=[]
    for c,items in sorted(by.items()):
        e={k:float(np.mean([abs(r[k]-r['actual_pic50']) for r in items])) for k in ['prediction','nearest_unchanged','medoid_unchanged']}
        cs.append(dict(component=c,queries=len(items),model_MAE=e['prediction'],nearest_MAE=e['nearest_unchanged'],medoid_MAE=e['medoid_unchanged'],gain=e['nearest_unchanged']-e['prediction']))
    if not cs: return dict(queries=0,components=0,acceptance=False)
    gains=np.array([r['gain'] for r in cs]); boot=np.mean(rng.choice(gains,size=(10000,len(gains))),axis=1)
    ci=np.quantile(boot,[.025,.975]).tolist(); mae=float(np.mean([r['model_MAE'] for r in cs])); gain=float(gains.mean())
    return dict(queries=len(rs),components=len(cs),component_mean_MAE=mae,nearest_component_MAE=float(np.mean([r['nearest_MAE'] for r in cs])),
        medoid_component_MAE=float(np.mean([r['medoid_MAE'] for r in cs])),gain_over_nearest=gain,conditional_component_bootstrap_gain_95=ci,
        acceptance=bool(len(rs)>=30 and len(cs)>=5 and mae<=.7 and gain>=.1 and ci[0]>0),component_details=cs)
summary=dict(series=len(series),all_queries=metrics(predictions),supported_queries=metrics([r for r in predictions if r['supported']]),
    independent_validation=False,biological_activity_validated=False,tool_changed=False)
save('predictions.json',predictions); save('series.json',summaries); save('summary.json',summary)
print(json.dumps({**summary,'all_queries':{k:v for k,v in summary['all_queries'].items() if k!='component_details'},'supported_queries':{k:v for k,v in summary['supported_queries'].items() if k!='component_details'}},indent=2))
