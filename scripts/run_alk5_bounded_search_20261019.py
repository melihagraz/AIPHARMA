"""Finite sequential search; numerical objectives are engineering targets, not safety gates."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('root', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
sys.path.insert(0, str(a.root/'backend'))
from app.services.alk5_analogue_design import design_analogues

a.output.mkdir(parents=True, exist_ok=False)
def save(name, data):
    (a.output/name).write_text(json.dumps(data, indent=2, allow_nan=False))

protocol = dict(max_candidates=100, batch_sizes=[30,30,40],
    objectives=dict(Solubility_AqSolDB_min_delta=.3, risk_max_delta=.02,
        one_risk_min_improvement=.05, Bioavailability_Ma_min_delta=-.05,
        Caco2_Wang_min_delta=-.2),
    risk_endpoints=['DILI','hERG','AMES'],
    reserved_documents=['doi:10.1016/j.ejmech.2020.112154',
        'doi:10.1016/j.bmcl.2016.07.030', 'doi:10.1016/j.ejmech.2020.112354'],
    stop='Finite library exhausted or 100 candidates; never change targets after outcomes',
    independence='Reserved-source engineering check only: baseline parents and model previously seen; not blind or independent biological validation',
    activity_gate='No validated prospective ALK5 activity model: all designs remain unmeasured',
    caveat='Same-model numerical optimization can exploit model error; no clinical interpretation of thresholds',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
save('protocol.json',protocol)
design=design_analogues(100)
save('input_designs.json',design)
reserved=set(protocol['reserved_documents'])
dev=[]
hold=[]
for c in design['candidates']:
    (hold if reserved.intersection(c['provenance']['parent_documents']) else dev).append(c)
save('partition.json',dict(development=[c['candidate_id'] for c in dev],reserved=[c['candidate_id'] for c in hold]))

import torch
torch.set_num_threads(4)
from admet_ai import ADMETModel
model=ADMETModel(num_workers=0)
cache={}
endpoints=protocol['risk_endpoints']+['Solubility_AqSolDB','Bioavailability_Ma','Caco2_Wang']
all_rows=[]
def evaluate(batch,name):
    inputs=sorted({s for c in batch for s in [c['smiles'],c['provenance']['parent_smiles']]}-cache.keys())
    if inputs:
        df=model.predict(smiles=inputs)
        assert df.index.is_unique and set(inputs)<=set(df.index)
        df.to_csv(a.output/(name+'_all_predictions.csv'),float_format='%.17g')
        for s in inputs:
            cache[s]={k:float(df.loc[s,k]) for k in endpoints}
    rows=[]
    for c in batch:
        child=cache[c['smiles']]
        parent=cache[c['provenance']['parent_smiles']]
        delta={k:child[k]-parent[k] for k in endpoints}
        checks=dict(solubility=delta['Solubility_AqSolDB']>=.3,
            risk_nonworsening=all(delta[k]<=.02 for k in protocol['risk_endpoints']),
            one_risk_improves=any(delta[k]<=-.05 for k in protocol['risk_endpoints']),
            bioavailability=delta['Bioavailability_Ma']>=-.05,
            permeability=delta['Caco2_Wang']>=-.2)
        rows.append(dict(candidate_id=c['candidate_id'],smiles=c['smiles'],
            parent_id=c['provenance']['parent_id'],parent=parent,prediction=child,
            delta=delta,checks=checks,engineering_objectives_met=all(checks.values()),
            activity_validated=False,safety_validated=False))
    save(name+'.json',rows)
    all_rows.extend(rows)
    print(json.dumps(dict(stage=name,count=len(rows),objectives_met=sum(r['engineering_objectives_met'] for r in rows))),flush=True)

for i in range(0,len(dev),30):
    evaluate(dev[i:i+30],f'development_{i//30+1}')
evaluate(hold,'reserved_final')
save('summary.json',dict(total=len(all_rows),development=len(dev),reserved=len(hold),
    engineering_hits=[r['candidate_id'] for r in all_rows if r['engineering_objectives_met']],
    biological_validation=False,safety_validation=False,
    conclusion='Numerical model-screen results only; target activity remains unverified'))
print((a.output/'summary.json').read_text(),flush=True)
