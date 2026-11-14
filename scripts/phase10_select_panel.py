#!/usr/bin/env python3
"""Frozen, descriptor-matched EGFR DUD-E selection (no docking scores)."""
import argparse, csv, hashlib, json, random
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from scipy.optimize import linear_sum_assignment

FIELDS = ['MW','logP','HBA','HBD','TPSA','RotB','charge']
SCALE = np.array([50.,0.8,2.,1.,30.,3.])
# Exact charge plus per-pair max differences; thresholds chosen before viewing outcomes.
CALIPER = np.array([65.,1.0,2.,1.,40.,4.])

def hash_file(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def desc(m):
    return np.array([Descriptors.MolWt(m), Crippen.MolLogP(m), Lipinski.NumHAcceptors(m),
                     Lipinski.NumHDonors(m), rdMolDescriptors.CalcTPSA(m),
                     Lipinski.NumRotatableBonds(m), Chem.GetFormalCharge(m)], dtype=float)
def parse(path):
    out = {}; bad = 0
    for line in Path(path).read_text(errors='replace').splitlines():
        toks = line.split()
        if not toks or toks[0].startswith('#'): continue
        mol=Chem.MolFromSmiles(toks[0]);
        if mol is None: bad+=1; continue
        canon=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if canon not in out:
            out[canon]={'smiles':canon,'source_id': toks[1] if len(toks)>1 else '', 'descs':desc(mol)}
    return out,bad

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--active',required=True,type=Path)
    p.add_argument('--decoy',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--n',type=int,default=50)
    p.add_argument('--seed',type=int,default=20261110)
    p.add_argument('--exclude',type=Path,default=None,help='CSV of previously docked molecules; excludes smiles')
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    act,bad_a=parse(a.active); dec,bad_d=parse(a.decoy)
    excluded=set();
    if a.exclude and a.exclude.exists():
        with a.exclude.open(newline='') as f:
            for row in csv.DictReader(f):
                smi=row.get('smiles','')
                mol=Chem.MolFromSmiles(smi) if smi else None
                if mol: excluded.add(Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True))
    act={s:v for s,v in act.items() if s not in excluded}
    dec={s:v for s,v in dec.items() if s not in excluded and s not in act}
    if len(act)<a.n: raise ValueError('Not enough actives')
    rng=random.Random(a.seed)
    active_rows=rng.sample([act[k] for k in sorted(act)],a.n)
    decoy_rows=[dec[k] for k in sorted(dec)]
    A=np.stack([r['descs'] for r in active_rows]); D=np.stack([r['descs'] for r in decoy_rows]);
    cost=np.empty((len(A),len(D)),dtype=np.float64)
    feasible=np.empty_like(cost,dtype=bool)
    for i,ar in enumerate(A):
        dif=np.abs(D[:,:6]-ar[:6]); exact_charge=D[:,6]==ar[6]
        valid=exact_charge & np.all(dif<=CALIPER,axis=1)
        c=np.sum(((D[:,:6]-ar[:6])/SCALE)**2,axis=1)
        cost[i]=np.where(valid,c,1e8)
        feasible[i]=valid
    row,col=linear_sum_assignment(cost)
    if any(cost[row[k],col[k]]>=1e8 for k in range(len(row))):
        available=feasible.sum(1).tolist()
        raise RuntimeError(f'Strict preset matching calipers cannot match all {a.n} actives; feasible counts per active: {available}. Do not silently relax; reduce n or pre-register revised calipers.')
    matches=[]
    for idx, (i,j) in enumerate(zip(row,col),1):
        ac=active_rows[i]; dc=decoy_rows[j]
        matches.append({'pair_id':f'pair_{idx:03d}', 'active':ac,'decoy':dc,'cost':round(float(cost[i,j]),6)})
    # Matched pairs retained for paired and unpaired sensitivity analyses.
    with (a.out/'matched_panel.csv').open('w',newline='') as f:
        w=csv.DictWriter(f, fieldnames=['pair_id','class','smiles','source_id'] ); w.writeheader()
        for match in matches:
            for c in ['active','decoy']:
                x=match[c]; w.writerow({'pair_id':match['pair_id'],'class':c,'smiles':x['smiles'],'source_id':x['source_id']})
    with (a.out/'matched_descriptors.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['pair_id','class','smiles']+FIELDS); w.writeheader()
        for match in matches:
            for c in ['active','decoy']:
                x=match[c]; w.writerow({'pair_id':match['pair_id'],'class':c,'smiles':x['smiles'],**dict(zip(FIELDS,np.round(x['descs'],4)))})
    diff=np.stack([m['active']['descs']-m['decoy']['descs'] for m in matches])
    AA=np.stack([m['active']['descs'] for m in matches]); DD=np.stack([m['decoy']['descs'] for m in matches])
    stats={}
    for j,k in enumerate(FIELDS):
        pooled=np.sqrt((AA[:,j].var(ddof=1)+DD[:,j].var(ddof=1))/2)
        stats[k]={'active_mean':round(float(AA[:,j].mean()),4),'decoy_mean':round(float(DD[:,j].mean()),4),
                  'smd':round(float((AA[:,j].mean()-DD[:,j].mean())/pooled),4) if pooled else 0,
                  'max_pair_abs_difference':round(float(np.abs(diff[:,j]).max()),4)}
    audit={'seed':a.seed,'n_per_class':a.n, 'selection_uses_docking_scores':False,
           'source_sha256':{'active':hash_file(a.active),'decoy':hash_file(a.decoy)},
           'invalid_smiles':{'active':bad_a,'decoy':bad_d},
           'previously_docked_smiles_excluded_count':len(excluded),
           'matching_features':FIELDS,'calipers_MW_logP_HBA_HBD_TPSA_RotB':CALIPER.tolist(),
           'exact_formal_charge_matching':True,
           'matching_method':'Minimum-total-cost 1:1 without replacement via linear_sum_assignment',
           'descriptor_balance':stats,'panel_sha256':hash_file(a.out/'matched_panel.csv')}
    (a.out/'selection_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({'selected_pairs':a.n,'max_abs_smd':max(abs(v['smd']) for v in stats.values()),
                     'balance':{k:v['smd'] for k,v in stats.items()},
                     'panel_sha256':audit['panel_sha256']},indent=2))
if __name__=='__main__': main()
