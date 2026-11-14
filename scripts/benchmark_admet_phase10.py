#!/usr/bin/env python3
"""AIPHARMA Phase 10: audit candidate third-party ADMET panels and frozen-model transfer.

NOT a proof of fully independent validation: ADMET-AI v2 was trained on a
potentially different TDC release; source identity overlap remains unresolved.

Usage:
 .venv/bin/python scripts/benchmark_admet_phase10.py --repo . \
   --candidate-repo reports/validation_phase1/phase10_candidate_source --mode all

Only reads inputs and writes new Phase-10 outputs; never changes previous analyses.
"""
from __future__ import annotations
import argparse, hashlib, json, math, subprocess, sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ENDPOINTS = ('hERG','DILI')
BASE_NOTE = ('Candidate external data; independence of ADMET-AI v2 training and model-selection data NOT established. '
             'ADMET-AI v2 reports training on PyTDC 1.1.15, while the Phase-8 reference audit used PyTDC 1.0.0. '
             'Exact/parent/connectivity identity screening against locally available TDC partitions cannot prove that '
             'the pretrained model never encountered these molecules. No compound safety or platform superiority claim.')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')

def standardize(s):
    if not isinstance(s,str) or not s.strip(): return None
    mol=Chem.MolFromSmiles(s)
    if mol is None: return None
    exact=Chem.MolToSmiles(mol,isomericSmiles=True)
    try:
        parent=rdMolStandardize.FragmentParent(mol)
        parent_s=Chem.MolToSmiles(parent,isomericSmiles=True)
    except Exception:
        parent_s=exact
        parent=mol
    try: connection=Chem.MolToInchiKey(parent).split('-')[0]
    except Exception: connection=None
    try: scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=parent)
    except Exception: scaffold=None
    return {'exact':exact,'parent':parent_s,'connectivity':connection,'scaffold':scaffold}

def read_sources(folder,endpoint):
    # Names in external repositories can evolve. Detect the expected columns, not an invented filename.
    picks=[]
    for f in sorted((folder/'data'/'ground_truth').rglob('*.csv')):
        if endpoint.lower() not in f.name.lower(): continue
        try: df=pd.read_csv(f,low_memory=False)
        except Exception: continue
        names={x.lower():x for x in df.columns}
        if not {'smiles','y_true'} <= set(names): continue
        d=df.rename(columns={names['smiles']:'smiles', names['y_true']:'y_true'}).copy()
        if 'endpoint' in names:
            e= d[names['endpoint']].fillna('').astype(str).str.lower()
            d=d[e.str.contains('herg' if endpoint=='hERG' else 'dili')]
        if len(d): picks.append((f,d))
    return picks

def load_reference(root,endpoint):
    p=root/'reports'/'validation_phase1'/'phase8_admet'/endpoint
    train=p/'train_val_curated.csv'; test=p/'test_curated.csv'
    if not train.is_file() or not test.is_file():
        raise RuntimeError(f'Missing Phase 8 curated data: {train} and/or {test}. Run Phase 8 audit first.')
    arr=[]; ids={'exact':set(),'parent':set(),'connectivity':set(),'scaffold':set()}
    for fn in (train,test):
        df=pd.read_csv(fn)
        if not {'Drug','Y'} <= set(df.columns): raise RuntimeError(f'Bad TDC schema at {fn}')
        for s in df.Drug:
            d=standardize(s)
            if d:
                for k,v in d.items():
                    if v: ids[k].add(v)
        arr.append({'path':str(fn),'sha256':digest(fn),'n':len(df)})
    return ids,arr,pd.read_csv(train)

def prepare(root,source,out,endpoint):
    picks=read_sources(source,endpoint)
    if not picks: return {'status':'NO_CANDIDATE_FILE','endpoint':endpoint,'glob':str(source/'data'/'ground_truth'/'*.csv')}
    if len(picks)>1:
        return {'status':'AMBIGUOUS_MULTIPLE_FILES','endpoint':endpoint,'files':[str(p) for p,_ in picks]}
    src,df=picks[0]
    if df['y_true'].isna().any(): raise RuntimeError(f'{src}: null labels, review source')
    labels=pd.to_numeric(df['y_true'],errors='coerce')
    if not labels.isin([0,1]).all(): raise RuntimeError(f'{src}: nonbinary labels, review source')
    df=df.copy();df['y_true']=labels.astype(int)
    fields=[standardize(x) for x in df.smiles]
    df['exact']=[z['exact'] if z else None for z in fields]
    df['parent']=[z['parent'] if z else None for z in fields]
    df['connectivity']=[z['connectivity'] if z else None for z in fields]
    df['scaffold']=[z['scaffold'] if z else None for z in fields]
    n_bad=int(df.exact.isna().sum());df=df[df.exact.notna()].copy()
    conflicts=df.groupby('exact').y_true.nunique()
    bad=set(conflicts[conflicts>1].index)
    n_conflict=int(df.exact.isin(bad).sum())
    df=df[~df.exact.isin(bad)].drop_duplicates('exact').copy()
    ref,refs,_=load_reference(root,endpoint)
    for k in ('exact','parent','connectivity','scaffold'):
        df['overlap_'+k] = df[k].isin(ref[k]) & df[k].notna() & (df[k]!='')
    # Cautious identity filter (not a scaffold filter): excludes exact, parent, and first InChIKey block.
    df['identity_overlap_tdc']=(df.overlap_exact | df.overlap_parent | df.overlap_connectivity)
    df['screened_for_scoring']=~df.identity_overlap_tdc
    output=out/endpoint;output.mkdir(parents=True,exist_ok=True)
    df.to_csv(output/'candidate_annotated.csv',index=False)
    score_df=df[df.screened_for_scoring].copy()
    score_df.to_csv(output/'identity_filtered_candidates.csv',index=False)
    sums={
        'status':'CANDIDATE_IDENTITY_AUDIT_COMPLETED_NOT_INDEPENDENCE_PROOF',
        'endpoint':endpoint,'source_file':str(src),'source_sha256':digest(src),
        'reference':refs,'rdkit_version':rdBase.rdkitVersion,
        'original_rows':len(pd.read_csv(src)),'invalid_structure_rows':n_bad,
        'conflicting_exact_identity_rows_removed':n_conflict,
        'deduplicated_n':int(len(df)),
        'deduplicated_positive':int(df.y_true.sum()),'deduplicated_negative':int((df.y_true==0).sum()),
        'overlap_counts':{k:int(df['overlap_'+k].sum()) for k in ('exact','parent','connectivity','scaffold')},
        'overlap_identity_union_n':int(df.identity_overlap_tdc.sum()),
        'filtered_n':int(len(score_df)),
        'filtered_positive':int(score_df.y_true.sum()),
        'filtered_negative':int((score_df.y_true==0).sum()),
        'filtered_scaffold_overlap_n':int(score_df.overlap_scaffold.sum()),
        'candidate_sha256':digest(output/'candidate_annotated.csv'),
        'filtered_sha256':digest(output/'identity_filtered_candidates.csv'),
        'warning':BASE_NOTE,
    }
    write(output/'audit.json',sums)
    return sums

def descriptors(smi):
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    m=Chem.MolFromSmiles(smi)
    return [Descriptors.MolWt(m),Descriptors.MolLogP(m),rdMolDescriptors.CalcTPSA(m),
            Lipinski.NumHDonors(m),Lipinski.NumHAcceptors(m),Lipinski.NumRotatableBonds(m),
            m.GetNumHeavyAtoms(),rdMolDescriptors.CalcFractionCSP3(m)]

def metrics(y,p):
    y=np.asarray(y,dtype=int);p=np.asarray(p,dtype=float)
    return {'n':len(y),'positive':int(y.sum()),'negative':int((y==0).sum()),
            'roc_auc':float(roc_auc_score(y,p)),
            'average_precision':float(average_precision_score(y,p)),
            'brier':float(brier_score_loss(y,p)) if np.min(p)>=0 and np.max(p)<=1 else None}

def paired_ci(y,a,b,B=2000,seed=20261112):
    rng=np.random.default_rng(seed);y=np.asarray(y,dtype=int);a=np.asarray(a);b=np.asarray(b)
    pos=np.where(y==1)[0];neg=np.where(y==0)[0];ds=[]
    for i in range(B):
        ids=np.concatenate([rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)])
        ds.append(float(roc_auc_score(y[ids],a[ids])-roc_auc_score(y[ids],b[ids])))
    return [float(x) for x in np.percentile(ds,[2.5,97.5])]

def evaluate(root,out,endpoint,batch):
    folder=out/endpoint;aud=folder/'audit.json';dat=folder/'identity_filtered_candidates.csv'
    if not aud.exists() or not dat.exists():
        return {'status':'MISSING_AUDIT_OR_CANDIDATES','endpoint':endpoint}
    a=json.loads(aud.read_text());df=pd.read_csv(dat)
    if digest(dat)!=a['filtered_sha256']:raise RuntimeError(f'{endpoint} filtered data changed')
    npos=int(df.y_true.sum());nneg=int((df.y_true==0).sum())
    if min(npos,nneg)<10:
        return {'status':'INSUFFICIENT_CLASSES_FOR_PREDEFINED_EVALUATION','endpoint':endpoint,
                'n':len(df),'positive':npos,'negative':nneg,'minimum_required_per_class':10,
                'warning':BASE_NOTE}
    sys.path.insert(0,str(root))
    from Admet_App.app.predictor import AdmetPredictor
    predictor=AdmetPredictor(num_workers=0)
    ref_train=pd.read_csv(root/'reports'/'validation_phase1'/'phase8_admet'/endpoint/'train_val_curated.csv')
    model=make_pipeline(StandardScaler(),LogisticRegression(max_iter=2000,C=1.0,class_weight='balanced',random_state=20261112))
    model.fit(np.asarray([descriptors(x) for x in ref_train.Drug]),ref_train.Y.to_numpy(dtype=int))
    baseline=model.predict_proba(np.asarray([descriptors(x) for x in df.smiles]))[:,1]
    preds=[]
    for start in range(0,len(df),batch):
        smiles=df.smiles.iloc[start:start+batch].tolist()
        records=predictor.predict(smiles)
        if len(records)!=len(smiles):raise RuntimeError(f'Bad ADMET-AI count for {endpoint}')
        for record in records:
            score=record.raw.get(endpoint)
            if score is None or not math.isfinite(float(score)):
                raise RuntimeError(f'Missing ADMET endpoint {endpoint}: {list(record.raw)}')
            preds.append(float(score))
    y=df.y_true.to_numpy(dtype=int)
    result=df[['smiles','y_true','exact','parent','connectivity','overlap_scaffold']].copy()
    result['admet_ai_score']=preds
    result['rdkit_baseline_score']=baseline
    result.to_csv(folder/'external_candidate_predictions.csv',index=False)
    ans={'status':'EXPLORATORY_CANDIDATE_TRANSFER_NOT_PROVEN_INDEPENDENT','endpoint':endpoint,
         'models':{'ADMET_AI':metrics(y,preds),'RDKit_8_feature_logistic':metrics(y,baseline)},
         'auc_difference':float(roc_auc_score(y,preds)-roc_auc_score(y,baseline)),
         'paired_stratified_bootstrap95_auc_difference':paired_ci(y,preds,baseline),
         'prediction_sha256':digest(folder/'external_candidate_predictions.csv'),
         'all_filtered_scaffold_overlap_n':int(df.overlap_scaffold.sum()),'warning':BASE_NOTE}
    write(folder/'evaluation.json',ans)
    return ans

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--candidate-repo',type=Path,required=True)
    p.add_argument('--mode',choices=['audit','evaluate','all'],default='audit')
    p.add_argument('--batch',type=int,default=32)
    opt=p.parse_args()
    root=opt.repo.resolve();source=opt.candidate_repo.resolve()
    if not source.is_dir():p.error(f'Candidate source directory missing: {source}')
    if opt.batch<1:p.error('--batch must be positive')
    out=root/'reports'/'validation_phase1'/'phase10_external_admet';out.mkdir(parents=True,exist_ok=True)

    try:
        commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True,stderr=subprocess.DEVNULL).strip()
    except Exception:
        commit=None
    all_results={'status':'CANDIDATE_SOURCE_UNVERIFIED','source_repo':str(source),'source_git_commit':commit,
                 'generated_utc':datetime.now(timezone.utc).isoformat(),'warning':BASE_NOTE,
                 'endpoints':{}}
    for ep in ENDPOINTS:
        if opt.mode in ('audit','all'):
            ans=prepare(root,source,out,ep)
            all_results['endpoints'][ep]={'audit':ans}
            print(f'{ep} AUDIT: n={ans.get("filtered_n")} positive={ans.get("filtered_positive")} '
                  f'negative={ans.get("filtered_negative")} identity_overlap={ans.get("overlap_identity_union_n")} '
                  f'scaffold_overlap={ans.get("filtered_scaffold_overlap_n")} status={ans["status"]}',flush=True)
        if opt.mode in ('evaluate','all'):
            ans=evaluate(root,out,ep,opt.batch)
            all_results['endpoints'].setdefault(ep,{})['evaluation']=ans
            print(f'{ep} EVALUATION: {ans["status"]} '
                  f'ADMET-AUC={ans.get("models",{}).get("ADMET_AI",{}).get("roc_auc")} '
                  f'RDKit-AUC={ans.get("models",{}).get("RDKit_8_feature_logistic",{}).get("roc_auc")}',flush=True)
    write(out/'phase10_summary.json',all_results)
    print('SAVED:',out/'phase10_summary.json')
    print('WARNING:',BASE_NOTE)
if __name__=='__main__':main()
