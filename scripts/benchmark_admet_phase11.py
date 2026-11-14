#!/usr/bin/env python3
"""AIPHARMA Phase 11: post hoc ADMET transfer quality/uncertainty audit.

Reads Phase 8 + Phase 10 v2 frozen predictions; DOES NOT refit models,
change labels, select score orientation, download data, or run inference.

Usage:
  .venv/bin/python scripts/benchmark_admet_phase11.py --repo . --bootstrap 2000

LIMITS: Original ADMET-AI training molecules are not fully known here;
reference similarity to locally installed PyTDC 1.0 partitions is a proxy,
NOT proof of chemical novelty to the pretrained v2 model.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

ENDPOINTS = ('DILI', 'hERG')
LABEL_NOTES = {
    'DILI': 'Candidate source: DILIrank/DILIst vMost+vLess=1, vNo=0; potentially distinct from the PyTDC DILI label definition.',
    'hERG': 'Candidate source: CardioTox blocker at IC50 10 uM threshold (1=blocker); potentially distinct from the PyTDC hERG label definition.',
}

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(1024*1024),b''):
            h.update(part)
    return h.hexdigest()

def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')

def good_metric(y, score):
    y=np.asarray(y, dtype=int)
    score=np.asarray(score, dtype=float)
    if len(y)==0 or not np.isfinite(score).all() or len(set(y.tolist()))<2:
        return None
    ans={'n':int(len(y)), 'positive':int(y.sum()),'negative':int((y==0).sum()),
         'positive_prevalence_no_skill_AP':float(y.mean()),
         'roc_auc':float(roc_auc_score(y,score)),
         'average_precision':float(average_precision_score(y,score))}
    if np.min(score)>=0 and np.max(score)<=1:
        ans['brier']=float(brier_score_loss(y,score))
    return ans

def ci(y, a, b, repeats, seed):
    rng=np.random.default_rng(seed)
    y=np.asarray(y,dtype=int);a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
    if min(len(pos),len(neg))<5:
        return None
    vals=np.empty((repeats,3),dtype=float)
    for i in range(repeats):
        idx=np.concatenate((rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)))
        x=roc_auc_score(y[idx],a[idx]); z=roc_auc_score(y[idx],b[idx]);vals[i]=(x,z,x-z)
    q=np.quantile(vals,[.025,.975],axis=0)
    return {'ADMET_AI_AUC_95CI':list(map(float,q[:,0])),
            'RDKit_AUC_95CI':list(map(float,q[:,1])),
            'delta_ADMET_minus_RDKit_AUC_95CI':list(map(float,q[:,2])),
            'method':'paired class-stratified compound bootstrap; ignores document/scaffold clusters'}

def fp(smi, generator):
    mol=Chem.MolFromSmiles(str(smi))
    return generator.GetFingerprint(mol) if mol is not None else None

def run(root, nboot, seed):
    source=root/'reports'/'validation_phase1'/'phase10_external_admet_v2'
    ref_root=root/'reports'/'validation_phase1'/'phase8_admet'
    out=root/'reports'/'validation_phase1'/'phase11_admet_robustness'
    if not source.is_dir():
        raise SystemExit(f'Missing Phase 10 v2 folder: {source}')
    generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    final={'status':'POST_HOC_TRANSFER_ROBUSTNESS_NOT_PROVEN_EXTERNAL',
           'date_utc':datetime.now(timezone.utc).isoformat(),
           'bootstrap_replicates':nboot,'bootstrap_seed':seed,
           'rdkit':rdBase.rdkitVersion,
           'warnings':[
                'ADMET-AI v2 original training/model-selection sources not completely reconstructed; candidate source independence unresolved.',
                'Candidates were filtered against local PyTDC 1.0 reference, which may differ from ADMET-AI v2 original training data.',
                'Candidate endpoint labels may not match PyTDC task definitions or assay conditions exactly.',
                'These are post hoc analyses of pre-existing third-party model outputs, NOT superiority of AIPHARMA prediction algorithm.',
                'Compound-level bootstrap does not account for publication, series or scaffold dependence.',
                'Scaffold and similarity strata are exploratory; no significance testing planned for subgroups.'
           ], 'endpoints':{}}
    md=['# AIPHARMA Phase 11 — ADMET transfer robustness', '',
        '**Exploratory. Not independent external validation. No platform superiority claim.**','']
    for j, ep in enumerate(ENDPOINTS):
        folder=source/ep
        a_path=folder/'audit.json'; e_path=folder/'evaluation.json'
        d_path=folder/'identity_filtered_candidates.csv'
        p_path=folder/'external_candidate_predictions.csv'
        for p in (a_path,e_path,d_path,p_path):
            if not p.is_file():raise SystemExit(f'Missing: {p}')
        aud=json.loads(a_path.read_text()); ev=json.loads(e_path.read_text())
        if sha256(d_path)!=aud['filtered_sha256']:
            raise SystemExit(f'{ep}: filtered candidate SHA mismatch!')
        if sha256(p_path)!=ev['prediction_sha256']:
            raise SystemExit(f'{ep}: predicted candidate SHA mismatch!')
        dat=pd.read_csv(d_path);pred=pd.read_csv(p_path)
        if len(dat)!=len(pred) or len(dat)!=aud['filtered_n']:
            raise SystemExit(f'{ep}: candidate/prediction length mismatch!')
        for field in ('smiles','y_true','connectivity'):
            if field not in dat or field not in pred or not dat[field].fillna('').astype(str).reset_index(drop=True).equals(pred[field].fillna('').astype(str).reset_index(drop=True)):
                raise SystemExit(f'{ep}: ordered {field} mismatch between frozen inputs and predictions!')
        y=pred.y_true.to_numpy(dtype=int)
        a=pred.admet_ai_score.to_numpy(dtype=float)
        b=pred.rdkit_baseline_score.to_numpy(dtype=float)
        if not set(y).issubset({0,1}) or not (np.isfinite(a).all() and np.isfinite(b).all()):
            raise SystemExit(f'{ep}: bad labels or nonfinite scores')
        met=good_metric(y,a); base=good_metric(y,b)
        if met is None or base is None:raise SystemExit(f'{ep}: not both classes')
        if abs(met['roc_auc']-ev['models']['ADMET_AI']['roc_auc'])>1e-8 or abs(base['roc_auc']-ev['models']['RDKit_8_feature_logistic']['roc_auc'])>1e-8:
            raise SystemExit(f'{ep}: AUC does not match frozen Phase 10 report!')
        bounds=ci(y,a,b,nboot,seed+j)
        section={'source':aud.get('source_file'), 'candidate_label_definition':LABEL_NOTES[ep],
                 'frozen_pred_sha256_verified':True,
                 'n':len(pred), 'ADMET_AI':met, 'RDKit_8_feature_logistic':base,
                 'delta_auc_admet_minus_rdkit':float(met['roc_auc']-base['roc_auc']),
                 'bootstrap':bounds,
                 'scaffold_overlap_to_local_TDC':{},'ecfp4_similarity_to_local_TDC':{}}
        # Scaffold proxy strata. This is never full original training set independence.
        if 'overlap_scaffold' not in pred:raise SystemExit('Missing overlap_scaffold')
        flags=pred.overlap_scaffold.astype(str).str.lower().isin(('true','1','yes')).to_numpy()
        for name, sel in [('overlap',flags),('no_overlap',~flags)]:
            yy=y[sel]
            section['scaffold_overlap_to_local_TDC'][name]={
                'n':int(sel.sum()),'positive':int(yy.sum()),'negative':int((yy==0).sum()),
                'ADMET_AI':good_metric(yy,a[sel]) if min(int(yy.sum()),int((yy==0).sum()))>=5 else None,
                'RDKit':good_metric(yy,b[sel]) if min(int(yy.sum()),int((yy==0).sum()))>=5 else None,
                'note':'Descriptive only; AUC suppressed for fewer than 5 members of either class'}
        refs=[]
        for fn in ('train_val_curated.csv','test_curated.csv'):
            p=ref_root/ep/fn
            if not p.is_file():raise SystemExit(f'Missing {p}')
            frame=pd.read_csv(p)
            if 'Drug' not in frame:raise SystemExit(f'Missing Drug column in {p}')
            for smi in frame.Drug:
                f=fp(smi,generator)
                if f is not None:refs.append(f)
        if not refs:raise SystemExit('No reference fingerprints')
        max_tanimoto=[]
        for smi in pred.smiles:
            f=fp(smi,generator)
            if f is None:raise SystemExit(f'Bad candidate SMILES: {smi}')
            max_tanimoto.append(max(DataStructs.BulkTanimotoSimilarity(f,refs)))
        similarity=np.asarray(max_tanimoto,dtype=float)
        section['ecfp4_similarity_to_local_TDC']['summary']={
            'median':float(np.median(similarity)),'min':float(similarity.min()),'max':float(similarity.max()),
            'fraction_ge_0_8':float((similarity>=.8).mean()),
            'fraction_lt_0_4':float((similarity<.4).mean()),
            'note':'Reference is locally downloaded TDC 1.0 train+test; not verified original ADMET-AI v2 training set'}
        for label,sel in [('lt_0.4',similarity<.4),('0.4_to_0.8',(similarity>=.4)&(similarity<.8)),('ge_0.8',similarity>=.8)]:
            yy=y[sel]
            section['ecfp4_similarity_to_local_TDC'][label]={
                'n':int(sel.sum()),'positive':int(yy.sum()),'negative':int((yy==0).sum()),
                'ADMET_AI':good_metric(yy,a[sel]) if min(int(yy.sum()),int((yy==0).sum()))>=5 else None,
                'RDKit':good_metric(yy,b[sel]) if min(int(yy.sum()),int((yy==0).sum()))>=5 else None}
        output=pred[['smiles','y_true','admet_ai_score','rdkit_baseline_score','overlap_scaffold']].copy()
        output['max_ECFP4_Tanimoto_local_TDC_1_0']=similarity
        out.mkdir(parents=True,exist_ok=True)
        output.to_csv(out/f'{ep}_scored_with_similarity.csv',index=False)
        final['endpoints'][ep]=section
        cistr=bounds['delta_ADMET_minus_RDKit_AUC_95CI'] if bounds else [None,None]
        print(f'{ep}: n={len(pred)} ADMET={met["roc_auc"]:.4f} RDKit={base["roc_auc"]:.4f} '
              f'delta={met["roc_auc"]-base["roc_auc"]:+.4f} 95%CI={cistr} '
              f'local-TDC similarity<0.4 n={int((similarity<.4).sum())}')
        md.extend([f'## {ep}', '',
                   f'- n={len(pred)}; positives={int(y.sum())}, negatives={int((y==0).sum())}.',
                   f'- ADMET-AI AUC={met["roc_auc"]:.4f}; RDKit AUC={base["roc_auc"]:.4f}.',
                   f'- Delta AUC={met["roc_auc"]-base["roc_auc"]:+.4f}; compound-bootstrap 95% CI={cistr}.',
                   f'- Median maximum ECFP4 similarity to *locally available* TDC 1.0 reference={np.median(similarity):.3f}.',
                   f'- Tanimoto<0.4: {(similarity<.4).sum()} compounds; Tanimoto>=0.8: {(similarity>=.8).sum()} compounds.',
                   f'- Label: {LABEL_NOTES[ep]}',''])
    write_json(out/'phase11_summary.json',final)
    md.extend(['## Interpretation limits', '', *['- '+w for w in final['warnings']]])
    (out/'phase11_summary.md').write_text('\n'.join(md)+'\n', encoding='utf-8')
    print('SAVED:',out/'phase11_summary.json')
    print('SAVED:',out/'phase11_summary.md')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--bootstrap',type=int,default=2000)
    p.add_argument('--seed',type=int,default=20261112)
    args=p.parse_args()
    if args.bootstrap<100: p.error('bootstrap must be >=100')
    run(args.repo.resolve(),args.bootstrap,args.seed)
