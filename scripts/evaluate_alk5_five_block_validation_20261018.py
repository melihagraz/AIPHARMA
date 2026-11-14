"""Evaluate every locked profile and candidate, including failed pose controls."""
from pathlib import Path
import argparse,csv,hashlib,itertools,json,re
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors,Crippen
from sklearn.metrics import roc_auc_score
from spyrmsd import rmsd
from run_alk5_five_block_validation_20261018 import coordinates


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def auc(y,s):
    y,s=np.asarray(y),np.asarray(s);delta=s[y==1,None]-s[y==0]
    return float(np.mean((delta>0)+.5*(delta==0)))


def block_metrics(rows,key):
    y=np.array([r['label'] for r in rows]);s=np.array([r[key] for r in rows]);observed=auc(y,s)
    assert abs(observed-roc_auc_score(y,s))<1e-12
    blocks=sorted({r['block'] for r in rows});indices=[[i for i,r in enumerate(rows) if r['block']==b] for b in blocks]
    assert len(blocks)==5 and all(len(ix)==2 and {int(y[i]) for i in ix}=={0,1} for ix in indices)
    pair_results=[]
    for ix in indices:
        a,w=sorted(ix,key=lambda i:-y[i]);pair_results.append(float(s[a]>s[w])+.5*float(s[a]==s[w]))
    null=[]
    for switches in itertools.product([False,True],repeat=5):
        yp=y.copy()
        for ix,switch in zip(indices,switches):
            if switch:yp[ix]=1-yp[ix]
        null.append(auc(yp,s))
    rng=np.random.default_rng(20261018);boot=[]
    for _ in range(2000):
        ix=np.concatenate([indices[k] for k in rng.integers(0,5,5)])
        boot.append(auc(y[ix],s[ix]))
    return {'auc':observed,'paired_block_bootstrap_95':np.quantile(boot,[.025,.975]).tolist(),
            'within_block_correct_fraction':float(np.mean(pair_results)),'within_block_results':pair_results,
            'unadjusted_one_sided_block_swap_p':float(np.mean(np.array(null)>=observed-1e-12)),
            'block_swaps':len(null),'bootstrap_samples':2000}


def evaluate(out):
    run=json.loads((out/'summary.json').read_text());protocol=run['protocol'];items=json.loads((out/'structures.json').read_text())
    labels=json.loads((out/'label_key.json').read_text());labelmap={r['id']:r for r in labels};raw=run['results']
    assert len(raw)==252 and len({(r['id'],r['pdb'],r['function'],r['seed']) for r in raw})==252
    assert sha(out/'structures.json')==protocol['structures_sha256'] and sha(out/'label_key.json')==protocol['label_key_sha256']
    o=out/'evaluation';o.mkdir(exist_ok=True)
    provenance=[]
    for label in sorted(labels,key=lambda r:(r['block'],-r['label'])):
        source=label['provenance']
        provenance.append({'id':label['id'],'readable_label':f"block_{label['block']}_"+('active' if label['label'] else 'weak'),
            'block':label['block'],'label':label['label'],'origin':source['origin'],
            'monomer':source.get('monomer'),'median_pchembl':source.get('median_pchembl'),
            'reported_IC50_interval_nM':source.get('reported_IC50_interval_nM'),
            'reported_weak_values_nM':source.get('weak_values_nM'),'documents':source['documents'],
            'original_documents':source.get('original_documents'),'assay_ids':source.get('assay_ids'),
            'activity_ids':source.get('activity_ids'),'smiles':source.get('original_smiles') or source['smiles']})
    with (o/'controls_provenance.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(provenance[0]));writer.writeheader();writer.writerows(provenance)
    # Check raw pose remarks independently of the Python energy array serialization.
    remark_checks=0
    for r in raw:
        if r.get('score') is None:continue
        path=out/'jobs'/f"{r['pdb']}_{r['function']}_{r['id']}_{r['seed']}"/'poses.pdbqt'
        match=re.search(r'VINA RESULT:\s*([-\d.]+)',path.read_text());assert match
        assert abs(float(match.group(1))-r['score'])<.0011;remark_checks+=1
    aggregate=[];profiles={};pose_checks=0
    for pdb in protocol['targets']:
        for sf in protocol['functions']:
            for item in items:
                rr=[r for r in raw if r['id']==item['id'] and r['pdb']==pdb and r['function']==sf]
                assert len(rr)==3
                row={'id':item['id'],'kind':item['kind'],'pdb':pdb,'function':sf,'completed_seeds':sum(r['status']=='complete' for r in rr),
                     'seeds':[{k:r.get(k) for k in ['seed','score','status','error']} for r in sorted(rr,key=lambda r:r['seed'])]}
                mol=Chem.MolFromSmiles(item['smiles']);row.update({'mw':Descriptors.MolWt(mol),'logp':Crippen.MolLogP(mol),'heavy_atoms':mol.GetNumHeavyAtoms()})
                if item['id'] in labelmap:row.update({k:labelmap[item['id']][k] for k in ['label','block']})
                if row['completed_seeds']==3:
                    scores=[r['score'] for r in rr];row.update({'mean_score':float(np.mean(scores)),'seed_score_range':float(np.ptp(scores)),
                        'higher_is_ranked_first':-float(np.mean(scores)),'score_per_heavy_atom':-float(np.mean(scores))/row['heavy_atoms']})
                    poses={}
                    template=(out/'ligands'/item['id']/'ligand.pdbqt').read_text()
                    for r in rr:
                        text=(out/'jobs'/f"{pdb}_{sf}_{item['id']}_{r['seed']}"/'poses.pdbqt').read_text()
                        m,x=coordinates(template,text);poses[r['seed']]=(m,x)
                    differences=[]
                    for a,b in itertools.combinations(sorted(poses),2):
                        m,x=poses[a];_,y=poses[b];numbers=np.array([a.GetAtomicNum() for a in m.GetAtoms()]);adj=Chem.GetAdjacencyMatrix(m)
                        value=float(rmsd.symmrmsd(x,y,numbers,numbers,adj,adj,minimize=False));differences.append({'seed_a':a,'seed_b':b,'RMSD_A':value});pose_checks+=1
                    row['seed_pose_RMSDs']=differences;row['max_seed_pose_RMSD_A']=max(r['RMSD_A'] for r in differences)
                    row['pose_variability_flag_gt2A']=row['max_seed_pose_RMSD_A']>2
                aggregate.append(row)
            part=[r for r in aggregate if r['pdb']==pdb and r['function']==sf]
            scored=[r for r in part if r.get('mean_score') is not None]
            for row in scored:row['rank_among_13']=1+sum(r['mean_score']<row['mean_score'] for r in scored);row['rank_denominator']=len(scored)
            native=[r for r in raw if r['pdb']==pdb and r['function']==sf and r['kind']=='native_control'];assert len(native)==3
            profile={'native_pose_passes':sum(bool(r.get('native_pose_pass')) for r in native),'native_pose_controls':native,
                     'all_native_pose_controls_pass':all(r.get('native_pose_pass') for r in native),
                     'pose_variable_molecules':sum(bool(r.get('pose_variability_flag_gt2A')) for r in part),
                     'candidate_results':[r for r in part if r['kind']=='candidate']}
            controls=[r for r in part if r['kind']=='control']
            profile['observed_control_property_ranges']={key:[min(r[key] for r in controls),max(r[key] for r in controls)] for key in ['mw','logp','heavy_atoms']}
            for candidate in profile['candidate_results']:
                candidate['outside_observed_control_ranges']=[key for key,(lo,hi) in profile['observed_control_property_ranges'].items() if candidate[key]<lo or candidate[key]>hi]
            profile['coverage_limit']='Observed minimum/maximum only, not a validated applicability domain; raw candidate ranks can be affected by property differences from controls.'
            if len(controls)==10 and all(r.get('mean_score') is not None for r in controls):
                methods={}
                for key in ['higher_is_ranked_first','score_per_heavy_atom','mw','logp','heavy_atoms']:
                    methods[key]=block_metrics(controls,key)
                # Report both size directions rather than choosing a favorable sign afterwards.
                for row in controls:row['smaller_mw']=-row['mw'];row['smaller_heavy_atoms']=-row['heavy_atoms']
                for key in ['smaller_mw','smaller_heavy_atoms']:methods[key]=block_metrics(controls,key)
                seed_aucs={}
                for seed in protocol['seeds']:
                    yy=[r['label'] for r in controls];ss=[-next(x['score'] for x in r['seeds'] if x['seed']==seed) for r in controls]
                    seed_aucs[str(seed)]=auc(yy,ss)
                    assert abs(seed_aucs[str(seed)]-roc_auc_score(yy,ss))<1e-12
                profile.update({'status':'complete','methods':methods,'per_seed_auc':seed_aucs})
            else:profile.update({'status':'incomplete','reason':'No complete-case control deletion; at least one locked control missing a seed'})
            profiles[pdb+'_'+sf]=profile
    results={'protocol':protocol,'completed_jobs':sum(r['status']=='complete' for r in raw),'scheduled_jobs':252,
             'failed_jobs':[r for r in raw if r['status']!='complete'],'profiles':profiles,'per_structure':aggregate,
             'limits':['Five blocks, only five weak controls; retrospective and partially previously exposed, not external blind validation.',
                       'Source/scaffold graph separation is not proof of biological independence; assay/ATP/substrate differences persist.',
                       'Confidence intervals condition on selected blocks, fixed poses and scores; they do not cover all curation or design uncertainty.',
                       'Block-swap p values are exploratory, one-sided and unadjusted across receptors, engines and baselines.',
                       'Failed native recovery and seed variability qualify interpretation; no best-receptor or best-seed selection.',
                       'Docking ranks do not establish candidate activity, affinity, efficacy or safety; no MD or ADMET validation in this experiment.']}
    (o/'metrics.json').write_text(json.dumps(results,indent=2))
    with (o/'scores.csv').open('w',newline='') as f:
        fields=['pdb','function','id','kind','label','block','completed_seeds','mean_score','seed_score_range','max_seed_pose_RMSD_A','pose_variability_flag_gt2A','rank_among_13','rank_denominator','mw','logp','heavy_atoms']
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(aggregate)
    verification={'passed':True,'raw_pose_score_remark_checks':remark_checks,'study_seed_pose_pair_checks':pose_checks,
                  'locked_structures_hash_passed':True,'locked_label_key_hash_passed':True,'record_count':len(raw),
                  'profiles':len(profiles),'all_three_candidates_reported_in_all_profiles':all(len(r['candidate_results'])==3 for r in profiles.values()),
                  'native_control_count':sum(r['kind']=='native_control' for r in raw),'native_pass_count':sum(bool(r.get('native_pose_pass')) for r in raw if r['kind']=='native_control'),
                  'script_sha256':sha(Path(__file__))}
    independent=json.loads((out/'native_rmsd_independent_audit.json').read_text())
    assert independent['native_records']==18 and independent['agreements_within_0_01A']==18 and independent['same_pass_fail_decisions']==18
    verification['independent_native_RMSD_agreements']=18
    (o/'verification.json').write_text(json.dumps(verification,indent=2))
    print(json.dumps({'profiles':{key:{'auc':v.get('methods',{}).get('higher_is_ranked_first',{}).get('auc'),'native_passes':v['native_pose_passes'],'candidate_ranks':{r['id']:r.get('rank_among_13') for r in v['candidate_results']}} for key,v in profiles.items()},'verification':verification},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);evaluate(p.parse_args().output)
