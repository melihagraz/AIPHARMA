"""Unblind only after scoring completes; evaluate fixed held-out structures.
Pair-bootstrap CIs preserve the matched design. Both primary engine ROC tests
are reported with Holm correction; consensus is exploratory.
"""
from pathlib import Path
import json,csv,hashlib,math,datetime
from collections import Counter
import numpy as np
from scipy.stats import rankdata,hypergeom,binom
from sklearn.metrics import roc_auc_score,average_precision_score,roc_curve
from rdkit import Chem,DataStructs
from rdkit.Chem import Descriptors,rdFingerprintGenerator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path.home()/'aipharma_holdout_20261016';d=p/'docking';locked=p/'locked_bounded';out=p/'evaluation';out.mkdir(exist_ok=True)
scorefile=d/'summary.json'
if not scorefile.exists():raise RuntimeError('Scoring summary absent; labels must remain withheld')
scoring=json.loads(scorefile.read_text());manifest=json.loads((locked/'manifest.json').read_text())
assert hashlib.sha256((locked/'structures_blinded.json').read_bytes()).hexdigest()==manifest['blinded_structures_sha256']==scoring['locked_selection_sha256']
assert hashlib.sha256((locked/'label_key.json').read_bytes()).hexdigest()==manifest['label_key_sha256']
assert scoring['completed']+scoring['failed']==manifest['actual_structures']
blinded=json.loads((locked/'structures_blinded.json').read_text())['structures']
for item in blinded:
 check=json.loads((d/item['case_id']/'result.json').read_text())
 if check.get('error')=='No DOCK6 scores' and not check.get('independent_vina_attempted'):raise RuntimeError('Complete independent Vina attempt before unblinding: '+item['case_id'])
labels=json.loads((locked/'label_key.json').read_text());inventory=json.loads((locked/'exclusion_inventory.json').read_text());parent_overlap={r['parent_identity'] for r in labels.values()}&set(inventory['normalized_parent_smiles']);scaffold_overlap={r['scaffold'] for r in labels.values()}&set(inventory['scaffolds']);assert not parent_overlap and not scaffold_overlap;assert len({r['scaffold'] for r in labels.values()})==len(labels);allids=sorted(labels);yall=np.array([labels[i]['label'] for i in allids],int)
assert len(set(yall))==2 and len(allids)==manifest['actual_structures']
(out/'unblinding_audit.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scoring_complete_before_label_key_read':True,'score_summary_sha256':hashlib.sha256(scorefile.read_bytes()).hexdigest(),'label_key_sha256':manifest['label_key_sha256'],'structures_sha256':manifest['blinded_structures_sha256'],'case_result_sha256':{i:hashlib.sha256((d/i/'result.json').read_bytes()).hexdigest() for i in allids}},indent=2))
results={}
for case in allids:
 resultfile=d/case/'result.json';results[case]=json.loads(resultfile.read_text())
 assert results[case]['smiles']==labels[case]['smiles']
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
known=[r for r in csv.DictReader((p/'alk5_active_control_benchmark.csv').open()) if r['label']=='1'];fps=[fpgen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in known]
base_mw={i:Descriptors.MolWt(Chem.MolFromSmiles(labels[i]['smiles'])) for i in allids};base_sim={i:max(DataStructs.BulkTanimotoSimilarity(fpgen.GetFingerprint(Chem.MolFromSmiles(labels[i]['smiles'])),fps)) for i in allids}
engine={name:{i:-float(results[i][field]) for i in allids if field in results[i] and math.isfinite(float(results[i][field]))} for name,field in [('Vina','vina_score'),('DOCK6','dock6_score')]}
common=sorted(set(engine['Vina'])&set(engine['DOCK6']));engine['Mean_rank']={}
if common:
 vr=rankdata([-engine['Vina'][i] for i in common]);dr=rankdata([-engine['DOCK6'][i] for i in common]);engine['Mean_rank']={i:-float((vr[j]+dr[j])/2) for j,i in enumerate(common)}
def auc_fast(y,values):
 n1=int(np.sum(y));n0=len(y)-n1
 return float((np.sum(rankdata(values)[y==1])-n1*(n1+1)/2)/(n1*n0))
def basic(ids,values):
 y=np.array([labels[i]['label'] for i in ids],int);values=np.asarray(values,float);n1=int(y.sum());n=len(y)
 if n1==0 or n1==n:raise RuntimeError('Only one class remains')
 auc=float(roc_auc_score(y,values));assert abs(auc-auc_fast(y,values))<1e-12
 order=sorted(range(n),key=lambda j:(-values[j],ids[j]));k=min(10,n);found=int(y[order[:k]].sum());threshold=values[order[k-1]];better=np.where(values>threshold)[0];tied=np.where(values==threshold)[0];slots=k-len(better);positive_better=int(y[better].sum());positive_tied=int(y[tied].sum());negative_tied=len(tied)-positive_tied;top_pair_counts=Counter(labels[ids[j]]['pair'] for j in order[:k]);both=sum(v==2 for v in top_pair_counts.values());single=sum(v==1 for v in top_pair_counts.values())
 return {'n':n,'active':n1,'weak':n-n1,'roc_auc':auc,'average_precision':float(average_precision_score(y,values)),'top_k':k,'top_k_active':found,'top_k_expected_random':k*n1/n,'top_k_active_tie_range':[positive_better+max(0,slots-negative_tied),positive_better+min(slots,positive_tied)],'top_k_case_ids':[ids[j] for j in order[:k]],'enrichment_top_k':(found/k)/(n1/n),'top_k_hypergeometric_p':float(hypergeom.sf(found-1,n,n1,k)),'top_k_matched_pair_null_p':float(binom.sf(found-both-1,single,0.5)),'top_k_pair_composition':{'both_selected':both,'one_selected':single}}
def evaluate(name,mapping):
 ids=sorted(mapping);values=np.array([mapping[i] for i in ids]);y=np.array([labels[i]['label'] for i in ids],int);result=basic(ids,values);rng=np.random.default_rng(20261016);pairids=sorted({labels[i]['pair'] for i in ids});by_pair={pair:[j for j,i in enumerate(ids) if labels[i]['pair']==pair] for pair in pairids};mw=np.array([base_mw[i] for i in ids]);sim=np.array([base_sim[i] for i in ids]);boot={'roc_auc':[],'average_precision':[],'auc_minus_MW':[],'auc_minus_similarity':[]}
 for _ in range(2000):
  sample=rng.choice(pairids,len(pairids),replace=True);ix=np.array([j for pair in sample for j in by_pair[pair]],int);yy=y[ix]
  if len(set(yy))<2:continue
  aa=auc_fast(yy,values[ix]);boot['roc_auc'].append(aa);boot['average_precision'].append(float(average_precision_score(yy,values[ix])));boot['auc_minus_MW'].append(aa-auc_fast(yy,mw[ix]));boot['auc_minus_similarity'].append(aa-auc_fast(yy,sim[ix]))
 result['pair_bootstrap_95pct_CI']={key:[float(x) for x in np.quantile(vals,[0.025,0.975])] for key,vals in boot.items()};result['baseline_same_cases']={'MW_roc_auc':float(roc_auc_score(y,mw)),'similarity_roc_auc':float(roc_auc_score(y,sim))}
 null_rng=np.random.default_rng(20261016);perms=5000;count=0
 allpairs=sorted({labels[i]['pair'] for i in allids});pair_position={pair:j for j,pair in enumerate(allpairs)}
 for _ in range(perms):
  flips=null_rng.integers(0,2,len(allpairs));permuted=np.array([labels[i]['label'] ^ int(flips[pair_position[labels[i]['pair']]]) for i in ids])
  if len(set(permuted))<2:continue
  if auc_fast(permuted,values)>=result['roc_auc']-1e-12:count+=1
 result['roc_permutation_p_one_sided']=(count+1)/(perms+1);result['permutations']=perms;result['permutation_unit']='label swaps within matched active/weak pairs';result['tie_break_rule']='fixed randomized case ID at cutoff; tie range also reported';result['excluded_case_ids']=sorted(set(allids)-set(ids));result['excluded_labels']=dict(Counter('active' if labels[i]['label'] else 'weak' for i in result['excluded_case_ids']))
 if result['excluded_case_ids']:
  missing_rank=min(values)-1;result['failure_at_bottom_sensitivity']=basic(allids,[mapping.get(i,missing_rank) for i in allids])
 return result
metrics={name:evaluate(name,mapping) for name,mapping in engine.items() if len(mapping)>=2 and len({labels[i]['label'] for i in mapping})==2}
ordered=sorted([n for n in ['Vina','DOCK6'] if n in metrics],key=lambda n:metrics[n]['roc_permutation_p_one_sided']);running=0
for j,name in enumerate(ordered):
 adjusted=max(running,min(1,(len(ordered)-j)*metrics[name]['roc_permutation_p_one_sided']));running=adjusted;metrics[name]['roc_permutation_p_Holm_two_engines']=adjusted
baselines={'MW_higher_is_better':basic(allids,[base_mw[i] for i in allids]),'prior_active_similarity':basic(allids,[base_sim[i] for i in allids])}
audit={'known_parent_overlap':len(parent_overlap),'known_scaffold_overlap':manifest['known_scaffold_overlap'],'holdout_scaffold_duplicates':manifest['within_holdout_scaffold_duplicates'],'labels_withheld_from_scoring':True,'class_counts':dict(Counter(int(x) for x in yall)),'endpoint':'IC50','active':'median exact-measurement pChEMBL >=7; threshold corresponds to the 100 nM IC50 scale','weak':'measured lower bound IC50 >=10,000 nM; exact value unknown','prior_training_overlap':'unknown','original_exact_only_viable':False,'bounded_measurement_protocol_amendment':True}
summary={'manifest':manifest,'audit':audit,'metrics':metrics,'baselines':baselines,'controls':scoring['controls'],'success_caveat':'Top10 >=8 is descriptive only; ROC uncertainty, baselines, missingness and assay-source bias must be considered. This tests retrospective docking recognition, not new-candidate activity or safety.'}
(out/'metrics.json').write_text(json.dumps(summary,indent=2))
fields=['case_id','label','label_evidence','median_exact_pchembl','weak_IC50_lower_bound_nM','vina_score','dock6_score','mean_rank','MW','prior_active_similarity','pair','scaffold','smiles','molecule_chembl_ids','document_chembl_ids','assay_chembl_ids','embedded_conformers','MMFF_converged','DOCK6_scored_conformers','workflow_status','workflow_error']
with (out/'unblinded_case_results.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(fields)
 for i in allids:
  row=labels[i];r=results[i];w.writerow([i,row['label'],row['label_evidence'],row['median_pchembl'],row['weak_ic50_lower_bound_nm'],r.get('vina_score'),r.get('dock6_score'),(-engine['Mean_rank'][i] if i in engine['Mean_rank'] else None),base_mw[i],base_sim[i],row['pair'],row['scaffold'],row['smiles'],';'.join(row['molecule_chembl_ids']),';'.join(row['document_chembl_ids']),';'.join(row['assay_chembl_ids']),r.get('embedded'),r.get('converged'),r.get('dock6_scored_conformers'),r.get('status'),r.get('error','')])
fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained');colors={'Vina':'#1768ac','DOCK6':'#e07826','Mean_rank':'#19836c'}
for name,mapping in engine.items():
 if name not in metrics:continue
 ids=sorted(mapping);y=[labels[i]['label'] for i in ids];vv=[mapping[i] for i in ids];fpr,tpr,_=roc_curve(y,vv);axes[0].plot(fpr,tpr,label=f"{name}: AUC {metrics[name]['roc_auc']:.2f} (n={metrics[name]['n']})",color=colors[name],linestyle='--' if name=='Mean_rank' else '-')
axes[0].plot([0,1],[0,1],':',color='gray',label='Random');axes[0].set(xlabel='False positive rate',ylabel='True positive rate',title='Label-withheld ALK5 holdout');axes[0].legend(loc='lower right')
names=list(metrics);found=[metrics[n]['top_k_active'] for n in names];axes[1].bar(names,found,color=[colors[n] for n in names]);axes[1].scatter(range(len(names)),[metrics[n]['top_k_expected_random'] for n in names],marker='D',color='black',s=30,label='Random expectation for scored subset');axes[1].set(ylim=(0,10.5),ylabel='Experimentally active compounds in top 10',title=f'{len(allids)//2} active + {len(allids)//2} weak; unseen recorded scaffolds');axes[1].legend(loc='upper right')
for i,value in enumerate(found):axes[1].text(i,value+0.15,str(value),ha='center')
fig.suptitle('Retrospective target-recognition test; weak IC50 values are right-censored',fontsize=11);fig.savefig(out/'holdout_results.png',dpi=170);fig.savefig(out/'holdout_results.svg');plt.close(fig)
document_counts={label:Counter(doc for row in labels.values() if row['label']==value for doc in row['document_chembl_ids']) for label,value in [('active',1),('weak',0)]}
(out/'assay_source_counts.json').write_text(json.dumps(document_counts,indent=2))
lines=['# ALK5 label-withheld holdout — 2026-10-16','',f"Final set: {manifest['actual_pairs']} experimentally active / weak pairs ({len(allids)} structures), one molecule per scaffold. Previously recorded parents and scaffolds are excluded.",'','## Protocol and data availability','','The original exact-IC50-only protocol found no new eligible weak controls after exclusions and was not scored. A separate, documented amendment admitted measured IC50 lower bounds of at least 10,000 nM as weak controls, while retaining the 100 nM active threshold, direct high-confidence biochemical target mapping, exclusions and matching. These censored observations are not assigned invented point IC50 values. All 18 selected weak controls use such lower-bound evidence. BindingDB was also audited; two additional exact weak structures survived structural exclusion, but were not included because comparable assay qualification was not established.','',f"Excluded history: {manifest['prior_parent_count']} normalized parents / {manifest['prior_scaffold_count']} exact Murcko scaffolds. No overlap remains. Pretrained-model training overlap and unrecorded generated structures remain unknown.",'','Scoring received only randomized case IDs and structures; labels remained in a separate hashed key. This is code-level label withholding, not a third-party administered blind trial. All structures retained representative experimental stereo/charge; no ADMET filter was applied to this binding-recognition benchmark.','', '## Results','','| Method | Scored | ROC-AUC (pair-bootstrap 95% CI) | Average precision | Actives in top 10 | Top10 chance p |','|:---|---:|:---|---:|---:|---:|']
for name,m in metrics.items():
 lo,hi=m['pair_bootstrap_95pct_CI']['roc_auc'];lines.append(f"| {name} | {m['n']}/{len(allids)} | {m['roc_auc']:.3f} ({lo:.3f}--{hi:.3f}) | {m['average_precision']:.3f} | {m['top_k_active']}/{m['top_k']} | {m['top_k_hypergeometric_p']:.4f} |")
lines.extend(['','Vina was also attempted independently on structures for which DOCK6 returned no score; failures in one engine did not prevent evaluation in the other. These attempts used the same conformers, box, seed and budget. DOCK6 failures and its original score threshold were preserved.','Consensus mean rank is secondary and uses the common scored subset; raw engine scores are never averaged. Primary engine ROC permutations swap labels within matched pairs and use Holm correction for two engines. Top10 also has a matched-pair null p-value in metrics.json alongside the descriptive hypergeometric reference. Full CIs, ties, failures and sensitivity results are in metrics.json.','', '![Holdout results](holdout_results.png)','', '## Simple baselines',''])
for name,b in baselines.items():lines.append(f"- {name}: ROC-AUC {b['roc_auc']:.3f}; top10 actives {b['top_k_active']}/10.")
lines.extend(['','Docking adds clear value over a baseline only if its paired-bootstrap AUC difference supports that interpretation; see auc_minus_MW and auc_minus_similarity intervals. A good docking AUC alone is insufficient to establish target specificity.','', '## Interpretation limits','','This is a small retrospective biochemical-IC50 recognition study on one ALK5 receptor conformation. It does not validate the generative prior, predicted ADMET safety, or biological activity of the new 20 candidates. Weak controls are censored and experimental conditions vary. The two pose controls informed protocol development, and exact scaffold exclusion does not guarantee low analogue similarity. Pair-bootstrap uncertainty does not account fully for study/assay-source correlations; assay_source_counts.json documents those sources.','', 'No parameters, membership, endpoint thresholds or rank weights were changed after scoring. The failure of the exact-only design and the bounded-value amendment remain archived.','', '## Files','','- metrics.json: complete evaluation and provenance.','- unblinded_case_results.csv: all observations and raw per-engine results; censored exact potencies remain blank.','- unblinding_audit.json: locked-file hashes and evaluation timestamp.','- ../locked_bounded/: scoring inputs, separate label key, exclusions and manifest.','- ../docking/: raw prepared molecules, scores and logs.','- [ChEMBL experimental data documentation](https://www.ebi.ac.uk/training/online/courses/chembl-quick-tour/accessing-chembl-data/viewing-results-from-the-web-interface/)'])
(out/'README.md').write_text('\n'.join(lines)+'\n');print(json.dumps({'metrics':{n:{k:m[k] for k in ['n','roc_auc','top_k_active','top_k_hypergeometric_p','pair_bootstrap_95pct_CI','baseline_same_cases']} for n,m in metrics.items()},'baselines':{n:b['roc_auc'] for n,b in baselines.items()}},indent=2),flush=True)
