"""Join fixed-set DOCK6 results to existing ADMET predictions and write audit tables."""
from pathlib import Path
import json,csv,re,hashlib
from rdkit import Chem
from scipy.stats import spearmanr,rankdata
root=Path.home();p=root/'.aipharma_dock6_candidates'
s=json.loads((p/'summary.json').read_text());ad=json.loads((p/'pah_admet_1PY5.json').read_text())
canon=lambda x:Chem.MolToSmiles(Chem.MolFromSmiles(x))
lookup={canon(r['smiles']):r for r in ad['rows']}
for r in s['ranked']:
 a=lookup[canon(r['smiles'])];r['admet']={k:a[k] for k in ['ames','herg','dili']};r['passes_all_three_at_0_5']=all(a[k]<0.5 for k in ['ames','herg','dili'])
 text=(p/f"candidate_{r['original_rank']:02d}"/'docked_scored.mol2').read_text();pose_rows=[]
 for e in re.split(r'(?=##########\s+Name:)',text):
  m=re.search(r'Name:\s*(\S+)',e);sc=re.search(r'Grid_Score:\s*([-\d.]+)',e)
  if not (m and sc):continue
  cid=int(m.group(1).rsplit('_',1)[1]);opt=r['mmff_results'][cid]
  pose_rows.append({'name':m.group(1),'score':float(sc.group(1)),'mmff_status':opt[0],'mmff_energy_kcal_mol':opt[1]})
 best=min(pose_rows,key=lambda x:x['score']);r['best_pose']=best
 converged=[x for x in pose_rows if x['mmff_status']==0]
 r['best_converged_conformer_dock6_score']=min(x['score'] for x in converged) if converged else None
 r['converged_conformers_scored']=len(converged)
 r['unassigned_stereocenter_count']=sum(x[1]=='?' for x in r['unassigned_stereocenters'])
s['admet_pass_all_three_count']=sum(r['passes_all_three_at_0_5'] for r in s['ranked'])
strict=[r for r in s['ranked'] if r['best_converged_conformer_dock6_score'] is not None]
if len(strict)>2:
 rho,pval=spearmanr([r['vina_score'] for r in strict],[r['best_converged_conformer_dock6_score'] for r in strict]);s['converged_only_sensitivity']={'n':len(strict),'spearman_rho':float(rho),'spearman_p':float(pval),'caveat':'Post hoc sensitivity analysis; does not replace the all-conformer primary analysis'}
s['pose_controls']={pdb:json.loads((root/('.aipharma_dock6_amber_'+pdb)/'ensemble_summary.json').read_text())['top_scored_pose'] for pdb in ['1PY5','1VJY']}
s['best_pose_mmff_nonconverged_candidates']=sum(r['best_pose']['mmff_status']!=0 for r in s['ranked'])
s['admet_source_sha256']=hashlib.sha256((p/'pah_admet_1PY5.json').read_bytes()).hexdigest()
(p/'integrated_summary.json').write_text(json.dumps(s,indent=2))
columns=['original_vina_order','vina_score','vina_rank_with_ties','dock6_score','dock6_rank','mean_rank','best_pose_mmff_status','embedded_conformers','scored_conformers','converged_conformers_scored','converged_only_dock6_score','ames','herg','dili','passes_all_three_at_0_5','unassigned_stereocenters','smiles']
with (p/'rank_comparison.csv').open('w',newline='') as f:
 writer=csv.writer(f);writer.writerow(columns)
 for r in s['ranked']:writer.writerow([r['original_rank'],r['vina_score'],r['vina_rank_with_ties'],r['dock6_score'],r['dock6_rank'],r['mean_rank'],r['best_pose']['mmff_status'],r['embedded_conformers'],r['scored_conformers'],r['converged_conformers_scored'],r['best_converged_conformer_dock6_score'],r['admet']['ames'],r['admet']['herg'],r['admet']['dili'],r['passes_all_three_at_0_5'],r['unassigned_stereocenter_count'],r['smiles']])
lines=['# Fixed 12-candidate orthogonal DOCK6 follow-up','',f"Completed: {s['completed']}/12; failed: {s['failed']}.",f"Vina/DOCK6 Spearman rho: {s['spearman_rho']:.4f}; descriptive p = {s['spearman_p']:.4f}.",'','Primary ranking uses all conformers under the fixed protocol. Mean rank is an exploratory summary, not a calibrated combined score. This is a Vina-preselected, range-restricted set; no enrichment or biological validation is inferred.','', '| Original Vina order | Vina score | DOCK6 score | Mean rank | Top-pose MMFF converged | Ames / hERG / DILI |','|---:|---:|---:|---:|:---:|:---|']
for r in s['ranked']:
 a=r['admet'];lines.append(f"| {r['original_rank']} | {r['vina_score']:.2f} | {r['dock6_score']:.2f} | {r['mean_rank']:.2f} | {'yes' if r['best_pose']['mmff_status']==0 else 'no'} | {a['ames']:.3f} / {a['herg']:.3f} / {a['dili']:.3f} |")
lines.extend(['',f"Candidates passing all three existing ADMET probability thresholds (<0.5): {s['admet_pass_all_three_count']}/12. No experimental compound nomination is made.",f"Top docking poses from nonconverged MMFF conformers: {s['best_pose_mmff_nonconverged_candidates']}/12.",'','MMFF status 0 means converged. Nonconverged conformers were retained consistently with the control protocol and flagged, not silently removed. Converged-only sensitivity is reported separately in integrated_summary.json. Protonation, tautomer and stereoisomer coverage remains incomplete. Scores from the two engines have different scales; compare ranks only.'])
if 'converged_only_sensitivity' in s:lines.append('Converged-only sensitivity: '+json.dumps(s['converged_only_sensitivity']))
(p/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:s[k] for k in ['completed','failed','spearman_rho','spearman_p','admet_pass_all_three_count','best_pose_mmff_nonconverged_candidates','converged_only_sensitivity']}))
