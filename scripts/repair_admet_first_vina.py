"""Repair Vina 1.1.2-incompatible Meeko macrocycle pseudoatoms, without reselection.
Only parse failures citing CG0 are eligible; preserve initial files/results.
Uses the same lowest-energy converged conformer, box, Vina seed and search.
"""
from pathlib import Path
import json,subprocess,re,shutil,hashlib
from rdkit import Chem
from meeko import MoleculePreparation,PDBQTWriterLegacy
from scipy.stats import rankdata,spearmanr
root=Path.home();run=root/'aipharma_generation_20261016';p=run/'docking';source=p/'summary.json'
if not source.exists():raise RuntimeError('Wait for the original batch to finish before repairs')
protocol=json.loads((p/'docking_protocol.json').read_text());selected=json.loads((run/'results/locked_top20.json').read_text())['candidates'];lookup={x['selection_order']:x for x in selected}
repairs=[]
for folder in sorted(p.glob('candidate_*')):
 result_path=folder/'result.json'
 if not result_path.exists():continue
 result=json.loads(result_path.read_text())
 if result['status']=='completed':continue
 log=(folder/'vina.log').read_text() if (folder/'vina.log').exists() else ''
 if not ('CG0' in log and 'not a valid AutoDock type' in log):continue
 initial=dict(result);(folder/'initial_failure.json').write_text(json.dumps(initial,indent=2));shutil.copyfile(folder/'vina.log',folder/'vina_initial.log');shutil.copyfile(folder/'ligand.pdbqt',folder/'ligand_initial.pdbqt')
 opts=result['MMFF_results'];good=[i for i,x in enumerate(opts) if x[0]==0];best=min(good,key=lambda i:opts[i][1]);mols=[m for m in Chem.SDMolSupplier(str(folder/'conformers.sdf'),removeHs=False) if m is not None];mol=next(m for m in mols if m.GetProp('_Name').endswith('_conf_'+str(best)))
 setups=MoleculePreparation(rigid_macrocycles=True).prepare(mol);text,ok,error=PDBQTWriterLegacy.write_string(setups[0])
 if not ok:raise RuntimeError(str(error))
 assert 'CG0' not in text
 (folder/'ligand.pdbqt').write_text(text);box=protocol['Vina']['pocket'];c=box['center'];s=box['size'];cmd=[str(root/'miniforge3/envs/aipharma/bin/vina'),'--receptor',protocol['Vina']['receptor'],'--ligand','ligand.pdbqt','--center_x',str(c[0]),'--center_y',str(c[1]),'--center_z',str(c[2]),'--size_x',str(s[0]),'--size_y',str(s[1]),'--size_z',str(s[2]),'--exhaustiveness','16','--seed','20261016','--cpu','1','--out','vina_pose.pdbqt']
 with (folder/'vina.log').open('w') as f:proc=subprocess.run(cmd,cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=900)
 if proc.returncode==0:
  match=re.search(r'VINA RESULT:\s*(-?\d+\.\d+)',(folder/'vina_pose.pdbqt').read_text())
  if not match:raise RuntimeError('No repaired Vina score')
  result.update({'status':'completed','vina_score':float(match.group(1)),'lowest_energy_input_conformer':best,'preparation_repair':'Meeko rigid_macrocycles=True; Vina 1.1.2 lacks CG0 support; ring conformation fixed to existing converged input','initial_failure_preserved':True});result.pop('error',None)
 else:result['repair_failed']=True
 result_path.write_text(json.dumps(result,indent=2));repairs.append({'candidate':result['name'],'status':result['status'],'vina_score':result.get('vina_score')});print(json.dumps(repairs[-1]),flush=True)
summary=json.loads(source.read_text());rows=[json.loads((folder/'result.json').read_text()) for folder in sorted(p.glob('candidate_*')) if (folder/'result.json').exists()];valid=[r for r in rows if r['status']=='completed'];failed=[r for r in rows if r['status']!='completed'];vr=rankdata([r['vina_score'] for r in valid]);dr=rankdata([r['dock6_score'] for r in valid])
for i,row in enumerate(valid):row.update({'vina_rank':float(vr[i]),'dock6_rank':float(dr[i]),'mean_rank':float((vr[i]+dr[i])/2),'selection_data':lookup[row['selection_order']]})
valid.sort(key=lambda r:(r['mean_rank'],r['selection_order']));rho,pval=spearmanr([r['vina_score'] for r in valid],[r['dock6_score'] for r in valid]);summary.update({'completed':len(valid),'failed':len(failed),'ranked':valid,'failures':failed,'spearman_rho':float(rho),'spearman_p':float(pval),'format_repairs':repairs})
assert summary['locked_selection_sha256']==hashlib.sha256((run/'results/locked_top20.json').read_bytes()).hexdigest()
(p/'repair_amendment.json').write_text(json.dumps({'reason':'Meeko macrocycle pseudoatom unsupported by installed Vina 1.1.2','rule':'Only failed CG0 parses are reprepared with rigid macrocycles; no changes to membership, conformers, score thresholds or Vina search','repairs':repairs,'ring_limitation':'Vina evaluates the selected ring conformer without macrocycle flexibility'},indent=2));(p/'summary_repaired.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'completed':len(valid),'failed':len(failed),'spearman_rho':float(rho),'format_repairs':repairs}),flush=True)
