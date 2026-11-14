"""Dock the locked ADMET-first set without changing membership after scores.
Includes two DOCK6 pose controls; only converged MMFF94s conformers retained.
"""
from pathlib import Path
import json,re,subprocess,hashlib,concurrent.futures,time,csv,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import rankdata,spearmanr
from spyrmsd import rmsd as srmsd
from meeko import MoleculePreparation,PDBQTWriterLegacy
root=Path.home();runroot=root/'aipharma_generation_20261016';selected_path=runroot/'results/locked_top20.json';lockhash=hashlib.sha256(selected_path.read_bytes()).hexdigest()
assert lockhash==(runroot/'results/locked_top20.sha256').read_text().strip()
selected=json.loads(selected_path.read_text())['candidates'];out=runroot/'docking';out.mkdir(exist_ok=True);D=root/'tools/dock6';env=root/'miniforge3/envs/aipharma/bin';source=root/'.aipharma_dock6_amber_1PY5'
pocket=json.loads((runroot/'pah_run_1PY5.json').read_text())['pocket']
protocol={'locked_selection_sha256':lockhash,'candidate_count':len(selected),'ETKDGv3':{'attempts':32,'pruneRmsThresh_A':0.25,'seed':20261015},'MMFF94s':{'max_iterations':1000,'retain_only_status_zero':True},'DOCK6':{'search':'rigid ensemble orientation + simplex minimization','orientations':1000,'grid_spacing_A':0.2,'receptor_charge':'OpenMM ff14SB','ligand_charge':'Gasteiger','score':'grid'},'Vina':{'exhaustiveness':16,'seed':20261016,'cpu_per_job':1,'ligand_prep':'Meeko from lowest-energy converged conformer','receptor':str(root/'.aipharma_redock/1PY5/receptor.pdbqt'),'pocket':pocket},'comparison':'per-engine tied ranks; unweighted mean rank; no averaging raw scores','pose_controls':['1PY5','1VJY'],'control_gate':'both DOCK6 top-score symmetry-corrected heavy RMSDs <2 A','limits':['Internal ALK5 pose controls, not experimental validation','No new candidates inserted after selection','No exhaustive protonation/tautomer enumeration']}
(out/'docking_protocol.json').write_text(json.dumps(protocol,indent=2))
base=dict(line.split(None,1) for line in (source/'ensemble.in').read_text().splitlines() if line.strip())

def process(item):
 start=time.time();name=item['name'];p=out/name;p.mkdir(exist_ok=True);pdb=item.get('pdb','1PY5');gridroot=root/('.aipharma_dock6_amber_'+pdb)
 result={'name':name,'type':item['type'],'smiles':item['smiles'],'selection_order':item.get('selection_order'),'docking_protocol_sha256':hashlib.sha256((out/'docking_protocol.json').read_bytes()).hexdigest()}
 try:
  mol=Chem.AddHs(Chem.MolFromSmiles(item['smiles']));settings=AllChem.ETKDGv3();settings.randomSeed=20261015;settings.pruneRmsThresh=0.25;settings.numThreads=1
  ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=32,params=settings));opts=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000,mmffVariant='MMFF94s');good=[cid for cid,opt in zip(ids,opts) if opt[0]==0]
  result.update({'embedded':len(ids),'converged':len(good),'MMFF_results':[list(x) for x in opts]})
  if not good:raise RuntimeError('No converged conformers; no fallback included')
  writer=Chem.SDWriter(str(p/'conformers.sdf'))
  for cid in good:mol.SetProp('_Name',f'{name}_conf_{cid}');writer.write(mol,confId=cid)
  writer.close()
  conv=subprocess.run([str(env/'obabel'),'conformers.sdf','-O','conformers_raw.mol2','--partialcharge','gasteiger'],cwd=p,capture_output=True,text=True,check=True);(p/'prep.log').write_text(conv.stderr)
  normal=[]
  for b in (p/'conformers_raw.mol2').read_text().split('@<TRIPOS>MOLECULE')[1:]:
   lines=b.strip('\n').splitlines();counts=lines[1].split();counts[2]='1';lines[1]=' '.join(counts);section=None
   for i,line in enumerate(lines):
    if line.startswith('@<TRIPOS>'):section=line
    elif section=='@<TRIPOS>ATOM' and line.strip():a=line.split();a[6]='1';lines[i]=' '.join(a)
   normal.append('@<TRIPOS>MOLECULE\n'+'\n'.join(lines).replace('@<TRIPOS>ATOM','\n@<TRIPOS>ATOM',1)+'\n@<TRIPOS>SUBSTRUCTURE\n1 UNL1 1 RESIDUE 1 A UNL1 0 ROOT\n')
  (p/'conformers.mol2').write_text(''.join(normal));params=base.copy();params.update({'ligand_atom_file':'conformers.mol2','receptor_site_file':str(gridroot/'site.sph'),'grid_score_grid_prefix':str(gridroot/'grid'),'ligand_outfile_prefix':'dock6','num_final_scored_poses':1})
  (p/'dock6.in').write_text(''.join(f'{k:60s} {v}\n' for k,v in params.items()))
  with (p/'dock6.log').open('w') as f:subprocess.run([str(D/'bin/dock6'),'-i','dock6.in','-o','dock6.out'],cwd=p,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=3600)
  entries=[]
  for e in re.split(r'(?=##########\s+Name:)',(p/'dock6_scored.mol2').read_text()):
   match=re.search(r'Grid_Score:\s*([-\d.]+)',e)
   if match:entries.append((float(match.group(1)),e))
  if not entries:raise RuntimeError('No DOCK6 scores')
  score,top=min(entries,key=lambda x:x[0]);result.update({'dock6_score':score,'dock6_scored_conformers':len(entries),'best_dock6_conformer':re.search(r'Name:\s*(\S+)',top).group(1)})
  if item['type']=='control':
   ref=Chem.SDMolSupplier(str(gridroot/'crystal.sdf'),removeHs=True)[0];refnums=np.array([a.GetAtomicNum() for a in ref.GetAtoms()]);refadj=Chem.GetAdjacencyMatrix(ref);refxyz=np.array(ref.GetConformer().GetPositions())
   gen=Chem.RemoveHs(mol);nums=np.array([a.GetAtomicNum() for a in gen.GetAtoms()]);adj=Chem.GetAdjacencyMatrix(gen)
   atoms=[a.split() for a in top.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].splitlines() if a.strip() and not a.split()[5].startswith('H')];xyz=np.array([a[2:5] for a in atoms],float)
   result['top_pose_rmsd_A']=float(srmsd.symmrmsd(refxyz,xyz,refnums,nums,refadj,adj,minimize=False));result['pose_recovery_pass']=result['top_pose_rmsd_A']<2
  else:
   bestid=min(good,key=lambda cid:opts[ids.index(cid)][1]);lig=Chem.Mol(mol);conf=Chem.Conformer(mol.GetConformer(bestid));lig.RemoveAllConformers();lig.AddConformer(conf,assignId=True)
   setups=MoleculePreparation().prepare(lig);pdbqt,ok,error=PDBQTWriterLegacy.write_string(setups[0]);
   if not ok:raise RuntimeError('Meeko preparation: '+str(error))
   (p/'ligand.pdbqt').write_text(pdbqt);c=pocket['center'];s=pocket['size'];cmd=[str(env/'vina'),'--receptor',str(root/'.aipharma_redock/1PY5/receptor.pdbqt'),'--ligand','ligand.pdbqt','--center_x',str(c[0]),'--center_y',str(c[1]),'--center_z',str(c[2]),'--size_x',str(s[0]),'--size_y',str(s[1]),'--size_z',str(s[2]),'--exhaustiveness','16','--seed','20261016','--cpu','1','--out','vina_pose.pdbqt']
   with (p/'vina.log').open('w') as f:subprocess.run(cmd,cwd=p,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=900)
   match=re.search(r'VINA RESULT:\s*(-?\d+\.\d+)',(p/'vina_pose.pdbqt').read_text())
   if not match:raise RuntimeError('No Vina result')
   result['vina_score']=float(match.group(1));result['lowest_energy_input_conformer']=bestid
  result['status']='completed'
 except Exception as e:result['status']='failed';result['error']=str(e)
 result['elapsed_seconds']=round(time.time()-start,2);(p/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['MMFF_results','smiles','docking_protocol_sha256']}),flush=True);return result

controls=[]
for pdb,code in [('1PY5','PY1'),('1VJY','460')]:
 desc=json.loads((root/('.aipharma_dock6_amber_'+pdb)/(code+'_ccd.json')).read_text())['pdbx_chem_comp_descriptor'];smi=next(d['descriptor'] for d in desc if d['type']=='SMILES_CANONICAL');controls.append({'type':'control','name':'control_'+pdb,'pdb':pdb,'smiles':smi})
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(process,i) for i in controls]
 for f in concurrent.futures.as_completed(futures):results.append(f.result())
 if not all(r.get('pose_recovery_pass',False) for r in results):
  (out/'control_failure.json').write_text(json.dumps(results,indent=2));raise RuntimeError('Pose control gate failed; candidate runs withheld')
 jobs=[{**r,'type':'candidate','name':f"candidate_{r['selection_order']:02d}"} for r in selected]
 futures=[pool.submit(process,i) for i in jobs]
 for f in concurrent.futures.as_completed(futures):results.append(f.result());(out/'checkpoint.json').write_text(json.dumps(results,indent=2))
assert lockhash==hashlib.sha256(selected_path.read_bytes()).hexdigest()
valid=[r for r in results if r['type']=='candidate' and r['status']=='completed'];failures=[r for r in results if r['type']=='candidate' and r['status']!='completed']
lookup={r['selection_order']:r for r in selected}
if valid:
 vr=rankdata([r['vina_score'] for r in valid]);dr=rankdata([r['dock6_score'] for r in valid])
 for i,r in enumerate(valid):r.update({'vina_rank':float(vr[i]),'dock6_rank':float(dr[i]),'mean_rank':float((vr[i]+dr[i])/2),'selection_data':lookup[r['selection_order']]})
 valid.sort(key=lambda r:(r['mean_rank'],r['selection_order']));rho,pval=spearmanr([r['vina_score'] for r in valid],[r['dock6_score'] for r in valid])
else:rho=pval=None
summary={'locked_selection_sha256':lockhash,'selected':len(selected),'completed':len(valid),'failed':len(failures),'controls':[r for r in results if r['type']=='control'],'spearman_rho':None if rho is None else float(rho),'spearman_p':None if pval is None else float(pval),'ranked':valid,'failures':failures,'limitations':protocol['limits']}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'completed':len(valid),'failed':len(failures),'spearman_rho':summary['spearman_rho']}),flush=True)
