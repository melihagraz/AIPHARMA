"""Exploratory orthogonal docking of the fixed 12 ALK5 Vina candidates.
Uses the same 32-attempt ETKDG/MMFF94s/rigid DOCK6 protocol as pose controls.
No new candidates are added. Scores are never averaged across engines.
"""
from pathlib import Path
import json,re,subprocess,hashlib,time,concurrent.futures,shutil,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import spearmanr,rankdata
root=Path.home();D=root/'tools/dock6';env=root/'miniforge3/envs/aipharma/bin'
source=root/'.aipharma_dock6_amber_1PY5';out=root/'.aipharma_dock6_candidates';out.mkdir(exist_ok=True)
input_file=out/'pah_run_1PY5.json';data=json.loads(input_file.read_text());candidates=data['ranked'];assert len(candidates)==12
manifest={'source_sha256':hashlib.sha256(input_file.read_bytes()).hexdigest(),'candidate_count':12,'target':'ALK5','pdb':'1PY5','seed':20261015,'conformer_attempts':32,'prune_rms_A':0.25,'grid_spacing_A':0.2,'orientation_limit':1000,'ranking':'mean rank of Vina and DOCK6; per-engine scores not combined','candidates':candidates}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
base=dict(line.split(None,1) for line in (source/'ensemble.in').read_text().splitlines() if line.strip())
base.update({'grid_score_grid_prefix':str(source/'grid'),'receptor_site_file':str(source/'site.sph')})

def process(candidate):
 started=time.time();rank=candidate['rank'];p=out/f'candidate_{rank:02d}';p.mkdir(exist_ok=True)
 result={'original_rank':rank,'smiles':candidate['smiles'],'vina_score':candidate['dock_score']}
 try:
  mol=Chem.MolFromSmiles(candidate['smiles'])
  if mol is None:raise RuntimeError('SMILES parse failure')
  result['unassigned_stereocenters']=Chem.FindMolChiralCenters(mol,includeUnassigned=True)
  mol=Chem.AddHs(mol);settings=AllChem.ETKDGv3();settings.randomSeed=20261015;settings.pruneRmsThresh=0.25;settings.numThreads=1
  ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=32,params=settings))
  if not ids:raise RuntimeError('No conformers embedded')
  opt=AllChem.MMFFOptimizeMoleculeConfs(mol,mmffVariant='MMFF94s',numThreads=1)
  writer=Chem.SDWriter(str(p/'ensemble.sdf'))
  for cid in ids:mol.SetProp('_Name',f'candidate_{rank:02d}_conf_{cid}');writer.write(mol,confId=cid)
  writer.close()
  conv=subprocess.run([str(env/'obabel'),'ensemble.sdf','-O','ensemble_raw.mol2','--partialcharge','gasteiger'],cwd=p,check=True,capture_output=True,text=True);(p/'prepare.log').write_text(conv.stderr)
  blocks=(p/'ensemble_raw.mol2').read_text().split('@<TRIPOS>MOLECULE')[1:];normal=[]
  for b in blocks:
   lines=b.strip('\n').splitlines();counts=lines[1].split();counts[2]='1';lines[1]=' '.join(counts);section=None
   for i,line in enumerate(lines):
    if line.startswith('@<TRIPOS>'):section=line
    elif section=='@<TRIPOS>ATOM' and line.strip():
     a=line.split();a[6]='1';lines[i]=' '.join(a)
   normal.append('@<TRIPOS>MOLECULE\n'+'\n'.join(lines).replace('@<TRIPOS>ATOM','\n@<TRIPOS>ATOM',1)+'\n@<TRIPOS>SUBSTRUCTURE\n1 UNL1 1 RESIDUE 1 A UNL1 0 ROOT\n')
  (p/'ensemble.mol2').write_text(''.join(normal));params=base.copy();params.update({'ligand_atom_file':'ensemble.mol2','ligand_outfile_prefix':'docked'})
  (p/'dock.in').write_text(''.join(f'{k:60s} {v}\n' for k,v in params.items()))
  with (p/'dock.log').open('w') as log:subprocess.run([str(D/'bin/dock6'),'-i','dock.in','-o','dock.out'],cwd=p,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=3600)
  text=(p/'docked_scored.mol2').read_text();scores=[float(x) for x in re.findall(r'Grid_Score:\s*([-\d.]+)',text)]
  if not scores:raise RuntimeError('DOCK6 returned no scored conformers')
  result.update({'status':'completed','dock6_score':min(scores),'embedded_conformers':len(ids),'scored_conformers':len(scores),'mmff_nonconverged':sum(x[0]!=0 for x in opt),'mmff_results':[list(x) for x in opt]})
 except Exception as exc:result.update({'status':'failed','error':str(exc)})
 result['elapsed_seconds']=round(time.time()-started,2);(p/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['mmff_results','smiles']}),flush=True);return result

results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 futures=[pool.submit(process,c) for c in candidates]
 for f in concurrent.futures.as_completed(futures):
  results.append(f.result());(out/'checkpoint.json').write_text(json.dumps(results,indent=2))
valid=[r for r in results if r['status']=='completed']
if valid:
 vina=np.array([r['vina_score'] for r in valid]);dock=np.array([r['dock6_score'] for r in valid]);vr=rankdata(vina);dr=rankdata(dock)
 for i,r in enumerate(valid):r.update({'vina_rank_with_ties':float(vr[i]),'dock6_rank':float(dr[i]),'mean_rank':float((vr[i]+dr[i])/2)})
 valid.sort(key=lambda r:r['mean_rank']);rho,pvalue=spearmanr(vina,dock)
else:rho=pvalue=None
summary={'target':'ALK5','pdb':'1PY5','scope':'exploratory orthogonal docking of fixed existing 12 candidates','manifest_sha256':hashlib.sha256((out/'manifest.json').read_bytes()).hexdigest(),'completed':len(valid),'failed':len(results)-len(valid),'spearman_rho':None if rho is None else float(rho),'spearman_p':None if pvalue is None else float(pvalue),'ranked':valid,'failures':[r for r in results if r['status']!='completed'],'limitations':['Protocol developed on two ALK5 co-crystals; not independent external validation','No biological affinity inference from scores','Conformer counts differ across compounds','Unspecified stereochemistry is not exhaustively enumerated','No new 20-member prospective set was constructed','ADMET results are not changed by this docking run']}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'completed':len(valid),'failed':summary['failed'],'spearman_rho':summary['spearman_rho']}),flush=True)
