"""Run Vina independently when the original sequence stopped at DOCK6 failure.
No labels are read. Inputs, Vina box/seed/budget and candidate membership stay fixed.
Do not retry Vina timeouts or change DOCK6's missing-score results.
"""
from pathlib import Path
import json,re,subprocess,concurrent.futures,shutil
from rdkit import Chem
from meeko import MoleculePreparation,PDBQTWriterLegacy
p=Path.home()/'aipharma_holdout_20261016';d=p/'docking';protocol=json.loads((d/'docking_protocol.json').read_text());folders=[]
for folder in sorted(d.glob('H*')):
 rp=folder/'result.json'
 if not rp.exists():continue
 result=json.loads(rp.read_text())
 if result.get('error')=='No DOCK6 scores' and 'vina_score' not in result and not result.get('independent_vina_attempted'):folders.append(folder)
def process(folder):
 rp=folder/'result.json';r=json.loads(rp.read_text());backup=folder/'before_independent_vina.json'
 if not backup.exists():shutil.copyfile(rp,backup)
 r['independent_vina_attempted']=True
 try:
  opts=r['MMFF_results'];good=[i for i,x in enumerate(opts) if x[0]==0];best=min(good,key=lambda i:opts[i][1]);mols=[m for m in Chem.SDMolSupplier(str(folder/'conformers.sdf'),removeHs=False) if m is not None];mol=next(m for m in mols if m.GetProp('_Name').endswith('_conf_'+str(best)))
  setup=MoleculePreparation(rigid_macrocycles=True).prepare(mol)[0];text,ok,error=PDBQTWriterLegacy.write_string(setup)
  if not ok:raise RuntimeError(str(error))
  (folder/'ligand.pdbqt').write_text(text);box=protocol['Vina']['pocket'];c=box['center'];s=box['size'];cmd=[str(Path.home()/'miniforge3/envs/aipharma/bin/vina'),'--receptor',protocol['Vina']['receptor'],'--ligand','ligand.pdbqt','--center_x',str(c[0]),'--center_y',str(c[1]),'--center_z',str(c[2]),'--size_x',str(s[0]),'--size_y',str(s[1]),'--size_z',str(s[2]),'--exhaustiveness','16','--seed','20261016','--cpu','1','--out','vina_pose.pdbqt']
  with (folder/'vina.log').open('w') as f:subprocess.run(cmd,cwd=folder,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=900)
  score=re.search(r'VINA RESULT:\s*(-?\d+\.\d+)',(folder/'vina_pose.pdbqt').read_text())
  if not score:raise RuntimeError('No Vina score')
  r['vina_score']=float(score.group(1));r['lowest_energy_input_conformer']=best;r['engine_status']={'DOCK6':'no_score','Vina':'completed'}
 except Exception as exc:r['independent_vina_error']=str(exc)
 rp.write_text(json.dumps(r,indent=2));print(json.dumps({'case_id':r['name'],'vina_score':r.get('vina_score'),'independent_vina_error':r.get('independent_vina_error'),'DOCK6':'failure preserved'}),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(process,folders))
print('Independent Vina attempts:',len(folders),flush=True)
