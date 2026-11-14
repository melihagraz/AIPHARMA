from pathlib import Path
import json,csv,hashlib,math,requests,shutil,os,re,subprocess,concurrent.futures,time
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem,Descriptors,Crippen,rdFingerprintGenerator
from meeko import MoleculePreparation,PDBQTWriterLegacy
p=Path.home()/'aipharma_ensemble_20261018';p.mkdir(exist_ok=True);env=Path.home()/'miniforge3/envs/aipharma/bin';os.environ['PATH']=str(env)+os.pathsep+os.environ['PATH']
if (p/'summary.json').exists():raise RuntimeError('Completed experiment is immutable')
locked=json.loads((Path.home()/'aipharma_generation_20261016/results/locked_top20.json').read_text())['candidates'];items=[]
for order in [19,18,2]:
 r=next(x for x in locked if x['selection_order']==order);m=Chem.MolFromSmiles(r['smiles']);assert not any(str(x.specified)=='Unspecified' for x in Chem.FindPotentialStereo(m));items.append({'id':f'candidate_{order:02d}','kind':'candidate','smiles':r['smiles']})
data=json.loads((Path.home()/'aipharma_holdout_20261016/endpoint_audit/datasets_all.json').read_text())['enzyme_catalytic__no_explicit_mutant'];used=set();fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
for label in [1,0]:
 for j,c in enumerate(items[:3],1):
  cm=Chem.MolFromSmiles(c['smiles']);pool=[]
  for r in data:
   m=Chem.MolFromSmiles(r['smiles'])
   if r['label']!=label or r['evidence']!='exact' or r['smiles'] in used or any(str(x.specified)=='Unspecified' for x in Chem.FindPotentialStereo(m)):continue
   if label==1 and r['median_pchembl']<8:continue
   dist=abs(Descriptors.MolWt(m)-Descriptors.MolWt(cm))/50+abs(Crippen.MolLogP(m)-Crippen.MolLogP(cm))
   pool.append((dist,r['smiles'],r))
  assert pool
  _,_,r=min(pool);used.add(r['smiles']);items.append({'id':('active_' if label else 'weak_')+str(j),'kind':'active_control' if label else 'weak_control','smiles':r['smiles'],'median_pchembl':r['median_pchembl'],'activity_ids':r['activity_ids'],'assay_ids':r['assay_ids'],'documents':r['documents'],'matched_candidate':c['id']})
protocol={'selection':'fixed candidate19,18,02 from historical review list; all stereo-resolved; controls property-nearest exact enzyme >=8 or <=5, chosen before new scores','targets':['1PY5','1VJY','2WOT'],'vina_seeds':[20261018,20261019],'exhaustiveness':8,'cpu_per_job':1,'ligands':'ETKDGv3 16 attempts seed20261018, MMFF94s1000 converged lowest-energy; Meeko rigid macrocycles','protein':'same OpenBabel pH7.4 rigid preparation for all three; local native-ligand pocket center, box24A each axis','controls':'native co-crystal redocking at seed20261018 gate symmetry-corrected heavy RMSD<2A. Failed controls remain reported; no tune','scope':'Exploratory receptor/seed sensitivity, three previously studied candidates. Not independent validation; absolute scores across receptors not combined as affinities','md_preselection':'candidate_02 plus property-matched active_3 and weak_3, 1PY5 seed20261018 poses; no switch based on scores'}
(p/'protocol.json').write_text(json.dumps(protocol,indent=2));(p/'locked_items.json').write_text(json.dumps(items,indent=2));print('Locked nine structures',flush=True)
for pdb,code in [('1PY5','PY1'),('1VJY','460'),('2WOT','ZZG')]:
 w=p/pdb;w.mkdir(exist_ok=True);raw=w/f'{pdb}.pdb'
 if not raw.exists():
  old=Path.home()/f'.aipharma_redock/{pdb}/{pdb}.pdb'
  if old.exists():shutil.copyfile(old,raw)
  else:r=requests.get(f'https://files.rcsb.org/download/{pdb}.pdb',timeout=90);r.raise_for_status();raw.write_text(r.text)
 lines=raw.read_text().splitlines();native=[x for x in lines if x.startswith('HETATM') and x[17:20].strip()==code];assert native
 first=(native[0][21:22],native[0][22:27]);native=[x for x in native if (x[21:22],x[22:27])==first];xyz=np.array([[float(x[30:38]),float(x[38:46]),float(x[46:54])] for x in native if x[76:78].strip()!='H']);center=(xyz.min(0)+xyz.max(0))/2
 (w/'native.pdb').write_text('\n'.join(native)+'\nEND\n');(w/'box.json').write_text(json.dumps({'center':center.tolist(),'size':[24,24,24],'ligand_code':code,'ligand_instance':first,'pdb_sha256':hashlib.sha256(raw.read_bytes()).hexdigest()},indent=2))
 clean=w/'protein.pdb';clean.write_text('\n'.join(x for x in lines if x.startswith(('ATOM','TER')))+'\nEND\n');done=subprocess.run([str(env/'obabel'),str(clean),'-O',str(w/'receptor.pdbqt'),'-xr','-p','7.4'],capture_output=True,text=True,check=True);(w/'receptor_prep.log').write_text(done.stderr)
 ccdpath=w/f'{code}_ccd.json'
 if not ccdpath.exists():r=requests.get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{code}',timeout=90);r.raise_for_status();ccdpath.write_text(r.text)
 ccd=json.loads(ccdpath.read_text());smi=next(r['descriptor'] for r in ccd['pdbx_chem_comp_descriptor'] if r['type']=='SMILES_CANONICAL');(w/'native_smiles.txt').write_text(smi)
# Prepare shared deterministic ligand files and native controls without starting from crystal coordinates.
allitems=items+[{'id':'native_'+pdb,'smiles':(p/pdb/'native_smiles.txt').read_text(),'kind':'native_control','pdb':pdb} for pdb in protocol['targets']]
for r in allitems:
 w=p/'ligands'/r['id'];w.mkdir(parents=True,exist_ok=True);m=Chem.AddHs(Chem.MolFromSmiles(r['smiles']));ep=AllChem.ETKDGv3();ep.randomSeed=20261018;ep.numThreads=1;ep.pruneRmsThresh=.25;ids=list(AllChem.EmbedMultipleConfs(m,16,ep));opt=AllChem.MMFFOptimizeMoleculeConfs(m,numThreads=1,maxIters=1000,mmffVariant='MMFF94s');good=[(o[1],i) for i,o in zip(ids,opt) if o[0]==0]
 if not good:(w/'error.json').write_text(json.dumps({'error':'no converged conformer'}));continue
 cid=min(good)[1];mm=Chem.Mol(m);conf=Chem.Conformer(m.GetConformer(cid));mm.RemoveAllConformers();mm.AddConformer(conf);writer=Chem.SDWriter(str(w/'prepared.sdf'));writer.write(mm);writer.close();setup=MoleculePreparation(rigid_macrocycles=True).prepare(mm)[0];txt,ok,err=PDBQTWriterLegacy.write_string(setup);assert ok,err;(w/'ligand.pdbqt').write_text(txt)
def dock(job):
 pdb,item,seed=job;w=p/'jobs'/f'{pdb}_{item["id"]}_{seed}';w.mkdir(parents=True,exist_ok=True);res={'pdb':pdb,'id':item['id'],'kind':item['kind'],'seed':seed}
 try:
  lig=p/'ligands'/item['id']/'ligand.pdbqt';box=json.loads((p/pdb/'box.json').read_text());cmd=[str(env/'vina'),'--receptor',str(p/pdb/'receptor.pdbqt'),'--ligand',str(lig),'--cpu','1','--seed',str(seed),'--exhaustiveness','8','--out',str(w/'pose.pdbqt')]
  for axis,c,sz in zip('xyz',box['center'],box['size']):cmd+=['--center_'+axis,str(c),'--size_'+axis,str(sz)]
  done=subprocess.run(cmd,capture_output=True,text=True,timeout=600);(w/'vina.log').write_text(done.stdout+'\n'+done.stderr);done.check_returncode();text=(w/'pose.pdbqt').read_text();res['score']=float(re.search(r'VINA RESULT:\s*(-?\d+\.\d+)',text).group(1));res['status']='complete'
  if item['kind']=='native_control':
   from spyrmsd import rmsd
   # Map crystal PDB atom names to CCD atoms; use CCD connectivity, not distance-inferred aromatic bonds.
   ccd=json.loads(next((p/pdb).glob('*_ccd.json')).read_text());atoms=ccd['chem_comp_atom'];name_to_idx={a['atom_id']:j for j,a in enumerate(atoms)};base=Chem.MolFromSmiles(item['smiles']);native=Chem.MolFromPDBFile(str(p/pdb/'native.pdb'),removeHs=True);gen=Chem.RemoveHs(Chem.SDMolSupplier(str(p/'ligands'/item['id']/'prepared.sdf'),removeHs=False)[0]);pose=[]
   for line in text.splitlines():
    if line.startswith(('ATOM','HETATM')) and line.split()[-1] not in ['H','HD','HS']:pose.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
   # PDB native connectivity may differ; compare symmetry against CCD native-order coordinates via OpenBabel SD conversion.
   conv=subprocess.run([str(env/'obabel'),str(p/pdb/'native.pdb'),'-O',str(w/'native.sdf')],capture_output=True,text=True,check=True);ref=Chem.SDMolSupplier(str(w/'native.sdf'),removeHs=True)[0]
   res['native_rmsd_A']=float(rmsd.symmrmsd(np.array(ref.GetConformer().GetPositions()),np.array(pose),np.array([a.GetAtomicNum() for a in ref.GetAtoms()]),np.array([a.GetAtomicNum() for a in gen.GetAtoms()]),Chem.GetAdjacencyMatrix(ref),Chem.GetAdjacencyMatrix(gen),minimize=False));res['pose_control_pass']=res['native_rmsd_A']<2
 except Exception as ex:res['error']=str(ex);res['status']='failed'
 (w/'result.json').write_text(json.dumps(res,indent=2));print(json.dumps(res),flush=True);return res
controls=[r for r in allitems if r['kind']=='native_control'];results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(dock,[(r['pdb'],r,20261018) for r in controls]))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results+=list(pool.map(dock,[(pdb,r,seed) for pdb in protocol['targets'] for r in items for seed in protocol['vina_seeds']]))
(p/'summary.json').write_text(json.dumps({'protocol':protocol,'results':results},indent=2));print('Completed ensemble study',flush=True)
