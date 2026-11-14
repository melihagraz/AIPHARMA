"""Small integrated CLI worked example. Not a web-service or biological validation."""
from pathlib import Path
import argparse,os,sys,json,time,hashlib,shutil,subprocess,platform,importlib.util,importlib.metadata

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pdb-file',type=Path,required=True);ap.add_argument('--reinvent-dir',type=Path,default=Path(__file__).resolve().parents[1]/'Reinvent_App');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--generate',type=int,default=32);ap.add_argument('--dock-top',type=int,default=2);ap.add_argument('--seed',type=int,default=7);ap.add_argument('--repeat',type=int,default=2);args=ap.parse_args()
 if args.output.exists():raise RuntimeError('Use a fresh output directory; never reuse prior artifacts')
 args.output.mkdir(parents=True);os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+os.environ['PATH'];os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
 for exe in ['fpocket','obabel','vina']:
  if not shutil.which(exe):raise RuntimeError('Missing executable '+exe)
 import torch;torch.set_num_threads(1)
 from rdkit import Chem
 from rdkit.Chem import AllChem
 from meeko import MoleculePreparation,PDBQTWriterLegacy
 spec=importlib.util.spec_from_file_location('pah',Path(__file__).with_name('run_pah_pipeline.py'));core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);core.REINVENT_DIR=args.reinvent_dir.resolve()
 sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
 protocol={'input_pdb_sha256':sha(args.pdb_file),'prior_sha256':sha(args.reinvent_dir/'prior/reinvent.prior'),'seed':args.seed,'requested':args.generate,'dock_top':args.dock_top,'repeats':args.repeat,'sampling':'unconditional base prior','selection':'QED descending then canonical SMILES; no ADMET gating in this smoke test','ligand_prep':'ETKDGv3 seed, MMFF94s maxIters1000 converged only; Meeko rigid_macrocycles','vina':{'seed':args.seed,'exhaustiveness':2,'cpu':1,'timeout_seconds':300},'reuse':'PDB input/checkpoint reused; pocket detection, receptor prep, generation, filtering, ADMET and docking recomputed per repeat','scope':'CLI integration example using existing pipeline core functions; not web UI/service validation or scientific affinity validation'}
 (args.output/'protocol.json').write_text(json.dumps(protocol,indent=2));versions={k:importlib.metadata.version(k) for k in ['torch','rdkit','admet-ai','chemprop','meeko','numpy']};hardware={'platform':platform.platform(),'cpu':platform.processor(),'logical_cpu_count':os.cpu_count(),'torch_threads':torch.get_num_threads(),'cuda_available':torch.cuda.is_available()}
 (args.output/'environment.json').write_text(json.dumps({'versions':versions,'hardware':hardware,'python':sys.version,'binaries':{x:shutil.which(x) for x in ['fpocket','obabel','vina']},'runner_sha256':sha(__file__),'core_sha256':sha(Path(__file__).with_name('run_pah_pipeline.py'))},indent=2));runs=[]
 for repeat in range(1,args.repeat+1):
  w=args.output/f'run_{repeat}';w.mkdir();timings={};t0=time.perf_counter()
  def stage(name,fn):
   t=time.perf_counter();v=fn();timings[name]=time.perf_counter()-t;print(f'run{repeat} {name}: {timings[name]:.3f}s',flush=True);return v
  pdb=w/'1PY5.pdb';shutil.copyfile(args.pdb_file,pdb)
  pocket=stage('pocket_detection',lambda:core.detect_pocket(pdb,w));rec=stage('receptor_preparation',lambda:core.prep_receptor(pdb,w,pocket));raw=stage('generation',lambda:core.generate(args.generate,args.seed))
  (w/'raw_samples.json').write_text(json.dumps(raw,indent=2))
  def filt():
   candidates={}
   for smi,nll in raw:
    mol=Chem.MolFromSmiles(smi)
    if mol is None:continue
    smi=Chem.MolToSmiles(mol);d=core.drug_like(smi)
    if d and d['pass']:candidates.setdefault(smi,{'smiles':smi,'nll':float(nll),'properties':d})
   return sorted(candidates.values(),key=lambda r:(-r['properties']['qed'],r['smiles']))
  cands=stage('physchem_filter',filt);top=cands[:args.dock_top]
  if len(top)!=args.dock_top:raise RuntimeError('Insufficient passing molecules; no threshold relaxation')
  def admet():
   from admet_ai import ADMETModel
   model=ADMETModel(num_workers=0);df=model.predict(smiles=[r['smiles'] for r in top]);df.to_csv(w/'admet.csv');return {s:{k:float(df.loc[s,k]) for k in ['AMES','hERG','DILI']} for s in df.index}
  predictions=stage('admet_model_load_and_prediction',admet)
  failures=[]
  def docking():
   result=[]
   for j,r in enumerate(top):
    row={**r,'toxicity_prediction':predictions[r['smiles']]};m=Chem.AddHs(Chem.MolFromSmiles(r['smiles']));params=AllChem.ETKDGv3();params.randomSeed=args.seed
    try:
     if AllChem.EmbedMolecule(m,params)!=0:raise RuntimeError('embedding failed')
     if AllChem.MMFFOptimizeMolecule(m,mmffVariant='MMFF94s',maxIters=1000)!=0:raise RuntimeError('MMFF did not converge')
     setup=MoleculePreparation(rigid_macrocycles=True).prepare(m)[0];pdbqt,ok,err=PDBQTWriterLegacy.write_string(setup)
     if not ok:raise RuntimeError(str(err))
     lig=w/f'ligand_{j}.pdbqt';lig.write_text(pdbqt);pose=w/f'pose_{j}.pdbqt';cmd=['vina','--receptor',str(rec),'--ligand',str(lig),'--out',str(pose),'--seed',str(args.seed),'--cpu','1','--exhaustiveness','2']
     for axis,c,size in zip('xyz',pocket['center'],pocket['size']):cmd+=['--center_'+axis,str(c),'--size_'+axis,str(size)]
     done=subprocess.run(cmd,capture_output=True,text=True,timeout=300);(w/f'vina_{j}.log').write_text(done.stdout+'\n'+done.stderr);done.check_returncode()
     import re
     match=re.search(r'VINA RESULT:\s*(-?\d+\.\d+)',pose.read_text())
     if not match:raise RuntimeError('No Vina score')
     row['vina_score']=float(match.group(1));row['ligand_pdbqt_sha256']=sha(lig)
    except Exception as ex:row['error']=str(ex);failures.append({'smiles':r['smiles'],'error':str(ex)})
    result.append(row)
   return result
  scored=stage('ligand_preparation_and_docking',docking);timings['total']=time.perf_counter()-t0
  result={'repeat':repeat,'generated_nonempty':len(raw),'unique_druglike':len(cands),'pocket':pocket,'receptor_sha256':sha(rec),'selected':scored,'failures':failures,'timings_seconds':timings};(w/'result.json').write_text(json.dumps(result,indent=2));runs.append(result)
 keys=['generated_nonempty','unique_druglike','pocket','receptor_sha256','selected','failures'];checks={k:all(r[k]==runs[0][k] for r in runs[1:]) for k in keys};checks['raw_sample_sequences_equal']=all((args.output/f'run_{i}/raw_samples.json').read_bytes()==(args.output/'run_1/raw_samples.json').read_bytes() for i in range(2,args.repeat+1));checks['all_docking_completed']=all(not r['failures'] for r in runs)
 summary={'checks':checks,'all_checks_pass':all(checks.values()),'runs':runs,'caveat':'Same host/env and reused input weights; not cross-machine reproducibility. Timings include model load/cache effects, not a performance benchmark.'};(args.output/'comparison.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'checks':checks,'timings':[r['timings_seconds'] for r in runs]},indent=2),flush=True)
if __name__=='__main__':main()
