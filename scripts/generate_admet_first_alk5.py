"""Prospectively specified ADMET-first ALK5 candidate selection.
No generator retraining or claim of experimentally established safety/activity.
"""
from pathlib import Path
import sys,json,csv,hashlib,importlib.metadata,math,time
from collections import Counter
root=Path.home()/'aipharma_generation_20261016';sys.path.insert(0,str(root/'Reinvent_App'));sys.path.insert(0,str(root))
import torch
from rdkit import Chem,RDLogger,DataStructs
from rdkit.Chem import Descriptors,Crippen,QED,Lipinski,rdMolDescriptors,rdFingerprintGenerator,FilterCatalog
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDConfig
sys.path.insert(0,str(Path(RDConfig.RDContribDir)/'SA_Score'))
import sascorer
from app.sampler import ReinventSampler
from pah_filter_profile import PAH_FILTER_PROFILE
RDLogger.DisableLog('rdApp.*');torch.set_num_threads(4)
seeds=[101,102,103,104,105];n=400;out=root/'results';out.mkdir(exist_ok=True)
if (out/'locked_top20.json').exists():raise RuntimeError('Selection already locked; use a new run directory rather than overwrite')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
versions={x:importlib.metadata.version(x) for x in ['torch','rdkit','admet-ai','chemprop']}
protocol={'run':'20261016_admet_first','seeds':seeds,'requested_per_seed':n,'generator':'unconditional REINVENT base prior; downstream target-informed selection','prior_sha256':sha(root/'Reinvent_App/prior/reinvent.prior'),'versions':versions,'base_filters':PAH_FILTER_PROFILE['base_filters'],'toxicity_thresholds':{'AMES':'<0.5','hERG':'<0.5','DILI':'<0.5'},'single_fragment':True,'exclude_unassigned_stereocenters':True,'PAINS':'hard reject','BRENK':'annotation only','exclude_exact_matches_to':'previous 12 candidates and 70 ALK5 benchmark actives','selection':'descending max Morgan radius-2/2048-bit Tanimoto to benchmark ALK5 actives, ascending max toxicity prediction, descending QED, canonical SMILES; max 2 per Murcko scaffold','selection_count_max':20,'threshold_relaxation_allowed':False,'benchmark_sha256':sha(root/'alk5_active_control_benchmark.csv'),'old_candidates_sha256':sha(root/'pah_run_1PY5.json'),'limitations':['Predictions are selection criteria, not experimental safety validation','No ALK5-conditioned generator training','Not exhaustive tautomer/protonation enumeration','Similarity is a selection heuristic, not binding evidence']}
(out/'protocol.json').write_text(json.dumps(protocol,indent=2));print('Protocol saved before generation: '+sha(out/'protocol.json'),flush=True)
canon=lambda smi:Chem.MolToSmiles(Chem.MolFromSmiles(smi))
old=json.loads((root/'pah_run_1PY5.json').read_text())['ranked'];excluded={canon(r['smiles']) for r in old}
active_rows=[r for r in csv.DictReader((root/'alk5_active_control_benchmark.csv').open()) if r['label']=='1'];excluded.update(canon(r['smiles']) for r in active_rows)
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);activefps=[fpgen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in active_rows]
def catalog(kind):
 p=FilterCatalog.FilterCatalogParams();p.AddCatalog(kind);return FilterCatalog.FilterCatalog(p)
pains=catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS);brenk=catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
base=PAH_FILTER_PROFILE['base_filters'];sampler=ReinventSampler();seen=set();records=[];eligible=[];counts=Counter();per_seed=[]
for seed in seeds:
 local=Counter();raw=sampler.sample(batch_size=n,seed=seed,target_prior='base')
 for sample in raw:
  row={'seed':seed,'raw_smiles':sample.smiles,'nll':float(sample.nll)};counts['returned']+=1;local['returned']+=1
  mol=Chem.MolFromSmiles(sample.smiles) if sample.smiles else None
  if mol is None:row['status']='invalid';records.append(row);counts['invalid']+=1;continue
  smi=Chem.MolToSmiles(mol);row['smiles']=smi;counts['valid']+=1;local['valid']+=1
  if smi in seen:row['status']='duplicate';records.append(row);counts['duplicate']+=1;continue
  seen.add(smi);counts['unique_valid']+=1
  if smi in excluded:row['status']='existing_exact_match';records.append(row);counts['existing_exact_match']+=1;continue
  props={'mw':Descriptors.MolWt(mol),'logp':Crippen.MolLogP(mol),'hbd':Lipinski.NumHDonors(mol),'hba':Lipinski.NumHAcceptors(mol),'rot':rdMolDescriptors.CalcNumRotatableBonds(mol),'tpsa':rdMolDescriptors.CalcTPSA(mol),'qed':QED.qed(mol),'sa':sascorer.calculateScore(mol)}
  row.update(props)
  ok=props['mw']<=base['mw_max'] and props['logp']<=base['logp_max'] and props['hbd']<=base['hbd_max'] and props['hba']<=base['hba_max'] and props['rot']<=base['rotatable_bonds_max'] and props['tpsa']<=base['tpsa_max'] and props['qed']>=base['qed_min'] and props['sa']<=base['sa_score_max']
  if not ok:row['status']='physchem_reject';records.append(row);counts['physchem_reject']+=1;continue
  if len(Chem.GetMolFrags(mol))!=1:row['status']='fragment_reject';records.append(row);counts['fragment_reject']+=1;continue
  if any(label=='?' for _,label in Chem.FindMolChiralCenters(mol,includeUnassigned=True)):row['status']='stereo_reject';records.append(row);counts['stereo_reject']+=1;continue
  if pains.HasMatch(mol):row['status']='PAINS_reject';records.append(row);counts['PAINS_reject']+=1;continue
  row['BRENK_alerts']=[match.GetDescription() for match in brenk.GetMatches(mol)]
  sims=DataStructs.BulkTanimotoSimilarity(fpgen.GetFingerprint(mol),activefps);idx=max(range(len(sims)),key=sims.__getitem__);row['max_alk5_active_tanimoto']=sims[idx];row['nearest_benchmark_active']=active_rows[idx]['molecule_chembl_id'];row['scaffold']=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)) or 'ACYCLIC'
  row['status']='eligible_for_ADMET';eligible.append(row);records.append(row);counts['eligible_for_ADMET']+=1;local['eligible_for_ADMET']+=1
 per_seed.append({'seed':seed,**dict(local)});(out/'generation_records.json').write_text(json.dumps(records,indent=2));print(json.dumps(per_seed[-1]),flush=True)
(out/'generation_counts.json').write_text(json.dumps({'counts':dict(counts),'per_seed':per_seed},indent=2))
print('Loading existing ADMET-AI '+versions['admet-ai']+'; eligible='+str(len(eligible)),flush=True)
from admet_ai import ADMETModel
model=ADMETModel(num_workers=0)
baseline=model.predict(smiles=[r['smiles'] for r in old]);baseline.to_csv(out/'baseline12_admet_recomputed.csv')
passing=[];predictions=[]
for start in range(0,len(eligible),128):
 batch=eligible[start:start+128];df=model.predict(smiles=[r['smiles'] for r in batch]);df.to_csv(out/f'admet_batch_{start//128:03d}.csv')
 if set(df.index)!=set(r['smiles'] for r in batch):raise RuntimeError('ADMET input/output identity mismatch')
 for row in batch:
  values={k:float(df.loc[row['smiles'],k]) for k in ['AMES','hERG','DILI']}
  if not all(math.isfinite(v) and 0<=v<=1 for v in values.values()):raise RuntimeError('Invalid toxicity output')
  result={**row,'toxicity':values,'max_toxicity':max(values.values()),'all_three_toxicity_pass':all(v<0.5 for v in values.values())};predictions.append(result)
  if result['all_three_toxicity_pass']:passing.append(result)
 (out/'admet_screened.json').write_text(json.dumps(predictions,indent=2));print(f'ADMET {len(predictions)}/{len(eligible)}; all-three pass={len(passing)}',flush=True)
passing.sort(key=lambda r:(-r['max_alk5_active_tanimoto'],r['max_toxicity'],-r['qed'],r['smiles']))
selected=[];scaffold_counts=Counter()
for r in passing:
 if scaffold_counts[r['scaffold']]>=2:continue
 r={**r,'selection_order':len(selected)+1};selected.append(r);scaffold_counts[r['scaffold']]+=1
 if len(selected)==20:break
baseline_count=int(((baseline[['AMES','hERG','DILI']]<0.5).all(axis=1)).sum())
summary={'protocol_sha256':sha(out/'protocol.json'),'requested':len(seeds)*n,'generation_counts':dict(counts),'per_seed':per_seed,'admet_screened':len(predictions),'all_three_toxicity_pass':len(passing),'all_three_pass_fraction_among_screened':len(passing)/len(predictions) if predictions else None,'recomputed_previous12_all_three_pass':baseline_count,'selected_for_docking':len(selected),'selected_scaffolds':len(scaffold_counts),'target_conditioned_model':False,'thresholds_relaxed':False,'docking_inspected_before_selection':False,'limitations':protocol['limitations']}
(out/'selection_summary.json').write_text(json.dumps(summary,indent=2));(out/'locked_top20.json').write_text(json.dumps({'protocol_sha256':sha(out/'protocol.json'),'selection_summary':summary,'candidates':selected},indent=2));(out/'locked_top20.sha256').write_text(sha(out/'locked_top20.json')+'\n')
print(json.dumps(summary,indent=2),flush=True)
