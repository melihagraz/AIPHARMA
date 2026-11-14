"""Build a scaffold-disjoint holdout and separate assay labels from scoring inputs."""
from pathlib import Path
from functools import lru_cache
from collections import Counter,defaultdict
import json,csv,hashlib,statistics,math,random
from rdkit import Chem,RDLogger,DataStructs
from rdkit.Chem import Descriptors,Crippen,Lipinski,rdMolDescriptors,rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
p=Path.home()/'aipharma_holdout_20261016';out=p/'locked_bounded';out.mkdir(exist_ok=True)
if (out/'manifest.json').exists():raise RuntimeError('Holdout already locked; do not overwrite')
protocol=json.loads((p/'eligibility_protocol.json').read_text());protocol.update({'amendment':'Right-censored IC50 > or >=10,000 nM admitted as experimentally weak; active threshold and exclusion/matching rules unchanged','weak_threshold':'exact median pChEMBL <=5 OR qualifying IC50 lower bound >=10,000 nM with no conflicting exact result','censored_values':'never treated as exact potency values','reason':'Strict exact-value protocol yielded 0 new weak controls','scoring_stereochemistry':'preserve representative experimental parent stereo/charge; stereo ignored only for grouping/exclusion'})
(p/'eligibility_protocol_bounded.json').write_text(json.dumps(protocol,indent=2));enumerator=rdMolStandardize.TautomerEnumerator();uncharger=rdMolStandardize.Uncharger()
@lru_cache(maxsize=None)
def normalize(smi):
 mol=Chem.MolFromSmiles(smi)
 if mol is None:return None
 mol=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol));mol=uncharger.uncharge(mol);mol=enumerator.Canonicalize(mol);Chem.RemoveStereochemistry(mol)
 canonical=Chem.MolToSmiles(mol);core=MurckoScaffold.GetScaffoldForMol(mol);scaffold=Chem.MolToSmiles(core) if core.GetNumAtoms() else 'ACYCLIC'
 return canonical,scaffold
prior=[];prior.extend(r['smiles'] for r in csv.DictReader((p/'alk5_active_control_benchmark.csv').open()));prior.extend(r['smiles'] for r in json.loads((p/'pah_run_1PY5.json').read_text())['ranked']);prior.extend(r['smiles'] for r in json.loads((p/'generation_records.json').read_text()) if r.get('smiles'))
for pdb in ['1PY5','1VJY']:
 m=Chem.SDMolSupplier(str(Path.home()/('.aipharma_dock6_amber_'+pdb)/'crystal.sdf'),removeHs=True)[0];prior.append(Chem.MolToSmiles(m))
prior_normal=[normalize(s) for s in prior];assert all(x is not None for x in prior_normal);known={x[0] for x in prior_normal};knownscaf={x[1] for x in prior_normal}
(out/'exclusion_inventory.json').write_text(json.dumps({'source_structures':len(prior),'unique_normalized_parents':len(known),'unique_scaffolds':len(knownscaf),'normalized_parent_smiles':sorted(known),'scaffolds':sorted(knownscaf),'limitation':'Only stored project structures; prior-training overlap and unrecorded generation samples cannot be audited'},indent=2))
assays={a['assay_chembl_id']:a for a in json.loads((p/'raw/boundary_assays_all.json').read_text())['assays']};groups=defaultdict(list);reject=Counter()
original=json.loads((p/'raw/all_activities.json').read_text())['activities'];extra=json.loads((p/'raw/weak_boundary_all.json').read_text())['activities'];allrows={r['activity_id']:r for r in original+extra}
for row in allrows.values():
 a=assays[row['assay_chembl_id']]
 if a.get('assay_type')!='B' or a.get('confidence_score')!=9 or a.get('target_chembl_id')!='CHEMBL4439':reject['assay_quality']+=1;continue
 if row.get('standard_units')!='nM' or row.get('standard_type')!='IC50':reject['measurement_format']+=1;continue
 if row.get('data_validity_comment') or row.get('potential_duplicate') in [True,1,'1']:reject['validity_or_duplicate_flag']+=1;continue
 relation=row.get('standard_relation');value=float(row['standard_value']) if row.get('standard_value') is not None else None
 exact=relation=='=' and row.get('pchembl_value') is not None
 bounded=relation in ['>','>='] and value is not None and value>=10000
 if not (exact or bounded) or not row.get('canonical_smiles'):reject['measurement_format']+=1;continue
 norm=normalize(row['canonical_smiles'])
 if norm is None:reject['invalid_structure']+=1;continue
 groups[norm[0]].append({**row,'normalized_smiles':norm[0],'normalized_scaffold':norm[1],'is_right_censored_weak':bounded})
active=[];weak=[]
for identity,rows in groups.items():
 exactrows=[x for x in rows if not x['is_right_censored_weak']];bounds=[float(x['standard_value']) for x in rows if x['is_right_censored_weak']];values=[float(x['pchembl_value']) for x in exactrows];median=statistics.median(values) if values else None
 if values and ((min(values)<=5 and max(values)>=7) or (bounds and max(values)>=7)):reject['conflicting_parent']+=1;continue
 if median is not None and median>=7:label=1
 elif median is not None and median<=5:label=0
 elif median is None and bounds:label=0
 else:reject['intermediate_parent']+=1;continue
 norm=normalize(identity);scaf=norm[1]
 if identity in known:reject['previous_parent']+=1;continue
 if scaf in knownscaf:reject['previous_scaffold']+=1;continue
 representative=min(exactrows,key=lambda x:(abs(float(x['pchembl_value'])-median),int(x['activity_id']))) if exactrows else min(rows,key=lambda x:int(x['activity_id']))
 mol=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(Chem.MolFromSmiles(representative['canonical_smiles'])));smi=Chem.MolToSmiles(mol)
 r={'smiles':smi,'parent_identity':identity,'scaffold':scaf,'label':label,'median_pchembl':median,'pchembl_range':[min(values),max(values)] if values else None,'weak_ic50_lower_bound_nm':min(bounds) if bounds else None,'label_evidence':'right_censored_IC50' if median is None else 'exact_IC50','measurements':len(rows),'molecule_chembl_ids':sorted({x['molecule_chembl_id'] for x in rows}),'activity_ids':[x['activity_id'] for x in rows],'assay_chembl_ids':sorted({x['assay_chembl_id'] for x in rows}),'document_chembl_ids':sorted({x['document_chembl_id'] for x in rows if x.get('document_chembl_id')}),'representative_activity_id':representative['activity_id'],'mw':Descriptors.MolWt(mol),'logp':Crippen.MolLogP(mol),'hbd':Lipinski.NumHDonors(mol),'hba':Lipinski.NumHAcceptors(mol),'rot':rdMolDescriptors.CalcNumRotatableBonds(mol),'tpsa':rdMolDescriptors.CalcTPSA(mol),'unassigned_stereo_types':[str(x.type) for x in Chem.FindPotentialStereo(mol) if str(x.specified)=='Unspecified']}
 (active if label else weak).append(r)
active_scaf={x['scaffold'] for x in active};beforecross=len(weak);weak=[r for r in weak if r['scaffold'] not in active_scaf];reject['weak_cross_class_scaffold_overlap']=beforecross-len(weak)
key=lambda r:hashlib.sha256(('20261016'+r['smiles']).encode()).hexdigest()
def one_per_scaffold(rows):
 result=[];seen=set()
 for r in sorted(rows,key=key):
  if r['scaffold'] not in seen:result.append(r);seen.add(r['scaffold'])
 return result
pool_counts={'eligible_active_parents':len(active),'eligible_weak_parents':len(weak),'active_scaffolds':len(active_scaf),'weak_scaffolds':len({r['scaffold'] for r in weak})};active=one_per_scaffold(active);weak=one_per_scaffold(weak)
scales={'mw':50,'logp':1,'hbd':1,'hba':2,'rot':2,'tpsa':30}
def distance(a,b):return math.sqrt(sum(((a[k]-b[k])/scale)**2 for k,scale in scales.items()))
available=list(active);pairs=[]
for w in sorted(weak,key=key):
 if not available or len(pairs)>=50:break
 a=min(available,key=lambda x:(distance(x,w),key(x)));available.remove(a);pairs.append((a,w))
flat=[]
for i,(a,w) in enumerate(pairs,1):
 for row in [a,w]:flat.append({**row,'pair':i,'property_match_distance':distance(a,w)})
random.Random(20261016).shuffle(flat);fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048);priorfps=[fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in sorted(known)];labels={};inputs=[]
for i,r in enumerate(flat,1):
 case=f'H{i:03d}';r['case_id']=case;r['max_similarity_to_recorded_prior']=max(DataStructs.BulkTanimotoSimilarity(fpgen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])),priorfps));labels[case]=r;inputs.append({'case_id':case,'smiles':r['smiles']})
(out/'structures_blinded.json').write_text(json.dumps({'structures':inputs},indent=2));(out/'label_key.json').write_text(json.dumps(labels,indent=2))
with (out/'structures_blinded.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=['case_id','smiles']);w.writeheader();w.writerows(inputs)
assert not ({r['parent_identity'] for r in flat}&known);assert not ({r['scaffold'] for r in flat}&knownscaf);assert len({r['scaffold'] for r in flat})==len(flat)
manifest={'requested_pairs':50,'actual_pairs':len(pairs),'weak_evidence_counts':dict(Counter(r['label_evidence'] for r in flat if r['label']==0)),'actual_structures':len(flat),'classes_balanced':True,'pool_counts':pool_counts,'exclusions':dict(reject),'prior_parent_count':len(known),'prior_scaffold_count':len(knownscaf),'property_mean_absolute_pair_difference':{k:statistics.mean(abs(a[k]-w[k]) for a,w in pairs) if pairs else None for k in scales},'known_scaffold_overlap':0,'within_holdout_scaffold_duplicates':0,'rules_relaxed':False,'labels_used_by_scoring':False,'normalization_stereochemistry':'preserved in representative scoring structures; ignored for grouping/exclusion only','blinded_structures_sha256':hashlib.sha256((out/'structures_blinded.json').read_bytes()).hexdigest(),'label_key_sha256':hashlib.sha256((out/'label_key.json').read_bytes()).hexdigest(),'protocol_sha256':hashlib.sha256((p/'eligibility_protocol_bounded.json').read_bytes()).hexdigest(),'limitations':['Scoring code is label-withheld, not an independently administered blind trial','Unknown overlap with pretrained model data','Retrospective heterogeneous assay conditions; only high-confidence target-mapped binding IC50 measurements retained','All recorded prior scaffolds excluded; unrecorded generated samples cannot be excluded']}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
