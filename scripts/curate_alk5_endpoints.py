"""Assay-description endpoint audit; fixed rules before outcome recomputation."""
from pathlib import Path
import json,re,hashlib,statistics,csv
from collections import Counter,defaultdict
from functools import lru_cache
from rdkit import Chem,RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
p=Path.home()/'aipharma_holdout_20261016';o=p/'endpoint_audit';o.mkdir(exist_ok=True)
protocol={'purpose':'Correct mixed endpoints; exploratory historical reanalysis, not new blind test','classification':'description evidence with recorded reasons; BAO and assay_type alone not decisive','groups':['enzyme_catalytic','biochemical_binding','biochemical_unspecified','cell_signaling','cell_target_engagement','cell_phenotype','uncertain'],'primary':'enzyme_catalytic; explicit mutants separated from no-explicit-mutant descriptions','mutations':'T204D and any other explicit mutant flagged; not-stated is NOT confirmed wild type','thresholds':'active median exact pChEMBL>=7; weak <=5 or right censored bound>=10000nM without contradictory exact measurements','pooling':'only within endpoint x mutation stratum; no cross-stratum medians','history':'all36 revealed holdout scaffolds excluded from model comparison','evaluation':'fixed Morgan logistic C1 balanced and max training-active similarity, same fivefold source+scaffold grouping; no parameter tuning','unknown':'excluded from model comparison; retained in audit','limitations':'Description-based classification, not experimental revalidation; truncated or inconsistent assay annotations remain. Confidence9 is target-mapping confidence only.'}
(o/'protocol.json').write_text(json.dumps(protocol,indent=2))
assays=json.loads((p/'raw/boundary_assays_all.json').read_text())['assays'];classes={}
def classify(d):
 t=d.lower()
 if 'nanobret' in t:return 'cell_target_engagement','NanoBRET; treated as cellular engagement, distinct from signal response'
 if any(k in t for k in ['cell growth','cell proliferation']):return 'cell_phenotype','Cell growth/proliferation endpoint, not direct kinase potency'
 if any(k in t for k in ['reporter gene','luciferase reporter','cell-based luciferase','cell-based luciferase','on-target cellular','signal transduction inhibitory','in human t cells','human primary t-cells','in nhlf','human nhlf','in human treg','in anti-cd3','in human a549','in human hepg2','in human mda','in human hacat','in hek293','in human hek293','in human hs-578t','in tgf-beta-stimulated human hepg2','in tgfbetar1 in human','in tgfbeta-stimulated human']):return 'cell_signaling','Explicit intact-cell signaling / reporter description'
 if 'cells' in t and any(k in t for k in ['luciferase','nuclear translocation','smad','tgfbeta-stimulated'] ) and not ('as substrate' in t or 'using casein' in t):return 'cell_signaling','Cell context with downstream signal readout'
 if 'expressed in hek293' in t:return 'uncertain','HEK293 expression alone does not specify purified-protein versus intact-cell readout'
 if any(k in t for k in ['binding affinity','displacement','binding assay','kinome scan','kinomescan','lathascreen']):return 'biochemical_binding','Explicit binding/tracer displacement; not catalytic IC50'
 if any(k in t for k in ['as substrate','peptide substrate','substrate and atp','presence of atp','presence of 1 microm atp','presence of 10 um atp','gamma atp','gammaatp','gamma-33p','33 p-gamma','[33p]','(33p)','kinase assay','kinase assays','kinase glo','adp-glo','enzyme kinetic','enzymatic activity','enzymatic assay','enyzmatic assay','autophosphorylation','radioactive kinase','radioisotopic protein kinase','biochemical assays (reaction biology)']):return 'enzyme_catalytic','Explicit ATP/substrate/enzyme activity/kinase assay evidence'
 if any(k in t for k in ['recombinant','kinase domain','htrf','tr-fret','gst-tagged','gst-fused','gst fused','tagged','expressed in sf9','expressed in baculovirus','fluorescen','phosphorylation expressed in sf9']):return 'biochemical_unspecified','Protein/preparation/readout evidence; catalytic versus binding mechanism insufficiently explicit'
 return 'uncertain','No sufficient endpoint evidence in description'
for a in assays:
 if a.get('target_chembl_id')!='CHEMBL4439' or a.get('assay_type')!='B' or a.get('confidence_score')!=9:continue
 group,reason=classify(a.get('description',''))
 if a['assay_chembl_id']=='CHEMBL1837693':group,reason='uncertain','Description names p38alpha despite ALK5 target mapping; exclude pending source review'
 if a['assay_chembl_id']=='CHEMBL4714986':group,reason='enzyme_catalytic','Manual description review: biotin-labelled peptide substrate and gamma33P ATP'
 if a['assay_chembl_id']=='CHEMBL814456':group,reason='uncertain','Autophosphorylation in Sf9 cells: intact-cell versus protein-production context ambiguous'
 t=a.get('description','').lower();mut='T204D' if 't204d' in t or 't204d mutation' in t else 'other_explicit_mutant' if 'mutant' in t or 'mutation' in t else 'no_explicit_mutant'
 classes[a['assay_chembl_id']]={'assay_id':a['assay_chembl_id'],'endpoint':group,'mutation':mut,'reason':reason,'description':a.get('description'),'bao_format':a.get('bao_format')}
# Sanity checks target known pitfalls, not keyword mirroring.
assert classes['CHEMBL4361670']['endpoint']=='enzyme_catalytic'
assert classes['CHEMBL4270419']['endpoint']=='cell_signaling'
assert classes['CHEMBL831709']['endpoint']=='biochemical_binding'
assert classes['CHEMBL4188411']['mutation']=='T204D'
assert classes['CHEMBL899823']['endpoint']=='uncertain'
(o/'assay_classification.json').write_text(json.dumps(list(classes.values()),indent=2))
with (o/'assay_classification.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(next(iter(classes.values()))));w.writeheader();w.writerows(classes.values())
raw={r['activity_id']:r for fn in ['all_activities.json','weak_boundary_all.json'] for r in json.loads((p/'raw'/fn).read_text())['activities']};e=rdMolStandardize.TautomerEnumerator();u=rdMolStandardize.Uncharger()
@lru_cache(maxsize=None)
def norm(s):
 m=Chem.MolFromSmiles(s)
 if m is None:return None
 m=e.Canonicalize(u.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))));Chem.RemoveStereochemistry(m);return Chem.MolToSmiles(m),Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) or 'ACYCLIC'
old=json.loads((p/'locked_bounded/label_key.json').read_text());oldscaf={r['scaffold'] for r in old.values()};groups=defaultdict(list);rejected=Counter();activityaudit=[]
for r in raw.values():
 if r['assay_chembl_id'] not in classes:continue
 if r.get('standard_type')!='IC50' or r.get('standard_units')!='nM' or r.get('data_validity_comment') or r.get('potential_duplicate') in [True,1,'1']:continue
 exact=r.get('standard_relation')=='=' and r.get('pchembl_value') is not None;bounded=r.get('standard_relation') in ['>','>='] and float(r.get('standard_value') or 0)>=10000
 if not (exact or bounded) or not r.get('canonical_smiles'):continue
 n=norm(r['canonical_smiles'])
 if not n:continue
 c=classes[r['assay_chembl_id']];key=c['endpoint']+'__'+c['mutation'];groups[(key,n[0])].append((r,exact,n[1]));activityaudit.append({'activity_id':r['activity_id'],'assay_id':r['assay_chembl_id'],'stratum':key,'smiles':n[0],'old_holdout_scaffold':n[1] in oldscaf})
datasets=defaultdict(list)
for (key,smi),rr in sorted(groups.items()):
 vals=[float(r['pchembl_value']) for r,x,s in rr if x];bounds=[float(r['standard_value']) for r,x,s in rr if not x];med=statistics.median(vals) if vals else None
 if vals and ((min(vals)<=5 and max(vals)>=7) or (bounds and max(vals)>=7)):rejected['conflicting_within_stratum']+=1;continue
 label=1 if med is not None and med>=7 else 0 if (med is not None and med<=5) or (med is None and bounds) else None
 if label is None:rejected['intermediate']+=1;continue
 row={'smiles':smi,'scaffold':rr[0][2],'label':label,'median_pchembl':med,'evidence':'exact' if vals else 'censored','documents':sorted({r['document_chembl_id'] for r,x,s in rr if r.get('document_chembl_id')}),'activity_ids':[r['activity_id'] for r,x,s in rr],'assay_ids':sorted({r['assay_chembl_id'] for r,x,s in rr}),'stratum':key,'old_holdout_scaffold':rr[0][2] in oldscaf};datasets[key].append(row)
(o/'activity_assignment.json').write_text(json.dumps(activityaudit,indent=2));(o/'datasets_all.json').write_text(json.dumps(datasets,indent=2))
summary={'assay_count':len(classes),'assay_endpoint_counts':dict(Counter(r['endpoint'] for r in classes.values())),'retained_measurement_count':len(activityaudit),'strata':{k:{'all_parents':len(v),'all_class_counts':dict(Counter(r['label'] for r in v)),'model_eligible_after_old_holdout_scaffold_exclusion':len([r for r in v if not r['old_holdout_scaffold']]),'model_class_counts':dict(Counter(r['label'] for r in v if not r['old_holdout_scaffold']))} for k,v in datasets.items()},'rejected':dict(rejected),'classification_limit':'Manual review of rule rationale needed before publication. BAO formats are not treated as ground truth.'}
(o/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
