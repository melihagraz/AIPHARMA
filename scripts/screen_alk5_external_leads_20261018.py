"""Screen external weak-IC50 leads for parent, scaffold and publication overlap.

This is an inventory, NOT an endpoint-validated control set.
"""
from pathlib import Path
from functools import lru_cache
from collections import Counter
import concurrent.futures
import csv
import hashlib
import json
import re
import requests
from rdkit import Chem,RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog('rdApp.*')
out=Path.home()/'aipharma_external_leads_20261018'
old=Path.home()/'aipharma_holdout_20261016'
cache=out/'known_documents';cache.mkdir(exist_ok=True)
raw={r['activity_id']:r for name in ['all_activities.json','weak_boundary_all.json'] for r in json.loads((old/'raw'/name).read_text())['activities']}
ids=sorted({r['document_chembl_id'] for r in raw.values() if r.get('document_chembl_id')})

def fetch_documents(start):
 dest=cache/f'{start:04d}.json'
 params={'document_chembl_id__in':','.join(ids[start:start+40]),'limit':1000}
 url='https://www.ebi.ac.uk/chembl/api/data/document.json'
 if dest.exists():return json.loads(dest.read_text())['documents']
 response=requests.get(url,params=params,timeout=45);response.raise_for_status();body=response.json();dest.write_bytes(response.content)
 return body['documents']

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:docs=[d for batch in pool.map(fetch_documents,range(0,len(ids),40)) for d in batch]
assert set(ids)=={d['document_chembl_id'] for d in docs}
(out/'known_document_registry.json').write_text(json.dumps(docs,indent=2))
doi_map={str(d['doi']).lower():d['document_chembl_id'] for d in docs if d.get('doi')}
pmid_map={str(d['pubmed_id']):d['document_chembl_id'] for d in docs if d.get('pubmed_id')}
enum=rdMolStandardize.TautomerEnumerator();uncharger=rdMolStandardize.Uncharger()

@lru_cache(maxsize=None)
def normalize(smi):
 m=Chem.MolFromSmiles(smi)
 if m is None:return None
 m=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m));m=uncharger.uncharge(m);m=enum.Canonicalize(m);Chem.RemoveStereochemistry(m)
 return Chem.MolToSmiles(m),MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False)

known_parents=set();known_scaffolds=set()
for r in raw.values():
 if r.get('canonical_smiles'):
  value=normalize(r['canonical_smiles'])
  if value:known_parents.add(value[0]);known_scaffolds.add(value[1])
print('normalized prior structures',len(known_parents),flush=True)
bdb=json.loads((out/'bindingdb_high_cutoff.json').read_text())
affinities=next(iter(bdb.values()))['affinities']
weak=[];strong_ids=set()
for r in affinities:
 if r['affinity_type']!='IC50':continue
 match=re.fullmatch(r'\s*([<>=~]*)\s*([\d.]+(?:[eE][+-]?\d+)?)\s*',r['affinity'])
 if not match:continue
 relation=match[1] or '=';value=float(match[2])
 if value<=100 and relation in ['=','<','<=']:strong_ids.add(r['monomerid'])
 if value<10000 or relation not in ['=','>','>=']:continue
 norm=normalize(r['smile'])
 if norm is None:continue
 weak.append({**r,'relation':relation,'value_nM':value,'normalized_parent':norm[0],'scaffold':norm[1],
              'parent_in_prior_raw':norm[0] in known_parents,'scaffold_in_prior_raw':bool(norm[1] and norm[1] in known_scaffolds),
              'source_in_prior_raw':r.get('doi','').lower() in doi_map or str(r.get('pmid','')) in pmid_map,
              'prior_document_chembl_id':doi_map.get(r.get('doi','').lower()) or pmid_map.get(str(r.get('pmid',''))),
              'source_is_bindingdb_dataset_doi':r.get('doi','').lower().startswith('10.7270/'),
              'endpoint_status':'unverified; no catalytic-versus-binding/cellular assumption from target ID'})
for r in weak:r['same_monomer_has_strong_IC50_in_external_query']=r['monomerid'] in strong_ids
with (out/'weak_lead_screen.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(weak[0]));w.writeheader();w.writerows(weak)
(out/'weak_lead_screen.json').write_text(json.dumps(weak,indent=2))
pc=list(csv.DictReader((out/'pubchem_1799383.csv').open()))
pc=[r for r in pc if r['Target Accession(s)']=='P36897' and r['Standard Type']=='IC50']
pc_records=[]
for r in pc:
 norm=normalize(r['PUBCHEM_EXT_DATASOURCE_SMILES']);value=float(r['IC50']);relation=r['IC50 Qualifier'] or '='
 pc_records.append({'cid':r['PUBCHEM_CID'],'bindingdb_id':r['Ligand'],'smiles':r['PUBCHEM_EXT_DATASOURCE_SMILES'],
                    'normalized_parent':norm[0],'value_nM':value,'relation':relation,'parent_in_prior_raw':norm[0] in known_parents,
                    'scaffold_in_prior_raw':bool(norm[1] and norm[1] in known_scaffolds),
                    'doi':'10.1021/cb9002865','source_in_prior_raw':'10.1021/cb9002865' in doi_map,
                    'conservative_status':'needs original table consistency review before acceptance'})
(out/'pubchem_target_only_review.json').write_text(json.dumps(pc_records,indent=2))
fresh=[r for r in weak if not r['parent_in_prior_raw'] and not r['source_in_prior_raw'] and not r['source_is_bindingdb_dataset_doi']]
summary={'bindingdb_total_measurement_rows':len(affinities),'prior_raw_activity_rows':len(raw),'prior_document_count':len(docs),
         'prior_parent_count':len(known_parents),'weak_IC50_rows':len(weak),'weak_distinct_parents':len({r['normalized_parent'] for r in weak}),
         'exact_weak_rows':sum(r['relation']=='=' for r in weak),'lower_bound_weak_rows':sum(r['relation']!='=' for r in weak),
         'new_parent_weak_rows':sum(not r['parent_in_prior_raw'] for r in weak),
         'new_parent_and_new_literature_source_rows':len(fresh),'new_parent_and_new_literature_source_distinct_parents':len({r['normalized_parent'] for r in fresh}),
         'new_literature_dois':sorted({r['doi'] for r in fresh}),'pubchem_target_rows':len(pc_records),
         'accepted_new_controls':0,'scope':'Public lead discovery only; no docking scores; not assay-validated negatives or independent blocks yet',
         'doi_caveat':'10.7270 DOI identifies BindingDB dataset; original publication/patent identity must be resolved before source novelty claims'}
(out/'screen_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
