"""Rebuild source-family/scaffold components after independent name mapping."""
from pathlib import Path
import concurrent.futures
import gzip
import hashlib
import json
import re
import requests
from html.parser import HTMLParser
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from audit_alk5_independent_blocks_20261018 import graph

home=Path.home();base=home/'aipharma_external_leads_20261018';out=home/'aipharma_structure_mapping_20261018'
registry={r['document_chembl_id']:r for r in json.loads((base/'known_document_registry.json').read_text())}
old=json.loads((home/'aipharma_block_audit_final_20261018/history_and_stereo_eligible_nodes.json').read_text())
mapped=json.loads((out/'mapping_results.json').read_text())
new=[r for r in mapped if r['mapping_status']=='parent_verified']
patents=sorted({registry[d]['patent_id'].replace('-','') for r in old for d in r['documents'] if registry[d].get('patent_id')})
cache=out/'old_patent_families';cache.mkdir(exist_ok=True)
class Text(HTMLParser):
 def __init__(self):super().__init__();self.parts=[];self.skip=0
 def handle_starttag(self,t,a):
  if t in ['script','style']:self.skip+=1
 def handle_endtag(self,t):
  if t in ['script','style']:self.skip=max(0,self.skip-1)
 def handle_data(self,d):
  if not self.skip and d.strip():self.parts.append(d.strip())
def fetch(ident):
 path=cache/(ident+'.html.gz');url=f'https://patents.google.com/patent/{ident}/en'
 try:
  if path.exists():raw=gzip.decompress(path.read_bytes())
  else:
   response=requests.get(url,timeout=25);response.raise_for_status();raw=response.content;path.write_bytes(gzip.compress(raw))
  parser=Text();parser.feed(raw.decode());families=sorted(set(re.findall(r'Family ID[=: ]+\s*(\d+)',' '.join(parser.parts))))
  return {'patent':ident,'family_ids':families,'status':'complete' if len(families)==1 else 'unresolved','url':url,'raw_sha256':hashlib.sha256(raw).hexdigest()}
 except Exception as e:return {'patent':ident,'status':'failed','error':str(e),'url':url}
print('Checking prior patent families',len(patents),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:prior=list(pool.map(fetch,patents))
(out/'old_patent_family_log.json').write_text(json.dumps(prior,indent=2))
families={r['patent']:r['family_ids'][0] for r in prior if r['status']=='complete'}
new_archive=json.loads((base/'patent_archive_log.json').read_text())
families.update({r['patent']:r['family_ids'][0] for r in new_archive if r['status']=='complete'})
def source(d):
 r=registry[d]
 if r.get('patent_id'):
  patent=r['patent_id'].replace('-','')
  return 'patent-family:'+families[patent] if patent in families else 'unresolved-patent:'+patent
 if r.get('doi'):return 'doi:'+r['doi'].lower()
 if r.get('pubmed_id'):return 'pmid:'+str(r['pubmed_id'])
 return 'chembl:'+d
old_rows=[{**r,'original_documents':r['documents'],'documents':sorted({source(d) for d in r['documents']}),'origin':'prior_exact_eligible'} for r in old]
new_rows=[]
for r in new:
 smi=r['bindingdb_parent'];scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(smi),includeChirality=False)
 new_rows.append({'smiles':smi,'scaffold':scaffold,'label':0,'documents':sorted({'patent-family:'+families[s] for s in r['patent_sources']}),
                  'origin':'source_parent_mapping_verified','monomer':r['monomer'],'original_smiles':r['original_smiles'],
                  'mapping_status':r['mapping_status'],'weak_values_nM':r['weak_values_nM'],'benchmark_accepted':False})
assert not ({r['smiles'] for r in old_rows}&{r['smiles'] for r in new_rows})
all_rows=old_rows+new_rows
(out/'expanded_nodes.json').write_text(json.dumps(all_rows,indent=2))
(out/'new_mapped_nodes.json').write_text(json.dumps(new_rows,indent=2))
summary={'prior_molecules':len(old_rows),'new_parent_verified_weak_leads':len(new_rows),'quarantined':len(mapped)-len(new_rows),
         'prior_patent_family_requests':len(prior),'prior_patent_family_complete':sum(r['status']=='complete' for r in prior),
         'new_patent_family_count':len({d for r in new_rows for d in r['documents']}),'mapping_alone_is_not_assay_validation':True,
         'score_selection_used':False,'benchmark_promotions':0}
for pool_name,rows in [('new_only',new_rows),('expanded',all_rows)]:
 for mode,threshold in [('source_family_exact_scaffold',None),('source_family_scaffold_similarity_040',.4)]:
  data=graph(rows,threshold)
  (out/f'{pool_name}_{mode}_graph.json').write_text(json.dumps(data,indent=2))
  summary[pool_name+'_'+mode]={k:v for k,v in data.items() if k not in ['components','edges']}
  if pool_name=='expanded':
   summary[pool_name+'_'+mode]['components_with_new_and_prior']=sum(any(i<len(old_rows) for i in b['indices']) and any(i>=len(old_rows) for i in b['indices']) for b in data['components'])
(out/'expanded_graph_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
