"""Fetch an independent activity pool and lock holdout eligibility rules first."""
from pathlib import Path
import json,urllib.request,urllib.parse,time,hashlib
p=Path.home()/'aipharma_holdout_20261016';raw=p/'raw';raw.mkdir(exist_ok=True)
protocol={'target_chembl_id':'CHEMBL4439','target_organism':'Homo sapiens','endpoint':'IC50 only','standard_relation':'=','units':'nM','pchembl_required':True,'active_threshold':'median exact IC50 pChEMBL >=7','weak_threshold':'median exact IC50 pChEMBL <=5','conflicts':'exclude molecules with measurements in both >=7 and <=5 ranges','measurements':'exclude flagged data_validity_comment and potential_duplicate records','exclusions':'all recorded molecules and exact Murcko scaffolds from previous 140-compound benchmark, previous 12 candidates, all 1961 valid structures in the recent 2000-sample run, and both ALK5 pose-control ligands','leakage_identity':'normalized parent, neutralization where possible, canonical tautomer, stereochemistry ignored for identity/scaffold exclusions','matching':'greedy nearest property match across mw/logP/HBD/HBA/rot/TPSA; no score-based selection','target_pairs':50,'threshold_relaxation':False,'primary_scoring':'Vina and DOCK6 separately; fixed unweighted mean rank secondary','labels':'retained in separate key; not supplied to docking jobs','success_interpretation':'ROC-AUC and 95% CI; top-10 active count against exact balanced null; 8/10 active is a descriptive target, not a guarantee or universal validation criterion','prior_training_overlap':'unknown; only overlap with recorded project data can be excluded'}
(p/'eligibility_protocol.json').write_text(json.dumps(protocol,indent=2))
base='https://www.ebi.ac.uk/chembl/api/data/'
def get(url,dest):
 if dest.exists():return json.loads(dest.read_text())
 last=None
 for attempt in range(3):
  try:
   request=urllib.request.Request(url,headers={'User-Agent':'AIPharma-retrospective-validation/1.0','Accept':'application/json'})
   with urllib.request.urlopen(request,timeout=30) as response:data=response.read()
   body=json.loads(data);dest.write_bytes(data);return body
  except Exception as exc:last=exc;print('Fetch retry: '+str(exc),flush=True)
 raise RuntimeError(str(last))
target=get(base+'target/CHEMBL4439.json',raw/'target.json');assert target['organism']=='Homo sapiens'
params={'target_chembl_id':'CHEMBL4439','standard_type':'IC50','pchembl_value__isnull':'false','limit':1000,'format':'json'}
url=base+'activity.json?'+urllib.parse.urlencode(params);page=0;rows=[]
while url:
 body=get(url,raw/f'activity_page_{page:03d}.json');rows.extend(body['activities']);print(f'page={page} cumulative={len(rows)} total={body["page_meta"]["total_count"]}',flush=True)
 nxt=body['page_meta'].get('next');url=urllib.parse.urljoin(base,nxt) if nxt else None;page+=1
(raw/'all_activities.json').write_text(json.dumps({'activities':rows},indent=2));print('Fetched '+str(len(rows))+' IC50 measurements',flush=True)
