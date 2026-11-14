"""Fetch public provenance-bearing ALK5 lead sources; no docking or label assignment."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import requests

out=Path.home()/'aipharma_external_leads_20261018';out.mkdir(exist_ok=True)
jobs={
 'bindingdb_high_cutoff.json':'https://bindingdb.org/rest/getLigandsByUniprots?uniprot=P36897&cutoff=100000000&response=application/json',
 'pubchem_1799383_description.json':'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1799383/description/JSON',
 'pubchem_1799383.csv':'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/1799383/CSV',
 'dorsomorphin_article_search.json':'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1021/cb9002865&format=json',
}
def fetch(job):
 name,url=job;dest=out/name
 if dest.exists():return {'name':name,'url':url,'status':'cached','sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}
 try:
  response=requests.get(url,timeout=40);response.raise_for_status()
  if name.endswith('.json'):response.json()
  dest.write_bytes(response.content)
  return {'name':name,'url':url,'status':'complete','bytes':len(response.content),'sha256':hashlib.sha256(response.content).hexdigest()}
 except Exception as e:return {'name':name,'url':url,'status':'failed','error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,jobs.items()))
(out/'fetch_log.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2),flush=True)
