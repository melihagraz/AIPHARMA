from pathlib import Path
import concurrent.futures
import json
import requests

p=Path.home()/'aipharma_external_leads_20261018'
ds=sorted({r['doi'] for r in json.loads((p/'weak_lead_screen.json').read_text()) if not r['parent_in_prior_raw']})
d=p/'datacite';d.mkdir(exist_ok=True)
def fetch(doi):
 try:
  response=requests.get('https://api.datacite.org/dois/'+doi,timeout=30);response.raise_for_status()
  (d/(doi.split('/')[-1]+'.json')).write_bytes(response.content)
  a=response.json()['data']['attributes']
  return {'doi':doi,'titles':a.get('titles'),'url':a.get('url'),'related':a.get('relatedIdentifiers'),'description':a.get('descriptions')}
 except Exception as e:return {'doi':doi,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:results=list(pool.map(fetch,ds))
(p/'source_resolution.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2),flush=True)
