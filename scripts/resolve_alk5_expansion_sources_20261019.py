"""Resolve all mixed and broad-range dataset sources before assay curation."""
import concurrent.futures,hashlib,json,sys,urllib.request
from pathlib import Path
p=Path(sys.argv[1]); dest=p/'datacite'; dest.mkdir(exist_ok=True)
sources=json.loads((p/'source_inventory.json').read_text())
selected=[r for r in sources if r['needs_original_source_resolution'] and ((r['exact_strong'] and r['exact_intermediate']+r['exact_weak']) or r['exact_weak'])]
def fetch(r):
    doi=r['doi']; target=dest/(doi.split('/')[-1]+'.json')
    try:
        if not target.exists():
            req=urllib.request.Request('https://api.datacite.org/dois/'+doi,headers={'User-Agent':'AIPHARMA-research-source-audit/1.0'})
            with urllib.request.urlopen(req,timeout=40) as response: content=response.read()
            json.loads(content); target.write_bytes(content)
        raw=target.read_bytes(); d=json.loads(raw)['data']['attributes']
        return dict(doi=doi,inventory=r,titles=d.get('titles'),url=d.get('url'),related=d.get('relatedIdentifiers'),sha256=hashlib.sha256(raw).hexdigest())
    except Exception as e: return dict(doi=doi,error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(fetch,selected))
(p/'source_resolution.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2),flush=True)
