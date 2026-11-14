"""Resolve public drug names against PubChem; archive responses and unresolved names."""
import concurrent.futures, hashlib, json, sys, time, urllib.request, urllib.parse
from pathlib import Path
p=Path(sys.argv[1]); out=p/'structure_resolution'; out.mkdir(exist_ok=True)
names=sorted({r['compound'] for r in json.loads((p/'reserved_alk5_single_dose.json').read_text())} |
             {r['compound'] for r in json.loads((p/'reserved_alk5_dose_series.json').read_text())})
def resolve(name):
    cache=out/(hashlib.sha256(name.encode()).hexdigest()[:16]+'.json')
    if cache.exists(): return json.loads(cache.read_text())
    url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(name,safe='')+'/property/IsomericSMILES,InChIKey/JSON'
    result=dict(name=name,query_name=name,url=url,source='PubChem PUG REST',status='unresolved')
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'AIPHARMA-public-data-curation/1.0'}),timeout=45) as response:
                raw=response.read().decode()
            result['raw_response']=json.loads(raw)
            props=result['raw_response'].get('PropertyTable',{}).get('Properties',[])
            result['status']='unique_name_match_pending_identity_audit' if len(props)==1 else 'ambiguous_or_empty'
            break
        except Exception as exc:
            result['error']=str(exc)
            if attempt<2: time.sleep(1+attempt)
    cache.write_text(json.dumps(result,indent=2)); return result
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    records=list(pool.map(resolve,names))
(out/'resolved_names.json').write_text(json.dumps(records,indent=2))
summary=dict(names=len(names),unique_name_matches=sum(r['status']=='unique_name_match_pending_identity_audit' for r in records),
 unresolved=[r['name'] for r in records if r['status']!='unique_name_match_pending_identity_audit'],
 structures_not_yet_chemically_audited=True,modeling_started=False)
(out/'summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
