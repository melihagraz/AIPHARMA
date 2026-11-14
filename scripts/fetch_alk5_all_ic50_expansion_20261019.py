"""Fetch all current IC50 records, without requiring precomputed pChEMBL."""
import concurrent.futures,hashlib,json,sys,urllib.request
from pathlib import Path
p=Path(sys.argv[1]); raw=p/'chembl_raw'; raw.mkdir(exist_ok=True)
def get(url,path):
    if not path.exists():
        req=urllib.request.Request(url,headers={'User-Agent':'AIPHARMA-public-data-audit/1.0'})
        with urllib.request.urlopen(req,timeout=60) as response: data=response.read()
        json.loads(data); path.write_bytes(data)
    return json.loads(path.read_bytes())
base='https://www.ebi.ac.uk/chembl/api/data'
protocol=dict(target='CHEMBL4439',standard_type='IC50',pchembl_required=False,stage='Raw collection only; endpoint, mutant, validity, stereo and source-family review required',holdout='No modeling before curation and source/scaffold partition locked')
(p/'chembl_fetch_protocol.json').write_text(json.dumps(protocol,indent=2))
status=get(base+'/status.json',raw/'status.json')
query=base+'/activity.json?target_chembl_id=CHEMBL4439&standard_type=IC50&limit=1000'
first=get(query,raw/'activities_00000.json'); total=first['page_meta']['total_count']
def page(offset): return get(query+f'&offset={offset}',raw/f'activities_{offset:05d}.json')['activities']
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    records=first['activities']+[r for rs in pool.map(page,range(1000,total,1000)) for r in rs]
assert len(records)==total and len({r['activity_id'] for r in records})==total
known=json.loads((Path(sys.argv[2])/'reports/alk5_holdout_20261016/raw/all_activities.json').read_text(encoding='utf-8-sig'))['activities']
known += json.loads((Path(sys.argv[2])/'reports/alk5_holdout_20261016/raw/weak_boundary_all.json').read_text(encoding='utf-8-sig'))['activities']
old={r['activity_id'] for r in known}; extra=[r for r in records if r['activity_id'] not in old]
(raw/'all_activities.json').write_text(json.dumps(dict(activities=records),indent=2))
(raw/'additional_activities.json').write_text(json.dumps(dict(activities=extra),indent=2))
assayids=sorted({r['assay_chembl_id'] for r in records})
def assaybatch(start): return get(base+'/assay.json?assay_chembl_id__in='+','.join(assayids[start:start+40])+'&limit=1000',raw/f'assays_{start:05d}.json')['assays']
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: assays=[r for rs in pool.map(assaybatch,range(0,len(assayids),40)) for r in rs]
assert set(assayids)=={r['assay_chembl_id'] for r in assays}
(raw/'all_assays.json').write_text(json.dumps(dict(assays=assays),indent=2))
summary=dict(chembl_status=status,raw_records=total,prior_record_ids=len(old),additional_record_ids=len(extra),all_assays=len(assayids),additional_assays=len({r['assay_chembl_id'] for r in extra}-{r['assay_chembl_id'] for r in known}),modeling_started=False,raw_sha256=hashlib.sha256((raw/'all_activities.json').read_bytes()).hexdigest())
(p/'chembl_fetch_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2),flush=True)
