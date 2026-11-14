from pathlib import Path
import json,requests,concurrent.futures
p=Path.home()/'aipharma_analogues_20261016';o=p/'safety_audit';o.mkdir(exist_ok=True)
parents=json.loads((p/'parents.json').read_text());docs=sorted({d for r in parents for d in r['documents']})
def fetch(t):
 kind,id=t;url=f'https://www.ebi.ac.uk/chembl/api/data/{kind}/{id}.json';r=requests.get(url,timeout=90);r.raise_for_status();data=r.json();(o/f'{id}.json').write_text(json.dumps(data,indent=2));return kind,id,data
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 data=list(ex.map(fetch,[('molecule',r['parent_id']) for r in parents]+[('document',d) for d in docs]))
summary=[]
for kind,id,d in data:
 if kind=='molecule':summary.append({'id':id,'name':d.get('pref_name'),'synonyms':[s['molecule_synonym'] for s in d.get('molecule_synonyms',[])],'max_phase':d.get('max_phase'),'inchi_key':d.get('molecule_structures',{}).get('standard_inchi_key')})
 else:print(json.dumps({'document':id,'title':d.get('title'),'doi':d.get('doi'),'pubmed_id':d.get('pubmed_id'),'year':d.get('year')}),flush=True)
(o/'identity_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
