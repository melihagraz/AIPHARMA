"""Archive patent sources and consolidate new-parent leads without claiming validation."""
from pathlib import Path
from collections import defaultdict
import concurrent.futures
import hashlib
import json
import re
import requests
from html.parser import HTMLParser

class PlainText(HTMLParser):
 def __init__(self):super().__init__();self.parts=[];self.skip=0
 def handle_starttag(self,tag,attrs):
  if tag in ['script','style']:self.skip+=1
 def handle_endtag(self,tag):
  if tag in ['script','style']:self.skip=max(0,self.skip-1)
 def handle_data(self,data):
  if not self.skip and data.strip():self.parts.append(data.strip())

p=Path.home()/'aipharma_external_leads_20261018';dest=p/'sources';dest.mkdir(exist_ok=True)
patents=['US9682983B2','US10017516B2','US20240018122A1','US20240116947A1','US20240116948A1']
def fetch(ident):
 path=dest/(ident+'.html');url=f'https://patents.google.com/patent/{ident}/en'
 try:
  if not path.exists():
   r=requests.get(url,timeout=40);r.raise_for_status();path.write_bytes(r.content)
  text=path.read_text()
  # Tables remain in original archived HTML; this is only a searchable text helper.
  parser=PlainText();parser.feed(text);plain=' '.join(parser.parts)
  (dest/(ident+'.txt')).write_text(plain,encoding='utf-8')
  families=re.findall(r'Family ID[=: ]+\s*(\d+)',plain)
  return {'patent':ident,'url':url,'status':'complete','family_ids':sorted(set(families)),
          'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
 except Exception as e:return {'patent':ident,'url':url,'status':'failed','error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:archives=list(pool.map(fetch,patents))
(p/'patent_archive_log.json').write_text(json.dumps(archives,indent=2));print(json.dumps(archives,indent=2),flush=True)
rows=json.loads((p/'weak_lead_screen.json').read_text())
new=[r for r in rows if not r['parent_in_prior_raw']]
mapping={'10.7270/Q28S4N3R':'US9682983B2','10.7270/Q20P12CN':'US10017516B2',
         '10.7270/Q2ZC873C':'US20240018122A1','10.7270/Q2BG2T8T':'US20240116948A1','10.7270/Q2G73K10':'US20240116947A1'}
groups=defaultdict(list)
for r in new:groups[r['normalized_parent']].append(r)
consolidated=[]
for parent,records in sorted(groups.items()):
 sources=sorted({mapping[r['doi']] for r in records})
 consolidated.append({'normalized_parent':parent,'original_smiles':records[0]['smile'],
                      'bindingdb_monomers':sorted({r['monomerid'] for r in records}),
                      'patent_sources':sources,'dataset_dois':sorted({r['doi'] for r in records}),
                      'weak_values_nM':sorted({r['value_nM'] for r in records}),
                      'relations':sorted({r['relation'] for r in records}),
                      'scaffold':records[0]['scaffold'],'scaffold_in_prior_raw':records[0]['scaffold_in_prior_raw'],
                      'strong_IC50_conflict_flag':any(r['same_monomer_has_strong_IC50_in_external_query'] for r in records),
                      'source_endpoint_review':'Catalytic enzyme methods found in patents; individual molecule-table mapping pending',
                      'accepted_for_benchmark':False})
(p/'new_parent_leads_consolidated.json').write_text(json.dumps(consolidated,indent=2))
print('consolidated_new_parents',len(consolidated),flush=True)
