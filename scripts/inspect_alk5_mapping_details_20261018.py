from pathlib import Path
from html.parser import HTMLParser
import json,re
class Text(HTMLParser):
 def __init__(self):super().__init__();self.parts=[]
 def handle_data(self,s):
  if s.strip():self.parts.append(s.strip())
p=Path.home()/'aipharma_external_leads_20261018'
s=(p/'sources/US9682983B2.txt').read_text()
for needle in ['Compound 58','Compound 67','Compound 48′','Example 1 Synthetic Protocols']:
 i=s.find(needle);print(needle,s[i:i+1300])
parser=Text();parser.feed((p/'monomer_118138.html').read_text());plain=' '.join(parser.parts)
for needle in ['US968','US100','Compound','Name','118138']:
 for m in list(re.finditer(re.escape(needle),plain))[:3]:print('BDB',needle,plain[max(0,m.start()-80):m.start()+280])
data=json.loads((p/'patent_table_review.json').read_text())
print('108945',[r['original_smiles'] for r in data if '108945' in r['bindingdb_monomers']])
