"""Check numeric ALK5-column support, keeping molecule mapping explicitly unresolved."""
from pathlib import Path
from html.parser import HTMLParser
from collections import defaultdict
import json
import re

class Tables(HTMLParser):
 def __init__(self):super().__init__();self.tables=[];self.current=None;self.row=None;self.cell=None;self.depth=0
 def handle_starttag(self,tag,attrs):
  if tag=='table':
   if self.depth==0:self.current=[]
   self.depth+=1
  if self.depth and tag=='tr':self.row=[]
  if self.depth and tag in ['td','th']:self.cell=[]
 def handle_data(self,data):
  if self.cell is not None:self.cell.append(data)
 def handle_endtag(self,tag):
  if self.depth and tag in ['td','th'] and self.cell is not None:
   if self.row is not None:self.row.append(' '.join(' '.join(self.cell).split()))
   self.cell=None
  if self.depth and tag=='tr' and self.row is not None:self.current.append(self.row);self.row=None
  if tag=='table':
   self.depth-=1
   if self.depth==0:self.tables.append(self.current);self.current=None

p=Path.home()/'aipharma_external_leads_20261018';alltables={}
for f in sorted((p/'sources').glob('*.html')):
 parser=Tables();parser.feed(f.read_text());alltables[f.stem]=parser.tables
(p/'patent_tables.json').write_text(json.dumps(alltables,indent=2))
rows=json.loads((p/'new_parent_leads_consolidated.json').read_text());reviews=[]
for r in rows:
 matches=[]
 for patent in r['patent_sources']:
  tables=alltables[patent]
  for ti,header in enumerate(tables):
   column=None
   for line in header:
    for i,cell in enumerate(line):
     if re.fullmatch(r'ALK\s*5',cell) or ('ALK5' in cell and 'IC50' in cell):column=i
   if column is None:continue
   # Google splits the printed header and its body into adjacent HTML tables.
   table=header+(tables[ti+1] if ti+1<len(tables) else [])
   for line in table:
    for i,cell in enumerate(line):
     if re.fullmatch(r'ALK\s*5',cell):column=i
    if column is not None and len(line)>column:
     cell=line[column].replace(',','')
     try:value=float(cell)
     except ValueError:continue
     if value in r['weak_values_nM']:matches.append({'patent':patent,'table_index':ti,'row_identifier':line[0],'ALK5_value_nM':value,'type':'ALK5 column numeric match; structure mapping not established'})
    # Comparative tables have two columns with a compound identifier and ALK5 value.
   if any('ALK5' in ' '.join(line) and 'IC' in ' '.join(line) for line in header):
    for line in table:
     if len(line)==2 and line[0] in ['C1','C2','C 1','C 2']:
      try:value=float(line[1].replace(',',''))
      except ValueError:continue
      if value in r['weak_values_nM']:matches.append({'patent':patent,'table_index':ti,'row_identifier':line[0],'ALK5_value_nM':value,'type':'comparative ALK5 table numeric match; structure mapping not established'})
 reviews.append({**r,'numeric_table_support':matches,'numeric_table_support_found':bool(matches),'accepted_for_benchmark':False})
(p/'patent_table_review.json').write_text(json.dumps(reviews,indent=2))
print(json.dumps({'new_parents':len(reviews),'with_numeric_ALK5_table_support':sum(r['numeric_table_support_found'] for r in reviews),'without_numeric_support':[{k:r[k] for k in ['bindingdb_monomers','weak_values_nM','patent_sources']} for r in reviews if not r['numeric_table_support_found']]},indent=2))
