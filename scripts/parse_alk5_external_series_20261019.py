"""Map archived official BindingDB HTML records to the primary patent IC50 table."""
import json,re,sys,html,hashlib
from pathlib import Path
from html.parser import HTMLParser
p=Path(sys.argv[1])
records=[]
for f in sorted(p.glob('bindingdb_entry12979*.html')):
    text=f.read_text(encoding='utf-8-sig')
    for block in text.split('<div><span class="header">Ligand</span>')[1:]:
        mid=re.search(r'monomerid=(\d+)',block)
        smi=re.search(r"setClipboard\('([^']+)'\).*?Copy&nbsp;SMILES",block)
        ic=re.search(r'IC50:&nbsp;([0-9.Ee+\-]+)nM',block)
        name=re.search(r'>BDBM\d+</a>\((.*?)\)<button',block,re.S)
        assert mid and smi and ic and name
        label=html.unescape(re.sub('<[^>]+>','',name[1]))
        num=re.search(r'^\s*(\d+)\s',label) or re.search(r'Compound\s+(\d+)',label)
        records.append(dict(monomer_id=mid[1],smiles=html.unescape(smi[1]),ic50_nM=float(ic[1]),compound_number=int(num[1]) if num else None,label=label,source_html=f.name))
assert len(records)==64, len(records)
assert len({r['monomer_id'] for r in records})==64
class Tables(HTMLParser):
    def __init__(self): super().__init__(); self.tables=[]; self.stack=[]; self.row=None; self.cell=None
    def handle_starttag(self,tag,attrs):
        if tag=='table': self.stack.append([])
        elif self.stack and tag=='tr': self.row=[]
        elif self.stack and tag in ('td','th'): self.cell=[]
    def handle_data(self,data):
        if self.cell is not None: self.cell.append(data)
    def handle_endtag(self,tag):
        if tag in ('td','th') and self.cell is not None:
            if self.row is not None: self.row.append(' '.join(''.join(self.cell).split()))
            self.cell=None
        elif tag=='tr' and self.row is not None:
            if self.stack: self.stack[-1].append(self.row)
            self.row=None
        elif tag=='table' and self.stack: self.tables.append(self.stack.pop())
t=Tables(); t.feed((p/'patent_US20250179094A1.html').read_text(encoding='utf-8-sig'))
candidates=[x for x in t.tables if any('LY364947' in ' '.join(row) for row in x) and any('20.9' in ' '.join(row) for row in x)]
assert len(candidates)==1, len(candidates)
(p/'primary_table_cells.json').write_text(json.dumps(candidates[0],indent=2))
clean_rows=[[cell for cell in row if cell.strip()] for row in candidates[0]]
table={row[0]:float(row[1]) for row in clean_rows if len(row)==2 and (row[0].isdigit() or row[0]=='LY364947')}
unmapped=[]
for r in records:
    key=str(r['compound_number']) if r['compound_number'] is not None else ('LY364947' if 'LY364947' in r['label'] else None)
    if key is None: unmapped.append(r); continue
    assert key in table and abs(table[key]-r['ic50_nM'])<1e-8,(key,r['ic50_nM'],table.get(key))
    r['primary_table_key']=key; r['primary_value_verified']=True
(p/'parsed_series.json').write_text(json.dumps(records,indent=2))
(p/'primary_table.json').write_text(json.dumps(table,indent=2))
summary=dict(records=len(records),primary_table_rows=len(table),verified=sum(r.get('primary_value_verified',False) for r in records),unmapped=unmapped,source='US20250179094A1',source_is_new_blind=False,all_values_active_range=True,structure_mapping='Official BindingDB SMILES; patent numeric table crosscheck does not independently verify drawn structure identity',hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in p.glob('*.html')})
(p/'parse_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
