"""Archive table cells and BindingDB compound labels for explicit numeric mapping."""
import hashlib,html,json,re,sys
from pathlib import Path
from html.parser import HTMLParser
p=Path(sys.argv[1]); out=p/'additional_local_series/primary_audit'
class Tables(HTMLParser):
    def __init__(self): super().__init__(); self.tables=[]; self.stack=[]; self.row=None; self.cell=None
    def handle_starttag(self,t,a):
        if t=='table': self.stack.append([])
        elif self.stack and t=='tr': self.row=[]
        elif self.stack and t in ('td','th'): self.cell=[]
    def handle_data(self,d):
        if self.cell is not None: self.cell.append(d)
    def handle_endtag(self,t):
        if t in ('td','th') and self.cell is not None:
            if self.row is not None: self.row.append(' '.join(''.join(self.cell).split()))
            self.cell=None
        elif t=='tr' and self.row is not None:
            if self.stack: self.stack[-1].append(self.row)
            self.row=None
        elif t=='table' and self.stack: self.tables.append(self.stack.pop())
result=[]
for eid,patent in [(10830,'US11447490B2'),(10870,'US11466003B2'),(11335,'US11691981B2')]:
    f=out/(patent+'.html'); parser=Tables(); parser.feed(f.read_text(encoding='utf-8-sig'))
    tables=[dict(index=i,rows=t) for i,t in enumerate(parser.tables) if re.search(r'IC\s*50|ALK5|TGF.?R1',str(t),re.I)]
    labels=[]
    for page in sorted((p/'bindingdb').glob(f'entry{eid}_page*.html')):
        for block in page.read_text(encoding='utf-8-sig').split('<div><span class="header">Ligand</span>')[1:]:
            mid=re.search(r'monomerid=(\d+)',block); label=re.search(r'>BDBM\d+</a>\((.*?)\)<button',block,re.S)
            if mid and label: labels.append(dict(monomer_id=mid[1],label=html.unescape(re.sub('<[^>]+>','',label[1])),source=page.name))
    data=dict(entry=eid,patent=patent,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),all_tables=parser.tables,candidate_tables=tables,compound_labels=labels)
    (out/(patent+'_tables.json')).write_text(json.dumps(data,indent=2))
    result.append(dict(entry=eid,tables=[dict(index=i,rows=len(t),preview=t[:3]) for i,t in enumerate(parser.tables[:8])],labels_preview=labels[:1]))
print(json.dumps(result,indent=2))
