from pathlib import Path
import json,re
p=Path.home()/'aipharma_external_leads_20261018'
data=json.loads((p/'patent_table_review.json').read_text())
print([(r['bindingdb_monomers'],[(m['patent'],m['row_identifier']) for m in r['numeric_table_support']]) for r in data])
s=(p/'sources/US9682983B2.txt').read_text();s=s[s.index('Example 1 Synthetic Protocols'):]
sections={}
for n in [1,4,6,7,31,36,37,41,43,45,48,54,55,58,66,67]:
 m=re.search(r'\bCompound '+str(n)+r'\b',s)
 if not m:continue
 rest=s[m.end():];end=re.search(r'\bCompound \d+\b|\bExample 2\b',rest)
 section=rest[:end.start()] if end else rest[:6000]
 sections[str(n)]=section
 if n in [4,43,48]:print('COMPOUND',n,section[-1800:])
(p/'mapping_synthesis_sections.json').write_text(json.dumps(sections,indent=2))
for ident in ['US20240018122A1','US20240116947A1','US20240116948A1']:
 t=(p/'sources'/f'{ident}.txt').read_text();matches=list(re.finditer(r'Comparative Example|Example C|COMPARATIVE EXAMPLE',t))
 print('PATENT',ident,[(t[max(0,m.start()-50):m.start()+500]) for m in matches[:5]])
