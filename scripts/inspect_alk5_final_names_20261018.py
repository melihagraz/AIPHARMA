from pathlib import Path
import json,re
p=Path.home()/'aipharma_external_leads_20261018';s=(p/'sources/US9682983B2.txt').read_text()
i=s.index('Compound 67 and Compound 68');t=s[i:i+10000]
for m in re.finditer(r'(?:give|obtain|afford)',t):
 if 'quinoline' in t[m.end():m.end()+400]:print(t[m.end():m.end()+420])
data=json.loads((p/'patent_table_review.json').read_text())
print([(r['bindingdb_monomers'],r['original_smiles']) for r in data if any(x in r['bindingdb_monomers'] for x in ['161935','161923','127883','161922','161934'])])
s=(p/'sources/US20240018122A1.txt').read_text();i=s.rfind('Example C3');print('C3',s[i:i+1800])
for m in list(re.finditer(r'piperidin-4-yl', (p/'sources/US9682983B2.txt').read_text()))[:5]:
 s=(p/'sources/US9682983B2.txt').read_text();print('48context',s[max(0,m.start()-180):m.start()+220])
