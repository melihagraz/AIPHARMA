from pathlib import Path
import json,re
p=Path.home()/'aipharma_external_leads_20261018'
sections=json.loads((p/'mapping_synthesis_sections.json').read_text())
for number,text in sections.items():
 # Read reaction-product contexts, retaining source wording for explicit review.
 text=text.split('1 H NMR')[0]
 hits=list(re.finditer(r'(?i)obtain(?:ed)?\s+|give\s+|afford(?:ed)?\s+',text))
 print('COMPOUND',number)
 for hit in hits[-2:]:print(text[hit.end():hit.end()+420])
 if not hits:print(text[:650])
s=(p/'sources/US9682983B2.txt').read_text()
print('ALL48',[(m.start(),s[m.start():m.start()+65]) for m in re.finditer(r'Compound 48',s)])
