from pathlib import Path
from html.parser import HTMLParser
import requests,re
p=Path.home()/'aipharma_structure_mapping_20261018'
url='https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid=667169'
r=requests.get(url,timeout=30);r.raise_for_status();(p/'monomer_667169.html').write_bytes(r.content)
title=re.search(r'<title>(.*?)</title>',r.text,re.S).group(1);print('TITLE',title)
s=(Path.home()/'aipharma_external_leads_20261018/sources/US20240116947A1.txt').read_text()
for needle in ['pyridin-4-yl','pyridin-3-yl','pyridin-2-yl']:
 for m in list(re.finditer(re.escape(needle),s))[-4:]:print('SOURCE',s[max(0,m.start()-170):m.start()+100])
