import json,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as E
p=Path(sys.argv[1]); ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(p/'41587_2026_3090_MOESM4_ESM.xlsx') as z:
    styles=E.fromstring(z.read('xl/styles.xml'))
    formats={int(f.get('numFmtId')):f.get('formatCode') for f in styles.findall('s:numFmts/s:numFmt',ns)}
    formats.update({0:'General',9:'0%',10:'0.00%'})
    xfs=styles.findall('s:cellXfs/s:xf',ns)
    checks=[]
    for name,column in [('xl/worksheets/sheet3.xml','K'),('xl/worksheets/sheet4.xml','L')]:
        for c in E.fromstring(z.read(name)).findall('.//s:sheetData/s:row/s:c',ns):
            address=c.get('r'); value=c.find('s:v',ns)
            if ''.join(x for x in address if x.isalpha())!=column or int(''.join(x for x in address if x.isdigit()))<=12 or value is None: continue
            style=int(c.get('s','0')); fmt=int(xfs[style].get('numFmtId','0'))
            checks.append(dict(sheet_xml=name,cell=address,raw_value=value.text,style=style,number_format_id=fmt,number_format=formats.get(fmt)))
            break
(p/'source_percent_style_check.json').write_text(json.dumps(checks,indent=2)); print(json.dumps(checks,indent=2))
