import json,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as E
p=Path(sys.argv[1]); ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(p/'41587_2026_3090_MOESM4_ESM.xlsx') as z:
    ss=[''.join(r.itertext()) for r in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    result=[]
    for name,row_ids in [('xl/worksheets/sheet3.xml',{5,6,7,8,9,10,11,12}),('xl/worksheets/sheet4.xml',{5,6,7,8,9})]:
        output=[]
        for r in E.fromstring(z.read(name)).findall('.//s:sheetData/s:row',ns):
            if int(r.get('r')) not in row_ids: continue
            cells={}
            for c in r.findall('s:c',ns):
                address=c.get('r'); col=''.join(x for x in address if x.isalpha())
                if len(col)>1 or col>'M': continue
                v=c.find('s:v',ns); text=v.text if v is not None else ''
                if c.get('t')=='s' and text: text=ss[int(text)]
                cells[address]=text
            output.append(cells)
        result.append(dict(sheet_xml=name,header_cells=output))
(p/'source_layout.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
