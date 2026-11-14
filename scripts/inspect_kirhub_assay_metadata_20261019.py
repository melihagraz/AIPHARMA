"""Read public assay metadata cells without executing workbook content."""
import json,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
p=Path(sys.argv[1]); ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(p/'Reaction_Biology_WT_Kinases_Info.xlsx') as z:
    shared=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        shared=[''.join(si.itertext()) for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    matched=[]
    for name in z.namelist():
        if not name.startswith('xl/worksheets/sheet') or not name.endswith('.xml'): continue
        rows=[]
        for row in ET.fromstring(z.read(name)).findall('.//s:sheetData/s:row',ns):
            cells={}
            for c in row.findall('s:c',ns):
                value=c.find('s:v',ns)
                text=value.text if value is not None else ''.join(c.find('s:is',ns).itertext()) if c.find('s:is',ns) is not None else ''
                if c.get('t')=='s' and text: text=shared[int(text)]
                cells[c.get('r')]=text
            rows.append(cells)
        hits=[r for r in rows if any('ALK5' in v.upper() or 'TGFBR1' in v.upper() for v in r.values())]
        if hits: matched.append(dict(sheet_xml=name,header_rows=rows[:3],target_rows=hits))
result=dict(source='Reaction Biology WT Kinases Info, official KIRHub download',matched_sheets=matched)
(p/'assay_metadata_ALK5.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
