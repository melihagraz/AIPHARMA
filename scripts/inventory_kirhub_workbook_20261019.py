"""Read-only workbook structure/target-cell inventory; no authoring or evaluation."""
import hashlib,json,posixpath,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
p=Path(sys.argv[1]); file=p/'41587_2026_3090_MOESM4_ESM.xlsx'
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
with zipfile.ZipFile(file) as z:
    shared=[]
    if 'xl/sharedStrings.xml' in z.namelist(): shared=[''.join(si.itertext()) for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    relationships={r.get('Id'):r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    result=[]
    for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',ns):
        target=relationships[sheet.get('{'+ns['r']+'}id')]
        name=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        xml=ET.fromstring(z.read(name)); values=[]; hits=[]; row_count=0; cells=0
        for row in xml.findall('.//s:sheetData/s:row',ns):
            row_count+=1; preview={}
            for c in row.findall('s:c',ns):
                cells+=1; value=c.find('s:v',ns); text=value.text if value is not None else ''
                if c.get('t')=='s' and text: text=shared[int(text)]
                elif c.get('t')=='inlineStr': text=''.join(c.find('s:is',ns).itertext())
                if row_count<=4 and len(preview)<10: preview[c.get('r')]=text
                if 'TGFBR1' in text.upper() or 'ALK5' in text.upper(): hits.append(dict(cell=c.get('r'),value=text))
            if preview: values.append(preview)
        result.append(dict(sheet=sheet.get('name'),sheet_xml=name,rows=row_count,cells=cells,header_preview=values,target_labels=hits))
report=dict(source_file_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),sheets=result,source_modified=False)
(p/'workbook_inventory.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))
