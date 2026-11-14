"""Extract complete ALK5 columns, preserving single-dose and dose-series endpoints separately."""
import hashlib,json,math,re,sys,zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as E
p=Path(sys.argv[1]); file=p/'41587_2026_3090_MOESM4_ESM.xlsx'; before=hashlib.sha256(file.read_bytes()).hexdigest()
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(file) as z:
    shared=[''.join(si.itertext()) for si in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
    def sheet(name):
        output={}
        for row in E.fromstring(z.read(name)).findall('.//s:sheetData/s:row',ns):
            cells={}
            for c in row.findall('s:c',ns):
                value=c.find('s:v',ns); text=value.text if value is not None else ''
                if c.get('t')=='s' and text: text=shared[int(text)]
                elif c.get('t')=='inlineStr': text=''.join(c.find('s:is',ns).itertext())
                col=''.join(x for x in c.get('r') if x.isalpha()); cells[col]=text
            output[int(row.get('r'))]=cells
        return output
    s4=sheet('xl/worksheets/sheet4.xml'); s3=sheet('xl/worksheets/sheet3.xml'); s2=sheet('xl/worksheets/sheet2.xml')
    assert s4[9]['A']=='Compound' and s4[9]['L']=='ALK5_TGFBR1'
    assert s3[12]['A']=='Compound' and s3[12]['B']=='Name' and s3[12]['C']=='Dose' and s3[12]['K']=='TGFBR1'
    single=[]; missing=[]
    for number,r in sorted(s4.items()):
        if number<=9 or not r.get('A','').strip(): continue
        raw=r.get('L','')
        try: value=float(raw); assert math.isfinite(value)
        except Exception:
            missing.append(dict(compound=r['A'],source_cell=f'L{number}',source_value=raw)); continue
        inhibition=100-value
        label='screen_inhibitor' if inhibition>=70 else ('low_screen_inhibition' if inhibition<=20 else 'intermediate')
        single.append(dict(compound=r['A'],dose_uM=1.,residual_activity_percent=value,inhibition_percent=inhibition,label=label,source_sheet='Table S4',source_cell=f'L{number}',source_name_cell=f'A{number}',outside_nominal_percent_range=not 0<=value<=100,IC50_nM=None))
    assert len({r['compound'] for r in single})==len(single)
    series=[]; series_missing=[]
    for number,r in sorted(s3.items()):
        if number<=12 or not r.get('A','').strip(): continue
        try:
            value=float(r['K'])
            unit_match=re.fullmatch(r'\s*([\d.Ee+\-]+)\s*[uµμ]M\s*',r['C'])
            assert unit_match, 'Dose unit not explicitly uM'
            dose=float(unit_match[1]); assert math.isfinite(value) and math.isfinite(dose) and dose>0
        except Exception:
            series_missing.append(dict(source_row=number,source_values=r)); continue
        series.append(dict(compound_id=r['A'],compound=r.get('B'),dose_source_text=r['C'],dose_uM=dose,response_raw=value,response_scale='Percent inhibition, as explicitly defined in Table S3 cell A9; numeric values retained without rescaling',source_sheet='Table S3',source_cell=f'K{number}',dose_cell=f'C{number}'))
    metadata=[dict(source_row=i,cells=r) for i,r in s2.items() if r.get('B')=='ALK5/TGFBR1' or r.get('C')=='TGFBR1']
assert hashlib.sha256(file.read_bytes()).hexdigest()==before
for name,data in [('reserved_alk5_single_dose.json',single),('reserved_alk5_dose_series.json',series),('reserved_missing_values.json',dict(single=missing,series=series_missing)),('official_ALK5_metadata.json',metadata)]: (p/name).write_text(json.dumps(data,indent=2,allow_nan=False))
summary=dict(source_doi='10.1038/s41587-026-03090-8',source_sha256=before,source_modified=False,whole_source_reserved=True,modeled=False,
    single_dose_records=len(single),single_dose_missing=len(missing),single_dose_labels=dict(Counter(r['label'] for r in single)),single_dose_percent_noise_flags=sum(r['outside_nominal_percent_range'] for r in single),
    dose_series_points=len(series),dose_series_compounds=len({r['compound_id'] for r in series}),dose_series_doses_uM=sorted({r['dose_uM'] for r in series}),dose_series_IC50_fitted=False,
    dose_series_raw_response_range=[min((r['response_raw'] for r in series),default=None),max((r['response_raw'] for r in series),default=None)],dose_series_response_unit='Percent inhibition per Table S3 A9; no multiplication by 100',dose_series_fit_blocker='No IC50 fitting performed: first audit dose coverage and curve shape; observed threshold crossings are not fitted IC50 measurements',
    ROC_AUC_permitted_by_predeclared_minimum=all(sum(r['label']==label for r in single)>=5 for label in ['screen_inhibitor','low_screen_inhibition']),
    primary_source_consistency='Article mentions 61 dose-response inhibitors, supplied TableS3 describes and supplies only the observed subset; use actual extracted counts, do not invent missing compounds',
    credit='Saifudeen et al. 2026; Reaction Biology Corporation biochemical data acknowledged')
(p/'extraction_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
