from pathlib import Path
import json,urllib.request,urllib.parse
p=Path.home()/'aipharma_holdout_20261016';raw=p/'raw'
protocol_path=p/'eligibility_protocol.json';protocol=json.loads(protocol_path.read_text());protocol.update({'assay_type':'B (binding/biochemical only)','assay_confidence_score':9,'assay_target_required':'CHEMBL4439','within_class_scaffold_cap':1,'cross_class_exact_scaffold_overlap':'exclude weak controls if scaffold occurs in eligible active pool','rule_freeze_note':'Assay-quality and within-set scaffold rules fixed before eligibility computation or any scoring'})
protocol_path.write_text(json.dumps(protocol,indent=2))
body=json.loads((raw/'all_activities.json').read_text());ids=sorted({r['assay_chembl_id'] for r in body['activities']});assays=[]
for start in range(0,len(ids),50):
 dest=raw/f'assays_{start//50:03d}.json'
 if dest.exists():b=json.loads(dest.read_text())
 else:
  url='https://www.ebi.ac.uk/chembl/api/data/assay.json?'+urllib.parse.urlencode({'assay_chembl_id__in':','.join(ids[start:start+50]),'limit':1000})
  last=None
  for attempt in range(3):
   try:
    with urllib.request.urlopen(url,timeout=30) as response:data=response.read()
    b=json.loads(data);dest.write_bytes(data);break
   except Exception as exc:last=exc
  else:raise RuntimeError(str(last))
 assays.extend(b['assays']);print(f'Assays {len(assays)}/{len(ids)}',flush=True)
assert set(ids)=={a['assay_chembl_id'] for a in assays}
(raw/'all_assays.json').write_text(json.dumps({'assays':assays},indent=2))
print('Assay metadata complete',flush=True)
