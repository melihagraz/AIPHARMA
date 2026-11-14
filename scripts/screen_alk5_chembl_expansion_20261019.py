"""Screen additional record IDs with explicit exclusions; never infer potency from missing values."""
import json,sys,math
from pathlib import Path
from collections import Counter,defaultdict
root=Path(sys.argv[1]); p=root/'reports/alk5_data_expansion_20261019'
rows=json.loads((p/'chembl_raw/additional_activities.json').read_text())['activities']
assays={r['assay_chembl_id']:r for r in json.loads((p/'chembl_raw/all_assays.json').read_text())['assays']}
classes={r['assay_id']:r for r in json.loads((root/'reports/alk5_endpoint_audit_20261016/assay_classification.json').read_text(encoding='utf-8-sig'))}
accepted=[]; rejected=[]; review=[]
for r in rows:
    a=assays[r['assay_chembl_id']]; reason=[]
    if a.get('target_chembl_id')!='CHEMBL4439' or a.get('confidence_score')!=9 or a.get('assay_type')!='B': reason.append('assay_quality')
    if r.get('standard_relation')!='=' or r.get('standard_units')!='nM': reason.append('not_exact_nM')
    if r.get('data_validity_comment') or r.get('potential_duplicate'): reason.append('validity_or_duplicate')
    try:
        value=float(r['standard_value']); assert value>0 and math.isfinite(value)
    except Exception: reason.append('no_positive_numeric_value'); value=None
    c=classes.get(r['assay_chembl_id'])
    if c and (c['endpoint']!='enzyme_catalytic' or c['mutation']!='no_explicit_mutant'): reason.append('not_eligible_catalytic_endpoint')
    if not c and not reason: reason.append('unreviewed_endpoint')
    output=dict(activity_id=r['activity_id'],assay_id=r['assay_chembl_id'],document_id=r.get('document_chembl_id'),smiles=r.get('canonical_smiles'),value_nM=value,reported_pchembl=r.get('pchembl_value'),reasons=reason)
    if not reason: accepted.append(output)
    else: rejected.append(output)
for assayid in sorted({r['assay_chembl_id'] for r in rows}):
    if assayid not in classes:
        a=assays[assayid]; review.append({k:a.get(k) for k in ['assay_chembl_id','assay_type','confidence_score','description','bao_format','document_chembl_id','assay_organism','variant_sequence']})
summary=dict(additional_records=len(rows),provisionally_eligible_exact_records=len(accepted),reason_counts=dict(Counter(reason for r in rejected for reason in r['reasons'])),unreviewed_assays=review,accepted_details=accepted)
(p/'chembl_additional_screen.json').write_text(json.dumps(summary,indent=2)); (p/'chembl_additional_exclusions.json').write_text(json.dumps(rejected,indent=2))
print(json.dumps(summary,indent=2),flush=True)
