"""Match patent-specific compound labels to specified ALK5 numeric table cells."""
import json,re,sys
from decimal import Decimal, ROUND_HALF_EVEN
from pathlib import Path
from collections import Counter
p=Path(sys.argv[1]); out=p/'additional_local_series/primary_audit'
configs={10830:('US11447490B2',2),10870:('US11466003B2',2),11335:('US11691981B2',4)}
tables={}; labels={}
for eid,(patent,index) in configs.items():
    data=json.loads((out/(patent+'_tables.json')).read_text()); table={}
    for ri,row in enumerate(data['all_tables'][index]):
        cells=[(ci,c.strip()) for ci,c in enumerate(row) if c.strip()]
        if len(cells)==1 and cells[0][1].startswith('*Compound A:'): continue  # Explicit reference-control footnote, not a numeric row.
        assert len(cells)%2==0,(eid,ri,cells)
        for j in range(0,len(cells),2):
            key=cells[j][1]; value=cells[j+1][1]
            assert key not in table,(eid,key)
            table[key]=dict(reported_value=value,row_zero_based=ri,value_column_zero_based=cells[j+1][0],table_index=index)
    tables[eid]=table
    for r in data['compound_labels']:
        numbers=re.findall(r'Compound\s+([0-9]+[a-z]?)\b',r['label'])
        if numbers: labels[eid,r['monomer_id']]=dict(compound_numbers=sorted(set(numbers)),bindingdb_label=r['label'],bindingdb_page=r['source'])
rows=json.loads((p/'additional_local_series/curated_complete_identity/records.json').read_text())
queries=json.loads((p/'additional_local_series/within_series/predictions.json').read_text())
keys={(r['assay_id'],r['normalized_parent']) for r in queries}; results=[]
for r in rows:
    if (r['assay_id'],r['normalized_parent']) not in keys: continue
    checks=[]
    for raw in r['source_records']:
        eid=raw['entry_id']; label=labels.get((eid,raw['monomer_id'])); check=dict(monomer_id=raw['monomer_id'],entry=eid,bindingdb_value=raw['reported_value'])
        if not label: check['status']='missing_compound_label'
        else:
            check.update(label); check['patent']=configs[eid][0]; candidates=[]
            for number in label['compound_numbers']:
                cell=tables[eid].get(number)
                if cell: candidates.append(dict(compound_number=number,**cell))
            check['primary_candidates']=candidates; matches=[]
            for cell in candidates:
                try:
                    b=Decimal(raw['reported_value'].removesuffix('nM')); a=Decimal(cell['reported_value'])
                    if a==b: matches.append(dict(**cell,match_type='exact'))
                    elif a.quantize(Decimal(1).scaleb(b.as_tuple().exponent),rounding=ROUND_HALF_EVEN)==b:
                        matches.append(dict(**cell,match_type='database_reported_precision_rounding'))
                except Exception: pass
            check['matching_candidates']=matches
            check['status']=('numeric_value_matches_nM' if any(c['match_type']=='exact' for c in matches) else 'primary_rounds_to_database_precision') if matches else 'unresolved_primary_mapping'
        checks.append(check)
    results.append(dict(assay_id=r['assay_id'],normalized_parent=r['normalized_parent'],checks=checks,
        all_alias_values_match=all(c['status'] in ['numeric_value_matches_nM','primary_rounds_to_database_precision'] for c in checks),
        drawn_structure_independently_verified=False))
summary=dict(tested_molecules=len(queries),matched_molecules=sum(r['all_alias_values_match'] for r in results),
    alias_status_counts=dict(Counter(c['status'] for r in results for c in r['checks'])),
    primary_value_units='nM, from explicitly identified ALK5/TGFBR1 table headers',
    structure_mapping='BindingDB source-specific compound-number label connects SMILES to patent row; original patent drawing identity not independently verified',
    full_structure_validation=False,model_or_thresholds_changed=False)
(out/'query_value_matches.json').write_text(json.dumps(results,indent=2)); (out/'match_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
