"""Inventory all exact IC50 source spectra, not only extreme negatives."""
import json,re,sys,statistics
from pathlib import Path
from collections import defaultdict
root=Path(sys.argv[1]); out=root/'reports/alk5_data_expansion_20261019'; out.mkdir(exist_ok=True)
data=json.loads((root/'reports/alk5_external_leads_20261018/bindingdb_high_cutoff.json').read_text(encoding='utf-8-sig'))
rows=next(iter(data.values()))['affinities']
docs=json.loads((root/'reports/alk5_external_leads_20261018/known_document_registry.json').read_text(encoding='utf-8-sig'))
dois={str(r['doi']).lower() for r in docs if r.get('doi')}; pmids={str(r['pubmed_id']) for r in docs if r.get('pubmed_id')}
sources=defaultdict(list)
for r in rows:
    if r['affinity_type']!='IC50': continue
    m=re.fullmatch(r'\s*([<>=~]*)\s*([\d.]+(?:[eE][+-]?\d+)?)\s*',r['affinity'])
    if not m: continue
    relation=m[1] or '='; value=float(m[2]); key=(r.get('doi','').lower(),str(r.get('pmid','')))
    sources[key].append(dict(**r,relation=relation,value_nM=value))
summary=[]
for (doi,pmid),rs in sources.items():
    exact=[r for r in rs if r['relation']=='=' and r['value_nM']>0]
    summary.append(dict(doi=doi,pmid=pmid,rows=len(rs),exact=len(exact),
        exact_strong=sum(r['value_nM']<=100 for r in exact),
        exact_intermediate=sum(100<r['value_nM']<10000 for r in exact),
        exact_weak=sum(r['value_nM']>=10000 for r in exact),
        bound_weak=sum(r['relation'] in ['>','>='] and r['value_nM']>=10000 for r in rs),
        min_exact=min((r['value_nM'] for r in exact),default=None),
        max_exact=max((r['value_nM'] for r in exact),default=None),
        source_matches_prior=doi in dois or pmid in pmids,
        needs_original_source_resolution=doi.startswith('10.7270/')))
summary.sort(key=lambda r:(-r['exact_weak'],-r['exact_intermediate'],-r['exact']))
(out/'source_inventory.json').write_text(json.dumps(summary,indent=2))
(out/'source_rows.json').write_text(json.dumps([dict(doi=doi,pmid=pmid,records=rs) for (doi,pmid),rs in sources.items()],indent=2))
print(json.dumps(dict(total_measurements=len(rows),IC50_sources=len(summary),new_resolved_sources=[r for r in summary if not r['source_matches_prior'] and not r['needs_original_source_resolution']],most_mixed_sources=summary[:15]),indent=2))
