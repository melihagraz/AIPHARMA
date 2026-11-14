"""Reserve an entire eligible, previously unfitted publication/scaffold component."""
import hashlib,json,sys
from pathlib import Path
from collections import defaultdict,Counter
p=Path(sys.argv[1]); source=p/'curated_complete_identity'; out=p/'locked_split'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(component_edges=['normalized identity','normalized exact Murcko scaffold','canonical publication DOI/document'],
    reserved_eligibility='>=10 exact assay/structure groups, contains newly reviewed catalytic endpoint, no identity from previous 541-query-model training data, all sources are verified publication documents',
    choice='First SHA256-sorted eligible component by canonical source keys; no outcomes, scores or class balance used',
    no_eligibility_fallback='If none qualify, stop before modeling and require another source; never relabel a previously fitted component as untouched',
    scope='Held out from upcoming model development and from previous 541-record regressors by identity; source was present in raw public archive and inspected, not blind or prospectively new measurements'))
rows=json.loads((source/'records.json').read_text()); docs=json.loads((p/'prior_document_registry.json').read_text())
docs={r['document_chembl_id']:r for r in docs}
prior=json.loads((p/'prior_measured_records.json').read_text()); known={r['normalized_parent'] for r in prior}
parents=list(range(len(rows)))
def find(i):
    while parents[i]!=i: parents[i]=parents[parents[i]]; i=parents[i]
    return i
def union(i,j): parents[find(j)]=find(i)
def sourcekey(d):
    item=docs.get(d,{})
    return 'doi:'+str(item['doi']).lower() if item.get('doi') else 'document:'+d
seen={}
for i,r in enumerate(rows):
    keys=[('identity',r['normalized_parent']),('scaffold',r['scaffold'])]+[('source',sourcekey(d)) for d in r['documents']]
    for k in keys:
        if k in seen: union(i,seen[k])
        else: seen[k]=i
groups=defaultdict(list)
for i in range(len(rows)): groups[find(i)].append(i)
components=[]
for indexes in groups.values():
    rs=[rows[i] for i in indexes]; documents=sorted({d for r in rs for d in r['documents']}); keys=sorted({sourcekey(d) for d in documents})
    complete_sources=all(d in docs and docs[d].get('doc_type')=='PUBLICATION' and docs[d].get('doi') for d in documents)
    overlap=sum(r['normalized_parent'] in known for r in rs)
    c=dict(key_sha256=hashlib.sha256('|'.join(keys).encode()).hexdigest(),indices=indexes,sources=keys,documents=documents,groups=len(rs),exact_groups=sum(r['median_pic50'] is not None for r in rs),prior_model_identity_overlap=overlap,publication_identity_complete=bool(complete_sources),has_new_endpoint=any(r['newly_reviewed_endpoint'] for r in rs))
    c['eligible_reserved']=bool(c['exact_groups']>=10 and c['has_new_endpoint'] and overlap==0 and complete_sources)
    components.append(c)
candidates=sorted([c for c in components if c['eligible_reserved']],key=lambda c:c['key_sha256'])
save('components.json',components)
if not candidates:
    save('summary.json',dict(status='no_eligible_unused_component',components=len(components),modeling_permitted=False)); print((out/'summary.json').read_text()); sys.exit()
reserved=set(candidates[0]['indices']); train=[r for i,r in enumerate(rows) if i not in reserved]; test=[r for i,r in enumerate(rows) if i in reserved]
for key in ['normalized_parent','scaffold']:
    assert {r[key] for r in train}.isdisjoint(r[key] for r in test)
assert {sourcekey(d) for r in train for d in r['documents']}.isdisjoint(sourcekey(d) for r in test for d in r['documents'])
save('training_records.json',train); save('reserved_records.json',test)
summary=dict(status='locked_before_modeling',components=len(components),eligible_unused_components=len(candidates),training_groups=len(train),training_exact=sum(r['median_pic50'] is not None for r in train),training_interval=sum(r['median_pic50'] is None for r in train),training_labels=dict(Counter(r['label'] for r in train)),reserved_groups=len(test),reserved_exact=sum(r['median_pic50'] is not None for r in test),reserved_sources=candidates[0]['sources'],reserved_labels_not_summarized=True,training_sha256=hashlib.sha256((out/'training_records.json').read_bytes()).hexdigest(),reserved_sha256=hashlib.sha256((out/'reserved_records.json').read_bytes()).hexdigest(),reserved_opened_for_model_selection=False,modeling_permitted=True)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
