"""Record explicit review scope and quarantine two newly identified isotope conflicts."""
import hashlib,json,sys
from pathlib import Path
from collections import Counter
p=Path(sys.argv[1]); out=p/'drawing_audit_current'; out.mkdir(exist_ok=False)
evidence=p/'primary_audit/structure_evidence'; remaining=evidence/'remaining_drawing_review'
manifest=json.loads((remaining/'manifest.json').read_text()); rows=json.loads((p/'mass_conflict_excluded/curated_complete_identity/records.json').read_text())
graph=json.loads((evidence/'drawing_review_11335/graph_crosschecks.json').read_text())
text_only={'1','6','7','8','10','20'}; review=[]; excluded=[]
for r in manifest:
    number=r['compound_numbers'][0]; status='visual_connectivity_consistent_not_atom_mapped'
    notes=['Original route and archived depiction visually compared; no full manual graph transcription', 'Alkene geometry is not validated'] if r['assay_id']=='BindingDB10830_assay1' else ['Original route and archived depiction visually compared; not a complete atom-mapped reconstruction']
    if r['assay_id']=='BindingDB10870_assay1' and number in {'50','51'}:
        status='confirmed_isotope_conflict_quarantined'
        notes=['Original shared route C00134 and synthesis text explicitly use deuterated methyl iodide',
               'Compound50: NH-CD3 versus archived NH-CH3; compound51: N(CD3)2 versus archived N(CH3)2',
               'Source MH masses364.2 and381.0; original 44-query evaluation is superseded']
        excluded.append((r['assay_id'],r['normalized_parent']))
    elif r['assay_id']=='BindingDB10870_assay1' and number in text_only:
        status='synthesis_text_consistent_drawing_not_located'
        notes=['Exact synthesis paragraph inspected; no separate final-product image located, not full graphical validation']
    if r['tag']=='BindingDB10830_assay1_22':
        notes+=['Initial panel mistakenly selected following comparative23 route; corrected source is C00085, not C00086; full corrected image reviewed']
    if r['tag']=='BindingDB10870_assay1_37':
        notes+=['Initial panel showed only intermediate37-a; final product verified against continuation C00121']
    if r['tag']=='BindingDB10870_assay1_12':
        notes+=['Patent hydroxamate and archived imidic-acid depiction are tautomer-related; do not require identical displayed bond/proton placement']
    review.append(dict(tag=r['tag'],assay_id=r['assay_id'],normalized_parent=r['normalized_parent'],status=status,notes=notes))
for r in graph:
    assert r['graphs_match']
    review.append(dict(tag='BindingDB11335_assay1_'+r['compound_number'],assay_id='BindingDB11335_assay1',normalized_parent=r['canonical_archive'],
        status='manual_transcription_canonical_graph_match',notes=['Same analyst saw references: not blind independent expert review']))
clean=[r for r in rows if (r['assay_id'],r['normalized_parent']) not in set(excluded)]
assert len(rows)==44 and len(clean)==42 and len(review)==44
(out/'curated_complete_identity').mkdir()
for name,data in [('curated_complete_identity/records.json',clean),('review_records.json',review),('quarantined_keys.json',excluded)]:
    (out/name).write_text(json.dumps(data,indent=2,allow_nan=False))
summary=dict(reviewed_records=len(review),status_counts=dict(Counter(r['status'] for r in review)),
 remaining_test_records=len(clean),confirmed_additional_isotope_conflicts=len(excluded),
 all_remaining_structures_fully_verified=False,biological_activity_or_safety_validated=False,
 model_parameters_changed=False,previous_44_query_run_superseded=True,
 interpretation='Complete triage of remaining entries, but six text-only and26 visual-only entries still lack manual graph transcription',
 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'review_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
