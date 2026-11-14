#!/usr/bin/env python3
"""AIPHARMA Phase 6: exact document x assay provenance for Phase 4 ChEMBL candidates.

READ-ONLY on original Phase 4/5 data. Fetches the same public ChEMBL activities
(or --input-json for offline QA), stores the raw snapshot, and reconstructs
per-document/per-assay candidate contributions from individual activity records.
Does NOT verify biological assay equivalence, source independence or AUC.

Why: Phase 5's multi-document/multi-assay summary is an inventory; its
per-document assay union may be contaminated by information from other
publications concerning the same molecule. Phase 6 uses actual activity links.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, Counter
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, keys):
    with path.open('w', newline='', encoding='utf-8') as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)


def phase4_module(path):
    if not path.exists():
        raise SystemExit(f'Cannot find Phase 4 code: {path}')
    spec = importlib.util.spec_from_file_location('aipharma_phase4_locked', str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def flags_for_description(desc):
    s = (desc or '').lower()
    flags = []
    for group, words in (
        ('CELL_OR_PHENOTYPE_CHECK', ('cell', 'proliferat', 'viab', 'reporter', 'luciferase', 'cytotox')),
        ('TARGET_CONSTRUCT_CHECK', ('mutant', 'mutation', 'fusion', 'domain', 'fragment', 'truncat')),
        ('SPECIES_CHECK', ('mouse', 'murine', 'rat ', 'bovine')),
        ('METHOD_DETAIL_CHECK', ('atp', 'radioactive', 'fluorescen', 'thermal', 'enzym')),
    ):
        if any(k in s for k in words):
            flags.append(group)
    if not s.strip(): flags.append('DESCRIPTION_UNAVAILABLE')
    return ';'.join(flags)


def run(repo: Path, input_json: Path | None, max_pages: int) -> None:
    phase4 = repo/'reports/validation_phase1/phase4_chembl_alk5'
    phase5 = repo/'reports/validation_phase1/phase5_source_audit'
    train_path = repo/'reports/alk5_endpoint_archive_20261016/dataset.json'
    cand_path = phase4/'candidate_external_alk5_NEEDS_REVIEW.csv'
    cur_path = phase4/'curation_audit.csv'
    sum_path = phase4/'summary.json'
    source_path = phase5/'summary.json'
    for p in (train_path,cand_path,cur_path,sum_path,source_path):
        if not p.is_file(): raise SystemExit(f'MISSING: {p}')
    old = json.loads(sum_path.read_text(encoding='utf-8'))
    required_sha = old.get('sha256', {})
    for key,path in [('training',train_path),('curation_audit',cur_path),('candidate_csv',cand_path)]:
        if sha256(path) != required_sha.get(key):
            raise SystemExit(f'Changed Phase 4 input: {key}; cannot compare provenance.')
    if old.get('source') != 'ChEMBL_API_chembl4439' or not old.get('download_complete'):
        if input_json is None:
            raise SystemExit('Phase 4 not recorded as complete official ChEMBL download; supply offline --input-json only for testing.')
    locked = {(r['smiles'], r['label']) for r in read_csv(cand_path)}
    if len(locked) != len(read_csv(cand_path)):
        raise SystemExit('Duplicate Phase 4 candidate SMILES+label.')
    mod = phase4_module(repo/'scripts/benchmark_alk5_phase4.py')
    if input_json:
        raw_obj=json.loads(input_json.read_text(encoding='utf-8'))
        raw=raw_obj['activities'] if isinstance(raw_obj,dict) else raw_obj
        source=f'offline fixture: {input_json.name}'
        complete=False
        pages=0
    else:
        raw,pages,complete=mod.download(max_pages=max_pages)
        source='ChEMBL API current live snapshot'
    if not isinstance(raw,list): raise SystemExit('Expected activity records array.')

    dst=repo/'reports/validation_phase1/phase6_assay_provenance'
    dst.mkdir(parents=True,exist_ok=True)
    snapshot_path=dst/'raw_chembl_activity_snapshot.json'
    snapshot_path.write_text(json.dumps({'activities':raw},ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')

    train=json.loads(train_path.read_text(encoding='utf-8'))
    reasons, meta, strict, _ = mod.build(train,raw)
    replay={(r['smiles'],r['label']) for r in strict}
    drift=bool(locked != replay)
    if drift:
        diff={
          'absent_from_new_snapshot':len(locked-replay),
          'new_in_snapshot':len(replay-locked),
          'example_absent_smiles':[s for s,l in sorted(locked-replay)[:5]],
          'example_added_smiles':[s for s,l in sorted(replay-locked)[:5]],
        }
    else: diff={}

    reasons2=Counter()
    groups=defaultdict(lambda: {'active':set(),'weak':set(),'activity_ids':set(),'descriptions':set(),'measurements':0})
    docs=defaultdict(lambda: {'active':set(),'weak':set(),'assays':set(),'links':set()})
    weak_links=[]
    candidate_activity_rows=0
    for row in raw:
        v=mod.filtered_activity(row,reasons2)
        if v is None or (v['smiles'],str(v['label'])) not in locked:
            continue
        candidate_activity_rows+=1
        key=(v['document'],v['assay'])
        group=groups[key]
        bucket='active' if v['label']==1 else 'weak'
        group[bucket].add(v['smiles'])
        if v['activity_id']:group['activity_ids'].add(v['activity_id'])
        group['measurements']+=1
        desc=v.get('assay_description','')
        if desc:group['descriptions'].add(desc)
        doc=docs[v['document']]
        doc[bucket].add(v['smiles']); doc['assays'].add(v['assay']); doc['links'].add(key)
        if v['label']==0:
            weak_links.append({'smiles':v['smiles'],'document_chembl_id':v['document'],
                               'assay_chembl_id':v['assay'],'activity_id':v['activity_id'],
                               'ic50_nm':v['ic50_nm'],'assay_description':desc[:500],
                               'heuristic_flags':flags_for_description(desc)})

    def item(g,k):
        d,a=k
        descriptions=sorted(g['descriptions'])
        example=descriptions[0] if descriptions else ''
        return {'document_chembl_id':d,'assay_chembl_id':a,'active_molecules':len(g['active']),
                'weak_molecules':len(g['weak']),'mixed_labels':int(bool(g['active'] and g['weak'])),
                'activity_records':g['measurements'],'distinct_activity_ids':len(g['activity_ids']),
                'assay_descriptions_n':len(descriptions),'assay_description_sample':example[:500],
                'heuristic_flags':flags_for_description(example),
                'manual_review_status':'PENDING_EXTERNAL_ASSAY_REVIEW'}
    group_rows=[item(v,k) for k,v in groups.items()]
    group_rows.sort(key=lambda r:(-r['mixed_labels'],-r['weak_molecules'],-r['active_molecules'],r['document_chembl_id'],r['assay_chembl_id']))
    doc_rows=[]
    for d,v in docs.items():
        doc_rows.append({'document_chembl_id':d,'active_molecules':len(v['active']),
                         'weak_molecules':len(v['weak']),'mixed_labels':int(bool(v['active'] and v['weak'])),
                         'n_distinct_assays':len(v['assays']),
                         'n_mixed_assays':sum(int(bool(groups[k]['active'] and groups[k]['weak'])) for k in v['links']),
                         'assay_ids':';'.join(sorted(v['assays'])),
                         'manual_review_status':'PENDING_EXTERNAL_SOURCE_REVIEW'})
    doc_rows.sort(key=lambda r:(-r['mixed_labels'],-r['weak_molecules'],-r['active_molecules'],r['document_chembl_id']))

    mixed_assay_keys={k for k,v in groups.items() if v['active'] and v['weak']}
    weak_in_mixed_assays=set().union(*(groups[k]['weak'] for k in mixed_assay_keys)) if mixed_assay_keys else set()
    assays_by_doc=defaultdict(set)
    for d,a in groups: assays_by_doc[a].add(d)
    shared={a:sorted(x) for a,x in assays_by_doc.items() if len(x)>1}
    status=('SOURCE_SNAPSHOT_DRIFT_NO_VALIDATION' if drift else
            'ASSAY_PROVENANCE_RECONSTRUCTED_PENDING_MANUAL_REVIEW')
    summary={
      'status':status,'source':source,'download_pages':pages,'download_complete':complete,
      'raw_records':len(raw),'raw_snapshot_sha256':sha256(snapshot_path),
      'phase4_candidate_n':len(locked),'replay_candidate_n':len(replay),
      'replay_matches_locked_phase4':not drift,'drift_details':diff,
      'candidate_activity_records':candidate_activity_rows,
      'document_assay_groups':len(groups),'mixed_document_assay_groups':len(mixed_assay_keys),
      'weak_molecules_in_at_least_one_mixed_assay':len(weak_in_mixed_assays),
      'weak_molecules_in_at_least_one_mixed_document':len(set().union(*(v['weak'] for v in docs.values() if v['active'] and v['weak'])) if docs else set()),
      'candidate_unique_documents':len(docs),'mixed_documents':sum(bool(v['active'] and v['weak']) for v in docs.values()),
      'assays_appearing_in_multiple_document_ids':len(shared),
      'shared_assay_examples':dict(list(sorted(shared.items()))[:10]),
      'counts_match_phase5_summary':{
        'documents':len(docs)==json.loads(source_path.read_text(encoding='utf-8'))['unique_source_documents_total'],
        'mixed_documents':sum(bool(v['active'] and v['weak']) for v in docs.values())==json.loads(source_path.read_text(encoding='utf-8'))['source_documents_with_both_classes'],
      },
      'interpretation':[
        'Document-assay links now come from individual ChEMBL activity records; Phase 5 per-document assay unions can over-assign assays to documents.',
        'A mixed-label assay is not automatically biochemically comparable or source-independent.',
        'Descriptions and document IDs require full original method/patent and cross-database review.',
        'This analysis does not run docking, train models, score predictions or estimate external AUC.',
      ],
      'generated_utc':datetime.now(timezone.utc).isoformat(),
    }
    (dst/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    write_csv(dst/'document_assay_groups.csv',group_rows,
              ['document_chembl_id','assay_chembl_id','active_molecules','weak_molecules','mixed_labels',
               'activity_records','distinct_activity_ids','assay_descriptions_n','assay_description_sample',
               'heuristic_flags','manual_review_status'])
    write_csv(dst/'document_review_queue.csv',doc_rows,
              ['document_chembl_id','active_molecules','weak_molecules','mixed_labels','n_distinct_assays',
               'n_mixed_assays','assay_ids','manual_review_status'])
    write_csv(dst/'weak_activity_links.csv',weak_links,
              ['smiles','document_chembl_id','assay_chembl_id','activity_id','ic50_nm',
               'assay_description','heuristic_flags'])
    print('PHASE 6 — ChEMBL ACTIVITY-LEVEL ASSAY PROVENANCE (not validation)')
    print('ChEMBL snapshot complete:',complete,'raw activity records:',len(raw))
    print('Exact Phase 4 candidate replay:', 'MATCH' if not drift else f'DRIFT {diff}')
    print('Candidate document-assay groups:',len(groups),'mixed-label:',len(mixed_assay_keys))
    print('Unique weak molecules in mixed assays:',len(weak_in_mixed_assays))
    print('Documents:',len(docs),'with mixed labels:',summary['mixed_documents'])
    print('Assays linked to >1 document:',len(shared))
    print('OUTPUT',dst)
    print('DO NOT claim external validation before manual assay/source audit.')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo',type=Path,default=Path('.'))
    ap.add_argument('--input-json',type=Path,help='Offline fixture (for testing; does not validate public data)')
    ap.add_argument('--max-pages',type=int,default=0,help='0=entire ChEMBL source')
    args=ap.parse_args()
    if args.max_pages<0: ap.error('max-pages cannot be negative')
    run(args.repo.resolve(),args.input_json,args.max_pages)
