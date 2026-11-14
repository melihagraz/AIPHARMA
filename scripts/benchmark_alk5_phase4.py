#!/usr/bin/env python3
"""AIPHARMA Phase 4: ChEMBL ALK5 independent-*candidate* dataset audit.

DOES NOT establish external independence or validate activity. Fetches measured
human TGFBR1/ALK5 IC50 records, constructs conservatively separated potency
classes, and flags historical training publication/scaffold/identity overlaps.

Default target CHEMBL4439 is the human TGF-beta receptor type I / ALK5.
Use the manual protocol review before feeding candidates to any evaluator.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path
from urllib.error import URLError, HTTPError
from urllib.parse import urljoin
from urllib.request import Request, urlopen

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

TARGET = 'CHEMBL4439'
API = ('https://www.ebi.ac.uk/chembl/api/data/activity.json'
       '?target_chembl_id=CHEMBL4439&standard_type=IC50&limit=500')


def canon_scaffold(smiles):
    if not smiles or not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms() < 3:
        return None
    csmi = Chem.MolToSmiles(mol, isomericSmiles=True)
    scaf = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
    return csmi, scaf


def get_json(url, tries=3):
    for n in range(tries):
        try:
            req = Request(url, headers={'User-Agent': 'AIPHARMA-validation-audit/1.0 (academic)'})
            with urlopen(req, timeout=40) as resp:
                return json.load(resp)
        except (URLError, HTTPError, TimeoutError, json.JSONDecodeError) as exc:
            if n == tries - 1:
                raise RuntimeError(f'ChEMBL download failed at {url}: {exc}') from exc
            time.sleep(2 * (n+1))


def download(max_pages=0):
    url = API
    all_rows = []
    pages = 0
    complete = False
    while url:
        block = get_json(url)
        if 'activities' not in block or 'page_meta' not in block:
            raise ValueError('Unexpected ChEMBL API response schema')
        all_rows.extend(block['activities'])
        pages += 1
        nxt = block['page_meta'].get('next')
        url = urljoin(url, nxt) if nxt else None
        if pages % 10 == 0:
            print(f'Fetched {pages} pages; {len(all_rows)} records', file=sys.stderr)
        if max_pages and pages >= max_pages and url:
            print('WARNING: max-pages limit reached; incomplete snapshot', file=sys.stderr)
            break
    else:
        complete = True
    if url is None:
        complete = True
    return all_rows, pages, complete


def filtered_activity(rec, reasons):
    # Exact, uncensored measurements only. Deliberately don't infer negatives
    # from missing records or presumed decoy status.
    if rec.get('target_chembl_id') != TARGET:
        reasons['wrong_target'] += 1; return None
    if rec.get('target_organism') != 'Homo sapiens':
        reasons['nonhuman_or_unknown'] += 1; return None
    if rec.get('assay_type') != 'B':
        reasons['not_binding_assay_type'] += 1; return None
    if rec.get('standard_type') != 'IC50' or rec.get('standard_units') != 'nM':
        reasons['wrong_endpoint_or_units'] += 1; return None
    if str(rec.get('standard_relation') or rec.get('relation')) != '=':
        reasons['censored_or_ambiguous_relation'] += 1; return None
    if rec.get('data_validity_comment'):
        reasons['data_validity_comment'] += 1; return None
    try:
        nm = float(rec.get('standard_value'))
        if not (math.isfinite(nm) and nm > 0):
            raise ValueError('not finite positive')
    except (TypeError, ValueError):
        reasons['invalid_value'] += 1; return None
    # Locked thresholds mirror the historical training archive exactly:
    # >= pIC50 7 : IC50 <= 100 nM; <= pIC50 5 : IC50 >= 10000 nM.
    if nm <= 100.0:
        label = 1
    elif nm >= 10000.0:
        label = 0
    else:
        reasons['grey_zone_100_to_10000_nM'] += 1; return None
    structural = canon_scaffold(rec.get('canonical_smiles'))
    if structural is None:
        reasons['invalid_or_missing_smiles'] += 1; return None
    doc = str(rec.get('document_chembl_id') or '').strip()
    assay = str(rec.get('assay_chembl_id') or '').strip()
    if not doc or not assay:
        reasons['missing_document_or_assay'] += 1; return None
    return {
        'smiles': structural[0], 'scaffold':structural[1], 'label':label,
        'ic50_nm':nm, 'document':doc, 'assay':assay,
        'activity_id':str(rec.get('activity_id') or ''),
        'assay_description':str(rec.get('assay_description') or ''),
    }


def build(train, records):
    known_smiles, known_scaffolds, known_docs, known_assays = set(),set(),set(),set()
    for ix, row in enumerate(train):
        result = canon_scaffold(row.get('smiles'))
        if result is None:
            raise ValueError(f'Invalid training molecule at row {ix}')
        csmi, scaf = result
        known_smiles.add(csmi)
        if scaf: known_scaffolds.add(scaf)
        known_docs.update(str(v) for v in row.get('documents', []) if v)
        known_assays.update(str(v) for v in row.get('assay_ids', []) if v)
    reasons = collections.Counter()
    mols = {}
    for r in records:
        obj = filtered_activity(r, reasons)
        if obj is not None:
            mols.setdefault(obj['smiles'], []).append(obj)
    meta_rows, export = [], []
    for smiles, obs in sorted(mols.items()):
        scaf = obs[0]['scaffold']
        labels = set(o['label'] for o in obs)
        docs = sorted(set(o['document'] for o in obs))
        assays = sorted(set(o['assay'] for o in obs))
        overlap_identity = smiles in known_smiles
        overlap_scaffold = (bool(scaf) and scaf in known_scaffolds)
        overlap_doc = bool(set(docs) & known_docs)
        overlap_assay = bool(set(assays) & known_assays)
        conflict = len(labels) > 1
        unresolved_scaffold = not bool(scaf)
        flags = [name for name, active in (
            ('conflicting_labels', conflict),
            ('identity', overlap_identity),
            ('scaffold', overlap_scaffold),
            ('training_publication', overlap_doc),
            ('training_assay', overlap_assay),
            ('unresolved_scaffold', unresolved_scaffold)
        ) if active]
        result = {
            'smiles':smiles,
            'label':str(obs[0]['label']) if not conflict else '',
            'source_id':docs[0],
            'all_documents':';'.join(docs),
            'assay_ids':';'.join(assays),
            'n_exact_ic50_records':len(obs),
            'min_ic50_nM':min(o['ic50_nm'] for o in obs),
            'max_ic50_nM':max(o['ic50_nm'] for o in obs),
            'exclusion_flags':';'.join(flags),
            'assay_descriptions_sample':obs[0]['assay_description'][:200]
        }
        meta_rows.append(result)
        if not flags:
            export.append({k:result[k] for k in ('smiles','label','source_id')})
    return reasons, meta_rows, export, {'train_n':len(train),'train_distinct_smiles':len(known_smiles),
                                        'train_scaffolds':len(known_scaffolds),'train_publications':len(known_docs)}


def write_csv(path, rows, columns):
    with path.open('w',encoding='utf-8',newline='') as fh:
        wr=csv.DictWriter(fh,fieldnames=columns,extrasaction='ignore')
        wr.writeheader();wr.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('.'))
    p.add_argument('--input-json',type=Path, help='Offline fixture with ChEMBL-style activities array; no download')
    p.add_argument('--max-pages',type=int,default=0,help='0=fetch all pages (default); positive may truncate')
    args=p.parse_args()
    if args.max_pages < 0:
        p.error('--max-pages must be >=0')
    repo=args.repo.resolve()
    train_path=repo/'reports/alk5_endpoint_archive_20261016/dataset.json'
    if not train_path.is_file():
        p.error(f'Training dataset missing: {train_path}')
    train=json.loads(train_path.read_text(encoding='utf-8'))
    if args.input_json:
        raw=json.loads(args.input_json.read_text(encoding='utf-8'))
        records=raw['activities'] if isinstance(raw,dict) else raw
        pages=0;complete=False;origin='local_fixture_not_external'
    else:
        records,pages,complete=download(args.max_pages)
        origin='ChEMBL_API_chembl4439'
    reasons,meta,strict,training=build(train,records)
    dst=repo/'reports/validation_phase1/phase4_chembl_alk5'
    dst.mkdir(parents=True,exist_ok=True)
    all_path=dst/'curation_audit.csv'
    cand_path=dst/'candidate_external_alk5_NEEDS_REVIEW.csv'
    write_csv(all_path,meta, ['smiles','label','source_id','all_documents','assay_ids',
                               'n_exact_ic50_records','min_ic50_nM','max_ic50_nM',
                               'exclusion_flags','assay_descriptions_sample'])
    write_csv(cand_path,strict,['smiles','label','source_id'])
    cnt=collections.Counter(r['label'] for r in strict)
    reject=collections.Counter()
    for r in meta:
        for flag in r['exclusion_flags'].split(';'):
            if flag: reject[flag]+=1
    summary={
      'status':'CANDIDATE_DATA_ONLY_REQUIRES_MANUAL_ASSAY_AND_CROSS_DB_SOURCE_REVIEW',
      'target':'human ALK5 / TGFBR1 (CHEMBL4439)',
      'source':origin,'download_pages':pages,'download_complete':complete,
      'retrieved_activity_records':len(records),
      'eligible_activity_measurements':sum(int(r['n_exact_ic50_records']) for r in meta),
      'unique_molecules_after_measurement_filters':len(meta),
      'candidate_molecules_after_overlap_screen':len(strict),
      'candidate_active':cnt.get('1',0),'candidate_weak':cnt.get('0',0),
      'rejections_at_measurement_level':dict(reasons),
      'flags_at_molecule_level':dict(reject),
      'training':training,
      'label_rule':'exact human ALK5 IC50 <=100 nM active; >=10000 nM weak; 100< IC50 <10000 excluded',
      'warning':'This is NOT an external-validation result. Check assay descriptions and protocol, duplicate papers across databases, patents, structures/salts/tautomers, and training overlap before scoring.',
      'raw_source_url':API,
      'sha256':{'training':sha(train_path),'curation_audit':sha(all_path),'candidate_csv':sha(cand_path)},
      'rdkit_version':rdBase.rdkitVersion,
    }
    (dst/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PHASE 4 CANDIDATE ALK5 AUDIT')
    print('OUTPUT',dst)
    print('Downloaded',len(records),'records in',pages,'pages; complete=',complete)
    print('Eligible measurement records',summary['eligible_activity_measurements'])
    print('Unique molecules after measurement filtering',len(meta))
    print('New candidate molecules after overlap screen',len(strict),
          '(active=',cnt.get('1',0),'weak=',cnt.get('0',0),')')
    print('Top flagged reasons',reject.most_common(6))
    print('IMPORTANT: not yet an independent external validation; manual assay and source audit required.')

if __name__=='__main__':
    main()
