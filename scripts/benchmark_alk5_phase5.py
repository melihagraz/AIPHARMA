#!/usr/bin/env python3
"""AIPHARMA Phase 5: source concentration and chemical-overlap audit.

READ-ONLY on Phase 4 inputs. This is NOT a validation score, independent
assay review, or proof of external dataset independence. Offline, RDKit only.

Outputs under reports/validation_phase1/phase5_source_audit/:
  summary.json, publication_breakdown.csv, candidate_review_queue.csv,
  conservative_candidates_NEED_MANUAL_REVIEW.csv
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized_parent(smiles, tautomer):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        # Largest fragment -> neutral form -> canonical tautomer. This handles
        # many salt/tautomer variants, not every conceivable chemical identity.
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        mol = tautomer.Canonicalize(mol)
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception:
        return None


def split_ids(value):
    return {s.strip() for s in (value or '').split(';') if s.strip()}


def run(repo):
    base = repo / 'reports/validation_phase1/phase4_chembl_alk5'
    trainfile = repo / 'reports/alk5_endpoint_archive_20261016/dataset.json'
    auditfile = base / 'curation_audit.csv'
    candidatefile = base / 'candidate_external_alk5_NEEDS_REVIEW.csv'
    sumfile = base / 'summary.json'
    for p in (trainfile, auditfile, candidatefile, sumfile):
        if not p.is_file():
            raise SystemExit(f'Missing required input: {p}')
    prior = json.loads(sumfile.read_text(encoding='utf-8'))
    expected = prior.get('sha256', {})
    check = {
        'training': trainfile,
        'curation_audit': auditfile,
        'candidate_csv': candidatefile,
    }
    for name, p in check.items():
        if name in expected and sha256(p) != expected[name]:
            raise SystemExit(f'INPUT CHANGED: {name} SHA256 does not match Phase 4 summary; abort.')

    training = json.loads(trainfile.read_text(encoding='utf-8'))
    audit = read_csv(auditfile)
    candidate = read_csv(candidatefile)
    if len(candidate) != int(prior['candidate_molecules_after_overlap_screen']):
        raise SystemExit('Candidate count differs from Phase 4 summary; abort.')
    metadata = {r['smiles']: r for r in audit}
    if len(metadata) != len(audit):
        raise SystemExit('Duplicate canonical SMILES in Phase 4 audit; abort.')
    for r in candidate:
        if r['smiles'] not in metadata:
            raise SystemExit('Candidate not present in original curation audit; abort.')
        if metadata[r['smiles']]['exclusion_flags'].strip():
            raise SystemExit('Flagged molecule in candidate CSV; abort.')
        if r['label'] not in ('0', '1'):
            raise SystemExit('Bad binary candidate label; abort.')

    tautomer = rdMolStandardize.TautomerEnumerator()
    train_parent = defaultdict(set)
    train_molecules = []
    for t in training:
        key = normalized_parent(t['smiles'], tautomer)
        if key:
            train_parent[key].add(int(t['label']))
            m = Chem.MolFromSmiles(key)
            if m:
                train_molecules.append(m)
    fingerprinter = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    train_fp = [fingerprinter.GetFingerprint(m) for m in train_molecules]

    candidate_parents = defaultdict(list)
    standardized = {}
    missing_parents = []
    for r in candidate:
        key = normalized_parent(r['smiles'], tautomer)
        standardized[r['smiles']] = key
        if key:
            candidate_parents[key].append(r)
        else:
            missing_parents.append(r['smiles'])
    conflicting_parents = {
        p for p, rows in candidate_parents.items()
        if len({int(r['label']) for r in rows}) > 1
    }

    source = defaultdict(lambda: {'active':set(), 'weak':set(), 'assays':set(), 'n_record_total':0})
    review = []
    for r in candidate:
        m = metadata[r['smiles']]
        key = standardized[r['smiles']]
        overlap = key in train_parent if key else False
        conflict = key in conflicting_parents if key else False
        siblings = len(candidate_parents.get(key, [])) if key else 0
        mol = Chem.MolFromSmiles(key or r['smiles'])
        if mol and train_fp:
            fp = fingerprinter.GetFingerprint(mol)
            nearest = max(DataStructs.BulkTanimotoSimilarity(fp, train_fp))
        else:
            nearest = None
        docs = split_ids(m['all_documents'])
        assays = split_ids(m['assay_ids'])
        for doc in docs:
            v = source[doc]
            v['active' if r['label'] == '1' else 'weak'].add(r['smiles'])
            v['assays'].update(assays)
            v['n_record_total'] += int(m['n_exact_ic50_records'])
        tags = []
        if key is None: tags.append('UNRESOLVED_STANDARDIZED_STRUCTURE')
        if overlap: tags.append('TRAINING_PARENT_MATCH')
        if conflict: tags.append('CONTRADICTORY_LABELS_AFTER_STANDARDIZATION')
        if siblings > 1: tags.append('MULTIPLE_CANDIDATES_SAME_PARENT')
        if nearest is not None and nearest >= 0.90: tags.append('HIGH_TRAINING_SIMILARITY_GE_0.90')
        if len(docs) > 1: tags.append('MULTIPLE_SOURCE_DOCUMENTS')
        desc = m['assay_descriptions_sample'] or ''
        desc_low = desc.lower()
        if any(k in desc_low for k in ('cell', 'cellular', 'proliferation', 'viability')):
            tags.append('ASSAY_DESCRIPTION_REQUIRES_CELLULAR_REVIEW')
        review.append({
            'smiles':r['smiles'], 'label':r['label'], 'standardized_parent':key or '',
            'parent_same_as_training':int(overlap),
            'parent_conflicting_candidate_labels':int(conflict),
            'same_parent_candidate_rows':siblings,
            'max_tanimoto_to_training':f'{nearest:.5f}' if nearest is not None else '',
            'source_documents':';'.join(sorted(docs)),
            'assay_ids':';'.join(sorted(assays)),
            'n_exact_ic50_records':m['n_exact_ic50_records'],
            'min_ic50_nM':m['min_ic50_nM'],
            'max_ic50_nM':m['max_ic50_nM'],
            'assay_description_sample':desc,
            'audit_flags':';'.join(tags),
            'automatic_candidate_status':(
                'EXCLUDE_POTENTIAL_LEAKAGE' if overlap or conflict or key is None
                else 'CANDIDATE_MANUAL_SOURCE_ASSAY_REVIEW_REQUIRED'),
        })

    review.sort(key=lambda r:(r['label']!='0',
                              r['automatic_candidate_status']!='EXCLUDE_POTENTIAL_LEAKAGE',
                              -float(r['max_tanimoto_to_training'] or 0), r['smiles']))
    conservative = [r for r in review if r['automatic_candidate_status'].startswith('CANDIDATE_')]
    weak_documents = [d for d,s in source.items() if s['weak']]
    mixed_documents = [d for d,s in source.items() if s['weak'] and s['active']]
    source_rows = [{
        'document_chembl_id':doc,
        'active_molecules':len(v['active']), 'weak_molecules':len(v['weak']),
        'all_unique_molecules':len(v['active'] | v['weak']),
        'mixed_labels':int(bool(v['active'] and v['weak'])),
        'distinct_assay_ids':len(v['assays']),
        'assay_ids':';'.join(sorted(v['assays'])),
        'sum_measurement_records_across_molecules':v['n_record_total'],
    } for doc,v in source.items()]
    source_rows.sort(key=lambda v:(-v['weak_molecules'], -v['active_molecules'], v['document_chembl_id']))

    counts = Counter(r['label'] for r in conservative)
    flags = Counter()
    for r in review:
        flags.update(r['audit_flags'].split(';') if r['audit_flags'] else [])
    flags.pop('', None)
    clean_prevalence = (counts['1']/len(conservative)) if conservative else None
    output = {
      'status': 'PRELIMINARY_SOURCE_AND_CHEMICAL_AUDIT_ONLY_NOT_EXTERNAL_VALIDATION',
      'source': 'Phase 4 ChEMBL 4439 candidates; no new download or docking',
      'input_hashes_verified': True,
      'phase4_candidate_n':len(candidate),
      'phase4_class_counts':dict(Counter(r['label'] for r in candidate)),
      'automated_parent_screen_remaining_n':len(conservative),
      'automated_parent_screen_remaining_class_counts':dict(counts),
      'automated_parent_screen_active_prevalence':clean_prevalence,
      'remaining_average_precision_no_skill_baseline':clean_prevalence,
      'overlap_with_training_after_parent_standardization':sum(r['parent_same_as_training'] for r in review),
      'conflicting_labels_after_parent_standardization':sum(r['parent_conflicting_candidate_labels'] for r in review),
      'unresolved_standardization':len(missing_parents),
      'unique_source_documents_total':len(source),
      'source_documents_with_weak_molecules':len(weak_documents),
      'source_documents_with_both_classes':len(mixed_documents),
      'unique_weak_molecules_with_any_mixed_document':len({r['smiles'] for r in review if r['label']=='0' and
                         any(d in mixed_documents for d in split_ids(r['source_documents']))}),
      'flag_counts_nonexclusive':dict(flags),
      'top_documents_by_weak_count':source_rows[:15],
      'warnings':[
        'ChEMBL source overlap was checked with exact ChEMBL publication IDs only; cross-database papers/patents remain unresolved.',
        'Fragment-parent/tautomer standardization does not prove identity equivalence for all salts, stereoisomers, and structural conventions.',
        'Near-identical analogues and source-specific label effects may remain. Similarity flags are diagnostic, not hard exclusions.',
        'Descriptions are truncated Phase 4 samples; assay targets, constructs, mutant status, and biochemical comparability require manual review.',
        'Molecule-based ROC-AUC confidence intervals can be overconfident when structures are clustered within publications.',
        'No model fitting, evaluation ROC-AUC, experimental validation or software superiority tested here.',
      ],
      'rdkit_version':rdBase.rdkitVersion,
      'input_sha256':{k:sha256(v) for k,v in check.items()},
    }
    out = repo/'reports/validation_phase1/phase5_source_audit'
    out.mkdir(parents=True,exist_ok=True)
    (out/'summary.json').write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    write_csv(out/'publication_breakdown.csv', source_rows, list(source_rows[0]) if source_rows else [
        'document_chembl_id','active_molecules','weak_molecules','all_unique_molecules','mixed_labels','distinct_assay_ids','assay_ids','sum_measurement_records_across_molecules'])
    columns = list(review[0]) if review else ['smiles','label','audit_flags']
    write_csv(out/'candidate_review_queue.csv', review, columns)
    write_csv(out/'conservative_candidates_NEED_MANUAL_REVIEW.csv', conservative, columns)
    print('PHASE 5 — SOURCE + CHEMICAL PARENT AUDIT (NOT external validation)')
    print('INPUT phase4 candidates:',len(candidate),'active',sum(r['label']=='1' for r in candidate),'weak',sum(r['label']=='0' for r in candidate))
    print('Training parent overlap:',output['overlap_with_training_after_parent_standardization'])
    print('Conflicts after parent standardization:',output['conflicting_labels_after_parent_standardization'])
    print('Remaining candidates, STILL manual-review only:',len(conservative),'active',counts['1'],'weak',counts['0'])
    print('Publication groups:',len(source),'with any weak:',len(weak_documents),'with both classes:',len(mixed_documents))
    print('Remaining positive prevalence / no-skill AP baseline:',
        f'{clean_prevalence:.4f}' if clean_prevalence is not None else 'N/A')
    print('OUTPUT',out)
    print('CAUTION: Do not report external validation metrics from these candidates before source/assay review.')

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path('.'))
    args = p.parse_args()
    run(args.repo.resolve())
