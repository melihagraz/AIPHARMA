#!/usr/bin/env python3
"""AIPHARMA EAVE evidence audit (read-only against scientific inputs)."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

# Work when executed by relative path from root, even if PYTHONPATH isn't set.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from aipharma_eave.audit import CASES, Audit


def markdown(report):
    lines = [f"# EAVE scientific audit: {report['case']}", '',
             f"**Overall integrity status:** `{report['audit_status']}`", '',
             '## Recorded evidence and independent calculations', '']
    for k,v in report['metrics'].items():
        lines.append(f'- **{k}:** {v}')
    lines.extend(['', '## Evidence-bound claims', ''])
    for k,v in report['claims'].items():
        lines.append(f"- **{k}** — `{v['status']}`: {v['reason']}")
    lines.extend(['', '## Integrity checks', ''])
    for x in report['checks']:
        lines.append(f"- `{x['status']}` **{x['check']}**: {x['detail']}")
    lines.extend(['','## Limits', ''])
    for x in report['limitations']:
        lines.append('- '+x)
    return '\n'.join(lines) + '\n'


def atomic_text(p, content):
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name+'.tmp')
    tmp.write_text(content, encoding='utf-8')
    tmp.replace(p)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--case', choices=sorted(CASES))
    g.add_argument('--all', action='store_true')
    p.add_argument('--root', type=Path, default=Path.cwd())
    p.add_argument('--out', type=Path, help='Output directory (default: <root>/reports/eave)')
    p.add_argument('--bootstrap', type=int, default=10000, help='0 skips interval recomputation; otherwise >= 100')
    p.add_argument('--reference-report', type=Path, help='Trusted previously exported EAVE JSON report to detect changes (single case only)')
    p.add_argument('--claim-type', choices=('provenance_verified','software_execution','pose_reproduction',
       'retrospective_discrimination','experimental_potency','clinical_efficacy','generalization','receptor_superiority'),
       help='Print a structured evidence decision for this claim in addition to the audit')
    args=p.parse_args(argv)
    if args.bootstrap != 0 and args.bootstrap < 100:
        p.error('--bootstrap must be 0 or >=100')
    if args.reference_report and args.all:
        p.error('--reference-report requires --case, not --all')
    root=args.root.resolve(); dest=(args.out or (root/'reports/eave')).resolve()
    case_list=list(CASES) if args.all else [args.case]
    print('EAVE read-only audit: reports and scores will not be modified.')
    results=[]
    for case in case_list:
        result=Audit(root,case,args.bootstrap).run()
        if args.reference_report:
            try:
                previous=json.loads(args.reference_report.read_text(encoding='utf-8'))
                if previous.get('case') != case or previous.get('audit_status') != 'verified':
                    raise ValueError('Reference must be a verified report for the SAME case')
                original=previous.get('artifact_sha256')
                if not isinstance(original,dict) or not original:
                    raise ValueError('No stored SHA256 artifact record in reference')
                changed={k: {'original':v, 'current':result['artifact_sha256'].get(k)}
                         for k,v in original.items() if result['artifact_sha256'].get(k) != v}
                result['checks'].append({'check':'reference_snapshot','status':'fail' if changed else 'verified',
                                         'detail': ('Snapshot differs: '+json.dumps(changed)) if changed else 'Reference report SHA256 fingerprints unchanged'})
                if changed:
                    result['audit_status']='fail'
                    result['counts']['fail']=result['counts'].get('fail',0)+1
                    for name in ('provenance_verified','software_execution','pose_reproduction','retrospective_discrimination','receptor_superiority'):
                        result['claims'][name]['status']='not_established'
                        result['claims'][name]['reason']='Stored reference snapshot differs; audit is not trustworthy until resolved'
                else:
                    result['counts']['verified']=result['counts'].get('verified',0)+1
            except (OSError,ValueError,json.JSONDecodeError) as e:
                result['checks'].append({'check':'reference_snapshot','status':'fail','detail':f'Invalid trusted reference: {e}'})
                result['audit_status']='fail'
                result['counts']['fail']=result['counts'].get('fail',0)+1
                for name in ('provenance_verified','software_execution','pose_reproduction','retrospective_discrimination','receptor_superiority'):
                    result['claims'][name]['status']='not_established'
                    result['claims'][name]['reason']='Cannot validate against requested reference snapshot'
        out=dest/case
        atomic_text(out/'evidence_report.json',json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+'\n')
        atomic_text(out/'evidence_report.md',markdown(result))
        print(f"{case}: {result['audit_status']}; AUC={result['metrics'].get('auc_recomputed', 'unavailable')}; "
              f"report={out/'evidence_report.md'}")
        if args.claim_type:
            print(f"  claim {args.claim_type}: {result['claims'][args.claim_type]['status']}: {result['claims'][args.claim_type]['reason']}")
        results.append(result)
    if args.all:
        atomic_text(dest/'index.json',json.dumps([{k:r[k] for k in ('case','audit_status','metrics')} for r in results],indent=2)+'\n')
    return 3 if any(r['audit_status']=='fail' for r in results) else (2 if any(r['audit_status']=='incomplete' for r in results) else 0)

if __name__=='__main__':
    raise SystemExit(main())
