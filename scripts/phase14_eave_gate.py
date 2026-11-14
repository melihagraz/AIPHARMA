#!/usr/bin/env python3
"""Opt-in AIPHARMA EAVE quality gate: audit, or launch-and-audit a frozen benchmark.

Examples:
  python scripts/phase14_eave_gate.py verify --all
  python scripts/phase14_eave_gate.py verify --case alk5-3hmm --claim experimental_potency
  python scripts/phase14_eave_gate.py run --case alk5-3hmm
  python scripts/phase14_eave_gate.py ablate
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from aipharma_eave.audit import CASES
from aipharma_eave.phase14 import gate, run


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--baseline-dir',type=Path,default=None)
    sub=p.add_subparsers(dest='action',required=True)
    verify=sub.add_parser('verify',help='Read-only EAVE gate for completed benchmark(s)')
    grp=verify.add_mutually_exclusive_group(required=True)
    grp.add_argument('--all',action='store_true')
    grp.add_argument('--case',choices=sorted(CASES))
    r=sub.add_parser('run',help='Launch selected existing docking script and audit upon completion')
    r.add_argument('--case',choices=sorted(CASES),required=True)
    r.add_argument('--limit',type=int,default=0,help='Dock at most N not-yet-scored molecules; 0=all remaining')
    ab=sub.add_parser('ablate',help='Synthetic-fixture detection ablation: no real scientific inputs modified')
    for cmd in (verify,r):
        cmd.add_argument('--bootstrap',type=int,default=1000)
        cmd.add_argument('--claim',choices=('provenance_verified','software_execution',
             'pose_reproduction','retrospective_discrimination','experimental_potency',
             'clinical_efficacy','generalization','receptor_superiority'))
        cmd.add_argument('--allow-limited',action='store_true',help='Accept limited_support for the explicitly requested claim')
    args=p.parse_args(argv)
    root=args.root.resolve()
    bd=args.baseline_dir.resolve() if args.baseline_dir else None
    if args.action=='ablate':
        from aipharma_eave.phase14_ablation import main as do_ablation
        return do_ablation(root)
    if args.bootstrap != 0 and args.bootstrap < 100:
        p.error('--bootstrap must be 0 or >=100')
    if args.action=='run':
        try: rc,out=run(root,args.case,bd,args.bootstrap,args.claim,args.allow_limited,args.limit)
        except (OSError,RuntimeError,ValueError) as exc:
            print(f'ERROR: {exc}',file=sys.stderr);return 3
        print(json.dumps({'case':args.case,'verdict':out['gate_verdict'],
             'auc':out.get('computed_auc'),'output':str(root/'reports/eave/phase14'/args.case/'gate_report.json')},indent=2))
        return rc
    codes=[]
    for case in CASES if args.all else [args.case]:
        try: rc,out=gate(root,case,bd,args.bootstrap,args.claim,args.allow_limited)
        except (OSError,RuntimeError,ValueError) as exc:
            print(f'ERROR {case}: {exc}',file=sys.stderr);return 3
        codes.append(rc)
        print(f"{case}: {out['gate_verdict']}; AUC={out['computed_auc']}; "
              f"baseline_changes={len(out['baseline_changed_artifacts'])}" +
              (f"; requested_claim_approved={out['claim_is_approved']}" if args.claim else ''))
        if rc:
            details = (out['baseline_error'] or
                       (str(out['baseline_changed_artifacts'])[:400] if out['baseline_changed_artifacts'] else None) or
                       f"EAVE={out['underlying_eave_status']}; claim={out['claim_decision']}")
            print('  DETAILS:',details)
    return max(codes)

if __name__=='__main__':raise SystemExit(main())
