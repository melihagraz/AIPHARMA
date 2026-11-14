#!/usr/bin/env python3
"""Phase 15: isolated real-docking replay with evidence-aware post-run checks.

This is an optional live integration test, not a new drug-discovery benchmark.
It calls the existing Phase 11 ALK5 runner on an isolated copy in a temporary
workspace and never rewrites the original scientific score files.

Commands (from AIPHARMA root, within its WSL Python environment):
  python scripts/phase15_live_replay.py check
  python scripts/phase15_live_replay.py run
  python scripts/phase15_live_replay.py run --keep-workspace

The Phase 14 snapshot gate pins *old outputs* and is therefore intended for
verifying already-completed archived runs. Phase 15 instead pins immutable
inputs before/after a fresh computation, audits the newly generated outputs,
and re-checks that original archived outputs were not modified.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from aipharma_eave.audit import Audit, CASES
from aipharma_eave.phase14 import gate, immutable_preflight

CASE = 'alk5-1vjy'
FIELDS = ('pair_id', 'class', 'smiles', 'source_id')
SCORES_REL = Path(CASES[CASE]['output']) / 'actual_vina_scores.csv'
REPORT_DIR = Path('reports/eave/phase15')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', encoding='utf-8', dir=path.parent,
                                     prefix=path.name+'.tmp-', delete=False) as tmp:
        json.dump(payload, tmp, indent=2, sort_keys=True, ensure_ascii=False)
        tmp.write('\n')
        tmpname = tmp.name
    os.replace(tmpname, path)


def read_csv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(newline='', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open('w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def row_key(row: dict) -> tuple[str, ...]:
    return tuple(row.get(f, '') for f in FIELDS)


def row_score(row: dict) -> float:
    result = float(row['score_kcal_mol'])
    if not math.isfinite(result):
        raise ValueError('Non-finite docking affinity')
    return result


def remove_one_row(path: Path, pair_id: str, klass: str) -> tuple[dict, dict]:
    fields, rows = read_csv(path)
    if len(rows) != 100:
        raise ValueError(f'Replay needs a completed 100-row baseline, found {len(rows)}')
    hits = [row for row in rows if row['pair_id'] == pair_id and row['class'] == klass]
    if len(hits) != 1:
        raise ValueError(f'Replay target {pair_id}/{klass} not uniquely present')
    removed = hits[0]
    original = {row_key(row): row for row in rows}
    if len(original) != 100:
        raise ValueError('Original score file contains duplicate rows')
    for row in rows:
        row_score(row)
    keep = [row for row in rows if row is not removed]
    write_csv(path, fields, keep)
    return removed, {row_key(row): row for row in keep}


def validate_completed_replay(path: Path, removed: dict, preserved: dict) -> float:
    _, final_rows = read_csv(path)
    if len(final_rows) != 100:
        raise ValueError(f'Expected 100 finite scores after replay, found {len(final_rows)}')
    mapping = {row_key(row): row for row in final_rows}
    if len(mapping) != 100 or set(mapping) != set(preserved) | {row_key(removed)}:
        raise ValueError('Replay result has missing, duplicate, or unexpected molecule identities')
    for key, old in preserved.items():
        new = mapping[key]
        if new != old:
            raise ValueError(f'Previously completed molecule was modified: {key}')
        row_score(new)
    return row_score(mapping[row_key(removed)])


def fingerprint_files(root: Path, files: list[str]) -> dict[str, str | None]:
    return {rel: sha256(root / rel) if (root / rel).is_file() else None for rel in files}


def copies_for_case(root: Path, target: Path) -> list[str]:
    """Minimal copy; intentionally omit all writable locks and EAVE reports."""
    c = CASES[CASE]
    paths = [
        f"{c['sources']}/active.ism", f"{c['sources']}/decoy.ism",
        f"{c['sources']}/source_manifest_hf.json",
        c['panel'], c['selection'], c['descriptors'],
        f"{c['output']}/actual_vina_scores.csv",
        f"{c['output']}/execution_manifest.json",
        f"{c['output']}/analysis.json",
    ]
    redock = root / c['redock']
    if not redock.is_dir():
        raise FileNotFoundError(f'Missing original redocking directory: {redock}')
    for original_file in sorted(redock.iterdir()):
        if original_file.is_file():
            paths.append(str(original_file.relative_to(root)))
    for rel in paths:
        original = root / rel
        if not original.is_file() or original.is_symlink():
            raise FileNotFoundError(f'Missing or unsafe replay input: {original}')
        target_file = target / rel
        target_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target_file)
        if sha256(original) != sha256(target_file):
            raise RuntimeError(f'Copy hash mismatch: {rel}')
    return paths


def check(root: Path) -> dict:
    requested = [
        'scripts/phase11_alk5.py', 'scripts/phase14_eave_gate.py',
        'aipharma_eave/audit.py', 'aipharma_eave/phase14.py',
        str(SCORES_REL)
    ]
    missing = [rel for rel in requested if not (root/rel).is_file()]
    current_gate = None
    error = None
    try:
        code, out = gate(root, CASE, bootstrap=0, write=False)
        current_gate = {'status': out['gate_verdict'], 'returncode': code,
                        'baseline_changes': len(out['baseline_changed_artifacts'])}
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    preflight = immutable_preflight(root, CASE)
    ready = not missing and current_gate is not None and current_gate['returncode'] == 0 and preflight['ok']
    return {'ready': bool(ready), 'case':CASE, 'missing_files':missing,
            'archive_gate':current_gate,'immutable_input_preflight':preflight,
            'error':error,
            'note':'Only existing local files were inspected; no docking executed.'}


def run_live(root: Path, pair_id: str, klass: str, bootstrap: int, keep_workspace: bool) -> tuple[int, dict]:
    report_dir = root / REPORT_DIR
    report_dir.mkdir(parents=True, exist_ok=True)
    lock_path = report_dir / '.phase15_live.lock'
    with lock_path.open('a+') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another Phase 15 replay is already running')
        start = datetime.now(timezone.utc).isoformat()
        report = {'schema':'aipharma-phase15-live-replay-v1', 'case':CASE,
                  'started_utc':start, 'mode':'isolated_single_molecule_real_runner',
                  'test_molecule':{'pair_id':pair_id,'class':klass},
                  'original_input_or_score_files_modified':None,
                  'experiment_recomputed':'one existing benchmark molecule, not a new independent benchmark',
                  'scientific_limitations':['Passing the replay does not validate experimental potency.',
                      'Local EAVE baseline snapshots are unsigned.',
                      'This is an opt-in integration test and does not replace direct-legacy-script controls.']}
        root_score = root / SCORES_REL
        source_snapshot = None
        tmpdir = None
        result_code = 3
        try:
            readiness = check(root)
            if not readiness['ready']:
                raise RuntimeError(f'Original project preflight failed: {readiness}')
            # Pin the full original EAVE baseline artifact inventory including output scores.
            old_pin_path = root/'phase14_eave_baselines'/f'{CASE}.json'
            pinned = json.loads(old_pin_path.read_text(encoding='utf-8'))['artifact_sha256']
            source_snapshot = fingerprint_files(root, list(pinned))
            if any(source_snapshot.get(path) != expected for path, expected in pinned.items()):
                raise RuntimeError('Original project does not match its Phase14 baseline')
            tmpdir = Path(tempfile.mkdtemp(prefix='replay-', dir=report_dir))
            copies_for_case(root, tmpdir)
            # Do not touch the original root score table.
            removed, preserved = remove_one_row(tmpdir/SCORES_REL,pair_id,klass)
            report['original_affinity_kcal_mol'] = row_score(removed)
            frozen_rel = [rel for rel in pinned if not rel.startswith(CASES[CASE]['output'].rstrip('/') + '/')]
            frozen_before = fingerprint_files(tmpdir, frozen_rel)
            cmd = [sys.executable,str(root/'scripts/phase11_alk5.py'),
                   '--root',str(tmpdir),'dock','--limit','1']
            env = os.environ.copy()
            env['PYTHONPATH'] = str(root) + (os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
            report['run_command'] = cmd
            print('Running original Phase11 docking script on an isolated copy (one molecule)...', flush=True)
            with ((report_dir/'live_runner_stdout.log').open('w', encoding='utf-8') as out,
                  (report_dir/'live_runner_stderr.log').open('w', encoding='utf-8') as err):
                proc = subprocess.Popen(cmd, cwd=tmpdir, env=env, stdout=out, stderr=err,
                                        start_new_session=True)
                try:
                    subcode = proc.wait()
                except KeyboardInterrupt:
                    try: os.killpg(proc.pid, signal.SIGINT)
                    except ProcessLookupError: pass
                    try: proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        try: os.killpg(proc.pid, signal.SIGTERM)
                        except ProcessLookupError: pass
                        try: proc.wait(timeout=10)
                        except subprocess.TimeoutExpired:
                            try: os.killpg(proc.pid, signal.SIGKILL)
                            except ProcessLookupError: pass
                            proc.wait()
                    raise
            report['runner_returncode'] = subcode
            if subcode != 0:
                raise RuntimeError(f'Phase11 docking runner returned {subcode}; check log files')
            replay_affinity = validate_completed_replay(tmpdir/SCORES_REL,removed,preserved)
            report['replayed_affinity_kcal_mol'] = replay_affinity
            report['score_difference_replayed_minus_original'] = round(replay_affinity-row_score(removed), 4)
            report['preserved_successful_rows'] = len(preserved)
            if fingerprint_files(tmpdir,frozen_rel) != frozen_before:
                raise RuntimeError('Frozen inputs changed inside the isolated replay workspace')
            # EAVE checks newly assembled 100-row completed run without comparing
            # its *new* score CSV checksum to the immutable old score CSV checksum.
            audit = Audit(tmpdir, CASE, bootstrap).run()
            report['postrun_eave_audit_status'] = audit['audit_status']
            report['postrun_auc_recomputed'] = audit.get('metrics',{}).get('auc_recomputed')
            report['postrun_check_counts'] = {
                state:sum(f['status']==state for f in audit['checks'])
                for state in ('verified','unverified','fail')
            }
            if audit['audit_status'] != 'verified':
                raise RuntimeError(f'Post-run EAVE not verified: {report["postrun_check_counts"]}')
            code, original_gate = gate(root, CASE, bootstrap=0, write=False)
            if code != 0 or original_gate['gate_verdict'] != 'verified':
                raise RuntimeError('Original archived Phase14 quality gate changed')
            if fingerprint_files(root,list(pinned)) != source_snapshot:
                raise RuntimeError('One or more original pinned artifacts changed during replay')
            report['original_input_or_score_files_modified'] = False
            report['verdict'] = 'verified_real_replay'
            result_code = 0
        except (OSError,ValueError,RuntimeError,KeyError,TypeError,KeyboardInterrupt) as exc:
            report['verdict'] = 'blocked'
            report['error'] = f'{type(exc).__name__}: {exc}'
            if source_snapshot is not None:
                report['original_input_or_score_files_modified'] = (
                    fingerprint_files(root,list(source_snapshot)) != source_snapshot)
        finally:
            report['ended_utc'] = datetime.now(timezone.utc).isoformat()
            if tmpdir is not None:
                if keep_workspace:
                    report['scratch_workspace'] = str(tmpdir)
                else:
                    shutil.rmtree(tmpdir, ignore_errors=True)
                    report['scratch_workspace'] = None
            atomic_json(report_dir/'phase15_live_report.json',report)
        return result_code, report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=Path.cwd(),help='Original local AIPHARMA project root')
    sub = ap.add_subparsers(dest='action',required=True)
    sub.add_parser('check',help='Read-only archive and immutable input checks')
    run = sub.add_parser('run',help='Actually replay one molecule with original Phase11/Vina runner')
    run.add_argument('--pair',default='pair_001')
    run.add_argument('--class',dest='klass',choices=['active','decoy'],default='active')
    run.add_argument('--bootstrap',type=int,default=1000)
    run.add_argument('--keep-workspace',action='store_true')
    args=ap.parse_args(argv)
    root=args.root.resolve()
    if args.action == 'check':
        status=check(root)
        print(json.dumps(status,indent=2))
        return 0 if status['ready'] else 2
    if args.bootstrap not in (0,) and args.bootstrap < 100:
        ap.error('--bootstrap must be 0 or >=100')
    if not args.pair.startswith('pair_') or len(args.pair)!=8 or not args.pair[-3:].isdigit():
        ap.error('--pair must look like pair_001')
    try:
        code, report = run_live(root,args.pair,args.klass,args.bootstrap,args.keep_workspace)
        print(json.dumps(report,indent=2,ensure_ascii=False))
        return code
    except (OSError,RuntimeError,ValueError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 3

if __name__ == '__main__':
    raise SystemExit(main())
