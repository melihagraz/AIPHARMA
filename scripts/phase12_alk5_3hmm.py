#!/usr/bin/env python3
"""AIPHARMA Phase 12: same frozen Phase 11 ALK5/TGFR1 panel, alternate 3HMM receptor.

This script performs REAL crystallographic redocking and optionally REAL Vina screening.
No fake/placeholder scores. Reuses, but NEVER modifies, the Phase 11 1VJY results.

Run from project root, using the AIPHARMA WSL2 virtual environment:
  PYTHONPATH="$PWD" python scripts/phase12_alk5_3hmm.py redock
  PYTHONPATH="$PWD" python scripts/phase12_alk5_3hmm.py check
  PYTHONPATH="$PWD" python scripts/phase12_alk5_3hmm.py dock --limit 10
  PYTHONPATH="$PWD" python scripts/phase12_alk5_3hmm.py dock
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import random
import signal
import sys
from datetime import datetime, timezone

PANEL_REL = Path('reports/phase11/alk5_matched50/matched_panel.csv')
BASELINE_REL = Path('reports/phase11/alk5_matched50_docking/actual_vina_scores.csv')
REDOCK_REL = Path('reports/phase12/redocking/3HMM')
OUT_REL = Path('reports/phase12/alk5_3hmm_matched50_docking')
FIELDS = ['pair_id', 'class', 'smiles', 'source_id', 'score_kcal_mol', 'error', 'num_modes']
KEYS = ['pair_id', 'class', 'smiles', 'source_id']
EXHAUSTIVENESS = 8
NUM_MODES = 3
SEED = 42


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def jload(p: Path):
    return json.loads(p.read_text())


def save_json(p: Path, data):
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name + '.tmp')
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    temp.replace(p)


def read_csv(p: Path):
    if not p.exists():
        return []
    with p.open(newline='') as f:
        return list(csv.DictReader(f))


def save_csv(p: Path, rows, fields=FIELDS):
    p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.name + '.tmp')
    with temp.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    temp.replace(p)


def prep_redocking(root: Path, exhaustiveness=16):
    import phase9_redock
    # Named alternate; never mutate existing EGFR/ALK5 1VJY target mappings.
    phase9_redock.TARGETS['ALK5_3HMM'] = {
        'pdb': '3HMM', 'ligand': '855',
        'source': 'https://www.rcsb.org/structure/3HMM',
    }
    redock = root / REDOCK_REL
    report_file = redock / 'redock_report.json'
    if report_file.exists():
        existing = jload(report_file)
        if existing.get('pdb_id') == '3HMM' and existing.get('ccd_ligand') == '855' \
           and existing.get('passed') and existing.get('pose_reproduced_within_2A'):
            rfile = redock / 'prepared_receptor.pdbqt'
            gfile = redock / 'vina_grid.txt'
            if rfile.exists() and gfile.exists() and all(
                digest(p) == existing.get('files_sha256', {}).get(p.name)
                for p in (rfile, gfile)
            ):
                print('3HMM redocking already passed and inputs unchanged; skipping.')
                print(f"RMSD={existing.get('pose_rmsd_angstrom')} Angstrom")
                return 0
    report = phase9_redock.run_target(
        'ALK5_3HMM', root / 'reports/phase12/redocking', root,
        offline=False, exhaustiveness=exhaustiveness,
    )
    print(f"3HMM: executed={report.get('passed')}, RMSD={report.get('pose_rmsd_angstrom')}, "
          f"pose<2A={report.get('pose_reproduced_within_2A')}, blocker={report.get('blocker')}")
    if not report.get('passed') or not report.get('pose_reproduced_within_2A'):
        print('STOP: redocking RMSD criterion not met. Do not run the 3HMM panel benchmark.', file=sys.stderr)
        return 2
    return 0


def verify(root: Path):
    # Revalidate frozen panel, original source hashes and completed 1VJY baseline.
    import phase11_alk5 as phase11
    rows, original, *_ = phase11.validate(root, exhaustiveness=EXHAUSTIVENESS)
    if len(rows) != 100 or len(original) != 100:
        raise RuntimeError('Phase 11 baseline must contain the same complete 100-row frozen panel.')
    if not all(math.isfinite(float(r['score_kcal_mol'])) for r in original):
        raise RuntimeError('Invalid Phase 11 baseline score')
    report_file = root / REDOCK_REL / 'redock_report.json'
    if not report_file.exists():
        raise FileNotFoundError('Missing 3HMM redocking report: run redock first')
    rep = jload(report_file)
    if rep.get('pdb_id') != '3HMM' or rep.get('ccd_ligand') != '855':
        raise RuntimeError('Wrong PDB or crystallographic ligand for 3HMM')
    if not rep.get('passed') or not rep.get('pose_reproduced_within_2A') \
       or not math.isfinite(float(rep.get('pose_rmsd_angstrom', float('nan')))) \
       or float(rep['pose_rmsd_angstrom']) >= 2:
        raise RuntimeError('3HMM redocking did not pass the <2 Angstrom criterion')
    receptor = root / REDOCK_REL / 'prepared_receptor.pdbqt'
    grid = root / REDOCK_REL / 'vina_grid.txt'
    for p in (receptor, grid):
        if not p.exists() or digest(p) != rep.get('files_sha256', {}).get(p.name):
            raise RuntimeError(f'3HMM prepared input hash mismatch: {p}')
    center, size = rep.get('center'), rep.get('size')
    if not isinstance(center, list) or not isinstance(size, list) or len(center) != 3 or len(size) != 3:
        raise RuntimeError('Missing cocrystal-based 3HMM grid')
    box = tuple(float(x) for x in center + size)
    if not all(map(math.isfinite, box)) or not all(x > 0 for x in box[3:]):
        raise RuntimeError('Bad docking box')
    baseline_path = root / BASELINE_REL
    manifest = {
        'experiment': 'phase12_3HMM_receptor_sensitivity',
        'panel_sha256': digest(root / PANEL_REL),
        'baseline_phase11_scores_sha256': digest(baseline_path),
        '3hmm_receptor_sha256': digest(receptor), '3hmm_grid_sha256': digest(grid),
        'redocking_rmsd_angstrom': float(rep['pose_rmsd_angstrom']),
        'receptor_pdb': '3HMM', 'ligand_ccd': '855',
        'exhaustiveness': EXHAUSTIVENESS, 'num_modes': NUM_MODES,
        'cpu': 1, 'vina_seed': SEED, 'timeout_per_molecule_sec': 360,
        'score_method': 'Docking_App.app.dock.dock_one',
    }
    out = root / OUT_REL
    mpath = out / 'execution_manifest.json'
    if mpath.exists() and jload(mpath) != manifest:
        raise RuntimeError('Frozen inputs or execution parameters changed since 3HMM run started')
    scores = read_csv(out / 'actual_vina_scores.csv')
    allowed = {tuple(r[k] for k in KEYS) for r in rows}
    observed = set()
    for r in scores:
        k = tuple(r.get(c) for c in KEYS)
        if k not in allowed or k in observed:
            raise RuntimeError('Unknown or duplicate 3HMM score row')
        observed.add(k)
        if not math.isfinite(float(r['score_kcal_mol'])):
            raise RuntimeError('Non-finite score in success-only 3HMM results')
    return rows, original, scores, out, receptor, box, manifest


def auc_for(pairs):
    pos = [-a for a, b in pairs]
    neg = [-b for a, b in pairs]
    n = len(pos)
    return sum((x > y) + 0.5 * (x == y) for x in pos for y in neg) / (n*n)


def paired_matrix(scores):
    pairs = {}
    for s in scores:
        row = pairs.setdefault(s['pair_id'], {})
        if s['class'] in row:
            raise ValueError('Duplicate pair/class')
        row[s['class']] = float(s['score_kcal_mol'])
    return pairs


def summarize(root: Path, baseline, observed, panel):
    base = paired_matrix(baseline)
    obs = paired_matrix(observed)
    expected_ids = {r['pair_id'] for r in panel}
    counts = {x: sum(r['class'] == x for r in observed) for x in ('active', 'decoy')}
    complete = len(observed) == len(panel) and set(obs) == expected_ids and all(
        set(d) == {'active', 'decoy'} for d in obs.values()
    )
    result = {
        'target': 'ALK5', 'receptor_pdb': '3HMM',
        'reference_receptor_pdb': '1VJY', 'panel_sha256': digest(root / PANEL_REL),
        'total': len(panel), 'completed_with_finite_scores': len(observed),
        'valid_docking_counts': counts, 'primary_result_ready': complete,
        'roc_auc_3HMM': None, 'roc_auc_1VJY': None,
        'delta_auc_3HMM_minus_1VJY': None, 'conditional_pair_bootstrap95_delta': None,
        '3hmm_within_pair_active_better': None,
        'comparative_caveat': 'Same compounds and scoring protocol but different crystallographic receptor and ligand-defined grid.',
        'created_utc': now(),
    }
    if complete:
        ids = sorted(expected_ids)
        pnew = [(obs[k]['active'], obs[k]['decoy']) for k in ids]
        pold = [(base[k]['active'], base[k]['decoy']) for k in ids]
        result['roc_auc_3HMM'] = round(auc_for(pnew), 4)
        result['roc_auc_1VJY'] = round(auc_for(pold), 4)
        result['delta_auc_3HMM_minus_1VJY'] = round(result['roc_auc_3HMM'] - result['roc_auc_1VJY'], 4)
        result['3hmm_within_pair_active_better'] = sum(a < b for a,b in pnew)
        rng = random.Random(20261110)
        diffs, aucs = [], []
        for _ in range(10000):
            ids_idx = rng.choices(range(len(ids)), k=len(ids))
            new_sample = [pnew[i] for i in ids_idx]
            old_sample = [pold[i] for i in ids_idx]
            aucs.append(auc_for(new_sample))
            diffs.append(aucs[-1] - auc_for(old_sample))
        diffs.sort(); aucs.sort()
        result['conditional_pair_bootstrap95_delta'] = [round(diffs[250], 4), round(diffs[9750], 4)]
        result['conditional_pair_bootstrap95_3HMM'] = [round(aucs[250], 4), round(aucs[9750], 4)]
        result['summary'] = 'Receptor-structure sensitivity test; not external experimental potency validation.'
    save_json(root / OUT_REL / 'analysis.json', result)
    print(json.dumps(result, indent=2))
    return result


def dock(root: Path, limit: int):
    rows, baseline, existing, out, receptor, box, manifest = verify(root)
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'phase12_single_process.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Phase 12 already running in another terminal. Do not start duplicates.')
        rows, baseline, existing, out, receptor, box, manifest = verify(root)
        if not (out / 'execution_manifest.json').exists():
            save_json(out / 'execution_manifest.json', manifest)
        observed = {tuple(r[k] for k in KEYS) for r in existing}
        todo = [r for r in rows if tuple(r[k] for k in KEYS) not in observed]
        if limit:
            todo = todo[:limit]
        interrupted = {'yes': False}
        old_handler = signal.getsignal(signal.SIGINT)
        def on_sigint(*_):
            interrupted['yes'] = True
            print('\nCtrl+C requested: saving completed molecules, stopping at safe checkpoint...', flush=True)
        signal.signal(signal.SIGINT, on_sigint)
        try:
            from Docking_App.app.dock import dock_one
            for r in todo:
                if interrupted['yes']:
                    break
                try:
                    outc = dock_one(r['smiles'], 'ALK5', exhaustiveness=EXHAUSTIVENESS,
                                    num_modes=NUM_MODES, cpu=1, timeout=360,
                                    pocket_box=box, custom_receptor_pdbqt=receptor)
                    score = float(outc.affinity) if outc.affinity is not None else None
                    error = outc.error
                    modes = outc.num_modes
                except KeyboardInterrupt:
                    interrupted['yes'] = True
                    break
                except Exception as exc:
                    score = None; error = f'{type(exc).__name__}: {str(exc)[:300]}'; modes = 0
                if interrupted['yes']:
                    break
                if score is None or not math.isfinite(score):
                    failpath = out / 'failed_attempts.csv'
                    failures = read_csv(failpath)
                    failures.append({'time_utc': now(), **r, 'error': str(error or 'no score')[:300]})
                    save_csv(failpath, failures, ['time_utc'] + KEYS + ['error'])
                    print(f"FAILED {r['pair_id']} {r['class']}: {error} (pending)", flush=True)
                    continue
                existing.append({**r, 'score_kcal_mol': str(score), 'error': '', 'num_modes': str(modes)})
                save_csv(out / 'actual_vina_scores.csv', existing)
                print(f"{len(existing)}/{len(rows)} {r['pair_id']} {r['class']} score={score:.3f}", flush=True)
        finally:
            signal.signal(signal.SIGINT, old_handler)
        result = summarize(root, baseline, existing, rows)
        if interrupted['yes']:
            print('INTERRUPTED. Rerun dock to resume without losing valid scores.')
        elif not result['primary_result_ready']:
            print('INCOMPLETE. Rerun dock to retry remaining molecules.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path.cwd())
    sub = p.add_subparsers(dest='cmd', required=True)
    sub.add_parser('redock', help='Crystallographic 3HMM/ligand 855 redocking')
    sub.add_parser('check', help='Verify existing Phase 11 panel and validated Phase 12 inputs')
    run = sub.add_parser('dock', help='Dock same 100 molecules against 3HMM, checkpointing successes')
    run.add_argument('--limit', type=int, default=0)
    sub.add_parser('report', help='Recalculate 3HMM vs 1VJY AUC after scoring')
    args = p.parse_args()
    root = args.root.resolve()
    try:
        if args.cmd == 'redock':
            return prep_redocking(root)
        rows, baseline, observed, out, receptor, box, manifest = verify(root)
        if args.cmd == 'check':
            print(json.dumps({
                'ready': True, 'receptor_pdb': '3HMM',
                'redocking_rmsd_angstrom': manifest['redocking_rmsd_angstrom'],
                'frozen_panel_sha256': manifest['panel_sha256'],
                'matched_panel_rows': len(rows), '3hmm_completed': len(observed),
                '3hmm_remaining': len(rows) - len(observed),
                'phase11_baseline_complete': len(baseline) == 100,
            }, indent=2))
        elif args.cmd == 'dock':
            if args.limit < 0:
                raise ValueError('--limit cannot be negative')
            dock(root, args.limit)
        elif args.cmd == 'report':
            summarize(root, baseline, observed, rows)
    except Exception as exc:
        print(f'ERROR {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
