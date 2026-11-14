#!/usr/bin/env python3
"""Score a pre-frozen EGFR matched panel with REAL AutoDock Vina; resumable and fail-closed.

Inputs: prior Phase 9 validated EGFR 1M17 receptor+grid, and a frozen Phase 10 panel.
Does not sample/replace any molecule, modify scores, or silently ignore failures.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, random
from datetime import datetime, timezone
from pathlib import Path

PANEL_PATH = Path('reports/phase10/egfr_matched50/matched_panel.csv')
AUDIT_PATH = Path('reports/phase10/egfr_matched50/selection_audit.json')
REPORT_DIR = Path('reports/phase10/egfr_matched50_docking')
REPO_REDOCK = Path('reports/phase9/redocking/1M17/redock_report.json')
RESULT_FIELDS = ['pair_id','class','smiles','source_id','score_kcal_mol','error','num_modes']

def sha256(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def readcsv(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))
def auc(pos, neg):
    if not pos or not neg:return None
    return sum((p>n)+0.5*(p==n) for p in pos for n in neg)/(len(pos)*len(neg))
def bootstrap_pairs(pairs, reps=2000, seed=20261110):
    rng = random.Random(seed)
    estimates=[]
    for _ in range(reps):
        draw=rng.choices(pairs,k=len(pairs))
        estimates.append(auc([x[0] for x in draw],[x[1] for x in draw]))
    estimates.sort()
    return [round(estimates[int(.025*reps)],4),round(estimates[min(reps-1,int(.975*reps))],4)]
def persist_scores(path,rows):
    tmp=path.with_suffix('.csv.tmp')
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=RESULT_FIELDS)
        w.writeheader();w.writerows(rows)
        f.flush(); os.fsync(f.fileno())
    os.replace(tmp,path)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path.cwd())
    parser.add_argument('--exhaustiveness',type=int,default=8)
    parser.add_argument('--limit',type=int,default=0,help='Maximum total scored panel rows this invocation (0: all)')
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args();root=args.root.resolve()
    panel=root/PANEL_PATH;audit=root/AUDIT_PATH
    if not panel.exists() or not audit.exists():parser.error('Phase 10 matched panel or audit missing')
    provenance=json.loads(audit.read_text())
    if sha256(panel)!=provenance.get('panel_sha256'):
        parser.error('Frozen panel SHA-256 mismatch. Refusing to run altered or reordered panel.')
    source_dir=root/'reports/phase9/egfr_benchmark/source'
    for klass in ('active','decoy'):
        source=source_dir/f'{klass}.ism'
        if source.exists() and sha256(source)!=provenance['source_sha256'][klass]:
            parser.error(f'{klass} source file hash changed since panel selection')
    panel_rows=readcsv(panel)
    keys=[(r['pair_id'],r['class'],r['smiles']) for r in panel_rows]
    if len(keys)!=len(set(keys)) or len(panel_rows)!=2*provenance['n_per_class']:
        parser.error('Duplicate or missing panel records')
    pairs={}
    for r in panel_rows:
        if r['class'] not in ('active','decoy'):parser.error('Invalid class label')
        pairs.setdefault(r['pair_id'],set()).add(r['class'])
    if len(pairs)!=provenance['n_per_class'] or any(x!={'active','decoy'} for x in pairs.values()):
        parser.error('Panel does not contain exactly 1 active and 1 decoy per pair')
    redock=root/REPO_REDOCK
    if not redock.exists():parser.error('No Phase 9 crystal redocking report')
    rdr=json.loads(redock.read_text())
    if not rdr.get('passed') or not rdr.get('pose_reproduced_within_2A'):
        parser.error('Crystal redocking RMSD acceptance criterion not passed')
    receptor=redock.parent/'prepared_receptor.pdbqt'
    grid=redock.parent/'vina_grid.txt'
    for f in (receptor,grid):
        if not f.exists() or sha256(f)!=rdr.get('files_sha256',{}).get(f.name):
            parser.error(f'Missing or changed Phase 9 validated input {f}')
    out=root/REPORT_DIR;out.mkdir(parents=True,exist_ok=True)
    scores_path=out/'actual_vina_scores.csv'
    manifest_path=out/'execution_manifest.json'
    manifest={'panel_sha256':sha256(panel),'receptor_sha256':sha256(receptor),
              'grid_sha256':sha256(grid),'exhaustiveness':args.exhaustiveness,
              'num_modes':3,'cpu':1,'vina_seed':42,'timeout_per_molecule_sec':360,
              'method':'Docking_App.app.dock.dock_one'}
    if manifest_path.exists():
        old=json.loads(manifest_path.read_text())
        if old!=manifest:parser.error('Existing results use different panel/receptor/parameters. Refusing to combine runs.')
    else:manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    previous=readcsv(scores_path) if scores_path.exists() else []
    bykey={(r['pair_id'],r['class'],r['smiles']):r for r in previous}
    if len(bykey)!=len(previous) or not set(bykey).issubset(set(keys)):
        parser.error('Existing score records contain duplicate/unknown molecules')
    if args.check_only:
        print(json.dumps({'ready':True,'total':len(panel_rows),'already_scored':len(previous),
                          'remaining':len(panel_rows)-len(previous),'panel_sha256':sha256(panel)},indent=2));return 0
    if args.exhaustiveness<1:parser.error('Invalid exhaustiveness')
    from Docking_App.app.dock import dock_one
    completed=list(previous)
    candidates=[r for r in panel_rows if (r['pair_id'],r['class'],r['smiles']) not in bykey]
    if args.limit>0:candidates=candidates[:args.limit]
    box=tuple(rdr['center']+rdr['size'])
    for r in candidates:
        try:
            outcome=dock_one(r['smiles'],'EGFR',exhaustiveness=args.exhaustiveness,
                             num_modes=3,cpu=1,timeout=360,pocket_box=box,
                             custom_receptor_pdbqt=receptor)
            v=outcome.affinity
            if v is not None and not math.isfinite(float(v)):
                v=None
            row={**r,'score_kcal_mol':v if v is not None else '',
                 'error':outcome.error or ('' if v is not None else 'no_finite_affinity'),
                 'num_modes':outcome.num_modes}
        except Exception as e:
            row={**r,'score_kcal_mol':'','error':f'{type(e).__name__}: {str(e)[:250]}','num_modes':0}
        completed.append(row);persist_scores(scores_path,completed)
        print(f"{len(completed)}/{len(panel_rows)} {r['pair_id']} {r['class']}  "
              f"score={row['score_kcal_mol'] or 'NA'}  {row['error']}",flush=True)
    # Always save status; only report primary AUC for full panel and adequate coverage.
    present={ (r['pair_id'],r['class']):r for r in completed }
    valid={c:0 for c in ('active','decoy')}
    vals={c:[] for c in ('active','decoy')}
    paired=[]
    for r in completed:
        if r['score_kcal_mol'] not in ('',None):
            score=-float(r['score_kcal_mol'])
            valid[r['class']]+=1;vals[r['class']].append(score)
    for pair in sorted(pairs):
        aa=present.get((pair,'active'));dd=present.get((pair,'decoy'))
        if aa and dd and aa['score_kcal_mol'] not in ('',None) and dd['score_kcal_mol'] not in ('',None):
            paired.append((-float(aa['score_kcal_mol']),-float(dd['score_kcal_mol'])))
    full=len(completed)==len(panel_rows)
    coverage_ok=all(valid[c]/provenance['n_per_class']>=0.9 for c in valid)
    success=full and coverage_ok
    result={'target':'EGFR','dataset':'DUD-E EGFR','validation':'matched descriptor pilot',
            'panel_sha256':sha256(panel),'docking_completed_count':len(completed),
            'total_panel_rows':len(panel_rows),'valid_docking_counts':valid,
            'matched_pairs_complete':len(paired),
            'primary_result_ready':success,
            'roc_auc':round(auc(vals['active'],vals['decoy']),4) if success else None,
            'conditional_pair_bootstrap95':bootstrap_pairs(paired) if success and len(paired)>3 else None,
            'within_matched_pair_fraction_active_better':round(sum(a>d for a,d in paired)/len(paired),4) if paired else None,
            'limitations':['DUD-E decoys are presumed, not experimentally confirmed, inactive',
                           'AUC is measured on matched subset and not general EGFR inhibitor populations',
                           'Docking scores are not experimental affinities'],
            'created_utc':datetime.now(timezone.utc).isoformat()}
    (out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'processed':len(completed),'total':len(panel_rows),
                      'valid':valid,'auc_if_complete':result['roc_auc'],
                      'matched_pair_success_fraction':result['within_matched_pair_fraction_active_better']},indent=2))
    return 0 if not full or success else 2
if __name__=='__main__':raise SystemExit(main())
