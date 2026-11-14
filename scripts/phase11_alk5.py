#!/usr/bin/env python3
"""Phase 11: DUD-E TGFR1/ALK5 frozen matched panel and resumable real-Vina docking.

Three commands from AIPHARMA root:
  python scripts/phase11_alk5.py prepare
  python scripts/phase11_alk5.py check
  PYTHONPATH="$PWD" python scripts/phase11_alk5.py dock --limit 10

Default receptor: ALK5 PDB 1VJY, verified by Phase 9 redocking.
This is a computational benchmark (presumed DUD-E decoys), not experimental activity validation.
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
import random
import signal
import sys

SOURCE = Path('reports/phase11/alk5_benchmark/source')
PANEL_DIR = Path('reports/phase11/alk5_matched50')
OUT_DIR = Path('reports/phase11/alk5_matched50_docking')
REDOCK_DIR = Path('reports/phase9/redocking/1VJY')
EXPECTED_SOURCE_HASHES = {
    'active': '6a1df16bfa5477c932057e574146c2d58c5854768525feb900aa8042c0930829',
    'decoy':  '503596358d8d72ed29f22595c8d1c5ab702768f4a30e9d1502a43c47130d8440',
}
FEATURES = ['MW', 'logP', 'HBA', 'HBD', 'TPSA', 'RotB', 'charge']
# Calipers and scales are fixed independently of docking outcomes.
CALIPERS = [65.0, 1.0, 2.0, 1.0, 40.0, 4.0]
SCALES = [50.0, 0.8, 2.0, 1.0, 30.0, 3.0]
N_PAIRS = 50
SEED = 20261110
FIELDS = ['pair_id', 'class', 'smiles', 'source_id', 'score_kcal_mol', 'error', 'num_modes']


def utc(): return datetime.now(timezone.utc).isoformat()
def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save_atomic(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    with temp.open('w', encoding='utf-8', newline='') as f:
        f.write(content)
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)

def json_read(path): return json.loads(path.read_text(encoding='utf-8'))
def json_save(path, obj): save_atomic(path, json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def csv_rows(path):
    if not path.exists(): return []
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def csv_save(path, rows, fields):
    from io import StringIO
    b = StringIO()
    w = csv.DictWriter(b, fieldnames=fields, lineterminator='\n', extrasaction='ignore')
    w.writeheader(); w.writerows(rows)
    save_atomic(path, b.getvalue())

def verify_sources(root):
    hashes = {}
    for klass in ('active', 'decoy'):
        f = root / SOURCE / (klass + '.ism')
        if not f.is_file(): raise FileNotFoundError(f'Missing TGFR1 source file: {f}')
        actual = sha256(f)
        if actual != EXPECTED_SOURCE_HASHES[klass]:
            raise RuntimeError(f'{klass}.ism SHA256 changed! expected {EXPECTED_SOURCE_HASHES[klass]}, observed {actual}. Stop and inspect source provenance.')
        hashes[klass] = actual
    return hashes

def chem_imports():
    try:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
        from rdkit import Chem, rdBase
        from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
    except ImportError as e:
        raise RuntimeError('Requires RDKit, numpy and scipy in the active WSL environment') from e
    return np, linear_sum_assignment, Chem, rdBase, Descriptors, Crippen, Lipinski, rdMolDescriptors

def parse_ism(path, chem):
    Chem, Descriptors, Crippen, Lipinski, rdMolDescriptors = chem
    parsed = {}; invalid = 0; duplicates = 0; lines = 0
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        tokens = line.split()
        if not tokens or tokens[0].startswith('#'): continue
        lines += 1
        mol = Chem.MolFromSmiles(tokens[0])
        if mol is None:
            invalid += 1; continue
        key = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        if key in parsed:
            duplicates += 1; continue
        vec = [float(Descriptors.MolWt(mol)), float(Crippen.MolLogP(mol)),
               float(Lipinski.NumHAcceptors(mol)), float(Lipinski.NumHDonors(mol)),
               float(rdMolDescriptors.CalcTPSA(mol)), float(Lipinski.NumRotatableBonds(mol)),
               float(Chem.GetFormalCharge(mol))]
        if any(not math.isfinite(x) for x in vec):
            invalid += 1; continue
        parsed[key] = {'smiles': key, 'source_id': tokens[1] if len(tokens)>1 else '', 'descriptor': vec}
    return parsed, {'input_rows': lines, 'unique_valid_smiles':len(parsed), 'duplicate_smiles':duplicates, 'invalid_smiles':invalid}

def prepare(root, allow_fixture=False, n_pairs=N_PAIRS, seed=SEED):
    np, assignment, Chem, rdBase, Descriptors, Crippen, Lipinski, rdMolDescriptors = chem_imports()
    hashes = verify_sources(root) if not allow_fixture else {cl:sha256(root/SOURCE/(cl+'.ism')) for cl in ('active','decoy')}
    chem = Chem, Descriptors, Crippen, Lipinski, rdMolDescriptors
    act, acount = parse_ism(root/SOURCE/'active.ism', chem)
    dec, dcount = parse_ism(root/SOURCE/'decoy.ism', chem)
    overlap = set(act).intersection(dec)
    for key in overlap: dec.pop(key)
    if len(act) < n_pairs or len(dec) < n_pairs:
        raise ValueError(f'Not enough unique compounds for {n_pairs} pairs (actives={len(act)}, decoys={len(dec)})')
    out = root / PANEL_DIR
    if (out/'matched_panel.csv').exists() or (out/'selection_audit.json').exists():
        if not (out/'matched_panel.csv').exists() or not (out/'selection_audit.json').exists():
            raise RuntimeError('Existing panel incomplete: stop, do not overwrite')
        audit = json_read(out/'selection_audit.json')
        if (sha256(out/'matched_panel.csv') != audit.get('panel_sha256') or
            audit.get('source_sha256') != hashes or audit.get('n_per_class') != n_pairs):
            raise RuntimeError('Existing ALK5 frozen panel differs from current inputs; refusing to overwrite')
        print(json.dumps({'already_prepared':True,'pairs':audit['n_per_class'],
                          'max_abs_smd':audit.get('max_abs_smd'),
                          'panel_sha256':audit['panel_sha256']},indent=2))
        return audit
    rng = random.Random(seed)
    actives = [act[k] for k in sorted(act)]
    decoys = [dec[k] for k in sorted(dec)]
    D = np.asarray([r['descriptor'] for r in decoys], dtype=np.float64)
    A = np.asarray([r['descriptor'] for r in actives], dtype=np.float64)
    # Pre-specified strict calipers; no docking information enters selection.
    feasible = (A[:,None,6] == D[None,:,6])
    # Construct iteratively to avoid large temporary multidimensional arrays.
    for j, c in enumerate(CALIPERS):
        feasible &= np.abs(A[:,None,j] - D[None,:,j]) <= c
    eligible = [i for i in range(len(actives)) if np.count_nonzero(feasible[i]) >= 1]
    rng.shuffle(eligible)
    if len(eligible) < n_pairs:
        raise RuntimeError(f'Only {len(eligible)} of {len(actives)} TGFR1 actives have a decoy within the pre-specified calipers; cannot select {n_pairs} pairs without changing protocol.')
    picked = None
    attempt = 0
    for attempt in range(1,65):
        if attempt == 1:
            chosen = eligible[:n_pairs]
        else:
            chosen = rng.sample(eligible, n_pairs)
        mat = np.empty((n_pairs, len(decoys)), dtype=np.float64)
        for k, idx in enumerate(chosen):
            diff = (D[:,:6] - A[idx,:6]) / np.asarray(SCALES)
            sq = np.einsum('ij,ij->i', diff, diff)
            mat[k,:] = np.where(feasible[idx], sq, 1e9)
        rows, cols = assignment(mat)
        if len(rows)==n_pairs and np.all(mat[rows, cols] < 1e9):
            picked = (chosen, rows, cols, mat)
            break
    if picked is None:
        raise RuntimeError('No 1:1 feasible assignment in 64 pre-specified descriptor-only attempts. Do not relax calipers without revising protocol.')
    chosen, rows, cols, mat = picked
    panel = []; descriptors=[]; aa=[]; dd=[]
    for idx, (i, j) in enumerate(zip(rows, cols), start=1):
        ar = actives[chosen[int(i)]]; dr = decoys[int(j)]
        aa.append(ar['descriptor']); dd.append(dr['descriptor'])
        for klass, obj in [('active',ar),('decoy',dr)]:
            pair_id = f'pair_{idx:03d}'
            panel.append({'pair_id':pair_id,'class':klass,'smiles':obj['smiles'],'source_id':obj['source_id']})
            descriptors.append({'pair_id':pair_id,'class':klass,'smiles':obj['smiles'],
                                **{FEATURES[k]: round(obj['descriptor'][k],4) for k in range(7)}})
    aa=np.asarray(aa); dd=np.asarray(dd)
    balance={}
    for j, name in enumerate(FEATURES):
        pooled=np.sqrt((aa[:,j].var(ddof=1) + dd[:,j].var(ddof=1))/2)
        mean_diff=float(aa[:,j].mean()-dd[:,j].mean())
        balance[name]={'active_mean':round(float(aa[:,j].mean()),4),
                       'decoy_mean':round(float(dd[:,j].mean()),4),
                       'smd':round(mean_diff/pooled,4) if pooled>0 else 0.0,
                       'max_pair_difference':round(float(abs(aa[:,j]-dd[:,j]).max()),4)}
    max_abs_smd = max(abs(v['smd']) for v in balance.values())
    panel_file=out/'matched_panel.csv'
    csv_save(panel_file,panel,['pair_id','class','smiles','source_id'])
    csv_save(out/'matched_descriptors.csv',descriptors,['pair_id','class','smiles']+FEATURES)
    audit={
        'target':'ALK5','dud_e_target':'TGFR1','source_dataset':'DUD-E TGFR1 via THU-ATOM DrugCLIP mirror',
        'benchmark_receptor_pdb':'1VJY','dud_e_reference_pdb':'3HMM',
        'n_per_class':n_pairs,'seed':seed,'matching_attempt':attempt,
        'source_sha256':hashes,'source_parse_counts':{'active':acount,'decoy':dcount},
        'active_decoy_canonical_smiles_overlap':len(overlap),
        'matching_features':FEATURES,'calipers_MW_logP_HBA_HBD_TPSA_RotB':CALIPERS,
        'scales_MW_logP_HBA_HBD_TPSA_RotB':SCALES,
        'exact_formal_charge_matching':True,
        'matching_method':'seeded eligible active sampling followed by one-to-one minimum-cost assignment without replacement',
        'selection_uses_docking_scores':False,'descriptor_balance':balance,
        'max_abs_smd':round(max_abs_smd,4),'balance_below_0_1':bool(max_abs_smd<0.1),
        'panel_sha256':sha256(panel_file),
        'rdkit_version':rdBase.rdkitVersion,'numpy_version':np.__version__,
        'caveats':['DUD-E 235 raw TGFR1 ligands vs 133 clustered in actives_final.ism',
                   'DUD-E decoys presumed inactive; not an experimentally confirmed inactive set',
                   '1VJY is not DUD-E TGFR1 benchmark reference 3HMM',
                   'Descriptor matching cannot remove all scaffold/analogue/decoy bias'],
        'created_utc':utc(),
    }
    json_save(out/'selection_audit.json',audit)
    print(json.dumps({'prepared':True,'matched_pairs':n_pairs,'eligible_actives':len(eligible),
                      'selection_attempt':attempt,'max_abs_smd':round(max_abs_smd,4),
                      'balance_below_0_1':bool(max_abs_smd<0.1),
                      'panel_sha256':audit['panel_sha256'],
                      'next':'python scripts/phase11_alk5.py check'},indent=2))
    return audit

def validate(root, exhaustiveness=8):
    source_hashes = verify_sources(root)
    p = root/PANEL_DIR/'matched_panel.csv'
    audit_file = root/PANEL_DIR/'selection_audit.json'
    if not p.exists() or not audit_file.exists():
        raise RuntimeError('No frozen ALK5 panel: run prepare first')
    audit = json_read(audit_file)
    if audit.get('source_sha256')!=source_hashes or sha256(p)!=audit.get('panel_sha256'):
        raise RuntimeError('Source or frozen panel SHA-256 mismatch')
    rows=csv_rows(p)
    expected=n= audit['n_per_class']
    if n != N_PAIRS or len(rows)!=2*n:
        raise RuntimeError(f'Unexpected panel size ({len(rows)} vs {2*N_PAIRS})')
    seen=set(); pairs={}
    for r in rows:
        key=(r['pair_id'],r['class'],r['smiles'])
        if key in seen: raise RuntimeError('Duplicate panel molecule')
        seen.add(key)
        if r['class'] not in ('active','decoy'):raise RuntimeError('Invalid class label')
        pairs.setdefault(r['pair_id'],set()).add(r['class'])
    if len(pairs)!=n or any(v != {'active','decoy'} for v in pairs.values()):
        raise RuntimeError('Unmatched or duplicate pairs in frozen panel')
    redock=root/REDOCK_DIR/'redock_report.json'
    if not redock.exists():raise FileNotFoundError('Missing Phase 9 ALK5 1VJY redocking report')
    rep=json_read(redock)
    if rep.get('target')!='ALK5' or rep.get('pdb_id')!='1VJY' or not rep.get('passed') or not rep.get('pose_reproduced_within_2A'):
        raise RuntimeError('ALK5 1VJY crystallographic redocking did not pass')
    rmsd=float(rep.get('pose_rmsd_angstrom',float('nan')))
    if not (math.isfinite(rmsd) and 0<=rmsd<2):
        raise RuntimeError('Invalid ALK5 redocking RMSD')
    receptor=root/REDOCK_DIR/'prepared_receptor.pdbqt'; grid=root/REDOCK_DIR/'vina_grid.txt'
    for f in (receptor,grid):
        if not f.exists() or sha256(f)!=rep.get('files_sha256',{}).get(f.name):
            raise RuntimeError(f'Missing or modified validated ALK5 docking input: {f}')
    center=rep.get('center'); size=rep.get('size')
    if not (isinstance(center,list) and isinstance(size,list) and len(center)==3 and len(size)==3):
        raise RuntimeError('Missing ALK5 1VJY validated grid coordinates')
    box=tuple(float(x) for x in center+size)
    if not all(map(math.isfinite, box)) or not all(x>0 for x in box[3:]):
        raise RuntimeError('Invalid validated grid geometry')
    manifest={
        'panel_sha256':sha256(p),'receptor_sha256':sha256(receptor),'grid_sha256':sha256(grid),
        'receptor_pdb':'1VJY','receptor_ligand_ccd':'460','exhaustiveness':exhaustiveness,
        'num_modes':3,'cpu':1,'vina_seed':42,'timeout_per_molecule_sec':360,
        'method':'Docking_App.app.dock.dock_one',
    }
    out=root/OUT_DIR
    existing=out/'execution_manifest.json'
    if existing.exists() and json_read(existing)!=manifest:
        raise RuntimeError('Existing ALK5 results have different panel/receptor/scoring parameters')
    scores=csv_rows(out/'actual_vina_scores.csv')
    known={(r['pair_id'],r['class'],r['smiles']) for r in rows}
    seen_scores=set()
    for r in scores:
        key=(r['pair_id'],r['class'],r['smiles'])
        if key not in known or key in seen_scores:
            raise RuntimeError('Unknown or duplicate scored molecule. Do not merge results.')
        seen_scores.add(key)
        try:
            x=float(r['score_kcal_mol'])
        except (ValueError,KeyError):
            raise RuntimeError('A non-finite/failed record was found in success-only score file; inspect before proceeding')
        if not math.isfinite(x):raise RuntimeError('Non-finite docking score')
    return rows, scores, out, receptor, box, manifest, rmsd


def roc_auc(pos, neg):
    if not pos or not neg: return None
    return sum((x>y)+0.5*(x==y) for x in pos for y in neg)/(len(pos)*len(neg))

def paired_bootstrap(paired, reps=6000, seed=SEED+3):
    rng=random.Random(seed); n=len(paired); vals=[]
    for _ in range(reps):
        sample=rng.choices(paired,k=n)
        vals.append(roc_auc([a for a,b in sample],[b for a,b in sample]))
    vals.sort()
    return [round(vals[int(0.025*reps)],4),round(vals[min(reps-1,int(0.975*reps))],4)]

def report(root, rows, scores, rmsd):
    by_pair={}; counts={'active':0,'decoy':0}; v={'active':[],'decoy':[]}
    for r in scores:
        cls=r['class']; value=-float(r['score_kcal_mol'])
        counts[cls]+=1;v[cls].append(value)
        by_pair.setdefault(r['pair_id'],{})[cls]=value
    pairs=[(x['active'],x['decoy']) for x in by_pair.values() if set(x)=={'active','decoy'}]
    wins=sum(a>d for a,d in pairs); ties=sum(a==d for a,d in pairs)
    full=len(scores)==len(rows)
    coverage_ok=all(counts[k]>=0.9*N_PAIRS for k in counts)
    ready=full and coverage_ok
    result={
        'target':'ALK5','dataset':'DUD-E TGFR1 (clustered actives + property-matched decoys)',
        'receptor_pdb':'1VJY', 'dud_e_reference_pdb':'3HMM',
        'redocking_rmsd_angstrom':rmsd,'panel_sha256':sha256(root/PANEL_DIR/'matched_panel.csv'),
        'completed_with_finite_scores':len(scores),'total':len(rows),'valid_docking_counts':counts,
        'matched_pairs_complete':len(pairs),'primary_result_ready':ready,
        'roc_auc':round(roc_auc(v['active'],v['decoy']),4) if ready else None,
        'conditional_pair_bootstrap95':paired_bootstrap(pairs) if ready else None,
        'pairwise_active_better':wins, 'pairwise_ties':ties,
        'within_matched_pair_fraction_active_better':round((wins+0.5*ties)/len(pairs),4) if pairs else None,
        'limitations':[
            'DUD-E presumed inactive decoys are not experimentally confirmed negatives',
            'Panel is a 50+50 subset; descriptor matching does not remove all chemical/scaffold biases',
            'Validation receptor 1VJY differs from DUD-E native TGFR1 3HMM target structure',
            'Docking score is not an experimental inhibition or binding affinity measurement',
        ],'created_utc':utc(),
    }
    out=root/OUT_DIR
    out.mkdir(parents=True,exist_ok=True)
    json_save(out/'analysis.json',result)
    print(json.dumps({'processed':len(scores),'total':len(rows),
                      'valid':counts,'auc_if_complete':result['roc_auc'],
                      'paired_success_fraction':result['within_matched_pair_fraction_active_better'],
                      'matched_pairs_complete':len(pairs)},indent=2))
    return result

def dock(root, limit, exhaustiveness):
    rows, current, out, receptor, box, manifest, rmsd=validate(root,exhaustiveness)
    out.mkdir(parents=True,exist_ok=True)
    lock_file=out/'phase11_single_process.lock'
    with lock_file.open('a+') as lock:
        try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Phase 11 docking is already running in another terminal. Stop this duplicate invocation.')
        # Race-safe validation under lock.
        rows, current, out, receptor, box, manifest, rmsd=validate(root,exhaustiveness)
        if not (out/'execution_manifest.json').exists():
            json_save(out/'execution_manifest.json', manifest)
        scored={(r['pair_id'],r['class'],r['smiles']) for r in current}
        todo=[r for r in rows if (r['pair_id'],r['class'],r['smiles']) not in scored]
        if limit>0:todo=todo[:limit]
        if not todo:
            print('No remaining rows; generating report from existing valid scores.')
            report(root,rows,current,rmsd)
            return
        # Ctrl+C cannot contaminate the score table with Vina rc=-2 failures.
        interrupted={'flag':False}
        old=signal.getsignal(signal.SIGINT)
        def interrupt_handler(sig, frame):
            interrupted['flag']=True
            print('\nInterrupt requested. Current molecule may fail; previous scores are preserved.',flush=True)
        signal.signal(signal.SIGINT,interrupt_handler)
        try:
            from Docking_App.app.dock import dock_one
            for r in todo:
                if interrupted['flag']: break
                try:
                    outcome=dock_one(r['smiles'],'ALK5',exhaustiveness=exhaustiveness,
                                    num_modes=3,cpu=1,timeout=360,pocket_box=box,
                                    custom_receptor_pdbqt=receptor)
                    value=outcome.affinity
                    err=outcome.error
                    modes=outcome.num_modes
                except KeyboardInterrupt:
                    interrupted['flag']=True; break
                except Exception as e:
                    value=None;err=f'{type(e).__name__}: {str(e)[:300]}';modes=0
                if interrupted['flag']: break
                if value is not None:
                    try:
                        value=float(value)
                        if not math.isfinite(value): value=None
                    except (ValueError,TypeError):value=None
                if value is None:
                    # Never persist NA as a completed score; failed rows remain resumable.
                    fail={'time_utc':utc(),'pair_id':r['pair_id'],'class':r['class'],
                          'smiles':r['smiles'],'source_id':r['source_id'],'error':str(err or 'missing_affinity')[:500]}
                    path=out/'failed_attempts.csv'; past=csv_rows(path);past.append(fail)
                    csv_save(path,past,['time_utc','pair_id','class','smiles','source_id','error'])
                    print(f"FAILED {r['pair_id']} {r['class']} {fail['error']} (will retry on next run)",flush=True)
                    continue
                item={**r,'score_kcal_mol':value,'error':'','num_modes':modes}
                current.append(item)
                csv_save(out/'actual_vina_scores.csv',current,FIELDS)
                print(f"{len(current)}/{len(rows)} {r['pair_id']} {r['class']} score={value}",flush=True)
        finally:
            signal.signal(signal.SIGINT,old)
        result=report(root,rows,current,rmsd)
        if interrupted['flag']: print('INTERRUPTED; valid scores saved, remaining compounds may be resumed safely.')
        if not result['primary_result_ready'] and not interrupted['flag']:
            print('Incomplete: retry with the same dock command; failed molecules remain pending.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd(),help='AIPHARMA project root')
    sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('prepare',help='Build frozen, descriptor-matched 50+50 panel from local .ism files')
    chk=sub.add_parser('check',help='Verify source/panel/redocking/hashes and remaining counts')
    chk.add_argument('--exhaustiveness',type=int,default=8)
    run=sub.add_parser('dock',help='Run actual Vina docking against validated ALK5 1VJY')
    run.add_argument('--limit',type=int,default=0)
    run.add_argument('--exhaustiveness',type=int,default=8)
    a=p.parse_args();root=a.root.resolve()
    try:
        if a.cmd=='prepare': prepare(root)
        elif a.cmd=='check':
            rows, scores, out, rec, box, manifest, rmsd=validate(root,a.exhaustiveness)
            print(json.dumps({'ready':True,'redocking_rmsd_angstrom':rmsd,
                              'panel_rows':len(rows),'successful_scores':len(scores),
                              'remaining':len(rows)-len(scores),
                              'max_abs_smd':json_read(root/PANEL_DIR/'selection_audit.json').get('max_abs_smd'),
                              'panel_sha256':sha256(root/PANEL_DIR/'matched_panel.csv')},indent=2))
        else:
            if a.exhaustiveness<1 or a.exhaustiveness>128: raise ValueError('exhaustiveness must be 1..128')
            dock(root,a.limit,a.exhaustiveness)
    except Exception as e:
        print(f'ERROR: {type(e).__name__}: {e}',file=sys.stderr)
        return 2
    return 0
if __name__=='__main__':raise SystemExit(main())
