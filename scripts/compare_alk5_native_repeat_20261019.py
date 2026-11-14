"""Describe actual native-repeat differences without replacing either run."""
import json,sys,hashlib
from pathlib import Path
old=Path(sys.argv[1]); new=Path(sys.argv[2])
rows=[]
historical={(r['pdb'],r['function'],r['id'],r['seed']):r for r in json.loads((old/'summary.json').read_text())['results']}
drift=[]
for folder in sorted((new/'jobs').glob('*_native_*')):
    cached=json.loads((old/'jobs'/folder.name/'result.json').read_text()); b=json.loads((folder/'result.json').read_text())
    a=historical[b['pdb'],b['function'],b['id'],b['seed']]
    if any(cached.get(k)!=a.get(k) for k in ['score','native_rmsd_A','native_pose_pass']):
        drift.append(dict(job=folder.name,warning='Historical aggregate differs from current historical job-cache file; aggregate retained, no overwrite',historical={k:a.get(k) for k in ['score','native_rmsd_A','native_pose_pass']},current_cache={k:cached.get(k) for k in ['score','native_rmsd_A','native_pose_pass']}))
    assert a['ligand_sha256']==b['ligand_sha256']
    for f in ['receptor.pdbqt','box.json','native_ref.sdf']:
        assert hashlib.sha256((old/'receptors'/b['pdb']/f).read_bytes()).digest()==hashlib.sha256((new/'receptors'/b['pdb']/f).read_bytes()).digest()
    rows.append(dict(pdb=b['pdb'],function=b['function'],seed=b['seed'],old_score=a['score'],new_score=b['score'],score_delta=b['score']-a['score'],old_native_rmsd_A=a['native_rmsd_A'],new_native_rmsd_A=b['native_rmsd_A'],same_pass_fail=a['native_pose_pass']==b['native_pose_pass']))
result=dict(records=len(rows),historical_source='Frozen summary.json aggregate, not mutable per-job cache',historical_cache_drift=drift,same_input_hashes=True,identical_pass_fail=sum(r['same_pass_fail'] for r in rows),max_abs_score_delta=max(abs(r['score_delta']) for r in rows),max_abs_rmsd_delta_A=max(abs(r['new_native_rmsd_A']-r['old_native_rmsd_A']) for r in rows),warning='Repetition does not imply bitwise reproducibility; differences retained. No cause established by this audit.',records_detail=rows)
(new/'native_repeat_comparison.json').write_text(json.dumps(result,indent=2)); print(json.dumps({k:v for k,v in result.items() if k!='records_detail'},indent=2))
