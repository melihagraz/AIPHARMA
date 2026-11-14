"""Independent arithmetic and split audit over saved results; no predictor calls."""
import json,math,statistics,sys
from pathlib import Path
from collections import defaultdict
p=Path(sys.argv[1])
def read(f): return json.loads((p/f).read_text(encoding='utf-8-sig'))
records=read('local_sar/measured_records.json'); pred=read('continuous_transfer/predictions.json'); summary=read('continuous_transfer/summary.json')
assignment={(r['assay_id'],r['normalized_parent']):r['component'] for r in pred}
assert len(assignment)==len(records)==541
owners=defaultdict(set); identities=defaultdict(set)
for r in records:
    c=assignment[r['assay_id'],r['normalized_parent']]
    for key in [('identity',r['normalized_parent']),('scaffold',r['scaffold'])]+[('document',m['document_id']) for m in r['measurements'] if m['document_id']]: owners[key].add(c)
    identities[r['normalized_parent']].add(c)
assert all(len(v)==1 for v in owners.values())
for r in pred: assert r['component'] not in identities[r['nearest_normalized_parent']]
components=defaultdict(list)
for r in pred: components[r['component']].append(r)
mae=statistics.mean(statistics.mean(abs(r['actual_pic50']-r['nearest_prediction']) for r in rs) for rs in components.values())
assert abs(mae-summary['component_mean_nearest_MAE'])<1e-12
external=read('external_evaluation/predictions.json'); training=read('external_evaluation/training_structure_medians.json')
assert not ({r['normalized_parent'] for r in external}&{r['normalized_parent'] for r in training})
extsum=read('external_evaluation/summary.json')
for label in ['nearest','median']:
    calc=statistics.mean(abs(r['pic50']-r[label+'_pic50']) for r in external)
    assert abs(calc-extsum['metrics'][label]['MAE'])<1e-12
def passing(r):
    d=r['delta']
    return (d['Solubility_AqSolDB']>=.3 and all(d[k]<=.02 for k in ['DILI','hERG','AMES']) and any(d[k]<=-.05 for k in ['DILI','hERG','AMES']) and d['Bioavailability_Ma']>=-.05 and d['Caco2_Wang']>=-.2)
rounds=[]
for f in sorted(p.glob('development_*.json'))+[p/'reserved_final.json',p/'refinement/results.json',p/'refinement_sar_guided/results.json']:
    rows=json.loads(f.read_text()); assert all(passing(r)==r['engineering_objectives_met'] for r in rows)
    rounds.append(dict(file=str(f.relative_to(p)),count=len(rows),numerical_hits=sum(passing(r) for r in rows)))
for stage in ['refinement','refinement_sar_guided']:
    rows=read(stage+'/results.json'); assert all(not r['activity_validated'] and not r['safety_validated'] for r in rows)
result=dict(status='passed',checks=['541 predictions accounted for','document/identity/scaffold groups never split across components','nearest references always outside held component','saved MAE arithmetic recomputed','external identity overlap absent','all numerical gates independently recomputed','no refined candidate marked activity/safety validated'],rounds=rounds,
    external_duplicate_group_sensitivity=dict(retained=sum(len(r['source_measurements'])==1 for r in external),nearest_MAE=statistics.mean(abs(r['pic50']-r['nearest_pic50']) for r in external if len(r['source_measurements'])==1),median_MAE=statistics.mean(abs(r['pic50']-r['median_pic50']) for r in external if len(r['source_measurements'])==1)),
    limitations='This verifier checks saved split labels/arithmetic, not the biological validity of source measurements or independently the drawn patent structures.')
if (p/'regression_challenger/summary.json').exists():
    challenger=read('regression_challenger/predictions.json'); cs=read('regression_challenger/summary.json')
    assert len(challenger)==541
    cg=defaultdict(list)
    actual={(r['assay_id'],r['normalized_parent']):r['actual_pic50'] for r in pred}
    for r in challenger:
        key=(r['assay_id'],r['normalized_parent'])
        assert r['component']==assignment[key] and abs(r['actual_pic50']-actual[key])<1e-12
        cg[r['component']].append(r)
    cm=statistics.mean(statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in rs) for rs in cg.values())
    assert abs(cm-cs['component_mean_MAE'])<1e-12
    accepted=cm<=.7 and cs['component_mean_gain']>=.1 and cs['fixed_prediction_component_bootstrap_gain_95'][0]>0
    assert accepted==cs['engineering_acceptance_met']==cs['candidate_predictions_permitted_by_protocol']
    result['checks'].append('regression challenger labels/splits/error and rejection gate recomputed')
if (p/'pose_check/geometry_summary.json').exists():
    pose=read('pose_check/summary.json')['results']; diag=read('pose_check/pose_diagnostics.json')
    audit=read('pose_check/native_rmsd_independent_audit.json')
    assert len(pose)==126 and len(diag)==126
    assert all(r['status']=='complete' for r in pose)
    assert sum(r.get('native_pose_pass',False) for r in pose if r['kind']=='native_control')==14
    assert audit['same_pass_fail_decisions']==audit['agreements_within_0_01A']==18
    for r in diag:
        assert r['distance_compatible_pose']==(bool(r['contacts']) and r['min_protein_heavy_atom_distance_A']>=1.6)
        assert all(2.4<=c['distance_A']<=3.5 for c in r['contacts'])
    for r in read('pose_check/geometry_summary.json'):
        actual_groups=0; paired=0
        for group in r['groups']:
            members=[d for d in diag if d['id']==r['id'] and d['pdb']==group['pdb'] and d['function']==group['function']]
            count=sum(d['distance_compatible_pose'] for d in members)
            assert count==group['compatible_seeds'] and len(members)==3
            actual_groups+=bool(group['interpretable'] and count>=2)
            paired+=bool(group['paired_parent_supports_interpretation'] and count>=2)
        assert actual_groups==r['compatible_interpretable_groups']
        assert paired==r['compatible_paired_parent_supported_groups']
    result['checks'].append('126 complete poses, 14/18 native recovery, independent 18/18 RMSD agreement and geometry-group arithmetic verified')
(p/'verification.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
