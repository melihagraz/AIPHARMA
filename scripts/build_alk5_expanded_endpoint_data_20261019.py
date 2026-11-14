"""Endpoint-based exact/censored dataset; no model scores or midpoint imputation."""
import hashlib,json,math,statistics,sys
from pathlib import Path
from collections import defaultdict,Counter
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
p=Path(sys.argv[1]); out=p/'curated_complete_identity'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
old=json.loads((p/'prior_assay_classification.json').read_text())
eligible={r['assay_id']:'previous_explicit_catalytic' for r in old if r['endpoint']=='enzyme_catalytic' and r['mutation']=='no_explicit_mutant'}
added={
 'CHEMBL3998849':'Human ALK5, casein substrate, 10uM ATP, radiometric kinase assay; database A code does not change explicit endpoint',
 'CHEMBL4406750':'Human ALK5, dephosphorylated casein, gamma33P ATP, scintillation',
 'CHEMBL4406752':'Human ALK5, casein substrate, gamma33P ATP, radioactive assay',
 'CHEMBL4626180':'Human ALK5, casein substrate, gamma33P ATP, radiometric hotspot assay',
 'CHEMBL5620859':'Human ALK5, casein substrate, gamma33P ATP, 120min enzyme assay'}
eligible.update(added)
save('protocol.json',dict(endpoint='Direct catalytic IC50 only, confidence9 targetCHEMBL4439, no explicit mutant; database A/B assay code accepted only with explicit catalytic description',
    rationale='Four explicit radiometric catalytic assays were excluded solely by prior B-code requirement. Addition defined before modeling, not after model outcomes.',
    eligible_assays=eligible,newly_reviewed_assays=added,
    numeric='Positive finite standard nM IC50, = or explicit inequality; pChEMBL not required',
    identity='Same prior normalization; exclude unresolved stereochemistry and incomplete tautomer enumeration; per assay/normalized structure grouping',
    exact='Median exact pIC50; retain all original measurements and bound constraints',
    censored='Keep pIC50 interval, never substitute boundary or midpoint as an exact response',
    definite_labels='Active if IC50<=100nM guaranteed, weak if IC50>=10000nM guaranteed, otherwise intermediate/uncertain',
    partition='No modeling until source/scaffold components and holdout allocation separately locked',
    limits='Not new experimental measurements, not wild-type confirmation from absent mutant annotation, not independently blind',
    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
raw=json.loads((p/'chembl_raw/all_activities.json').read_text())['activities']
assays={r['assay_chembl_id']:r for r in json.loads((p/'chembl_raw/all_assays.json').read_text())['assays']}
enum=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger()
cache={}; rejected=[]; grouped=defaultdict(list)
for r in raw:
    aid=r['assay_chembl_id']; a=assays[aid]; reasons=[]
    if aid not in eligible: reasons.append('endpoint_not_eligible')
    if a.get('target_chembl_id')!='CHEMBL4439' or a.get('confidence_score')!=9: reasons.append('target_confidence')
    if r.get('data_validity_comment') or r.get('potential_duplicate'): reasons.append('flagged_validity_or_duplicate')
    if r.get('standard_type')!='IC50' or r.get('standard_units')!='nM' or r.get('standard_relation') not in ['=','<','<=','>','>=']: reasons.append('unsupported_numeric_relation_units')
    try:
        value=float(r['standard_value']); assert math.isfinite(value) and value>0
    except Exception: reasons.append('nonpositive_or_missing_value'); value=None
    if reasons:
        rejected.append(dict(activity_id=r['activity_id'],assay_id=aid,reasons=reasons)); continue
    try:
        text=r['canonical_smiles']
        if text not in cache:
            mol=Chem.MolFromSmiles(text); assert mol is not None
            assert not any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(mol))
            standardized=un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol)))
            enumeration=enum.Enumerate(standardized)
            assert str(enumeration.status)=='Completed', 'Incomplete tautomer enumeration'
            normalized=enum.PickCanonical(enumeration)
            Chem.RemoveStereochemistry(normalized)
            cache[text]=(Chem.MolToSmiles(normalized),Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(normalized)))
        identity,scaffold=cache[text]
    except Exception as exc:
        reason='incomplete_tautomer_enumeration' if 'Incomplete tautomer' in str(exc) else 'structure_or_unresolved_stereo'
        rejected.append(dict(activity_id=r['activity_id'],assay_id=aid,reasons=[reason])); continue
    grouped[aid,identity].append(dict(activity_id=r['activity_id'],document_id=r['document_chembl_id'],molecule_id=r.get('molecule_chembl_id'),smiles=text,scaffold=scaffold,relation=r['standard_relation'],ic50_nM=value,pic50_boundary=9-math.log10(value)))
records=[]; conflicts=[]
for (aid,identity),rs in sorted(grouped.items()):
    exact=[r['pic50_boundary'] for r in rs if r['relation']=='=']
    lower=max([r['pic50_boundary'] for r in rs if r['relation'] in ['<','<=']],default=-math.inf)
    upper=min([r['pic50_boundary'] for r in rs if r['relation'] in ['>','>=']],default=math.inf)
    median=statistics.median(exact) if exact else None
    if exact and min(exact)<=5 and max(exact)>=7:
        conflicts.append(dict(assay_id=aid,normalized_parent=identity,measurements=rs,reason='exact_activity_cliff_conflict_within_same_assay_identity')); continue
    if lower>upper or (median is not None and not lower-1e-9<=median<=upper+1e-9):
        conflicts.append(dict(assay_id=aid,normalized_parent=identity,measurements=rs,reason='inconsistent_exact_and_interval_constraints')); continue
    lo=median if median is not None else lower; hi=median if median is not None else upper
    label='active' if lo>=7 else ('weak' if hi<=5 else ('intermediate' if median is not None else 'uncertain_interval'))
    records.append(dict(assay_id=aid,normalized_parent=identity,scaffold=rs[0]['scaffold'],smiles=rs[0]['smiles'],documents=sorted({r['document_id'] for r in rs}),median_pic50=median,lower_pic50=None if not math.isfinite(lo) else lo,upper_pic50=None if not math.isfinite(hi) else hi,label=label,measurement_kind='exact_median' if exact else 'censored_interval',newly_reviewed_endpoint=aid in added,measurements=rs))
save('records.json',records); save('exclusions.json',rejected); save('constraint_conflicts.json',conflicts)
summary=dict(raw_records=len(raw),eligible_assays=len(eligible),represented_assays=len({r['assay_id'] for r in records}),assay_structure_groups=len(records),unique_normalized_structures=len({r['normalized_parent'] for r in records}),exact_groups=sum(r['measurement_kind']=='exact_median' for r in records),interval_groups=sum(r['measurement_kind']=='censored_interval' for r in records),labels=dict(Counter(r['label'] for r in records)),newly_reviewed_endpoint_groups=sum(r['newly_reviewed_endpoint'] for r in records),newly_reviewed_endpoint_labels=dict(Counter(r['label'] for r in records if r['newly_reviewed_endpoint'])),constraint_conflict_groups=len(conflicts),exclusion_reasons=dict(Counter(reason for r in rejected for reason in r['reasons'])),modeling_started=False)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
