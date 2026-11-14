"""Audit measured local SAR without assigning activity to an unmeasured design."""
import argparse, hashlib, json, math, statistics
from collections import defaultdict, Counter
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize

ap=argparse.ArgumentParser(); ap.add_argument('root',type=Path); ap.add_argument('output',type=Path)
a=ap.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
def save(n,x): (a.output/n).write_text(json.dumps(x,indent=2,allow_nan=False))
def read(p): return json.loads((a.root/p).read_text(encoding='utf-8-sig'))
classifications=read('reports/alk5_endpoint_audit_20261016/assay_classification.json')
eligible={r['assay_id'] for r in classifications if r['endpoint']=='enzyme_catalytic' and r['mutation']=='no_explicit_mutant'}
save('protocol.json',dict(assays=sorted(eligible),endpoint='IC50 exact nM; catalytic assay, no explicit mutant',
    exclusions=['nonexact or nonpositive IC50','data validity comment','potential duplicate','unresolved stereochemistry'],
    identity='FragmentParent Cleanup Uncharger canonical tautomer, stereo removed for identity only',
    support='Measured analogues, not predictions or independent validation of current designs',
    neighbor_similarity_min=.6, pair_similarity_min=.65, pair_rule='Same exact Murcko scaffold and same assay',
    uncertainty='No explicit mutant annotation is not confirmation of wild type; chemical neighborhood can contain activity cliffs',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
enum=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger()
def identity(m):
    m=enum.Canonicalize(un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))))
    Chem.RemoveStereochemistry(m); return Chem.MolToSmiles(m)
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def scaffold(m): return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
raw={}
for file in ['all_activities.json','weak_boundary_all.json']:
    for r in read('reports/alk5_holdout_20261016/raw/'+file)['activities']: raw[r['activity_id']]=r
excluded=Counter(); groups=defaultdict(list)
for r in raw.values():
    if r['assay_chembl_id'] not in eligible: excluded['endpoint_or_mutant']+=1; continue
    if r.get('standard_type')!='IC50' or r.get('standard_units')!='nM' or r.get('standard_relation')!='=':
        excluded['nonexact_endpoint']+=1; continue
    if r.get('data_validity_comment') or r.get('potential_duplicate'):
        excluded['validity_or_duplicate']+=1; continue
    try:
        value=float(r['standard_value']); assert math.isfinite(value) and value>0
        m=Chem.MolFromSmiles(r['canonical_smiles']); assert m is not None
        assert not any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(m))
        key=identity(m)
    except Exception: excluded['value_structure_stereo']+=1; continue
    groups[(r['assay_chembl_id'],key)].append(dict(activity_id=r['activity_id'],molecule_id=r.get('molecule_chembl_id'),document_id=r.get('document_chembl_id'),smiles=r['canonical_smiles'],ic50_nM=value,pic50=9-math.log10(value)))
records=[]; mols=[]; fps=[]
for (assay,key), rs in sorted(groups.items()):
    mol=Chem.MolFromSmiles(rs[0]['smiles'])
    record=dict(assay_id=assay,normalized_parent=key,scaffold=scaffold(mol),smiles=rs[0]['smiles'],median_pic50=statistics.median(r['pic50'] for r in rs),measurements=rs)
    records.append(record); mols.append(mol); fps.append(fpgen.GetFingerprint(mol))
save('measured_records.json',records)
designs=read('reports/alk5_bounded_search_20261019/input_designs.json')['candidates']
refined=read('reports/alk5_bounded_search_20261019/refinement/results.json')
refined += [dict(r,candidate_id='sar_guided/'+r['candidate_id']) for r in read('reports/alk5_bounded_search_20261019/refinement_sar_guided/results.json') if r['engineering_objectives_met']]
selected=[dict(candidate_id=c['candidate_id'],smiles=c['smiles'],parent_smiles=c['provenance']['parent_smiles']) for c in designs if c['candidate_id'] in ['ALK5_design_011','ALK5_design_020']]+refined
results=[]
for c in selected:
    cm=Chem.MolFromSmiles(c['smiles']); pm=Chem.MolFromSmiles(c['parent_smiles'])
    cf=fpgen.GetFingerprint(cm); pf=fpgen.GetFingerprint(pm); pk=identity(pm); ck=identity(cm)
    neighbors=[]
    for i,r in enumerate(records):
        sim=DataStructs.TanimotoSimilarity(cf,fps[i])
        if sim>=.6:
            neighbors.append(dict(**r,candidate_similarity=sim,parent_similarity=DataStructs.TanimotoSimilarity(pf,fps[i]),is_measured_parent=r['normalized_parent']==pk,is_exact_candidate=r['normalized_parent']==ck))
    neighbors.sort(key=lambda r:(-r['candidate_similarity'],r['assay_id'],r['normalized_parent']))
    parent_rows=[(i,r) for i,r in enumerate(records) if r['normalized_parent']==pk]
    pairs=[]
    for i,pr in parent_rows:
        for j,r in enumerate(records):
            if r['assay_id']!=pr['assay_id'] or r['normalized_parent']==pk or r['scaffold']!=pr['scaffold']: continue
            sim=DataStructs.TanimotoSimilarity(pf,fps[j])
            if sim<.65: continue
            pairs.append(dict(assay_id=r['assay_id'],parent_measurements=pr['measurements'],analogue_measurements=r['measurements'],analogue_smiles=r['smiles'],parent_similarity=sim,candidate_similarity=DataStructs.TanimotoSimilarity(cf,fps[j]),delta_pic50=r['median_pic50']-pr['median_pic50'],activity_cliff=abs(r['median_pic50']-pr['median_pic50'])>=1))
    results.append(dict(candidate_id=c['candidate_id'],smiles=c['smiles'],measured_exact_candidate=any(r['normalized_parent']==ck for r in records),neighbors=neighbors,same_assay_parent_pairs=pairs,activity_prediction=None,activity_validated=False))
save('local_sar.json',results)
save('summary.json',dict(raw_deduplicated=len(raw),eligible_assays=len(eligible),measured_assay_structure_groups=len(records),exclusions=dict(excluded),candidates=[dict(candidate_id=r['candidate_id'],neighbors=len(r['neighbors']),same_assay_parent_pairs=len(r['same_assay_parent_pairs']),activity_cliffs=sum(p['activity_cliff'] for p in r['same_assay_parent_pairs']),measured_exact_candidate=r['measured_exact_candidate']) for r in results]))
print((a.output/'summary.json').read_text(),flush=True)
