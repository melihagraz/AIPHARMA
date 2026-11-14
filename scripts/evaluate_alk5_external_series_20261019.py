"""Diagnostic external-source check of frozen potency baselines, no retuning."""
import json,sys,math,statistics,hashlib
from pathlib import Path
from collections import defaultdict
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import spearmanr
p=Path(sys.argv[1]); out=p/'external_evaluation_deduplicated'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('protocol.json',dict(methods=['nearest Morgan2/2048 reference','training median'],
    no_retuning=True,source='US20250179094A1 / BindingDB entry12979',
    test='63 patent examples before normalized-identity/stereo exclusions; LY364947 positive control excluded; median pIC50 per normalized identity, all aliases and potency ranges retained',
    train='Previously extracted exact catalytic IC50; compound medians across eligible assays',
    exclusion='Remove all training structures sharing normalized identity or exact Murcko scaffold with any eligible external test structure',
    endpoint='Continuous pIC50; active-only series cannot evaluate active versus weak discrimination',
    metrics=['MAE','RMSE','Spearman correlation'],
    independence='New source evaluation after model fixed, but source potency values already inspected. Not blind. Registry source-family check is limited, not an exhaustive external audit.',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
enum=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger()
def norm(m):
    m=enum.Canonicalize(un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))))
    Chem.RemoveStereochemistry(m); return Chem.MolToSmiles(m)
def scaf(m): return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
external=json.loads((p/'external_source/parsed_series.json').read_text()); test=[]; excluded=[]
for r in external:
    if r['primary_table_key']=='LY364947': excluded.append(dict(monomer_id=r['monomer_id'],reason='positive control')); continue
    m=Chem.MolFromSmiles(r['smiles'])
    if m is None or any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(m)):
        excluded.append(dict(monomer_id=r['monomer_id'],reason='structure/stereo')); continue
    test.append(dict(**r,normalized_parent=norm(m),scaffold=scaf(m),pic50=9-math.log10(r['ic50_nM'])))
test_groups=defaultdict(list)
for r in test: test_groups[r['normalized_parent']].append(r)
test=[dict(**{k:v for k,v in rs[0].items() if k!='pic50'},pic50=statistics.median(r['pic50'] for r in rs),source_measurements=[dict(monomer_id=r['monomer_id'],compound_number=r['compound_number'],ic50_nM=r['ic50_nM']) for r in rs],source_pic50_range=[min(r['pic50'] for r in rs),max(r['pic50'] for r in rs)]) for rs in test_groups.values()]
keys={r['normalized_parent'] for r in test}; scaffolds={r['scaffold'] for r in test}
records=json.loads((p/'local_sar/measured_records.json').read_text()); grouped=defaultdict(list); train_exclusions=[]
for r in records:
    if r['normalized_parent'] in keys or r['scaffold'] in scaffolds:
        train_exclusions.append(dict(normalized_parent=r['normalized_parent'],assay_id=r['assay_id'],reason='external identity/scaffold')); continue
    grouped[r['normalized_parent']].append(r)
train=[dict(normalized_parent=k,smiles=rs[0]['smiles'],pic50=statistics.median(r['median_pic50'] for r in rs)) for k,rs in sorted(grouped.items())]
assert train
fps=[gen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in train]
baseline=statistics.median(r['pic50'] for r in train); predictions=[]
for r in test:
    similarities=DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])),fps)
    j=max(range(len(train)),key=lambda j:similarities[j])
    predictions.append(dict(**r,nearest_pic50=train[j]['pic50'],median_pic50=baseline,nearest_similarity=similarities[j],nearest_training_identity=train[j]['normalized_parent']))
metrics={}
for label in ['nearest','median']:
    actual=[r['pic50'] for r in predictions]; pred=[r[label+'_pic50'] for r in predictions]
    rho=float(spearmanr(actual,pred).statistic) if len(set(pred))>1 else None
    metrics[label]=dict(MAE=statistics.mean(abs(x-y) for x,y in zip(actual,pred)),RMSE=math.sqrt(statistics.mean((x-y)**2 for x,y in zip(actual,pred))),Spearman=rho if rho is not None and math.isfinite(rho) else None)
save('predictions.json',predictions); save('training_structure_medians.json',train)
save('exclusions.json',dict(test=excluded,training=train_exclusions))
summary=dict(test_count=len(test),training_unique_compounds=len(train),excluded_training_records=len(train_exclusions),metrics=metrics,nearest_similarity_range=[min(r['nearest_similarity'] for r in predictions),max(r['nearest_similarity'] for r in predictions)],candidate_activity_validated=False,classification_validation=False)
save('summary.json',summary); print(json.dumps(summary,indent=2),flush=True)
