"""One fixed model on a source-separated screen; no hyperparameter selection."""
import hashlib,json,sys
from pathlib import Path
from collections import Counter
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
p=Path(sys.argv[1]); k=p/'kirhub'; out=k/'frozen_rf_screen'; out.mkdir(exist_ok=False)
parameters=dict(n_estimators=300,min_samples_leaf=3,max_features=.5,random_state=20261019,n_jobs=4)
protocol=dict(model='Previous fixed Morgan2048 radius2 RandomForest configuration; no search',parameters=parameters,
 training='Exact-only records from scaffold-excluded expanded training; equal total weight per normalized identity',
 censored='Retained in curated data but excluded from this exact-response estimator; never imputed',
 external='Whole KIRHub source excluded from fit; retrospective and target values already inspected, not blind',
 metrics='Spearman prediction versus percent inhibition, fixed-prediction bootstrap and permutation; no IC50 MAE across endpoints',
 AUC='Only if >=5 structures per prespecified extreme class; otherwise not calculated',
 promotion='Diagnostic only: does not supersede failed development acceptance or authorize generated-candidate ranking',
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'protocol.json').write_text(json.dumps(protocol,indent=2))
training=json.loads((k/'overlap_audit/scaffold_excluded_training.json').read_text())
train=[r for r in training if r['median_pic50'] is not None]
test=[r for r in json.loads((k/'overlap_audit/structure_records.json').read_text()) if r['status']=='normalized_for_overlap_audit']
assert {r['normalized_parent'] for r in train}.isdisjoint(r['normalized_parent'] for r in test)
assert {r['scaffold'] for r in train}.isdisjoint(r['scaffold'] for r in test)
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def fp(s):
    a=np.zeros(2048,dtype=np.uint8); DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(s)),a); return a
X=np.array([fp(r['smiles']) for r in train]); Xt=np.array([fp(r['source_smiles']) for r in test])
counts=Counter(r['normalized_parent'] for r in train)
model=RandomForestRegressor(**parameters).fit(X,[r['median_pic50'] for r in train],sample_weight=[1/counts[r['normalized_parent']] for r in train])
pred=model.predict(Xt)
# Save structure-only predictions before merging measured external target values.
predictions=[dict(compound=r['compound'],prediction_pic50=float(v)) for r,v in zip(test,pred)]
(out/'predictions_before_endpoint_merge.json').write_text(json.dumps(predictions,indent=2))
measured={r['compound']:r for r in json.loads((k/'reserved_alk5_single_dose.json').read_text())}
observed=np.array([measured[r['compound']]['inhibition_percent'] for r in test])
rho=float(spearmanr(pred,observed).statistic); rng=np.random.default_rng(20261019)
boot=[]
for _ in range(5000):
    ix=rng.integers(0,len(test),len(test)); value=float(spearmanr(pred[ix],observed[ix]).statistic)
    if np.isfinite(value): boot.append(value)
null=np.array([spearmanr(pred,rng.permutation(observed)).statistic for _ in range(5000)])
labels=Counter(measured[r['compound']]['label'] for r in test)
for r in predictions:
    r.update(inhibition_percent=measured[r['compound']]['inhibition_percent'],screen_label=measured[r['compound']]['label'])
predictions.sort(key=lambda r:-r['prediction_pic50'])
for rank,r in enumerate(predictions,1): r['predicted_rank']=rank
(out/'screen_results.json').write_text(json.dumps(predictions,indent=2))
summary=dict(training_exact_groups=len(train),training_unique_identities=len(counts),excluded_censored_groups=len(training)-len(train),
 evaluated_drugs=len(test),labels=dict(labels),spearman_rho=rho,
 fixed_prediction_drug_bootstrap_95=[float(x) for x in np.quantile(boot,[.025,.975])],
 two_sided_permutation_p=float((1+sum(np.abs(null)>=abs(rho)))/(len(null)+1)),
 AUC_calculated=False,IC50_MAE_calculated=False,strong_screen_inhibitor_ranks=[dict(compound=r['compound'],rank=r['predicted_rank']) for r in predictions if r['screen_label']=='screen_inhibitor'],
 caveats=['One external source, not multiple source-level replications','Mostly low inhibition with measurement floor ties','Bootstrap conditional on fitted predictions, not refitting uncertainty','Screen at 1uM and training IC50 are different endpoints'],
 candidate_ranking_authorized=False,activity_or_safety_validated=False)
(out/'summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
