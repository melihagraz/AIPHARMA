"""Independent inverse-matrix LOO identity and label perturbation audit."""
import json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
p=Path(sys.argv[1]); rows=json.loads((p/'predictions.json').read_text()); groups=defaultdict(list)
for r in rows: groups[r['series']].append(r)
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
maxdiff=0.; maxleak=0.
for rs in groups.values():
    fp=[gen.GetFingerprint(Chem.MolFromSmiles(r['normalized_parent'])) for r in rs]
    K=np.array([DataStructs.BulkTanimotoSimilarity(f,fp) for f in fp]); inv=np.linalg.inv(K+.25*np.eye(len(rs)))
    y=np.array([r['actual_pic50'] for r in rs])
    for i,r in enumerate(rs):
        center=float(np.median(np.delete(y,i))); a=inv@(y-center)
        independent=float(y[i]-a[i]/inv[i,i]); maxdiff=max(maxdiff,abs(independent-r['prediction']))
        perturbed=y.copy(); perturbed[i]+=100
        other=float(perturbed[i]-(inv@(perturbed-center))[i]/inv[i,i]); maxleak=max(maxleak,abs(independent-other))
        assert r['nearest_identity']!=r['normalized_parent']
        assert r['training_analogues']==len(rs)-1
assert maxdiff<1e-9 and maxleak<1e-9
summary=dict(verified_queries=len(rows),series=len(groups),inverse_LOO_max_abs_difference=maxdiff,
    held_label_plus100_prediction_max_change=maxleak,passed=True,
    limits='Numerical and exclusion audit, not independent biological validation')
(p/'verification.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
