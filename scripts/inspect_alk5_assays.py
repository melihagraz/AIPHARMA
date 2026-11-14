from pathlib import Path
import json
from collections import Counter
p=Path.home()/'aipharma_holdout_20261016';a=json.loads((p/'raw/boundary_assays_all.json').read_text())['assays']
for r in a:print(json.dumps({k:r.get(k) for k in ['assay_chembl_id','assay_type','confidence_score','description','bao_format']},ensure_ascii=False))
