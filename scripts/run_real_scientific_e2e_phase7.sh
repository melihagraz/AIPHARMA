#!/usr/bin/env bash
# RUN ONLY ON A TRUSTED LOCAL SCIENTIFIC MACHINE / WSL2 WITH INTERNET + MODELS.
# This script is not a package installer. It fails closed if engines are absent.
set -euo pipefail
cd "$(dirname "$0")/.."
TARGET="${1:-EGFR}"
COUNT="${2:-4}"
python scripts/live_e2e_preflight_phase7.py --target "$TARGET"
python scripts/run_true_three_stage_e2e_20261109.py --target "$TARGET" --count "$COUNT" --out "reports/phase7/real_e2e_${TARGET}.json"
python - "$TARGET" <<'PY'
import json,sys
from pathlib import Path
p=Path(f'reports/phase7/real_e2e_{sys.argv[1]}.json')
r=json.loads(p.read_text())
assert r['full_three_stage_success'] is True, 'Complete real model integration was not verified'
assert all(r['stages'][k]['status'] in {'real_inference_passed','real_model_inference_passed','real_docking_passed'} for k in ('reinvent','admet','vina'))
print('REAL three-stage software integration PASSED:', p)
print('NOT a biological target-validation or experimental drug-discovery claim')
PY
