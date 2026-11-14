#!/usr/bin/env python3
"""Evidence-based preflight for real AIPHARMA REINVENT -> ADMET-AI -> Vina.

Does not install packages, simulate predictions, or claim assay validation.
Returns 0 only when runtime prerequisites exist. A successful preflight is
NOT proof the complete model run or experimental validation succeeded.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from importlib import metadata, util
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def checksum(path: Path):
    h = sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def version(name):
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def preflight(target: str) -> dict:
    from Docking_App.app.dock import _vina_binary, list_targets, _resolve_target

    prior = ROOT / 'Reinvent_App/prior/reinvent.prior'
    available_targets = list_targets()
    report = {
        'audit_kind': 'real_engine_preflight_not_an_experiment',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'target': target,
        'python': sys.version.split()[0],
        'platform': platform.platform(),
        'dependencies': {name: {'installed': util.find_spec(module) is not None, 'version': version(dist)}
                         for name, module, dist in [
                             ('admet_ai', 'admet_ai', 'admet-ai'),
                             ('meeko', 'meeko', 'meeko'),
                             ('rdkit', 'rdkit', 'rdkit'),
                             ('torch', 'torch', 'torch')]},
        'model_prior_sha256': checksum(prior) if prior.is_file() else None,
        'available_target_names': available_targets,
        'scientific_scope': 'workflow smoke test; NOT ALK5 activity or safety validation',
        'blockers': [],
    }
    if not prior.is_file():
        report['blockers'].append('REINVENT prior weights missing')
    if not report['dependencies']['admet_ai']['installed']:
        report['blockers'].append('ADMET-AI real package/models missing')
    if not report['dependencies']['meeko']['installed']:
        report['blockers'].append('Meeko ligand-preparation package missing')
    try:
        receptor, conf = _resolve_target(target)
        report['receptor_sha256'] = checksum(receptor)
        report['grid_config_sha256'] = checksum(conf)
    except ValueError:
        report['blockers'].append(f'No PDBQT receptor and Vina grid config for target {target!r}')
    try:
        binary = Path(_vina_binary())
        proc = subprocess.run([str(binary), '--version'], capture_output=True,
                              text=True, timeout=12, check=False)
        report['vina_binary'] = str(binary)
        report['vina_binary_sha256'] = checksum(binary) if binary.is_file() else None
        report['vina_version_output'] = (proc.stdout + proc.stderr).strip()[:400]
        if proc.returncode != 0:
            report['blockers'].append('Vina --version failed')
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired) as exc:
        report['blockers'].append('Native Vina CLI missing/non-executable: ' + str(exc)[:180])
    report['prerequisites_present'] = not report['blockers']
    report['live_three_stage_run_performed'] = False
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', default='EGFR')
    parser.add_argument('--output', type=Path, default=ROOT / 'reports/phase7/engine_preflight.json')
    args = parser.parse_args()
    report = preflight(args.target)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'prerequisites_present': report['prerequisites_present'],
                      'blockers': report['blockers'], 'output': str(args.output)}, indent=2))
    return 0 if report['prerequisites_present'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
