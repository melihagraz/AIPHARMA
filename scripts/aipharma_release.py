#!/usr/bin/env python3
"""AIPHARMA local research-software release preparation, stdlib only.

Does not upload, tag, license, or publish anything; does not alter experiments.
"""
from __future__ import annotations
import argparse
import datetime as dt
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import zipfile

VERSION = "0.1.0-rc1"
TAG = "v0.1.0-rc1"
EXCLUDED_PARTS = {
    '.git', '.venv', 'venv', 'env', 'node_modules', '__pycache__', '.pytest_cache',
    '.mypy_cache', '.tox', '.next', '.cache', 'dist', 'build', 'dist_release',
    'reports', 'datasets', 'downloads', 'checkpoints', 'weights', 'secrets',
    '.aws', '.ssh', 'private', 'scratch', 'tmp', 'logs', 'replay_workspaces', 'local',
}
CODE_DIRS = {
    'scripts','tests','aipharma_eave','backend','Reinvent_App','Admet_App',
    'Docking_App','Pocket_App','server','src','frontend','yeni_motor',
    'public','brand','examples','.github','reproducibility','docs',
}
KEEP_SUFFIXES = {
    '.py','.sh','.ps1','.js','.jsx','.mjs','.cjs','.ts','.tsx','.css','.html',
    '.md','.rst','.txt','.toml','.yml','.yaml','.json','.lock','.svg',
    '.png','.jpeg','.jpg','.webp','.ipynb','.cfg','.ini',
}
SKIP_SUFFIXES = {
    '.csv','.tsv','.pdb','.pdbqt','.sdf','.mol2','.ism','.zip','.gz','.tar',
    '.sqlite','.sqlite3','.db','.pem','.key','.p12','.pfx','.pt','.pth','.ckpt',
    '.safetensors','.bin','.npy','.npz','.parquet','.xlsx','.pdf','.docx',
}
# Source assets are included only below max size and after manual review.
MAX_FILE_BYTES = 2_000_000
SECRET_RE = re.compile(
    rb'(-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}|sk-[A-Za-z0-9_-]{30,}|AKIA[0-9A-Z]{16})',
    re.IGNORECASE,
)
ROOT_FILES = {
    'VERSION','CITATION.cff','LICENSE','LICENSE.txt','LICENSE.md',
    'LICENSE_PENDING.md','README.md','SETUP.md','package.json','package-lock.json',
    'pnpm-lock.yaml','yarn.lock','Dockerfile','docker-compose.yml',
    '.gitignore','.dockerignore','.env.example','requirements.txt',
    'requirements-delivery.txt','requirements-real-scientific-validation.txt',
    'requirements-delivery-linux-py311.lock.txt','REPRODUCIBILITY.md',
    'DATA_AVAILABILITY.md','THIRD_PARTY_NOTICES.md','CHANGELOG.md',
    'eslint.config.js','vite.config.js','vitest.config.js','tsconfig.json',
    'tsconfig.node.json','index.html',
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def cmd(argv: list[str], cwd: Path, timeout=12) -> dict:
    try:
        p = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return {'ok': p.returncode == 0, 'returncode': p.returncode,
                'stdout': p.stdout.strip()[:1500], 'stderr': p.stderr.strip()[:700]}
    except (OSError, subprocess.TimeoutExpired) as e:
        return {'ok': False, 'error': type(e).__name__}


def runcheck(root: Path) -> dict:
    must = [
        'scripts/phase9_redock.py','scripts/phase10_run_matched_egfr.py',
        'scripts/phase11_alk5.py','scripts/phase12_alk5_3hmm.py',
        'scripts/eave_audit.py','scripts/phase14_eave_gate.py',
        'scripts/phase15_live_replay.py',
        'aipharma_eave/audit.py','aipharma_eave/phase14.py',
        'tests/test_eave_engine.py','tests/test_phase14_gate.py',
        'tests/test_phase15_live_replay.py',
        'requirements-real-scientific-validation.txt','package-lock.json',
    ]
    missing = [f for f in must if not (root / f).is_file()]
    required_evidence = [
        'reports/phase10/egfr_matched50_docking/analysis.json',
        'reports/phase11/alk5_matched50_docking/analysis.json',
    ]
    evidence = {k: (root / k).is_file() for k in required_evidence}
    phase12 = list((root / 'reports/phase12').rglob('analysis.json')) if (root / 'reports/phase12').is_dir() else []
    has_license = any((root / n).is_file() for n in ('LICENSE','LICENSE.txt','LICENSE.md'))
    license_clearance = (root / 'reproducibility' / 'DISTRIBUTION_APPROVED.txt').is_file()
    has_lock = bool(list((root / 'reproducibility' / 'local').glob('pip-freeze-*.txt')))
    has_capture = (root / 'reproducibility' / 'local' / 'environment.json').is_file()
    has_public_lock = bool(list((root / 'reproducibility' / 'locks').glob('pip-freeze-*.txt')))
    has_citation = (root / 'CITATION.cff').is_file()
    scientific_conflicts = []
    rein = root / 'Reinvent_App' / 'requirements.txt'
    if rein.is_file() and 'numpy<2' in rein.read_text(encoding='utf8', errors='replace').replace(' ', ''):
        for local_lock in (root / 'reproducibility' / 'local').glob('pip-freeze-*.txt'):
            if re.search(r'(?im)^numpy==2(?:\.|$)', local_lock.read_text(encoding='utf8',errors='replace')):
                scientific_conflicts.append('observed_numpy_2_vs_Reinvent_App_numpy_lt_2')
                break
    readme_text = (root / 'README.md').read_text(encoding='utf8',errors='replace') if (root / 'README.md').exists() else ''
    legacy_readme = ('Phase 6 update' in readme_text and 'Phase 15' not in readme_text)
    has_git = (root / '.git').exists()
    remote = cmd(['git','remote','-v'], root) if has_git else {'ok': False, 'error': 'not_a_git_repo'}
    tag = cmd(['git', 'tag', '-l', TAG], root) if has_git else {'ok': False}
    # Results are not promised as distributable even if local files exist.
    d = {
      'version': VERSION,'intended_tag': TAG, 'local_code_present': not missing,
      'missing_components': missing,
      'evidence_files_present': evidence, 'phase12_analysis_files_found': len(phase12),
      'root_license_file_present': has_license,
      'distribution_rights_explicitly_approved': license_clearance,
      'citation_file_present': has_citation,
      'legacy_root_README_may_need_update': legacy_readme,
      'local_python_lock_present': has_lock,
      'observed_scientific_dependency_conflicts': scientific_conflicts,
      'environment_capture_present': has_capture,
      'reviewed_public_dependency_lock_present': has_public_lock,
      'git_repository_present': has_git,
      'git_tag_exists_locally': bool(tag.get('stdout')),
      'git_remote_configured': bool(remote.get('stdout')),
      'public_repository_verified': False,
      'doi_verified': False,
      'dataset_distribution_rights_verified': False,
    }
    d['publication_blockers'] = [key for key, value in [
        ('missing_current_pipeline_code', not missing),
        ('root_LICENSE_missing', has_license),
        ('third_party_rights_not_approved', license_clearance),
        ('CITATION_missing', has_citation),
        ('legacy_readme_needs_review', not legacy_readme),
        ('local_pip_freeze_missing', has_lock),
        ('environment_capture_missing', has_capture),
        ('reviewed_public_dependency_lock_missing', has_public_lock),
        ('scientific_dependency_conflicts_detected', not scientific_conflicts),
        ('git_tag_not_verified', bool(tag.get('stdout'))),
        ('public_repository_not_verified', False),
        ('dataset_distribution_not_verified', False),
    ] if not value]
    return d


def snapshot(root: Path) -> Path:
    out = root / 'reproducibility' / 'local'
    out.mkdir(parents=True, exist_ok=True)
    pyv = f'py{sys.version_info.major}{sys.version_info.minor}'
    lock = out / f'pip-freeze-{pyv}.txt'
    # Do not use cmd() here: its diagnostic output is intentionally truncated.
    # The lock file must contain every line of the actual pip freeze.
    try:
        frozen = subprocess.run([sys.executable, '-m', 'pip', 'freeze'], cwd=root,
                                capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise RuntimeError(f'pip freeze failed: {type(e).__name__}') from e
    if frozen.returncode != 0:
        raise RuntimeError('pip freeze failed; active environment cannot be captured')
    lock.write_text(frozen.stdout.rstrip('\n') + '\n', encoding='utf-8')
    vina_bin = shutil.which('vina')
    vina_version = cmd([vina_bin, '--version'], root) if vina_bin else {'ok': False,'error':'vina_not_on_PATH'}
    # Never store environment variables or private credentials.
    spec = {
        'schema': 'aipharma-release-environment-v1',
        'captured_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'release_candidate': VERSION,
        'os': platform.platform(),
        'machine': platform.machine(),
        'python_version': platform.python_version(),
        'python_implementation': platform.python_implementation(),
        'python_lock_path': str(lock.relative_to(root)),
        'python_lock_sha256': sha(lock),
        'node_version': cmd(['node','--version'], root),
        'npm_version': cmd(['npm','--version'], root),
        'npm_lock_sha256': sha(root / 'package-lock.json') if (root / 'package-lock.json').is_file() else None,
        'vina_version': vina_version,
        'pip_check': cmd([sys.executable, '-m', 'pip', 'check'], root, timeout=35),
        'vina_binary_sha256': sha(Path(vina_bin)) if vina_bin and Path(vina_bin).is_file() else None,
        'git_head': cmd(['git','rev-parse','HEAD'], root),
        'note': 'Freeze is an observed environment, not a clean-install success guarantee. Model weights, datasets and Vina executable are not bundled.',
    }
    (out / 'environment.json').write_text(json.dumps(spec, indent=2) + '\n', encoding='utf8')
    return out / 'environment.json'


def included_file(root: Path, p: Path) -> bool:
    rel = p.relative_to(root)
    parts = rel.parts
    if any(part in EXCLUDED_PARTS or part.startswith('.env.') for part in parts):
        return False
    if p.is_symlink() or not p.is_file():
        return False
    if p.stat().st_size > MAX_FILE_BYTES:
        return False
    if p.suffix.lower() in SKIP_SUFFIXES:
        return False
    if p.name in ('.env', 'credentials.json', 'service-account.json', 'id_rsa'):
        return False
    # Explicit roots only; do not accidentally publish confidential user files.
    if len(parts) == 1:
        return p.name in ROOT_FILES or (p.name.startswith('requirements-') and p.suffix=='.txt')
    if parts[0] not in CODE_DIRS:
        return False
    if 'data' in (x.lower() for x in parts) or 'inputs' in (x.lower() for x in parts):
        return False
    return p.suffix.lower() in KEEP_SUFFIXES or p.name in ('Dockerfile','.env.example','.gitignore')


def collect(root: Path) -> list[Path]:
    paths = []
    for dirname in sorted(CODE_DIRS):
        folder=root/dirname
        if folder.is_dir():
            for d, dirs, files in os.walk(folder, followlinks=False):
                dirs[:] = sorted(x for x in dirs if x not in EXCLUDED_PARTS and not x.startswith('.env.') and not (Path(d)/x).is_symlink())
                for filename in files:
                    p = Path(d)/filename
                    if included_file(root,p): paths.append(p)
    for item in ROOT_FILES:
        p = root/item
        if included_file(root,p): paths.append(p)
    return sorted(set(paths), key=lambda p: str(p.relative_to(root)))


def assemble(root: Path, *, strict=False) -> dict:
    status = runcheck(root)
    if strict and status['publication_blockers']:
        raise RuntimeError('PUBLIC RELEASE BLOCKED: '+', '.join(status['publication_blockers']))
    candidates = collect(root)
    if len(candidates) < 15:
        raise RuntimeError('Suspiciously few source files: '+str(len(candidates)))
    suspect=[]
    for p in candidates:
        data=p.read_bytes()
        if SECRET_RE.search(data):suspect.append(str(p.relative_to(root)))
    if suspect:
        raise RuntimeError('Potential secrets in source files; archive blocked: '+', '.join(suspect[:12]))
    outdir=root/'dist_release'
    outdir.mkdir(exist_ok=True)
    output=outdir/f'AIPHARMA_{VERSION}_source_candidate.zip'
    manifest={'schema':'aipharma-source-candidate-v1','version':VERSION,
        'status':'LOCAL_CANDIDATE_NOT_PUBLIC_RELEASE',
        'source_count':len(candidates), 'publication_blockers': status['publication_blockers'],
        'excluded_data_and_weights':True,
        'files':{str(p.relative_to(root)):sha(p) for p in candidates}}
    # Build reproducible deterministic candidate bytes (file metadata independent of mtime).
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=8) as z:
        for p in candidates:
            arc=f'AIPHARMA-{VERSION}/'+p.relative_to(root).as_posix()
            info=zipfile.ZipInfo(arc,date_time=(2026,11,10,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=(0o644 << 16)
            z.writestr(info,p.read_bytes())
        info=zipfile.ZipInfo(f'AIPHARMA-{VERSION}/RELEASE_MANIFEST.json',date_time=(2026,11,10,0,0,0))
        info.compress_type=zipfile.ZIP_DEFLATED
        z.writestr(info,json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    checksum = outdir/(output.name+'.sha256')
    checksum.write_text(f'{sha(output)}  {output.name}\n',encoding='utf8')
    (outdir/'release_audit.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf8')
    return {'archive':str(output),'archive_sha256':sha(output),'files':len(candidates),
        'publishable':not bool(status['publication_blockers']),
        'publication_blockers':status['publication_blockers']}


def test(root: Path):
    specs=['test_eave_engine.py','test_phase14_gate.py','test_phase15_live_replay.py','test_release_packaging.py']
    results={}
    for name in specs:
        found=(root/'tests'/name).is_file()
        if not found:
            results[name]={'ok':False,'error':'missing'}
            continue
        results[name]=cmd([sys.executable,'-m','unittest','discover','-s','tests','-p',name,'-q'],root,timeout=120)
    outdir=root/'reproducibility'/'local';outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'test_results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8')
    return results


def main():
    ap=argparse.ArgumentParser(description='Local, non-publishing AIPHARMA reproducibility release tooling')
    ap.add_argument('action',choices=('inspect','freeze','test','bundle','strict-check'))
    ap.add_argument('--root',default=None,help='Project root; defaults to parent of scripts/')
    args=ap.parse_args()
    root=Path(args.root).expanduser().resolve() if args.root else Path(__file__).resolve().parents[1]
    if not (root/'scripts').is_dir():
        ap.error('Not an AIPHARMA project root')
    try:
        if args.action in ('inspect','strict-check'):
            data=runcheck(root)
            print(json.dumps(data,indent=2))
            return 2 if args.action=='strict-check' and data['publication_blockers'] else 0
        if args.action=='freeze':
            p=snapshot(root); print('Environment captured:',p);return 0
        if args.action=='test':
            data=test(root);print(json.dumps(data,indent=2))
            return 0 if all(x.get('ok') for x in data.values()) else 2
        if args.action=='bundle':
            print(json.dumps(assemble(root),indent=2));return 0
    except Exception as e:
        print(f'ERROR: {e}',file=sys.stderr)
        return 2
    return 1

if __name__=='__main__':
    raise SystemExit(main())
