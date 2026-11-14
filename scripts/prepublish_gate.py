#!/usr/bin/env python3
"""Read-only release-preflight checker for a proposed public AIPHARMA Git repository.
Does not read passwords, upload, push, alter, delete or publish files.
Check output is only a local preliminary review; manual IP/legal/data audit remains essential.
"""
from __future__ import annotations
import argparse
import re
import subprocess
from pathlib import Path

RED = re.compile(r'(^|/)(?:\.env(?:\..*)?|\.venv|venv|__pycache__|node_modules|\.git|private|secrets?|credentials?|tokens?|patient|phi)(/|$)|\.(?:pem|key|p12|pfx|sqlite3?|db)$', re.I)
BAD_LICENSE = re.compile(r'pending|to be determined|tbd|placeholder|TODO', re.I)
DOI = re.compile(r'10\.\d{4,9}/[^\s,;\]\)]+', re.I)

def git(root: Path, *args: str):
    try:
        cp = subprocess.run(['git', '-C', str(root), *args], capture_output=True, text=True, timeout=12)
        return cp.stdout.strip() if cp.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        return None

def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--root', default='.', help='Local repository root; nothing is modified')
    p.add_argument('--require-doi', action='store_true', help='Require DOI in CITATION.cff and a real GitHub remote')
    args = p.parse_args()
    root = Path(args.root).expanduser().resolve()
    blockers, warnings, passed = [], [], []
    def require(name: str, pred, why: str):
        (passed if pred else blockers).append(f'{name}: {"OK" if pred else why}')
    if not root.is_dir():
        print(f'BLOCK: no directory {root}'); return 2
    for path in ['README.md','VERSION','CITATION.cff','CHANGELOG.md','DATA_AVAILABILITY.md','REPRODUCIBILITY.md','THIRD_PARTY_NOTICES.md']:
        require(path, (root/path).is_file(), 'MISSING')
    license_path = root/'LICENSE'
    if not license_path.is_file():
        blockers.append('LICENSE: MISSING (LICENSE_PENDING.md is NOT a license)')
    else:
        try:
            lic = license_path.read_text('utf-8', errors='ignore')
            require('LICENSE', len(lic.strip()) >= 100 and not BAD_LICENSE.search(lic[:1000]), 'EMPTY / PLACEHOLDER; legal review needed')
        except OSError:
            blockers.append('LICENSE: UNREADABLE')
    if (root/'LICENSE_PENDING.md').is_file():
        warnings.append('LICENSE_PENDING.md exists; remove when actual root LICENSE is approved')
    if (root/'VERSION').exists():
        version = (root/'VERSION').read_text('utf-8', errors='ignore').strip()
        require('VERSION syntax', bool(re.fullmatch(r'\d+\.\d+\.\d+(?:[-+][0-9A-Za-z.-]+)?',version)), 'INVALID')
        if 'rc' in version.lower():
            warnings.append(f'Version {version} is a release candidate; do not label a final release without review')
    else: version = '<missing>'
    has_locks = any((root / f).exists() for f in ['requirements.lock','uv.lock','poetry.lock','conda-lock.yml', 'reproducibility/locks/pip-freeze.txt','reproducibility/locks/pip-freeze-WSL.txt'])
    require('Python dependency lock', has_locks, 'MISSING (capture reviewed working environment)')
    has_tests = (root/'tests').is_dir() and any((root/'tests').glob('test*.py'))
    require('Test suite', has_tests, 'MISSING tests/test*.py')
    has_source = any((root/x).is_dir() for x in ('backend','frontend','scripts','eave','src'))
    require('Source tree', has_source, 'MISSING')
    provenance = root/'DATA_AVAILABILITY.md'
    if provenance.is_file() and re.search(r'pending|not complete|not publicly|todo',provenance.read_text('utf-8',errors='ignore'),re.I):
        warnings.append('DATA_AVAILABILITY.md appears unfinished; manually verify redistributable inputs, panel hashes and source URLs')
    tracked_raw = git(root,'ls-files','--cached','-z')
    if tracked_raw is None:
        blockers.append('Git repository: NOT INITIALIZED or git unavailable')
    else:
        # -z output strips trailing NUL in stdout.strip(); names are not otherwise modified.
        tracked = [x for x in tracked_raw.split('\0') if x]
        require('Git tracked content',len(tracked)>0,'NO STAGED/TRACKED FILES')
        risks = [name for name in tracked if not name.endswith('.env.example') and RED.search(name.replace('\\','/'))]
        if risks:
            blockers.append('SENSITIVE/NON-DISTRIBUTABLE tracked paths ('+str(len(risks))+'): '+', '.join(risks[:12])+(f' ... ({len(risks)-12} more)' if len(risks)>12 else ''))
        large = []
        for name in tracked:
            f=root/name
            try:
                if f.is_file() and f.stat().st_size > 25*1024*1024: large.append(name)
            except OSError: pass
        if large:
            warnings.append('Tracked large files (>25MiB), review for models/raw data: '+', '.join(large[:10]))
        # Presence of a .gitignore alone cannot prevent already tracked data from publishing.
    remote = git(root,'remote','get-url','origin')
    if remote:
        passed.append('Git origin: CONFIGURED (public reachability NOT verified)')
        if 'github.com' not in remote.lower(): warnings.append('Git origin is not GitHub; ensure selected archive integration works')
    else:
        blockers.append('Git origin: MISSING (create repo only after ownership review)')
    citation = (root/'CITATION.cff').read_text('utf-8',errors='ignore') if (root/'CITATION.cff').is_file() else ''
    if not re.search(r'^authors\s*:',citation,re.M): blockers.append('CITATION.cff: authors metadata MISSING')
    if re.search(r'example\.com|example-user|YOUR_|TODO|pending', citation, re.I):
        blockers.append('CITATION.cff: contains placeholders')
    doi = DOI.search(citation)
    if doi: passed.append('CITATION.cff: DOI PRESENT (not checked against Zenodo or release)')
    elif args.require_doi: blockers.append('CITATION.cff: RELEASE DOI MISSING')
    else: warnings.append('CITATION.cff: release DOI missing; REQUIRED before final submission')
    tag=git(root,'tag','--list',f'v{version}')
    if tag: passed.append(f'Local tag v{version}: PRESENT (remote/release not checked)')
    else: warnings.append(f'Local tag v{version} not created yet')
    print(f'AIPHARMA PUBLICATION GATE — {root}\n')
    for x in passed: print('PASS', x)
    for x in warnings: print('WARN', x)
    for x in blockers: print('BLOCK', x)
    print(f'\nSUMMARY: {len(blockers)} blocking local conditions, {len(warnings)} warnings. This is NOT a public-link/rights audit.')
    return 1 if blockers else 0
if __name__=='__main__': raise SystemExit(main())
