"""Create explicit source bundles; never include personal env, DB or model weights."""
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

root = Path(__file__).resolve().parents[1]
report = root/'reports/alk5_delivery_20261019'
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
download = json.loads((report/'downloaded_designs_with_profiles.json').read_text())
assert len(download['candidates']) == len(download['admet_comparison']['pairs']) == 30
assert all(a['smiles'] == b['smiles'] and a['candidate_id'] == b['candidate_id'] for a, b in zip(download['candidates'], download['admet_comparison']['pairs']))
browser = dict(status='passed_for_defined_local_workflow', actual_browser_sign_in=True, empty_database_access=True, real_design_generation=True, real_admet_service=True, displayed_designs=30, actual_browser_download=True, downloaded_pairs=30, downloaded_design_profile_alignment=True, download_sha256=digest(report/'downloaded_designs_with_profiles.json'), screenshot='design_admet_ui.jpg', frontend_targeted_tests=6, fresh_backend_tests=8, ranking_regression_tests=3, general_usability_study=False, clean_scientific_installation=False, source_models_environment_reused=True)
(report/'browser_verification.json').write_text(json.dumps(browser, indent=2))
latex = root/'reports/latex'
tex = latex/'PharmaProject_arXiv_manuscript_UPDATED.tex'
assert tex.is_file()
latex_zip = root/'reports/PharmaProject_LaTeX_DELIVERY_20261019.zip'
with ZipFile(latex_zip, 'w', ZIP_DEFLATED) as z:
    z.write(tex, tex.name)
    for folder in ['media', 'figures']:
        for f in sorted((latex/folder).rglob('*')):
            if f.is_file():
                z.write(f, str(f.relative_to(latex)))
    z.writestr('README.txt', 'Main file: PharmaProject_arXiv_manuscript_UPDATED.tex\nAll original and added figures are included.\nPDF compilation NOT verified: built-in Windows preview failed to start (CreateProcessWithLogonW 267).\nThis is a source bundle, not a verified new PDF.\n')
    z.write(report/'manuscript_source_check.json', 'source_check.json')
source_zip = root/'reports/AIPHARMA_delivery_20261019.zip'
with ZipFile(source_zip, 'w', ZIP_DEFLATED) as z:
    for folder in ['backend/app', 'Admet_App/app', 'dist', 'src', 'public']:
        for f in sorted((root/folder).rglob('*')):
            if not f.is_file() or '__pycache__' in f.parts or f.suffix in {'.pyc', '.db', '.sqlite3'} or f.name.startswith('.env'):
                continue
            z.write(f, str(f.relative_to(root)))
    for name in ['requirements.txt', 'requirements-delivery.txt', 'requirements-delivery-linux-py311.lock.txt', 'SETUP.md', 'package.json', 'package-lock.json', 'index.html', 'vite.config.js']:
        z.write(root/name, name)
    for name in ['build-delivery.mjs', 'run_delivery_http_20261019.py', 'vite-denovo-local.mjs', 'dev-orchestrator.mjs']:
        z.write(root/'scripts'/name, 'scripts/'+name)
    for folder in ['Admet_App']:
        for name in ['NOTICE', 'LICENSE', 'README.md', 'requirements.txt']:
            f = root/folder/name
            if f.is_file():
                z.write(f, str(f.relative_to(root)))
    z.write(report/'downloaded_designs_with_profiles.json', 'examples/ALK5_designs_with_profiles.json')
    z.write(report/'README.md', 'DELIVERY_EVIDENCE.md')
    z.writestr('README.txt', 'AIPHARMA selected ALK5 design/profiling delivery path.\nSee SETUP.md and DELIVERY_EVIDENCE.md.\nBuilt UI is in dist; start the main API with its delivery environment.\nAdmet_App requires its separate scientific runtime and model weights, which are NOT bundled.\nREINVENT/pocket/docking/LLM services and public deployment are outside this checked package.\nNo personal .env, existing database, venv or model weights are included.\nNo molecular activity or safety validation is claimed.\n')
for p in [latex_zip, source_zip]:
    with ZipFile(p) as z:
        assert z.testzip() is None
        assert not any(Path(n).name.startswith('.env') or Path(n).suffix in {'.db', '.sqlite3'} for n in z.namelist())
files = {str(f.relative_to(report)).replace('\\', '/'): digest(f) for f in sorted(report.rglob('*')) if f.is_file() and f.name != 'MANIFEST_SHA256.json'}
(report/'MANIFEST_SHA256.json').write_text(json.dumps(files, indent=2))
print(json.dumps(dict(packages=[dict(path=str(p), bytes=p.stat().st_size, sha256=digest(p)) for p in [source_zip, latex_zip]], report_files=len(files)), indent=2))
