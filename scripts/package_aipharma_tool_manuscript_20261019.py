"""Package the existing manuscript and every referenced figure; no PDF compilation."""
from pathlib import Path
import hashlib,json,re,zipfile
root=Path.cwd(); tex=root/'reports/latex/PharmaProject_arXiv_manuscript_UPDATED.tex'; text=tex.read_text(encoding='utf-8')
destination=tex.parent/'PharmaProject_arXiv_manuscript_UPDATED_toolpaper_20261019.zip'
assert not destination.exists(), 'Preserve existing package; select a new version name'
items={tex.name:tex}
for name in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',text): items[name]=tex.parent/name
reports=[
 'reports/manuscript_static_check_20261019.json',
 'reports/ALK5_CURRENT_STATUS.md',
 'reports/alk5_bounded_search_20261019/decision.json',
 'reports/alk5_data_expansion_20261019/curated/summary.json',
 'reports/alk5_data_expansion_20261019/within_series/summary.json',
 'reports/alk5_data_expansion_20261019/within_series/protocol.json',
 'reports/alk5_data_expansion_20261019/external_local_transfer/interpretation.json',
 'reports/alk5_data_expansion_20261019/kirhub/frozen_rf_screen/summary.json',
 'reports/alk5_data_expansion_20261019/kirhub/frozen_rf_screen/protocol.json',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/README.md',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/review_summary.json',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/review_records.json',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/within_series/summary.json',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/within_series/protocol.json',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/within_series/predictions.json',
 'reports/alk5_data_expansion_20261019/additional_local_series/drawing_audit_current/within_series/verification.json',
 'scripts/evaluate_alk5_within_series_20261019.py',
 'scripts/verify_alk5_within_series_20261019.py',
 'scripts/check_aipharma_manuscript_source_20261019.py']
for name in reports: items['evidence/'+name]=root/name
assert all(f.is_file() for f in items.values())
manifest={name:dict(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for name,f in items.items()}
readme='''AIPHARMA tool-paper source package — 19 October 2026

Main source: PharmaProject_arXiv_manuscript_UPDATED.tex. This is the full existing manuscript revised in place; all 11 referenced illustrations are included at their original relative paths. No replacement condensed manuscript or newly compiled PDF is included.

Scope: integrated molecular design, property comparison, provenance and export. ALK5 activity-ranking methods did not meet recorded acceptance criteria; generated molecules have no measured activity/safety validation. Current additional-source result:42 queries,3 series,2 dependence components; model MAE0.3921 versus nearest analogue0.4305. Prior50/48/44-query results are superseded.

Included evidence is a selected snapshot, not the complete executable tool or full research archive. Analysis scripts retain original data-root assumptions. Full third-party raw workbook/patent images are not redistributed in this ZIP; source links and local archive details appear in the manuscript/reports. Audit source/reuse permissions before public deposit.

Compilation status: built-in preview attempted twice, blocked before LaTeX startup by Windows sandbox CreateProcessWithLogonW267. Static checks passed:11 assets present, no missing references, duplicate labels/manual section numbers or unbalanced environments. These checks do not establish successful compilation or PDF layout.

Before submission: supply verified author/affiliation/disclosure information (author field currently empty); verify retained literature references and public data/code access; compile the full source and review layout. The package is a reviewable manuscript source, not a claim of publication readiness.
'''
with zipfile.ZipFile(destination,'x',compression=zipfile.ZIP_DEFLATED) as z:
    for name,f in items.items(): z.write(f,name)
    z.writestr('README_PACKAGE.txt',readme)
    z.writestr('MANIFEST_SHA256.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(destination) as z:
    assert z.testzip() is None
    assert z.read(tex.name)==tex.read_bytes()
    for name,data in manifest.items(): assert hashlib.sha256(z.read(name)).hexdigest()==data['sha256']
summary=dict(package=str(destination),bytes=destination.stat().st_size,included_figures=11,included_files=len(items)+2,
    zip_integrity_passed=True,source_and_assets_hash_verified=True,compile_verified=False)
(root/'reports/manuscript_package_20261019.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
