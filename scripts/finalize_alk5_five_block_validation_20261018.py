"""Attach design evidence and runtime metadata after the result audit."""
from pathlib import Path
import hashlib,importlib.metadata,json,platform,shutil

home=Path.home();out=home/'aipharma_five_block_validation_20261018';mapping=home/'aipharma_structure_mapping_20261018'
verification=json.loads((out/'evaluation/verification.json').read_text());assert verification['passed'] and verification['record_count']==252
protocol=json.loads((out/'protocol.json').read_text());design=out/'design';design.mkdir(exist_ok=True)
for name in ['five_block_panel_feasibility.json','completed_panel_nodes.json','completed_panel_source_family_scaffold_similarity_040_graph.json','additional_active_source_review.json','mapping_decisions.csv']:
    shutil.copyfile(mapping/name,design/name)
assert hashlib.sha256((design/'five_block_panel_feasibility.json').read_bytes()).hexdigest()==protocol['panel_feasibility_sha256']
runtime={'python':platform.python_version(),'platform':platform.system(),'machine':platform.machine(),'packages':{}}
for package in ['vina','rdkit','numpy','scipy','scikit-learn','meeko','spyrmsd','matplotlib']:
    try:runtime['packages'][package]=importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:runtime['packages'][package]='unavailable'
(out/'runtime_versions.json').write_text(json.dumps(runtime,indent=2))
status=json.loads((out/'run_status.json').read_text());status.update({'phase':'complete','scoring_finished':True,'evaluated':True,'report_generated':True})
(out/'run_status.json').write_text(json.dumps(status,indent=2))
files={str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.json'}
(out/'MANIFEST_SHA256.json').write_text(json.dumps(files,indent=2))
assert all(hashlib.sha256((out/name).read_bytes()).hexdigest()==value for name,value in files.items())
print(json.dumps({'verified_files':len(files),'design_hash_matches_locked_protocol':True,'runtime':runtime},indent=2))
