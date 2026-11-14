"""Prepare original linked drawings and source-independent RDKit depictions for review."""
import hashlib,json,sys,urllib.request
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
p=Path(sys.argv[1]); evidence=p/'primary_audit/structure_evidence'
out=evidence/'drawing_review_11335'; out.mkdir(exist_ok=False)
valid=json.loads((p/'mass_conflict_excluded/curated_complete_identity/records.json').read_text())
lookup={(r['assay_id'],r['normalized_parent']):r for r in json.loads((evidence/'primary_structure_evidence.json').read_text())}
manifest=[]
for r in valid:
    if r['assay_id']!='BindingDB11335_assay1': continue
    e=lookup[r['assay_id'],r['normalized_parent']]['primary_evidence'][0]
    number=e['compound_number']; urls=e.get('characterization_image_urls',[])
    kind='characterization_table'
    if not urls:
        urls=[url for b in e['synthesis_blocks'] for url in b['image_urls']]; kind='synthesis_route'
    row=dict(compound_number=number,normalized_parent=r['normalized_parent'],source_kind=kind,source_urls=urls,
             visual_reviewed=False,full_structure_verified=False)
    if urls:
        source=out/(number+'_patent.png')
        try:
            with urllib.request.urlopen(urls[0],timeout=45) as response: raw=response.read()
            assert raw.startswith(b'\x89PNG\r\n\x1a\n')
            source.write_bytes(raw); row.update(source_file=source.name,source_sha256=hashlib.sha256(raw).hexdigest())
        except Exception as exc: row['download_error']=str(exc)
    Draw.MolToFile(Chem.MolFromSmiles(r['normalized_parent']),str(out/(number+'_reference.png')),size=(600,350))
    row['reference_file']=number+'_reference.png'; manifest.append(row)
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)); print(json.dumps([dict(compound=r['compound_number'],kind=r['source_kind'],downloaded='source_file' in r) for r in manifest],indent=2))
