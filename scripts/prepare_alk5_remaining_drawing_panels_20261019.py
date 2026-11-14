"""Read-only scientific comparison panels; original full route images preserved."""
import hashlib,json,sys,urllib.request,concurrent.futures
from pathlib import Path
from PIL import Image,ImageDraw,ImageOps
from rdkit import Chem
from rdkit.Chem import Draw
p=Path(sys.argv[1]); evidence=p/'primary_audit/structure_evidence'; out=evidence/'remaining_drawing_review'; out.mkdir(exist_ok=False)
valid=json.loads((p/'mass_conflict_excluded/curated_complete_identity/records.json').read_text())
lookup={(r['assay_id'],r['normalized_parent']):r for r in json.loads((evidence/'primary_structure_evidence.json').read_text())}
manifest=[]
for r in valid:
    if r['assay_id']=='BindingDB11335_assay1': continue
    es=lookup[r['assay_id'],r['normalized_parent']]['primary_evidence']; urls=[u for e in es for b in e['synthesis_blocks'] for u in b['image_urls']]
    number=es[0]['compound_number']; tag=r['assay_id']+'_'+number
    row=dict(tag=tag,assay_id=r['assay_id'],compound_numbers=[e['compound_number'] for e in es],normalized_parent=r['normalized_parent'],source_urls=sorted(set(urls)),
             source_file=tag+'_source.png',reference_file=tag+'_reference.png',visual_status='pending')
    Draw.MolToFile(Chem.MolFromSmiles(r['normalized_parent']),str(out/row['reference_file']),size=(600,350)); manifest.append(row)
def fetch(r):
    if not r['source_urls']: r['download_status']='no_route_image_located'; return r
    try:
        with urllib.request.urlopen(r['source_urls'][0],timeout=40) as response: raw=response.read()
        assert raw.startswith(b'\x89PNG\r\n\x1a\n'); (out/r['source_file']).write_bytes(raw)
        r.update(download_status='downloaded',sha256=hashlib.sha256(raw).hexdigest())
    except Exception as e: r.update(download_status='failed',error=str(e))
    return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool: manifest=list(pool.map(fetch,manifest))
for page,start in enumerate(range(0,len(manifest),4),1):
    subset=manifest[start:start+4]; canvas=Image.new('RGB',(1200,400*len(subset)),'white'); pen=ImageDraw.Draw(canvas)
    for j,r in enumerate(subset):
        y=j*400; pen.text((15,y+8),r['tag']+' | patent final-route viewport (left), archived structure (right)',fill='black')
        if r['download_status']=='downloaded':
            original=Image.open(out/r['source_file']).convert('RGB'); w,h=original.size
            crop=(0,max(0,h-max(500,int(w*.75))),w,h); r['panel_source_crop_box']=crop; r['original_dimensions']=[w,h]
            view=ImageOps.contain(original.crop(crop),(580,365)); canvas.paste(view,((600-view.width)//2,y+30))
        else: pen.text((15,y+45),r['download_status'],fill='red')
        canvas.paste(Image.open(out/r['reference_file']).convert('RGB'),(600,y+30))
        r['panel']='panel_'+str(page)+'.png'
    canvas.save(out/('panel_'+str(page)+'.png'))
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)); print(json.dumps(dict(records=len(manifest),downloaded=sum(r['download_status']=='downloaded' for r in manifest),missing=[r['tag'] for r in manifest if r['download_status']!='downloaded'],panels=(len(manifest)+3)//4),indent=2))
