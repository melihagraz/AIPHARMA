"""Mass consistency check; not atom/bond reconstruction."""
import json,re,sys
from pathlib import Path
from collections import Counter
from rdkit import Chem
from rdkit.Chem import Descriptors,Draw
p=Path(sys.argv[1]); rows=json.loads((p/'primary_structure_evidence.json').read_text()); output=[]; representatives=[]; seen=set()
for r in rows:
    mass=Descriptors.ExactMolWt(Chem.MolFromSmiles(r['normalized_parent']))+1.007276466621
    checks=[]; urls=[]
    for e in r['primary_evidence']:
        values=[]
        for block in e['synthesis_blocks']:
            text=block['text']; urls+=block['image_urls']
            found=re.findall(r'm/z\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*\[M\s*\+\s*H\]',text,re.I)
            if not found: found=re.findall(r'\[M\s*\+\s*H\]\s*\+?\s*=\s*([0-9]+(?:\.[0-9]+)?)',text,re.I)
            # Only first product MS after the exact synthesis heading; later text can describe subsequent examples.
            if found: values.append(float(found[0]))
        for row in e['characterization_table_rows']:
            if len(row)>2 and re.fullmatch(r'\d+(?:\.\d+)?',row[2]): values.append(float(row[2]))
        urls+=e.get('characterization_image_urls',[])
        checks.append(dict(patent=e['patent'],compound_number=e['compound_number'],observed_MH_masses=values,
            status='no_unambiguous_MH_mass_extracted' if not values else ('consistent_within_0.5_Da' if all(abs(v-mass)<=.5 for v in values) else 'mass_discrepancy')))
    output.append(dict(assay_id=r['assay_id'],normalized_parent=r['normalized_parent'],calculated_MH=mass,checks=checks,full_structure_verified=False))
    if r['assay_id'] not in seen and urls:
        seen.add(r['assay_id']); tag=r['assay_id']; Draw.MolToFile(Chem.MolFromSmiles(r['normalized_parent']),str(p/(tag+'_reference.png')),size=(700,420))
        representatives.append(dict(assay_id=tag,normalized_parent=r['normalized_parent'],evidence=r['primary_evidence'],drawing_url=urls[0],reference_file=tag+'_reference.png'))
(p/'mass_checks.json').write_text(json.dumps(output,indent=2)); (p/'representatives.json').write_text(json.dumps(representatives,indent=2))
summary=dict(molecules=len(rows),alias_checks=dict(Counter(c['status'] for r in output for c in r['checks'])),representatives=len(representatives),full_structure_verified=False)
(p/'mass_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
