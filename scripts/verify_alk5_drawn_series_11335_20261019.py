"""Manual source-drawing transcription compared with archived molecular graphs.

The analyst has seen references: this is an auditable crosscheck, not blinded
independent reconstruction. Only the ten explicitly reviewed drawings apply.
"""
import hashlib,json,sys
from pathlib import Path
from rdkit import Chem
p=Path(sys.argv[1]); manifest=json.loads((p/'manifest.json').read_text())
# Source drawings: common tertiary-alcohol/aminopyridine/oxypyrazole-CF3 core.
# Phenyl positions are counted from its bond to the pyrazole nitrogen.
aryl={
 '17':('c1ccccc1','unsubstituted phenyl'),
 '34':('c1ccc(C)cc1','4-methylphenyl'),
 '35':('c1ccc(F)cc1','4-fluorophenyl'),
 '36':('c1cc(Cl)ccc1','3-chlorophenyl'),
 '38':('c1ccccc1F','2-fluorophenyl'),
 '40':('c1ccccc1Cl','2-chlorophenyl'),
 '43':('c1ccc(Cl)cc1','4-chlorophenyl'),
 '44':('c1ccc(C#N)cc1','4-cyanophenyl, final product44 in route, not precursor43'),
 '45':('c1c(F)cccc1F','2,6-difluorophenyl'),
 '51':('c1ccccc1C#N','2-cyanophenyl')}
results=[]
for r in manifest:
    number=r['compound_number']; group,note=aryl[number]
    group=group.replace('1','4')  # Outer pyridine ring1 is still open at this branch.
    transcribed='OC(C)(C)c1nccc(Nc2nccc(Oc3cc(C(F)(F)F)nn3-'+group+')c2)c1'
    mol=Chem.MolFromSmiles(transcribed); reference=Chem.MolFromSmiles(r['normalized_parent'])
    assert mol is not None and reference is not None
    drawn=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
    archived=Chem.MolToSmiles(reference,canonical=True,isomericSmiles=True)
    results.append(dict(compound_number=number,source_url=r['source_urls'][0],source_sha256=r['source_sha256'],
        manually_transcribed_smiles=transcribed,source_substituent_observation=note,
        canonical_transcription=drawn,canonical_archive=archived,graphs_match=drawn==archived,
        visual_reviewed=True,scope='Constitutional graph, explicit isotope/stereo labels only; no experimental sample identity verification',
        source_kind=r['source_kind']))
summary=dict(reviewed=len(results),graph_matches=sum(r['graphs_match'] for r in results),
    blind_or_independent_analyst=False,automatic_image_recognition=False,
    assay='BindingDB11335_assay1 / US11691981B2',other_series_review_complete=False,
    biology_or_safety_validated=False,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(p/'graph_crosschecks.json').write_text(json.dumps(results,indent=2)); (p/'graph_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
assert all(r['graphs_match'] for r in results), 'Source-transcribed graph mismatch'
