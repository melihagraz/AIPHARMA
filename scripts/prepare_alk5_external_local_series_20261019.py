"""Curate an archived separate patent for the frozen multiple-reference task."""
import hashlib,json,math,statistics,sys
from pathlib import Path
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
p=Path(sys.argv[1]); source=p/'external_source'; out=p/'external_local_transfer'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
save('source_protocol.json',dict(source='US20250179094A1, BindingDB entry12979 assay1',
 endpoint='Direct catalytic ADP-Glo TGFBR1tide phosphorylation IC50; no explicit mutation in archived method, not sequence confirmation',
 task='Freeze prior within-series Morgan2/2048 Tanimoto KRR lambda0.25; same-assay same-Murcko >=8 identities; leave one out, no tuning',
 independence='Separate from the 684-row development data by identity/scaffold audit; previously inspected and used for different single-reference model tests, not untouched or blind',
 exclusions='Positive control LY364947, unspecified stereo, incomplete normalization; normalized aliases median with ranges preserved; inconsistent >=2 log aliases quarantined',
 interpretation='One external source cannot satisfy five-component acceptance requirement; report diagnostic generalization only',
 files_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in source.glob('*.html')}))
enum=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger(); grouped=defaultdict(list); exclusions=[]
for r in json.loads((source/'parsed_series.json').read_text()):
    try:
        assert r.get('primary_value_verified'), 'Primary value not verified'
        assert r['primary_table_key']!='LY364947', 'Positive control'
        mol=Chem.MolFromSmiles(r['smiles']); assert mol is not None
        assert not any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(mol)), 'Unspecified stereo'
        mol=un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol)))
        options=enum.Enumerate(mol); assert str(options.status)=='Completed', 'Incomplete tautomer enumeration'
        mol=enum.PickCanonical(options); Chem.RemoveStereochemistry(mol); ident=Chem.MolToSmiles(mol)
        grouped[ident].append(dict(**r,scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol)),pic50=9-math.log10(r['ic50_nM'])))
    except Exception as exc: exclusions.append(dict(monomer_id=r['monomer_id'],reason=str(exc)))
rows=[]
for ident,items in sorted(grouped.items()):
    values=[r['pic50'] for r in items]
    if max(values)-min(values)>=2:
        exclusions.append(dict(normalized_parent=ident,reason='Alias potency spread >=2 logs',source_records=items)); continue
    rows.append(dict(assay_id='BindingDB12979_assay1',normalized_parent=ident,smiles=ident,scaffold=items[0]['scaffold'],
        documents=['US20250179094A1'],median_pic50=statistics.median(values),source_pic50_range=[min(values),max(values)],source_records=items))
development=json.loads((p/'development_records.json').read_text())
id_overlap=sorted({r['normalized_parent'] for r in rows}&{r['normalized_parent'] for r in development})
scaf_overlap=sorted({r['scaffold'] for r in rows}&{r['scaffold'] for r in development})
save('overlap.json',dict(identity_overlap=id_overlap,scaffold_overlap=scaf_overlap))
assert not id_overlap and not scaf_overlap, 'External source overlaps current development chemistry'
counts=defaultdict(int)
for r in rows: counts[r['scaffold']]+=1
save('exclusions.json',exclusions); (out/'curated_complete_identity').mkdir(); save('curated_complete_identity/records.json',rows)
summary=dict(input_records=64,normalized_identities=len(rows),excluded_records_or_alias_groups=len(exclusions),
    scaffold_sizes=sorted(counts.values(),reverse=True),eligible_series=sum(n>=8 for n in counts.values()),eligible_queries=sum(n for n in counts.values() if n>=8),
    development_identity_overlap=len(id_overlap),development_scaffold_overlap=len(scaf_overlap),previously_inspected=True,source_not_blind=True)
save('curation_summary.json',summary); print(json.dumps(summary,indent=2))
