"""Prepare three archived catalytic BindingDB entries for provisional diagnostics."""
import hashlib,json,math,re,statistics,sys
from pathlib import Path
from collections import Counter,defaultdict
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
p=Path(sys.argv[1]); out=p/'additional_local_series_complete'; out.mkdir(exist_ok=False)
def save(n,x): (out/n).write_text(json.dumps(x,indent=2,allow_nan=False))
entries={10830:'US11447490',10870:'US11466003',11335:'US11691981'}
save('protocol_source.json',dict(entries=entries,endpoint='Archived assay1 descriptions explicitly catalytic ALK5 ATP/peptide ADP-Glo',
 exclusions='Only positive exact nM IC50, complete normalization and stereo; remove development identity and exact-scaffold overlaps before eligibility; retain exclusions',
 method='Existing evaluator unchanged, lambda0.25, same assay/scaffold >=8; no tuning or outcome-based series selection',
 sources='Previously archived/inspected and possibly prior source families; not untouched confirmation. BindingDB numbers not independently crosschecked to these patent tables.',
 conservative_dependence='Group US11447490 and US11466003 together as same protocol/source family pending exhaustive genealogy; do not count as two independent sources',
 acceptance='Provisional diagnostics cannot promote a model regardless of numerical performance',
 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
dev=json.loads((p/'curated_complete_identity/records.json').read_text()); ids={r['normalized_parent'] for r in dev}; scaffolds={r['scaffold'] for r in dev}
archive=json.loads((p/'entry_archive.json').read_text()); en=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger()
grouped=defaultdict(list); excluded=[]
for entry in archive:
    eid=entry['entry_id']
    if eid not in entries: continue
    method=entry['assays']['1']['text']; assert 'ATP' in method and ('ADP-Glo' in method or 'ADP-GLO' in method)
    assert 'T204D' not in method
    for r in entry['alk5_records']:
        try:
            assert r['measurement_type']=='IC50'
            m=re.fullmatch(r'([0-9.Ee+\-]+)nM',r['reported_value']); assert m, 'Nonexact/unsupported value'
            val=float(m[1]); assert val>0 and math.isfinite(val)
            mol=Chem.MolFromSmiles(r['smiles']); assert mol is not None
            assert not any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(mol)), 'Unspecified stereo'
            mol=un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol)))
            ts=en.Enumerate(mol); assert str(ts.status)=='Completed', 'Incomplete tautomer enumeration'
            mol=en.PickCanonical(ts); Chem.RemoveStereochemistry(mol); identity=Chem.MolToSmiles(mol); scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
            assert identity not in ids and scaffold not in scaffolds, 'Development identity/scaffold overlap'
            grouped[eid,identity].append(dict(**r,normalized_parent=identity,scaffold=scaffold,pic50=9-math.log10(val)))
        except Exception as exc: excluded.append(dict(entry=eid,monomer_id=r['monomer_id'],reported_value=r['reported_value'],reason=str(exc)))
rows=[]
for (eid,identity),rs in sorted(grouped.items()):
    ys=[r['pic50'] for r in rs]
    if max(ys)-min(ys)>=2:
        excluded.append(dict(entry=eid,reason='Alias potency spread >=2 logs',normalized_parent=identity)); continue
    family='conservative_family_10830_10870' if eid in [10830,10870] else entries[eid]
    rows.append(dict(assay_id=f'BindingDB{eid}_assay1',normalized_parent=identity,scaffold=rs[0]['scaffold'],smiles=identity,
        documents=[family],median_pic50=statistics.median(ys),source_records=rs,source_patent=entries[eid],numeric_primary_crosscheck=False))
(out/'curated_complete_identity').mkdir(); save('curated_complete_identity/records.json',rows); save('exclusions.json',excluded)
counts=Counter((r['assay_id'],r['scaffold']) for r in rows)
summary=dict(normalized_assay_identity_groups=len(rows),unique_identities=len({r['normalized_parent'] for r in rows}),
    entry_counts=dict(Counter(r['assay_id'] for r in rows)),exclusion_counts=dict(Counter(r['reason'] for r in excluded)),
    eligible_series=sum(n>=8 for n in counts.values()),eligible_queries=sum(n for n in counts.values() if n>=8),provisional_not_confirmatory=True)
save('curation_summary.json',summary); print(json.dumps(summary,indent=2))
