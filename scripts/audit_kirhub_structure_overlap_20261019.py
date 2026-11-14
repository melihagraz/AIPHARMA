"""Audit name-resolved public structures; never use endpoint values to select structures."""
import json,sys,hashlib
from pathlib import Path
from collections import Counter
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
p=Path(sys.argv[1]); k=p/'kirhub'; resolved=json.loads((k/'structure_resolution/resolved_names.json').read_text())
byname={r['name']:r for r in resolved}
prior=json.loads((p/'prior_measured_records.json').read_text())
training=json.loads((p/'curated_complete_identity/records.json').read_text())
old_ids={r['normalized_parent'] for r in prior}; old_scaf={r['scaffold'] for r in prior}
new_ids={r['normalized_parent'] for r in training}; new_scaf={r['scaffold'] for r in training}
enum=rdMolStandardize.TautomerEnumerator(); un=rdMolStandardize.Uncharger()
records=[]
for item in json.loads((k/'reserved_alk5_single_dose.json').read_text()):
    name=item['compound']; match=byname[name]; row=dict(compound=name,source_name_cell=item['source_name_cell'],resolution_url=match['url'])
    try:
        assert match['status']=='unique_name_match_pending_identity_audit'
        prop=match['raw_response']['PropertyTable']['Properties'][0]
        smiles=prop.get('SMILES') or prop.get('IsomericSMILES')
        mol=Chem.MolFromSmiles(smiles); assert mol is not None
        row.update(pubchem_cid=prop['CID'],pubchem_inchikey=prop['InChIKey'],source_smiles=smiles)
        assert not any(str(s.specified)=='Unspecified' for s in Chem.FindPotentialStereo(mol)), 'Unresolved stereochemistry'
        standardized=un.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol)))
        options=enum.Enumerate(standardized); assert str(options.status)=='Completed', 'Incomplete tautomer enumeration'
        normalized=enum.PickCanonical(options); Chem.RemoveStereochemistry(normalized)
        identity=Chem.MolToSmiles(normalized); scaffold=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(normalized))
        row.update(status='normalized_for_overlap_audit',normalized_parent=identity,scaffold=scaffold,
            prior_model_identity_overlap=identity in old_ids,prior_model_scaffold_overlap=scaffold in old_scaf,
            expanded_data_identity_overlap=identity in new_ids,expanded_data_scaffold_overlap=scaffold in new_scaf)
    except Exception as exc: row.update(status='quarantined',reason=str(exc))
    records.append(row)
valid=[r for r in records if r['status']=='normalized_for_overlap_audit']
testids={r['normalized_parent'] for r in valid}; testscaf={r['scaffold'] for r in valid}
identity_clean=[r for r in training if r['normalized_parent'] not in testids]
scaffold_clean=[r for r in training if r['normalized_parent'] not in testids and r['scaffold'] not in testscaf]
out=k/'overlap_audit'; out.mkdir(exist_ok=True)
for name,data in [('structure_records.json',records),('identity_excluded_training.json',identity_clean),('scaffold_excluded_training.json',scaffold_clean)]:
    (out/name).write_text(json.dumps(data,indent=2,allow_nan=False))
summary=dict(screen_drugs=len(records),normalized=len(valid),quarantined=[dict(compound=r['compound'],reason=r['reason']) for r in records if r['status']=='quarantined'],
    overlaps={key:sum(r[key] for r in valid) for key in ['prior_model_identity_overlap','prior_model_scaffold_overlap','expanded_data_identity_overlap','expanded_data_scaffold_overlap']},
    identity_excluded_training_groups=len(identity_clean),scaffold_excluded_training_groups=len(scaffold_clean),
    scaffold_excluded_training_exact=sum(r['median_pic50'] is not None for r in scaffold_clean),
    interpretation='PubChem unique name matches with chemical normalization, not independent confirmation of vial identity. Earlier models with overlapping compounds cannot be called independent on those compounds.',
    models_fitted=False,source_endpoint_values_not_used_for_structure_selection=True,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(out/'summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
