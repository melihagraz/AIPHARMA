from pathlib import Path
import json,re
from rdkit import Chem,DataStructs
from rdkit.Chem import Descriptors,Crippen,rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

base=Path.home()/'aipharma_external_leads_20261018';out=Path.home()/'aipharma_structure_mapping_20261018'
rows=json.loads((out/'expanded_nodes.json').read_text());weak=next(r for r in rows if r.get('monomer')=='667180')
def props(m):return {'mw':Descriptors.MolWt(m),'logp':Crippen.MolLogP(m),'heavy_atoms':m.GetNumHeavyAtoms()}
wp=props(Chem.MolFromSmiles(weak['smiles']))
g=json.loads((out/'expanded_source_family_scaffold_similarity_040_graph.json').read_text())
weakblock=next(b for b in g['components'] if any(rows[i].get('monomer')=='667180' for i in b['indices']))
generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
outside={r['scaffold'] for i,r in enumerate(rows) if i not in weakblock['indices'] and r['scaffold']}
outsidefp=[generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in sorted(outside)]
aff=next(iter(json.loads((base/'bindingdb_high_cutoff.json').read_text()).values()))['affinities']
enumerator=rdMolStandardize.TautomerEnumerator();uncharger=rdMolStandardize.Uncharger();candidates=[];total=0
for r in aff:
 if r['doi']!='10.7270/Q2G73K10' or r['affinity_type']!='IC50':continue
 match=re.fullmatch(r'\s*([<>=]*)\s*([\d.]+)\s*',r['affinity'])
 if not match or float(match[2])>100 or match[1] not in ['','<','<=','=']:continue
 total+=1;m=Chem.MolFromSmiles(r['smile']);pr=props(m)
 delta={k:abs(pr[k]-wp[k]) for k in pr}
 if delta['mw']>30 or delta['logp']>.75 or delta['heavy_atoms']>2:continue
 m=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m));m=uncharger.uncharge(m);m=enumerator.Canonicalize(m)
 Chem.RemoveStereochemistry(m);smiles=Chem.MolToSmiles(m);scaf=MurckoScaffold.MurckoScaffoldSmiles(mol=m)
 fp=generator.GetFingerprint(Chem.MolFromSmiles(scaf));maxsim=max(DataStructs.BulkTanimotoSimilarity(fp,outsidefp))
 candidates.append({**r,'normalized_parent':smiles,'properties':pr,'differences':delta,'scaffold':scaf,
                    'maximum_scaffold_similarity_to_other_blocks':maxsim,'preserves_040_separation':maxsim<.4,
                    'score':delta['mw']/30+delta['logp']/.75+delta['heavy_atoms']/2})
result={'strong_records_same_source':total,'fixed_caliper_candidates':len(candidates),'candidates':sorted(candidates,key=lambda r:(not r['preserves_040_separation'],r['score'],r['monomerid'])),
        'rules':'Same source family and original fixed MW/logP/heavy atom calipers; no docking scores read; all nominal/upper bound active values <=100nM need source-table validation'}
(out/'missing_active_probe.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
