"""Source-verify a missing active and test a stricter, five-block panel design."""
from pathlib import Path
from urllib.parse import quote
import hashlib,json,re,requests
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem
from rdkit.Chem import Descriptors,Crippen
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from audit_alk5_independent_blocks_20261018 import graph

home=Path.home();base=home/'aipharma_external_leads_20261018';out=home/'aipharma_structure_mapping_20261018'
name='4-[7-(5-chloro-2-fluorophenyl)-1H,2H,3H-pyrido[3,4-b][1,4]oxazin-1-yl]pyridine'
text=(base/'sources/US20240116947A1.txt').read_text()
assert re.sub(r'\s+','',name).lower() in re.sub(r'\s+','',text).lower()
url='https://www.ebi.ac.uk/opsin/ws/'+quote(name,safe='')+'.json'
dest=out/'opsin/667169.json'
if not dest.exists():
 response=requests.get(url,timeout=30);response.raise_for_status();dest.write_bytes(response.content)
body=json.loads(dest.read_text());assert body['status']=='SUCCESS'
enumerator=rdMolStandardize.TautomerEnumerator();uncharger=rdMolStandardize.Uncharger()
def normalize(smi):
 m=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(Chem.MolFromSmiles(smi)));m=uncharger.uncharge(m);m=enumerator.Canonicalize(m);Chem.RemoveStereochemistry(m)
 return Chem.MolToSmiles(m)
probe=json.loads((out/'missing_active_probe.json').read_text())
candidate=next(r for r in probe['candidates'] if r['monomerid']=='667169')
assert normalize(body['smiles'])==candidate['normalized_parent']
# BindingDB's 7.50 is a midpoint of a patent category, NOT a measured exact IC50.
tables=json.loads((base/'patent_tables.json').read_text())['US20240116947A1']
category_rows=[]
for table in tables:
 for row in table:
  cells=[c for c in row if c]
  if len(cells)>=2 and cells[-1]=='5-10' and '1' in re.findall(r'\d+',cells[-2]):category_rows.append(row)
assert category_rows
active={'smiles':candidate['normalized_parent'],'original_smiles':candidate['smile'],'label':1,'documents':['patent-family:71620262'],
        'scaffold':candidate['scaffold'],'origin':'new_source_verified_active_interval','monomer':'667169',
        'reported_IC50_interval_nM':[5,10],'point_value_is_measured':False,'source_compound':'1',
        'source_patent':'US20240116947A1','source_name':name,'opsin_smiles':body['smiles'],'category_row':category_rows[0]}
(out/'additional_active_source_review.json').write_text(json.dumps(active,indent=2))
rows=json.loads((out/'expanded_nodes.json').read_text());assert active['smiles'] not in {r['smiles'] for r in rows};rows.append(active)
mapping=json.loads((out/'mapping_results.json').read_text())
relative={r['monomer'] for r in mapping if r.get('bindingdb_relative_stereo_annotation')}
for r in rows:
 m=Chem.MolFromSmiles(r.get('original_smiles') or r['smiles']);r['mw']=Descriptors.MolWt(m);r['logp']=Crippen.MolLogP(m);r['heavy_atoms']=m.GetNumHeavyAtoms()
 r['pair_eligible']=r.get('monomer') not in relative
g=graph(rows,.4);(out/'completed_panel_source_family_scaffold_similarity_040_graph.json').write_text(json.dumps(g,indent=2))
(out/'completed_panel_nodes.json').write_text(json.dumps(rows,indent=2))
weakblocks=[b for b in g['components'] if b['weak']];activeblocks=[b for b in g['components'] if b['active']]
cost=np.full((len(weakblocks),len(activeblocks)),1e6);options={}
for i,wb in enumerate(weakblocks):
 for j,ab in enumerate(activeblocks):
  if ab['weak'] and ab['block']!=wb['block']:continue
  choices=[]
  for wi in wb['indices']:
   w=rows[wi]
   if w['label']!=0 or not w['pair_eligible']:continue
   for ai in ab['indices']:
    a=rows[ai]
    if a['label']!=1 or not a['pair_eligible']:continue
    delta={k:abs(w[k]-a[k]) for k in ['mw','logp','heavy_atoms']}
    if delta['mw']<=30 and delta['logp']<=.75 and delta['heavy_atoms']<=2:
     choices.append((delta['mw']/30+delta['logp']/.75+delta['heavy_atoms']/2,w['smiles'],a['smiles'],wi,ai,delta))
  if choices:best=min(choices);cost[i,j]=best[0];options[i,j]=best
ii,jj=linear_sum_assignment(cost);pairs=[]
for i,j in zip(ii,jj):
 if cost[i,j]>=1e6:continue
 _,_,_,wi,ai,delta=options[i,j]
 pairs.append({'block':len(pairs)+1,'source_components':sorted({weakblocks[i]['block'],activeblocks[j]['block']}),
               'weak':rows[wi],'active':rows[ai],'property_differences':delta})
summary={'source_name_parent_verified':True,'active_interval_nM':[5,10],'discarded_midpoint_as_exact_measurement':True,
         'maximum_matched_blocks_after_additional_active':len(pairs),'threshold_scaffold_tanimoto':.4,
         'relative_stereo_new_records_excluded_from_pair_selection':sorted(relative),'pairs':pairs,
         'limits':'Previously examined old data; properties-only matching and connected-component separation are not biological independence proof. No docking scores used.'}
assert len({c for r in pairs for c in r['source_components']})==sum(len(r['source_components']) for r in pairs)
(out/'five_block_panel_feasibility.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'matched_blocks':len(pairs),'pairs':[{'block':r['block'],'weak_monomer':r['weak'].get('monomer'),'active_monomer':r['active'].get('monomer'),'weak_origin':r['weak']['origin'],'active_origin':r['active']['origin']} for r in pairs]},indent=2))
