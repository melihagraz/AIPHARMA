"""Feasibility only: can full dependency components supply five matched blocks?"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem
from rdkit.Chem import Descriptors,Crippen

p=Path.home()/'aipharma_structure_mapping_20261018'
rows=json.loads((p/'expanded_nodes.json').read_text())
for r in rows:
 m=Chem.MolFromSmiles(r['smiles']);r['mw']=Descriptors.MolWt(m);r['logp']=Crippen.MolLogP(m);r['heavy_atoms']=m.GetNumHeavyAtoms()
results={}
for mode in ['source_family_exact_scaffold','source_family_scaffold_similarity_040']:
 g=json.loads((p/f'expanded_{mode}_graph.json').read_text());blocks=g['components']
 weakblocks=[b for b in blocks if b['weak']];activeblocks=[b for b in blocks if b['active']]
 cost=np.full((len(weakblocks),len(activeblocks)),1e6);pairs={}
 for i,wb in enumerate(weakblocks):
  for j,ab in enumerate(activeblocks):
   # Never merge two weak-containing components or split a full component.
   if ab['weak'] and ab['block']!=wb['block']:continue
   choices=[]
   for wi in wb['indices']:
    if rows[wi]['label']!=0:continue
    for ai in ab['indices']:
     if rows[ai]['label']!=1:continue
     w,a=rows[wi],rows[ai]
     delta={k:abs(w[k]-a[k]) for k in ['mw','logp','heavy_atoms']}
     if delta['mw']<=30 and delta['logp']<=.75 and delta['heavy_atoms']<=2:
      score=delta['mw']/30+delta['logp']/.75+delta['heavy_atoms']/2
      choices.append((score,w['smiles'],a['smiles'],wi,ai,delta))
   if choices:
    best=min(choices);cost[i,j]=best[0];pairs[i,j]=best
 wi,ai=linear_sum_assignment(cost);selected=[]
 for i,j in zip(wi,ai):
  if cost[i,j]>=1e6:continue
  score,_,_,wn,an,delta=pairs[i,j]
  selected.append({'weak_source_component':weakblocks[i]['block'],'active_source_component':activeblocks[j]['block'],
                   'weak_node':wn,'active_node':an,'property_differences':delta,
                   'weak_origin':rows[wn]['origin'],'weak_monomer':rows[wn].get('monomer'),
                   'weak_smiles':rows[wn]['smiles'],'active_smiles':rows[an]['smiles']})
 results[mode]={'weak_components':len(weakblocks),'active_components':len(activeblocks),
                'feasible_disjoint_full_component_blocks':len(selected),'five_block_design_feasible':len(selected)>=5,
                'selected_feasibility_examples':selected,'calipers':{'mw':30,'logp':.75,'heavy_atoms':2},
                'limits':'Maximum matching under full-component preservation and fixed property calipers, not a proof over all possible molecular subsetting schemes; no docking or prospective validation.'}
(p/'block_matching_feasibility.json').write_text(json.dumps(results,indent=2))
print(json.dumps({k:{a:b for a,b in v.items() if a!='selected_feasibility_examples'} for k,v in results.items()},indent=2))
