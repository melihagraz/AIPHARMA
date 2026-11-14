"""Fixed, property-matched enzyme control challenge; exploratory, not blind validation."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def properties(r):
    m = Chem.MolFromSmiles(r['smiles'])
    return {'mw': Descriptors.MolWt(m), 'logp': Crippen.MolLogP(m), 'heavy_atoms': m.GetNumHeavyAtoms()}


def prepare(source, audit, out):
    if (out/'protocol.json').exists(): raise RuntimeError('Protocol already locked')
    old = {r['smiles'] for r in json.loads((source/'locked_items.json').read_text())}
    dataset = json.loads((audit/'datasets_all.json').read_text())['enzyme_catalytic__no_explicit_mutant']
    pool = []
    for r in dataset:
        m = Chem.MolFromSmiles(r['smiles'])
        if r['evidence'] != 'exact' or r['smiles'] in old or r['old_holdout_scaffold']: continue
        if any(str(x.specified) == 'Unspecified' for x in Chem.FindPotentialStereo(m)): continue
        pool.append({**r, **properties(r)})
    active = sorted([r for r in pool if r['label'] == 1], key=lambda r:r['smiles'])
    weak = sorted([r for r in pool if r['label'] == 0], key=lambda r:r['smiles'])
    cost = np.full((len(weak),len(active)), 1e6)
    for i,w in enumerate(weak):
        for j,a in enumerate(active):
            dm, dl, dh = abs(w['mw']-a['mw']), abs(w['logp']-a['logp']), abs(w['heavy_atoms']-a['heavy_atoms'])
            if dm<=30 and dl<=.75 and dh<=2: cost[i,j]=dm/30+dl/.75+dh/2
    wi, ai = linear_sum_assignment(cost)
    pairs = [(weak[i],active[j]) for i,j in zip(wi,ai) if cost[i,j]<1e6]
    out.mkdir(parents=True,exist_ok=True)
    key, structures = [], []
    for k,(w,a) in enumerate(pairs,1):
        for r in [w,a]:
            ident = f'item_{len(key)+1:02d}'
            key.append({**r,'id':ident,'pair':k})
            structures.append({'id':ident,'smiles':r['smiles']})
    # Opaque ordering prevents label-based behavior in docking code; labels are not truly blinded.
    structures.sort(key=lambda r:hashlib.sha256(r['smiles'].encode()).hexdigest())
    (out/'structures.json').write_text(json.dumps(structures,indent=2))
    (out/'label_key.json').write_text(json.dumps(key,indent=2))
    protocol={'scope':'Exploratory new scored subset from previously examined data, not independent external or third-party blind validation',
              'selection':'All eligible exact weak controls, one-to-one global minimum property distance; no score selection',
              'calipers':{'MW_absolute_max':30,'logP_absolute_max':.75,'heavy_atoms_absolute_max':2},
              'exclusions':'Previous nine ensemble molecules, old holdout scaffolds, unspecified stereochemistry, censored measurements',
              'labels':'Median exact catalytic pChEMBL >=7 active; <=5 weak; no explicit mutant; descriptions only',
              'weak_pool':len(weak),'active_pool':len(active),'matched_pairs':len(pairs),
              'scoring':{'receptor':'1PY5','functions':['vina','vinardo'],'seed':20261018,'exhaustiveness':8,'cpu':1,'num_poses':9},
              'preparation':{'ETKDGv3_seed':20261018,'conformers':16,'MMFF94s_iterations':1000,'converged_only':True,'Meeko_rigid_macrocycles':True},
              'metrics':['ROC_AUC','within_pair_correct_fraction','paired_label_swap_test','heavy_atom_normalization','smaller_size_baseline'],
              'promotion_rule':'No production promotion from this small exploratory set',
              'structures_sha256':sha(out/'structures.json'),'label_key_sha256':sha(out/'label_key.json'),
              'dataset_sha256':sha(audit/'datasets_all.json'),'receptor_sha256':sha(source/'1PY5/receptor.pdbqt')}
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2))
    print(json.dumps(protocol),flush=True)


def dock(task):
    out, source, item, sf = task
    w=out/'jobs'/f"{item['id']}_{sf}";w.mkdir(parents=True,exist_ok=True)
    result={'id':item['id'],'function':sf}
    try:
        mol=Chem.AddHs(Chem.MolFromSmiles(item['smiles']))
        ep=AllChem.ETKDGv3();ep.randomSeed=20261018;ep.numThreads=1;ep.pruneRmsThresh=.25
        ids=list(AllChem.EmbedMultipleConfs(mol,16,ep))
        opts=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000,mmffVariant='MMFF94s')
        good=[(o[1],cid) for cid,o in zip(ids,opts) if o[0]==0]
        if not good: raise RuntimeError('No converged conformer')
        cid=min(good)[1];m=Chem.Mol(mol);conf=Chem.Conformer(mol.GetConformer(cid));m.RemoveAllConformers();m.AddConformer(conf)
        writer=Chem.SDWriter(str(w/'prepared.sdf'));writer.write(m);writer.close()
        setup=MoleculePreparation(rigid_macrocycles=True).prepare(m)[0]
        text,ok,err=PDBQTWriterLegacy.write_string(setup)
        if not ok: raise RuntimeError(str(err))
        (w/'ligand.pdbqt').write_text(text)
        box=json.loads((source/'1PY5/box.json').read_text())
        v=Vina(sf_name=sf,cpu=1,seed=20261018,verbosity=0)
        v.set_receptor(str(source/'1PY5/receptor.pdbqt'));v.set_ligand_from_string(text)
        v.compute_vina_maps(center=box['center'],box_size=box['size'])
        v.dock(exhaustiveness=8,n_poses=9)
        energy=v.energies(n_poses=9)
        v.write_poses(str(w/'poses.pdbqt'),n_poses=9,overwrite=True)
        result.update({'status':'complete','score':float(energy[0,0]),'energies':energy.tolist()})
    except Exception as e: result.update({'status':'failed','error':str(e)})
    (w/'result.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
    return result


def run(source,out):
    protocol=json.loads((out/'protocol.json').read_text())
    assert sha(out/'structures.json')==protocol['structures_sha256']
    assert sha(source/'1PY5/receptor.pdbqt')==protocol['receptor_sha256']
    if (out/'results.json').exists(): raise RuntimeError('Completed run is immutable')
    structures=json.loads((out/'structures.json').read_text())
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(dock,[(out,source,r,sf) for r in structures for sf in ['vina','vinardo']]))
    (out/'results.json').write_text(json.dumps(results,indent=2))
    print('COMPLETED',len(results),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','run']);p.add_argument('--source',type=Path,required=True);p.add_argument('--audit',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.mode=='prepare':prepare(args.source,args.audit,args.output)
    else:run(args.source,args.output)
