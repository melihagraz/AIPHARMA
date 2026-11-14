"""Execute the immutable 13-molecule, 252-job follow-up protocol.

Scoring reads structures, not the label key. Labels were previously exposed;
this separation is an implementation safeguard, not independent blinding.
"""
from pathlib import Path
import argparse
import concurrent.futures
import hashlib
import json
import shutil
import time
import traceback
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation,PDBQTWriterLegacy
from spyrmsd import rmsd
from vina import Vina
import vina


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def coordinates(template,pose):
    smi=next(line[len('REMARK SMILES '):] for line in template.splitlines() if line.startswith('REMARK SMILES ') and not line.startswith('REMARK SMILES IDX'))
    mol=Chem.MolFromSmiles(smi);mapping=[]
    for line in template.splitlines():
        if line.startswith('REMARK SMILES IDX '):mapping.extend(map(int,line[len('REMARK SMILES IDX '):].split()))
    positions={int(line[6:11]):[float(line[30:38]),float(line[38:46]),float(line[46:54])] for line in pose.split('ENDMDL',1)[0].splitlines() if line.startswith(('ATOM','HETATM'))}
    xyz=np.full((mol.GetNumAtoms(),3),np.nan)
    for i,j in zip(mapping[::2],mapping[1::2]):xyz[i-1]=positions[j]
    assert np.isfinite(xyz).all()
    return mol,xyz


def prepare(item,out):
    folder=out/'ligands'/item['id'];folder.mkdir(parents=True,exist_ok=True)
    result={'id':item['id'],'smiles':item['smiles']}
    try:
        mol=Chem.AddHs(Chem.MolFromSmiles(item['smiles']))
        ep=AllChem.ETKDGv3();ep.randomSeed=20261018;ep.numThreads=1;ep.pruneRmsThresh=.25
        ids=list(AllChem.EmbedMultipleConfs(mol,16,ep))
        opts=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=1000,mmffVariant='MMFF94s')
        good=[(opt[1],cid) for cid,opt in zip(ids,opts) if opt[0]==0]
        result.update({'embedded':len(ids),'converged':len(good),'MMFF_results':[list(o) for o in opts]})
        if not good:raise RuntimeError('No converged conformer; no fallback')
        energy,cid=min(good);prepared=Chem.Mol(mol);conf=Chem.Conformer(mol.GetConformer(cid));prepared.RemoveAllConformers();prepared.AddConformer(conf)
        writer=Chem.SDWriter(str(folder/'prepared.sdf'));writer.write(prepared);writer.close()
        text,ok,error=PDBQTWriterLegacy.write_string(MoleculePreparation(rigid_macrocycles=True).prepare(prepared)[0])
        if not ok:raise RuntimeError(str(error))
        (folder/'ligand.pdbqt').write_text(text)
        result.update({'status':'complete','selected_conformer':cid,'MMFF_energy':energy,'ligand_sha256':sha(folder/'ligand.pdbqt')})
    except Exception as exc:
        result.update({'status':'failed','error':str(exc)});(folder/'error.log').write_text(traceback.format_exc())
    (folder/'preparation.json').write_text(json.dumps(result,indent=2))
    return result


def dock(task):
    out,item,pdb,sf,seed=task;folder=out/'jobs'/f"{pdb}_{sf}_{item['id']}_{seed}";folder.mkdir(parents=True,exist_ok=True)
    record={'id':item['id'],'kind':item['kind'],'pdb':pdb,'function':sf,'seed':seed}
    path=folder/'result.json'
    if path.exists():return json.loads(path.read_text())
    start=time.monotonic()
    try:
        prep=json.loads((out/'ligands'/item['id']/'preparation.json').read_text())
        if prep['status']!='complete':raise RuntimeError('Ligand preparation failed: '+prep.get('error','unknown'))
        template=out/'ligands'/item['id']/'ligand.pdbqt';assert sha(template)==prep['ligand_sha256']
        box=json.loads((out/'receptors'/pdb/'box.json').read_text())
        v=Vina(sf_name=sf,cpu=1,seed=seed,verbosity=0)
        v.set_receptor(str(out/'receptors'/pdb/'receptor.pdbqt'));v.set_ligand_from_file(str(template))
        v.compute_vina_maps(center=box['center'],box_size=box['size'])
        v.dock(exhaustiveness=8,n_poses=9)
        energies=v.energies(n_poses=9)
        v.write_poses(str(folder/'poses.pdbqt'),n_poses=9,overwrite=True)
        record.update({'status':'complete','score':float(energies[0,0]),'energies':energies.tolist(),'ligand_sha256':prep['ligand_sha256']})
        if item['kind']=='native_control':
            mol,xyz=coordinates(template.read_text(),(folder/'poses.pdbqt').read_text())
            ref=Chem.SDMolSupplier(str(out/'receptors'/pdb/'native_ref.sdf'),removeHs=True)[0];assert ref is not None
            value=float(rmsd.symmrmsd(np.asarray(ref.GetConformer().GetPositions()),xyz,
                np.array([a.GetAtomicNum() for a in ref.GetAtoms()]),np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),
                Chem.GetAdjacencyMatrix(ref),Chem.GetAdjacencyMatrix(mol),minimize=False))
            record.update({'native_rmsd_A':value,'native_pose_pass':value<2})
    except Exception as exc:
        # A scored native with an RMSD audit error is retained as an audit failure.
        record.update({'status':'failed','error':str(exc)});(folder/'error.log').write_text(traceback.format_exc())
    record['elapsed_seconds']=round(time.monotonic()-start,3);path.write_text(json.dumps(record,indent=2))
    return record


def main(args):
    source,locked,out=args.source,args.locked,args.output
    protocol=json.loads((locked/'protocol.json').read_text())
    assert vina.__version__==protocol['vina_package_version']
    assert sha(locked/'structures.json')==protocol['structures_sha256']
    assert sha(locked/'label_key.json')==protocol['label_key_sha256']
    out.mkdir(parents=True,exist_ok=True)
    if (out/'summary.json').exists():raise RuntimeError('Completed run immutable')
    for name in ['protocol.json','structures.json','label_key.json']:
        dest=out/name
        if dest.exists():assert sha(dest)==sha(locked/name)
        else:shutil.copyfile(locked/name,dest)
    items=json.loads((out/'structures.json').read_text());native=[];input_hashes={}
    for pdb in protocol['targets']:
        folder=out/'receptors'/pdb;folder.mkdir(parents=True,exist_ok=True)
        for name in ['receptor.pdbqt','box.json','native.pdb']:
            src=source/pdb/name;dest=folder/name
            if dest.exists():assert sha(src)==sha(dest)
            else:shutil.copyfile(src,dest)
            input_hashes[str(src)]=sha(src)
        ref=source/'jobs'/f'{pdb}_native_{pdb}_20261018'/'native_ref.sdf'
        shutil.copyfile(ref,folder/'native_ref.sdf');input_hashes[str(ref)]=sha(ref)
        native.append({'id':'native_'+pdb,'kind':'native_control','smiles':(source/pdb/'native_smiles.txt').read_text().strip(),'pdb':pdb})
    preparations=[]
    for item in items+native:
        dest=out/'ligands'/item['id']/'preparation.json'
        preparations.append(json.loads(dest.read_text()) if dest.exists() else prepare(item,out))
    (out/'preparation_summary.json').write_text(json.dumps(preparations,indent=2))
    (out/'input_hashes.json').write_text(json.dumps(input_hashes,indent=2))
    tasks_native=[(out,item,item['pdb'],sf,seed) for item in native for sf in protocol['functions'] for seed in protocol['seeds']]
    tasks_study=[(out,item,pdb,sf,seed) for pdb in protocol['targets'] for sf in protocol['functions'] for item in items for seed in protocol['seeds']]
    assert len(tasks_native)==18 and len(tasks_study)==234
    results=[]
    for phase,tasks in [('native',tasks_native),('study',tasks_study)]:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures=[pool.submit(dock,task) for task in tasks]
            for future in concurrent.futures.as_completed(futures):
                result=future.result();results.append(result)
                state={'phase':phase,'finished':len(results),'scheduled':252,'complete':sum(r['status']=='complete' for r in results),'failed':sum(r['status']!='complete' for r in results)}
                (out/'run_status.json').write_text(json.dumps(state,indent=2))
                print(json.dumps({**state,'last':{k:result.get(k) for k in ['id','pdb','function','seed','score','native_rmsd_A','error']}}),flush=True)
    assert len(results)==252
    assert all(sha(Path(path))==value for path,value in input_hashes.items())
    (out/'summary.json').write_text(json.dumps({'protocol':protocol,'vina_version':vina.__version__,'results':results,'input_files_unchanged':True},indent=2))
    print('DONE 252',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--locked',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=6);main(p.parse_args())
