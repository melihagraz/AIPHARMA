"""DOCK6 conformer-ensemble control generated from SMILES, never crystal coordinates.
Fixed development protocol: ETKDGv3 seed 20261015, 32 attempts, prune 0.25 A,
MMFF94s conformer optimization, rigid orientation search + simplex minimization.
"""
from pathlib import Path
import json,re,subprocess,time,numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from spyrmsd import rmsd as srmsd
root=Path.home();D=root/'tools/dock6';ob=root/'miniforge3/envs/aipharma/bin/obabel'
all_results={}
for pdb,code in [('1PY5','PY1'),('1VJY','460')]:
 p=root/('.aipharma_dock6_amber_'+pdb)
 desc=json.loads((p/(code+'_ccd.json')).read_text())['pdbx_chem_comp_descriptor']
 smiles=next(d['descriptor'] for d in desc if d['type']=='SMILES_CANONICAL')
 mol=Chem.AddHs(Chem.MolFromSmiles(smiles));settings=AllChem.ETKDGv3();settings.randomSeed=20261015;settings.pruneRmsThresh=0.25;settings.numThreads=2
 ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=32,params=settings))
 if not ids:raise RuntimeError('No conformers embedded')
 opt=AllChem.MMFFOptimizeMoleculeConfs(mol,mmffVariant='MMFF94s',numThreads=2)
 writer=Chem.SDWriter(str(p/'ensemble.sdf'))
 for cid in ids:
  mol.SetProp('_Name',f'{pdb}_conf_{cid}');writer.write(mol,confId=cid)
 writer.close()
 subprocess.run([str(ob),'ensemble.sdf','-O','ensemble_raw.mol2','--partialcharge','gasteiger'],cwd=p,check=True,capture_output=True)
 blocks=(p/'ensemble_raw.mol2').read_text().split('@<TRIPOS>MOLECULE')[1:];normal=[]
 for b in blocks:
  lines=b.strip('\n').splitlines();counts=lines[1].split();counts[2]='1';lines[1]=' '.join(counts);section=None
  for i,line in enumerate(lines):
   if line.startswith('@<TRIPOS>'):section=line
   elif section=='@<TRIPOS>ATOM' and line.strip():
    a=line.split();a[6]='1';lines[i]=' '.join(a)
  normal.append('@<TRIPOS>MOLECULE\n'+'\n'.join(lines).replace('@<TRIPOS>ATOM','\n@<TRIPOS>ATOM',1)+'\n@<TRIPOS>SUBSTRUCTURE\n1 UNL1 1 RESIDUE 1 A UNL1 0 ROOT\n')
 (p/'ensemble.mol2').write_text(''.join(normal))
 params=dict(line.split(None,1) for line in (p/'rigid_control.in').read_text().splitlines() if line.strip())
 params.update({'ligand_atom_file':'ensemble.mol2','ligand_outfile_prefix':'ensemble_docked','num_final_scored_poses':1})
 (p/'ensemble.in').write_text(''.join(f'{k:60s} {v}\n' for k,v in params.items()))
 subprocess.run([str(D/'bin/dock6'),'-i','ensemble.in','-o','ensemble.out'],cwd=p,check=True)
 output=(p/'ensemble_docked_scored.mol2').read_text()
 entries=re.split(r'(?=##########\s+Name:)',output)
 ref=Chem.SDMolSupplier(str(p/'crystal.sdf'),removeHs=True)[0]
 gen=Chem.RemoveHs(mol)
 refxyz=np.array(ref.GetConformer().GetPositions());refnums=np.array([a.GetAtomicNum() for a in ref.GetAtoms()]);refadj=Chem.GetAdjacencyMatrix(ref)
 nums=np.array([a.GetAtomicNum() for a in gen.GetAtoms()]);adj=Chem.GetAdjacencyMatrix(gen)
 rows=[]
 for entry in entries:
  if '@<TRIPOS>ATOM' not in entry:continue
  atoms=[r.split() for r in entry.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].splitlines() if r.strip() and not r.split()[5].startswith('H')]
  xyz=np.array([r[2:5] for r in atoms],float)
  score=float(re.search(r'Grid_Score:\s*([-\d.]+)',entry).group(1))
  rmsd=float(srmsd.symmrmsd(refxyz,xyz,refnums,nums,refadj,adj,minimize=False))
  rows.append({'name':re.search(r'Name:\s*(\S+)',entry).group(1),'score':score,'rmsd_A':rmsd})
 rows.sort(key=lambda r:r['score'])
 if not rows:raise RuntimeError('No docked conformers')
 result={'pdb':pdb,'protocol':'SMILES ETKDGv3/MMFF94s ensemble + DOCK6 rigid orientation/minimization','seed':20261015,'attempted_conformers':32,'embedded_conformers':len(ids),'docked_conformers':len(rows),'mmff_optimization_results':[list(x) for x in opt],'top_scored_pose':rows[0],'top_scored_pose_recovery_pass':rows[0]['rmsd_A']<2,'best_sampled_rmsd_A':min(r['rmsd_A'] for r in rows),'rows':rows,'limitations':['Protocol development on these two ALK5 structures','Ligand-derived site matching centers','Not biological or prospective external validation']}
 (p/'ensemble_summary.json').write_text(json.dumps(result,indent=2));all_results[pdb]=result;print(json.dumps({pdb:result}),flush=True)
