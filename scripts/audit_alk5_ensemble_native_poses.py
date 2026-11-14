from pathlib import Path
import json,re,subprocess
import numpy as np
from rdkit import Chem
from spyrmsd import rmsd
p=Path.home()/'aipharma_ensemble_20261018';env=Path.home()/'miniforge3/envs/aipharma/bin';results=[]
for pdb in ['1PY5','1VJY','2WOT']:
 w=p/'jobs'/f'{pdb}_native_{pdb}_20261018';old=json.loads((w/'result.json').read_text());(w/'before_rmsd_audit.json').write_text(json.dumps(old,indent=2)) if not (w/'before_rmsd_audit.json').exists() else None;text=(w/'pose.pdbqt').read_text().split('ENDMDL',1)[0];template=(p/'ligands'/f'native_{pdb}'/'ligand.pdbqt').read_text();smiles=next(l[len('REMARK SMILES '):] for l in template.splitlines() if l.startswith('REMARK SMILES ') and not l.startswith('REMARK SMILES IDX'))
 mapping=[]
 for l in template.splitlines():
  if l.startswith('REMARK SMILES IDX '):mapping.extend(map(int,l[len('REMARK SMILES IDX '):].split()))
 positions={int(l[6:11]):[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in text.splitlines() if l.startswith(('ATOM','HETATM'))};mol=Chem.MolFromSmiles(smiles);coords=np.full((mol.GetNumAtoms(),3),np.nan)
 for i,j in zip(mapping[::2],mapping[1::2]):coords[i-1]=positions[j]
 assert np.isfinite(coords).all()
 native=w/'native_ref.sdf';subprocess.run([str(env/'obabel'),str(p/pdb/'native.pdb'),'-O',str(native)],capture_output=True,check=True);ref=Chem.SDMolSupplier(str(native),removeHs=True)[0];assert ref is not None
 value=float(rmsd.symmrmsd(np.array(ref.GetConformer().GetPositions()),coords,np.array([a.GetAtomicNum() for a in ref.GetAtoms()]),np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),Chem.GetAdjacencyMatrix(ref),Chem.GetAdjacencyMatrix(mol),minimize=False));out={**old,'status':'complete','error':None,'native_rmsd_A':value,'pose_control_pass':value<2,'rmsd_audit':'Corrected nonexistent CCD field access and used Meeko SMILES IDX mapping; existing docking not rerun; original result preserved'};(w/'result.json').write_text(json.dumps(out,indent=2));results.append(out)
(p/'native_pose_audit.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
