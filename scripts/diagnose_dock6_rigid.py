from pathlib import Path
import subprocess,json,re,numpy as np
from rdkit import Chem
from spyrmsd import rmsd as srmsd
root=Path.home();D=root/'tools/dock6';results={}
for pdb in ['1PY5','1VJY']:
 p=root/('.aipharma_dock6_amber_'+pdb)
 params=dict(line.split(None,1) for line in (p/'crystal.in').read_text().splitlines() if line.strip())
 params.update({'conformer_search_type':'rigid','ligand_outfile_prefix':'rigid_control'})
 (p/'rigid_control.in').write_text(''.join(f'{k:60s} {v}\n' for k,v in params.items()))
 subprocess.run([str(D/'bin/dock6'),'-i','rigid_control.in','-o','rigid_control.out'],cwd=p,check=True)
 def atoms(file):return [r.split() for r in (p/file).read_text().split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].splitlines() if r.strip() and not r.split()[5].startswith('H')]
 a=atoms('crystal.mol2');b=atoms('rigid_control_scored.mol2');x=np.array([r[2:5] for r in a],float);y=np.array([r[2:5] for r in b],float)
 mol=Chem.SDMolSupplier(str(p/'crystal.sdf'),removeHs=True)[0];nums=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]);adj=Chem.GetAdjacencyMatrix(mol)
 result={'rmsd_A':float(srmsd.symmrmsd(x,y,nums,nums,adj,adj,minimize=False)),'grid_score':float(re.search(r'Grid_Score:\s*([-\d.]+)',(p/'rigid_control_scored.mol2').read_text()).group(1))}
 (p/'rigid_control_summary.json').write_text(json.dumps(result,indent=2));results[pdb]=result
print(json.dumps(results))
