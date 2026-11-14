"""Validate reported RMSD using connectivity parsed independently from output MOL2."""
from pathlib import Path
from rdkit import Chem
import numpy as np,json,re
from spyrmsd import rmsd as srmsd
out={}
def parse(block):
 atoms=[r.split() for r in block.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].splitlines() if r.strip()]
 bonds=[r.split() for r in block.split('@<TRIPOS>BOND\n')[1].split('@<TRIPOS>')[0].splitlines() if r.strip()]
 heavy=[r for r in atoms if not r[5].startswith('H')];mapping={int(r[0]):i for i,r in enumerate(heavy)};adj=np.zeros((len(heavy),len(heavy)),dtype=int)
 for b in bonds:
  i,j=map(int,b[1:3])
  if i in mapping and j in mapping:adj[mapping[i],mapping[j]]=adj[mapping[j],mapping[i]]=1
 xyz=np.array([r[2:5] for r in heavy],float);nums=np.array([Chem.GetPeriodicTable().GetAtomicNumber(r[5].split('.')[0]) for r in heavy])
 return xyz,nums,adj,bonds
for pdb in ['1PY5','1VJY']:
 p=Path.home()/('.aipharma_dock6_amber_'+pdb);summary=json.loads((p/'ensemble_summary.json').read_text())
 ref=Chem.SDMolSupplier(str(p/'crystal.sdf'),removeHs=True)[0];nums=np.array([a.GetAtomicNum() for a in ref.GetAtoms()]);adj=Chem.GetAdjacencyMatrix(ref);xyz=np.array(ref.GetConformer().GetPositions())
 entries=re.split(r'(?=##########\s+Name:)',(p/'ensemble_docked_scored.mol2').read_text())
 top=next(e for e in entries if re.search(r'Name:\s*'+re.escape(summary['top_scored_pose']['name'])+r'\s',e))
 ox,on,oa,bonds=parse(top)
 actual=float(srmsd.symmrmsd(xyz,ox,nums,on,adj,oa,minimize=False));assert abs(actual-summary['top_scored_pose']['rmsd_A'])<1e-6
 prepared=(p/'ensemble.mol2').read_text();prepblock=next(b for b in prepared.split('@<TRIPOS>MOLECULE') if summary['top_scored_pose']['name'] in b)
 _,_,_,pb=parse(prepblock);assert bonds==pb
 perm=np.arange(len(nums))[::-1];selfrms=float(srmsd.symmrmsd(xyz,xyz[perm],nums,nums[perm],adj,adj[perm][:,perm],minimize=False));assert selfrms<1e-6
 assert summary['docked_conformers']==len(summary['rows']);assert summary['top_scored_pose']==min(summary['rows'],key=lambda r:r['score'])
 out[pdb]={'mol2_connectivity_and_bond_orders_preserved':True,'independently_recomputed_top_rmsd_A':actual,'symmetry_mapping_self_test_rmsd_A':selfrms,'pose_count_and_top_score_verified':True}
print(json.dumps(out));(Path.home()/'.aipharma_dock6_alk5/quality_checks.json').write_text(json.dumps(out,indent=2))
