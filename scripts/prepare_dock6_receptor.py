"""Prepare rigid receptor with AMBER14 charges and optimized hydrogens.
Existing heavy atom coordinates are fixed; missing residues are not modeled.
Run in WSL aipharma Python, with OpenMM/PDBFixer in ~/tools/aipharma-prep-python.
"""
from pathlib import Path
import sys,json,argparse,subprocess,random
import numpy as np
sys.path.insert(0,str(Path.home()/'tools/aipharma-prep-python'))
from pdbfixer import PDBFixer
from openmm import app,unit,NonbondedForce,VerletIntegrator,Context,LocalEnergyMinimizer,Platform
parser=argparse.ArgumentParser();parser.add_argument('--pdb',default='1PY5');args=parser.parse_args()
random.seed(20261015);np.random.seed(20261015)
root=Path.home(); out=root/('.aipharma_amberprep_'+args.pdb);out.mkdir(exist_ok=True)
source=root/'.aipharma_redock'/args.pdb/'receptor.pdb'
fixer=PDBFixer(filename=str(source))
fixer.findMissingResidues(); skipped={str(k):v for k,v in fixer.missingResidues.items()};fixer.missingResidues={}
fixer.findMissingAtoms();missing={r.name+r.id:[a.name for a in atoms] for r,atoms in fixer.missingAtoms.items()}
fixer.addMissingAtoms(seed=20261015)
ff=app.ForceField('amber14-all.xml')
model=app.Modeller(fixer.topology,fixer.positions)
variants=model.addHydrogens(ff,pH=7.4,platform=Platform.getPlatformByName('CPU'))
system=ff.createSystem(model.topology,nonbondedMethod=app.NoCutoff,constraints=None)
initial=np.array(model.positions.value_in_unit(unit.angstrom))
for atom in model.topology.atoms():
 if atom.element.symbol!='H': system.setParticleMass(atom.index,0)
integrator=VerletIntegrator(0.001)
context=Context(system,integrator,Platform.getPlatformByName('CPU'),{'Threads':'2'})
context.setPositions(model.positions)
LocalEnergyMinimizer.minimize(context,10*unit.kilojoule_per_mole/unit.nanometer,500)
positions=context.getState(getPositions=True).getPositions()
coords=np.array(positions.value_in_unit(unit.angstrom))
heavy=[a.index for a in model.topology.atoms() if a.element.symbol!='H']
assert np.max(np.abs(coords[heavy]-initial[heavy]))<1e-6
with (out/'receptor.pdb').open('w') as f: app.PDBFile.writeFile(model.topology,positions,f,keepIds=True)
nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
charges=[nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles())]
subprocess.run([str(root/'miniforge3/envs/aipharma/bin/obabel'),str(out/'receptor.pdb'),'-O',str(out/'receptor_raw.mol2')],check=True,capture_output=True)
lines=(out/'receptor_raw.mol2').read_text().splitlines();section=None;residues={};used=set()
for i,line in enumerate(lines):
 if line.startswith('@<TRIPOS>'):section=line
 elif section=='@<TRIPOS>ATOM' and line.strip():
  a=line.split();delta=np.linalg.norm(coords-np.array(a[2:5],float),axis=1);idx=int(np.argmin(delta));assert delta[idx]<0.002 and idx not in used;used.add(idx)
  key=(a[6],a[7]);residues.setdefault(key,(len(residues)+1,a[0]));a[6]=str(residues[key][0]);a[8]=f'{charges[idx]:.6f}';lines[i]=' '.join(a)
assert len(used)==len(charges)
counts=lines[2].split();counts[2]=str(len(residues));lines[2]=' '.join(counts);lines[4]='USER_CHARGES'
lines.extend(['@<TRIPOS>SUBSTRUCTURE']+[f'{idx} {key[1]} {first} RESIDUE 1 A {key[1]} 0 ROOT' for key,(idx,first) in residues.items()])
(out/'receptor.mol2').write_text('\n'.join(lines).replace('@<TRIPOS>ATOM','\n@<TRIPOS>ATOM',1)+'\n')
meta={'pdb':args.pdb,'charge_method':'OpenMM amber14-all.xml (protein ff14SB)','pH':7.4,'hydrogens_minimized':True,'heavy_atom_max_shift_A':float(np.max(np.abs(coords[heavy]-initial[heavy]))),'net_charge':sum(charges),'atoms':len(charges),'missing_atoms_added':missing,'missing_residues_not_modeled':skipped,'hydrogen_variants':variants,'seed':20261015}
(out/'preparation.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))
