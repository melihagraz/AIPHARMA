"""Reproducible preliminary DOCK6 ALK5 run, executed with WSL Python.
Defaults to Gasteiger charges; accepts a separately prepared AMBER receptor.
Ligand-derived matching spheres define the known crystallographic pocket.
"""
from pathlib import Path
import subprocess, json, re, hashlib, time, argparse, shutil
parser=argparse.ArgumentParser()
parser.add_argument('--pdb',default='1PY5')
parser.add_argument('--receptor-mol2',type=Path)
parser.add_argument('--output',type=Path)
parser.add_argument('--grid-spacing',type=float,default=0.5)
args=parser.parse_args()
ROOT=Path.home()
DOCK=ROOT/'tools/dock6'
ENV=ROOT/'miniforge3/envs/aipharma/bin'
OUT=args.output or ROOT/'.aipharma_dock6_alk5'
OUT.mkdir(exist_ok=True)

def run(args, log, stdin=None):
    with (OUT/log).open('w') as f:
        p=subprocess.run([str(a) for a in args],cwd=OUT,input=stdin,text=True,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
    if p.returncode: raise RuntimeError(f'{log}: exit {p.returncode}')

def config(template, changes):
    params=dict(line.split(None,1) for line in template.read_text().splitlines() if line.strip())
    params.update({k:str(v) for k,v in changes.items()})
    return ''.join(f'{k:60s} {v}\n' for k,v in params.items())

def mol2_atoms(path):
    text=path.read_text().split('@<TRIPOS>ATOM\n',1)[1].split('@<TRIPOS>',1)[0]
    return [line.split() for line in text.splitlines() if line.strip()]

started=time.time()
source=ROOT/'.aipharma_redock'/args.pdb
ligand_code=next(line[17:20].strip() for line in (source/'crystal_lig.pdb').read_text().splitlines() if line.startswith('HETATM'))
if args.receptor_mol2:
    shutil.copyfile(args.receptor_mol2,OUT/'receptor.mol2')
else:
    run([ENV/'obabel',source/'receptor.pdb','-O','receptor.mol2','-h','--partialcharge','gasteiger'],'prepare_receptor.log')
from rdkit import Chem
from rdkit.Chem import AllChem
import urllib.request
ccd=OUT/(ligand_code+'_ccd.json')
if not ccd.exists():
    with urllib.request.urlopen('https://data.rcsb.org/rest/v1/core/chemcomp/'+ligand_code,timeout=30) as response:
        ccd.write_bytes(response.read())
descriptors=json.loads(ccd.read_text())['pdbx_chem_comp_descriptor']
smiles=next(d['descriptor'] for d in descriptors if d['type']=='SMILES_CANONICAL')
reference=Chem.MolFromPDBFile(str(source/'crystal_lig.pdb'),removeHs=True)
reference=AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles),reference)
reference=Chem.AddHs(reference,addCoords=True)
writer=Chem.SDWriter(str(OUT/'crystal.sdf')); writer.write(reference); writer.close()
run([ENV/'obabel','crystal.sdf','-O','crystal.mol2','--partialcharge','gasteiger'],'prepare_crystal.log')
for name in (('crystal.mol2',) if args.receptor_mol2 else ('receptor.mol2','crystal.mol2')):
    path=OUT/name
    content=path.read_text().replace('@<TRIPOS>ATOM', '\n@<TRIPOS>ATOM', 1)
    lines=content.splitlines()
    section=None
    residues={}
    for i,line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):
            section=line
        elif section=='@<TRIPOS>ATOM' and line.strip():
            fields=line.split()
            key=(fields[6],fields[7])
            if key not in residues: residues[key]=(len(residues)+1, fields[0])
            fields[6]=str(residues[key][0])
            lines[i]=' '.join(fields)
    counts=lines[2].split(); counts[2]=str(len(residues)); lines[2]=' '.join(counts)
    lines.extend(['@<TRIPOS>SUBSTRUCTURE']+[f'{idx} {key[1]} {root} RESIDUE 1 A {key[1]} 0 ROOT' for key,(idx,root) in residues.items()])
    path.write_text('\n'.join(lines)+'\n')
atoms=mol2_atoms(OUT/'crystal.mol2')
heavy=[a for a in atoms if not a[5].startswith('H')]
with (OUT/'site.sph').open('w') as f:
    f.write(f'Ligand-derived matching centers; {args.pdb} known pocket\ncluster     1   number of spheres in cluster {len(heavy):5d}\n')
    for i,a in enumerate(heavy,1):
        x,y,z=map(float,a[2:5]); f.write(f'{i:5d}{x:10.5f}{y:10.5f}{z:10.5f}{1.4:8.3f}{i:5d} 0  0\n')
run([DOCK/'bin/showbox'],'showbox.log','Y\n6\nsite.sph\n1\nbox.pdb\n')
(OUT/'grid.in').write_text(config(DOCK/'install/test/grid_generation/grid.in',{'receptor_file':'receptor.mol2','box_file':'box.pdb','vdw_definition_file':DOCK/'parameters/vdw_AMBER_parm99.defn','score_grid_prefix':'grid','contact_score':'no','grid_spacing':args.grid_spacing}))
run([DOCK/'bin/grid','-i','grid.in','-o','grid.out'],'grid.log')
base={'ligand_atom_file':'crystal.mol2','receptor_site_file':'site.sph','grid_score_grid_prefix':'grid','vdw_defn_file':DOCK/'parameters/vdw_AMBER_parm99.defn','flex_defn_file':DOCK/'parameters/flex.defn','flex_drive_file':DOCK/'parameters/flex_drive.tbl','ligand_outfile_prefix':'crystal_redock','min_anchor_size':5,'pruning_max_orients':100,'pruning_clustering_cutoff':2,'max_orientations':1000,'num_final_scored_poses':10,'score_threshold':100000,'rank_ligands':'no'}
base.update({'minimize_ligand':'yes','simplex_max_iterations':500,'simplex_tors_premin_iterations':0,'simplex_max_cycles':1,'simplex_score_converge':0.1,'simplex_cycle_converge':1.0,'simplex_trans_step':1.0,'simplex_rot_step':0.1,'simplex_tors_step':10.0,'simplex_final_min':'no','simplex_random_seed':0,'simplex_restraint_min':'no'})
(OUT/'crystal.in').write_text(config(DOCK/'install/test/grid_score/orientflexiblemulti.dockin',base))
run([DOCK/'bin/dock6','-i','crystal.in','-o','crystal.out'],'crystal.log')
report={'stage':'preliminary_crystal_redocking','target':'ALK5','pdb':args.pdb,'charge_method':('OpenMM AMBER14 receptor; Gasteiger ligand' if args.receptor_mol2 else 'Open Babel Gasteiger (receptor and ligand)'),'site_method':'crystal heavy-atom matching centers; 6 A grid padding','elapsed_seconds':round(time.time()-started,2),'limitations':['Protocol under development; not external validation','Known-ligand-derived matching spheres','No prospective validation claim'],'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.in')}}
scored=OUT/'crystal_redock_scored.mol2'
if scored.exists(): report['score_headers']=re.findall(r'^##########.*$',scored.read_text(),re.M)
if scored.exists():
    import numpy as np
    from spyrmsd import rmsd as srmsd
    pose=mol2_atoms(scored)
    pose_heavy=[a for a in pose if not a[5].startswith('H')]
    original=Chem.RemoveHs(reference)
    # DOCK preserves MOL2 atom identifiers through search.
    xyz={int(a[0]):[float(v) for v in a[2:5]] for a in pose_heavy}
    refheavy=[a for a in mol2_atoms(OUT/'crystal.mol2') if not a[5].startswith('H')]
    refcoords=np.array([[float(v) for v in a[2:5]] for a in refheavy])
    posecoords=np.array([xyz[int(a[0])] for a in refheavy])
    numbers=np.array([a.GetAtomicNum() for a in original.GetAtoms()])
    adjacency=Chem.GetAdjacencyMatrix(original)
    report['symmetry_corrected_rmsd_A']=float(srmsd.symmrmsd(refcoords,posecoords,numbers,numbers,adjacency,adjacency,minimize=False))
    report['pose_recovery_pass']=report['symmetry_corrected_rmsd_A']<2.0
native=base.copy()
native.update({'conformer_search_type':'rigid','orient_ligand':'no','minimize_ligand':'no','ligand_outfile_prefix':'native_score'})
(OUT/'native.in').write_text(config(DOCK/'install/test/grid_score/orientflexiblemulti.dockin',native))
run([DOCK/'bin/dock6','-i','native.in','-o','native.out'],'native.log')
native_text=(OUT/'native_score_scored.mol2').read_text()
report['native_crystal_grid_score']=float(re.search(r'Grid_Score:\s*([-\d.]+)',native_text).group(1))
report['redocked_grid_score']=float(re.search(r'Grid_Score:\s*([-\d.]+)',scored.read_text()).group(1))
report['grid_spacing_A']=args.grid_spacing
report['candidate_ranking_status']='withheld_pending_protocol_pose_recovery'
report['bond_order_source']='RCSB CCD '+ligand_code+' canonical SMILES assigned with RDKit'
report['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.in')}
(OUT/'summary.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2),flush=True)
