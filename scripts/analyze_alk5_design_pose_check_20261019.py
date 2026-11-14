"""Auditable distance and retained-core pose diagnostics, without affinity ranking."""
import importlib.util,json,math,statistics,sys
from pathlib import Path
import numpy as np
from rdkit import Chem,RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Scaffolds import MurckoScaffold
p=Path(sys.argv[1]); data=json.loads((p/'summary.json').read_text()); rows=data['results']
assert len(rows)==126
spec=importlib.util.spec_from_file_location('pose_helper',Path('/home/cobre/run_alk5_five_block_validation_20261018.py'))
helper=importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
factory=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
def atom(line): return dict(name=line[12:16].strip(),resname=line[17:20].strip(),chain=line[21],resnum=int(line[22:26]),xyz=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
protein={}; hinge={}; mapping=[]
for pdb in ['1PY5','1VJY','2WOT']:
    source=Path('/home/cobre/aipharma_ensemble_20261018')/pdb/'protein.pdb'
    original=[atom(l) for l in source.read_text().splitlines() if l.startswith('ATOM')]
    prepared=[atom(l) for l in (p/'receptors'/pdb/'receptor.pdbqt').read_text().splitlines() if l.startswith('ATOM') and not l.split()[-1].startswith('H')]
    protein[pdb]=np.array([r['xyz'] for r in prepared]); hinge[pdb]={}
    for name in ['N','O']:
        matches=[r for r in original if r['resnum']==283 and r['resname']=='HIS' and r['name']==name]
        assert len(matches)==1,(pdb,name,len(matches))
        orig=matches[0]; nearest=min(prepared,key=lambda r:float(np.linalg.norm(r['xyz']-orig['xyz'])))
        distance=float(np.linalg.norm(nearest['xyz']-orig['xyz']))
        assert distance<.02 and nearest['name']==name and nearest['resname']=='HIS'
        hinge[pdb][name]=nearest['xyz']; mapping.append(dict(pdb=pdb,original_chain=orig['chain'],original_residue=283,prepared_residue=nearest['resnum'],atom=name,coordinate_difference_A=distance))
native={}
for pdb in protein:
    for sf in ['vina','vinardo']:
        ns=[r for r in rows if r['kind']=='native_control' and r['pdb']==pdb and r['function']==sf]
        assert len(ns)==3
        native[pdb,sf]=sum(r.get('native_pose_pass',False) for r in ns)>=2
cache={}; diagnostics=[]
for r in rows:
    if r['status']!='complete':
        diagnostics.append(dict(**r,geometry_status='calculation_failed')); continue
    template=(p/'ligands'/r['id']/'ligand.pdbqt').read_text()
    pose=(p/'jobs'/f"{r['pdb']}_{r['function']}_{r['id']}_{r['seed']}"/'poses.pdbqt').read_text()
    m,xyz=helper.coordinates(template,pose); cache[r['id'],r['pdb'],r['function'],r['seed']]=(m,xyz)
    contacts=[]
    for f in factory.GetFeaturesForMol(m):
        target='N' if f.GetFamily()=='Acceptor' else ('O' if f.GetFamily()=='Donor' else None)
        if target:
            for i in f.GetAtomIds():
                d=float(np.linalg.norm(xyz[i]-hinge[r['pdb']][target]))
                if 2.4<=d<=3.5: contacts.append(dict(ligand_atom=i,feature=f.GetFamily(),hinge_atom=target,distance_A=d))
    min_distance=float(np.sqrt(((xyz[:,None,:]-protein[r['pdb']][None,:,:])**2).sum(axis=2)).min())
    diagnostics.append(dict(id=r['id'],kind=r['kind'],pdb=r['pdb'],function=r['function'],seed=r['seed'],geometry_status='calculated',native_group_interpretable=native[r['pdb'],r['function']],contacts=contacts,min_protein_heavy_atom_distance_A=min_distance,severe_clash=min_distance<1.6,distance_compatible_pose=bool(contacts) and min_distance>=1.6))
comparisons=[]
for cid,pid in [('design_11','parent_11'),('design_20','parent_20'),('refine_hydroxymethyl','parent_20'),('refine_carbamate_methyl','parent_20')]:
    for pdb in protein:
        for sf in ['vina','vinardo']:
            rmsds=[]
            for seed in [20261018,20261019,20261020]:
                if (cid,pdb,sf,seed) not in cache or (pid,pdb,sf,seed) not in cache: continue
                pm,px=cache[pid,pdb,sf,seed]; cm,cx=cache[cid,pdb,sf,seed]
                core=MurckoScaffold.GetScaffoldForMol(pm); pa=pm.GetSubstructMatches(core,uniquify=False); ca=cm.GetSubstructMatches(core,uniquify=False)
                assert pa and ca
                rmsds.append(dict(seed=seed,core_rmsd_A=min(float(np.sqrt(((px[list(i)]-cx[list(j)])**2).sum(axis=1).mean())) for i in pa for j in ca)))
            comparisons.append(dict(candidate=cid,parent=pid,pdb=pdb,function=sf,paired_seed_core_rmsds=rmsds,median_core_rmsd_A=statistics.median(r['core_rmsd_A'] for r in rmsds) if rmsds else None,native_group_interpretable=native[pdb,sf]))
summary=[]
for cid in ['parent_11','parent_20','design_11','design_20','refine_hydroxymethyl','refine_carbamate_methyl']:
    groups=[]
    for pdb in protein:
        for sf in ['vina','vinardo']:
            ds=[r for r in diagnostics if r['id']==cid and r['pdb']==pdb and r['function']==sf]
            groups.append(dict(pdb=pdb,function=sf,interpretable=native[pdb,sf],compatible_seeds=sum(r.get('distance_compatible_pose',False) for r in ds),scheduled_seeds=3))
    pid={'design_11':'parent_11','design_20':'parent_20','refine_hydroxymethyl':'parent_20','refine_carbamate_methyl':'parent_20'}.get(cid)
    for g in groups:
        parent_rows=[r for r in diagnostics if r['id']==pid and r['pdb']==g['pdb'] and r['function']==g['function']]
        g['paired_parent_compatible_seeds']=sum(r.get('distance_compatible_pose',False) for r in parent_rows) if pid else None
        g['paired_parent_supports_interpretation']=bool(g['interpretable'] and g['paired_parent_compatible_seeds'] is not None and g['paired_parent_compatible_seeds']>=2)
    summary.append(dict(id=cid,groups=groups,interpretable_groups=sum(g['interpretable'] for g in groups),compatible_interpretable_groups=sum(g['interpretable'] and g['compatible_seeds']>=2 for g in groups),paired_parent_supported_groups=sum(g['paired_parent_supports_interpretation'] for g in groups),compatible_paired_parent_supported_groups=sum(g['paired_parent_supports_interpretation'] and g['compatible_seeds']>=2 for g in groups),activity_validated=False))
for name,obj in [('hinge_atom_mapping.json',mapping),('pose_diagnostics.json',diagnostics),('core_pose_comparison.json',comparisons),('geometry_summary.json',summary)]: (p/name).write_text(json.dumps(obj,indent=2,allow_nan=False))
print(json.dumps(summary,indent=2),flush=True)
