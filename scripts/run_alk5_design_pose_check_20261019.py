"""Paired parent/design pose diagnostic; docking scores never rank these designs."""
import concurrent.futures,hashlib,importlib.util,json,shutil,sys,time
from pathlib import Path

root=Path(sys.argv[1]); source=Path(sys.argv[2]); out=root/'pose_check'
out.mkdir(exist_ok=True)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
helper=Path('/home/cobre/run_alk5_five_block_validation_20261018.py')
spec=importlib.util.spec_from_file_location('pose_helper',helper); mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
protocol=dict(purpose='Exploratory pose geometry only; no activity ranking',
    structures=['parent_11','parent_20','design_11','design_20','refine_hydroxymethyl','refine_carbamate_methyl'],
    receptors=['1PY5','1VJY','2WOT'],functions=['vina','vinardo'],seeds=[20261018,20261019,20261020],
    preparation='Unchanged archived helper ETKDG16/MMFF94s/Meeko; no crystal pose initialization',
    controls='Repeat all 18 native control runs; require >=2/3 RMSD<2A within receptor/scoring group for interpretable geometry',
    geometry='Backbone His283 N-donor to ligand acceptor or O-acceptor to ligand donor heavy-atom distance 2.4-3.5A; severe protein heavy atom clash <1.6A',
    interpretation='Distance-compatible contact is not a confirmed H-bond; no angle, protonation ensemble, water network or free-energy inference',
    reporting='All seeds/receptors/functions, including native failures; no best-score choice',
    acceptance='No biological activity or safety acceptance; pose compatibility only',
    helper_sha256=sha(helper),script_sha256=sha(Path(__file__)))
pp=out/'protocol.json'
if pp.exists(): assert json.loads(pp.read_text())==protocol
else: pp.write_text(json.dumps(protocol,indent=2))
initial=json.loads((root/'input_designs.json').read_text())['candidates']; lookup={r['candidate_id']:r for r in initial}
items=[]
for suffix,cid in [('11','ALK5_design_011'),('20','ALK5_design_020')]:
    c=lookup[cid]
    items.extend([dict(id='parent_'+suffix,kind='measured_parent',smiles=c['provenance']['parent_smiles']),dict(id='design_'+suffix,kind='unmeasured_design',smiles=c['smiles'])])
for stage,cid,name in [('refinement','ALK5_refine_001','refine_hydroxymethyl'),('refinement_sar_guided','ALK5_refine_003','refine_carbamate_methyl')]:
    c=next(r for r in json.loads((root/stage/'results.json').read_text()) if r['candidate_id']==cid)
    items.append(dict(id=name,kind='unmeasured_design',smiles=c['smiles']))
native=[]; hashes={}
for pdb in protocol['receptors']:
    dest=out/'receptors'/pdb; dest.mkdir(parents=True,exist_ok=True)
    for name in ['receptor.pdbqt','box.json','native.pdb','native_ref.sdf']:
        orig=source/'receptors'/pdb/name; target=dest/name
        if target.exists(): assert sha(target)==sha(orig)
        else: shutil.copyfile(orig,target)
        hashes[str(orig)]=sha(orig)
    smi=(Path('/home/cobre/aipharma_ensemble_20261018')/pdb/'native_smiles.txt').read_text().strip()
    native.append(dict(id='native_'+pdb,kind='native_control',smiles=smi,pdb=pdb))
(out/'structures.json').write_text(json.dumps(items+native,indent=2)); (out/'input_hashes.json').write_text(json.dumps(hashes,indent=2))
for item in items+native:
    if not (out/'ligands'/item['id']/'preparation.json').exists(): mod.prepare(item,out)
tasks=[(out,r,r['pdb'],sf,seed) for r in native for sf in protocol['functions'] for seed in protocol['seeds']]
tasks += [(out,r,pdb,sf,seed) for r in items for pdb in protocol['receptors'] for sf in protocol['functions'] for seed in protocol['seeds']]
assert len(tasks)==126
def execute(task): return mod.dock(task)
results=[]
with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
    for f in concurrent.futures.as_completed([pool.submit(execute,t) for t in tasks]):
        r=f.result(); results.append(r)
        status=dict(finished=len(results),scheduled=len(tasks),complete=sum(x['status']=='complete' for x in results),failed=sum(x['status']!='complete' for x in results))
        (out/'run_status.json').write_text(json.dumps(status,indent=2))
        print(json.dumps(dict(**status,last=r['id'])),flush=True)
assert all(sha(Path(k))==v for k,v in hashes.items())
(out/'summary.json').write_text(json.dumps(dict(results=results,input_files_unchanged=True,activity_validated=False),indent=2))
