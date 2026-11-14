"""Compare exact and grid native-pose scores without changing scoring physics."""
import runpy,json,re
ns=runpy.run_path(str(__import__('pathlib').Path.home()/'validate_dock6_alk5.py'))
p=ns['OUT']; D=ns['DOCK']; run=ns['run']; config=ns['config']
base=ns['native'].copy()
results={}
cont=base.copy(); cont.update({'grid_score_primary':'no','continuous_score_primary':'yes','cont_score_rec_filename':'receptor.mol2','cont_score_att_exp':6,'cont_score_rep_exp':12,'cont_score_rep_rad_scale':1,'cont_score_use_dist_dep_dielectric':'yes','cont_score_dielectric':4,'cont_score_vdw_scale':1,'cont_score_es_scale':1,'ligand_outfile_prefix':'native_continuous','score_threshold':1000000})
(p/'native_continuous.in').write_text(config(D/'install/test/grid_score/orientflexiblemulti.dockin',cont))
run([D/'bin/dock6','-i','native_continuous.in','-o','native_continuous.out'],'native_continuous.log')
results['continuous']=(p/'native_continuous_scored.mol2').read_text().split('@<TRIPOS>MOLECULE')[0]
for spacing in [0.25,0.2]:
 tag=str(spacing).replace('.','_'); prefix='grid_'+tag
 gridparams={'receptor_file':'receptor.mol2','box_file':'box.pdb','vdw_definition_file':D/'parameters/vdw_AMBER_parm99.defn','score_grid_prefix':prefix,'contact_score':'no','grid_spacing':spacing}
 (p/(prefix+'.in')).write_text(config(D/'install/test/grid_generation/grid.in',gridparams))
 run([D/'bin/grid','-i',prefix+'.in','-o',prefix+'.out'],prefix+'.log')
 b=base.copy(); b.update({'grid_score_grid_prefix':prefix,'ligand_outfile_prefix':'native_'+tag})
 (p/('native_'+tag+'.in')).write_text(config(D/'install/test/grid_score/orientflexiblemulti.dockin',b))
 run([D/'bin/dock6','-i','native_'+tag+'.in','-o','native_'+tag+'.out'],'native_'+tag+'.log')
 results[str(spacing)]=(p/('native_'+tag+'_scored.mol2')).read_text().split('@<TRIPOS>MOLECULE')[0]
(p/'scoring_diagnostics.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
