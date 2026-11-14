from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
p=Path.home()/'aipharma_holdout_20261016/endpoint_audit';o=p/'publication_figures';o.mkdir(exist_ok=True)
s=json.loads((p/'summary.json').read_text());m=json.loads((p/'comparison_summary.json').read_text());ci=json.loads((p/'uncertainty_supplement.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
fig,ax=plt.subplots(figsize=(11,5));ax.axis('off');ax.set(xlim=(0,11),ylim=(0,5))
boxes=[(0.2,3.3,3.1,1.25,'Archived ChEMBL records\nIC50 + assay descriptions\n153 eligible assays'),(3.95,3.3,3.1,1.25,'Separate endpoints\nCatalysis / binding / cells\nUncertain records retained'),(7.7,3.3,3.1,1.25,'Separate explicit mutants\nStandardize parent identity\nLabel within each stratum'),(7.7,.65,3.1,1.35,'Exclude old-test scaffolds\nGroup by document + scaffold\nFixed five-fold evaluation'),(3.95,.65,3.1,1.35,'Compare fixed methods\nFingerprint logistic model\nTraining-active similarity'),(.2,.65,3.1,1.35,'Report limitations\nSparse weak controls\nNo new-candidate validation')]
for x,y,w,h,t in boxes:ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.06',facecolor='#edf4fa',edgecolor='#2d5977'));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=10)
for start,end in [((3.35,3.9),(3.85,3.9)),((7.1,3.9),(7.6,3.9)),((9.25,3.2),(9.25,2.1)),((7.6,1.3),(7.15,1.3)),((3.85,1.3),(3.4,1.3))]:ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','lw':1.8,'color':'#2d5977'})
fig.tight_layout();fig.savefig(o/'figure1_workflow.png',bbox_inches='tight');plt.close(fig)
labels={'enzyme_catalytic':'Catalytic enzyme','biochemical_binding':'Biochemical binding','biochemical_unspecified':'Protein-based, unspecified','cell_signaling':'Cellular signaling','cell_target_engagement':'Cellular engagement','cell_phenotype':'Cell proliferation','uncertain':'Uncertain'}
names=list(labels);values=[s['assay_endpoint_counts'][x] for x in names];fig,ax=plt.subplots(figsize=(9,4.4));ax.barh(list(labels.values())[::-1],values[::-1],color=['#aab2bb' if x=='uncertain' else '#357ca5' for x in names[::-1]])
for j,v in enumerate(values[::-1]):ax.text(v+.5,j,str(v),va='center')
ax.set(xlim=(0,58),xlabel='Number of assay descriptions',title='Endpoint classification of 153 ALK5 assays');fig.tight_layout();fig.savefig(o/'figure2_endpoint_counts.png',bbox_inches='tight');plt.close(fig)
key='enzyme_catalytic__no_explicit_mutant';r=m[key]['evaluations']['document_and_scaffold'];c=ci['strata'][key]['document_and_scaffold']['fixed_OOF_group_bootstrap_AUC_95CI'];fig,axs=plt.subplots(1,2,figsize=(11,4.5));colors=['#1768ac','#d07a23']
for j,name in enumerate(['logistic','similarity']):
 v=r['methods'][name]['roc_auc'];lo,hi=c[name];axs[0].errorbar(j,v,yerr=[[v-lo],[hi-v]],fmt='o',capsize=6,color=colors[j],ms=8);axs[0].text(j-0.12,v,f'{v:.3f}',ha='right',color=colors[j])
axs[0].axhline(.5,color='gray',ls=':');axs[0].set(xticks=[0,1],xticklabels=['Logistic model','Similarity'],xlim=(-.6,1.6),ylim=(0,1.06),ylabel='ROC-AUC',title='Pooled out-of-fold estimates')
x=np.arange(1,6)
for j,name in enumerate(['logistic','similarity']):axs[1].plot(x,r['methods'][name]['fold_AUC'],'o-',label=['Logistic model','Similarity'][j],color=colors[j])
axs[1].axhline(.5,color='gray',ls=':');axs[1].set(xticks=x,xticklabels=[f'{k}\nweak n={f["test_weak"]}' for k,f in zip(x,r['folds'])],ylim=(0,1.07),ylabel='ROC-AUC',xlabel='Test fold',title='Performance varies across folds');axs[1].legend(loc='lower left',fontsize=9)
fig.suptitle('Catalytic enzyme stratum: 345 active / 11 weak; document + scaffold split',fontsize=11);fig.text(.5,.005,'Intervals resample fixed out-of-fold predictions by group; they do not include model retraining uncertainty.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.035,1,.96]);fig.savefig(o/'figure3_enzyme_performance.png',bbox_inches='tight');plt.close(fig)
(o/'figure_provenance.json').write_text(json.dumps({'source_sha256':{f:hashlib.sha256((p/f).read_bytes()).hexdigest() for f in ['summary.json','comparison_summary.json','uncertainty_supplement.json']},'figures':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in o.glob('*.png')},'note':'Plots use archived numeric results; no new fitted results or illustrative molecule activity claims.'},indent=2));print(str(o))
