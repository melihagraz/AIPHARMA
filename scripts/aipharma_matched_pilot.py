#!/usr/bin/env python3
"""Lock and audit a single matched AIPHARMA-vs-DockStream Vina pilot.

Does NOT silently run DockStream or infer compatible configurations.
A 'matched comparison' is reported only after an externally executed
DockStream output and protocol are supplied.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,subprocess
from datetime import datetime,timezone
from pathlib import Path

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
 return h.hexdigest()
def readcsv(p):
 with p.open('r',newline='',encoding='utf-8-sig') as f:
  r=csv.DictReader(f);return list(r),r.fieldnames or []
def writecsv(p,rows,cols):
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
def roc_auc_simple(labels,scores):
 positives=[x for x,y in zip(scores,labels) if y==1]
 negatives=[x for x,y in zip(scores,labels) if y==0]
 if not positives or not negatives:return None
 return sum(1 if a>b else 0.5 if a==b else 0 for a in positives for b in negatives)/(len(positives)*len(negatives))
def spearman(a,b):
 if len(a)<3:return None
 def ranks(arr):
  order=sorted(range(len(arr)),key=lambda i:arr[i]);r=[0]*len(arr);j=0
  while j<len(arr):
   k=j+1
   while k<len(arr) and arr[order[k]]==arr[order[j]]:k+=1
   v=(j+k+1)/2
   for m in range(j,k):r[order[m]]=v
   j=k
  return r
 x,y=ranks(a),ranks(b);xm=sum(x)/len(x);ym=sum(y)/len(y)
 dx=math.sqrt(sum((v-xm)**2 for v in x));dy=math.sqrt(sum((v-ym)**2 for v in y))
 return sum((u-xm)*(v-ym) for u,v in zip(x,y))/dx/dy if dx and dy else None

def lock(args):
 panel=args.panel.resolve();receptor=args.receptor.resolve();out=args.out.resolve()
 rows,cols=readcsv(panel)
 needed=[args.id_col,args.smiles_col,args.label_col]
 if any(x not in cols for x in needed):raise ValueError(f'Panel lacks required columns {needed}; has {cols}')
 if not receptor.is_file():raise FileNotFoundError(receptor)
 seen=set();standard=[]
 for r in rows:
  ident=r[args.id_col].strip(); smi=r[args.smiles_col].strip();v=r[args.label_col].strip()
  if not ident or not smi or v not in {'0','1'}:raise ValueError(f'Invalid id/smiles/binary class at row {len(standard)+1}')
  if ident in seen:raise ValueError(f'Duplicate id: {ident}')
  seen.add(ident);standard.append({'compound_id':ident,'smiles':smi,'label':int(v)})
 npos=sum(r['label'] for r in standard)
 if args.expected_n and len(standard)!=args.expected_n:raise ValueError(f'Expected {args.expected_n} compounds; got {len(standard)}')
 if npos==0 or npos==len(standard):raise ValueError('Both binary classes required')
 out.mkdir(parents=True,exist_ok=True)
 writecsv(out/'panel_locked.csv',standard,['compound_id','smiles','label'])
 # DockStream input is label-free; label file must not be passed to scoring process.
 (out/'input_for_dockstream.smi').write_text(''.join(f'{r["smiles"]} {r["compound_id"]}\n' for r in standard))
 protocol={'task':'Single-target, same-compound-set ALK5 AutoDock Vina workflow comparison','locked_utc':datetime.now(timezone.utc).isoformat(),'panel_source':str(panel),'panel_source_sha256':digest(panel),'locked_panel_sha256':digest(out/'panel_locked.csv'),'docking_input_sha256':digest(out/'input_for_dockstream.smi'),'n':len(standard),'positive':npos,'negative':len(standard)-npos,'receptor':str(receptor),'receptor_sha256':digest(receptor),'box_center_A':list(args.center),'box_size_A':list(args.size),'vina_exhaustiveness':args.exhaustiveness,'seed':args.seed,'cpu_per_ligand':args.cpu,'docking_backend':'AutoDock Vina','protocol_equivalence_status':'NOT_YET_VALIDATED','notes':['Same input molecules do not guarantee same protonation, 3D embedding, receptor/charge preparation, Vina version or mode.','DockStream workflow completion and data-reconciliation benchmark; not an algorithm superiority or biological validation test.','All failures and missing molecules must remain in denominators.']}
 (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
 print(f'LOCKED: n={len(standard)} positive={npos} negative={len(standard)-npos}; receptor SHA256={protocol["receptor_sha256"][:12]}')
 print('SAVED:',out/'protocol.json')

def read_scoremap(path,idcol,scorecol):
 rows,cols=readcsv(path)
 if idcol not in cols or scorecol not in cols:raise ValueError(f'{path.name}: required {idcol},{scorecol}; available={cols}')
 d={};bad=[]
 for r in rows:
  i=r[idcol].strip()
  if not i:bad.append('blank_id');continue
  if i in d:raise ValueError(f'Duplicate score row for {i} in {path}')
  try: v=float(r[scorecol]);v=v if math.isfinite(v) else None
  except (ValueError,TypeError):v=None
  d[i]=v
 return d,bad

def compare(args):
 root=args.lock.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
 prot=json.loads((root/'protocol.json').read_text()); panel,cols=readcsv(root/'panel_locked.csv')
 for row in panel: row['label']=int(row['label'])
 if digest(root/'panel_locked.csv')!=prot['locked_panel_sha256']:raise ValueError('Locked panel hash changed')
 if not args.dockstream_config.exists():raise FileNotFoundError('DockStream executed configuration JSON required')
 cfg=args.dockstream_config.read_text(encoding='utf-8')
 if 'vina' not in cfg.lower(): raise ValueError('DockStream executed config must explicitly contain Vina backend')
 try: json.loads(cfg)
 except json.JSONDecodeError as exc: raise ValueError(f'Invalid JSON in DockStream config: {exc}') from exc
 config_warning='Could not verify actual DockStream Vina protocol equivalence automatically; manually compare prep, box, receptor and engine.'
 a,ab=read_scoremap(args.aipharma_scores,args.aipharma_id,args.aipharma_score)
 b,bb=read_scoremap(args.dockstream_scores,args.dockstream_id,args.dockstream_score)
 ids={r['compound_id'] for r in panel}
 outsiders={'aipharma':sorted(set(a)-ids),'dockstream':sorted(set(b)-ids)}
 if any(outsiders.values()):raise ValueError(f'Scores contain IDs outside locked panel: {outsiders}')
 rows=[];completeA=completeB=0;both=[]
 for r in panel:
  ident=r['compound_id']; x=a.get(ident);y=b.get(ident)
  completeA+=x is not None; completeB+=y is not None
  if x is not None and y is not None:both.append((r['label'],x,y))
  rows.append({'compound_id':ident,'label':r['label'],'aipharma_vina':x if x is not None else '', 'dockstream_vina':y if y is not None else '', 'both_success':int(x is not None and y is not None)})
 writecsv(out/'paired_compounds.csv',rows,['compound_id','label','aipharma_vina','dockstream_vina','both_success'])
 def auc_for(d):
  z=[(r['label'],d[r['compound_id']]) for r in panel if d.get(r['compound_id']) is not None]
  return roc_auc_simple([y for y,_ in z],[-x for _,x in z]) if z else None
 x=[v[1] for v in both];y=[v[2] for v in both];dy=[v[0] for v in both]
 result={'status':'EXECUTED_OUTPUTS_AUDITED_PROTOCOL_EQUIVALENCE_UNVERIFIED','locked_panel_sha256':prot['locked_panel_sha256'],'n_locked':len(panel),'aipharma_completed':completeA,'dockstream_completed':completeB,'both_completed':len(both),'aipharma_completion_fraction':completeA/len(panel),'dockstream_completion_fraction':completeB/len(panel),'aipharma_auc_on_completed':auc_for(a),'dockstream_auc_on_completed':auc_for(b),'aipharma_auc_on_joint_complete':roc_auc_simple(dy,[-v for v in x]) if both else None,'dockstream_auc_on_joint_complete':roc_auc_simple(dy,[-v for v in y]) if both else None,'vina_score_spearman_joint_complete':spearman(x,y) if both else None,'aipharma_score_file_sha256':digest(args.aipharma_scores),'dockstream_score_file_sha256':digest(args.dockstream_scores),'dockstream_config_sha256':digest(args.dockstream_config),'protocol_equivalence_warning':config_warning,'limits':['Differences in Vina scores need not imply differences in docking algorithms: ligand preparation and versions may differ.','Completed-subset AUC is subject to differential failure-selection bias.','No speed comparison unless separate timed runs and machine metadata recorded.','No biological or integrated-platform superiority claim.']}
 (out/'matched_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 md=['# Matched ALK5 DockStream comparison — descriptive audit','',f'Locked: {len(panel)} molecules; AIPHARMA scored {completeA}; DockStream scored {completeB}; paired complete {len(both)}.','',f'- AIPHARMA AUC (successful rows): {result["aipharma_auc_on_completed"]}',f'- DockStream AUC (successful rows): {result["dockstream_auc_on_completed"]}',f'- Spearman (paired complete): {result["vina_score_spearman_joint_complete"]}','','**Protocol equivalence not established by this script.** Verify Vina binary version, ligand preparation, receptor identity, grid, pH, seed, exhaustive search and CPU before using this as a controlled workflow comparison. Do not claim superiority.']
 (out/'matched_summary.md').write_text('\n'.join(md)+'\n')
 print(f'COMPARISON AUDIT: AIPHARMA={completeA}/{len(panel)} DockStream={completeB}/{len(panel)} paired={len(both)}')
 print('SAVED:',out/'matched_summary.md')

def inventory(args):
 repo=args.repo.resolve(); print('Repo:',repo)
 print('git_sha:',subprocess.run(['git','-C',str(repo),'rev-parse','HEAD'],capture_output=True,text=True).stdout.strip() if (repo/'.git').exists() else 'NOT GIT')
 print('Docker/DockStream installed?:', 'YES' if (repo/'docker.py').exists() and (repo/'examples').exists() else 'NO')
 for p in sorted((repo/'examples').rglob('*.json')) if (repo/'examples').exists() else []:
  try:txt=p.read_text(errors='replace').lower()
  except OSError:continue
  if any(x in txt for x in ('vina','autodock')): print('VINA_EXAMPLE',p.relative_to(repo))

def local_inventory(args):
 root=args.root.resolve()
 print('=== Candidate ALK5 CSV panels and score files ===')
 prefixes=['reports/phase11','reports/phase12','reports/validation_phase1','reports/alk5_']
 count=0
 for p in sorted((root/'reports').rglob('*.csv')) if (root/'reports').exists() else []:
  rel=p.relative_to(root).as_posix()
  if not (any(rel.startswith(q) for q in prefixes) or 'matched' in rel.lower() or 'dude' in rel.lower()):continue
  try:
   with p.open(encoding='utf-8-sig',newline='') as f: cols=next(csv.reader(f))
  except (OSError,StopIteration,UnicodeError):continue
  print(rel, 'COLUMNS:',','.join(cols[:20]))
  count+=1
  if count>=70:print('... remaining candidates not shown');break
 print('=== Receptor PDBQT files ===')
 for p in sorted(root.rglob('*.pdbqt')):
  if any(x in p.parts for x in ('.venv','node_modules','.git')):continue
  print(p.relative_to(root))
  count+=1
  if count>100:break

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 subs=parser.add_subparsers(dest='action',required=True)
 l=subs.add_parser('lock');l.add_argument('--panel',type=Path,required=True);l.add_argument('--receptor',type=Path,required=True);l.add_argument('--id-col',default='compound_id');l.add_argument('--smiles-col',default='smiles');l.add_argument('--label-col',default='label');l.add_argument('--center',type=float,nargs=3,required=True);l.add_argument('--size',type=float,nargs=3,required=True);l.add_argument('--seed',type=int,required=True);l.add_argument('--exhaustiveness',type=int,required=True);l.add_argument('--cpu',type=int,default=1);l.add_argument('--expected-n',type=int,default=100);l.add_argument('--out',type=Path,required=True)
 c=subs.add_parser('compare');c.add_argument('--lock',type=Path,required=True);c.add_argument('--aipharma-scores',type=Path,required=True);c.add_argument('--dockstream-scores',type=Path,required=True);c.add_argument('--dockstream-config',type=Path,required=True);c.add_argument('--aipharma-id',default='compound_id');c.add_argument('--dockstream-id',default='compound_id');c.add_argument('--aipharma-score',default='score');c.add_argument('--dockstream-score',default='score');c.add_argument('--out',type=Path,required=True)
 i=subs.add_parser('dockstream-inventory');i.add_argument('--repo',type=Path,required=True)
 v=subs.add_parser('local-inventory');v.add_argument('--root',type=Path,default=Path('.'))
 args=parser.parse_args()
 try:
  if args.action=='lock':lock(args)
  if args.action=='compare':compare(args)
  if args.action=='dockstream-inventory':inventory(args)
  if args.action=='local-inventory':local_inventory(args)
 except (ValueError,FileNotFoundError) as e:parser.error(str(e))
