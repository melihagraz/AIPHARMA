#!/usr/bin/env python3
"""Read-only AIPHARMA release readiness audit. Does not license or publish anything."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, platform, re, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

SKIP_DIRS={'.git','.venv','.venv_tdc','venv','node_modules','.pytest_cache','__pycache__','dist','build','.mypy_cache'}
SKIP_FILES={'.env','.env.local','.env.production','.env.development','id_rsa','id_ed25519'}
WEIGHT_SUFFIX={'.pt','.pth','.ckpt','.safetensors','.onnx','.h5','.hdf5','.joblib','.pkl'}
CREDENTIAL_MARKERS=[
 ('github_token',re.compile(r'gh[pousr]_[A-Za-z0-9_]{25,}')),
 ('openai_key',re.compile(r'sk-[A-Za-z0-9_-]{30,}')),
 ('aws_access_key',re.compile(r'AKIA[0-9A-Z]{16}')),
 ('private_key',re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC |DSA )?PRIVATE KEY-----')),
 ('generic_secret',re.compile(r'(?im)^\s*(?:export\s+)?(?:API_KEY|SECRET_KEY|ADMIN_PASSWORD|SESSION_SIGNING_KEY|DATABASE_PASSWORD|ACCESS_TOKEN)\s*[:=]\s*["\']?([^\s"\'#]{12,})')),
]

def run(args,cwd):
 try:
  p=subprocess.run(args,cwd=cwd,capture_output=True,text=True,timeout=10)
  return p.stdout.strip() if p.returncode==0 else None
 except (OSError,subprocess.TimeoutExpired):return None

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def inventory(root):
 entries=[]
 for d,dirs,files in os.walk(root):
  dirs[:]=[x for x in dirs if x not in SKIP_DIRS]
  for fn in files:
   p=Path(d)/fn
   if not p.is_symlink() and p.is_file():
    try: entries.append((p.relative_to(root).as_posix(),p.stat().st_size,p))
    except OSError:pass
 return sorted(entries)

def audit(root,out):
 out.mkdir(parents=True,exist_ok=True)
 files=inventory(root);names={n for n,_,_ in files}
 is_repo=run(['git','rev-parse','--is-inside-work-tree'],root)=='true'
 commit=run(['git','rev-parse','HEAD'],root) if is_repo else None
 tracked=(run(['git','ls-files','-z'],root) or '').split('\x00') if is_repo else []
 tracked={x for x in tracked if x}
 status=run(['git','status','--porcelain'],root) if is_repo else None
 root_licenses=[n for n in names if '/' not in n and n.lower() in {'license','license.txt','license.md','copying'}]
 third_licenses=[n for n in names if n.lower().endswith(('/license','/license.txt','/license.md','/notice'))]
 weights=[{'path':n,'bytes':size,'sha256':sha(p)} for n,size,p in files if p.suffix.lower() in WEIGHT_SUFFIX]
 personal_files=[n for n in names if Path(n).name in SKIP_FILES or n.endswith(('.pem','.key','.p12','.pfx'))]
 suspicious=[]
 candidate_scan=[(n,p) for n,size,p in files if size<=500_000 and p.suffix.lower() in {'.py','.js','.ts','.tsx','.jsx','.json','.yml','.yaml','.toml','.sh','.ps1','.txt','.md','.env','.example'} and n not in personal_files and not n.startswith('reports/')]
 for n,p in candidate_scan:
  try: content=p.read_text(encoding='utf-8',errors='replace')
  except OSError:continue
  for label,pat in CREDENTIAL_MARKERS:
   matches=list(pat.finditer(content))
   if matches:
    # Do not reveal values; these may be dummy/example placeholders.
    suspicious.append({'path':n,'marker':label,'count':len(matches)})
 explicit_paths=[
  ('ALK5_1VJY_matched_panel','reports/phase11/alk5_matched50_docking'),
  ('ALK5_3HMM_matched_panel','reports/phase12/alk5_3hmm_matched50_docking'),
  ('ADMET_TDC','reports/validation_phase1/phase8_admet'),
  ('ADMET_transfer','reports/validation_phase1/phase10_external_admet_v2'),
  ('ADMET_robustness','reports/validation_phase1/phase11_admet_robustness'),
  ('EAVE_phase12','reports/eave/phase12_holdout'),
 ]
 availability=[]
 for label,rel in explicit_paths:
  p=root/rel
  availability.append({'artifact':label,'path':rel,'exists_in_checked_tree':p.exists(), 'tracked_files_count':sum(1 for t in tracked if t==rel or t.startswith(rel+'/')) if is_repo else None})
 blockers=[]
 if not root_licenses:blockers.append('No root-level AIPHARMA LICENSE detected: ownership and intended project license must be decided by authorized rights holder.')
 if not is_repo:blockers.append('Not an accessible Git checkout: no verifiable tracked commit for the checked snapshot.')
 if is_repo and status:blockers.append('Working tree has changes/untracked files; clean-clone evidence may differ from local execution.')
 if personal_files:blockers.append(f'{len(personal_files)} potential credential/private-key file(s) found; inspect locally before publishing.')
 if suspicious:blockers.append(f'{len(suspicious)} potential embedded secret pattern(s); manually review without sharing the values.')
 blockers.append('Third-party model weight rights and fetchable pinned versions are NOT verified by static file inventory.')
 blockers.append('Independent public clean-clone scientific replay is NOT verified by this audit.')
 blockers.append('Reported benchmark artifacts must be approved and deposited with checksums; local paths alone are not public links.')
 summary={'created_utc':datetime.now(timezone.utc).isoformat(),'root':str(root),'git':{'accessible':is_repo,'commit':commit,'clean_worktree':status=='' if is_repo else None,'tracked_count':len(tracked) if is_repo else None},'project_root_license':root_licenses,'third_party_license_files':third_licenses,'weight_files':weights,'private_path_candidates':personal_files,'potential_embedded_secrets':suspicious,'artifact_checks':availability,'blockers':blockers,'status':'NOT_PUBLICATION_READY','scan_scope':'Static scan only, not a full security or legal review. No keys printed or files published.'}
 (out/'release_gate_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
 with (out/'local_inventory.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['relative_path','bytes','sha256','hash_status'])
  for n,size,p in files:
   # Exclude potential private files from CSV, do not hash personal secrets.
   if n in personal_files:continue
   w.writerow([n,size,sha(p) if size<=32*1024*1024 else '', 'computed' if size<=32*1024*1024 else 'deferred_larger_than_32MiB'])
 md=['# AIPHARMA release audit','',f'Checked snapshot: `{root}`',f'Git commit: `{commit or "unavailable"}`',f'Root project LICENSE: `{", ".join(root_licenses) or "MISSING"}`','', '## Blocking conditions','']
 md += [f'- {x}' for x in blockers]
 md += ['','## Manuscript-associated artifact presence','', '| Artifact | Checked tree | Git-tracked file count |','|---|---|---:|']
 for x in availability:md.append(f'| {x["artifact"]} | {"present" if x["exists_in_checked_tree"] else "missing"} | {x["tracked_files_count"] if x["tracked_files_count"] is not None else "unknown"} |')
 md += ['','## Hash policy','Inventory hashes are computed for regular files at most 32 MiB; larger files are marked deferred and require explicit artifact hashing before release.','', '## Important scope','This audit neither uploads code nor makes a legal decision about licensing or weights. Running a test on an existing local environment is not equivalent to a public clean-clone replay. Potential embedded-secret matches require local manual review.']
 (out/'release_gate_summary.md').write_text('\n'.join(md)+'\n')
 print('AUDIT COMPLETE: NOT_PUBLICATION_READY')
 print('Root LICENSE:',root_licenses or 'MISSING')
 print('Git commit:',commit or 'UNAVAILABLE')
 print('Candidate private files:',len(personal_files),'embedded-secret patterns:',len(suspicious))
 print('Benchmark artifact directories present:',sum(x['exists_in_checked_tree'] for x in availability),'/',len(availability))
 print('SAVED:',out/'release_gate_summary.md')

if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--root',type=Path,default=Path('.'))
 ap.add_argument('--out',type=Path,default=Path('reports/release_gate'))
 a=ap.parse_args();root=a.root.expanduser().resolve();out=(root/a.out).resolve() if not a.out.is_absolute() else a.out.resolve()
 if not root.is_dir():ap.error('Root not found')
 if root==out or root in out.parents and out==root:ap.error('Output cannot equal root')
 audit(root,out)
