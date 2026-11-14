"""Independently reconstruct source-named compounds; preserve unresolved cases."""
from pathlib import Path
from functools import lru_cache
import concurrent.futures
import hashlib
import json
import re
import requests
from urllib.parse import quote
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

base=Path.home()/'aipharma_external_leads_20261018'
out=Path.home()/'aipharma_structure_mapping_20261018';out.mkdir(exist_ok=True)
specs={
 '111122':('US9682983B2','1','4-(6-(6-(piperazin-1-yl)pyridin-3-yl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '118138':('US9682983B2','4','4-(6-(3-(piperazin-1-yl)prop-1-ynyl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '127883':('US9682983B2','6','4-(6-(4-((1S,4S)-2,5-diazabicyclo[2.2.1]heptan-2-yl)phenyl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '143207':('US9682983B2','7','5-(piperazin-1-yl)-2-(3-(quinolin-4-yl)pyrazolo[1,5-a]pyrimidin-6-yl)thiazole'),
 '148370':('US9682983B2','31','5-(6-(4-(2-(piperidin-1-yl)ethoxy)phenyl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '148387':('US9682983B2','36','5-(6-(3-chloro-4-(piperazin-1-yl)phenyl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '148388':('US9682983B2','37','2-(piperazin-1-yl)-5-(3-(quinolin-4-yl)pyrazolo[1,5-a]pyrimidin-6-yl)benzonitrile'),
 '160478':('US9682983B2','41','4-(6-(6-(2-(piperazin-1-yl)ethoxy)pyridin-3-yl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '161911':('US9682983B2','43','4-(6-(6-(2-(4-methylpiperazin-1-yl)ethoxy)pyridin-3-yl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '161913':('US9682983B2','45','4-(6-(3-chloro-4-(2-(piperidin-1-yl)ethoxy)phenyl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '161922':('US9682983B2','54','5-(6-(4-((2R,6S)-2,6-dimethylpiperazin-1-yl)phenyl)pyrazolo[1,5-a]pyrimidin-3-yl)-2-methylquinoline'),
 '161926':('US9682983B2','58','2-methyl-5-(6-(5-((2-methyl-1-(piperidin-1-yl)propan-2-yl)oxy)pyridin-2-yl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '161934':('US9682983B2','66','(S)-2-methyl-5-(6-(5-(pyrrolidin-3-yloxy)pyridin-2-yl)pyrazolo[1,5-a]pyrimidin-3-yl)quinoline'),
 '161935':('US9682983B2','67','5-(6-(4-((2S,6R)-2,6-dimethylpiperazin-1-yl)phenyl)pyrazolo[1,5-a]pyrimidin-3-yl)-2-methoxyquinoline'),
 '667163':('US20240116948A1','C1','6-(5-chloro-2-fluorophenyl)-8-(4-methoxyphenyl)-2H,3H,4H-pyrido[3,2-b][1,4]oxazine'),
 '667180':('US20240116947A1','C2','5-[7-(5-Chloro-2-fluorophenyl)-1H,2H,3H-pyrido[3,4-b][1,4]oxazin-1-yl]-1,3-thiazole'),
}
quarantine={
 '108945':'Numeric table row 48 is not automatically synthetic Compound 48-prime. Its reconstructed name differs from the BindingDB piperidinyl-phenyl structure; original row image/MOL identity unresolved.',
 '161923':'Source Compound 55 explicitly names 2-deuterated quinoline; BindingDB SMILES has no deuterium isotope. Do not silently erase isotope identity.',
 '645599':'Source C3 heading names a Boc-protected compound, but HCl deprotection procedure, measured mass and BindingDB free-amine structure point to deprotected material. Source inconsistency requires image/MOL review; do not rewrite the source name.'}
leads=json.loads((base/'patent_table_review.json').read_text())
byid={r['bindingdb_monomers'][0]:r for r in leads}
assert set(byid)==set(specs)|set(quarantine)
def squash(s):return re.sub(r'\s+','',s).lower()
for ident,(patent,number,name) in specs.items():
 text=(base/'sources'/f'{patent}.txt').read_text()
 assert squash(name) in squash(text),(ident,'name not found verbatim in source')
 assert any(m['patent']==patent and (m['row_identifier']==number or not m['row_identifier']) for m in byid[ident]['numeric_table_support'])
protocol={'scope':'Source name -> OPSIN structure -> normalized-parent identity versus BindingDB; no docking scores',
          'normalization':'RDKit Cleanup, FragmentParent, Uncharger, canonical tautomer; stereochemistry retained for stereo comparison and removed only for parent identity',
          'name_rule':'Manually transcribed final-product name must occur verbatim modulo whitespace/case in archived original patent; counterion TFA removed explicitly',
          'quarantine_rule':'No inferred correction of source numbering, isotopes or protection-state inconsistencies',
          'input_sha256':hashlib.sha256((base/'patent_table_review.json').read_bytes()).hexdigest(),
          'services':['https://www.ebi.ac.uk/opsin/','https://www.ebi.ac.uk/opsin/ws/'],
          'benchmark_promotion':'Mapping alone does not establish independent source blocks or suitability for docking'}
(out/'protocol.json').write_text(json.dumps(protocol,indent=2))
(out/'transcribed_names.json').write_text(json.dumps({k:{'patent':v[0],'compound':v[1],'source_name':v[2]} for k,v in specs.items()},indent=2))
responses=out/'opsin';responses.mkdir(exist_ok=True)
def fetch(job):
 ident,(patent,number,name)=job;dest=responses/(ident+'.json');url='https://www.ebi.ac.uk/opsin/ws/'+quote(name,safe='')+'.json'
 try:
  if dest.exists():body=json.loads(dest.read_text())
  else:
   response=requests.get(url,timeout=30);body=response.json();dest.write_bytes(response.content)
  return {'id':ident,'url':url,'response':body}
 except Exception as e:return {'id':ident,'url':url,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:downloads=list(pool.map(fetch,specs.items()))
(out/'fetch_log.json').write_text(json.dumps(downloads,indent=2))
enum=rdMolStandardize.TautomerEnumerator();uncharger=rdMolStandardize.Uncharger()
def normalized(smi,stereo):
 m=Chem.MolFromSmiles(smi)
 if m is None:return None
 m=rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m));m=uncharger.uncharge(m);m=enum.Canonicalize(m)
 if not stereo:Chem.RemoveStereochemistry(m)
 return Chem.MolToSmiles(m,isomericSmiles=True)
results=[]
for download in downloads:
 ident=download['id'];lead=byid[ident];patent,number,name=specs[ident]
 record={**lead,'monomer':ident,'source_patent':patent,'source_compound':number,'source_name':name,'accepted_for_benchmark':False}
 body=download.get('response',{});smi=body.get('smiles')
 if smi and body.get('status') in ['SUCCESS','WARNING']:
  parent=normalized(smi,False);record.update({'opsin_smiles':smi,'opsin_status':body['status'],'opsin_warnings':body.get('warnings'),
    'source_parent':parent,'bindingdb_parent':normalized(lead['original_smiles'],False),
    'parent_identity_match':parent==normalized(lead['original_smiles'],False),
    'stereo_smiles_match':normalized(smi,True)==normalized(lead['original_smiles'],True),
    'mapping_status':'parent_verified' if parent==normalized(lead['original_smiles'],False) else 'structure_mismatch',
    'bindingdb_relative_stereo_annotation':'|r|' in lead['original_smiles']})
 else:record.update({'mapping_status':'parser_unresolved','error':download.get('error') or body.get('message') or body.get('status')})
 results.append(record)
for ident,reason in quarantine.items():results.append({**byid[ident],'monomer':ident,'mapping_status':'quarantined','reason':reason,'accepted_for_benchmark':False})
(out/'mapping_results.json').write_text(json.dumps(results,indent=2))
print(json.dumps([{'id':r['monomer'],'status':r['mapping_status'],'stereo_match':r.get('stereo_smiles_match'),'error':r.get('error')} for r in results],indent=2),flush=True)
