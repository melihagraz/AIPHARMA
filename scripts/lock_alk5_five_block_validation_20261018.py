"""Lock the reviewable follow-up panel; do not run or claim docking results."""
from pathlib import Path
import json,hashlib,re,csv
from rdkit import Chem

home=Path.home();out=home/'aipharma_structure_mapping_20261018'
panel=json.loads((out/'five_block_panel_feasibility.json').read_text());assert panel['maximum_matched_blocks_after_additional_active']==5
dest=out/'locked_followup';dest.mkdir(exist_ok=True)
if (dest/'protocol.json').exists():raise RuntimeError('Panel already locked')
structures=[];labels=[]
for pair in panel['pairs']:
 for kind in ['active','weak']:
  row=pair[kind];smiles=row.get('original_smiles') or row['smiles'];mol=Chem.MolFromSmiles(smiles)
  assert not any(str(x.specified)=='Unspecified' for x in Chem.FindPotentialStereo(mol))
  ident='control_'+hashlib.sha256(smiles.encode()).hexdigest()[:12]
  structures.append({'id':ident,'kind':'control','smiles':smiles})
  labels.append({'id':ident,'block':pair['block'],'label':int(kind=='active'),'provenance':row,'source_components':pair['source_components']})
old=json.loads((home/'aipharma_ensemble_20261018/locked_items.json').read_text())
structures.extend({k:r[k] for k in ['id','kind','smiles']} for r in old if r['kind']=='candidate')
assert len(structures)==13 and len({r['smiles'] for r in structures})==13
structures.sort(key=lambda r:r['id'])
(dest/'structures.json').write_text(json.dumps(structures,indent=2));(dest/'label_key.json').write_text(json.dumps(labels,indent=2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol={
 'purpose':'Exploratory five-source/scaffold-block follow-up; not independent external blinded validation',
 'status':'locked_not_started','controls':10,'blocks':5,'fixed_candidates':['candidate_19','candidate_18','candidate_02'],
 'history':'Some control structures/labels and scores were exposed previously; new paired panel selected only by source structure and physicochemical rules. No third-party blinding.',
 'eligibility':'Source-family plus exact Murcko or scaffold Morgan radius2/2048 Tanimoto >=0.4 connected components kept intact; 0.4 is a sensitivity threshold, not biological independence proof.',
 'matching':{'max_MW_difference':30,'max_logP_difference':.75,'max_heavy_atom_difference':2,'new_relative_stereo_records_excluded':True},
 'activity_labels':'Primary controls: old exact median pChEMBL >=7 or <=5; new weak reported IC50 >=10000nM; new active monomer667169 source interval5-10nM, not midpoint7.5 as an exact measurement.',
 'targets':['1PY5','1VJY','2WOT'],'functions':['vina','vinardo'],'vina_package_version':'1.2.7',
 'seeds':[20261018,20261019,20261020],'exhaustiveness':8,'cpu_per_job':1,
 'ligand_preparation':{'ETKDGv3_conformers':16,'seed':20261018,'MMFF94s_iterations':1000,'converged_only':True,'Meeko_rigid_macrocycles':True,'no_crystal_pose_as_start':True},
 'receptor_preparation':'Reuse archived ensemble receptors and24A boxes unchanged, preserve hashes; no score-based pocket or protonation tuning',
 'native_controls':'One native ligand per receptor, both engines, all three seeds; symmetry-corrected heavy-atom RMSD with atom mapping, no superposition, <2A pose recovery criterion. Failures retained and interpretations qualified.',
 'pose_sensitivity':'All pairwise seed pose RMSDs reported per molecule/receptor/engine; >2A flagged as pose variability, not proof of inactivity or a universal convergence failure.',
 'expected_jobs':{'study':234,'native':18,'total':252},
 'evaluation':'Mean seed score within each receptor/engine; do not select best seed/receptor; ROC-AUC, within-block ranking, 5-block resampling/permutation, MW/logP/HA and score-per-HA baselines; no cross-receptor affinity averaging.',
 'candidate_reporting':'All three fixed candidates shown beside all ten controls for every receptor/engine, including poor ranks; no post-score substitutions',
 'claim_limit':'No experimental activity/safety, no general performance or Vinardo superiority from this small exposed set; pose recovery and discrimination separate',
 'structures_sha256':sha(dest/'structures.json'),'label_key_sha256':sha(dest/'label_key.json'),
 'panel_feasibility_sha256':sha(out/'five_block_panel_feasibility.json')}
(dest/'protocol.json').write_text(json.dumps(protocol,indent=2))
readme=out/'README.md';text=readme.read_text()
text=text.replace('## Kaynak ve iskelet grafiği','## Ek aktif bulunmadan önceki kaynak ve iskelet grafiği',1)
text=text.replace('Daha sıkı tasarımın beş bloğa ulaşması için mevcut gruplara bağlanmayan ve sabit özellik sınırlarında aktif–zayıf karşılaştırması sağlayan en az bir ek uygun blok gerekiyor. Sadece daha fazla aynı-seri molekül eklemek bu açığı kapatmaz. Bu matematiksel fizibilite, yeterli istatistiksel güç veya dış doğrulama garantisi değildir.',
 'Bu ilk aşamada eksik olan, beşinci zayıf kaynak grubuna uygun bir aktifti. Aşağıdaki ek kaynak incelemesiyle bu açık kapatıldı; ilk aşamanın dosyaları karşılaştırma için korundu.')
prefix='''# Güncel karar: beş gruplu takip paneli hazır; henüz docking sonucu yok

16/19 dış zayıf yapı ipucunun kaynak adı–yapı–ALK5 satırı bağlantısı doğrulandı; üçü belirsizlik nedeniyle ayrı tutuldu. Beşinci kaynak grubunda eksik kalan aktif ayrıca özgün patentten doğrulandı: BindingDB 667169, US20240116947A1 Example 1. Patent gerçek bir kesin IC50 yerine **5–10 nM aralığı** raporlar. BindingDB 7,50 nM değerinin bu aralığın ortası olduğu görüldü; tam ölçüm gibi kullanılmadı.

Bu aktif eklenince 344 yapı üzerinde daha sıkı kaynak ailesi + iskelet benzerliği ≥0.4 grafiği yeniden hesaplandı. Aynı özellik sınırlarında beş ayrık tam bileşen bloğu kurulabildi. Göreli stereo notasyonu taşıyan dört yeni kayıt, çift seçiminden çıkarıldı. Her seçilen çift MW≤30 g/mol, logP≤0.75 ve ağır atom farkı≤2 sınırlarını sağlıyor.

**Hazır panel: 5 aktif + 5 zayıf kontrol ve önceki üç aday.** Aday 19/18/02 değiştirilmedi. Bazı eski kontrollerin geçmiş sonuçları bilindiği için bu panel yeni bağımsız kör doğrulama diye sunulamaz. Grafik ayrıklığı da biyolojik bağımsızlık ispatı değildir.

`locked_followup/` yapı listesini, kaynak anahtarını ve skorlar hesaplanmadan belirlenen takip protokolünü içerir. Plan: üç reseptör, Vina/Vinardo, üç tohum; 234 çalışma hesabı ve 18 native poz kontrolü, toplam 252 iş. **Bu dosya oluşturulurken işler başlatılmadı.** Yeni AUC, aday sırası veya yöntem üstünlüğü sonucu yoktur.

Kayıtlar ve değişiklikler aşağıdadır; ilk dört-blok fizibilite sonucu korunarak ek aktifle beşe ulaşılması ayrı belgelenmiştir.

---

'''
readme.write_text(prefix+text,encoding='utf-8')
verification=json.loads((out/'verification.json').read_text())
verification.update({'additional_active_source_parent_verified':True,'additional_active_reported_interval_nM':[5,10],
                     'matched_disjoint_blocks_after_active':5,'locked_controls':10,'locked_candidates':3,'locked_expected_jobs':252,'jobs_started':False})
(out/'verification.json').write_text(json.dumps(verification,indent=2))
print(json.dumps({'controls':10,'candidates':3,'blocks':5,'protocol_sha256':sha(dest/'protocol.json'),'jobs_started':False},indent=2))
