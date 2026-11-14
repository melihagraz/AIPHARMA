"""Report a source-checked lead inventory without promoting it to a benchmark."""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json

p=Path.home()/'aipharma_external_leads_20261018'
screen=json.loads((p/'screen_summary.json').read_text())
leads=json.loads((p/'patent_table_review.json').read_text())
archives=json.loads((p/'patent_archive_log.json').read_text())
family={r['patent']:r['family_ids'][0] for r in archives if r['status']=='complete' and len(r['family_ids'])==1}
assert len(leads)==19 and len({r['normalized_parent'] for r in leads})==19
assert all(r['numeric_table_support_found'] and not r['accepted_for_benchmark'] for r in leads)
assert len(family)==5 and len(set(family.values()))==4
flat=[]
for i,r in enumerate(leads,1):
 flat.append({'lead_id':f'lead_{i:02d}','bindingdb_monomers':';'.join(r['bindingdb_monomers']),
              'reported_IC50_nM':';'.join(map(str,r['weak_values_nM'])),'relation':';'.join(r['relations']),
              'patents':';'.join(r['patent_sources']),'patent_family_ids':';'.join(sorted({family[s] for s in r['patent_sources']})),
              'parent_in_old_raw':False,'scaffold_in_old_raw':r['scaffold_in_prior_raw'],
              'numeric_ALK5_table_support':True,'structure_to_compound_mapping_independently_verified':False,
              'accepted_for_benchmark':False,'smiles':r['original_smiles'],'normalized_parent':r['normalized_parent']})
with (p/'new_leads_19.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
summary={**screen,'resolved_new_parent_patent_publications':5,'patent_family_ids':sorted(set(family.values())),
         'patent_family_count':4,'new_distinct_parent_leads':19,'numeric_ALK5_table_support_count':19,
         'structure_to_patent_compound_mapping_verified_count':0,'new_docking_jobs':0,'accepted_new_controls':0,
         'source_novelty_limit':'Literal patent IDs checked against old ChEMBL document registry; complete old/new patent-family overlap and source/scaffold block audit pending'}
(p/'final_summary.json').write_text(json.dumps(summary,indent=2))
lines=['# ALK5: dış kaynaklardan zayıf kontrol araması — 18 Ekim 2026','',
       '## Sonuç','',
       '**Önceki ChEMBL ham yapı listesinde bulunmayan 19 normalize ebeveyn yapı ipucu bulundu. Bunlar henüz kabul edilmiş 19 kontrol veya 19 bağımsız karşılaştırma değildir.** Sayısal zayıf IC50 değerlerinin tümü ilgili patentin ALK5 sütununda bulunuyor. Yapının patentteki bileşik numarasıyla birebir eşleşmesi ayrıca doğrulanmadığından bu kayıtlarla yeni doğrulama testi başlatılmadı.','',
       '## Aramanın kapsamı','',
       f"BindingDB P36897 sorgusu IC50 üst eşiği 100.000.000 nM ile yeniden yapıldı. Yanıtta {screen['bindingdb_total_measurement_rows']} ölçüm satırı, IC50 ≥10.000 nM veya uygun alt sınır olarak {screen['weak_IC50_rows']} satır vardı. Bunlar {screen['weak_distinct_parents']} normalize ebeveyne karşılık geldi; aynı molekülün tekrar ölçümü ayrı kontrol sayılmadı.",'',
       f"Eski ChEMBL havuzunda {screen['prior_raw_activity_rows']} ham aktivite, {screen['prior_parent_count']} normalize ebeveyn ve {screen['prior_document_count']} kaynak belgenin metaverisi karşılaştırıldı. Zayıf dış kayıtların 35 satırı eski ham ebeveyn listesinde bulunmuyordu; tekrarlar birleştirildiğinde 19 yapı kaldı.",'',
       'Eski BindingDB dosyası 481 ölçüm içerirken yeni geniş sorgu 5384 ölçüm döndürdü. Önceki dosyanın sorgu URL’si saklanmadığından farkın tamamının yalnızca eşikten kaynaklandığı söylenemez; eski API biçimi/sürümü ve kapsamı da farklıdır. Bununla birlikte düşük üst eşik, aranan zayıf kayıtları dışarıda bırakabileceği için yeni aramada geniş eşik açıkça kaydedildi.','',
       'Bu arama tamlık garantisi vermez; Ki/Kd, tek konsantrasyon inhibisyonu ve hücresel endpointler zayıf enzim IC50 yerine kullanılmadı. Büyük ölçüm sayısı bağımsız veri sayısı değildir.','',
       '## Kaynaklar ve tekrarlar','',
       '| Kaynak | Yeni yapı ipucu | Kaynak incelemesi |','|---|---:|---|',
       '| [US9682983B2](https://patents.google.com/patent/US9682983B2/en) ve [US10017516B2](https://patents.google.com/patent/US10017516B2/en) | 16 | Aynı patent ailesi (51625632); 32 kayıt tek 16 yapı olarak sayıldı. ALK1–6 enzim tablosunun ALK5 sütunu denetlendi. |',
       '| [US20240018122A1](https://patents.google.com/patent/US20240018122A1/en) | 1 | ALK5 ADP-Glo enzim yöntemi; sayısal tablo desteği var. |',
       '| [US20240116948A1](https://patents.google.com/patent/US20240116948A1/en) | 1 | ALK5 ADP-Glo; karşılaştırma tablosunda 10.000 nM. |',
       '| [US20240116947A1](https://patents.google.com/patent/US20240116947A1/en) | 1 | ALK5 ADP-Glo; karşılaştırma tablosunda 10.200 nM. |','',
       'Beş patent yayını dört patent ailesine karşılık geliyor. Patent ailesi sayısı, kimyasal kaynak/iskelet bağımsızlığı sayısı değildir. Benzer seriler ve ortak araştırma programları daha sıkı gruplamayı gerektirebilir. Eski kaynaklarla tam patent ailesi örtüşmesi henüz denetlenmedi; yalnızca birebir patent numarası eşleşmesi kontrol edildi.','',
       '10.7270 ile başlayan DOI’ler BindingDB veri kümesi kimlikleridir. Bunlar yeni makale gibi sayılmadı; DataCite IsDerivedFrom bağlantılarıyla özgün patentlere çözümlendi. DOI saymak yerine asıl kaynak birleştirildi.','',
       '## Yayın üzerinden çapraz kontrol','',
       '[Hao ve arkadaşlarının özgün yayını](https://pmc.ncbi.nlm.nih.gov/articles/PMC2825290/) ve ona bağlı [PubChem AID 1799383](https://pubchem.ncbi.nlm.nih.gov/bioassay/1799383/) ayrıca incelendi. Kaynakta saflaştırılmış insan kinazlarıyla enzim ölçümleri var. ALK5 hedef satırları diğer kinazlardan ayrıldı; PubChem standart değerlerinin µM, orijinal IC50 alanının nM olduğu kontrol edildi.','',
       'Bu assay’deki beş ALK5 yapı kaydı eski ham havuzda zaten bulunuyordu; yeni molekül diye eklenmedi. Kaynakta “No” veya “n.t.” yazan hücrelere sayısal IC50 uydurulmadı.','',
       '## Neden henüz yeni benchmark değil?','',
       '1. ALK5 sütununda aynı sayının bulunması, o sayının doğru kimyasal yapıya bağlandığını tek başına doğrulamaz. 19 yapının patent bileşik numarası ve yapı çizimi/MOL kaydıyla eşleştirilmesi gerekiyor.','2. Karşılaştırma değerleri ve yuvarlak 10.000 nM kayıtları için ölçümün kesin değer mi, sınır/üst test konsantrasyonu mu olduğu korunmalı; sınıflama otomatik genişletilmedi.','3. Yeni ve eski veriler birlikte kaynak/patent ailesi/iskelet grafiğine alınmalı. Dört patent ailesini otomatik dört bağımsız blok saymıyoruz.','4. Aynı assay içinde karşılaştırılabilir aktifler bulunup özellik eşlemesi skor görülmeden sabitlenmeli. ADP-Glo ile başka ATP/substrat koşulları açıklamasız birleştirilmemeli.','',
       'Bu adımlar tamamlanmadan yeni docking, yeni AUC veya aday etkinliği iddiası üretilmedi. Önceki 0.833/0.861 sonuçları bağımsız doğrulanmış performans gibi yeniden adlandırılmadı.','',
       '## Dosyalar','',
       '- `new_leads_19.csv`: 19 yapı, kaynaklar, değerler ve açık inceleme durumları.','- `patent_table_review.json`: sayısal tablo eşleşmeleri, aday bileşik numaraları; birebir yapı kanıtı değildir.','- `weak_lead_screen.json/.csv`: bütün 266 zayıf ölçüm ipucunun eski yapı/kaynak örtüşmesi.','- `sources/`: beş özgün patent HTML arşivi ve aranabilir metinleri.','- `datacite/`, `source_resolution.json`: veri kümesi–patent kaynak bağı.','- `known_document_registry.json`: eski 139 kaynak metaverisi.','- `fetch_log.json`: URL, başarı ve SHA256 özetleri.','',
       'Normalleştirme: RDKit Cleanup, FragmentParent, Uncharger, kanonik tautomer, stereo kaldırılarak kimlik karşılaştırması. Orijinal stereo içeren SMILES ayrıca korundu. Kimlik karşılaştırmasında stereo kaldırılması, belirsiz stereokimyayı docking için onaylamak anlamına gelmez.']
(p/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(p/'verification.json').write_text(json.dumps({'passed':True,'unique_new_parents':19,'numeric_table_support':19,
                                             'five_patent_archives_complete':True,'patent_family_count':4,
                                             'benchmark_promotions':0,'new_docking_jobs':0,
                                             'limits':'Numerical source support only; individual structure mapping and full independent-block audit pending'},indent=2))
print(json.dumps(summary,indent=2),flush=True)
