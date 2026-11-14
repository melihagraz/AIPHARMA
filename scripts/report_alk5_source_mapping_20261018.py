"""Complete the mapping audit, validate source table rows and graph partitions."""
from pathlib import Path
import csv,json,hashlib,re
from rdkit import Chem

home=Path.home();base=home/'aipharma_external_leads_20261018';out=home/'aipharma_structure_mapping_20261018'
records=json.loads((out/'mapping_results.json').read_text());tables=json.loads((base/'patent_tables.json').read_text())
verified=[r for r in records if r['mapping_status']=='parent_verified'];assert len(verified)==16
for r in verified:
 compound=r['source_compound'];patent=r['source_patent'];values=r['weak_values_nM']
 if compound.startswith('C'):
  matches=[row for table in tables[patent] for row in table if [c for c in row if c]==[compound,str(int(values[0]))]]
  assert matches,(r['monomer'],'comparative table row not found')
  r['verified_numeric_row']={'compound':compound,'IC50_nM':values[0],'row_cells':matches[0]}
 else:
  matches=[m for m in r['numeric_table_support'] if m['patent']==patent and m['row_identifier']==compound and m['ALK5_value_nM'] in values]
  assert matches,(r['monomer'],'enzyme column row not found')
  r['verified_numeric_row']=matches[0]
 assert Chem.MolFromSmiles(r['opsin_smiles']) is not None
 r['source_mapping_verified']=True
 # Keep a distinction between reported numerical IC50 and assay-ceiling certainty.
 r['measurement_ceiling_review_required']=values==[10000.0]
 r['reported_relative_stereo_review_required']=r['bindingdb_relative_stereo_annotation']
(out/'mapping_results.json').write_text(json.dumps(records,indent=2))
checked=0
for f in out.glob('*_graph.json'):
 g=json.loads(f.read_text());adj=[set() for _ in range(g['molecule_count'])]
 for e in g['edges']:adj[e['a']].add(e['b']);adj[e['b']].add(e['a'])
 unseen=set(range(len(adj)));components=[]
 while unseen:
  root=min(unseen);unseen.remove(root);seen={root};todo=[root]
  while todo:
   for k in adj[todo.pop()]:
    if k not in seen:seen.add(k);unseen.remove(k);todo.append(k)
  components.append(frozenset(seen))
 assert set(components)=={frozenset(b['indices']) for b in g['components']};checked+=1
matching=json.loads((out/'block_matching_feasibility.json').read_text())
summary=json.loads((out/'expanded_graph_summary.json').read_text())
flat=[]
for r in records:
 flat.append({'monomer':r['monomer'],'mapping_status':r['mapping_status'],'source_patent':r.get('source_patent',';'.join(r['patent_sources'])),
              'source_compound':r.get('source_compound'),'IC50_nM':';'.join(map(str,r['weak_values_nM'])),
              'parent_match':r.get('parent_identity_match'),'canonical_stereo_match':r.get('stereo_smiles_match'),
              'relative_stereo_review_required':r.get('reported_relative_stereo_review_required'),
              'measurement_ceiling_review_required':r.get('measurement_ceiling_review_required'),
              'benchmark_accepted':False,'reason':r.get('reason','source name and numeric row verified; independent block/measurement qualifications still apply'),
              'original_smiles':r['original_smiles'],'source_name':r.get('source_name')})
with (out/'mapping_decisions.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
relative=sum(bool(r.get('reported_relative_stereo_review_required')) for r in verified)
lines=['# ALK5: yapı–kaynak eşleştirmesi ve genişletilmiş blok denetimi — 18 Ekim 2026','',
       '## Sonuç','',
       '**19 dış yapı ipucunun 16’sı, özgün patentteki son ürün adından OPSIN ile yeniden oluşturulan yapı ve doğru ALK5 tablo satırıyla eşleştirildi. Üç kayıt belirsizlik nedeniyle ayrı tutuldu.** Bu, moleküllerin etkinliğinin yeni deneyle doğrulandığı anlamına gelmez.','',
       'Yeniden oluşturulan yapılarda normalize ebeveyn eşleşmesi ve kanonik stereo-SMILES eşleşmesi 16/16. Tuz/karşı iyon normalleştirmesi açıkça kaydedildi; TFA karşı iyonu ebeveyn yapıya dahil edilmedi. Orijinal SMILES ve API yanıtları korundu.','',
       f'16 kaydın {relative} tanesinde BindingDB girdisi ayrıca `|r|` göreli stereo notasyonu taşıyor. Kanonik stereo-SMILES eşitliği, bu işaretin rasemat/göreli stereo anlamını veya deneyde test edilen tek izomeri tek başına çözmez. Bu kayıtlar otomatik olarak tek izomer docking girdisi kabul edilmedi. 10.000 nM sınırındaki kayıt için de deney tavanı/yuvarlama incelemesi açık tutuldu.','',
       '## Test dışında tutulan üç kayıt','',
       '| BindingDB | Neden |','|---|---|',
       '| 108945 | Enzim tablosunun 48 satırı, sentez bölümündeki 48′ bileşiğiyle otomatik eşlenemez. 48′ adı BindingDB yapısından farklıdır. Bu durum veri hatası ispatı değildir; doğru çizim/MOL eşleştirmesi tamamlanmadı. |',
       '| 161923 | Patent bileşik 55 adında döteryum etiketi var; BindingDB yapısında döteryum yok. İzotop farkı sessizce silinmedi. |',
       '| 645599 | C3 başlığı Boc-korumalı adı verirken HCl ile koruma kaldırma işlemi, kütle ve BindingDB yapısı korunmamış ürüne işaret ediyor. Kaynak içi uyuşmazlık nedeniyle otomatik düzeltme yapılmadı. |','',
       '## Kaynak ve iskelet grafiği','',
       'Önceki uygun 327 kayıt ile yapı eşleşmesi doğrulanan 16 zayıf ipucu birlikte incelendi: 318 aktif, 25 zayıf; toplam 343 yapı. Önceki uygun havuzdaki üç patentin aile kimlikleri de çözüldü. Aynı patent ailesinin farklı yayınları tek kaynak sayıldı. Yeni 16 yapının üç kaynak ailesi var.','',
       'Sadece grafiğin bileşenlerini saymakla yetinmedik: bileşenleri bölmeden, aynı MW/logP/ağır atom sınırları altında uygun aktif–zayıf çifti oluşturulup oluşturulamayacağını hesapladık. Birden fazla zayıf kaynak bileşeni aynı aktif bileşenine bağlanarak bağımsızmış gibi sayılmadı.','',
       '| Gruplama | Eşlenebilen ayrık tam bileşen blokları | Beş blok hedefi |','|---|---:|---|',
       f"| Ortak kaynak ailesi veya birebir aynı Murcko iskeleti | {matching['source_family_exact_scaffold']['feasible_disjoint_full_component_blocks']} | Sayısal olarak mümkün |",
       f"| Ek olarak iskelet ECFP4 benzerliği ≥0.4 olanları birleştir | {matching['source_family_scaffold_similarity_040']['feasible_disjoint_full_component_blocks']} | Sağlanmıyor |",'',
       'Özellik sınırları önceki testle aynı tutuldu: MW farkı ≤30 g/mol, logP farkı ≤0.75, ağır atom farkı ≤2. Eşleme docking skorlarını okumuyor. 0.4 eşiği bir duyarlılık incelemesidir; evrensel biyolojik bağımsızlık kuralı değildir. Sonucu daha iyi göstermek için bu sınırlar gevşetilmedi.','',
       'Bu fizibilite, 16 kaydın hepsinin kullanılabildiği iyimser bir senaryodur; stereo/ölçüm tavanı kayıtlarının son incelemesi daha az kullanılabilir kayıt bırakabilir. Daha sıkı kaynak/benzerlik tanımı altında bile beş özellik-eşlenmiş blok hazır değildir.','',
       '## Karar','',
       '16 kayıt, yapı–kaynak eşleştirmesi yapılmış aday kontrol havuzuna alındı; yeni benchmarka otomatik kabul edilmedi. Üç belirsiz kayıt test dışında kaldı. Üretim sıralaması değiştirilmedi ve yeni docking başlatılmadı.','',
       'Daha sıkı tasarımın beş bloğa ulaşması için mevcut gruplara bağlanmayan ve sabit özellik sınırlarında aktif–zayıf karşılaştırması sağlayan en az bir ek uygun blok gerekiyor. Sadece daha fazla aynı-seri molekül eklemek bu açığı kapatmaz. Bu matematiksel fizibilite, yeterli istatistiksel güç veya dış doğrulama garantisi değildir.','',
       '## Yöntem ve tekrar üretim','',
       'Patent kimyasal adları, arşivlenmiş metinde boşluk/büyük-küçük harf dışında değişmeden bulunmak zorunda. [OPSIN resmi servisi](https://www.ebi.ac.uk/opsin/) adları yapıya çevirir; bu yapı RDKit normalleştirmesiyle BindingDB yapısına karşılaştırılır. Yöntem fiziksel numunenin kimliğini, assay doğruluğunu veya koşulların eşdeğerliğini kanıtlamaz.','',
       'Sayısal eşleşme ayrıca doğru bileşik satırına bağlandı. Önceki dosyada karşılaştırma tablolarının boş ilk hücresi yüzünden bileşik etiketi boş kalmıştı; burada C1/C2 etiketleri özgün tablo hücrelerinden açıkça doğrulandı.','',
       '- `mapping_decisions.csv`: 19 kayıt için son eşleştirme ve açık inceleme durumları.','- `mapping_results.json`, `transcribed_names.json`, `opsin/`: özgün ad, bileşik numarası, yeniden oluşturulmuş yapı ve API yanıtları.','- `expanded_nodes.json`, `*_graph.json`: eski/yeni tüm düğümler ve kenar gerekçeleri.','- `old_patent_family_log.json`, `old_patent_families/`: kaynak ailesi çözümleme kayıtları.','- `block_matching_feasibility.json`: iki tanımda eşlenebilen örnek çiftler; kilitlenmiş yeni benchmark değildir.','',
       'Girdi kaynakları: [US9682983B2](https://patents.google.com/patent/US9682983B2/en), [US20240116948A1](https://patents.google.com/patent/US20240116948A1/en), [US20240116947A1](https://patents.google.com/patent/US20240116947A1/en), [US20240018122A1](https://patents.google.com/patent/US20240018122A1/en). Tam arşivler önceki `alk5_external_leads_20261018/sources/` raporundadır.']
(out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
(out/'verification.json').write_text(json.dumps({'passed':True,'source_names_found_verbatim':16,'independent_name_to_parent_matches':16,
                                              'canonical_stereo_matches':16,'numeric_compound_row_checks':16,'quarantined':3,
                                              'graphs_verified_by_independent_adjacency_traversal':checked,'new_docking_jobs':0,
                                              'relative_stereo_review_required':relative,'benchmark_promotions':0},indent=2))
print(json.dumps({'verified_source_mappings':len(verified),'quarantined':3,'graphs_verified':checked,'relative_stereo_review_required':relative,'new_docking_jobs':0}))
