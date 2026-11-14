"""Make the final report from saved, verified complete results only."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np


def main(out):
    evaluation=out/'evaluation';data=json.loads((evaluation/'metrics.json').read_text());v=json.loads((evaluation/'verification.json').read_text())
    assert v['passed'];profiles=data['profiles'];protocol=data['protocol']
    vina_auc=[profiles[p+'_vina']['methods']['higher_is_ranked_first']['auc'] for p in protocol['targets']]
    vinardo_auc=[profiles[p+'_vinardo']['methods']['higher_is_ranked_first']['auc'] for p in protocol['targets']]
    all_intervals_include_chance=all(p['methods']['higher_is_ranked_first']['paired_block_bootstrap_95'][0]<=.5<=p['methods']['higher_is_ranked_first']['paired_block_bootstrap_95'][1] for p in profiles.values() if p['status']=='complete')
    lines=['# AIPHARMA: beş bloklu kontrol ve aday değerlendirmesi — protokol 20261018','',
       '## Ana bulgu','',
       '**Bu test, mevcut docking skorunun ALK5 için güvenilir bir etkinlik sıralaması sağladığını doğrulamadı.** Vina ROC-AUC değerleri '+ ' / '.join(f'{a:.3f}' for a in vina_auc)+', Vinardo değerleri '+ ' / '.join(f'{a:.3f}' for a in vinardo_auc)+' (sırasıyla 1PY5 / 1VJY / 2WOT). Vina bütün yapılarda 0.5 altındadır; Vinardo ilk yapıda daha olumlu görünse de sonuç diğer yapılarda tekrarlanmamıştır. Bu küçük panel, skorlama fonksiyonunun hedefte hiçbir bilgi taşımadığını genel olarak kanıtlamaz.','',
       ('Altı koşulun tümünde beş-blok bootstrap %95 aralığı 0.5 değerini içeriyor. ' if all_intervals_include_chance else '')+'Örneklem yalnızca beş aktif ve beş zayıftır; yöntem üstünlüğü veya genellenebilir doğrulama iddiası yapılmamalıdır.','',
       'Aday 19 Vina ile üç yapıda 11/13, Vinardo ile sırasıyla 2/13, 3/13 ve 5/13: önceliği yönteme duyarlıdır. Aday 18 ve 02 bütün koşullarda alt sıralardadır. Üç aday da kontrol grubunun MW, logP ve ağır atom min–maks aralıklarının dışındadır; bu karşılaştırma adayların etkisizliğini veya etkinliğini kanıtlamaz.','',
       '## Araç ne yapıyor; bu tur neyi test ettik?','',
       'AIPHARMA molekül üretimi, yapı/cep bilgisi, docking ve ADMET tahminlerini birleştiren tasarım ve önceliklendirme aracıdır. Bu tur yalnızca docking tabanlı sıralama ve poz duyarlılığını inceler. ADMET güvenliliği, yeni moleküllerin biyolojik etkinliği veya klinik etki test edilmedi.','',
       'Yapı listesi ve değerlendirme kuralları hesaplamadan önce kilitlendi: beş kaynak/iskelet bloğunda 5 aktif + 5 zayıf kontrol ve önceki aday 19, 18, 02. Bazı eski kontrollerin verileri/sonuçları daha önce görüldüğü için çalışma bağımsız dış/kör doğrulama değildir. Yeni kontrol seçimi kaynak, yapı ve fizikokimyasal özelliklere dayanır; bu turun docking skorlarına dayanmaz.','',
       'Üç reseptör (1PY5, 1VJY, 2WOT), iki skorlama fonksiyonu (Vina/Vinardo 1.2.7), üç tohum (20261018–20261020), exhaustiveness=8 ve CPU/iş=1. Reseptör dosyaları ve 24 Å kutular önceki arşivden değişmeden alındı. Ligand hazırlığı ETKDGv3 16 konformer, MMFF94s en fazla 1000 iterasyon, yalnızca yakınsayan en düşük enerjili konformer, Meeko rijit makrosiklus seçeneği.','',
       f"**{data['completed_jobs']}/{data['scheduled_jobs']} hesap tamamlandı.** 234 çalışma hesabı, 18 native kontrol. Native pozların **{v['native_pass_count']}/18** tanesi kristal pozu geri bulma ölçütünü (simetri düzeltilmiş ağır atom RMSD <2 Å, süperpozisyon yok) geçti. Başarısız kontroller saklandı; hiçbir tohum/reseptör başarıya göre seçilip çıkarılmadı.",'',
       '## Kontrollerde ayırt edicilik','',
       'Her reseptör ve yöntem için üç tohumun ortalama skoru kullanıldı. Reseptörler veya skorlama fonksiyonlarının ham değerleri tek afinite gibi birleştirilmedi. AUC doğruluk yüzdesi değildir; 0.5 rastgele sıralama referansıdır.','',
       '| Reseptör | Yöntem | ROC-AUC | Beş-blok bootstrap %95 | Doğru çift / 5 | Native poz geçen / 3 |',
       '|---|---|---:|---|---:|---:|']
    for key,p in profiles.items():
        pdb,sf=key.split('_',1)
        if p['status']!='complete':lines.append(f"| {pdb} | {sf} | eksik | — | — | {p['native_pose_passes']} |");continue
        m=p['methods']['higher_is_ranked_first'];lo,hi=m['paired_block_bootstrap_95'];correct=sum(m['within_block_results'])
        lines.append(f"| {pdb} | {sf} | {m['auc']:.3f} | {lo:.3f}–{hi:.3f} | {correct:g} | {p['native_pose_passes']} |")
    lines+=['','Aralıklar yalnızca seçilmiş beş blok ve sabit hesaplar üzerinde koşullu bootstrap aralıklarıdır; veri kürasyonu, eski sonuçların bilinmesi ve assay farklılıklarının tüm belirsizliğini kapsamaz. Tek yönlü blok-değişim p değerleri metrics.json içinde keşifsel ve çoklu karşılaştırma düzeltmesi yapılmadan saklanmıştır; anlamlılık/üstünlük iddiası kurulmadı.','',
            '## Boyut ve lipofilisite referansları','',
            '| Reseptör / yöntem | Ham docking AUC | Skor / ağır atom AUC | MW yüksek | MW düşük | logP yüksek | Ağır atom yüksek | Ağır atom düşük |',
            '|---|---:|---:|---:|---:|---:|---:|---:|']
    for key,p in profiles.items():
        if p['status']!='complete':continue
        methods=p['methods'];keys=['higher_is_ranked_first','score_per_heavy_atom','mw','smaller_mw','logp','heavy_atoms','smaller_heavy_atoms']
        lines.append('| '+key+' | '+' | '.join(f"{methods[k]['auc']:.3f}" for k in keys)+' |')
    lines+=['','Özellik yönleri açıkça raporlandı; iyi görünen yön seçilerek referans performansı saklanmadı. Kontroller kendi aralarında özellik-eşlenmiştir, üç adayın kontrol dağılımında bulunması garanti edilmez.','',
            '## Üç adayın tamamı','',
            'Sıra, aynı reseptör ve yöntem içinde 13 molekülün üç-tohum ortalama docking skoruna göredir. 1 en negatif skordur; en etkili ilaç anlamına gelmez. Eşit ortalamalarda ortak sıra kullanılır. Kontrollerin karşılaştırıldığı deney koşulları, adayların aktivitesine atanmaz.','',
            '| Reseptör / yöntem | Aday | Ortalama skor (kcal/mol) | 13 içinde sıra | Tohum skor aralığı | En büyük tohumlar arası RMSD (Å) | Kontrol min–maks dışında |',
            '|---|---|---:|---:|---:|---:|---|']
    for key,p in profiles.items():
        for c in sorted(p['candidate_results'],key=lambda r:r['id']):
            if c.get('mean_score') is None:lines.append(f"| {key} | {c['id']} | eksik | — | — | — | — |");continue
            outside=', '.join(c['outside_observed_control_ranges']) or 'yok'
            lines.append(f"| {key} | {c['id']} | {c['mean_score']:.3f} | {c['rank_among_13']} | {c['seed_score_range']:.3f} | {c['max_seed_pose_RMSD_A']:.3f} | {outside} |")
    lines+=['','Adayların kontrol özellik aralıkları dışında kalması bir biyolojik başarısızlık kararı değildir; bu kontrol setinden adaylara doğrudan performans genellemesini sınırlar. Kaydedilen min–maks aralığı bir doğrulanmış applicability-domain modeli değildir. Adaylar olumsuz sıra aldıysa saklanmadı; olumlu sıra da deneysel bağlanma kanıtı sayılmadı.','',
            '## Native poz ve tohum duyarlılığı','',
            '| Reseptör / yöntem | Native RMSD, üç tohum (Å) | 13 çalışma molekülünde tohum poz farkı >2 Å | Tohum başına kontrol AUC |',
            '|---|---|---:|---|']
    for key,p in profiles.items():
        native=sorted(p['native_pose_controls'],key=lambda r:r['seed'])
        values=', '.join(f"{r['native_rmsd_A']:.3f}" if r.get('native_rmsd_A') is not None else 'eksik' for r in native)
        aucs=', '.join(f"{seed}: {value:.3f}" for seed,value in p.get('per_seed_auc',{}).items()) or 'eksik'
        lines.append(f"| {key} | {values} | {p['pose_variable_molecules']} | {aucs} |")
    lines+=['','Tohum poz değişkenliği, aynı başlangıç konformerinden farklı docking aramalarının duyarlılığıdır; moleküler dinamik değildir. >2 Å işareti poz değişkenliği olarak kullanıldı; otomatik etkisizlik veya evrensel yakınsamama kararı değildir. Native poz kurtarma arama/poz makullüğünü inceler; aktif–zayıf ayrımını tek başına doğrulamaz.','',
            '## Kaynak ve ölçüm sınırları','',
            'Beş blok, kaynak ailesi ve Murcko iskeleti/Morgan radius2-2048 benzerliği ≥0.4 bağlantılarını koruyan grafikte ayrı bileşenlerden oluşturuldu. Bu, biyolojik veya tam istatistiksel bağımsızlık ispatı değildir. Her çiftte MW farkı≤30 g/mol, logP farkı≤0.75, ağır atom farkı≤2 korundu. Kaynak/ATP/substrat farkları bütün bloklarda aynı değildir.','',
            'Ek aktif 667169 için gerçek kaynak 5–10 nM aralığıdır; BindingDB orta noktası 7.50 nM kesin ölçüm gibi kullanılmadı. Yeni zayıf kontrolün 10.000 nM sınır değerinde olma belirsizliği korunur. Deney kaynakları patent/literatür kayıtlarıdır; yeni laboratuvar ölçümü yoktur. Göreli stereo notasyonu taşıyan dört yeni zayıf kayıt ve üç yapı-kaynak belirsizliği olan kayıt seçilmedi.','',
            '## Araç ve makale için karar','',
            'Bu turun bütün reseptör/yöntem sonuçları birlikte raporlandı. Bir başarılı reseptör veya AUC seçilerek genel doğrulama iddiasına dönüştürülmedi. Üretim sıralaması bu küçük keşifsel karşılaştırma üzerinden değiştirilmedi. Mevcut docking skoru doğrulanmış bir ALK5 aktivite filtresi gibi kullanılmamalıdır. Adayların biyolojik etkinliği veya güvenliliği onaylanmadı.','',
            'İncelenen koşullarda ayrım, poz kurtarma ve tohum duyarlılığı ayrı kanıt katmanları olarak sunulabilir. Yazılımın hesapları yürütmesi ayrı bir işlevsellik bulgusudur; bunun için aday etkinliği iddiasına gerek yoktur.','',
            '## Dosyalar ve yeniden üretim','',
            '- `protocol.json`, `structures.json`, `label_key.json`: skor öncesi kilitli seçim ve kurallar. Protokolün özgün `locked_not_started` alanı değiştirilmedi; gerçek tamamlanma `run_status.json` ve `summary.json` içindedir.','- `input_hashes.json`, `receptors/`, `ligands/`: dosya özetleri, kutular, hazırlık durumları ve yapılar.','- `jobs/`: 252 işin pozları, ham enerjileri ve varsa hataları.','- `evaluation/metrics.json`, `scores.csv`, `verification.json`: tam metrikler ve denetimler.','- `code/`: çalıştırma, değerlendirme ve raporlama betikleri.','',
            'Çalıştırma: ayrı kurulmuş Vina1.2.7 Python klasörü PYTHONPATH üzerinde, mevcut RDKit/Meeko/spyrmsd ortamı içinde run betiğine --source önceki ensemble, --locked kilitli panel, --output yeni klasör verilir. Tamamlanmış çalışma üstüne yeni docking yazılmaz. Değerlendirme yalnızca summary.json oluştuğunda çalışır. Skorlama etiket anahtarını ayrıştırmaz; buna rağmen geçmiş etiketlerin bilinmesi nedeniyle bağımsız körlük iddiası yapılmaz.']
    lines += ['', 'Ek denetim dosyaları: `evaluation/controls_provenance.csv` kontrolleri okunabilir blok adları ve kaynaklarıyla eşler. `native_rmsd_independent_audit.json` tüm 18 native RMSD için RDKit in-place kontrolünü içerir; 0.01 Å içinde sayısal uyum ve aynı geçer/kalır kararları elde edilmiştir. İlk otomatik bağ-türü eşleştirmesinin hataları initial dosyasında korunur; açık element/bağlantı eşlemeleriyle karşılaştırma tamamlanmıştır. `design/` kilitli tasarımın kaynak grafiğini, `runtime_versions.json` ortam sürümlerini içerir.']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    # Standard scientific plotting; all points are saved-data results, no generated artwork.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figures=out/'figures';figures.mkdir(exist_ok=True)
    names=list(profiles);x=np.arange(len(names));fig,ax=plt.subplots(figsize=(10,4.5))
    for j,seed in enumerate(protocol['seeds']):
        values=[next((r.get('native_rmsd_A',np.nan) for r in profiles[k]['native_pose_controls'] if r['seed']==seed),np.nan) for k in names]
        ax.scatter(x+(j-1)*.16,values,label=str(seed),s=45)
    ax.axhline(2,color='#b91c1c',linestyle='--',label='2 A pose criterion');ax.set_xticks(x,names,rotation=20,ha='right');ax.set_ylabel('Native heavy-atom RMSD (A; no alignment)');ax.set_title('Native pose recovery: all 18 controls');ax.legend(fontsize=8,ncol=4);fig.tight_layout()
    fig.savefig(figures/'native_pose_controls.png',dpi=220);fig.savefig(figures/'native_pose_controls.svg');plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.5))
    for candidate in ['candidate_19','candidate_18','candidate_02']:
        ranks=[next(r.get('rank_among_13',np.nan) for r in profiles[k]['candidate_results'] if r['id']==candidate) for k in names]
        ax.plot(x,ranks,marker='o',label=candidate)
    ax.set_xticks(x,names,rotation=20,ha='right');ax.set_yticks(range(1,14));ax.set_ylim(13.5,.5);ax.set_ylabel('Rank of mean docking score among 13');ax.set_title('All fixed candidates; rank 1 = most negative score, not measured activity');ax.legend();fig.tight_layout()
    fig.savefig(figures/'candidate_ranks.png',dpi=220);fig.savefig(figures/'candidate_ranks.svg');plt.close(fig)
    with (out/'README.md').open('a',encoding='utf-8') as f:f.write('\n![Native poz kontrolleri](figures/native_pose_controls.png)\n\n![Üç adayın bütün sıraları](figures/candidate_ranks.png)\n')
    print('Report and two scientific figures saved',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
