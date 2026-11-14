"""Independent selection/arithmetic audit and readable method-comparison table."""
import json,statistics,sys
from collections import Counter,defaultdict
from pathlib import Path
p=Path(sys.argv[1]); out=p/'nested_method_search'
def read(f): return json.loads(f.read_text())
protocol=read(out/'protocol.json'); summary=read(out/'summary.json'); predictions=read(out/'nested_predictions.json')
references=read(p/'reference_delta/predictions.json'); lookup={r['index']:r for r in references}
assert len(predictions)==len(references)==501
assert {r['index'] for r in predictions}==set(lookup)
options=[c['name'] for c in protocol['options']]; assert len(options)==14
method_errors=defaultdict(list); selection_counts=Counter()
for f in sorted(out.glob('component_*.json')):
    c=read(f); chosen=min(options,key=lambda name:(c['inner_scores'][name],options.index(name)))
    assert c['selected']==chosen; selection_counts[chosen]+=1
    ids=c['outer_query_indices']; assert all(lookup[i]['component']==c['held_component'] for i in ids)
    assert set(c['outer_all_method_predictions'])==set(options)
    for method,z in c['outer_all_method_predictions'].items():
        assert len(z)==len(ids)
        method_errors[method].append(statistics.mean(abs(value-lookup[i]['actual_pic50']) for value,i in zip(z,ids)))
    for r in [r for r in predictions if r['component']==c['held_component']]:
        j=ids.index(r['index'])
        assert r['selected_method']==chosen and r['prediction']==c['outer_all_method_predictions'][chosen][j]
        assert r['actual_pic50']==lookup[r['index']]['actual_pic50'] and r['anchor_pic50']==lookup[r['index']]['anchor_pic50']
        assert r['anchor_index']==lookup[r['index']]['anchor_index']
groups=defaultdict(list)
for r in predictions: groups[r['component']].append(r)
mae=statistics.mean(statistics.mean(abs(r['prediction']-r['actual_pic50']) for r in rs) for rs in groups.values())
assert abs(mae-summary['nested_component_MAE'])<1e-12
assert dict(selection_counts)==summary['selection_counts']
ranked=sorted([(m,statistics.mean(v)) for m,v in method_errors.items()],key=lambda r:r[1])
for m,v in ranked:
    saved=next(r for r in summary['exploratory_configuration_leaderboard'] if r['method']==m)
    assert abs(saved['component_mean_MAE']-v)<1e-12
passed=mae<=.7 and summary['nested_gain']>=.1 and summary['conditional_component_bootstrap_gain_95'][0]>0
assert passed==summary['engineering_acceptance_met']
verify=dict(status='passed',queries=501,options=14,checks=['every outer choice minimizes its recorded inner scores','outer targets/anchors match the frozen 501-query task','selected outer predictions match the chosen method','14 configuration mean errors independently recomputed','nested error and acceptance gate recomputed'],scope='Saved output arithmetic and choices checked; fitted training separation enforced by runtime assertions and source-group input. Not biological validation.')
(out/'verification.json').write_text(json.dumps(verify,indent=2))
lines=['# ALK5 yontem karsilastirmasi','',
 '14 secenek: degisim-yok referansi, 4 Ridge ayari, 9 deney-farkli Gaussian process ayari. Ayni 501 olculmus hedef/referans cifti kullanildi.',
 '','Her dis kaynak grubunun etiketleri yontem secimine katilmadi. Secim, kalan kaynaklarin uc ic bolumundeki hataya gore yapildi. Baslangic referansinin olculmus etkinligi hem modele hem karsilastirmaya tahmin girdisi olarak verildi.',
 '',f"Secim algoritmasinin ayri dis kaynaklardaki grup-ortalama hatasi: **{mae:.5f} pIC50**. Degisim-yok referansi: **{summary['baseline_component_MAE']:.5f}**.",
 f"Kazanc: {summary['nested_gain']:.5f}; sabit tahminlere kosullu kaynak-bootstrap %95 araligi: {summary['conditional_component_bootstrap_gain_95']}.",
 f"Onceden belirlenmis muhendislik kabul kosullari: **{'gecti' if passed else 'gecmedi'}**.",
 '', '## Seceneklerin kesif amacli tablosu','',
 'Bu tablodan sonradan secilen en iyi satir icin bagimsiz dogrulama iddia edilmez. Daha guvenilir karsilastirma yukaridaki icte secilip dista sinanan politikanin sonucudur.',
 '', '| Secenek | Grup-ortalama mutlak hata (pIC50) |', '| --- | ---: |']
lines += [f'| {m} | {v:.5f} |' for m,v in ranked]
lines += ['', '## Ic secimler','']+[f'- {m}: {n} dis grupta secildi.' for m,n in sorted(selection_counts.items())]
lines += ['', 'Ayni geriye donuk veri daha once incelendi. Bu bir gelistirme karsilastirmasidir; kor, prospektif veya yeni molekul biyolojik dogrulamasi degildir. Hicbir aday etkinligi ya da guvenliligi bu arama ile dogrulanmadi. Tum basarisiz seceneklerin ic hata puanlari ve dis tahminleri dosyalarda korunur.']
(out/'README.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(dict(verification=verify,status={k:v for k,v in summary.items() if k!='exploratory_configuration_leaderboard'},best_exploratory_options=ranked[:3]),indent=2))
