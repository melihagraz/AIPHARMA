"""Evaluate all locked matched controls, retaining failures and size-only baselines."""
import argparse
import csv
import hashlib
import itertools
import json
import re
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def auc(y,s):
    y,s=np.asarray(y),np.asarray(s)
    delta=s[y==1,None]-s[y==0]
    return float(np.mean((delta>0)+.5*(delta==0)))


def main(out,audit,rescoring):
    protocol=json.loads((out/'protocol.json').read_text())
    assert sha(out/'label_key.json')==protocol['label_key_sha256']
    assert sha(out/'structures.json')==protocol['structures_sha256']
    key=json.loads((out/'label_key.json').read_text());raw=json.loads((out/'results.json').read_text())
    assays={r['assay_id']:r for r in json.loads((audit/'assay_classification.json').read_text())}
    rows=[];review=[]
    for item in key:
        assert item['evidence']=='exact'
        assert (item['label']==1 and item['median_pchembl']>=7) or (item['label']==0 and item['median_pchembl']<=5)
        records=[r for r in raw if r['id']==item['id']]
        assert len(records)==2
        descriptions=[assays[a] for a in item['assay_ids']]
        assert all(a['endpoint']=='enzyme_catalytic' and a['mutation']=='no_explicit_mutant' for a in descriptions)
        review.append({'id':item['id'],'label':item['label'],'assays':descriptions})
        row={k:item[k] for k in ['id','pair','smiles','label','median_pchembl','mw','logp','heavy_atoms','documents','scaffold']}
        for sf in ['vina','vinardo']:
            result=next(r for r in records if r['function']==sf)
            if result['status']=='complete':
                pose=(out/'jobs'/f"{item['id']}_{sf}"/'poses.pdbqt').read_text()
                written_score=float(re.search(r'VINA RESULT:\s*([-\d.]+)',pose).group(1))
                assert abs(written_score-result['score'])<.0011
            row[sf+'_status']=result['status']
            row[sf]=(-result['score']) if result['status']=='complete' else None
            row[sf+'_per_heavy_atom']=row[sf]/item['heavy_atoms'] if row[sf] is not None else None
        row['smaller_heavy_atoms']=-item['heavy_atoms'];row['smaller_MW']=-item['mw'];row['logP']=item['logp']
        rows.append(row)
        if all(r['status']=='complete' for r in records):
            assert sha(out/'jobs'/f"{item['id']}_vina"/'ligand.pdbqt')==sha(out/'jobs'/f"{item['id']}_vinardo"/'ligand.pdbqt')
    pairs=sorted({r['pair'] for r in rows});pair_indices=[];balance=[]
    for pair in pairs:
        indices=[i for i,r in enumerate(rows) if r['pair']==pair]
        assert len(indices)==2 and {rows[i]['label'] for i in indices}=={0,1}
        pair_indices.append(indices)
        a,w=sorted([rows[i] for i in indices],key=lambda r:-r['label'])
        delta={k:abs(a[k]-w[k]) for k in ['mw','logp','heavy_atoms']}
        assert delta['mw']<=30 and delta['logp']<=.75 and delta['heavy_atoms']<=2
        balance.append({'pair':pair,'active':a['id'],'weak':w['id'],'differences':delta})
    y=np.array([r['label'] for r in rows]);metrics={}
    # A pair is not independent when it shares a source or scaffold with another pair.
    parent=list(range(len(pairs)))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i in range(len(pairs)):
        for j in range(i):
            ri=[rows[k] for k in pair_indices[i]];rj=[rows[k] for k in pair_indices[j]]
            di={d for r in ri for d in r['documents']};dj={d for r in rj for d in r['documents']}
            si={r['scaffold'] for r in ri};sj={r['scaffold'] for r in rj}
            if di&dj or si&sj:parent[root(i)]=root(j)
    components={}
    for i in range(len(pairs)):components.setdefault(root(i),[]).extend(pair_indices[i])
    groups=list(components.values())
    methods=['vina','vinardo','vina_per_heavy_atom','vinardo_per_heavy_atom','smaller_heavy_atoms','smaller_MW','logP']
    for method in methods:
        if any(r[method] is None for r in rows):
            metrics[method]={'status':'incomplete','reason':'No complete-case deletion allowed'};continue
        s=np.array([r[method] for r in rows]);observed=auc(y,s)
        assert abs(observed-roc_auc_score(y,s))<1e-12
        pair_wins=[]
        for indices in pair_indices:
            i,j=sorted(indices,key=lambda i:-y[i]);pair_wins.append(float(s[i]>s[j])+.5*float(s[i]==s[j]))
        null=[]
        for swaps in itertools.product([False,True],repeat=len(pairs)):
            yp=y.copy()
            for indices,swap in zip(pair_indices,swaps):
                if swap:yp[indices]=yp[indices[::-1]]
            null.append(auc(yp,s))
        rng=np.random.default_rng(20261018);boot=[]
        for _ in range(2000):
            ix=np.concatenate([pair_indices[k] for k in rng.integers(0,len(pairs),len(pairs))])
            boot.append(auc(y[ix],s[ix]))
        metrics[method]={'status':'complete','auc':observed,'within_pair_correct_fraction':float(np.mean(pair_wins)),
                         'pair_results':pair_wins,'paired_bootstrap_95':np.quantile(boot,[.025,.975]).tolist(),
                         'paired_label_swap_p_one_sided':float(np.mean(np.array(null)>=observed-1e-12)),
                         'paired_permutations':len(null)}
        grouped_null=[]
        for swaps in itertools.product([False,True],repeat=len(groups)):
            yp=y.copy()
            for indices,swap in zip(groups,swaps):
                if swap:yp[indices]=1-yp[indices]
            grouped_null.append(auc(yp,s))
        metrics[method]['source_scaffold_block_swap_p_one_sided']=float(np.mean(np.array(grouped_null)>=observed-1e-12))
        metrics[method]['source_scaffold_block_swaps']=len(grouped_null)
    summary={'protocol':protocol,'metrics':metrics,'matching':balance,'completed_jobs':sum(r['status']=='complete' for r in raw),
             'scheduled_jobs':len(key)*2,'document_count':len({d for r in rows for d in r['documents']}),
             'source_scaffold_pair_component_count':len(groups),'source_scaffold_pair_components':[[rows[i]['id'] for i in g] for g in groups],
             'limits':['Previously examined dataset; only structures new to this docking challenge, not an external independent validation.',
                       'No original-publication protocol verification; no-explicit-mutant does not prove wild type.',
                       'Only one receptor and one search seed in new matched test; no Vinardo native pose-recovery control in this comparison.',
                       'Bootstrap treats matched pairs as sampling units; shared source documents/scaffolds limit independence.',
                       'Degenerate bootstrap intervals in perfectly separated small samples do not establish certainty.',
                       'Exploratory multiple-method comparisons, no multiple-testing adjustment; no production promotion or candidate efficacy/safety claims.']}
    (out/'metrics.json').write_text(json.dumps(summary,indent=2));(out/'control_audit.json').write_text(json.dumps(review,indent=2))
    with (out/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    # Explicitly show the confounding in the original, unmatched six-control set.
    old=list(csv.DictReader((rescoring/'comparison.csv').open()))
    old_checks={}
    for pdb in ['1PY5','1VJY','2WOT']:
        controls=[r for r in old if r['pdb']==pdb and r['kind']!='candidate']
        yy=[int(r['kind']=='active_control') for r in controls]
        old_checks[pdb]={k:auc(yy,[-float(r[k]) for r in controls]) for k in ['heavy_atoms','mw']}
    (out/'original_size_confounding.json').write_text(json.dumps({'smaller_is_ranked_first_AUC':old_checks,'interpretation':'Size alone perfectly separates the original labels; perfect score-per-atom AUC is not evidence of target-specific information.'},indent=2))
    lines=['# ALK5: eşlenmiş kontrol ve alternatif skor karşılaştırması','',
           '## Ne yaptık?',
           'AIPHARMA molekül üretir, docking ve ADMET hesaplarını bir araya getirir. Bu test, hedefe göre sıralamanın bilinen aktifleri zayıflardan ayırıp ayırmadığını inceler; yeni adayların etkinlik veya güvenliliğini ölçmez.','',
           'Önceki üç aktif/üç zayıfta ağır atom başına Vina AUC üç reseptörde de 1.000 çıktı. Ancak yalnızca küçük molekülü öne almak da AUC=1.000 verdi. Bu yüzden bu sonuç bir düzeltme başarısı sayılmadı. Aynı pozlarda Vinardo AUC: 0.333 / 0.000 / 0.000.','',
           f"Sonra skorları hesaplamadan önce {len(pairs)} yeni aktif–zayıf çiftini kilitledik. Çift başına MW farkı ≤30 g/mol, logP farkı ≤0.75 ve ağır atom farkı ≤2. Önceki dokuz molekül, eski holdout iskeletleri, sansürlü ve stereokimyası belirsiz kayıtlar dışlandı. {len(key)} molekülde, 1PY5 üzerinde iki fonksiyonla {len(raw)} docking işi yapıldı.",'',
           '## Yeni eşlenmiş test sonuçları','',
           '| Yöntem | ROC-AUC | Çift içinde doğru sıralama | Çift-bootstrap %95 |',
           '|---|---:|---:|---|']
    names={'vina':'Vina 1.2.7','vinardo':'Vinardo','vina_per_heavy_atom':'Vina / ağır atom','vinardo_per_heavy_atom':'Vinardo / ağır atom','smaller_heavy_atoms':'Yalnızca daha az ağır atom','smaller_MW':'Yalnızca daha düşük MW','logP':'Yalnızca daha yüksek logP'}
    for method in methods:
        m=metrics[method]
        if m['status']!='complete':lines.append(f'| {names[method]} | eksik | — | — |');continue
        lo,hi=m['paired_bootstrap_95'];lines.append(f"| {names[method]} | {m['auc']:.3f} | {m['within_pair_correct_fraction']:.3f} | {lo:.3f}–{hi:.3f} |")
    lines+=['','AUC doğruluk yüzdesi değildir; 0.5 rastgele ayırma referansıdır. Aynı moleküller üzerinde birden çok yöntem incelendi; aralarından iyi görüneni seçmek bağımsız doğrulama değildir.',
            '',f"Kontroller {summary['document_count']} kaynak belgeye dayanır; ortak belge veya iskeletle bağlanan çiftler yalnızca {len(groups)} bağımsız bağlantı bileşeni oluşturur. Çift-bootstrap aralıkları bu bağımlılığı kapsamaz. Kaynak/iskelet blokları birlikte çevrildiğinde Vina için keşifsel tek yönlü p={metrics['vina']['source_scaffold_block_swap_p_one_sided']:.3f} olur; güçlü istatistiksel doğrulama iddiası yapılamaz. Bu test tek reseptör/tek tohumla yapılmıştır; Vinardo için yeni native poz kontrolü yapılmamıştır.",'',
            '## Karar','',
            'Bu küçük, keşifsel çalışma üzerinden üretim sıralamasını değiştirmedik. Boyut normalizasyonu tek başına çözüm kabul edilmedi. Yeni adaylara etkili/güvenli etiketi verilmedi. Eşlenmiş test, hangi yaklaşımın genişletilmeye değer olduğunu değerlendirmek içindir.','',
            'Skorlar comparison.csv içinde sıralama yönüne çevrilmiştir: daha yüksek daha önce gelir; poz dosyaları ve raw job sonuçlarında orijinal enerji işaretleri korunur.','',
            '## Tekrar üretim ve sınırlar','',
            'protocol.json skor öncesi kuralları ve girdi SHA256 özetlerini; label_key.json deney kimliklerini; control_audit.json assay açıklamalarını; jobs/ hazırlanan ligandları, pozları ve hataları; metrics.json tüm metrikleri içerir. Başarısız işler varsa saklanır ve analizden sessizce çıkarılmaz.','',
            'Bu çalışma yeni dış veri, laboratuvar testi veya bağımsız kör doğrulama değildir. Native pozun geri bulunması ayrı bir ölçüttür; ayırt ediciliği kanıtlamaz.','',
            'Vinardo kullanımı: https://autodock-vina.readthedocs.io/en/stable/docking_basic.html',
            'Sabit poz skoru: https://autodock-vina.readthedocs.io/en/stable/vina.html']
    (out/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    verification={'passed':True,'AUC_verified_against_sklearn':len(methods),'matching_all_calipers_passed':True,
                  'endpoint_review_count':len(review),'locked_input_hashes_passed':True,'completed_jobs':summary['completed_jobs'],
                  'script_sha256':sha(Path(__file__))}
    (out/'verification.json').write_text(json.dumps(verification,indent=2))
    print(json.dumps({k:round(v['auc'],3) for k,v in metrics.items() if v['status']=='complete'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--audit',type=Path,required=True);p.add_argument('--rescoring',type=Path,required=True);a=p.parse_args();main(a.output,a.audit,a.rescoring)
