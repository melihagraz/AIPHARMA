"""Audit source/scaffold connected components without reading docking scores."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def graph(rows, threshold=None):
    parent=list(range(len(rows)));edges=[]
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    fps={r['scaffold']:gen.GetFingerprint(Chem.MolFromSmiles(r['scaffold'])) for r in rows if r['scaffold']}
    for i,a in enumerate(rows):
        for j in range(i):
            b=rows[j];docs=sorted(set(a['documents'])&set(b['documents']))
            same=bool(a['scaffold'] and a['scaffold']==b['scaffold'])
            similarity=None
            if threshold is not None and a['scaffold'] and b['scaffold']:
                similarity=DataStructs.TanimotoSimilarity(fps[a['scaffold']],fps[b['scaffold']])
            if docs or same or (similarity is not None and similarity>=threshold):
                parent[root(i)]=root(j)
                edges.append({'a':i,'b':j,'shared_documents':docs,'identical_scaffold':same,'scaffold_tanimoto':similarity})
    components={}
    for i in range(len(rows)): components.setdefault(root(i),[]).append(i)
    blocks=[]
    for index,indices in enumerate(sorted(components.values(),key=lambda v:(-sum(rows[i]['label']==0 for i in v),-len(v),v[0])),1):
        members=[rows[i] for i in indices];counts=Counter(r['label'] for r in members)
        blocks.append({'block':index,'molecules':len(indices),'active':counts[1],'weak':counts[0],
                       'both_classes':bool(counts[1] and counts[0]),'documents':sorted({d for r in members for d in r['documents']}),
                       'scaffold_count':len({r['scaffold'] for r in members}),'indices':indices})
    mixed=sum(b['both_classes'] for b in blocks)
    return {'molecule_count':len(rows),'active':sum(r['label']==1 for r in rows),'weak':sum(r['label']==0 for r in rows),
            'component_count':len(blocks),'mixed_class_components':mixed,'weak_containing_components':sum(b['weak']>0 for b in blocks),
            'components':blocks,'edges':edges}


def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    if (a.output/'summary.json').exists():raise RuntimeError('Completed audit is immutable')
    data=json.loads((a.audit/'datasets_all.json').read_text())['enzyme_catalytic__no_explicit_mutant']
    # Recompute actual scaffolds: the historical ACYCLIC sentinel is not a molecule.
    data=[{**r,'source_scaffold':r['scaffold'],'scaffold':MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(r['smiles']),includeChirality=False)} for r in data]
    previous=json.loads((a.ensemble/'locked_items.json').read_text())
    matched=json.loads((a.matched/'label_key.json').read_text())
    exposed={r['smiles'] for r in previous+matched}
    exact=[r for r in data if r['evidence']=='exact']
    eligible=[];exclusions=[]
    for r in exact:
        reasons=[]
        if r['old_holdout_scaffold']:reasons.append('old_holdout_scaffold')
        if any(str(x.specified)=='Unspecified' for x in Chem.FindPotentialStereo(Chem.MolFromSmiles(r['smiles']))):reasons.append('unspecified_stereochemistry')
        if not reasons:eligible.append(r)
        else:exclusions.append({'smiles':r['smiles'],'label':r['label'],'reasons':reasons})
    unseen=[r for r in eligible if r['smiles'] not in exposed]
    pools={'all_exact':exact,'history_and_stereo_eligible':eligible,'not_in_last_two_docking_sets':unseen}
    protocol={'scope':'Source/scaffold dependency audit, no docking score access or threshold tuning',
              'primary_edges':'shared source document OR identical nonempty Murcko scaffold; transitive connected components',
              'sensitivity_edges':'primary edges OR scaffold Morgan radius2/2048 Tanimoto >=0.4; exploratory cutoff, not a universal independence criterion',
              'empty_scaffolds':'Not joined just because both scaffolds are empty',
              'pools':list(pools),'five_blocks':'Design target only; five independent blocks do not guarantee adequate power or assay comparability',
              'minimum_p':'2^-k only for a one-sided sign/label-swap test with k independent binary blocks; not a general AUC significance formula',
              'exposure':'Unseen means absent from the last two named docking sets and old holdout scaffolds, not never examined or truly external',
              'input_hashes':{str(p):sha(p) for p in [a.audit/'datasets_all.json',a.ensemble/'locked_items.json',a.matched/'label_key.json']}}
    (a.output/'protocol.json').write_text(json.dumps(protocol,indent=2))
    results={}
    for name,rows in pools.items():
        results[name]={}
        (a.output/f'{name}_nodes.json').write_text(json.dumps(rows,indent=2))
        for mode,threshold in [('source_exact_scaffold',None),('source_scaffold_similarity_040',.4)]:
            result=graph(rows,threshold);results[name][mode]={k:v for k,v in result.items() if k not in ['components','edges']}
            (a.output/f'{name}_{mode}_graph.json').write_text(json.dumps(result,indent=2))
            with (a.output/f'{name}_{mode}_blocks.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=['block','molecules','active','weak','both_classes','documents','scaffold_count']);w.writeheader()
                for b in result['components']:w.writerow({k:v for k,v in b.items() if k!='indices'})
            assert sum(b['molecules'] for b in result['components'])==len(rows)
            # Every recorded relation must stay inside its assigned component.
            assignment={i:b['block'] for b in result['components'] for i in b['indices']}
            assert all(assignment[e['a']]==assignment[e['b']] for e in result['edges'])
    (a.output/'exclusions.json').write_text(json.dumps(exclusions,indent=2))
    (a.output/'summary.json').write_text(json.dumps(results,indent=2))
    lines=['# ALK5: bağımsız kontrol blokları denetimi','',
           'Amaç: yeni docking çalıştırmadan önce mevcut katalitik enzim havuzunda kaynak ve iskelet bakımından ayrık, hem aktif hem zayıf içeren grupların bulunup bulunmadığını belirlemek. Docking skorları okunmadı.','',
           '## Sonuçlar','',
           '| Havuz | Aktif | Zayıf | Tüm bileşenler | Zayıf içeren | Her iki sınıfı içeren |',
           '|---|---:|---:|---:|---:|---:|']
    titles={'all_exact':'Tüm kesin ölçümler','history_and_stereo_eligible':'Eski holdout iskeletleri ve belirsiz stereo dışarıda','not_in_last_two_docking_sets':'Son iki docking kümesinde bulunmayan uygun yapılar'}
    for name in pools:
        r=results[name]['source_exact_scaffold'];lines.append(f"| {titles[name]} | {r['active']} | {r['weak']} | {r['component_count']} | {r['weak_containing_components']} | {r['mixed_class_components']} |")
    lines+=['','Bir bileşen, ortak kaynak VEYA aynı Murcko iskeleti üzerinden doğrudan veya dolaylı bağlanan moleküllerdir. Molekül sayısı, bağımsız karşılaştırma sayısı değildir. Büyük bir bileşeni çizimde parçalamak, aradaki bağımlılığı kaldırmaz.','',
            '## Daha sıkı benzerlik duyarlılığı','',
            'Scaffold Morgan/ECFP4 benzerliği ≥0.4 olanlar da bağlandığında sonuçlar:']
    for name in pools:
        r=results[name]['source_scaffold_similarity_040'];lines.append(f"- {titles[name]}: {r['component_count']} bileşen; {r['mixed_class_components']} iki sınıflı bileşen.")
    lines+=['','0.4 sınırı bir duyarlılık incelemesidir; biyolojik bağımsızlığı kanıtlayan evrensel bir eşik değildir. Molekül ECFP benzerliği ile Murcko iskeleti benzerliği aynı kavram değildir.','',
            '## Tasarım kararı','',
            'Yeni docking için en az beş kaynak/iskelet ayrık ve iki sınıfı içeren blok hedefleniyorsa, bu sayının mevcut havuzdan sağlanıp sağlanmadığı yukarıdaki tablodan değerlendirilmelidir. Beş blok tek başına yeterli güç veya yayınlanabilir doğrulama garantisi değildir.','',
            'Yalnızca aktif içeren bir kaynak, zayıf içeren başka bir kaynağa eşlenince ikisi tek karşılaştırma bloğu olur. Böyle eşlemeler ek bağımsız kaynak yaratmaz. Aynı grafiği kullanarak alt küme seçmek bazı köprüleri kaldırabilir; bunun bağımsızlık sağladığına karar vermeden önce ortak kaynak/seri ve assay farkları ayrıca incelenmelidir.','',
            'Bu denetim, mevcut verilerden elde edilebilecek en fazla alt küme sayısı için optimizasyon veya kesin imkânsızlık ispatı değildir. Raporlanan sayılar tüm havuzun bağlantı bileşenleridir.','',
            'Son iki docking kümesinde olmayan yapı, yeni dış veri demek değildir. Tüm havuz önceki veri kürasyonu/model değerlendirmesinde görülmüştür.','',
            'Yeni kontrol gerekirse özgün yayın/protokolü doğrulanmış, aynı enzim endpointinde ölçülmüş aktif ve zayıf kayıtlar aranmalı; hücresel ölçümler veya yapay çeldiriciler deneyle doğrulanmış enzim negatifleri gibi eklenmemelidir.','',
            '## Tek yönlü işaret testinin teorik alt sınırı','', '| Bağımsız ikili blok (k) | En küçük p = 2^-k |','|---|---:|']
    for k in range(1,9):lines.append(f'| {k} | {2**(-k):.6f} |')
    lines+=['','Bu tablo, tüm blokların aynı yönde sonuç verdiği bağımsız ikili test içindir; önceki AUC veya bootstrap sonuçlarına otomatik p değeri sağlamaz.','',
            '## Dosyalar','',
            'Her havuz için nodes.json molekülleri; graph.json tüm kenar gerekçelerini ve bileşen üyelerini; blocks.csv kaynak ve sınıf sayılarını içerir. protocol.json kuralları ve girdi özetlerini; verification.json denetimleri; exclusions.json dışlamaları saklar.']
    (a.output/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (a.output/'verification.json').write_text(json.dumps({'passed':True,'no_docking_scores_read':True,'all_nodes_assigned_once':True,'all_edges_within_components':True,'script_sha256':sha(Path(__file__))},indent=2))
    print(json.dumps(results,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--audit',type=Path,required=True);p.add_argument('--ensemble',type=Path,required=True);p.add_argument('--matched',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
