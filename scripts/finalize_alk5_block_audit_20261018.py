"""Independently verify graph partitions and create a component overview."""
import hashlib
import json
from pathlib import Path

base=Path(__file__).resolve().parents[1]/'reports/alk5_block_audit_20261018'
checked=0
for path in sorted(base.glob('*_graph.json')):
    data=json.loads(path.read_text());n=data['molecule_count'];adj=[set() for _ in range(n)]
    for edge in data['edges']:
        adj[edge['a']].add(edge['b']);adj[edge['b']].add(edge['a'])
    remaining=set(range(n));found=[]
    while remaining:
        start=min(remaining);todo=[start];seen={start};remaining.remove(start)
        while todo:
            node=todo.pop()
            for other in adj[node]:
                if other not in seen:seen.add(other);remaining.remove(other);todo.append(other)
        found.append(frozenset(seen))
    assert set(found)=={frozenset(c['indices']) for c in data['components']}
    checked+=1
data=json.loads((base/'history_and_stereo_eligible_source_exact_scaffold_graph.json').read_text())
blocks=data['components'];height=130+90*((len(blocks)+2)//3)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1020" height="{height}" viewBox="0 0 1020 {height}" role="img" aria-labelledby="title desc">',
     '<title id="title">ALK5 source and exact-scaffold components</title>',
     '<desc id="desc">Eligible pool: 318 active and 9 weak molecules. 21 components; only 2 contain both classes.</desc>',
     f'<rect width="1020" height="{height}" fill="#f8fafc"/>',
     '<g font-family="Arial, sans-serif" fill="#0f172a">',
     '<text x="24" y="32" font-size="22">ALK5: source / exact-scaffold components</text>',
     '<text x="24" y="59" font-size="15">Eligible pool: 318 active + 9 weak; edges = shared document OR identical scaffold.</text>',
     '<text x="24" y="82" font-size="14">Teal: both classes | Amber: weak only | Grey: active only. Boxes are disconnected components.</text>']
for i,b in enumerate(blocks):
    x=24+(i%3)*332;y=100+(i//3)*90
    color='#ccfbf1' if b['both_classes'] else '#fef3c7' if b['weak'] else '#e2e8f0'
    svg += [f'<rect x="{x}" y="{y}" width="310" height="72" rx="7" fill="{color}" stroke="#64748b"/>',
            f'<text x="{x+12}" y="{y+25}" font-size="16">Component {b["block"]}</text>',
            f'<text x="{x+12}" y="{y+49}" font-size="14">Active {b["active"]} | Weak {b["weak"]} | Documents {len(b["documents"])}</text>']
svg+=['</g></svg>'];(base/'components.svg').write_text('\n'.join(svg),encoding='utf-8')
readme=base/'README.md';s=readme.read_text(encoding='utf-8')
s+='''
## Uygulama kararı

Önceki holdout iskeletleri ve belirsiz stereokimya dışlandıktan sonra uygun havuzda 318 aktif ve 9 zayıf var. Son iki docking seti çıkarılınca 309 aktif, **0 zayıf** kalıyor. Bu nedenle aynı havuzdan yeni bir kontrol testi seçip bağımsız doğrulama diye sunmuyoruz; bu aşamada yeni docking başlatılmadı.

Tüm uygun havuzda 21 kaynak/aynı-iskelet bileşeni var, fakat bunların yalnızca 2'sinde iki sınıf birlikte bulunuyor. Bu sayı, beş bağımsız aktif–zayıf bloğunun hazır bulunduğunu göstermiyor; alt küme tasarımının kesin matematiksel üst sınırı olarak da yorumlanmamalıdır.

Sonraki veri toplama hedefi daha fazla aktif molekül değil, farklı kaynak ve kimyasal serilerden **deneysel enzim zayıf kontrolleri** ve bunlarla karşılaştırılabilir aktiflerdir. Yeni kayıtların önceki havuzla molekül, kaynak, iskelet ve assay tekrarları denetlenmelidir. Aynı yayın ChEMBL ve BindingDB'de bulunuyorsa iki bağımsız kaynak sayılmamalıdır. Protokol koşulları ve ara aktivite kayıtları ayrıca korunmalıdır.

İskeletler molekül SMILES'larından yeniden hesaplandı. Tarihsel `ACYCLIC` yer tutucusu geçerli bir iskelet SMILES'ı sayılmadı; bütün halkasız moleküller sırf halkasız oldukları için birleştirilmedi. Orijinal alan `source_scaffold` olarak korundu.

![Uygun havuzun bağlantı bileşenleri](components.svg)
'''
readme.write_text(s,encoding='utf-8')
verification=json.loads((base/'verification.json').read_text())
verification.update({'independent_adjacency_traversal_verified_graphs':checked,'new_docking_jobs':0})
(base/'verification.json').write_text(json.dumps(verification,indent=2))
manifest={str(p.relative_to(base)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(base.rglob('*')) if p.is_file() and p.name!='MANIFEST_SHA256.json'}
(base/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2))
assert all(hashlib.sha256((base/p).read_bytes()).hexdigest()==h for p,h in manifest.items())
print(json.dumps({'graphs_verified':checked,'files_hashed':len(manifest),'new_docking_jobs':0}))
