"""Archive target-specific source/assay records without treating target ID as endpoint evidence."""
import concurrent.futures,hashlib,html,json,re,sys,urllib.request
from pathlib import Path
p=Path(sys.argv[1]); folder=p/'bindingdb'; folder.mkdir(exist_ok=True)
entries=[9244,9757,2984,11765,11969,11968,8608,10870,11335,10175,10830]
(p/'fetch_protocol.json').write_text(json.dumps(dict(entries=entries,purpose='Endpoint/source curation only; no model fitting or scores',rules='Retain raw target and assay labels; explicit mutants excluded later; do not infer enzyme assay from target name',holdout='Allocate entire verified source-family/scaffold components before any new model run; no claim of analyst blinding',max_hits_per_entry=3000),indent=2))
base='https://www.bindingdb.org'
def fetch(url,dest):
    if not dest.exists():
        request=urllib.request.Request(url,headers={'User-Agent':'AIPHARMA-source-audit/1.0'})
        with urllib.request.urlopen(request,timeout=45) as response: raw=response.read()
        dest.write_bytes(raw)
    return dest.read_text(encoding='utf-8-sig')
def stripped(text): return ' '.join(html.unescape(re.sub('<[^>]+>',' ',text)).split())
def process(entry):
    try:
        prefix=f'{base}/rwd/jsp/dbsearch/PrimarySearch_ki.jsp?entryid={entry}&tag=entry'
        first=fetch(prefix,folder/f'entry{entry}_page0.html')
        match=re.search(r'Found\s*<span class="red">(\d+)</span>',first); assert match
        count=int(match[1]); assert count<=3000
        pages=[first]+[fetch(prefix+f'&startPg={start}&Increment=50&submit=Search',folder/f'entry{entry}_page{start}.html') for start in range(50,count,50)]
        records=[]; assays={}
        for page in pages:
            for block in page.split('<div><div><span class="header">Target</span>')[1:]:
                target=re.search(r'<a class="big"[^>]*>(.*?)</a>',block,re.S)
                if not target: continue
                name=stripped(target[1])
                if not re.search(r'TGF.*(?:receptor.*(?:type-?1|type I)|beta.*receptor.*1)',name,re.I): continue
                mid=re.search(r'monomerid=(\d+)',block)
                smi=re.search(r"setClipboard\('([^']+)'\).*?Copy&nbsp;SMILES",block)
                affinity=re.search(r'<span class="big">(IC50|Ki|Kd|EC50):&nbsp;([^<]+)</span>',block)
                details=re.search(r'href="([^"]*Summary_ki.jsp\?[^"]+)"',block)
                assay_links=re.findall(r'href="([^"]*assay.jsp\?[^"]+)"',block)
                assay_snip=re.search(r'Assay Description:</(?:a|span)><span>(.*?)</span>',block,re.S)
                assert mid and smi and affinity and details
                record=dict(entry_id=entry,target=name,monomer_id=mid[1],smiles=html.unescape(smi[1]),measurement_type=affinity[1],reported_value=stripped(affinity[2]),details_url=base+html.unescape(details[1]),assay_snippet=stripped(assay_snip[1]) if assay_snip else None,assay_links=[base+html.unescape(a) for a in assay_links])
                records.append(record)
                for link in record['assay_links']:
                    aid=re.search(r'assayid=(\d+)',link)[1]
                    text=fetch(link,folder/f'entry{entry}_assay{aid}.html')
                    assays[aid]=dict(url=link,text=stripped(text),sha256=hashlib.sha256(text.encode()).hexdigest())
        return dict(entry_id=entry,all_target_hits=count,alk5_records=records,assays=assays,status='complete')
    except Exception as exc: return dict(entry_id=entry,status='failed',error=str(exc))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    results=[]
    for r in pool.map(process,entries):
        results.append(r); (p/'entry_archive.json').write_text(json.dumps(results,indent=2))
        print(json.dumps(dict(entry=r['entry_id'],status=r['status'],ALK5_records=len(r.get('alk5_records',[])),assays=list(r.get('assays',{})),error=r.get('error'))),flush=True)
