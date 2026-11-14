"""Retrieve the publisher's publicly linked supplementary-table workbook."""
import concurrent.futures,hashlib,json,sys,urllib.request,zipfile
from pathlib import Path
p=Path(sys.argv[1]); p.mkdir(exist_ok=True)
base='https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41587-026-03090-8/MediaObjects/'
def fetch(i):
    name=f'41587_2026_3090_MOESM{i}_ESM.xlsx'; url=base+name; dest=p/name
    try:
        if not dest.exists():
            request=urllib.request.Request(url,headers={'User-Agent':'AIPHARMA-public-research/1.0'})
            with urllib.request.urlopen(request,timeout=40) as response: data=response.read()
            assert data[:2]==b'PK', 'Not an Excel package'
            dest.write_bytes(data)
        with zipfile.ZipFile(dest) as z: assert 'xl/workbook.xml' in z.namelist()
        return dict(url=url,file=name,status='complete',sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),bytes=dest.stat().st_size)
    except Exception as exc: return dict(url=url,status='unavailable',error=str(exc))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool: results=list(pool.map(fetch,range(1,7)))
(p/'supplement_fetch.json').write_text(json.dumps(results,indent=2)); print(json.dumps(results,indent=2),flush=True)
