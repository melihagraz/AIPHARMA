"""Collect primary synthesis/table text and linked drawings, without claiming structure validation."""
import html,json,re,sys
from pathlib import Path
p=Path(sys.argv[1]); root=p/'additional_local_series/primary_audit'; out=root/'structure_evidence'; out.mkdir(exist_ok=True)
matches=json.loads((root/'query_value_matches.json').read_text()); output=[]
for r in matches:
    if not r['all_alias_values_match']: continue
    evidence=[]
    for check in r['checks']:
        patent=check['patent']; text=(root/(patent+'.html')).read_text(encoding='utf-8-sig')
        for candidate in check['matching_candidates']:
            number=candidate['compound_number']; escaped=re.escape(number)
            headings=list(re.finditer(r'<(?:div|heading)\b[^>]*>\s*(?:Example\s+\d+\s+)?(?:Synthetic Route of Compound|Synthesis of Compound)\s+'+escaped+r'\s*</(?:div|heading)>',text,re.I))
            blocks=[]
            for match in headings:
                following=text[match.end():]; end=re.search(r'<(?:div|heading)\b[^>]*>\s*(?:Example\s+\d+\s+)?(?:Synthetic Route of Compound|Synthesis of Compound)\s+\S+',following,re.I)
                block=following[:end.start() if end else min(15000,len(following))]
                plain=' '.join(html.unescape(re.sub('<[^>]+>',' ',block)).split())
                urls=re.findall(r'<img\b[^>]*src="([^"]+)"',block)
                blocks.append(dict(heading_offset=match.start(),text=plain,image_urls=urls))
            # Some patents give selected structures in the characterization table only.
            data=json.loads((root/(patent+'_tables.json')).read_text()); table_rows=[]
            if patent=='US11691981B2':
                for row in data['all_tables'][1]:
                    cells=[c.strip() for c in row if c.strip()]
                    if cells and cells[0]==number: table_rows.append(row)
            table_images=[]
            if patent=='US11691981B2':
                for tr in re.findall(r'<tr\b[^>]*>.*?</tr>',text,re.S):
                    cells=re.findall(r'<td\b[^>]*>(.*?)</td>',tr,re.S)
                    if cells and re.sub('<[^>]+>','',cells[0]).strip()==number:
                        table_images.extend(re.findall(r'<img\b[^>]*src="([^"]+)"',tr))
            evidence.append(dict(patent=patent,compound_number=number,source_label=check['bindingdb_label'],synthesis_blocks=blocks,characterization_table_rows=table_rows,characterization_image_urls=table_images))
    output.append(dict(assay_id=r['assay_id'],normalized_parent=r['normalized_parent'],primary_evidence=evidence,
                       independently_reconstructed_structure=False))
(out/'primary_structure_evidence.json').write_text(json.dumps(output,indent=2))
summary=dict(molecules=len(output),with_synthesis_text=sum(any(e['synthesis_blocks'] for e in r['primary_evidence']) for r in output),
 with_characterization_rows=sum(any(e['characterization_table_rows'] for e in r['primary_evidence']) for r in output),
 independently_reconstructed_structures=0,conclusion='Evidence locations only; molecular mass or a linked drawing is not a full atom/bond identity check')
(out/'summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
