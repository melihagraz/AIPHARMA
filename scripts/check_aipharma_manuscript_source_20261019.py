"""Static source checks only; never substitutes for LaTeX compilation/layout review."""
import collections,hashlib,json,re
from pathlib import Path
p=Path('reports/latex/PharmaProject_arXiv_manuscript_UPDATED.tex'); text=p.read_text(encoding='utf-8')
assets=re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',text)
missing=[a for a in assets if not (p.parent/a).is_file()]
labels=re.findall(r'\\label\{([^}]+)\}',text); refs=re.findall(r'\\(?:ref|eqref|pageref)\{([^}]+)\}',text)
duplicates=[k for k,v in collections.Counter(labels).items() if v>1]
undefined=sorted(set(refs)-set(labels)); stack=[]; env_errors=[]
for match in re.finditer(r'\\(begin|end)\{([^}]+)\}',text):
    kind,env=match.groups()
    if kind=='begin': stack.append(env)
    elif not stack or stack[-1]!=env: env_errors.append(dict(position=match.start(),environment=env))
    else: stack.pop()
sections=re.findall(r'\\(?:section|subsection|subsubsection)\{(\d+(?:\.\d+)*)[ .]',text)
dup_sections=[k for k,v in collections.Counter(sections).items() if v>1]
summary=dict(source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),figures=len(assets),missing_assets=missing,
    duplicate_labels=duplicates,undefined_references=undefined,environment_errors=env_errors,unclosed_environments=stack,
    duplicate_manual_section_numbers=dup_sections,empty_author_field=bool(re.search(r'\\author\{\s*\}',text)),
    static_checks_passed=not any([missing,duplicates,undefined,env_errors,stack,dup_sections]),
    compile_attempts=2,compile_success=False,compile_blocker='Windows sandbox CreateProcessWithLogonW failed:267 before compiler startup',
    pdf_visual_review_performed=False)
Path('reports/manuscript_static_check_20261019.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
assert summary['static_checks_passed']
