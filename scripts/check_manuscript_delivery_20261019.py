"""Static source/asset checks; explicitly not a substitute for TeX compilation."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
tex = root/'reports/latex/PharmaProject_arXiv_manuscript_UPDATED.tex'
text = tex.read_text(encoding='utf-8')
before = tex.with_name('PharmaProject_arXiv_manuscript_UPDATED.before_delivery_20261019.tex').read_text(encoding='utf-8')
images = re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', text)
images = [s for s in images if not s.startswith('#')]
assert all((tex.parent/s).is_file() for s in images), 'Missing image'
original_images = [s for s in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', before) if not s.startswith('#')]
assert set(original_images) <= set(images), 'An original illustration was removed'
labels = re.findall(r'\\label\{([^}]+)\}', text)
duplicates = [label for label, n in Counter(labels).items() if n > 1]
assert not duplicates, duplicates
references = re.findall(r'\\(?:ref|eqref|pageref)\{([^}]+)\}', text)
assert set(references) <= set(labels), 'Missing referenced label'
stack = []
for action, name in re.findall(r'\\(begin|end)\{([^}]+)\}', text):
    if action == 'begin':
        stack.append(name)
    else:
        assert stack and stack.pop() == name, f'Environment mismatch: {name}'
assert not stack
assert len(text) > len(before), 'Unexpected manuscript shortening'
result = dict(status='static_checks_passed', characters_before=len(before), characters_after=len(text), original_image_references=len(original_images), current_image_references=len(images), all_original_illustrations_retained=True, all_image_files_present=True, unique_labels=True, references_resolved=True, environments_balanced=True, source_sha256=hashlib.sha256(tex.read_bytes()).hexdigest(), pdf_compilation_verified=False, compiler_attempts=3, compiler_error='Could not prepare LaTeX preview: windows sandbox CreateProcessWithLogonW failed: 267')
out = root/'reports/alk5_delivery_20261019'
out.mkdir(exist_ok=True)
(out/'manuscript_source_check.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
