"""Resolve standalone evidence links and retain the two rendered figure views."""
import base64
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urlsplit

leaf = Path(__file__).resolve().parents[1]
out = leaf.parents[2] / 'build/issue-729/render'
html = out / 'notebook.html'
text = html.read_text()


class Evidence(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images, self.links = [], []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'img' and a.get('src', '').startswith('data:image/png;base64,'):
            self.images.append(a['src'].split(',', 1)[1])
        if tag == 'a' and a.get('href'):
            url = urlsplit(a['href'])
            if url.path and not url.scheme and not url.netloc:
                self.links.append(a['href'])


p = Evidence()
p.feed(text)
for link in set(p.links):
    target = (out / link).resolve() if link.startswith('../../../') else (leaf / link).resolve()
    if not target.is_relative_to(leaf) or not target.exists():
        raise ValueError(f'local evidence link does not resolve within the leaf: {link}')
    text = text.replace(f'href="{link}"', f'href="../../../analyses/reproduction/2025-perdomo/{target.relative_to(leaf)}"')
html.write_text(text)
assert len(p.images) == 2, 'expected both numerical comparison figures'
for name, data in zip(['nacl-pressure.png', 'licl-endpoints.png'], p.images):
    (out / name).write_bytes(base64.b64decode(data))
paths = [leaf / 'notebook.qmd', *sorted((leaf / 'inputs/current').glob('*')), *sorted((leaf / 'inputs/source').glob('*.csv')), *sorted((leaf / 'results/current').glob('*')), html, out / 'nacl-pressure.png', out / 'licl-endpoints.png']
(out / 'presentation-inputs.json').write_text(json.dumps({str(p.relative_to(leaf.parents[2])): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}, indent=2)+'\n')
print(f'Checked {len(set(p.links))} local evidence links; retained HTML, both figures and presentation hashes in {out}.')
