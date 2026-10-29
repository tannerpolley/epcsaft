"""Transport frozen input bytes and the complete /3 parameter declaration."""
import csv
import hashlib
import json
from pathlib import Path
import shutil
import re
from epcsaft import Parameters

ROOT = Path(__file__).resolve().parents[3]
MAIN = Path('/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT')
OUT = ROOT / 'analyses/2026-mixed-chloride-water-activity/inputs/issue-902-p12-r1'
paths = [*[(MAIN / 'build/issue-899' / group / name) for group, names in [
    ('sources-r2', ['el-guendouzi2001-mg-ca.csv', 'el-guendouzi2005-table1.csv',
                    'azougen2010-mixtures.csv', 'azougen2010-binary-Ca.csv',
                    'goldberg-high-m.csv', 'rard1981-tableIX.csv']),
    ('sources-r4', ['rard-reconstructed.csv', 'rard-reference-paths.csv',
                    'robinson-reconstructed.csv', 'robinson-reference-reconstructed.csv'])]
    for name in names],
    *sorted((ROOT / 'analyses/2026-hydrochloric-acid-chloride-basis/inputs').glob('*hcl*.csv')),
    ROOT / 'data/packets/yao-1992-licl-mgcl2-observations/1/observations/yao1992-licl-mgcl2-table1.csv',
    *[(MAIN.parent / 'Lithium_Extraction/data/reference/extraction_models/top_fecl3/wu2025' / n)
      for n in ['printed_values.csv', 'states.csv', 'sources.json', 'source_conflicts.csv']]]
paths += [ROOT / 'analyses/2026-lithium-chloride-refit/inputs/frozen-realizations.csv']
records = []
design = (MAIN / 'build/issue-902/design-p12-r1.md').read_text()
expected = dict(re.findall(r"\| `([^`]+)` \| `([a-f0-9]{64})`", design))
for path in paths:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    for key, frozen in expected.items():
        if key.endswith('/' + path.name) or key == path.name:
            assert digest == frozen, path
    records.append(dict(path=str(path), sha256=digest, retained=path.name))
    if not path.is_relative_to(ROOT):
        shutil.copyfile(path, OUT / path.name)
packet = ROOT / 'data/packets/figiel-2025-chloride-mixture-refit/3'
for path in [packet / 'packet.toml', packet / 'parameters/parameters.json', packet / 'provenance/source-origin.toml']:
    records.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
p = Parameters.from_json(packet / 'parameters/parameters.json')._document
with (OUT / 'parameters.csv').open('w') as f:
    w = csv.writer(f)
    w.writerow(['field', 'identity', 'indices', 'value', 'unit'])
    for attr in ['component_ids', 'site_component', 'site_multiplicity']:
        for v in getattr(p, attr):
            w.writerow([attr, '', '', v, ''])
    for attr in ['base_choice', 'association_choice', 'electrolyte_choice', 'ionic_dispersion_choice', 'permittivity_choice', 'polar_choice']:
        w.writerow([attr, '', '', getattr(p, attr), ''])
    for e in p.entries:
        w.writerow([e.field, e.identity, ' '.join(map(str, e.index)), e.magnitude, e.unit])
(OUT / 'input-identities.json').write_text(json.dumps(records, indent=2) + '\n')
roles = json.loads((ROOT / 'analyses/2026-lithium-chloride-refit/inputs/fit-declaration.json').read_text())['row_lists']['Yao']
with (OUT / 'yao-roles.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=['csv_row', 'locator', 'role'])
    w.writeheader()
    w.writerows(roles)
replay = ROOT / 'analyses/2026-hydrochloric-acid-chloride-basis/results/packet-readback-r4/numerical/training.csv'
# The native reader uses the retained acid owner directly.
