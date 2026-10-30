"""Freeze and recheck declaration, method, comparator and preparation bytes; no Engine import."""
import ast
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
SOURCE = ROOT / 'build/issue-929/main'
OUT = HERE / 'results/preparation.json'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def check():
    record = json.loads(OUT.read_text())
    for name, expected in record['sha256'].items():
        assert digest(Path(name)) == expected, f'changed image/method: {name}'
    return record


def main():
    if sys.argv[1] == 'check':
        check()
        print('Preparation bytes unchanged')
        return
    assert not OUT.exists()
    image = json.loads((HERE / 'results/image.json').read_text())
    assert image['commit'] == (ROOT / 'build/issue-929/main-sha.txt').read_text().strip()
    files = [SOURCE / 'README.md', SOURCE / 'uv.lock', Path(image['loader']),
             HERE / 'results/image.json', ROOT / 'build/issue-929/setup.log']
    for parent in (ROOT / 'analyses/2026-audit-timing-campaign/scripts', HERE / 'scripts', HERE / 'inputs',
                   SOURCE / 'data', SOURCE / 'analyses/2026-competitor-cross-check',
                   SOURCE / 'analyses/2026-head-to-head-comparisons', SOURCE / 'analyses/2026-density-root-benchmark',
                   ROOT / 'analyses/2026-competitor-cross-check', ROOT / 'analyses/2026-density-root-benchmark'):
        files += [p for p in parent.rglob('*') if p.is_file() and p.suffix in ('.py', '.jl', '.json', '.csv', '.toml')]
    build = Path(image['build_dir'])
    files += [build / n for n in image['libraries_sha256']]
    files += [p for p in (SOURCE / 'engine/src/epcsaft').rglob('*.py')]
    files += [build / 'build.ninja', build / 'meson-info/intro-install_plan.json', build / 'meson-info/intro-buildoptions.json']
    cache = Path.home() / '.cache/epcsaft-comparators/222'
    files += [cache / 'julia-1.13.1/bin/julia', cache / 'feos/bin/python']
    files += list((cache / 'feos/lib').glob('python*/site-packages/feos/*.so'))
    for parent in (cache / 'depot/packages/Clapeyron',):
        files += [p for p in parent.rglob('*') if p.is_file() and p.suffix in ('.jl', '.toml')]
    files += list((cache / 'julia-1.13.1/lib').glob('libjulia.so*'))
    feeds_path = SOURCE / 'analyses/reproduction/2026-khudaida/figures/figure-02/source/source-points.csv'
    files.append(feeds_path)
    # The fixed physical snapshots must equal the captured source's declarations.
    declared = json.loads((HERE / 'inputs/declared-cases.json').read_text())
    current = json.loads((SOURCE / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    assert all(v == current['models'][k] for k, v in declared['models'].items())
    assert all(r in current['requests'] for r in declared['requests'])
    hh = SOURCE / 'analyses/2026-head-to-head-comparisons/scripts/h2h_engine.py'
    tree = ast.parse(hh.read_text())
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in ('coefficient', 'khudaida_document', 'khudaida_feeds')]
    namespace = dict(json=json, csv=csv, KH=SOURCE / 'data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json', FEEDS=feeds_path, KH_IDS=('water', 'ethanol', 'isobutanol', 'sodium-cation', 'chloride-anion'))
    exec(compile(tree, str(hh), 'exec'), namespace)
    assert namespace['khudaida_document']('constant-permittivity') == json.loads((HERE / 'inputs/electrolyte-model.json').read_text())
    feeds = {int(r.pop('point')): [float(v) for v in r.values()] for r in csv.DictReader((HERE / 'inputs/electrolyte-feeds.csv').open())}
    assert namespace['khudaida_feeds']() == feeds
    record = dict(implementation_base='da73ef241844efb01fba3249e369568dd1e99608', main_sha=image['commit'], build_jobs=6,
        compiler=subprocess.check_output(['c++', '--version'], text=True),
        sha256={str(p): digest(p) for p in sorted(set(files))},
        readme=(SOURCE / 'README.md').read_text(), build_options=json.loads((build / 'meson-info/intro-buildoptions.json').read_text()))
    OUT.write_text(json.dumps(record, indent=2) + '\n')
    check()


if __name__ == '__main__':
    main()
