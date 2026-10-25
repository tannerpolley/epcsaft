#!/usr/bin/env python3
"""Construct current analysis inputs with the installed Engine; never solve."""
import csv
from pathlib import Path
import re
import runpy
import sys
from types import SimpleNamespace

import epcsaft
from epcsaft import Mixture, Pure, equilibrium as eq, regression as reg

ROOT = Path(__file__).resolve().parents[1]


def forbidden(*args, **kwargs):
    raise RuntimeError('declaration check attempted numerical execution')


# Each covered script defines declarations(); a script whose top level exits under
# --declare-only declares there.
DECLARED = re.compile(r'^def declarations\(\):', re.M)
# These method records require their historical wheel or application checkout.
HISTORICAL = {r['analysis']: r['reason'] for r in csv.DictReader((ROOT / 'tools/analysis-declarations-excluded.tsv').open(), delimiter='\t')}


def leaf(path):
    parts = path.relative_to(ROOT / 'analyses').parts
    return '/'.join(parts[:2] if parts[0] == 'reproduction' else parts[:1])


def check(path):
    old_path, old_argv, old_modules = sys.path[:], sys.argv[:], set(sys.modules)
    try:
        # Inputs and nested folders may import their leaf's scripts/ modules.
        scripts = ROOT / 'analyses' / leaf(path) / 'scripts' if path.is_relative_to(ROOT / 'analyses') else path.parent
        sys.path[:0] = [str(path.parent), str(scripts), str(ROOT / 'engine')]
        sys.argv = [str(path), 'initial', '--declare-only'] if path.stem == 'replay' else [str(path), '--declare-only', '--output', str(ROOT / 'build/analysis-declarations/unused.json')]
        try:
            module = SimpleNamespace(**runpy.run_path(str(path)))
            module.declarations()
        except SystemExit as error:
            if error.code != 0:
                raise
    finally:
        sys.path[:], sys.argv[:] = old_path, old_argv
        for name in set(sys.modules) - old_modules:
            if str(getattr(sys.modules[name], '__file__', '')).startswith(str(ROOT / 'analyses')):
                del sys.modules[name]


def main():
    if 'site-packages' not in epcsaft._core.__file__ and ROOT not in Path(epcsaft._core.__file__).parents:
        raise RuntimeError("declarations require this checkout's installed Engine")
    for owner in (epcsaft._core, eq, reg, Mixture, Pure):
        for name in dir(owner):
            if name.startswith(('solve_', 'continue_', 'stability_', 'evaluate_')) or name in ('fit', 'predict', 'differentiate', 'state', 'reference', 'compile_problem', 'observations_at'):
                setattr(owner, name, forbidden)
    for row, reason in HISTORICAL.items():
        print(f'excluded {row}: {reason}')
    sources = {path: path.read_text() for path in sorted((ROOT / 'analyses').rglob('*.py'))}
    built = [path for path, text in sources.items() if DECLARED.search(text)]
    covered = {leaf(path) for path in built}
    uncovered = []
    for path, text in sources.items():
        name = str(path.relative_to(ROOT / 'analyses'))
        # A leaf's declarations() declares the inputs its sibling scripts/*.py share.
        shared = leaf(path) in covered and path.parent == ROOT / 'analyses' / leaf(path) / 'scripts'
        if path in built or shared or any(name == row or name.startswith(row + '/') for row in HISTORICAL):
            continue
        if re.search(r'\b(?:import|from) epcsaft\b', text):
            uncovered.append(name)
    if uncovered:
        raise RuntimeError('uncovered Engine declarations; define declarations() or add an exclusion row: ' + ', '.join(uncovered))
    for path in built:
        name = str(path.relative_to(ROOT / 'analyses'))
        try:
            check(path)
        except Exception as error:
            raise RuntimeError(f'{path.relative_to(ROOT)}: {error}') from error
        print(f'declared {name}')


if __name__ == '__main__':
    main()
