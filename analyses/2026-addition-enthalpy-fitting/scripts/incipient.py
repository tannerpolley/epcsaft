"""Declare #269 boundaries; retain installed-wheel values and native diagnostics."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys

import epcsaft
from epcsaft import equilibrium as eq

ROOT = Path(__file__).resolve().parents[3]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--wheel', type=Path, required=True)
parser.add_argument('--source-commit', help='resolved source snapshot for a controlled replay')
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
sys.argv = ['zero_inventory.py', '--output', str(args.output / 'fresh-solvent.json'),
            '--native-input', str(args.output / 'mea-input.txt')]
source = runpy.run_path(str(Path(__file__).with_name('zero_inventory.py')))
def mea_boundary(feed=None, support=source['ids'][:3]):
    original = source['problem']
    return eq.Problem(T=original.T, P=eq.Free(), feed=original.feed if feed is None else feed,
        phases=[original.phases[0], eq.Phase('V', kind='vapor', amount=eq.Pinned(0), support=support)],
        reactions=original.reactions, neutral_reference=original.neutral_reference,
        standard_gibbs_j_per_mol=original.standard_gibbs_j_per_mol,
        topology_declared=original.topology_declared, row_tolerance=original.row_tolerance)
mea = mea_boundary()
sys.path.insert(0, str(ROOT / 'engine'))
from tests.native_inputs import parameters  # declaration helper; runtime remains installed
propane = epcsaft.Mixture(parameters([(44.096, 2.0020, 3.6184, 208.11)] * 3,
    pairs=[('k_ij', (0, 1), 0.), ('k_ij', (0, 2), 0.), ('k_ij', (1, 2), 0.)]))
cases = []
for dew in (False, True):
    problem = eq.Problem(T=300., P=eq.Free(1e6), feed=eq.Amounts({'s0': 1., 's1': 0., 's2': 0.}),
        phases=[eq.Phase('L', kind='liquid', amount=eq.Pinned(0) if dew else eq.Free(), composition_guess=[.8, .1, .1]),
                eq.Phase('V', kind='vapor', amount=eq.Free() if dew else eq.Pinned(0), composition_guess=[.8, .1, .1])],
        reactions=[eq.Reaction({'s1': -1., 's2': 1.}, equilibrium_constant=1., name='D=E')], row_tolerance=1e-10)
    cases.append(('propane-dew' if dew else 'propane-bubble', propane, problem))
cases.append(('mea-zero', source['model'], mea))
control = mea_boundary(support=source['ids'][1:3])
cases.append(('mea-control', source['model'], control))
for alpha in (1e-3, 1e-4, 1e-6, 1e-8):
    loaded = mea_boundary(feed=eq.Amounts(dict(zip(source['ids'], (alpha, 1., 7.911062165006594, 0., 0., 0., 0., 0., 0.)))))
    cases.append((f'mea-alpha-{alpha}', source['model'], loaded))
rows, failures = [], {}
for name, model, problem in cases:
    try:
        result = eq.solve_equilibrium(model, problem)
        rows.append(dict(case=name, pressure_Pa=result.pressure, amounts_mol=json.dumps(result.amounts),
            fractions=json.dumps(result.mole_fractions), densities_mol_m3=json.dumps(result.molar_densities),
            raw_rows=json.dumps(result.residuals), classified_rows=json.dumps(result.classified_residuals),
            success=result.success, message=result.message))
        if not result.success: failures[name] = result.message
    except Exception as error:
        failures[name] = str(error)
if rows:
    with (args.output / 'incipient-wheel.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
(args.output / 'incipient-wheel.json').write_text(json.dumps(dict(source_commit=args.source_commit or subprocess.check_output(
    ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(), wheel_sha256=hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
    module=epcsaft.__file__, inputs=source['out']['inputs'], failures=failures,
    limits='Numerical verification; no physical validation, global stability, or boundary caloric use.'), indent=2) + '\n')
print(json.dumps(failures), flush=True)
