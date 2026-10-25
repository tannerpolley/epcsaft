"""Compare native Engine calcium laws with pinned native PHREEQC, without EOS solves."""
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import tarfile
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
LEAF = Path(__file__).resolve().parents[1]
BUILD = ROOT / 'build/issue-789/2018-iapws'
COMMIT = 'b0b3be767158ccc3322d2c816625cf470045e67e'
DATABASE_SHA = '93547b0343d9f151e73fb48e7927aa9e9c777399fedcb8c7497d00371af4d0ae'
TEMPERATURES = ((25, 298.15), (50, 323.15), (75, 348.15), (100, 373.15), (125, 398.15), (150, 423.15))
REACTIONS = {'calcium_bicarbonate': 'CaHCO3+', 'calcium_carbonate': 'CaCO3'}
RESULTS = LEAF / 'results'
SOURCE = LEAF / 'inputs/cases/carbonate-complexation/source/phreeqc-v3.9.0-carbonate-reactions.csv'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(args, **kwargs):
    return subprocess.run([str(a) for a in args], check=True, text=True, **kwargs)

BUILD.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(exist_ok=True)
os.environ.pop('EPCSAFT_HSL_LIBRARY', None)
source = BUILD / f'phreeqc3-{COMMIT}'
archive = BUILD / 'phreeqc-source.tar.gz'
if not archive.exists():
    urllib.request.urlretrieve(f'https://codeload.github.com/phreeqc-dev/phreeqc3/tar.gz/{COMMIT}', archive)
if not source.exists():
    with tarfile.open(archive) as bundle:
        bundle.extractall(BUILD, filter='data')
phreeqc = BUILD / 'phreeqc-build/phreeqc'
if not phreeqc.exists():
    run(['nice', '-n', '10', 'cmake', '-S', source, '-B', phreeqc.parent, '-DCMAKE_BUILD_TYPE=Release'])
    run(['nice', '-n', '10', 'cmake', '--build', phreeqc.parent, '-j', os.cpu_count()])
database = source / 'database/wateq4f.dat'
assert digest(database) == DATABASE_SHA, 'wrong pinned PHREEQC database'
meson = ROOT / 'build/engine-meson/shared'
library = meson / 'libepcsaft_native.so'
driver = BUILD / 'calcium_laws'
run(['nice', '-n', '10', 'c++', '-std=c++17', '-I' + str(ROOT / 'engine/native'),
     '-I' + str(ROOT / 'engine/build/subprojects/CppAD-20250000.3/include'),
     '-I' + str(ROOT / 'engine/build/subprojects/CppAD-20250000.3'),
     '-I' + str(meson / 'build/subprojects/CppAD-20250000.3'),
     LEAF / 'scripts/calcium_laws.cpp', '-L' + str(meson), '-Wl,-rpath,' + str(meson),
     '-lepcsaft_native', '-o', driver])
laws = {r['reaction_id']: r for r in csv.DictReader(SOURCE.open()) if r['reaction_id'] in REACTIONS}
assert set(laws) == set(REACTIONS)
assert all(float(r[c]) == 0 for r in laws.values() for c in ('analytic_A5', 'analytic_A6'))
request = ''.join(' '.join([name, str(t)] + [row[f'analytic_A{i}'] for i in range(1, 5)]) + '\n'
                  for name, row in laws.items() for _, t in TEMPERATURES)
raw = run([driver], input=request, capture_output=True).stdout
(RESULTS / 'engine-calcium-laws.csv').write_text(raw)
pqi = 'SELECTED_OUTPUT\n-file calcium.sel\n-reset false\n-high_precision true\nUSER_PUNCH\n-headings T_K pressure_atm CaHCO3 CaCO3\n-start\n10 PUNCH TK, PRESSURE, LK_SPECIES("CaHCO3+"), LK_SPECIES("CaCO3")\n-end\n'
pqi += ''.join(f'SOLUTION {i}\n-temp {c}\n-pressure 1\npH 7\nEND\n' for i, (c, _) in enumerate(TEMPERATURES, 1))
(RESULTS / 'phreeqc-calcium-laws.pqi').write_text(pqi)
completed = run([phreeqc, RESULTS / 'phreeqc-calcium-laws.pqi', BUILD / 'calcium.out', database], cwd=BUILD, capture_output=True)
(RESULTS / 'phreeqc-console.txt').write_text(completed.stdout + completed.stderr)
(RESULTS / 'phreeqc-selected.txt').write_bytes((BUILD / 'calcium.sel').read_bytes())
(RESULTS / 'phreeqc-output.txt').write_bytes((BUILD / 'calcium.out').read_bytes())
reference = [list(map(float, line.split())) for line in (BUILD / 'calcium.sel').read_text().splitlines()[1:] if line.strip()]
assert len(reference) == len(TEMPERATURES) and all(math.isfinite(r[1]) and r[1] > 0 for r in reference)
assert all(abs(row[0] - pair[1]) < 1e-9 for row, pair in zip(reference, TEMPERATURES))
engine = list(csv.DictReader(io.StringIO(raw)))
expected = {(name, t) for name in REACTIONS for _, t in TEMPERATURES}
assert len(engine) == len(expected) and {(r['reaction_id'], float(r['T_K'])) for r in engine} == expected
rows = []
for row in engine:
    name, t = row['reaction_id'], float(row['T_K'])
    ref = next(r for r in reference if abs(r[0] - t) < 1e-9)[2 + list(REACTIONS).index(name)]
    value = float(row['engine_log10_K'])
    difference = value - ref
    assert all(map(math.isfinite, (value, ref, difference)))
    rows.append(dict(reaction_id=name, T_K=t, engine_log10_K=value, phreeqc_log10_K=ref,
                     signed_difference=difference, absolute_difference=abs(difference), tolerance=1e-10,
                     outcome='pass' if abs(difference) <= 1e-10 else 'fail',
                     evidence_label='numerical verification against the PHREEQC source model'))
with (RESULTS / 'calcium-comparison.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
identity = {'engine_source_commit': run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True).stdout.strip(),
            'engine_tree_clean': not run(['git', 'status', '--porcelain', '--', 'engine'], cwd=ROOT, capture_output=True).stdout,
            'native_library': str(library), 'native_library_sha256': digest(library),
            'driver_sha256': digest(driver), 'source_input_sha256': digest(SOURCE),
            'phreeqc_source_commit': COMMIT, 'phreeqc_executable': str(phreeqc),
            'phreeqc_executable_sha256': digest(phreeqc), 'database_sha256': digest(database),
            'requested_pressure_atm': 1, 'actual_reference_pressures_atm': [r[1] for r in reference],
            'pressure_correction': 'zero reaction volumes in pinned WATEQ4F; PHREEQC raises hot solution pressure; no liquid-state qualification',
            'calcium_298K_nominal_log10_K': {k: v['log10_K_record'] for k, v in laws.items()},
            'nominal_record_diagnostic_tolerance': 0.0005, 'source_archive_sha256': digest(archive)}
(RESULTS / 'run-identity.json').write_text(json.dumps(identity, indent=2) + '\n')
assert all(r['outcome'] == 'pass' for r in rows), 'calcium source-law comparison failed; retain operands'
print(f'{len(rows)} calcium-law comparisons; maximum absolute log10 K difference {max(r["absolute_difference"] for r in rows):.3g}')
