"""A missing fit work sample must not hide a missing fitted coordinate (#940)."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(dir=ROOT / 'build') as directory:
    root = Path(directory)
    for label in 'AB':
        (root / label).mkdir()
        for section in ('declared', 'continuation', 'fitting', 'asl'):
            with (root / label / f'{section}.csv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=('case', 'quantity', 'kind', 'value', 'binary64'))
                writer.writeheader()
                writer.writerow(dict(case='fitting/219-held-NaCl', quantity='physical/0', kind='float', value='-0.25'))
                if label == 'A':
                    writer.writerow(dict(case='fitting/219-held-NaCl', quantity='iteration_costs/0', kind='float', value='3.196e-23'))
    command = [sys.executable, str(ROOT / 'analyses/2026-tpd-phase-discovery/scripts/replay_compare.py'),
               str(root / 'A'), str(root / 'B'), str(root / 'comparison.json')]
    assert subprocess.run(command, stdout=subprocess.DEVNULL).returncode == 0
    report = json.loads((root / 'comparison.json').read_text())
    assert all(len(section['work']) == 1 and not section['failures'] for section in report.values())
    path = root / 'B/fitting.csv'
    original = path.read_text()
    path.write_text(original.splitlines()[0] + '\n')
    assert subprocess.run(command, stdout=subprocess.DEVNULL).returncode == 0
    assert json.loads((root / 'comparison.json').read_text())['fitting']['failures']
    path.write_text(original)
    path = root / 'A/fitting.csv'
    path.write_text(path.read_text().replace('3.196e-23', 'nan'))
    assert subprocess.run(command, stdout=subprocess.DEVNULL).returncode == 0
    assert any(r['reason'] == 'non-finite value' for r in json.loads((root / 'comparison.json').read_text())['fitting']['failures'])
print('Variable work history accepted; missing fitted coordinate refused: pass')
