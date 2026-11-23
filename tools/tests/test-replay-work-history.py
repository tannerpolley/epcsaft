"""Fit history may vary while fitted coordinates and scalar counters stay present."""
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
                if section == 'fitting':
                    writer.writerow(dict(case='fitting/219-held-NaCl', quantity='trial_failures', kind='int',
                                         value=0 if label == 'A' else 1))
                if label == 'A':
                    writer.writerow(dict(case='fitting/219-held-NaCl', quantity='iteration_costs/0', kind='float', value='3.196e-23'))
    command = [sys.executable, str(ROOT / 'engine/tests/cases/tpd-phase-discovery/scripts/replay_compare.py'),
               str(root / 'A'), str(root / 'B'), str(root / 'comparison.json')]
    assert subprocess.run(command, stdout=subprocess.DEVNULL).returncode == 0
    report = json.loads((root / 'comparison.json').read_text())
    assert all(not section['failures'] for section in report.values())
    assert [r['key'][1] for r in report['fitting']['work']] == ['iteration_costs/0', 'trial_failures']
    path = root / 'B/fitting.csv'
    original = path.read_text()
    for quantity in ('physical/0', 'trial_failures'):
        path.write_text(''.join(line for line in original.splitlines(keepends=True) if f',{quantity},' not in line))
        assert subprocess.run(command, stdout=subprocess.DEVNULL).returncode == 0
        failure, = json.loads((root / 'comparison.json').read_text())['fitting']['failures']
        assert failure['key'] == ['fitting/219-held-NaCl', quantity] and failure['reason'] == 'row missing'
    assert report['fitting']['rows'] == 2
    path.write_text(original)
    path = root / 'A/fitting.csv'
    path.write_text(path.read_text().replace('3.196e-23', 'nan'))
    assert subprocess.run(command, stdout=subprocess.DEVNULL).returncode == 0
    assert any(r['reason'] == 'non-finite value' for r in json.loads((root / 'comparison.json').read_text())['fitting']['failures'])
print('Variable work accepted; missing fitted coordinate and scalar counter refused: pass')
