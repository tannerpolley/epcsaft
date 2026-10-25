"""Retain the lead's charged-trial diagnostic outputs without rerunning them."""
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import struct

HERE = Path(__file__).resolve().parents[1]
SOURCES = [Path('/tmp/charged_trial_probe.out'), Path('/tmp/charged_trial_raw.out')]
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


def main():
    output = HERE/'results/diagnostics/charged-trial'
    output.mkdir(parents=True, exist_ok=True)
    assert not (output/'probe.csv').exists(), 'retained output already exists'
    rows, messages, empty, occurrences = [], {}, [], {}

    def record(case, name, value):
        rows.append(dict(case=case, quantity=name, kind='float', value=repr(value),
                         binary64=struct.pack('>d', value).hex()))

    for source in SOURCES:
        case = 'charged/artificial-zero-plane/raw-driver'
        for line in source.read_text().splitlines():
            if line.startswith('RESULT '):
                case = 'charged/'+line.split()[1]
                messages[case] = re.search(r'message=(.*?) objective=', line).group(1)
            elif line.startswith(('POLISH ', 'DIRECT ')):
                role = line.split()[0].lower()
                index = occurrences.get(role, 0)
                occurrences[role] = index+1
                case = f'charged/artificial-zero-plane/{role}-{index}'
            for name, value in re.findall(r'(\w+)=('+NUMBER+r')(?=\s|$)', line):
                record(case, name, float(value))
            fields = line.split()
            if not fields or not re.fullmatch(r'[A-Za-z_]+', fields[0]):
                continue
            if len(fields) > 1 and all(re.fullmatch(NUMBER, item) for item in fields[1:]):
                key = (case, fields[0])
                index = occurrences.get(key, 0)
                occurrences[key] = index+1
                name = fields[0]+(f'/{index}' if fields[0] == 'B' else '')
                for i, value in enumerate(fields[1:]):
                    record(case, f'{name}/{i}', float(value))
            elif len(fields) == 1 and fields[0] in ('coordinates', 'multipliers', 'composition', 'density'):
                empty.append(dict(case=case, quantity=fields[0]))
    with (output/'probe.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=('case', 'quantity', 'kind', 'value', 'binary64'))
        writer.writeheader()
        writer.writerows(rows)
    code = Path('/tmp/charged_trial_probe.cpp')
    retained_code = HERE/'scripts/charged_trial_probe.cpp'
    assert not retained_code.exists()
    shutil.copyfile(code, retained_code)
    (output/'probe.json').write_text(json.dumps(dict(
        sources={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in [*SOURCES, code]},
        retained_source=str(retained_code.relative_to(HERE)), messages=messages, unavailable_fields=empty,
        conditions=dict(temperature_K=298.15, pressure_Pa=100000., component_permutation=[1, 2, 0],
                        fixture='properties_support::aqueous_electrolyte'),
        units=dict(composition='mol/mol', density='mol/m3', raw_volumes='m3/mol',
                   objective='G/RT on the native trial basis', dual='G/RT on the native trial basis',
                   classified_rows='dimensionless row-scaled residuals',
                   coordinates='native mixed coordinate layout: log composition, amount in mol, log packing and scaled multipliers'),
        interpretation='The artificial zero plane drives the solvent toward the composition cutoff; the actual declared-candidate plane yields a stationary trial.',
        limits=['lead disposable native diagnosis, not physical validation',
                'compile/link identity was not supplied with these probe files',
                'read-only retention; no diagnostic is recompiled or rerun']), indent=2)+'\n')
    print(f'Retained {len(rows)} numeric diagnostic rows')


if __name__ == '__main__':
    main()
