"""Check the accepted transcriptions, including decimal strings and page locators."""
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_HASHES = {
    'hamer-wu-1972-table-values.csv': '9a672e7b2f5f31dbb282312efe0ea0916ea1113a37cd07353ed788efcb4b2f89',
    'haslam-2020-table8.csv': '715e00790905fc417955ccb4b9a526bd1aa60af840bbe090697d692a48dc4d90',
    'haslam-2020-parameters.csv': 'cec59ee1c0b6c9807774e0486c6357ac43a8c60da2e0b39b459727744a18d173',
}


def validate():
    for name, expected in SOURCE_HASHES.items():
        assert hashlib.sha256((ROOT / 'inputs/source' / name).read_bytes()).hexdigest() == expected, name
    with (ROOT / 'inputs/source/hamer-wu-1972-table-values.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len({(r['salt'], r['molality_mol_kg']) for r in rows}) == len(rows) == 207
    assert all(r['haslam_exact_row_selection'] == 'not independently established' for r in rows)
    return rows


if __name__ == '__main__':
    print(f'Source fidelity passed: {len(validate())} unique Hamer–Wu rows; all three accepted CSV hashes match.')
