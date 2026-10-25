"""Retain and read timing/work samples without calling the thermodynamic model."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import struct


def retain(path, metadata):
    samples = metadata['costs']
    absent = []
    for sample in samples:
        if 'native_work' in sample:
            for name in ('iterations', 'evaluations', 'hessian_calls'):
                if name not in sample['native_work']:
                    sample.pop(name, None)
                    absent.append({key: sample[key] for key in ('case', 'role', 'sample')} | {'field': name})
    metadata['unavailable_native_work'] = absent
    target = path.with_name(path.stem+'-costs.csv')
    if target.exists():
        raise RuntimeError(f'refusing to overwrite {target}')
    source = path.parent/metadata.get('numerical_rows_path', path.stem+'.csv')
    reported = {(row['case'], row['quantity']): row for row in csv.DictReader(source.open())}
    rows = []
    for sample in samples:
        prefix = {key: sample[key] for key in ('case', 'role', 'sample')}
        for name, value in sample.items():
            if name in prefix:
                continue
            values = value.items() if isinstance(value, dict) else [('', value)]
            for field, number in values:
                quantity = name+'/'+field if field else name
                kind = 'float' if isinstance(number, float) else 'int'
                row = dict(**prefix, quantity=quantity, kind=kind, value=repr(number),
                           binary64=struct.pack('>d', number).hex() if kind == 'float' else '',
                           reference_csv='', reference_case='', reference_quantity='')
                if name in ('iterations', 'evaluations', 'hessian_calls') and 'native_work' in sample:
                    assert number == sample['native_work'][name]
                    row.update(kind='alias', value='native_work/'+name, binary64='')
                report_case = sample['case']+'/'+sample['role']+'/work'
                operand = reported.get((report_case, field))
                if sample['sample'] == 0 and name == 'native_work' and operand is not None:
                    assert operand['kind'] == 'float' and operand['binary64'] == struct.pack('>d', number).hex()
                    row.update(kind='reference', value='', binary64='', reference_csv=os.path.relpath(source, path.parent),
                               reference_case=report_case, reference_quantity=field)
                rows.append(row)
    with target.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=('case', 'role', 'sample', 'quantity', 'kind', 'value', 'binary64', 'reference_csv', 'reference_case', 'reference_quantity'))
        writer.writeheader()
        writer.writerows(rows)
    metadata['costs'] = dict(csv=target.name, sha256=hashlib.sha256(target.read_bytes()).hexdigest(), samples=len(samples),
                             first_report_work='sample 0 native work references the identical reported CSV operands',
                             count_aliases='iterations/evaluations/Hessian counts alias native work where that map is present')
    metadata.pop('median3', None)
    assert read(path, metadata) == samples, 'cost extraction changed an observed operand'
    return metadata


def read(path, metadata):
    if isinstance(metadata['costs'], list):
        return metadata['costs']
    target = path.parent/metadata['costs']['csv']
    assert hashlib.sha256(target.read_bytes()).hexdigest() == metadata['costs']['sha256']
    groups = {}
    for row in csv.DictReader(target.open()):
        groups.setdefault((row['case'], row['role'], int(row['sample'])), {})[row['quantity']] = row
    references = {}
    result = []
    for (case, role, sample), fields in groups.items():
        def value(row):
            if row['kind'] == 'alias':
                return value(fields[row['value']])
            if row['kind'] == 'reference':
                source = path.parent/row['reference_csv']
                if source not in references:
                    references[source] = {(r['case'], r['quantity']): r for r in csv.DictReader(source.open())}
                return value(references[source][(row['reference_case'], row['reference_quantity'])])
            number = float(row['value']) if row['kind'] == 'float' else int(row['value'])
            if row['kind'] == 'float':
                assert struct.pack('>d', number).hex() == row['binary64']
            return number
        record = dict(case=case, role=role, sample=sample)
        for quantity, row in fields.items():
            if '/' in quantity:
                name, field = quantity.split('/', 1)
                record.setdefault(name, {})[field] = value(row)
            else:
                record[quantity] = value(row)
        result.append(record)
    return result


def write(path, metadata):
    if path.exists():
        raise RuntimeError(f'refusing to overwrite {path}')
    path.write_text(json.dumps(retain(path, metadata), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('metadata', type=Path, nargs='+')
    for path in parser.parse_args().metadata:
        original = path.read_bytes()
        metadata = json.loads(original)
        metadata = retain(path, metadata)
        metadata['retention'] = dict(original_json_sha256=hashlib.sha256(original).hexdigest(),
                                     transformation='timing/work samples moved to CSV without recalculation; medians are derived by present.py')
        path.write_text(json.dumps(metadata, indent=2)+'\n')
        print(path, metadata['costs']['samples'], 'samples preserved exactly')
