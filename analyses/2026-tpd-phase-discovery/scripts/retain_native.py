"""Retain native numerical comparison operands without another calculation."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import struct


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--wheel-sha256', required=True)
    parser.add_argument('--compare', type=Path)
    parser.add_argument('--preparation', type=Path)
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--installation-log', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError('refusing to overwrite retained native operands or provenance')
    rows, logs, versions, comparisons = [], [], set(), set()

    def retain(prefix, field, value):
        if isinstance(value, dict):
            for name, item in value.items():
                retain(prefix, f'{field}/{name}' if field else name, item)
            return
        if isinstance(value, list):
            for index, item in enumerate(value):
                retain(prefix, f'{field}/{index}', item)
            return
        kind = 'bool' if isinstance(value, bool) else 'float' if isinstance(value, float) else 'none' if value is None else 'text'
        rows.append(dict(**prefix, field=field, kind=kind, value=repr(value) if kind == 'float' else str(value),
                         binary64=struct.pack('>d', value).hex() if kind == 'float' else ''))

    for path in args.log:
        contents = path.read_text()
        cases = re.findall(r'^EPCSAFT_CASE .*test_id=(\S+) execution=(\S+)', contents, re.MULTILINE)
        ids = {case for case, _ in cases}
        if len(ids) != 1:
            raise RuntimeError(f'{path}: expected one identified native test')
        test_id = ids.pop()
        count, asserted, observations = 0, 0, 0
        for line in contents.splitlines():
            if not line.startswith('EPCSAFT_COMPARE_JSON '):
                continue
            # Parse native numeric tokens as doubles, also preserving emitted -0.
            record = json.loads(line.split(' ', 1)[1], parse_int=float)
            versions.add(record.pop('schema_version'))
            comparison_id = record.pop('comparison_id')
            key = (test_id, comparison_id)
            if key in comparisons:
                raise RuntimeError(f'duplicate native comparison {key}')
            comparisons.add(key)
            prefix = dict(test_id=test_id, comparison_id=comparison_id,
                          quantity=record.pop('quantity'), comparison_unit=record.pop('unit'))
            retain(prefix, '', record)
            count += 1
            asserted += 'verdict' in record
            observations += 'verdict' not in record
        if not count:
            raise RuntimeError(f'{path}: no numerical operands')
        logs.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(), test_id=test_id,
                         execution=cases[-1][1], comparisons=count, asserted=asserted, observations=observations))
    args.output.mkdir(parents=True)
    target = args.output/'comparisons.csv'
    with target.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=('test_id', 'comparison_id', 'quantity', 'comparison_unit', 'field', 'kind', 'value', 'binary64'))
        writer.writeheader()
        writer.writerows(rows)
    comparison = None
    if args.compare is not None:
        identical = target.read_bytes() == args.compare.read_bytes()
        comparison = dict(identical=identical, reference=str(args.compare),
                          reference_sha256=hashlib.sha256(args.compare.read_bytes()).hexdigest())
        if identical:
            target.unlink()
            target = args.compare
    execution_identity = {}
    for name in ('preparation', 'executable', 'installation_log'):
        path = getattr(args, name)
        if path is not None:
            before = path.stat()
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            after = path.stat()
            assert (before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_ino, after.st_size, after.st_mtime_ns), f'{path} changed while recording its identity'
            execution_identity[name] = dict(path=str(path), sha256=digest)
    if args.preparation is not None:
        preparation = json.loads(args.preparation.read_text())
        assert preparation['wheel_sha256'] == args.wheel_sha256
        execution_identity['build_dir'] = preparation['build_dir']
    (args.output/'provenance.json').write_text(json.dumps(dict(logs=logs, native_format_versions=sorted(versions),
        associated_wheel_sha256=args.wheel_sha256, numerical_rows=str(target.name),
        numerical_rows_path=os.path.relpath(target, args.output), comparison=comparison,
        execution_identity=execution_identity,
        numerical_rows_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        columns='comparison_unit applies to model, reference and absolute tolerance; input units follow their emitted names and source definitions',
        limits=['read-only extraction of supplied native numerical logs',
                'wheel/run association supplied by the lead; optional binary and preparation hashes identify the inspected files',
                'no Engine calculation, build or installation; no physical validation']), indent=2)+'\n')
    print(f'Retained {len(comparisons)} comparisons and {len(rows)} operand/diagnostic rows')


if __name__ == '__main__':
    main()
