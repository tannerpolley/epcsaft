"""Run native executables in the isolated prefix and retain failed checks too."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess


def parse(output, internal):
    def vector(label, n):
        return [float(re.search(re.escape(label) + r'\[' + str(i) + r'\] = (\S+)', output)[1]) for i in range(n)]
    row = {key: vector(key, 4) for key in ('x', 'z_L', 'z_U')}
    row['g'] = [float(re.search(r'g\(' + str(i) + r'\) = (\S+)', output)[1]) for i in range(2)]
    row['objective'] = float(re.search(r'f\(x\*\) = (\S+)', output)[1])
    if internal:
        for key in ('diagnostics', 'scaled_diagnostics', 'work'):
            row[key] = [float(v) for v in re.search(r'^'+key+r'=(.*)$', output, re.M)[1].split(',')]
        row['lambda'] = vector('lambda', 2)
        row['bound_relax_factor'] = float(re.search(r'^bound_relax_factor=(.*)$', output, re.M)[1])
        row['internal'] = [[float(v) for v in line.split(',')[1:]] for line in output.splitlines() if line.startswith('internal,')]
        row['modules'] = re.findall(r'^module=(.*)$', output, re.M)
        assert 'version=3.11.9' in output and 'linear_solver=mumps' in output and 'status=0' in output
    return row


def check(row, reference, relaxation):
    for key in ('x', 'z_L', 'z_U', 'g', 'lambda', 'diagnostics', 'scaled_diagnostics', 'work'):
        if key in row:
            assert all(math.isfinite(v) for v in row[key]), key
    assert math.isfinite(row['objective'])
    # Lead amendment: original-bound feasibility includes Ipopt's unchanged
    # relative bound relaxation plus four ulps, not an altered solve option.
    bounds = [(1., 5.)]*4 + [(25., None), (40., 40.)]
    row['original_bound_checks'] = []
    for value, (lower, upper) in zip(row['x']+row['g'], bounds):
        for bound, sign in ((lower, -1), (upper, 1)):
            if bound is None:
                continue
            allowance = relaxation*max(1., abs(bound)) + 4*math.ulp(bound)
            violation = max(0., sign*(value-bound))
            row['original_bound_checks'].append(dict(value=value, bound=bound, violation=violation, allowance=allowance))
            assert violation <= allowance, row['original_bound_checks'][-1]
    # Compare the full primal/dual solution at the requested solver tolerance.
    row['linux_difference'] = {}
    for key in ('x', 'z_L', 'z_U', 'g', 'lambda'):
        if key in row:
            diff = [a-b for a, b in zip(row[key], reference[key])]
            row['linux_difference'][key] = diff
            assert all(abs(d) <= 1e-7*max(1., abs(b)) for d, b in zip(diff, reference[key])), (key, diff)
    row['linux_objective_difference'] = row['objective']-reference['objective']
    assert abs(row['linux_objective_difference']) <= 1e-7*max(1., abs(reference['objective']))
    if 'diagnostics' in row:
        assert row['diagnostics'][1] <= 1e-7 and row['diagnostics'][3] <= 1e-7
        assert row['scaled_diagnostics'][3] <= 1e-7
        assert row['internal'] and all(math.isfinite(v) for line in row['internal'] for v in line)
        assert all(v >= 0 and v.is_integer() for v in row['work'])


def main():
    prefix = Path(os.environ['PREFIX'])
    evidence = Path(os.environ.get('IPOPT_PROBE_EVIDENCE', 'ipopt-probe-evidence')).resolve()
    evidence.mkdir(parents=True, exist_ok=True)
    reference = json.loads(Path('linux-reference.json').read_text())
    results = {}
    codes = []
    for name in ('hs071', 'internal_caller'):
        exe = prefix / 'Library/bin' / (name + '.exe')
        run = subprocess.run([str(exe)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (evidence / (name + '.log')).write_text(run.stdout)
        print(run.stdout, flush=True)
        codes.append(run.returncode)
        results[name] = parse(run.stdout, name == 'internal_caller')
    # Write raw numbers before any numerical acceptance assertion.
    def retain():
        (evidence / 'numbers.json').write_text(json.dumps(results, indent=2, allow_nan=False))
        with (evidence/'solution.csv').open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['program', 'quantity', 'index', 'value'])
            for name, row in results.items():
                for key in ('x', 'z_L', 'z_U', 'g', 'lambda', 'diagnostics', 'scaled_diagnostics', 'work'):
                    for i, value in enumerate(row.get(key, [])):
                        writer.writerow([name, key, i, value])
                writer.writerow([name, 'objective', 0, row['objective']])
    retain()
    try:
        assert codes == [0, 0], codes
        factor = results['internal_caller']['bound_relax_factor']
        assert factor == reference['bound_relax_factor']
        for row in results.values():
            check(row, reference, factor)
        modules = results['internal_caller']['modules']
        private = [p for p in modules if Path(p).name.lower() == 'epcsaft_ipopt_3119.dll']
        assert len(private) == 1 and Path(private[0]).resolve() == (prefix/'Library/bin/epcsaft_ipopt_3119.dll').resolve(), modules
        assert not any(re.search(r'(^|[\\/])(lib)?ipopt[^\\/]*\.dll$', p, re.I) for p in modules)
    finally:
        retain()
        print(json.dumps(results, indent=2, allow_nan=False))
    hashes = {}
    for name in ('Library/bin/epcsaft_ipopt_3119.dll', 'Library/lib/epcsaft_ipopt_3119.lib'):
        hashes[name] = hashlib.sha256((prefix / name).read_bytes()).hexdigest()
    (evidence/'sha256.json').write_text(json.dumps(hashes, indent=2))
    print(hashes)
    print('Dependency numerical verification passed; Engine Windows qualification remains pending.')


if __name__ == '__main__':
    main()
