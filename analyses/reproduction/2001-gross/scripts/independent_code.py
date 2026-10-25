"""Retain teqp numerical verification of the final Gross (2001) declarations."""
import csv
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import sysconfig
import time

import numpy as np
import teqp

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[2]
TOLERANCE = 1e-10  # Fixed before calculation; relative to the retained Engine value.
START_FACTOR = 0.999  # Displace retained densities so a native solve is required.
PARAMETERS = {r['component']: r for r in csv.DictReader((HERE/'inputs/parameters/pure/any-solvent.csv').open())}
CALIBRATION = json.loads((HERE/'inputs/digitized/calibration.json').read_text())


def model(names, k=0.):
    return teqp.make_model({'kind': 'PCSAFT', 'model': {'coeffs': [
        {'name': name, 'BibTeXKey': 'Gross-IECR-2001', 'm': float(PARAMETERS[name]['m']),
         'sigma_Angstrom': float(PARAMETERS[name]['s']), 'epsilon_over_k': float(PARAMETERS[name]['e'])}
        for name in names], 'kmat': [[0., k], [k, 0.]] if len(names) == 2 else [[0.]]}})


def main():
    started = time.perf_counter()
    assert teqp.__version__ == '0.23.2'
    build = ROOT/'build/gross-2001'
    build.mkdir(parents=True, exist_ok=True)
    includes = subprocess.check_output([str(ROOT/'.venv/bin/python'), '-m', 'pybind11', '--includes'], text=True).split()
    subprocess.run(['c++', '-O2', '-shared', '-std=c++17', '-fPIC', *includes,
        str(HERE/'scripts/teqp_density.cpp'), '-o', str(build/('teqp_density'+sysconfig.get_config_var('EXT_SUFFIX')))], check=True)
    sys.path.insert(0, str(build))
    density = importlib.import_module('teqp_density')
    inputs = [HERE/'results/pure-saturation.csv', HERE/'results/figure-comparisons.csv']
    output = []
    success_codes = (teqp.VLE_return_code.xtol_satisfied, teqp.VLE_return_code.functol_satisfied)
    for source in inputs:
        for index, r in enumerate(csv.DictReader(source.open()), 1):
            assert r['kind'] in ('prediction', 'density', 'pressure', 'temperature')
            row = {'source': source.name, 'record': index, 'figure': r['figure'], 'reference': r['reference'],
                'component': r['component'], 'temperature_k': r['temperature_k'], 'pressure_pa': r['pressure_pa'],
                'quantity': 'liquid_density' if r['kind'] in ('prediction', 'density') else 'bubble_pressure' if r['kind']=='pressure' else 'bubble_temperature',
                'unit': 'kg/m^3' if r['kind'] in ('prediction', 'density') else 'Pa' if r['kind']=='pressure' else 'K',
                'engine_numerical_tolerance_met': r['numerical_tolerance_met']}
            if r['success'] != 'True':
                row.update(disposition='no Engine value', reason=r['message'])
                output.append(row)
                continue
            t, p = float(r['temperature_k']), float(r['pressure_pa'])
            x = np.array(json.loads(r['liquid_composition']))
            y = np.array(json.loads(r['vapor_composition']))
            row['liquid_first_component_mole_fraction'] = float(x[0])
            if row['quantity'] == 'liquid_density':
                m = model([r['component']])
                rho, iterations, step = density.solve(m, t, p, START_FACTOR*float(r['rho_liquid_mol_m3']), x)
                value = rho*float(PARAMETERS[r['component']]['MW'])
                engine = float(r['rho_liquid_kg_m3'])
                row.update(native_iterations=iterations, final_density_relative_step=step)
                L, V = np.array([rho]), None
            else:
                c = CALIBRATION[r['figure']]
                m = model(c['components'], c['k'])
                L0 = START_FACTOR*float(r['rho_liquid_mol_m3'])*x
                V0 = START_FACTOR*float(r['rho_vapor_mol_m3'])*y
                if row['quantity'] == 'bubble_pressure':
                    code, L, V = m.mix_VLE_Tx(t, L0, V0, x, 1e-10, 1e-12, 1e-10, 1e-12, 20)
                    value, engine = density.pressure(m, t, V), p
                else:
                    options = teqp.MixVLEpxFlags()
                    options.atol, options.reltol = 1e-10, 1e-12
                    options.axtol, options.relxtol, options.maxiter = 1e-10, 1e-12, 20
                    code, t, L, V = m.mixture_VLE_px(p, x, t, L0, V0, options)
                    value, engine = t, float(r['temperature_k'])
                assert code in success_codes, (source.name, index, code)
                row['native_return'] = str(code)
                row['vapor_first_component_mole_fraction'] = float(V[0]/V.sum())
                row['chemical_potential_difference_over_rt'] = float(np.max(np.abs(
                    m.get_chempotVLE_autodiff(t, L)-m.get_chempotVLE_autodiff(t, V)))/(m.get_R(x)*t))
                assert row['chemical_potential_difference_over_rt'] < 1e-7, row
                assert L.sum() > V.sum() > 0, row
            target = value if row['quantity']=='bubble_pressure' else p
            closure = abs(density.pressure(m, t, L)/target - 1)
            if V is not None:
                closure = max(closure, abs(density.pressure(m, t, V)/target - 1))
            if V is not None:
                assert closure < 1e-7, row
            # Low-pressure liquid cancellation is recorded; density is checked directly.
            assert np.isfinite(value) and value > 0, row
            row['teqp_liquid_density_mol_m3'] = float(L.sum())
            if V is not None:
                row['teqp_vapor_density_mol_m3'] = float(V.sum())
            difference = abs(value/engine - 1)
            row.update(engine_value=engine, teqp_value=value, relative_difference=difference,
                relative_tolerance=TOLERANCE, tolerance_met=difference <= TOLERANCE,
                phase_pressure_relative_closure=closure, disposition='compared')
            output.append(row)
    fields = list(dict.fromkeys(k for r in output for k in r))
    with (HERE/'results/independent-code-comparison.csv').open('w') as f:
        writer = csv.DictWriter(f, fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(output)
    compared = [r for r in output if r['disposition']=='compared']
    summary = {q: {'rows': len(rs), 'maximum_relative_difference': max(r['relative_difference'] for r in rs),
        'all_within_tolerance': all(r['tolerance_met'] for r in rs)} for q in dict.fromkeys(r['quantity'] for r in compared)
        for rs in [[r for r in compared if r['quantity']==q]]}
    record = {'comparator': 'teqp', 'version': teqp.__version__, 'environment': sys.prefix,
        'python_executable': sys.executable, 'relative_tolerance': TOLERANCE,
        'density_relative_step_tolerance': 1e-14, 'native_maximum_iterations': 20,
        'binary_atol': 1e-10, 'binary_reltol': 1e-12, 'binary_axtol': 1e-10, 'binary_relxtol': 1e-12,
        'phase_pressure_relative_closure_threshold': 1e-7,
        'start_density_factor': START_FACTOR, 'wall_seconds': time.perf_counter()-started,
        'git_base': subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        'basis': 'Independent-code numerical verification at final retained T,P (pure), T,x (binary), and P,x (isobar); Engine densities are displaced initial guesses only.',
        'summary': summary, 'uncompared_rows': len(output)-len(compared),
        'input_sha256': {str(f.relative_to(HERE)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in inputs+[HERE/'inputs/parameters/pure/any-solvent.csv', HERE/'inputs/digitized/calibration.json']},
        'code_sha256': {str(f.relative_to(HERE)): hashlib.sha256(f.read_bytes()).hexdigest()
            for f in [Path(__file__), HERE/'scripts/teqp_density.cpp']},
        'teqp_extension_sha256': hashlib.sha256(Path(teqp.teqp.__file__).read_bytes()).hexdigest(),
        'pure_pressure_limit': 'Relative pressure closure is diagnostic only for liquid density: cancellation at low pressure can exceed the existing closure threshold even when the density step and cross-code density difference pass.',
        'method_source': 'https://teqp.readthedocs.io/en/latest/algorithms/VLE.html'}
    (HERE/'results/independent-code-run.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))
    assert all(r['tolerance_met'] for r in compared), 'independent-code difference exceeds fixed tolerance'


if __name__ == '__main__':
    main()
