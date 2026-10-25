"""Native fixed-temperature, fixed-pressure flashes for digitized ternary tie lines."""
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from PIL import Image, ImageDraw

from run import HERE, TOL, eq, model, write_rows


def main():
    start = time.perf_counter()
    cases = json.loads((HERE / 'inputs/digitized/ternary-pixels.json').read_text())
    inputs, outputs, diagnostics = [], [], []
    for figure, case in cases.items():
        vertices = np.array(case['vertices'], dtype=float)
        # Barycentric calibration maps the three plotted vertices to mole fractions.
        inverse = np.linalg.inv(np.vstack([vertices.T, np.ones(3)]))
        source = HERE / f'figures/figure-{figure}/source/paper-source-01-gross-2001-figure-0{figure}.png'
        image = Image.open(source).convert('RGB')
        m = model(case['components'], interactions=case['pairs'])
        for index, point in enumerate(case['points']):
            x = inverse @ np.array([*point['liquid_pixel'], 1.])
            y = inverse @ np.array([*point['vapor_pixel'], 1.])
            assert min(x) > 0 and min(y) > 0
            assert abs(sum(x) - 1.) < 1e-12 and abs(sum(y) - 1.) < 1e-12
            record = {'figure': figure, 'point': index, 'reference': point['reference'],
                'temperature_k': case['temperature_k'], 'pressure_pa': point['pressure_pa']}
            for phase in ('liquid', 'vapor'):
                px, py = point[phase + '_pixel']
                ImageDraw.Draw(image).ellipse((px-5, py-5, px+5, py+5), width=2,
                    outline='red' if point['reference'] == 'experiment' else 'blue')
            for i, component in enumerate(case['components']):
                inputs.append({**record, 'component': component, 'liquid_mole_fraction': x[i], 'vapor_mole_fraction': y[i]})
            problem = eq.Problem(T=case['temperature_k'], P=point['pressure_pa'],
                phases=[eq.Phase('L', kind='liquid', composition_guess=x),
                        eq.Phase('V', kind='vapor', composition_guess=y)],
                feed=eq.MoleFractions(dict(zip(case['components'], (x+y)/2))),
                row_tolerance=TOL['equilibrium_maximum_absolute_row_residual'])
            r = eq.solve_equilibrium(m, problem)
            diagnostic = {**record, 'success': r.success, 'message': r.message,
                'max_abs_residual': max(map(abs, r.residuals), default=float('inf'))}
            diagnostics.append(diagnostic)
            if not r.success:
                continue
            n = len(case['components'])
            xp, yp = np.array(r.mole_fractions[:n]), np.array(r.mole_fractions[n:])
            closure = max(abs(m.state(case['temperature_k'], rho=r.molar_densities[i],
                x=r.mole_fractions[i*n:(i+1)*n]).pressure / point['pressure_pa'] - 1.) for i in (0, 1))
            diagnostic.update(phase_pressure_relative_closure=closure,
                numerical_tolerance_met=closure <= TOL['phase_pressure_relative_closure'] and
                diagnostic['max_abs_residual'] <= TOL['equilibrium_maximum_absolute_row_residual'])
            for i, component in enumerate(case['components']):
                observed_k, predicted_k = y[i]/x[i], yp[i]/xp[i]
                outputs.append({**record, 'component': component,
                    'observed_liquid': x[i], 'observed_vapor': y[i],
                    'predicted_liquid': xp[i], 'predicted_vapor': yp[i],
                    'maximum_composition_absolute_error': max(abs(xp[i]-x[i]), abs(yp[i]-y[i])),
                    'observed_k': observed_k, 'predicted_k': predicted_k,
                    'k_absolute_relative_error_percent': 100*abs(predicted_k/observed_k-1.),
                    'model_tolerance_met': max(abs(xp[i]-x[i]), abs(yp[i]-y[i])) <= TOL['ternary_composition_absolute']
                    if point['reference'] == 'model' else '',
                    'numerical_tolerance_met': diagnostic['numerical_tolerance_met']})
        case['source_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
        image.save(HERE / f'figures/figure-{figure}/digitization-overlay.png')
    with (HERE / 'inputs/digitized/ternary-values.csv').open('w') as out:
        writer = csv.DictWriter(out, fieldnames=inputs[0], lineterminator="\n"); writer.writeheader(); writer.writerows(inputs)
    write_rows('ternary-comparisons.csv', outputs)
    write_rows('ternary-diagnostics.csv', diagnostics)
    (HERE / 'results/ternary-run.json').write_text(json.dumps({'wall_seconds': time.perf_counter()-start,
        'source_hashes': {f:c['source_sha256'] for f,c in cases.items()},
        'qa': 'review: approximate manual endpoint pairing; no independent extraction'}, indent=2)+'\n')
    print([(r['figure'], r['reference'], r['success']) for r in diagnostics])


if __name__ == '__main__':
    main()
