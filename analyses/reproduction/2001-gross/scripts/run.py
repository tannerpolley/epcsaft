"""Declare published PC-SAFT inputs; independent states are solved by the Engine."""
import csv
import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import subprocess
import time

import epcsaft
from epcsaft import Mixture, Parameters, Pure, equilibrium as eq, _core

HERE = Path(__file__).resolve().parents[1]
TOL = json.loads((HERE / 'inputs/tolerances.json').read_text())
PARAMETERS = {r['component']: r for r in csv.DictReader((HERE / 'inputs/parameters/pure/any-solvent.csv').open())}
CALIBRATION = json.loads((HERE / 'inputs/digitized/calibration.json').read_text())
STARTS = json.loads((HERE / 'inputs/phase-starts.json').read_text())


def model(components, k=0., interactions=None):
    records = []
    for name in components:
        r = PARAMETERS[name]
        records.append({'component_id': name, 'fixed': {
            'molar_mass': {'value': {'magnitude': float(r['MW']), 'unit': 'kilogram / mole'}},
            'charge_number': {'value': {'magnitude': 0., 'unit': 'elementary-charge'}}},
            'coefficients': [{'identity': f'{name}/{family}', 'family': family,
                'value': {'magnitude': float(r[column]), 'unit': unit}}
                for family, column, unit in [('segment_count', 'm', 'dimensionless'),
                ('segment_diameter', 's', 'angstrom'), ('dispersion_energy_over_k', 'e', 'kelvin')]]})
    interactions = interactions if interactions is not None else ([] if len(components) == 1 else [(0, 1, k)])
    pairs = [{'component_id_a': components[i], 'component_id_b': components[j],
        'coefficients': [{'identity': f'pair/{i}/{j}/k_ij', 'family': 'k_ij',
        'value': {'magnitude': value, 'unit': 'dimensionless'}}]} for i, j, value in interactions]
    p = Parameters.from_mapping({'schema': 'epcsaft.parameters', 'schema_version': 1,
        'components': records, 'pairs': pairs, 'model_families': [
        {'kind': kind, 'choice': choice} for kind, choice in
        [('base', 'pc-saft'), ('association', 'none'), ('electrolyte', 'none'),
         ('permittivity', 'none'), ('polar', 'none')]],
        'topology': {'sites': [], 'edges': []}})
    return (Pure if len(components) == 1 else Mixture)(p)


def pure_rows():
    rows = []
    for table in ('002', '003'):
        for r in list(csv.reader((HERE / f'results/tables/table-{table}/source/table-{table}.csv').open()))[2:]:
            if not r[1]:
                continue
            name = 'CO2' if r[0] == 'carbon dioxide' else r[0]
            p = PARAMETERS[name]
            assert [float(r[i]) for i in (2, 3, 4)] == [float(p[c]) for c in ('m', 's', 'e')]
            assert math.isclose(float(r[1]) / 1000, float(p['MW']), rel_tol=1e-14)
            low, high = map(float, r[7].split('-'))
            rows.append({'figure': 'table-2', 'component': name, 'kind': 'prediction',
                'temperature_k': low + .3 * (high - low), 'reference': 'none',
                'paper_psat_aad_percent': r[5], 'paper_volume_aad_percent': r[6],
                'paper_fit_temperature_range_k': r[7]})
    assert len(rows) == len(PARAMETERS)
    return rows


def declarations():
    return [model((r['component'],)) for r in pure_rows()] + [
        model(tuple(c['components']), c['k']) for c in CALIBRATION.values() if 'components' in c] + [
        model(c['components'], interactions=c['pairs']) for c in
        json.loads((HERE / 'inputs/digitized/ternary-pixels.json').read_text()).values()]


def calculate(row):
    row = dict(row)
    try:
        if row['kind'] in ('density', 'prediction'):
            m = model((row['component'],))
            t = float(row['y']) if row['kind'] == 'density' else row['temperature_k']
            row['temperature_k'] = t
            problem = m.saturation(t)
            use_coordinates = False
            if row['kind'] == 'density':
                sample = next(r for r in pure_rows() if r['component'] == row['component'])
                use_coordinates = t >= sample['temperature_k']
                if use_coordinates:
                    initial = eq.solve_equilibrium(m, m.saturation(sample['temperature_k']))
                    assert initial.success, initial.message
                row['initialization'] = 'native coordinates at prediction temperature' if use_coordinates else 'native cold seed'
        else:
            c = CALIBRATION[row['figure']]
            m = model(tuple(c['components']), c['k'])
            x = float(row['x'])
            feed = eq.MoleFractions({c['components'][0]: x, c['components'][1]: 1 - x})
            if row['kind'] == 'temperature':
                problem = eq.Problem(phases=[eq.Phase('L', kind='liquid'),
                    eq.Phase('V', kind='vapor', amount=eq.Pinned(0.))], T=eq.Free(),
                    P=c['pressure_pa'], feed=feed)
            else:
                y = STARTS['vapor_first_component_mole_fraction'][row['figure']][row['temperature_c']]
                problem = eq.Problem(phases=[eq.Phase('L', kind='liquid', composition_guess=[x, 1-x]),
                    eq.Phase('V', kind='vapor', amount=eq.Pinned(0.), composition_guess=[y, 1-y])],
                    T=float(row['temperature_c']) + 273.15, P=eq.Free(float(row['y']) * 1e5), feed=feed)
        problem.row_tolerance = TOL['equilibrium_maximum_absolute_row_residual']
        if row['kind'] == 'density' and use_coordinates:
            attempt = eq.solve_equilibrium_starts(m, problem,
                [_core.EquilibriumStart(coordinates=initial.coordinates)],
                _core.EquilibriumWorkLimits(attempts=1))
            r = attempt.results[0]
            if r is None:
                raise RuntimeError(attempt.causes[0])
        else:
            r = eq.solve_equilibrium(m, problem)
        row.update(success=r.success, message=r.message, solver_status=r.solver_status,
            iterations=r.iterations, evaluations=r.evaluations,
            max_abs_residual=max(map(abs, r.residuals), default=math.inf))
        if not r.success:
            return row
        t, p = r.temperature, r.pressure
        nc = len(m.component_ids)
        phases = [m.state(t, rho=r.molar_densities[i], x=r.mole_fractions[i*nc:(i+1)*nc]) for i in (0, 1)]
        row.update(temperature_k=t, pressure_pa=p, rho_liquid_mol_m3=r.molar_densities[0],
            rho_vapor_mol_m3=r.molar_densities[1],
            phase_pressure_relative_closure=max(abs(s.pressure - p) / p for s in phases),
            phase_pressure_absolute_closure_pa=max(abs(s.pressure - p) for s in phases),
            # Diagnostic only: pressure change of half a represented density step, (1/2)|P_rho| ulp(rho).
            phase_pressure_density_ulp_pa=max(.5 * abs(s.pressure_density_derivative) * math.ulp(s.molar_density)
                for s in phases),
            liquid_composition=json.dumps(r.mole_fractions[:nc]),
            vapor_composition=json.dumps(r.mole_fractions[nc:]),
            log_fugacity_difference=max(abs(phases[0].log_fugacity_coefficient[i]
                + math.log(r.mole_fractions[i]) - phases[1].log_fugacity_coefficient[i]
                - math.log(r.mole_fractions[nc+i])) for i in range(nc)))
        row['numerical_tolerance_met'] = (row['max_abs_residual'] <= TOL['equilibrium_maximum_absolute_row_residual']
            and row['phase_pressure_relative_closure'] <= TOL['phase_pressure_relative_closure'])
        if row['kind'] in ('density', 'prediction'):
            row['rho_liquid_kg_m3'] = r.molar_densities[0] * float(PARAMETERS[row['component']]['MW'])
            row['rho_vapor_kg_m3'] = r.molar_densities[1] * float(PARAMETERS[row['component']]['MW'])
        if row['kind'] != 'prediction':
            observed = float(row['x'] if row['kind'] == 'density' else row['y'])
            prediction = row['rho_liquid_kg_m3'] if row['kind'] == 'density' else t - 273.15 if row['kind'] == 'temperature' else p / 1e5
            error = abs(prediction - observed)
            relative, absolute = (TOL['figure_density_relative'], TOL['figure_density_absolute_kg_m3']) if row['kind'] == 'density' else (0., TOL['figure_temperature_absolute_k']) if row['kind'] == 'temperature' else (TOL['figure_pressure_relative'], TOL['figure_pressure_absolute_bar'])
            row.update(observed=observed, predicted=prediction, absolute_error=error,
                absolute_relative_error_percent=100 * error / abs(observed + 273.15 if row['kind'] == 'temperature' else observed),
                model_tolerance_met=error <= max(absolute, relative * abs(observed)) if row['reference'] == 'model' else '')
            if row['kind'] == 'density':
                row['volume_absolute_relative_error_percent'] = 100 * abs(observed / prediction - 1.)
    except Exception as error:
        row.update(success=False, message=f'{type(error).__name__}: {error}')
    return row


def write_rows(name, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE / 'results' / name).open('w') as out:
        writer = csv.DictWriter(out, fieldnames=fields, lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def main():
    start = time.perf_counter()
    inputs = pure_rows() + list(csv.DictReader((HERE / 'inputs/digitized/figure-values.csv').open()))
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
        rows = list(pool.map(calculate, inputs))
    write_rows('pure-saturation.csv', [r for r in rows if r['kind'] == 'prediction'])
    write_rows('figure-comparisons.csv', [r for r in rows if r['kind'] != 'prediction'])
    metadata = {'wall_seconds': time.perf_counter() - start, 'processes': os.cpu_count(),
        'engine_provenance': epcsaft.provenance(), 'source_parent': 'KXZN4KLQ',
        'source_pdf_attachment': 'ZAETGZK9',
        'pdf_sha256': '373b1c6224407933ee8e123bdd78b186b8bf49cabe1f17dffc9bd631d660bc72',
        'source_doi': '10.1021/ie0003887', 'git_base': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'total_rows': len(rows), 'successful_rows': sum(r.get('success', False) for r in rows),
        'input_hashes': {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (HERE / 'inputs').rglob('*') if p.is_file()}}
    (HERE / 'results/run.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps({k: metadata[k] for k in ('wall_seconds', 'total_rows', 'successful_rows')}), flush=True)


if __name__ == '__main__':
    main()
