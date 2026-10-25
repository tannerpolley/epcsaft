"""Evaluate native Figiel declarations; compare retained quantities on explicit bases."""
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

from epcsaft import GAS_CONSTANT_J_PER_MOL_K, Mixture, Parameters, provenance, tolerances
from epcsaft import _core

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
OUT = LEAF / 'results/data'
C = json.loads((LEAF / 'inputs/comparison.json').read_text())
T, P, B = C['temperature_K'], C['pressure_Pa'], C['tolerances']
TERMS = ['hard_chain', 'dispersion', 'association', 'born', 'debye_huckel']
FIELDS = TERMS + ['residual_helmholtz', 'bulk_relative_permittivity', 'molar_density', 'pressure',
                  'pressure_density_derivative', 'compressibility_factor', 'log_fugacity_coefficient',
                  'component_property_status', 'density_convergence', 'solver_options']

def plain(v):
    if v is None or isinstance(v, (str, int, bool)):
        return v
    if isinstance(v, float):
        return v if math.isfinite(v) else None
    if isinstance(v, dict):
        return {k: plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [plain(x) for x in v]
    if hasattr(v, 'name') and not hasattr(v, 'value'):
        return v.name
    return {k: plain(getattr(v, k)) for k in dir(v) if not k.startswith('_') and not callable(getattr(v, k))}
def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

historical = {k: json.loads((ROOT / p).read_text()) for k, p in C['historical'].items()}
indexed = {}
for choice, data in historical.items():
    assert set(data) == set(C['salts']), (choice, 'salt set')
    for salt, source in data.items():
        assert source['components'] == C['salts'][salt], (choice, salt, 'component order')
        rows = {r['molality']: r for r in source['rows']}
        assert len(rows) == len(source['rows']) and set(rows) == set(C['molalities_mol_kg'])
        indexed[choice, salt] = rows
for salt in C['salts']:
    for m in C['molalities_mol_kg']:
        assert indexed['exclusion', salt][m]['mole_fractions'] == indexed['all-pairs', salt][m]['mole_fractions']
packet = json.loads((ROOT / C['packet']).read_text())
alternative = json.loads((LEAF / 'inputs/all-pairs.json').read_text())
expected = json.loads(json.dumps(packet))
next(f for f in expected['model_families'] if f['kind'] == 'electrolyte')['ionic_dispersion'] = 'all-ion-pairs'
assert alternative == expected, 'control must change only ionic_dispersion'
curve = list(csv.DictReader((ROOT / C['figure']).open()))
paths = [ROOT / C[k] for k in ['packet', 'figure', 'figure_method']] + [ROOT / p for p in C['historical'].values()]
paths += [LEAF / 'inputs/comparison.json', LEAF / 'inputs/all-pairs.json', Path(__file__)]
def declarations():
    return {(choice, salt): Mixture(Parameters.from_mapping(d, components=components))
            for choice, d in [('exclusion', packet), ('all-pairs', alternative)] for salt, components in C['salts'].items()}
def main():
    meta = {'gas_constant_J_per_mol_K': GAS_CONSTANT_J_PER_MOL_K, 'engine': plain(provenance()), 'extension': _core.__file__, 'extension_sha256': digest(Path(_core.__file__)),
            'command': sys.argv, 'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'git_status': subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True),
            'input_hashes': {str(p.relative_to(ROOT)): digest(p) for p in paths}, 'native_tolerances': plain(tolerances())}
    write('run.json', meta)  # Record declarations and effective policy before the campaign.
    checks, states, activities, references, failures = [], [], [], [], []

    def check(group, context, quantity, native, reference, tolerance, unit='dimensionless'):
        delta = native - reference if native is not None and reference is not None else None
        checks.append(dict(group=group, **context, quantity=quantity, native=native, reference=reference, delta=delta,
                           absolute_error=abs(delta) if delta is not None else None, tolerance=tolerance, unit=unit,
                           passed=abs(delta) <= tolerance if delta is not None else None))

    def attempt(context, operation):
        try:
            return operation()
        except Exception as error:
            failures.append(dict(context, cause=f'{type(error).__name__}: {error}'))
            return None

    started = time.perf_counter()
    models = declarations()
    for choice in historical:
        for salt, components in C['salts'].items():
            base = dict(convention=choice, salt=salt, molality='', chart='reference')
            model = models[choice, salt]
            meta['solver_options'] = plain(model._model.solver_options)
            write('run.json', meta)
            ref = attempt(base, lambda: model.reference(T, P, solvent_fractions=(1., 0., 0.), ions=components[1:], stoichiometry=C['stoichiometry']))
            if ref is not None:
                references.append(dict(base, **{k: plain(getattr(ref, k)) for k in ['log_fugacity', 'molality', 'convergence_error']}))
                for i, value in enumerate(ref.log_fugacity):
                    check('reference', base, f'log_phi_ion_{i}', value, historical[choice][salt]['reference_log_phi_inf'][i], B['log_fugacity_absolute'])
            for m, old in indexed[choice, salt].items():
                context = dict(base, molality=m, chart='activity')
                activity = attempt(context, lambda: ref.activity(m)) if ref is not None else None
                if activity is not None:
                    gamma = activity.mean_ionic_molality_scale
                    activities.append(dict(context, gamma=gamma, log_gamma=activity.mean_ionic_log_molality_scale, solver_options=plain(activity.solver_options)))
                    check('gamma', context, 'gamma_m', gamma, old['gamma_pm_m'], B['gamma_relative'] * abs(old['gamma_pm_m']))
                for chart in ['pressure', 'exclusion-density', 'all-pairs-density']:
                    context = dict(base, molality=m, chart=chart)
                    source = indexed[chart.removesuffix('-density'), salt][m] if chart != 'pressure' else old
                    args = dict(P=P, phase='liquid') if chart == 'pressure' else dict(rho=source['molar_density'])
                    state = attempt(context, lambda: model.state(T, x=source['mole_fractions'], **args))
                    if state is None:
                        continue
                    values = {k: plain(getattr(state, k)) for k in FIELDS}
                    states.append(dict(context, mole_fractions=source['mole_fractions'], **values))
                    states[-1]['log_fugacity_unavailable_cause'] = 'nonpositive Z: ln(Z) undefined' if state.compressibility_factor <= 0 else None
                    if any(values[k] is None for k in TERMS + ['residual_helmholtz']):
                        failures.append(dict(context, cause='nonfinite required reduced term'))
                    z, phi = state.compressibility_factor, state.log_fugacity_coefficient
                    if z > 0 and phi is not None and all(v is not None for v in phi):
                        check('euler', context, 'sum_x_log_phi', sum(x*v for x, v in zip(source['mole_fractions'], phi)), state.residual_helmholtz + z - 1 - math.log(z), B['euler_absolute'])
                    if chart != 'pressure':
                        continue
                    conv = state.density_convergence
                    check('root', context, 'reduced_pressure_residual', conv.residual, 0., conv.criterion.value)
                    check('root', context, 'positive_slope', float(state.pressure_density_derivative > 0), 1., 0.)
                    states[-1]['pressure_error_Pa'] = state.pressure - P
                    for field in TERMS + ['residual_helmholtz', 'bulk_relative_permittivity', 'molar_density']:
                        target = sum(old[k] for k in TERMS) if field == 'residual_helmholtz' else old[{'bulk_relative_permittivity': 'permittivity'}.get(field, field)]
                        tol = B['density_relative']*abs(target) if field == 'molar_density' else B['permittivity_absolute'] if field == 'bulk_relative_permittivity' else B['reduced_energy_scaled_absolute']*max(1, abs(target))
                        check('history', context, field, values[field], target, tol, 'mol/m^3' if field == 'molar_density' else 'dimensionless')
                    for i, target in enumerate(old['log_phi']):
                        check('history', context, f'log_phi_{i}', phi[i] if phi is not None else None, target, B['log_fugacity_absolute'])
            if salt == 'LiBr' and ref is not None:
                for row in curve:
                    m, target, uncertainty = (float(row[k]) for k in ['molality_mol_kg', 'gamma_pm_m', 'digitization_uncertainty'])
                    context = dict(base, molality=m, chart='figure-5')
                    activity = attempt(context, lambda: ref.activity(m))
                    gamma = activity.mean_ionic_molality_scale if activity else None
                    activities.append(dict(context, gamma=gamma, source_gamma=target, uncertainty=uncertainty,
                                           signed_relative_error=gamma/target-1 if gamma is not None else None,
                                           absolute_relative_error=abs(gamma/target-1) if gamma is not None else None,
                                           uncertainty_ratio=abs(gamma-target)/uncertainty if gamma is not None else None))
            print(choice, salt, 'complete', flush=True)
    for a in [s for s in states if s['convention'] == 'exclusion' and s['chart'] == 'exclusion-density']:
        b = next((s for s in states if s['convention'] == 'all-pairs' and all(s[k] == a[k] for k in ['salt', 'molality', 'chart'])), None)
        if b is None:
            continue
        context = {k: a[k] for k in ['convention', 'salt', 'molality', 'chart']}
        for field in ['hard_chain', 'association', 'born', 'debye_huckel', 'bulk_relative_permittivity']:
            check('isolation', context, field, b[field], a[field], B['matched_scaled_absolute']*max(1, abs(a[field])))
        check('isolation', context, 'total_minus_dispersion_change', b['residual_helmholtz']-a['residual_helmholtz'], b['dispersion']-a['dispersion'], B['matched_scaled_absolute']*max(1, abs(a['residual_helmholtz'])))
        if a['salt'] == 'KCl' and a['molality'] == 2:
            check('witness', context, 'resolved_dispersion_change', float(abs(b['dispersion']-a['dispersion']) > B['witness_minimum_absolute']), 1., 0.)
    meta.update(wall_seconds=time.perf_counter()-started, state_calls=len(states), activity_calls=len(activities), reference_calls=len(references), failed_calls=len(failures))
    for name, data in [('run.json', meta), ('states.json', states), ('activities.json', activities), ('references.json', references), ('failures.json', failures)]:
        write(name, data)
    with (OUT / 'crosswalk.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(checks[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(checks)
    print(json.dumps(meta | {'input_hashes': 'see run.json', 'native_tolerances': 'see run.json'}, indent=2))
if __name__ == '__main__':
    main()
