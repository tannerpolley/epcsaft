"""NaCl(s)–water: native equilibrium, source uncertainty and closure at measured temperatures."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, Parameters, _core
from epcsaft.equilibrium import (
    Molalities, NeutralReference, Phase, Problem, ReactionLogPolynomial,
    ReactionReference, ReferenceSourceBasis, Solid, solve_equilibrium, equilibrium_numerical_estimates,
)

HERE = Path(__file__).resolve().parents[1]
SPECIES = ['water', 'sodium-cation', 'chloride-anion']
P, P0, MW, MS, VS = 101300., 101325., .01801528, .05844277, 27.1e-6


def inputs(name):
    with (HERE / 'inputs' / name).open() as stream:
        return [{k: float(v) if k != 'source_locator' and k != 'salt' else v
                 for k, v in row.items()} for row in csv.DictReader(stream)]


def model(version=2, constant=False):
    doc = json.loads((HERE / f'inputs/held-{version}-malmberg-maryott.json').read_text())
    if constant:
        eps = next(c for c in doc['correlations'] if c['family'] == 'relative_permittivity')
        eps.update(form='constant', independent_variables=[])
        eps['constant']['value']['magnitude'] = 78.09
        eps.pop('terms')
        eps.pop('reference_temperature')
    return Mixture(Parameters.from_mapping(doc, components=SPECIES))


def problem(row, temperature=None, feed=12., dg_shift=0.):
    t0 = row['temperature_K']
    g, h, cp = row['g_J_mol'] + dg_shift, row['h_J_mol'], row['cp_J_mol_K']
    b, c = -(h - t0 * cp) / R, cp / R
    law = ReactionLogPolynomial(-g / (R * t0) - b / t0, b, c, 0., t0,
                                t0 - .1, t0 + .1, 'molality-infinite-dilution')
    solid = Solid('halite', formula={SPECIES[1]: 1., SPECIES[2]: 1.}, correlation=law,
                  reference=ReactionReference(ReferenceSourceBasis.CommonMolalityInfiniteDilution, P0),
                  molar_volume_m3_per_mol=VS)
    return Problem(T=t0 if temperature is None else temperature, P=P,
                   phases=[Phase('L', kind='liquid')], solids=[solid],
                   feed=Molalities({SPECIES[1]: feed, SPECIES[2]: feed},
                                  solvent=('water',), solvent_mass_kg=1.),
                   neutral_reference=NeutralReference(solvent='water', standard_molality_mol_per_kg=1.,
                                                     solute_molality_ratios={SPECIES[1]: 1., SPECIES[2]: 1.}))


def declarations():
    for version in (1, 2):
        model(version)
    for row in inputs('clarke-glew-1985-table34.csv'):
        problem(row)


def solve(mix, row, **options):
    request = problem(row, **options)
    result = solve_equilibrium(mix, request)
    record = {'temperature_K': request.T, 'pressure_Pa': P, 'reference_pressure_Pa': P0, 'feed_mol': options.get('feed', 12.),
              'dg_shift_J_mol': options.get('dg_shift', 0.),
              'success': result.success, 'message': result.message}
    if not result.amounts:
        return record | {'numerically_closed': False}
    n = [result.amounts[0] * x for x in result.mole_fractions]
    water = n[0] * MW
    m = n[1] / water
    ref = mix.reference(request.T, P, solvent_fractions=[1., 0., 0.], ions=SPECIES[1:], stoichiometry=(1, 1))
    ref0 = mix.reference(request.T, P0, solvent_fractions=[1., 0., 0.], ions=SPECIES[1:], stoichiometry=(1, 1))
    correction = sum(ref.log_fugacity) - sum(ref0.log_fugacity) + 2 * math.log(P / P0)
    gamma = ref.activity(m).mean_ionic_molality_scale
    state = mix.state(request.T, rho=result.molar_densities[0], x=result.mole_fractions)
    law = request.solids[0].correlation
    ln_k = law.a + law.b / request.T + law.c * math.log(request.T / law.reference_temperature)
    identity = 2 * math.log(m * gamma) - ln_k + correction - VS * (P - P0) / (R * request.T)
    identity -= math.log(10) * result.solid_affinity_decades[0]
    balance = [water - 1., (n[1] + result.solid_amounts[0] - record['feed_mol']) / record['feed_mol'],
               (n[2] + result.solid_amounts[0] - record['feed_mol']) / record['feed_mol']]
    charge = (n[1] - n[2]) / record['feed_mol']
    residual = max(map(abs, result.classified_residuals))
    record.update(molality_mol_kg=m, solid_mol=result.solid_amounts[0], water_kg=water,
                  si_decade=result.solid_affinity_decades[0], classified_max=residual,
                  residual_max=max(map(abs, result.residuals)), inventory_scaled_max=max(map(abs, balance)),
                  charge_scaled=charge, identity_ln=identity, reference_pressure_correction_ln=correction,
                  gamma_mean=gamma, ln_k=ln_k, density_mol_m3=state.molar_density,
                  pressure_error_Pa=state.pressure - P, dp_drho_Pa_m3_mol=state.pressure_density_derivative,
                  reference_convergence_error=ref.convergence_error,
                  reference_error_note='Analytic insertion limit: no terminal sequence estimate',
                  classified_residuals=list(result.classified_residuals), residuals=list(result.residuals))
    estimate = equilibrium_numerical_estimates(mix, result)
    record.update(source_reference_method=estimate.source_reference_method,
                  numerical_estimate_cause=estimate.cause, row_inner_error=estimate.row_inner,
                  row_floor=estimate.row_floor)
    record['numerically_closed'] = (result.success and math.isfinite(m) and m > 0
        and math.isfinite(state.molar_density) and state.molar_density > 0
        and math.isfinite(state.pressure_density_derivative) and state.pressure_density_derivative > 0
        and residual <= 1e-10 and max(map(abs, balance)) <= 1e-10 and abs(charge) <= 1e-10
        and result.solid_amounts[0] > 0 and abs(record['si_decade']) <= 1e-9 and abs(identity) < 1e-8)
    return record


def main():
    thermo = inputs('clarke-glew-1985-table34.csv')
    measured = inputs('pinho-macedo-2005.csv')
    mix = model()
    records, comparisons = [], []
    def attempt(row, variant='primary', **options):
        try:
            result = solve(mix if variant == 'primary' else model(1 if variant == 'packet-1' else 2,
                                                                 constant=variant == 'constant-epsilon'), row, **options)
        except Exception as error:
            result = {'temperature_K': options.get('temperature', row['temperature_K']),
                      'numerically_closed': False, 'message': str(error)}
        records.append(result | {'variant': variant, 'options': options})
        return result
    for row, obs in zip(thermo, measured, strict=True):
        assert row['temperature_K'] == obs['temperature_K']
        center = attempt(row)
        endpoints = [attempt(row, dg_shift=sign * row['g_se_J_mol']) for sign in (-1, 1)]
        thermal = [attempt(row, temperature=row['temperature_K'] + sign * .1) for sign in (-1, 1)]
        constant = attempt(row, 'constant-epsilon')
        w, u = obs['mass_fraction_salt'], obs['expanded_uncertainty_mass_fraction_95_percent']
        observed = w / (MS * (1 - w))
        lo, hi = [(w + sign * u) / (MS * (1 - w - sign * u)) for sign in (-1, 1)]
        comparison = {'temperature_K': row['temperature_K'], 'measured_mol_kg': observed,
                      'measured_low_95': lo, 'measured_high_95': hi,
                      'measured_u95_linear': u / (MS * (1 - w)**2),
                      'u_ln_k_standard': row['g_se_J_mol'] / (R * row['temperature_K']),
                      'numerically_closed': all(r['numerically_closed'] for r in [center, *endpoints, *thermal])}
        if comparison['numerically_closed']:
            low, high = sorted(r['molality_mol_kg'] for r in endpoints)
            m = center['molality_mol_kg']
            comparison.update(predicted_mol_kg=m, signed_error=m - observed,
                              relative_error_percent=100 * (m - observed) / observed,
                              predicted_low_se=low, predicted_high_se=high,
                              interval_disjoint=high < lo or low > hi,
                              temperature_low=min(r['molality_mol_kg'] for r in thermal),
                              temperature_high=max(r['molality_mol_kg'] for r in thermal),
                              constant_epsilon_mol_kg=constant.get('molality_mol_kg'))
        comparisons.append(comparison)
    feed = attempt(thermo[0], feed=24.)
    attempt(thermo[0], 'packet-1')
    primary = [r for r in records if r['variant'] == 'primary' and not r['options']]
    available = [r for r in comparisons if r['numerically_closed']]
    summary = {'attempted_temperatures': len(thermo), 'closed_temperatures': len(available),
               'native_attempts': len(records), 'native_failed': sum(not r['numerically_closed'] for r in records),
               'aard_percent': sum(abs(r['relative_error_percent']) for r in available) / len(available) if len(available) == len(thermo) else None,
               'max_absolute_error_mol_kg': max(abs(r['signed_error']) for r in available) if available else None,
               'mean_absolute_error_mol_kg': sum(abs(r['signed_error']) for r in available) / len(available) if available else None,
               'disjoint_count': sum(r.get('interval_disjoint', False) for r in comparisons),
               'feed_invariance_relative': (abs(feed['molality_mol_kg'] / primary[0]['molality_mol_kg'] - 1)
                 if feed.get('numerically_closed') and primary[0].get('numerically_closed') else None)}
    (HERE / 'results/native-states.json').write_text(json.dumps(records, indent=2) + '\n')
    (HERE / 'results/summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (HERE / 'results/comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, sorted({k for r in comparisons for k in r}), lineterminator='\n')
        writer.writeheader()
        writer.writerows(comparisons)
    paths = sorted(p for p in (HERE / 'inputs').rglob('*') if p.is_file())
    identity = {'native_module': _core.__file__, 'native_sha256': hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
                'input_sha256': {str(p.relative_to(HERE / 'inputs')): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (HERE / 'results/run-identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    assert len(primary) == 8 and all(r['numerically_closed'] for r in records)
    assert summary['feed_invariance_relative'] < 1e-8


if __name__ == '__main__':
    declarations() if '--declare-only' in sys.argv else main()
