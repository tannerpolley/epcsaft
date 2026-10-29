"""NaCl/KCl co-saturation: seven native solves and the sampled interval falsifier."""
import hashlib
from itertools import product
import subprocess
from load_engine import ROOT, LIB
from run import HERE, MW, R, Mixture, Parameters, json, math, sys
from epcsaft import equilibrium as eq
from co_saturation_reference import SOURCE, thermochemistry

SPECIES = ['water', 'sodium-cation', 'potassium-cation', 'chloride-anion']
OUT = HERE / 'results/co-saturation'


def declarations():
    doc = json.loads((HERE / 'inputs/held-2-malmberg-maryott.json').read_text())
    sodium = next(p for p in doc['pairs'] if {p['component_id_a'], p['component_id_b']} == {SPECIES[0], SPECIES[1]})
    potassium = next(p for p in doc['pairs'] if {p['component_id_a'], p['component_id_b']} == {SPECIES[0], SPECIES[2]})
    law = json.loads(json.dumps(sodium['coefficients'][0]['temperature_law']).replace('sodium', 'potassium'))
    terms = law['terms']
    terms[0]['coefficient']['value']['magnitude'] = SOURCE['K_water_law']['intercept'] + SOURCE['K_water_law']['slope_per_K'] * 298.15
    terms[1]['coefficient']['value']['magnitude'] = SOURCE['K_water_law']['slope_per_K']
    potassium['coefficients'][0].pop('value')
    potassium['coefficients'][0]['temperature_law'] = law
    return Mixture(Parameters.from_mapping(doc, components=SPECIES))


def problem(thermo, pressure, feed, shifts):
    solids = [eq.Solid(name, formula={SPECIES[i]: 1., SPECIES[3]: 1.},
        equilibrium_constant=math.exp(-(g + dg) / (R * thermo['T_K'])),
        reference=eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution, SOURCE[salt]['source_pressure_Pa']),
        molar_volume_m3_per_mol=volume) for i, name, salt, g, dg, volume in zip(
            (1, 2), ('halite', 'sylvite'), ('NaCl', 'KCl'), thermo['G_J_mol'], shifts, SOURCE['solid_volume_m3_mol'])]
    return eq.Problem(T=thermo['T_K'], P=pressure, phases=[eq.Phase('L', kind='liquid')], solids=solids,
        feed=eq.Molalities(dict(zip(SPECIES[1:], [*feed, sum(feed)])), solvent=('water',), solvent_mass_kg=1.),
        neutral_reference=eq.NeutralReference(solvent='water', standard_molality_mol_per_kg=1.,
            solute_molality_ratios=dict(zip(SPECIES[1:], (1., 1., 2.)))))


def main():
    OUT.mkdir(exist_ok=True)
    mix = declarations()
    thermo = thermochemistry(R)
    t, p = thermo['T_K'], SOURCE['assumed_pressure_Pa']
    cases = [('central', p, (10., 8.), (0., 0.)), ('excess', p, (20., 16.), (0., 0.)),
             ('pressure', SOURCE['sensitivity_pressure_Pa'], (10., 8.), (0., 0.))]
    cases += [('corner', p, (10., 8.), shifts) for shifts in product(*[(-u, u) for u in thermo['G_halfwidth_J_mol']])]
    records = []
    for name, pressure, feed, shifts in cases:
        request = problem(thermo, pressure, feed, shifts)
        result = eq.solve_equilibrium(mix, request)
        record = dict(name=name, T_K=t, P_Pa=pressure, feed_mol=list(feed), shifts_J_mol=list(shifts))
        fields = 'amounts mole_fractions molar_densities solid_amounts solid_affinity_decades classified_residuals residuals'
        record.update({field: list(getattr(result, field)) for field in fields.split()})
        record.update(success=result.success, message=result.message)
        records.append(record)
        (OUT / 'native-states.json').write_text(json.dumps(records, indent=2) + '\n')
        assert result.success and result.amounts, record
        n = [result.amounts[0] * x for x in result.mole_fractions]
        water = n[0] * MW
        state = mix.state(t, rho=result.molar_densities[0], x=result.mole_fractions)
        errors = []
        for i, solid, g, shift in zip((1, 2), request.solids, thermo['G_J_mol'], shifts):
            p0, vs = solid.reference.reference_pressure_pa, solid.molar_volume_m3_per_mol
            r0 = mix.reference(t, p0, solvent_fractions=[1., 0., 0., 0.], ions=(SPECIES[i], SPECIES[3]), stoichiometry=(1, 1))
            affinity = (math.log(result.mole_fractions[i] * result.mole_fractions[3] / MW**2)
                + state.log_fugacity_coefficient[i] + state.log_fugacity_coefficient[3] - sum(r0.log_fugacity)
                + 2 * math.log(pressure / p0) + (g + shift) / (R * t) - vs * (pressure - p0) / (R * t))
            errors.append(affinity - math.log(10) * result.solid_affinity_decades[i - 1])
        balance = [n[0] - 1 / MW, n[1] + result.solid_amounts[0] - feed[0],
            n[2] + result.solid_amounts[1] - feed[1], n[3] + sum(result.solid_amounts) - sum(feed)]
        record.update(molality_mol_kg=[v / water for v in n[1:]], water_kg=water, balance_mol=balance,
            charge_mol=n[1] + n[2] - n[3], identity_error_ln=errors,
            dp_drho=state.pressure_density_derivative, pressure_error_Pa=state.pressure - pressure)
        positive = [*record['molality_mol_kg'], *result.solid_amounts, state.molar_density, state.pressure_density_derivative]
        record['numerically_closed'] = (all(math.isfinite(v) and v > 0 for v in positive)
            and max(map(abs, result.classified_residuals)) <= 1e-10 and max(map(abs, balance)) <= 1e-10
            and abs(record['charge_mol']) <= 1e-10 and max(map(abs, result.solid_affinity_decades)) <= 1e-9
            and max(map(abs, errors)) < 1e-8)
        (OUT / 'native-states.json').write_text(json.dumps(records, indent=2) + '\n')
        assert record['numerically_closed'], record
        print(name, record['molality_mol_kg'], flush=True)
    comparisons = []
    for i, salt in enumerate(SOURCE['measurement']['salts']):
        central, measured = records[0]['molality_mol_kg'][i], SOURCE['measurement']['molality_mol_kg'][i]
        calculated = [r['molality_mol_kg'][i] for r in records if r['name'] == 'corner']
        row = dict(salt=salt, measured=measured, predicted=central, signed_error=central - measured,
            relative_error_percent=100 * (central / measured - 1), sampled_range=[min(calculated), max(calculated)])
        for key in ('direct_intervals', 'mass_fraction_intervals'):
            lo, hi = thermo[key][i]
            separation = max(lo - max(calculated), min(calculated) - hi, 0.)
            row[key] = dict(measured=[lo, hi], separation_mol_kg=separation, verdict='disjoint' if separation > 0 else 'overlapping')
        comparisons.append(row)
    feed_error = [abs(b / a - 1) for a, b in zip(records[0]['molality_mol_kg'][:2], records[1]['molality_mol_kg'][:2])]
    assert max(feed_error) < 1e-8
    summary = dict(comparison=comparisons, excess_feed_relative_error=feed_error, thermochemistry=thermo)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    paths = [HERE / name for name in ('inputs/co-saturation.json', 'inputs/held-2-malmberg-maryott.json')]
    paths += list((HERE / 'scripts').glob('*.py')) + [LIB, *LIB.parent.glob('libepcsaft*.so')]
    identity = dict(source_HEAD=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        engine_tree=subprocess.check_output(['git', 'rev-parse', 'HEAD:engine'], cwd=ROOT, text=True).strip(),
        sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths})
    (OUT / 'run-identity.json').write_text(json.dumps(identity, indent=2) + '\n')


if __name__ == '__main__':
    declarations() if '--declare-only' in sys.argv else main()
