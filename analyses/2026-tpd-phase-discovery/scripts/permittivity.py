"""#223 permittivity-form admission and defining low-density relations."""
import csv
import json
import math

import baseline as b

SPECIES = ('water', 'sodium-cation', 'chloride-anion')
SELECTIONS = b.ROOT/'analyses/2026-dielectric-alternatives/inputs/selections'
ASCANI = b.ROOT/'data/packets/ascani-2022-case-study-2-electrolyte-lle/3/parameters/parameters.json'
ZUBER = b.ROOT/'data/packets/zuber-2014-electrolyte-permittivity/1/parameters/ion-suppression.csv'


def qualify(q):
    from scipy.constants import Avogadro, Boltzmann
    x = (.98, .01, .01)
    for form in ('density-temperature', 'dipolar-perturbation', 'ion-fraction-suppression', 'ion-specific-suppression', 'component-mixing'):
        path = ASCANI if form == 'component-mixing' else SELECTIONS/('ion-fraction-suppression.json' if form == 'ion-specific-suppression' else form+'.json')
        document = json.loads(path.read_text())
        b.INPUTS['permittivity/'+form] = dict(path=str(path.relative_to(b.ROOT)), sha256=b.digest(path), components=SPECIES,
            defining_equations='full-eos sections 6.1-6.4; density_correlated_bulk_permittivity / effective_dipolar_density / dipolar_correlation_integral / dipolar_perturbation_bulk_permittivity')
        if form == 'ion-specific-suppression':
            next(f for f in document['model_families'] if f['kind'] == 'permittivity')['choice'] = form
            document['model_coefficients'] = [c for c in document['model_coefficients'] if c['family'] != 'ion_fraction_suppression_coefficient']
            with ZUBER.open(newline='') as stream:
                entries = {r['component_id']: r for r in csv.DictReader(stream) if r['solvent'] == 'water' and r['component_id'] in SPECIES}
            for component in document['components']:
                name = component['component_id']
                if name in entries:
                    row = entries[name]
                    component['coefficients'].append(dict(identity=f'component/{name}/ion_specific_suppression_coefficient',
                        family='ion_specific_suppression_coefficient', value=dict(magnitude=float(row['ion_specific_suppression_coefficient']), unit='dimensionless'),
                        provenance=dict(source_id='zuber-2014', locator=row['source_locator'], domain_id='aqueous-298-15-k')))
            b.INPUTS['permittivity/'+form].update(suppression_source=str(ZUBER.relative_to(b.ROOT)), suppression_sha256=b.digest(ZUBER),
                transformation='Explicit Zuber water ion coefficients replace the Figiel universal coefficient; other selected base inputs unchanged')
        model = b.Mixture(b.Parameters.from_mapping(document, components=SPECIES))
        problem = q.eq.Problem(T=298.15, P=1e5, row_tolerance=1e-10, feed=q.eq.Amounts(dict(zip(SPECIES, x))),
                               phases=[q.eq.Phase('L', kind='liquid', composition_guess=x)])
        case = 'permittivity/'+form
        if form in ('density-temperature', 'dipolar-perturbation'):
            components = {c['component_id']: {r['family']: r['value'] for r in c['coefficients']} for c in document['components']}
            if form == 'density-temperature':
                volume = components['water']['dielectric_volume']
                q.require(volume['unit'] == 'liter / mole', 'density dielectric volume requires source liter/mole conversion')
                dv = volume['magnitude']*1e-3
                dt = components['water']['dielectric_temperature']['magnitude']
                q.quantities(case+'/defining-inputs', dict(dielectric_volume_m3_per_mol=dv, dielectric_temperature_K=dt))
            else:
                mu2 = [components[name]['dipole_scaling']['magnitude'] for name in SPECIES]
                alpha = [components[name]['polarizability_scaling']['magnitude'] for name in SPECIES]
                a2 = [components[name]['correlation_integral_parameter']['magnitude'] for name in SPECIES]
                q.quantities(case+'/defining-inputs', dict(mu2_debye2=mu2, alpha_angstrom3=alpha, a2=a2,
                    Avogadro_per_mol=Avogadro, Boltzmann_J_per_K=Boltzmann))
            distances = []
            for rho in (1., 1e-3, 1e-6):
                if form == 'density-temperature':
                    expected = 1+rho*x[0]*dv*(dt/298.15-1)
                else:
                    debye2 = (1e-21/299792458.0)**2/(4*math.pi*8.8541878128e-12*1e-30)  # Debye^2 -> J A^3, SI Debye and CODATA 2018 epsilon_0
                    y = 4*math.pi/9*sum(rho*Avogadro*x[i]*1e-30*(mu2[i]*debye2/(Boltzmann*298.15)+3*alpha[i]) for i in range(3))
                    integral = 1+sum(a2[i]*x[i] for i in range(3))*math.expm1(-y)
                    expected = 1+3*y*(1+y+(17/16*integral-1)*y*y)
                actual = model.state(298.15, rho=rho, x=x).bulk_relative_permittivity
                error = abs(actual/expected-1)
                q.quantities(f'{case}/rho-{rho:g}', dict(rho_mol_m3=rho, epsilon_actual=actual, epsilon_reference=expected,
                    epsilon_relative_error=error, epsilon_distance_from_vacuum=abs(actual-1)))
                q.require(error <= 1e-13, f'{case}: defining permittivity equation differs by {error}')
                distances.append(abs(actual-1))
            q.require(distances[2] < distances[1] < distances[0], f'{case}: low-density sequence does not approach epsilon=1')
            central = q.timed(case, 'declared', lambda: q.eq.solve_equilibrium(model, problem))
            q.require(central.success and central.requested_tolerance_met, f'{case}: declared solve failed: {central.message}')
            compiled = q.eq.compile_problem(model, problem)
            report = q.timed(case, 'explicit-test', lambda: q.eq.stability_test(compiled, central, q.eq.StabilityOptions()))
            q.check_report(case, report, charges=(0, 1, -1))
            q.require(any(any(i in (1, 2) for i in t.phase.support) for t in report.trials), f'{case}: admitted charged family ran no charged trials')
        else:
            q.condensed_case(case, model, problem)
