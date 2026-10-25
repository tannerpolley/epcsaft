"""Engine side of the 2026-10-08 head-to-head comparisons; Python declares inputs and reads native results.

Writes OUT (JSON): Khudaida electrolyte LLE (packet v7 and the Clapeyron-matched variant), exact C2/C4
solved-state actions with central-difference checks, Born/permittivity fixed states, and the
declared inputs the FeOs and Clapeyron scripts read.
"""
import copy, csv, json, statistics, sys, time
from pathlib import Path
from epcsaft import Mixture, Parameters, provenance, regression as r
from epcsaft import equilibrium as q
HERE = Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]; CC = ROOT / 'analyses/2026-competitor-cross-check'
sys.path.insert(0, str(CC / 'scripts'))
import engine as cc  # the cross-check's C2/C4 parameter documents

KH = ROOT / 'data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json'
KH_IDS = ('water', 'ethanol', 'isobutanol', 'sodium-cation', 'chloride-anion')
FEEDS = ROOT / 'analyses/reproduction/2026-khudaida/figures/figure-02/source/source-points.csv'
SELECTIONS = ROOT / 'analyses/2026-dielectric-alternatives/inputs/selections'
NACL = ('water', 'sodium-cation', 'chloride-anion')
MOLALITIES = (0.1, 1.0, 3.0)


def timed(f, warm=3, samples=10):
    first = time.perf_counter_ns(); value = f(); first = time.perf_counter_ns() - first
    for _ in range(warm): f()
    ns = []
    for _ in range(samples):
        t = time.perf_counter_ns(); f(); ns.append(time.perf_counter_ns() - t)
    return value, dict(first_ms=first / 1e6, median_ms=statistics.median(ns) / 1e6, min_ms=min(ns) / 1e6, max_ms=max(ns) / 1e6, samples=samples)


def coefficient(component, family):
    return next(c for c in component['coefficients'] if c['family'] == family)


# ---- 2. Khudaida electrolyte LLE ----
def khudaida_document(variant):
    d = json.loads(KH.read_text())
    if variant in ('matched', 'constant-permittivity'):  # direct Born and DH at sigma: the forms Clapeyron's DHBorn offers
        next(f for f in d['model_families'] if f['kind'] == 'electrolyte').update(c_shell=0.0, c_dielectric=0.0)
        for c in d['components']:
            if c['component_id'] in KH_IDS[3:]:
                coefficient(c, 'debye_huckel_diameter')['value']['magnitude'] = coefficient(c, 'segment_diameter')['value']['magnitude']
    if variant == 'constant-permittivity':  # every solvent at water's 78.09: a composition-independent permittivity (Clapeyron ConstRSP)
        next(f for f in d['model_families'] if f['kind'] == 'permittivity')['choice'] = 'solvent-only'
        d['model_coefficients'] = [m for m in d['model_coefficients'] if m['family'] != 'ion_fraction_suppression_coefficient']
        for c in d['correlations']:
            if c['family'] == 'relative_permittivity': c['constant']['value']['magnitude'] = 78.09
    return d


def khudaida_feeds():
    out = {}
    for row in csv.DictReader(FEEDS.open()):
        if row['dataset'] == 'feed':
            f = [float(row[k]) for k in ('x_water', 'x_ethanol', 'x_isobutanol', 'x_nacl')]
            f = f + [f[3]]; s = sum(f); out[int(row['point_id'])] = [v / s for v in f]
    return dict(sorted(out.items()))


def khudaida(variant):
    model = Mixture(Parameters.from_mapping(khudaida_document(variant), components=KH_IDS)); rows = []
    for point, feed in khudaida_feeds().items():
        problem = q.Problem(T=293.15, P=1e5, feed=q.Amounts(dict(zip(KH_IDS, feed))), row_tolerance=1e-10, phases=[
            q.Phase('phase-0', kind='liquid', support=KH_IDS, composition_guess=(0.80, 0.05, 0.05, 0.05, 0.05)),
            q.Phase('phase-1', kind='liquid', support=KH_IDS, composition_guess=(0.20, 0.10, 0.6998, 0.0001, 0.0001))])
        res, cost = timed(lambda: q.solve_equilibrium(model, problem))
        n = len(KH_IDS); xs = [list(res.mole_fractions[i * n:(i + 1) * n]) for i in range(res.phase_count)]
        order = sorted(range(res.phase_count), key=lambda i: -xs[i][0])
        rows.append(dict(point=point, feed=feed, success=res.success, message=res.message, phase_count=res.phase_count,
                         x=[xs[i] for i in order], rho=[res.molar_densities[i] for i in order],
                         beta=res.amounts[order[0]] / sum(res.amounts), iterations=res.iterations, evaluations=res.evaluations,
                         max_residual=max(map(abs, res.residuals)), cost=cost))
    return rows


def clapeyron_khudaida_spec():
    """Matched-variant inputs for Clapeyron, read from the same Engine document (no re-typed values)."""
    d = khudaida_document('matched'); comps = {c['component_id']: c for c in d['components']}
    corr = {(c['component_id'], c['family']): c for c in d['correlations']}
    def val(cid, fam):
        c = comps[cid]
        if any(x['family'] == fam for x in c['coefficients']): return coefficient(c, fam)['value']['magnitude']
        return corr[(cid, fam)]['constant']['value']['magnitude']
    edges = {(e['endpoint_a']['component_id'], e['endpoint_b']['component_id']): e for e in d['topology']['edges'] if e['source']['kind'] == 'explicit'}
    out = []
    for cid in KH_IDS:
        z = comps[cid]['fixed']['charge_number']['value']['magnitude']; e = edges.get((cid, cid))
        out.append(dict(id=cid, z=z, Mw=comps[cid]['fixed']['molar_mass']['value']['magnitude'], m=val(cid, 'segment_count'),
                        sigma=val(cid, 'segment_diameter'), epsilon=val(cid, 'dispersion_energy_over_k'),
                        epsilon_ab=e['energy_over_k']['value']['magnitude'] if e else 0., kappa_ab=e['volume']['value']['magnitude'] if e else 0.,
                        relative_permittivity=val(cid, 'relative_permittivity') if z == 0 else None,
                        born_diameter=val(cid, 'born_diameter') if z else None))
    k = [[0.] * 5 for _ in range(5)]
    for p in d['pairs']:
        i, j = KH_IDS.index(p['component_id_a']), KH_IDS.index(p['component_id_b']); k[i][j] = k[j][i] = coefficient(p, 'k_ij')['value']['magnitude']
    suppression = next(m for m in d['model_coefficients'] if m['family'] == 'ion_fraction_suppression_coefficient')['value']['magnitude']
    return dict(components=out, k=k, ion_fraction_suppression=suppression, T=293.15, P=1e5, feeds=khudaida_feeds())


# ---- 3. Exact sensitivities ----
def sensitivity_case(inputs, ident):
    req = next(x for x in inputs['requests'] if x['id'] == ident); spec = copy.deepcopy(inputs['models'][req['case']])
    ids = [c['id'] for c in spec['components']]; c0 = spec['components'][0]
    coords = [('k_ij', tuple(ids), 'k')] + ([('association_energy_over_k', ((ids[0], 'a'), (ids[0], 'b')), 'epsilon_ab'),
                                             ('association_volume', ((ids[0], 'a'), (ids[0], 'b')), 'kappa_ab')] if req['case'] == 'C2' else [])
    def value_of(s, key): return s['k'] if key == 'k' else s['components'][0][key]
    def build(s, T, active=False):
        p = Parameters.from_mapping(cc.document(s, req), components=ids)
        kw = {}
        if active:
            chart = r.CoordinateSet(); chart.coordinates = [r.coordinate(p, [(f, t, 1)], origin=value_of(s, key), scale=1) for f, t, key in coords]
            kw['active_parameters'] = chart.active()
        model = Mixture(p, **kw)
        if req['case'] == 'C2':
            prob = q.bubble_point(T=T, liquid=q.Phase('L', kind='liquid'), feed=q.MoleFractions(dict(zip(ids, req['x']))))
        else:
            prob = q.Problem(T=T, P=req['P'], feed=q.MoleFractions(dict(zip(ids, req['x']))), phases=[
                q.Phase('water-rich', kind='liquid', composition_guess=[.99, .01]), q.Phase('water-poor', kind='liquid', composition_guess=[.01, .99])])
        prob.row_tolerance = 1e-12
        return model, prob
    def observe(res):  # observable values in a fixed order, phases sorted by first-component fraction
        xs = [res.mole_fractions[2 * i:2 * i + 2] for i in range(res.phase_count)]
        if req['case'] == 'C2': return [res.pressure, xs[1][0]]
        rich = max(range(2), key=lambda i: xs[i][0])
        return [xs[rich][0], xs[1 - rich][0], res.amounts[rich] / sum(res.amounts)]
    names = ['P_Pa', 'y_first'] if req['case'] == 'C2' else ['x_water_rich', 'x_water_poor', 'beta_water_rich']
    model, prob = build(spec, req['T'], active=True)
    central = q.solve_equilibrium(model, prob)
    xs = [central.mole_fractions[2 * i:2 * i + 2] for i in range(central.phase_count)]
    K = q.SolvedStateObservableKind
    if req['case'] == 'C2':
        observables = [q.SolvedStateObservable(K.Pressure), q.SolvedStateObservable(K.PhaseMoleFraction, 1, 0)]
    else:
        rich = max(range(2), key=lambda i: xs[i][0])
        observables = [q.SolvedStateObservable(K.PhaseMoleFraction, rich, 0), q.SolvedStateObservable(K.PhaseMoleFraction, 1 - rich, 0),
                       q.SolvedStateObservable(K.PhaseAmount, rich, 0)]
    nc = len(coords)
    directions = [('T', q.SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [0.0] * nc))] + \
        [(f, q.SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [float(i == j) for j in range(nc)])) for i, (f, _, _) in enumerate(coords)]
    def exact():
        res = q.solve_equilibrium(model, prob)
        return q.solved_state_first_actions(q.compile_problem(model, prob), res, observables, [d for _, d in directions])
    batches, exact_cost = timed(exact)
    def actions_only(): return q.solved_state_first_actions(q.compile_problem(model, prob), central, observables, [d for _, d in directions])
    _, action_cost = timed(actions_only)
    rows = []
    for (name, _), batch in zip(directions, batches):
        for obs_name, item in zip(names, batch.results):
            rows.append(dict(case=ident, observable=obs_name, direction=name, value=item.value, exact=item.action, status=item.status.name))
    # Central-difference check: re-solve the declared problem at shifted inputs.
    steps = {'T': 1e-3, 'k_ij': 1e-5, 'association_energy_over_k': 1e-5 * c0.get('epsilon_ab', 1), 'association_volume': 1e-5 * c0.get('kappa_ab', 1)}
    def shifted(name, sign):
        s = copy.deepcopy(spec); T = req['T']
        if name == 'T': T += sign * steps['T']
        else:
            key = next(k for f, _, k in coords if f == name)
            if key == 'k': s['k'] += sign * steps[name]
            else: s['components'][0][key] += sign * steps[name]
        m, p = build(s, T); res = q.solve_equilibrium(m, p)
        if not res.success: raise RuntimeError(res.message)
        return observe(res)
    def finite_differences():
        return {name: [(a - b) / (2 * steps[name]) for a, b in zip(shifted(name, 1), shifted(name, -1))] for name, _ in directions}
    fd, fd_cost = timed(finite_differences, warm=0, samples=3)
    for row in rows:
        row['fd'] = fd[row['direction']][names.index(row['observable'])]; row['fd_step'] = steps[row['direction']]
    return dict(rows=rows, central=observe(central), names=names, coordinates=[f for f, _, _ in coords],
                cost=dict(solve_plus_actions=exact_cost, actions_on_solved_state=action_cost, central_fd_all_directions=fd_cost),
                iterations=central.iterations, evaluations=central.evaluations)


# ---- 4. Born term and permittivity at fixed states ----
def born_document(form, born):
    """Figiel water/NaCl selection with one permittivity form and one Born form (matched or source SSM+DS)."""
    base = {'ion-fraction-suppression': 'ion-fraction-suppression', 'mole-fraction-mixing': 'ion-fraction-suppression',
            'schreckenberg': 'density-temperature', 'perturbation': 'dipolar-perturbation'}[form]
    d = json.loads((SELECTIONS / f'{base}.json').read_text())
    if form == 'mole-fraction-mixing':  # Ascani Eq. 13 = sum x_i eps_i for one solvent; ions eps = 8 (FeOs/Clapeyron linear rule)
        next(f for f in d['model_families'] if f['kind'] == 'permittivity')['choice'] = 'component-permittivity-mixing'
        d['model_coefficients'] = [m for m in d['model_coefficients'] if m['family'] != 'ion_fraction_suppression_coefficient']
        for c in d['components']:
            if c['component_id'] != 'water':
                c['coefficients'].append(dict(identity=f"component/{c['component_id']}/relative_permittivity", family='relative_permittivity', value=dict(magnitude=8.0, unit='dimensionless')))
    if born == 'direct-hs-diameter':  # FeOs advanced Born: direct form at the ionic hard-sphere diameter 0.88 sigma
        next(f for f in d['model_families'] if f['kind'] == 'electrolyte').update(c_shell=0.0, c_dielectric=0.0)
        for c in d['components']:
            if c['component_id'] != 'water':
                coefficient(c, 'born_diameter')['value']['magnitude'] = coefficient(c, 'packing_diameter')['value']['magnitude']
    return d


def born_states():
    water_mw = 0.0180153  # kg/mol, the selections' water molar mass
    model = Mixture(Parameters.from_mapping(born_document('ion-fraction-suppression', 'ssm-ds'), components=NACL)); states = []
    for m in MOLALITIES:
        x = [1.0, m * water_mw, m * water_mw]; x = [v / sum(x) for v in x]
        s = model.state(298.15, P=1e5, x=x, phase='liquid')
        states.append(dict(m=m, T=298.15, x=x, rho=s.molar_density, density_source='Engine liquid root at 1e5 Pa, ion-fraction-suppression + SSM+DS'))
    rows = []
    for form in ('ion-fraction-suppression', 'mole-fraction-mixing', 'schreckenberg', 'perturbation'):
        for born in ('ssm-ds', 'direct-hs-diameter'):
            model = Mixture(Parameters.from_mapping(born_document(form, born), components=NACL))
            for st in states:
                s = model.state(st['T'], rho=st['rho'], x=st['x'])
                rows.append(dict(form=form, born=born, m=st['m'], rho=st['rho'], x=st['x'], permittivity=s.bulk_relative_permittivity,
                                 a_born=s.born, a_debye_huckel=s.debye_huckel, a_residual=s.residual_helmholtz))
    return states, rows


def born_comparator_inputs():
    d = born_document('perturbation', 'direct-hs-diameter'); s = json.loads((SELECTIONS / 'density-temperature.json').read_text())
    comps = {c['component_id']: c for c in d['components']}; sc = {c['component_id']: c for c in s['components']}
    def v(c, f): return coefficient(c, f)['value']['magnitude']
    return dict(components=[dict(id=i, z=comps[i]['fixed']['charge_number']['value']['magnitude'], m=v(comps[i], 'segment_count'),
                                 sigma=2.7927 if i == 'water' else v(comps[i], 'segment_diameter'), epsilon=v(comps[i], 'dispersion_energy_over_k'),
                                 born_diameter=None if i == 'water' else v(comps[i], 'packing_diameter'),
                                 dipole_scaling=v(comps[i], 'dipole_scaling'), polarizability_scaling=v(comps[i], 'polarizability_scaling'),
                                 correlation_integral_parameter=v(comps[i], 'correlation_integral_parameter'),
                                 relative_permittivity=78.09 if i == 'water' else 8.0,
                                 dielectric_volume_m3_per_mol=v(sc[i], 'dielectric_volume') * 1e-3 if i == 'water' else None,
                                 dielectric_temperature=v(sc[i], 'dielectric_temperature') if i == 'water' else None) for i in NACL])


def main():
    out = Path(sys.argv[1]); inputs = json.loads((CC / 'inputs/cases.json').read_text()); started = time.time()
    result = dict(provenance=[str(x) for x in provenance()], clapeyron_khudaida=clapeyron_khudaida_spec())
    t = time.time(); result['khudaida'] = {v: khudaida(v) for v in ('matched', 'constant-permittivity', 'packet')}; result['khudaida_wall_s'] = time.time() - t
    t = time.time(); result['sensitivities'] = {i: sensitivity_case(inputs, i) for i in ('C2-0.5', 'C4-298.15')}; result['sensitivities_wall_s'] = time.time() - t
    t = time.time(); result['born_states'], result['born'] = born_states(); result['born_comparators'] = born_comparator_inputs(); result['born_wall_s'] = time.time() - t
    result['wall_s'] = time.time() - started
    out.write_text(json.dumps(result, indent=1, default=str))


def declarations():
    return (
        [
            Mixture(Parameters.from_mapping(khudaida_document(v), components=KH_IDS))
            for v in ('matched', 'constant-permittivity', 'packet')
        ],
        [
            Mixture(Parameters.from_mapping(born_document(f, b), components=NACL))
            for f in ('ion-fraction-suppression', 'mole-fraction-mixing', 'schreckenberg', 'perturbation')
            for b in ('ssm-ds', 'direct-hs-diameter')
        ],
        clapeyron_khudaida_spec(),
        born_comparator_inputs(),
    )


if __name__ == '__main__':
    main()
