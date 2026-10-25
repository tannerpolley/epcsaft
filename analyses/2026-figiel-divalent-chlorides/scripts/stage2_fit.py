"""#209 binary divalent fits on frozen /1; all thermodynamics and optimization are native."""
import argparse
import copy
import hashlib
import itertools
import json
import math
import re
import subprocess
import time
from pathlib import Path

from epcsaft import Mixture, Parameters, _core, regression as r
import stage1_fit as F
import stage1_minimax as M

SOURCE = F.REPO / 'data/packets/figiel-2025-chloride-activity-refit/1/parameters/parameters.json'
HELD = F.REPO / 'data/packets/held-2014-strategy2/1/parameters/parameters.json'
DATA = F.REPO / 'analyses/2026-held-2014-chlorides/inputs/evaluated.csv'
BASE = F.declared(SOURCE)
HELD_DOC = json.loads(HELD.read_text())
SALTS = {'MgCl2': ('magnesium-cation', -1812100., (3.1327, -.25, .817)),
         'CaCl2': ('calcium-cation', -1487100., (3.2648, .0041, 1.))}
FAMILIES = ('segment_diameter', 'packing_diameter', 'debye_huckel_diameter')
OUT = F.OUT


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, record):
    arrays = []
    for i, result in enumerate(record.get('results', [])):
        for owner, keys in ((result, ('iteration_seconds', 'iteration_costs')),
                            (result['covariance'], ('left_vectors', 'right_vectors', 'optimizer', 'physical'))):
            for key in keys:
                values = owner.get(key, [])
                arrays += [dict(start=i, quantity=('covariance/' if owner is result['covariance'] else '')+key,
                                element=j, value=value) for j,value in enumerate(values)]
                if values:
                    owner[key] = dict(csv=f'stage2-{name}-numeric-arrays.csv', elements=len(values))
    if arrays:
        write(f'{name}-numeric-arrays', arrays)
    (OUT / f'stage2-{name}.json').write_text(json.dumps(M.finite(record), indent=2, allow_nan=False)+'\n')


def write(name, rows):
    prior = F.STEM
    F.STEM = 'stage2'
    F.write(name, rows)
    F.STEM = prior


def extended():
    doc = copy.deepcopy(BASE)
    for ion, _, _ in SALTS.values():
        component = copy.deepcopy(next(c for c in HELD_DOC['components'] if c['component_id'] == ion))
        for family, value in (('born_diameter', 3.5), ('solvation_factor', 1.)):
            component['coefficients'].append({'identity': f'component/{ion}/{family}', 'family': family,
                'value': {'magnitude': value, 'unit': 'angstrom' if family == 'born_diameter' else 'dimensionless'}})
        doc['components'].append(component)
        for partner in ('water', 'chloride-anion'):
            a, b = sorted((ion, partner))
            doc['pairs'].append({'component_id_a': a, 'component_id_b': b, 'coefficients': [
                {'identity': f'pair/{a}/{b}/k_ij', 'family': 'k_ij', 'value': {'magnitude': 0., 'unit': 'dimensionless'}}]})
    return doc


def declarations(salt, document=None):
    ion, target, held_start = SALTS[salt]
    p = Parameters.from_mapping(extended() if document is None else document, components=('water', ion, 'chloride-anion'))
    cs = [r.coordinate(p, [('segment_diameter', ion, 1), ('packing_diameter', ion, .88),
                           ('debye_huckel_diameter', ion, .88)], origin=0, scale=1, bounds=(2.4, 4.4))]
    cs += [r.coordinate(p, [('k_ij', (ion, partner), 1)], origin=0, scale=1, bounds=(-1, 1))
           for partner in ('water', 'chloride-anion')]
    if document is None:
        chart = r.evaluate_coordinates(p, cs, held_start, [])
        p = Parameters(chart.parameters, p.component_ids, p.site_ids)
    energy = r.electrolyte(p, 'solvation_gibbs_energy', F.T, F.P, {'water': 1}, ions=(ion, 'chloride-anion'),
        stoichiometry=(1, 2), molality=0, species=ion, observed=target, form='difference', scale=1)
    hs = [r.hydration_hold(energy, bracket=(1, 8), energy_tolerance=1e-6,
                           diameter_tolerance=1e-10, maximum_iterations=50)]
    records = [v for v in F.read(DATA) if v['salt'] == salt and float(v['molality_mol_kg']) <= 4.5]
    assert len(records) == 42 and all(float(v['temperature_K']) == F.T for v in records)
    obs, meta = [], []
    for row in records:
        for quantity, kind in (('gamma', 'mean_ionic_activity_coefficient'), ('osmotic', 'osmotic_coefficient')):
            ob = r.electrolyte(p, kind, F.T, F.P, {'water': 1}, ions=(ion, 'chloride-anion'), stoichiometry=(1, 2),
                molality=float(row['molality_mol_kg']), observed=float(row[f'{quantity}_evaluated']), form='log_ratio')
            obs.append(ob)
            meta.append(dict(row=len(meta), salt=salt, quantity=quantity, molality=ob.molality, observed=ob.observed))
    obs.append(energy)
    checks = [F.group(quantity, list(range(j, 84, 2)), r.CheckMetric.AbsoluteLog10Ratio, .02, report=True)
              for j, quantity in enumerate(('gamma', 'osmotic'))]
    checks += [F.group('hydration', [84], r.CheckMetric.AbsoluteDifference, 10)]
    return p, cs, hs, obs, meta, checks


def native_start(p, cs, physical):
    point = r.evaluate_coordinates(p, cs, physical, [])
    chart = r.CoordinateSet()
    chart.coordinates = cs
    return chart.start(point.parameters)


def candidate(salt, physical, document=None):
    doc = extended() if document is None else copy.deepcopy(document)
    ion = SALTS[salt][0]
    values = dict(zip((*FAMILIES, 'water', 'chloride-anion', 'born_diameter'), physical))
    for component in doc['components']:
        if component['component_id'] == ion:
            for c in component['coefficients']:
                if c['family'] in values:
                    c['value']['magnitude'] = values[c['family']]
    for pair in doc['pairs']:
        names = {pair['component_id_a'], pair['component_id_b']}
        if ion in names:
            partner = next(v for v in names if v != ion)
            pair['coefficients'][0]['value']['magnitude'] = values[partner]
    return doc


def guards(document, label):
    assert document['components'][:len(BASE['components'])] == BASE['components']
    assert document['pairs'][:len(BASE['pairs'])] == BASE['pairs']
    assert all(document[key] == BASE[key] for key in ('model_coefficients', 'model_families', 'correlations', 'topology'))
    source, fixed = M.assembled(), M.assembled(document=BASE)
    stored = M.assembled(document=document)
    _, _, obs, meta, groups = F.declarations(source, 'all')
    old = _core.observations_at(fixed._native, [], [], [], obs, [0.]*len(obs))
    new = _core.observations_at(stored._native, [], [], [], obs, [0.]*len(obs))
    checks = r.evaluate_checks(stored, obs, groups, source)
    records = [v | dict(reference=a.observed, stored=b.observed, error=abs(b.observed-a.observed),
        tolerance=1e-12*abs(a.observed) if a.observed else 1e-12,
        passed=abs(b.observed-a.observed) <= (1e-12*abs(a.observed) if a.observed else 1e-12))
        for v, a, b in zip(meta, old, new)]
    if label == 'combined':
        write('combined-stage1-preservation', records)
    group_records = F.checks_record(checks, label)
    if label in ('baseline', 'combined'):
        write(f'{label}-stage1-groups', group_records)
    if label != 'baseline':
        operands = lambda rows: [{k:str(v) if v is not None else '' for k,v in row.items() if k!='start'} for row in rows]
        assert operands(group_records) == operands(F.read(OUT/'stage2-baseline-stage1-groups.csv'))
    invariants = []
    for solvent in ('water', 'methanol', 'ethanol'):
        x = [float(v==solvent) for v in stored.component_ids]
        a,b = [Mixture(p).state(F.T,P=F.P,x=x,phase='liquid') for p in (fixed,stored)]
        i = stored.component_ids.index(solvent)
        invariants += [dict(solvent=solvent, quantity=q, reference=av, stored=bv, error=abs(bv-av)) for q,av,bv in (
            ('density_mol_m3',a.molar_density,b.molar_density), ('ln_phi',a.log_fugacity_coefficient[i],b.log_fugacity_coefficient[i]))]
    if label == 'combined':
        write('combined-stage1-invariance', invariants)
    assert all(v['error'] <= (1e-12*abs(v['reference']) if v['reference'] else 1e-12) for v in invariants)
    assert all(v['passed'] for v in records) and all(c.status == 'pass' for c in checks if not c.report_only)
    return dict(rows=len(records), maximum_change=max(v['error'] for v in records),
                enforced_groups=sum(not c.report_only for c in checks), passed=True)


def retain(salt, phase, results, meta, obs):
    summaries, values, columns, constraints = [], [], [], []
    for j, result in enumerate(results):
        record = {key: getattr(result, key) for key in ('status', 'message', 'usable', 'optimizer', 'physical', 'iterations',
            'initial_cost', 'final_cost', 'residual_evaluations', 'jacobian_evaluations', 'trial_failures', 'iteration_seconds', 'iteration_costs')}
        record['holds'] = [{key: getattr(h, key) for key in ('energy_error', 'correction', 'uncertainty', 'derivative', 'parameter_derivative')}
                           for h in result.chart.holds]
        record['covariance'] = {key: getattr(result.covariance, key) for key in ('status', 'rank', 'singular_values', 'active_bounds',
            'optimizer', 'physical', 'variance_factor', 'projected_gradient_norm', 'left_vectors', 'right_vectors', 'assumption')}
        record['failures'] = [{key: (f.status.name if key == 'status' else getattr(f, key)) for key in ('trial', 'observation', 'optimizer', 'status', 'message', 'jacobian')}
                              for f in result.failures]
        for v, y in zip(meta, result.predictions):
            values.append(v | dict(phase=phase, start=j, prediction=y, error_decade=math.log10(y/v['observed'])))
        for i, value in enumerate(result.optimizer_jacobian):
            columns.append(dict(start=j, row=i//3, coordinate=i%3, weighted_residual_column=value))
        if result.minimax:
            m = result.minimax
            record['minimax'] = {key: getattr(m, key) for key in ('availability', 'solver_status', 'feasible', 'stationary', 't',
                'target_maximum', 'primal', 'stationarity', 'dual_sign', 'complementarity', 'active_lower', 'active_upper',
                'lower_distances', 'upper_distances', 'hold_multipliers', 'hold_causes', 'full_stationarity')}
            for i, row in enumerate(m.rows):
                constraints.append(dict(start=j, constraint=i, kind=row.kind, row=row.row, sign=row.sign,
                    threshold=row.threshold, scaled_value=m.values[i], physical_value=m.physical_values[i],
                    physical_multiplier=m.physical_multipliers[i], root_error=m.root_errors[i]))
        summaries.append(record)
    for name, rows in (('rows', values), ('columns', columns), ('constraints', constraints)):
        if rows:
            write(f'{salt}-{phase}-{name}', rows)
    return summaries


def fit(salt, phase):
    p, cs, hs, obs, meta, checks = declarations(salt)
    controls = r.FitControls(maximum_iterations=200 if phase == 'ceres' else 600,
                             maximum_elapsed_time_seconds=120 if phase == 'ceres' else 900,
                             function_tolerance=1e-12, gradient_tolerance=1e-10, parameter_tolerance=1e-12)
    spec = None
    if phase == 'ceres':
        physical = [SALTS[salt][2], *itertools.product((2.9, 3.4, 3.9), (-.5, 0., .5), (-.5, 0., .5))]
        starts = [native_start(p, cs, v) for v in physical]
    else:
        previous = json.loads((OUT/f'stage2-{salt}-ceres.json').read_text())
        starts = []
        for v in sorted((v for v in previous['results'] if v['usable']), key=lambda v: v['final_cost']):
            if all(max(abs(x-y)/2 for x, y in zip(v['optimizer'], u)) > .01 for u in starts):
                starts.append(v['optimizer'])
        spec = r.MinimaxRequest()
        spec.targets, spec.group_scales = list(range(84)), [.02, .02, 1000.]
        spec.q, spec.tolerance, spec.activity_distance, spec.margin = .02, 1e-10, 1e-6, .1
    assert starts, 'No usable Ceres point; inspect retained failures before minimax.'
    write(f'{salt}-{phase}-starts', [dict(start=i, sigma_angstrom=u[0], k_water=u[1], k_chloride=u[2]) for i,u in enumerate(starts)])
    script_sha256 = digest(Path(__file__))
    print(salt, phase, len(starts), 'declared native starts', flush=True)
    began = time.perf_counter()
    results = r.fit(p, cs, obs, holds=hs, checks=checks, controls=controls, minimax=spec, singular_vectors=True,
        weights=[1/(.02*math.log(10))**2]*84+[0.], starts=starts,
        shared_iterations=5600 if phase == 'ceres' else 600, shared_seconds=1800 if phase == 'ceres' else 900)
    summaries = retain(salt, phase, results, meta, obs)
    save(f'{salt}-{phase}', dict(results=summaries, declared_starts=starts, completed_starts=len(results),
        seconds=time.perf_counter()-began, wheel_sha256=digest(next((F.REPO/'build/environment-wheel').glob('*.whl'))),
        source_sha256=digest(SOURCE), evaluated_sha256=digest(DATA), script_sha256=script_sha256,
        shared_seconds=1800 if phase=='ceres' else 900, shared_iterations=5600 if phase=='ceres' else 600,
        source_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=F.REPO,text=True).strip(), controls={key:getattr(controls,key) for key in
        ('maximum_iterations', 'maximum_elapsed_time_seconds', 'function_tolerance', 'gradient_tolerance', 'parameter_tolerance')}))
    print(salt, phase, [(v.status, v.final_cost, v.minimax.target_maximum if v.minimax else None) for v in results], flush=True)


def qualify(salt, document=None):
    result = json.loads((OUT/f'stage2-{salt}-minimax.json').read_text())
    usable = [(j, v) for j, v in enumerate(result['results']) if v['usable']]
    assert usable, 'No usable native minimax result; inspect the retained stopping point.'
    j, selected = min(usable, key=lambda pair: pair[1]['minimax']['target_maximum'])
    doc = candidate(salt, selected['physical']) if document is None else document
    save(f'{salt}-candidate', doc)
    p, cs, hs, obs, meta, checks = declarations(salt, doc)
    chart_set = r.CoordinateSet()
    chart_set.coordinates = cs
    u = chart_set.start(p._native)
    values = _core.observations_at(p._native, [], [], [], obs, [0.]*len(obs))
    native_checks = r.evaluate_checks(p, obs, checks, p)
    rows = [v | dict(prediction=y.observed, error_decade=math.log10(y.observed/v['observed'])) for v, y in zip(meta, values)]
    write(f'{salt}-stored-rows', rows)
    model, states = Mixture(p), []
    for ob in obs[::2][:-1]:
        x = _core.formula_composition(ob.composition, p._native, ob.ion_indices, (1, 2), ob.molality)
        s = model.state(F.T, P=F.P, x=x, phase='liquid')
        charge = 2*x[1]-x[2]
        states.append(dict(molality=ob.molality, density_mol_m3=s.molar_density, pressure_error_Pa=abs(s.pressure-F.P),
            pressure_density_derivative=s.pressure_density_derivative, charge=charge,
            passed=abs(s.pressure-F.P)<.01 and s.pressure_density_derivative>0 and charge==0))
    write(f'{salt}-stored-states', states)
    bands = []
    for quantity in ('gamma', 'osmotic'):
        for lo, hi in ((0, .1), (.1, 1), (1, 3), (3, 4.5)):
            band = [v for v in rows if v['quantity']==quantity and lo<v['molality']<=hi]
            worst = max(band, key=lambda v: abs(v['error_decade']))
            bands.append(dict(quantity=quantity, lower_molality=lo, upper_molality=hi, rows=len(band),
                maximum_decade=abs(worst['error_decade']), argmax_molality=worst['molality']))
    write(f'{salt}-bands', bands)
    born = r.coordinate(p, [('born_diameter', SALTS[salt][0], 1)], origin=0, scale=1, bounds=(1, 8))
    endpoints = [_core.observations_at(p._native, [born], [b], [], [obs[-1]], [0.])[0].observed for b in (1., 8.)]
    exact = _core.evaluate_observations(p._native, cs, u, hs, obs[:84])
    sigma_map = exact.chart.map[:9]
    sigma_map_passed = sigma_map == [1., 0., 0., .88, 0., 0., .88, 0., 0.]
    h = .01  # Symmetric fourth-order differences resolve dilute osmotic cancellation.
    estimates = []
    for step in (h, h/2):
        sides = []
        for offset in (-2*step, -step, step, 2*step):
            point = list(u)
            point[0] += offset
            y = _core.observations_at(p._native, cs, point, hs, obs[:84], [0.]*84)
            sides.append([math.log(v.observed/ob.observed) for v, ob in zip(y, obs)])
        estimates.append([(a-8*b+8*c-d)/(12*step) for a,b,c,d in zip(*sides)])
    columns = []
    for i, (coarse, fine) in enumerate(zip(*estimates)):
        value = exact.residual_optimizer_gradient[3*i]
        reference, uncertainty = (16*fine-coarse)/15, abs(fine-coarse)/15
        tolerance = 1e-7+2e-5*abs(value)
        columns.append(dict(row=i, exact_log_column=value, Richardson=reference, error=abs(value-reference),
            uncertainty=uncertainty, tolerance=tolerance, passed=abs(value-reference)<=tolerance and uncertainty<=tolerance/4))
    write(f'{salt}-sigma-columns', columns)
    frozen = guards(doc, salt)
    accepted = (sigma_map_passed and all(v['passed'] for v in states+columns) and max(abs(v['error_decade']) for v in rows)<=.02
        and abs(values[-1].observed-SALTS[salt][1])<=10 and selected['minimax']['feasible'] and selected['minimax']['stationary']
        and result['completed_starts']==len(result['declared_starts']))
    report = dict(selected=j, parameters=dict(sigma_angstrom=selected['physical'][0], packing_angstrom=selected['physical'][1],
        DH_angstrom=selected['physical'][2], k_water=selected['physical'][3], k_chloride=selected['physical'][4], born_angstrom=selected['physical'][5]),
        hydration_target_J_mol=SALTS[salt][1], stored_hydration_J_mol=values[-1].observed,
        bracket_endpoint_energies_J_mol=endpoints, holds=selected['holds'], minimax=selected['minimax'], sigma_map=sigma_map,
        maximum_gamma_decade=max(abs(v['error_decade']) for v in rows if v['quantity']=='gamma'),
        maximum_osmotic_decade=max(abs(v['error_decade']) for v in rows if v['quantity']=='osmotic'),
        maximum_pressure_error_Pa=max(v['pressure_error_Pa'] for v in states), stage1_preservation=frozen,
        failed_sigma_columns=sum(not v['passed'] for v in columns), accepted=accepted,
        claim_limits=['Binary aqueous fit to evaluated reference data at 298.15 K and 1 bar.',
                      'Finite declared multistart search; no global optimum, parameter uncertainty or physical validation.', 'Density is unqualified.'])
    save(f'{salt}-qualification', report)
    print(salt, 'qualification', report, flush=True)
    return report


def packet_document(doc):
    doc = copy.deepcopy(doc)
    doc['document_version'] = 2
    doc['purpose'] = 'Aqueous binary MgCl2/CaCl2 fit on frozen stage-1 /1, issue 209. Evaluated-reference fitting only; no mixed-salt, density, alcohol or parameter-uncertainty qualification.'
    sources = [HELD_DOC['sources'][0], {'source_id': 'marcus-1991', 'citation': 'Marcus, Y. Thermodynamics of solvation of ions. Part 5. Gibbs free energy of hydration at 298.15 K. J. Chem. Soc., Faraday Trans. 1991, 87, 2995-2999.',
        'doi': '10.1039/FT9918702995', 'use_basis': 'Table 1 p. 2996; G* plus the owner-verified Figiel Eq. 18 shift 17.908 kJ/mol.'},
        {'source_id': 'issue-209-stage2-native-fit', 'citation': 'ePC-SAFT issue 209 native stage-2 fit, 2026-10-05',
         'use_basis': 'Goldberg-Nuttall 1978 Tables 17/20; native hydration holds to Marcus 1991; frozen stage-1 packet /1.'}]
    doc['sources'] += [s for s in sources if s['source_id'] not in {v['source_id'] for v in doc['sources']}]
    doc['domains'] += [v for v in HELD_DOC['domains'] if v['domain_id'] not in {d['domain_id'] for d in doc['domains']}]
    ions = {v[0] for v in SALTS.values()}
    for owner in doc['components']+doc['pairs']:
        new = owner.get('component_id') in ions or bool(ions & {owner.get('component_id_a'), owner.get('component_id_b')})
        if not new:
            continue
        for c in owner['coefficients']:
            if c['family'] in (*FAMILIES, 'k_ij', 'born_diameter', 'solvation_factor'):
                c['provenance'] = {'source_id': 'figiel-yu-held-2025' if c['family']=='solvation_factor' else 'issue-209-stage2-native-fit',
                    'locator': 'Figiel solvation factor 1' if c['family']=='solvation_factor' else 'stage2-summary.json; sigma/packing/DH affine chart and Marcus hydration hold', 'domain_id':'figiel-298-k-1-bar'}
    return doc


def publish(doc):
    root = SOURCE.parents[1].parent/'2'
    assert not root.exists(), 'Immutable packet /2 already exists.'
    packet = root/'parameters/parameters.json'
    packet.parent.mkdir(parents=True)
    packet.write_text(json.dumps(doc, indent=2)+'\n')
    origin = root/'provenance/source-origin.toml'
    origin.parent.mkdir()
    text = 'provenance_kind = "derived-record-provenance-v1"\nschema_version = 1\ntransformations = ["Native binary Mg/Ca fits and Marcus hydration holds on frozen /1; preserve every stage-1 parameter."]\n'
    text += '\n[[records]]\nretained_path = "parameters/parameters.json"\nsha256 = '+json.dumps(digest(packet))+'\n'
    for i, path in enumerate((SOURCE, HELD, DATA, OUT/'stage2-summary.json', OUT/'stage2-MgCl2-stored-rows.csv', OUT/'stage2-CaCl2-stored-rows.csv')):
        text += f'\n[[sources]]\nsource_id = "input-{i}"\nrole = "Retained fit input or result"\nretained_path = '+json.dumps(str(path.relative_to(F.REPO)))+'\nsha256 = '+json.dumps(digest(path))+'\n'
    origin.write_text(text)
    manifest = (SOURCE.parents[1]/'packet.toml').read_text().split('[[files]]')[0]
    updates = {'packet_version':'2', 'title':'Figiel aqueous chloride stage-2 Mg/Ca fit', 'description':doc['purpose'],
               'tolerance_rationale':'All 168 binary Mg/Ca rows meet 0.02 decade; stored hydration within 10 J/mol; frozen stage-1 rows unchanged within 1e-12. Numerical fitting, not physical validation.'}
    for key, value in updates.items():
        manifest = re.sub(r'^'+key+r' = .*$', key+' = '+json.dumps(value), manifest, flags=re.MULTILINE)
    for source in doc['sources']:
        if source['source_id'] not in manifest:
            manifest += '\n[[sources]]\n'+'\n'.join(k+' = '+json.dumps(v) for k,v in source.items())+'\n'
    for path, role, media in ((packet,'parameters','application/json'), (origin,'provenance','application/toml')):
        manifest += '\n[[files]]\npath = '+json.dumps(str(path.relative_to(root)))+'\nrole = '+json.dumps(role)+'\nmedia_type = '+json.dumps(media)+'\nsha256 = '+json.dumps(digest(path))+'\n'
        manifest += 'source_refs = '+json.dumps([v['source_id'] for v in doc['sources']])+'\ndomain_refs = ["figiel-298-k-1-bar"]\nunits_basis = "Explicit coefficient units in the parameter document."\n'
    (root/'packet.toml').write_text(manifest)


def finish():
    summaries = {salt: json.loads((OUT/f'stage2-{salt}-qualification.json').read_text()) for salt in SALTS}
    accepted = all(v['accepted'] for v in summaries.values())
    doc = extended()
    for salt in SALTS:
        selected = json.loads((OUT/f'stage2-{salt}-minimax.json').read_text())['results'][summaries[salt]['selected']]
        doc = candidate(salt, selected['physical'], doc)
    doc = packet_document(doc)
    save('combined-candidate', doc)
    doc = json.loads((OUT/'stage2-combined-candidate.json').read_text())
    if accepted:
        summaries = {salt: qualify(salt, doc) for salt in SALTS}
        accepted = all(v['accepted'] for v in summaries.values())
    frozen = guards(doc, 'combined')
    save('summary', dict(salts=summaries, accepted=accepted, stage1_preservation=frozen, source_sha256=digest(SOURCE),
        wheel_sha256=digest(next((F.REPO/'build/environment-wheel').glob('*.whl'))), packet_written=accepted,
        input_gate=dict(source='Marcus 1991 Table 1 p. 2996; Zotero 9J34P8R4/T2JYYU7D', G_star_kJ_mol=[-1830,-1505],
            conversion_kJ_mol=17.908, mole_fraction_targets_kJ_mol=[-1812.1,-1487.1]),
        accounting=dict(engine_runtime_added=0, analysis_gross_ceiling=450)))
    if not accepted:
        print('Acceptance failed: retained bands, bounds, hydration bracket, columns and closure; packet /2 not written.', flush=True)
        return
    publish(doc)


def diagnose():
    p,cs,hs,obs,_,_ = declarations('CaCl2')
    chart = r.evaluate_coordinates(p,cs,native_start(p,cs,(2.9,.5,-.5)),hs)
    model = Mixture(Parameters(chart.parameters,p.component_ids,p.site_ids))
    bad = obs[81]
    bad.form, bad.scale = r.ResidualForm.Difference, 1.
    prediction = r.predict(model,bad)
    save('diagnosis', dict(ca_initial_start=dict(sigma_angstrom=2.9,k_water=.5,k_chloride=-.5,
        molality=bad.molality,raw_osmotic=prediction.value,status=prediction.status.name),
        numerical_fix_attempt='28 declared Ceres starts per salt followed by native minimax from each distinct usable minimum; native sigma-column, hydration and pressure checks.',
        conclusion='The accepted finite search returns stationary bounded points outside the row criterion. Examine retained bands and active bounds; this is not a global impossibility proof.',
        coordinate_or_tolerance_change=False,packet_written=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('preflight', 'ceres', 'minimax', 'qualify', 'finish', 'diagnose'))
    parser.add_argument('--salt', choices=tuple(SALTS))
    args = parser.parse_args()
    if args.phase == 'preflight':
        guards(BASE, 'baseline')
        for salt in SALTS:
            p, cs, hs, obs, _, _ = declarations(salt)
            point = native_start(p, cs, SALTS[salt][2])
            chart = r.evaluate_coordinates(p, cs, point, hs)
            print(salt, 'Held start', chart.physical, [(h.energy_error, h.correction) for h in chart.holds], flush=True)
    elif args.phase in ('ceres', 'minimax'):
        fit(args.salt, args.phase)
    elif args.phase == 'qualify':
        qualify(args.salt)
    elif args.phase == 'diagnose':
        diagnose()
    else:
        finish()
