"""Declare the accepted cases; all thermodynamic calculations run in the Engine."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from decimal import Decimal
import time
import sys
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

LEAF = Path(__file__).resolve().parents[1]
INPUT = LEAF / 'inputs/current'
OUTPUT = LEAF / 'results/current'
CRITERIA = json.loads((INPUT / 'criteria.json').read_text())


def declarations():
    return {name: epcsaft.Mixture(epcsaft.Parameters.from_json(LEAF.parents[2] / spec['path'], components=spec['components'])) for name, spec in json.loads((INPUT / 'models.json').read_text()).items()}


def plain(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    if hasattr(type(value), '__members__'):
        return str(value)
    return {k: plain(getattr(value, k)) for k in dir(value) if not k.startswith('_') and not callable(getattr(value, k))}


def state(model, T, P, x, phase, rho=None):
    s = model.state(T, P=P, x=x, phase=phase) if rho is None else model.state(T, rho=rho, x=x)
    return {k: plain(getattr(s, k)) for k in ['pressure', 'molar_density', 'log_fugacity_coefficient', 'bulk_relative_permittivity']}


def checked(model, result, feed):
    flat, n, P = result.mole_fractions, result.amounts, result.pressure
    x = [flat[p*len(feed):(p+1)*len(feed)] for p in range(result.phase_count)]
    if not n or not all(math.isfinite(v) for v in flat+n+result.molar_densities) or len(flat) != len(n)*len(feed):
        return {'usable': False, 'cause': 'native result has no complete finite phase inventory'}
    checks = {'normalization': max(abs(sum(v)-1) for v in x),
              'charge': max(abs(v[-2]-v[-1]) for v in x),
              'balance': max(abs(sum(n[p]*x[p][i] for p in range(len(n)))-v) for i, v in enumerate(feed)),
              'minimum_amount': min(n), 'feed_total': sum(feed)}
    states, errors = [], []
    for v, rho in zip(x, result.molar_densities):
        try:
            states.append(state(model, result.temperature, P, v, 'vapor' if v[0] == 1 else 'liquid', rho))
        except Exception as error:
            states.append(None)
            errors.append(str(error))
    checks['states'], checks['state_errors'] = states, errors
    checks['pressure_relative'] = max((abs(s['pressure']-P)/max(P, 1.) for s in states if s), default=None)
    if len(states) == 2 and all(states):
        mu = [[math.log(v[i])+s['log_fugacity_coefficient'][i] if v[i] > 0 else None for i in range(len(v))] for v, s in zip(x, states)]
        checks['water_gap_RT'] = mu[0][0]-mu[1][0]
        if len(x[0]) == 4:
            checks['butanol_gap_RT'] = mu[0][1]-mu[1][1]
            checks['neutral_salt_gap_RT'] = sum(mu[0][-2:])-sum(mu[1][-2:])
    checks['usable'] = bool(result.requested_tolerance_met and checks['normalization'] <= CRITERIA['normalization'] and checks['charge'] <= CRITERIA['charge'] and checks['balance'] <= CRITERIA['balance']*max(1., sum(feed)) and min(n) >= 0 and len(errors) == 0 and checks['pressure_relative'] <= CRITERIA['pressure_relative'] and all(abs(checks[k]) <= CRITERIA['chemical_potential'] for k in ['water_gap_RT', 'butanol_gap_RT', 'neutral_salt_gap_RT'] if k in checks))
    return checks


def run(case):
    started = time.monotonic()
    ids = ('water', '1-butanol', 'lithium-cation', 'chloride-anion') if case['model'].startswith('licl') else ('water', 'sodium-cation', 'chloride-anion')
    model = MODELS[case['model']]
    output = {'case_id': case['id'], 'limitation': case.get('limitation'), 'engine': plain(epcsaft.provenance())}
    try:
        if case['kind'] == 'state':
            output['state'] = state(model, case['T'], case['P'], case['x'], case['phase'])
            try:
                packing_problem = eq.Problem(T=case['T'], P=case['P'], feed=eq.Amounts(dict(zip(ids,case['x']))), phases=[eq.Phase('source', kind=case['phase'], support=ids)], row_tolerance=CRITERIA['row'])
                output['packing_result'] = plain(eq.solve_equilibrium(model, packing_problem))
            except Exception as error:
                output['packing_request_error'] = f'{type(error).__name__}: {error}'
            return output
        if 'molality' in case:
            feed = eq.Molalities({ids[-2]: case['molality'], ids[-1]: case['molality']}, solvent=('water',), solvent_mass_kg=1.)
            inventory = [1/.0180153, case['molality'], case['molality']]
        else:
            inventory = case['feed']
            feed = eq.Amounts(dict(zip(ids, inventory)))
        if case['kind'] in ('bubble', 'charged'):
            problem = eq.bubble_point(liquid=eq.Phase('L', kind='liquid'), T=case['T'], feed=feed)
            problem.phases[1].support = case['vapor_support']
            output['vapor_support'] = case['vapor_support']
        else:
            phases = eq.Discover() if case['kind'] == 'discover' else [eq.Phase(p, kind='liquid', composition_guess=v) for p, v in zip(['aqueous', 'organic'], case['endpoints'])]
            problem = eq.Problem(T=case['T'], P=case['P'], feed=feed, phases=phases)
        problem.row_tolerance = CRITERIA['row']
        result = eq.solve_equilibrium(model, problem)
        output['result'] = plain(result)
        output['checks'] = checked(model, result, inventory)
        if not output['checks']['usable'] and case['kind'] != 'charged':
            limits = eq.EquilibriumWorkLimits(attempts=3, iterations=600, attempt_iterations=300, seconds=180., attempt_seconds=90.)
            starts = [eq.EquilibriumStart(), eq.EquilibriumStart(compositions=case.get('endpoints', []))]
            if len(ids) == 4:
                seed = json.loads((INPUT / 'diagnostic-start.json').read_text())
                starts.append(eq.EquilibriumStart(compositions=seed['compositions'], packing_fractions=seed['packing_fractions']))
            elif result.coordinates:
                starts.append(eq.EquilibriumStart(previous=result))
            probe_problem = eq.Problem(T=case['T'], P=case['P'], feed=feed, phases=[eq.Phase(p.name, kind=p.kind) for p in problem.phases], row_tolerance=problem.row_tolerance) if case['kind']=='split' else problem
            probe = eq.solve_equilibrium_starts(model, probe_problem, starts, limits)
            output['diagnosis_starts'] = plain(probe)
            if probe.selected is not None:
                result = probe.results[probe.selected]
                output['corrected_result'], output['corrected_checks'] = plain(result), checked(model, result, inventory)
        usable = output.get('corrected_checks', output['checks'])['usable']
        if usable and case['kind'] != 'discover':
            output['stability'] = plain(eq.stability_test(eq.compile_problem(model, problem), result, eq.StabilityOptions()))
        if case.get('refine'):
            problem.row_tolerance = CRITERIA['refined_row']
            refined = eq.solve_equilibrium(model, problem, result)
            output['refined_result'], output['refined_checks'] = plain(refined), checked(model, refined, inventory)
        if case['kind'] == 'split':
            discovery = eq.Problem(T=case['T'], P=case['P'], feed=feed, phases=eq.Discover(), row_tolerance=CRITERIA['row'])
            discovered = eq.solve_equilibrium(model, discovery)
            output['discovered_result'], output['discovered_checks'] = plain(discovered), checked(model, discovered, inventory)
    except Exception as error:
        output['error'] = f'{type(error).__name__}: {error}'
    finally:
        output['wall_seconds'] = time.monotonic()-started
    return output


def execution_inputs(case):
    spec = json.loads((INPUT / 'models.json').read_text())[case['model']]
    files = [INPUT / 'criteria.json', LEAF.parents[2] / spec['path']]
    if case['model'] == 'licl-conditional' and case['kind'] != 'state':
        files.append(INPUT / 'diagnostic-start.json')
    hashes = {str(p.relative_to(LEAF.parents[2])): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    hashes.update({name: hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest() for name, value in [('declaration', case), ('model_selection', spec)]})
    return hashes


def campaign_identity(cases):
    identities = {}
    for case in cases:
        with gzip.open(OUTPUT / (case['id']+'.json.gz'), 'rt') as stream:
            record = json.load(stream)
        identity = record.get('execution_identity', {})
        if not identity or identity['inputs'] != execution_inputs(case):
            raise ValueError(f"execution-input mismatch for {case['id']}; recalculate this case before comparing")
        identities[case['id']] = dict(identity, engine=record['engine'])
    return {'cases': identities}


def comparisons(cases):
    summaries, endpoints, scalars = [], [], []
    sources = list(csv.DictReader((LEAF / 'inputs/source/tables-3-5-source.csv').open()))
    for c in cases:
        with gzip.open(OUTPUT / (c['id']+'.json.gz'), 'rt') as stream:
            d = json.load(stream)
        r, checks = d.get('corrected_result', d.get('result', {})), d.get('corrected_checks', d.get('checks', {}))
        row = dict(case_id=c['id'], kind=c['kind'], model=c['model'], limitation=c.get('limitation'), wall_seconds=d['wall_seconds'], usable=checks.get('usable'), phase_count=r.get('phase_count'), error=d.get('error', r.get('message')), pressure_Pa=r.get('pressure'), source_pressure_Pa=c.get('source_pressure'), molality=c.get('molality'))
        row.update({k: checks.get(k) for k in ['normalization', 'charge', 'balance', 'pressure_relative', 'water_gap_RT', 'butanol_gap_RT', 'neutral_salt_gap_RT']})
        if c['kind'] == 'bubble' and checks.get('usable'):
            row.update(pressure_error_Pa=r['pressure']-c['source_pressure'], display_error_Pa=r['pressure']-round(c['source_pressure'], 2), source_half_unit_Pa=c['half_unit'])
        n = len(c.get('feed', [0, 0, 0]))
        ordered = lambda v: sorted([v[i:i+n] for i in range(0, len(v), n)], key=lambda x: -x[0])
        if 'refined_result' in d:
            f = d['refined_result']
            row.update(refined_usable=d['refined_checks']['usable'], refinement_same_phases=f['phase_count']==r['phase_count'], refinement_pressure_Pa=abs(f['pressure']-r['pressure']), refinement_x=max(abs(a-b) for a,b in zip(sum(ordered(f['mole_fractions']), []), sum(ordered(r['mole_fractions']), []))), refinement_pressure_limit_Pa=c.get('refinement_limit'))
        if 'discovered_result' in d:
            v, vc = d['discovered_result'], d['discovered_checks']
            row.update(discovered_usable=vc['usable'], discovered_phase_count=v['phase_count'])
            if checks.get('usable') and vc['usable'] and r['phase_count'] == v['phase_count']:
                row.update(objective_difference_per_feed=abs(r['objective']-v['objective'])/sum(c['feed']), declared_discovered_x=max(abs(a-b) for a,b in zip(sum(ordered(r['mole_fractions']), []),sum(ordered(v['mole_fractions']), []))))
            for method, result in [('declared', r), ('discovered', v)]:
                for phase, x in zip(['aqueous', 'organic'] if result['phase_count']==2 else ['single'], ordered(result['mole_fractions'])):
                    mass = [x[0]*.01801528, x[1]*.0741216, x[2]*.00694+x[3]*.03545]
                    reference = c['endpoints'][0 if phase=='aqueous' else 1] if phase!='single' else [math.nan]*n
                    endpoints.append(dict(source_x_error_max=None if phase=='single' else max(abs(a-b) for a,b in zip(x, reference)), case_id=c['id'], method=method, phase=phase, water_x=x[0], butanol_x=x[1], lithium_x=x[2], chloride_x=x[3], water_w=mass[0]/sum(mass), butanol_w=mass[1]/sum(mass), licl_w=mass[2]/sum(mass), limitation=c['limitation']))
        if 'stability' in d:
            row.update(stability_completed=d['stability']['completed'], stability_message=d['stability']['message'], resolved_descent=any(t['resolved_descent'] for t in d['stability']['trials']), trial_count=len(d['stability']['trials']))
        if c['kind']=='state' and 'state' in d:
            row['state_pressure_error_Pa']=d['state']['pressure']-c['P']
        summaries.append(row)
        if c['kind']=='split' and c['id'].endswith('101300') and checks.get('usable'):
            order=sorted(range(r['phase_count']), key=lambda p:-r['mole_fractions'][p*n])
            for phase, p in zip(['aqueous','organic'],order):
                source=[s for s in sources if s['case_id']==c['source_case'] and s['phase']==phase]
                xv=[s for s in source if s['quantity'].startswith('x_')]
                for s in source:
                    value=float(s['printed_value']); half=0.5*10**Decimal(s['printed_value']).as_tuple().exponent
                    if s in xv:
                        vals=[float(t['printed_value']) for t in xv]; hs=[.5*10**Decimal(t['printed_value']).as_tuple().exponent for t in xv]; lo=(value-half)/(sum(vals)+sum(hs)-2*half); hi=(value+half)/(sum(vals)-sum(hs)+2*half)
                        actual=r['mole_fractions'][p*n+{'x_water':0,'x_1_butanol':1,'x_lithium':2,'x_chloride':3}[s['quantity']]]
                    else:
                        lo,hi=value-half,value+half; actual=r['packing_fractions'][p] if s['quantity']=='packing_fraction' else 1/r['molar_densities'][p]
                    scalars.append(dict(case_id=c['id'],phase=phase,quantity=s['quantity'],source_printed=s['printed_value'],source_lower=lo,source_upper=hi,current=actual,error=actual-value,within_source_rounding=lo<=actual<=hi,limitation=c['limitation']))
    for name, rows in [('summary',summaries),('endpoints',endpoints),('table5-comparison',scalars)]:
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (OUTPUT / (name+'.csv')).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader();writer.writerows(rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--read-results', action='store_true', help='compare retained results without solving')
    parser.add_argument('--case', help='run just this exact retained case ID')
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    cases = json.loads((INPUT / 'cases.json').read_text())
    if args.case and args.case not in {case['id'] for case in cases}:
        parser.error(f'unknown case: {args.case}')
    if not args.read_results:
        MODELS = declarations()
    for case in ([] if args.read_results else cases):
        if args.case and case['id'] != args.case:
            continue
        identity = dict(producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), python=sys.executable, inputs=execution_inputs(case))
        record = run(case)
        record['execution_identity'] = identity
        with gzip.open(OUTPUT / (case['id']+'.json.gz'), 'wt') as stream:
            json.dump(record, stream, allow_nan=False, separators=(',', ':'))
        print(case['id'], record['wall_seconds'], record.get('error', record.get('checks', {}).get('usable')), flush=True)

    identity = campaign_identity(cases)
    comparisons(cases)
    (OUTPUT / 'run-identity.json').write_text(json.dumps(identity, indent=2)+'\n')
