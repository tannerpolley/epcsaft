"""Evaluate the frozen #640 falsifiers; no refit, packet adoption, or stage B."""
from __future__ import annotations

from functools import lru_cache
import json
import math
import os
import time

from epcsaft import Mixture, PropertyActionStatus

import fit
from trace_activity import Activities, COMPONENTS as TRACE_COMPONENTS

HERE, ROOT, D, T, P = fit.HERE, fit.ROOT, fit.D, fit.T, fit.P
RESULTS = HERE / 'results'
ION = {'Li': ('lithium-cation', 1), 'Na': ('sodium-cation', 1), 'K': ('potassium-cation', 1),
       'Mg': ('magnesium-cation', 2), 'Ca': ('calcium-cation', 2), 'Cl': ('chloride-anion', -1)}


def metrics(values):
    assert values and all(math.isfinite(v) for v in values)
    return {'rows': len(values), 'bias': sum(values)/len(values),
            'RMS': math.sqrt(sum(v*v for v in values)/len(values)), 'maximum_absolute_error': max(map(abs, values))}


def passes(m, limits):
    return m['RMS'] <= limits['RMS'] and m['maximum_absolute_error'] <= limits['maximum_absolute_error']


class Water:
    def __init__(self, p):
        self.p = p
        self.reference = Mixture(p).reference(T, P, solvent_fractions=[float(c == 'water') for c in p.component_ids])

    @lru_cache(maxsize=None)
    def predict(self, li=0.0, na=0.0, mg=0.0, ca=0.0, k=0.0):
        x = fit.composition(self.p, {ION['Li'][0]: li, ION['Na'][0]: na, ION['Mg'][0]: mg,
                                   ION['Ca'][0]: ca, ION['K'][0]: k, ION['Cl'][0]: li+na+2*(mg+ca)+k})
        a = self.reference.solvent_activity(x=[x[c] for c in self.p.component_ids])
        assert a.status == PropertyActionStatus.Available, (a.status, a.cause, a.osmotic_cause, x)
        assert isinstance(a.osmotic_coefficient, float) and math.isfinite(a.osmotic_coefficient)
        s = a.state
        return {'phi': a.osmotic_coefficient, 'aw': a.solvent_activity, 'density_mol_m3': s.molar_density,
                'pressure_error_Pa': s.pressure-P, 'density_residual': s.density_convergence.residual,
                'density_stop': str(s.density_convergence.stop)}


def mixture_rows(source):
    rows = fit.read(ROOT / D['held_out_paths'][source])
    for i,row in enumerate(rows, 1):
        row['csv_row'] = i
        row['phi'] = None if row.get('flag') == 'illegible' else float(row['phi'])
        if source == 'Yao':
            row['a'], row['b'] = float(row['m_LiCl']), float(row['m_MgCl2'])
            role = D['row_lists']['Yao'][i-1]
            row['locator'], row['bucket'] = role['locator'], role['role']
        else:
            row['a'], row['b'] = float(row['m_KCl']), float(row['m_CaCl2'])
            row['bucket'] = 'endpoint' if not row['a'] else 'I<=3' if row['a']+3*row['b'] <= 3 else 'I>3'
        row['ionic_strength'] = row['a']+3*row['b']
    return rows


def held_out(water, source, optimum):
    rows = mixture_rows(source)
    scenarios = [s for s in fit.read(HERE / 'inputs/frozen-realizations.csv') if s['source'] in ('all', source)]
    fields = ('m_LiCl', 'm_MgCl2', 'phi') if source == 'Yao' else ('m_KCl', 'm_CaCl2', 'phi')
    nominal, buckets, rounds = [], [], []
    for scenario in scenarios:
        residuals, details = {}, []
        for row in rows:
            if row['bucket'] == 'excluded':
                continue
            values = [row['a'], row['b'], row['phi']]
            selected = str(scenario['row']) in ('all', str(row['csv_row']))
            if selected:
                values = [v+int(scenario['sign'])*0.00005 if v and scenario['field'] in ('all', field) else v
                          for v,field in zip(values, fields, strict=True)]
            a,b,phi = values
            predicted = water.predict(li=a, mg=b) if source == 'Yao' else water.predict(k=a, ca=b)
            r = predicted['phi']-phi
            residuals[row['csv_row']] = r
            details.append({'optimum': optimum, 'source': source, 'realization': scenario['id'], 'row': row['csv_row'],
                            'locator': row['locator'], 'set': row['set'], 'bucket': row['bucket'],
                            'm_first_salt': a, 'm_divalent': b, 'I': a+3*b, 'observed': phi,
                            'model': predicted['phi'], 'raw_residual': r, **{k:v for k,v in predicted.items() if k not in ('phi', 'aw')}})
        for out in details:
            out['mixing_residual'] = None
            if source == 'Yao' and out['bucket'] in ('I<=3', '3<I<=8.3'):
                ends = [r for r in rows if r['set'] == out['set'] and r['bucket'] == 'endpoint']
                (li,) = [r for r in ends if r['b'] == 0]
                (mg,) = [r for r in ends if r['a'] == 0]
                a,b = out['m_first_salt'],out['m_divalent']
                x = 3*b/(2*a+3*b)
                out['mixing_residual'] = out['raw_residual'] - (x*residuals[mg['csv_row']] + (1-x)*residuals[li['csv_row']])
        if scenario['id'] == 'nominal':
            nominal = details
        for bucket in sorted({r['bucket'] for r in details} - {'endpoint'}):
            selected = [r for r in details if r['bucket'] == bucket]
            raw = metrics([r['raw_residual'] for r in selected])
            mixing = metrics([r['mixing_residual'] for r in selected]) if selected[0]['mixing_residual'] is not None else None
            rounds.append({'optimum': optimum, 'source': source, 'realization': scenario['id'], 'bucket': bucket,
                           **raw, 'passed': passes(raw, D['limits']['Phi']) if bucket != 'report-only' else None,
                           'mixing_RMS': mixing['RMS'] if mixing else None,
                           'mixing_maximum': mixing['maximum_absolute_error'] if mixing else None,
                           'mixing_passed': passes(mixing, D['limits']['Phi']) if mixing else None})
    for bucket in sorted({r['bucket'] for r in rounds}):
        group = [r for r in rounds if r['bucket'] == bucket]
        first = next(r for r in group if r['realization'] == 'nominal')
        buckets.append({**first, 'scored': bucket != 'report-only', 'realizations': len(group),
                        'passing_realizations': sum(r['passed'] for r in group) if bucket != 'report-only' else None,
                        'all_passed': all(r['passed'] for r in group) if bucket != 'report-only' else None,
                        'worst_RMS': max(r['RMS'] for r in group), 'worst_maximum': max(r['maximum_absolute_error'] for r in group),
                        'mixing_all_passed': all(r['mixing_passed'] for r in group) if first['mixing_passed'] is not None else None,
                        'worst_mixing_RMS': max(r['mixing_RMS'] for r in group) if first['mixing_passed'] is not None else None,
                        'worst_mixing_maximum': max(r['mixing_maximum'] for r in group) if first['mixing_passed'] is not None else None})
    return nominal, rounds, buckets


def transfers(p):
    mixture = Mixture(p)
    values, diagnostics = {}, []
    for row in fit.read(HERE / D['lithium']['path']):
        m = {e: float(row[e])/1000 for e in ION if e != 'Cl'}
        m['Ca'] += float(row['Sr'])/1000
        m['Cl'] = math.fsum(n*ION[e][1] for e,n in m.items())
        amounts = {'water': float(row['water_mol']), **{ION[e][0]:n for e,n in m.items()}}
        total = math.fsum(amounts.values())
        x = [amounts[c]/total for c in p.component_ids]
        state = mixture.state(T, P=P, x=x, phase='liquid')
        a = {e: (math.log(x[p.component_ids.index(s)]*state.molar_density) + state.residual_chemical_potential_over_rt[p.component_ids.index(s)])/math.log(10)
             for e,(s,_) in ION.items()}
        values[row['state']] = combinations(a)
        diagnostics.append({'state': row['state'], 'density_mol_m3': state.molar_density,
                            'pressure_error_Pa': state.pressure-P, 'charge_mol': sum(m[e]*ION[e][1] for e in m)})
    return {f'{target}/{ratio}': values[target][ratio]-values['fig7'][ratio]
            for target in ('table1','ms2') for ratio in values['fig7']}, diagnostics


def combinations(a):
    # Lithium #188's neutral activity combinations; species reference offsets cancel on transfer.
    return {'F2': a['Li']-a['Na'], 'F2-prime': a['Li']+a['Cl'], 'F4-Mg': 2*a['Na']-a['Mg'],
            'F4-Ca': 2*a['Na']-a['Ca'], 'F5': a['K']+a['Cl']}


def main():
    began = time.perf_counter()
    fit.verify()
    assert 'EPCSAFT_HSL_LIBRARY' not in os.environ
    record = json.loads((RESULTS / 'fit.json').read_text())
    baseline = fit.parameters()
    control, control_diag = transfers(baseline)
    pitzer = json.loads((RESULTS / 'pitzer.json').read_text())
    pv = {r['state']: combinations(r['log10_activity']) for r in pitzer['states']}
    pt = {f'{target}/{ratio}': pv[target][ratio]-pv['fig7'][ratio] for target in ('table1','ms2') for ratio in pv['fig7']}
    f6_nominal, f6_rounds, f6 = held_out(Water(baseline), 'Robinson-Covington', 'unchanged-K-Ca')
    nominal, rounds, verdicts, traces, lithium_rows, all_training, summaries, diagnostics = [], [], f6[:], [], [], [], [], []
    parameter_diffs = []
    for index,o in enumerate(v for v in record['optima'] if v['within_band']):
        label = f'optimum-{index+1}'
        p = fit.parameters(o['physical'])
        water = Water(p)
        before = [{k: getattr(e,k) for k in ('identity','field','index','magnitude','unit')} for e in baseline._document.entries]
        after = [{k: getattr(e,k) for k in ('identity','field','index','magnitude','unit')} for e in p._document.entries]
        changes = [{'optimum': label, 'before': a, 'after': b} for a,b in zip(before,after,strict=True) if a != b]
        assert len(changes) == len(fit.LABELS)
        assert all('lithium-cation' in c['before']['identity'] for c in changes)
        parameter_diffs += changes
        f6_change = max(abs(water.predict(k=r['m_first_salt'], ca=r['m_divalent'])['phi']-r['model']) for r in f6_nominal)
        assert f6_change <= 1e-10
        yn,yr,yb = held_out(water, 'Yao', label)
        nominal += yn
        rounds += yr
        verdicts += yb
        trace_p = fit.parameters(o['physical'], components=TRACE_COMPONENTS)
        activity = Activities(trace_p)
        f3 = [{**activity.row(m), 'optimum': label, 'Song_theta_eff': D['F3']['Song_2010_report_only']['theta_Li_Na'] + 0.5*m*D['F3']['Song_2010_report_only']['psi_Li_Na_Cl'],
               'binary_route_error_decade': activity.binary_check(m)} for m in D['F3']['molality_NaCl']]
        assert max(abs(r['binary_route_error_decade']) for r in f3) <= 1e-10
        traces += f3
        engine,diag = transfers(p)
        diagnostics += [{'optimum': label, **r} for r in diag]
        lr = [{'optimum': label, 'state': key.split('/')[0], 'check': key.split('/')[1], 'engine_transfer': value,
               'pitzer_transfer': pt[key], 'difference_decade': value-pt[key], 'baseline_difference': control[key]-pt[key],
               'change_from_baseline': value-control[key], 'target_met': abs(value-pt[key]) <= D['limits']['Lithium_abs_decade']}
              for key,value in engine.items()]
        for row in lr:
            row['F4_preserved'] = abs(row['change_from_baseline']) <= D['limits']['F4_change_decade'] if row['check'].startswith('F4') else None
        lithium_rows += lr
        training = []
        for row in fit.training(p):
            predicted = water.predict(li=row['m_LiCl'], na=row['m_NaCl'])
            observed_aw = math.exp(-p.molar_masses[0]*2*(row['m_LiCl']+row['m_NaCl'])*row['observed'])
            training.append({**{k:v for k,v in row.items() if k != 'observation'}, 'optimum': label,
                             'model': predicted['phi'], 'residual': predicted['phi']-row['observed'],
                             'observed_aw': observed_aw, 'model_aw': predicted['aw'], 'aw_residual': predicted['aw']-observed_aw})
        all_training += training
        guard = metrics([r['residual'] for r in training if r['source'] == 'Hamer-Wu'])
        eg = {b: metrics([r['aw_residual'] for r in training if r['source'] == 'El Guendouzi' and ('I<=3' if r['m_LiCl']+r['m_NaCl'] <= 3 else 'I>3') == b]) for b in ('I<=3','I>3')}
        exceeds = any(not passes(eg[b], D['reference_model_aw_residual'][b]) for b in eg)
        f3_ok = all(abs(r['theta_eff']) <= D['limits']['F3_abs_theta_kg_mol'] for r in f3)
        f1_ok = all(b['all_passed'] for b in yb if b['scored'])
        mixing_failure = any(not b['mixing_all_passed'] for b in yb if b['scored'])
        f4_ok = all(r['target_met'] and r['F4_preserved'] for r in lr if r['check'].startswith('F4'))
        guard_ok = passes(guard, D['Hamer_Wu_guard'])
        f6_ok = all(b['all_passed'] for b in f6)
        trigger = not f3_ok or (mixing_failure and exceeds)
        summaries.append({'optimum': label, 'starts': o['starts'], 'cost': o['final_cost'], 'F1_passed': f1_ok,
                          'F3_passed': f3_ok, 'F4_passed': f4_ok, 'F6_passed': f6_ok,
                          'Hamer_Wu': guard, 'Hamer_Wu_guard_passed': guard_ok, 'El_Guendouzi_aw': eg,
                          'El_Guendouzi_exceeds_reference': exceeds, 'Yao_mixing_failure': mixing_failure,
                          'stage_B_trigger': trigger, 'stop': not (f1_ok and f3_ok and f4_ok and f6_ok and guard_ok),
                          'packet_eligible': f1_ok and f3_ok and f4_ok and f6_ok and guard_ok and not exceeds,
                          'active_bounds': o['active_bounds'], 'physical': o['physical'], 'F6_max_change': f6_change,
                          'unchanged_native_parameter_entries': len(before)-len(changes)})
    for name,rows in [('held-out-nominal.csv', nominal+f6_nominal), ('rounding-metrics.csv', rounds+f6_rounds),
                      ('buckets.csv', verdicts), ('trace.csv', traces), ('lithium-transfers.csv', lithium_rows),
                      ('training.csv', all_training), ('lithium-density-checks.csv', diagnostics)]:
        fit.write(RESULTS / name, rows)
    elapsed = time.perf_counter()-began
    summary = {'declaration_sha256': fit.sha(fit.DECL), 'selected': summaries[0], 'optima': summaries,
               'optimum_dependent': len({(s['F1_passed'],s['F3_passed'],s['Hamer_Wu_guard_passed']) for s in summaries}) > 1,
               'evaluation_wall_seconds': elapsed, 'scientific_wall_seconds': record['wall_seconds']+elapsed+pitzer['wall_seconds'],
               'native_failure_count': 0, 'packet_created': False, 'stage_B_run': False,
               'missing_evidence': D['missing'], 'baseline_lithium_diagnostics': control_diag}
    (RESULTS / 'evaluation.json').write_text(json.dumps(summary, indent=2) + '\n')
    (RESULTS / 'parameter-changes.json').write_text(json.dumps(parameter_diffs, indent=2) + '\n')
    print(json.dumps(summary['selected'], indent=2), flush=True)


if __name__ == '__main__':
    main()
