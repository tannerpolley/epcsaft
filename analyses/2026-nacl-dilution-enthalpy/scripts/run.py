"""Native NaCl dilution readouts and the amended empirical sampling procedure.

The extensive balances and their reporting-weight derivatives follow the
accepted design. All EOS values, actions and pressure solves stay native.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

from epcsaft import (Mixture, Pure, Parameters, PressureActionDirection,
                    PropertyObservable as Obs, IdealPolynomial, IdealCorrelation,
                    IdealInterval, ThermochemistryRecord, _core,
                    GAS_CONSTANT_J_PER_MOL_K as R)

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
MW, P = .01801528, 101325.
GAMMA16 = 16 * 2.**-53 / (1 - 16 * 2.**-53)
SPECIES = ['water', 'sodium-cation', 'chloride-anion']
DOCUMENT = json.loads((HERE / 'inputs/parameters.json').read_text())
CRITERIA = dict(pressure_closure_Pa=1.e-3, density_relative_spread=1.e-8,
    benchmark_coefficient=.003, refinement_relative_to_benchmark=1.e-4,
    numerical_relative_to_benchmark=.01, H_identity_J_mol=1.e-7,
    H_T_Cp_identity_J_mol_K=1.e-7, ideal_Cp_identity_J_mol_K=1.e-8, FD_J_mol_K=1.e-5)
SETTINGS = {'ordinary': {'ASSOCIATION_RESIDUAL': 2.e-13, 'DENSITY_ROOT_RESIDUAL': 1.e-12},
            'tight': {'ASSOCIATION_RESIDUAL': 2.e-14, 'DENSITY_ROOT_RESIDUAL': 1.e-13}}


def thermal(species, alternate=False):
    values = ({'water': ([40., .1], 100.), 'sodium-cation': ([30., -.1], 10000.),
               'chloride-anion': ([25., .05], -5000.)} if alternate else
              {'water': ([33.6], 0.), 'sodium-cation': ([20.8], 0.), 'chloride-anion': ([20.8], 0.)})
    return ThermochemistryRecord(P, [[IdealInterval(270., 375., True, True,
        IdealCorrelation(IdealPolynomial(values[s][0], 298.15, values[s][1], 0.), P))] for s in species])


def models(alternate=False):
    return ({k: Mixture(Parameters.from_mapping(DOCUMENT, components=SPECIES),
                       thermochemistry=thermal(SPECIES, alternate), solver_options=v) for k, v in SETTINGS.items()},
            {k: Pure(Parameters.from_mapping(DOCUMENT, components=['water']),
                     thermochemistry=thermal(['water'], alternate), solver_options=v) for k, v in SETTINGS.items()})


def observations():
    rows = list(csv.DictReader((HERE / 'inputs/measurement-candidates.csv').open()))
    # Independent re-extraction of the three measured columns preserves signs,
    # trailing zeros, replicates and half-last-digit radii from the source pages.
    text = (HERE / 'inputs/sources/partanen2023-tables-s19-s21.txt').read_text()
    for table, count in [('S19', 19), ('S20', 23), ('S21', 18)]:
        section = text.split('Table ' + table + ' Heats of dilution')[1].split('\f')[0]
        triples = [line.split()[:3] for line in section.splitlines()
                   if line.split() and line.split()[0][0].isdigit() and len(line.split()) >= 7]
        group = [r for r in rows if r['table'] == table]
        assert len(group) == len(triples) == count
        for r, triple in zip(group, triples, strict=True):
            for key, radius, value in zip(['mi_mol_kg', 'mf_mol_kg', 'H_J_mol'],
                    ['half_printed_digit_mi_mol_kg', 'half_printed_digit_mf_mol_kg', 'half_printed_digit_H_J_mol'], triple, strict=True):
                value = value.replace('−', '-').replace('–', '-')
                assert float(r[key]) == float(value), (table, r['source_row'], key)
                assert r['source_line'].split()[['mi_mol_kg', 'mf_mol_kg', 'H_J_mol'].index(key)].replace('−', '-').replace('–', '-') == value
                assert float(r[radius]) == .5 * 10**(-len(value.split('.')[1]) if '.' in value else 0)
            assert (r['included_primary'] == 'True') == (float(r['mi_mol_kg']) <= 5.)
    assert len(rows) == 60 and sum(r['included_primary'] == 'True' for r in rows) == 54
    return rows


def declarations():
    models()
    models(True)
    observations()
    PressureActionDirection(1., 0., [0., 0., 0.])


def write_json(name, value):
    (HERE / 'results' / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def write_csv(name, rows):
    with (HERE / 'results' / name).open('w') as stream:
        writer = csv.DictWriter(stream, sorted({k for row in rows for k in row}), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check_state(st):
    assert all(math.isfinite(v) for v in [st.residual_enthalpy, st.residual_isobaric_heat_capacity,
        st.total_enthalpy, st.total_isobaric_heat_capacity, st.pressure, st.molar_density,
        st.pressure_density_derivative]), 'nonfinite state readout'
    assert st.molar_density > 0 and st.pressure_density_derivative > 0, 'nonpositive density or pressure slope'
    assert abs(st.pressure - P) <= CRITERIA['pressure_closure_Pa'], 'pressure closure exceeds 0.001 Pa'


def evaluate(model, T, m, setting, start, records, temperature_actions):
    nw = 1./(MW*m) if m is not None else 0.
    N = nw + 2.
    x = [nw/N, 1./N, 1./N] if m is not None else [1.]
    kw = {'x': x} if m is not None else {}
    row = dict(T_K=T, m_mol_kg=m, setting=setting, start=start, x=x, nw_mol=nw, N_mol=N)
    records.append(row)
    try:
        seed = None
        if start == 'anchored':
            seed = model.state(T - (.1 if m is not None else 1.), P=P, phase='liquid', **kw)
            row['seed'] = dict(T_K=T - (.1 if m is not None else 1.), rho=seed.molar_density, pressure_error_Pa=seed.pressure-P,
                               pressure_density_slope=seed.pressure_density_derivative)
            check_state(seed)
        st = model.state(T, P=P, phase='liquid', anchor=seed, **kw)
        row.update(rho=st.molar_density, pressure_error_Pa=st.pressure-P,
                   pressure_density_slope=st.pressure_density_derivative,
                   h_residual_J_mol=st.residual_enthalpy, cp_residual_J_mol_K=st.residual_isobaric_heat_capacity,
                   h_total_J_mol=st.total_enthalpy, cp_total_J_mol_K=st.total_isobaric_heat_capacity)
        check_state(st)
        a = MW/(1.+2.*MW*m)**2 if m is not None else 0.
        dx = [-2*a, a, a] if m is not None else [0.]
        batches = {'molality': model.pressure_actions(T, P, phase='liquid', anchor=st, **kw,
            observables=[Obs.ResidualEnthalpy], directions=[PressureActionDirection(0., 0., dx)])}
        if temperature_actions:
            batches['temperature'] = model.pressure_actions(T, P, phase='liquid', anchor=st, **kw,
                observables=[Obs.ResidualEnthalpy, Obs.ResidualIsobaricHeatCapacity],
                directions=[PressureActionDirection(1., 0., [0.]*len(x))])
        row['actions'] = []
        for direction, batch in batches.items():
            row['actions'].extend(dict(direction=direction, observable=str(r.observable), status=r.status.name,
                                      value=r.value, action=r.action) for r in batch.results)
            diag = dict(association_residual=batch.diagnostics.association_residual,
                        pressure_error_Pa=batch.diagnostics.pressure_residual,
                        pressure_density_slope=batch.diagnostics.pressure_density_derivative)
            row[direction+'_diagnostics'] = diag
            assert all(r.status.name == 'Available' and math.isfinite(r.value) and math.isfinite(r.action)
                       for r in batch.results), 'unavailable or nonfinite action'
            assert diag['association_residual'] <= SETTINGS[setting]['ASSOCIATION_RESIDUAL'], 'association diagnostic exceeds requested criterion'
            assert abs(diag['pressure_error_Pa']) <= CRITERIA['pressure_closure_Pa']
            assert diag['pressure_density_slope'] > 0
        row['dh_dm'] = batches['molality'].results[0].action
        if temperature_actions:
            row.update(dh_dT=batches['temperature'].results[0].action, dcp_dT=batches['temperature'].results[1].action)
    except Exception as exc:
        row['failure'] = repr(exc)
        raise
    return row


def endpoint(mixes, waters, T, setting, start, m, cache, records, temperature_actions=False):
    key = (T, setting, start, m, temperature_actions)
    if key in cache:
        row = cache[key]
        if 'failure' in row:
            raise RuntimeError(row['failure'])
        return row
    try:
        row = evaluate(mixes[setting] if m is not None else waters[setting], T, m, setting, start, records, temperature_actions)
        cache[key] = row
        if m is not None:
            water = endpoint(mixes, waters, T, setting, start, None, cache, records, temperature_actions)
            nw, N = row['nw_mol'], row['N_mol']
            for quantity, field in [('H', 'h_residual_J_mol'), ('Cp', 'cp_residual_J_mol_K')]:
                A, C = N*row[field], nw*water[field]
                row.update({f'A_{quantity}': A, f'C_{quantity}': C, f'B_{quantity}': A-C})
            row.update(B_H_total=N*row['h_total_J_mol']-nw*water['h_total_J_mol'],
                       B_Cp_total=N*row['cp_total_J_mol_K']-nw*water['cp_total_J_mol_K'],
                       dB_H_dm=(-nw/m)*(row['h_residual_J_mol']-water['h_residual_J_mol'])+N*row['dh_dm'])
            if temperature_actions:
                row.update(dB_H_dT=N*row['dh_dT']-nw*water['dh_dT'],
                           dB_Cp_dT=N*row['dcp_dT']-nw*water['dcp_dT'])
    except Exception as exc:
        row = cache.get(key, {'T_K': T, 'm_mol_kg': m, 'setting': setting, 'start': start})
        row['failure'] = repr(exc)
        cache[key] = row
        raise
    return row


def sample_chord(T, mi, mf, rows, mixes, waters, cache, records, samples):
    U = CRITERIA['benchmark_coefficient']*R*T
    tau = CRITERIA['refinement_relative_to_benchmark']*U
    ri, rf = [float(rows[0][k]) for k in ['half_printed_digit_mi_mol_kg', 'half_printed_digit_mf_mol_kg']]
    history, previous, stable, cause = [], None, 0, None
    nominal = None
    try:
        nominal = endpoint(mixes, waters, T, 'tight', 'cold', mf, cache, records)['B_H'] - endpoint(
            mixes, waters, T, 'tight', 'cold', mi, cache, records)['B_H']
        # Pure water must agree across the same four variants, just like endpoints.
        ds = [endpoint(mixes, waters, T, s, v, None, cache, records)['rho']
              for s in SETTINGS for v in ['cold', 'anchored']]
        assert (max(ds)-min(ds))/min(ds) <= CRITERIA['density_relative_spread'], 'water root disagreement'
        for n in [2, 4, 8, 16, 32, 64]:
            level = []
            for i in [mi-ri+2*ri*k/n for k in range(n+1)]:
                for f in [mf-rf+2*rf*k/n for k in range(n+1)]:
                    assert 0 < f < i
                    variants, fp, gradients, densities = {}, [], {}, {'i': [], 'f': []}
                    for setting in SETTINGS:
                        for start in ['cold', 'anchored']:
                            ai = endpoint(mixes, waters, T, setting, start, i, cache, records)
                            af = endpoint(mixes, waters, T, setting, start, f, cache, records)
                            name = setting+'_'+start
                            variants[name] = af['B_H']-ai['B_H']
                            fp.append(GAMMA16*(abs(af['A_H'])+abs(af['C_H'])+abs(ai['A_H'])+abs(ai['C_H'])))
                            gradients[name] = [-ai['dB_H_dm'], af['dB_H_dm']]
                            densities['i'].append(ai['rho'])
                            densities['f'].append(af['rho'])
                    for ds in densities.values():
                        assert (max(ds)-min(ds))/min(ds) <= CRITERIA['density_relative_spread'], 'liquid-root density disagreement'
                    ref = variants['tight_cold']
                    row = dict(T_K=T, mi_nominal=mi, mf_nominal=mf, n=n, mi=i, mf=f,
                        reference_J_mol=ref, variants_J_mol=variants,
                        E_J_mol=max(abs(v-ref) for v in variants.values()), F_J_mol=max(fp), gradients=gradients)
                    samples.append(row)
                    level.append(row)
            lo, hi = min(r['reference_J_mol'] for r in level), max(r['reference_J_mol'] for r in level)
            E, F = max(r['E_J_mol'] for r in level), max(r['F_J_mol'] for r in level)
            bnum = 2*E+F
            gs = [[v[d] for r in level for v in r['gradients'].values()] for d in [0, 1]]
            sign_resolved = all(min(v) > 0 or max(v) < 0 for v in gs)
            change = max(abs(lo-previous[0]), abs(hi-previous[1])) if previous else None
            stable = stable+1 if change is not None and change <= tau else 0
            history.append(dict(n=n, points=len(level), sampled_min_J_mol=lo, sampled_max_J_mol=hi,
                E_J_mol=E, F_J_mol=F, b_num_J_mol=bnum, endpoint_change_J_mol=change,
                consecutive_small_changes=stable, derivative_ranges=[[min(v), max(v)] for v in gs],
                sampled_signs_resolved=sign_resolved))
            previous = (lo, hi)
            if stable >= 2 and sign_resolved and bnum <= CRITERIA['numerical_relative_to_benchmark']*U:
                break
    except Exception as exc:
        cause = repr(exc)
    ok = cause is None and stable >= 2 and sign_resolved and bnum <= CRITERIA['numerical_relative_to_benchmark']*U
    item = dict(T_K=T, mi_mol_kg=mi, mf_mol_kg=mf, mi_radius=ri, mf_radius=rf,
                uniform_benchmark_J_mol=U, tau_J_mol=tau, levels=history,
                cause=cause if cause else (None if ok else 'sampling stop not reached'),
                sampling_resolved=ok, nominal_J_mol=nominal)
    if ok:
        bsample = 2*max(h['endpoint_change_J_mol'] for h in history[-2:])
        item.update(b_sample_J_mol=bsample, b_num_J_mol=bnum, interval_J_mol=[lo-bnum-bsample, hi+bnum+bsample])
    comparisons = []
    for r in rows:
        obs, radius = float(r['H_J_mol']), float(r['half_printed_digit_H_J_mol'])
        bench = [obs-U-radius, obs+U+radius]
        result = dict(table=r['table'], source_row=r['source_row'], T_K=T, mi_mol_kg=mi, mf_mol_kg=mf,
            observed_J_mol=obs, rounding_J_mol=radius, benchmark_J_mol=U, benchmark_low=bench[0], benchmark_high=bench[1],
            predicted_J_mol=nominal, signed_nominal_error_J_mol=None if nominal is None else nominal-obs,
            benchmark_normalized_error=None if nominal is None else (nominal-obs)/U,
            nonnegative_interval_gap_J_mol=None, classification='unresolved', cause=item['cause'])
        if ok:
            interval = item['interval_J_mol']
            gap = max(interval[0]-bench[1], bench[0]-interval[1], 0.)
            depth = min(interval[1], bench[1])-max(interval[0], bench[0])
            result.update(calculated_low=interval[0], calculated_high=interval[1], b_num_J_mol=bnum,
                b_sample_J_mol=bsample, nonnegative_interval_gap_J_mol=gap, overlap_depth_J_mol=depth,
                classification='disjoint' if gap > 2*tau else 'overlap' if depth > 2*tau else 'unresolved')
            if result['classification'] == 'unresolved':
                result['cause'] = 'gap and overlap depth do not exceed two tau'
        comparisons.append(result)
    print(json.dumps({k: v for k, v in item.items() if k != 'levels'}), flush=True)
    item['classification'] = ('unresolved' if any(r['classification'] == 'unresolved' for r in comparisons)
        else 'disjoint' if any(r['classification'] == 'disjoint' for r in comparisons) else 'overlap')
    return item, comparisons


def qualification(mixes, waters, cache, records):
    alternate, alternate_water = models(True)
    altcache, checks, ladders, heat, alt_records = {}, [], [], [], []
    def ep(T, m, alt=False):
        return endpoint(alternate if alt else mixes, alternate_water if alt else waters,
                        T, 'ordinary', 'cold', m, altcache if alt else cache, alt_records if alt else records, True)
    for T in [298.15, 323.15, 348.15]:
        for mi, mf in [(.0744, .0371), (3.222, 2.1)]:
            row = dict(T_K=T, mi_mol_kg=mi, mf_mol_kg=mf)
            try:
                a, b, aa, bb = ep(T, mi), ep(T, mf), ep(T, mi, True), ep(T, mf, True)
                D, C = b['B_H']-a['B_H'], b['B_Cp']-a['B_Cp']
                row.update(total_minus_residual_J_mol=b['B_H_total']-a['B_H_total']-D,
                    changed_record_H_J_mol=bb['B_H_total']-aa['B_H_total']-D,
                    changed_record_Cp_J_mol_K=bb['B_Cp_total']-aa['B_Cp_total']-C,
                    thermal_identity_J_mol_K=b['dB_H_dT']-a['dB_H_dT']-C,
                    delta_Cp_J_mol_K=C, exact_dH_dT_J_mol_K=b['dB_H_dT']-a['dB_H_dT'])
                row['passed'] = (abs(row['total_minus_residual_J_mol']) <= CRITERIA['H_identity_J_mol']
                    and abs(row['changed_record_H_J_mol']) <= CRITERIA['H_identity_J_mol']
                    and abs(row['changed_record_Cp_J_mol_K']) <= CRITERIA['ideal_Cp_identity_J_mol_K']
                    and abs(row['thermal_identity_J_mol_K']) <= CRITERIA['H_T_Cp_identity_J_mol_K'])
            except Exception as exc:
                row.update(passed=False, cause=repr(exc))
            checks.append(row)
        for m in sorted({float(r[k]) for r in observations() if float(r['T_K']) == T and r['included_primary'] == 'True'
                         for k in ['mi_mol_kg', 'mf_mol_kg']} | {.1, 1., 3.}):
            row = dict(T_K=T, m_mol_kg=m, reference_m_mol_kg=1.)
            try:
                a, b = ep(T, 1.), ep(T, m)
                row.update(delta_Cp_from_1m_J_mol_K=b['B_Cp']-a['B_Cp'],
                    exact_d_delta_Cp_dT_J_mol_K2=b['dB_Cp_dT']-a['dB_Cp_dT'],
                    H_T_minus_Cp_J_mol_K=b['dB_H_dT']-b['B_Cp'])
            except Exception as exc:
                row['cause'] = repr(exc)
            heat.append(row)
    for mi, mf in [(1.5, .0204), (3.222, 2.1)]:
        row = dict(T_K=298.15, mi_mol_kg=mi, mf_mol_kg=mf)
        try:
            exact = ep(298.15, mf)['dB_H_dT']-ep(298.15, mi)['dB_H_dT']
            values = []
            for h in [.1, .05]:
                plus = ep(298.15+h, mf)['B_H']-ep(298.15+h, mi)['B_H']
                minus = ep(298.15-h, mf)['B_H']-ep(298.15-h, mi)['B_H']
                values.append((plus-minus)/(2*h))
                row[f'plus_h_{h}_J_mol'], row[f'minus_h_{h}_J_mol'] = plus, minus
            rich = (4*values[1]-values[0])/3
            row.update(exact_J_mol_K=exact, FD_h_0_1_K=values[0], FD_h_0_05_K=values[1],
                Richardson_J_mol_K=rich, fine_minus_coarse_J_mol_K=values[1]-values[0],
                exact_minus_Richardson_J_mol_K=exact-rich, passed=abs(exact-rich) <= CRITERIA['FD_J_mol_K'])
        except Exception as exc:
            row.update(passed=False, cause=repr(exc))
        ladders.append(row)
    records.extend(r | {'ideal_record': 'alternate'} for r in alt_records)
    return checks, ladders, heat


def main():
    started = datetime.now(timezone.utc).isoformat()
    measured = observations()
    mixes, waters = models()
    records, samples, chords, comparisons, cache = [], [], [], [], {}
    keys = sorted({(float(r['T_K']), float(r['mi_mol_kg']), float(r['mf_mol_kg']))
                   for r in measured if r['included_primary'] == 'True'})
    assert len(keys) == 24
    for T, mi, mf in keys:
        rows = [r for r in measured if r['included_primary'] == 'True' and
                (float(r['T_K']), float(r['mi_mol_kg']), float(r['mf_mol_kg'])) == (T, mi, mf)]
        chord, scores = sample_chord(T, mi, mf, rows, mixes, waters, cache, records, samples)
        chords.append(chord)
        comparisons.extend(scores)
    checks, ladders, heat = qualification(mixes, waters, cache, records)
    for chord in chords:
        try:
            a, b = [endpoint(mixes, waters, chord['T_K'], 'ordinary', 'cold', chord[k], cache, records, True)
                    for k in ['mi_mol_kg', 'mf_mol_kg']]
            chord.update(delta_Cp_J_mol_K=b['B_Cp']-a['B_Cp'], exact_dH_dT_J_mol_K=b['dB_H_dT']-a['dB_H_dT'])
        except Exception as exc:
            chord['thermal_cause'] = repr(exc)
    failures = [r for r in records if 'failure' in r]
    summaries = []
    for T in [298.15, 323.15, 348.15, None]:
        rows = [r for r in comparisons if T is None or r['T_K'] == T]
        cs = [c for c in chords if T is None or c['T_K'] == T]
        errors = [r['signed_nominal_error_J_mol'] for r in rows if r['signed_nominal_error_J_mol'] is not None]
        item = dict(T_K=T, rows=len(rows), chords=len(cs), nominal_available=len(errors),
            rows_disjoint=sum(r['classification'] == 'disjoint' for r in rows),
            rows_overlap=sum(r['classification'] == 'overlap' for r in rows),
            rows_unresolved=sum(r['classification'] == 'unresolved' for r in rows),
            **{'chords_'+k: sum(c['classification'] == k for c in cs) for k in ['disjoint', 'overlap', 'unresolved']},
            positive_errors=sum(e > 0 for e in errors), negative_errors=sum(e < 0 for e in errors))
        if len(errors) == len(rows):
            item.update(mean_signed_error_J_mol=sum(errors)/len(errors),
                mean_absolute_error_J_mol=sum(map(abs, errors))/len(errors),
                RMSE_J_mol=math.sqrt(sum(e*e for e in errors)/len(errors)), max_absolute_error_J_mol=max(map(abs, errors)))
        summaries.append(item)
    write_json('native-states.json', records)
    write_json('samples.json', samples)
    write_json('chords.json', chords)
    write_json('failures.json', failures)
    write_json('qualification.json', {'identities': checks, 'thermal_ladder': ladders})
    write_csv('comparison.csv', comparisons)
    write_csv('heat-capacity.csv', heat)
    write_json('summary.json', {'groups': summaries, 'excluded_rows': [r for r in measured if r['included_primary'] != 'True'],
        'native_failures': len(failures), 'state_records': len(records), 'samples': len(samples),
        'numerical_qualification_passed': all(r['passed'] for r in checks+ladders)})
    paths = sorted(p for p in (HERE/'inputs').rglob('*') if p.is_file()) + [Path(__file__)]
    lib = Path(_core.__file__)
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    identity = dict(started_UTC=started, finished_UTC=datetime.now(timezone.utc).isoformat(),
        checkout=str(ROOT), checkout_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        python=sys.executable, native_module=str(lib), native_sha256=sha(lib),
        native_library_sha256=sha(lib.parent/'libepcsaft_native.so'),
        packet_sha256=sha(ROOT/'data/packets/held-2014-strategy2/2/parameters/parameters.json'),
        input_sha256={str(p.relative_to(HERE)): sha(p) for p in paths}, R_J_mol_K=R, MW_kg_mol=MW, P_Pa=P,
        gamma16=GAMMA16, criteria=CRITERIA, settings=SETTINGS,
        grids=[2, 4, 8, 16, 32, 64], solution_seed_delta_K=-.1, water_seed_delta_K=-1.,
        heat_capacity_reference_m_mol_kg=1., finite_difference_steps_K=[.1, .05],
        manufactured_thermal_interval_K=[270., 375.], effective_solver_options={k: m._model.solver_options for k, m in mixes.items()},
        native_defaults={v.name: v.value for v in _core.tolerances()})
    write_json('run-identity.json', identity)
    print(json.dumps(summaries, indent=2))
    assert len(comparisons) == 54
    assert all(r['passed'] for r in checks+ladders), 'numerical qualification failed; retained'


if __name__ == '__main__':
    declarations() if '--declare-only' in sys.argv else main()
