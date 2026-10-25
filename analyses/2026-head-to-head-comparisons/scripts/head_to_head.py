"""Read native head-to-head outputs (2026-10-08) and write comparison CSVs; no thermodynamics here."""
import csv, gzip, json, math, os, re, shutil, statistics, sys
from datetime import datetime, timezone
from pathlib import Path
HERE = Path(__file__).resolve().parents[1]; RES = Path(os.environ.get('EPCSAFT_H2H_OUTPUT', HERE / 'results')); CC = HERE.parent / '2026-competitor-cross-check'

def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, 'w', newline='') as f: w = csv.DictWriter(f, keys, lineterminator='\n'); w.writeheader(); w.writerows(rows)

def median_of_process_medians(rows, impl, ident, selection, kind):
    m = [float(r['median_ns']) for r in rows if (r['implementation'], r['id'], r['selection'], r['kind']) == (impl, ident, selection, kind)]
    return statistics.median(m) / 1e6 if m else None

def engine_timing(raw):
    """Fresh main-wheel timing: per-process statistics, block summary, returned-value differences."""
    runs = [json.loads((raw / f'Engine-main-timing-{k}.json').read_text()) for k in (1, 2, 3)]
    stats = []
    for d in runs:
        groups = {}
        for t in d['timings']: groups.setdefault((t['id'], t['selection'], t['kind'], t['rep']), []).append(t['elapsed_ns'])
        for (i, s, k, rep), ns in groups.items():
            p = statistics.quantiles(ns, n=10, method='inclusive') if len(ns) > 1 else [ns[0]] * 9
            stats.append(dict(id=i, implementation='Engine', selection=s, kind=k, rep=rep, samples=len(ns), median_ns=statistics.median(ns), p10_ns=p[0], p90_ns=p[8],
                              wheel_sha256=d['host']['preparation']['wheel_sha256'], load_start=d['host']['load_start'][0], niceness=d['host']['niceness']))
    write_csv(RES / 'h2h-engine-timing-statistics.csv', stats)
    old = list(csv.DictReader(open(CC / 'results/timing-statistics.csv'))); prev = list(csv.DictReader(open(CC / 'results/density-closure-final-engine-replay-statistics.csv')))
    blocks = [('C1-300', 'complete-block', 'C1', 'warm'), ('C2-0.5', 'complete-block', 'C2', 'warm'), ('C4-298.15', 'complete-block', 'C4', 'warm'), ('C4-323.15', 'complete-block', 'C4', 'warm'),
              ('C3-0.88sigma-1', 'reused-reference', '0.88sigma', 'warm'), ('C3-0.88sigma-1', 'fresh-reference', None, None),
              ('C1-gradient', 'complete-block', None, 'warm'), ('C2-gradient', 'complete-block', None, 'warm'),
              ('C1-300-common-liquid', 'complete-block', 'C1', 'warm'), ('C2-0.5-common-liquid', 'complete-block', 'C2', 'warm'), ('C3-0.88sigma-1-common-liquid', 'complete-block', '0.88sigma', 'warm')]
    summary = []
    for ident, sel, csel, fsel in blocks:
        new = median_of_process_medians(stats_as_str(stats), 'Engine', ident, sel, 'warm')
        row = dict(id=ident, engine_selection=sel, engine_main_026b3031_ms=new,
                   engine_original_ms=median_of_process_medians(old, 'Engine', ident, sel, 'warm'),
                   engine_835d9bc2_ms=median_of_process_medians(prev, 'Engine', ident, sel, 'warm'),
                   feos_ms=median_of_process_medians(old, 'FeOs', ident, 'warm', 'warm') if fsel else None,
                   clapeyron_ms=median_of_process_medians(old, 'Clapeyron', ident, csel, 'warm') if csel else None)
        for k in ('feos', 'clapeyron'):
            row[f'engine_over_{k}'] = new / row[f'{k}_ms'] if row[f'{k}_ms'] else None
        pm = [float(s['median_ns']) / 1e6 for s in stats if (s['id'], s['selection'], s['kind']) == (ident, sel, 'warm')]
        row['process_medians_ms'] = '/'.join(f'{v:.3f}' for v in pm)
        # Scatter of the three process medians about their median: largest deviation, and full range, each over the median.
        row['process_max_deviation_over_median'] = max(abs(v - new) for v in pm) / new
        row['process_range_over_median'] = (max(pm) - min(pm)) / new
        summary.append(row)
    write_csv(RES / 'h2h-engine-timing-summary.csv', summary)
    # Returned values of the timed calls, against the retained Engine and comparator values; the Engine
    # comparison uses the retained quantity-specific metric and target of density-closure-final-values.csv.
    retained = {(r['id'], r['implementation'], r['phase'], r['component'], r['quantity']): r for r in csv.DictReader(open(CC / 'results/density-closure-final-values.csv'))}
    vals = []
    for d in runs:
        for rec in d['records']:
            q = [('', '', 'P', rec.get('P'))] if 'P' in rec else []
            q += [(p['phase'], '', 'rho', p['rho']) for p in rec.get('phases', [])]
            q += [(p['phase'], str(i), 'x', v) for p in rec.get('phases', []) for i, v in enumerate(p['x'])]
            q += [(p['phase'], '', k, p[k]) for p in rec.get('phases', []) for k in ('ar', 'Z') if k in p]
            q += [(p['phase'], str(i), 'lnphi', v) for p in rec.get('phases', []) for i, v in enumerate(p.get('lnphi', []))]
            q += [('', '', k, rec[k]) for k in ('gamma', 'phi', 'beta') if k in rec]
            q += [('', str(i), 'gradient', g) for i, g in enumerate(rec.get('gradients', []))]
            for ph, comp, name, v in q:
                row = dict(id=rec['id'], selection=rec['selection'], rep=runs.index(d) + 1, phase=ph, component=comp, quantity=name, engine_main=v)
                for impl in ('Engine', 'FeOs', 'Clapeyron'):
                    r = retained.get((rec['id'], impl, ph, comp, name)); ref = None if r is None else float(r['value'])
                    row[f'{impl}_retained'] = ref
                    row[f'main_minus_{impl}'] = None if ref is None else v - ref
                    row[f'main_minus_{impl}_relative'] = None if ref in (None, 0) or name == 'x' else (v - ref) / abs(ref)
                r = retained.get((rec['id'], 'Engine', ph, comp, name))
                if r and r['target']:
                    ref = float(r['value'])
                    diff = {'relative': lambda: abs(v - ref) / abs(ref), 'absolute': lambda: abs(v - ref), 'absolute_log': lambda: abs(math.log(v) - math.log(ref))}[r['metric']]()
                    row.update(engine_metric=r['metric'], engine_difference=diff, engine_target=float(r['target']), engine_within_target=diff <= float(r['target']))
                vals.append(row)
    write_csv(RES / 'h2h-engine-timed-values.csv', vals)
    (RES / 'h2h-engine-raw').mkdir(exist_ok=True)
    for f in sorted(raw.glob('Engine-main*.json')):
        with open(f, 'rb') as a, gzip.GzipFile(RES / 'h2h-engine-raw' / (f.name + '.gz'), 'wb', mtime=0) as b: shutil.copyfileobj(a, b)
    for r in summary: print(r)

def stats_as_str(stats): return [{k: str(v) for k, v in s.items()} for s in stats]

def rel(a, b): return None if a is None or b is None or b == 0 else (a - b) / abs(b)

VALID_PRESSURE = 1e-6  # relative pressure error above which a returned flash state is not an equilibrium
SEPARATION = 1e-8  # largest mole-fraction difference below which two returned phases are the same phase

def flash_status(d):
    """A Clapeyron flash return is an equilibrium only with a finite phase fraction on the specified pressure;
    composition separation alone does not count. Invalid returns stay in the rows with their diagnostics."""
    if 'error' in d: return 'error'
    if not math.isfinite(d['beta']) or d['residuals']['pressure_relative'] > VALID_PRESSURE: return 'invalid: no valid equilibrium returned'
    return 'two-phase' if max(abs(a - b) for a, b in zip(*d['x'])) > SEPARATION else 'single-phase'

def comparators(raw):
    """Electrolyte LLE, solved-state sensitivities and Born/permittivity rows from the three native outputs."""
    e = json.loads((raw / 'engine-h2h.json').read_text()); f = json.loads((raw / 'feos-h2h.json').read_text()); c = json.loads((raw / 'clapeyron-h2h.json').read_text())
    ids = ('water', 'ethanol', 'isobutanol', 'sodium-cation', 'chloride-anion'); lle, conv = [], []
    packet = {r['point']: r for r in e['khudaida']['packet']}
    for variant in ('constant-permittivity', 'matched'):
        crow = {r['point']: r for r in c['khudaida'][variant]['rows']}
        for m in e['khudaida'][variant]:
            p = packet[m['point']]; cr = crow[m['point']]; cd, cg = cr['default-initialization'], cr['engine-guesses']
            def got(d, key, k, j=None):
                if flash_status(d) != 'two-phase': return None
                return d[key][k] if j is None else d[key][k][j]
            for k, phase in enumerate(('water-rich', 'organic')):
                for j, comp in enumerate(ids + ('molar_density_mol_m3',)):
                    key, jj = ('rho', None) if j == 5 else ('x', j)
                    ev = m[key][k] if jj is None else m[key][k][jj]
                    row = dict(variant=variant, point=m['point'], phase=phase, quantity=comp, engine=ev, engine_packet_v7=p[key][k] if jj is None else p[key][k][jj])
                    for tag, d in (('clapeyron_default_init', cd), ('clapeyron_engine_guesses', cg)):
                        v = got(d, key, k, jj); row[tag] = v; row[f'{tag}_minus_engine'] = None if v is None else v - ev
                    lle.append(row)
            row = dict(variant=variant, point=m['point'], phase='water-rich', quantity='phase_fraction_beta', engine=m['beta'], engine_packet_v7=p['beta'])
            for tag, d in (('clapeyron_default_init', cd), ('clapeyron_engine_guesses', cg)):
                v = d['beta'] if flash_status(d) == 'two-phase' else None; row[tag] = v; row[f'{tag}_minus_engine'] = None if v is None else v - m['beta']
            lle.append(row)
            out = dict(variant=variant, point=m['point'], engine_success=m['success'], engine_phase_count=m['phase_count'], engine_iterations=m['iterations'], engine_evaluations=m['evaluations'],
                       engine_max_residual=m['max_residual'], engine_first_ms=m['cost']['first_ms'], engine_median_ms=m['cost']['median_ms'])
            for tag, d in (('clapeyron_default_init', cd), ('clapeyron_engine_guesses', cg)):
                out[f'{tag}_status'] = flash_status(d); out[f'{tag}_error'] = d.get('error')
                out[f'{tag}_beta'] = d.get('beta'); out[f'{tag}_rho'] = '/'.join(map(str, d.get('rho', [])))
                out[f'{tag}_max_abs_x_minus_engine'] = max(abs(a - b) for k in range(2) for a, b in zip(d['x'][k], m['x'][k])) if flash_status(d) == 'two-phase' else None
                for k, v in d.get('residuals', {}).items(): out[f'{tag}_{k}'] = v
                if 'cost' in d: out[f'{tag}_first_ms'] = d['cost']['first_ms']; out[f'{tag}_median_ms'] = d['cost']['median_ms']
            conv.append(out)
        for i, chk in enumerate(c['khudaida'][variant]['dh_mu_check_point1']):
            for j, comp in enumerate(ids):
                conv.append(dict(variant=variant, point=1, diagnostic=f"Clapeyron DH mu/RT at Engine phase {('water-rich', 'organic')[i]}: {comp}", automatic=chk['automatic'][j], central_difference=chk['central_difference'][j]))
    write_csv(RES / 'h2h-khudaida-lle.csv', lle); write_csv(RES / 'h2h-khudaida-convergence.csv', conv)
    sens, costs = [], []
    cs = c.get('sensitivities', {}); fg = f['C2_gradient']['value']
    for case, s in e['sensitivities'].items():
        cc = {r['direction']: r['derivative'] for r in cs.get(case, {}).get('rows', [])}
        for r in s['rows']:
            i = s['names'].index(r['observable']); cv = cc.get(r['direction'], [None] * 3)[i]
            fv = fg['gradients'][s['coordinates'].index(r['direction'])] if case.startswith('C2') and r['observable'] == 'P_Pa' and r['direction'] != 'T' else None
            sens.append(dict(case=case, observable=r['observable'], direction=r['direction'], central_value=r['value'], engine_exact=r['exact'], engine_central_fd=r['fd'], fd_step=r['fd_step'],
                             engine_exact_vs_fd_relative=rel(r['exact'], r['fd']), clapeyron_ift=cv, clapeyron_vs_engine_relative=rel(cv, r['exact']),
                             feos_ad=fv, feos_vs_engine_relative=rel(fv, r['exact']), engine_status=r['status']))
        for k, v in s['cost'].items(): costs.append(dict(case=case, implementation='Engine', work=k, **v))
        for r in cs.get(case, {}).get('rows', []): costs.append(dict(case=case, implementation='Clapeyron', work=f"solve + IFT derivative, direction {r['direction']}", **r['cost']))
        if case in cs:
            costs.append(dict(case=case, implementation='Clapeyron', work='primal solve only', **cs[case]['primal_cost']))
            costs.append(dict(case=case, implementation='Clapeyron', work='solve + joint IFT Jacobian, k_ij and T', **cs[case]['two_direction_jacobian_cost']))
    costs.append(dict(case='C2-0.5', implementation='FeOs', work='bubble_point_pressure_derivatives (k_ij, epsilon_k_ab, kappa_ab), incl. its pressure-guess bubble solve', **f['C2_gradient']['cost']))
    write_csv(RES / 'h2h-sensitivities.csv', sens); write_csv(RES / 'h2h-sensitivity-costs.csv', costs)
    born = []; F = {(r['form'], r['m']): r for r in f['born']}; C = {(r['form'], r['m']): r for r in c.get('born', [])}
    ce, cc_ = e_constants(), c.get('constants')
    # Born prefactor e^2/(eps0 kB) ratios to the Engine's; FeOs uses eps0 = 8.85416e-12 F/m with the Engine's e and kB.
    ratio = dict(feos=ce['epsilon0'] / 8.85416e-12, clapeyron=(cc_['e'] / ce['e'])**2 * ce['epsilon0'] / cc_['epsilon0'] * ce['kB'] / cc_['kB'])
    for r in e['born']:
        row = dict(form=r['form'], born=r['born'], m=r['m'], rho=r['rho'], engine_permittivity=r['permittivity'], engine_a_born=r['a_born'])
        if r['born'] == 'direct-hs-diameter':
            fr, cr = F.get((r['form'], r['m'])), C.get(('schreckenberg' if r['form'] == 'schreckenberg' else r['form'], r['m']))
            row.update(feos_a_born=fr and fr['a_born'], feos_a_born_relative=fr and rel(fr['a_born'], r['a_born']),
                       clapeyron_permittivity=cr and cr['permittivity'], clapeyron_permittivity_relative=cr and rel(cr['permittivity'], r['permittivity']),
                       clapeyron_a_born=cr and cr['a_born'], clapeyron_a_born_relative=cr and rel(cr['a_born'], r['a_born']),
                       feos_constant_ratio=ratio['feos'], feos_a_born_relative_after_ratio=fr and rel(fr['a_born'] / ratio['feos'], r['a_born']),
                       clapeyron_constant_ratio=ratio['clapeyron'], clapeyron_a_born_relative_after_ratio=cr and rel(cr['a_born'] / ratio['clapeyron'], r['a_born']))
        born.append(row)
    write_csv(RES / 'h2h-born-permittivity.csv', born)
    (RES / 'h2h-raw').mkdir(exist_ok=True)
    for name in ('engine-h2h.json', 'feos-h2h.json', 'clapeyron-h2h.json', 'clapeyron-h2h-process1.json', 'julia.log', 'julia-process1.log'):
        if not (raw / name).is_file():
            continue
        with open(raw / name, 'rb') as a, gzip.GzipFile(RES / 'h2h-raw' / (name + '.gz'), 'wb', mtime=0) as b: shutil.copyfileobj(a, b)
    # Session time stamps written by the run commands (README), and the end of each Julia process log.
    stamps = [dict(event=n, utc=(raw / n).read_text().strip()) for n in ('c1-start', 'c1-end', 'julia-start', 'julia-tasks-start', 'julia-khudaida2-start', 'julia2-start') if (raw / n).is_file()]
    stamps += [dict(event=f'{n} last write', utc=datetime.fromtimestamp((raw / n).stat().st_mtime, timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')) for n in ('julia-process1.log', 'julia.log') if (raw / n).is_file()]
    write_csv(RES / 'h2h-session-times.csv', stamps)

def probes(raw):
    """Point-1 matched-variant term checks at the Engine's solved phases (water-rich, organic), from the retained
    probe outputs in results/h2h-probes/ and the Clapeyron probe lines of results/h2h-raw/julia-process1.log.gz."""
    P = RES / 'h2h-probes'; load = lambda n: json.loads((P / n).read_text())
    terms, mu, pert, cmu = load('probe-engine-terms.json'), load('probe-engine-mu-terms.json'), load('probe-perturbed.json'), load('probe-clapeyron-mu-terms.json')
    log = gzip.open(RES / 'h2h-raw' / 'julia-process1.log.gz', 'rt').read().splitlines()
    phases, ids, R, T = ('water-rich', 'organic'), ('water', 'ethanol', 'isobutanol', 'sodium-cation', 'chloride-anion'), 8.31446261815324, 293.15
    rows = []
    def add(check, phase, comp, qty, a, b): rows.append(dict(check=check, phase=phase, component=comp, quantity=qty, a=a, b=b, a_minus_b=a - b, relative=rel(a, b)))
    p3 = [re.match(r'probe3 (\w+) clap=(\S+) engine=(\S+)', l) for l in log if l.startswith('probe3 ') and 'clap=' in l]
    for k, m in enumerate(p3): add('Helmholtz term: Clapeyron (a) vs Engine (b)', phases[k // 8], '', m[1], float(m[2]), float(m[3]))
    p6 = [l for l in log if l.startswith('probe6 ')]
    for line, st, label in zip(p6, pert, ('solved state', '+0.05 mol water', '+0.05 mol ethanol', '+0.05 mol isobutanol')):
        v = re.findall(r'(\w+) (\S+) vs (\S+)', line)
        for name, c, e in v: add('perturbed water-rich composition: Clapeyron (a) vs Engine (b)', label, '', name, float(c), float(e))
    for k, (e, c) in enumerate(zip(mu, cmu)):
        for t in e:
            for i in range(5): add('term mu_i/RT at constant V,T: Clapeyron automatic (a) vs Engine central difference (b)', phases[k], ids[i], t, c[t][i], e[t][i])
    # Engine exact lnphi against central differences of its own Helmholtz terms, and cross-phase equality of mu/RT.
    full = {'exact': [], 'central difference': []}
    for k, (s, d) in enumerate(zip(terms, mu)):
        lnZ = math.log(s['pressure'] / (s['rho'] * R * T)); cd = [sum(d[t][i] for t in d) - lnZ for i in range(5)]
        for i in range(5): add('Engine lnphi: exact (a) vs central difference of its Helmholtz terms (b)', phases[k], ids[i], 'lnphi', s['lnphi'][i], cd[i])
        full['exact'].append([s['lnphi'][i] + math.log(s['x'][i] * s['pressure']) for i in range(5)])  # mu/RT up to a pure-component constant
        full['central difference'].append([cd[i] + lnZ + math.log(s['rho'] * s['x'][i]) for i in range(5)])
    for how, (a, b) in full.items():
        for i in range(3): add(f'Engine mu/RT equality between phases, {how}: water-rich (a) vs organic (b)', '', ids[i], 'mu/RT', a[i], b[i])
        add(f'Engine mu/RT equality between phases, {how}: water-rich (a) vs organic (b)', '', 'NaCl (cation + anion)', 'mu/RT', a[3] + a[4], b[3] + b[4])
    for name in ('engine mu diff', 'clap mu diff'):
        line = next(l for l in log if l.startswith(f'probe4 {name}'))
        for i, v in enumerate(json.loads(line.split(name, 1)[1])):
            add(f"Clapeyron mu/RT difference between phases at {'the Engine' if name.startswith('engine') else 'its own'} solved states (a)", '', ids[i], 'mu/RT', v, 0.)
    ad, cd = (json.loads(next(l for l in log if l.startswith(f'probe7 eps {k}')).split(k, 1)[1]) for k in ('AD', 'FD'))
    for i in range(5): add('Clapeyron permittivity gradient: automatic (a) vs central difference (b)', 'water-rich', ids[i], 'd eps_r / d n_i', ad[i], cd[i])
    write_csv(RES / 'h2h-khudaida-term-checks.csv', rows)

def e_constants(): return json.loads((CC / 'inputs/cases.json').read_text())['constants']['Engine']

if __name__ == '__main__':
    {'engine-timing': engine_timing, 'comparators': comparators, 'probes': probes}[sys.argv[1]](Path(sys.argv[2]))
