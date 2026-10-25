"""Reduce retained raw blocks to timings, host, work, checks and statistics CSVs (no Engine import).

Host label per block: disturbed (outside fit or compiler process at either snapshot), quiet (load1 < 2 at both),
otherwise loaded. Ratios are reported only for value-matched equal work and are ties inside the A=A factor.
"""
import csv
import gzip
import json
import math
from pathlib import Path
import statistics as st

HERE = Path(__file__).resolve().parents[1]
RAW = HERE / 'results/raw'
OUT = HERE / 'results'
INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
TOL = INPUTS['tolerances']
T_OF = {'T1': 'T1', 'T2': 'T2', 'T3': 'T3', 'T4': 'T4', 'T6': 'T6', 'T7': 'T7', 'T8': 'T8', 'T9-serial': 'T9', 'T9-fits': 'T9',
        'T10': 'T10', 'T10-cohorts': 'T10', 'T11': 'T11'}


def blocks(name):
    path = RAW / f'{name}.jsonl.gz'
    if not path.exists():
        return []
    out = []
    try:  # a stopped section leaves an unterminated stream; keep every complete line
        for line in gzip.open(path, 'rt'):
            r = json.loads(line)
            if r['kind'] == 'block':
                out.append(r)
    except (EOFError, json.JSONDecodeError):
        pass
    return out


def host(r):
    snaps = [r['before'], r['after']]
    if any(s.get('outside_fits') or s.get('compilers') for s in snaps):
        return 'disturbed'
    return 'quiet' if all(s['load'][0] < 2.0 for s in snaps) else 'loaded'


def tool(arm):
    return 'Engine' if arm.startswith('Engine') else arm


def write(name, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (gzip.open(OUT / name, 'wt', newline='') if name.endswith('.gz') else (OUT / name).open('w', newline='')) as f:
        w = csv.DictWriter(f, fields, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def close(a, b, rel, floor):
    return abs(a - b) <= rel * max(abs(a), abs(b), floor)


# ---------------------------------------------------------------- value extraction for parity checks

def t1_values(r):
    v = r['reply']['samples'][0]['value']
    if 'P' not in v:  # neutral-property or derivative rows
        return {k: v[k] for k in ('value', 'gradient', 'hessian') if k in v}
    return dict(P=v['P'], h_res=v['h_res'], cv_res=v['cv_res'], a_res=v['a_res'], lnphi=v['lnphi'])


def flat(d, prefix=''):
    out = {}
    items = d.items() if isinstance(d, dict) else ((f'[{i}]', v) for i, v in enumerate(d))
    for k, v in items:
        key = f'{prefix}{k}' if k.startswith('[') else (f'{prefix}.{k}' if prefix else k)
        if isinstance(v, (dict, list)):
            out |= flat(v, key)
        else:
            out[key] = v
    return out


FLOORS = dict(P=1., h_res=1., cv_res=1., a_res=1., lnphi=1., value=1., gradient=1., hessian=1.)


def compare(section, operation, reference, other, label, rel):
    rows = []
    a, b = flat(reference), flat(other)
    for key in a:
        if key not in b:
            continue
        x, y = a[key], b[key]
        if x is None or y is None:
            rows.append(dict(section=section, operation=operation, comparison=label, quantity=key, reference=x, other=y,
                             passed=x is None and y is None, cause='unavailable value on one side' if (x is None) != (y is None) else ''))
            continue
        floor = FLOORS.get(key.split('[')[0], 1.)
        err = abs(x - y) / max(abs(x), abs(y), floor)
        rows.append(dict(section=section, operation=operation, comparison=label, quantity=key, reference=x, other=y,
                         scaled_error=err, tolerance=rel, passed=err <= rel, cause=''))
    return rows


def t3_values(r):
    v = r['reply']['samples'][0]['value']
    arm = tool(r['arm'])
    if arm == 'Engine':
        row = v[0]
        n = len(row['mole_fractions']) // 2
        phases = sorted(range(2), key=lambda j: -row['molar_densities'][j])
        return dict(P=row['pressure'][0], rho=[row['molar_densities'][j] for j in phases],
                    x=[row['mole_fractions'][j * n:(j + 1) * n] for j in phases])
    if arm == 'FeOs':
        ph = sorted(v['phases'], key=lambda p: -p['rho'])
        return dict(P=v['P'], rho=[p['rho'] for p in ph], x=[p['x'] for p in ph])
    order = sorted(range(2), key=lambda j: -v['densities'][j])
    return dict(P=v['P'], rho=[v['densities'][j] for j in order], x=[v['compositions'][j] for j in order])


def eq_compare(section, operation, e, o, label):
    rows = []
    for key, rel in (('P', TOL['equilibrium_rel']),):
        rows.append(dict(section=section, operation=operation, comparison=label, quantity=key, reference=e[key], other=o[key],
                         scaled_error=abs(e[key] - o[key]) / abs(e[key]), tolerance=rel, passed=close(e[key], o[key], rel, 0)))
    for j in range(2):
        rows.append(dict(section=section, operation=operation, comparison=label, quantity=f'rho[{j}]', reference=e['rho'][j],
                         other=o['rho'][j], scaled_error=abs(e['rho'][j] - o['rho'][j]) / e['rho'][j], tolerance=TOL['equilibrium_rel'],
                         passed=close(e['rho'][j], o['rho'][j], TOL['equilibrium_rel'], 0)))
        for i, (a, b) in enumerate(zip(e['x'][j], o['x'][j])):
            rows.append(dict(section=section, operation=operation, comparison=label, quantity=f'x[{j}][{i}]', reference=a, other=b,
                             scaled_error=abs(a - b), tolerance=TOL['composition_abs'], passed=abs(a - b) <= TOL['composition_abs']))
    return rows


def check_row(section, operation, comparison, quantity, reference, other, error, tolerance, cause=''):
    return dict(section=section, operation=operation, comparison=comparison, quantity=quantity, reference=reference, other=other,
                scaled_error=error, tolerance=tolerance, passed=error is not None and error <= tolerance, cause=cause)


def t2_checks(section, operation, ref, r):
    """Engine roots against the retained density oracle; comparators against it or against the Engine root."""
    out, arm = [], tool(r['arm'])
    value = r['reply']['samples'][0]['value']
    case = next((c for c in INPUTS['t2']['ordinary'] if operation.startswith(c['id'] + '/')), None)
    if case and operation.endswith('/pressure-cold'):
        for st, v in zip(case['states'], value):
            rho = v.get('rho') if arm == 'Engine' else v.get('rho_mol_m3')
            err = None if rho is None else abs(rho / st['oracle_rho'] - 1)
            out.append(check_row(section, operation, f'{arm} vs density oracle', f'rho[{st["index"]}]', st['oracle_rho'], rho, err,
                           TOL['density_rel'], '' if rho is not None else v.get('note') or v.get('status', 'no root')))
    elif operation.endswith('/pressure-cold') and arm != 'Engine' and ref is not None:
        for j, (a, b) in enumerate(zip(ref['reply']['samples'][0]['value'], value)):
            rho = b.get('rho')
            err = None if rho is None else abs(rho / a['rho'] - 1)
            out.append(check_row(section, operation, f'Engine vs {arm}', f'rho[{j}]', a['rho'], rho, err, TOL['density_rel'], b.get('error', '')))
    return out


def phases_of(r):
    v, arm = r['reply']['samples'][0]['value'], tool(r['arm'])
    if arm == 'Engine':
        e = v[0]
        if 'error' in e or not e.get('success'):
            return 'failed'
        n = len(e['mole_fractions']) // e['phases']
        return sorted(((e['molar_densities'][j], e['mole_fractions'][n * j:n * j + n]) for j in range(e['phases'])), reverse=True)
    if arm == 'FeOs':
        if v.get('phases') is None:
            return 'failed'
        return None if v.get('phases') == 1 else sorted(zip(v['densities'], v['compositions']), reverse=True)
    present = [j for j, f in enumerate(v['fractions']) if f > 1e-12]
    return sorted(((1 / v['volumes'][j], v['compositions'][j]) for j in present), reverse=True)


def t4_checks(section, operation, ref, r):
    if not operation.endswith('/discover') or 'sequence' in operation:
        return []
    a, b = phases_of(ref), phases_of(r)
    label = f'Engine vs {tool(r["arm"])}' if not r['arm'].startswith('Engine') else 'same image'
    if a == 'failed':
        return [check_row(section, operation, label, 'phase count', None, None, None, 0.0, 'Engine discovery returned no accepted state')]
    if b == 'failed':
        return [check_row(section, operation, label, 'phase count', None, None, None, 0.0, r['reply']['samples'][0]['value']['error'][:200])]
    na, nb = len(a) if a else 1, len(b) if b else 1
    out = [check_row(section, operation, label, 'phase count', na, nb, float(na != nb), 0.0,
               '' if na == nb else 'different number of phases returned')]
    if a and b and na == nb:
        for j, ((ra, xa), (rb, xb)) in enumerate(zip(a, b)):
            out.append(check_row(section, operation, label, f'rho[{j}]', ra, rb, abs(rb / ra - 1), TOL['equilibrium_rel']))
            out.append(check_row(section, operation, label, f'x[{j}]', xa, xb, max(abs(p - q) for p, q in zip(xa, xb)), TOL['composition_abs']))
    return out


def t6_checks(section, ref_rows, other_rows, label, operation):
    out = []
    for a, b in zip(ref_rows, other_rows):
        out.append(check_row(section, operation, label, f'gamma(m={a["m"]:.6g})', a['gamma'], b['gamma'],
                       abs(a['gamma'] - b['gamma']) / max(abs(a['gamma']), 1.), TOL['property_rel']))
    return out


def node_error(a, b):
    worst = 0.0
    for ra, rb in zip(a, b):
        for key in ('value', 'gradient', 'hessian'):
            x, y = ra.get(key), rb.get(key)
            if x is None or y is None:
                continue
            x, y = (x, y) if isinstance(x, list) else ([x], [y])
            scale = max([abs(v) for v in x + y] + [1e-300])
            worst = max(worst, max(abs(p - q) for p, q in zip(x, y)) / scale)
    return worst


def main():
    timing, work, checks = [], [], []
    for name in ('T7', 'T3', 'T1', 'T2', 'T4', 'T6', 'T8', 'T9-serial', 'T9-fits', 'T10', 'T10-cohorts', 'T11',
                 'T7-disturbed-amine-fit'):
        rows = blocks(name)
        first = {}
        for r in rows:
            label = host(r)
            reply = r['reply']
            base = dict(section=name, T=T_OF.get(name, name.split('-')[0]), operation=r['operation'], tool=tool(r['arm']), arm=r['arm'],
                        repetition=r['repetition'], ordinal=r.get('ordinal'), host=label, load1_before=r['before']['load'][0],
                        load1_after=r['after']['load'][0], outside_fits=len(r['before'].get('outside_fits', [])) + len(r['after'].get('outside_fits', [])),
                        mem_available_kb=r['before']['mem_available_kb'])
            if 'error' in reply:
                timing.append(base | dict(sample=None, error=reply['error'][:300]))
                continue
            for i, s in enumerate(reply['samples']):
                err = s.get('error') or (s['value'].get('error') if isinstance(s.get('value'), dict) else None)
                calls = s.get('per_call') or 1
                timing.append(base | dict(sample=i, seconds=s['seconds'], per_call=calls, seconds_per_call=s['seconds'] / calls,
                                          cpu_seconds=s.get('cpu_seconds'), child_cpu_seconds=s.get('child_cpu_seconds'),
                                          reference_seconds=s.get('reference_seconds'), curve_seconds=s.get('curve_seconds'),
                                          nodes_per_second=s.get('nodes_per_second'), error=err or ''))
                counts = s.get('work') or s.get('counters')
                if counts:
                    work.append(base | dict(sample=i) | {f'work_{k}': v for k, v in counts.items()})
            # same-image and cross-tool value parity, from the first sample of each arm's first repetition
            if r['repetition'] == 0 and reply.get('samples'):
                first[(r['operation'], r['arm'])] = r
        for (operation, arm), r in first.items():
            if arm == 'Engine-A':
                continue
            ref = first.get((operation, 'Engine-A'))
            if ref is None:
                continue
            if name == 'T1':
                label = 'same image' if arm == 'Engine-control' else f'Engine vs {arm}'
                rel = TOL['same_image'] if arm == 'Engine-control' else TOL['property_rel']
                checks += compare(name, operation, t1_values(ref), t1_values(r), label, rel)
            elif name == 'T3':
                if arm == 'Engine-control':
                    checks += compare(name, operation, t3_values(ref), t3_values(r), 'same image', TOL['same_image'])
                else:
                    checks += eq_compare(name, operation, t3_values(ref), t3_values(r), f'Engine vs {arm}')
            elif name == 'T2' and (arm != 'Engine-control' or True):
                checks += t2_checks(name, operation, ref, r) if arm != 'Engine-control' else []
                if arm == 'Engine-control':
                    a, b = json.dumps(ref['reply']['samples'][0].get('value'), sort_keys=True), json.dumps(r['reply']['samples'][0].get('value'), sort_keys=True)
                    checks.append(dict(section=name, operation=operation, comparison='same image', quantity='complete returned value',
                                       passed=a == b, cause='' if a == b else 'returned values differ between same-image workers'))
            elif name == 'T4':
                checks += t4_checks(name, operation, ref, r)
            elif arm == 'Engine-control':
                a, b = json.dumps(ref['reply']['samples'][0].get('value'), sort_keys=True), json.dumps(r['reply']['samples'][0].get('value'), sort_keys=True)
                checks.append(dict(section=name, operation=operation, comparison='same image', quantity='complete returned value',
                                   passed=a == b, cause='' if a == b else 'returned values differ between same-image workers'))
        if name == 'T2':  # the Engine-A oracle rows are not reached by the pairwise loop above
            for (operation, arm), r in first.items():
                if arm == 'Engine-A':
                    checks += t2_checks(name, operation, None, r)
        if name == 'T6':
            e = first.get(('nacl/reused', 'Engine-A'))
            f = first.get(('nacl/fresh', 'Engine-A'))
            c = first.get(('nacl/clapeyron-curve', 'Clapeyron'))
            if e and f:
                checks += t6_checks(name, e['reply']['samples'][0]['value']['points'], f['reply']['samples'][0]['value']['points'], 'reused vs fresh reference', 'nacl/fresh')
            if e and c:
                rows_c = [dict(m=m, gamma=v['gamma']) for m, v in zip(INPUTS['t6']['molalities'], c['reply']['samples'][0]['value'])]
                checks += t6_checks(name, e['reply']['samples'][0]['value']['points'], rows_c, 'Engine vs Clapeyron', 'nacl/clapeyron-curve')
        if name == 'T10':
            for (operation, arm), r in first.items():
                if arm == 'Engine-A' and not operation.endswith('/native') and (operation + '/native', 'Engine-A') in first:
                    a = r['reply']['samples'][0]['value']
                    b = first[(operation + '/native', 'Engine-A')]['reply']['samples'][0]['value']
                    err = max(node_error(x, y) for x, y in zip(a, b))
                    checks.append(check_row(name, operation, 'ASL vs native', 'value, gradient and Hessian entries', None, None, err, TOL['same_image']))
    verification = []
    for name in ('T1-verify', 'T8-verify'):
        for r in blocks(name):
            for row in r['reply']['samples'][0]['value']:
                levels = row.get('levels') or []
                errs = [lv.get('error', lv.get('worst_excess')) for lv in levels]
                verification.append(dict(
                    section=name, state=r['operation'], order=row.get('order', len(row.get('combo', []))),
                    directions='+'.join(row['names']) if 'names' in row else '+'.join(map(str, row.get('combo', []))),
                    check=row.get('check', 'stencil' if 'names' in row else ''), verdict=row.get('verdict', 'passed' if row.get('passed', True) else 'failed'),
                    cause=row.get('cause', ''), action=row.get('action'),
                    coarse_error=errs[0] if errs else row.get('error', row.get('scaled_error', row.get('scaled_difference'))),
                    fine_error=errs[1] if len(errs) > 1 else None,
                    allowed=levels[0].get('allowed') if levels else row.get('allowed'),
                    halving_ratio=(errs[0] / errs[1]) if len(errs) > 1 and errs[1] not in (None, 0) and errs[0] is not None
                    and errs[0] > 0 and errs[1] > 0 else None))
                ratio = verification[-1]['halving_ratio']
                if verification[-1]['verdict'] == 'unresolved verification':
                    verification[-1]['diagnosis'] = ('step truncation: the error falls about fourfold when the step halves'
                                                     if ratio and 3.5 <= ratio <= 4.5 else 'no convergence on halving: cancellation or round-off')
    write('verification.csv', verification)
    write('timings.csv.gz', timing)
    write('work.csv.gz', work)
    write('checks.csv', checks)
    stats = []
    groups = {}
    for t in timing:
        if t.get('sample') is None or t['section'] == 'T7-disturbed-amine-fit':
            continue
        groups.setdefault((t['section'], t['T'], t['operation'], t['tool']), {}).setdefault((t['arm'], t['repetition']), []).append(t)
    failed = {(c['section'], c['operation'], c['comparison'].removeprefix('Engine vs ')) for c in checks
              if not c['passed'] and c['comparison'] != 'same image'}
    for (section, T, operation, tl), by_block in groups.items():
        medians = [st.median(x['seconds_per_call'] for x in v) for v in by_block.values()]
        samples = sorted(x['seconds_per_call'] for v in by_block.values() for x in v)
        labels = sorted({x['host'] for v in by_block.values() for x in v})
        row = dict(section=section, T=T, operation=operation, tool=tl, blocks=len(medians), samples=len(samples),
                   median_s=st.median(medians), min_block_s=min(medians), max_block_s=max(medians),
                   spread_factor=max(medians) / min(medians), mad_s=st.median(abs(m - st.median(medians)) for m in medians),
                   p10_s=samples[int(0.1 * (len(samples) - 1))], p90_s=samples[int(0.9 * (len(samples) - 1))], hosts='+'.join(labels),
                   load1_min=min(x['load1_before'] for v in by_block.values() for x in v),
                   load1_max=max(x['load1_before'] for v in by_block.values() for x in v),
                   errors=sum(1 for v in by_block.values() for x in v if x.get('error')))
        stats.append(row)
    engine = {(s['section'], s['operation']): s for s in stats if s['tool'] == 'Engine'}
    for s in stats:
        e = engine.get((s['section'], s['operation']))
        if s['tool'] != 'Engine' and e:
            ratio = s['median_s'] / e['median_s']
            s.update(engine_median_s=e['median_s'], ratio_to_engine=ratio, engine_aa_factor=e['spread_factor'],
                     tie=1 / e['spread_factor'] <= ratio <= e['spread_factor'],
                     value_matched=(s['section'], s['operation'], s['tool']) not in failed)
    write('statistics.csv', stats)
    print(len(timing), 'timing rows;', len(work), 'work rows;', len(checks), 'checks,', sum(1 for c in checks if not c['passed']), 'failed;',
          len(stats), 'statistics rows')


if __name__ == '__main__':
    main()
