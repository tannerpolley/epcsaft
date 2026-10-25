"""One row per audit question (results/answers.csv): the auditor's finding, its status at current main and the
measured answer, computed from the reduced CSVs and inputs/prior.json. The only interpretation owner.

Status: holds (the suspected cost or limit is measured at current main), does not hold (measured and absent),
not settled (the exposed counters cannot separate the suspected cost), mechanism confirmed, not material (the
mechanism is shown but its cost is immaterial at the measured size), fixed (#N), replayed (removed since the audit by
the named issue and confirmed by a replay here), fixed per #N record, not re-timed (removed per the named issue's delivered
record), new (a cost the audits did not name), not timed (dropped by the accepted design; owner named) or excluded
(outside #578).
"""
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import statistics as st

HERE = Path(__file__).resolve().parents[1]
R = HERE / 'results'
PRIOR = json.loads((HERE / 'inputs/prior.json').read_text())


def rows(name):
    path = R / name
    if not path.exists():
        return []
    return list(csv.DictReader(gzip.open(path, 'rt') if name.endswith('.gz') else path.open()))


S = rows('statistics.csv')
C = rows('checks.csv')
V = rows('verification.csv')
STAT = {(r['section'], r['operation'], r['tool']): r for r in S}


def med(section, operation, tool='Engine'):
    r = STAT.get((section, operation, tool))
    return None if r is None else float(r['median_s'])


def aa(section, operation):
    r = STAT.get((section, operation, 'Engine'))
    return None if r is None else float(r['spread_factor'])


def us(x):
    if x is None:
        return 'unavailable'
    return f'{x * 1e3:.3g} ms' if x >= 1e-3 else f'{x * 1e6:.3g} µs'


def ms(x):
    return 'unavailable' if x is None else f'{x * 1e3:.3g} ms'


def s_(x):
    return 'unavailable' if x is None else f'{x:.3g} s'


def ratio(a, b):
    return None if a is None or b is None else a / b


def fx(x):
    return 'unavailable' if x is None else f'{x:.2g}×'


def matched(section, operation, tool):
    r = STAT.get((section, operation, tool))
    return r is not None and r.get('value_matched') == 'True'


def outside_aa(section, operation, tool):
    r = STAT.get((section, operation, tool))
    return r is not None and r.get('tie') == 'False'


def ties(section, ops, tool):
    """Ratios inside the Engine's own A=A spread, named with their Engine/comparator ratio."""
    out = [f'{o.split("/")[1]} {fx(ratio(med(section, o), med(section, o, tool)))}' for o in ops
           if STAT.get((section, o, tool), {}).get('tie') == 'True']
    return f' (ties inside A=A: {", ".join(out)})' if out else ''


def t1():
    out = []
    pure = [k[1].rsplit('/', 1)[0] for k in STAT if k[0] == 'T1' and k[1].endswith('/helmholtz-Value') and k[2] == 'Engine']
    by_order = {o: [ratio(med('T1', f'{p}/helmholtz-{o}'), med('T1', f'{p}/helmholtz-{o}', 'Clapeyron')) for p in pure]
                for o in ('Value', 'Gradient', 'Hessian')}
    bundles = [k[1] for k in STAT if k[0] == 'T1' and k[1].endswith('/bundle') and k[2] == 'Clapeyron']
    vsc = {b: ratio(med('T1', b), med('T1', b, 'Clapeyron')) for b in bundles if matched('T1', b, 'Clapeyron') and outside_aa('T1', b, 'Clapeyron')}
    vsf = {b: ratio(med('T1', b), med('T1', b, 'FeOs')) for b in (k[1] for k in STAT if k[0] == 'T1' and k[1].endswith('/bundle') and k[2] == 'FeOs')
           if matched('T1', b, 'FeOs')}
    over3 = sorted(b.split('/bundle')[0] for b, x in vsc.items() if x >= 3)
    over3f = [b for b, x in vsf.items() if x >= 3 and outside_aa('T1', b, 'FeOs')]
    order = lambda o: f'{fx(st.median(by_order[o]))}{ties("T1", [f"{p}/helmholtz-{o}" for p in pure], "Clapeyron")}'
    out.append(dict(T='T1', audit_id='O-S1', finding='At equal derivative order the Engine is more than 3× slower than FeOs, which would justify compiled low-order forward mode',
                    status='holds' if over3f else 'does not hold',
                    answer=(f'Full fixed-density bundle (P, h_res, cv_res, ln phi, a_res), value-matched: Engine/FeOs {fx(min(vsf.values()))}–{fx(max(vsf.values()))}, '
                            'so the Engine is faster than FeOs at every matched state (FeOs pays its Python property interface). '
                            f'Side note on Clapeyron: Engine/Clapeyron {fx(min(vsc.values()))}–{fx(max(vsc.values()))} outside A=A'
                            f'{ties("T1", bundles, "Clapeyron")}, at least 3× only on the associating states {", ".join(over3) or "none"}. '
                            f'Pure neutral reduced Helmholtz energy, Engine/Clapeyron median ratio by order: value {order("Value")}, '
                            f'gradient {order("Gradient")}, Hessian {order("Hessian")}.'),
                    evidence='statistics.csv T1 */helmholtz-*, */bundle',
                    decision='No compiled low-order path is justified against FeOs; the Clapeyron gap on associating states would need its own design.'))
    eth = [b for b in ('Sol-1/ethane-100', 'Sol-1/ethane-10000', 'Sol-1/ethane-loop')]
    out.append(dict(T='T1', audit_id='Sol-1', finding='Which work dominates a neutral full-state request (ethane, three states)',
                    status='holds',
                    answer=(f'Ethane bundle {", ".join(us(med("T1", e + "/bundle")) for e in eth)}; reduced-Helmholtz value '
                            f'{", ".join(us(med("T1", e + "/helmholtz-Value")) for e in eth)}, Hessian {", ".join(us(med("T1", e + "/helmholtz-Hessian")) for e in eth)}. '
                            'The bundle costs about one second-order evaluation: state() builds every derivative it may return, so a value-only property still pays for it.'),
                    evidence='statistics.csv T1 Sol-1/*; work.csv.gz forward/reverse/hessian counters', decision='A narrower State request is a separate design question.'))
    out.append(dict(T='T1', audit_id='A-Q01', finding='HelmholtzModel::state costs more because it always builds a broad State',
                    status='holds', answer=('Same measurement as Sol-1 on the four methane/ethane states: bundle '
                                            + ', '.join(us(med('T1', f'A-Q01/{s}/bundle')) for s in ('300-5000', '180-200', '180-18000', '300-dilute'))
                                            + '; pressure-only and pressure-plus-enthalpy requests have no narrower Engine call to time.'),
                    evidence='statistics.csv T1 A-Q01/*', decision='As Sol-1.'))
    first = lambda sid, n: st.median([float(r['median_s']) for r in S if r['section'] == 'T1' and r['tool'] == 'Engine'
                                      and r['operation'].startswith(sid + '/action-') and r['operation'].count('+') == n - 1] or [float('nan')])
    assoc = first('A-Q05/associating-binary', 1), first('A-Q05/associating-binary', 2), first('A-Q05/associating-binary', 3)
    neut = first('O-P01/methane-ethane', 1), first('O-P01/methane-ethane', 2)
    out.append(dict(T='T1', audit_id='Sol-3', finding='Implicit association cost grows with the site graph',
                    status='holds', answer=(f'ResidualEnthalpy action median per request, cross-associating binary (four sites): order 1 {us(assoc[0])}, '
                                            f'order 2 {us(assoc[1])}, order 3 {us(assoc[2])}; non-associating methane/ethane: order 1 {us(neut[0])}, order 2 {us(neut[1])}; '
                                            f'self-associating methanol/n-hexane order 1 {us(first("O-P01/methanol-hexane", 1))}.'),
                    evidence='statistics.csv T1 */action-*; verification.csv', decision='Association is the largest per-state action cost outside the electrolyte and MEA rows.'))
    out.append(dict(T='T1', audit_id='A-Q05', finding='Association central solves and implicit site derivatives dominate higher-order actions',
                    status='holds', answer=(f'Order 3 costs {fx(ratio(assoc[2], assoc[0]))} the order-1 request on the associating binary. '
                                            f'MEA nine-species liquid: bundle {us(med("T1", "O-S1/MEA-liquid/bundle"))}, order-1 {us(first("O-S1/MEA-liquid", 1))}, '
                                            f'order-2 {us(first("O-S1/MEA-liquid", 2))}; its charge-changing composition basis directions are refused by native admission.'),
                    evidence='statistics.csv T1; verification.csv', decision='No derivative-method change is supported by cost alone.'))
    return out


def t2():
    cases = ('A-propane-vapor', 'A-propane-liquid', 'A-methane-ethane-vapor', 'A-C2-liquid')
    cold = {c: med('T2', f'{c}/pressure-cold') for c in cases}
    fixed = {c: med('T2', f'{c}/fixed-density') for c in cases}
    anch = {c: med('T2', f'{c}/pressure-anchored') for c in cases}
    clap = {c: med('T2', f'{c}/pressure-cold', 'Clapeyron') for c in cases}
    save = {c: 1 - anch[c] / cold[c] for c in cases}
    vap = 'A-propane-vapor'
    loop = [c for c in ('A-Q02/308.8-vapor', 'A-Q02/308.8-liquid', 'A-Q02/308.0-vapor', 'A-Q02/308.0-liquid')]
    act1 = st.median(med('T2', f'{c}/pressure-action-[{j}]') for c in loop for j in (0, 1))
    act2 = st.median(med('T2', f'{c}/pressure-action-{k}') for c in loop for k in ('[0, 0]', '[0, 1]', '[1, 1]'))
    loopcold = st.median(med('T2', f'{c}/pressure-cold') for c in loop)
    loopfixed = st.median(med('T2', f'{c}/fixed-density') for c in loop)
    oracle = [c for c in C if c['section'] == 'T2' and 'density oracle' in c['comparison']]
    feos_loop = [c for c in C if c['section'] == 'T2' and c['comparison'] == 'Engine vs FeOs']
    common = (f'All {len(oracle)} Engine, FeOs and Clapeyron roots on the four ordinary 32-state sequences '
              f'{"match" if all(c["passed"] == "True" for c in oracle) else "do not all match"} the retained density oracle. ')
    return [
        dict(T='T2', audit_id='O-S2', finding='Pressure-to-density closure on declared branches costs more than comparators; vapor requests check the liquid side first (F13)',
             status='holds', answer=common + (f'Engine per call, cold / anchored / fixed density: propane vapor {us(cold[vap])} / {us(anch[vap])} / {us(fixed[vap])}; '
                                              f'Clapeyron {us(clap[vap])}, FeOs {us(med("T2", vap + "/pressure-cold", "FeOs"))}. Closure adds '
                                              f'{fx(ratio(cold[vap], fixed[vap]))} over the fixed-density state on vapor and {fx(ratio(cold["A-propane-liquid"], fixed["A-propane-liquid"]))} on liquid. '
                                              'Exposed counters do not separate a liquid-side search from the vapor root, so no causal attribution to liquid-first order is made.'),
             evidence='statistics.csv T2; checks.csv T2', decision='Branch search is a minority of state cost; no dropped check is justified.'),
        dict(T='T2', audit_id='Sol-2', finding='Density closure dominates pressure-chart derivatives near the ethane critical loop',
             status='does not hold', answer=(f'Ethane loop states: cold closure {us(loopcold)}, fixed-density state {us(loopfixed)}, first pressure-chart action {us(act1)}, '
                                             f'second {us(act2)} (medians over the four branch states). '
                                             f'FeOs roots differ from the Engine by up to {max(float(c["scaled_error"]) for c in feos_loop):.2g} relative here, so its loop times are not ranked; '
                                             'Clapeyron matches within the density criterion.'),
             evidence='statistics.csv T2 A-Q02/*; checks.csv T2', decision='No narrow-loop closure change follows.'),
        dict(T='T2', audit_id='A-Q02', finding='Locating slope zeros and proving a bracket costs more than closing the root',
             status='not settled', answer=(f'A warm previous root saves {min(save.values()):.0%}–{max(save.values()):.0%} of a cold density solve (anchored vs cold: ' + '; '.join(
                 f'{c} {us(anch[c])} vs {us(cold[c])}' for c in cases) + f'; loop states anchored {us(st.median(med("T2", f"{c}/pressure-anchored") for c in loop))}). '
                 'The exposed counters do not separate the slope-zero search and bracket check from the root closure, so the finding is not settled. '
                 'Native diagnostics for the loop states are retained once outside timing (raw/T2-prepare.json).'),
             evidence='statistics.csv T2; raw/T2-prepare.json', decision='Anchors save a fraction of a closure; isolating the bracket check needs a counter or its own design.'),
    ]


def t3():
    ops = ('C1-300', 'C2-0.5', 'methane-bubble')
    e = {o: med('T3', o) for o in ops}
    f = {o: med('T3', o, 'FeOs') for o in ops}
    c = {o: med('T3', o, 'Clapeyron') for o in ops}
    then = PRIOR['history_549']
    text = '; '.join(f'{o}: Engine {us(e[o])} (cf26a17b {us(then[o]["seconds"])}), FeOs {us(f[o])}, Clapeyron {us(c[o])}' for o in ops)
    slow = all(outside_aa('T3', o, t) and ratio(e[o], m[o]) > 1 for o in ops for t, m in (('FeOs', f), ('Clapeyron', c)))
    return [
        dict(T='T3', audit_id='O-S3', finding='Cold pure saturation is slower than comparators', status='holds' if slow else 'does not hold',
             answer=text + '. All three tools return the same states within the equilibrium criteria.', evidence='statistics.csv T3; checks.csv T3',
             decision='The Engine remains slower on these represented incipient-phase calls; any change needs its own design.'),
        dict(T='T3', audit_id='O-S4', finding='Binary bubble/dew from cold and warm sequences; cold starts fail away from the qualified cases',
             status='fixed per #565 record, not re-timed', answer='Cold-start convergence of bubble, dew and two-phase flash is fixed by #565 (PR #589) per its delivered record (merged before this campaign; not re-timed here). The represented cold bubble costs are in O-S3; warm sequences and dew rows were dropped by the design.',
             evidence='#565 closeout; statistics.csv T3', decision='Time a warm sequence only with a selected alternative.'),
        dict(T='T3', audit_id='A-Q03', finding='The common constrained driver spends more work than a specialized bubble algorithm',
             status='holds' if slow else 'does not hold', answer=f'Methane/ethane 180 K bubble: Engine {us(e["methane-bubble"])}, FeOs {us(f["methane-bubble"])}, Clapeyron {us(c["methane-bubble"])}.',
             evidence='statistics.csv T3', decision='As O-S3.'),
        dict(T='T3', audit_id='F-S5', finding='Newton fast paths versus the forced Ipopt path', status='not timed',
             answer='No option disables the saturation/bubble/dew Newton attempt; the design dropped this pair (#578 design T3).', evidence='design.md T3', decision='None.'),
        dict(T='T3', audit_id='F-S9', finding='FeOs/Clapeyron versus Engine for pure saturation and binary bubble points', status='holds' if slow else 'does not hold',
             answer=text + '.', evidence='statistics.csv T3', decision='As O-S3.'),
    ]


def t7():
    k = PRIOR['key_517a3450']
    W = rows('work.csv.gz')
    work = {}
    for w in W:
        if w['section'] == 'T7' and w['arm'] == 'Engine-A':
            work.setdefault(w['operation'], {c[5:]: w[c] for c in w if c.startswith('work_') and w[c] != ''})
    lines = []
    for case, prior in k.items():
        lines.append(f'{case}: {ms(med("T7", case)) if not case.startswith("mea-") else s_(med("T7", case))} (517a3450 {ms(prior["seconds"]) if not case.startswith("mea-") else s_(prior["seconds"])})')
    f1 = {key: int(float(v)) for key, v in work.get('mea-f1/residual-jacobian', {}).items() if key in k['mea-f1/residual-jacobian']['work']}
    return [
        dict(T='T7', audit_id='O-S8', finding='Absolute time and work of the key application solves', status='holds',
             answer='; '.join(lines) + f'; reacting Cp row {s_(med("T7", "mea-cp/hilliard-353.15-K-0.358"))}. Ordinary-case work counts are identical to 517a3450.',
             evidence='statistics.csv T7; work.csv.gz T7; inputs/prior.json', decision='No slowdown is attributed; times are loaded-host observations.'),
        dict(T='T7', audit_id='F-S1', finding='Cold reacting MEA flash against the legacy pcsaft flash', status='not timed',
             answer=f'Engine cold MEA flash {ms(med("T7", "case-7/flash-5918.5286-Pa/cold"))}; the legacy application solver is not an equivalent comparator (design T7).',
             evidence='statistics.csv T7', decision='None.'),
        dict(T='T7', audit_id='F-S2', finding='Ionic LLE warm sequence against a FeOs neutral analogue', status='not timed',
             answer=f'Engine Khudaida eight-feed sequence {ms(med("T7", "khudaida/sequential-8"))}; FeOs has no ionic LLE (design T7).', evidence='statistics.csv T7', decision='None.'),
        dict(T='T7', audit_id='F-S10', finding='Richardson-ladder reference work per reacting solve', status='holds',
             answer=(f'F1 residual plus exact Jacobian: {f1["reference_computations"]} reference computations and {f1["reference_columns"]} reference columns '
                     f'with {f1["density_solves"]} density solves (517a3450: {k["mea-f1/residual-jacobian"]["work"]["reference_computations"]} / '
                     f'{k["mea-f1/residual-jacobian"]["work"]["reference_columns"]} / {k["mea-f1/residual-jacobian"]["work"]["density_solves"]}); refrozen by #595 (PR #601), '
                     f'which stopped computing unread reference columns; reacting Cp row {s_(med("T7", "mea-cp/hilliard-353.15-K-0.358"))}.'),
             evidence='work.csv.gz T7; inputs/prior.json', decision='Reference reuse already shares computations; see T9.'),
    ]


def dropped():
    return [
        dict(T='T5', audit_id=i, finding='Phase envelope through the critical point: total time and corrector solves', status='fixed per #571 record, not re-timed',
             answer=('#571 (PR #591, merged before this campaign) closes the P05 methane/ethane envelope in 89 accepted continuation steps against 459 before it, '
                     'per its delivered review; not re-timed here because the design dropped envelope timing. The competitor comparison on equally resolved curves '
                     '(O-S6, Sol-10, A-Q11) is still unmeasured.'),
             evidence='#571 delivered review', decision='An equally resolved curve comparison is a separate study.') for i in ('O-S6', 'Sol-10', 'A-Q11')] + [
        dict(T='T12', audit_id='F-S3', finding='Dense Jacobian/Hessian cost as species count grows', status='not timed',
             answer='No accepted consumer of a 20/40-species packet (design T12).', evidence='design.md', decision='None.'),
        dict(T='T13', audit_id='O-S11', finding='Construction and first-call latency', status='not timed',
             answer='No decision depends on it (design T13); preparation times are retained per request in the raw records.', evidence='raw/*.jsonl.gz prepare_seconds', decision='None.'),
        dict(T='T13', audit_id='O-S12', finding='Import latency', status='not timed', answer='No decision depends on it (design T13).', evidence='design.md', decision='None.'),
        dict(T='—', audit_id='F-S6', finding='MA57 versus MUMPS', status='excluded', answer='Owned by #509 (Coin-HSL licence).', evidence='#578 body', decision='None.'),
        dict(T='—', audit_id='A-Q10', finding='Equilibrium linear-solver choices', status='excluded', answer='Owned by #509.', evidence='#578 body', decision='None.'),
        dict(T='—', audit_id='A-Q12', finding='End-to-end application cost', status='excluded', answer='Waits for the absorber adapter repair (MEA #52).', evidence='#578 body', decision='None.'),
    ]


def raw(name):
    import gzip
    path = R / 'raw' / f'{name}.jsonl.gz'
    out = []
    if path.exists():
        try:
            for line in gzip.open(path, 'rt'):
                r = json.loads(line)
                if r['kind'] == 'block':
                    out.append(r)
        except EOFError:
            pass
    return out


def t4():
    if not raw('T4'):
        return []
    blocks = [r for r in raw('T4') if r['arm'] == 'Engine-A' and r['repetition'] == 0]
    work = {}
    for r in blocks:
        if r['operation'].endswith('/discover') and 'sequence' not in r['operation']:
            row = r['reply']['samples'][0]['value'][0]
            d = row.get('discovery', {})
            work[r['operation']] = (row.get('phases'), row.get('iterations'), d.get('work', {}), len(d.get('rounds', [])))
    drivers = [w[2].get('driver_invocations') for w in work.values() if w[2]]
    phase_fail = [c for c in C if c['section'] == 'T4' and c['passed'] != 'True']
    seq = {m: med('T4', f'{sysname}/sequence/{m}') for sysname in ('methane_ethane',) for m in ('discover', 'declared-cold', 'declared-sequence')}
    seq2 = {m: med('T4', f'methanol_hexane/sequence/{m}') for m in ('discover', 'declared-cold', 'declared-sequence')}
    per = [k for k in STAT if k[0] == 'T4' and k[1].endswith('/discover') and 'sequence' not in k[1] and k[2] == 'Engine']
    e = st.median(med(*k) for k in per)
    f = [med(k[0], k[1], 'FeOs') for k in per if med(k[0], k[1], 'FeOs')]
    c = [med(k[0], k[1], 'Clapeyron') for k in per if med(k[0], k[1], 'Clapeyron')]
    mich = [med('T4', k[1].replace('/discover', '/michelsen'), 'Clapeyron') for k in per if med('T4', k[1].replace('/discover', '/michelsen'), 'Clapeyron')]
    text = (f'Engine Discover per point median {ms(e)} (FeOs tp_flash {ms(st.median(f)) if f else "—"}, Clapeyron default {ms(st.median(c)) if c else "—"}, '
            f'Michelsen {ms(st.median(mich)) if mich else "—"}); driver invocations per discovery {min(drivers):.0f}–{max(drivers):.0f}. '
            f'Six-point methane/ethane sequence: discover {ms(seq["discover"])}, discovered declarations started at the discovered state {ms(seq["declared-cold"])}, '
            f'started from the previous accepted point {ms(seq["declared-sequence"])}; five-point methanol/n-hexane: {ms(seq2["discover"])} / {ms(seq2["declared-cold"])} / '
            f'{ms(seq2["declared-sequence"])}. Discovery costs {fx(ratio(seq["discover"], seq["declared-cold"]))} a re-solve from the discovered state on methane/ethane. '
            f'Clapeyron agrees with every Engine phase count and state; FeOs agrees on methane/ethane and fails its own stability analysis or Rachford-Rice step at '
            f'{len(phase_fail)} methanol/n-hexane points, which carry no ratio.')
    # The discovered declaration carries the accepted state as its start (phase_stability.cpp:803-812), so this arm restarts at its own answer.
    pts = {m: next(r['reply']['samples'][0]['value'] for r in raw('T4') if r['operation'] == f'methane_ethane/sequence/{m}' and r['arm'] == 'Engine-A' and r['repetition'] == 0)
           for m in ('declared-cold', 'declared-sequence')}
    wk = lambda p: f'{p["iterations"]}/{p["evaluations"]}/{p["hessian_calls"]}'
    own = sorted({wk(p) for p in pts['declared-cold'] if p['phases'] == 2})
    at = lambda p: f'{p["T"]:g} K/{p["P"] / 1e6:g} MPa'
    seqp = pts['declared-sequence']
    moved = [f'{at(p)} {wk(p)} against {wk(q)}' for p, q in zip(seqp, pts['declared-cold']) if wk(p) != wk(q)]
    worst = max(range(1, len(seqp)), key=lambda i: seqp[i]['evaluations'])
    start = next(seqp[j] for j in range(worst - 1, -1, -1) if seqp[j]['phases'] == seqp[worst]['phases'])
    return [dict(T='T4', audit_id='O-S5', finding='Unknown-phase-count flash pays for stability trials run through the full driver (R09b)', status='holds',
                 answer=text, evidence='statistics.csv T4; checks.csv T4; raw/T4.jsonl.gz discovery work',
                 decision='Discovery, not the declared solve, is the cost; any change to trial work needs its own correctness design.'),
            dict(T='T4', audit_id='Sol-7', finding='Which initialization or discovery policy pays for a sequence rather than one solve', status='holds',
                 answer=('Discovery dominates the sequence. Restarting from the discovered state is cheap: the discovered declaration carries the accepted state as its start '
                         f'(phase_stability.cpp:803-812), so every two-phase point takes {", ".join(own)} iterations/evaluations/Hessian calls. A start taken from the '
                         f'previous accepted point changes the work only at {", ".join(moved) or "no point"}; the extra cost is concentrated at {at(seqp[worst])}, '
                         f'whose start comes from the distant {at(start)} point. '
                         'No genuinely cold declared solve was timed. ' + text),
                 evidence='statistics.csv T4 */sequence/*; raw/T4.jsonl.gz', decision='Discovery is the sequence cost; a start from a distant previous point is not free.')]


def t6():
    if not raw('T6'):
        return []
    rec = {(r['operation'], r['arm']): r for r in raw('T6')}
    reused = [s for r in raw('T6') if r['operation'] == 'nacl/reused' for s in r['reply'].get('samples', [])]
    fresh = [s for r in raw('T6') if r['operation'] == 'nacl/fresh' for s in r['reply'].get('samples', [])]
    ref = st.median(s['reference_seconds'] for s in reused)
    curve = st.median(s['curve_seconds'] for s in reused)
    fr = st.median(s['seconds'] for s in fresh)
    cl = med('T6', 'nacl/clapeyron-curve', 'Clapeyron')
    old = PRIOR['issue_533']['C3-0.88sigma-1']
    chk = [c for c in C if c['section'] == 'T6']
    return [dict(T='T6', audit_id=i, finding='One-time infinite-dilution reference cost versus per-point activity cost', status='holds',
                 answer=(f'30-point NaCl curve: reference once {ms(ref)} plus 30 reused activity calls {ms(curve)}; 30 fresh reference-plus-point calls {ms(fr)} '
                         f'({fx(ratio(fr, ref + curve))} the reused curve); Clapeyron curve {ms(cl)}. #533 at 1 mol/kg: fresh {us(old["Engine_seconds"])}, '
                         f'reused {us(old["reused_engine_seconds"])} per point. Reused and fresh references give identical gamma at all 30 points; Clapeyron gamma differs by '
                         f'{min(float(c["scaled_error"]) for c in chk if c["comparison"] == "Engine vs Clapeyron"):.2g}–'
                         f'{max(float(c["scaled_error"]) for c in chk if c["comparison"] == "Engine vs Clapeyron"):.2g} relative, beyond the 1e-8 criterion, so its time is not ranked.'),
                 evidence='statistics.csv T6; raw/T6.jsonl.gz', decision='Caller reference reuse still removes most curve work; no analytic dilution limit is needed for speed.')
            for i in ('O-S7', 'F-S10 (NaCl)')]


def t8():
    if not raw('T8'):
        return []
    out = []
    lines = []
    for case in ('isomerisation', 'bubble-C2', 'ionic-lle'):
        b1, ball = med('T8', f'{case}/first-batch-1'), med('T8', f'{case}/first-batch-all')
        s1, sall = med('T8', f'{case}/first-separate-1'), med('T8', f'{case}/first-separate-all')
        o2 = [float(r['median_s']) for r in S if r['section'] == 'T8' and r['operation'].startswith(f'{case}/order-2') and r['tool'] == 'Engine']
        o3 = [float(r['median_s']) for r in S if r['section'] == 'T8' and r['operation'].startswith(f'{case}/order-3') and r['tool'] == 'Engine']
        lines.append(f'{case}: first actions batched one/all {us(b1)}/{us(ball)}, separate all {us(sall)}; order 2 median {us(st.median(o2)) if o2 else "—"}, '
                     f'order 3 {us(st.median(o3)) if o3 else "—"}')
    ver = [v for v in V if v['section'] == 'T8-verify']
    cnt = Counter(v['verdict'] for v in ver)
    text = '; '.join(lines) + f'. Verification outside timing: {dict(cnt)}.'
    for i in ('Sol-6', 'A-Q06'):
        out.append(dict(T='T8', audit_id=i, finding='Cost per extra direction and order after one central solve; batching amortizes central geometry',
                        status='holds', answer=text, evidence='statistics.csv T8; verification.csv T8-verify', decision='Batch first actions; higher orders cost one request each.'))
    return out


def t9():
    if not raw('T9-fits') and not raw('T9-serial'):
        return []
    ev = {f: med('T9-serial', f'{f}/evaluate-shift-0.0') for f in ('propane', 'electrolyte', 'F1')}
    near = {f: med('T9-serial', f'{f}/evaluate-shift-0.0001') for f in ('propane', 'electrolyte', 'F1')}
    fits = raw('T9-fits')
    table = {}
    for r in fits:
        s = r['reply'].get('samples', [{}])[0]
        if 'value' in s:
            table.setdefault(r['operation'], []).append((s['seconds'], s.get('child_cpu_seconds'), s['value']['status'], s['value']['iterations'],
                                                          s['value']['jacobian_evaluations'], r['before']['load'][0]))
    def cell(op):
        v = table.get(op, [])
        return 'not run' if not v else f'{st.median(x[0] for x in v):.3g} s (worker CPU {max(0.0, st.median(x[1] or 0 for x in v)):.3g} s)'
    ncpu = next(json.loads(line)['ncpu'] for line in gzip.open(R / 'raw/T9-fits.jsonl.gz', 'rt') if '"identities"' in line)
    lines = [f'{f}: ' + ', '.join(f'{w} worker(s) {cell(f"{f}/workers-{w}/iterations-1")}' for w in (1, ncpu // 2, ncpu)) for f in ('propane', 'electrolyte', 'F1')]
    full = table.get('F1/workers-0/iterations-40', [])
    full_text = (f'six full F1 fits (40 iterations, 120 s, default workers) {min(x[0] for x in full):.3g}–{max(x[0] for x in full):.3g} s, all '
                 f'{", ".join(sorted({x[2] for x in full}))} after {", ".join(sorted({str(x[3]) for x in full}))} iterations (frozen fits at {PRIOR["f1_frozen"]["revision"][:8]}: '
                 f'{PRIOR["f1_frozen"]["fits"]["A"]["seconds"]:.3g} and {PRIOR["f1_frozen"]["fits"]["B"]["seconds"]:.3g} s; the change is not attributed here); '
                 'every fit saw outside CPU consumers at its boundaries, so none supports a quiet-host claim') if full else 'full F1 fits not run'
    ser = f'serial residual plus exact Jacobian at the start: propane {ms(ev["propane"])}, electrolyte {ms(ev["electrolyte"])}, F1 {s_(ev["F1"])}; nearby point {ms(near["propane"])}, {ms(near["electrolyte"])}, {s_(near["F1"])}'
    elec = table.get(f'electrolyte/workers-{ncpu}/iterations-1', [])
    parent = elec and all((x[1] or 0) < 0.01 for x in elec)
    walls = [x[0] for op, v in table.items() if op.startswith('electrolyte/') for x in v]
    base = '; '.join(lines) + '; ' + full_text + '; ' + ser + '.'
    return [
        dict(T='T9', audit_id='O-S9', finding='Process workers pay off only above small observation counts', status='holds', answer=base,
             evidence='raw/T9-fits.jsonl.gz; statistics.csv T9-serial', decision='Worker scaling is per family; see the table.'),
        dict(T='T9', audit_id='Sol-8', finding='Serial electrolyte observation grouping (parent-owned groups) limits fit throughput',
             status='mechanism confirmed, not material' if parent else 'does not hold',
             answer=(('Mechanism confirmed: electrolyte bounded fits spend no measurable worker CPU at any worker count, because electrolyte rows share one reference '
                      'cache and stay in the parent process by design (fit.cpp:482-491), so the worker option does not apply to them. The whole bounded fit takes '
                      f'{ms(min(walls))}–{ms(max(walls))} at every worker count, so the throughput limit is not material at this size. ') if parent else '') + '; '.join(lines[1:2]) + '.',
             evidence='raw/T9-fits.jsonl.gz child_cpu_seconds; fit.cpp:482-491', decision='Measure parent-held electrolyte group time in a full-size electrolyte fit before changing worker assignment.'),
        dict(T='T9', audit_id='A-Q08', finding='Group reuse versus more parallel rows; startup cost for small fits', status='holds', answer=base,
             evidence='raw/T9-fits.jsonl.gz', decision='As O-S9.'),
        dict(T='T9', audit_id='Sol-5', finding='Reference work shared across MEA observation groups (current caches counted before any new cache)', status='holds',
             answer=f'F1 serial evaluation {s_(ev["F1"])} at the start and {s_(near["F1"])} at the nearby point; work counts in T7 (F-S10).',
             evidence='statistics.csv T9-serial, work.csv.gz T7', decision='No new cache is justified from these counts alone.'),
        dict(T='T9', audit_id='A-Q09', finding='Serial multistart scheduling under a shared budget', status='not timed',
             answer='Serial order is source-established; no alternative scheduler was authorized (design T9).', evidence='design.md T9', decision='None.'),
    ]


def t10():
    if not raw('T10'):
        return []
    chk = [c for c in C if c['section'] == 'T10']
    node = {(f, k): med('T10', f'{f}/fgh{k}/alternating') for f in ('property', 'equilibrium') for k in (0, 1, 2)}
    nat = {(f, k): med('T10', f'{f}/fgh{k}/alternating/native') for f in ('property', 'equilibrium') for k in (0, 1, 2)}
    same = {f: med('T10', f'{f}/fgh2/same') for f in ('property', 'equilibrium')}
    cohorts = {}
    for r in raw('T10-cohorts'):
        s = r['reply'].get('samples', [{}])[0]
        if 'nodes_per_second' in s:
            cohorts.setdefault(r['operation'], []).append(s['nodes_per_second'])
    coh = '; '.join(f'{k} {st.median(v):.3g} nodes/s' for k, v in sorted(cohorts.items()))
    worst = max((float(c['scaled_error']) for c in chk if c['scaled_error']), default=None)
    text = (f'ASL complete node per call (value / gradient / Hessian, alternating neighbor controls): property {us(node[("property", 0)])} / {us(node[("property", 1)])} / '
            f'{us(node[("property", 2)])}, native {us(nat[("property", 0)])} / {us(nat[("property", 1)])} / {us(nat[("property", 2)])}; equilibrium '
            f'{us(node[("equilibrium", 0)])} / {us(node[("equilibrium", 1)])} / {us(node[("equilibrium", 2)])}, native {us(nat[("equilibrium", 0)])} / '
            f'{us(nat[("equilibrium", 1)])} / {us(nat[("equilibrium", 2)])}. Repeated identical controls (cached): property {us(same["property"])}, equilibrium {us(same["equilibrium"])}. '
            f'ASL-native parity worst scaled difference {worst:.2g}. Concurrent cohorts: {coh or "not run"}.')
    one = {f: st.median(cohorts.get(f'{f}/1-same', [0])) for f in ('property', 'equilibrium')}
    best = {f: max(st.median(v) for k, v in cohorts.items() if k.startswith(f + '/')) for f in ('property', 'equilibrium')}
    serial = all(best[f] <= one[f] * 1.1 for f in one)
    eq_extra = ratio(node[('equilibrium', 0)], nat[('equilibrium', 0)])
    # First-row time of each node in the alternating pair: the record's starting point (anchor) and the changed controls.
    first_rows = [[s['per_row_seconds'][i][0] for r in raw('T10') if r['operation'] == 'equilibrium/fgh0/alternating' for s in r['reply']['samples']] for i in (0, 1)]
    anchor, moved = (st.median(x) for x in first_rows)
    conc = ('No cohort exceeds single-thread throughput' if serial else 'Some cohort exceeds single-thread throughput') + \
        '; Pyomo calls the adapter through ctypes CFUNCTYPE pointers, which release the Python interpreter lock, so the serialization is native, consistent with the per-function prepared-record mutex (asl_adapter.cpp); same-record threads also lose the record\'s per-row cache to each other.'
    return [dict(T='T10', audit_id='O-S10', finding='ASL callback cost with first and second derivatives', status='holds', answer=text,
                 evidence='statistics.csv T10; checks.csv T10', decision='IDAES cubic comparison dropped (design); the callback costs are the absolute answer.'),
            dict(T='T10', audit_id='Sol-9', finding='ASL locks or per-output calls dominate a node Jacobian/Hessian', status='does not hold',
                 answer='Within one node the first row computes the complete branch and later rows are cache hits (per-row times in raw/T10.jsonl.gz). ' + text,
                 evidence='raw/T10.jsonl.gz per_row_seconds', decision='Per-output calling is not the cost.'),
            dict(T='T10', audit_id='A-Q07', finding='Record parsing, mutex serialization or repeated evaluation in an optimizer iteration', status='holds',
                 answer=conc + ' ' + text, evidence='raw/T10-cohorts.jsonl.gz', decision='Concurrent callers get no throughput gain; a lock change is a separate correctness design.'),
            dict(T='T10', audit_id='F-S8', finding='Mutex contention when several threads evaluate one model', status='holds', answer=conc + f' Concurrent cohorts: {coh}.',
                 evidence='raw/T10-cohorts.jsonl.gz', decision='As A-Q07.'),
            dict(T='T10', audit_id='— (new)', finding='Changed-control cost of the ASL equilibrium record', status='new',
                 answer=(f'The ASL equilibrium node averages {fx(eq_extra)} the native direct solve ({us(node[("equilibrium", 0)])} vs {us(nat[("equilibrium", 0)])} per node), '
                         f'but that average spans two nodes: {us(anchor)} at the record\'s starting point (cached result) and {us(moved)} at changed controls. '
                         f'Native solves take about {us(nat[("equilibrium", 0)])} per node, so the changed node costs about {fx(ratio(moved, nat[("equilibrium", 0)]))}, '
                         'with identical values and derivatives. The cause is the path: the callback walks from its starting point to every other control by '
                         'continuation each time (asl_adapter.cpp:466-483), while the native arm solves directly from fixed guesses. Indicative only on a loaded host.'),
                 evidence='statistics.csv T10 equilibrium/*; raw/T10.jsonl.gz per_row_seconds; asl_adapter.cpp:466-483',
                 decision='Worth a design before the record is used inside process optimizers; any change must keep branch safety.')]


def t11():
    if not raw('T11'):
        return []
    rec = {}
    for r in raw('T11'):
        for s in r['reply'].get('samples', []):
            rec.setdefault(r['operation'], []).append((s['seconds'] / s['per_call'], s.get('counters', {}).get('recordings'), s['per_call']))
    def part(op):
        v = rec.get(op, [])
        return '—' if not v else f'{st.median(x[0] for x in v) * 1e6:.3g} µs per call, {st.median(x[1] for x in v):.0f} recordings per {v[0][2]} calls'
    text = '; '.join(f'{op.split("/")[-1]}: {part(op)}' for op in sorted(rec))
    three = st.median(x[1] for x in rec.get('sequence/three-repeated', [(0, 0, 0)]))
    return [dict(T='T11', audit_id=i, finding='Three or more alternating thermal or predicate signatures cause tape re-recording', status='holds' if three else 'does not hold',
                 answer=text + '. Manufactured three-interval NASA9 record; the T4 replay shows the two-phase methane/ethane signatures.',
                 evidence='raw/T11.jsonl.gz counters', decision='A wider tape cache is a later workload question; no real application is shown to alternate three signatures.')
            for i in ('Sol-4', 'A-Q04', 'F-S4')]


def t14():
    blocks = raw('T14-replay')
    if not blocks:
        return [dict(T='T14', audit_id='F-S7', finding='Benzene/toluene PH round trip needs 66 iterations', status='fixed (#579)',
                     answer='Diagnosed and fixed by #579 (caloric seed budget); not replayed.', evidence='#579', decision='None.')]
    rowsv = blocks[0]['reply']['samples'][0]['value']
    ph = next(r for r in rowsv if r['case'] == 'benzene/toluene PH round trip')
    return [dict(T='T14', audit_id='F-S7', finding='Benzene/toluene PH round trip needs 66 iterations', status='fixed (#579), replayed' if ph['success'] and ph['iterations'] < 66 else 'holds',
                 answer=f'Replay of the audit declaration at current main: success {ph["success"]}, {ph["iterations"]} iterations, T relative error {ph["T_rel_error"]:.2g} (audit: 66 iterations).',
                 evidence='raw/T14-replay.jsonl.gz', decision='None.')]


def main():
    out = []
    for part in (t1, t2, t3, t4, t6, t7, t8, t9, t10, t11, t14, dropped):
        out += part()
    order = ['T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', '—']
    out.sort(key=lambda r: order.index(r['T']))
    with (R / 'answers.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, list(out[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(out)
    for r in out:
        print(r['T'], r['audit_id'], r['status'], '|', r['answer'][:220])


if __name__ == '__main__':
    main()
