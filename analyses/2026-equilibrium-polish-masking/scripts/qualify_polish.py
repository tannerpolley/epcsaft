"""Installed-wheel comparisons of physical state and fresh polish corrections (#272)."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'analyses/2026-equilibrium-polish-masking/results/resumed'
SMALL = ROOT / 'analyses/2026-equilibrium-polish-masking/results/amendment-2/small-phase-input.json'
sys.path.insert(0, str(ROOT / 'engine'))
sys.path.insert(0, str(ROOT / 'analyses/2026-tpd-phase-discovery/scripts'))
import baseline as b
from unloaded_input import declaration
from epcsaft import Mixture, Parameters, _core
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, Reaction, ReactionLogPolynomial, EOS_STANDARD_STATE_ID, solve_equilibrium
from tests.native_inputs import ETHANE, parameters
from tests.regression.test_fit import bubble

PHYSICAL = ('temperature', 'pressure', 'amounts', 'mole_fractions', 'molar_densities', 'packing_fractions')
FIELDS = PHYSICAL + ('success', 'requested_tolerance_met', 'iterations', 'evaluations', 'hessian_calls', 'residuals', 'classified_residuals', 'coordinates', 'multipliers')

def cases():
    model = Mixture(parameters((ETHANE,) * 3, [('k_ij', p, 0.) for p in ((0, 1), (0, 2), (1, 2))]))
    T = 350.
    law = ReactionLogPolynomial(math.log(3.)-1200/T-.4*math.log(T/300)-.002*T, 1200, .4, .002, 300, 250, 450, EOS_STANDARD_STATE_ID)
    yield 'masking', model, Problem(phases=[Phase('V', kind='vapor')], T=T, P=2e5, row_tolerance=1e-10,
        feed=Amounts({'s0': 1., 's1': 1e-12, 's2': 0.}), reactions=[Reaction({'s1': -1., 's2': 1.}, correlation=law)])
    for name, model, problem in b.declared_cases():
        if name in ('case-6/figure-2-point-7', 'case-7/flash-5918.5286-Pa'):
            yield name, model, problem
    model, problem = declaration(ROOT)
    yield 'mea-nominally-unloaded', model, problem
    inputs = json.loads((ROOT / 'analyses/2026-reaction-constant-fitting/results/cold-start.json').read_text())['inputs']
    model = Mixture(parameters(list(zip(inputs['molar_mass_g_mol'], inputs['m'], inputs['sigma_angstrom'], inputs['epsilon_K'])), pairs=[('k_ij', (0, 1), inputs['kij'])]))
    problem = bubble(inputs['T_guess_K'])
    problem.T, problem.P = Free(inputs['T_guess_K']), inputs['P_fixed_Pa']
    problem.reactions = (Reaction(dict(zip(model.component_ids, inputs['stoichiometry'])), name=inputs['reaction'], correlation=ReactionLogPolynomial(inputs['a'], inputs['b_K'], inputs['c'], inputs['d_per_K'], inputs['Tref_K'], inputs['Tmin_K'], inputs['Tmax_K'], inputs['standard_state_id'])),)
    yield 'free-T', model, problem
    from tests.equilibrium.test_continuation import mixture, bubble as gross_bubble
    yield 'free-P', mixture(), gross_bubble(180.)
    small = json.loads(SMALL.read_text())
    if small:
        yield 'small-positive-phase', mixture(), Problem(T=180., P=Free(small['P_Pa']), feed=Amounts(small['feed_mol']),
            phases=[Phase('L', kind='liquid', amount=Pinned(1.-small['beta']), composition_guess=list(small['xL'].values()), packing_guess=small['etaL']),
                    Phase('V', kind='vapor', composition_guess=list(small['xV'].values()), packing_guess=small['etaV'])])

def component_names(model, problem):
    result = []
    for p in problem.phases:
        support = p.support or model.component_ids
        result += [f'log x/{p.name}/{s}' for s in support]
        if not isinstance(p.amount, Pinned):
            result += [f'amount/{p.name}']
        result += [f'log eta/{p.name}']
    if isinstance(problem.T, Free): result += ['T/T']
    if isinstance(problem.P, Free): result += ['P/P']
    return result

def save_csv(path, rows):
    if rows:
        with path.open('w') as f:
            w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n"); w.writeheader(); w.writerows(rows)

def run(label, wheel):
    states, traces, inputs = [], [], {}
    def solve(name, model, problem, previous=None):
        with _core.Capture(diagnostics=True) as capture:
            result = solve_equilibrium(model, problem, previous)
        state = dict(case=name, **{k: getattr(result, k) for k in FIELDS})
        states.append(state)
        names = component_names(model, problem)
        for event in capture.events:
            if not event.name.startswith('polish'): continue
            row = dict(case=name, event=event.name, step=event.iteration)
            row.update({f'v{i}': v for i, v in enumerate(event.values)})
            if event.name == 'polish_correction':
                index = int(event.values[1]); row['component'] = names[index] if index < len(names) else f'multiplier/RT/{index-len(names)}'
            traces.append(row)
        inputs[name] = b.problem_input(problem)
        inputs[name]['components'] = list(model.component_ids)
        print(label, name, result.success, result.iterations, result.evaluations, result.hessian_calls, flush=True)
        assert result.success and result.requested_tolerance_met, (name, result.message)
        return result
    for name, model, problem in cases(): solve(name, model, problem)
    fixed = b.load('run_fixed_settler', b.KH / 'scripts/run_fixed_settler.py')
    model = Mixture(Parameters.from_json(b.KH_PACKET, components=fixed.COMPONENTS))
    previous = None
    for i, feed in sorted(fixed.feed_rows().items()):
        previous = solve(f'khudaida/point-{i}', model, fixed.engine_problem(feed), previous)
    # Homogeneous numerical rows have one CSV owner.
    fields = sorted({k for row in traces for k in row})
    with (OUT / f'{label}-trace.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n"); w.writeheader(); w.writerows(traces)
    numerical = []
    for state in states:
        for quantity in FIELDS:
            values = state[quantity] if isinstance(state[quantity], list) else [state[quantity]]
            for i, value in enumerate(values): numerical.append(dict(case=state['case'], quantity=quantity, index=i, value=value))
    state_path = OUT / f'{label}-states.csv'
    save_csv(state_path, numerical)
    retained = OUT.parent / 'amendment-2' / state_path.name
    if retained.exists() and state_path.read_bytes() == retained.read_bytes():
        state_path.unlink()
        state_path = retained
    input_path = OUT / 'inputs.json'
    declarations = json.dumps(inputs, indent=2)+'\n'
    if input_path.exists():
        assert input_path.read_text() == declarations, 'candidate and reference declarations differ'
    else:
        input_path.write_text(declarations)
    (OUT / f'{label}-policy.json').write_text(json.dumps({t.name: t.value for t in _core.tolerances()
        if t.name in ('EQUILIBRIUM_POLISH_CORRECTION', 'PHASE_DISAPPEARANCE_AMOUNT')}, indent=2)+'\n')
    (OUT / f'{label}-wheel.json').write_text(json.dumps(dict(path=str(wheel), sha256=hashlib.file_digest(wheel.open('rb'), 'sha256').hexdigest(),
        states_path=str(state_path.relative_to(ROOT)), states_sha256=hashlib.file_digest(state_path.open('rb'), 'sha256').hexdigest(),
        inputs_path=str(input_path.relative_to(ROOT)), inputs_sha256=hashlib.file_digest(input_path.open('rb'), 'sha256').hexdigest()), indent=2)+'\n')


def declarations():
    return list(cases())


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('label'); p.add_argument('wheel', type=Path); a = p.parse_args()
    run(a.label, a.wheel)
