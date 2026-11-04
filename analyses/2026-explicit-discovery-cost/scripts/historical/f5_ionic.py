"""#647 falsifier 5 records for ionic systems: #568 F1 (Khudaida Figure 2 point 7 Discover) and the #223 case-6 and
case-7 declarations (baseline.declared_cases, unchanged): declared solve, stability_test with default options on it,
and Discover from the feed, with every result, round, trial and work record.
usage: <env>/bin/python f5_ionic.py OUT.json"""
import copy
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import falsify as f  # noqa: E402
from falsify import b, eq  # noqa: E402


def report_record(rep):
    return dict(message=rep.message, completed=rep.completed, plane=list(rep.plane_over_rt), work=dict(rep.work),
                trials=[f.trial_record(t) for t in rep.trials])


def cases():
    names = ('water', 'ethanol', 'isobutanol', 'sodium-cation', 'chloride-anion')
    model = b.Mixture(b.Parameters.from_json(b.KH_PACKET, components=names))
    feed = eq.Amounts(dict(zip(names, (0.646925133827633, 0.07045345803743744, 0.22564721508764596,
                                       0.028487096523641797, 0.028487096523641797))))
    yield '568-F1/khudaida-figure-2-point-7', model, eq.Problem(T=293.15, P=1e5, feed=feed, row_tolerance=1e-10, phases=[
        eq.Phase('Lw', kind='liquid', composition_guess=(0.80, 0.05, 0.05, 0.05, 0.05)),
        eq.Phase('Lo', kind='liquid', composition_guess=(0.20, 0.10, 0.6998, 0.0001, 0.0001))])
    for case, model, problem in b.declared_cases():
        if case.startswith(('case-6/', 'case-7/')):
            yield case, model, problem


rows = []
for case, model, problem in cases():
    declared = eq.solve_equilibrium(model, problem)
    row = dict(case=case, declared=f.result_record(declared))
    try:
        row['test'] = report_record(eq.stability_test(eq.compile_problem(model, problem), declared, eq.StabilityOptions()))
    except Exception as error:
        row['test'] = f'{type(error).__name__}: {error}'
    request = copy.copy(problem)
    request.phases, request.topology_declared = eq.Discover(), False
    try:
        found = eq.solve_equilibrium(model, request)
        d = found.discovery
        row['discover'] = dict(result=f.result_record(found), work=dict(d.work),
                               declarations=[dict(role=role, **f.result_record(r)) for role, r in zip(d.declaration_roles, d.declarations)],
                               rounds=[report_record(rep) for rep in d.rounds])
    except Exception as error:
        row['discover'] = f'{type(error).__name__}: {error}'
    rows.append(row)
    print(case, declared.success, row['discover']['result']['phases'] if isinstance(row['discover'], dict) else row['discover'], flush=True)
Path(sys.argv[1]).write_text(json.dumps(dict(engine=list(f.epcsaft.provenance()), rows=rows), indent=1, default=str) + '\n')
