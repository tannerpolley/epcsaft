"""Capture pre-build ASL work while replaying the accepted numerical operands."""
import csv
import hashlib
import json
from pathlib import Path

import epcsaft
from pyomo.environ import ExternalFunction
import baseline as b

OUT = b.HERE/'results/asl-trace-before'
IDENTITY = b.HERE/'results/installed-baseline-identity.json'
EXPECTED = '298a360e532f6f4c5e7ef70c5be45e8b9852310d1c388d6225c504714c3dc410'


def main():
    identity = json.loads(IDENTITY.read_text())
    assert identity['wheel_sha256'] == EXPECTED
    assert b.digest(identity['wheel']) == EXPECTED
    assert b.digest(epcsaft._core.__file__) == identity['installed_core_sha256']
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT/'summary.json').exists():
        raise RuntimeError('refusing to overwrite the ASL trace supplement')
    original = ExternalFunction.evaluate_fgh
    original_load = b.load
    library = Path(epcsaft.__file__).parent/'libepcsaft_asl.so'

    def load(name, path):
        module = original_load(name, path)
        if name == 'khudaida_live':
            # Build preparation may remove direct_url's copy; the verified immutable cache owns this identity.
            module.neutral.installed_identity = lambda: (Path(identity['wheel']), EXPECTED, library)
        return module
    captures, work, trace_values = [], [], []
    labels = ('anchor', 'start-0.9', 'start-1.1')

    def evaluate(self, args, *positional, **keywords):
        case = f'asl/{labels[len(captures)//10]}/ASL/phase-{args[5]}/component-{args[6]}'
        with epcsaft.trace() as captured:
            returned = original(self, args, *positional, **keywords)
        captures.append(captured)
        events = [e for e in captured.events if e.name == 'hessian_calls']
        solves = [s for s in captured.spans if s.name == 'solve_equilibrium' and s.convergence is not None]
        work.append(dict(case=case, feed=[float(x) for x in args[:4]], fgh=keywords.get('fgh'),
                         driver_finalizations=len(events), iterations=sum(s.convergence.iterations for s in solves),
                         evaluations=sum(s.convergence.evaluations for s in solves), hessian_calls=sum(sum(e.values) for e in events),
                         tangent_plane_trials=sum(s.name == 'tangent_plane_trial' for s in captured.spans), spans=len(captured.spans)))
        if not captured.spans:
            work[-1].update(driver_finalizations_visible=0, trial_spans_visible=0, driver_finalizations=None,
                            iterations=None, evaluations=None, hessian_calls=None, tangent_plane_trials=None,
                            availability='ASL shared-library trace unavailable')
        trace_values.extend({**row, 'case': case+'/'+row['case']} for row in b.trace_rows(captured))
        return returned

    failure = None
    control = None
    try:
        ExternalFunction.evaluate_fgh = evaluate
        b.load = load
        with epcsaft.trace() as control:
            b.asl()
        prior = list(csv.DictReader((b.HERE/'results/before/asl.csv').open()))
        assert prior == b.ROWS, 'accepted ASL numerical operands changed'
        assert len(work) == 30, f'expected 30 ASL calls, captured {len(work)}'
        assert all(row['tangent_plane_trials'] in (None, 0) for row in work), 'implicit tangent-plane trial detected'
        if not sum(row['driver_finalizations'] or 0 for row in work):
            raise RuntimeError('ASL shared-library calls expose no underlying driver trace to the installed epcsaft Capture')
    except Exception as error:
        failure = f'{type(error).__name__}: {error}'
        raise
    finally:
        ExternalFunction.evaluate_fgh = original
        b.load = original_load
        if trace_values:
            b.csv_write(OUT/'trace.csv', trace_values)
        public_control = None
        if control is not None:
            selected = {f'span-{i}' for i, s in enumerate(control.spans) if s.name == 'solve_equilibrium' and s.convergence is not None}
            selected.update(f'event-{i}' for i, e in enumerate(control.events) if e.name == 'hessian_calls')
            b.csv_write(OUT/'public-control-trace.csv', [r for r in b.trace_rows(control) if r['case'].split('/')[0] in selected])
            public_control = dict(spans=len(control.spans), driver_finalizations=sum(e.name == 'hessian_calls' for e in control.events),
                iterations=sum(s.convergence.iterations for s in control.spans if s.convergence and s.name == 'solve_equilibrium'),
                evaluations=sum(s.convergence.evaluations for s in control.spans if s.convergence and s.name == 'solve_equilibrium'),
                hessian_calls=sum(sum(e.values) for e in control.events if e.name == 'hessian_calls'))
        (OUT/'summary.json').write_text(json.dumps(dict(identity=identity, work=work, failure=failure,
            numerical_operand_rows=len(b.ROWS), numerical_operands_sha256=b.digest(b.HERE/'results/before/asl.csv'),
            numerical_operands_identical=prior == b.ROWS if 'prior' in locals() else None, public_control=public_control,
            limits=['native trace measurements on the original installed wheel', 'zero visible trial spans is sufficient only when the ASL driver trace is visible']), indent=2)+'\n')


if __name__ == '__main__':
    main()
