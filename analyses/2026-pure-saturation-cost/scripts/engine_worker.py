"""Select a guarded prepared image; reuse release declarations and whole-call timers."""
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('release_engine', ROOT / 'analyses/2026-release-computational-cost/scripts/engine_worker.py')
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)
image = Path(os.environ['EPCSAFT_COST_IMAGE'])
e.RELEASE = Path(os.environ['EPCSAFT_COST_SOURCE'])
e.w.ROOT = e.RELEASE
e.w.HERE = image.parents[1]


def prepared(req):
    from epcsaft import equilibrium as q
    cc = e.cross()
    inputs = json.loads((e.RELEASE / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    request = next(r for r in inputs['requests'] if r['id'] == req['case'])
    model, problem, _, _ = cc.prepare(inputs, request)
    compiled = q.compile_problem(model, problem)

    dedicated = os.environ.get('EPCSAFT_COST_ARM') == 'B'
    if not dedicated:
        spec = importlib.util.spec_from_file_location('_phase_pressures', os.environ['EPCSAFT_COST_PRESSURES'])
        pressures = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(pressures)

    def call():
        r = compiled.solve_saturation() if dedicated else compiled.solve()
        if not r.success or not r.requested_tolerance_met:
            raise RuntimeError(r.message)
        if dedicated:
            value = r.saturation_values()
            del value['diagnostics']
        else:
            phase_pressures = pressures.phase_pressures(compiled, r)
            d = r.linear_solve_diagnostics
            linear = None if d is None else dict(zip(
                ('method', 'reason', 'residual', 'evaluation_floor', 'inverse_norm', 'condition_estimate',
                 'forward_error_estimate', 'working_precision_bits', 'refinement_iterations'),
                (getattr(d, key) for key in ('system', 'reason', 'residual_norm', 'evaluation_floor',
                 'jacobian_inverse_norm', 'condition_estimate', 'forward_error_estimate', 'precision_bits', 'refinement_iterations'))))
            k = r.last_rank
            value = dict(success=r.success, requested_tolerance_met=r.requested_tolerance_met, message=r.message,
                solver_status=r.solver_status, P=r.pressure, T=r.temperature,
                phases=[dict(phase=name, rho=r.molar_densities[i], x=[r.mole_fractions[i]], P=phase_pressures[i])
                        for i, name in enumerate(('liquid', 'vapor'))], convergence=cc.convergence(r.convergence),
                max_residual=max(map(abs, r.residuals)), iterations=r.iterations, evaluations=r.evaluations,
                hessian_calls=r.hessian_calls, linear=linear,
                rank=None if k is None else {key: getattr(k, key) for key in ('basis', 'rank', 'dimension', 'threshold')},
                first_numerical_failure=r.first_numerical_failure, terminal_numerical_failure=None, rejected_reduction=None,
                last_full_evaluation=cc.convergence(r.convergence) | dict(cause='', rank=None, backward_error=None))
        return dict(value=value)
    return call


e.w.OPS['prepared-equilibrium'] = prepared
if __name__ == '__main__':
    os.environ.pop('CAMPAIGN_FIT_LAUNCHER', None)
    e.w.main()
