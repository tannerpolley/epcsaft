"""ALG-DENSITY: finite crossings, slope zeros and dilute ethane limits."""
import math

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Pure, _core

from ..native_inputs import ETHANE, parameters


def test_density_crossings_and_finite_domain() -> None:
    model = Pure(parameters((ETHANE,)))
    T, P = 200.0, 217800.69494692178
    cap = model.molar_density_upper_bound(T)
    nodes = [math.exp(math.log(1e-3) + k / 512 * math.log(cap * .999 / 1e-3)) for k in range(513)]
    diagnostic = model.density_diagnostics(
        T, P, nodes=nodes, density_tolerance=1e-9, log_tolerance=1e-14,
        max_evaluations=100000, max_seconds=60.0, adaptive=True,
        finite_anchor=1e-3, finite_upper=2 * cap, evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    expected = (138.1005458522085, 3218.7215247759887, 17541.25129083695)
    assert len(diagnostic.pressure_roots) == 3
    for root, reference in zip(diagnostic.pressure_roots, expected):
        assert abs(root.molar_density - reference) <= 1e-9
        assert root.stop == 'exact_numerical_zero' or math.log(root.upper / root.lower) <= 1e-14
    assert [root.pressure_density_derivative > 0 for root in diagnostic.pressure_roots] == [True, False, True]
    assert len(diagnostic.spinodals) == 2
    assert max(abs(fold.residual) for fold in diagnostic.spinodals) <= 1e-10
    boundary = diagnostic.finite_domain
    assert boundary is not None
    assert boundary.lower <= cap <= boundary.upper
    assert (boundary.upper - boundary.lower) / cap <= 1e-13
    errors = []
    for rho in (1e-6, 5e-7):
        state = model.state(T, rho=rho)
        errors.append((abs(state.pressure / (rho * R * T) - 1), abs(state.pressure_density_derivative / (R * T) - 1)))
    assert max(errors[0]) <= 1e-6
    assert all(fine < coarse for coarse, fine in zip(errors[0], errors[1]))
    partial = model.density_diagnostics(
        T, P, nodes=(100.0, 200.0), density_tolerance=1e-9, log_tolerance=1e-14,
        max_evaluations=4, max_seconds=60.0, adaptive=True,
        finite_upper=0.0, evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    assert len(partial.unresolved) == 1
    assert partial.unresolved[0].lower <= expected[0] <= partial.unresolved[0].upper
    assert partial.evaluations == 4
