"""ALG-DENSITY: finite crossings, slope zeros and dilute ethane limits."""
import math

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, Pure, _core
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, solve_equilibrium

from tests.native_inputs import ETHANE, parameters


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
    for budget in (650, 780):
        retained = model.density_diagnostics(
            T, P, nodes=nodes, density_tolerance=1e-9, log_tolerance=1e-14,
            max_evaluations=budget, max_seconds=10.0, adaptive=True,
            finite_upper=0.0, evaluation_order=_core.DerivativeOrder.PressureSlope,
        )
        assert len(retained.pressure_roots) == 3
        assert all(abs(root.molar_density - reference) <= 1e-9
                   for root, reference in zip(retained.pressure_roots, expected))
    exact = model.density_diagnostics(
        T, model.state(T, rho=100.0).pressure, nodes=(50.0, 100.0, 200.0, 300.0),
        density_tolerance=1e-9, log_tolerance=1e-14, max_evaluations=3, max_seconds=10.0,
        adaptive=True, finite_upper=0.0, evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    assert [root.molar_density for root in exact.pressure_roots] == [100.0]
    slope_search = model.density_diagnostics(
        308.8, 1e7, nodes=(5000.0, 5500.0, 6000.0, 7000.0, 8000.0),
        density_tolerance=1e-9, log_tolerance=1e-14, max_evaluations=5, max_seconds=10.0,
        adaptive=True, finite_upper=0.0, evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    assert len(slope_search.unresolved) == 1
    bracket = slope_search.unresolved[0]
    assert abs(bracket.lower - 5500.0) <= 1e-9 and abs(bracket.upper - 7000.0) <= 1e-9
    assert bracket.cause == "density diagnostic budget exhausted"


def test_declared_liquid_stays_on_its_metastable_root_inside_the_loop() -> None:
    """#314: an Ipopt barrier on ln eta pushed the stationary liquid start across the spinodal."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))  # identical species: pure ethane
    for P in (4.75e6, 4.6735e6):  # saturation 4.7802 MPa, liquid fold 4.673 MPa at 305 K
        result = solve_equilibrium(model, Problem(T=305.0, P=P, feed=Amounts({"s0": 1.0, "s1": 1.0}),
                                                  phases=[Phase("L", kind="liquid")]))
        assert result.success, result.message
        root = model.state(305.0, P=P, x=(0.5, 0.5), phase="liquid").molar_density
        assert math.isclose(result.molar_densities[0], root, rel_tol=1e-10, abs_tol=0.0), (P, root)


def test_declared_liquid_less_densely_packed_than_the_vapor_is_refused() -> None:
    """#599: methane/ethane at 266.5 K, above the z = 0.5 critical temperature, reaches the upper dew point with L and V swapped."""
    model = Mixture(parameters(((16.043, 1.0, 3.7039, 150.03), ETHANE), [("k_ij", (0, 1), 0.0)]))
    result = solve_equilibrium(model, Problem(T=266.5, P=Free(5.0e6), feed=Amounts({"s0": 0.5, "s1": 0.5}),
                                              phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))]))
    assert not result.success
    assert result.message.startswith("declared liquid L (packing fraction 0.152977") and "declared vapor V (0.160234" in result.message
    assert (result.iterations, result.evaluations, result.hessian_calls) == {"mumps": (24, 49, 25), "ma57": (24, 48, 25)}[result.linear_solver]
