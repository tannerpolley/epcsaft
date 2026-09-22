"""#87/#90 trace-state reliability: the analytic trace solution below the former composition box
and the intensive accuracy of trace rows after the amount-weighted optimizer stop."""
import math
from typing import Any

import pytest

from epcsaft import GAS_CONSTANT_J_PER_MOL_K, Mixture
from epcsaft.equilibrium import (
    Amounts, EOS_STANDARD_STATE_ID, Phase, Problem, Reaction, ReactionLogPolynomial,
    SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateActionStatus,
    SolvedStateObservable, SolvedStateObservableKind as Kind,
    compile_problem, solve_equilibrium, solved_state_actions,
)
from tests.native_inputs import ETHANE, parameters, same


@pytest.mark.parametrize("total", [1.0, 10.0])
@pytest.mark.parametrize("equilibrium_constant", [1e-4, 1e-8, 1e-12, 1e-16])
def test_trace_isomerisation_meets_intensive_screens(evidence: Any, equilibrium_constant: float, total: float) -> None:
    """Identical-ethane A <=> B at 350 K and 200000 Pa has x_B = K/(1+K) exactly, with the logistic
    A1-A3 references of analyses/2026-trace-action-diagnosis. At K = 1e-8 the solution lies below the
    former 1e-8 log-composition box of the driver, which rejected both feeds (classified residuals
    2.98e-10 and 1.09e-10 against the 1e-10 request despite Ipopt SUCCESS). Screens are the retained
    diagnosis screens: composition 1e-6 relative, actions 5e-6, 5e-5 and 5e-4 relative, no absolute floor."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    temperature = 350.0
    l1 = -1200 / temperature**2 + 0.4 / temperature + 0.002
    l2 = 2400 / temperature**3 - 0.4 / temperature**2
    l3 = -7200 / temperature**4 + 0.8 / temperature**3
    p = equilibrium_constant / (1 + equilibrium_constant)
    references = [p, p * (1 - p) * l1, p * (1 - p) * (l2 + (1 - 2 * p) * l1**2),
                  p * (1 - p) * (l3 + 3 * (1 - 2 * p) * l1 * l2 + (1 - 6 * p + 6 * p * p) * l1**3)]
    a = math.log(equilibrium_constant) - 1200 / temperature - 0.4 * math.log(temperature / 300) - 0.002 * temperature
    law = ReactionLogPolynomial(a, 1200, 0.4, 0.002, 300, 250, 450, EOS_STANDARD_STATE_ID)
    problem = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5,
                      feed=Amounts({"s0": total / 2, "s1": total / 2}), row_tolerance=1e-10,
                      reactions=[Reaction({"s0": -1.0, "s1": 1.0}, correlation=law)])
    central = solve_equilibrium(model, problem)
    assert central.success and central.requested_tolerance_met and central.convergence.residual <= 1e-10, central.message
    label = f"trace/{equilibrium_constant:g}/{total:g}"
    inputs = {"K": equilibrium_constant, "feed_mol": total, "T_K": temperature, "P_Pa": 2e5,
              "classified_residual": central.convergence.residual,
              "unscaled_residual_max": max(abs(value) for value in central.residuals),
              "solver_status": central.solver_status, "iterations": central.iterations}
    same(evidence, f"{label}/x_B", central.mole_fractions[1], references[0], atol=0.0, rtol=1e-6, unit="1",
         inputs=inputs, basis="Identity", method="logistic-analytic",
         claim="trace mole fraction equals K/(1+K) without a composition box")
    if equilibrium_constant <= 1e-8:
        assert central.mole_fractions[1] < 1e-8  # the analytic solution sits below the former box
    compiled = compile_problem(model, problem)
    observable = SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)
    direction = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [])
    for order, tolerance in enumerate((5e-6, 5e-5, 5e-4), 1):
        item = solved_state_actions(compiled, central, SolvedStateActionRequest([observable], [direction] * order)).results[0]
        assert item.status == SolvedStateActionStatus.Available, item.status
        same(evidence, f"{label}/A{order}", item.action, references[order], atol=0.0, rtol=tolerance, unit=f"K^-{order}",
             inputs=inputs, basis="Identity", method="logistic-analytic",
             claim=f"order-{order} temperature action of the trace mole fraction matches the logistic chain")
    # First-order optimality of the extensive objective with the returned constraint multipliers, from the
    # compiled objective gradient and row Jacobian rather than the driver's own reconstruction: a wrong
    # feed scaling of the material multiplier or a nonzero normalization multiplier leaves a residual of
    # order |grad F| (0.4 to 142 here) instead of rounding.
    n = len(central.coordinates)
    gradient = compiled.objective_gradient(list(central.coordinates))
    jacobian = compiled.jacobian(list(central.coordinates))
    constraint_rows = [k for k, item in enumerate(compiled.rows) if str(item.kind).split(".")[-1] not in ("Chemical", "Pressure")]
    assert len(constraint_rows) == len(central.constraint_multipliers)
    kkt = max(abs(gradient[i] + sum(mu * jacobian[row * n + i] for mu, row in zip(central.constraint_multipliers, constraint_rows)))
              for i in range(compiled.primal_size))
    same(evidence, f"{label}/kkt_residual", kkt, 0.0, atol=1e-12, rtol=0.0, unit="1",
         inputs={**inputs, "objective_gradient_max": max(abs(value) for value in gradient)}, basis="Identity",
         method="objective-gradient-plus-jacobian",
         claim="the returned constraint multipliers make the extensive objective stationary on the constraint rows")
    # Stationarity and curvature at the returned point, with the standard-state gauge made explicit so
    # the identity is independent of the driver's internal gauge: g0 = (0, -RT ln K) with the consistent
    # scalar K reproduces the same state. The balance multipliers act through the compiled orthonormal
    # balance basis, read from the row Jacobian's multiplier columns (-B_ji on a chemical row), so
    # ln(x_i rho) + mu_res_i/RT + g0_i/RT = sum_j lambda_j B_ji for both species from public state quantities.
    explicit = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5,
                       feed=Amounts({"s0": total / 2, "s1": total / 2}), row_tolerance=1e-10,
                       reactions=[Reaction({"s0": -1.0, "s1": 1.0}, equilibrium_constant=equilibrium_constant)],
                       standard_gibbs_j_per_mol=(0.0, -GAS_CONSTANT_J_PER_MOL_K * temperature * math.log(equilibrium_constant)))
    gauged = solve_equilibrium(model, explicit)
    assert gauged.success and math.isclose(gauged.mole_fractions[1], central.mole_fractions[1], rel_tol=1e-12)
    compiled_explicit = compile_problem(model, explicit)
    n = len(gauged.coordinates)
    jacobian = compiled_explicit.jacobian(list(gauged.coordinates))
    state = model.state(temperature, rho=gauged.molar_densities[0], x=tuple(gauged.mole_fractions))
    for species, standard in ((0, 0.0), (1, -math.log(equilibrium_constant))):
        row = next(k for k, item in enumerate(compiled_explicit.rows) if "Chemical" in str(item.kind) and item.index == species)
        basis = [-jacobian[row * n + compiled_explicit.primal_size + j] for j in range(compiled_explicit.rank)]
        potential = (math.log(gauged.mole_fractions[species] * state.molar_density)
                     + state.residual_chemical_potential_over_rt[species] + standard)
        same(evidence, f"{label}/stationarity/{species}", sum(m * b for m, b in zip(gauged.multipliers, basis)), potential,
             atol=1e-12, rtol=1e-12, unit="1", inputs={**inputs, "balance_basis": basis}, basis="Identity",
             method="public-state-chemical-potential",
             claim="the returned balance multipliers reproduce the species chemical potential over RT at the returned point")
    # Curvature certificate at the returned point: positive, reproduced by a fresh solve started there to 1e-9
    # (the pre-correction certificate with optimizer duals differed by 5.9e-8 at K = 1e-16), and on the
    # physical scale of the trace amount N x_B, the log-chart curvature of the trace coordinate. The
    # reproduction row checks consistency only; the scale row falsifies the feed scaling of the multipliers
    # at N = 10 (mu/N or mu N would give 4.8e-7 or -3.7e-6 against 1e-7), while at N = 1 the scalings coincide.
    warm = solve_equilibrium(model, explicit, start=gauged)
    assert gauged.reduced_hessian_minimum > 0.0 and warm.success
    same(evidence, f"{label}/curvature", warm.reduced_hessian_minimum, gauged.reduced_hessian_minimum, atol=0.0, rtol=1e-9,
         unit="1", inputs=inputs, basis="Relationship", method="warm-resolve",
         claim="the reduced-Hessian minimum certified at the returned point is reproduced from that point")
    same(evidence, f"{label}/curvature_scale", gauged.reduced_hessian_minimum, total * gauged.mole_fractions[1], atol=0.0,
         rtol=2.0 * equilibrium_constant + 1e-9, unit="1", inputs=inputs, basis="Relationship", method="trace-amount-scale",
         claim="the minimum reduced curvature of the isomerisation is the trace amount to first order in K")
