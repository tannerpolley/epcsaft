"""Independent concentration and molal activity comparisons on a manufactured electrolyte LLE."""
import math
from collections.abc import Sequence
from itertools import combinations
from typing import Any

import pytest

from epcsaft import Mixture, Parameters, _core, regression
from epcsaft import equilibrium as eq
from epcsaft.equilibrium import Amounts, Phase, Problem, Reaction, solve_equilibrium
from tests.native_inputs import parameters, same

TEMPERATURE, PRESSURE = 298.15, 5e6
BASE_LOG_K = (-7.0, -7.15)


def _parameters() -> Parameters:
    aqueous = (18.0153, 1.2047, 2.7927, 353.95)
    organic = (86.18, 3.0576, 3.7983, 236.77)
    cation, anion = (1.008, 1.0, 2.82, 230.0), (35.45, 1.0, 2.76, 170.0)
    species = (aqueous, organic, cation, anion, organic, organic)
    return parameters(species, [("k_ij", pair, 0.0) for pair in combinations(range(len(species)), 2)],
                      charges=(0, 0, 1, -1, 0, 0),
                      extra=[("relative_permittivity", (i,), 78.09 if i == 0 else 4.0, "dimensionless")
                             for i in (0, 1, 4, 5)]
                      + [(family, (i,), diameter, "angstrom")
                         for family in ("packing_diameter", "debye_huckel_diameter")
                         for i, diameter in ((2, 2.48), (3, 2.43))],
                      choices=("fully-dissociated-debye-huckel", "exclude-same-sign-ion-pairs", "solvent-only"))


def _problem(log_k: Sequence[float], salt: float = 0.05, *, temperature: float = TEMPERATURE) -> Problem:
    return Problem(T=temperature, P=PRESSURE,
                   phases=[Phase("aqueous", kind="liquid", support=("s0", "s2", "s3"),
                                 composition_guess=(0.94, 0.03, 0.03)),
                           Phase("organic", kind="liquid", support=("s1", "s4", "s5"),
                                 composition_guess=(0.96, 0.02, 0.02))],
                   feed=Amounts({"s0": 1.0, "s1": 1.0, "s2": salt, "s3": salt}),
                   reactions=[Reaction({"s2": -1, "s3": -1, product: 1}, name=name,
                                       equilibrium_constant=math.exp(value))
                              for product, name, value in zip(("s4", "s5"), ("extract_a", "extract_b"), log_k)])


def _reference(model: Mixture, problem: Problem) -> tuple[float, float, float, float, float]:
    solved = solve_equilibrium(model, problem)
    assert solved.success, solved.message
    density_aq, density_org = solved.molar_densities
    x = solved.mole_fractions
    temperature = solved.temperature
    assert temperature is not None
    aqueous = model.state(temperature, rho=density_aq, x=x[:6])
    reference = model.reference(temperature, PRESSURE, solvent_fractions=(1, 0, 0, 0, 0, 0), ions=("s2", "s3"), stoichiometry=(1, 1))
    fugacity = aqueous.log_fugacity_coefficient
    assert fugacity is not None and fugacity[2] is not None
    log_activity = math.log(x[2] / 0.0180153) + fugacity[2] - reference.log_fugacity[0]
    c_aq, c_a, c_b = x[2] * density_aq, x[10] * density_org, x[11] * density_org
    return log_activity, -log_activity / math.log(10.0), c_a, c_a / c_aq, (c_b / c_aq) / (c_a / c_aq)**2


def _rows(p: Parameters, problem: Problem, values: Sequence[float]) -> tuple[regression.Observation, ...]:
    ln_activity, ph, concentration, distribution, exchange = values
    observation = regression.observation
    raw = observation(p, "solute_activity", problem, phase="aqueous", species=("s2",), solvent="s0",
                      observed=ln_activity, form="difference", scale=1.0)
    acidity = observation(p, "solute_activity", problem, phase="aqueous", species=("s2",), solvent="s0",
                          negative_log10=True, observed=ph, form="difference", scale=1.0)
    c_a = observation(p, "concentration", problem, phase="organic", species=("s4",),
                      observed=concentration, form="difference", scale=1.0)
    c_b = observation(p, "concentration", problem, phase="organic", species=("s5",),
                      observed=1.0, form="difference", scale=1.0)
    c_aq = observation(p, "concentration", problem, phase="aqueous", species=("s2",),
                       observed=1.0, form="difference", scale=1.0)
    d = regression.log_linear([(1, c_a), (-1, c_aq)], observed=distribution, form="log_ratio")
    ratio = regression.log_linear([(1, c_b), (-2, c_a), (1, c_aq)],
                                  observed=math.log(exchange), form="difference", scale=1.0)
    return raw, acidity, c_a, d, ratio


def _coordinates() -> list[regression.Coordinate]:
    return [regression.reaction_coordinate("extract_a", "a", origin=BASE_LOG_K[0], scale=0.5),
            regression.reaction_log_shift({"extract_a": 1.0, "extract_b": -0.5},
                                         origin=0.0, scale=0.5, bounds=(-1.0, 1.0))]


def _shifted(point: Sequence[float]) -> tuple[float, float]:
    return BASE_LOG_K[0] + 0.5 * (point[0] + point[1]), BASE_LOG_K[1] - 0.25 * point[1]


def _countercurrent(p: Parameters, point: float, *, salt: float = 0.05, observed: float = 1.0,
                    tolerance: float = 1e-11, count: int = 2, max_sweeps: int = 80) -> regression.Observation:
    solvent = {"s0": 1.0, "s1": 1.0}
    return regression.countercurrent_concentration(
        p, [_problem(_shifted((0.0, point)), salt) for _ in range(count)],
        fixed_feeds=[{**solvent, "s2": salt, "s3": salt}, *([solvent] * (count - 1))],
        aqueous_transfer={"s2": {"s2": 1.0, "s3": 1.0}},
        organic_transfer={"s4": {"s4": 1.0}, "s5": {"s5": 1.0}},
        species="s2", volume_m3=0.001, link_tolerance=tolerance, max_sweeps=max_sweeps,
        observed=observed, form="difference", scale=observed,
    )


def _countercurrent_value(p: Parameters, point: float, tolerance: float = 1e-11) -> float:
    prediction = regression.predict(Mixture(p), _countercurrent(p, point, tolerance=tolerance))
    assert prediction.status == regression.ObservationStatus.Available, prediction.message
    assert prediction.value is not None and prediction.convergence is not None
    assert prediction.convergence.residual <= tolerance
    entering, terminal = prediction.countercurrent_stages
    assert entering is not None and entering.success and terminal is not None and terminal.success
    raffinate = terminal.amounts[0] * terminal.mole_fractions[2]
    extract = entering.amounts[1] * sum(entering.mole_fractions[10:12])
    assert abs(raffinate + extract - 0.05) <= 1e-9
    for stage in (entering, terminal):
        assert abs(stage.amounts[0] * (stage.mole_fractions[2] - stage.mole_fractions[3])) <= 1e-9
    assert prediction.value == pytest.approx(raffinate / 0.001)
    return prediction.value


def test_countercurrent_concentration_one_stage(evidence: Any) -> None:
    p = _parameters()
    model = Mixture(p)
    point = 0.0
    predicted = regression.predict(model, _countercurrent(p, point, count=1))
    direct = solve_equilibrium(model, _problem(BASE_LOG_K))
    assert direct.success and predicted.status == regression.ObservationStatus.Available
    expected = direct.amounts[0] * direct.mole_fractions[2] / 0.001
    same(evidence, "countercurrent/one-stage", predicted.value, expected, 1e-9, 0.0)
    compiled = eq.compile_problem(model, _problem(BASE_LOG_K))
    observable = eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, 0, 2)
    direction = eq.SolvedStateActionDirection(0.0, 0.0, [0.0] * 6, [],
                                              reaction_coefficients=[[0.5, 0, 0, 0], [-0.25, 0, 0, 0]])
    batch, = eq.solved_state_first_actions(compiled, direct, [observable], [direction])
    assert batch.results[0].action is not None
    action, = regression.differentiate(p, [_coordinates()[1]], [_countercurrent(p, point, count=1)],
                                       optimizer=[point]).optimizer_gradient
    same(evidence, "countercurrent/one-stage-action", action, batch.results[0].action / 0.001, 1e-9, 0.0)


def test_countercurrent_concentration_linked_actions(evidence: Any) -> None:
    p = _parameters()
    coordinate = _coordinates()[1]
    row = _countercurrent(p, 0.0)
    result = regression.differentiate(p, [coordinate], [row], optimizer=[0.0])
    action, = result.optimizer_gradient
    center = regression.predict(Mixture(p), row)
    entering = center.countercurrent_stages[0]
    assert entering is not None and entering.success
    linked_salt = entering.amounts[0] * entering.mole_fractions[2]
    fixed = {"s0": 1.0, "s1": 1.0, "s2": linked_salt, "s3": linked_salt}
    terminal = regression.countercurrent_concentration(
        p, [_problem(BASE_LOG_K)], fixed_feeds=[fixed], aqueous_transfer={}, organic_transfer={},
        species="s2", volume_m3=0.001, link_tolerance=1e-11, max_sweeps=1,
        observed=1.0, form="difference", scale=1.0,
    )
    fixed_action, = regression.differentiate(p, [coordinate], [terminal], optimizer=[0.0]).optimizer_gradient
    for step in (0.001, 0.0005):
        plus, minus = [_countercurrent_value(p, q) for q in (step, -step)]
        refined = [_countercurrent_value(p, q, 2.5e-12) for q in (step, -step)]
        floor = sum(abs(value - tight) + math.ulp(value) for value, tight in zip((plus, minus), refined)) / (2 * step)
        finite = (plus - minus) / (2 * step)
        limit = 0.0002 * abs(action) + floor
        same(evidence, f"countercurrent/linked-action/{step}", action, finite, limit, 0.0,
             fixed_terminal_feed=fixed_action)
        assert abs(fixed_action - finite) > limit


def test_countercurrent_concentration_fit_and_refusal(evidence: Any) -> None:
    p, target = _parameters(), 0.12
    model = Mixture(p)
    rows = []
    for salt in (0.04, 0.06):
        prediction = regression.predict(model, _countercurrent(p, target, salt=salt))
        assert prediction.status == regression.ObservationStatus.Available and prediction.value is not None, prediction.message
        rows.append(_countercurrent(p, 0.0, salt=salt, observed=prediction.value))
    controls = regression.FitControls(maximum_iterations=20, maximum_elapsed_time_seconds=90.0, workers=2)
    fitted = regression.fit(p, [_coordinates()[1]], rows, controls=controls)
    assert fitted.usable and fitted.status == "converged", fitted.message
    same(evidence, "countercurrent/fit-shift", fitted.optimizer[0], target, 1e-6, 0.0, residuals=fitted.residuals)

    refused = regression.predict(model, _countercurrent(p, 0.0, max_sweeps=1))
    assert refused.status == regression.ObservationStatus.CentralSolveUnavailable
    assert refused.value is None and refused.refusal_convergence is not None
    assert refused.refusal_convergence.residual > 1e-11


def test_reactive_lle_log_observations_and_actions(evidence: Any) -> None:
    p = _parameters()
    model = Mixture(p)
    problem = _problem(BASE_LOG_K)
    ln_activity, ph, c_a, distribution, exchange = values = _reference(model, problem)
    rows = _rows(p, problem, values)
    expected = (ln_activity, ph, c_a, math.log(distribution), math.log(exchange))
    for name, row, value in zip(("ln_activity", "pH", "concentration", "ln_D", "ln_exchange"), rows, expected):
        prediction = regression.predict(model, row)
        assert prediction.status == regression.ObservationStatus.Available, prediction.message
        same(evidence, f"reactive-lle/{name}/prediction", prediction.value, value, 1e-10, 1e-10)
        assert prediction.residual is not None
        assert abs(prediction.residual) <= 1e-10 * max(1.0, abs(value))

    derivatives = regression.differentiate(p, _coordinates(), rows, optimizer=[0, 0])
    for column in range(2):
        step = 1e-4  # centered O(h²) truncation, with the existing E5 first-action comparison limit
        points = [[0.0, 0.0], [0.0, 0.0]]
        points[0][column], points[1][column] = step, -step
        plus, minus = [_reference(model, _problem(_shifted(point))) for point in points]
        for i, (a, b) in enumerate(zip(plus, minus)):
            finite = (math.log(a / b) if i >= 3 else a - b) / (2 * step)
            same(evidence, f"reactive-lle/action/{i}/{column}", derivatives.optimizer_gradient[2 * i + column],
                 finite, 5e-6, 5e-6)


def test_reactive_lle_fit_recovers_reaction_constants(evidence: Any) -> None:
    p = _parameters()
    model = Mixture(p)
    known = (0.4, -0.3)
    rows: list[regression.Observation] = []
    for salt in (0.04, 0.06):
        observed = _reference(model, _problem(_shifted(known), salt))
        _, ph, concentration, distribution, exchange = _rows(p, _problem(BASE_LOG_K, salt), observed)
        concentration.scale = concentration.observed
        rows.extend((ph, concentration, distribution, exchange))
    fitted = regression.fit(p, _coordinates(), rows, controls=regression.FitControls(
        maximum_iterations=40, maximum_elapsed_time_seconds=30.0, workers=2))
    assert fitted.usable and fitted.status == "converged", fitted.message
    assert fitted.covariance.rank == 2
    assert fitted.optimizer == pytest.approx(known, abs=1e-6, rel=0)
    for i, (predicted, actual) in enumerate(zip(_shifted(fitted.optimizer), _shifted(known))):
        same(evidence, f"reactive-lle/recovered-ln-k/{i}", predicted, actual, 1e-6, 0.0,
             residuals=fitted.residuals, rank=fitted.covariance.rank)


@pytest.mark.parametrize("use_starts", [False, True], ids=["final-refresh", "perturbed-starts"])
def test_log_linear_fit_refreshes_zero_weight_row(use_starts: bool, evidence: Any) -> None:
    p = _parameters()
    known = (0.4, -0.3)
    observed = _reference(Mixture(p), _problem(_shifted(known)))
    _, ph, _, distribution, exchange = _rows(p, _problem(BASE_LOG_K), observed)
    rows = [ph, distribution, exchange]
    controls = regression.FitControls(maximum_iterations=40, maximum_elapsed_time_seconds=30.0, workers=1)
    if use_starts:
        fits = regression.fit(p, _coordinates(), rows, weights=[1, 1, 0], controls=controls,
                              starts=[[0.2, 0.0], [-0.1, 0.2]], shared_iterations=80, shared_seconds=30.0)
        assert len(fits) == 2
    else:
        fits = [regression.fit(p, _coordinates(), rows, weights=[1, 1, 0], controls=controls)]
    for i, fitted in enumerate(fits):
        assert fitted.usable and fitted.status == "converged", fitted.message
        assert fitted.optimizer == pytest.approx(known, abs=1e-6, rel=0)
        assert fitted.statuses[-1] == regression.ObservationStatus.Available
        # This held-out row comes from the final solve at the fitted constants.
        same(evidence, f"log-linear/final-refresh/{use_starts}/{i}", fitted.predictions[-1], math.log(observed[-1]),
             1e-8, 0.0)


def test_log_linear_fit_refreshes_reconciled_zero_weight_row(evidence: Any) -> None:
    p = _parameters()
    model = Mixture(p)
    problem = _problem(BASE_LOG_K)
    observed = _reference(model, problem)
    _, _, _, distribution, exchange = _rows(p, problem, observed)
    # A held reaction constant and an offset distribution target move the measured temperature.
    target = regression.log_linear(distribution.log_terms, observed=math.log(observed[3]) + 0.02,
                                   form="difference", scale=1.0)
    uncertainty = regression.MeasurementUncertainty()
    uncertainty.output_sigma = [0.005, 1.0]
    uncertainty.groups = [regression.measurement_group(p, [0, 1], [
        regression.measured_input(p, "temperature", sigma=2.0, bounds=(285.0, 310.0), unit="K")])]
    ln_k = BASE_LOG_K[0]
    coordinate = regression.reaction_coordinate("extract_a", "a", origin=ln_k, scale=0.5, bounds=(ln_k, ln_k))
    fitted = regression.fit(p, [coordinate], [target, exchange], weights=[1, 0], uncertainty=uncertainty,
                            controls=regression.FitControls(maximum_iterations=40,
                                                           maximum_elapsed_time_seconds=30.0, workers=1))
    assert fitted.usable and fitted.status == "converged", fitted.message
    assert fitted.reconciliation is not None
    temperature, = fitted.reconciliation.values
    assert temperature != TEMPERATURE
    expected = math.log(_reference(model, _problem(BASE_LOG_K, temperature=temperature))[-1])
    same(evidence, "log-linear/reconciled-refresh", fitted.predictions[-1], expected, 1e-10, 0.0,
         temperature_K=temperature)


@pytest.mark.parametrize("decades", [0.0, 0.04], ids=["matched", "offset"])
def test_log_linear_minimax_uses_decade_error_and_slope(decades: float, evidence: Any) -> None:
    p = _parameters()
    model = Mixture(p)
    problem = _problem(BASE_LOG_K)
    observed = _reference(model, problem)
    distribution = _rows(p, problem, observed)[3]
    distribution.observed /= 10.0**decades
    spec = regression.MinimaxRequest()
    spec.targets, spec.group_scales = [0], []
    spec.q, spec.tolerance, spec.activity_distance, spec.margin = 1.0, 1e-10, 1e-6, 0.1
    result = _core.evaluate_minimax(p._native, _coordinates(), [0.0, 0.0], [], [distribution], [], None, spec)
    m = result.minimax
    assert m is not None and m.availability == "values_only", result.message
    same(evidence, f"log-linear/minimax-error/{decades}", m.values[0] + m.point[-1], decades, 1e-10, 0.0)
    step = 1e-4
    dq = [_reference(model, _problem(_shifted([sign * step, 0])))[3] for sign in (1, -1)]
    finite = math.log10(dq[0] / dq[1]) / (2 * step)
    same(evidence, f"log-linear/minimax-action/{decades}", m.jacobian[0], finite, 5e-6, 0.0)


@pytest.mark.parametrize("decades", [0.0, 0.04], ids=["matched", "offset"])
def test_observations_at_manufactures_log_linear_target_in_linear_units(decades: float, evidence: Any) -> None:
    p = _parameters()
    model = Mixture(p)
    problem = _problem(BASE_LOG_K)
    template = _rows(p, problem, _reference(model, problem))[3]
    point = (0.2, -0.1)
    expected_q = _reference(model, _problem(_shifted(point)))[3]
    row, = _core.observations_at(p._native, _coordinates(), point, [], [template], [decades])
    same(evidence, f"log-linear/manufactured-target/{decades}", row.observed, expected_q / 10.0**decades,
         1e-10, 1e-10)
    prediction = regression.predict(model, row)
    same(evidence, f"log-linear/manufactured-residual/{decades}", prediction.residual, decades * math.log(10.0),
         1e-10, 0.0)
