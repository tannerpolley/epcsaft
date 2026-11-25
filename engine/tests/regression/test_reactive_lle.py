"""Independent concentration and molal activity comparisons on a manufactured electrolyte LLE."""
import math
from collections.abc import Sequence
from itertools import combinations
from typing import Any

import pytest

from epcsaft import Mixture, Parameters, regression
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


def _problem(log_k: Sequence[float], salt: float = 0.05) -> Problem:
    return Problem(T=TEMPERATURE, P=PRESSURE,
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
    aqueous = model.state(TEMPERATURE, rho=density_aq, x=x[:6])
    reference = model.reference(TEMPERATURE, PRESSURE, solvent_fractions=(1, 0, 0, 0, 0, 0), ions=("s2", "s3"), stoichiometry=(1, 1))
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
        _, ph, _, distribution, exchange = _rows(p, _problem(BASE_LOG_K, salt), observed)
        rows.extend((ph, distribution, exchange))
    fitted = regression.fit(p, _coordinates(), rows, controls=regression.FitControls(
        maximum_iterations=40, maximum_elapsed_time_seconds=30.0, workers=2))
    assert fitted.usable and fitted.status == "converged", fitted.message
    assert fitted.covariance.rank == 2
    assert fitted.optimizer == pytest.approx(known, abs=1e-6, rel=0)
    for i, (predicted, actual) in enumerate(zip(_shifted(fitted.optimizer), _shifted(known))):
        same(evidence, f"reactive-lle/recovered-ln-k/{i}", predicted, actual, 1e-6, 0.0,
             residuals=fitted.residuals, rank=fitted.covariance.rank)
