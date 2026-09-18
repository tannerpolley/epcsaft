"""Public fitting transport (#49): the native R2-J2/D86 rows through epcsaft.regression, name mapping and typed failures."""
import math
from typing import Any

import pytest

from epcsaft import (
    DensityBranch, Mixture, ParameterActionDirection, ParameterError, ParameterFamily, PressureActionDirection, Pure, _core, regression,
)
from epcsaft.equilibrium import Amounts, Discover, Free, Phase, Pinned, Problem, solve_equilibrium

from ..native_inputs import ETHANE, ROWS, associating_binary, parameters, same

METHANE = (16.043, 1.0, 3.7039, 150.03)
BOUNDS = {"dispersion_energy_over_k": (140.0, 170.0), "k_ij": (-0.2, 0.2)}


def neutral_binary(epsilon: float = 150.03, kij: float = 0.0) -> Any:
    """properties_support::neutral_binary() with the R2-J2 truth or start values."""
    return parameters(((METHANE[0], METHANE[1], METHANE[2], epsilon), ETHANE), pairs=[("k_ij", (0, 1), kij)])


def bubble(T: float) -> Problem:
    """equilibrium_support::boundary(): L free, V pinned 0, P free from 1e5 Pa, equimolar amounts."""
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                   T=T, P=Free(1.0e5), feed=Amounts({"s0": 0.5, "s1": 0.5}))


def manufactured() -> tuple[Any, list[Any]]:
    """The R2-J2 truth (epsilon 152 K, k_ij -0.02) and its noise-free observations."""
    truth = neutral_binary(152.0, -0.02)
    model = Mixture(truth)
    observations = [regression.observation(truth, "pressure", bubble(T), observed=1.0, form="log_ratio", mass_basis=False) for T in (170.0, 180.0, 190.0)]
    observations.append(regression.observation(truth, "mole_fraction", bubble(180.0), observed=1.0, phase="V", species=("s0",),
                                               form="difference", scale=0.01, mass_basis=False))
    P_bubble = solve_equilibrium(model, bubble(180.0)).pressure[0]
    observations.append(regression.density(truth, 180.0, 2.0 * P_bubble, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=1.0,
                                           form="log_ratio", mass_basis=True))
    for item in observations:
        prediction = regression.predict(model, item)
        assert prediction.status == regression.ObservationStatus.Available and prediction.value is not None
        item.observed = prediction.value
    return truth, observations


def coordinates(start: Any) -> list[Any]:
    return [regression.coordinate(start, "dispersion_energy_over_k", "s0", origin=150.0, scale=10.0, bounds=BOUNDS["dispersion_energy_over_k"]),
            regression.coordinate(start, "k_ij", ("s0", "s1"), origin=0.0, scale=0.1, bounds=BOUNDS["k_ij"])]


def test_public_recovery_matches_the_native_rows(evidence: Any) -> None:
    """The same manufactured request through epcsaft.regression reproduces regression_checks R2-J2 (r4/*)."""
    _, observations = manufactured()
    start = neutral_binary(148.0, 0.0)
    result = regression.fit(start, coordinates(start), observations)
    assert result.status == "converged" and result.usable and result.trial_failures == 0
    same(evidence, "r4/recovery/epsilon", result.physical[0])
    same(evidence, "r4/recovery/kij", result.physical[1])
    same(evidence, "r4/recovery/initial_cost", result.initial_cost)
    assert result.final_cost <= 1.0e-20 and ROWS["r4/recovery/final_cost"] <= 1.0e-20
    for i, value in enumerate(result.predictions):
        same(evidence, f"r4/recovery/prediction/{i}", value)
    assert result.covariance.status == "available" and result.covariance.rank == 2
    for k in range(4):
        same(evidence, f"r4/recovery/covariance/{k}", result.covariance.optimizer[k])
        same(evidence, f"r4/recovery/physical_covariance/{k}", result.covariance.physical[k])
    same(evidence, "r4/recovery/variance_factor", result.covariance.variance_factor)
    for k, value in enumerate(result.optimizer_jacobian):
        same(evidence, f"r4/recovery/optimizer_jacobian/{k}", value)
    assert math.isclose(result.physical[0], 152.0, rel_tol=1.0e-6) and abs(result.physical[1] + 0.02) <= 1.0e-6
    # The fitted document, renamed for the public constructors, reproduces the fitted predictions; the fit's
    # final evaluation was warm-started, a cold re-solve lands within the 1e-10 row criterion (observed 1e-13).
    fitted = Mixture(regression.fitted_parameters(start, result))
    for item, value in zip(observations, result.predictions):
        assert math.isclose(regression.predict(fitted, item).value or 0.0, value, rel_tol=1.0e-9)


def test_public_density_actions_match_the_native_rows(evidence: Any) -> None:
    """Mixture/Pure.pressure_actions with no observables is the #86 density-only request (d86/*)."""
    active = [(ParameterFamily.DispersionEnergyOverK, "epsilon"), (ParameterFamily.SegmentCount, "m"), (ParameterFamily.SigmaConstant, "sigma")]
    model = Pure(parameters((ETHANE,)), active_parameters=[_core.ActiveParameter(family, (0,)) for family, _ in active])
    directions = [PressureActionDirection(parameter=[ParameterActionDirection(family, (0,), 1.0)]) for family, _ in active]
    batch = model.pressure_actions(200.0, 217800.6949469248, phase="liquid", observables=(), directions=directions)
    assert batch.density_status == _core.PropertyActionStatus.Available and batch.results == []
    for k, (_, name) in enumerate(active):
        same(evidence, f"d86/parameter/{name}/fine", batch.density_action[k])
    thermal = model.pressure_actions(200.0, 217800.6949469248, phase="liquid", observables=(),
                                     directions=[PressureActionDirection(temperature=1.0), PressureActionDirection(pressure=1.0)])
    same(evidence, "d86/temperature/fine", thermal.density_action[0])
    same(evidence, "d86/pressure", thermal.density_action[1])
    absent = model.pressure_actions(200.0, 1.0e10, phase="liquid", observables=(), directions=directions)
    assert absent.density_status == _core.PropertyActionStatus.DensityUnavailable and absent.density_action == []
    with pytest.raises(ParameterError):
        model.density_actions(200.0, 1000.0, observables=(), directions=[_core.DensityActionDirection()])


def test_name_mapping_and_tie_map() -> None:
    """Names resolve to the owner's slot convention; ties map physical rows to optimizer columns."""
    p = associating_binary()
    pair = regression.coordinate(p, "k_ij", ("s1", "s0"), origin=0.0, scale=1.0)
    assert pair.slots == [[1, 2]] and pair.family == ParameterFamily.Kij
    site = regression.coordinate(p, "association_energy_over_k", (("s0", "site0"), ("s0", "site1")), origin=0.0, scale=1.0)
    assert site.slots == [[1, 4]]
    tied = regression.coordinate(p, "sigma_constant", "s0", "s1", origin=2.5, scale=0.1, bounds=(2.0, 4.0))
    assert tied.slots == [[0], [1]] and (tied.lower, tied.upper) == (2.0, 4.0)
    a_set = regression.CoordinateSet()
    a_set.coordinates = [tied, pair]
    assert a_set.tie() == [0, 0, 1] and [(item.family, item.indices) for item in a_set.active()] == [
        (ParameterFamily.SigmaConstant, [0]), (ParameterFamily.SigmaConstant, [1]), (ParameterFamily.Kij, [1, 2])]
    with pytest.raises(ParameterError):
        regression.density(p, 300.0, 1.0e5, {"s0": 1.0}, phase="liquid", observed=1.0, form="difference", mass_basis=False)
    for bad in (lambda: regression.coordinate(p, "k_ij", ("s0", "s0"), origin=0.0, scale=1.0),
                lambda: regression.coordinate(p, "sigma_constant", "s9", origin=0.0, scale=1.0),
                lambda: regression.coordinate(p, "dielectric", "s0", origin=0.0, scale=1.0),
                lambda: regression.observation(p, "pressure", Problem(phases=Discover(), T=300.0, P=1.0e5, feed=Amounts({"s0": 1.0})), observed=1.0, form="log_ratio", mass_basis=False),
                lambda: regression.observation(p, "mole_fraction", bubble(300.0), observed=1.0, phase="X", form="log_ratio", mass_basis=False),
                lambda: regression.density(p, 300.0, 1.0e5, {"s9": 1.0}, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False)):
        with pytest.raises(ParameterError):
            bad()


def test_predictions_and_failures_keep_their_meaning() -> None:
    """Values come from the state owner; a failing observation is typed, and an initial failure names its row."""
    p = neutral_binary()
    model = Mixture(p)
    item = regression.density(p, 180.0, 5.0e6, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=15000.0, form="log_ratio", mass_basis=False)
    prediction = regression.predict(model, item)
    state = model.state(180.0, P=5.0e6, x=(0.5, 0.5), phase="liquid")
    assert prediction.value == state.molar_density and prediction.unit == "mol/m^3"
    assert prediction.residual == math.log(state.molar_density / 15000.0)
    mass = regression.density(p, 180.0, 5.0e6, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=15000.0, form="log_ratio", mass_basis=True)
    mass_prediction = regression.predict(model, mass)
    assert mass_prediction.unit == "kg/m^3" and mass_prediction.value is not None
    assert math.isclose(mass_prediction.value, state.molar_density * (0.5 * 0.016043 + 0.5 * 0.03007), rel_tol=1.0e-12)
    none = regression.predict(model, regression.density(p, 200.0, 1.0e10, {"s0": 1.0}, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False))
    assert none.status == regression.ObservationStatus.DensityUnavailable and none.value is None and none.message
    assert prediction.branch == state.density_branch and prediction.stable_root_count == state.stable_root_count
    assert state.density_branch == DensityBranch.Unspecified and state.stable_root_count == 1  # one stable root at 5 MPa
    _, observations = manufactured()
    observations.append(regression.density(p, 200.0, 1.0e10, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False))
    result = regression.fit(p, coordinates(p), observations)
    assert result.status == "initial_evaluation_failed" and not result.usable and result.predictions == []
    assert [(item.observation, item.status) for item in result.failures] == [(5, regression.ObservationStatus.DensityUnavailable)]
    single = regression.fit(p, coordinates(p), observations[:1])
    assert single.usable and single.covariance.status == "rank_deficient" and single.covariance.optimizer == []
