"""B1-B3 (#56-#59): public State values, actions and diagnostics against fresh C++ rows."""
import math
from typing import Any, cast

import pytest

from epcsaft import (
    ActiveParameter,
    DensityActionDirection,
    DensityBranch,
    IdealPolynomial,
    Mixture,
    ParameterActionDirection,
    ParameterError,
    ParameterFamily,
    PressureActionDirection,
    PropertyActionStatus,
    PropertyObservable,
    Pure,
    ThermochemistryStatus,
)
from tests.native_inputs import ETHANE, ROWS, associating_binary, p3_forms, parameters, record, same, surface_model

IDEAL = ("ideal_enthalpy", "ideal_internal_energy", "ideal_entropy",
         "ideal_isobaric_heat_capacity", "ideal_isochoric_heat_capacity")
TOTAL = tuple(name.replace("ideal", "total") for name in IDEAL)
KEYS = ("h", "u", "s", "cp", "cv")


def test_thermal_records_reach_the_native_state(evidence: Any) -> None:
    ethane = parameters((ETHANE,))
    for name, form in p3_forms().items():
        model = Pure(ethane, thermochemistry=record(form))
        interior = model.state(320.0, rho=1000.0)
        for key, ideal, total in zip(KEYS, IDEAL, TOTAL):
            same(evidence, f"p3-form-value/{name}/{key}", getattr(interior, ideal))
            same(evidence, f"p3-total/{name}/{key}", getattr(interior, total))
        for endpoint in (250.0, 350.0):  # owned closed endpoints keep their values
            owned = model.state(endpoint, rho=1000.0)
            for key, ideal in zip(KEYS, IDEAL):
                same(evidence, f"p3-endpoint/{name}/{endpoint:g}/{key}", getattr(owned, ideal))

    state = surface_model().state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid")
    for index, field in enumerate(IDEAL + TOTAL):
        same(evidence, f"state/calorics/{index}", getattr(state, field))
    # Same-state assembly (P3's criterion): Total - Residual is the supplied Ideal contribution.
    residual = (state.residual_enthalpy, state.residual_internal_energy, state.residual_entropy,
                state.residual_isobaric_heat_capacity, state.residual_isochoric_heat_capacity)
    for key, ideal, total, part in zip(KEYS, IDEAL, TOTAL, residual):
        same(evidence, f"decomposition/{key}", getattr(state, total) - cast(float, part), getattr(state, ideal), 2e-12, 2e-12)


def test_isobaric_enthalpy_slope_is_the_state_heat_capacity(evidence: Any) -> None:
    """(dh/dT) at fixed P, x and liquid branch equals that State's cp; P3's 1e-8 ladder criterion."""
    model = surface_model()
    cp = cast(float, model.state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid").total_isobaric_heat_capacity)
    slopes = []
    for step in (0.01, 0.003):
        high, low = (cast(float, model.state(353.15 + sign * step, P=1.0e5, x=(0.25, 0.75),
                                             phase="liquid").total_enthalpy) for sign in (1.0, -1.0))
        slopes.append((high - low) / (2.0 * step))
        same(evidence, f"isobaric-slope/{step}", slopes[-1], cp, 0.0, 1e-8, T=353.15, P=1.0e5)
    assert math.isclose(slopes[0], slopes[1], rel_tol=1e-8)
    # A lost thermal record would leave the residual slope: the two differ by the ideal cp.
    assert cp - cast(float, Mixture(associating_binary()).state(
        353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid").residual_isobaric_heat_capacity) > 40.0


def test_the_model_copies_its_thermal_record(evidence: Any) -> None:
    thermal = record(p3_forms()["polynomial"])
    model = Pure(parameters((ETHANE,)), thermochemistry=thermal)
    before = cast(float, model.state(320.0, rho=1000.0).ideal_enthalpy)
    thermal.component_intervals = record(IdealPolynomial([30.0, 0.0], 300.0, 3000.0, 50.0)).component_intervals
    same(evidence, "copy/unchanged", model.state(320.0, rho=1000.0).ideal_enthalpy, before)
    rebuilt = cast(float, Pure(parameters((ETHANE,)), thermochemistry=thermal).state(320.0, rho=1000.0).ideal_enthalpy)
    same(evidence, "copy/rebuilt-anchor-shift", rebuilt - before, 2000.0)
    with pytest.raises(ParameterError, match="ideal anchor outside"):  # native admission stays the owner
        Pure(parameters((ETHANE,)), thermochemistry=record(IdealPolynomial([30.0], 200.0, 0.0, 0.0)))


ENERGY, VOLUME = ParameterFamily.AssociationEnergyOverK, ParameterFamily.AssociationVolume
PAIR = (1, 4)  # both native slots of the reciprocal site pair (0, 1) in the 4 x 4 site matrix
P4_STATE = (313.15, 42441.09427799721, (0.25, 0.75))
OBSERVABLES = tuple(getattr(PropertyObservable, f"{group}{name}") for group in ("Ideal", "Residual", "Total") for name in (
    "Enthalpy", "InternalEnergy", "Entropy", "IsobaricHeatCapacity", "IsochoricHeatCapacity"))


def p4_model() -> Mixture:
    forms = (p3_forms()["polynomial"], IdealPolynomial([40.0, 0.0], 300.0, 2000.0, 70.0))
    return Mixture(associating_binary(), thermochemistry=record(*forms),
                   active_parameters=[ActiveParameter(ENERGY, PAIR), ActiveParameter(VOLUME, PAIR)])


def test_mixed_state_parameter_actions_match_native_in_both_charts(evidence: Any) -> None:
    """P4's three mixed T, rho or P, x and reciprocal-pair parameter directions, orders one to three."""
    model = p4_model()
    temperature, density, x = P4_STATE
    raw = ((0.2, 5.0, 15000.0, (-0.01, 0.01), 0.03, 0.001), (-0.1, -3.0, -10000.0, (0.005, -0.005), -0.02, -0.0005),
           (0.05, 2.0, 6000.0, (0.002, -0.002), 0.01, 0.0002))
    moves = [[ParameterActionDirection(ENERGY, PAIR, item[4]), ParameterActionDirection(VOLUME, PAIR, item[5])] for item in raw]
    direct = [DensityActionDirection(item[0], item[1], item[3], move) for item, move in zip(raw, moves)]
    isobaric = [PressureActionDirection(item[0], item[2], item[3], move) for item, move in zip(raw, moves)]
    pressure = model.state(temperature, rho=density, x=x).pressure
    for order in (1, 2, 3):
        batches = (("direct", model.density_actions(temperature, density, x=x, observables=OBSERVABLES, directions=direct[:order])),
                   ("pressure", model.pressure_actions(temperature, pressure, x=x, phase="liquid", observables=OBSERVABLES,
                                                       directions=isobaric[:order])))
        for chart, batch in batches:
            for index, result in enumerate(batch.results):
                assert result.status == PropertyActionStatus.Available
                same(evidence, f"{chart}/action/{index}/{order}", result.action)
        assert batches[1][1].diagnostics.stable_root_count is not None and batches[0][1].diagnostics.stable_root_count is None
    # Losing the parameter part moves the first-order residual-enthalpy action ten orders above the criterion.
    bare = model.density_actions(temperature, density, x=x, observables=OBSERVABLES[5:6],
                                 directions=[DensityActionDirection(*raw[0][:2], raw[0][3])]).results[0].action
    assert abs(cast(float, bare) - ROWS["direct/action/5/1"]) > 0.1


def test_parameter_action_matches_a_refining_parameter_perturbation(evidence: Any) -> None:
    """d h_res / d(eps_AB / k) against rebuilt models; P4's 5e-6 first-order criterion and 2x refinement."""
    temperature, density, x = P4_STATE
    unit = [DensityActionDirection(parameter=[ParameterActionDirection(ENERGY, PAIR, 1.0)], composition=(0.0, 0.0))]
    action = cast(float, p4_model().density_actions(temperature, density, x=x, observables=OBSERVABLES[5:6],
                                                    directions=unit).results[0].action)
    slopes = []
    for step in (1.0, 0.3):
        high, low = (Mixture(associating_binary(sign * step)).state(temperature, rho=density, x=x).residual_enthalpy
                     for sign in (1.0, -1.0))
        slopes.append((high - low) / (2.0 * step))
        same(evidence, f"parameter-perturbation/{step}", action, slopes[-1], 5e-6, 5e-6, step_K=step)
    assert abs(action) > 0.5
    assert abs(slopes[0] - slopes[1]) <= 1e-5 * max(1.0, abs(slopes[1]))
    model = p4_model()
    fixed = [DensityActionDirection(parameter=[ParameterActionDirection(VOLUME, (0,), 1.0)], composition=(0.0, 0.0))]
    zero = model.density_actions(temperature, density, x=x, observables=OBSERVABLES[5:6], directions=fixed).results[0]
    assert zero.status == PropertyActionStatus.Available and abs(cast(float, zero.action)) <= 1e-12  # known fixed slot
    unknown = [DensityActionDirection(parameter=[ParameterActionDirection(VOLUME, (999,), 1.0)], composition=(0.0, 0.0))]
    with pytest.raises(ParameterError):
        model.density_actions(temperature, density, x=x, observables=OBSERVABLES[5:6], directions=unknown)


def test_pure_forwards_its_parameter_direction(evidence: Any) -> None:
    """A Pure scalar-parameter action is nonzero and equals a perturbed-epsilon difference in both charts."""
    epsilon = ParameterFamily.DispersionEnergyOverK
    model = Pure(parameters((ETHANE,)), active_parameters=[ActiveParameter(epsilon, (0,))])
    move = [ParameterActionDirection(epsilon, (0,), 1.0)]
    observable = (PropertyObservable.ResidualEnthalpy,)
    direct = model.density_actions(350.0, 1000.0, observables=observable,
                                   directions=[DensityActionDirection(parameter=move)]).results[0].action
    pressure = model.state(350.0, rho=1000.0).pressure
    isobaric = model.pressure_actions(350.0, pressure, phase="vapor", observables=observable,
                                      directions=[PressureActionDirection(parameter=move)]).results[0].action
    for step in (0.1, 0.03):
        shifted = [Pure(parameters(((*ETHANE[:3], ETHANE[3] + sign * step),))) for sign in (1.0, -1.0)]
        for chart, action, states in (("density", direct, [item.state(350.0, rho=1000.0) for item in shifted]),
                                      ("pressure", isobaric, [item.state(350.0, P=pressure, phase="vapor") for item in shifted])):
            slope = (states[0].residual_enthalpy - states[1].residual_enthalpy) / (2.0 * step)
            assert abs(slope) > 1.0
            same(evidence, f"pure-parameter/{chart}/{step}", action, slope, 5e-6, 5e-6)


def test_root_and_thermal_diagnostics_accompany_their_states(evidence: Any) -> None:
    model = surface_model()
    rooted = model.state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid")
    closure = cast(Any, rooted.density_convergence)
    for index, value in enumerate((rooted.molar_density, rooted.pressure, rooted.pressure_density_derivative,
                                   closure.residual, float(cast(int, rooted.stable_root_count)))):
        same(evidence, f"state/root/{index}", value)  # rho, achieved P, dP/drho, root residual, stable roots
    assert rooted.density_branch == DensityBranch.Liquid and closure.residual <= closure.criterion.value
    direct = model.state(353.15, rho=rooted.molar_density, x=(0.25, 0.75))  # no root search, so no root evidence
    assert direct.density_branch == DensityBranch.Unspecified and direct.stable_root_count is None
    assert direct.requested_density_branch == direct.returned_density_branch == DensityBranch.Unspecified
    assert not direct.density_branch_mismatch
    assert direct.density_convergence is None and direct.thermochemistry_status == ThermochemistryStatus.Available
    same(evidence, "charts/h_res", direct.residual_enthalpy, rooted.residual_enthalpy)
    # Missing and out-of-interval records are distinct causes; the residual sibling is the same number.
    outside = model.state(460.0, rho=rooted.molar_density, x=(0.25, 0.75))
    bare = Mixture(associating_binary()).state(460.0, rho=rooted.molar_density, x=(0.25, 0.75))
    assert outside.thermochemistry_status == ThermochemistryStatus.OutsideInterval
    assert all(getattr(outside, field) is None for field in IDEAL + TOTAL)
    assert bare.thermochemistry_status == ThermochemistryStatus.MissingRecord and bare.total_enthalpy is None
    same(evidence, "outside/h_res", outside.residual_enthalpy, bare.residual_enthalpy)
    # At an owned endpoint the value survives and the thermal action is typed, beside an available residual sibling.
    endpoint = Pure(parameters((ETHANE,)), thermochemistry=record(p3_forms()["polynomial"])).density_actions(
        250.0, 1000.0, observables=(PropertyObservable.IdealEnthalpy, PropertyObservable.ResidualEnthalpy),
        directions=[DensityActionDirection(1.0)]).results
    same(evidence, "p3-endpoint/polynomial/250/h", endpoint[0].value)
    assert endpoint[0].status == PropertyActionStatus.ThermalBoundary and endpoint[0].action is None
    assert endpoint[1].status == PropertyActionStatus.Available and endpoint[1].action is not None


def test_requested_liquid_returns_the_sole_vapor_root() -> None:
    # #824, ALG-DENSITY: the retained negative-slope sample establishes the returned side.
    model = Pure(parameters((ETHANE,)))
    liquid = model.state(300.0, P=100.0, phase="liquid")
    vapor = model.state(300.0, P=100.0, phase="vapor")
    assert liquid.molar_density == vapor.molar_density == 0.04009108203959122
    assert liquid.pressure == vapor.pressure
    assert liquid.density_branch == vapor.density_branch == DensityBranch.Unspecified
    assert liquid.requested_density_branch == DensityBranch.Liquid
    assert liquid.returned_density_branch == vapor.returned_density_branch == DensityBranch.Vapor
    assert vapor.requested_density_branch == DensityBranch.Vapor
    assert liquid.density_branch_mismatch and not vapor.density_branch_mismatch
