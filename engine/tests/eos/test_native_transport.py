"""B1-B3 (#56-#59): public State values, actions and diagnostics against fresh C++ rows."""
import math
from typing import Any, cast

import pytest

from epcsaft import IdealPolynomial, Mixture, ParameterError, Pure
from tests.native_inputs import ETHANE, ROWS, associating_binary, p3_forms, parameters, record, same

IDEAL = ("ideal_enthalpy", "ideal_internal_energy", "ideal_entropy",
         "ideal_isobaric_heat_capacity", "ideal_isochoric_heat_capacity")
TOTAL = tuple(name.replace("ideal", "total") for name in IDEAL)
KEYS = ("h", "u", "s", "cp", "cv")


def surface_model(**inputs: Any) -> Mixture:
    """The C++ example's model: constant cp with anchors at 300 K and 1 bar on [250, 450] K."""
    forms = (IdealPolynomial([80.0, 0.0], 300.0, 1000.0, 50.0), IdealPolynomial([35.0, 0.0], 300.0, 2000.0, 70.0))
    return Mixture(associating_binary(), thermochemistry=record(*forms, upper=450.0), **inputs)


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
        outside = model.state(400.0, rho=1000.0)  # residual work survives a missing ideal part
        bare = Pure(ethane).state(400.0, rho=1000.0)
        assert all(getattr(outside, field) is None for field in IDEAL + TOTAL)
        assert bare.total_enthalpy is None
        same(evidence, f"outside/{name}/h_res", outside.residual_enthalpy, bare.residual_enthalpy)

    state = surface_model().state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid")
    fields = (IDEAL[0], IDEAL[1], TOTAL[0], TOTAL[1])
    for index, field in enumerate(fields):
        same(evidence, f"state/energy/{index}", getattr(state, field))
    for index, field in enumerate(IDEAL[2:] + TOTAL[2:]):
        same(evidence, f"state/entropy/{index}", getattr(state, field))
    # Same-state assembly (P3's criterion): Total - Residual is the supplied Ideal contribution.
    residual = (state.residual_enthalpy, state.residual_internal_energy, state.residual_entropy,
                state.residual_isobaric_heat_capacity, state.residual_isochoric_heat_capacity)
    for key, ideal, total, part in zip(KEYS, IDEAL, TOTAL, residual):
        evidence.compare(f"decomposition/{key}", model=getattr(state, total) - cast(float, part),
                         reference=getattr(state, ideal),
                         criterion={"type": "scalar", "rule": "symmetric_close", "atol": 2e-12, "rtol": 2e-12})
        assert math.isclose(getattr(state, total) - cast(float, part), getattr(state, ideal), rel_tol=2e-12, abs_tol=2e-12)


def test_isobaric_enthalpy_slope_is_the_state_heat_capacity(evidence: Any) -> None:
    """(dh/dT) at fixed P, x and liquid branch equals that State's cp; P3's 1e-8 ladder criterion."""
    model = surface_model()
    cp = cast(float, model.state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid").total_isobaric_heat_capacity)
    slopes = []
    for step in (0.01, 0.003):
        high, low = (cast(float, model.state(353.15 + sign * step, P=1.0e5, x=(0.25, 0.75),
                                             phase="liquid").total_enthalpy) for sign in (1.0, -1.0))
        slopes.append((high - low) / (2.0 * step))
        evidence.compare(f"isobaric-slope/{step}", model=slopes[-1], reference=cp, inputs={"T": 353.15, "P": 1.0e5},
                         criterion={"type": "scalar", "rule": "symmetric_close", "atol": 1e-8, "rtol": 1e-8})
        assert math.isclose(slopes[-1], cp, rel_tol=1e-8)
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
    assert ROWS["p3-form-value/polynomial/h"] == before
