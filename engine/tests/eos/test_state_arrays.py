"""Independent homogeneous rows preserve scalar values and failure meaning."""
import copy
import math
from collections.abc import Mapping
from typing import Any, cast

import pytest

from epcsaft import Mixture, Parameters, Pure, State, StateStatus
from epcsaft.parameters import document_mapping
from tests.eos.test_association_diagnostics import water_parameters
from tests.eos.test_state import _basic_model
from tests.native_inputs import parameters
from tests.packets import parameters_json

FIELDS = tuple(name for name in dir(State) if not name.startswith("_"))


def plain(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, Mapping):
        return {k: plain(v) for k, v in value.items()}
    if hasattr(value, "name"):
        return value.name
    return {k: plain(getattr(value, k)) for k in dir(value) if not k.startswith("_") and not callable(getattr(value, k))}


def equal(actual: Any, expected: Any, field: str = "") -> None:
    if isinstance(expected, float):
        if not math.isfinite(expected):
            assert (math.isnan(actual) and math.isnan(expected)) or actual == expected
        else:
            scale = {"isothermal_compressibility": 1e-12, "isentropic_compressibility": 1e-12,
                     "isobaric_expansivity": 1e-6, "joule_thomson_coefficient": 1e-6}.get(field, 1.0)
            assert abs(actual - expected) <= 1e-10 * max(abs(expected), scale)
    elif isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            equal(actual[key], expected[key], field)
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected, strict=True):
            equal(a, b, field)
    else:
        assert actual == expected


def model_for(name: str) -> Any:
    if name in ("propane", "methane-ethane"):
        return _basic_model(name)
    if name == "salt":
        return Mixture(Parameters.from_json(parameters_json("figiel-yu-held-2025", 3),
                       components=("water", "sodium-cation", "chloride-anion")))
    pairs = [("k_ij", (0, 1), 0.0)]
    volume = math.sqrt(.032384 * .034868) * (math.sqrt(3.1771 * 3.0007) / ((3.1771 + 3.0007) / 2)) ** 3
    for edge, energy, size in [((0, 1), 2653.4, .032384), ((2, 3), 2500.7, .034868),
                               ((0, 3), 2577.05, volume), ((1, 2), 2577.05, volume)]:
        pairs.extend([("association_energy_over_k", edge, energy), ("association_volume", edge, size)])
    return Mixture(parameters(((46.068, 2.3827, 3.1771, 198.24), (18.015, 1.0656, 3.0007, 366.51)),
                   pairs=pairs, sites=((0, 0, 1, 1), (1, 1, 1, 1))))


def inputs(name: str, count: int) -> tuple[list[float], list[float], list[float], list[list[float]], str]:
    s = [k / (count - 1) for k in range(count)]
    if name == "salt":
        return ([298.15 + 5*t for t in s], [54000 + 1000*t for t in s], [1e5 + 1e5*t for t in s],
                [[.98-.02*t, .01+.01*t, .01+.01*t] for t in s], "liquid")
    if name == "association":
        return ([300 + 20*t for t in s], [18000 + 1000*t for t in s], [1e5 + 1e5*t for t in s],
                [[.4+.2*t, .6-.2*t] for t in s], "liquid")
    return ([280+40*t for t in s], [100+200*t for t in s],
            [1e5+(9e5 if name == "methane-ethane" else 1e5)*t for t in s],
            [[.3+.2*t, .7-.2*t] if name == "methane-ethane" else [1.] for t in s], "vapor")


@pytest.mark.parametrize("name", ["methane-ethane", "propane", "salt", "association"])
@pytest.mark.parametrize("chart", ["rho", "P"])
def test_state_arrays_preserve_scalar_values(name: str, chart: str) -> None:
    scalar, batch = model_for(name), model_for(name)
    temperatures, densities, pressures, compositions, phase = inputs(name, 100)
    values = densities if chart == "rho" else pressures
    rows = [scalar.state(t, **{chart: v}, x=x, phase=phase, derived_properties=True)
            for t, v, x in zip(temperatures, values, compositions, strict=True)]
    result = batch.states(temperatures, **{chart: values}, x=compositions, phase=phase, derived_properties=True)
    assert result.status == [StateStatus.Available] * 100
    for field in FIELDS:
        equal(plain(getattr(result, field)), plain([getattr(r, field) for r in rows]), field)
    equal(plain(batch.counters()), plain(scalar.counters()))


def test_state_arrays_isolate_failures() -> None:
    model = model_for("salt")
    result = model.states([298.15] * 3, rho=55000., x=[[.98,.01,.01], [.98,.012,.008], [.98,.01,.01]])
    assert result.status == [StateStatus.Available, StateStatus.ParameterError, StateStatus.Available]
    assert result.pressure[1] is None and "charge sum" in result.cause[1]
    equal(result.pressure[2], model_for("salt").state(298.15, rho=55000., x=[.98,.01,.01]).pressure)
    doc = cast(dict[str, Any], copy.deepcopy(document_mapping(water_parameters())))
    doc["topology"]["edges"][0]["energy_over_k"]["value"]["magnitude"] = 100000.
    water = Pure(Parameters.from_mapping(doc, components=["water"]))
    sequence = water.states([100000., 298.15, 100000.], P=1e5, phase="liquid")
    assert sequence.status == [StateStatus.Available, StateStatus.EosError, StateStatus.Available]
    assert "water/a, water/b" in sequence.cause[1]
    equal(sequence.pressure[2], Pure(water.parameters).state(100000., P=1e5, phase="liquid").pressure)
    doc["topology"]["edges"][0]["energy_over_k"]["value"]["magnitude"] = 10000.
    interrupted = Pure(Parameters.from_mapping(doc, components=["water"]))
    fresh = Pure(interrupted.parameters).state(100000., rho=1.)
    recovered = interrupted.states([100000., 100., 100000.], rho=[1., 1e-5, 1.])
    assert recovered.status == [StateStatus.Available, StateStatus.EosError, StateStatus.Available]
    for field in FIELDS:
        equal(plain(getattr(recovered, field)[2]), plain(getattr(fresh, field)), field)
    ideal_doc = cast(dict[str, Any], copy.deepcopy(document_mapping(model_for("propane").parameters)))
    next(f for f in ideal_doc["model_families"] if f["kind"] == "base")["choice"] = "ideal-gas"
    ideal_doc["components"][0]["coefficients"] = []
    ideal_doc["components"][0]["fixed"]["charge_number"]["value"]["unit"] = "dimensionless"
    ideal = Pure(Parameters.from_mapping(ideal_doc))
    overflow = ideal.states([300., 1e308, 300.], P=1e5)
    assert overflow.status == [StateStatus.Available, StateStatus.DomainError, StateStatus.Available]
    with pytest.raises(ValueError) as error:
        ideal.state(1e308, P=1e5)
    assert overflow.cause[1] == str(error.value)
    equal(overflow.pressure[2], Pure(ideal.parameters).state(300., P=1e5).pressure)
    assert ideal.states([], rho=200.).rows == []
    with pytest.raises(ValueError, match="array length"):
        ideal.states([300., 310.], rho=[200., 300., 400.])


def test_state_arrays_keep_optional_and_active_outputs() -> None:
    from epcsaft import ActiveParameter, ParameterFamily
    from tests.eos.test_permittivity_forms import DIELECTRIC, RHO, T, X, _dielectric

    neutral = Mixture(model_for("methane-ethane").parameters,
                      active_parameters=[ActiveParameter(ParameterFamily.DispersionEnergyOverK, [0])])
    dielectric = Mixture(_dielectric(), active_parameters=DIELECTRIC)
    for model, temperature, density, x in [(neutral, 300., 200., [.4,.6]), (dielectric, T, RHO, list(X))]:
        result = model.states([temperature, temperature+1.], rho=density, x=x)
        for k, t in enumerate([temperature, temperature+1.]):
            scalar = model.state(t, rho=density, x=x)
            for field in FIELDS:
                equal(plain(getattr(result, field)[k]), plain(getattr(scalar, field)), field)
        assert result.total_enthalpy == [None, None]
    thermal = model_for("methane-ethane")
    result = thermal.states([300., 1200.], rho=200., x={"methane": .4, "ethane": .6}, derived_properties=True)
    assert result.total_enthalpy[0] is not None and result.total_enthalpy[1] is None
    assert result.residual_enthalpy[1] is not None
