from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from math import isinf, log
from typing import Any, cast

import pytest

from epcsaft import ActiveParameterSet, _core
from epcsaft.eos import Mixture, Pure
from epcsaft.parameters import Parameters
from epcsaft.regression._model import ParameterFit
from epcsaft.regression._result import FitResult
from epcsaft.units import unit_registry
from tests.eos._raw_cases import NEUTRAL_BINARY

_IDENTITY = "component/ethane/dispersion_energy_over_k"
_TRUE_VALUE = 191.42
_PAIR_IDENTITY = "pair/ethane/methane/k_ij"
_PAIR_TRUE_VALUE = 0.03
_TEMPERATURES = (190.0, 200.0, 210.0)


def _model(value: float = _TRUE_VALUE) -> Pure:
    mapping = deepcopy(NEUTRAL_BINARY)
    components = mapping["components"]
    assert isinstance(components, list)
    for component in components:
        if component["component_id"] != "ethane":
            continue
        for coefficient in component["coefficients"]:
            if coefficient["identity"] == _IDENTITY:
                coefficient["value"]["magnitude"] = value
    return Pure(Parameters.from_mapping(mapping, components=("ethane",)))


def _observations(model: Pure, residual: str = "log_ratio") -> list[dict[str, Any]]:
    observations = []
    for temperature in _TEMPERATURES:
        saturation = model.saturation(T=temperature * unit_registry.kelvin)
        observations.append(
            {
                "family": "pure_saturation",
                "temperature_k": temperature,
                "liquid_density_upper_mol_m3": model._native.molar_density_upper_bound(
                    temperature, [1.0]
                ),
                "multiplier": 1.0,
                "targets": [
                    {
                        "label": "saturation_pressure_pa",
                        "observed": saturation.saturation_pressure_pa,
                        "residual": residual,
                        "scale": saturation.saturation_pressure_pa,
                        "value_factor": 1.0,
                    },
                    {
                        "label": "liquid_molar_density_mol_per_m3",
                        "observed": saturation.liquid.molar_density_mol_m3,
                        "residual": residual,
                        "scale": saturation.liquid.molar_density_mol_m3,
                        "value_factor": 1.0,
                    },
                ],
            }
        )
    return observations


def _fit(
    model: Pure,
    observations: list[dict[str, Any]],
    *,
    workers: int,
    lower: float = 185.0,
    upper: float = 198.0,
    input_value: float = 187.0,
    starts: int = 4,
) -> dict[str, Any]:
    return _core._fit_resolved(
        model._native,
        model._active_coordinates(ActiveParameterSet(model._parameters, (_IDENTITY,))),
        [0],
        [{"scale": 100.0, "lower": lower, "upper": upper, "input": input_value}],
        observations,
        {
            "start_count": starts,
            "include_input": True,
            "workers": workers,
            "maximum_iterations": 40,
            "function_tolerance": 1.0e-12,
            "gradient_tolerance": 1.0e-12,
            "parameter_tolerance": 1.0e-12,
        },
    )


def _log_residuals(value: float, observed: list[dict[str, Any]]) -> list[float]:
    model = _model(value)
    residuals: list[float] = []
    for temperature, observation in zip(_TEMPERATURES, observed, strict=True):
        saturation = model.saturation(T=temperature * unit_registry.kelvin)
        targets = observation["targets"]
        assert isinstance(targets, list)
        residuals.extend(
            (
                log(saturation.saturation_pressure_pa / targets[0]["observed"]),
                log(saturation.liquid.molar_density_mol_m3 / targets[1]["observed"]),
            )
        )
    return residuals


def _starts(result: dict[str, Any]) -> list[dict[str, Any]]:
    return cast(list[dict[str, Any]], result["starts"])


def _parameter_value(parameters: Parameters, identity: str) -> float:
    spec = next(item for item in parameters.parameter_specs if item.identity == identity)
    return float(spec.value.to(spec.unit).magnitude)


def _resolved_fit(tmp_path: Any, parameters: Parameters, observations: list[dict[str, Any]]) -> Any:
    parameter_path = tmp_path / "parameters.json"
    parameter_bytes = json.dumps(
        parameters.to_mapping(), separators=(",", ":"), sort_keys=True
    ).encode()
    parameter_path.write_bytes(parameter_bytes)
    rows = []
    for observation in observations:
        targets = observation["targets"]
        rows.append(
            {
                "temperature_k": observation["temperature_k"],
                "saturation_pressure": targets[0]["observed"],
                "saturated_liquid_molar_density": targets[1]["observed"],
            }
        )
    row_digest = (
        "sha256:"
        + hashlib.sha256(
            json.dumps(rows, separators=(",", ":"), sort_keys=True).encode()
        ).hexdigest()
    )
    fit = ParameterFit.from_mapping(
        {
            "schema_version": 1,
            "parameters": {
                "path": str(parameter_path),
                "sha256": "sha256:" + hashlib.sha256(parameter_bytes).hexdigest(),
            },
            "model": {"kind": "pure"},
            "active": {
                "identities": [_IDENTITY],
                "coordinates": {
                    _IDENTITY: {"bounds": [185.0, 198.0], "scale": 100.0, "unit": "kelvin"}
                },
                "ties": [],
            },
            "observations": [
                {
                    "family": "pure_saturation",
                    "source": {
                        "kind": "table",
                        "path": "inline-saturation",
                        "sha256": row_digest,
                        "citation": "Project-authored compact numerical fixture.",
                        "domain": {"temperature_k": [190.0, 210.0]},
                        "rows": rows,
                    },
                    "targets": [
                        "saturation_pressure",
                        "saturated_liquid_molar_density",
                    ],
                    "residuals": {
                        "saturation_pressure": "log_ratio",
                        "saturated_liquid_molar_density": "log_ratio",
                    },
                    "units": {
                        "saturation_pressure": "pascal",
                        "saturated_liquid_molar_density": "mole / meter**3",
                    },
                    "basis": {
                        "saturation_pressure": "absolute",
                        "saturated_liquid_molar_density": "molar",
                    },
                    "reference_convention": {
                        "saturation_pressure": "absolute",
                        "saturated_liquid_molar_density": "saturated_liquid",
                    },
                    "transformation": {
                        "saturation_pressure": {"input": "identity", "output": "identity"},
                        "saturated_liquid_molar_density": {
                            "input": "identity",
                            "output": "identity",
                        },
                    },
                    "multiplier": 1.0,
                }
            ],
            "solver": {
                "multistart": {
                    "count": 4,
                    "include_input": True,
                },
            },
        }
    )
    return fit.resolve(parameters, base_path=tmp_path)


def test_exact_generic_fit_recovers_parameter_and_is_schedule_deterministic() -> None:
    model = _model()
    observations = _observations(model)
    first = _fit(model, observations, workers=2)
    second = _fit(model, observations, workers=2)
    first_starts = _starts(first)
    second_starts = _starts(second)
    assert [item["usable"] for item in first_starts] == [item["usable"] for item in second_starts]
    assert [item["final_physical"][0] for item in first_starts] == pytest.approx(
        [item["final_physical"][0] for item in second_starts], rel=0.0, abs=1.0e-12
    )
    usable = [item for item in first_starts if item["usable"]]
    assert usable
    assert all(
        item["final_physical"][0] == pytest.approx(_TRUE_VALUE, abs=1.0e-7) for item in usable
    )
    assert {item["basin"] for item in usable} == {0}
    assert all(item["diagnostics"]["rank"] == 1 for item in usable)
    assert all(item["active_bounds"] == [] for item in usable)

    winner = min(usable, key=lambda item: item["cost"])
    step = 1.0e-4
    below = _log_residuals(_TRUE_VALUE - step, observations)
    above = _log_residuals(_TRUE_VALUE + step, observations)
    finite_difference = [(high - low) / (2.0 * step) for low, high in zip(below, above)]
    assert winner["physical_jacobian"] == pytest.approx(finite_difference, rel=2.0e-5, abs=1.0e-9)


def test_fixed_state_pressure_recovers_generic_binary_parameter_with_exact_jacobian() -> None:
    base = Parameters.from_mapping(NEUTRAL_BINARY)

    def model(value: float) -> Mixture:
        return Mixture(
            base.with_values(
                {_PAIR_IDENTITY: value},
                source_id="pair-fit-fixture",
                locator="project-authored numerical fixture",
            )
        )

    states = ((240.0, 800.0), (260.0, 1500.0), (280.0, 3000.0), (300.0, 5000.0))
    composition = (0.4, 0.6)
    target_model = model(_PAIR_TRUE_VALUE)
    observations: list[dict[str, Any]] = []
    for temperature, density in states:
        pressure = float(
            target_model.state(
                T=temperature * unit_registry.kelvin,
                rho=density * unit_registry.mole / unit_registry.meter**3,
                x=composition,
            )
            .pressure.to("pascal")
            .magnitude
        )
        observations.append(
            {
                "family": "fixed_state_pressure",
                "temperature_k": temperature,
                "molar_density_mol_m3": density,
                "composition": composition,
                "multiplier": 1.0,
                "targets": [
                    {
                        "label": "pressure_pa",
                        "observed": pressure,
                        "residual": "log_ratio",
                        "scale": pressure,
                        "value_factor": 1.0,
                    }
                ],
            }
        )

    input_model = model(0.0)
    native = _core._fit_resolved(
        input_model._native,
        input_model._active_coordinates(
            ActiveParameterSet(input_model._parameters, (_PAIR_IDENTITY,))
        ),
        [0],
        [{"scale": 0.05, "lower": -0.1, "upper": 0.1, "input": 0.0}],
        observations,
        {
            "start_count": 4,
            "include_input": True,
            "workers": 2,
            "maximum_iterations": 40,
            "function_tolerance": 1.0e-12,
            "gradient_tolerance": 1.0e-12,
            "parameter_tolerance": 1.0e-12,
        },
    )
    usable = [item for item in _starts(native) if item["usable"]]
    assert usable
    assert all(
        item["final_physical"][0] == pytest.approx(_PAIR_TRUE_VALUE, abs=1.0e-10) for item in usable
    )

    def residuals(value: float) -> list[float]:
        candidate = model(value)
        values: list[float] = []
        for (temperature, density), observation in zip(states, observations, strict=True):
            pressure = float(
                candidate.state(
                    T=temperature * unit_registry.kelvin,
                    rho=density * unit_registry.mole / unit_registry.meter**3,
                    x=composition,
                )
                .pressure.to("pascal")
                .magnitude
            )
            values.append(log(pressure / observation["targets"][0]["observed"]))
        return values

    step = 1.0e-5
    below = residuals(_PAIR_TRUE_VALUE - step)
    above = residuals(_PAIR_TRUE_VALUE + step)
    finite_difference = [(high - low) / (2.0 * step) for low, high in zip(below, above)]
    winner = min(usable, key=lambda item: item["cost"])
    assert winner["physical_jacobian"] == pytest.approx(finite_difference, rel=3.0e-5, abs=1.0e-9)


def test_scaled_difference_uses_declared_scale_and_failed_start_has_no_penalty() -> None:
    model = _model()
    scaled = _starts(
        _fit(
            model,
            _observations(model, "scaled_difference"),
            workers=1,
            lower=_TRUE_VALUE - 0.01,
            upper=_TRUE_VALUE + 0.01,
            input_value=_TRUE_VALUE,
            starts=1,
        )
    )[0]
    assert scaled["usable"] is True
    assert scaled["residuals"] == pytest.approx([0.0] * 6, abs=1.0e-12)

    failed = _starts(
        _fit(
            model,
            _observations(model),
            workers=1,
            lower=20.0,
            upper=210.0,
            input_value=20.0,
            starts=1,
        )
    )[0]
    assert failed["usable"] is False
    assert failed["non_evaluable_count"] >= 1
    assert failed["failure_reason"]
    assert isinf(failed["cost"])


def test_fit_result_transfers_the_accepted_parameter_without_mutating_input(tmp_path: Any) -> None:
    input_model = _model(187.0)
    observations = _observations(_model())
    native = _fit(input_model, observations, workers=2)
    resolved = _resolved_fit(tmp_path, input_model._parameters, observations)
    result = FitResult.from_native(
        resolved,
        native,
        artifact_fingerprint="sha256:" + "a" * 64,
        package_fingerprint="sha256:" + "b" * 64,
    )
    assert result.accepted_start is not None
    fitted = result.to_parameters()
    assert _parameter_value(input_model._parameters, _IDENTITY) == 187.0
    assert _parameter_value(fitted, _IDENTITY) == pytest.approx(_TRUE_VALUE, abs=1.0e-7)
