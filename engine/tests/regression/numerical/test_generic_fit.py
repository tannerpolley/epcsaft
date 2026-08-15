from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any, SupportsFloat, cast

import pytest

from epcsaft.eos import Mixture, Pure
from epcsaft.parameters import Parameters
from epcsaft.regression import FitPreflightError, FitResult, ParameterFit, fit
from epcsaft.units import unit_registry
from tests.eos._raw_cases import NEUTRAL_BINARY
from tests.eos.numerical.test_parameter_derivatives import _neutral_ternary_pairs

_PURE_IDENTITY = "component/ethane/dispersion_energy_over_k"
_PURE_TRUE_VALUE = 191.42
_MIXED_IDENTITIES = (
    "pair/beta/gamma/k_ij",
    "component/gamma/dispersion_energy_over_k",
    "pair/alpha/beta/k_ij",
)
_MIXED_TRUE_VALUES = (-0.05, 290.0, 0.03)


def _canonical(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode()


def _sha256(value: object) -> str:
    return "sha256:" + hashlib.sha256(_canonical(value)).hexdigest()


def _write_parameters(tmp_path: Any, parameters: Parameters, name: str) -> tuple[str, str]:
    path = tmp_path / name
    payload = _canonical(parameters.to_mapping())
    path.write_bytes(payload)
    return str(path), "sha256:" + hashlib.sha256(payload).hexdigest()


def _parameter_value(parameters: Parameters, identity: str) -> float:
    specification = next(item for item in parameters.parameter_specs if item.identity == identity)
    value = specification.value
    return (
        float(value.to(specification.unit).magnitude)
        if hasattr(value, "to")
        else float(cast(SupportsFloat, value))
    )


def _pressure(
    model: Mixture, temperature: float, density: float, composition: Sequence[float]
) -> float:
    return float(
        model.state(
            T=temperature * unit_registry.kelvin,
            rho=density * unit_registry.mole / unit_registry.meter**3,
            x=composition,
        )
        .pressure.to("pascal")
        .magnitude
    )


def _fixed_state_observation(
    rows: list[dict[str, float]], composition: tuple[float, ...], index: int
) -> dict[str, object]:
    temperatures = tuple(row["temperature_k"] for row in rows)
    densities = tuple(row["molar_density_mol_per_m3"] for row in rows)
    return {
        "family": "fixed_state_pressure",
        "composition": list(composition),
        "source": {
            "kind": "table",
            "path": f"inline-fixed-state-{index}",
            "sha256": _sha256(rows),
            "citation": "Project-authored compact numerical fixture.",
            "domain": {
                "temperature_k": [min(temperatures) - 1.0, max(temperatures) + 1.0],
                "molar_density_mol_per_m3": [min(densities) - 1.0, max(densities) + 1.0],
            },
            "rows": rows,
        },
        "targets": ["pressure"],
        "residuals": {"pressure": "log_ratio"},
        "units": {"pressure": "pascal"},
        "basis": {"pressure": "absolute"},
        "reference_convention": {"pressure": "absolute"},
        "transformation": {"pressure": {"input": "identity", "output": "identity"}},
        "multiplier": 1.0,
    }


def _mixed_fit_mapping(
    tmp_path: Any,
) -> tuple[
    dict[str, object], Parameters, tuple[tuple[float, float, tuple[float, ...], float], ...]
]:
    true_parameters = Parameters.from_mapping(_neutral_ternary_pairs())
    true_model = Mixture(true_parameters)
    state_groups = (
        ((305.0, 300.0), (325.0, 1800.0), (0.15, 0.25, 0.60)),
        ((310.0, 600.0), (335.0, 2000.0), (0.55, 0.30, 0.15)),
        ((300.0, 1000.0), (330.0, 1500.0), (0.25, 0.60, 0.15)),
    )
    observations: list[dict[str, object]] = []
    states: list[tuple[float, float, tuple[float, ...], float]] = []
    for index, (first, second, composition) in enumerate(state_groups):
        rows: list[dict[str, float]] = []
        for temperature, density in (first, second):
            pressure = _pressure(true_model, temperature, density, composition)
            rows.append(
                {
                    "temperature_k": temperature,
                    "molar_density_mol_per_m3": density,
                    "pressure": pressure,
                }
            )
            states.append((temperature, density, composition, pressure))
        observations.append(_fixed_state_observation(rows, composition, index))

    input_parameters = true_parameters.with_values(
        {
            _MIXED_IDENTITIES[0]: -0.01,
            _MIXED_IDENTITIES[1]: 270.0 * unit_registry.kelvin,
            _MIXED_IDENTITIES[2]: 0.0,
        },
        source_id="regression-mixed-input",
        locator="project-authored numerical fixture",
    )
    path, digest = _write_parameters(tmp_path, input_parameters, "mixed-parameters.json")
    mapping: dict[str, object] = {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "mixture"},
        "active": {
            "identities": list(_MIXED_IDENTITIES),
            "coordinates": {
                _MIXED_IDENTITIES[0]: {
                    "bounds": [-0.12, 0.08],
                    "origin": -0.01,
                    "scale": 0.05,
                },
                _MIXED_IDENTITIES[1]: {
                    "bounds": [250.0, 330.0],
                    "origin": 270.0,
                    "scale": 300.0,
                    "unit": "kelvin",
                },
                _MIXED_IDENTITIES[2]: {
                    "bounds": [-0.08, 0.12],
                    "origin": 0.0,
                    "scale": 0.05,
                },
            },
            "ties": [],
        },
        "observations": observations,
        "solver": {
            "multistart": {
                "starts": [
                    {"identity": "input", "values": [-0.01, 270.0, 0.0]},
                    {"identity": "confirmation-a", "values": [-0.06, 300.0, 0.05]},
                    {"identity": "confirmation-b", "values": [0.04, 285.0, -0.05]},
                ],
                "affine_bound_margin": 1.0e-7,
            }
        },
    }
    return mapping, input_parameters, tuple(states)


def _complete_residuals(
    parameters: Parameters,
    states: Sequence[tuple[float, float, tuple[float, ...], float]],
    values: Mapping[str, float],
    *,
    residual: str = "log_ratio",
) -> tuple[float, ...]:
    candidate = parameters.with_values(
        {
            identity: (
                value * unit_registry.kelvin
                if identity == "component/gamma/dispersion_energy_over_k"
                else value
            )
            for identity, value in values.items()
        },
        source_id="regression-centered-difference",
        locator="independent complete-residual oracle",
    )
    model = Mixture(candidate)
    predictions = tuple(
        (_pressure(model, temperature, density, composition), observed)
        for temperature, density, composition, observed in states
    )
    if residual == "log_ratio":
        return tuple(math.log(predicted / observed) for predicted, observed in predictions)
    if residual == "scaled_difference":
        return tuple((predicted - observed) / observed for predicted, observed in predictions)
    raise ValueError(f"unknown residual: {residual}")


def test_owner_blocks_preserve_mixed_active_order_and_complete_residual_jacobian(
    tmp_path: Any,
) -> None:
    """Falsify identity/order loss and residual scaling at the EOS--Ceres seam."""

    mapping, input_parameters, states = _mixed_fit_mapping(tmp_path)
    result = fit(ParameterFit.from_mapping(mapping), base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert (
        tuple(
            identity
            for coordinate in result.specification.coordinates
            for identity in coordinate.parameter_identities
        )
        == _MIXED_IDENTITIES
    )
    assert accepted.final_physical == pytest.approx(_MIXED_TRUE_VALUES, rel=2.0e-8, abs=2.0e-8)
    assert accepted.predictions == pytest.approx(
        [state[3] for state in states], rel=2.0e-9, abs=1.0e-5
    )
    assert accepted.rank == 3
    assert math.isfinite(accepted.condition_number)

    evaluated = result.evaluated_rows
    assert len(evaluated) == len(states)
    assert tuple(row.predicted for row in evaluated) == pytest.approx(
        accepted.predictions, rel=0.0, abs=0.0
    )
    assert tuple(row.scaled_residual for row in evaluated) == pytest.approx(
        accepted.residuals, rel=0.0, abs=0.0
    )
    assert all(row.evaluation_status == "evaluated" for row in evaluated)
    assert all(row.selector == "pressure" for row in evaluated)
    assert all(row.unit == "pascal" and row.basis == "absolute" for row in evaluated)
    assert result.evaluated_row_table_hash.startswith("sha256:")

    candidate = result.export_candidate_parameters()
    assert input_parameters.fingerprint == result.specification.parameters.fingerprint
    assert candidate.fingerprint != input_parameters.fingerprint
    for identity, expected in zip(_MIXED_IDENTITIES, accepted.final_physical, strict=True):
        assert _parameter_value(candidate, identity) == pytest.approx(expected)
    assert all(
        row.candidate_parameter_identity == result.candidate_parameter_identity
        for row in evaluated
    )

    point = dict(zip(_MIXED_IDENTITIES, accepted.final_physical, strict=True))
    steps = (1.0e-5, 1.0e-3, 1.0e-5)
    columns: list[tuple[float, ...]] = []
    for identity, step in zip(_MIXED_IDENTITIES, steps, strict=True):
        below = _complete_residuals(
            input_parameters, states, {**point, identity: point[identity] - step}
        )
        above = _complete_residuals(
            input_parameters, states, {**point, identity: point[identity] + step}
        )
        columns.append(
            tuple((high - low) / (2.0 * step) for low, high in zip(below, above, strict=True))
        )
    centered = tuple(
        columns[column][row]
        for row in range(len(states))
        for column in range(len(_MIXED_IDENTITIES))
    )
    # Pair steps of 1e-5 and the 1e-3 K component step balance O(h^2)
    # truncation against pressure cancellation. Failure means Regression
    # changed owner column order or scaled values and derivatives differently.
    assert accepted.physical_jacobian == pytest.approx(centered, rel=5.0e-5, abs=2.0e-9)


def test_tied_owner_coordinate_reports_interior_covariance_and_active_bound(
    tmp_path: Any,
) -> None:
    """Falsify tied derivatives, free covariance, and bound ambiguity."""

    pair_beta_gamma, gamma_energy, pair_alpha_beta = _MIXED_IDENTITIES
    true_parameters = Parameters.from_mapping(_neutral_ternary_pairs()).with_values(
        {pair_beta_gamma: 0.05, pair_alpha_beta: 0.05},
        source_id="regression-tied-truth",
        locator="project-authored numerical fixture",
    )
    composition = (0.45, 0.35, 0.2)
    states: list[tuple[float, float, tuple[float, ...], float]] = []
    rows: list[dict[str, float]] = []
    for temperature, density in ((305.0, 500.0), (318.0, 1100.0), (330.0, 1700.0)):
        pressure = _pressure(Mixture(true_parameters), temperature, density, composition)
        states.append((temperature, density, composition, pressure))
        rows.append(
            {
                "temperature_k": temperature,
                "molar_density_mol_per_m3": density,
                "pressure": pressure,
            }
        )
    observation = _fixed_state_observation(rows, composition, 0)
    cast(dict[str, object], observation)["residuals"] = {"pressure": "scaled_difference"}
    input_parameters = true_parameters.with_values(
        {
            pair_beta_gamma: 0.0,
            gamma_energy: 270.0 * unit_registry.kelvin,
            pair_alpha_beta: 0.0,
        },
        source_id="regression-tied-input",
        locator="project-authored numerical fixture",
    )
    path, digest = _write_parameters(tmp_path, input_parameters, "tied-parameters.json")
    mapping = {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "mixture"},
        "active": {
            "identities": list(_MIXED_IDENTITIES),
            "coordinates": {
                pair_beta_gamma: {"bounds": [-0.08, 0.08], "origin": 0.0, "scale": 0.05},
                gamma_energy: {
                    "bounds": [250.0, 330.0],
                    "origin": 270.0,
                    "scale": 300.0,
                    "unit": "kelvin",
                },
            },
            "ties": [{"identities": [pair_beta_gamma, pair_alpha_beta]}],
        },
        "observations": [observation],
        "solver": {
            "multistart": {
                "starts": [{"identity": "input", "values": [0.0, 270.0]}],
                "affine_bound_margin": 1.0e-7,
            }
        },
    }
    interior = fit(mapping, base_path=tmp_path)
    interior_start = interior.best_usable_start
    assert interior_start is not None
    assert interior_start.final_physical[0] == pytest.approx(0.05, abs=2.0e-8)
    assert interior_start.rank == 2
    assert interior_start.active_bounds == ()
    assert interior_start.covariance_status == "available"
    assert interior.covariance_coordinate_identities == (pair_beta_gamma, gamma_energy)
    assert len(interior.physical_covariance) == 4

    pair_coordinate = cast(
        dict[str, object],
        cast(dict[str, object], cast(dict[str, object], mapping["active"])["coordinates"])[
            pair_beta_gamma
        ],
    )
    pair_coordinate["bounds"] = [-0.08, 0.04]
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical[0] == pytest.approx(0.04, abs=1.0e-12)
    assert accepted.active_bounds == (0,)
    assert accepted.rank == 2
    assert math.isfinite(accepted.condition_number)
    assert accepted.covariance_status == "unavailable_active_bounds"

    point = accepted.final_physical
    coordinate_values = (
        {pair_beta_gamma: point[0], pair_alpha_beta: point[0]},
        {gamma_energy: point[1]},
    )
    steps = (1.0e-5, 1.0e-3)
    columns: list[tuple[float, ...]] = []
    for values, step in zip(coordinate_values, steps, strict=True):
        baseline = {
            pair_beta_gamma: point[0],
            pair_alpha_beta: point[0],
            gamma_energy: point[1],
        }
        below = _complete_residuals(
            input_parameters,
            states,
            {
                **baseline,
                **{identity: value - step for identity, value in values.items()},
            },
            residual="scaled_difference",
        )
        above = _complete_residuals(
            input_parameters,
            states,
            {
                **baseline,
                **{identity: value + step for identity, value in values.items()},
            },
            residual="scaled_difference",
        )
        columns.append(
            tuple((high - low) / (2.0 * step) for low, high in zip(below, above, strict=True))
        )
    centered = tuple(
        columns[column][row]
        for row in range(len(states))
        for column in range(len(coordinate_values))
    )
    assert accepted.physical_jacobian == pytest.approx(centered, rel=5.0e-5, abs=2.0e-9)
    fitted = {identity: value for identity, value, _ in result.fitted_values()}
    assert fitted[pair_beta_gamma] == fitted[pair_alpha_beta] == pytest.approx(0.04)


def _pure_parameters(value: float) -> Parameters:
    mapping = deepcopy(NEUTRAL_BINARY)
    components = cast(list[dict[str, Any]], mapping["components"])
    mapping["components"] = [
        component for component in components if component["component_id"] == "ethane"
    ]
    mapping["pairs"] = []
    mapping["document_id"] = "pure-ethane-regression-fixture"
    components = cast(list[dict[str, Any]], mapping["components"])
    coefficient = next(
        coefficient
        for component in components
        if component["component_id"] == "ethane"
        for coefficient in cast(list[dict[str, Any]], component["coefficients"])
        if coefficient["identity"] == _PURE_IDENTITY
    )
    cast(dict[str, object], coefficient["value"])["magnitude"] = value
    return Parameters.from_mapping(mapping)


def _pure_fit_mapping(tmp_path: Any, parameters: Parameters, name: str) -> dict[str, object]:
    true_model = Pure(_pure_parameters(_PURE_TRUE_VALUE))
    rows = []
    for temperature in (190.0, 200.0, 210.0):
        saturation = true_model.saturation(T=temperature * unit_registry.kelvin)
        rows.append(
            {
                "temperature_k": temperature,
                "saturation_pressure": saturation.saturation_pressure_pa,
                "saturated_liquid_molar_density": saturation.liquid.molar_density_mol_m3,
            }
        )
    path, digest = _write_parameters(tmp_path, parameters, f"{name}-parameters.json")
    return {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "pure"},
        "active": {
            "identities": [_PURE_IDENTITY],
            "coordinates": {
                _PURE_IDENTITY: {
                    "bounds": [20.0, 210.0],
                    "origin": float(
                        next(
                            item.value.to("kelvin").magnitude
                            for item in parameters.parameter_specs
                            if item.identity == _PURE_IDENTITY
                        )
                    ),
                    "scale": 100.0,
                    "unit": "kelvin",
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "pure_saturation",
                "source": {
                    "kind": "table",
                    "path": f"inline-{name}",
                    "sha256": _sha256(rows),
                    "citation": "Project-authored compact numerical fixture.",
                    "domain": {"temperature_k": [190.0, 210.0]},
                    "rows": rows,
                },
                "targets": ["saturation_pressure", "saturated_liquid_molar_density"],
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
                "starts": [
                    {
                        "identity": "input",
                        "values": [
                            float(
                                next(
                                    item.value.to("kelvin").magnitude
                                    for item in parameters.parameter_specs
                                    if item.identity == _PURE_IDENTITY
                                )
                            )
                        ],
                    }
                ],
                "affine_bound_margin": 1.0e-7,
            }
        },
    }


def test_pure_saturation_owner_transfer_and_non_evaluable_trial_diagnostics(
    tmp_path: Any,
) -> None:
    """Falsify duplicated saturation evaluation and fabricated failed-trial penalties."""

    input_parameters = _pure_parameters(187.0)
    successful = fit(_pure_fit_mapping(tmp_path, input_parameters, "successful"))
    assert isinstance(successful, FitResult)
    assert successful.best_usable_start is not None
    assert successful.best_usable_start.owner_receipts == ()
    assert successful.to_mapping()["best_usable_start"] is not None
    assert successful.best_usable_start.final_physical[0] == pytest.approx(
        _PURE_TRUE_VALUE, abs=1.0e-7
    )
    fitted = successful.export_candidate_parameters()
    assert _parameter_value(input_parameters, _PURE_IDENTITY) == 187.0
    assert _parameter_value(fitted, _PURE_IDENTITY) == pytest.approx(_PURE_TRUE_VALUE, abs=1.0e-7)

    with pytest.raises(FitPreflightError) as error:
        fit(_pure_fit_mapping(tmp_path, _pure_parameters(20.0), "failed"))
    assert error.value.code == "physical_rejection"
