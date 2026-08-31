from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from typing import Any, SupportsFloat, cast

import epcsaft
from epcsaft import equilibrium
from epcsaft.eos import Mixture, Pure
from epcsaft.parameters import Parameters
from tests.eos._cases import MULTINEUTRAL_REACTIVE, NEUTRAL_BINARY, neutral_ternary_pairs
from tests.equilibrium._cases import (
    multineutral_reaction,
    provider_boundary_case,
    solved_pressure_case,
)

PAIR_IDENTITY = "pair/ethane/methane/k_ij"
PURE_IDENTITY = "component/ethane/dispersion_energy_over_k"
PURE_TRUTH = 191.42


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode()


def sha256(value: object) -> str:
    return "sha256:" + hashlib.sha256(canonical_bytes(value)).hexdigest()


def write_parameters(tmp_path: Path, parameters: Parameters, name: str) -> tuple[str, str]:
    path = tmp_path / name
    payload = canonical_bytes(parameters.to_mapping())
    path.write_bytes(payload)
    return str(path), "sha256:" + hashlib.sha256(payload).hexdigest()


def parameter_value(parameters: Parameters, identity: str) -> float:
    specification = next(item for item in parameters.parameter_specs if item.identity == identity)
    value = specification.value
    return (
        float(value.to(specification.unit).magnitude)
        if hasattr(value, "to")
        else float(cast(SupportsFloat, value))
    )


def with_parameter(parameters: Parameters, identity: str, value: float, source: str) -> Parameters:
    specification = next(item for item in parameters.parameter_specs if item.identity == identity)
    return parameters.with_values(
        {identity: epcsaft.unit_registry.Quantity(value, specification.unit)},
        source_id=source,
        locator="project-authored manufactured fit case",
    )


def pressure(
    parameters: Parameters,
    temperature_k: float,
    density_mol_m3: float,
    composition: Sequence[float],
) -> float:
    return float(
        Mixture(parameters)
        .state(
            T=temperature_k * epcsaft.unit_registry.kelvin,
            rho=density_mol_m3 * epcsaft.unit_registry.mole / epcsaft.unit_registry.meter**3,
            x=composition,
        )
        .pressure.to("pascal")
        .magnitude
    )


def fixed_state_observation(
    rows: Sequence[Mapping[str, float]],
    composition: Sequence[float],
    *,
    identity: str,
    residual: str = "scaled_difference",
    multiplier: float = 1.0,
) -> dict[str, object]:
    temperatures = tuple(row["temperature_k"] for row in rows)
    densities = tuple(row["molar_density_mol_per_m3"] for row in rows)
    frozen_rows = [dict(row) for row in rows]
    return {
        "family": "fixed_state_pressure",
        "composition": list(composition),
        "source": {
            "kind": "table",
            "path": f"inline-{identity}",
            "sha256": sha256(frozen_rows),
            "citation": "Project-authored manufactured fixed-state case.",
            "domain": {
                "temperature_k": [min(temperatures) - 1.0, max(temperatures) + 1.0],
                "molar_density_mol_per_m3": [min(densities) - 1.0, max(densities) + 1.0],
            },
            "rows": frozen_rows,
        },
        "targets": ["pressure"],
        "residuals": {"pressure": residual},
        "units": {"pressure": "pascal"},
        "basis": {"pressure": "absolute"},
        "reference_convention": {"pressure": "absolute"},
        "transformation": {"pressure": {"input": "identity", "output": "identity"}},
        "multiplier": multiplier,
    }


def fixed_state_fit_case(
    tmp_path: Path,
    *,
    truth: float = 0.04,
    initial: float = -0.02,
    starts: Sequence[tuple[str, float]] | None = None,
) -> tuple[
    dict[str, object],
    Parameters,
    tuple[tuple[float, float, tuple[float, float], float], ...],
]:
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    truth_parameters = with_parameter(parameters, PAIR_IDENTITY, truth, "fixed-state-truth")
    states = tuple(
        (
            temperature,
            density,
            composition,
            pressure(truth_parameters, temperature, density, composition),
        )
        for temperature, density, composition in (
            (285.0, 450.0, (0.40, 0.60)),
            (310.0, 1050.0, (0.40, 0.60)),
            (335.0, 1800.0, (0.40, 0.60)),
        )
    )
    rows = [
        {
            "temperature_k": temperature,
            "molar_density_mol_per_m3": density,
            "pressure": observed,
        }
        for temperature, density, _, observed in states
    ]
    input_parameters = with_parameter(parameters, PAIR_IDENTITY, initial, "fixed-state-input")
    path, digest = write_parameters(tmp_path, input_parameters, "fixed-state-parameters.json")
    declared_starts = starts or (("input", initial), ("confirmation", 0.08))
    mapping: dict[str, object] = {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "mixture"},
        "active": {
            "identities": [PAIR_IDENTITY],
            "coordinates": {
                PAIR_IDENTITY: {"bounds": [-0.10, 0.10], "origin": initial, "scale": 0.05}
            },
            "ties": [],
        },
        "observations": [
            fixed_state_observation(rows, states[0][2], identity="fixed-state-recovery")
        ],
        "solver": {
            "multistart": {
                "starts": [
                    {"identity": identity, "values": [value]} for identity, value in declared_starts
                ],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {"maximum_iterations": 25, "optimizer_tolerance": 1.0e-12},
        },
    }
    return mapping, input_parameters, states


def pure_parameters(value: float) -> Parameters:
    mapping = deepcopy(NEUTRAL_BINARY)
    components = cast(list[dict[str, Any]], mapping["components"])
    mapping["components"] = [
        component for component in components if component["component_id"] == "ethane"
    ]
    mapping["pairs"] = []
    mapping["document_id"] = "pure-saturation-fit-fixture"
    coefficient = next(
        coefficient
        for component in cast(list[dict[str, Any]], mapping["components"])
        for coefficient in cast(list[dict[str, Any]], component["coefficients"])
        if coefficient["identity"] == PURE_IDENTITY
    )
    cast(dict[str, object], coefficient["value"])["magnitude"] = value
    return Parameters.from_mapping(mapping)


def pure_saturation_fit_case(
    tmp_path: Path,
    *,
    initial: float = 187.0,
    starts: Sequence[tuple[str, float]] | None = None,
) -> tuple[dict[str, object], Parameters, tuple[float, ...]]:
    truth_model = Pure(pure_parameters(PURE_TRUTH))
    rows: list[dict[str, float]] = []
    predictions: list[float] = []
    for temperature in (195.0, 210.0):
        saturation = truth_model.saturation(T=temperature * epcsaft.unit_registry.kelvin)
        rows.append(
            {
                "temperature_k": temperature,
                "saturation_pressure": saturation.saturation_pressure_pa,
                "saturated_liquid_molar_density": saturation.liquid.molar_density_mol_m3,
            }
        )
        predictions.extend(
            (saturation.saturation_pressure_pa, saturation.liquid.molar_density_mol_m3)
        )
    parameters = pure_parameters(initial)
    path, digest = write_parameters(tmp_path, parameters, "pure-saturation-parameters.json")
    declared_starts = starts or (("truth", PURE_TRUTH),)
    mapping: dict[str, object] = {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "pure"},
        "active": {
            "identities": [PURE_IDENTITY],
            "coordinates": {
                PURE_IDENTITY: {
                    "bounds": [170.0, 205.0],
                    "origin": initial,
                    "scale": 20.0,
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
                    "path": "inline-pure-saturation",
                    "sha256": sha256(rows),
                    "citation": "Project-authored manufactured saturation case.",
                    "domain": {"temperature_k": [195.0, 210.0]},
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
                    {"identity": identity, "values": [value]} for identity, value in declared_starts
                ],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {"maximum_iterations": 20, "optimizer_tolerance": 1.0e-10},
        },
    }
    return mapping, parameters, tuple(predictions)


def equilibrium_fit_mapping(
    tmp_path: Path,
    parameters: Parameters,
    problem: equilibrium.GeneralReactiveEquilibriumProblem,
    parameter_identity: str,
    observed: float,
    *,
    bounds: tuple[float, float],
    origin: float,
    scale: float,
    starts: Sequence[tuple[str, float]],
    target_scale: float,
) -> dict[str, object]:
    output = problem.outputs[0]
    payload_path, digest = write_parameters(tmp_path, parameters, "equilibrium-parameters.json")
    return {
        "schema_version": 2,
        "parameters": {"path": payload_path, "sha256": digest},
        "model": {"kind": "mixture"},
        "active": {
            "identities": [parameter_identity],
            "coordinates": {
                parameter_identity: {
                    "bounds": list(bounds),
                    "origin": origin,
                    "scale": scale,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "equilibrium",
                "identity": "manufactured-equilibrium-batch",
                "request": equilibrium.general_reactive_equilibrium_problem_to_mapping(problem),
                "targets": [
                    {
                        "identity": "manufactured-equilibrium-row",
                        "prediction_identity": output.identity,
                        "observed": observed,
                        "unit": output.unit,
                        "residual": "scaled_difference",
                        "scale": target_scale,
                        "multiplier": 1.0,
                        "role": "verification",
                        "measurement_classification": "manufactured",
                        "basis": output.basis,
                        "source_identity": "manufactured-equilibrium-source",
                        "source_hash": f"sha256:{'5' * 64}",
                        "included": True,
                        "drop_reason": None,
                        "aggregate_identity": output.aggregate_identity,
                        "covariance_identity": output.covariance_identity,
                    }
                ],
                "multiplier": 1.0,
            }
        ],
        "solver": {
            "multistart": {
                "starts": [{"identity": identity, "values": [value]} for identity, value in starts],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {"maximum_iterations": 25, "optimizer_tolerance": 1.0e-8},
        },
    }


def full_space_fit_case(
    tmp_path: Path,
) -> tuple[
    dict[str, object],
    Parameters,
    Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
    str,
    float,
    equilibrium.GeneralReactiveEquilibriumResult,
]:
    parameters = Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = Mixture(parameters)
    reaction = multineutral_reaction()
    phase = equilibrium.ReactivePhase(
        "finite-liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
    )
    output = equilibrium.EquilibriumOutput(
        "cation-fraction",
        "phase.mole_fraction",
        "dimensionless",
        "true-species-mole-fraction",
        phase.identity,
        (0.0, 0.0, 1.0, 0.0),
        support="positive",
    )
    unanchored = equilibrium.GeneralReactiveEquilibriumProblem(
        "finite-phase-fit",
        equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(119_199.78746707023 * epcsaft.unit_registry.pascal),
        (phase,),
        reaction,
        (phase.identity,),
        (output,),
    )
    anchor = equilibrium.solve(model, unanchored)
    assert anchor.status == "evaluated" and anchor.continuation_state is not None
    problem = replace(
        unanchored,
        continuation_identity="finite-liquid-branch",
        continuation_state=anchor.continuation_state,
    )
    parameter_identity = "component/solvent/segment_diameter"
    truth_value = 2.9
    truth = equilibrium.solve(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters, (parameter_identity,), values=(truth_value,)
        ),
    )
    assert truth.status == "evaluated" and truth.values[0] is not None
    mapping = equilibrium_fit_mapping(
        tmp_path,
        parameters,
        problem,
        parameter_identity,
        float(truth.values[0]),
        bounds=(2.5, 3.1),
        origin=2.8,
        scale=0.2,
        starts=(("off-truth", 2.8),),
        target_scale=float(truth.values[0]),
    )
    return mapping, parameters, model, problem, parameter_identity, truth_value, truth


def solved_pressure_fit_case(
    tmp_path: Path,
    *,
    ideal_vapor: bool,
    temperature_k: float = 313.15,
    reaction_correlation: bool = False,
    reaction_b_upper: float = 0.0,
    parameter_identity: str = "component/solvent/segment_diameter",
    truth_value: float = 2.9,
    bounds: tuple[float, float] = (2.7, 3.05),
    origin: float = 2.8,
    scale: float = 0.1,
) -> tuple[
    dict[str, object],
    Parameters,
    Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
    str,
    float,
    equilibrium.GeneralReactiveEquilibriumResult,
]:
    parameters, model, unanchored = solved_pressure_case(ideal_vapor=ideal_vapor)
    unanchored = replace(
        unanchored,
        identity=f"{unanchored.identity}-{temperature_k:g}-kelvin",
        temperature=equilibrium.Fixed(temperature_k * epcsaft.unit_registry.kelvin),
    )
    if reaction_correlation:
        parameter_mapping = deepcopy(MULTINEUTRAL_REACTIVE)
        parameter_mapping["domains"].append(
            {
                "domain_id": "engine-candidate-313-15-to-393-15-k",
                "kind": "fit-range",
                "temperature_min": {"magnitude": 313.15, "unit": "kelvin"},
                "temperature_max": {"magnitude": 393.15, "unit": "kelvin"},
                "pressure_min": {"magnitude": 1.0, "unit": "pascal"},
                "pressure_max": {"magnitude": 100_000_000.0, "unit": "pascal"},
            }
        )
        parameter_mapping["reaction_correlations"] = [
            {
                "reaction_id": "R1",
                "kind": "ln-k-a-plus-b-over-t",
                "source": {
                    "source_id": "engine-numerical-fixture",
                    "locator": "manufactured constant-plus-reciprocal-temperature reaction",
                    "domain_id": "engine-numerical-fixture-domain",
                },
                "candidate_domain_id": "engine-candidate-313-15-to-393-15-k",
                "qualification": "candidate_extrapolation",
                "source_sha256": f"sha256:{'6' * 64}",
                "coefficients": [
                    {
                        "name": "a",
                        "identity": "reaction:R1:correlation:a",
                        "value": {"magnitude": -240.0, "unit": "dimensionless"},
                        "optimizer_domain": {"lower": -250.0, "upper": -230.0},
                    },
                    {
                        "name": "b_k",
                        "identity": "reaction:R1:correlation:b_k",
                        "value": {"magnitude": -100.0, "unit": "kelvin"},
                        "optimizer_domain": {"lower": -200.0, "upper": reaction_b_upper},
                    },
                ],
            }
        ]
        candidate = {
            "candidate_domain_id": "engine-candidate-313-15-to-393-15-k",
            "qualification": "candidate_extrapolation",
            "source_sha256": f"sha256:{'6' * 64}",
        }
        for family in parameter_mapping["model_families"]:
            family.update(candidate)
        for component in parameter_mapping["components"]:
            for record in (*component["fixed"].values(), *component["coefficients"]):
                record.update(candidate)
        for pair in parameter_mapping["pairs"]:
            for record in pair["coefficients"]:
                record.update(candidate)
        for record in parameter_mapping["model_coefficients"]:
            record.update(candidate)
        for correlation in parameter_mapping["correlations"]:
            correlation.update(candidate)
        for site in parameter_mapping["topology"]["sites"]:
            site.update(candidate)
        for edge in parameter_mapping["topology"]["edges"]:
            edge.update(candidate)
        parameters = Parameters.from_mapping(parameter_mapping)
        model = Mixture(parameters)
        reaction = unanchored.reaction_system
        unanchored = replace(
            unanchored,
            pressure=equilibrium.Solved(
                3.2e6 * epcsaft.unit_registry.pascal,
                (1.0e4 * epcsaft.unit_registry.pascal, 1.0e8 * epcsaft.unit_registry.pascal),
                tuple(
                    value * epcsaft.unit_registry.pascal
                    for value in (1.0e5, 1.0e6, 3.2e6, 1.0e7, 5.0e7)
                ),
            ),
            reaction_system=replace(
                reaction,
                equilibrium_constants=(
                    replace(
                        reaction.equilibrium_constants[0],
                        reaction_id="R1",
                        correlation_kind="ln-k-a-plus-b-over-t",
                        coefficient_identities=(
                            "reaction:R1:correlation:a",
                            "reaction:R1:correlation:b_k",
                        ),
                        coefficient_values=(-240.0, -100.0),
                    ),
                ),
            ),
            outputs=(
                *unanchored.outputs,
                equilibrium.EquilibriumOutput(
                    "liquid-cation-fraction",
                    "phase.mole_fraction",
                    "dimensionless",
                    "true-species-mole-fraction",
                    unanchored.phases[0].identity,
                    (0.0, 0.0, 1.0, 0.0),
                    support="positive",
                ),
            ),
        )
    truth_active = epcsaft.ActiveParameterSet(
        parameters, (parameter_identity,), values=(truth_value,)
    )
    anchor = equilibrium.solve(model, unanchored, active_parameters=truth_active)
    assert anchor.status == "evaluated" and anchor.continuation_state is not None, anchor.failure
    problem = replace(unanchored, continuation_state=anchor.continuation_state)
    truth = equilibrium.solve(model, problem, active_parameters=truth_active)
    assert truth.status == "evaluated" and all(value is not None for value in truth.values)
    path, digest = write_parameters(tmp_path, parameters, "solved-pressure-parameters.json")
    residuals = ("log_ratio", "scaled_difference", "log_ratio") + (
        ("scaled_difference",) if reaction_correlation else ()
    )
    scales = (0.05, 0.01, 0.05) + ((0.01,) if reaction_correlation else ())
    mapping: dict[str, object] = {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "mixture"},
        "active": {
            "identities": [parameter_identity],
            "coordinates": {
                parameter_identity: {
                    "bounds": list(bounds),
                    "origin": origin,
                    "scale": scale,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "equilibrium",
                "identity": problem.identity,
                "request": equilibrium.general_reactive_equilibrium_problem_to_mapping(problem),
                "targets": [
                    {
                        "identity": f"manufactured-{output.identity}",
                        "prediction_identity": output.identity,
                        "observed": float(observed),
                        "unit": output.unit,
                        "residual": residual,
                        "scale": scale,
                        "multiplier": 1.0,
                        "role": "verification",
                        "measurement_classification": "manufactured",
                        "basis": output.basis,
                        "source_identity": "manufactured-solved-pressure",
                        "source_hash": f"sha256:{'4' * 64}",
                        "included": True,
                        "drop_reason": None,
                        "aggregate_identity": None,
                        "covariance_identity": None,
                    }
                    for output, observed, residual, scale in zip(
                        problem.outputs, truth.values, residuals, scales, strict=True
                    )
                ],
                "multiplier": 1.0,
            }
        ],
        "solver": {
            "multistart": {
                "starts": [{"identity": "off-truth", "values": [origin]}],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {"maximum_iterations": 80, "optimizer_tolerance": 1.0e-8},
        },
    }
    return mapping, parameters, model, problem, parameter_identity, truth_value, truth


def boundary_fit_case(
    tmp_path: Path,
    *,
    parameter_identity: str | None = None,
    truth_value: float = 0.08,
    bounds: tuple[float, float] = (-0.05, 0.15),
    origin: float = 0.0,
    scale: float = 0.1,
    starts: Sequence[tuple[str, float]] = (("truth", 0.08),),
) -> tuple[
    dict[str, object],
    Parameters,
    Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
    str,
    equilibrium.GeneralReactiveEquilibriumResult,
]:
    parameters, model, problem, default_identity = provider_boundary_case()
    output = replace(
        problem.outputs[0],
        aggregate_identity="boundary-output-group",
        covariance_identity="diagonal-whitening-group",
    )
    problem = replace(problem, outputs=(output,))
    identity = default_identity if parameter_identity is None else parameter_identity
    truth = equilibrium.solve(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters, (identity,), values=(truth_value,)
        ),
    )
    assert truth.status == "evaluated" and truth.values[0] is not None
    mapping = equilibrium_fit_mapping(
        tmp_path,
        parameters,
        problem,
        identity,
        float(truth.values[0]),
        bounds=bounds,
        origin=origin,
        scale=scale,
        starts=starts,
        target_scale=float(truth.values[0]),
    )
    return mapping, parameters, model, problem, identity, truth


def tied_coordinate_fit_case(
    tmp_path: Path,
) -> tuple[
    dict[str, object],
    Parameters,
    tuple[tuple[float, float, tuple[float, ...], float], ...],
    tuple[str, str, str],
]:
    identities = (
        "pair/beta/gamma/k_ij",
        "component/gamma/dispersion_energy_over_k",
        "pair/alpha/beta/k_ij",
    )
    pair_a, energy, pair_b = identities
    base = Parameters.from_mapping(neutral_ternary_pairs())
    truth = base.with_values(
        {pair_a: 0.05, pair_b: 0.05},
        source_id="tied-coordinate-truth",
        locator="project-authored manufactured fit case",
    )
    composition = (0.45, 0.35, 0.20)
    states = tuple(
        (
            temperature,
            density,
            composition,
            pressure(truth, temperature, density, composition),
        )
        for temperature, density in ((305.0, 500.0), (318.0, 1100.0), (330.0, 1700.0))
    )
    rows = [
        {
            "temperature_k": temperature,
            "molar_density_mol_per_m3": density,
            "pressure": observed,
        }
        for temperature, density, _, observed in states
    ]
    input_parameters = base.with_values(
        {
            pair_a: 0.0,
            energy: 270.0 * epcsaft.unit_registry.kelvin,
            pair_b: 0.0,
        },
        source_id="tied-coordinate-input",
        locator="project-authored manufactured fit case",
    )
    path, digest = write_parameters(tmp_path, input_parameters, "tied-parameters.json")
    mapping: dict[str, object] = {
        "schema_version": 2,
        "parameters": {"path": path, "sha256": digest},
        "model": {"kind": "mixture"},
        "active": {
            "identities": list(identities),
            "coordinates": {
                pair_a: {"bounds": [-0.08, 0.08], "origin": 0.0, "scale": 0.05},
                energy: {
                    "bounds": [250.0, 330.0],
                    "origin": 270.0,
                    "scale": 300.0,
                    "unit": "kelvin",
                },
            },
            "ties": [{"identities": [pair_a, pair_b]}],
        },
        "observations": [fixed_state_observation(rows, composition, identity="tied-coordinates")],
        "solver": {
            "multistart": {
                "starts": [{"identity": "input", "values": [0.0, 270.0]}],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {"maximum_iterations": 40, "optimizer_tolerance": 1.0e-12},
        },
    }
    return mapping, input_parameters, states, identities


def scaled_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(2.0e-9, abs(actual), abs(expected))
