from __future__ import annotations

import copy
import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import cast

import pytest

import epcsaft
from epcsaft import equilibrium, regression
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE, MULTISOLVENT_MULTIION

from .chemical_equilibrium_cases import typed_problem as _typed_problem
from .test_chemical_provider import (
    _held_parameters,
    _held_water_ionization_problem,
    _multineutral_problem,
    _typed_standard_state,
)


def _reaction_mapping(problem: equilibrium.ChemicalEquilibriumProblem) -> dict[str, object]:
    standard = problem.source_standard_state
    return {
        "species_ids": list(problem.species_ids),
        "charges": list(problem.charges),
        "molar_masses_kg_per_mol": list(problem.molar_masses_kg_per_mol),
        "balance_matrix": [list(row) for row in problem.balance_matrix],
        "conserved_totals": list(problem.conserved_totals),
        "reaction_matrix": [list(row) for row in problem.reaction_matrix],
        "feed_amounts_mol": list(problem.feed_amounts_mol),
        "equilibrium_constants": [
            [
                item.ln_value,
                item.source_id,
                item.reference_id,
                item.reaction_orientation,
                item.conversion_id,
                item.dimensionless,
            ]
            for item in problem.equilibrium_constants
        ],
        "strict_interior_amount_floor_mol": problem.strict_interior_amount_floor_mol,
        "source_standard_state": None
        if standard is None
        else {
            "id": standard.id,
            "activity_scale_id": standard.activity_scale_id,
            "log_activity_scale_factors": list(standard.log_activity_scale_factors),
            "reference_pressure_pa": standard.reference_pressure_pa,
            "component_ids": list(standard.source_reference_component_ids),
            "solvent_composition": list(standard.source_reference_solvent_composition),
            "ion_pairs": [
                [list(pair.component_ids), list(pair.stoichiometry)]
                for pair in standard.source_reference_ion_pairs
            ],
            "phase": standard.source_reference_phase,
            "reference_convention": standard.source_reference_convention,
            "activity_convention_id": standard.source_reference_activity_convention_id,
            "standard_molality_mol_per_kg": (
                standard.source_reference_standard_molality_mol_per_kg
            ),
        },
    }


def _row_mapping(row: equilibrium.HomogeneousReactiveObservationRow) -> dict[str, object]:
    return {
        "identity": row.identity,
        "quantity": row.quantity,
        "coefficients": list(row.coefficients),
        "solvent_mass_coefficients_kg_per_mol": list(
            row.solvent_mass_coefficients_kg_per_mol
        ),
        "support": row.support,
        "censor_limit": row.censor_limit,
        "aggregate_identity": row.aggregate_identity,
        "covariance_identity": row.covariance_identity,
    }


def _target_mapping(
    identity: str,
    prediction_identity: str,
    observed: float,
    unit: str,
    residual: str,
    scale: float,
    *,
    included: bool = True,
    drop_reason: str | None = None,
    aggregate_identity: str | None = None,
    covariance_identity: str | None = None,
) -> dict[str, object]:
    return {
        "identity": identity,
        "prediction_identity": prediction_identity,
        "observed": observed,
        "unit": unit,
        "residual": residual,
        "scale": scale,
        "multiplier": 1.0,
        "role": "diagnostic",
        "measurement_classification": "direct",
        "source_identity": "manufactured-reactive-observation",
        "source_hash": f"sha256:{'2' * 64}",
        "included": included,
        "drop_reason": drop_reason,
        "aggregate_identity": aggregate_identity,
        "covariance_identity": covariance_identity,
    }


def _reactive_fit_mapping(
    tmp_path: Path,
    parameters: epcsaft.Parameters,
    observation_problem: equilibrium.HomogeneousReactiveObservationProblem,
    rows: tuple[equilibrium.HomogeneousReactiveObservationRow, ...],
    targets: tuple[dict[str, object], ...],
    parameter_identity: str,
    *,
    lower: float,
    upper: float,
    origin: float,
    scale: float,
    starts: tuple[tuple[str, float], ...],
) -> dict[str, object]:
    payload = json.dumps(
        parameters.to_mapping(), separators=(",", ":"), sort_keys=True
    ).encode()
    parameter_path = tmp_path / "reactive-parameters.json"
    parameter_path.write_bytes(payload)
    reference = observation_problem.continuation_reference
    return {
        "schema_version": 2,
        "parameters": {
            "path": str(parameter_path),
            "sha256": f"sha256:{hashlib.sha256(payload).hexdigest()}",
        },
        "model": {"kind": "mixture"},
        "active": {
            "identities": [parameter_identity],
            "coordinates": {
                parameter_identity: {
                    "bounds": [lower, upper],
                    "origin": origin,
                    "scale": scale,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "homogeneous_reactive",
                "identity": observation_problem.identity,
                "phase": {
                    "identity": observation_problem.phase_identity,
                    "role": observation_problem.phase_role,
                    "model": "provider",
                    "admissible_packing_fraction_interval": list(
                        observation_problem.phase.admissible_packing_fraction_interval
                    ),
                },
                "reaction_system": _reaction_mapping(observation_problem.reaction_system),
                "continuation": {
                    "identity": observation_problem.continuation_identity,
                    "branch_policy": observation_problem.branch_policy,
                    "reference": None
                    if reference is None
                    else reference.to_mapping(),
                },
                "rows": [_row_mapping(row) for row in rows],
                "temperature_k": 313.15,
                "pressure_pa": 119_199.78746707023,
                "targets": list(targets),
                "multiplier": 1.0,
            }
        ],
        "solver": {
            "multistart": {
                "starts": [
                    {"identity": identity, "values": [value]} for identity, value in starts
                ],
                "affine_bound_margin": 1.0e-7,
            }
        },
    }


def _continuation_reference(
    result: equilibrium.HomogeneousReactiveObservationResult,
    *,
    maximum_log_composition_distance: float,
    maximum_log_volume_distance: float,
) -> equilibrium.HomogeneousContinuationReference:
    assert result.status == "evaluated", result.failure
    assert result.state is not None
    assert result.trial is not None
    state = result.state
    return equilibrium.HomogeneousContinuationReference(
        state.mole_fractions,
        state.volume_m3 / math.fsum(state.amounts_mol),
        maximum_log_composition_distance,
        maximum_log_volume_distance,
        state.phase_identity,
        state.phase_role,
        state.state_identity,
        state.phase_role_receipt.fingerprint,
        result.trial.parameter_fingerprint,
        state.reaction_fingerprint,
        state.topology_fingerprint,
        ((state.mole_fractions, state.volume_m3 / math.fsum(state.amounts_mol)),),
    )


def _observation_problem(
    phase: equilibrium.ProviderPhase,
    problem: equilibrium.ChemicalEquilibriumProblem,
    *,
    identity: str = "multineutral-fixed-tp",
    phase_identity: str = "multineutral-liquid",
    continuation_identity: str = "multineutral-liquid-branch",
    temperature_k: float = 313.15,
    pressure_pa: float = 119_199.78746707023,
) -> equilibrium.HomogeneousReactiveObservationProblem:
    unanchored = equilibrium.HomogeneousReactiveObservationProblem(
        identity=identity,
        phase_identity=phase_identity,
        phase_role="liquid",
        phase=phase,
        reaction_system=problem,
        continuation_identity=continuation_identity,
        branch_policy="local_certified_role_selected_state",
    )
    reference = equilibrium.certify_homogeneous_continuation_reference(
        unanchored,
        temperature_k * epcsaft.unit_registry.kelvin,
        pressure_pa * epcsaft.unit_registry.pascal,
        maximum_log_composition_distance=1.0,
        maximum_log_volume_distance=1.0,
    )
    return replace(
        unanchored,
        continuation_reference=reference,
    )


def test_non_evaluable_trial_vocabulary_is_closed() -> None:
    for code in epcsaft.NON_EVALUABLE_TRIAL_CODES:
        assert epcsaft.NonEvaluableTrial(code, "manufactured diagnostic").code == code
    with pytest.raises(ValueError, match="outside the public vocabulary"):
        epcsaft.NonEvaluableTrial("unregistered_failure", "manufactured diagnostic")  # type: ignore[arg-type]


@pytest.mark.slow
def test_source_reference_pressure_and_parameter_derivatives_match_resolves() -> None:
    """Falsify frozen-pressure or frozen-parameter source-reference transport."""

    parameters = epcsaft.Parameters.from_mapping(MULTISOLVENT_MULTIION)
    model = epcsaft.Mixture(parameters)
    identity = "pair/anion/solvent/k_ij"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    declaration = epcsaft.SourceReferenceDeclaration(
        reference_state_id="multisolvent-multiion-derivative-fixture",
        component_ids=parameters.component_ids,
        solvent_mole_fractions=(0.7, 0.3, 0.0, 0.0, 0.0, 0.0),
        ion_pairs=(
            epcsaft.SourceReferenceIonPair(("cation", "anion"), (1, 1)),
            epcsaft.SourceReferenceIonPair(("cation", "iodide-anion"), (1, 1)),
            epcsaft.SourceReferenceIonPair(("divalent-cation", "anion"), (1, 2)),
        ),
        phase="liquid",
        reference_convention="molality-infinite-dilution",
        activity_convention_id="molality-infinite-dilution-v1",
        standard_molality_mol_per_kg=1.0,
        reference_pressure_pa=100_000.0,
        required_derivatives=("pressure", "parameters"),
    )
    temperature = 313.15
    pressure = 1_000_000.0
    canonical = epcsaft.source_reference_transfer(
        model,
        replace(declaration, required_derivatives=()),
        T=temperature * epcsaft.unit_registry.kelvin,
        P=pressure * epcsaft.unit_registry.pascal,
    )
    exact = epcsaft.source_reference_transfer(
        model,
        declaration,
        T=temperature * epcsaft.unit_registry.kelvin,
        P=pressure * epcsaft.unit_registry.pascal,
        active_parameters=active,
    )
    assert exact.transfer_log_contractions == pytest.approx(
        canonical.transfer_log_contractions,
        rel=2.0e-13,
        abs=2.0e-13,
    )
    assert exact.log_fugacity_contractions == pytest.approx(
        canonical.log_fugacity_contractions,
        rel=2.0e-13,
        abs=2.0e-13,
    )

    step_pa = 10.0
    pressure_resolves = tuple(
        epcsaft.source_reference_transfer(
            model,
            declaration,
            T=temperature * epcsaft.unit_registry.kelvin,
            P=(pressure + direction * step_pa) * epcsaft.unit_registry.pascal,
            active_parameters=active,
        ).transfer_log_contractions
        for direction in (-1.0, 1.0)
    )
    pressure_difference = tuple(
        (above - below) / (2.0 * step_pa) for below, above in zip(*pressure_resolves, strict=True)
    )
    assert exact.pressure_derivatives_per_pa == pytest.approx(
        pressure_difference, rel=3.0e-6, abs=3.0e-12
    )

    step = 1.0e-5
    parameter_resolves = tuple(
        epcsaft.source_reference_transfer(
            model,
            declaration,
            T=temperature * epcsaft.unit_registry.kelvin,
            P=pressure * epcsaft.unit_registry.pascal,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters,
                (identity,),
                values=(active.values[0] + direction * step,),
            ),
        ).transfer_log_contractions
        for direction in (-1.0, 1.0)
    )
    parameter_difference = tuple(
        (above - below) / (2.0 * step) for below, above in zip(*parameter_resolves, strict=True)
    )
    assert tuple(row[0] for row in exact.parameter_derivatives) == pytest.approx(
        parameter_difference, rel=2.0e-6, abs=2.0e-9
    )
    assert exact.active_parameter_identities == (identity,)
    assert exact.derivative_availability == ("pressure", "parameters")


@pytest.mark.slow
def test_value_only_provider_potential_observation_has_zero_column_jacobian() -> None:
    """Falsify requiring a dummy active parameter for a Provider-owned value."""

    model, specification, source = _held_water_ionization_problem()
    reaction_problem = _typed_problem(
        copy.deepcopy(specification),
        source_standard_state=_typed_standard_state(copy.deepcopy(source)),
    )
    observation_problem = _observation_problem(
        equilibrium.ProviderPhase(model, (1.0e-6, 0.74)),
        reaction_problem,
        identity="water-ionization-value-only",
        phase_identity="water-ionization-liquid",
        continuation_identity="water-ionization-liquid-branch",
        temperature_k=float(specification["temperature_k"]),
        pressure_pa=float(specification["pressure_pa"]),
    )
    coefficients = (-2.0, 1.0, 1.0)
    row = equilibrium.HomogeneousReactiveObservationRow(
        "water-ionization-log-fugacity",
        "log_fugacity_contraction",
        coefficients,
    )
    descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        observation_problem,
        (row,),
    )
    result = equilibrium.evaluate_homogeneous_reactive_observations(
        observation_problem,
        float(specification["temperature_k"]) * epcsaft.unit_registry.kelvin,
        float(specification["pressure_pa"]) * epcsaft.unit_registry.pascal,
        descriptor,
    )

    assert result.status == "evaluated"
    assert result.trial is not None
    assert result.trial.identities == result.trial.values == result.trial.units == ()
    assert result.state is not None
    assert result.jacobian == ((),)
    fixed_state = model.state(
        T=result.state.temperature_k * epcsaft.unit_registry.kelvin,
        rho=(
            result.state.molar_density_mol_per_m3
            * epcsaft.unit_registry.mole
            / epcsaft.unit_registry.meter**3
        ),
        x=result.state.mole_fractions,
    )
    assert fixed_state.fugacity is not None
    expected = math.fsum(
        coefficient * math.log(value.to("pascal").magnitude)
        for coefficient, value in zip(
            coefficients,
            fixed_state.fugacity.value,
            strict=True,
        )
    )
    assert result.values == pytest.approx((expected,), abs=2.0e-10)

    assert result.equilibrium is not None
    assert observation_problem.continuation_reference is not None
    assert len(observation_problem.continuation_reference.anchor_basin_states) == 2
    ambiguous_problem = replace(
        observation_problem,
        continuation_reference=replace(
            observation_problem.continuation_reference,
            maximum_log_composition_distance=100.0,
            maximum_log_volume_distance=100.0,
        ),
    )
    ambiguous_descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        ambiguous_problem,
        (row,),
    )
    ambiguous = equilibrium.evaluate_homogeneous_reactive_observations(
        ambiguous_problem,
        float(specification["temperature_k"]) * epcsaft.unit_registry.kelvin,
        float(specification["pressure_pa"]) * epcsaft.unit_registry.pascal,
        ambiguous_descriptor,
    )
    assert ambiguous.status == "non_evaluable"
    assert ambiguous.failure is not None
    assert ambiguous.failure.code == "multiple_unresolved_branches"


@pytest.mark.slow
def test_homogeneous_observation_total_jacobian_matches_independent_resolves(
    tmp_path: Path,
) -> None:
    """Falsify partial-only Jacobians and individual-ion potential exposure."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    problem = _multineutral_problem()
    phase = equilibrium.ProviderPhase(model, (1.0e-6, 0.74))
    observation_problem = _observation_problem(phase, problem)
    identities = (
        "pair/anion/solvent/k_ij",
        "component/solvent/segment_diameter",
    )
    active = epcsaft.ActiveParameterSet(parameters, identities)
    rows = (
        equilibrium.HomogeneousReactiveObservationRow(
            "cation-mole-fraction",
            "mole_fraction",
            (0.0, 0.0, 1.0, 0.0),
            support="positive",
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "molar-density",
            "molar_density_mol_per_m3",
            support="positive",
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "reaction-log-fugacity",
            "log_fugacity_contraction",
            (-1.0, 0.0, 1.0, 1.0),
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "solvent-fugacity",
            "fugacity_pa",
            (1.0, 0.0, 0.0, 0.0),
            support="positive",
        ),
    )
    descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        observation_problem,
        rows,
        active_parameters=active,
    )
    evaluated = equilibrium.evaluate_homogeneous_reactive_observations(
        observation_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        descriptor,
        active_parameters=active,
    )
    assert evaluated.status == "evaluated"
    assert evaluated.physical_branch is not None
    assert evaluated.state is not None
    assert evaluated.trial is not None
    assert evaluated.equilibrium is not None
    assert not evaluated.equilibrium.diagnostics.derivative_audits_performed
    assert evaluated.descriptor.branch_policy == "local_certified_role_selected_state"
    assert evaluated.state.phase_role == "liquid"
    assert evaluated.state.phase_role_receipt.scope == (
        "local_reference_anchored_fixed_composition_role"
    )
    assert evaluated.state.phase_role_receipt.mechanical_class == "strict_stable"
    assert evaluated.state.phase_role_receipt.globality_claim == "not_established"
    assert evaluated.state.phase_role_receipt.phase_selection_claim == (
        "reference_anchored_locally_certified_declared_role"
    )
    assert (
        "exact_state_negative_pressure_derivative_log_volume"
        in evaluated.state.phase_role_receipt.role_evidence
    )
    assert (
        evaluated.state.phase_role_receipt.completeness_criterion
        == "one_certified_chemical_basin_inside_declared_reference_region;"
        "global_phase_stability_not_established"
    )
    assert evaluated.state.phase_role_receipt.provider_evaluation_count == 1
    assert (
        evaluated.state.phase_role_receipt.active_trial_fingerprint == evaluated.trial.fingerprint
    )
    assert evaluated.state.phase_role_receipt.fingerprint.startswith("sha256:")
    assert evaluated.physical_branch.phase_role_certificate_fingerprint == (
        evaluated.state.phase_role_receipt.fingerprint
    )
    assert evaluated.trial.values == active.values
    assert evaluated.trial.identities == active.identities
    assert evaluated.trial.active_parameter_fingerprint == active.fingerprint
    perturbed_active = epcsaft.ActiveParameterSet(
        parameters,
        identities,
        values=(active.values[0] + 0.01, active.values[1] + 0.001),
    )
    perturbed = equilibrium.evaluate_homogeneous_reactive_observations(
        observation_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        descriptor,
        active_parameters=perturbed_active,
    )
    assert perturbed.status == "evaluated"
    assert perturbed.trial is not None
    assert perturbed.trial.values == perturbed_active.values
    assert perturbed.trial.values != active.values
    assert perturbed.trial.active_parameter_fingerprint == perturbed_active.fingerprint
    repeated = equilibrium.evaluate_homogeneous_reactive_observations(
        observation_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        descriptor,
        active_parameters=active,
    )
    assert repeated.status == "evaluated"
    assert repeated.state is not None
    assert repeated.physical_branch is not None
    assert repeated.state.state_identity == evaluated.state.state_identity
    assert repeated.physical_branch.identity == evaluated.physical_branch.identity
    assert repeated.values == evaluated.values
    assert repeated.jacobian == evaluated.jacobian
    fixed_density_state = model.state(
        T=evaluated.state.temperature_k * epcsaft.unit_registry.kelvin,
        rho=(
            evaluated.state.molar_density_mol_per_m3
            * epcsaft.unit_registry.mole
            / epcsaft.unit_registry.meter**3
        ),
        x=evaluated.state.mole_fractions,
    )
    assert fixed_density_state.pressure.to("pascal").magnitude == pytest.approx(
        evaluated.state.eos_pressure_pa, rel=5.0e-12
    )
    assert math.fsum(evaluated.state.amounts_mol) / fixed_density_state.molar_density.to(
        "mole / meter**3"
    ).magnitude == pytest.approx(evaluated.state.volume_m3, rel=2.0e-12)
    assert fixed_density_state.fugacity is not None
    expected_log_fugacity = math.fsum(
        coefficient * math.log(value.to("pascal").magnitude)
        for coefficient, value in zip(
            rows[2].coefficients,
            fixed_density_state.fugacity.value,
            strict=True,
        )
    )
    assert evaluated.values[2] == pytest.approx(expected_log_fugacity, abs=2.0e-10)
    assert evaluated.values[3] == pytest.approx(
        fixed_density_state.fugacity.value[0].to("pascal").magnitude,
        rel=2.0e-10,
    )
    assert all(value is not None and math.isfinite(value) for value in evaluated.values)

    columns: list[tuple[float, ...]] = []
    for parameter, step in enumerate((1.0e-5, 1.0e-5)):
        resolves = tuple(
            equilibrium.evaluate_homogeneous_reactive_observations(
                observation_problem,
                313.15 * epcsaft.unit_registry.kelvin,
                119_199.78746707023 * epcsaft.unit_registry.pascal,
                descriptor,
                active_parameters=epcsaft.ActiveParameterSet(
                    parameters,
                    identities,
                    values=tuple(
                        value + direction * step if index == parameter else value
                        for index, value in enumerate(active.values)
                    ),
                ),
            )
            for direction in (-1.0, 1.0)
        )
        assert all(result.status == "evaluated" for result in resolves)
        columns.append(
            tuple(
                (float(above) - float(below)) / (2.0 * step)
                for below, above in zip(resolves[0].values, resolves[1].values, strict=True)
            )
        )
    differences = tuple(
        columns[column][row] for row in range(len(rows)) for column in range(len(identities))
    )
    exact = tuple(value for row in evaluated.jacobian if row is not None for value in row)
    assert exact == pytest.approx(differences, rel=3.0e-5, abs=3.0e-7)

    strict_continuation_problem = replace(
        observation_problem,
        continuation_reference=_continuation_reference(
            evaluated,
            maximum_log_composition_distance=1.0e-10,
            maximum_log_volume_distance=1.0e-10,
        ),
    )
    branch_active = epcsaft.ActiveParameterSet(parameters, (identities[0],))
    strict_descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        strict_continuation_problem,
        rows,
        active_parameters=branch_active,
    )
    anchor_before = equilibrium.evaluate_homogeneous_reactive_observations(
        strict_continuation_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        strict_descriptor,
        active_parameters=branch_active,
    )
    switched = equilibrium.evaluate_homogeneous_reactive_observations(
        strict_continuation_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        strict_descriptor,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters,
            (identities[0],),
            values=(branch_active.values[0] + 0.01,),
        ),
    )
    assert switched.status == "non_evaluable"
    assert switched.failure is not None
    assert switched.failure.code == "branch_loss"
    anchor_after = equilibrium.evaluate_homogeneous_reactive_observations(
        strict_continuation_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        strict_descriptor,
        active_parameters=branch_active,
    )
    assert anchor_before.status == anchor_after.status == "evaluated"
    assert anchor_before.values == anchor_after.values
    assert anchor_before.jacobian == anchor_after.jacobian
    assert anchor_before.physical_branch == anchor_after.physical_branch
    continuation_target = _target_mapping(
        "synthetic-row/branch-tracer",
        rows[0].identity,
        float(evaluated.values[0]),
        rows[0].unit,
        "scaled_difference",
        float(evaluated.values[0]),
    )
    mixed = regression.fit(
        _reactive_fit_mapping(
            tmp_path,
            parameters,
            strict_continuation_problem,
            rows,
            (continuation_target,),
            identities[0],
            lower=branch_active.values[0] - 0.02,
            upper=branch_active.values[0] + 0.02,
            origin=branch_active.values[0],
            scale=0.01,
            starts=(
                ("A", branch_active.values[0]),
                ("B", branch_active.values[0] + 0.01),
            ),
        ),
        base_path=tmp_path,
    )
    assert tuple(start.start_identity for start in mixed.starts) == ("A", "B")
    assert mixed.starts[0].usable
    assert not mixed.starts[1].usable
    assert mixed.starts[1].trial_failure_code == "branch_loss"
    assert mixed.starts[1].batch_outcomes[0].status == "rejected"
    assert mixed.accounting.attempted_starts == 2
    assert mixed.accounting.accepted_starts == 1
    assert mixed.accounting.failed_starts == 1
    assert mixed.accounting.evaluated_row_trials == 1
    assert mixed.accounting.rejected_row_trials == 1

    selected_basin = evaluated.equilibrium.diagnostics.search.selected_basin_ordinal
    selected_starts = {
        attempt.start_identity
        for attempt in evaluated.equilibrium.diagnostics.search.attempts
        if attempt.basin_ordinal == selected_basin
    }
    assert selected_basin is not None
    assert len(selected_starts) >= 2

    mislabeled_problem = replace(
        observation_problem,
        identity="multineutral-fixed-tp-mislabeled-vapor",
        phase_role="vapor",
    )
    mislabeled_descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        mislabeled_problem,
        rows,
        active_parameters=active,
    )
    mislabeled = equilibrium.evaluate_homogeneous_reactive_observations(
        mislabeled_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        119_199.78746707023 * epcsaft.unit_registry.pascal,
        mislabeled_descriptor,
        active_parameters=active,
    )
    assert mislabeled.status == "non_evaluable"
    assert mislabeled.failure is not None
    assert mislabeled.failure.code == "physical_rejection"

    with pytest.raises(ValueError, match="charge neutral"):
        equilibrium.homogeneous_reactive_observation_descriptor(
            observation_problem,
            (
                equilibrium.HomogeneousReactiveObservationRow(
                    "individual-ion-fugacity",
                    "log_fugacity_contraction",
                    (0.0, 0.0, 1.0, 0.0),
                ),
            ),
            active_parameters=active,
        )


@pytest.mark.slow
def test_physical_branch_receipt_binds_problem_and_exact_state() -> None:
    """Falsify a branch receipt that omits the exact problem or certified state."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    phase = equilibrium.ProviderPhase(model, (1.0e-6, 0.74))
    first_reaction = _multineutral_problem()
    second_reaction = replace(
        first_reaction,
        feed_amounts_mol=(0.75, 0.2, 0.15, 0.15),
    )
    identity = "pair/anion/solvent/k_ij"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    row = equilibrium.HomogeneousReactiveObservationRow(
        "cation-mole-fraction",
        "mole_fraction",
        (0.0, 0.0, 1.0, 0.0),
        support="positive",
    )
    results = []
    for problem_identity, reaction in (
        ("feed-start-a", first_reaction),
        ("feed-start-b", second_reaction),
    ):
        observation_problem = _observation_problem(
            phase,
            reaction,
            identity=problem_identity,
        )
        descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
            observation_problem,
            (row,),
            active_parameters=active,
        )
        results.append(
            equilibrium.evaluate_homogeneous_reactive_observations(
                observation_problem,
                313.15 * epcsaft.unit_registry.kelvin,
                119_199.78746707023 * epcsaft.unit_registry.pascal,
                descriptor,
                active_parameters=active,
            )
        )
    assert all(result.status == "evaluated" for result in results)
    assert results[0].physical_branch is not None
    assert results[1].physical_branch is not None
    assert results[0].physical_branch.identity != results[1].physical_branch.identity
    assert (
        results[0].physical_branch.problem_fingerprint
        != results[1].physical_branch.problem_fingerprint
    )
    assert results[0].physical_branch.state_identity != results[1].physical_branch.state_identity
    assert results[0].physical_branch.amounts_mol == pytest.approx(
        results[1].physical_branch.amounts_mol,
        rel=2.0e-9,
        abs=2.0e-11,
    )
    assert results[0].physical_branch.volume_m3 == pytest.approx(
        results[1].physical_branch.volume_m3,
        rel=2.0e-9,
    )


@pytest.mark.slow
def test_source_reference_observation_is_coordinate_invariant_and_total() -> None:
    """Falsify order-, orientation-, or partial-only source-reference derivatives."""

    def evaluate_variant(
        variant: str,
    ) -> tuple[
        equilibrium.HomogeneousReactiveObservationResult,
        epcsaft.Parameters,
        str,
    ]:
        base_model, base_spec, base_source = _held_water_ionization_problem()
        spec = copy.deepcopy(base_spec)
        source = copy.deepcopy(base_source)
        components = tuple(base_model.component_ids)
        coefficients = (-2.0, 1.0, 1.0)
        if variant == "component_reorder":
            permutation = (2, 0, 1)
            components = tuple(components[index] for index in permutation)
            for field in (
                "species_ids",
                "charges",
                "molar_masses_kg_per_mol",
                "feed_amounts",
            ):
                values = spec[field]
                spec[field] = tuple(values[index] for index in permutation)
            spec["balance_matrix"] = tuple(
                tuple(row[index] for index in permutation) for row in spec["balance_matrix"]
            )
            spec["reaction_matrix"] = tuple(
                tuple(row[index] for index in permutation) for row in spec["reaction_matrix"]
            )
            for field in (
                "log_activity_scale_factors",
                "source_reference_component_ids",
                "source_reference_solvent_composition",
            ):
                values = source[field]
                source[field] = tuple(values[index] for index in permutation)
            coefficients = tuple(coefficients[index] for index in permutation)
        elif variant == "reaction_reversal":
            spec["reaction_matrix"] = tuple(
                tuple(-value for value in row) for row in spec["reaction_matrix"]
            )
            spec["ln_k"] = tuple(-value for value in spec["ln_k"])
        elif variant != "baseline":
            raise AssertionError(f"unknown source-reference variant {variant}")

        parameters = _held_parameters(components)
        model = epcsaft.Mixture(parameters)
        spec["provider_fingerprint"] = model.parameter_fingerprint
        reaction_problem = _typed_problem(
            spec,
            source_standard_state=_typed_standard_state(source),
        )
        observation_problem = _observation_problem(
            equilibrium.ProviderPhase(model, (1.0e-6, 0.74)),
            reaction_problem,
            identity=f"water-ionization-{variant}",
            phase_identity="water-ionization-liquid",
            continuation_identity="water-ionization-liquid-branch",
            temperature_k=float(spec["temperature_k"]),
            pressure_pa=float(spec["pressure_pa"]),
        )
        parameter_identity = "pair/hydronium-cation/water/k_ij"
        active = epcsaft.ActiveParameterSet(parameters, (parameter_identity,))
        descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
            observation_problem,
            (
                equilibrium.HomogeneousReactiveObservationRow(
                    "water-ionization-log-activity",
                    "log_activity_contraction",
                    coefficients,
                ),
            ),
            active_parameters=active,
        )
        result = equilibrium.evaluate_homogeneous_reactive_observations(
            observation_problem,
            float(spec["temperature_k"]) * epcsaft.unit_registry.kelvin,
            float(spec["pressure_pa"]) * epcsaft.unit_registry.pascal,
            descriptor,
            active_parameters=active,
        )
        return result, parameters, parameter_identity

    baseline, parameters, parameter_identity = evaluate_variant("baseline")
    reordered, _, _ = evaluate_variant("component_reorder")
    reversed_reaction, _, _ = evaluate_variant("reaction_reversal")
    assert baseline.status == reordered.status == reversed_reaction.status == "evaluated"
    assert float(baseline.values[0]) == pytest.approx(
        float(reordered.values[0]), rel=2.0e-9, abs=2.0e-11
    )
    assert float(baseline.values[0]) == pytest.approx(
        float(reversed_reaction.values[0]), rel=2.0e-9, abs=2.0e-11
    )
    assert baseline.jacobian[0] is not None
    assert reordered.jacobian[0] is not None
    assert reversed_reaction.jacobian[0] is not None
    assert baseline.jacobian[0][0] == pytest.approx(
        reordered.jacobian[0][0], rel=2.0e-6, abs=2.0e-8
    )
    assert baseline.jacobian[0][0] == pytest.approx(
        reversed_reaction.jacobian[0][0], rel=2.0e-6, abs=2.0e-8
    )
    assert (
        baseline.descriptor.reaction_fingerprint
        != reversed_reaction.descriptor.reaction_fingerprint
    )

    step = 1.0e-5
    resolved_values = []
    base_model, base_spec, base_source = _held_water_ionization_problem()
    reaction_problem = _typed_problem(
        copy.deepcopy(base_spec),
        source_standard_state=_typed_standard_state(copy.deepcopy(base_source)),
    )
    observation_problem = _observation_problem(
        equilibrium.ProviderPhase(base_model, (1.0e-6, 0.74)),
        reaction_problem,
        identity="water-ionization-baseline",
        phase_identity="water-ionization-liquid",
        continuation_identity="water-ionization-liquid-branch",
        temperature_k=float(base_spec["temperature_k"]),
        pressure_pa=float(base_spec["pressure_pa"]),
    )
    active = epcsaft.ActiveParameterSet(parameters, (parameter_identity,))
    row = equilibrium.HomogeneousReactiveObservationRow(
        "water-ionization-log-activity",
        "log_activity_contraction",
        (-2.0, 1.0, 1.0),
    )
    descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        observation_problem,
        (row,),
        active_parameters=active,
    )
    for direction in (-1.0, 1.0):
        resolved = equilibrium.evaluate_homogeneous_reactive_observations(
            observation_problem,
            float(base_spec["temperature_k"]) * epcsaft.unit_registry.kelvin,
            float(base_spec["pressure_pa"]) * epcsaft.unit_registry.pascal,
            descriptor,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters,
                (parameter_identity,),
                values=(active.values[0] + direction * step,),
            ),
        )
        assert resolved.status == "evaluated"
        resolved_values.append(float(resolved.values[0]))
    difference = (resolved_values[1] - resolved_values[0]) / (2.0 * step)
    assert baseline.jacobian[0][0] == pytest.approx(difference, rel=3.0e-5, abs=3.0e-7)


@pytest.mark.slow
def test_equilibrium_regression_batches_exact_rows_and_reports_accounting(
    tmp_path: Path,
) -> None:
    """Falsify row loss, residual penalties, and opaque identifiability diagnostics."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    phase = equilibrium.ProviderPhase(model, (1.0e-6, 0.74))
    problem = _multineutral_problem()
    observation_problem = _observation_problem(phase, problem)
    identity = "pair/anion/solvent/k_ij"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    rows = (
        equilibrium.HomogeneousReactiveObservationRow(
            "cation-mole-fraction", "mole_fraction", (0.0, 0.0, 1.0, 0.0), support="positive"
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "total-ion-mole-fraction",
            "mole_fraction",
            (0.0, 0.0, 1.0, 1.0),
            support="positive",
            aggregate_identity="manufactured-total-ions",
            covariance_identity="manufactured-prepared-stock",
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "molar-density", "molar_density_mol_per_m3", support="positive"
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "cation-detection-limit",
            "mole_fraction",
            (0.0, 0.0, 1.0, 0.0),
            support="left_censored",
            censor_limit=0.0018,
        ),
        equilibrium.HomogeneousReactiveObservationRow(
            "cation-lower-limit",
            "mole_fraction",
            (0.0, 0.0, 1.0, 0.0),
            support="right_censored",
            censor_limit=0.0,
        ),
        equilibrium.HomogeneousReactiveObservationRow("state-pressure", "pressure_pa"),
    )
    descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        observation_problem,
        rows,
        active_parameters=active,
    )
    temperature = 313.15 * epcsaft.unit_registry.kelvin
    pressure = 119_199.78746707023 * epcsaft.unit_registry.pascal
    truth = equilibrium.evaluate_homogeneous_reactive_observations(
        observation_problem,
        temperature,
        pressure,
        descriptor,
        active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(0.08,)),
    )
    assert truth.status == "evaluated"
    assert truth.state is not None
    observation_problem = replace(
        observation_problem,
        continuation_reference=_continuation_reference(
            truth,
            maximum_log_composition_distance=1.0,
            maximum_log_volume_distance=0.5,
        ),
    )
    descriptor = equilibrium.homogeneous_reactive_observation_descriptor(
        observation_problem,
        rows,
        active_parameters=active,
    )
    truth = equilibrium.evaluate_homogeneous_reactive_observations(
        observation_problem,
        temperature,
        pressure,
        descriptor,
        active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(0.08,)),
    )
    assert truth.status == "evaluated"
    cation, total_ions, density, censored_cation, right_censored_cation, state_pressure = (
        float(value) for value in truth.values
    )
    assert censored_cation < 0.0018
    assert right_censored_cation > 0.0
    targets = (
        _target_mapping(
            "synthetic-row/cation-lower-limit",
            "cation-lower-limit",
            0.0,
            "dimensionless",
            "right_censored",
            0.0001,
        ),
        _target_mapping(
            "synthetic-row/cation",
            "cation-mole-fraction",
            cation,
            "dimensionless",
            "log_ratio",
            cation,
        ),
        _target_mapping(
            "synthetic-row/total-ions",
            "total-ion-mole-fraction",
            total_ions,
            "dimensionless",
            "log_ratio",
            total_ions,
            aggregate_identity="manufactured-total-ions",
            covariance_identity="manufactured-prepared-stock",
        ),
        _target_mapping(
            "synthetic-row/density",
            "molar-density",
            density,
            "mole / meter**3",
            "scaled_difference",
            density,
        ),
        _target_mapping(
            "synthetic-row/cation-detection-limit",
            "cation-detection-limit",
            0.0018,
            "dimensionless",
            "left_censored",
            0.0001,
        ),
        _target_mapping(
            "synthetic-row/pressure-diagnostic",
            "state-pressure",
            state_pressure,
            "pascal",
            "scaled_difference",
            state_pressure,
            included=False,
            drop_reason="model-derived pressure is not a fit target",
        ),
    )
    fit = regression.fit(
        _reactive_fit_mapping(
            tmp_path,
            parameters,
            observation_problem,
            rows,
            targets,
            identity,
            lower=-0.05,
            upper=0.2,
            origin=0.0,
            scale=0.1,
            starts=(("origin", 0.0), ("confirmation", 0.05)),
        ),
        base_path=tmp_path,
    )
    accepted = fit.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((0.08,), abs=2.0e-8)
    assert accepted.predictions == pytest.approx(
        (right_censored_cation, cation, total_ions, density, censored_cation), rel=2.0e-8
    )
    assert accepted.rank == 1
    assert accepted.free_subspace_rank == 1
    assert accepted.information_row_count == 3
    assert accepted.degrees_of_freedom == 2
    assert accepted.covariance_status == "available"
    assert len(fit.physical_covariance) == 1
    assert fit.accounting.input_rows == 6
    assert fit.accounting.evaluated_row_trials == 10
    assert fit.accounting.rejected_row_trials == 0
    assert fit.accounting.failed_row_trials == 0
    assert fit.accounting.dropped_rows == 1
    assert fit.accounting.attempted_starts == 2
    assert fit.accounting.accepted_starts == 2
    candidate = fit.export_candidate_parameters()
    assert candidate.fingerprint != parameters.fingerprint
    assert fit.evaluated_row_table_hash.startswith("sha256:")
    assert len(fit.evaluated_rows) == 5
    fit_mapping = fit.to_mapping()
    accepted_receipts = fit_mapping["starts"][accepted.index]["owner_receipts"]
    assert isinstance(accepted_receipts, list)
    accepted_receipt = json.loads(accepted_receipts[0]["canonical_bytes_utf8"])
    assert accepted_receipts[0]["digest"].startswith("sha256:")
    accepted_observations = accepted_receipt["observations"]
    aggregate_receipt = next(
        observation
        for observation in accepted_observations
        if observation["prediction_identity"] == "total-ion-mole-fraction"
    )
    assert aggregate_receipt["aggregate_identity"] == "manufactured-total-ions"
    assert aggregate_receipt["covariance_identity"] == "manufactured-prepared-stock"
    assert accepted.optimizer_basin is not None
    assert accepted_receipt["physical_branch"]["identity"].startswith("sha256:")
    search_certificate = accepted_receipt["certified_state"]["chemical_search_certificate"]
    assert search_certificate["basin_count"] == 1
    continuation_certificate = accepted_receipt["certified_state"][
        "continuation_certificate"
    ]
    assert continuation_certificate["in_region_basin_count"] == 1
    assert continuation_certificate["selected_basin_ordinal"] == search_certificate[
        "selected_basin_ordinal"
    ]
    assert search_certificate["budget_prefix_selection_history"]

    tampered_receipt = json.loads(accepted.owner_receipts[0].canonical_bytes_utf8)
    tampered_receipt["certified_state"]["chemical_search_certificate"][
        "selected_basin_representative_attempt_ordinal"
    ] += 1
    tampered_json = json.dumps(
        tampered_receipt,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    with pytest.raises(ValueError, match="receipt|certificate"):
        equilibrium.validate_homogeneous_reactive_owner_receipt(
            tampered_json,
            descriptor=descriptor,
            batch_identity=observation_problem.identity,
            requested_temperature=temperature,
            requested_pressure=pressure,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (identity,), values=accepted.final_physical
            ),
            observations=tuple(
                equilibrium.HomogeneousReactiveOwnerObservation(
                    str(target["identity"]),
                    str(target["prediction_identity"]),
                    cast(str | None, target["aggregate_identity"]),
                    cast(str | None, target["covariance_identity"]),
                )
                for target in targets
                if target["included"]
            ),
        )

    bounded_fit = regression.fit(
        _reactive_fit_mapping(
            tmp_path,
            parameters,
            observation_problem,
            rows,
            targets,
            identity,
            lower=-0.05,
            upper=0.05,
            origin=0.0,
            scale=0.1,
            starts=(("origin", 0.0),),
        ),
        base_path=tmp_path,
    )
    bounded = bounded_fit.best_usable_start
    assert bounded is not None
    assert bounded.final_physical == pytest.approx((0.05,), abs=2.0e-9)
    assert bounded.rank == 1
    assert bounded.active_bounds == (0,)
    assert bounded.covariance_status == "unavailable_active_bounds"
    assert bounded.affine_covariance == ()
