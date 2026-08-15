from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium, regression
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE

from .test_chemical_provider import _multineutral_problem
from .test_homogeneous_reactive_observations import (
    _observation_problem,
    _reaction_mapping,
    _row_mapping,
    _target_mapping,
)


def _fixture() -> tuple[
    epcsaft.Parameters,
    equilibrium.ReactiveBubbleVLEProblem,
    tuple[equilibrium.ReactiveBubbleObservationRow, ...],
]:
    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    liquid = _observation_problem(
        equilibrium.ProviderPhase(model, (1.0e-6, 0.70)),
        _multineutral_problem(),
    )
    problem = equilibrium.ReactiveBubbleVLEProblem(
        identity="source-neutral-reactive-bubble",
        liquid_problem=liquid,
        vapor_phase_identity="source-neutral-vapor",
        vapor_component_ids=("solvent", "cosolvent"),
        vapor_model=equilibrium.ProviderNonidealVapor("installed-provider-eos"),
        pressure_interval_pa=(2.5e6, 4.0e6),
        pressure_starts_pa=(2.8e6, 3.2e6, 3.6e6),
        continuation_identity="source-neutral-reactive-bubble-branch",
    )
    rows = (
        equilibrium.ReactiveBubbleObservationRow("bubble-pressure", "bubble_pressure_pa"),
        equilibrium.ReactiveBubbleObservationRow(
            "first-vapor-fraction", "vapor_mole_fraction", (1.0, 0.0)
        ),
        equilibrium.ReactiveBubbleObservationRow(
            "first-partial-pressure", "vapor_partial_pressure_pa", (1.0, 0.0)
        ),
        equilibrium.ReactiveBubbleObservationRow(
            "first-liquid-fraction", "liquid_mole_fraction", (1.0, 0.0, 0.0, 0.0)
        ),
    )
    return parameters, problem, rows


def _multireaction_fixture() -> tuple[
    epcsaft.Parameters,
    equilibrium.ReactiveBubbleVLEProblem,
    tuple[equilibrium.ReactiveBubbleObservationRow, ...],
]:
    """Extend the generic fixture with one neutral isomerization reaction."""

    mapping = deepcopy(MULTINEUTRAL_REACTIVE)
    mapping["document_id"] = "multireaction-reactive-bubble-fixture"
    components = {item["component_id"]: item for item in mapping["components"]}
    isomer = deepcopy(components["cosolvent"])
    isomer["component_id"] = "isomer"
    for coefficient in isomer["coefficients"]:
        coefficient["identity"] = coefficient["identity"].replace("/cosolvent/", "/isomer/")
    provenance = deepcopy(components["solvent"]["fixed"]["molar_mass"]["provenance"])
    for other in ("solvent", "cosolvent", "cation", "anion"):
        left, right = sorted(("isomer", other))
        mapping["pairs"].append(
            {
                "component_id_a": left,
                "component_id_b": right,
                "coefficients": [
                    {
                        "identity": f"pair/{left}/{right}/k_ij",
                        "family": "k_ij",
                        "value": {"magnitude": 0.0, "unit": "dimensionless"},
                        "provenance": deepcopy(provenance),
                    }
                ],
            }
        )
    for correlation in tuple(mapping["correlations"]):
        if correlation["component_id"] == "cosolvent":
            isomer_correlation = deepcopy(correlation)
            isomer_correlation["component_id"] = "isomer"
            isomer_correlation["correlation_id"] = isomer_correlation["correlation_id"].replace(
                "/cosolvent/", "/isomer/"
            )
            isomer_correlation["constant"]["identity"] = isomer_correlation["constant"][
                "identity"
            ].replace("/cosolvent/", "/isomer/")
            mapping["correlations"].append(isomer_correlation)
    components["isomer"] = isomer
    component_ids = ("solvent", "cosolvent", "isomer", "cation", "anion")
    mapping["components"] = [components[identity] for identity in component_ids]
    parameters = epcsaft.Parameters.from_mapping(mapping)
    reference_id = "provider-helmholtz-coordinate-basis"
    reaction_system = equilibrium.ChemicalEquilibriumProblem(
        species_ids=component_ids,
        charges=(0, 0, 0, 1, -1),
        molar_masses_kg_per_mol=(0.018, 0.030, 0.030, 0.009, 0.009),
        balance_matrix=(
            (1.0, 0.0, 0.0, 1.0, 0.0),
            (1.0, 0.0, 0.0, 0.0, 1.0),
            (0.0, 1.0, 1.0, 0.0, 0.0),
        ),
        conserved_totals=(0.9, 0.9, 0.2),
        reaction_matrix=((-1.0, 0.0, 0.0, 1.0, 1.0), (0.0, -1.0, 1.0, 0.0, 0.0)),
        feed_amounts_mol=(0.8, 0.15, 0.05, 0.1, 0.1),
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                -269.6037501925047,
                "synthetic-multireaction-source",
                reference_id,
                "products_positive",
                "already-provider-basis",
                True,
            ),
            equilibrium.ChemicalEquilibriumConstant(
                0.0,
                "synthetic-multireaction-source",
                reference_id,
                "products_positive",
                "already-provider-basis",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-30,
    )
    liquid = _observation_problem(
        equilibrium.ProviderPhase(epcsaft.Mixture(parameters), (1.0e-6, 0.70)),
        reaction_system,
        identity="multireaction-fixed-tp",
        phase_identity="multireaction-liquid",
        continuation_identity="multireaction-liquid-branch",
    )
    problem = equilibrium.ReactiveBubbleVLEProblem(
        "multireaction-reactive-bubble",
        liquid,
        "multireaction-vapor",
        ("solvent", "cosolvent", "isomer"),
        equilibrium.ProviderNonidealVapor("installed-provider-eos"),
        (2.0e6, 5.0e6),
        (2.5e6, 3.2e6, 4.0e6),
        "multireaction-reactive-bubble-branch",
    )
    rows = (
        equilibrium.ReactiveBubbleObservationRow("bubble-pressure", "bubble_pressure_pa"),
        equilibrium.ReactiveBubbleObservationRow(
            "first-vapor-fraction", "vapor_mole_fraction", (1.0, 0.0, 0.0)
        ),
        equilibrium.ReactiveBubbleObservationRow(
            "cation-liquid-fraction", "liquid_mole_fraction", (0.0, 0.0, 0.0, 1.0, 0.0)
        ),
    )
    return parameters, problem, rows


def _reordered_reaction_system(
    problem: equilibrium.ChemicalEquilibriumProblem,
    permutation: tuple[int, ...],
) -> equilibrium.ChemicalEquilibriumProblem:
    return replace(
        problem,
        species_ids=tuple(problem.species_ids[index] for index in permutation),
        charges=tuple(problem.charges[index] for index in permutation),
        molar_masses_kg_per_mol=tuple(
            problem.molar_masses_kg_per_mol[index] for index in permutation
        ),
        balance_matrix=tuple(
            tuple(row[index] for index in permutation) for row in problem.balance_matrix
        ),
        reaction_matrix=tuple(
            tuple(row[index] for index in permutation) for row in problem.reaction_matrix
        ),
        feed_amounts_mol=tuple(problem.feed_amounts_mol[index] for index in permutation),
    )


@pytest.mark.slow
def test_nonideal_reactive_bubble_total_jacobian_and_measurement_semantics() -> None:
    """Falsify frozen speciation, ideal-vapor fallback, or p_i=f_i measurement drift."""

    parameters, problem, rows = _fixture()
    identity = "component/solvent/dispersion_energy_over_k"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    descriptor = equilibrium.reactive_bubble_observation_descriptor(
        problem, rows, active_parameters=active
    )
    result = equilibrium.evaluate_reactive_bubble_observations(
        problem,
        313.15 * epcsaft.unit_registry.kelvin,
        descriptor,
        active_parameters=active,
    )
    assert result.status == "evaluated", result.failure
    assert result.pressure_pa is not None
    assert result.physical_branch is not None
    assert result.certification is not None
    assert result.certification.status == "certified_local_equilibrium"
    assert result.certification.second_order_status == "passed"
    assert result.certification.globality_claim == "not_established"
    assert result.physical_branch.globality_claim == "not_established"
    assert result.search.distinct_branch_count == 1
    assert result.search.accepted_candidate_count == len(problem.pressure_starts_pa)
    assert tuple(attempt.start_pressure_pa for attempt in result.search.attempts) == (
        problem.pressure_starts_pa
    )
    assert result.phases[1].support == "declared_neutral_components"
    assert result.phases[1].supported_component_ids == problem.vapor_component_ids
    assert result.liquid_state is not None
    assert result.phases[0].mole_fractions == pytest.approx(
        result.liquid_state.mole_fractions, rel=0.0, abs=2.0e-13
    )
    assert result.phases[0].molar_volume_m3_per_mol == pytest.approx(
        result.liquid_state.volume_m3 / sum(result.liquid_state.amounts_mol),
        rel=2.0e-9,
    )
    pressure, vapor_fraction, partial_pressure, _ = result.values
    assert pressure is not None and vapor_fraction is not None and partial_pressure is not None
    assert partial_pressure == pytest.approx(pressure * vapor_fraction, rel=2.0e-13)
    vapor_phase = result.phases[1]
    vapor_state = epcsaft.Mixture(parameters.select(problem.vapor_component_ids))._state(
        T=313.15 * epcsaft.unit_registry.kelvin,
        rho=(
            vapor_phase.molar_density_mol_m3
            * epcsaft.unit_registry.mole
            / epcsaft.unit_registry.meter**3
        ),
        x=vapor_phase.mole_fractions,
    )
    assert vapor_state.fugacity is not None
    equilibrium_fugacity = float(vapor_state.fugacity.value[0].to("pascal").magnitude)
    assert equilibrium_fugacity != pytest.approx(partial_pressure, rel=1.0e-5)
    assert result.phases[0].chemical_potential_over_rt[:2] == pytest.approx(
        result.phases[1].chemical_potential_over_rt,
        rel=2.0e-8,
        abs=2.0e-8,
    )

    value = active.values[0]
    # The bubble root is certified to an absolute reduced-potential residual.
    # A 1e-2 K parameter step keeps root stopping error below the derivative
    # comparison while remaining a 3e-5 relative perturbation here.
    step = 1.0e-2
    resolved = []
    for sign in (-1.0, 1.0):
        trial = epcsaft.ActiveParameterSet(parameters, (identity,), values=(value + sign * step,))
        trial_descriptor = equilibrium.reactive_bubble_observation_descriptor(
            problem, rows, active_parameters=trial
        )
        trial_result = equilibrium.evaluate_reactive_bubble_observations(
            problem,
            313.15 * epcsaft.unit_registry.kelvin,
            trial_descriptor,
            active_parameters=trial,
        )
        assert trial_result.status == "evaluated", trial_result.failure
        resolved.append(trial_result)
    for row in range(len(rows)):
        lower = resolved[0].values[row]
        upper = resolved[1].values[row]
        analytic = result.jacobian[row]
        assert lower is not None and upper is not None and analytic is not None
        central = (upper - lower) / (2.0 * step)
        scale = max(1.0, abs(central), abs(analytic[0]))
        assert abs(analytic[0] - central) / scale < 2.0e-7


@pytest.mark.slow
def test_reactive_bubble_start_order_is_pure_and_receipt_is_self_validating() -> None:
    """Falsify mutable continuation, search-ordinal branch identity, or receipt drift."""

    parameters, problem, rows = _fixture()
    identity = "component/solvent/dispersion_energy_over_k"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))

    def solve(candidate: equilibrium.ReactiveBubbleVLEProblem):
        descriptor = equilibrium.reactive_bubble_observation_descriptor(
            candidate, rows, active_parameters=active
        )
        return descriptor, equilibrium.evaluate_reactive_bubble_observations(
            candidate,
            313.15 * epcsaft.unit_registry.kelvin,
            descriptor,
            active_parameters=active,
        )

    first_descriptor, first = solve(problem)
    reversed_problem = replace(
        problem, pressure_starts_pa=tuple(reversed(problem.pressure_starts_pa))
    )
    second_descriptor, second = solve(reversed_problem)
    _third_descriptor, third = solve(problem)
    for result in (first, second, third):
        assert result.status == "evaluated", result.failure
        assert result.physical_branch is not None
    assert first.values == pytest.approx(third.values, rel=2.0e-12, abs=2.0e-12)
    for first_row, third_row in zip(first.jacobian, third.jacobian, strict=True):
        assert first_row == pytest.approx(third_row, rel=2.0e-12, abs=2.0e-12)
    assert first.physical_branch == third.physical_branch
    assert first.physical_branch.branch_family_identity != ""
    # Search order is accounting, not physical identity.
    assert first_descriptor.fingerprint == second_descriptor.fingerprint
    assert first.physical_branch.branch_family_identity == (
        second.physical_branch.branch_family_identity
    )
    assert first.physical_branch.state_identity == second.physical_branch.state_identity

    receipt = equilibrium.serialize_reactive_bubble_owner_receipt(
        first, batch_identity="source-neutral-batch"
    )
    digest = equilibrium.validate_reactive_bubble_owner_receipt(
        receipt,
        descriptor=first_descriptor,
        batch_identity="source-neutral-batch",
        requested_temperature=313.15 * epcsaft.unit_registry.kelvin,
        active_parameters=active,
    )
    assert digest.startswith("sha256:")
    with pytest.raises(ValueError, match="physical branch is inconsistent"):
        equilibrium.validate_reactive_bubble_owner_receipt(
            receipt.replace('"pressure_pa":3348', '"pressure_pa":3349', 1),
            descriptor=first_descriptor,
            batch_identity="source-neutral-batch",
            requested_temperature=313.15 * epcsaft.unit_registry.kelvin,
            active_parameters=active,
        )


@pytest.mark.slow
def test_reactive_bubble_rejects_charged_vapor_and_evaluates_explicit_ideal_model() -> None:
    """Falsify component-named support, dummy coordinates, or silent vapor fallback."""

    _, problem, rows = _fixture()
    with pytest.raises(ValueError, match="topology"):
        replace(problem, vapor_component_ids=("solvent", "cation"))
    comparison = replace(
        problem,
        vapor_model=equilibrium.IdealGasVaporComparison("explicit-comparison"),
    )
    descriptor = equilibrium.reactive_bubble_observation_descriptor(comparison, rows)
    result = equilibrium.evaluate_reactive_bubble_observations(
        comparison, 313.15 * epcsaft.unit_registry.kelvin, descriptor
    )
    assert result.status == "evaluated", result.failure
    assert result.jacobian == ((),) * len(rows)
    assert result.phases[1].mechanical_class == "strict_stable_ideal_gas"
    assert result.phases[0].chemical_potential_over_rt[:2] == pytest.approx(
        result.phases[1].chemical_potential_over_rt, rel=2.0e-12, abs=2.0e-12
    )


@pytest.mark.slow
def test_reactive_bubble_reference_fails_closed_on_branch_loss() -> None:
    """Falsify permissive continuation outside the immutable coupled neighborhood."""

    _, problem, rows = _fixture()
    descriptor = equilibrium.reactive_bubble_observation_descriptor(problem, rows)
    result = equilibrium.evaluate_reactive_bubble_observations(
        problem, 313.15 * epcsaft.unit_registry.kelvin, descriptor
    )
    assert result.status == "evaluated", result.failure
    assert result.physical_branch is not None
    assert result.pressure_pa is not None
    liquid, vapor = result.phases
    reference = equilibrium.ReactiveBubbleContinuationReference(
        pressure_pa=0.5 * result.pressure_pa,
        vapor_mole_fractions=vapor.mole_fractions,
        liquid_molar_volume_m3_per_mol=liquid.molar_volume_m3_per_mol,
        vapor_molar_volume_m3_per_mol=vapor.molar_volume_m3_per_mol,
        maximum_log_pressure_distance=1.0e-3,
        maximum_log_composition_distance=1.0e-3,
        maximum_log_volume_distance=1.0e-3,
        anchor_state_identity=result.physical_branch.state_identity,
        anchor_problem_fingerprint=descriptor.problem_fingerprint,
        anchor_topology_fingerprint=descriptor.topology_fingerprint,
        anchor_parameter_fingerprint=result.trial.parameter_fingerprint,
    )
    rejected_problem = replace(problem, continuation_reference=reference)
    rejected_descriptor = equilibrium.reactive_bubble_observation_descriptor(rejected_problem, rows)
    rejected = equilibrium.evaluate_reactive_bubble_observations(
        rejected_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        rejected_descriptor,
    )
    assert rejected.status == "non_evaluable"
    assert rejected.failure is not None
    assert rejected.failure.code == "branch_loss"


@pytest.mark.slow
def test_reactive_bubble_is_invariant_to_component_order_and_reaction_orientation() -> None:
    """Falsify positional component coupling or reaction-row orientation dependence."""

    parameters, base_problem, rows = _fixture()
    base_problem = replace(base_problem, pressure_starts_pa=(3.2e6,))
    base_reaction = base_problem.liquid_problem.reaction_system
    constant = base_reaction.equilibrium_constants[0]
    reversed_reaction = replace(
        base_reaction,
        reaction_matrix=(tuple(-value for value in base_reaction.reaction_matrix[0]),),
        equilibrium_constants=(replace(constant, ln_value=-constant.ln_value),),
    )
    reversed_liquid = _observation_problem(
        equilibrium.ProviderPhase(epcsaft.Mixture(parameters), (1.0e-6, 0.70)),
        reversed_reaction,
        identity="source-neutral-reversed-reaction",
        phase_identity="multineutral-liquid",
        continuation_identity="source-neutral-reversed-liquid-branch",
    )
    reversed_problem = replace(base_problem, liquid_problem=reversed_liquid)

    permutation = (1, 0, 3, 2)
    reordered_parameters = parameters.select(
        tuple(parameters.component_ids[index] for index in permutation)
    )
    reordered_liquid = _observation_problem(
        equilibrium.ProviderPhase(epcsaft.Mixture(reordered_parameters), (1.0e-6, 0.70)),
        _reordered_reaction_system(base_reaction, permutation),
        identity="source-neutral-reordered-components",
        phase_identity="multineutral-liquid",
        continuation_identity="source-neutral-reordered-liquid-branch",
    )
    reordered_problem = replace(base_problem, liquid_problem=reordered_liquid)
    pressure_row = (rows[0],)
    results = []
    for candidate in (base_problem, reversed_problem, reordered_problem):
        descriptor = equilibrium.reactive_bubble_observation_descriptor(candidate, pressure_row)
        results.append(
            equilibrium.evaluate_reactive_bubble_observations(
                candidate,
                313.15 * epcsaft.unit_registry.kelvin,
                descriptor,
            )
        )
    assert all(result.status == "evaluated" for result in results), tuple(
        result.failure for result in results
    )
    assert float(results[1].values[0]) == pytest.approx(float(results[0].values[0]), rel=3.0e-9)
    assert float(results[2].values[0]) == pytest.approx(float(results[0].values[0]), rel=3.0e-9)


@pytest.mark.slow
def test_reactive_bubble_failures_preserve_no_incipient_domain_and_support_codes() -> None:
    """Falsify generic failure collapse or derivative-chart reuse after support loss."""

    parameters, problem, rows = _fixture()
    row = (rows[0],)

    no_incipient_problem = replace(
        problem,
        pressure_interval_pa=(1.0e5, 2.0e5),
        pressure_starts_pa=(1.5e5,),
    )
    no_incipient = equilibrium.evaluate_reactive_bubble_observations(
        no_incipient_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        equilibrium.reactive_bubble_observation_descriptor(no_incipient_problem, row),
    )
    assert no_incipient.failure is not None
    assert no_incipient.failure.code == "no_incipient_phase"

    domain_identity = "component/solvent/segment_diameter"
    domain_trial = epcsaft.ActiveParameterSet(
        parameters,
        (domain_identity,),
        values=(100.0 * epcsaft.unit_registry.angstrom,),
    )
    reduced_problem = replace(problem, pressure_starts_pa=(3.2e6,))
    domain = equilibrium.evaluate_reactive_bubble_observations(
        reduced_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        equilibrium.reactive_bubble_observation_descriptor(
            reduced_problem, row, active_parameters=domain_trial
        ),
        active_parameters=domain_trial,
    )
    assert domain.failure is not None
    assert domain.failure.code == "provider_domain_rejection"

    floor_problem = replace(
        reduced_problem,
        liquid_problem=replace(
            reduced_problem.liquid_problem,
            reaction_system=replace(
                reduced_problem.liquid_problem.reaction_system,
                strict_interior_amount_floor_mol=0.11,
            ),
        ),
    )
    support = equilibrium.evaluate_reactive_bubble_observations(
        floor_problem,
        313.15 * epcsaft.unit_registry.kelvin,
        equilibrium.reactive_bubble_observation_descriptor(floor_problem, row),
    )
    assert support.failure is not None
    assert support.failure.code == "active_bound_ambiguity"


@pytest.mark.slow
def test_multireaction_nonideal_bubble_total_jacobian_matches_resolves() -> None:
    """Falsify single-reaction assumptions in the coupled Schur derivative."""

    parameters, problem, rows = _multireaction_fixture()
    identity = "component/solvent/dispersion_energy_over_k"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))

    def solve(value: float):
        trial = epcsaft.ActiveParameterSet(parameters, (identity,), values=(value,))
        descriptor = equilibrium.reactive_bubble_observation_descriptor(
            problem, rows, active_parameters=trial
        )
        result = equilibrium.evaluate_reactive_bubble_observations(
            problem,
            313.15 * epcsaft.unit_registry.kelvin,
            descriptor,
            active_parameters=trial,
        )
        assert result.status == "evaluated", result.failure
        return result

    center = solve(active.values[0])
    step = 1.0e-3
    below = solve(active.values[0] - step)
    above = solve(active.values[0] + step)
    assert center.certification is not None
    assert center.certification.status == "certified_local_equilibrium"
    assert center.certification.globality_claim == "not_established"
    assert center.search.distinct_branch_count == 1
    assert center.phases[1].supported_component_ids == problem.vapor_component_ids
    for row, exact in enumerate(center.jacobian):
        assert exact is not None
        low = below.values[row]
        high = above.values[row]
        assert low is not None and high is not None
        numerical = (high - low) / (2.0 * step)
        scale = max(1.0, abs(exact[0]), abs(numerical))
        assert abs(exact[0] - numerical) / scale < 2.0e-7


@pytest.mark.slow
def test_reactive_bubble_uses_canonical_regression_fit_and_opaque_receipt(
    tmp_path: Path,
) -> None:
    """Falsify a parallel CAP-13 fit surface or non-final owner receipt."""

    parameters, problem, _ = _fixture()
    problem = replace(problem, pressure_starts_pa=(3.2e6,))
    identity = "component/solvent/dispersion_energy_over_k"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    problem = replace(
        problem,
        continuation_reference=equilibrium.certify_reactive_bubble_continuation_reference(
            problem,
            313.15 * epcsaft.unit_registry.kelvin,
            active_parameters=active,
            maximum_log_pressure_distance=0.5,
            maximum_log_composition_distance=0.5,
            maximum_log_volume_distance=0.5,
        ),
    )
    row = equilibrium.ReactiveBubbleObservationRow(
        "regression-bubble-pressure", "bubble_pressure_pa", support="positive"
    )
    descriptor = equilibrium.reactive_bubble_observation_descriptor(
        problem, (row,), active_parameters=active
    )
    truth = equilibrium.evaluate_reactive_bubble_observations(
        problem,
        313.15 * epcsaft.unit_registry.kelvin,
        descriptor,
        active_parameters=active,
    )
    assert truth.status == "evaluated", truth.failure
    pressure = float(truth.values[0])
    parameter_payload = json.dumps(
        parameters.to_mapping(), separators=(",", ":"), sort_keys=True
    ).encode()
    parameter_path = tmp_path / "reactive-bubble-parameters.json"
    parameter_path.write_bytes(parameter_payload)
    liquid_reference = problem.liquid_problem.continuation_reference
    assert liquid_reference is not None
    bubble_reference = problem.continuation_reference
    assert bubble_reference is not None
    value = active.values[0]
    mapping = {
        "schema_version": 2,
        "parameters": {
            "path": str(parameter_path),
            "sha256": f"sha256:{hashlib.sha256(parameter_payload).hexdigest()}",
        },
        "model": {"kind": "mixture"},
        "active": {
            "identities": [identity],
            "coordinates": {
                identity: {
                    "bounds": [value - 5.0, value + 5.0],
                    "origin": value,
                    "scale": 2.0,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "reactive_bubble_vle",
                "identity": problem.identity,
                "liquid_problem_identity": problem.liquid_problem.identity,
                "liquid_phase": {
                    "identity": problem.liquid_problem.phase_identity,
                    "role": "liquid",
                    "model": "provider",
                    "admissible_packing_fraction_interval": list(
                        problem.liquid_problem.phase.admissible_packing_fraction_interval
                    ),
                },
                "reaction_system": _reaction_mapping(problem.liquid_problem.reaction_system),
                "liquid_continuation": {
                    "identity": problem.liquid_problem.continuation_identity,
                    "branch_policy": problem.liquid_problem.branch_policy,
                    "reference": liquid_reference.to_mapping(),
                },
                "vapor": {
                    "identity": problem.vapor_phase_identity,
                    "component_ids": list(problem.vapor_component_ids),
                    "model": "provider_nonideal",
                    "reference_id": problem.vapor_model.reference_id,
                },
                "bubble_continuation": {
                    "identity": problem.continuation_identity,
                    "reference": bubble_reference.to_mapping(),
                },
                "pressure_interval_pa": list(problem.pressure_interval_pa),
                "pressure_starts_pa": list(problem.pressure_starts_pa),
                "rows": [_row_mapping(row)],
                "temperature_k": 313.15,
                "targets": [
                    _target_mapping(
                        "source-neutral/bubble-pressure",
                        row.identity,
                        pressure,
                        "pascal",
                        "scaled_difference",
                        pressure,
                    )
                ],
                "multiplier": 1.0,
            }
        ],
        "solver": {
            "multistart": {
                "starts": [{"identity": "exact-source-neutral", "values": [value]}],
                "affine_bound_margin": 1.0e-7,
            }
        },
    }
    fitted = regression.fit(mapping, base_path=tmp_path)
    accepted = fitted.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((value,), abs=2.0e-9)
    assert accepted.predictions == pytest.approx((pressure,), rel=2.0e-10)
    assert len(accepted.owner_receipts) == 1
    assert accepted.owner_receipts[0].digest.startswith("sha256:")
