from __future__ import annotations

import hashlib
import json
import math
from copy import deepcopy

import pytest

import epcsaft
from epcsaft import equilibrium
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE, NEUTRAL_BINARY

from .test_chemical_provider import _held_parameters, _multineutral_problem
from .test_homogeneous_reactive_observations import _observation_problem


def test_imposed_log_fugacity_recovers_analytic_open_component_uptake() -> None:
    """Falsify hidden feed sweeps and incorrect open-component balance coupling."""

    component_ids = ("water", "hydronium-cation", "hydroxide-anion")
    model = epcsaft.Mixture(_held_parameters(component_ids))
    pressure_pa = 200_000.0
    closed_amount = 1.5
    expected_uptake = 0.375
    target_log_fugacity = math.log(
        pressure_pa * expected_uptake / (closed_amount + expected_uptake)
    )
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=component_ids,
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(0.01801528, 0.01902322, 0.01700734),
        balance_matrix=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        conserved_totals=(1.0, 0.25, 0.5),
        reaction_matrix=(),
        feed_amounts_mol=(1.0, 0.25, 0.5),
        equilibrium_constants=(),
        strict_interior_amount_floor_mol=1.0e-14,
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        identity="manufactured-open-component-ideal-gas",
        temperature=equilibrium.Fixed(350.0 * epcsaft.unit_registry.kelvin),
        pressure=equilibrium.Fixed(pressure_pa * epcsaft.unit_registry.pascal),
        phases=(
            equilibrium.ReactivePhase(
                "configured-vapor",
                "vapor",
                "finite",
                equilibrium.AllComponents(),
                equilibrium.IdealGasModel("provider-helmholtz-coordinate-basis"),
                "manufactured-open-component-branch",
                "local_certified_role_selected_state",
            ),
        ),
        reaction_system=reaction,
        reaction_phase_ids=(),
        outputs=(
            equilibrium.EquilibriumOutput(
                "open-component-uptake",
                "phase.amount",
                "mole",
                "true-species-amount",
                "configured-vapor",
                (0.0, 1.0, 0.0),
                support="positive",
            ),
        ),
        feed=equilibrium.AffineFeed(
            component_ids,
            (1.0, 0.0, 0.5),
            ((0.0,), (1.0,), (0.0,)),
            (equilibrium.FeedCoordinate("open-component", (0.05, 1.5)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "reservoir-potential",
                "open-component",
                "configured-vapor",
                "phase.log_fugacity",
                (0.0, 1.0, 0.0),
                target_log_fugacity,
                0.2,
            ),
        ),
    )

    result = equilibrium.phase_equilibrium(model, problem)

    assert result.status == "evaluated", result.failure
    assert result.values == pytest.approx((expected_uptake,), rel=0.0, abs=2.0e-11)
    assert result.jacobian == ((),)
    assert result.certificate_scope == "fixed_TP_affine_feed_imposed_intensive_boundaries"
    assert result.branch_signature is not None
    mapping = equilibrium.general_reactive_equilibrium_problem_to_mapping(problem)
    assert equilibrium.general_reactive_equilibrium_problem_from_mapping(mapping) == problem


def test_imposed_log_fugacity_accepts_neutral_provider_mixture() -> None:
    """Falsify an electrolyte-only callback gate in generic homogeneous GREPE."""

    parameter_mapping = deepcopy(NEUTRAL_BINARY)
    parameter_mapping["domains"][0] = {
        "domain_id": "engine-numerical-fixture-domain",
        "kind": "fit-range",
        "temperature_min": {"magnitude": 150.0, "unit": "kelvin"},
        "temperature_max": {"magnitude": 300.0, "unit": "kelvin"},
    }
    parameters = epcsaft.Parameters.from_mapping(parameter_mapping)
    model = epcsaft.Mixture(parameters)
    component_ids = parameters.component_ids
    expected_x = 0.2
    amounts = (expected_x, 1.0 - expected_x)
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=component_ids,
        charges=(0, 0),
        molar_masses_kg_per_mol=(0.016043, 0.03007),
        balance_matrix=((1.0, 0.0), (0.0, 1.0)),
        conserved_totals=amounts,
        reaction_matrix=(),
        feed_amounts_mol=amounts,
        equilibrium_constants=(),
        strict_interior_amount_floor_mol=1.0e-14,
    )
    provider_phase = equilibrium.ProviderPhase(model, (1.0e-6, 0.74))
    anchor = equilibrium.HomogeneousReactiveObservationProblem(
        "neutral-provider-anchor",
        "configured-liquid",
        "liquid",
        provider_phase,
        reaction,
        "neutral-provider-branch",
        "local_certified_role_selected_state",
    )
    temperature = 170.0 * epcsaft.unit_registry.kelvin
    pressure = 1.0e6 * epcsaft.unit_registry.pascal
    reference = equilibrium.certify_homogeneous_continuation_reference(
        anchor,
        temperature,
        pressure,
        maximum_log_composition_distance=1.0,
        maximum_log_volume_distance=1.0,
    )
    phase = equilibrium.ReactivePhase(
        "configured-liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
        "neutral-provider-branch",
        "local_certified_role_selected_state",
        reference,
    )
    potential = equilibrium.EquilibriumOutput(
        "methane-log-fugacity",
        "phase.log_fugacity",
        "dimensionless",
        "natural-log-pascal-true-species",
        phase.identity,
        (1.0, 0.0),
    )
    truth = equilibrium.phase_equilibrium(
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "neutral-provider-truth",
            equilibrium.Fixed(temperature),
            equilibrium.Fixed(pressure),
            (phase,),
            reaction,
            (),
            (potential,),
            feed=equilibrium.direct_amounts(component_ids, amounts),
        ),
    )
    assert truth.status == "evaluated", truth.failure
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "neutral-provider-uptake",
        equilibrium.Fixed(temperature),
        equilibrium.Fixed(pressure),
        (phase,),
        reaction,
        (),
        (
            equilibrium.EquilibriumOutput(
                "methane-x",
                "phase.mole_fraction",
                "dimensionless",
                "true-species-liquid-mole-fraction",
                phase.identity,
                (1.0, 0.0),
                support="positive",
            ),
        ),
        feed=equilibrium.AffineFeed(
            component_ids,
            (0.0, 0.8),
            ((1.0,), (0.0,)),
            (equilibrium.FeedCoordinate("open-methane", (0.02, 0.5)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "methane-reservoir",
                "open-methane",
                phase.identity,
                "phase.log_fugacity",
                (1.0, 0.0),
                float(truth.values[0]),
                amounts[0],
            ),
        ),
    )

    evaluated = equilibrium.phase_equilibrium(model, problem)

    assert evaluated.status == "evaluated", evaluated.failure
    assert evaluated.values == pytest.approx((expected_x,), abs=2.0e-10)


@pytest.mark.slow
def test_imposed_potential_output_jacobian_includes_uptake_response() -> None:
    """Falsify fixed-feed partials masquerading as total regression derivatives."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    reaction = _multineutral_problem()
    anchor = _observation_problem(
        equilibrium.ProviderPhase(model, (1.0e-6, 0.74)), reaction
    )
    phase = equilibrium.ReactivePhase(
        anchor.phase_identity,
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
        anchor.continuation_identity,
        anchor.branch_policy,
        anchor.continuation_reference,
    )
    temperature = equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin)
    pressure = equilibrium.Fixed(
        119_199.78746707023 * epcsaft.unit_registry.pascal
    )
    potential_output = equilibrium.EquilibriumOutput(
        "cosolvent-log-fugacity",
        "phase.log_fugacity",
        "dimensionless",
        "natural-log-pascal-true-species",
        phase.identity,
        (0.0, 1.0, 0.0, 0.0),
    )
    fixed_truth = equilibrium.GeneralReactiveEquilibriumProblem(
        "open-component-provider-truth",
        temperature,
        pressure,
        (phase,),
        reaction,
        (phase.identity,),
        (potential_output,),
        feed=equilibrium.direct_amounts(
            reaction.species_ids, (0.8, 0.2, 0.1, 0.1)
        ),
    )
    identity = "pair/anion/solvent/k_ij"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    truth = equilibrium.phase_equilibrium(model, fixed_truth, active_parameters=active)
    assert truth.status == "evaluated", truth.failure
    target = float(truth.values[0])
    output = equilibrium.EquilibriumOutput(
        "cation-mole-fraction",
        "phase.mole_fraction",
        "dimensionless",
        "true-species-liquid-mole-fraction",
        phase.identity,
        (0.0, 0.0, 1.0, 0.0),
        support="positive",
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "open-component-provider",
        temperature,
        pressure,
        (phase,),
        reaction,
        (phase.identity,),
        (output,),
        feed=equilibrium.AffineFeed(
            reaction.species_ids,
            (0.8, 0.1, 0.1, 0.1),
            ((0.0,), (1.0,), (0.0,), (0.0,)),
            (equilibrium.FeedCoordinate("open-neutral", (0.02, 0.5)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "imposed-neutral-potential",
                "open-neutral",
                phase.identity,
                "phase.log_fugacity",
                potential_output.coefficients,
                target,
                0.16,
            ),
        ),
    )

    evaluated = equilibrium.phase_equilibrium(model, problem, active_parameters=active)

    assert evaluated.status == "evaluated", evaluated.failure
    assert evaluated.values[0] is not None
    assert 0.0 < evaluated.values[0] < 1.0
    step = 2.0e-5
    perturbed_values = []
    for sign in (-1.0, 1.0):
        trial = epcsaft.ActiveParameterSet(
            parameters,
            (identity,),
            values=(active.values[0] + sign * step,),
        )
        result = equilibrium.phase_equilibrium(model, problem, active_parameters=trial)
        assert result.status == "evaluated", result.failure
        perturbed_values.append(float(result.values[0]))
    independent = (perturbed_values[1] - perturbed_values[0]) / (2.0 * step)
    assert evaluated.jacobian[0] is not None
    assert evaluated.jacobian[0][0] == pytest.approx(independent, rel=3.0e-5, abs=2.0e-8)
    descriptor = equilibrium.general_reactive_equilibrium_descriptor(
        model, problem, active_parameters=active
    )
    observations = (
        equilibrium.GeneralEquilibriumOwnerObservation(
            "synthetic-cation-target", output.identity
        ),
    )
    receipt = equilibrium.serialize_general_reactive_equilibrium_receipt(
        evaluated,
        batch_identity="open-component-provider-batch",
        observations=observations,
    )
    assert equilibrium.validate_general_reactive_equilibrium_receipt(
        receipt,
        descriptor=descriptor,
        problem=problem,
        batch_identity="open-component-provider-batch",
        active_parameters=active,
        observations=observations,
    ).startswith("sha256:")

    tampered = json.loads(receipt)
    owner_receipt = tampered["owner_receipt"]
    owner_receipt["problem_fingerprint"] = "sha256:tampered"
    owner_payload = dict(owner_receipt)
    owner_payload.pop("fingerprint")
    fingerprint_payload = json.dumps(
        {
            "schema": "epcsaft.intensive-boundary-owner-receipt",
            "value": owner_payload,
        },
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    owner_receipt["fingerprint"] = f"sha256:{hashlib.sha256(fingerprint_payload).hexdigest()}"
    tampered_json = json.dumps(
        tampered,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    with pytest.raises(ValueError, match="identity"):
        equilibrium.validate_general_reactive_equilibrium_receipt(
            tampered_json,
            descriptor=descriptor,
            problem=problem,
            batch_identity="open-component-provider-batch",
            active_parameters=active,
            observations=observations,
        )
