from __future__ import annotations

import json
import math
from copy import deepcopy
from dataclasses import replace

import pytest

import epcsaft
from epcsaft import equilibrium
from tests.eos._raw_cases import MULTINEUTRAL_REACTIVE, NEUTRAL_BINARY

from .test_chemical_provider import _held_parameters, _multineutral_problem


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

    result = equilibrium.solve(model, problem)

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
    temperature = 170.0 * epcsaft.unit_registry.kelvin
    pressure = 1.0e6 * epcsaft.unit_registry.pascal
    state = model.state(T=temperature, P=pressure, x=amounts, phase="liquid")
    phase = equilibrium.ReactivePhase(
        "configured-liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
        start=equilibrium.FinitePhaseStart(
            amounts,
            1.0 / float(state.molar_density.to("mole / meter**3").magnitude),
        ),
    )
    potential = equilibrium.EquilibriumOutput(
        "methane-log-fugacity",
        "phase.log_fugacity",
        "dimensionless",
        "natural-log-pascal-true-species",
        phase.identity,
        (1.0, 0.0),
    )
    truth = equilibrium.solve(
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

    evaluated = equilibrium.solve(model, problem)

    assert evaluated.status == "evaluated", evaluated.failure
    assert evaluated.values == pytest.approx((expected_x,), abs=2.0e-10)


@pytest.mark.slow
def test_imposed_potential_output_jacobian_includes_uptake_response() -> None:
    """Falsify fixed-feed partials masquerading as total regression derivatives."""

    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    model = epcsaft.Mixture(parameters)
    reaction = _multineutral_problem()
    phase = equilibrium.ReactivePhase(
        "multineutral-liquid",
        "liquid",
        "finite",
        equilibrium.AllComponents(),
        equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
    )
    temperature = equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin)
    pressure = equilibrium.Fixed(119_199.78746707023 * epcsaft.unit_registry.pascal)
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
        feed=equilibrium.direct_amounts(reaction.species_ids, (0.8, 0.2, 0.1, 0.1)),
    )
    identity = "pair/anion/solvent/k_ij"
    active = epcsaft.ActiveParameterSet(parameters, (identity,))
    truth = equilibrium.solve(model, fixed_truth, active_parameters=active)
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

    evaluated = equilibrium.solve(model, problem, active_parameters=active)

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
        result = equilibrium.solve(model, problem, active_parameters=trial)
        assert result.status == "evaluated", result.failure
        perturbed_values.append(float(result.values[0]))
    independent = (perturbed_values[1] - perturbed_values[0]) / (2.0 * step)
    assert evaluated.jacobian[0] is not None
    assert evaluated.jacobian[0][0] == pytest.approx(independent, rel=3.0e-5, abs=2.0e-8)
    descriptor = equilibrium.general_reactive_equilibrium_descriptor(
        model, problem, active_parameters=active
    )
    observations = (
        equilibrium.GeneralEquilibriumOwnerObservation("synthetic-cation-target", output.identity),
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
    tampered["problem_fingerprint"] = "sha256:tampered"
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


@pytest.mark.slow
def test_solved_pressure_uses_generic_grepe_route() -> None:
    parameters = epcsaft.Parameters.from_mapping(MULTINEUTRAL_REACTIVE)
    reaction = _multineutral_problem()
    liquid_id = "reactive-liquid"
    vapor_id = "incipient-vapor"
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "generic-incipient-pressure",
        equilibrium.Fixed(313.15 * epcsaft.unit_registry.kelvin),
        equilibrium.Solved(
            3.2e6 * epcsaft.unit_registry.pascal,
            (2.5e6 * epcsaft.unit_registry.pascal, 4.0e6 * epcsaft.unit_registry.pascal),
            (2.8e6 * epcsaft.unit_registry.pascal, 3.2e6 * epcsaft.unit_registry.pascal),
        ),
        (
            equilibrium.ReactivePhase(
                liquid_id,
                "liquid",
                "finite",
                equilibrium.AllComponents(),
                equilibrium.ProviderModel((1.0e-6, 0.70), "installed-provider-eos"),
            ),
            equilibrium.ReactivePhase(
                vapor_id,
                "vapor",
                "incipient",
                equilibrium.DeclaredNeutralComponents(("solvent", "cosolvent")),
                equilibrium.ProviderModel((1.0e-6, 0.70), "installed-provider-eos"),
            ),
        ),
        reaction,
        (liquid_id,),
        (
            equilibrium.EquilibriumOutput(
                "pressure",
                "system.pressure",
                "pascal",
                "true-species-pressure",
                support="positive",
            ),
            equilibrium.EquilibriumOutput(
                "solvent-vapor-fraction",
                "phase.mole_fraction",
                "dimensionless",
                "true-species-vapor-mole-fraction",
                vapor_id,
                (1.0, 0.0),
                support="positive",
            ),
            equilibrium.EquilibriumOutput(
                "solvent-partial-pressure",
                "phase.partial_pressure",
                "pascal",
                "true-species-vapor-partial-pressure",
                vapor_id,
                (1.0, 0.0),
                support="positive",
            ),
        ),
        continuation_identity="generic-incipient-branch",
    )

    result = equilibrium.solve(epcsaft.Mixture(parameters), problem)

    assert result.status == "evaluated", result.failure
    assert 2.5e6 < result.values[0] < 4.0e6
    assert 0.0 < result.values[1] < 1.0
    assert result.values[2] == pytest.approx(result.values[0] * result.values[1], rel=2.0e-13)
    assert len(result.phases) == 2
    assert result.phases[1].supported_component_ids == ("solvent", "cosolvent")
    assert result.phases[0].chemical_potential_over_rt[:2] == pytest.approx(
        result.phases[1].chemical_potential_over_rt, rel=2.0e-8, abs=2.0e-8
    )
    vapor_oracle = epcsaft.Mixture(parameters.select(("solvent", "cosolvent"))).state(
        T=313.15 * epcsaft.unit_registry.kelvin,
        P=result.values[0] * epcsaft.unit_registry.pascal,
        x=result.phases[1].mole_fractions,
        phase="vapor",
    )
    assert result.phases[1].molar_volume_m3_per_mol == pytest.approx(
        1.0 / float(vapor_oracle.molar_density.to("mole / meter**3").magnitude),
        rel=2.0e-9,
    )
    assert result.continuation_state is not None
    assert dict(result.evidence)["bubble_log_closure_abs"] < 2.0e-10

    pure = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            identity="generic-pure-incipient-pressure",
            pressure=equilibrium.Solved(
                10000.0 * epcsaft.unit_registry.pascal,
                (
                    1000.0 * epcsaft.unit_registry.pascal,
                    100000.0 * epcsaft.unit_registry.pascal,
                ),
                (10000.0 * epcsaft.unit_registry.pascal,),
            ),
            phases=(
                problem.phases[0],
                replace(
                    problem.phases[1],
                    support=equilibrium.DeclaredNeutralComponents(("solvent",)),
                ),
            ),
            outputs=problem.outputs[:1],
            continuation_identity="generic-pure-incipient-branch",
        ),
    )
    assert pure.status == "evaluated", pure.failure
    assert pure.phases[1].supported_component_ids == ("solvent",)
    assert pure.phases[1].mole_fractions == (1.0,)

    ideal = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            identity="generic-ideal-incipient-pressure",
            phases=(
                problem.phases[0],
                replace(
                    problem.phases[1],
                    model=equilibrium.IdealGasModel("ideal-gas"),
                ),
            ),
        ),
    )
    assert ideal.status == "evaluated", ideal.failure
    assert 2.5e6 < ideal.values[0] < 4.0e6
    assert ideal.phases[1].mechanical_class == "analytical_ideal_gas"

    boundary = replace(
        problem,
        identity="fixed-pressure-incipient-composition-boundary",
        feed=equilibrium.AffineFeed(
            reaction.species_ids,
            reaction.feed_amounts_mol,
            ((1.0,), (-1.0,), (0.0,), (0.0,)),
            (equilibrium.FeedCoordinate("composition-shift", (-0.02, 0.02)),),
        ),
        intensive_boundaries=(
            equilibrium.IntensiveBoundary(
                "imposed-bubble-pressure",
                "composition-shift",
                None,
                "system.pressure",
                (),
                float(result.values[0]),
                0.0,
            ),
        ),
    )
    boundary_result = equilibrium.solve(epcsaft.Mixture(parameters), boundary)
    assert boundary_result.status == "evaluated", boundary_result.failure
    assert boundary_result.values == pytest.approx(result.values, rel=2.0e-9, abs=2.0e-8)
    assert boundary_result.branch_signature is not None
    assert boundary_result.branch_signature[1][0][1] == pytest.approx(0.0, abs=2.0e-10)

    continued = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(problem, continuation_state=result.continuation_state),
    )
    assert continued.values == pytest.approx(result.values, rel=2.0e-12, abs=2.0e-12)
    reversed_starts = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            pressure=replace(problem.pressure, starts=tuple(reversed(problem.pressure.starts))),
        ),
    )
    assert reversed_starts.values == pytest.approx(result.values, rel=2.0e-12, abs=2.0e-12)

    finite = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            pressure=equilibrium.Fixed(result.values[0] * epcsaft.unit_registry.pascal),
            phases=(problem.phases[0],),
            outputs=problem.outputs[:1],
            continuation_identity=None,
        ),
    )
    assert finite.status == "evaluated", finite.failure
    assert result.phases[0].mole_fractions == pytest.approx(
        finite.phases[0].mole_fractions, abs=2.0e-13
    )
    assert result.phases[0].molar_volume_m3_per_mol == pytest.approx(
        finite.phases[0].molar_volume_m3_per_mol, rel=2.0e-9
    )

    parameter_id = "component/solvent/segment_diameter"
    center = epcsaft.ActiveParameterSet(parameters, (parameter_id,))
    exact = equilibrium.solve(
        epcsaft.Mixture(parameters), problem, active_parameters=center
    )
    assert exact.status == "evaluated"
    step = 1.0e-5
    sides = tuple(
        equilibrium.solve(
            epcsaft.Mixture(parameters),
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (parameter_id,), values=(center.values[0] + sign * step,)
            ),
        )
        for sign in (-1.0, 1.0)
    )
    continued_neighbor = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(problem, continuation_state=result.continuation_state),
        active_parameters=epcsaft.ActiveParameterSet(
            parameters, (parameter_id,), values=(center.values[0] + step,)
        ),
    )
    assert continued_neighbor.values == pytest.approx(
        sides[1].values, rel=2.0e-10, abs=2.0e-9
    )
    for row, lower, upper in zip(exact.rows, sides[0].values, sides[1].values, strict=True):
        assert row.jacobian is not None and lower is not None and upper is not None
        finite_difference = (float(upper) - float(lower)) / (2.0 * step)
        derivative = row.jacobian[0]
        assert abs(derivative - finite_difference) / max(
            1.0, abs(derivative), abs(finite_difference)
        ) < 2.0e-7

    boundary_exact = equilibrium.solve(
        epcsaft.Mixture(parameters), boundary, active_parameters=center
    )
    boundary_sides = tuple(
        equilibrium.solve(
            epcsaft.Mixture(parameters),
            boundary,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (parameter_id,), values=(center.values[0] + sign * step,)
            ),
        )
        for sign in (-1.0, 1.0)
    )
    assert boundary_exact.status == "evaluated", boundary_exact.failure
    assert all(side.status == "evaluated" for side in boundary_sides)
    for row, lower, upper in zip(
        boundary_exact.rows,
        boundary_sides[0].values,
        boundary_sides[1].values,
        strict=True,
    ):
        assert row.jacobian is not None and lower is not None and upper is not None
        finite_difference = (float(upper) - float(lower)) / (2.0 * step)
        assert row.jacobian[0] == pytest.approx(
            finite_difference, rel=5.0e-5, abs=2.0e-9
        )

    constant = reaction.equilibrium_constants[0]
    reversed_reaction = replace(
        reaction,
        reaction_matrix=(tuple(-value for value in reaction.reaction_matrix[0]),),
        equilibrium_constants=(replace(constant, ln_value=-constant.ln_value),),
    )
    permutation = (1, 0, 3, 2)
    permuted_reaction = replace(
        reaction,
        species_ids=tuple(reaction.species_ids[index] for index in permutation),
        charges=tuple(reaction.charges[index] for index in permutation),
        molar_masses_kg_per_mol=tuple(
            reaction.molar_masses_kg_per_mol[index] for index in permutation
        ),
        balance_matrix=tuple(
            tuple(row[index] for index in permutation) for row in reaction.balance_matrix
        ),
        reaction_matrix=tuple(
            tuple(row[index] for index in permutation) for row in reaction.reaction_matrix
        ),
        feed_amounts_mol=tuple(reaction.feed_amounts_mol[index] for index in permutation),
    )
    for model, candidate in (
        (epcsaft.Mixture(parameters), replace(problem, reaction_system=reversed_reaction)),
        (
            epcsaft.Mixture(
                parameters.select(tuple(parameters.component_ids[index] for index in permutation))
            ),
            replace(problem, reaction_system=permuted_reaction),
        ),
    ):
        invariant = equilibrium.solve(model, candidate)
        assert invariant.status == "evaluated", (model.component_ids, invariant.failure)
        assert invariant.values[0] == pytest.approx(result.values[0], rel=3.0e-9)

    no_root = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            pressure=equilibrium.Solved(
                1.5e5 * epcsaft.unit_registry.pascal,
                (1.0e5 * epcsaft.unit_registry.pascal, 2.0e5 * epcsaft.unit_registry.pascal),
                (1.5e5 * epcsaft.unit_registry.pascal,),
            ),
        ),
    )
    assert no_root.status == "non_evaluable"
    assert no_root.failure is not None
    assert no_root.failure.code == "no_incipient_phase"

    nonconverged = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            pressure=equilibrium.Solved(
                1.0e-3 * epcsaft.unit_registry.pascal,
                (1.0e-3 * epcsaft.unit_registry.pascal, 4.0e6 * epcsaft.unit_registry.pascal),
                (1.0e-3 * epcsaft.unit_registry.pascal,),
            ),
        ),
    )
    assert nonconverged.status == "non_evaluable"
    assert nonconverged.failure is not None
    assert nonconverged.failure.code == "numerical_convergence_failure"

    assert result.continuation_state is not None
    vapor_index = result.continuation_state.phase_identities.index(vapor_id)
    phases = list(result.continuation_state.phases)
    phases[vapor_index] = replace(
        phases[vapor_index], mole_fractions=(1.0e-300, 1.0 - 1.0e-300)
    )
    lost = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            continuation_state=replace(result.continuation_state, phases=tuple(phases)),
        ),
    )
    assert lost.status == "non_evaluable"
    assert lost.failure is not None
    assert lost.failure.code == "branch_loss"

    support_loss = equilibrium.solve(
        epcsaft.Mixture(parameters),
        replace(
            problem,
            reaction_system=replace(reaction, strict_interior_amount_floor_mol=0.11),
        ),
    )
    assert support_loss.status == "non_evaluable"
    assert support_loss.failure is not None
    assert support_loss.failure.code == "active_bound_ambiguity", support_loss.failure
