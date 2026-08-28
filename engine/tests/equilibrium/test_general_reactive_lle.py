from __future__ import annotations

import hashlib
import json
import math
from copy import deepcopy
from dataclasses import replace
from pathlib import Path

import pytest

import epcsaft
from epcsaft import _core, equilibrium, regression
from tests.eos._cases import (
    ASCANI_2023_HYPOTHETICAL,
    CHARGED_TERNARY_IMMISCIBLE,
    MULTINEUTRAL_REACTIVE,
    NEUTRAL_TERNARY_IMMISCIBLE,
)

_BASIS = "manufactured-reactive-lle-basis"
_PARAMETER_IDENTITY = "pair/alpha/beta/k_ij"


def _fixture() -> tuple[
    epcsaft.Parameters,
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
]:
    parameters = epcsaft.Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    model = epcsaft.Mixture(parameters)
    seed = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTPFeed(
                300.0 * epcsaft.unit_registry.kelvin,
                1.0e5 * epcsaft.unit_registry.pascal,
                (1.0, 1.0, 1.0),
                equilibrium.DeclaredTopology(
                    tuple(
                        equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents())
                        for _ in range(3)
                    )
                ),
            )
        ),
    )
    aqueous, organic = seed.phases[:2]
    starts = (aqueous.mole_fractions, organic.mole_fractions)
    feed = tuple(starts[0][index] + starts[1][index] for index in range(3))
    ln_k = aqueous.chemical_potential_over_rt[1] - aqueous.chemical_potential_over_rt[0]
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=model.component_ids,
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(0.03, 0.03, 0.03),
        balance_matrix=((1.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        conserved_totals=(feed[0] + feed[1], feed[2]),
        reaction_matrix=((-1.0, 1.0, 0.0),),
        feed_amounts_mol=feed,
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                ln_k,
                "manufactured-reactive-lle-oracle",
                "provider-helmholtz-coordinate-basis",
                "products_positive",
                "already-provider-basis",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-12,
    )
    phases = tuple(
        equilibrium.ReactivePhase(
            identity,
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-8, 0.74), "installed-provider-eos"),
            start=equilibrium.FinitePhaseStart(
                amounts,
                state.molar_volume_m3_per_mol,
            ),
        )
        for identity, amounts, state in zip(
            ("aqueous", "organic"),
            starts,
            (aqueous, organic),
            strict=True,
        )
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"amount-{phase.identity}-{component_id}",
            "phase.amount",
            "mole",
            _BASIS,
            phase.identity,
            tuple(float(index == component) for index in range(3)),
            support="positive",
        )
        for phase in phases
        for component, component_id in enumerate(model.component_ids)
    ) + (
        equilibrium.EquilibriumOutput(
            "aqueous-alpha-mole-fraction",
            "phase.mole_fraction",
            "dimensionless",
            _BASIS,
            "aqueous",
            (1.0, 0.0, 0.0),
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "organic-molar-density",
            "phase.molar_density",
            "mole / meter**3",
            _BASIS,
            "organic",
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "aqueous-pressure",
            "phase.pressure",
            "pascal",
            _BASIS,
            "aqueous",
            support="positive",
        ),
        equilibrium.EquilibriumOutput(
            "aqueous-beta-chemical-potential",
            "phase.chemical_potential_over_rt",
            "dimensionless",
            _BASIS,
            "aqueous",
            (0.0, 1.0, 0.0),
        ),
    )
    return (
        parameters,
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "manufactured-reactive-lle",
            equilibrium.Fixed(300.0 * epcsaft.unit_registry.kelvin),
            equilibrium.Fixed(1.0e5 * epcsaft.unit_registry.pascal),
            phases,
            reaction,
            ("aqueous",),
            outputs,
        ),
    )


def _explicit_phase_direction_problem(
    problem: equilibrium.GeneralReactiveEquilibriumProblem,
) -> equilibrium.GeneralReactiveEquilibriumProblem:
    transfer_rows = tuple(
        equilibrium.PhaseReactionRow(
            f"transfer-{species_id}",
            (
                equilibrium.PhaseSpeciesTerm("aqueous", species_id, -1.0),
                equilibrium.PhaseSpeciesTerm("organic", species_id, 1.0),
            ),
        )
        for species_id in problem.reaction_system.species_ids
    )
    reaction_rows = tuple(
        equilibrium.PhaseReactionRow(
            f"aqueous-reaction-{index}",
            tuple(
                equilibrium.PhaseSpeciesTerm("aqueous", species_id, coefficient)
                for species_id, coefficient in zip(
                    problem.reaction_system.species_ids,
                    stoichiometry,
                    strict=True,
                )
                if coefficient != 0.0
            ),
            ln_equilibrium_constant=record.ln_value,
        )
        for index, (stoichiometry, record) in enumerate(
            zip(
                problem.reaction_system.reaction_matrix,
                problem.reaction_system.equilibrium_constants,
                strict=True,
            )
        )
    )
    return replace(
        problem,
        phase_reactions=equilibrium.PhaseReactionSystem(
            tuple(phase.identity for phase in problem.phases),
            problem.reaction_system.species_ids,
            transfer_rows + reaction_rows,
            tuple(float(value) for value in problem.reaction_system.charges),
            problem.reaction_system.balance_matrix,
        ),
    )


def test_finite_phase_is_one_certified_coupled_minimum_and_label_covariant() -> None:
    """Falsify staged speciation, missing transfer equilibrium, or label ordering."""

    parameters, model, problem = _fixture()
    active = epcsaft.ActiveParameterSet(parameters, (_PARAMETER_IDENTITY,))
    result = equilibrium.solve(model, problem, active_parameters=active)

    assert result.status == "evaluated", result.failure
    assert result.numerical_status == "passed"
    assert result.physical_status == "passed"
    assert result.provider_domain_status == "passed"
    assert result.descriptor is not None
    assert result.descriptor.formulation == "finite"
    assert result.certificate_scope == ("fixed_TP_configured_finite_phases_local_minimum")
    assert result.globality_claim == "not_established"
    assert len(result.phases) == 2
    phase_amounts = {
        phase.identity: tuple(
            result.rows[phase_index * 3 + component].value for component in range(3)
        )
        for phase_index, phase in enumerate(problem.phases)
    }
    assert all(value is not None for values in phase_amounts.values() for value in values)
    reconstructed_species = tuple(
        sum(float(values[component]) for values in phase_amounts.values()) for component in range(3)
    )
    reconstructed_conserved = (
        reconstructed_species[0] + reconstructed_species[1],
        reconstructed_species[2],
    )
    assert reconstructed_conserved == pytest.approx(
        problem.reaction_system.conserved_totals, abs=2.0e-10
    )
    aqueous, organic = result.phases
    for left, right in zip(
        aqueous.chemical_potential_over_rt,
        organic.chemical_potential_over_rt,
        strict=True,
    ):
        assert left == pytest.approx(right, abs=3.0e-8)
    for phase in result.phases:
        assert (
            -phase.chemical_potential_over_rt[0] + phase.chemical_potential_over_rt[1]
        ) == pytest.approx(0.0, abs=3.0e-8)
        assert phase.pressure_pa == pytest.approx(1.0e5, rel=7.0e-9)
    assert (
        max(
            abs(left - right)
            for left, right in zip(aqueous.mole_fractions, organic.mole_fractions, strict=True)
        )
        > 0.5
    )

    observations = tuple(
        equilibrium.GeneralEquilibriumOwnerObservation(output.identity, output.identity)
        for output in problem.outputs
    )
    receipt = equilibrium.serialize_general_reactive_equilibrium_receipt(
        result,
        batch_identity="manufactured/reactive-lle",
        observations=observations,
    )
    assert equilibrium.validate_general_reactive_equilibrium_receipt(
        receipt,
        descriptor=result.descriptor,
        problem=problem,
        batch_identity="manufactured/reactive-lle",
        active_parameters=active,
        observations=observations,
    ).startswith("sha256:")
    evidence = dict(json.loads(receipt)["evidence"])
    assert evidence["declared_direction_rank"] == evidence["admissible_nullity"] == 4
    assert evidence["reduced_hessian_min_eigenvalue"] > 0.0

    permuted = replace(
        problem,
        phases=tuple(reversed(problem.phases)),
        reaction_phase_ids=tuple(reversed(problem.reaction_phase_ids)),
    )
    swapped = equilibrium.solve(model, permuted, active_parameters=active)
    assert swapped.status == "evaluated", swapped.failure
    states = {
        phase.identity: state for phase, state in zip(problem.phases, result.phases, strict=True)
    }
    swapped_states = {
        phase.identity: state for phase, state in zip(permuted.phases, swapped.phases, strict=True)
    }
    for identity, state in states.items():
        assert swapped_states[identity].mole_fractions == pytest.approx(
            state.mole_fractions, abs=2.0e-10
        )


def test_explicit_phase_directions_reproduce_legacy_admissible_space() -> None:
    """Explicit transfer/reaction rows must drive the same coupled Gibbs problem."""

    parameters, model, problem = _fixture()
    explicit = _explicit_phase_direction_problem(problem)

    active = epcsaft.ActiveParameterSet(parameters, (_PARAMETER_IDENTITY,))
    legacy = equilibrium.solve(model, problem, active_parameters=active)
    directed = equilibrium.solve(model, explicit, active_parameters=active)

    assert legacy.status == directed.status == "evaluated"
    for legacy_phase, directed_phase in zip(legacy.phases, directed.phases, strict=True):
        assert directed_phase.mole_fractions == pytest.approx(
            legacy_phase.mole_fractions, rel=2.0e-7, abs=2.0e-10
        )
        assert directed_phase.molar_volume_m3_per_mol == pytest.approx(
            legacy_phase.molar_volume_m3_per_mol, rel=2.0e-7
        )
    for legacy_row, directed_row in zip(legacy.rows, directed.rows, strict=True):
        assert directed_row.jacobian == pytest.approx(legacy_row.jacobian, rel=2.0e-6, abs=2.0e-8)


def test_declared_finite_phase_start_recovers_from_collapsed_candidate() -> None:
    """Keep a valid caller start first, then recover through canonical partitions."""

    _, model, problem = _fixture()
    feed = tuple(
        math.fsum(phase.start.amounts_mol[component] for phase in problem.phases)
        for component in range(3)
    )
    collapsed_starts = tuple(
        equilibrium.FinitePhaseStart(
            tuple(amount / 2.0 for amount in feed),
            phase.start.molar_volume_m3_per_mol,
        )
        for phase in problem.phases
    )
    collapsed_problem = replace(
        problem,
        identity="manufactured-reactive-lle-collapsed-start",
        phases=tuple(
            replace(phase, start=start)
            for phase, start in zip(problem.phases, collapsed_starts, strict=True)
        ),
    )
    collapsed = _core._solve_grepe(
        model._native, model._native, collapsed_problem, None, "", False, True
    ).equilibrium
    generated = _core._solve_grepe(
        model._native,
        model._native,
        replace(
            collapsed_problem,
            phases=tuple(replace(phase, start=None) for phase in collapsed_problem.phases),
        ),
        None,
        "",
    ).equilibrium
    result = equilibrium.solve(model, collapsed_problem)

    assert result.status == "evaluated", result.failure
    assert not collapsed.accepted
    assert generated.accepted
    assert collapsed.reduced_hessian_min_eigenvalue < 0.0
    assert result.physical_status == "passed"
    evidence = dict(result.evidence)
    assert evidence["provider_evaluation_count"] == (
        collapsed.provider_evaluation_count + generated.provider_evaluation_count
    )
    assert evidence["phase_distinctness"] > 1.0e-6
    assert evidence["reduced_hessian_min_eigenvalue"] > 0.0
    for phase in result.phases:
        assert phase.pressure_pa == pytest.approx(1.0e5, rel=7.0e-9)
    for left, right in zip(
        result.phases[0].chemical_potential_over_rt,
        result.phases[1].chemical_potential_over_rt,
        strict=True,
    ):
        assert left == pytest.approx(right, abs=3.0e-8)


def test_reduced_tangent_hessian_rejects_genuine_negative_curvature() -> None:
    """Reject an analytic Provider saddle independently of the HELD campaign."""

    _, model, problem = _fixture()
    feed = (1.0, 1.0, 0.01)
    mole_fractions = tuple(amount / math.fsum(feed) for amount in feed)
    state = model.state(
        T=300.0 * epcsaft.unit_registry.kelvin,
        P=1.0e5 * epcsaft.unit_registry.pascal,
        x=mole_fractions,
        phase="liquid",
    )
    molar_volume = 1.0 / float(state.molar_density.to("mole / meter**3").magnitude)
    reaction = equilibrium.ChemicalEquilibriumProblem(
        model.component_ids,
        (0, 0, 0),
        (0.03, 0.03, 0.03),
        (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        ),
        feed,
        (),
        feed,
        (),
        1.0e-12,
    )
    problem = replace(
        problem,
        identity="analytic-provider-saddle",
        reaction_system=reaction,
        reaction_phase_ids=(),
        phases=tuple(
            replace(
                phase,
                model=equilibrium.ProviderModel(
                    (1.0e-6, 0.74), "installed-provider-eos"
                ),
                start=equilibrium.FinitePhaseStart(
                    tuple(amount / 2.0 for amount in feed), molar_volume
                ),
            )
            for phase in problem.phases
        ),
    )
    native = _core._solve_grepe(
        model._native,
        model._native,
        problem,
        None,
        "",
        False,
        True,
    ).equilibrium

    assert not native.accepted
    assert native.balance_inf_norm <= 1.0e-9
    assert native.pressure_relative_inf_norm <= 1.0e-8
    assert native.transfer_chemical_potential_inf_norm <= 1.0e-6
    assert native.kkt_stationarity_inf_norm <= 1.0e-7
    assert native.minimum_bound_margin > 1.0e-7
    assert native.reduced_hessian_min_eigenvalue < -1.0
    assert native.reduced_hessian_threshold == pytest.approx(
        1.0e-10 * native.reduced_hessian_scale, rel=2.0e-15
    )


def test_phase_directions_cannot_create_mass_beyond_the_amount_cap() -> None:
    """Reject a count-conserving direction that violates the cap's positive mass invariant."""

    _, model, _ = _fixture()
    feed = (1.0, 1.0, 1.0)
    reaction = equilibrium.ChemicalEquilibriumProblem(
        model.component_ids,
        (0, 0, 0),
        (0.02, 0.03, 0.04),
        ((1.0, 1.0, 1.0),),
        (3.0,),
        (),
        feed,
        (),
        1.0e-12,
    )
    phases = tuple(
        equilibrium.ReactivePhase(
            identity,
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-6, 0.74), "installed-provider-eos"),
            equilibrium.FinitePhaseStart((0.5, 0.5, 0.5), 1.0e-4),
        )
        for identity in ("left", "right")
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "mass-creating-direction",
        equilibrium.Fixed(300.0 * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(1.0e5 * epcsaft.unit_registry.pascal),
        phases,
        reaction,
        (),
        (
            equilibrium.EquilibriumOutput(
                "left-A", "phase.amount", "mole", _BASIS, "left", (1.0, 0.0, 0.0)
            ),
        ),
        phase_reactions=equilibrium.PhaseReactionSystem(
            ("left", "right"),
            model.component_ids,
            (
                equilibrium.PhaseReactionRow(
                    "mass-creating",
                    (
                        equilibrium.PhaseSpeciesTerm("left", model.component_ids[0], -1.0),
                        equilibrium.PhaseSpeciesTerm("left", model.component_ids[1], 1.0),
                    ),
                ),
            ),
            (0.0, 0.0, 0.0),
            ((1.0, 1.0, 1.0),),
        ),
    )
    with pytest.raises(ValueError, match="global molar-mass conservation"):
        _core._solve_grepe(model._native, model._native, problem, None, "")


def test_explicit_phase_directions_are_the_reactive_authority() -> None:
    """A configured direction space must not require a duplicate legacy reaction basis."""

    parameters, model, problem = _fixture()
    explicit = _explicit_phase_direction_problem(problem)
    direction_only = replace(
        explicit,
        reaction_system=replace(
            explicit.reaction_system,
            reaction_matrix=(),
            equilibrium_constants=(),
        ),
        reaction_phase_ids=(),
    )
    active = epcsaft.ActiveParameterSet(parameters, (_PARAMETER_IDENTITY,))

    expected = equilibrium.solve(model, explicit, active_parameters=active)
    actual = equilibrium.solve(model, direction_only, active_parameters=active)

    assert expected.status == actual.status == "evaluated"
    actual_evidence = dict(actual.evidence)
    assert actual_evidence["declared_direction_rank"] == 4
    assert actual_evidence["admissible_nullity"] == 4
    for expected_phase, actual_phase in zip(expected.phases, actual.phases, strict=True):
        assert actual_phase.mole_fractions == pytest.approx(
            expected_phase.mole_fractions, rel=2.0e-7, abs=2.0e-10
        )
    for expected_row, actual_row in zip(expected.rows, actual.rows, strict=True):
        assert actual_row.jacobian == pytest.approx(expected_row.jacobian, rel=2.0e-6, abs=2.0e-8)


def test_finite_phase_total_parameter_jacobian_matches_independent_resolves() -> None:
    """Falsify a partial-only LLE Jacobian or fixed-pressure derivative leakage."""

    parameters, model, problem = _fixture()
    # Keep the finite-difference stencil inside the certified two-finite-phase
    # basin.  At 0.2 the corrected scale-invariant chart is on a local basin
    # boundary, so the neighboring solves truthfully collapse one phase.
    center_value = 0.205
    step = 1.0e-5

    def solve(value: float) -> equilibrium.GeneralReactiveEquilibriumResult:
        result = equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters,
                (_PARAMETER_IDENTITY,),
                values=(value,),
            ),
        )
        assert result.status == "evaluated", result.failure
        return result

    center = solve(center_value)
    below = solve(center_value - step)
    above = solve(center_value + step)
    for row, low, high in zip(center.rows, below.values, above.values, strict=True):
        assert row.jacobian is not None and low is not None and high is not None
        numerical = (high - low) / (2.0 * step)
        scale = max(1.0, abs(row.jacobian[0]), abs(numerical))
        assert abs(row.jacobian[0] - numerical) / scale < 3.0e-6


@pytest.mark.slow
def test_finite_phase_exact_nlp_derivatives_match_centered_oracles() -> None:
    """Falsify wrong objective, constraint, or Lagrangian derivatives supplied to IPOPT."""

    parameters, model, problem = _fixture()
    active = epcsaft.ActiveParameterSet(parameters, (_PARAMETER_IDENTITY,))
    native = _core._solve_grepe(
        model._native, model._native, problem, active, "", True
    ).equilibrium

    assert native.accepted
    assert native.derivative_audits_performed
    assert native.derivative_audit_status == "passed"
    assert native.objective_gradient_audit_max_relative_error < 2.0e-6
    assert native.constraint_jacobian_audit_max_relative_error < 2.0e-8
    assert native.lagrangian_hessian_audit_max_relative_error < 2.0e-5
    local_width = len(native.local_variables)
    assert local_width == len(native.local_lower_bounds) == len(native.local_upper_bounds)
    assert native.kkt_dimension == len(native.kkt_root_values)
    assert len(native.kkt_root_jacobian) == native.kkt_dimension * local_width
    assert len(native.kkt_input_jacobian) == native.kkt_dimension * len(native.input_order)
    assert len(native.output_local_derivatives) == len(problem.outputs) * local_width
    assert max(abs(value) for value in native.kkt_root_values) <= 1.0e-7
    assert native.minimum_bound_margin > 1.0e-7
    assert native.reduced_hessian_threshold == pytest.approx(
        1.0e-10 * native.reduced_hessian_scale, rel=2.0e-15
    )
    total_feed = math.fsum(problem.reaction_system.feed_amounts_mol)
    conserved_mass = math.fsum(
        amount * molar_mass
        for amount, molar_mass in zip(
            problem.reaction_system.feed_amounts_mol,
            problem.reaction_system.molar_masses_kg_per_mol,
            strict=True,
        )
    )
    phase_width = len(model.component_ids) + 1
    for phase in range(len(problem.phases)):
        for species, molar_mass in enumerate(problem.reaction_system.molar_masses_kg_per_mol):
            assert native.local_upper_bounds[phase * phase_width + species] == pytest.approx(
                math.log(2.0 * conserved_mass / (molar_mass * total_feed))
            )
            assert native.phase_amounts_mol[phase][species] <= (conserved_mass / molar_mass)


def test_finite_phase_provider_domain_failure_preserves_failure_axis() -> None:
    """Falsify collapsing a Provider source-domain rejection into nonconvergence."""

    _, model, problem = _fixture()
    first_phase = problem.phases[0]
    assert isinstance(first_phase.model, equilibrium.ProviderModel)
    outside_source_domain = replace(
        problem,
        phases=(
            replace(
                first_phase,
                model=replace(
                    first_phase.model,
                    admissible_packing_fraction_interval=(1.0e-8, 1.0e-7),
                ),
            ),
            problem.phases[1],
        ),
    )
    result = equilibrium.solve(model, outside_source_domain)

    assert result.status == "non_evaluable"
    assert result.failure is not None
    assert result.failure.code == "provider_domain_rejection"
    assert result.solver_status == "not_started"
    assert result.numerical_status == "not_adjudicated"
    assert result.physical_status == "not_adjudicated"
    assert result.provider_domain_status == "failed"


def test_finite_phase_amount_chart_is_scale_invariant_with_trace_species() -> None:
    """Preserve phase composition while scaling extensive feed and amounts."""

    _, model, problem = _fixture()
    trace_floor = 0.3
    base_amount_scale = 10.0
    scale = 1.0e6

    def configured(scale_factor: float) -> equilibrium.GeneralReactiveEquilibriumProblem:
        phases = []
        for phase in problem.phases:
            assert isinstance(phase.start, equilibrium.FinitePhaseStart)
            amounts = tuple(
                amount * base_amount_scale * scale_factor for amount in phase.start.amounts_mol
            )
            phases.append(
                replace(
                    phase,
                    start=equilibrium.FinitePhaseStart(
                        amounts, phase.start.molar_volume_m3_per_mol
                    ),
                )
            )
        feed = tuple(
            math.fsum(phase.start.amounts_mol[component] for phase in phases)
            for component in range(3)
        )
        reaction = replace(
            problem.reaction_system,
            feed_amounts_mol=feed,
            conserved_totals=(feed[0] + feed[1], feed[2]),
            strict_interior_amount_floor_mol=trace_floor,
        )
        return replace(
            problem,
            identity=f"{problem.identity}-scale-{scale_factor:g}",
            phases=tuple(phases),
            reaction_system=reaction,
        )

    reference = equilibrium.solve(model, configured(1.0))
    scaled = equilibrium.solve(model, configured(scale))
    assert reference.status == "evaluated", reference.failure
    assert scaled.status == "evaluated", scaled.failure
    assert reference.physical_status == scaled.physical_status == "passed"
    assert reference.evidence and scaled.evidence

    for reference_phase, scaled_phase in zip(reference.phases, scaled.phases, strict=True):
        assert scaled_phase.mole_fractions == pytest.approx(
            reference_phase.mole_fractions, rel=2.0e-7, abs=2.0e-12
        )
        assert scaled_phase.molar_volume_m3_per_mol == pytest.approx(
            reference_phase.molar_volume_m3_per_mol, rel=2.0e-7
        )
    for reference_row, scaled_row in zip(reference.rows, scaled.rows, strict=True):
        if reference_row.selector == "phase.amount":
            assert scaled_row.value == pytest.approx(
                float(reference_row.value) * scale, rel=2.0e-6, abs=2.0e-10
            )


def test_finite_phase_strict_mapping_is_closed() -> None:
    _, _, problem = _fixture()
    mapping = equilibrium.general_reactive_equilibrium_problem_to_mapping(problem)
    loaded = equilibrium.general_reactive_equilibrium_problem_from_mapping(mapping)
    assert loaded == problem
    resolved_problem = replace(
        problem,
        feed=equilibrium.direct_amounts(
            problem.reaction_system.species_ids,
            problem.reaction_system.feed_amounts_mol,
        ),
    )
    resolved_mapping = equilibrium.general_reactive_equilibrium_problem_to_mapping(resolved_problem)
    assert (
        equilibrium.general_reactive_equilibrium_problem_from_mapping(resolved_mapping)
        == resolved_problem
    )

    phase_direction_problem = replace(
        resolved_problem,
        phase_reactions=equilibrium.PhaseReactionSystem(
            ("aqueous", "organic"),
            problem.reaction_system.species_ids,
            tuple(
                equilibrium.PhaseReactionRow(
                    f"transfer-{species_id}",
                    (
                        equilibrium.PhaseSpeciesTerm("aqueous", species_id, -1.0),
                        equilibrium.PhaseSpeciesTerm("organic", species_id, 1.0),
                    ),
                )
                for species_id in problem.reaction_system.species_ids
            ),
            tuple(float(value) for value in problem.reaction_system.charges),
            problem.reaction_system.balance_matrix,
        ),
    )
    phase_direction_mapping = equilibrium.general_reactive_equilibrium_problem_to_mapping(
        phase_direction_problem
    )
    assert (
        equilibrium.general_reactive_equilibrium_problem_from_mapping(phase_direction_mapping)
        == phase_direction_problem
    )

    mapping["phases"][0]["start"]["amounts_mol"][0] = -1.0
    with pytest.raises(ValueError, match="nonnegative amounts"):
        equilibrium.general_reactive_equilibrium_problem_from_mapping(mapping)


def test_unsupported_finite_solved_pressure_keeps_typed_refusal() -> None:
    _, model, problem = _fixture()
    result = equilibrium.solve(
        model,
        replace(
            problem,
            pressure=equilibrium.Solved(
                1.0e5 * epcsaft.unit_registry.pascal,
                (
                    5.0e4 * epcsaft.unit_registry.pascal,
                    2.0e5 * epcsaft.unit_registry.pascal,
                ),
                (1.0e5 * epcsaft.unit_registry.pascal,),
            ),
        ),
    )

    assert result.status == "non_evaluable"
    assert result.failure is not None
    assert result.failure.code == "unsupported_topology"
    assert result.solver_status == "not_started"
    assert result.numerical_status == "unavailable"


def _ascani_2023_problem() -> tuple[
    epcsaft.Mixture,
    equilibrium.GeneralReactiveEquilibriumProblem,
    tuple[float, ...],
]:
    """Configure Ascani's hypothetical Figure 3 through generic GREPE inputs."""

    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(ASCANI_2023_HYPOTHETICAL))
    temperature = 300.0 * epcsaft.unit_registry.kelvin
    pressure = 1.0e5 * epcsaft.unit_registry.pascal
    pure_liquid_log_fugacity = []
    pure_liquid_chemical_potential = []
    for component in range(3):
        composition = tuple(float(index == component) for index in range(3))
        state = model.state(T=temperature, P=pressure, x=composition, phase="liquid")
        assert state.fugacity is not None
        pure_liquid_log_fugacity.append(state.fugacity.ln_coefficient[component])
        trace = 1.0e-11
        trace_composition = tuple(
            1.0 - 2.0 * trace if index == component else trace for index in range(3)
        )
        pure_result = equilibrium.solve(
            model,
            equilibrium.PhaseEquilibriumProblem(
                equilibrium.FixedTPFeed(
                    temperature,
                    pressure,
                    trace_composition,
                    equilibrium.DeclaredTopology(
                        (equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),)
                    ),
                )
            ),
        )
        pure_liquid_chemical_potential.append(
            pure_result.phases[0].chemical_potential_over_rt[component]
        )
    source_ln_k = math.log(2.25)
    provider_ln_k = (
        source_ln_k
        - pure_liquid_chemical_potential[0]
        - pure_liquid_chemical_potential[1]
        + pure_liquid_chemical_potential[2]
    )

    # Figure 3 publishes the two endpoints graphically, not as a numerical table.
    # These digitized tie-line-I coordinates are starts only; every asserted result
    # below is independently certified by the configured total-Gibbs solve.
    starts = (
        (0.496692, 0.202820, 0.300488),
        (0.255113, 0.094788, 0.650099),
    )
    first_a_total = starts[0][0] + starts[0][2]
    second_a_total = starts[1][0] + starts[1][2]
    first_b_total = starts[0][1] + starts[0][2]
    second_b_total = starts[1][1] + starts[1][2]
    determinant = first_a_total * second_b_total - second_a_total * first_b_total
    phase_amounts = (
        (0.7 * second_b_total - 0.5 * second_a_total) / determinant,
        (0.5 * first_a_total - 0.7 * first_b_total) / determinant,
    )
    phases = tuple(
        equilibrium.ReactivePhase(
            identity,
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-8, 0.74), "installed-provider-eos"),
            start=equilibrium.FinitePhaseStart(
                tuple(amount * fraction for fraction in composition),
                1.0
                / float(
                    model.state(
                        T=temperature,
                        P=pressure,
                        x=composition,
                        phase="liquid",
                    )
                    .molar_density.to("mole / meter**3")
                    .magnitude
                ),
            ),
        )
        for identity, composition, amount in zip(
            ("a-rich-liquid", "c-rich-liquid"), starts, phase_amounts, strict=True
        )
    )
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=model.component_ids,
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(0.05, 0.05, 0.10),
        balance_matrix=((1.0, 0.0, 1.0), (0.0, 1.0, 1.0)),
        conserved_totals=(0.7, 0.5),
        reaction_matrix=((-1.0, -1.0, 1.0),),
        feed_amounts_mol=(0.5, 0.3, 0.2),
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                provider_ln_k,
                "ascani-2023-eq10-pure-liquid-standard-converted",
                "provider-helmholtz-coordinate-basis",
                "products_positive",
                "already-provider-basis",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-12,
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"{phase.identity}-{component_id}-mole-fraction",
            "phase.mole_fraction",
            "dimensionless",
            "mole-fraction",
            phase.identity,
            tuple(float(index == component) for index in range(3)),
            support="positive",
        )
        for phase in phases
        for component, component_id in enumerate(model.component_ids)
    )
    return (
        model,
        equilibrium.GeneralReactiveEquilibriumProblem(
            "ascani-2023-hypothetical-figure-3",
            equilibrium.Fixed(temperature),
            equilibrium.Fixed(pressure),
            phases,
            reaction,
            (phases[0].identity,),
            outputs,
        ),
        tuple(pure_liquid_log_fugacity),
    )


@pytest.mark.slow
def test_generic_two_phase_grepe_solves_ascani_2023_hypothetical_case() -> None:
    """Prove the generic two-phase solve on Ascani's source-configured system."""

    model, problem, pure_liquid_log_fugacity = _ascani_2023_problem()
    result = equilibrium.solve(model, problem)

    assert result.status == "evaluated", result.failure
    assert result.numerical_status == "passed"
    assert result.physical_status == "passed"
    assert result.provider_domain_status == "passed"
    assert len(result.phases) == 2
    evidence = dict(result.evidence)
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-7
    assert evidence["transfer_chemical_potential_inf_norm"] <= 1.0e-7
    assert evidence["phase_distinctness"] >= 0.35
    digitized_figure_3_tie_line = (
        (0.496692, 0.202820, 0.300488),
        (0.255113, 0.094788, 0.650099),
    )
    for phase, digitized in zip(result.phases, digitized_figure_3_tie_line, strict=True):
        assert phase.mole_fractions == pytest.approx(digitized, abs=0.12)
    for phase in result.phases:
        state = model.state(
            T=300.0 * epcsaft.unit_registry.kelvin,
            P=1.0e5 * epcsaft.unit_registry.pascal,
            x=phase.mole_fractions,
            phase="liquid",
        )
        assert state.fugacity is not None
        source_activity_reaction = (
            -math.log(phase.mole_fractions[0])
            - state.fugacity.ln_coefficient[0]
            + pure_liquid_log_fugacity[0]
            - math.log(phase.mole_fractions[1])
            - state.fugacity.ln_coefficient[1]
            + pure_liquid_log_fugacity[1]
            + math.log(phase.mole_fractions[2])
            + state.fugacity.ln_coefficient[2]
            - pure_liquid_log_fugacity[2]
        )
        assert source_activity_reaction == pytest.approx(math.log(2.25), abs=1.0e-7)


def test_generic_finite_phase_grepe_accepts_an_empty_reaction_set() -> None:
    """Prove that configured phases do not require a fictitious reaction."""

    model, reactive_problem, _ = _ascani_2023_problem()
    feed = tuple(
        math.fsum(
            phase.start.amounts_mol[component]
            for phase in reactive_problem.phases
            if phase.start is not None
        )
        for component in range(len(model.component_ids))
    )
    nonreactive = equilibrium.ChemicalEquilibriumProblem(
        species_ids=model.component_ids,
        charges=(0, 0, 0),
        molar_masses_kg_per_mol=(0.05, 0.05, 0.10),
        balance_matrix=tuple(
            tuple(float(row == column) for column in range(3)) for row in range(3)
        ),
        conserved_totals=feed,
        reaction_matrix=(),
        feed_amounts_mol=feed,
        equilibrium_constants=(),
        strict_interior_amount_floor_mol=1.0e-12,
    )
    problem = replace(
        reactive_problem,
        identity="ascani-2023-nonreactive-configured-phases",
        reaction_system=nonreactive,
        reaction_phase_ids=(),
        outputs=reactive_problem.outputs
        + tuple(
            equilibrium.EquilibriumOutput(
                f"{phase.identity}-{component_id}-amount",
                "phase.amount",
                "mole",
                "amount",
                phase.identity,
                tuple(float(index == component) for index in range(3)),
                support="positive",
            )
            for phase in reactive_problem.phases
            for component, component_id in enumerate(model.component_ids)
        ),
    )

    result = equilibrium.solve(model, problem)

    assert result.status == "evaluated", result.failure
    evidence = dict(result.evidence)
    assert evidence["reaction_affinity_inf_norm"] == 0.0
    assert evidence["transfer_chemical_potential_inf_norm"] <= 1.0e-7
    phase_component_amounts = tuple(row.value for row in result.rows[-6:])
    assert all(amount is not None for amount in phase_component_amounts)
    for component, expected in enumerate(feed):
        assert math.fsum(
            float(phase_component_amounts[phase * 3 + component]) for phase in range(2)
        ) == pytest.approx(expected, abs=2.0e-10)


@pytest.mark.slow
def test_generic_two_phase_grepe_generates_ascani_start_candidates() -> None:
    """Prove known phases do not require source-derived composition starts."""

    model, problem, _ = _ascani_2023_problem()
    generated_start_problem = replace(
        problem,
        identity="ascani-2023-hypothetical-generated-starts",
        phases=tuple(replace(phase, start=None) for phase in problem.phases),
    )
    result = equilibrium.solve(model, generated_start_problem)

    assert result.status == "evaluated", result.failure
    evidence = dict(result.evidence)
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-7
    assert evidence["transfer_chemical_potential_inf_norm"] <= 1.0e-7
    expected = (
        (0.5444273618587753, 0.21016731227664495, 0.2454053258645798),
        (0.17715297890104614, 0.06383171565566206, 0.7590153054432918),
    )
    actual = tuple(
        phase.mole_fractions
        for phase in sorted(result.phases, key=lambda phase: phase.mole_fractions[0], reverse=True)
    )
    for phase, reference in zip(actual, expected, strict=True):
        assert phase == pytest.approx(reference, abs=2.0e-9)


@pytest.mark.slow
def test_configured_charged_two_phase_reactive_minimum() -> None:
    """Prove charged two-phase GREPE with source-state parameter derivatives."""

    parameter_mapping = deepcopy(CHARGED_TERNARY_IMMISCIBLE)
    parameter_mapping["domains"] = deepcopy(MULTINEUTRAL_REACTIVE["domains"])
    parameters = epcsaft.Parameters.from_mapping(parameter_mapping)
    model = epcsaft.Mixture(parameters)
    seed_feed = (0.32, 0.32, 0.32, 0.02, 0.02)
    seed = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTPFeed(
                300.0 * epcsaft.unit_registry.kelvin,
                1.0e5 * epcsaft.unit_registry.pascal,
                seed_feed,
                equilibrium.DeclaredTopology(
                    tuple(
                        equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents())
                        for _ in range(2)
                    )
                ),
            )
        ),
    )
    assert isinstance(seed.extent, equilibrium.FeedSplit)
    assert len(seed.phases) == 2
    starts = tuple(phase.mole_fractions for phase in seed.phases)
    feed = tuple(
        sum(phase.mole_fractions[species] for phase in seed.phases)
        for species in range(len(model.component_ids))
    )
    ln_k = (
        seed.phases[0].chemical_potential_over_rt[1] - seed.phases[0].chemical_potential_over_rt[0]
    )
    standard_state = equilibrium.ChemicalStandardState(
        id="manufactured-charged-two-phase-source",
        activity_scale_id="mole-fraction-to-molality",
        log_activity_scale_factors=(0.0, 0.0, 0.0, 0.0, 0.0),
        reference_pressure_pa=1.0e5,
        source_reference_component_ids=model.component_ids,
        source_reference_solvent_composition=(1.0, 0.0, 0.0, 0.0, 0.0),
        source_reference_ion_pairs=(epcsaft.SourceReferenceIonPair(("cation", "anion"), (1, 1)),),
        source_reference_phase="liquid",
        source_reference_convention="molality-infinite-dilution",
        source_reference_activity_convention_id="molality-infinite-dilution-v1",
        source_reference_standard_molality_mol_per_kg=1.0,
    )
    transfer = epcsaft.source_reference_transfer(
        model,
        epcsaft.SourceReferenceDeclaration(
            reference_state_id=standard_state.id,
            component_ids=model.component_ids,
            solvent_mole_fractions=standard_state.source_reference_solvent_composition,
            ion_pairs=standard_state.source_reference_ion_pairs,
            phase="liquid",
            reference_convention="molality-infinite-dilution",
            activity_convention_id=standard_state.source_reference_activity_convention_id,
            standard_molality_mol_per_kg=1.0,
            reference_pressure_pa=1.0e5,
            required_derivatives=(),
        ),
        T=300.0 * epcsaft.unit_registry.kelvin,
        P=1.0e5 * epcsaft.unit_registry.pascal,
    )
    standard_state = replace(
        standard_state,
        log_activity_scale_factors=(
            transfer.activity_scale_log_contractions[0],
            transfer.activity_scale_log_contractions[1],
            transfer.activity_scale_log_contractions[2],
            0.5 * transfer.activity_scale_log_contractions[3],
            0.5 * transfer.activity_scale_log_contractions[3],
        ),
    )
    source_ln_k = ln_k - (
        transfer.transfer_log_contractions[1] - transfer.transfer_log_contractions[0]
    )
    reaction = equilibrium.ChemicalEquilibriumProblem(
        species_ids=model.component_ids,
        charges=(0, 0, 0, 1, -1),
        molar_masses_kg_per_mol=(0.03, 0.03, 0.03, 0.01, 0.01),
        balance_matrix=(
            (1.0, 1.0, 0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, 1.0, 0.0),
        ),
        conserved_totals=(feed[0] + feed[1], feed[2], feed[3]),
        reaction_matrix=((-1.0, 1.0, 0.0, 0.0, 0.0),),
        feed_amounts_mol=feed,
        equilibrium_constants=(
            equilibrium.ChemicalEquilibriumConstant(
                source_ln_k,
                "manufactured-charged-two-phase-oracle",
                standard_state.id,
                "products_positive",
                "source-standard-state-to-provider-neutral-reference",
                True,
            ),
        ),
        strict_interior_amount_floor_mol=1.0e-14,
        source_standard_state=standard_state,
    )
    phases = tuple(
        equilibrium.ReactivePhase(
            f"liquid-{index}",
            "liquid",
            "finite",
            equilibrium.AllComponents(),
            equilibrium.ProviderModel((1.0e-8, 0.74), "installed-provider-eos"),
            start=equilibrium.FinitePhaseStart(amounts, state.molar_volume_m3_per_mol),
        )
        for index, (amounts, state) in enumerate(zip(starts, seed.phases, strict=True))
    )
    outputs = tuple(
        equilibrium.EquilibriumOutput(
            f"{phase.identity}-alpha",
            "phase.mole_fraction",
            "dimensionless",
            _BASIS,
            phase.identity,
            (1.0, 0.0, 0.0, 0.0, 0.0),
            support="positive",
        )
        for phase in phases
    )
    problem = equilibrium.GeneralReactiveEquilibriumProblem(
        "configured-charged-two-phase",
        equilibrium.Fixed(300.0 * epcsaft.unit_registry.kelvin),
        equilibrium.Fixed(1.0e5 * epcsaft.unit_registry.pascal),
        phases,
        reaction,
        (phases[0].identity,),
        outputs,
    )

    result = equilibrium.solve(model, problem)

    assert result.status == "evaluated", result.failure
    assert len(result.phases) == 2
    for phase in result.phases:
        assert abs(phase.mole_fractions[3] - phase.mole_fractions[4]) <= 2.0e-10
    evidence = dict(result.evidence)
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-7
    assert evidence["transfer_chemical_potential_inf_norm"] <= 1.0e-7

    generated = equilibrium.solve(
        model,
        replace(
            problem,
            identity="configured-charged-two-phase-generated-starts",
            phases=tuple(replace(phase, start=None) for phase in phases),
        ),
    )
    assert generated.status == "evaluated", generated.failure
    assert len(generated.phases) == 2
    for phase in generated.phases:
        assert abs(phase.mole_fractions[3] - phase.mole_fractions[4]) <= 2.0e-10
    generated_evidence = dict(generated.evidence)
    assert generated_evidence["reaction_affinity_inf_norm"] <= 1.0e-7
    assert generated_evidence["transfer_chemical_potential_inf_norm"] <= 1.0e-7
    assert generated_evidence["phase_distinctness"] > 0.1

    step = 1.0e-7

    def active_solve(value: float) -> equilibrium.GeneralReactiveEquilibriumResult:
        active_result = equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters,
                (_PARAMETER_IDENTITY,),
                values=(value,),
            ),
        )
        assert active_result.status == "evaluated", active_result.failure
        return active_result

    center = active_solve(0.2)
    below = active_solve(0.2 - step)
    above = active_solve(0.2 + step)
    for row, low, high in zip(center.rows, below.values, above.values, strict=True):
        assert row.jacobian is not None and low is not None and high is not None
        numerical = (high - low) / (2.0 * step)
        scale = max(1.0, abs(row.jacobian[0]), abs(numerical))
        assert abs(row.jacobian[0] - numerical) / scale < 2.0e-5


def test_generic_regression_recovers_finite_phase_parameter(tmp_path: Path) -> None:
    """Falsify any hidden LLE-specific Regression objective or owner receipt."""

    parameters, model, problem = _fixture()
    problem = _explicit_phase_direction_problem(problem)
    total_feed = math.fsum(problem.reaction_system.feed_amounts_mol)
    problem = replace(
        problem,
        feed=equilibrium.mole_fractions(
            problem.reaction_system.species_ids,
            tuple(amount / total_feed for amount in problem.reaction_system.feed_amounts_mol),
            total_feed,
        ),
    )
    true_value = 0.202
    truth = equilibrium.solve(
        model,
        problem,
        active_parameters=epcsaft.ActiveParameterSet(
            parameters,
            (_PARAMETER_IDENTITY,),
            values=(true_value,),
        ),
    )
    assert truth.status == "evaluated", truth.failure
    observed = truth.values[0]
    assert observed is not None
    parameter_payload = json.dumps(
        parameters.to_mapping(), separators=(",", ":"), sort_keys=True
    ).encode()
    parameter_path = tmp_path / "reactive-lle-parameters.json"
    parameter_path.write_bytes(parameter_payload)
    fit_mapping = {
        "schema_version": 2,
        "parameters": {
            "path": str(parameter_path),
            "sha256": f"sha256:{hashlib.sha256(parameter_payload).hexdigest()}",
        },
        "model": {"kind": "mixture"},
        "active": {
            "identities": [_PARAMETER_IDENTITY],
            "coordinates": {
                _PARAMETER_IDENTITY: {
                    "bounds": [0.19, 0.21],
                    "origin": 0.2,
                    "scale": 0.005,
                }
            },
            "ties": [],
        },
        "observations": [
            {
                "family": "equilibrium",
                "identity": "manufactured/reactive-lle-regression",
                "request": equilibrium.general_reactive_equilibrium_problem_to_mapping(problem),
                "targets": [
                    {
                        "identity": "manufactured/reactive-lle-alpha-amount",
                        "prediction_identity": problem.outputs[0].identity,
                        "observed": observed,
                        "unit": "mole",
                        "residual": "scaled_difference",
                        "scale": observed,
                        "multiplier": 1.0,
                        "role": "diagnostic",
                        "measurement_classification": "manufactured",
                        "basis": _BASIS,
                        "source_identity": "manufactured-reactive-lle-oracle",
                        "source_hash": f"sha256:{'3' * 64}",
                        "included": True,
                        "drop_reason": None,
                        "aggregate_identity": None,
                        "covariance_identity": None,
                    }
                ],
                "multiplier": 1.0,
            }
        ],
        "solver": {
            "multistart": {
                "starts": [{"identity": "offset-start", "values": [0.198]}],
                "affine_bound_margin": 1.0e-7,
            },
            "controls": {
                "optimizer_tolerance": 1.0e-12,
            },
        },
    }
    fitted = regression.fit(fit_mapping, base_path=tmp_path)
    accepted = fitted.best_usable_start
    assert accepted is not None, fitted.to_mapping()
    assert accepted.final_physical == pytest.approx((true_value,), abs=2.0e-8)
    assert accepted.predictions == pytest.approx((observed,), rel=2.0e-9)
    assert len(accepted.owner_receipts) == 1
    assert accepted.owner_receipts[0].digest.startswith("sha256:")
