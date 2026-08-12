from __future__ import annotations

import pytest

from epcsaft import ActiveParameterSet, Mixture, Parameters
from epcsaft import unit_registry as u
from epcsaft.equilibrium import (
    AllComponents,
    DeclaredTopology,
    DiscoverPhases,
    EquilibriumPhaseState,
    FeedSplit,
    FixedTCoexistence,
    FixedTPFeed,
    IncipientBoundary,
    NeutralComponents,
    PhaseEquilibriumError,
    PhaseEquilibriumProblem,
    PhaseEquilibriumResult,
    PhaseHint,
    PhaseRole,
    PhaseSpecification,
    phase_equilibrium,
    tp_flash,
)
from tests.eos._raw_cases import (
    AQUEOUS_ELECTROLYTE,
    CHARGED_TERNARY_IMMISCIBLE,
    MULTINEUTRAL_REACTIVE,
    MULTISOLVENT_MULTIION,
    NEUTRAL_BINARY,
    NEUTRAL_TERNARY_IMMISCIBLE,
)

_TEMPERATURE = 300.0 * u.kelvin
_PRESSURE = 1.0e5 * u.pascal
_FEED = (0.5, 0.5)
_TOTAL_AMOUNT_MOL = 2.0
_AMOUNT_FEED = tuple(_TOTAL_AMOUNT_MOL * value for value in _FEED)


def test_discovery_preserves_protected_held_state_and_parameter_order() -> None:
    """Falsify state or exact selected-parameter drift from protected HELD discovery.

    Oracle: the retained public tp_flash result at identical Parameters,
    temperature, pressure, and normalized feed. Phase intensive states use
    2e-12 relative/absolute agreement; the amount-feed sensitivity carries the
    exact incumbent parameter identity/order and its established extensive
    coordinate scaling. Failure means the new typed discovery seam changed the
    protected numerical workflow.
    """
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    active = ActiveParameterSet(
        parameters,
        (
            "pair/ethane/methane/k_ij",
            "component/methane/dispersion_energy_over_k",
        ),
    )
    model = Mixture(parameters)
    retained = tp_flash(
        model,
        _TEMPERATURE,
        _PRESSURE,
        _FEED,
        active_parameters=active,
    )
    problem = PhaseEquilibriumProblem(
        FixedTPFeed(
            _TEMPERATURE,
            _PRESSURE,
            _AMOUNT_FEED,
            DiscoverPhases(),
        )
    )
    discovered = phase_equilibrium(model, problem, active_parameters=active)

    assert discovered.diagnostics.route == "held_discovery"
    assert discovered.diagnostics.globality_claim == "not_established"
    assert discovered.diagnostics.provider_evaluation_count > 0
    assert isinstance(discovered.extent, FeedSplit)
    assert discovered.extent.total_free_energy_over_rt == pytest.approx(
        _TOTAL_AMOUNT_MOL * retained.total_free_energy_over_rt,
        rel=2.0e-12,
        abs=2.0e-12,
    )
    assert discovered.extent.phase_amounts_mol == pytest.approx(
        tuple(_TOTAL_AMOUNT_MOL * value for value in retained.phase_fractions),
        rel=2.0e-12,
        abs=2.0e-12,
    )
    for phase, expected in zip(discovered.phases, retained.phases, strict=True):
        assert phase.mole_fractions == pytest.approx(
            expected.mole_fractions, rel=2.0e-12, abs=2.0e-12
        )
        assert phase.molar_volume_m3_per_mol == pytest.approx(
            expected.volume_m3 / expected.amount_mol,
            rel=2.0e-12,
            abs=2.0e-12,
        )
        assert phase.pressure_pa == pytest.approx(
            expected.pressure_pa, rel=2.0e-12, abs=2.0e-12
        )

    selected = discovered.sensitivity
    incumbent = retained.sensitivity
    assert selected is not None
    assert incumbent is not None
    assert selected.parameters == incumbent.parameters
    assert (
        selected.residual_parameter_jacobian
        == incumbent.residual_parameter_jacobian
    )
    for selected_row, incumbent_row in zip(
        selected.residual_state_jacobian,
        incumbent.residual_state_jacobian,
        strict=True,
    ):
        assert selected_row == pytest.approx(
            tuple(value / _TOTAL_AMOUNT_MOL for value in incumbent_row),
            rel=2.0e-12,
            abs=2.0e-12,
        )


def test_declared_single_phase_solves_role_and_hints_do_not_change_root() -> None:
    """Falsify admission-only role handling or hint-dependent scientific results.

    Oracle: the protected discovery state identifies the physical one-phase
    role. The declared HELD pressure-envelope solve must independently recover
    that intensive root, certify a strict-stable mechanical branch, reject the
    opposite role, and return the same state for two distinct finite hints.
    The 2e-8 volume tolerance covers pressure-root bisection, not a changed
    thermodynamic claim.
    """
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    model = Mixture(parameters)
    discovered = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(
            FixedTPFeed(
                _TEMPERATURE,
                _PRESSURE,
                _AMOUNT_FEED,
                DiscoverPhases(),
            )
        ),
    )
    assert len(discovered.phases) == 1
    role = discovered.phases[0].role
    topology = DeclaredTopology(
        (PhaseSpecification(role, AllComponents()),)
    )
    specifications = (
        PhaseEquilibriumProblem(
            FixedTPFeed(
                _TEMPERATURE,
                _PRESSURE,
                _AMOUNT_FEED,
                topology,
            ),
            (
                PhaseHint(
                    0,
                    discovered.phases[0].mole_fractions,
                    discovered.phases[0].molar_density_mol_m3,
                ),
            ),
        ),
        PhaseEquilibriumProblem(
            FixedTPFeed(
                _TEMPERATURE,
                _PRESSURE,
                _AMOUNT_FEED,
                topology,
            ),
            (
                PhaseHint(
                    0,
                    tuple(reversed(discovered.phases[0].mole_fractions)),
                    1.2 * discovered.phases[0].molar_density_mol_m3,
                ),
            ),
        ),
    )
    solved = tuple(
        phase_equilibrium(model, specification)
        for specification in specifications
    )
    for result in solved:
        assert result.diagnostics.certificate_scope == "fixed_tp_declared_topology"
        assert result.phases[0].mechanical_class == "strict_stable"
        assert result.phases[0].role == role
        assert result.phases[0].molar_volume_m3_per_mol == pytest.approx(
            discovered.phases[0].molar_volume_m3_per_mol,
            rel=2.0e-8,
            abs=2.0e-12,
        )
    assert solved[0].phases == solved[1].phases

    opposite = "liquid" if role == "vapor" else "vapor"
    with pytest.raises(PhaseEquilibriumError) as rejected:
        phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                FixedTPFeed(
                    _TEMPERATURE,
                    _PRESSURE,
                    _AMOUNT_FEED,
                    DeclaredTopology(
                        (PhaseSpecification(opposite, AllComponents()),)
                    ),
                )
            ),
        )
    assert rejected.value.diagnostics.failure_stage == "phase_role_certificate"


def test_declared_three_liquid_phases_satisfy_independent_equilibrium_claims() -> None:
    """Falsify a length-only three-phase implementation.

    Oracle: an interior symmetric ternary split has three positive phase
    extents, reconstructs the equimolar feed, and independently satisfies
    common pressure and every ordinary chemical-potential equality. The
    fixture keeps every minority mole fraction well away from the coordinate
    boundary, so equality rather than active-bound complementarity is the
    applicable thermodynamic claim.
    """
    parameters = Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    identity = "pair/alpha/beta/k_ij"
    selected = ActiveParameterSet(parameters, (identity,))
    model = Mixture(parameters)
    specification = FixedTPFeed(
        300.0 * u.kelvin,
        _PRESSURE,
        (1.0, 1.0, 1.0),
        DeclaredTopology(
            tuple(
                PhaseSpecification("liquid", AllComponents())
                for _ in range(3)
            )
        ),
    )
    result = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(specification),
        active_parameters=selected,
    )

    assert isinstance(result.extent, FeedSplit)
    assert len(result.phases) == 3
    assert result.diagnostics.route == "held_step8_declared_topology"
    assert result.diagnostics.material_balance_inf <= 2.0e-10
    assert result.diagnostics.pressure_residual_inf <= 1.0e-8
    assert result.diagnostics.chemical_potential_residual_inf <= 2.0e-8
    assert result.diagnostics.kkt_stationarity_inf <= 1.0e-7
    assert result.diagnostics.phase_distinctness >= 1.0e-3
    assert min(result.extent.phase_fractions) > 0.0
    for component in range(3):
        reconstructed = sum(
            fraction * phase.mole_fractions[component]
            for fraction, phase in zip(
                result.extent.phase_fractions, result.phases, strict=True
            )
        )
        assert reconstructed == pytest.approx(1.0 / 3.0, abs=2.0e-10)
    for phase in result.phases:
        assert phase.role == "liquid"
        assert phase.mechanical_class == "strict_stable"
        assert min(phase.mole_fractions) > 0.03
        assert phase.pressure_pa == pytest.approx(1.0e5, rel=1.0e-8)
    for phase in result.phases[1:]:
        assert phase.chemical_potential_over_rt == pytest.approx(
            result.phases[0].chemical_potential_over_rt,
            rel=2.0e-8,
            abs=2.0e-8,
        )

    hinted = PhaseEquilibriumProblem(
        specification,
        tuple(
            PhaseHint(
                index,
                phase.mole_fractions,
                phase.molar_density_mol_m3,
            )
            for index, phase in enumerate(result.phases)
        ),
    )
    perturbed = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(
            specification,
            tuple(
                PhaseHint(
                    index,
                    tuple(reversed(phase.mole_fractions)),
                    (1.2 if index % 2 == 0 else 0.8)
                    * phase.molar_density_mol_m3,
                )
                for index, phase in enumerate(result.phases)
            ),
        ),
        active_parameters=selected,
    )
    expected_set = sorted(
        (phase.mole_fractions for phase in result.phases),
        key=lambda values: max(range(3), key=values.__getitem__),
    )
    perturbed_set = sorted(
        (phase.mole_fractions for phase in perturbed.phases),
        key=lambda values: max(range(3), key=values.__getitem__),
    )
    for candidate, expected in zip(perturbed_set, expected_set, strict=True):
        assert candidate == pytest.approx(expected, abs=2.0e-12)

    def solve(value: float) -> PhaseEquilibriumResult:
        return phase_equilibrium(
            Mixture(parameters),
            hinted,
            active_parameters=ActiveParameterSet(
                parameters, (identity,), values=(value,)
            ),
        )

    def outputs(candidate: PhaseEquilibriumResult) -> tuple[float, ...]:
        assert isinstance(candidate.extent, FeedSplit)
        values: list[float] = []
        for amount, phase in zip(
            candidate.extent.phase_amounts_mol,
            candidate.phases,
            strict=True,
        ):
            values.extend(amount * fraction for fraction in phase.mole_fractions)
            values.append(amount * phase.molar_volume_m3_per_mol)
        return tuple(values)

    step = 1.0e-5
    center = solve(selected.values[0])
    lower = outputs(solve(selected.values[0] - step))
    upper = outputs(solve(selected.values[0] + step))
    sensitivity = center.sensitivity
    assert sensitivity is not None
    assert sensitivity.status == "available"
    exact = tuple(row[0] for row in sensitivity.output_parameter_jacobian)
    centered = tuple(
        (right - left) / (2.0 * step)
        for left, right in zip(lower, upper, strict=True)
    )
    assert exact == pytest.approx(centered, rel=2.0e-6, abs=3.0e-10)


def test_declared_multisolvent_multiion_roles_select_distinct_stable_roots() -> None:
    """Falsify binary-only or compressibility-only charged role selection.

    The six-component feed is exactly electroneutral and admits distinct
    mechanically strict-stable liquid and vapor pressure roots at the same
    fixed T,P. Each declared one-phase solve must select its requested root
    through the native pressure envelope while retaining every component.
    """
    feed = (0.60, 0.33, 0.02, 0.02, 0.01, 0.02)
    model = Mixture(Parameters.from_mapping(MULTISOLVENT_MULTIION))

    def solve(
        candidate: Mixture,
        candidate_feed: tuple[float, ...],
        role: PhaseRole,
    ) -> EquilibriumPhaseState:
        result = phase_equilibrium(
            candidate,
            PhaseEquilibriumProblem(
                FixedTPFeed(
                    298.15 * u.kelvin,
                    _PRESSURE,
                    candidate_feed,
                    DeclaredTopology(
                        (PhaseSpecification(role, AllComponents()),)
                    ),
                )
            ),
        )
        assert result.diagnostics.pressure_residual_inf <= 1.0e-8
        return result.phases[0]

    liquid = solve(model, feed, "liquid")
    vapor = solve(model, feed, "vapor")
    assert liquid.mechanical_class == "strict_stable"
    assert vapor.mechanical_class == "strict_stable"
    assert liquid.mole_fractions == pytest.approx(feed, abs=2.0e-12)
    assert vapor.mole_fractions == pytest.approx(feed, abs=2.0e-12)
    assert liquid.molar_density_mol_m3 > 100.0 * vapor.molar_density_mol_m3
    assert abs(
        liquid.mole_fractions[2]
        - liquid.mole_fractions[3]
        + 2.0 * liquid.mole_fractions[4]
        - liquid.mole_fractions[5]
    ) <= 1.0e-12

    feed_by_id = dict(zip(model.component_ids, feed, strict=True))
    permuted_ids = tuple(reversed(model.component_ids))
    permuted = Mixture(
        Parameters.from_mapping(
            MULTISOLVENT_MULTIION, components=permuted_ids
        )
    )
    permuted_liquid = solve(
        permuted,
        tuple(feed_by_id[component] for component in permuted_ids),
        "liquid",
    )
    permuted_composition = dict(
        zip(permuted_ids, permuted_liquid.mole_fractions, strict=True)
    )
    assert tuple(
        permuted_composition[component] for component in model.component_ids
    ) == pytest.approx(liquid.mole_fractions, abs=2.0e-12)
    assert permuted_liquid.molar_volume_m3_per_mol == pytest.approx(
        liquid.molar_volume_m3_per_mol, rel=2.0e-11
    )

    def solve_split(roles: tuple[PhaseRole, PhaseRole]) -> PhaseEquilibriumResult:
        return phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                FixedTPFeed(
                    298.15 * u.kelvin,
                    _PRESSURE,
                    feed,
                    DeclaredTopology(
                        tuple(
                            PhaseSpecification(role, AllComponents())
                            for role in roles
                        )
                    ),
                )
            ),
        )

    forward = solve_split(("liquid", "vapor"))
    reversed_order = solve_split(("vapor", "liquid"))
    assert forward.diagnostics.provider_evaluation_count > 0
    assert tuple(phase.role for phase in forward.phases) == ("liquid", "vapor")
    assert tuple(phase.role for phase in reversed_order.phases) == (
        "vapor",
        "liquid",
    )
    for selected_phase, reversed_phase in zip(
        forward.phases,
        reversed(reversed_order.phases),
        strict=True,
    ):
        assert selected_phase.mole_fractions == pytest.approx(
            reversed_phase.mole_fractions, rel=2.0e-11, abs=2.0e-12
        )
        assert selected_phase.molar_volume_m3_per_mol == pytest.approx(
            reversed_phase.molar_volume_m3_per_mol, rel=2.0e-11
        )
    assert isinstance(forward.extent, FeedSplit)
    assert isinstance(reversed_order.extent, FeedSplit)
    assert forward.extent.phase_amounts_mol == pytest.approx(
        tuple(reversed(reversed_order.extent.phase_amounts_mol)),
        rel=2.0e-11,
        abs=2.0e-12,
    )


def test_coexistence_structurally_excludes_ions_but_keeps_every_neutral() -> None:
    """Falsify a water-specific or epsilon-ion vapor support implementation.

    The model contains two neutral solvents and a balanced ion pair. The
    partner phase must be built on exactly the two-neutral support, satisfy
    both admitted chemical-potential equalities and common pressure, and
    retain a mechanically stable liquid--vapor volume ordering. No published
    value is used as an oracle; this compact case isolates the structural
    support and defining-equation claim.
    """
    model = Mixture(Parameters.from_mapping(MULTINEUTRAL_REACTIVE))
    result = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(
            FixedTCoexistence(
                298.15 * u.kelvin,
                DeclaredTopology(
                    (
                        PhaseSpecification("liquid", AllComponents()),
                        PhaseSpecification("vapor", NeutralComponents()),
                    )
                ),
                0,
                (0.70, 0.20, 0.05, 0.05),
            )
        ),
    )

    assert isinstance(result.extent, IncipientBoundary)
    liquid, vapor = result.phases
    assert liquid.supported_component_ids == (
        "solvent",
        "cosolvent",
        "cation",
        "anion",
    )
    assert vapor.supported_component_ids == ("solvent", "cosolvent")
    assert vapor.chemical_potential_component_ids == (
        "solvent",
        "cosolvent",
    )
    assert sum(vapor.mole_fractions) == pytest.approx(1.0, abs=2.0e-12)
    assert liquid.chemical_potential_over_rt[:2] == pytest.approx(
        vapor.chemical_potential_over_rt,
        abs=2.0e-8,
    )
    assert result.diagnostics.chemical_potential_residual_inf <= 2.0e-8
    assert result.diagnostics.pressure_residual_inf <= 1.0e-8
    assert liquid.mechanical_class == vapor.mechanical_class == "strict_stable"
    assert liquid.molar_volume_m3_per_mol < vapor.molar_volume_m3_per_mol

    reference_by_id = dict(zip(model.component_ids, (0.70, 0.20, 0.05, 0.05), strict=True))
    reversed_ids = tuple(reversed(model.component_ids))
    reversed_result = phase_equilibrium(
        Mixture(
            Parameters.from_mapping(
                MULTINEUTRAL_REACTIVE,
                components=reversed_ids,
            )
        ),
        PhaseEquilibriumProblem(
            FixedTCoexistence(
                298.15 * u.kelvin,
                DeclaredTopology(
                    (
                        PhaseSpecification("liquid", AllComponents()),
                        PhaseSpecification("vapor", NeutralComponents()),
                    )
                ),
                0,
                tuple(reference_by_id[component] for component in reversed_ids),
            )
        ),
    )
    assert isinstance(reversed_result.extent, IncipientBoundary)
    assert reversed_result.extent.pressure_pa == pytest.approx(
        result.extent.pressure_pa,
        rel=2.0e-10,
    )
    for candidate, expected in zip(
        reversed_result.phases,
        result.phases,
        strict=True,
    ):
        candidate_by_id = dict(
            zip(candidate.supported_component_ids, candidate.mole_fractions, strict=True)
        )
        expected_by_id = dict(
            zip(expected.supported_component_ids, expected.mole_fractions, strict=True)
        )
        assert candidate_by_id == pytest.approx(expected_by_id, abs=2.0e-10)


def test_coexistence_model_parameter_sensitivity_uses_both_phase_partials() -> None:
    """Falsify omission of model-wide partner-phase parameter derivatives.

    The dielectric suppression coefficient applies to both all-component
    phases. The implicit boundary-pressure derivative must therefore match a
    centered full-model re-solve; a reference-phase-only derivative has the
    wrong sign for this state.
    """
    parameters = Parameters.from_mapping(AQUEOUS_ELECTROLYTE)
    model = Mixture(parameters)
    identity = "model/relative_permittivity/ion_fraction_suppression_coefficient"
    selected = ActiveParameterSet(parameters, (identity,))
    problem = PhaseEquilibriumProblem(
        FixedTCoexistence(
            298.15 * u.kelvin,
            DeclaredTopology(
                (
                    PhaseSpecification("liquid", AllComponents()),
                    PhaseSpecification("vapor", AllComponents()),
                )
            ),
            0,
            (0.96, 0.02, 0.02),
        )
    )

    def solve(value: float) -> PhaseEquilibriumResult:
        return phase_equilibrium(
            model,
            problem,
            active_parameters=ActiveParameterSet(
                parameters,
                (identity,),
                values=(value,),
            ),
        )

    step = 1.0e-3
    center = solve(selected.values[0])
    lower = solve(selected.values[0] - step)
    upper = solve(selected.values[0] + step)
    assert isinstance(center.extent, IncipientBoundary)
    assert isinstance(lower.extent, IncipientBoundary)
    assert isinstance(upper.extent, IncipientBoundary)
    sensitivity = center.sensitivity
    assert sensitivity is not None
    assert sensitivity.status == "available"
    exact = sensitivity.output_parameter_jacobian[0][0]
    centered = (
        upper.extent.pressure_pa - lower.extent.pressure_pa
    ) / (2.0 * step)
    assert exact == pytest.approx(centered, rel=2.0e-6, abs=2.0e-3)


def test_charged_three_phase_split_preserves_charge_order_and_sensitivity() -> None:
    """Falsify a neutral-only three-phase or fixed-position HELD2 extension.

    Three immiscible neutral solvents share a balanced ion pair. Every returned
    phase must remain electroneutral, reconstruct the feed, satisfy the HELD2
    modified-potential equations, and preserve the same physical phase set
    when component order is reversed. The exact phase-indexed derivative is
    independently checked by centered parameter re-solves.
    """
    parameters = Parameters.from_mapping(CHARGED_TERNARY_IMMISCIBLE)
    identity = "pair/alpha/beta/k_ij"
    selected = ActiveParameterSet(parameters, (identity,))
    feed = (0.32, 0.32, 0.32, 0.02, 0.02)
    topology = DeclaredTopology(
        tuple(PhaseSpecification("liquid", AllComponents()) for _ in range(3))
    )
    specification = FixedTPFeed(
        300.0 * u.kelvin,
        _PRESSURE,
        feed,
        topology,
    )

    def solve(value: float) -> PhaseEquilibriumResult:
        return phase_equilibrium(
            Mixture(parameters),
            PhaseEquilibriumProblem(specification),
            active_parameters=ActiveParameterSet(
                parameters,
                (identity,),
                values=(value,),
            ),
        )

    center = solve(selected.values[0])
    assert isinstance(center.extent, FeedSplit)
    assert len(center.phases) == 3
    assert center.diagnostics.electroneutrality_inf <= 1.0e-9
    assert center.diagnostics.pressure_residual_inf <= 1.0e-8
    assert center.diagnostics.chemical_potential_residual_inf <= 1.0e-6
    assert center.diagnostics.kkt_stationarity_inf <= 1.0e-7
    for component, expected in enumerate(feed):
        reconstructed = sum(
            fraction * phase.mole_fractions[component]
            for fraction, phase in zip(
                center.extent.phase_fractions,
                center.phases,
                strict=True,
            )
        )
        assert reconstructed == pytest.approx(expected, abs=2.0e-10)
    for phase in center.phases:
        assert phase.mechanical_class == "strict_stable"
        assert abs(phase.mole_fractions[3] - phase.mole_fractions[4]) <= 1.0e-12

    def outputs(result: PhaseEquilibriumResult) -> tuple[float, ...]:
        assert isinstance(result.extent, FeedSplit)
        values: list[float] = []
        for amount, phase in zip(
            result.extent.phase_amounts_mol,
            result.phases,
            strict=True,
        ):
            values.extend(amount * fraction for fraction in phase.mole_fractions)
            values.append(amount * phase.molar_volume_m3_per_mol)
        return tuple(values)

    step = 1.0e-5
    lower = outputs(solve(selected.values[0] - step))
    upper = outputs(solve(selected.values[0] + step))
    sensitivity = center.sensitivity
    assert sensitivity is not None
    assert sensitivity.status == "available"
    exact = tuple(row[0] for row in sensitivity.output_parameter_jacobian)
    centered = tuple(
        (right - left) / (2.0 * step)
        for left, right in zip(lower, upper, strict=True)
    )
    assert exact == pytest.approx(centered, rel=3.0e-6, abs=4.0e-10)

    canonical_ids = parameters.component_ids
    reversed_ids = tuple(reversed(canonical_ids))
    feed_by_id = dict(zip(canonical_ids, feed, strict=True))
    reversed_parameters = Parameters.from_mapping(
        CHARGED_TERNARY_IMMISCIBLE,
        components=reversed_ids,
    )
    reversed_result = phase_equilibrium(
        Mixture(reversed_parameters),
        PhaseEquilibriumProblem(
            FixedTPFeed(
                300.0 * u.kelvin,
                _PRESSURE,
                tuple(feed_by_id[component] for component in reversed_ids),
                topology,
            )
        ),
    )
    canonical_by_dominant = {
        canonical_ids[max(range(3), key=phase.mole_fractions.__getitem__)]: dict(
            zip(canonical_ids, phase.mole_fractions, strict=True)
        )
        for phase in center.phases
    }
    reversed_by_dominant = {
        reversed_ids[
            max(
                (reversed_ids.index(name) for name in ("alpha", "beta", "gamma")),
                key=phase.mole_fractions.__getitem__,
            )
        ]: dict(zip(reversed_ids, phase.mole_fractions, strict=True))
        for phase in reversed_result.phases
    }
    assert canonical_by_dominant.keys() == reversed_by_dominant.keys()
    for dominant, expected in canonical_by_dominant.items():
        assert tuple(
            reversed_by_dominant[dominant][component] for component in canonical_ids
        ) == pytest.approx(
            tuple(expected[component] for component in canonical_ids),
            abs=2.0e-9,
        )
