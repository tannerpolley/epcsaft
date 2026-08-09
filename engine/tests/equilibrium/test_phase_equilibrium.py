from __future__ import annotations

from collections.abc import Sequence

import pytest

from epcsaft import ActiveParameterSet, Mixture, Parameters
from epcsaft import unit_registry as u
from epcsaft.equilibrium import (
    PhaseEquilibriumError,
    PhaseEquilibriumProblem,
    PhaseHint,
    PhaseState,
    phase_equilibrium,
)
from tests.eos._raw_cases import (
    AQUEOUS_ELECTROLYTE,
    MULTISOLVENT_MULTIION,
    NEUTRAL_ASSOCIATING_MIXTURE,
    NEUTRAL_BINARY,
    _component,
    _document,
    _pair,
)

_NEUTRAL_TEMPERATURE = 243.58 * u.kelvin
_NEUTRAL_PRESSURE = 3.949e6 * u.pascal
_NEUTRAL_FEED = (0.48815, 0.51185)
_THREE_LIQUID = _document(
    "symmetric-three-liquid-fixture",
    [
        _component(
            component,
            molar_mass=0.03,
            segment_count=1.5,
            segment_diameter=3.5,
            dispersion_energy_over_k=250.0,
        )
        for component in ("a", "b", "c")
    ],
    pairs=[_pair(left, right, k_ij=0.4) for left, right in (("a", "b"), ("a", "c"), ("b", "c"))],
)


def _outputs(phases: Sequence[PhaseState]) -> tuple[float, ...]:
    return tuple(value for phase in phases for value in (*phase.amounts_mol, phase.volume_m3))


def _assert_balances(
    phases: Sequence[PhaseState],
    feed: tuple[float, ...],
    charges: tuple[int, ...],
) -> None:
    reconstructed = tuple(
        sum(phase.amounts_mol[component] for phase in phases) for component in range(len(feed))
    )
    assert reconstructed == pytest.approx(feed, abs=2.0e-9)
    for phase in phases:
        assert sum(phase.mole_fractions) == pytest.approx(1.0, abs=2.0e-12)
        assert sum(
            charge * amount for charge, amount in zip(charges, phase.amounts_mol, strict=True)
        ) == pytest.approx(0.0, abs=2.0e-9)


def test_discovery_hints_roles_and_order_reach_the_same_neutral_state() -> None:
    """Falsify discovery scope, role, hint, order, and fail-closed claims.

    Oracle: methane/ethane one-vapor and two-phase Gibbs minima plus independently
    permuted constrained/component representations. Amount and phase-state
    invariants use 2e-9 absolute and 2e-10 relative tolerances; pressure,
    chemical-potential, and KKT residuals must remain within the solver's
    declared certificate limits. A wrong one-liquid declaration at a dilute
    vapor state must fail closed with its best candidate. Failure means one-phase
    discovery is absent, hints changed the finite mathematical search, phase
    labels changed the physical result, component positions selected a route,
    or an uncertified topology escaped the public boundary.
    """
    parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    model = Mixture(parameters)
    problem = PhaseEquilibriumProblem(
        _NEUTRAL_TEMPERATURE,
        _NEUTRAL_PRESSURE,
        _NEUTRAL_FEED,
    )
    discovered = phase_equilibrium(model, problem)
    assert tuple(phase.role for phase in discovered.phases) == ("liquid", "vapor")
    assert discovered.diagnostics.searched_topologies == (
        "L",
        "V",
        "LL",
        "LV",
        "LLL",
        "LLV",
    )
    assert discovered.diagnostics.topology_count == discovered.diagnostics.topology_budget == 6
    assert discovered.diagnostics.globality_claim == "not_established"
    assert discovered.diagnostics.phase_selection_claim == (
        "lowest_certified_local_candidate_in_recorded_finite_search"
    )
    assert discovered.diagnostics.completeness_criterion
    assert set(discovered.diagnostics.excluded_topologies) == {
        "more_than_one_vapor_phase",
        "solid_phases",
        "reactive_multiphase_equilibrium",
        "more_than_three_phases",
    }
    _assert_balances(discovered.phases, _NEUTRAL_FEED, (0, 0))

    hinted_problem = PhaseEquilibriumProblem(
        _NEUTRAL_TEMPERATURE,
        _NEUTRAL_PRESSURE,
        _NEUTRAL_FEED,
        hints=tuple(
            PhaseHint(phase.role, phase.mole_fractions, phase.molar_density_mol_m3)
            for phase in reversed(discovered.phases)
        ),
    )
    hinted = phase_equilibrium(model, hinted_problem)
    assert hinted.diagnostics.hint_count == 2
    assert hinted.diagnostics.hint_reordered_start_count > 0
    assert _outputs(hinted.phases) == pytest.approx(
        _outputs(discovered.phases), rel=2.0e-10, abs=2.0e-10
    )
    assert hinted.diagnostics.phase_selection_claim == discovered.diagnostics.phase_selection_claim

    constrained = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(
            _NEUTRAL_TEMPERATURE,
            _NEUTRAL_PRESSURE,
            _NEUTRAL_FEED,
            mode="constrained",
            phase_roles=("vapor", "liquid"),
        ),
    )
    assert constrained.diagnostics.certificate_scope == "declared_topology_only"
    assert constrained.diagnostics.globality_claim == "not_established"
    assert tuple(phase.role for phase in constrained.phases) == ("liquid", "vapor")
    assert _outputs(constrained.phases) == pytest.approx(
        _outputs(discovered.phases), rel=2.0e-10, abs=2.0e-10
    )

    reversed_parameters = Parameters.from_mapping(NEUTRAL_BINARY, components=("ethane", "methane"))
    reordered = phase_equilibrium(
        Mixture(reversed_parameters),
        PhaseEquilibriumProblem(
            _NEUTRAL_TEMPERATURE,
            _NEUTRAL_PRESSURE,
            tuple(reversed(_NEUTRAL_FEED)),
            mode="constrained",
            phase_roles=("liquid", "vapor"),
        ),
    )
    for canonical, permuted in zip(discovered.phases, reordered.phases, strict=True):
        assert permuted.amounts_mol == pytest.approx(
            tuple(reversed(canonical.amounts_mol)), rel=2.0e-10, abs=2.0e-10
        )
        assert permuted.volume_m3 == pytest.approx(canonical.volume_m3, rel=2.0e-10)

    with pytest.raises(PhaseEquilibriumError) as rejected:
        phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                300.0 * u.kelvin,
                1.0e5 * u.pascal,
                (0.5, 0.5),
                mode="constrained",
                phase_roles=("liquid",),
            ),
        )
    assert rejected.value.diagnostics.failure_stage == "certificate"
    assert rejected.value.diagnostics.certified_candidate_count == 0
    assert len(rejected.value.best_candidate) == 1

    one_phase = phase_equilibrium(
        model,
        PhaseEquilibriumProblem(
            300.0 * u.kelvin,
            1.0e5 * u.pascal,
            (0.5, 0.5),
        ),
    )
    assert tuple(phase.role for phase in one_phase.phases) == ("vapor",)


def test_three_liquid_and_multisolvent_multiion_charts_are_certified() -> None:
    """Falsify real three-phase and general exact-electroneutral charts.

    Oracle: a symmetric repulsive ternary PC-SAFT system with three distinct
    strict liquid minima; a six-component two-solvent/two-cation/two-anion
    topology evaluated in canonical and permuted order; charged discovery; and
    an associating neutral binary. Conservation and charge use 2e-9 mol absolute
    tolerance; pressure and KKT checks use the solver's reported limits. Invalid
    model temperature and ion-fraction inputs must fail before any start. Failure
    means discovery cannot select three phases or a charged phase, the owner
    assumes one solvent/ion pair, drops neutral association, or launches a solve
    outside the EOS domain.
    """
    three_model = Mixture(Parameters.from_mapping(_THREE_LIQUID))
    three = phase_equilibrium(
        three_model,
        PhaseEquilibriumProblem(
            250.0 * u.kelvin,
            1.0e5 * u.pascal,
            (1.0 / 3.0,) * 3,
            mode="constrained",
            phase_roles=("liquid",) * 3,
        ),
    )
    assert len(three.phases) == 3
    _assert_balances(three.phases, (1.0 / 3.0,) * 3, (0, 0, 0))
    assert all(max(phase.mole_fractions) > 0.99 for phase in three.phases)
    assert (
        len({phase.mole_fractions.index(max(phase.mole_fractions)) for phase in three.phases}) == 3
    )
    assert three.diagnostics.minimum_reduced_hessian_eigenvalue is not None
    assert three.diagnostics.minimum_reduced_hessian_eigenvalue > 0.0

    discovered_three = phase_equilibrium(
        three_model,
        PhaseEquilibriumProblem(
            250.0 * u.kelvin,
            1.0e6 * u.pascal,
            (1.0 / 3.0,) * 3,
        ),
    )
    assert tuple(phase.role for phase in discovered_three.phases) == ("liquid",) * 3
    assert all(max(phase.mole_fractions) > 0.99 for phase in discovered_three.phases)
    _assert_balances(discovered_three.phases, (1.0 / 3.0,) * 3, (0, 0, 0))

    canonical_ids = (
        "solvent",
        "cosolvent",
        "cation",
        "anion",
        "divalent-cation",
        "iodide-anion",
    )
    canonical_feed = (0.60, 0.33, 0.02, 0.02, 0.01, 0.02)
    feed_by_id = dict(zip(canonical_ids, canonical_feed, strict=True))
    permuted_ids = (
        "iodide-anion",
        "cosolvent",
        "divalent-cation",
        "solvent",
        "anion",
        "cation",
    )
    results = []
    for component_ids in (canonical_ids, permuted_ids):
        model = Mixture(Parameters.from_mapping(MULTISOLVENT_MULTIION, components=component_ids))
        feed = tuple(feed_by_id[component] for component in component_ids)
        pressure = model.state(
            T=298.15 * u.kelvin,
            rho=2000.0 * u.mole / u.meter**3,
            x=feed,
        ).pressure
        result = phase_equilibrium(
            model,
            PhaseEquilibriumProblem(
                298.15 * u.kelvin,
                pressure,
                feed,
                mode="constrained",
                phase_roles=("liquid",),
            ),
        )
        charges = tuple(
            {"cation": 1, "anion": -1, "divalent-cation": 2, "iodide-anion": -1}.get(component, 0)
            for component in component_ids
        )
        _assert_balances(result.phases, feed, charges)
        assert result.diagnostics.electroneutrality_inf <= (
            result.diagnostics.electroneutrality_tolerance_mol
        )
        assert "charge_equivalents" in " ".join(result.diagnostics.equation_state_identities)
        results.append(result)

    canonical_by_id = dict(zip(canonical_ids, results[0].phases[0].mole_fractions, strict=True))
    assert results[1].phases[0].mole_fractions == pytest.approx(
        tuple(canonical_by_id[component] for component in permuted_ids),
        rel=2.0e-10,
        abs=2.0e-10,
    )

    associating = phase_equilibrium(
        Mixture(Parameters.from_mapping(NEUTRAL_ASSOCIATING_MIXTURE)),
        PhaseEquilibriumProblem(
            300.0 * u.kelvin,
            1.0e5 * u.pascal,
            (0.6, 0.4),
            mode="constrained",
            phase_roles=("vapor",),
        ),
    )
    _assert_balances(associating.phases, (0.6, 0.4), (0, 0))

    aqueous_model = Mixture(Parameters.from_mapping(AQUEOUS_ELECTROLYTE))
    discovered_charged = phase_equilibrium(
        aqueous_model,
        PhaseEquilibriumProblem(
            298.15 * u.kelvin,
            1.0e5 * u.pascal,
            (0.98, 0.01, 0.01),
        ),
    )
    assert tuple(phase.role for phase in discovered_charged.phases) == ("liquid",)
    _assert_balances(discovered_charged.phases, (0.98, 0.01, 0.01), (0, 1, -1))
    assert discovered_charged.diagnostics.kkt_stationarity_inf <= (
        discovered_charged.diagnostics.kkt_stationarity_tolerance
    )

    for invalid_problem in (
        PhaseEquilibriumProblem(
            300.0 * u.kelvin,
            1.0e5 * u.pascal,
            (0.98, 0.01, 0.01),
        ),
        PhaseEquilibriumProblem(
            298.15 * u.kelvin,
            1.0e5 * u.pascal,
            (0.5, 0.25, 0.25),
        ),
    ):
        with pytest.raises(PhaseEquilibriumError) as domain_rejected:
            phase_equilibrium(
                aqueous_model,
                invalid_problem,
            )
        assert domain_rejected.value.diagnostics.failure_stage == "input"
        assert domain_rejected.value.diagnostics.attempted_start_count == 0


@pytest.mark.parametrize(
    ("mapping", "feed", "temperature", "pressure", "roles", "identities", "steps"),
    (
        (
            NEUTRAL_BINARY,
            _NEUTRAL_FEED,
            243.58,
            3.949e6,
            ("liquid", "vapor"),
            ("pair/ethane/methane/k_ij",),
            (1.0e-4,),
        ),
        (
            AQUEOUS_ELECTROLYTE,
            (0.98, 0.01, 0.01),
            298.15,
            1.0e5,
            ("liquid",),
            (
                "association/water/a/water/b/energy_over_k",
                "model/relative_permittivity/ion_fraction_suppression_coefficient",
            ),
            (1.0e-3, 1.0e-4),
        ),
    ),
    ids=("neutral-two-phase", "charged-one-phase"),
)
def test_public_kkt_sensitivity_matches_independent_centered_resolves(
    mapping: dict[str, object],
    feed: tuple[float, ...],
    temperature: float,
    pressure: float,
    roles: tuple[str, ...],
    identities: tuple[str, ...],
    steps: tuple[float, ...],
) -> None:
    """Falsify public KKT F_y/F_p layout and selected-parameter totals.

    Oracle: independent centered re-solves of the same fixed topology with each
    Parameters-owned coordinate perturbed at the explicit step. Amount and
    volume derivatives use 3e-6 relative plus 3e-11 absolute tolerance, which
    covers double-precision nonlinear solve error and centered O(h^2)
    truncation. The implicit backward error must be below 1e-13. Failure means
    the general owner changed active identity/order, omitted charged-chart or
    packing geometry terms, or differentiated solver history instead of the
    declared KKT equations.
    """
    parameters = Parameters.from_mapping(mapping)
    selected = ActiveParameterSet(parameters, identities)
    problem = PhaseEquilibriumProblem(
        temperature * u.kelvin,
        pressure * u.pascal,
        feed,
        mode="constrained",
        phase_roles=roles,  # type: ignore[arg-type]
    )

    def solve(values: tuple[float, ...]) -> tuple[float, ...]:
        active = ActiveParameterSet(parameters, identities, values=values)
        return _outputs(
            phase_equilibrium(Mixture(parameters), problem, active_parameters=active).phases
        )

    result = phase_equilibrium(Mixture(parameters), problem, active_parameters=selected)
    sensitivity = result.sensitivity
    assert sensitivity is not None
    assert sensitivity.status == "available"
    assert sensitivity.parameters == selected.coordinates
    assert sensitivity.rank == sensitivity.dimension
    assert sensitivity.linear_solve_backward_error is not None
    assert sensitivity.linear_solve_backward_error <= 1.0e-13
    exact = sensitivity.output_parameter_jacobian
    assert len(exact) == len(_outputs(result.phases))

    for parameter, step in enumerate(steps):
        lower = list(selected.values)
        upper = list(selected.values)
        lower[parameter] -= step
        upper[parameter] += step
        lower_output = solve(tuple(lower))
        upper_output = solve(tuple(upper))
        centered = tuple(
            (upper_value - lower_value) / (2.0 * step)
            for lower_value, upper_value in zip(lower_output, upper_output, strict=True)
        )
        assert tuple(row[parameter] for row in exact) == pytest.approx(
            centered, rel=3.0e-6, abs=3.0e-11
        )
