from __future__ import annotations

import math

import pytest

from epcsaft import Mixture, Parameters, equilibrium
from epcsaft import unit_registry as u
from tests.eos._cases import NEUTRAL_BINARY, NEUTRAL_TERNARY_IMMISCIBLE


def test_declared_one_phase_certifies_pressure_feed_stability_and_role() -> None:
    feed = (1.0, 1.0)
    result = equilibrium.solve(
        Mixture(Parameters.from_mapping(NEUTRAL_BINARY)),
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTPFeed(
                300.0 * u.kelvin,
                1.0e5 * u.pascal,
                feed,
                equilibrium.DeclaredTopology(
                    (equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),)
                ),
            )
        ),
    )

    assert isinstance(result.extent, equilibrium.FeedSplit)
    assert result.extent.phase_amounts_mol == pytest.approx((math.fsum(feed),), abs=1.0e-9)
    phase = result.phases[0]
    assert tuple(
        result.extent.phase_amounts_mol[0] * x for x in phase.mole_fractions
    ) == pytest.approx(feed, abs=1.0e-9)
    assert phase.role == "vapor"
    assert phase.pressure_pa == pytest.approx(1.0e5, rel=1.0e-8)
    assert phase.mechanical_class == "strict_stable"


def test_declared_hint_reorders_grepe_generated_starts_without_discovery() -> None:
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))
    state = model.state(
        T=300.0 * u.kelvin,
        P=1.0e5 * u.pascal,
        x=(0.5, 0.5),
        phase="vapor",
    )
    result = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            equilibrium.FixedTPFeed(
                300.0 * u.kelvin,
                1.0e5 * u.pascal,
                (1.0, 1.0),
                equilibrium.DeclaredTopology(
                    (equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),)
                ),
            ),
            hints=(
                equilibrium.PhaseHint(
                    0,
                    (0.5, 0.5),
                    float(state.molar_density.to("mole / meter**3").magnitude),
                ),
            ),
        ),
    )

    assert result.diagnostics.route == "grepe_declared_topology"
    assert result.diagnostics.hint_count == 1
    assert result.diagnostics.hint_reordered_start_count == 1


def test_declared_multiphase_certifies_balance_equilibrium_curvature_and_permutation() -> None:
    parameters = Parameters.from_mapping(NEUTRAL_TERNARY_IMMISCIBLE)
    model = Mixture(parameters)
    starts = ((0.90, 0.06, 0.04), (0.11, 0.83, 0.06), (0.14, 0.12, 0.74))
    start_fractions = (0.26, 0.33, 0.41)

    def initial_guesses(
        order: tuple[int, ...],
    ) -> tuple[equilibrium.DeclaredPhaseInitialGuess, ...]:
        guesses = []
        for phase_index, start_index in enumerate(order):
            composition = starts[start_index]
            state = model.state(
                T=300.0 * u.kelvin,
                P=1.0e5 * u.pascal,
                x=composition,
                phase="liquid",
            )
            guesses.append(
                equilibrium.DeclaredPhaseInitialGuess(
                    phase_index,
                    composition,
                    float(state.molar_density.to("mole / meter**3").magnitude),
                    start_fractions[start_index],
                )
            )
        return tuple(guesses)

    specification = equilibrium.FixedTPFeed(
        300.0 * u.kelvin,
        1.0e5 * u.pascal,
        (1.0, 1.0, 1.0),
        equilibrium.DeclaredTopology(
            tuple(
                equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents())
                for _ in range(3)
            )
        ),
    )
    forward = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            specification, initial_guesses=initial_guesses((0, 1, 2))
        ),
    )
    permuted = equilibrium.solve(
        model,
        equilibrium.PhaseEquilibriumProblem(
            specification, initial_guesses=initial_guesses((2, 0, 1))
        ),
    )

    for result in (forward, permuted):
        assert isinstance(result.extent, equilibrium.FeedSplit)
        assert len(result.phases) == 3
        assert min(result.extent.phase_amounts_mol) > 0.0
        assert result.diagnostics.material_balance_inf <= 1.0e-9
        assert result.diagnostics.electroneutrality_inf <= 1.0e-9
        assert result.diagnostics.pressure_residual_inf <= 1.0e-8
        assert result.diagnostics.chemical_potential_residual_inf <= 1.0e-6
        assert result.diagnostics.kkt_stationarity_inf <= 1.0e-7
        assert result.diagnostics.phase_distinctness >= 1.0e-3
        for component in range(3):
            reconstructed = math.fsum(
                fraction * phase.mole_fractions[component]
                for fraction, phase in zip(
                    result.extent.phase_fractions, result.phases, strict=True
                )
            )
            assert reconstructed == pytest.approx(1.0 / 3.0, abs=1.0e-9)
        for phase in result.phases:
            assert phase.mechanical_class == "strict_stable"
            assert phase.pressure_pa == pytest.approx(1.0e5, rel=1.0e-8)
            assert min(phase.mole_fractions) > 0.03
        for phase in result.phases[1:]:
            assert phase.chemical_potential_over_rt == pytest.approx(
                result.phases[0].chemical_potential_over_rt, abs=1.0e-6
            )

    def physical_set(
        result: equilibrium.PhaseEquilibriumResult,
    ) -> list[tuple[float, ...]]:
        return sorted(
            (phase.mole_fractions for phase in result.phases),
            key=lambda values: max(range(3), key=values.__getitem__),
        )

    for actual, expected in zip(physical_set(permuted), physical_set(forward), strict=True):
        assert actual == pytest.approx(expected, abs=2.0e-10)
