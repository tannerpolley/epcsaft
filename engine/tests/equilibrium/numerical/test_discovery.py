from __future__ import annotations

import math

import pytest

from epcsaft import Mixture, Parameters, equilibrium
from epcsaft import unit_registry as u
from tests.eos._cases import AQUEOUS_ELECTROLYTE, NEUTRAL_BINARY


def _reconstructed_feed(
    result: equilibrium.PhaseEquilibriumResult,
) -> tuple[float, ...]:
    assert isinstance(result.extent, equilibrium.FeedSplit)
    return tuple(
        math.fsum(
            amount * phase.mole_fractions[component]
            for amount, phase in zip(result.extent.phase_amounts_mol, result.phases, strict=True)
        )
        for component in range(len(result.extent.feed_amounts_mol))
    )


@pytest.mark.slow
def test_neutral_discovery_certifies_feed_equilibrium_stability_and_scaling() -> None:
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))

    def solve(feed: tuple[float, ...]) -> equilibrium.PhaseEquilibriumResult:
        return equilibrium.solve(
            model,
            equilibrium.PhaseEquilibriumProblem(
                equilibrium.FixedTPFeed(300.0 * u.kelvin, 1.0e5 * u.pascal, feed)
            ),
        )

    unit = solve((0.5, 0.5))
    scaled = solve((1.0, 1.0))

    for result, feed in ((unit, (0.5, 0.5)), (scaled, (1.0, 1.0))):
        assert isinstance(result.extent, equilibrium.FeedSplit)
        assert _reconstructed_feed(result) == pytest.approx(feed, abs=1.0e-9)
        assert result.diagnostics.material_balance_inf <= 1.0e-9
        assert result.diagnostics.pressure_residual_inf <= 1.0e-8
        assert result.diagnostics.chemical_potential_residual_inf <= 1.0e-6
        assert result.diagnostics.kkt_stationarity_inf <= 1.0e-7
        for phase in result.phases:
            assert phase.pressure_pa == pytest.approx(1.0e5, rel=1.0e-8)
            assert phase.mechanical_class == "strict_stable"
        for phase in result.phases[1:]:
            assert phase.chemical_potential_over_rt == pytest.approx(
                result.phases[0].chemical_potential_over_rt, abs=1.0e-6
            )

    assert scaled.extent.total_free_energy_over_rt == pytest.approx(
        2.0 * unit.extent.total_free_energy_over_rt, rel=2.0e-12, abs=2.0e-12
    )
    assert scaled.extent.phase_amounts_mol == pytest.approx(
        tuple(2.0 * amount for amount in unit.extent.phase_amounts_mol),
        rel=2.0e-12,
        abs=2.0e-12,
    )


def test_charged_discovery_certifies_charge_support_and_component_permutation() -> None:
    base_parameters = Parameters.from_mapping(AQUEOUS_ELECTROLYTE)
    canonical_ids = base_parameters.component_ids
    feed_by_id = dict(zip(canonical_ids, (0.98, 0.01, 0.01), strict=True))
    charge_by_id = {"water": 0, "sodium-cation": 1, "chloride-anion": -1}

    def solve(order: tuple[str, ...]) -> equilibrium.PhaseEquilibriumResult:
        model = Mixture(base_parameters.select(order))
        return equilibrium.solve(
            model,
            equilibrium.PhaseEquilibriumProblem(
                equilibrium.FixedTPFeed(
                    298.15 * u.kelvin,
                    1.0e5 * u.pascal,
                    tuple(feed_by_id[identity] for identity in order),
                )
            ),
        )

    canonical = solve(canonical_ids)
    permuted_ids = tuple(reversed(canonical_ids))
    permuted = solve(permuted_ids)

    for result, order in ((canonical, canonical_ids), (permuted, permuted_ids)):
        assert isinstance(result.extent, equilibrium.FeedSplit)
        assert _reconstructed_feed(result) == pytest.approx(
            tuple(feed_by_id[identity] for identity in order), abs=1.0e-9
        )
        assert result.diagnostics.material_balance_inf <= 1.0e-9
        assert result.diagnostics.electroneutrality_inf <= 1.0e-9
        assert result.diagnostics.chemical_potential_residual_inf <= 1.0e-6
        assert result.diagnostics.kkt_stationarity_inf <= 1.0e-7
        for phase in result.phases:
            assert phase.supported_component_ids == order
            assert phase.mechanical_class == "strict_stable"
            assert (
                abs(
                    math.fsum(
                        charge_by_id[identity] * fraction
                        for identity, fraction in zip(
                            phase.supported_component_ids, phase.mole_fractions, strict=True
                        )
                    )
                )
                <= 1.0e-9
            )

    canonical_composition = dict(
        zip(canonical_ids, canonical.phases[0].mole_fractions, strict=True)
    )
    permuted_composition = dict(zip(permuted_ids, permuted.phases[0].mole_fractions, strict=True))
    assert permuted_composition == pytest.approx(canonical_composition, abs=2.0e-10)
