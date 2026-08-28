from __future__ import annotations

import pytest

from epcsaft import Mixture, Parameters, equilibrium
from epcsaft import unit_registry as u
from tests.eos._cases import MULTINEUTRAL_REACTIVE, STRONGLY_ASSOCIATING_BINARY


def _two_phase_topology() -> equilibrium.DeclaredTopology:
    return equilibrium.DeclaredTopology(
        (
            equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
            equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),
        )
    )


def test_fixed_temperature_boundary_certifies_support_stability_and_component_permutation() -> None:
    reference_by_id = {
        "solvent": 0.70,
        "cosolvent": 0.20,
        "cation": 0.05,
        "anion": 0.05,
    }

    def solve(order: tuple[str, ...]) -> equilibrium.PhaseEquilibriumResult:
        model = Mixture(Parameters.from_mapping(MULTINEUTRAL_REACTIVE, components=order))
        return equilibrium.solve(
            model,
            equilibrium.PhaseEquilibriumProblem(
                equilibrium.FixedTCoexistence(
                    298.15 * u.kelvin,
                    equilibrium.DeclaredTopology(
                        (
                            equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
                            equilibrium.PhaseSpecification(
                                "vapor", equilibrium.NeutralComponents()
                            ),
                        )
                    ),
                    0,
                    tuple(reference_by_id[identity] for identity in order),
                )
            ),
        )

    canonical_ids = tuple(reference_by_id)
    canonical = solve(canonical_ids)
    permuted_ids = tuple(reversed(canonical_ids))
    permuted = solve(permuted_ids)

    for result, order in ((canonical, canonical_ids), (permuted, permuted_ids)):
        assert isinstance(result.extent, equilibrium.IncipientBoundary)
        liquid, vapor = result.phases
        neutral_ids = tuple(identity for identity in order if identity in {"solvent", "cosolvent"})
        assert result.extent.pressure_pa > 0.0
        assert result.diagnostics.pressure_residual_inf <= 1.0e-8
        assert result.diagnostics.chemical_potential_residual_inf <= 1.0e-6
        assert liquid.mechanical_class == vapor.mechanical_class == "strict_stable"
        assert liquid.molar_volume_m3_per_mol < vapor.molar_volume_m3_per_mol
        assert liquid.supported_component_ids == order
        assert vapor.supported_component_ids == neutral_ids
        assert vapor.chemical_potential_component_ids == neutral_ids
        liquid_mu = dict(
            zip(
                liquid.chemical_potential_component_ids,
                liquid.chemical_potential_over_rt,
                strict=True,
            )
        )
        for identity, potential in zip(
            vapor.chemical_potential_component_ids,
            vapor.chemical_potential_over_rt,
            strict=True,
        ):
            assert potential == pytest.approx(liquid_mu[identity], abs=1.0e-6)

    assert permuted.extent.pressure_pa == pytest.approx(canonical.extent.pressure_pa, rel=2.0e-10)
    for actual, expected in zip(permuted.phases, canonical.phases, strict=True):
        actual_by_id = dict(zip(actual.supported_component_ids, actual.mole_fractions, strict=True))
        expected_by_id = dict(
            zip(expected.supported_component_ids, expected.mole_fractions, strict=True)
        )
        assert actual_by_id == pytest.approx(expected_by_id, abs=2.0e-10)


@pytest.mark.slow
def test_fixed_pressure_boundary_certifies_bubble_dew_closure_and_component_permutation() -> None:
    def solve(
        order: tuple[str, ...], reference_index: int, reference: tuple[float, ...]
    ) -> equilibrium.PhaseEquilibriumResult:
        model = Mixture(Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY, components=order))
        return equilibrium.solve(
            model,
            equilibrium.PhaseEquilibriumProblem(
                equilibrium.FixedPCoexistence(
                    1200.0 * u.pascal,
                    (295.0 * u.kelvin, 305.0 * u.kelvin),
                    _two_phase_topology(),
                    reference_index,
                    reference,
                )
            ),
        )

    canonical_ids = ("solute", "solvent")
    bubble = solve(canonical_ids, 0, (0.5, 0.5))
    dew = solve(canonical_ids, 1, (0.01, 0.99))
    permuted_ids = tuple(reversed(canonical_ids))
    permuted = solve(permuted_ids, 0, (0.5, 0.5))

    for result, reference_index in ((bubble, 0), (dew, 1), (permuted, 0)):
        assert isinstance(result.extent, equilibrium.IncipientTemperatureBoundary)
        assert result.extent.reference_phase_index == reference_index
        assert 295.0 < result.extent.temperature_k < 305.0
        assert result.extent.pressure_pa == pytest.approx(1200.0, rel=1.0e-8)
        assert result.diagnostics.pressure_residual_inf <= 1.0e-8
        assert result.diagnostics.chemical_potential_residual_inf <= 1.0e-6
        assert result.phases[1].chemical_potential_over_rt == pytest.approx(
            result.phases[0].chemical_potential_over_rt, abs=1.0e-6
        )

    assert permuted.extent.temperature_k == pytest.approx(bubble.extent.temperature_k, rel=2.0e-10)
    for actual, expected in zip(permuted.phases, bubble.phases, strict=True):
        assert dict(
            zip(actual.supported_component_ids, actual.mole_fractions, strict=True)
        ) == pytest.approx(
            dict(zip(expected.supported_component_ids, expected.mole_fractions, strict=True)),
            abs=2.0e-10,
        )
