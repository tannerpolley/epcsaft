from __future__ import annotations

import math

import pytest

from epcsaft import equilibrium
from tests.equilibrium._cases import ionic_reaction_case


def test_reaction_system_certifies_analytic_equilibrium_charge_and_component_permutation() -> None:
    canonical_model, canonical_problem = ionic_reaction_case()
    canonical = equilibrium.solve(canonical_model, canonical_problem)
    permuted_ids = tuple(reversed(canonical_model.component_ids))
    permuted_model, permuted_problem = ionic_reaction_case(permuted_ids)
    permuted = equilibrium.solve(permuted_model, permuted_problem)

    expected = {"water": 1.0, "sodium-cation": 0.5, "chloride-anion": 0.5}
    charges = {"water": 0, "sodium-cation": 1, "chloride-anion": -1}
    for model, problem, result in (
        (canonical_model, canonical_problem, canonical),
        (permuted_model, permuted_problem, permuted),
    ):
        assert result.status == "evaluated", result.failure
        amounts = dict(zip(model.component_ids, result.values, strict=True))
        assert amounts == pytest.approx(expected, abs=1.0e-9)
        assert (
            abs(math.fsum(charges[identity] * float(amounts[identity]) for identity in amounts))
            <= 1.0e-9
        )
        reaction = problem.reaction_system
        for row, target in zip(reaction.balance_matrix, reaction.conserved_totals, strict=True):
            assert math.fsum(
                coefficient * float(value)
                for coefficient, value in zip(row, result.values, strict=True)
            ) == pytest.approx(target, abs=1.0e-9)
        evidence = dict(result.evidence)
        assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
        assert evidence["kkt_stationarity_inf_norm"] <= 1.0e-7
        assert evidence["declared_direction_rank"] == evidence["admissible_nullity"] == 1
