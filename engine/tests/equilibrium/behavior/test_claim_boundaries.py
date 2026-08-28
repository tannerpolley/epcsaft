from __future__ import annotations

from dataclasses import replace

from epcsaft import equilibrium
from epcsaft import unit_registry as u
from tests.equilibrium._cases import finite_phase_case, ionic_reaction_case


def test_evaluated_result_limits_its_claim_to_the_configured_local_certificate() -> None:
    model, problem = ionic_reaction_case()
    result = equilibrium.solve(model, problem)

    assert result.status == "evaluated"
    assert result.certificate_scope
    assert result.globality_claim == "not_established"


def test_unsupported_solved_state_returns_typed_failure_without_fallback() -> None:
    _, model, problem = finite_phase_case()
    result = equilibrium.solve(
        model,
        replace(
            problem,
            pressure=equilibrium.Solved(
                1.0e5 * u.pascal,
                (5.0e4 * u.pascal, 2.0e5 * u.pascal),
                (1.0e5 * u.pascal,),
            ),
        ),
    )

    assert result.status == "non_evaluable"
    assert result.failure is not None
    assert result.failure.code == "unsupported_topology"
    assert result.phases == ()
    assert all(value is None for value in result.values)
