from __future__ import annotations

from dataclasses import replace

import pytest

from epcsaft import Mixture, Parameters, equilibrium
from epcsaft import unit_registry as u
from tests.eos._cases import STRONGLY_ASSOCIATING_BINARY
from tests.equilibrium._cases import ionic_reaction_case


def test_public_solve_rejects_invalid_support_rank_conservation_and_charge() -> None:
    model, problem = ionic_reaction_case()
    reaction = problem.reaction_system
    row = reaction.reaction_matrix[0]
    constant = reaction.equilibrium_constants[0]
    invalid_reactions = (
        replace(
            reaction,
            reaction_matrix=(row, row),
            equilibrium_constants=(constant, constant),
        ),
        replace(reaction, reaction_matrix=((1.0, -1.0, -1.0),)),
        replace(reaction, reaction_matrix=((0.0, -1.0, 1.0),)),
    )
    invalid_problems = tuple(
        replace(problem, reaction_system=candidate) for candidate in invalid_reactions
    ) + (
        replace(
            problem,
            phases=(
                replace(
                    problem.phases[0],
                    support=equilibrium.DeclaredComponents(("water",)),
                ),
            ),
        ),
    )

    for invalid in invalid_problems:
        with pytest.raises(ValueError):
            equilibrium.solve(model, invalid)


def test_public_boundary_rejects_an_impossible_temperature_bracket() -> None:
    model = Mixture(Parameters.from_mapping(STRONGLY_ASSOCIATING_BINARY))
    topology = equilibrium.DeclaredTopology(
        (
            equilibrium.PhaseSpecification("liquid", equilibrium.AllComponents()),
            equilibrium.PhaseSpecification("vapor", equilibrium.AllComponents()),
        )
    )

    with pytest.raises(equilibrium.PhaseEquilibriumError):
        equilibrium.solve(
            model,
            equilibrium.PhaseEquilibriumProblem(
                equilibrium.FixedPCoexistence(
                    1.0e6 * u.pascal,
                    (295.0 * u.kelvin, 305.0 * u.kelvin),
                    topology,
                    0,
                    (0.5, 0.5),
                )
            ),
        )
