"""Native transport, phase/reaction order and unavailable-state causes."""
import pytest
from epcsaft.equilibrium import Amounts, IdealMoleFraction, Phase, Problem, Reaction, reaction_activities, solve_equilibrium


def test_reaction_activities() -> None:
    model = IdealMoleFraction(["A", "B", "C"], [0, 0, 0])
    problem = Problem(T=300, P=100000, phases=[Phase("L", kind="liquid")], feed=Amounts({"A": 1}),
                      reactions=[Reaction({"A": -1, "B": 1}, equilibrium_constant=2)])
    result = solve_equilibrium(model, problem)
    row, = reaction_activities(model, result)
    assert (row.phase_index, row.reaction_index, row.basis) == (0, 0, "ideal-mole-fraction")
    assert row.species.activities == result.mole_fractions
    assert row.species.activities[2] == 0 and row.species.log_activities[2] == -float("inf")
    assert row.species.causes == ["", "", ""] and row.log_quotient is not None and not row.cause
    with pytest.raises(ValueError, match="model/result mismatch"):
        reaction_activities(IdealMoleFraction(["A", "B", "C"], [0, 0, 0]), result)
