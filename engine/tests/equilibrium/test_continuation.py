from epcsaft import Mixture, Parameters, trace
from epcsaft.equilibrium import (
    Amounts,
    Continuation,
    ContinuationStop,
    Free,
    Phase,
    Pinned,
    Problem,
    continue_equilibrium,
    solve_equilibrium,
)

from tests.packets import parameters_json


def mixture() -> Mixture:
    return Mixture(Parameters.from_json(
        parameters_json("gross-2001-methane-ethane", 3),
        components=("methane", "ethane"),
    ))


def bubble(temperature: float) -> Problem:
    return Problem(
        phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
        T=temperature, P=Free(1.0e5), feed=Amounts({"methane": 0.5, "ethane": 0.5}),
    )


def test_bubble_branch_follows_the_cold_solve() -> None:
    model = mixture()
    # A 0.5 K step moves the incipient vapor density beyond the 1% anchored-root window
    # (branch_lost, then step_too_small); the quarter step follows the branch.
    # minimum_step = initial_step = 0.5 forbids halving, so this refusal is the branch-identity safeguard working.
    lost = continue_equilibrium(
        model, bubble(180.0), bubble(181.0),
        Continuation(initial_step=0.5, minimum_step=0.5, max_steps=2, max_retries_per_step=1),
    )
    assert not lost.success and lost.stop == ContinuationStop.step_too_small
    assert lost.cause_stop == ContinuationStop.branch_lost and lost.path_values == [0.0]
    result = continue_equilibrium(
        model, bubble(180.0), bubble(181.0),
        Continuation(initial_step=0.25, minimum_step=0.25, max_steps=4, max_retries_per_step=1),
    )
    assert result.success and result.stop == ContinuationStop.reached_target and result.path_values[-1] == 1.0
    for q, state in zip(result.path_values, result.states):
        cold = solve_equilibrium(model, bubble(180.0 + q))
        assert cold.success and state.requested_tolerance_met is True and state.pressure is not None and cold.pressure is not None
        assert abs(state.pressure / cold.pressure - 1.0) <= 1.0e-8
        assert max(abs(a - b) for a, b in zip(state.mole_fractions, cold.mole_fractions)) <= 1.0e-8


def test_changed_declaration_is_an_invalid_path() -> None:
    model = mixture()
    target = Problem(
        phases=[Phase("V", kind="vapor", amount=Pinned(0.0)), Phase("L", kind="liquid")],
        T=181.0, P=Free(1.0e5), feed=Amounts({"methane": 0.5, "ethane": 0.5}),
    )
    result = continue_equilibrium(
        model, bubble(180.0), target,
        Continuation(initial_step=0.5, minimum_step=0.5, max_steps=2, max_retries_per_step=1),
    )
    assert not result.success and result.stop == ContinuationStop.invalid_path
    assert result.states == [] and result.attempts == 0


def test_branch_rejection_discards_the_optimizer_start() -> None:
    model = mixture()
    point = solve_equilibrium(model, bubble(180.0))
    assert point.success and point.pressure is not None
    pressure = point.pressure
    feed = {
        "methane": 0.5 * (point.mole_fractions[0] + point.mole_fractions[2]),
        "ethane": 0.5 * (point.mole_fractions[1] + point.mole_fractions[3]),
    }

    def flash(pressure: float) -> Problem:
        return Problem(
            phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor")],
            T=180.0, P=pressure, feed=Amounts(feed),
        )

    start = flash(0.5 * pressure)
    path = continue_equilibrium(
        model, start, flash(2.0 * pressure),
        Continuation(initial_step=0.01, minimum_step=0.01, max_steps=2, max_retries_per_step=1),
        start_result=solve_equilibrium(model, start),
    )
    assert path.stop == ContinuationStop.step_too_small and path.cause_stop == ContinuationStop.branch_lost
    assert path.attempted_path == 0.01 and path.reached_path == 0.0
    rejected = path.last_result
    assert rejected is not None and rejected.success and rejected.requested_tolerance_met
    with trace(True) as capture:
        solved = solve_equilibrium(model, flash((0.5 + 1.5 * path.attempted_path) * pressure), rejected)
    assert solved.success and solved.requested_tolerance_met
    assert not any(event.name == "warm_initialized" for event in capture.events)
