import math
from pathlib import Path

import pytest

from epcsaft import Mixture, ParameterError, Parameters, _core, equilibrium
from epcsaft.equilibrium import (
    Amounts,
    Free,
    Molalities,
    Phase,
    Problem,
    Reaction,
    problem_from_mapping,
    solve_equilibrium,
)

from tests.packets import parameters_json


def mixture() -> Mixture:
    return Mixture(Parameters.from_json(
        parameters_json("gross-2001-methane-ethane", 3),
        components=("methane", "ethane"),
    ))


def problem(feed: Amounts, **phase: object) -> Problem:
    return Problem(phases=[Phase("L", kind="liquid", **phase)], T=300.0, P=1.0e5, feed=feed)  # type: ignore[arg-type]


def test_unknown_species_and_amount_status_are_rejected_by_name() -> None:
    model = mixture()
    with pytest.raises(ParameterError, match="unknown species methnae"):
        problem(Amounts({"methnae": 1.0}))._native(model)
    with pytest.raises(ParameterError, match="unknown amount status fixed"):
        problem_from_mapping({
            "phases": [{"name": "L", "kind": "liquid", "amount": {"status": "fixed"}}],
            "T": 300.0,
            "P": 1.0e5,
            "feed": {"basis": "amounts", "values": {"methane": 1.0}},
        })


def test_physical_arguments_are_required() -> None:
    with pytest.raises(TypeError, match="equilibrium_constant"):
        Reaction({"methane": -1.0})
    with pytest.raises(ParameterError, match="one of equilibrium_constant or affinity"):
        Reaction({"methane": -1.0}, equilibrium_constant=1.0, affinity_j_per_mol=2.0)
    with pytest.raises(TypeError, match="solvent_mass_kg"):
        Molalities({"methane": 1.0}, solvent=("methane",))  # type: ignore[call-arg]


def test_generic_phase_requires_an_explicit_kind() -> None:
    model = mixture()
    feed = Amounts({"methane": 0.5, "ethane": 0.5})
    with pytest.raises(ParameterError, match="phase V requires an explicit kind"):
        Problem(phases=[Phase("V")], T=300.0, P=1.0e5, feed=feed)._native(model)
    named = Problem(phases=[Phase("V", kind="liquid")], T=300.0, P=1.0e5, feed=feed)
    assert named._native(model).phases[0].kind == _core.PhaseKind.Liquid
    # Name-based inference is gone: a phase named V declared liquid solves as the liquid branch.
    result = solve_equilibrium(model, named)
    liquid = solve_equilibrium(model, problem(feed))
    assert result.success and liquid.success
    assert result.molar_densities == liquid.molar_densities


def test_omitted_and_explicit_empty_support_are_distinct() -> None:
    model = mixture()
    feed = Amounts({"methane": 0.5, "ethane": 0.0})
    omitted = problem(feed)._native(model)
    assert omitted.phases[0].support == [] and omitted.phases[0].support_declared is False
    compiled = _core.compile_problem(model._model, omitted)
    assert compiled.removed_species == [1]
    explicit = problem(feed, support=("methane", "ethane"))._native(model)
    assert explicit.phases[0].support == [0, 1] and explicit.phases[0].support_declared is True
    assert _core.compile_problem(model._model, explicit).removed_species == [1]
    empty = problem(feed, support=())._native(model)
    assert empty.phases[0].support == [] and empty.phases[0].support_declared is True
    with pytest.raises(ParameterError, match="empty support in phase L"):
        _core.compile_problem(model._model, empty)


def test_caller_row_tolerance_maps_and_reports_its_outcome() -> None:
    model = mixture()
    feed = Amounts({"methane": 0.5, "ethane": 0.5})
    default = solve_equilibrium(model, problem(feed))
    assert default.success and default.convergence.criterion is not None
    assert default.convergence.criterion.name == "EQUILIBRIUM_ROW_RESIDUAL"
    assert default.requested_tolerance_met is True
    requested = Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1.0e5, feed=feed, row_tolerance=1.0e-8)
    assert requested._native(model).row_tolerance == 1.0e-8
    result = solve_equilibrium(model, requested)
    criterion = result.convergence.criterion
    del result
    assert criterion is not None
    assert criterion.name == "EQUILIBRIUM_REQUESTED_ROW_RESIDUAL" and criterion.value == 1.0e-8
    with pytest.raises(ParameterError, match="invalid row tolerance"):
        solve_equilibrium(model, Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1.0e5, feed=feed, row_tolerance=-1.0))
    mapped = problem_from_mapping({
        "phases": [{"name": "L", "kind": "liquid"}], "T": 300.0, "P": 1.0e5,
        "feed": {"basis": "amounts", "values": {"methane": 1.0}}, "row_tolerance": 1.0e-9,
    })
    assert mapped.row_tolerance == 1.0e-9


def test_all_pinned_rank_deficient_declaration_is_not_run() -> None:
    model = mixture()
    feed = Amounts({"methane": 0.5, "ethane": 0.5})
    pinned = Problem(
        phases=[Phase("L", kind="liquid", amount=equilibrium.Pinned(1.0))], T=Free(300.0), P=1.0e5, feed=feed
    )
    result = solve_equilibrium(model, pinned)
    assert not result.success and result.solver_status == "NOT_RUN"
    assert result.generic_nullity > 0 and result.unknown_count == result.jacobian_rank + result.generic_nullity
    assert result.residuals == [] and math.isnan(result.convergence.residual)
    assert result.requested_tolerance_met is None


def test_malformed_explicit_mapping_fields_are_rejected() -> None:
    def mapped(**changes: object) -> Problem:
        return problem_from_mapping({"phases": [{"name": "L", "kind": "liquid"}], "T": 300.0, "P": 1.0e5,
                                     "feed": {"basis": "amounts", "values": {"methane": 1.0}}, **changes})
    # Omissions and defaults keep their meaning: a free amount, declared topology, a fixed mapped value.
    default = mapped()
    assert isinstance(default.phases, tuple) and isinstance(default.phases[0].amount, Free)
    assert default.topology_declared is True and mapped(T={"value": 310.0}).T == 310.0
    # Each malformed field used to become a free phase, a zero amount, T = 0 K or declared topology.
    for changes, message in (
        ({"phases": [{"name": "L", "kind": "liquid", "amount": 1.0}]}, "phase L amount must be a mapping"),
        ({"phases": [{"name": "L", "kind": "liquid", "amount": {"status": "pinned"}}]}, "pinned amount of phase L requires a value"),
        ({"T": {"free": False}}, "missing T value"),
        ({"topology_declared": "false"}, "topology_declared must be a bool"),
    ):
        with pytest.raises(ParameterError, match=message):
            mapped(**changes)


def test_ipopt_initialization_or_fatal_return_publishes_no_state(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Ipopt reads ipopt.opt from the working directory. An unknown option fails Initialize(); an unavailable
    # linear solver fails OptimizeTNLP() before any iterate, leaving the fallback start, whose classified
    # rows (1.50) meet this deliberately loose criterion. A zero iteration budget is an eligible
    # non-error stop: its polished start meets the rows and still publishes.
    model = mixture()
    loose = Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1.0e5,
                    feed=Amounts({"methane": 0.5, "ethane": 0.5}), row_tolerance=2.0)
    monkeypatch.chdir(tmp_path)
    for option, failure in (("bogus_option 1", "Ipopt Initialize returned -12"),
                            ("linear_solver ma57", "Ipopt OptimizeTNLP returned"), ("max_iter 0", None)):
        (tmp_path / "ipopt.opt").write_text(option + "\n")
        result = solve_equilibrium(model, loose)
        if failure is None:
            assert result.success and result.solver_status == "MAXITER_EXCEEDED" and math.isclose(result.amounts[0], 1.0)
            continue
        assert not result.success and result.amounts == [] and result.mole_fractions == []
        assert result.convergence.stop == _core.Stop.infeasible and failure in result.message
