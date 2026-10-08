"""Ordered declared starts: analytic exchange amounts, a shared work quota (#238) and the shared cold start (#350)."""
import importlib.util
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from epcsaft import Mixture, Parameters, Pure
from epcsaft import equilibrium as eq
from epcsaft.equilibrium import EOS_STANDARD_STATE_ID, Reaction, ReactionLogPolynomial
from tests.equilibrium.test_trace_reliability import _MEA, _MEA_REACTIONS
from tests.native_inputs import ETHANE, parameters
from tests.packets import parameters_json


def test_declared_starts_close_analytic_exchange_and_bound_work(evidence: Any) -> None:
    p = parameters((ETHANE,) * 6, pairs=[("k_ij", (i, j), 0.) for i in range(6) for j in range(i+1, 6)])
    model = Mixture(p)
    problem = eq.Problem(T=200., P=5e6,
        phases=[eq.Phase("aqueous", kind="liquid", support=("s0", "s2", "s4")),
                eq.Phase("organic", kind="liquid", support=("s1", "s3", "s5"))],
        feed=eq.Amounts({"s0": 1., "s1": 1., "s2": .3, "s4": .3}), row_tolerance=1e-10,
        reactions=[eq.Reaction({"s2": -1., "s3": 1.}, equilibrium_constant=2., name="exchange_a"),
                   eq.Reaction({"s4": -1., "s5": 1.}, equilibrium_constant=.5, name="exchange_b")])
    exact = ([10/13, 1/13, 2/13], [10/13, 2/13, 1/13])
    packing = []
    for support, x in zip(((0, 2, 4), (1, 3, 5)), exact, strict=True):
        full = [0.] * 6
        for i, value in zip(support, x, strict=True):
            full[i] = value
        packing.append(model._model.packing_fraction_coefficient(200., full)
                       * model.state(200., P=5e6, x=full, phase="liquid").molar_density)

    def seed(x: Any, eta: Any) -> eq.EquilibriumStart:
        return eq.EquilibriumStart(compositions=x, packing_fractions=eta, amounts=[1.3, 1.3])

    invalid = seed(exact, [-1., packing[1]])
    positive = seed(exact, packing)
    later = seed(([.8, .08, .12], [.8, .12, .08]), packing)
    solved = eq.solve_equilibrium_starts(model, problem, [invalid, positive, later], eq.EquilibriumWorkLimits())
    assert solved.selected == 1, solved.causes
    assert solved.results[0] is None and solved.results[2] is None
    result = solved.results[1]
    assert result is not None
    for phase, support in enumerate(((0, 2, 4), (1, 3, 5))):
        amount_error = abs(result.amounts[phase] - 1.3)
        fraction_error = max(abs(result.mole_fractions[6*phase+i] - x)
                             for i, x in zip(support, exact[phase], strict=True))
        assert amount_error <= 1e-9 and fraction_error <= 1e-9
        for quantity, value, unit in (("amount", amount_error, "mol"), ("composition", fraction_error, "1")):
            evidence.compare(f"declared-starts/exchange/{phase}/{quantity}", model=value, reference=0., unit=unit,
                             criterion={"rule": "symmetric_close", "atol": 1e-9})
    raw = max(abs(r) for r in result.residuals)
    assert raw <= 1e-10
    for a, b, target in ((2, 3, 2.), (4, 5, .5)):
        ratio = result.amounts[1]*result.mole_fractions[6+b] / (result.amounts[0]*result.mole_fractions[a])
        assert abs(ratio-target) <= 1e-8
        evidence.compare(f"declared-starts/ratio/{a}", model=ratio, reference=target, unit="1",
                         criterion={"rule": "symmetric_close", "atol": 1e-8})
    for replay in (eq.EquilibriumStart(coordinates=result.coordinates), eq.EquilibriumStart(previous=result)):
        restored = eq.solve_equilibrium_starts(model, problem, [replay], eq.EquilibriumWorkLimits())
        assert restored.selected == 0 and restored.results[0] is not None
        assert max(map(abs, restored.results[0].residuals)) <= 1e-10
    previous = eq.EquilibriumStart(previous=result)
    for field, changed_value in (("P", 4e6), ("T", 201.), ("feed", eq.Amounts({"s0": 2., "s1": 1., "s2": .3, "s4": .3}))):
        different = deepcopy(problem)
        setattr(different, field, changed_value)
        with pytest.raises(ValueError, match="fixed pressure|fixed temperature|feed"):
            eq.solve_equilibrium_starts(model, different, [positive, previous], eq.EquilibriumWorkLimits())
    partial = seed(exact, [packing[0], None])
    partial_result = eq.solve_equilibrium_starts(model, problem, [partial], eq.EquilibriumWorkLimits())
    assert partial_result.selected == 0 and partial_result.results[0] is not None
    assert max(map(abs, partial_result.results[0].residuals)) <= 1e-10
    pure = Pure(parameters((ETHANE,)))
    bubble = eq.Problem(T=200., P=eq.Free(1e5), feed=eq.Amounts({"s0": 1.}),
                        phases=[eq.Phase("L", kind="liquid"), eq.Phase("V", kind="vapor", amount=eq.Pinned(0.))],
                        row_tolerance=1e-10)
    reference = eq.solve_equilibrium(pure, bubble)
    assert reference.success
    pinned_seed = eq.EquilibriumStart()
    pinned_seed.compositions = [[1.], [1.]]
    pinned_seed.packing_fractions = [*reference.packing_fractions]
    pinned_seed.pressure = reference.pressure[0]
    for amount in (float("nan"), 1.):
        pinned_seed.amounts = [1., amount]
        refused = eq.solve_equilibrium_starts(pure, bubble, [pinned_seed], eq.EquilibriumWorkLimits())
        assert refused.selected is None and refused.results[0] is None and refused.iterations == 0
        assert "amount for phase V" in refused.causes[0]
    pinned_seed.amounts = [1., 0.]
    admitted = eq.solve_equilibrium_starts(pure, bubble, [pinned_seed], eq.EquilibriumWorkLimits())
    assert admitted.selected == 0 and admitted.results[0] is not None
    assert abs(admitted.results[0].pressure[0]/reference.pressure[0]-1) <= 1e-9
    previous.previous = reference
    bubble.P = eq.Free(2e5)
    reused = eq.solve_equilibrium_starts(pure, bubble, [previous], eq.EquilibriumWorkLimits())
    assert reused.selected == 0 and reused.results[0] is not None
    assert abs(reused.results[0].pressure[0]/reference.pressure[0]-1) <= 1e-9
    only = eq.solve_equilibrium_starts(model, problem, [invalid], eq.EquilibriumWorkLimits())
    assert only.selected is None and only.iterations == 0
    for name in ("iterations", "seconds"):
        stopped = eq.solve_equilibrium_starts(model, problem, [positive], eq.EquilibriumWorkLimits(**{name: 0}))
        assert stopped.selected is None and stopped.budget_exhausted and stopped.results[0] is None
    limits = eq.EquilibriumWorkLimits(iterations=1)
    far = seed(([.8, .08, .12], [.8, .12, .08]), [.001, .001])
    stopped = eq.solve_equilibrium_starts(model, problem, [invalid, far, positive], limits)
    assert stopped.selected is None and stopped.budget_exhausted and stopped.iterations <= 1
    assert stopped.results[1] is not None and stopped.results[2] is None
    unfinished = max(abs(r) for r in stopped.results[1].residuals)
    assert unfinished > 1e-10
    evidence.observe("declared-starts/work", raw_rows=raw, limited_raw_rows=unfinished,
                     limited_iterations=stopped.iterations, selected=solved.selected)


def test_declared_compositions_on_zero_exchange_inventory(evidence: Any) -> None:
    model = Mixture(parameters((ETHANE,) * 6,
        pairs=[("k_ij", (i, j), 0.) for i in range(6) for j in range(i+1, 6)]))
    problem = eq.Problem(T=200., P=5e6, feed=eq.Amounts({"s0": 1., "s1": 1.}),
        phases=[eq.Phase("aqueous", kind="liquid", support=("s0", "s2", "s4"), composition_guess=(1., 0., 0.)),
                eq.Phase("organic", kind="liquid", support=("s1", "s3", "s5"), composition_guess=(1., 0., 0.))],
        reactions=[eq.Reaction({"s2": -1., "s3": 1.}, equilibrium_constant=2.),
                   eq.Reaction({"s4": -1., "s5": 1.}, equilibrium_constant=.5)],
        row_tolerance=1e-10)
    reference = eq.solve_equilibrium(model, problem)
    assert reference.success, reference.message
    start = eq.EquilibriumStart(compositions=[[1., 0., 0.], [1., 0., 0.]],
                                packing_fractions=list(reference.packing_fractions), amounts=[1., 1.])
    solved = eq.solve_equilibrium_starts(model, problem, [start], eq.EquilibriumWorkLimits())
    assert solved.selected == 0, solved.causes
    result = solved.results[0]
    assert result is not None
    for phase, solvent in enumerate((0, 1)):
        amount_error = abs(result.amounts[phase] - 1.)
        fraction_error = abs(result.mole_fractions[6*phase+solvent] - 1.)
        assert amount_error <= 1e-9 and fraction_error <= 1e-9
        assert all(result.mole_fractions[6*phase+i] == 0. for i in (2, 3, 4, 5))
        evidence.compare(f"declared-starts/zero-exchange/{phase}/amount", model=amount_error,
                         reference=0., unit="mol", criterion={"rule": "symmetric_close", "atol": 1e-9})
    raw = max(map(abs, result.residuals))
    assert raw <= 1e-10
    evidence.observe("declared-starts/zero-exchange", amounts=result.amounts,
                     compositions=result.mole_fractions, raw_rows=raw)
    for invalid_fraction in (-.1, 1.1, float("nan")):
        start.compositions = [[1., invalid_fraction, 0.], [1., 0., 0.]]
        refused = eq.solve_equilibrium_starts(model, problem, [start], eq.EquilibriumWorkLimits())
        assert refused.selected is None and refused.results[0] is None and refused.iterations == 0


def test_default_start_replays_the_cold_start_and_counts_the_bubble_presolve() -> None:
    """#350: packet-3 MEA at loading 0.2556 and 318.15 K. A default start runs the plain solve's cold start (the
    two-phase K-factor estimate, the single-liquid speciation, the bubble's liquid pre-solve at its pressure seed and
    fugacity pressure step) under its work limits, so it returns the same state and summed work. The bubble reports
    that pre-solve's work, and a budget below it stops before the bubble solve."""
    model = Mixture(Parameters.from_json(parameters_json("mea-co2-h2o-nine-species-estimation", 3), components=_MEA))
    reactions = [Reaction(nu, correlation=ReactionLogPolynomial(a, 0., 0., 0., 318.15, 250., 450., EOS_STANDARD_STATE_ID))
                 for nu, a in zip(_MEA_REACTIONS, (-605.348, -587.8, -948.2312, -25.1014, -298.6415), strict=True)]
    feed = eq.Amounts({"carbon-dioxide": .25557798094241313 * 9.7292821, "monoethanolamine": 9.7292821,
                       "water": 76.9751925})

    def problem(pressure: Any, *vapor: eq.Phase) -> eq.Problem:
        return eq.Problem(T=318.15, P=pressure, phases=[eq.Phase("L", kind="liquid"), *vapor], feed=feed,
                          reactions=reactions, row_tolerance=1e-10)

    bubble = problem(eq.Free(7398.16), eq.Phase("V", kind="vapor", support=_MEA[:3], amount=eq.Pinned(0.)))
    work = ("iterations", "evaluations", "hessian_calls")
    for declared in (bubble, problem(5918.5286, eq.Phase("V", kind="vapor", support=_MEA[:3])), problem(8877.7929)):
        plain = eq.solve_equilibrium(model, declared)
        started = eq.solve_equilibrium_starts(model, declared, [eq.EquilibriumStart()], eq.EquilibriumWorkLimits())
        result = started.results[0]
        assert plain.success and started.selected == 0 and result is not None, started.causes
        assert [getattr(result, k) for k in work] == [getattr(plain, k) for k in work]
        assert started.iterations == plain.iterations
        for key in ("pressure", "amounts", "mole_fractions", "molar_densities"):
            for x, y in zip(getattr(plain, key), getattr(result, key), strict=True):
                assert abs(x - y) <= 1e-10 * max(abs(x), abs(y), 1e-14), key
    presolve = eq.solve_equilibrium(model, problem(7398.16))
    plain = eq.solve_equilibrium(model, bubble)
    assert all(getattr(plain, k) > getattr(presolve, k) for k in work)
    stopped = eq.solve_equilibrium_starts(model, bubble, [eq.EquilibriumStart()],
                                          eq.EquilibriumWorkLimits(iterations=presolve.iterations - 1))
    assert stopped.selected is None and stopped.results[0] is None and stopped.iterations == presolve.iterations - 1
    assert "iteration budget" in stopped.causes[0], stopped.causes


def test_supplied_packings_need_no_boundary_density_root() -> None:
    """#350: Born case 1 of analyses/2026-born-unavailable-states, the stable liquid at c_dielectric 0.33 moved by a
    water-reaction extent of 0.2. The start supplies both compositions, both packing fractions and the pressure but
    no amounts. The liquid has no density root at the boundary pressure seed; the supplied packing replaces it, so
    the start reaches the driver instead of being refused at initialization."""
    script = Path(__file__).resolve().parents[3] / "analyses/2026-born-unavailable-states/scripts/reproduce.py"
    spec = importlib.util.spec_from_file_location("born_reproduce", script)
    assert spec is not None and spec.loader is not None
    born = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(born)
    mapping, request, anchor = born.load_case(1)
    problem = born.problem_from_request(request, anchor)
    stable = [4.3156821899964646e-13, .06744004970829823, .8601794757184709, .04476857189213988, 1.02472867796729e-05,
              2.0895783143159884e-06, .017157150074218773, 2.404313692039684e-07, .010442175309977538]
    water = request["reaction_system"]["reaction_matrix"][0]
    start = eq.EquilibriumStart()
    start.compositions = [[x + .2 * nu for x, nu in zip(stable, water, strict=True)],
                          list(problem.phases[1].composition_guess)]
    start.packing_fractions = [.4468151338619424, 2.0279177672656488e-05]
    start.pressure = anchor["pressure_pa"]
    model = Mixture(Parameters.from_mapping(mapping))
    tried = eq.solve_equilibrium_starts(model, problem, [start], eq.EquilibriumWorkLimits(iterations=1))
    assert not tried.causes[0].startswith("initialization"), tried.causes
    assert tried.results[0] is not None and tried.iterations == 1
