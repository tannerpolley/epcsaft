"""Ordered declared starts: analytic exchange amounts, a shared work quota (#238) and the shared cold start (#350)."""
import importlib.util
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest

from epcsaft import Mixture, Parameters, Pure
from epcsaft import equilibrium as eq
from epcsaft.equilibrium import EOS_STANDARD_STATE_ID, Reaction, ReactionLogPolynomial
from tests.equilibrium.test_trace_reliability import _MEA, _MEA_REACTIONS, _mea_mapping
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
    assert reference.success and reference.pressure is not None
    pinned_seed = eq.EquilibriumStart()
    pinned_seed.compositions = [[1.], [1.]]
    pinned_seed.packing_fractions = [*reference.packing_fractions]
    pinned_seed.pressure = reference.pressure
    for amount in (float("nan"), 1.):
        pinned_seed.amounts = [1., amount]
        refused = eq.solve_equilibrium_starts(pure, bubble, [pinned_seed], eq.EquilibriumWorkLimits())
        assert refused.selected is None and refused.results[0] is None and refused.iterations == 0
        assert "amount for phase V" in refused.causes[0]
    pinned_seed.amounts = [1., 0.]
    admitted = eq.solve_equilibrium_starts(pure, bubble, [pinned_seed], eq.EquilibriumWorkLimits())
    assert admitted.selected == 0 and admitted.results[0] is not None and admitted.results[0].pressure is not None
    assert abs(admitted.results[0].pressure/reference.pressure-1) <= 1e-9
    previous.previous = reference
    bubble.P = eq.Free(2e5)
    reused = eq.solve_equilibrium_starts(pure, bubble, [previous], eq.EquilibriumWorkLimits())
    assert reused.selected == 0 and reused.results[0] is not None and reused.results[0].pressure is not None
    assert abs(reused.results[0].pressure/reference.pressure-1) <= 1e-9
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
    binary = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3),
                                         components=("methane", "ethane")))
    declared = eq.Problem(T=300., P=1e5, phases=[eq.Phase("L", kind="liquid")],
                          feed=eq.Amounts({"methane": .5, "ethane": .5}), row_tolerance=1e-8)
    reference = eq.solve_equilibrium(binary, declared)
    assert reference.success and reference.solver_settings
    timed = eq.solve_equilibrium_starts(binary, declared, [eq.EquilibriumStart(coordinates=reference.coordinates)],
                                        eq.EquilibriumWorkLimits(attempt_seconds=.001))
    result = timed.results[0]
    if result is not None and result.evaluations > 0:
        assert result.solver_settings == reference.solver_settings
        assert result.diagnostic_options_file == reference.diagnostic_options_file
        if "elapsed" in timed.causes[0]:
            assert not result.success and timed.selected is None
            assert result.message == timed.causes[0] and timed.iterations == result.iterations
        evidence.observe("declared-starts/elapsed-work", cause=timed.causes[0], iterations=result.iterations,
                         evaluations=result.evaluations, solver_settings=result.solver_settings,
                         diagnostic_options_file=result.diagnostic_options_file)


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


def _mea_flash(pressure: Any, *vapor: eq.Phase) -> eq.Problem:
    """Packet-3 MEA at loading 0.2556 and 318.15 K: a liquid and the declared vapor phases."""
    reactions = [Reaction(nu, correlation=ReactionLogPolynomial(a, 0., 0., 0., 318.15, 250., 450., EOS_STANDARD_STATE_ID))
                 for nu, a in zip(_MEA_REACTIONS, (-605.348, -587.8, -948.2312, -25.1014, -298.6415), strict=True)]
    feed = eq.Amounts({"carbon-dioxide": .25557798094241313 * 9.7292821, "monoethanolamine": 9.7292821,
                       "water": 76.9751925})
    return eq.Problem(T=318.15, P=pressure, phases=[eq.Phase("L", kind="liquid"), *vapor], feed=feed,
                      reactions=reactions, row_tolerance=1e-10)


def test_default_start_replays_the_cold_start_and_counts_the_bubble_presolve() -> None:
    """#350: packet-3 MEA at loading 0.2556 and 318.15 K. A default start runs the plain solve's cold start (the
    true-species split, the single-liquid speciation, the bubble's liquid pre-solve at its pressure seed and
    fugacity pressure step) under its work limits, so it returns the same state and summed work. The bubble reports
    that pre-solve's work, and a budget below it stops before the bubble solve."""
    model = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 3)), components=_MEA))
    bubble = _mea_flash(eq.Free(7398.16), eq.Phase("V", kind="vapor", support=_MEA[:3], amount=eq.Pinned(0.)))
    work = ("iterations", "evaluations", "hessian_calls")
    for declared in (bubble, _mea_flash(5918.5286, eq.Phase("V", kind="vapor", support=_MEA[:3])), _mea_flash(8877.7929)):
        plain = eq.solve_equilibrium(model, declared)
        started = eq.solve_equilibrium_starts(model, declared, [eq.EquilibriumStart()], eq.EquilibriumWorkLimits())
        result = started.results[0]
        assert plain.success and started.selected == 0 and result is not None, started.causes
        assert [getattr(result, k) for k in work] == [getattr(plain, k) for k in work]
        assert started.iterations == plain.iterations
        assert plain.pressure is not None and result.pressure is not None
        assert abs(plain.pressure - result.pressure) <= 1e-10 * plain.pressure
        for key in ("amounts", "mole_fractions", "molar_densities"):
            for x, y in zip(getattr(plain, key), getattr(result, key), strict=True):
                assert abs(x - y) <= 1e-10 * max(abs(x), abs(y), 1e-14), key
    presolve = eq.solve_equilibrium(model, _mea_flash(7398.16))
    plain = eq.solve_equilibrium(model, bubble)
    assert all(getattr(plain, k) > getattr(presolve, k) for k in work)
    stopped = eq.solve_equilibrium_starts(model, bubble, [eq.EquilibriumStart()],
                                          eq.EquilibriumWorkLimits(iterations=presolve.iterations - 1))
    assert stopped.selected is None and stopped.results[0] is None and stopped.iterations == presolve.iterations - 1
    assert "iteration budget" in stopped.causes[0], stopped.causes


def test_feed_start_regrows_the_vapor_at_the_356_point() -> None:
    """#356: the MEA flash at 3699.0803800848512 Pa. From the feed start (both phases at the feed composition, unfed
    liquid species at 1/9 of the feed before normalization, equal amounts) one step nearly emptied the liquid on main
    7fa511b6, and the solve ended "phase V vanished"; the amount step rule refuses that step. That start and the
    generated true-species split both reach main's apparent-K state within 1e-10."""
    model = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 3)), components=_MEA))
    problem = _mea_flash(3699.0803800848512, eq.Phase("V", kind="vapor", support=_MEA[:3]))
    z = [.25557798094241313 * 9.7292821, 9.7292821, 76.9751925]
    liquid = z + [sum(z) / 9] * 6
    feed = eq.EquilibriumStart(compositions=[[v / sum(liquid) for v in liquid], [v / sum(z) for v in z]],
                               amounts=[sum(z) / 2] * 2)
    started = eq.solve_equilibrium_starts(model, problem, [feed], eq.EquilibriumWorkLimits()).results[0]
    reference = [27.345065157988326, 59.35940946171173, 2.0609303345668776e-14, .1313557341458005, .5828610304006816,
                 .18833784813341686, 4.0587181258156585e-05, 5.054849860062213e-07, .0908926915858919,
                 1.1538627774106781e-07, .00651148768166661, 3.318679391441496e-10, .01661287421506598,
                 .9833871254530661, 0., 0., 0., 0., 0., 0., 42120.183730573895, 1.4010589629129806]
    for result in (eq.solve_equilibrium(model, problem), started):
        assert result is not None and result.success, result and result.message
        state = [*result.amounts, *result.mole_fractions, *result.molar_densities]
        assert all(abs(x - y) <= 1e-10 * max(abs(x), abs(y), 1e-14) for x, y in zip(state, reference, strict=True))


def test_true_species_split_ignores_an_unfed_spectator() -> None:
    """#356: a reacting flash s0 -> s1 (K = 2) at 200 K and 50 kPa with an unfed, nonreacting s3. Default supports and
    supports declared without s3 compile to the same effective supports, so both take the true-species split with the
    same work and state; the supports without s3 took the apparent-K start (18/50/19 vs 9/14/10)."""
    heavy = (*ETHANE[:3], 2 * ETHANE[3])
    model = Mixture(parameters((ETHANE, ETHANE, heavy, ETHANE), [("k_ij", (i, j), 0.) for i in range(4) for j in range(i + 1, 4)]))

    def solve(support: Any) -> Any:
        return eq.solve_equilibrium(model, eq.Problem(
            T=200., P=5e4, feed=eq.Amounts({"s0": 1., "s2": 1.}),
            phases=[eq.Phase("L", kind="liquid", support=support), eq.Phase("V", kind="vapor", support=support)],
            reactions=[eq.Reaction({"s0": -1., "s1": 1.}, equilibrium_constant=2.)]))

    default, declared = solve(None), solve(("s0", "s1", "s2"))
    for result in (default, declared):
        assert result.success, result.message
    work = ("iterations", "evaluations", "hessian_calls")
    assert [getattr(declared, k) for k in work] == [getattr(default, k) for k in work]
    for key in ("amounts", "mole_fractions", "molar_densities"):
        for x, y in zip(getattr(default, key), getattr(declared, key), strict=True):
            assert abs(x - y) <= 1e-10 * max(abs(x), abs(y), 1e-14), key


def test_supplied_packings_need_no_boundary_density_root() -> None:
    """#350: Born case 1 of analyses/2026-born-unavailable-states, the stable liquid at c_dielectric 0.33 (archive
    member fine/case1-continuation-0.330000-state.csv of results/retired-charts.tar.gz) moved by a water-reaction
    extent of 0.2. The start supplies both compositions, both packing fractions and the pressure but
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


@pytest.mark.parametrize(("species", "association", "T", "reference"), [
    ((32.042, 1.5255, 3.2300, 188.90), (2899.5, .035176), 505., (7098526.457937359, 14223.498914825788, 3264.4582146128782)),
    ((88.15, 3.6260, 3.4508, 247.28), (2252.1, .010319), 585., (3753197.3286717464, 4685.426417139512, 1511.9545691073597))])
def test_pure_saturation_default_start_takes_the_loop_pressure(species: Any, association: Any, T: float, reference: Any) -> None:
    """#391: Gross and Sadowski (2002) Table 1 methanol and 1-pentanol (2B) saturation, from the default start with
    the 0.1 MPa guess. That pressure has only a vapor root; main 3b307095 refused it as "no coexisting density
    roots". The pressure at the density scan's negative-slope node seeds both roots, and the solve reaches the
    retained declared-start state of analyses/reproduction/2002-gross (state-comparisons.csv lines 21 and 41)."""
    model = Pure(parameters((species,), [("association_energy_over_k", (0, 1), association[0]),
                                         ("association_volume", (0, 1), association[1])], ((0, 0), (1, 1))))
    result = eq.solve_equilibrium(model, eq.Problem(
        T=T, P=eq.Free(1e5), feed=eq.Amounts({"s0": 1.}), row_tolerance=1e-8,
        phases=[eq.Phase("L", kind="liquid"), eq.Phase("V", kind="vapor", amount=eq.Pinned(0.))]))
    assert result.success and result.requested_tolerance_met, result.message
    for x, y in zip([result.pressure, *result.molar_densities], reference, strict=True):
        assert abs(x - y) <= 1e-10 * y
    assert result.molar_densities[0] > 2 * result.molar_densities[1]
    assert all(model.state(T, rho=rho, x=[1.]).pressure_density_derivative > 0 for rho in result.molar_densities)
