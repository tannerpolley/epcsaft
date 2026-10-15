import math
from collections.abc import Sequence
from pathlib import Path

import pytest

from epcsaft import DeclarationError, Mixture, ParameterError, Parameters, _core, equilibrium
from epcsaft.equilibrium import (
    Amounts,
    Free,
    Molalities,
    Phase,
    Problem,
    problem_from_mapping,
    solve_equilibrium,
)

from tests.equilibrium.test_declared_starts import _mea_flash
from tests.equilibrium.test_trace_reliability import _MEA
from tests.native_inputs import aqueous_electrolyte
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


def test_named_guess_follows_the_declared_support() -> None:
    # A supported phase's guess has one entry per support component; a name map must start the same solve.
    def solve(first: object, second: object) -> _core.EquilibriumResult:
        return solve_equilibrium(mixture(), Problem(T=200.0, P=2e6, feed=equilibrium.MoleFractions({"methane": 0.5, "ethane": 0.5}),
            phases=[Phase("Le", kind="liquid", support=["ethane"], composition_guess=first),  # type: ignore[arg-type]
                    Phase("Lm", kind="liquid", support=["methane"], composition_guess=second)]))  # type: ignore[arg-type]
    positional, named = solve((1.0,), (1.0,)), solve({"ethane": 1.0}, {"methane": 1.0})
    assert positional.success and named.mole_fractions == positional.mole_fractions
    assert named.molar_densities == positional.molar_densities and named.amounts == positional.amounts


def test_named_guess_follows_native_support_order() -> None:
    # Native sorts a declared support into component order; a name map must not follow the caller's order.
    model = Mixture(Parameters.from_json(
        parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7), components=("water", "isobutanol")))

    def solve(first: object, second: object) -> _core.EquilibriumResult:
        return solve_equilibrium(model, Problem(T=293.15, P=1e5, feed=equilibrium.Amounts({"water": .6, "isobutanol": .4}),
            phases=[Phase("Lw", kind="liquid", support=["isobutanol", "water"], composition_guess=first),  # type: ignore[arg-type]
                    Phase("Lo", kind="liquid", support=["isobutanol", "water"], composition_guess=second)]))  # type: ignore[arg-type]
    positional = solve((.999, .001), (.3, .7))
    named = solve({"isobutanol": .001, "water": .999}, {"isobutanol": .7, "water": .3})
    assert positional.success and named.mole_fractions == positional.mole_fractions and named.amounts == positional.amounts


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
    assert default.requested_tolerance_met is True
    requested = Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1.0e5, feed=feed, row_tolerance=1.0e-8)
    assert requested._native(model).row_tolerance == 1.0e-8
    result = solve_equilibrium(model, requested)
    criterion = result.convergence.criterion
    del result
    assert criterion is not None
    assert criterion.value == 1.0e-8
    with pytest.raises(ParameterError):
        solve_equilibrium(model, Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1.0e5, feed=feed, row_tolerance=-1.0))

    defaults = {entry.name: entry.value for entry in _core.tolerances()}
    inherited = Mixture(model.parameters, solver_options={"EQUILIBRIUM_ROW_RESIDUAL": 1e-7})
    inherited_result = solve_equilibrium(inherited, problem(feed))
    assert inherited_result.convergence.criterion is not None
    assert inherited_result.convergence.criterion.value == 1e-7
    action = equilibrium.SolvedStateActionRequest(
        [equilibrium.SolvedStateObservable(equilibrium.SolvedStateObservableKind.Pressure)],
        [equilibrium.SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [])])
    with pytest.raises(ParameterError, match="central solver options"):
        equilibrium.solved_state_actions(equilibrium.compile_problem(inherited, problem(feed)), default, action)
    for tolerance in (defaults["EQUILIBRIUM_ROW_RESIDUAL"], 1e-8, 1e-12):
        ordinary = Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1e5,
            feed=Amounts({"methane": 50.0, "ethane": 50.0}), solver_options={"EQUILIBRIUM_ROW_RESIDUAL": tolerance})
        solved = solve_equilibrium(inherited, ordinary)
        assert solved.convergence.criterion is not None
        assert solved.convergence.criterion.value == tolerance
        inner = tolerance * defaults["EQUILIBRIUM_INNER_TOLERANCE_FACTOR"]
        for name in ("tol", "acceptable_tol", "constr_viol_tol"):
            assert math.isclose(solved.solver_settings[name], inner, rel_tol=1e-15)
        assert math.isclose(solved.solver_settings["dual_inf_tol"], 100.0 * inner, rel_tol=1e-15)
        # The declared-start route always runs Ipopt even on a square bubble cell.
        bubble = Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=equilibrium.Pinned(0.0))],
            T=180.0, P=Free(1e5), feed=feed, solver_options={"EQUILIBRIUM_ROW_RESIDUAL": tolerance})
        starts = equilibrium.solve_equilibrium_starts(inherited, bubble, [equilibrium.EquilibriumStart()], equilibrium.EquilibriumWorkLimits())
        candidate = starts.results[0]
        assert candidate is not None, starts.causes
        for name in ("tol", "acceptable_tol", "constr_viol_tol", "dual_inf_tol"):
            assert math.isclose(candidate.solver_settings[name], tolerance, rel_tol=1e-15)
    for options in ({"unknown": 1.0}, {"CLOSE_PACKING_FRACTION": 0.7}, {"DENSITY_ROOT_RESIDUAL": 0.0},
                    {"DENSITY_ROOT_RESIDUAL": -1.0}, {"DENSITY_ROOT_RESIDUAL": math.inf},
                    {"ASSOCIATION_RESIDUAL": math.nan}, {"ASSOCIATION_RESIDUAL": True},
                    {"ASSOCIATION_RESIDUAL": "bad"}, {"FIT_FUNCTION_TOLERANCE": 1e-8}):
        with pytest.raises(ParameterError, match=next(iter(options))):
            Mixture(model.parameters, solver_options=options)  # type: ignore[arg-type]
    with pytest.raises(ParameterError, match="outside owner scope"):
        Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1e5, feed=feed,
            solver_options={"DENSITY_ROOT_RESIDUAL": 1e-8})._native(model)
    with pytest.raises(ParameterError, match="duplicate"):
        Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1e5, feed=feed,
            row_tolerance=1e-8, solver_options={"EQUILIBRIUM_ROW_RESIDUAL": 1e-8})._native(model)
    with pytest.raises(ParameterError, match="EQUILIBRIUM_ROW_RESIDUAL"):
        Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1e5, feed=feed, row_tolerance=True)._native(model)
    with pytest.raises(ParameterError, match="FIT_FUNCTION_TOLERANCE"):
        _core.FitControls(maximum_iterations=1, maximum_elapsed_time_seconds=10, function_tolerance=True)
    for tolerance in (math.ulp(0.0), float.fromhex("0x1.fffffffffffffp+1023")):
        with pytest.raises(ParameterError, match="EQUILIBRIUM_ROW_RESIDUAL.*dual_inf_tol"):
            solve_equilibrium(model, Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1e5,
                feed=Amounts({"methane": 50.0, "ethane": 50.0}), solver_options={"EQUILIBRIUM_ROW_RESIDUAL": tolerance}))
    with pytest.raises(ParameterError, match="duplicate"):
        _core.FitControls(maximum_iterations=1, maximum_elapsed_time_seconds=10,
            function_tolerance=1e-8, solver_options={"FIT_FUNCTION_TOLERANCE": 1e-8})
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
    assert not result.success
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
    for changes in (
        {"phases": [{"name": "L", "kind": "liquid", "amount": 1.0}]},
        {"phases": [{"name": "L", "kind": "liquid", "amount": {"status": "pinned"}}]},
        {"T": {"free": False}},
        {"topology_declared": "false"},
    ):
        with pytest.raises(ParameterError):
            mapped(**changes)


def test_incidental_ipopt_options_preserve_equilibrium_state(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("EPCSAFT_IPOPT_OPTIONS_FILE", raising=False)
    monkeypatch.chdir(tmp_path)
    model = Mixture(Parameters.from_json(parameters_json("mea-co2-h2o-nine-species-estimation", 3), components=_MEA))
    declaration = _mea_flash(5918.5286, Phase("V", kind="vapor", support=_MEA[:3]))
    clean = solve_equilibrium(model, declaration)
    (tmp_path / "ipopt.opt").write_text("mu_strategy adaptive\nmax_iter 0\nbogus_option 1\n")
    stray = solve_equilibrium(model, declaration)
    assert clean.success and stray.success
    for key in ("amounts", "mole_fractions", "molar_densities"):
        assert all(math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-12)
                   for x, y in zip(getattr(stray, key), getattr(clean, key), strict=True))


def test_ipopt_initialization_or_fatal_return_publishes_no_state(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # An explicitly selected options file with an unknown option fails Initialize(); an unavailable
    # linear solver fails OptimizeTNLP() before any iterate, leaving the fallback start, whose classified
    # rows (1.50) meet this deliberately loose criterion. A zero iteration budget is an eligible
    # non-error stop: its polished start meets the rows and still publishes.
    model = mixture()
    loose = Problem(phases=[Phase("L", kind="liquid")], T=300.0, P=1.0e5,
                    feed=Amounts({"methane": 0.5, "ethane": 0.5}), row_tolerance=2.0)
    monkeypatch.chdir(tmp_path)
    missing = tmp_path / "missing.opt"
    monkeypatch.setenv("EPCSAFT_IPOPT_OPTIONS_FILE", str(missing))
    with pytest.raises(ParameterError, match="EPCSAFT_IPOPT_OPTIONS_FILE") as failure:
        solve_equilibrium(model, loose)
    assert str(missing) in str(failure.value)
    monkeypatch.setenv("EPCSAFT_IPOPT_OPTIONS_FILE", str(tmp_path / "diagnostic.opt"))
    for option, fatal in (("bogus_option 1", True), ("linear_solver ma57", True), ("max_iter 0", False)):
        (tmp_path / "diagnostic.opt").write_text(option + "\n")
        result = solve_equilibrium(model, loose)
        if not fatal:
            assert result.success and math.isclose(result.amounts[0], 1.0)
            continue
        assert not result.success and result.amounts == [] and result.mole_fractions == []
        # The fallback rows meet the criterion: the fatal return alone withholds the state.
        assert result.requested_tolerance_met is True and result.convergence.stop == _core.Stop.infeasible


def test_supplied_charged_family_refuses_before_trial_work() -> None:
    # ALG-STABILITY: density-independent charged admission, including packing replay.
    species = ("water", "sodium-cation", "chloride-anion")
    model = Mixture(Parameters.from_json(
        parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7), components=species,
    ))
    declaration = Problem(T=298.15, P=1e5, row_tolerance=1e-10,
        feed=Amounts(dict(zip(species, (0.98, 0.01, 0.01)))), phases=[Phase("L", kind="liquid")])
    central = solve_equilibrium(model, declaration)
    assert central.success
    compiled = equilibrium.compile_problem(model, declaration)
    charged = declaration._native(model).phases[0]
    charged.support = [0, 1, 2]
    charged.composition_guess = central.mole_fractions
    charged.packing_guess = central.packing_fractions[0]
    neutral = _core.PhaseSpec()
    neutral.support, neutral.composition_guess = [0], [1.0, 0.0, 0.0]
    neutral.packing_guess = central.packing_fractions[0]
    options = equilibrium.StabilityOptions()
    options.evaluate_supplied_trials, options.trials = True, [neutral, charged, neutral]
    with _core.Capture() as captured:
        with pytest.raises(ParameterError, match="^permittivity model has no low-density limit; declare phases$"):
            equilibrium.stability_test(compiled, central, options)
    assert all(s.name not in ("tangent_plane_trial", "solve_density") for s in captured.spans)
    assert all(e.name != "hessian_calls" for e in captured.events)
    options.trials = [neutral]
    assert "tpd_over_rt" in equilibrium.stability_test(compiled, central, options).trials[0].quantities


def test_omitted_vapor_support_is_the_zero_charge_components() -> None:
    """#438: the default bubble_point vapor holds only zero-charge components. NaCl 1 mol/kg at 373.15 K (Held 2008
    packet 4) and packet-3 MEA at 318.15 K equal the explicit neutral vapor and the pressures retained from it on
    the unchanged Engine (numerical verification, not validation); the ions stay in the balanced liquid. Full-width
    guesses and starts map onto the neutral vapor, an omitted vapor trial resolves to it, and an explicit charged
    vapor is refused before the reference calculation that would otherwise fail at 400 K."""
    nacl = ("water", "sodium-cation", "chloride-anion")
    salt = Mixture(Parameters.from_json(parameters_json("held-cameretti-sadowski-2008", 4), components=nacl))
    mea = Mixture(Parameters.from_json(parameters_json("mea-co2-h2o-nine-species-estimation", 3), components=_MEA))
    reactions = _mea_flash(0.).reactions
    cases: list[tuple[Mixture, float, Amounts | Molalities, Sequence[equilibrium.Reaction], float, int]]
    cases = [(salt, 373.15, Molalities({"sodium-cation": 1., "chloride-anion": 1.}, solvent=("water",), solvent_mass_kg=1.),
              (), 97924.5272181854, 1)]
    cases += [(mea, 318.15, Amounts({"carbon-dioxide": a * 9.7292821, "monoethanolamine": 9.7292821, "water": 76.9751925}),
               reactions, p, 3) for a, p in ((.1, 7286.707555004345), (.3, 7427.553577760993), (.5, 7522.2862426425445))]
    for model, T, feed, chemistry, pressure, neutral in cases:
        n = len(model.component_ids)
        default, explicit = (equilibrium.bubble_point(liquid=Phase("L", kind="liquid"), T=T, feed=feed, reactions=chemistry)
                             for _ in range(2))
        assert isinstance(explicit.phases, tuple)
        explicit.phases[1].support = model.component_ids[:neutral]
        compiled = equilibrium.compile_problem(model, default)
        assert [r.index for r in compiled.rows if r.kind == _core.RowKind.Chemical and r.phase == 1] == list(range(neutral))
        result, reference = solve_equilibrium(model, default), solve_equilibrium(model, explicit)
        assert result.success and reference.success, result.message
        assert abs(result.pressure[0] - pressure) <= 1e-10 * pressure
        for key in ("pressure", "amounts", "mole_fractions", "molar_densities"):
            for x, y in zip(getattr(result, key), getattr(reference, key), strict=True):
                assert abs(x - y) <= 1e-10 * max(abs(x), abs(y), 1e-14), key
        assert result.amounts[1] == 0. and result.mole_fractions[n + neutral:] == [0.] * (n - neutral)
        x = dict(zip(model.component_ids, (result.amounts[0] * v for v in result.mole_fractions[:n]), strict=True))
        balances = [(x["sodium-cation"], 1.), (x["chloride-anion"], 1.)] if model is salt else [
            (x["carbon-dioxide"] + x["carbamate-anion"] + x["bicarbonate-anion"] + x["carbonate-anion"],
             feed.values["carbon-dioxide"]),
            (x["monoethanolamine"] + x["protonated-monoethanolamine"] + x["carbamate-anion"], 9.7292821),
            (x["protonated-monoethanolamine"] + x["hydronium-cation"], x["carbamate-anion"] + x["bicarbonate-anion"]
             + 2. * x["carbonate-anion"] + x["hydroxide-anion"])]
        for amount, target in balances:
            assert abs(amount - target) <= 1e-10 * target
    guess, vapor = ([*head, *[0.] * (n - neutral)] for head in ([.3, .3, .4], reference.mole_fractions[n:n + neutral]))
    start = equilibrium.EquilibriumStart(compositions=[result.mole_fractions[:n], vapor])
    narrow = equilibrium.EquilibriumStart(compositions=[result.mole_fractions[:n], vapor[:neutral]])
    pairs = ((solve_equilibrium(mea, _mea_flash(5918.5286, Phase("V", kind="vapor", composition_guess=guess))),
              solve_equilibrium(mea, _mea_flash(5918.5286, Phase("V", kind="vapor", support=_MEA[:3],
                                                                 composition_guess=guess[:neutral])))),
             *(equilibrium.solve_equilibrium_starts(mea, problem, [seed], equilibrium.EquilibriumWorkLimits()).results
               for problem, seed in ((default, start), (explicit, narrow))))
    for full, reduced in (pairs[0], (pairs[1][0], pairs[2][0])):
        assert full is not None and reduced is not None and full.success
        assert [full.iterations, full.evaluations] == [reduced.iterations, reduced.evaluations]
        for x, y in zip([*full.pressure, *full.amounts], [*reduced.pressure, *reduced.amounts], strict=True):
            assert abs(x - y) <= 1e-10 * max(abs(x), abs(y))
    options = equilibrium.StabilityOptions()
    options.trials = [_core.PhaseSpec()]
    options.trials[0].kind = _core.PhaseKind.Vapor
    salt_default = equilibrium.bubble_point(liquid=Phase("L", kind="liquid"), T=373.15, feed=cases[0][2])
    salt_result = solve_equilibrium(salt, salt_default)
    report = equilibrium.stability_test(equilibrium.compile_problem(salt, salt_default), salt_result, options)
    assert report.trials[0].phase.support == [0]
    law = equilibrium.ReactionLogPolynomial(math.log(1e-8), 0., 0., 0., 298.15, 250., 350., "analytic-source")
    ionic = Mixture(aqueous_electrolyte())
    for T, support, error, message in ((373.15, ("water", "sodium-cation"), DeclarationError, "charged component index 1"),
                                       (400., ("s0", "s1"), DeclarationError, "charged component index 1"),
                                       (400., ("s0",), ParameterError, "reference_unavailable")):
        model, values = (salt, {"water": 1.}) if T == 373.15 else (ionic, {"s0": 1., "s1": 1e-6, "s2": 1e-6})
        refused = equilibrium.bubble_point(liquid=Phase("L", kind="liquid"), T=T, feed=Amounts(values), reactions=() if model is salt else [
            equilibrium.Reaction({"s0": -1., "s1": 1., "s2": 1.}, correlation=law, reference=equilibrium.ReactionReference(
                equilibrium.ReferenceSourceBasis.RawMoleFractionInfiniteDilution, None))])
        assert isinstance(refused.phases, tuple)
        refused.phases[1].support = support
        if model is ionic:
            refused.neutral_reference = equilibrium.NeutralReference(
                solvent="s0", solute_molality_ratios={"s1": 1., "s2": 1.}, standard_molality_mol_per_kg=1.)
        with pytest.raises(error, match=message):
            equilibrium.compile_problem(model, refused)


def test_mid_dome_flash_returns_the_split() -> None:
    """#314 review: without ln eta bounds this two-phase flash fell to the trivial point."""
    result = solve_equilibrium(mixture(), Problem(
        phases=[Phase("V", kind="vapor"), Phase("L", kind="liquid")],
        T=250.0, P=3.75e6, feed=Amounts({"methane": 0.3, "ethane": 0.7})))
    assert result.success, result.message
    vapor, liquid = result.mole_fractions[0], result.mole_fractions[2]  # methane, phase-major
    assert vapor - liquid > 0.1 and 0.0 < result.amounts[0] < 1.0
