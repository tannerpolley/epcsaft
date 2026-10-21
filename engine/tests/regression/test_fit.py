"""Public fitting transport (#49): the native R2-J2/D86 rows through epcsaft.regression, name mapping and typed failures."""
import math
from typing import Any

import pytest

from epcsaft import (
    Mixture, ParameterActionDirection, ParameterError, ParameterFamily, Parameters, PressureActionDirection, Pure, _core,
    regression,
)
from epcsaft.equilibrium import Amounts, Continuation, Discover, Phase, Problem, asl_record, solve_equilibrium

from tests.native_inputs import ETHANE, ROWS, bubble, budget, aqueous_electrolyte, associating_binary, p3_forms, parameters, record, same

METHANE = (16.043, 1.0, 3.7039, 150.03)
BOUNDS = {"dispersion_energy_over_k": (140.0, 170.0), "k_ij": (-0.2, 0.2)}


def neutral_binary(epsilon: float = 150.03, kij: float = 0.0) -> Any:
    """properties_support::neutral_binary() with the R2-J2 truth or start values."""
    return parameters(((METHANE[0], METHANE[1], METHANE[2], epsilon), ETHANE), pairs=[("k_ij", (0, 1), kij)])


def manufactured() -> tuple[Any, list[Any]]:
    """The R2-J2 truth (epsilon 152 K, k_ij -0.02) and its noise-free observations."""
    truth = neutral_binary(152.0, -0.02)
    model = Mixture(truth)
    observations = [regression.observation(truth, "pressure", bubble(T), observed=1.0, form="log_ratio", mass_basis=False) for T in (170.0, 180.0, 190.0)]
    observations.append(regression.observation(truth, "mole_fraction", bubble(180.0), observed=1.0, phase="V", species=("s0",),
                                               form="difference", scale=0.01, mass_basis=False))
    P_bubble = solve_equilibrium(model, bubble(180.0)).pressure
    assert P_bubble is not None
    observations.append(regression.density(truth, 180.0, 2.0 * P_bubble, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=1.0,
                                           form="log_ratio", mass_basis=True))
    for item in observations:
        prediction = regression.predict(model, item)
        assert prediction.status == regression.ObservationStatus.Available and prediction.value is not None
        item.observed = prediction.value
    return truth, observations


def coordinates(start: Any) -> list[Any]:
    return [regression.coordinate(start, [("dispersion_energy_over_k", "s0", 1)], origin=150.0, scale=10.0, bounds=BOUNDS["dispersion_energy_over_k"]),
            regression.coordinate(start, [("k_ij", ("s0", "s1"), 1)], origin=0.0, scale=0.1, bounds=BOUNDS["k_ij"])]


def test_public_recovery_matches_the_native_rows(evidence: Any) -> None:
    """The same manufactured request through epcsaft.regression reproduces regression_checks R2-J2 (r4/*)."""
    _, observations = manufactured()
    start = neutral_binary(148.0, 0.0)
    result = regression.fit(start, coordinates(start), observations, controls=budget())
    assert result.status == "converged" and result.usable and result.trial_failures == 0
    assert result.residuals == result.weighted_residuals  # unit weights: one accepted objective evaluation
    same(evidence, "r4/recovery/epsilon", result.physical[0])
    same(evidence, "r4/recovery/kij", result.physical[1])
    same(evidence, "r4/recovery/initial_cost", result.initial_cost)
    assert result.final_cost <= 1.0e-20 and ROWS["r4/recovery/final_cost"] <= 1.0e-20
    for i, value in enumerate(result.predictions):
        same(evidence, f"r4/recovery/prediction/{i}", value)
    assert result.covariance.status == "available" and result.covariance.rank == 2
    for k in range(4):
        same(evidence, f"r4/recovery/covariance/{k}", result.covariance.optimizer[k])
        same(evidence, f"r4/recovery/physical_covariance/{k}", result.covariance.physical[k])
    # equivalence standard, CODING_STANDARDS.md, PR #311: variance uses the squared 1e-12
    # normalized-residual absolute floor at the noise-free fit, rather than relative error at zero.
    variance_reference = ROWS["r4/recovery/variance_factor"]
    same(evidence, "r4/recovery/variance_factor", result.covariance.variance_factor, variance_reference,
         atol=max(1e-24, 1e-10 * abs(variance_reference)), rtol=0.0)
    for k, value in enumerate(result.optimizer_jacobian):
        same(evidence, f"r4/recovery/optimizer_jacobian/{k}", value)
    assert math.isclose(result.physical[0], 152.0, rel_tol=1.0e-6) and abs(result.physical[1] + 0.02) <= 1.0e-6
    # The fitted document, renamed for the public constructors, reproduces the fitted predictions; the fit's
    # final evaluation was warm-started, a cold re-solve lands within the 1e-10 row criterion (observed 1e-13).
    fitted = Mixture(regression.fitted_parameters(start, result))
    for item, value in zip(observations, result.predictions):
        assert math.isclose(regression.predict(fitted, item).value or 0.0, value, rel_tol=1.0e-9)


def test_public_density_actions_match_the_native_rows(evidence: Any) -> None:
    """Mixture/Pure.pressure_actions with no observables is the #86 density-only request (d86/*)."""
    active = [(ParameterFamily.DispersionEnergyOverK, "epsilon"), (ParameterFamily.SegmentCount, "m"), (ParameterFamily.SigmaConstant, "sigma")]
    model = Pure(parameters((ETHANE,)), active_parameters=[_core.ActiveParameter(family, (0,)) for family, _ in active])
    directions = [PressureActionDirection(parameter=[ParameterActionDirection(family, (0,), 1.0)]) for family, _ in active]
    batch = model.pressure_actions(200.0, 217800.6949469248, phase="liquid", observables=(), directions=directions)
    assert batch.density_status == _core.PropertyActionStatus.Available and batch.results == []
    for k, (_, name) in enumerate(active):
        same(evidence, f"d86/parameter/{name}/fine", batch.density_action[k])
    thermal = model.pressure_actions(200.0, 217800.6949469248, phase="liquid", observables=(),
                                     directions=[PressureActionDirection(temperature=1.0), PressureActionDirection(pressure=1.0)])
    same(evidence, "d86/density_value", batch.molar_density)
    same(evidence, "d86/temperature/fine", thermal.density_action[0])
    same(evidence, "d86/pressure", thermal.density_action[1])
    absent = model.pressure_actions(200.0, 1.0e10, phase="liquid", observables=(), directions=directions)
    assert absent.density_status == _core.PropertyActionStatus.DensityUnavailable and absent.density_action == []
    assert absent.molar_density is None and "no stable density root" in absent.density_message and absent.density_convergence is None
    assert batch.density_message == "" and batch.density_convergence is None
    with pytest.raises(ParameterError):
        model.density_actions(200.0, 1000.0, observables=(), directions=[_core.DensityActionDirection()])


def test_global_born_scalar_chart_values_and_exact_map(evidence: Any) -> None:
    """Global dimensionless targets follow theta = coefficient*(origin + scale*u)."""
    p = aqueous_electrolyte(c_shell=1.0)
    coordinates = [regression.coordinate(p, [("born_c_dielectric", None, 2)], origin=0.25, scale=0.05),
                   regression.coordinate(p, [("ionic_region_relative_permittivity", None, 2)], origin=4.0, scale=0.5)]
    chart = regression.evaluate_coordinates(p, coordinates, (1.0, -5.0), ())
    for i, expected in enumerate((0.6, 3.0)):
        same(evidence, f"born-global-chart/physical/{i}", chart.physical[i], expected, atol=1e-15, rtol=0.0)
    for i, expected in enumerate((0.1, 0.0, 0.0, 1.0)):
        same(evidence, f"born-global-chart/map/{i}", chart.map[i], expected, atol=1e-15, rtol=0.0)
    reference = regression.evaluate_coordinates(p, coordinates, (5.0, 0.0), ())  # c=1, epsilon=8
    baseline = regression.evaluate_coordinates(p, coordinates, (-5.0, 0.0), ())  # c=0
    def born(native: Any) -> float:
        return Mixture(Parameters(native, p.component_ids, p.site_ids)).state(313.15, rho=45000.0, x=(0.98, 0.01, 0.01)).born
    # At fixed shell geometry/bulk permittivity, the correction is proportional to c*(1-1/epsilon).
    base = born(baseline.parameters)
    same(evidence, "born-global-chart/radial-ratio", (born(chart.parameters) - base) / (born(reference.parameters) - base), 16/35, atol=0.0, rtol=1e-12)


def test_point_reader_takes_the_fit_thermochemistry(evidence: Any) -> None:
    """One total c_p row read through evaluate_observations equals its fit prediction at the moved point."""
    p, thermal = parameters((ETHANE,)), record(p3_forms()["polynomial"])
    row = regression.isobaric_heat_capacity(p, 300.0, 1.0e5, {"s0": 1.0}, phase="vapor", observed=30.5, form="log_ratio")
    chart = [regression.coordinate(p, [("dispersion_energy_over_k", "s0", 1)], origin=191.42, scale=10.0)]
    controls = regression.FitControls(maximum_iterations=2, maximum_elapsed_time_seconds=60.0)
    fitted = regression.fit(p, chart, [row], thermochemistry=thermal, controls=controls)
    assert fitted.usable
    read = _core.evaluate_observations(p._native, chart, fitted.optimizer, [], [row], thermochemistry=thermal)
    same(evidence, "point-reader/total-cp", read.predictions[0], fitted.predictions[0], atol=0.0, rtol=1e-10, T_K=300.0, P_Pa=1.0e5)


def test_minimax_point_without_fit_budgets(evidence: Any) -> None:
    """The aqueous gamma point in native minimax-consistency needs no optimizer budget (#753)."""
    p = aqueous_electrolyte()
    chart = [regression.coordinate(p, [("k_ij", ("s0", "s1"), 1)],
                                  origin=-0.3, scale=0.2, bounds=(-1.0, 1.0))]
    row = regression.electrolyte(p, "mean_ionic_activity_coefficient", 298.15, 1e5, {"s0": 1.0},
        ions=("s1", "s2"), stoichiometry=(1, 1), molality=.1, observed=1.0, form="log_ratio")
    guard = regression.CheckGroup()
    guard.name, guard.rows, guard.maximum = "decades", [0], 1.0
    guard.metric = regression.CheckMetric.AbsoluteLog10Ratio
    spec = regression.MinimaxRequest()
    spec.targets, spec.group_scales = [0], [1.0]
    spec.q, spec.tolerance, spec.activity_distance, spec.margin = .02, 1e-10, 1e-6, .1
    refused = regression.fit(p, chart, [row], checks=[guard], minimax=spec,
                             controls=regression.FitControls(maximum_iterations=0, maximum_elapsed_time_seconds=0.0))
    assert refused.status == "invalid_setup"
    assert refused.message == "fit controls: maximum_iterations and maximum_elapsed_time_seconds must be set and positive"
    result = _core.evaluate_minimax(p._native, chart, [0.0], [], [row], [guard], None, spec)
    assert result.minimax is not None, result.message
    m = result.minimax
    assert m.availability == "values_only"
    assert m.point[0] == 0.0 and len(m.point) == 2
    assert m.lower == pytest.approx([-3.5, 0.0]) and m.upper == [6.5, math.inf]
    assert len(m.rows) == len(m.values) == 4 and len(m.jacobian) == 8
    assert [(item.kind, item.row) for item in m.rows] == [("target", 0)] * 2 + [("maximum", 0)] * 2
    assert all(math.isfinite(value) for values in (m.point, m.values, m.jacobian) for value in values)
    assert len(m.frozen_checks) == 1
    frozen = m.frozen_checks[0]
    assert (frozen.name, frozen.rows, frozen.maximum_threshold) == ("decades", [0], 1.0)
    assert math.isfinite(frozen.maximum) and frozen.status == "pass"
    evidence.observe("minimax/point-without-fit-budgets", point=m.point, lower=m.lower, upper=m.upper,
        scaled_values=m.values, jacobian=m.jacobian, frozen_maximum=frozen.maximum)


def test_addition_enthalpy_composes_its_endpoint_states(evidence: Any) -> None:
    """A nonreacting addition row equals (N_f h_f - N_i h_i)/n_add - h_feed read from the public states."""
    p = parameters((ETHANE, (30.07, 1.6069, 3.5206, 220.0)), pairs=[("k_ij", (0, 1), 0.0)])
    model = Mixture(p, thermochemistry=record(p3_forms()["polynomial"], p3_forms()["polynomial"]))
    T, P, n = 260.0, 3.0e6, 0.2
    initial, final = {"s0": 2.0, "s1": 0.2}, {"s0": 2.0 + n, "s1": 0.2}
    start, end = (Problem(phases=[Phase("L", kind="liquid")], T=T, P=P, feed=Amounts(feed)) for feed in (initial, final))
    row = regression.addition_enthalpy(p, start, end, {"s0": 1.0}, n, T, 1.0e5, phase="vapor", observed=0.0, form="difference", scale=1.0)

    def enthalpy(feed: dict[str, float], pressure: float, phase: str) -> float:
        total = sum(feed.values())
        molar = model.state(T, P=pressure, x=[feed["s0"] / total, feed["s1"] / total], phase=phase).total_enthalpy
        assert molar is not None
        return total * molar

    reference = (enthalpy(final, P, "liquid") - enthalpy(initial, P, "liquid")) / n - enthalpy({"s0": 1.0, "s1": 0.0}, 1.0e5, "vapor")
    same(evidence, "addition-enthalpy/composed", regression.predict(model, row).value, reference, atol=0.0, rtol=1e-10,
         T_K=T, P_Pa=P, n_add_mol=n)


def test_name_mapping_and_tie_map() -> None:
    """Names resolve to the owner's slot convention; ties map physical rows to optimizer columns."""
    p = associating_binary()
    kij = regression.active_parameter(p, "k_ij", ("s1", "s0"))
    assert (kij.family, list(kij.indices)) == (ParameterFamily.Kij, [1, 2])
    with pytest.raises(ParameterError):
        regression.density(p, 300.0, 1.0e5, {"s0": 1.0}, phase="liquid", observed=1.0, form="difference", mass_basis=False)
    for bad in (lambda: regression.coordinate(p, [("k_ij", ("s0", "s0"), 1)], origin=0.0, scale=1.0),
                lambda: regression.coordinate(p, [("segment_diameter", "s9", 1)], origin=0.0, scale=1.0),
                lambda: regression.coordinate(p, [("dielectric", "s0", 1)], origin=0.0, scale=1.0),
                lambda: regression.active_parameter(p, "dielectric", "s0"),
                lambda: regression.observation(p, "pressure", Problem(phases=Discover(), T=300.0, P=1.0e5, feed=Amounts({"s0": 1.0})), observed=1.0, form="log_ratio", mass_basis=False),
                lambda: regression.observation(p, "mole_fraction", bubble(300.0), observed=1.0, phase="X", form="log_ratio", mass_basis=False),
                lambda: regression.observation(p, "phase_fraction", bubble(300.0), observed=1.0, form="log_ratio", mass_basis=False),
                lambda: regression.density(p, 300.0, 1.0e5, {"s9": 1.0}, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False)):
        with pytest.raises(ParameterError):
            bad()


def test_predictions_and_failures_keep_their_meaning() -> None:
    """Values come from the state owner; a failing observation is typed, and an initial failure names its row."""
    p = neutral_binary()
    model = Mixture(p)
    item = regression.density(p, 180.0, 5.0e6, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=15000.0, form="log_ratio", mass_basis=False)
    prediction = regression.predict(model, item)
    state = model.state(180.0, P=5.0e6, x=(0.5, 0.5), phase="liquid")
    assert prediction.value == state.molar_density and prediction.unit == "mol/m^3"
    assert prediction.residual is not None and math.isclose(  # log(a/b) and log a - log b differ by an ulp
        prediction.residual, math.log(state.molar_density / 15000.0), rel_tol=1e-15, abs_tol=0.0)
    mass = regression.density(p, 180.0, 5.0e6, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=15000.0, form="log_ratio", mass_basis=True)
    mass_prediction = regression.predict(model, mass)
    assert mass_prediction.unit == "kg/m^3" and mass_prediction.value is not None
    assert math.isclose(mass_prediction.value, state.molar_density * (0.5 * 0.016043 + 0.5 * 0.03007), rel_tol=1.0e-12)
    # Public predict admits by name: an out-of-range phase never becomes a finite value.
    broken = regression.observation(p, "mole_fraction", bubble(180.0), observed=1.0, phase="L", species=("s0",), form="log_ratio", mass_basis=False)
    broken.phase = 2
    with pytest.raises(ParameterError):
        regression.predict(model, broken)
    # A value survives a nonfinite residual.
    overflow = regression.predict(model, regression.density(p, 180.0, 5.0e6, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=-1.0e308, form="difference", scale=1.0e-308, mass_basis=False))
    assert overflow.status == regression.ObservationStatus.Nonfinite and overflow.value == state.molar_density and overflow.residual is None
    none = regression.predict(model, regression.density(p, 200.0, 1.0e10, {"s0": 1.0}, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False))
    assert none.status == regression.ObservationStatus.DensityUnavailable and none.value is None and none.message
    assert prediction.branch == state.density_branch and prediction.stable_root_count == state.stable_root_count
    _, observations = manufactured()
    observations.append(regression.density(p, 200.0, 1.0e10, {"s0": 0.5, "s1": 0.5}, phase="liquid", observed=1.0, form="log_ratio", mass_basis=False))
    result = regression.fit(p, coordinates(p), observations, controls=budget())
    assert result.status == "initial_evaluation_failed" and not result.usable and result.predictions == []
    assert [(item.trial, item.observation, item.status) for item in result.failures] == [(1, 5, regression.ObservationStatus.DensityUnavailable)]
    refused = regression.fit(p, coordinates(p), [observations[0], observations[-1], observations[0]], weights=[0, 1, 1], controls=budget())
    assert refused.statuses == [regression.ObservationStatus.NotEvaluated, regression.ObservationStatus.DensityUnavailable, regression.ObservationStatus.NotEvaluated]
    single = regression.fit(p, coordinates(p), observations[:1], controls=budget())
    assert single.usable and single.covariance.status == "rank_deficient" and single.covariance.optimizer == []
    one = regression.FitControls(maximum_iterations=1, maximum_elapsed_time_seconds=600.0)
    unfinished = regression.fit(p, coordinates(p), observations[:5], controls=one)
    assert unfinished.status == "limit_reached" and unfinished.usable and unfinished.covariance.status == "not_converged"


def test_asl_record_rejects_document_less_parameters_early() -> None:
    """Parameters built without a parameter document have nothing to export: a precise error, not a cast failure."""
    start = neutral_binary()
    built = Mixture(Parameters(start._native, start.component_ids, start.site_ids))
    with pytest.raises(ValueError, match="ASL record: parameters built without a parameter document"):
        asl_record(built, bubble(180.0), Continuation(initial_step=0.1, minimum_step=0.1, max_steps=1, max_retries_per_step=0), [], [], [], 0, 0.0)


def test_native_held_electrolyte_recovery_checks_and_columns() -> None:
    """#219 compact F2/F5/F6: one salt, an identifiable held chart, final native metrics."""
    from tests.packets import parameters_json

    p = Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 4),
        components=("water", "sodium-cation", "chloride-anion"),
    )
    cs = [
        regression.coordinate(p, [("k_ij", pair, 1)], origin=0, scale=0.2, bounds=(-1, 1))
        for pair in (("sodium-cation", "water"), ("sodium-cation", "chloride-anion"))
    ]
    cs += [
        regression.coordinate(
            p,
            [("ion_suppression_coefficient", "sodium-cation", 1), ("ion_suppression_coefficient", "chloride-anion", 1)],
            origin=0,
            scale=2,
            bounds=(3, 9),
        )
    ]
    energy = regression.electrolyte(
        p,
        "solvation_gibbs_energy",
        298.15,
        1e5,
        {"water": 1},
        ions=("sodium-cation", "chloride-anion"),
        stoichiometry=(1, 1),
        molality=0,
        species="sodium-cation",
        observed=-381930,
        form="difference",
        scale=1,
    )
    holds = [
        regression.hydration_hold(
            energy,
            bracket=(1, 8),
            energy_tolerance=1e-6,
            diameter_tolerance=1e-10,
            maximum_iterations=50,
        )
    ]
    physical_inputs = [-0.25, 0.9, 6.01]
    physical_cs = []
    for coordinate in cs:
        unit = regression.Coordinate()
        unit.targets = coordinate.targets
        unit.origin, unit.scale = 0, 1
        physical_cs.append(unit)
    declared = regression.evaluate_coordinates(p, physical_cs, physical_inputs, [])
    native = regression.CoordinateSet()
    native.coordinates = cs
    truth = native.start(declared.parameters)
    chart = regression.evaluate_coordinates(p, cs, truth, holds)
    templates = [
        regression.electrolyte(
            p,
            kind,
            298.15,
            1e5,
            {"water": 1},
            ions=("sodium-cation", "chloride-anion"),
            stoichiometry=(1, 1),
            molality=m,
            observed=1,
            form="log_ratio",
        )
        for m in (0.1, 1, 4.5)
        for kind in ("mean_ionic_activity_coefficient", "osmotic_coefficient")
    ]
    # Manufacture values from the chart already evaluated above; its hold is solved.
    observations = _core.observations_at(chart.parameters, [], [], [], templates, [0] * 6)
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0,
                                      function_tolerance=1e-15, gradient_tolerance=1e-14, parameter_tolerance=1e-14)
    result = regression.fit(
        p,
        cs,
        observations,
        holds=holds,
        controls=controls,
        weights=[1 / (0.02 * math.log(10)) ** 2] * 6,
    )
    assert result.status == "converged", result.message
    assert max(abs(v) / math.log(10) for v in result.residuals) <= 1e-10
    for value, expected in zip(result.optimizer, truth):
        assert abs(value - expected) <= 1e-6
    assert abs(result.physical[-1] - chart.physical[-1]) <= 1e-7
    assert (
        abs(result.chart.holds[0].energy_error) <= 1e-6
        and result.chart.holds[0].correction <= 1e-10
    )
    stored = Parameters(chart.parameters, p.component_ids, p.site_ids)
    for ob in templates:
        ob.form = regression.ResidualForm.Difference
        ob.scale = 1
    exact_data = _core.observations_at(chart.parameters, [], [], [], templates, [0] * 6)
    exact = regression.fit(stored, cs, exact_data, holds=holds, controls=budget())
    assert exact.final_cost <= 1e-24
    for c in range(3):
        step = 1e-4 * max(abs(truth[c]), 1)
        plus, minus = list(truth), list(truth)
        plus[c] += step
        minus[c] -= step
        high = regression.evaluate_coordinates(p, cs, plus, holds)
        low = regression.evaluate_coordinates(p, cs, minus, holds)
        cd = (high.physical[-1] - low.physical[-1]) / (2 * step)
        assert abs(chart.map[-3 + c] - cd) <= 2e-5 * max(abs(v) for v in chart.map[-3:]) + 1e-7
        yp = _core.observations_at(high.parameters, [], [], [], templates, [0] * 6)
        ym = _core.observations_at(low.parameters, [], [], [], templates, [0] * 6)
        column = exact.optimizer_jacobian[c::3]
        for i, value in enumerate(column):
            assert (
                abs(value - (yp[i].observed - ym[i].observed) / (2 * step))
                <= 2e-5 * max(abs(v) for v in column) + 1e-7
            )
    checks_data = _core.observations_at(
        chart.parameters, [], [], [], observations[:2], [0.01, -0.03]
    )
    group = regression.CheckGroup()
    group.name = "known-decades"
    group.rows = [6, 7]
    group.metric = regression.CheckMetric.AbsoluteLog10Ratio
    group.maximum = 0.02
    group.mean = 0.025
    final = regression.fit(
        stored,
        cs,
        observations + checks_data,
        holds=holds,
        checks=[group],
        controls=budget(),
        weights=[1] * 6 + [0, 0],
    )
    check = final.checks[0]
    assert final.status == "converged" and final.acceptance == "fail"
    assert abs(check.maximum - 0.03) <= 1e-10 and abs(check.mean - 0.02) <= 1e-10
    assert check.argmax == 7 and check.failed_rows == [7]
    sequence = regression.fit(
        stored,
        cs,
        exact_data,
        holds=holds,
        controls=budget(),
        starts=[truth, truth],
        shared_iterations=20,
        shared_seconds=20,
    )
    assert len(sequence) == 2 and all(candidate.status == "converged" for candidate in sequence)
    assert all(candidate.final_cost <= 1e-24 for candidate in sequence)


def test_mixture_osmotic_coefficient_matches_one_salt_and_solvent_activity() -> None:
    """#575: the full-composition osmotic row equals the one-salt row and the direct solvent activity, value and column."""
    from tests.packets import parameters_json

    names = ("water", "sodium-cation", "magnesium-cation", "chloride-anion")
    p = Parameters.from_json(parameters_json("figiel-2025-chloride-activity-refit", 2), components=names)
    chart = regression.CoordinateSet()
    chart.coordinates = [regression.coordinate(p, [target], origin=0, scale=1)
                         for target in (("k_ij", ("magnesium-cation", "water"), 1), ("debye_huckel_diameter", "sodium-cation", 1))]
    model = Mixture(p, active_parameters=chart.active())
    salt = regression.electrolyte(p, "osmotic_coefficient", 298.15, 1e5, {"water": 1}, ions=("magnesium-cation", "chloride-anion"),
                                  stoichiometry=(1, 2), molality=2.0, observed=1, form="difference", scale=1)
    x = _core.formula_composition([1, 0, 0, 0], p._native, (2, 3), (1, 2), 2.0)
    amounts = (1 / p.molar_masses[0], 1.5, 0.75, 3.0)  # Rard-Miller-like NaCl-MgCl2 composition, I = 3.75 mol/kg
    mixed = [dict(zip(names, x)), dict(zip(names, (n / sum(amounts) for n in amounts)))]
    rows = [regression.mixture_osmotic_coefficient(p, 298.15, 1e5, y, observed=1, scale=1) for y in mixed]
    direct = model.reference(298.15, 1e5, solvent_fractions=(1, 0, 0, 0)).solvent_activity(x=list(mixed[1].values()))
    expected = [(regression.predict(model, salt).value, regression.differentiate(model, salt).value),
                (direct.osmotic_coefficient, direct.osmotic_parameter_derivative)]
    for row, (value, column) in zip(rows, expected):
        assert regression.predict(model, row).value == pytest.approx(value, rel=1e-12)
        assert regression.differentiate(model, row).value == pytest.approx(column, rel=1e-10, abs=1e-14)
    with pytest.raises(ParameterError, match="non-electroneutral"):
        regression.predict(model, regression.mixture_osmotic_coefficient(p, 298.15, 1e5, {"water": 0.9, "sodium-cation": 0.1},
                                                                         observed=1, scale=1))


def test_density_options_reach_fit_worker(evidence: Any) -> None:
    """Two EOS state groups distinguish lower-model options from restored defaults (#486 R3)."""
    source = parameters((ETHANE,))
    model = Pure(source)
    observations = [regression.density(source, 180.0, P, {"s0": 1.0}, phase="liquid", observed=1.0,
        form="log_ratio", mass_basis=False) for P in (1e5, 1e6)]
    for observation in observations:
        value = regression.predict(model, observation).value
        assert value is not None
        observation.observed = value
    permittivity = regression.Observation()
    permittivity.kind = regression.ObservationKind.BulkRelativePermittivity
    permittivity.temperature, permittivity.pressure = 180.0, 1e5
    permittivity.composition, permittivity.branch = [1.0], _core.DensityBranch.Liquid
    permittivity.observed, permittivity.scale = 0.0, 1.0
    observations.append(permittivity)
    check = regression.CheckGroup()
    check.name, check.rows, check.maximum = "permittivity", [2], 1e-12
    check.baseline = True
    chart = [regression.coordinate(source, [("dispersion_energy_over_k", "s0", 1.0)], origin=191.42,
        scale=10.0, bounds=(180.0, 210.0))]
    def run(workers: int, options: dict[str, float]) -> Any:
        return regression.fit(source, chart, observations, weights=[1.0, 1.0, 0.0], checks=[check], baseline=source, controls=regression.FitControls(
            maximum_iterations=1, maximum_elapsed_time_seconds=60, workers=workers,
            gradient_tolerance=1e100, solver_options=options))
    default = run(1, {})
    single, parallel = (run(workers, {"DENSITY_ROOT_RESIDUAL": 1e-5}) for workers in (1, 2))
    assert default.usable and single.usable and parallel.usable, (default.message, single.message, parallel.message)
    # Row 1 belongs to the child; the default restoration exceeds the 1e-10 equivalence limit.
    assert abs(single.predictions[1] / default.predictions[1] - 1.0) > 1e-10
    evidence.observe("worker/default-control", default=default.predictions, nondefault=single.predictions,
        solver_options=single.solver_options, workers=2)
    assert single.optimizer == parallel.optimizer and single.solver_options == parallel.solver_options
    assert single.status == parallel.status and single.statuses == parallel.statuses
    assert [(f.status, f.message) for f in single.failures] == [(f.status, f.message) for f in parallel.failures] == []
    for field in ("predictions", "residuals", "weighted_residuals", "physical_jacobian", "optimizer_jacobian"):
        for i, (a, b) in enumerate(zip(getattr(single, field), getattr(parallel, field), strict=True)):
            assert math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12), (field, i, a, b)
            evidence.compare(f"worker/{field}/{i}", model=b, reference=a, unit="field units",
                criterion={"rule": "symmetric_close", "rtol": 1e-10, "atol": 1e-12})
    for a, b in zip(single.checks, parallel.checks, strict=True):
        assert a.status == b.status == "pass" and a.errors == b.errors
        assert a.baseline_maximum == b.baseline_maximum and a.failed_rows == b.failed_rows
        evidence.observe("worker/final-check", maximum=b.maximum, mean=b.mean, baseline=b.baseline_maximum, status=b.status)
    point = regression.differentiate(source, chart, observations[:2], optimizer=single.optimizer,
        solver_options={"DENSITY_ROOT_RESIDUAL": 1e-5})
    assert point.predictions == single.predictions[:2]
    with pytest.raises(ParameterError, match="outside owner scope"):
        regression.differentiate(source, chart, observations[:2], optimizer=single.optimizer,
            solver_options={"FIT_FUNCTION_TOLERANCE": 1e-8})
    assert point.solver_options == single.solver_options or all(
        entry["name"].startswith("FIT_") or entry in point.solver_options for entry in single.solver_options)


def test_solver_option_defaults_and_model_isolation() -> None:
    source = parameters((ETHANE,))
    defaults = {entry.name: entry.value for entry in _core.tolerances()}
    names = ("ASSOCIATION_RESIDUAL", "DENSITY_ROOT_RESIDUAL", "EQUILIBRIUM_ROW_RESIDUAL")
    def state(model: Pure) -> Any:
        return model.state(180.0, P=1e5, phase="liquid")
    records = [state(Pure(source, solver_options=options)) for options in (None, {}, {name: defaults[name] for name in names})]
    for result in records[1:]:
        assert result.molar_density == records[0].molar_density and result.pressure == records[0].pressure
        assert result.residual_helmholtz == records[0].residual_helmholtz
        assert result.density_convergence.evaluations == records[0].density_convergence.evaluations
        assert result.density_convergence.iterations == records[0].density_convergence.iterations
    ordinary = Pure(source)
    loose = Pure(source, solver_options={"DENSITY_ROOT_RESIDUAL": 1e-5})
    saved = state(loose).density_convergence
    assert state(ordinary).molar_density == records[0].molar_density
    assert state(loose).density_convergence.criterion.value == 1e-5
    del loose
    assert saved.criterion.value == 1e-5 and saved.criterion.upstream is None
    assert {entry.name: entry.value for entry in _core.tolerances()} == defaults
