"""Native B3 consistency and B4 guarded-minimax qualification."""
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path

from epcsaft import Parameters, _core, regression as r
from declarations import CATIONS, PACKET, ROOT, TARGETS, coordinates, holds, parameters, row
from run import manufactured, native_start
from checks import group

OUT = Path(__file__).resolve().parents[1] / "results"
COMPONENTS_B = ("water", "methanol", *CATIONS, "chloride-anion")
MOLALITIES = (0.1, 1.0, 4.5)
EPSILON_MOLALITIES = (0.1, 1.0)
Q = 0.02
TOLERANCE = 1e-10
ACTIVITY_DISTANCE = 1e-6
MARGIN = 0.1
COMPARE_RELATIVE = 1e-12
COMPARE_ABSOLUTE_AT_ZERO = 1e-12
KKT_LIMIT = 1e-8
records = []
failures = []


def record(case, quantity, exact, reference=None, tolerance=None, unit="", **labels):
    if exact is None or reference is None:
        error, passed = "", ""
    else:
        error = abs(exact - reference)
        passed = math.isfinite(error) and error <= tolerance
        if not passed:
            failures.append({"case": case, "quantity": quantity,
                             "message": "numerical comparison exceeded its declared tolerance"})
    records.append({"case": case, "quantity": quantity, "exact": exact, "reference": reference,
                    "error": error, "tolerance": tolerance, "unit": unit, "passed": passed, **labels})
    if passed is False:
        raise AssertionError((case, quantity, exact, reference, error, tolerance))


def compare(case, quantity, value_b, value_a, unit=""):
    if value_b is None or value_a is None:
        if value_b is not None or value_a is not None:
            raise AssertionError((case, quantity, value_b, value_a))
        return
    b, a = float(value_b), float(value_a)
    if math.isnan(b) and math.isnan(a):
        records.append({"case": case, "quantity": quantity, "exact": "nan", "reference": "nan",
                        "error": 0.0, "tolerance": COMPARE_ABSOLUTE_AT_ZERO, "unit": unit, "passed": True})
        return
    tolerance = COMPARE_RELATIVE * abs(a) if a != 0.0 else COMPARE_ABSOLUTE_AT_ZERO
    record(case, quantity, b, a, tolerance, unit)


def require(condition, message, **details):
    if not condition:
        failures.append({"case": details.pop("case", "qualification"),
                         "quantity": details.pop("quantity", "condition"), "message": message})
        raise AssertionError((message, details))


def minimax_spec(targets, scales):
    spec = r.MinimaxRequest()
    spec.targets = list(targets)
    spec.group_scales = list(scales)
    spec.q, spec.tolerance = Q, TOLERANCE
    spec.activity_distance, spec.margin = ACTIVITY_DISTANCE, MARGIN
    return spec


def transfer_observation(p, solvent_from="water", solvent_to="methanol"):
    return row(p, "transfer_gibbs_energy", "sodium-cation", x={solvent_to: 1.0},
               reference={solvent_from: 1.0}, form="difference", observed=0.0)


def composition_at_molality(p, molality):
    solvent = [1.0 if name == "water" else 0.0 for name in p.component_ids]
    ions = (p.component_ids.index("sodium-cation"), p.component_ids.index("chloride-anion"))
    x = _core.formula_composition(solvent, p._native, ions, (1, 1), molality)
    return dict(zip(p.component_ids, x))


def source_and_targets(p, cs, hs, source):
    templates = [row(p, kind, ion, m=m)
                 for ion in CATIONS
                 for kind in ("mean_ionic_activity_coefficient", "osmotic_coefficient")
                 for m in MOLALITIES]
    errors = [0.005 if item.kind == r.ObservationKind.MeanIonicActivityCoefficient
              and item.ion_indices[0] == p.component_ids.index("sodium-cation") else 0.0
              for item in templates]
    target_data = manufactured(p, cs, source, hs, templates, errors)
    epsilon_templates = [row(p, "bulk_relative_permittivity", x=composition_at_molality(p, m), form="difference")
                         for m in EPSILON_MOLALITIES]
    epsilon_data = manufactured(p, cs, source, hs, epsilon_templates, [0.1, 0.1])
    energy_rows = [hold.reference for hold in hs]

    forward = transfer_observation(p, "water", "methanol")
    reverse = transfer_observation(p, "methanol", "water")
    source_values = _core.evaluate_observations(p._native, cs, source, hs, [forward, reverse])
    candidates = [(source_values.predictions[i], item) for i, item in enumerate((forward, reverse))
                  if source_values.statuses[i] == r.ObservationStatus.Available and source_values.predictions[i] > 0.0]
    require(bool(candidates), "neither solvent transfer direction produced a positive native energy",
            case="B3", quantity="transfer direction")
    transfer_prediction, transfer = candidates[0]
    # A tenfold synthetic datum stays positive and creates a report-only ratio
    # failure; the shared native Difference row remains available for columns.
    transfer_data = manufactured(p, cs, source, hs, [transfer], [0.9 * transfer_prediction])
    transfer_detail = {"solvent_from": "water" if transfer.composition[p.component_ids.index("methanol")] > 0 else "methanol",
                       "solvent_to": "methanol" if transfer.composition[p.component_ids.index("methanol")] > 0 else "water",
                       "composition": {name: value for name, value in zip(p.component_ids, transfer.composition) if value > 0.0},
                       "reference_composition": {name: value for name, value in zip(p.component_ids, transfer.reference_composition) if value > 0.0}}
    return target_data + epsilon_data + energy_rows + transfer_data, transfer_prediction, transfer_detail


def indices_for(data):
    ntarget = len(CATIONS) * 2 * len(MOLALITIES)
    nepsilon = len(EPSILON_MOLALITIES)
    nenergy = len(TARGETS)
    epsilon = list(range(ntarget, ntarget + nepsilon))
    energy = list(range(ntarget + nepsilon, ntarget + nepsilon + nenergy))
    transfer = ntarget + nepsilon + nenergy
    require(len(data) == transfer + 1, "unexpected B4 observation layout", case="B4", quantity="row count")
    return list(range(ntarget)), epsilon, energy, transfer


def b4_checks(indices):
    _, epsilon, energy, transfer = indices
    diff = r.CheckMetric.AbsoluteDifference
    log = r.CheckMetric.AbsoluteLog10Ratio
    return [group("aqueous_permittivity", epsilon, diff, 0.1, 0.1, baseline=True),
            group("hydration_energy_report_only", energy, diff, 10.0, report_only=True),
            group("transfer_energy_log_report_only", [transfer], log, 0.05, report_only=True)]


def metric_name(metric):
    return "absolute_log10_ratio" if metric == r.CheckMetric.AbsoluteLog10Ratio else "absolute_difference"


def b3_scales(p):
    cs = coordinates(p)
    source_theta = native_start(p, cs)
    scales = (0.2, 0.4, 0.5, 0.25, 0.3, 0.6, 2.0)
    for coordinate, scale in zip(cs, scales):
        coordinate.scale = scale
    return cs, scales, source_theta


def point_checks_equal(case, name, b, a):
    require((b.name, b.unit, b.status, b.rows, b.failed_rows, b.argmax, b.report_only)
            == (a.name, a.unit, a.status, a.rows, a.failed_rows, a.argmax, a.report_only),
            "native check identity/status differs", case=case, quantity=name)
    for field in ("maximum", "mean", "maximum_threshold", "mean_threshold", "baseline_maximum", "baseline_mean"):
        compare(case, f"{name}/{field}", getattr(b, field), getattr(a, field), b.unit)
    require(len(b.errors) == len(a.errors), "native check row count differs", case=case, quantity=name)
    for row_index, (b_error, a_error) in enumerate(zip(b.errors, a.errors)):
        compare(case, f"{name}/row-{b.rows[row_index]}/error", b_error, a_error, b.unit)


def b3_point(name, p, cs, u, hs, data, baseline, checks, target_indices, scales):
    print(f"B3 {name}: native one-point values", flush=True)
    chart = r.evaluate_coordinates(p, cs, u, hs)
    stored = Parameters(chart.parameters, p.component_ids, p.site_ids)
    spec = minimax_spec(target_indices, scales)
    b = _core.evaluate_minimax(p._native, cs, list(u), hs, data, checks, baseline._native, spec)
    require(b.minimax is not None and b.minimax.availability == "values_only",
            "native minimax point columns unavailable", case=name, quantity="availability",
            status=b.status, native_message=b.message)
    a = _core.evaluate_observations(p._native, cs, list(u), hs, data)
    active = sorted(set(target_indices + [i for check in checks if not check.report_only for i in check.rows]))
    physical_width = len(chart.physical)
    optimizer_width = len(cs)
    for i in active:
        require(a.statuses[i] == r.ObservationStatus.Available and b.statuses[i] == r.ObservationStatus.Available,
                "an enforced B3 row is unavailable", case=name, quantity=f"row-{i}-status")
        compare(name, f"row-{i}/prediction", b.predictions[i], a.predictions[i], "native value")
        for c in range(physical_width):
            compare(name, f"row-{i}/physical-column-{c}", b.physical_jacobian[i * physical_width + c],
                    a.gradient[i * physical_width + c], "value per physical parameter unit")
        for c in range(optimizer_width):
            compare(name, f"row-{i}/optimizer-column-{c}", b.optimizer_jacobian[i * optimizer_width + c],
                    a.optimizer_gradient[i * optimizer_width + c], "value per optimizer unit")
    if name == "source":
        target_errors = []
        for i in target_indices:
            error = math.log10(a.predictions[i] / data[i].observed)
            expected = 0.005 if data[i].kind == r.ObservationKind.MeanIonicActivityCoefficient \
                and data[i].ion_indices[0] == p.component_ids.index("sodium-cation") else 0.0
            compare(name, f"target-row-{i}/manufactured-error", error, expected, "decade")
            target_errors.append(abs(error))
        record("B4-source-witness", "target-maximum", max(target_errors), 0.005, 1e-12, "decade")

    frozen_a = r.evaluate_checks(baseline, data, checks, baseline)
    checks_a = r.evaluate_checks(stored, data, checks, baseline)
    require(len(b.minimax.frozen_checks) == len(frozen_a) and len(b.checks) == len(checks_a),
            "native check group count differs", case=name, quantity="check count")
    for i, (b_frozen, a_frozen) in enumerate(zip(b.minimax.frozen_checks, frozen_a)):
        point_checks_equal(name, f"frozen-{i}", b_frozen, a_frozen)
    for i, (b_check, a_check) in enumerate(zip(b.checks, checks_a)):
        point_checks_equal(name, f"final-{i}", b_check, a_check)

    for c, (value, expected) in enumerate(zip(b.chart.physical[:optimizer_width], chart.physical[:optimizer_width])):
        compare(name, f"chart/physical-{c}", value, expected, "native parameter unit")
    for i, hold in enumerate(b.chart.holds):
        a_hold = chart.holds[i]
        for field, unit in (("energy_error", "J/mol"), ("correction", "Angstrom"), ("uncertainty", "Angstrom")):
            compare(name, f"hold-{i}/{field}", getattr(hold, field), getattr(a_hold, field), unit)
    for i, value in enumerate(b.minimax.values):
        record(name, f"constraint-{i}/scaled-value", value, unit="scaled constraint")
    return b


def b3_fitted_point(name, p, cs, hs, result, data, baseline, checks, active):
    print(f"B3 {name}: replay fitted coordinates with common cold initialization", flush=True)
    replay = b3_point(name + "-cold", p, cs, result.optimizer, hs, data, baseline,
                      checks, indices_for(data)[0], [1.0, 1000.0, 0.02])
    a = _core.evaluate_observations(result.parameters, cs, result.optimizer, hs, data)
    physical_width = len(result.physical)
    optimizer_width = len(cs)
    for i in active:
        require(a.statuses[i] == r.ObservationStatus.Available
                and result.statuses[i] == r.ObservationStatus.Available,
                "a returned target or enforced guard is unavailable", case=name, quantity=f"row-{i}-status")
        record(name, f"row-{i}/fitted-minus-cold", result.predictions[i] - replay.predictions[i],
               unit="native value; density initialization diagnostic")
        for c in range(physical_width):
            compare(name, f"row-{i}/physical-column-{c}",
                    result.physical_jacobian[i * physical_width + c],
                    a.gradient[i * physical_width + c], "value per physical parameter unit")
        for c in range(optimizer_width):
            compare(name, f"row-{i}/optimizer-column-{c}",
                    result.optimizer_jacobian[i * optimizer_width + c],
                    a.optimizer_gradient[i * optimizer_width + c], "value per optimizer unit")
    stored = Parameters(result.parameters, p.component_ids, p.site_ids)
    frozen = r.evaluate_checks(baseline, data, checks, baseline)
    final = r.evaluate_checks(stored, data, checks, baseline)
    require(len(result.minimax.frozen_checks) == len(frozen) and len(result.checks) == len(final),
            "returned-point native check group count differs", case=name, quantity="check count")
    for i, (b, ref) in enumerate(zip(result.minimax.frozen_checks, frozen)):
        point_checks_equal(name, f"frozen-{i}", b, ref)
    for i, (b, ref) in enumerate(zip(result.checks, final)):
        point_checks_equal(name, f"final-{i}", b, ref)
    for i, (b, ref) in enumerate(zip(result.physical, a.chart.physical)):
        compare(name, f"chart/physical-{i}", b, ref, "native parameter unit")
    for i, (b, ref) in enumerate(zip(result.chart.holds, a.chart.holds)):
        for field, unit in (("energy_error", "J/mol"), ("correction", "Angstrom"), ("uncertainty", "Angstrom")):
            compare(name, f"hold-{i}/{field}", getattr(b, field), getattr(ref, field), unit)


def record_constraints(case, result, start):
    m = result.minimax
    for c, value in enumerate(start):
        record(case, f"start/{c}", value, unit="optimizer coordinate")
    for c, value in enumerate(result.optimizer):
        record(case, f"fit/coordinate-{c}", value, unit="optimizer coordinate")
    for h, hold in enumerate(result.chart.holds):
        for field, unit in (("energy_error", "J/mol"), ("correction", "Angstrom"), ("uncertainty", "Angstrom"), ("derivative", "J/(mol Angstrom)")):
            record(case, f"hold-{h}/{field}", getattr(hold, field), unit=unit)
    for h, value in enumerate(m.hold_multipliers):
        record(case, f"hold-{h}/multiplier", value, unit="decade/(J/mol)", status=m.hold_causes[h])
        for c, derivative in enumerate(result.chart.holds[h].parameter_derivative):
            record(case, f"hold-{h}/partial-{c}", derivative, unit="J/mol per physical parameter unit")
    for c, value in enumerate(m.full_stationarity):
        unit = "decade per optimizer unit" if c < len(result.optimizer) else "decade/Angstrom"
        record(case, f"full-stationarity/{c}", value, 0.0, KKT_LIMIT, unit)
    for c, value in enumerate(result.physical):
        unit = "Angstrom" if result.chart.active[c].family.name == "BornDiameter" else "1"
        record(case, f"fit/physical-{c}", value, unit=unit)
    for c, value in enumerate(result.chart.map):
        record(case, f"chart-map/{c}", value, unit="physical parameter unit per optimizer unit")
    for field, unit in (("t", "decade"), ("target_maximum", "decade"), ("gap", "decade"),
                        ("primal", "scaled"), ("stationarity", "scaled"), ("dual_sign", "scaled"),
                        ("complementarity", "scaled")):
        record(case, f"minimax/{field}", getattr(m, field), unit=unit)
    for c, _ in enumerate(m.point):
        for field, array, unit in (("point", m.point, "optimizer/auxiliary"), ("lower", m.lower, "optimizer/auxiliary"),
                                   ("upper", m.upper, "optimizer/auxiliary"), ("z_lower", m.z_lower, "scaled multiplier"),
                                   ("z_upper", m.z_upper, "scaled multiplier"),
                                   ("physical_z_lower", m.physical_z_lower, "physical multiplier"),
                                   ("physical_z_upper", m.physical_z_upper, "physical multiplier"),
                                   ("lower_distance", m.lower_distances, "scaled variable"),
                                   ("upper_distance", m.upper_distances, "scaled variable")):
            record(case, f"variable-{c}/{field}", array[c], unit=unit)
    for i, row_data in enumerate(m.rows):
        for field, array, unit in (("scaled_value", m.values, "scaled constraint"),
                                   ("physical_value", m.physical_values, "native constraint unit"),
                                   ("scaled_multiplier", m.multipliers, "scaled multiplier"),
                                   ("physical_multiplier", m.physical_multipliers, "physical multiplier"),
                                   ("root_error", m.root_errors, "scaled constraint")):
            if i < len(array):
                record(case, f"constraint-{i}/{field}", array[i], unit=unit,
                       kind=row_data.kind, row=row_data.row, group=row_data.group, metric=metric_name(row_data.metric),
                       sign=row_data.sign, threshold=row_data.threshold, scale=row_data.scale,
                       auxiliary=row_data.auxiliary, members=":".join(map(str, row_data.members)))
    for i, check in enumerate(result.checks):
        record(case, f"check-{i}/maximum", check.maximum, unit=check.unit, status=check.status)
        record(case, f"check-{i}/mean", check.mean, unit=check.unit, status=check.status)
        if check.maximum_threshold is not None:
            record(case, f"check-{i}/maximum_threshold", check.maximum_threshold, unit=check.unit)
        if check.mean_threshold is not None:
            record(case, f"check-{i}/mean_threshold", check.mean_threshold, unit=check.unit)
        for row_index, error in zip(check.rows, check.errors):
            record(case, f"check-{i}/row-{row_index}/error", error, unit=check.unit, status=check.status)


def make_case():
    p = parameters(COMPONENTS_B)
    cs = coordinates(p)
    hs = holds(p)
    source = native_start(p, cs)
    data, transfer_prediction, transfer_detail = source_and_targets(p, cs, hs, source)
    target_indices, epsilon_indices, energy_indices, transfer_index = indices_for(data)
    baseline_chart = r.evaluate_coordinates(p, cs, source, hs)
    baseline = Parameters(baseline_chart.parameters, p.component_ids, p.site_ids)
    return p, cs, hs, source, data, baseline, (target_indices, epsilon_indices, energy_indices, transfer_index), transfer_prediction, transfer_detail


def run_b3(p, hs, data, baseline, indices):
    b3_cs, scales, source_theta = b3_scales(p)
    source_u = native_start(p, b3_cs)
    truth_theta = [value + (0.05 if abs(value + 0.05) <= 0.95 else -0.05) if i % 2 == 0
                   else value + (0.1 if abs(value + 0.1) <= 0.95 else -0.1)
                   for i, value in enumerate(source_theta[:6])] + [source_theta[6] - 1.0]
    truth_chart = r.evaluate_coordinates(p, coordinates(p), truth_theta, hs)
    truth_u = native_start(Parameters(truth_chart.parameters,p.component_ids,p.site_ids),b3_cs)
    targets, epsilon, energy, transfer = indices
    diff, log = r.CheckMetric.AbsoluteDifference, r.CheckMetric.AbsoluteLog10Ratio
    checks = [group("aqueous_permittivity", epsilon, diff, 0.1, 0.1, baseline=True),
              group("hydration_energy", energy, diff, 10.0),
              group("transfer_energy", [transfer], diff, 0.0, baseline=True),
              group("transfer_energy_log_report_only", [transfer], log, 0.05, report_only=True)]
    for index, value in enumerate(source_theta):
        record("B3-input", f"physical-source-{index}", value, unit="native parameter unit")
    for index, value in enumerate(scales):
        record("B3-input", f"affine-scale-{index}", value, unit="native parameter unit per optimizer unit")
    point_results = []
    for name, u in (("source", source_u), ("F2-manufactured-truth", truth_u)):
        result = b3_point(name, p, b3_cs, u, hs, data, baseline, checks, targets,
                          [1.0, 1000.0, 1000.0, 0.02])
        point_results.append(result)
    return b3_cs, scales, checks, point_results


def run_b4(p, cs, hs, source, data, baseline, indices):
    targets, epsilon, _, _ = indices
    checks = b4_checks(indices)
    weights = [1.0 / (Q * math.log(10.0)) ** 2] * len(targets)
    weights += [0.0] * (len(data) - len(targets))
    for index in epsilon:
        require(weights[index] == 0.0, "aqueous guard unexpectedly has positive Ceres weight",
                case="B4", quantity=f"weight-row-{index}")
    for index, value in enumerate(weights):
        record("B4-input", f"ceres-weight-row-{index}", value, unit="Ceres weight")
    controls = r.FitControls(maximum_iterations=300, maximum_elapsed_time_seconds=120.0)
    shifted = list(source)
    shifted[2] = math.copysign(max(0.0, abs(shifted[2]) - 0.01), shifted[2]) if shifted[2] else 0.01
    shifted[6] += 0.1
    spec = minimax_spec(targets, [1.0, 1000.0, 0.02])
    print("B4: native minimax, two declared starts", flush=True)
    started = time.perf_counter()
    results = r.fit(p, cs, data, weights=weights, controls=controls, holds=hs, checks=checks, baseline=baseline,
                    starts=[source, shifted], shared_iterations=600, shared_seconds=240.0, minimax=spec)
    elapsed = time.perf_counter() - started
    require(isinstance(results, list) and len(results) == 2, "both declared starts were not attempted",
            case="B4", quantity="starts", count=len(results) if isinstance(results, list) else 1)
    for index, (result, start) in enumerate(zip(results, (source, shifted))):
        case = f"B4-start-{index}"
        minimax_result = result.minimax
        print(f"{case}: status={result.status} solver={minimax_result.solver_status if minimax_result else 'unavailable'} "
              f"iterations={result.iterations} message={result.message}", flush=True)
        require(minimax_result is not None, "minimax record unavailable", case=case, quantity="minimax")
        record_constraints(case, result, start)
        m = minimax_result
        require(m.availability == "available", "minimax candidate or KKT operands unavailable",
                case=case, quantity="minimax availability", status=result.status, native_message=result.message)
        require(result.status == "converged" and m.feasible and m.stationary,
                "native minimax did not converge to a feasible stationary point", case=case,
                quantity="termination", status=result.status, feasible=m.feasible, stationary=m.stationary)
        require(m.t <= 0.005 + 1e-8 and abs(m.gap) <= 1e-8,
                "B4 target epigraph does not match the recomputed maximum", case=case,
                quantity="target maximum", t=m.t, maximum=m.target_maximum, gap=m.gap)
        require(max(m.primal, m.stationarity, m.dual_sign, m.complementarity) <= KKT_LIMIT,
                "B4 scaled KKT norm exceeded its accepted limit", case=case, quantity="KKT")
        require(max((abs(value) for value in m.root_errors), default=0.0) <= 1e-9,
                "B4 estimated held-root action error exceeded its limit", case=case, quantity="root estimate")
        require(result.checks[0].status == "pass", "aqueous permittivity guard failed", case=case,
                quantity="aqueous permittivity", status=result.checks[0].status)
        require(result.checks[0].maximum <= 0.1 and result.checks[0].mean <= 0.1,
                "aqueous permittivity check exceeds its frozen source thresholds", case=case,
                quantity="aqueous permittivity")
        hydration = result.checks[1]
        require(hydration.maximum <= 10.0, "stored hydration energy check exceeds 10 J/mol",
                case=case, quantity="hydration energy", maximum=hydration.maximum)
        for h, hold in enumerate(result.chart.holds):
            require(abs(hold.energy_error) <= 1e-6 and max(abs(hold.correction), abs(hold.uncertainty)) <= 1e-10,
                    "stored dependent Born root misses the D5 criteria", case=case,
                    quantity=f"hold-{h}", energy=hold.energy_error, correction=hold.correction,
                    uncertainty=hold.uncertainty)
    return results, checks, elapsed


def root_refinement(p, cs, hs, result, data, baseline, targets, checks, case):
    refined = [r.hydration_hold(hold.reference, bracket=(hold.lower, hold.upper), energy_tolerance=1e-7,
                                diameter_tolerance=1e-11, maximum_iterations=hold.maximum_iterations)
               for hold in hs]
    spec = minimax_spec(targets, [1.0, 1000.0, 0.02])
    original = _core.evaluate_minimax(p._native, cs, result.optimizer, hs, data, checks, baseline._native, spec)
    fine = _core.evaluate_minimax(p._native, cs, result.optimizer, refined, data, checks, baseline._native, spec)
    require(original.minimax is not None and fine.minimax is not None
            and original.minimax.availability == fine.minimax.availability == "values_only",
            "held-root refinement point evaluation unavailable", case="B4-refinement", quantity="availability")
    require(len(original.minimax.values) == len(fine.minimax.values)
            and len(original.minimax.jacobian) == len(fine.minimax.jacobian),
            "held-root refinement changed constraint dimensions", case="B4-refinement", quantity="dimensions")
    for i, (a, b) in enumerate(zip(original.minimax.values, fine.minimax.values)):
        record(case, f"constraint-{i}/scaled-value-change", b, a, 1e-9, "scaled constraint")
    for i, (a, b) in enumerate(zip(original.minimax.jacobian, fine.minimax.jacobian)):
        record(case, f"constraint-column-{i}/change", b, a, 1e-9, "scaled derivative")
    return len(original.minimax.values), len(original.minimax.jacobian)


def no_report_only_effect(p, cs, hs, u, data, baseline, targets, checks, scales):
    with_report = _core.evaluate_minimax(p._native, cs, list(u), hs, data, checks, baseline._native,
                                         minimax_spec(targets, scales))
    without_report = [item for item in checks if not item.report_only]
    without_scales = [scale for scale, item in zip(scales, checks) if not item.report_only]
    no_report = _core.evaluate_minimax(p._native, cs, list(u), hs, data, without_report, baseline._native,
                                       minimax_spec(targets, without_scales))
    require(with_report.minimax is not None and no_report.minimax is not None
            and with_report.minimax.availability == no_report.minimax.availability == "values_only",
            "report-only constraint comparison unavailable", case="B3-report-only", quantity="availability")
    a, b = with_report.minimax, no_report.minimax
    require(len(a.rows) == len(b.rows) and len(a.values) == len(b.values) and len(a.jacobian) == len(b.jacobian),
            "report-only rows changed the enforced minimax structure", case="B3-report-only",
            quantity="row dimensions")
    for i, (row_a, row_b) in enumerate(zip(a.rows, b.rows)):
        group_a = None if row_a.kind == "target" else row_a.group
        group_b = None if row_b.kind == "target" else row_b.group
        require((row_a.kind, row_a.row, group_a, row_a.metric, row_a.sign, row_a.threshold, row_a.scale,
                 row_a.auxiliary, row_a.members) ==
                (row_b.kind, row_b.row, group_b, row_b.metric, row_b.sign, row_b.threshold, row_b.scale,
                 row_b.auxiliary, row_b.members), "report-only check changed an enforced row",
                case="B3-report-only", quantity=f"constraint-{i}/identity")
    for i, (x, y) in enumerate(zip(a.values, b.values)):
        compare("B3-report-only", f"constraint-{i}/value", x, y, "scaled constraint")
    for i, (x, y) in enumerate(zip(a.jacobian, b.jacobian)):
        compare("B3-report-only", f"constraint-{i}/column", x, y, "scaled derivative")
    for i, (x, y) in enumerate(zip(a.point, b.point)):
        compare("B3-report-only", f"point-{i}/value", x, y, "optimizer/auxiliary")
    with_feasible = all(value <= TOLERANCE for value in a.values)
    without_feasible = all(value <= TOLERANCE for value in b.values)
    require(with_feasible == without_feasible, "report-only result changed the point feasibility classification",
            case="B3-report-only", quantity="feasibility")
    return with_report


def write_outputs(summary):
    OUT.mkdir(exist_ok=True)
    columns = ("case", "quantity", "exact", "reference", "error", "tolerance", "atol", "molality", "unit", "passed",
               "kind", "row", "group", "metric", "sign", "threshold", "scale", "auxiliary", "members", "status")
    with (OUT / "minimax.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    (OUT / "minimax.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


def main():
    started = time.perf_counter()
    summary = {"status": "running",
               "conditions": {"temperature_K": 298.15, "pressure_Pa": 100000.0,
                              "branch": "liquid", "source_packet": "figiel-2025-reference-electrolytes/4"},
               "declarations": {"targets": "LiCl/NaCl/KCl gamma and Phi at 0.1, 1, 4.5 mol/kg",
                                "guard": "NaCl aqueous relative permittivity at 0.1, 1 mol/kg; source max and mean 0.1",
                                "q_decade": Q, "solver_tolerance": TOLERANCE,
                                "activity_distance": ACTIVITY_DISTANCE, "initial_margin": MARGIN,
                                "hold_energy_tolerance_J_mol": 1e-6, "hold_diameter_tolerance_angstrom": 1e-10,
                                "per_solve_iterations": 300, "per_solve_seconds": 120,
                                "shared_iterations": 600, "shared_seconds": 240,
                                "root_refinement_energy_J_mol": 1e-7,
                                "root_refinement_diameter_angstrom": 1e-11},
               "thread_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "EPCSAFT_BUILD_JOBS", "EPCSAFT_HSL_LIBRARY")},
               "commands": ["env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 agent-heavy .venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/minimax.py"],
               "claim_limits": ["synthetic rows verify solver integration, not physical validity",
                                "local stationarity does not establish global optimality or infeasibility",
                                "no #209 parameter adoption or physical validation"]}
    try:
        p, cs, hs, source, data, baseline, indices, transfer_prediction, transfer_detail = make_case()
        targets, epsilon, _, transfer_index = indices
        summary["inputs"] = {"packet_sha256": hashlib.sha256(PACKET.read_bytes()).hexdigest(),
                             "water_data_sha256": hashlib.sha256((ROOT / "data/packets/hamer-wu-1972-aqueous-alkali-halides/1/observations/hamer-wu-1972-aqueous-alkali-halides.csv").read_bytes()).hexdigest(),
                             "component_ids": list(p.component_ids), "coordinate_count": len(cs),
                             "target_rows": targets, "epsilon_rows": epsilon, "transfer_row": transfer_index,
                             "aqueous_solvent_fraction_of_methanol": 0.0,
                             "transfer": {**transfer_detail, "positive_source_value_J_mol": transfer_prediction,
                                          "report_only_datum_ratio_to_source": 0.1}}
        summary["inputs"]["observations"] = [{"row": i, "kind": item.kind.name, "observed": item.observed,
            "unit": "J/mol" if item.kind in (r.ObservationKind.SolvationGibbsEnergy, r.ObservationKind.TransferGibbsEnergy) else "1",
            "composition": list(item.composition), "reference_composition": list(item.reference_composition),
            "molality_mol_kg": item.molality, "ion_indices": list(item.ion_indices), "stoichiometry": list(item.stoichiometry)}
            for i, item in enumerate(data)]
        chart = r.evaluate_coordinates(p, cs, source, hs)
        summary["inputs"]["physical_parameter_order"] = [{"family": item.family.name, "indices": list(item.indices),
            "unit": "Angstrom" if item.family.name == "BornDiameter" else "1"} for item in chart.active]
        summary["inputs"]["optimizer_coordinates"] = [{"targets": [{"family": t.family.name, "indices": t.indices, "coefficient": t.coefficient} for t in c.targets],
            "origin": c.origin, "scale": c.scale, "lower": c.lower, "upper": c.upper} for c in cs]
        summary["inputs"]["hydration_targets_J_mol"] = [h.reference.observed for h in hs]
        b3_cs, scales, b3checks, point_results = run_b3(p, hs, data, baseline, indices)
        require(point_results[0].checks[-1].status == "fail",
                "source report-only energy LogRatio check did not fail",
                case="B3", quantity="report-only energy LogRatio",
                status=point_results[0].checks[-1].status)
        summary["B3"] = {"points": ["source", "F2-manufactured-truth",
                                    "B4-start-0-return", "B4-start-1-return"],
                         "affine_scales": list(scales), "check_groups": [item.name for item in b3checks],
                         "raw_rows_compared": len(data), "report_only_energy_log_ratio_failure_at_source": True}
        results, b4checks, elapsed = run_b4(p, cs, hs, source, data, baseline, indices)
        summary["B4"] = {"start_count": len(results), "elapsed_seconds": elapsed,
                         "check_groups": [item.name for item in b4checks],
                         "solves": [{"status": result.status,
                                     "termination": result.minimax.solver_status if result.minimax else "unavailable",
                                     "application_status": result.minimax.application_status if result.minimax else "unavailable",
                                     "availability": result.minimax.availability if result.minimax else "unavailable",
                                     "iterations": result.iterations, "failures": len(result.failures),
                                     "acceptance": result.acceptance,
                                     "active_lower": list(result.minimax.active_lower),
                                     "active_upper": list(result.minimax.active_upper),
                                     "hold_causes": list(result.minimax.hold_causes)} for result in results]}
        active_b4_rows = sorted(set(targets + epsilon))
        for index, result in enumerate(results):
            b3_fitted_point(f"B4-start-{index}-native", p, cs, hs, result, data, baseline, b4checks, active_b4_rows)
            b3_u = native_start(Parameters(result.parameters,p.component_ids,p.site_ids),b3_cs)
            b4_point = b3_point(f"B4-start-{index}-B3", p, b3_cs, b3_u, hs, data, baseline,
                                b3checks, targets, [1.0, 1000.0, 1000.0, 0.02])
            point_results.append(b4_point)
        # Removing report-only groups at the same point must leave every enforced
        # TNLP row, value and exact action unchanged.
        with_report = no_report_only_effect(p, cs, hs, source, data, baseline, targets,
                                            b4checks, [1.0, 1000.0, 0.02])
        require(with_report.checks[-1].status == "fail",
                "B4 report-only energy LogRatio did not fail at the source point",
                case="B3-report-only", quantity="energy LogRatio", status=with_report.checks[-1].status)
        guard_row_ids = {item.row for item in results[0].minimax.rows if item.kind in ("maximum", "link")}
        require(set(epsilon).issubset(guard_row_ids),
                "zero-Ceres-weight aqueous guard rows did not enter B",
                case="B3", quantity="zero-weight guard membership")
        summary["B3"]["report_only_constraint_system_unchanged"] = True
        summary["B3"]["zero_weight_enforced_rows"] = sorted(epsilon)
        summary["B3"]["onepoint_constraint_availability"] = with_report.minimax.availability

        targets, _, _, _ = indices
        refinement_shapes = []
        for index, result in enumerate(results):
            refinement_shapes.append(root_refinement(p, cs, hs, result, data, baseline, targets, b4checks, f"B4-start-{index}-refinement"))
        summary["root_refinement"] = {"constraint_shapes": refinement_shapes,
                                      "maximum_original_root_error": max((max(result.minimax.root_errors, default=0.0)
                                                                         for result in results if result.minimax), default=0.0)}
        wheels = list((ROOT / "build/environment-wheel").glob("*.whl"))
        require(len(wheels) == 1, "expected one prepared noneditable wheel",
                case="provenance", quantity="wheel count")
        summary["wheel_sha256"] = hashlib.sha256(wheels[0].read_bytes()).hexdigest()
        summary["status"] = "qualified"
    except Exception as error:
        summary["status"] = "failed"
        summary["diagnostic"] = {"type": type(error).__name__, "message": str(error)}
        summary["failures"] = failures
        summary["seconds"] = time.perf_counter() - started
        write_outputs(summary)
        raise
    summary["failures"] = failures
    summary["seconds"] = time.perf_counter() - started
    write_outputs(summary)
    print(json.dumps({"status": summary["status"], "seconds": summary["seconds"],
                      "csv": str(OUT / "minimax.csv"), "json": str(OUT / "minimax.json")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
