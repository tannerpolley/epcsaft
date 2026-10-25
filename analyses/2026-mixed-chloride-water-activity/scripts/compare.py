"""Apply the committed source realizations; obtain every water property natively."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
from pathlib import Path

from epcsaft import Mixture, Parameters, PropertyActionStatus, _core

HERE = Path(__file__).resolve().parents[1]
INPUTS, RESULTS = HERE / "inputs", HERE / "results"
OBS = INPUTS / "mixed-chloride-water-observations/1/observations"
FREEZE = "f594304279ed61e5d8c0e13f094815a9dfc48b1e"
PACKETS = {
    "Held": "held-2014-strategy2/1",
    "Figiel": "figiel-2025-chloride-activity-refit/2",
}
COMPONENTS = (
    "water",
    "sodium-cation",
    "magnesium-cation",
    "calcium-cation",
    "chloride-anion",
)
SOURCES = ("robinson-bower-1966", "rard-miller-1987", "dinane-mounir-2003")
MW_SOURCE = 0.0180153
T, P = 298.15, 100000.0


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(name, rows):
    path = RESULTS / name
    with (
        gzip.open(path, "wt", newline="")
        if name.endswith(".gz")
        else path.open("w", newline="")
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def half(token):
    return float(Decimal(5).scaleb(Decimal(token).as_tuple().exponent - 1))


def matches(row, selector):
    return all(
        (
            row["reference_salt"] == "NaCl"
            and float(row["reference_molality_mol_kg"]) <= 1.2590
        )
        if key == "low_NaCl_reference"
        else row[key] == value
        for key, value in selector.items()
    )


def interpolate(xs, ys, x):
    i = min(max(bisect.bisect_right(xs, x) - 1, 0), len(xs) - 2)
    return ys[i] + (x - xs[i]) * (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])


def realize(source, scenario):
    """Source arithmetic only; retain original row labels for selectors and buckets."""
    op, field, amount = (
        scenario["operation"],
        scenario["field"],
        float(scenario["amount"]),
    )
    selector = json.loads(scenario["selector"])
    filename = {
        SOURCES[0]: "robinson-bower-1966-table-1.csv",
        SOURCES[1]: "rard-miller-1987-tables-i-ii-mixtures.csv",
        SOURCES[2]: "dinane-mounir-2003-table-1.csv",
    }[source]
    raw = read(OBS / filename)
    anchors = (
        read(OBS / "robinson-bower-1966-table-2.csv") if source == SOURCES[0] else []
    )
    converted = (
        read(OBS / "robinson-bower-1966-source-converted.csv") if anchors else []
    )
    anchor_x = [float(r["mean_MB_mol_kg"]) for r in anchors]
    anchor_phi = [float(r["phi_NaCl_reference_printed"]) for r in anchors]
    if anchors:
        groups = {(r["set"], r["subseries"]): r for r in raw}
        for i, anchor in enumerate(anchors):
            subseries = [
                r for (group, _), r in groups.items() if group == anchor["set"]
            ]
            shifts = [
                amount
                if op == "round_shared_reference" and matches(r, selector)
                else amount * half(r["MB_reference_mol_kg"])
                if op == "joint_primary_rounding"
                else 0.0
                for r in subseries
            ]
            anchor_x[i] += sum(shifts) / len(shifts)
            if op == "round_reference_anchor" and matches(anchor, selector):
                anchor_phi[i] += amount
            if op == "joint_primary_rounding":
                anchor_phi[i] += amount * half(anchor["phi_NaCl_reference_printed"])
    rows = []
    for at, original in enumerate(raw):
        r = dict(original)
        selected = matches(original, selector)
        if (
            op in ("round_input", "round_shared_reference", "round_shared_fraction")
            and selected
        ):
            r[field] = str(float(r[field]) + amount)
        if op == "joint_primary_rounding":
            fields = {
                SOURCES[0]: ("m_NaCl_mol_kg", "m_CaCl2_mol_kg", "MB_reference_mol_kg"),
                SOURCES[1]: (
                    "total_salt_molality_mT_mol_kg",
                    "y_NaCl_printed",
                    "reference_molality_mol_kg",
                    "reference_osmotic_coefficient",
                ),
                SOURCES[2]: ("m_nacl_mol_kg_h2o", "m_mgcl2_mol_kg_h2o", "observed_aw"),
            }[source]
            for key in fields:
                if float(r[key]) != 0.0:
                    r[key] = str(float(r[key]) + amount * half(original[key]))
        if source == SOURCES[0]:
            mna, md, mr = (
                float(r[k])
                for k in ("m_NaCl_mol_kg", "m_CaCl2_mol_kg", "MB_reference_mol_kg")
            )
            q, ionic = 2 * mna + 3 * md, mna + 3 * md
            pref = interpolate(anchor_x, anchor_phi, mr)
            near = sorted(range(len(anchor_x)), key=lambda k: abs(anchor_x[k] - mr))[:3]
            quadratic = sum(
                anchor_phi[j]
                * math.prod(
                    (mr - anchor_x[k]) / (anchor_x[j] - anchor_x[k])
                    for k in near
                    if k != j
                )
                for j in near
            )
            reconstruction = abs(2 * mr * (quadratic - pref) / q)
            observed = 2 * mr * pref / q
            if op == "observation_alternative":
                old = converted[at]
                observed = float(old["source_converted_phi"])
                if field in ("printed_R", "quadratic_reference_and_printed_R"):
                    observed = float(old["phi_using_printed_R_sensitivity"])
                if field in (
                    "quadratic_reference",
                    "quadratic_reference_and_printed_R",
                ):
                    observed *= float(
                        old["reference_phi_quadratic_sensitivity"]
                    ) / float(old["reference_phi_reconstructed"])
            identity = "/".join(
                original[k] for k in ("set", "subseries", "row_in_subseries")
            )
            fraction = original["y_Ca_printed"]
            role = (
                "CaCl2 endpoint"
                if original["inferred_zero_endpoint"] == "true"
                else "mixture"
            )
            quantity, salt, reference = "Phi", "CaCl2", "NaCl"
            ratio, printed_ratio, flagged = (
                2 * mr / q,
                float(original["R_observed"]),
                False,
            )
        elif source == SOURCES[1]:
            y, mt = (
                float(r["y_NaCl_printed"]),
                float(r["total_salt_molality_mT_mol_kg"]),
            )
            mna, md = 3 * y * mt / (1 + 2 * y), (1 - y) * mt / (1 + 2 * y)
            mr, pref, reference = (
                float(r["reference_molality_mol_kg"]),
                float(r["reference_osmotic_coefficient"]),
                r["reference_salt"],
            )
            if selected:
                if op == "common_concentration_scale":
                    mna, md, mr = (
                        mna * (1 + amount),
                        md * (1 + amount),
                        mr * (1 + amount),
                    )
                elif op == "Mg_stock_scale":
                    md *= 1 + amount
                elif op == "Na_stock_scale":
                    mna *= 1 + amount
                    if reference == "NaCl":
                        mr *= 1 + amount
                elif op == "Ca_reference_Phi_scale" and reference == "CaCl2":
                    pref *= 1 + amount
            q, ionic = 2 * mna + 3 * md, mna + 3 * md
            ratio = (2 if reference == "NaCl" else 3) * mr / q
            observed = ratio * pref
            if (
                op == "reference_change_Phi_offset"
                and selected
                and reference == "CaCl2"
            ):
                observed += amount
            identity = "/".join(original[k] for k in ("table", "row", "y_NaCl_printed"))
            fraction = original["y_NaCl_printed"]
            flagged = original["source_excluded_fit_point"] == "true"
            role, quantity, salt, printed_ratio = "mixture", "Phi", "MgCl2", None
            reconstruction = None
        else:
            mna, md = float(r["m_nacl_mol_kg_h2o"]), float(r["m_mgcl2_mol_kg_h2o"])
            q, ionic = 2 * mna + 3 * md, mna + 3 * md
            observed = float(r["observed_aw"])
            nominal_na, nominal_md = (
                float(original["m_nacl_mol_kg_h2o"]),
                float(original["m_mgcl2_mol_kg_h2o"]),
            )
            fraction = str(round(nominal_na / (nominal_na + 3 * nominal_md), 12))
            identity, role, quantity, salt = original["row"], "mixture", "aw", "MgCl2"
            reference, mr, pref, ratio, printed_ratio, flagged = (
                "",
                None,
                None,
                None,
                None,
                False,
            )
            reconstruction = None
        rows.append(
            {
                "source": source,
                "row_id": identity,
                "role": role,
                "salt": salt,
                "quantity": quantity,
                "fraction_label": fraction,
                "m_NaCl": mna,
                "m_divalent": md,
                "I_mol_kg": ionic,
                "q_mol_kg": q,
                "observed": observed,
                "reference_salt": reference,
                "reference_molality": mr,
                "reference_Phi": pref,
                "source_ratio": ratio,
                "printed_ratio": printed_ratio,
                "source_flagged": flagged,
                "reference_reconstruction_sensitivity": reconstruction,
                "included": not (op == "omit_source_flagged" and flagged),
            }
        )
    return rows


def declarations():
    return {
        name: Mixture(
            Parameters.from_json(
                INPUTS / packet / "parameters/parameters.json", components=COMPONENTS
            )
        )
        for name, packet in PACKETS.items()
    }


def initialize():
    global MODELS, REFERENCES, MASSES
    MODELS = declarations()
    REFERENCES = {
        name: model.reference(T, P, solvent_fractions=(1, 0, 0, 0, 0))
        for name, model in MODELS.items()
    }
    MASSES = {
        name: Parameters.from_json(
            INPUTS / packet / "parameters/parameters.json", components=COMPONENTS
        ).molar_masses[0]
        for name, packet in PACKETS.items()
    }


def predict(task):
    key, packet, salt, mna, md, repeat = task
    amounts = (
        1 / MASSES[packet],
        mna,
        md if salt == "MgCl2" else 0.0,
        md if salt == "CaCl2" else 0.0,
        mna + 2 * md,
    )
    x = tuple(n / sum(amounts) for n in amounts)
    result = REFERENCES[packet].solvent_activity(x=x)
    state = result.state
    scalar_phi = (
        result.osmotic_coefficient
        if isinstance(result.osmotic_coefficient, float)
        else None
    )
    valid = (
        result.status == PropertyActionStatus.Available
        and result.solvent_activity is not None
        and scalar_phi is not None
        and math.isfinite(scalar_phi)
        and math.isfinite(result.solvent_activity)
    )
    convergence = state.density_convergence
    row = {
        "prediction_id": key,
        "declaration": packet,
        "salt": salt,
        "m_NaCl": mna,
        "m_divalent": md,
        "water_molar_mass_kg_mol": MASSES[packet],
        "aw": result.solvent_activity,
        "Phi": scalar_phi,
        "available": valid,
        "status": str(result.status),
        "cause": result.cause,
        "osmotic_cause": result.osmotic_cause,
        "density_mol_m3": state.molar_density,
        "pressure_Pa": state.pressure,
        "pressure_error_Pa": state.pressure - P,
        "pressure_density_derivative": state.pressure_density_derivative,
        "branch": str(state.density_branch),
        "stable_roots": state.stable_root_count,
        "relative_permittivity": state.bulk_relative_permittivity,
        "charge_residual": x[1] + 2 * x[2] + 2 * x[3] - x[4],
        "density_residual": convergence.residual if convergence else None,
        "density_stop": str(convergence.stop) if convergence else "",
        "density_iterations": convergence.iterations if convergence else None,
        "density_evaluations": convergence.evaluations if convergence else None,
        "repeat_aw": None,
        "repeat_Phi": None,
    }
    if repeat:
        second = REFERENCES[packet].solvent_activity(x=x)
        row.update(
            repeat_aw=second.solvent_activity,
            repeat_Phi=second.osmotic_coefficient
            if isinstance(second.osmotic_coefficient, float)
            else None,
        )
    return row


def prediction_key(packet, row):
    values = (packet, row["salt"], row["m_NaCl"], row["m_divalent"])
    return hashlib.sha256(json.dumps(values).encode()).hexdigest()[:20]


def metrics(rows, quantity):
    valid = [r for r in rows if r["available"]]
    values = [r["residual"] for r in valid]
    return {
        "total_count": len(rows),
        "valid_count": len(valid),
        "native_failure_count": len(rows) - len(valid),
        "bias": sum(values) / len(values) if values else None,
        "MAE": sum(abs(v) for v in values) / len(values) if values else None,
        "RMS": math.sqrt(sum(v * v for v in values) / len(values)) if values else None,
        "maximum_absolute_error": max(map(abs, values)) if values else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--declare-only", action="store_true")
    parser.add_argument("--workers", type=int, default=os.cpu_count())
    args = parser.parse_args()
    if args.declare_only:
        declarations()
        return
    RESULTS.mkdir(exist_ok=True)
    manifest = json.loads((INPUTS / "input-files.json").read_text())
    for item in manifest["files"]:
        path = INPUTS / item["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], path
        committed = subprocess.check_output(
            ["git", "show", f"{FREEZE}:{path.relative_to(HERE.parents[1])}"]
        )
        assert hashlib.sha256(committed).hexdigest() == item["sha256"], path
    specification = json.loads((INPUTS / "comparison-specification.json").read_text())
    scenarios = read(INPUTS / "frozen-realizations.csv")
    assert len(scenarios) == 1467 and scenarios[0]["operation"] == "nominal"
    nominal = {source: realize(source, scenarios[0]) for source in SOURCES}
    nominal_by_id = {
        (r["source"], r["row_id"]): r for rows in nominal.values() for r in rows
    }
    for source, filename in [
        (SOURCES[0], "robinson-bower-1966-source-converted.csv"),
        (SOURCES[1], "rard-miller-1987-source-converted.csv"),
    ]:
        for row, converted in zip(nominal[source], read(OBS / filename), strict=True):
            assert (
                abs(row["observed"] - float(converted["source_converted_phi"])) < 1e-14
            )
    realized = {}
    for scenario in scenarios:
        selected_sources = (
            SOURCES if scenario["operation"] == "nominal" else (scenario["source"],)
        )
        for source in selected_sources:
            realized[scenario["realization_id"], source] = realize(source, scenario)
    # Separate binary reference checks from every mixture metric.
    references = []
    for raw in read(OBS / "robinson-bower-1966-table-2.csv"):
        references.append(
            {
                "source": SOURCES[0],
                "row_id": "anchor/" + raw["set"],
                "role": "NaCl anchor",
                "salt": "CaCl2",
                "m_NaCl": float(raw["mean_MB_mol_kg"]),
                "m_divalent": 0.0,
                "observed": float(raw["phi_NaCl_reference_printed"]),
                "quantity": "Phi",
            }
        )
    for raw in read(OBS / "rard-miller-1987-table-iii-mgcl2.csv"):
        references.append(
            {
                "source": SOURCES[1],
                "row_id": "binary/" + raw["row"],
                "role": "MgCl2 binary",
                "salt": "MgCl2",
                "m_NaCl": 0.0,
                "m_divalent": float(raw["molality_MgCl2_mol_kg"]),
                "observed": float(raw["phi_MgCl2"]),
                "quantity": "Phi",
            }
        )
    tasks = {}
    for (scenario_id, _), rows in realized.items():
        for row in rows:
            for packet in PACKETS:
                key = prediction_key(packet, row)
                if key not in tasks or scenario_id == "r0000":
                    tasks[key] = (
                        key,
                        packet,
                        row["salt"],
                        row["m_NaCl"],
                        row["m_divalent"],
                        scenario_id == "r0000",
                    )
    for row in references:
        for packet in PACKETS:
            key = prediction_key(packet, row)
            tasks[key] = (
                key,
                packet,
                row["salt"],
                row["m_NaCl"],
                row["m_divalent"],
                True,
            )
    print(
        f"{len(scenarios)} frozen realizations; {len(tasks)} unique native states; {args.workers} processes",
        flush=True,
    )
    initialize()
    pure_rows = []
    for packet, reference in REFERENCES.items():
        pure = reference.solvent_activity(x=(1, 0, 0, 0, 0))
        assert (
            pure.solvent_activity == 1.0
            and pure.osmotic_cause == "OsmoticUndefinedPureSolvent"
        )
        pure_rows.append(
            {
                "declaration": packet,
                "aw": pure.solvent_activity,
                "Phi": str(pure.osmotic_coefficient),
                "osmotic_cause": pure.osmotic_cause,
                "density_mol_m3": pure.state.molar_density,
                "pressure_Pa": pure.state.pressure,
                "pressure_error_Pa": pure.state.pressure - P,
                "pressure_density_derivative": pure.state.pressure_density_derivative,
                "relative_permittivity": pure.state.bulk_relative_permittivity,
                "branch": str(pure.state.density_branch),
            }
        )
    write("pure-reference-checks.csv", pure_rows)
    predictions = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize) as pool:
        for at, row in enumerate(pool.map(predict, tasks.values(), chunksize=16), 1):
            predictions.append(row)
            if at % 1000 == 0:
                print(f"native states {at}/{len(tasks)}", flush=True)
    write("native-predictions.csv", predictions)
    by_prediction = {r["prediction_id"]: r for r in predictions}
    repeats = [r for r in predictions if r["repeat_aw"] is not None]
    assert all(
        r["available"]
        and abs(r["aw"] - r["repeat_aw"]) <= 1e-10
        and abs(r["Phi"] - r["repeat_Phi"]) <= 1e-8
        for r in repeats
    )
    residual_rows, summaries, fraction_summaries = [], [], []
    for (scenario_id, source), rows in realized.items():
        for packet in PACKETS:
            mixed = []
            for row in rows:
                prediction = by_prediction[prediction_key(packet, row)]
                old = nominal_by_id[source, row["row_id"]]
                primary = prediction[row["quantity"]]
                phi_observed = (
                    row["observed"]
                    if row["quantity"] == "Phi"
                    else -math.log(row["observed"]) / (MW_SOURCE * row["q_mol_kg"])
                )
                phi_residual = (
                    prediction["Phi"] - phi_observed
                    if prediction["available"]
                    else None
                )
                result = dict(
                    realization_id=scenario_id,
                    declaration=packet,
                    prediction_id=prediction["prediction_id"],
                    **row,
                    nominal_I_mol_kg=old["I_mol_kg"],
                    I_bucket="I<=3" if old["I_mol_kg"] <= 3 else "I>3",
                    crossed_I_boundary=(old["I_mol_kg"] <= 3) != (row["I_mol_kg"] <= 3),
                    model=primary,
                    available=prediction["available"],
                    residual=primary - row["observed"]
                    if prediction["available"]
                    else None,
                    diagnostic_Phi_observed=phi_observed,
                    diagnostic_Phi_residual=phi_residual,
                    residual_over_aw_dispersion=(primary - row["observed"]) / 0.0006
                    if row["quantity"] == "aw" and prediction["available"]
                    else None,
                    residual_over_Phi_dispersion=phi_residual / 0.028
                    if row["quantity"] == "aw" and phi_residual is not None
                    else None,
                    relative_Phi_residual=phi_residual / phi_observed
                    if phi_residual is not None
                    else None,
                    relative_residual_over_scatter_0p1pct=phi_residual
                    / phi_observed
                    / 0.001
                    if source == SOURCES[1] and phi_residual is not None
                    else None,
                    relative_residual_over_scatter_0p25pct=phi_residual
                    / phi_observed
                    / 0.0025
                    if source == SOURCES[1] and phi_residual is not None
                    else None,
                    relative_residual_over_concentration_0p13pct=phi_residual
                    / phi_observed
                    / 0.0013
                    if source == SOURCES[1] and phi_residual is not None
                    else None,
                    relative_residual_over_Ca_reference_0p2pct=phi_residual
                    / phi_observed
                    / 0.002
                    if source == SOURCES[1]
                    and row["reference_salt"] == "CaCl2"
                    and phi_residual is not None
                    else None,
                    relative_residual_over_Ca_reference_0p4pct=phi_residual
                    / phi_observed
                    / 0.004
                    if source == SOURCES[1]
                    and row["reference_salt"] == "CaCl2"
                    and phi_residual is not None
                    else None,
                    water_mass_conversion_delta_aw=math.exp(
                        -prediction["water_molar_mass_kg_mol"]
                        * row["q_mol_kg"]
                        * phi_observed
                    )
                    - math.exp(-MW_SOURCE * row["q_mol_kg"] * phi_observed),
                    R_residual=(
                        prediction["Phi"] / row["reference_Phi"] - row["source_ratio"]
                    )
                    if row["reference_Phi"] and prediction["available"]
                    else None,
                    printed_R_residual=(
                        prediction["Phi"] / row["reference_Phi"] - row["printed_ratio"]
                    )
                    if row["printed_ratio"] and prediction["available"]
                    else None,
                )
                residual_rows.append(result)
                if row["role"] == "mixture" and (
                    row["included"] or not prediction["available"]
                ):
                    mixed.append(result)
            quantity = rows[0]["quantity"]
            for bucket in ("I<=3", "I>3", "all"):
                selected = [
                    r for r in mixed if bucket == "all" or r["I_bucket"] == bucket
                ]
                result = dict(
                    realization_id=scenario_id,
                    declaration=packet,
                    source=source,
                    system="NaCl-" + rows[0]["salt"] + "-water",
                    quantity=quantity,
                    I_bucket=bucket,
                    **metrics(selected, quantity),
                )
                limits = specification["primary_limits"][quantity]
                result["passed"] = (
                    result["native_failure_count"] == 0
                    and result["RMS"] is not None
                    and result["RMS"] <= limits["RMS"]
                    and result["maximum_absolute_error"]
                    <= limits["maximum_absolute_error"]
                )
                summaries.append(result)
            for fraction in sorted({r["fraction_label"] for r in mixed}):
                fraction_summaries.append(
                    dict(
                        realization_id=scenario_id,
                        declaration=packet,
                        source=source,
                        quantity=quantity,
                        fraction_label=fraction,
                        **metrics(
                            [r for r in mixed if r["fraction_label"] == fraction],
                            quantity,
                        ),
                    )
                )
    # Keep repeated source rows adjacent so the complete retained record compresses.
    residual_rows.sort(
        key=lambda r: (r["source"], r["row_id"], r["declaration"], r["realization_id"])
    )
    write("residuals.csv.gz", residual_rows)
    write("bucket-metrics.csv", summaries)
    write("fraction-metrics.csv", fraction_summaries)
    write(
        "nominal-residuals.csv",
        [r for r in residual_rows if r["realization_id"] == "r0000"],
    )
    binary_rows = []
    for row in references:
        for packet in PACKETS:
            prediction = by_prediction[prediction_key(packet, row)]
            binary_rows.append(
                dict(
                    declaration=packet,
                    prediction_id=prediction["prediction_id"],
                    **row,
                    model=prediction["Phi"],
                    available=prediction["available"],
                    residual=prediction["Phi"] - row["observed"]
                    if prediction["available"]
                    else None,
                )
            )
    write("binary-reference-checks.csv", binary_rows)
    indexed = {
        (r["realization_id"], r["source"], r["I_bucket"], r["declaration"]): r
        for r in summaries
    }
    transfers, coverage = [], []
    for scenario in scenarios:
        for bucket in specification["primary_buckets"]:
            source = bucket["source"]
            effective = (
                scenario["realization_id"]
                if scenario["source"] in ("all", source)
                else "r0000"
            )
            held, figiel = [
                indexed[effective, source, bucket["I_bucket"], packet]
                for packet in PACKETS
            ]
            resolution = specification["transfer_test"][
                bucket["quantity"] + "_resolution"
            ]
            improvement = (
                held["RMS"] - figiel["RMS"]
                if held["RMS"] is not None and figiel["RMS"] is not None
                else None
            )
            growth = (
                figiel["maximum_absolute_error"] - held["maximum_absolute_error"]
                if improvement is not None
                else None
            )
            passed = (
                not held["native_failure_count"]
                and not figiel["native_failure_count"]
                and improvement is not None
                and improvement > resolution
                and growth <= resolution
            )
            transfers.append(
                dict(
                    realization_id=scenario["realization_id"],
                    effective_source_realization=effective,
                    **bucket,
                    RMS_Held=held["RMS"],
                    RMS_Figiel=figiel["RMS"],
                    RMS_improvement=improvement,
                    maximum_error_growth=growth,
                    resolution=resolution,
                    passed=passed,
                )
            )
            coverage.append(
                {
                    "realization_id": scenario["realization_id"],
                    "source": source,
                    "I_bucket": bucket["I_bucket"],
                    "effective_source_realization": effective,
                    "Held_passed": held["passed"],
                    "Figiel_passed": figiel["passed"],
                }
            )
    write("transfer-tests.csv", transfers)
    write("realization-coverage.csv", coverage)
    reconstruction_max = max(
        r["reference_reconstruction_sensitivity"]
        for r in residual_rows
        if r["reference_reconstruction_sensitivity"] is not None
    )
    reconstruction_passed = reconstruction_max <= 0.001
    verdicts = []
    for bucket in specification["primary_buckets"]:
        for packet in PACKETS:
            selected = [
                r
                for r in coverage
                if r["source"] == bucket["source"]
                and r["I_bucket"] == bucket["I_bucket"]
            ]
            central = indexed["r0000", bucket["source"], bucket["I_bucket"], packet]
            verdicts.append(
                dict(
                    declaration=packet,
                    **bucket,
                    **{
                        k: central[k]
                        for k in (
                            "total_count",
                            "bias",
                            "MAE",
                            "RMS",
                            "maximum_absolute_error",
                            "native_failure_count",
                        )
                    },
                    nominal_passed=central["passed"],
                    frozen_realizations=len(scenarios),
                    passing_realizations=sum(r[packet + "_passed"] for r in selected),
                    robust_passed=all(r[packet + "_passed"] for r in selected)
                    and (reconstruction_passed or bucket["source"] != SOURCES[0]),
                )
            )
    write("verdicts.csv", verdicts)
    # Formal inventory conversion on the frozen one-liter water-mass basis.
    application = json.loads((INPUTS / "application-basis.json").read_text())
    feed = {
        r["analyte"]: float(r["value"])
        for r in read(INPUTS / "smackover-selected-base-feed-ms2.csv")
        if r["value"] != "not_reported"
    }
    water_mass = feed["specific_gravity"] - feed["TDS"] * 1e-6
    application_rows = [
        {
            "ion": ion,
            "concentration_mg_L": feed[ion],
            "molar_mass_g_mol": mass,
            "charge": application["charges"][ion],
            "molality_mol_kg": feed[ion] / 1000 / mass / water_mass,
        }
        for ion, mass in application["molar_masses_g_mol"].items()
    ]
    write("ms2-formal-molalities.csv", application_rows)
    m = {r["ion"]: r["molality_mol_kg"] for r in application_rows}
    failures = [r for r in predictions if not r["available"]]
    summary = {
        "specification_commit": FREEZE,
        "engine_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "specification_sha256": hashlib.sha256(
            (INPUTS / "comparison-specification.json").read_bytes()
        ).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "native_binaries": {
            str(path.relative_to(HERE.parents[1])): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in (
                Path(_core.__file__).resolve(),
                HERE.parents[1] / "build/engine-meson/shared/libepcsaft_native.so",
            )
        },
        "conditions": {"T_K": T, "P_Pa": P, "component_order": COMPONENTS},
        "realizations": len(scenarios),
        "unique_native_states": len(predictions),
        "native_failure_count": len(failures),
        "repeated_states": len(repeats),
        "repeat_aw_max": max(abs(r["aw"] - r["repeat_aw"]) for r in repeats),
        "repeat_Phi_max": max(abs(r["Phi"] - r["repeat_Phi"]) for r in repeats),
        "maximum_pressure_error_Pa": max(
            abs(r["pressure_error_Pa"]) for r in predictions if r["available"]
        ),
        "minimum_pressure_density_derivative": min(
            r["pressure_density_derivative"] for r in predictions if r["available"]
        ),
        "maximum_charge_residual": max(abs(r["charge_residual"]) for r in predictions),
        "Robinson_maximum_reconstruction_sensitivity": reconstruction_max,
        "Robinson_reconstruction_passed": reconstruction_passed,
        "transfer_passed": all(r["passed"] for r in transfers),
        "transfer_rows": len(transfers),
        "failed_transfer_rows": sum(not r["passed"] for r in transfers),
        "declarations_qualified": {
            packet: all(
                r["robust_passed"] for r in verdicts if r["declaration"] == packet
            )
            for packet in PACKETS
        },
        "ms2": {
            "water_kg_per_L": water_mass,
            "major_ion_I_mol_kg": sum(
                r["molality_mol_kg"] * r["charge"] ** 2 / 2 for r in application_rows
            ),
            "charge_imbalance_mol_kg": sum(
                r["molality_mol_kg"] * r["charge"] for r in application_rows
            ),
            "chloride_approximation_I_mol_kg": m["Na"] + 3 * m["Mg"] + 3 * m["Ca"],
            "extrapolation": "Composition to simultaneous Na/Ca/Mg feed; I lies inside the measured ranges.",
        },
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    if failures:
        print(
            "Native failures retained; diagnose before interpreting or changing approach.",
            flush=True,
        )


if __name__ == "__main__":
    main()
