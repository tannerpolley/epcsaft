"""Freeze main-derived states and independent fixed-density verification roots."""

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

from declarations import HERE, ROOT, cross_check, key_solves, models
from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R
from epcsaft import Parameters, _core, provenance
from epcsaft import equilibrium as eq

INTERVALS = 16384
LOG_TOLERANCE = 1e-14


def write_csv(path, rows):
    assert rows, f"empty result: {path}"
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=list(dict.fromkeys(k for row in rows for k in row)),
            lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows)


def diagnostic_point(point):
    fields = (
        "molar_density",
        "packing_fraction",
        "pressure",
        "pressure_density_derivative",
        "lower",
        "upper",
        "residual",
        "mechanical_sign",
        "linked",
        "declared",
        "cause",
        "stop",
    )
    return {name: getattr(point, name) for name in fields}


def compare_roots(previous, densities):
    assert len(previous) == len(densities), "frozen oracle root count"
    for before, density in zip(previous, densities, strict=True):
        delta = abs(density / float(before["rho_mol_m3"]) - 1)
        assert delta <= 1e-10, "frozen oracle density"
        yield before, density, delta


def enumerate_roots(model, case, diagnostic_output=None, *, intervals=INTERVALS, native_result=None):
    """Nonadaptive native diagnostic crossings/refine_crossing; never solve_density."""
    T, P, x = case["T_K"], case["P_Pa"], case["x"]
    upper = model.molar_density_upper_bound(T, x=x)
    lo = P / (R * T) * 1e-8
    hi = upper * (1 - 1e-12)
    loglo, loghi = math.log(lo), math.log(hi)
    nodes = [
        math.exp(loglo + (loghi - loglo) * (i / intervals)) for i in range(intervals + 1)
    ]
    result = model.density_diagnostics(
        T,
        P,
        x=x,
        nodes=nodes,
        density_tolerance=upper * LOG_TOLERANCE,
        log_tolerance=LOG_TOLERANCE,
        max_evaluations=intervals * 3,
        max_seconds=120,
        adaptive=False,
        finite_upper=hi,
        evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    if native_result is not None:
        native_result.append((result, lo, hi, upper * LOG_TOLERANCE))
        return
    if diagnostic_output is not None:
        diagnostic_output.write_text(json.dumps(dict(
            case=case, stop=result.stop, evaluations=result.evaluations,
            roots=[diagnostic_point(point) for point in result.pressure_roots],
            spinodals=[diagnostic_point(point) for point in result.spinodals],
            unresolved=[diagnostic_point(point) for point in result.unresolved],
        ), indent=2) + "\n")
    assert result.stop == "resolved_crossings" and not result.unresolved, (
        case["id"],
        result.stop,
    )
    assert not any(sample.cause for sample in result.samples), case["id"]
    roots = []
    for i, root in enumerate(result.pressure_roots):
        roots.append(
            dict(
                id=case["id"],
                state_index=case.get("index", 0),
                root=i,
                rho_mol_m3=root.molar_density,
                slope_Pa_m3_mol=root.pressure_density_derivative,
                slope_sign=root.mechanical_sign,
                pressure_Pa=root.pressure,
                pressure_relative_error=abs(root.pressure - P) / P,
                branch="unstable" if root.mechanical_sign < 0 else "",
                scan_intervals=intervals,
                lower_mol_m3=lo,
                upper_mol_m3=hi,
                close_packing_mol_m3=upper,
                log_tolerance=LOG_TOLERANCE,
                bracket_log_width=math.log(root.upper / root.lower),
                stop=root.stop,
            )
        )
    stable = [root for root in roots if root["slope_sign"] > 0]
    negative = [
        sample.molar_density for sample in result.samples if sample.mechanical_sign < 0
    ]
    for i, root in enumerate(stable):
        root["branch"] = (
            ("vapor" if i == 0 else "liquid")
            if len(stable) > 1
            else ("vapor" if root["rho_mol_m3"] < negative[0] else "liquid")
            if negative
            else "single"
        )
    assert roots, case["id"]
    return roots


def water_sigma(T):
    packet = json.loads(
        (
            ROOT / "data/packets/held-2014-strategy2/1/parameters/parameters.json"
        ).read_text()
    )
    correlation = next(
        c
        for c in packet["correlations"]
        if c.get("component_id") == "water" and c["family"] == "segment_diameter"
    )
    values = [
        T,
        correlation["constant"]["value"]["magnitude"],
        len(correlation["terms"]),
    ]
    for term in correlation["terms"]:
        values.extend(
            [
                term["amplitude"]["value"]["magnitude"],
                term["exponent_coefficient"]["value"]["magnitude"],
            ]
        )
    executable = ROOT / "build/density-benchmark/water-sigma"
    if not executable.exists():
        available = int(
            next(
                line.split()[1]
                for line in Path("/proc/meminfo").read_text().splitlines()
                if line.startswith("MemAvailable:")
            )
        )
        assert available >= 12 * 1024 * 1024, (
            "water-sigma helper compile requires MemAvailable >=12 GiB"
        )
        subprocess.run(
            [
                "nice",
                "-n",
                "10",
                "c++",
                "-std=c++17",
                "-O2",
                "-I" + str(ROOT / "engine/native"),
                "-I" + str(ROOT / "engine/build/subprojects/CppAD-20250000.3/include"),
                str(HERE / "scripts/water_sigma.cpp"),
                "-o",
                str(executable),
            ],
            check=True,
        )
    return float(
        subprocess.check_output(
            [str(executable)], input=" ".join(map(str, values)) + "\n", text=True
        )
    )


def refresh(output=None, release_metadata=None):
    assert os.getpriority(os.PRIO_PROCESS, 0) == 10, (
        "oracle generation must use nice 10"
    )
    document = json.loads((HERE / "inputs/cases.json").read_text())
    wheel, wheel_sha256 = provenance()
    wheel = Path(wheel)
    assert hashlib.sha256(wheel.read_bytes()).hexdigest() == wheel_sha256
    release = json.loads(release_metadata.read_text()) if release_metadata else None
    if release is not None:
        assert release["sha256"] == wheel_sha256
    source_revision = release["commit"] if release else document["main_sha"]
    output = output or HERE
    output.mkdir(parents=True, exist_ok=True)

    def tree(revision):
        return subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", revision + ":engine"], text=True
        ).strip()

    assert tree("HEAD") == tree(source_revision), (
        "oracle Engine sources differ from pinned main"
    )
    document["oracle_engine_provenance"] = [str(wheel), wheel_sha256]
    document["wheel_sha256"] = wheel_sha256
    document["reference_engine_tree"] = tree(source_revision)
    ms = models()
    old = list(csv.DictReader((HERE / "results/oracle.csv").open()))
    legacy = [row for row in old if int(row.get("state_index", 0)) == 0]
    rows = []
    comparisons = []
    quantities = document["model_quantities"]
    gap = quantities["vapor_spinodal_P_Pa"] - quantities["liquid_spinodal_P_Pa"]
    for case in document["cases"]:
        if case["tier"] == "C":
            changing = case["operation"] in (
                "C1",
                "C2",
                "bubble",
                "flash",
                "actions-C2",
                "actions-C3",
            )
            baseT = {
                "C1": 300.0,
                "C2": 373.15,
                "bubble": 180.0,
                "flash": 180.0,
                "actions-C2": 373.15,
                "actions-C3": 298.15,
            }.get(case["operation"])
            case["states"] = (
                [dict(index=i, T_K=baseT + i * 0.001) for i in range(32)]
                if changing
                else [dict(index=0)]
            )
            case["mode"] = (
                "changing temperature"
                if changing
                else "repeated identical input (rewards caching)"
            )
            if changing:
                case["N"] = 32
            continue
        case["mode"] = "changing pressure"
        case["N"] = 32
        increment = case["P_Pa"] * 1e-4
        if "metastable-liquid" in case["id"]:
            increment = 0.0005 * gap
        if "metastable-vapor" in case["id"]:
            increment = -0.0005 * gap
        if "critical-liquid" in case["id"] or "critical-vapor" in case["id"]:
            increment = 0.0001 * gap
        case["sequence_source"] = (
            "32 distinct pressures from base to base + "
            + repr(increment)
            + " Pa; fixed T and x"
        )
        case["states"] = []
        expected = None
        for i in range(32):
            state = dict(
                index=i,
                T_K=case["T_K"],
                P_Pa=case["P_Pa"] + increment * i / 31,
                x=case["x"],
                branch=case["branch"],
            )
            if i % 8 == 0:
                print("native oracle", case["id"], i, flush=True)
            roots = enumerate_roots(ms[case["model"]], dict(state, id=case["id"]),
                                    output / f"{case['id']}-state-{i}-diagnostics.json" if output != HERE else None)
            if expected is None:
                expected = [r["branch"] for r in roots]
            assert [r["branch"] for r in roots] == expected, (
                "sequence left the root regime"
            )
            stable = [root for root in roots if root["slope_sign"] > 0]
            selected = stable[0] if case["branch"] == "vapor" else stable[-1]
            state.update(
                oracle_rho_mol_m3=selected["rho_mol_m3"],
                oracle_branch=selected["branch"],
                oracle_stable_root_count=len(stable),
            )
            if case["model"] == "water":
                state["sigma_angstrom"] = water_sigma(state["T_K"])
            case["states"].append(state)
            rows.extend(roots)
            if i == 0:
                previous = [r for r in legacy if r["id"] == case["id"]]
                for before, _, delta in compare_roots(previous, [r["rho_mol_m3"] for r in roots]):
                    comparisons.append(
                        dict(
                            id=case["id"], root=int(before["root"]), relative_difference=delta
                        )
                    )
        if "metastable" in case["id"]:
            for state in case["states"]:
                assert (
                    quantities["liquid_spinodal_P_Pa"]
                    < state["P_Pa"]
                    < quantities["P_sat_0_99_Tc_Pa"]
                    if "liquid" in case["id"]
                    else quantities["P_sat_0_99_Tc_Pa"]
                    < state["P_Pa"]
                    < quantities["vapor_spinodal_P_Pa"]
                )
    assert len(comparisons) == 56
    document["oracle_method"] = (
        "Engine density_diagnostics, adaptive=False: separate crossing detection and refine_crossing; shares at_density/EOS/bounds, never scan_roots/refine/anchor"
    )
    document.pop("pressure_floor_source", None)
    (output / "cases.json" if output != HERE else HERE / "inputs/cases.json").write_text(json.dumps(document, indent=2) + "\n")
    write_csv(output / "oracle.csv" if output != HERE else HERE / "results/oracle.csv", rows)
    comparison_path = output / "oracle-reproduction.csv" if output != HERE else HERE / "results/oracle-reproduction.csv"
    if not comparison_path.exists():
        write_csv(comparison_path, comparisons)


def saturation(model, T):
    result = eq.solve_equilibrium(model, model.saturation(T))
    assert result.success, result.message
    return result.pressure


def main():
    ms = models()
    propane = ms["C1"]
    critical = eq.CriticalDeclaration()
    critical.density_guess_mol_per_m3 = 5000
    problem = eq.Problem(
        T=eq.Free(375),
        P=eq.Free(1e6),
        feed=eq.Amounts({"propane": 1}),
        phases=[eq.Phase("F", kind="liquid")],
        row_tolerance=1e-11,
        critical=critical,
    )
    result = eq.solve_equilibrium(propane, problem)
    assert result.success, result.message
    Tc, Pc = result.temperature, result.pressure
    Ts = 0.99 * Tc
    upper = propane.molar_density_upper_bound(Ts)
    diagnostics = propane.density_diagnostics(
        Ts,
        Pc,
        nodes=[upper * i / 1024 for i in range(1, 1024)],
        density_tolerance=1e-10,
        log_tolerance=LOG_TOLERANCE,
        max_evaluations=20000,
        max_seconds=60,
        adaptive=False,
        finite_upper=upper * (1 - 1e-12),
        evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    assert len(diagnostics.spinodals) == 2, diagnostics.stop
    folds = sorted(diagnostics.spinodals, key=lambda s: s.molar_density)
    pv, pl = folds[0].pressure, folds[1].pressure
    assert pl > 0 and pv > pl
    ps300, psnear, psacetone = (
        saturation(propane, 300),
        saturation(propane, Ts),
        saturation(ms["C5"], 300),
    )
    _, keys = key_solves()
    mea_model, mea_problem = keys["case-7/flash-5918.5286-Pa"]
    mea = eq.solve_equilibrium(mea_model, mea_problem)
    assert mea.success, mea.message
    module, inputs = cross_check()
    c4q = next(q for q in inputs["requests"] if q["id"] == "C4-298.15")
    c4 = eq.solve_equilibrium(ms["C4"], module.prepare(inputs, c4q)[1])
    assert c4.success, c4.message
    n = len(ms["C4"].component_ids)
    c4x = max(
        [c4.mole_fractions[i * n : (i + 1) * n] for i in range(c4.phase_count)],
        key=lambda x: x[0],
    )
    cases = []

    def add(id, tier, model, T, P, phase, x, source, N=100):
        cases.append(
            {
                "id": id,
                "tier": tier,
                "model": model,
                "T_K": T,
                "P_Pa": P,
                "branch": phase,
                "x": x,
                "N": N,
                "warmups": 5,
                "source": source,
            }
        )

    add(
        "A-propane-vapor",
        "A",
        "C1",
        300,
        0.5 * ps300,
        "vapor",
        [1],
        "C1-300; 0.5 model P_sat",
    )
    add(
        "A-propane-liquid",
        "A",
        "C1",
        300,
        2 * ps300,
        "liquid",
        [1],
        "C1-300; 2 model P_sat",
    )
    add("A-propane-supercritical", "A", "C1", 400, 5e6, "vapor", [1], "#478")
    add(
        "A-methane-ethane-vapor",
        "A",
        "methane-ethane",
        250,
        1e6,
        "vapor",
        [0.5, 0.5],
        "bubble-point collect_engine(declare_only=True)",
    )
    add("A-C2-liquid", "A", "C2", 373.15, 2e6, "liquid", [0.5, 0.5], "C2-0.5")
    for branch in ("liquid", "vapor"):
        add(
            "B-propane-critical-" + branch,
            "B",
            "C1",
            Ts,
            psnear,
            branch,
            [1],
            "0.99 model T_c; model P_sat",
        )
    add(
        "B-propane-above-critical",
        "B",
        "C1",
        1.01 * Tc,
        Pc,
        "vapor",
        [1],
        "1.01 model T_c; model P_c",
    )
    # One thousandth of the spinodal pressure gap stays inside each metastable interval.
    add(
        "B-propane-metastable-liquid",
        "B",
        "C1",
        Ts,
        pl + 0.001 * (pv - pl),
        "liquid",
        [1],
        "liquid spinodal +0.001 spinodal pressure gap",
    )
    add(
        "B-propane-metastable-vapor",
        "B",
        "C1",
        Ts,
        pv - 0.001 * (pv - pl),
        "vapor",
        [1],
        "vapor spinodal -0.001 spinodal pressure gap",
    )
    add("B-propane-1-Pa", "B", "C1", 200, 1, "vapor", [1], "#478")
    add(
        "B-propane-100-MPa",
        "B",
        "C1",
        200,
        1e8,
        "liquid",
        [1],
        "#478; same 200 K extreme isotherm",
    )
    for branch, P in [("liquid", 1e5), ("vapor", 3000)]:
        add(
            "B-water-" + branch,
            "B",
            "water",
            298.15,
            P,
            branch,
            [1],
            "neutral water from C3 packet",
        )
    add(
        "B-C4-water-rich",
        "B",
        "C4",
        298.15,
        c4q["P"],
        "liquid",
        c4x,
        "C4-298.15 main water-rich equilibrium composition",
    )
    for molality in (1, 3, 6):
        request = next(q for q in inputs["requests"] if q["id"] == "C3-sigma-1")
        parameters = Parameters.from_mapping(
            module.document(inputs["models"]["C3"], request),
            components=ms["C3"].component_ids,
        )
        x = _core.formula_composition(
            [1, 0, 0], parameters._native, (1, 2), (1, 1), molality
        )
        add(
            f"B-C3-{molality}-mol-kg",
            "B",
            "C3",
            298.15,
            1e5,
            "liquid",
            x,
            "C3-sigma matched electrostatics; Engine formula_composition",
            N=30,
        )
    for branch in ("liquid", "vapor"):
        add(
            "B-C5-" + branch,
            "B",
            "C5",
            300,
            psacetone,
            branch,
            [1],
            "C5-300 model P_sat; matched comparator dipole conversion",
        )
    add(
        "B-MEA-liquid",
        "B",
        "MEA",
        318.15,
        5918.5286,
        "liquid",
        mea.mole_fractions[: len(mea_model.component_ids)],
        "time-key-solves case-7 main flash liquid",
        N=20,
    )
    upper300 = propane.molar_density_upper_bound(300)
    diagnostic300 = propane.density_diagnostics(
        300,
        1e5,
        nodes=[upper300 * i / 1024 for i in range(1, 1024)],
        density_tolerance=1e-10,
        log_tolerance=LOG_TOLERANCE,
        max_evaluations=20000,
        max_seconds=60,
        adaptive=False,
        finite_upper=upper300 * (1 - 1e-12),
        evaluation_order=_core.DerivativeOrder.PressureSlope,
    )
    vapor300 = max(s.pressure for s in diagnostic300.spinodals)
    add(
        "B-propane-absent-vapor",
        "B",
        "C1",
        300,
        1.001 * vapor300,
        "vapor",
        [1],
        "vapor spinodal at 300 K times 1.001; sole liquid-side root; state branch Unspecified and stable count 1",
    )
    applications = [
        {
            "id": id,
            "tier": "C",
            "operation": op,
            "N": N,
            "warmups": 1 if op == "fit" else 5,
            "source": source,
        }
        for id, op, N, source in [
            ("C1-saturation", "C1", 20, "C1-300"),
            ("C2-bubble", "C2", 10, "C2-0.5"),
            (
                "methane-ethane-bubble",
                "bubble",
                20,
                "collect_engine(declare_only=True)",
            ),
            *[
                (name, "key", 1, "tools/time-key-solves")
                for name in key_solves()[0].CASES
            ],
            (
                "methane-ethane-flash",
                "flash",
                10,
                "tp-flash performance.declarations()",
            ),
            (
                "phase-stability",
                "stability",
                1,
                "methane/ethane TP flash; default native TPD options",
            ),
            ("C2-pressure-actions", "actions-C2", 20, "C2-0.5; ln phi, d/dP"),
            ("C3-pressure-actions", "actions-C3", 10, "C3-sigma-1; ln phi, d/dP"),
            (
                "pure-saturation-objective-gradient",
                "objective",
                2,
                "pure-saturation-regression replay setup; packet optimizer point",
            ),
            (
                "pure-saturation-fit",
                "fit",
                1,
                "pure-saturation-regression retained origin start; full controls and observations",
            ),
        ]
    ]
    quantities = {
        "T_c_K": Tc,
        "P_c_Pa": Pc,
        "P_sat_300_Pa": ps300,
        "P_sat_0_99_Tc_Pa": psnear,
        "acetone_P_sat_300_Pa": psacetone,
        "spinodal_T_K": Ts,
        "vapor_spinodal_P_Pa": pv,
        "liquid_spinodal_P_Pa": pl,
        "vapor_spinodal_300_P_Pa": vapor300,
        "method": {
            "critical": "Engine CriticalDeclaration; same declaration as test_critical_candidate_matches_teqp",
            "saturation": "Engine solve_equilibrium(Pure.saturation(T))",
            "spinodal": "Engine density_diagnostics; 1024 linear density nodes; PressureSlope; log tolerance 1e-14; no adaptive scan",
        },
    }
    comparator_inputs = {k: inputs[k] for k in ("models", "constants")}
    comparator_inputs["models"]["water"] = {
        "components": [inputs["models"]["C3"]["components"][0]]
    }
    packet = json.loads(
        (
            ROOT / "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json"
        ).read_text()
    )
    components = []
    for cid in ms["methane-ethane"].component_ids:
        c = next(c for c in packet["components"] if c["component_id"] == cid)
        coefficients = {v["family"]: v["value"]["magnitude"] for v in c["coefficients"]}
        components.append(
            {
                "id": cid,
                "m": coefficients["segment_count"],
                "sigma": coefficients["segment_diameter"],
                "epsilon": coefficients["dispersion_energy_over_k"],
                "molar_mass": c["fixed"]["molar_mass"]["value"]["magnitude"],
            }
        )
    comparator_inputs["models"]["methane-ethane"] = {
        "components": components,
        "k": packet["pairs"][0]["coefficients"][0]["value"]["magnitude"],
        "source": "gross-2001-methane-ethane/3 packet in Engine component order",
    }
    document = {
        "gas_constant_J_mol_K": R,
        "comparator_inputs": comparator_inputs,
        "main_sha": json.loads((HERE / "inputs/cases.json").read_text())["main_sha"],
        "model_quantities": quantities,
        "cases": cases + applications,
        "accepted_starts": {
            name: eq.solve_equilibrium(*keys[name]).coordinates
            for name in key_solves()[0].ORDINARY
        },
    }
    (HERE / "inputs/cases.json").write_text(json.dumps(document, indent=2) + "\n")
    rows = []
    for case in cases:
        print("oracle", case["id"], flush=True)
        rows.extend(enumerate_roots(ms[case["model"]], case))
    write_csv(HERE / "results/oracle.csv", rows)


if __name__ == "__main__":
    if "--initialize" in sys.argv:
        main()
    parser = argparse.ArgumentParser()
    parser.add_argument("--initialize", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--release-metadata", type=Path)
    args = parser.parse_args()
    refresh(args.output, args.release_metadata)
