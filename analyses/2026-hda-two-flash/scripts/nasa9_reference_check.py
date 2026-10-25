#!/usr/bin/env python3
"""Independent NASA9 evaluation of the HDA species records against the issue #73 table.

Evaluates cp, h, s of nasa-glenn-2002-hda-species/1 with NASA/TP-2002-211556 Eqs. (1)-(3) at
R = 8.314510, checks Delta_r H(298.15 K) of C7H8 + H2 -> C6H6 + CH4, then reads the same records
through the Engine ideal-gas path, epcsaft.Pure(...).state(T, rho=p0/(R T)), with every
coefficient rescaled once by 8.314510 / epcsaft.GAS_CONSTANT_J_PER_MOL_K in epcsaft.ThermochemistryPacket.
"""

import csv
import hashlib
import json
import math
import sys

import epcsaft

from hda_inputs import ANALYSIS, NASA9, ROOT, parameters, thermochemistry

R_SOURCE = 8.314510  # J/(mol K), NASA/TP-2002-211556; the independent evaluator's own gas constant

TOL = {"cp": 1e-6, "h": 1e-3, "s": 1e-6}
# Issue #73 reference table: (T_K, cp J/(mol K), h J/mol, s J/(mol K)).
EXPECTED = {
    "hydrogen": [(298.15, 28.836346808, -0.000002717, 130.681014278),
                 (600.0, 29.318148662, 8810.697863106, 151.078576537),
                 (999.999, 30.206308334, 20678.985281782, 166.216666798),
                 (1000.001, 30.206315551, 20679.045711203, 166.216727227)],
    "methane": [(298.15, 35.691325971, -74599.999986132, 186.371072175),
                (600.0, 52.690762017, -61409.585225929, 216.230772386),
                (999.999, 73.676417107, -35915.490589879, 248.330854028),
                (1000.001, 73.676504486, -35915.343437105, 248.331001094)],
    "benzene": [(298.15, 81.934906370, 82879.999913057, 269.159049522),
                (600.0, 159.386095767, 120371.868728919, 352.837937461),
                (999.999, 211.033740544, 195714.942294757, 447.936824030),
                (1000.001, 211.033918424, 195715.363324160, 447.937247761)],
    "toluene": [(298.15, 103.279830279, 50169.999997183, 320.188413034),
                (600.0, 195.752836731, 96427.557022348, 423.559518225),
                (999.999, 260.911511092, 189271.387719278, 540.713963974),
                (1000.001, 260.911736254, 189271.909555163, 540.714485488)],
}
DELTA_R_H_EXPECTED, DELTA_R_H_TOL = -41890.0, 1e-2  # J/mol, toluene + H2 -> benzene + CH4 at 298.15 K


def declarations():
    return {sid: epcsaft.Pure(parameters([sid]), thermochemistry=thermochemistry([sid])) for sid in EXPECTED}


def interval(record, T):
    """Low interval below a shared bound, high interval above; exact bounds are never queried."""
    return next(item for item in record["intervals"] if item["lower_K"] <= T < item["upper_K"])


def nasa9(a, T):
    """Return cp [J/(mol K)], h [J/mol], s at 1 bar [J/(mol K)]."""
    lnT = math.log(T)
    cp = a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4
    h = (-a[0] / T**2 + a[1] * lnT / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3
         + a[5] * T**3 / 4 + a[6] * T**4 / 5 + a[7] / T)
    s = (-a[0] / T**2 / 2 - a[1] / T + a[2] * lnT + a[3] * T + a[4] * T**2 / 2
         + a[5] * T**3 / 3 + a[6] * T**4 / 4 + a[8])
    return R_SOURCE * cp, R_SOURCE * T * h, R_SOURCE * s


def main():
    doc = json.loads(NASA9.read_text())
    records = {r["species_id"]: r for r in doc["species"]}
    models, r_engine, p0 = declarations(), epcsaft.GAS_CONSTANT_J_PER_MOL_K, doc["standard_pressure_Pa"]
    rows = []
    for sid, table in EXPECTED.items():
        for T, *expected in table:
            item = interval(records[sid], T)
            computed = nasa9(item["coefficients"], T)
            state = models[sid].state(T, rho=p0 / (r_engine * T))
            engine = (state.ideal_isobaric_heat_capacity, state.ideal_enthalpy, state.ideal_entropy)
            row = {"species": sid, "T_K": T, "interval_lower_K": item["lower_K"], "interval_upper_K": item["upper_K"]}
            for key, e, c, g in zip(TOL, expected, computed, engine):
                row |= {f"{key}_expected": e, f"{key}_computed": c, f"{key}_abs_error": abs(c - e),
                        f"{key}_engine": g, f"{key}_engine_abs_error": abs(g - e)}
            row["pass"] = all(row[f"{k}_abs_error"] <= t for k, t in TOL.items())
            row["engine_pass"] = all(row[f"{k}_engine_abs_error"] <= t for k, t in TOL.items())
            rows.append(row)
    h298 = {sid: nasa9(interval(records[sid], 298.15)["coefficients"], 298.15)[1] for sid in records}
    delta_r_h = h298["benzene"] + h298["methane"] - h298["toluene"] - h298["hydrogen"]
    with (ANALYSIS / "results/nasa9_reference_check.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "records": str(NASA9.relative_to(ROOT)), "records_sha256": hashlib.sha256(NASA9.read_bytes()).hexdigest(),
        "source_gas_constant": R_SOURCE, "engine_gas_constant": r_engine, "coefficient_scale": R_SOURCE / r_engine,
        "standard_pressure_Pa": p0, "tolerances": TOL, "rows": len(rows),
        "max_abs_error": {k: max(r[f"{k}_abs_error"] for r in rows) for k in TOL},
        "max_engine_abs_error": {k: max(r[f"{k}_engine_abs_error"] for r in rows) for k in TOL},
        "all_pass": all(r["pass"] for r in rows), "engine_all_pass": all(r["engine_pass"] for r in rows),
        "delta_r_h_298_J_per_mol": delta_r_h, "delta_r_h_expected": DELTA_R_H_EXPECTED,
        "delta_r_h_pass": abs(delta_r_h - DELTA_R_H_EXPECTED) <= DELTA_R_H_TOL,
        "engine": epcsaft.provenance(),
        "claim_limit": "Checks record transcription and the NASA9 forms only; no comparison with calorimetric data.",
    }
    (ANALYSIS / "results/nasa9_reference_check.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    sys.exit(0 if summary["all_pass"] and summary["engine_all_pass"] and summary["delta_r_h_pass"] else 1)


if __name__ == "__main__":
    main()
