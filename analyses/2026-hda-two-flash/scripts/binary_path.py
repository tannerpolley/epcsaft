#!/usr/bin/env python3
"""Equimolar benzene/toluene at 105000 Pa (#73): both density branches on a 5 K grid from 300 to 450 K
plus the native bubble and dew temperatures, with root counts, branch tags, the recomputed
equal-fugacity residual at the two crossings, and a central difference of h that only verifies the
exact cp. Fails unless every row has two stable roots, both crossing residuals are <= 1e-8 and the
crossings are within 1e-6 K of the retained 367.43107107297914 / 373.50696978716286 K."""

import csv
import json
import math
import sys

import epcsaft
from epcsaft import Mixture
from epcsaft.equilibrium import solve_equilibrium

from hda_inputs import ANALYSIS, PAIR, P, X, boundary_problem, parameters, thermochemistry

RETAINED_K = {"bubble": 367.43107107297914, "dew": 373.50696978716286}


def declarations():
    return (Mixture(parameters(PAIR), thermochemistry=thermochemistry(PAIR)),
            {kind: boundary_problem(kind, guess) for kind, guess in (("bubble", 365.0), ("dew", 372.0))})


def main():
    model, problems = declarations()
    crossings = {}
    for kind, problem in problems.items():
        result = solve_equilibrium(model, problem)
        T = result.temperature
        x, y = result.mole_fractions[:2], result.mole_fractions[2:]
        liquid, vapor = model.state(T, P=P, x=x, phase="liquid"), model.state(T, P=P, x=y, phase="vapor")
        crossings[kind] = {"success": result.success, "T_K": T, "retained_T_K": RETAINED_K[kind],
                           "abs_T_change_K": abs(T - RETAINED_K[kind]), "x": list(x), "y": list(y),
                           "rho_liquid": liquid.molar_density, "rho_vapor": vapor.molar_density,
                           "roots": [liquid.stable_root_count, vapor.stable_root_count],
                           "fugacity_residual": max(abs(math.log(x[i]) + liquid.log_fugacity_coefficient[i]
                                                        - math.log(y[i]) - vapor.log_fugacity_coefficient[i]) for i in range(2))}
    rows = []
    for T in sorted({300.0 + 5.0 * i for i in range(31)} | {c["T_K"] for c in crossings.values()}):
        for phase in ("liquid", "vapor"):
            try:
                state = model.state(T, P=P, x=X, phase=phase)
            except Exception as error:  # an unavailable branch is retained, never substituted
                rows.append({"T_K": T, "phase": phase, "available": False, "reason": str(error)})
                continue
            h = lambda t: model.state(t, P=P, x=X, phase=phase).total_enthalpy  # noqa: E731
            rows.append({"T_K": T, "phase": phase, "available": True, "branch": str(state.density_branch),
                         "roots": state.stable_root_count, "rho_mol_m3": state.molar_density, "h_J_mol": state.total_enthalpy,
                         "cp_J_molK": state.total_isobaric_heat_capacity,
                         "dh_dT_fd_minus_cp": (h(T + 1e-3) - h(T - 1e-3)) / 2e-3 - state.total_isobaric_heat_capacity,
                         "ln_phi_benzene": state.log_fugacity_coefficient[0], "ln_phi_toluene": state.log_fugacity_coefficient[1]})
    with (ANALYSIS / "results/binary_path.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    available = [r for r in rows if r["available"]]
    summary = {"pressure_Pa": P, "composition": X, "k_ij": -0.003, "crossings": crossings, "rows": len(rows),
               "unavailable": [(r["T_K"], r["phase"]) for r in rows if not r["available"]],
               "rows_without_two_stable_roots": [(r["T_K"], r["phase"]) for r in available if r["roots"] != 2],
               "max_abs_dh_dT_fd_minus_cp": max(abs(r["dh_dT_fd_minus_cp"]) for r in available),
               "engine": epcsaft.provenance()}
    summary["pass"] = (not summary["unavailable"] and not summary["rows_without_two_stable_roots"]
                       and all(c["success"] and c["fugacity_residual"] <= 1e-8 and c["abs_T_change_K"] <= 1e-6
                               for c in crossings.values()))
    (ANALYSIS / "results/binary_path.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    sys.exit(0 if summary["pass"] else 1)


if __name__ == "__main__":
    main()
