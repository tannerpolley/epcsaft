#!/usr/bin/env python3
"""Benzene/toluene rounding sensitivity (#73): Yu (2026) printed pure set vs the packet's full-precision Gross set.

Equimolar benzene/toluene, 105000 Pa, 350/370/400 K, k_ij = -0.003 in both. The limits are relative
differences < 1e-4 (density, total enthalpy), |d ln(phi_i)| < 1e-4 and bubble/dew differences < 0.01 K.
The liquid ln(phi_i) limit fails; the owner adopted the Gross set on 2026-09-26 with that deviation
accepted and recorded, so this check now fails only if the deviation leaves [4.75e-4, 4.85e-4).
"""

import csv
import json
import sys

import epcsaft
from epcsaft import Mixture, PressureActionDirection, PropertyObservable
from epcsaft.equilibrium import solve_equilibrium

from hda_inputs import ANALYSIS, PAIR, P, X, boundary_problem, parameters, thermochemistry

ACCEPTED_LN_PHI = (4.75e-4, 4.85e-4)


def yu_override():
    with (ANALYSIS / "inputs/yu_rounded_pure.csv").open() as handle:
        return {row["component"]: {"segment_count": float(row["m"]), "segment_diameter": float(row["sigma_A"]),
                                   "dispersion_energy_over_k": float(row["epsilon_over_k_K"])}
                for row in csv.DictReader(handle)}


def declarations():
    models = [Mixture(parameters(PAIR, override), thermochemistry=thermochemistry(PAIR)) for override in (None, yu_override())]
    return models, [boundary_problem("bubble", 365.0), boundary_problem("dew", 372.0)]


def state_values(m, T, phase):
    batch = m.pressure_actions(T, P, x=X, phase=phase, observables=[
        PropertyObservable.TotalEnthalpy, PropertyObservable.LogFugacityCoefficient],
        directions=[PressureActionDirection(temperature=1.0, composition=[0.0, 0.0])])
    if batch.molar_density is None:
        return {"available": False, "message": batch.density_message}
    values = [item.value for item in batch.results]
    return {"available": True, "branch": str(batch.diagnostics.branch), "roots": batch.diagnostics.stable_root_count,
            "rho": batch.molar_density, "h": values[0], "ln_phi_benzene": values[1], "ln_phi_toluene": values[2]}


def main():
    (gross, yu), problems = declarations()
    rows, within = [], True
    for T in (350.0, 370.0, 400.0):
        for phase in ("liquid", "vapor"):
            a, b = state_values(gross, T, phase), state_values(yu, T, phase)
            row = {"T_K": T, "phase": phase, "available": a["available"] and b["available"]}
            within &= row["available"]
            if row["available"]:
                row.update({"branch_gross": a["branch"], "branch_yu": b["branch"], "roots": a["roots"]})
                for key in ("rho", "h"):
                    row[f"{key}_gross"], row[f"{key}_yu"] = a[key], b[key]
                    row[f"rel_{key}"] = abs(b[key] - a[key]) / abs(a[key])
                    within &= row[f"rel_{key}"] < 1e-4
                for key in ("ln_phi_benzene", "ln_phi_toluene"):
                    row[f"{key}_gross"], row[f"{key}_yu"] = a[key], b[key]
                    row[f"abs_{key}"] = abs(b[key] - a[key])
            rows.append(row)
    temperatures = {}
    for problem, kind in zip(problems, ("bubble", "dew")):
        a, b = (solve_equilibrium(m, problem) for m in (gross, yu))
        a, b = (r.temperature if r.success else None for r in (a, b))
        temperatures[kind] = {"gross_K": a, "yu_K": b, "abs_K": None if a is None or b is None else abs(a - b)}
        within &= temperatures[kind]["abs_K"] is not None and temperatures[kind]["abs_K"] < 0.01
    with (ANALYSIS / "results/rounding_check.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    ln_phi = max(max(r["abs_ln_phi_benzene"], r["abs_ln_phi_toluene"]) for r in rows if r["available"])
    summary = {"pressure_Pa": P, "composition": X, "k_ij": -0.003, "boundary_temperatures": temperatures,
               "limits": {"rel_rho": 1e-4, "rel_h": 1e-4, "abs_ln_phi": 1e-4, "abs_T_K": 0.01},
               "max_rel_rho": max(r["rel_rho"] for r in rows if r["available"]),
               "max_rel_h": max(r["rel_h"] for r in rows if r["available"]),
               "max_abs_ln_phi": ln_phi, "accepted_ln_phi_band": ACCEPTED_LN_PHI,
               "density_enthalpy_boundaries_within_limits": within,
               "ln_phi_in_accepted_band": ACCEPTED_LN_PHI[0] <= ln_phi < ACCEPTED_LN_PHI[1],
               "decision": "Owner adoption 2026-09-26: the full-precision Gross pure set with k_ij = -0.003; the "
                           "liquid ln(phi_i) deviation from Yu's printed set is accepted and recorded; no refit.",
               "engine": epcsaft.provenance()}
    (ANALYSIS / "results/rounding_check.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    sys.exit(0 if within and summary["ln_phi_in_accepted_band"] else 1)


if __name__ == "__main__":
    main()
