"""Cameretti, Sadowski and Mollerup (2005) with the Engine: every property below is a native solve.

Writes results/*.csv and results/run-identity.json, or the same files to the directory given as the
only argument (tools/check-equilibrium-corpora). Four model variants are evaluated: ion table
(corrected 8944 / original 3359) x ionic packing diameter (0.88 sigma / sigma). criterion-violations.csv
lists every bubble result whose native criterion is not finite and positive, or that succeeds with a
residual outside it; it has no rows when every result holds.
"""
import csv
import json
import math
import os
import sys
import time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import epcsaft
import epcsaft.equilibrium as q

import model as M
import water_reference as W

HERE = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "results"
sys.path.insert(0, str(HERE.parents[2] / "tools"))
from engine_identity import loaded_identity

MW = M.WATER["M"] * 1e-3
P_DENSITY = 101325.0  # Pa; Figures 1, 3 and 4 densities at ambient pressure
VARIANTS = {"corrected/0.88": ("IONS", "0.88"), "corrected/sigma": ("IONS", "sigma"),
            "original/0.88": ("IONS_ORIGINAL", "0.88"), "original/sigma": ("IONS_ORIGINAL", "sigma")}


@lru_cache(maxsize=None)
def mixture(ions, variant):
    table, packing = VARIANTS[variant]
    return M.mixture(list(ions), packing, getattr(M, table))


def composition(amounts):
    """Mole fractions of water and ions for formula-unit molalities {ion: mol/kg water}."""
    total = 1 / MW + sum(amounts.values())
    return [1 / MW / total] + [v / total for v in amounts.values()]


def bubble(variant, T, amounts):
    """Native bubble pressure with a water-only vapor: the Engine default start, then a declared ideal-solution start."""
    ions = tuple(amounts)
    x = composition(amounts)
    feed = q.Molalities(dict(amounts), solvent=["water"], solvent_mass_kg=1.0)
    problem = q.Problem(phases=[q.Phase("L", kind="liquid"), q.Phase("V", kind="vapor", support=("water",), amount=q.Pinned(0.0))],
                        T=T, P=q.Free(), feed=feed)
    model = mixture(ions, variant)
    default = q.solve_equilibrium(model, problem)
    start = q.EquilibriumStart(pressure=W.saturation_pressure(T) * x[0], compositions=[x, [1.0]], packing_fractions=[0.45, None])
    result = q.solve_equilibrium_starts(model, problem, [start], q.EquilibriumWorkLimits(attempts=1)).results[0]
    return dict(default_success=default.success, default_status=default.solver_status, default_iterations=default.iterations,
                default_residual=default.convergence.residual, default_criterion=default.convergence.criterion.value,
                declared_success=result.success, declared_status=result.solver_status,
                declared_iterations=result.iterations, declared_residual=result.convergence.residual,
                declared_criterion=result.convergence.criterion.value, p_kPa=result.pressure / 1e3,
                liquid_density_mol_m3=result.molar_densities[0], vapor_density_mol_m3=result.molar_densities[1],
                liquid_amount_mol=result.amounts[0], vapor_amount_mol=result.amounts[1])


def density(variant, T, amounts):
    ions = tuple(amounts)
    x = composition(amounts)
    state = mixture(ions, variant).state(T, P=P_DENSITY, x=x, phase="liquid")
    mass = x[0] * M.WATER["M"] + sum(xi * M.IONS[i][3] for xi, i in zip(x[1:], ions))
    return state.molar_density * mass * 1e-6  # g/cm3


def salt_amounts(salt, x_salt):
    cation, anion, nc, na = M.SALTS[salt]
    m = x_salt / (1 - x_salt) / MW
    return {cation: nc * m, anion: na * m}


def evaluate(task):
    kind, variant, T, amounts = task
    return bubble(variant, T, dict(amounts)) if kind == "p" else density(variant, T, dict(amounts))


def unmet(record):
    """Sides of a bubble result whose criterion is not finite and positive, or that succeed with a residual outside it."""
    return [side for side in ("default", "declared") if not 0 < record[f"{side}_criterion"] < math.inf
            or record[f"{side}_success"] and not 0 <= record[f"{side}_residual"] <= record[f"{side}_criterion"]]


def pressure(record):
    return record["p_kPa"] if record["declared_success"] else math.nan


def write(name, rows):
    with (OUT / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else [], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def water_rows():
    pure = M.water()
    rows = []
    for T in [278.15 + 5 * i for i in range(24)]:
        result = q.solve_equilibrium(pure, pure.saturation(T))
        liquid = pure.state(T, P=P_DENSITY, phase="liquid")
        rows.append(dict(temperature_K=T, engine_saturation_pressure_kPa=result.pressure / 1e3,
                         iapws_saturation_pressure_kPa=W.saturation_pressure(T) / 1e3,
                         engine_saturated_liquid_density_g_cm3=result.molar_densities[0] * MW * 1e-3,
                         iapws_saturated_liquid_density_g_cm3=W.saturated_liquid_density(T) / 1e3,
                         engine_liquid_density_101kPa_g_cm3=liquid.molar_density * MW * 1e-3, converged=result.success))
    return rows


def table5_rows():
    """Both printed versions of Table 5; the 8944 correction prints pressures in units of 100 kPa."""
    original = list(csv.reader((HERE / "results/tables/table-005/source/table-005.csv").open()))[2:]
    corrected = list(csv.DictReader((HERE / "inputs/source/table-005-corrected.csv").open()))
    rows = []
    for k, (line, fix) in enumerate(zip(original, corrected)):
        T = float(line[0])
        for system, cols, ions in (("NaCl+KBr", (1, 2, 3, 4), ("Na+", "Cl-", "K+", "Br-")), ("NaBr+KCl", (6, 7, 8, 9), ("Na+", "Br-", "K+", "Cl-"))):
            _, _, p_exp, p_mod = (float(line[c]) for c in cols)
            prefix = "nacl_kbr" if system == "NaCl+KBr" else "nabr_kcl"
            # Molalities from the correction, which undoes the original's swapped NaCl+KBr rows at 328.15 K.
            m1, m2 = (float(fix[k]) for k in (("nacl_m", "kbr_m") if system == "NaCl+KBr" else ("nabr_m", "kcl_m")))
            amounts = {"Na+": m1, "K+": m2, ions[1]: m1, ions[3]: m2}
            rows.append(dict(row=k + 1, system=system, temperature_K=T, molality_1=m1, molality_2=m2, p_exp_kPa=p_exp,
                             p_model_original_kPa=p_mod, p_model_corrected_kPa=100 * float(fix[f"{prefix}_p_mod"]),
                             p_exp_corrected_kPa=100 * float(fix[f"{prefix}_p_exp"]), amounts=amounts))
    return rows


def run():
    start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    tasks, keys = [], []
    water = water_rows()
    t5 = table5_rows()
    for v in VARIANTS:
        for r in t5:
            tasks.append(("p", v, r["temperature_K"], tuple(r["amounts"].items())))
            keys.append(("t5", v, r["row"], r["system"]))
    lines = list(csv.DictReader((HERE / "results/digitized-lines.csv").open()))
    curves = {}
    for r in lines:
        curves.setdefault((r["figure"], r["salt"], r["quantity"], float(r["temperature_K"])), []).append(r)
    samples = []
    for (figure, salt, quantity, T), pts in curves.items():
        pts.sort(key=lambda r: float(r["x_salt"]))
        lo, hi = float(pts[0]["x_salt"]), float(pts[-1]["x_salt"])
        for i in range(8):  # eight evenly spaced positions along each recovered line
            target = lo + (hi - lo) * i / 7
            r = min(pts, key=lambda r: abs(float(r["x_salt"]) - target))
            samples.append(r)
    for v in VARIANTS:
        for r in samples:
            tasks.append((r["quantity"], v, float(r["temperature_K"]), tuple(salt_amounts(r["salt"], float(r["x_salt"])).items())))
            keys.append(("line", v, r))
    points = list(csv.DictReader((HERE / "results/digitized-points.csv").open()))
    for v in VARIANTS:
        for r in points:
            tasks.append((r["quantity"], v, float(r["temperature_K"]), tuple(salt_amounts(r["salt"], float(r["x_salt"])).items())))
            keys.append(("point", v, r))
    unique = list(dict.fromkeys(tasks))
    with Pool(os.cpu_count()) as pool:
        solved = dict(zip(unique, pool.map(evaluate, unique, chunksize=4)))
    solve_seconds = time.time() - start
    values = [solved[t] if t[0] == "rho" else pressure(solved[t]) for t in tasks]

    t5_out, line_out, point_out = [], [], []
    for key, value in zip(keys, values):
        if key[0] == "t5":
            r = next(x for x in t5 if x["row"] == key[2] and x["system"] == key[3])
            t5_out.append({**{k: v for k, v in r.items() if k != "amounts"}, "variant": key[1], "p_engine_kPa": value})
        elif key[0] == "line":
            r = key[2]
            line_out.append(dict(variant=key[1], figure=r["figure"], salt=r["salt"], quantity=r["quantity"], temperature_K=r["temperature_K"],
                                 x_salt=r["x_salt"], paper_model=r["value"], engine=value, value_per_pixel=r["value_per_pixel"], pixel_x=r["pixel_x"], pixel_y=r["pixel_y"]))
        else:
            r = key[2]
            point_out.append(dict(variant=key[1], figure=r["figure"], salt=r["salt"], quantity=r["quantity"], temperature_K=r["temperature_K"],
                                  x_salt=r["x_salt"], experiment=r["value"], engine=value, value_per_pixel=r["value_per_pixel"], nearest_line_px=r["nearest_line_px"]))
    write("water.csv", water)
    write("table-5.csv", t5_out)
    write("figure-lines.csv", line_out)
    write("figure-points.csv", point_out)
    write("activity-298K.csv", activity_rows())
    write_start_records(solved)
    kcl_rows, kcl = kcl_misprint(t5)
    write("table-5-kcl-misprint.csv", kcl_rows)
    write("criterion-violations.csv", [dict(variant=v, temperature_K=repr(T), ion_molalities_mol_kg=";".join(f"{ion}={m!r}" for ion, m in a),
                                            sides=";".join(sides), **record)
                                       for (kind, v, T, a), record in {**solved, **kcl}.items() if kind == "p" and (sides := unmet(record))])
    (OUT / "run-identity.json").write_text(json.dumps(dict(
        engine=loaded_identity(), processes=os.cpu_count(), native_solves=len(unique) + sum(t[0] == "p" for t in unique) + 2 * len(kcl),
        wall_seconds=round(time.time() - start, 1), parallel_solve_seconds=round(solve_seconds, 1)), indent=1) + "\n")


def kcl_misprint(t5):
    """Hypothesis only: Table 5 prints KCl = 0.009 mol/kg in one NaBr+KCl row; rerun it at 0.999 and 1.009 mol/kg.
    Returns the rows and the bubble results keyed like run()'s tasks."""
    r = next(r for r in t5 if r["system"] == "NaBr+KCl" and r["molality_2"] == 0.009)
    m1, T = r["molality_1"], r["temperature_K"]
    solved = {("p", v, T, (("Na+", m1), ("K+", m), ("Br-", m1), ("Cl-", m))): bubble(v, T, {"Na+": m1, "K+": m, "Br-": m1, "Cl-": m})
              for v in VARIANTS for m in (0.009, 0.999, 1.009)}
    return [dict(variant=v, temperature_K=T, nabr_molality=m1, kcl_molality=dict(a)["K+"], p_exp_kPa=r["p_exp_kPa"],
                 p_model_original_kPa=r["p_model_original_kPa"], p_engine_kPa=pressure(record)) for (_, v, _, a), record in solved.items()], solved


BANDS = [(0, 1), (1, 2), (2, 4), (4, 8), (8, math.inf)]  # total salt molality, formula units per kg water


def write_start_records(solved):
    """Every bubble state the default start fails, with its exact inputs and the declared-start result, and the failure fraction by band."""
    failures, counts = [], {}
    for (kind, variant, T, amounts), record in solved.items():
        if kind != "p":
            continue
        m_salt = sum(m for ion, m in amounts if ion.endswith("-"))  # anion molality = salt formula units for 1:1, 2:1 and mixed 1:1
        band = next(b for b in BANDS if b[0] <= m_salt < b[1])
        n = counts.setdefault((variant, band), [0, 0, 0.0, True, 0.0])
        n[0] += 1
        n[2] = max(n[2], record["declared_residual"])
        n[3] &= record["declared_success"]
        n[4] = max(n[4], record["declared_criterion"])
        if not record["default_success"]:
            n[1] += 1
            failures.append(dict(variant=variant, temperature_K=repr(T), salt_molality_mol_kg=m_salt,
                                 ion_molalities_mol_kg=";".join(f"{ion}={m!r}" for ion, m in amounts), **record))
    failures.sort(key=lambda r: (r["variant"], r["salt_molality_mol_kg"], float(r["temperature_K"])))
    write("default-start-failures.csv", failures)
    write("default-start-summary.csv", [dict(variant=v, salt_molality_band_mol_kg=f"{lo}-{hi}", states=n[0], default_failures=n[1],
                                             declared_all_succeeded=n[3], declared_max_residual=n[2],
                                             declared_max_criterion=n[4])
                                        for (v, (lo, hi)), n in sorted(counts.items())])


def activity_rows():
    """Molality-scale mean ionic activity and osmotic coefficients at 298.15 K, 1 bar against Hamer and Wu (1972)."""
    path = HERE.parents[2] / "data/packets/hamer-wu-1972-aqueous-alkali-halides/1/observations/hamer-wu-1972-aqueous-alkali-halides.csv"
    rows = []
    for obs in csv.DictReader(path.open()):
        cation, anion, _, _ = M.SALTS[obs["salt"]]
        m = float(obs["molality_mol_kg"])
        for v in VARIANTS:
            reference = reference_for(cation, anion, v)
            activity = reference.activity(m)
            rows.append(dict(variant=v, salt=obs["salt"], molality_mol_kg=m, gamma_hamer_wu=float(obs["gamma_pm_m"]),
                             gamma_engine=activity.mean_ionic_molality_scale, osmotic_engine=activity.osmotic_coefficient))
    return rows


@lru_cache(maxsize=None)
def reference_for(cation, anion, variant):
    return mixture((cation, anion), variant).reference(298.15, 1e5, solvent_fractions=[1, 0, 0], ions=(cation, anion), stoichiometry=(1, 1))


if __name__ == "__main__":
    run()
