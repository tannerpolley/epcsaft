"""Engine rows for the Gross (2005) quadrupole term against the retained Clapeyron rows and reference data (#574).

Run from the repository root after clapeyron_reference.jl: .venv/bin/python analyses/2026-polar-quadrupole/scripts/compare.py
Writes results/engine-rows.csv (every compared value) and results/summary.json (maxima and AADs). The nonpolar
Gross-Sadowski 2001 CO2 and benzene sets are the comparison models with their own parameters (ADR 0011).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import epcsaft
from epcsaft import Mixture, Parameters, Pure
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, solve_equilibrium

ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).resolve().parents[1] / "results"
GROSS = ROOT / "data/packets/gross-2005/1/parameters/parameters.json"
GS = ROOT / "data/packets/gross-sadowski-2001/3/parameters/parameters.json"


def model(names: tuple[str, ...], polar: bool = True, path: Path = GROSS) -> Mixture | Pure:
    document = json.loads(path.read_text())
    document["components"] = [c for c in document["components"] if c["component_id"] in names]
    document["pairs"] = [{"component_id_a": a, "component_id_b": b, "coefficients": [
        {"identity": f"{a}/{b}/k_ij", "family": "k_ij", "value": {"magnitude": 0.0, "unit": "dimensionless"}}]}
        for i, a in enumerate(names) for b in names[i + 1:]]
    if not polar:
        document["model_families"] = [f for f in document["model_families"] if f["kind"] != "polar"]
        for c in document["components"]:
            c["coefficients"] = [k for k in c["coefficients"] if not k["family"].startswith("quadrupole")]
    parameters = Parameters.from_mapping(document, components=names)
    return Pure(parameters) if len(names) == 1 else Mixture(parameters)


def problem(T: float, guess: float, feed: dict[str, float]) -> Problem:
    guess_x = [1.0] if len(feed) == 1 else None
    return Problem(phases=[Phase("L", kind="liquid", composition_guess=guess_x),
                           Phase("V", kind="vapor", amount=Pinned(0.0), composition_guess=guess_x)],
                   T=T, P=Free(guess), feed=Amounts(feed), row_tolerance=1.0e-12)


def solve(m: Mixture | Pure, T: float, guess: float, feed: dict[str, float]) -> tuple[float, list[float]]:
    result = solve_equilibrium(m, problem(T, guess, feed))
    if not result.success:
        raise RuntimeError(f"{feed} at {T} K: {result.message}")
    return result.pressure, result.phase_mole_fractions(1)


def dippr(T: float) -> float:
    # DIPPR-101 form, Pa; benzene coefficients from the #574 readiness probe, not checked against the DIPPR
    # compilation (Gross 2005 data ref. 2).
    return math.exp(83.107 - 6486.2 / T - 9.2194 * math.log(T) + 6.9844e-6 * T**2)


def read(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open() as handle:
        return list(csv.DictReader(handle))


def max_rel(rows: list[dict], kind: str) -> float:
    return max(abs(r["engine"] / r["comparator"] - 1.0) for r in rows if r["kind"] == kind)


def aad(rows: list[dict], kind: str, species: str) -> float:
    group = [r for r in rows if r["kind"] == kind and r["case"].startswith(species + " ")]
    return 100.0 * sum(abs(r["engine"] / r["reference"] - 1.0) for r in group) / len(group)


def declarations():
    saturation = read("clapeyron-saturation-rows.csv")
    result: list[tuple[Mixture | Pure, list[Problem]]] = [
        (m, [problem(float(r["T_K"]), float(r["p_sat_clapeyron_Pa"]), {c: 1.0}) for r in saturation if r["species"] == c])
        for c in ("carbon-dioxide", "benzene") for m in (model((c,)), model((c,), False, GS))]
    bubbles = [problem(float(r["T_K"]), float(r["p_bubble_clapeyron_Pa"]), {"carbon-dioxide": 1.0 - float(r["x_nitrogen"]),
               "nitrogen": float(r["x_nitrogen"])}) for r in read("clapeyron-bubble-rows.csv")]
    return result + [(model(("carbon-dioxide", "nitrogen")), bubbles)]


def main() -> None:
    rows = []
    for r in read("clapeyron-term-rows.csv"):
        names, x = tuple(r["species"].split("+")), tuple(float(v) for v in r["x"].split())
        T, rho = float(r["T_K"]), float(r["rho_mol_m3"])
        a = [(model(names, p).state(T, rho=rho) if len(names) == 1 else model(names, p).state(T, rho=rho, x=x)).residual_helmholtz
             for p in (True, False)]
        rows.append({"kind": "a_qq", "case": f"{r['species']} {r['x']} {T:g} K {rho:g} mol/m3", "engine": a[0] - a[1],
                     "comparator": float(r["a_qq_clapeyron"]), "reference": ""})
    for r in read("feos-term-rows.csv"):
        names, x = tuple(r["species"].split("+")), tuple(float(v) for v in r["x"].split())
        T, rho = float(r["T_K"]), float(r["rho_mol_m3"])
        a = [(model(names, p).state(T, rho=rho) if len(names) == 1 else model(names, p).state(T, rho=rho, x=x)).residual_helmholtz
             for p in (True, False)]
        rows.append({"kind": "a_qq_feos", "case": f"{r['species']} {r['x']} {T:g} K {rho:g} mol/m3", "engine": a[0] - a[1],
                     "comparator": float(r["a_qq_feos"]), "reference": ""})
    nonpolar = {c: model((c,), False, GS) for c in ("carbon-dioxide", "benzene")}
    for r in read("clapeyron-saturation-rows.csv"):
        c, T, p_c = r["species"], float(r["T_K"]), float(r["p_sat_clapeyron_Pa"])
        reference = float(r["p_sat_span_wagner_Pa"]) if c == "carbon-dioxide" else dippr(T)
        rows.append({"kind": "p_sat", "case": f"{c} {T:g} K", "engine": solve(model((c,)), T, p_c, {c: 1.0})[0],
                     "comparator": p_c, "reference": reference})
        rows.append({"kind": "p_sat_nonpolar_gs2001", "case": f"{c} {T:g} K", "engine": solve(nonpolar[c], T, p_c, {c: 1.0})[0],
                     "comparator": "", "reference": reference})
    for r in read("clapeyron-bubble-rows.csv"):
        T, x, p_c = float(r["T_K"]), float(r["x_nitrogen"]), float(r["p_bubble_clapeyron_Pa"])
        p, y = solve(model(("carbon-dioxide", "nitrogen")), T, p_c, {"carbon-dioxide": 1.0 - x, "nitrogen": x})
        rows.append({"kind": "p_bubble", "case": f"x_N2 {x:g} {T:g} K", "engine": p, "comparator": p_c, "reference": ""})
        rows.append({"kind": "y_nitrogen", "case": f"x_N2 {x:g} {T:g} K", "engine": y[1], "comparator": float(r["y_nitrogen_clapeyron"]),
                     "reference": ""})
    with (RESULTS / "engine-rows.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["kind", "case", "engine", "comparator", "reference"], lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()} for row in rows)

    summary = {"engine_provenance": list(epcsaft.provenance()),
               "max_relative_vs_clapeyron": {k: max_rel(rows, k) for k in ("a_qq", "p_sat", "p_bubble", "y_nitrogen")},
               "max_relative_vs_feos": {"a_qq": max_rel(rows, "a_qq_feos")},
               "aad_percent": {"carbon-dioxide quadrupole vs Span-Wagner 220-300 K": aad(rows, "p_sat", "carbon-dioxide"),
                               "carbon-dioxide nonpolar GS2001 vs Span-Wagner 220-300 K": aad(rows, "p_sat_nonpolar_gs2001", "carbon-dioxide"),
                               "benzene quadrupole vs DIPPR-101 290-540 K": aad(rows, "p_sat", "benzene"),
                               "benzene nonpolar GS2001 vs DIPPR-101 290-540 K": aad(rows, "p_sat_nonpolar_gs2001", "benzene")}}
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
