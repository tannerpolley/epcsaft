"""Evaluate pure-MEA saturation with the Baygi 2015 2B parameters under the sigma^3 and d^3 conventions (#165).

Fakouri Baygi & Pahlavanzadeh 2015 Eq. 4 (p. 791) puts d_ij^3 in the association strength; the Engine
puts sigma_ij^3. Delta is linear in kappa, so at fixed T the d^3 model equals the Engine model with
kappa(T) = kappa (d/sigma)^3: exact for saturation at a stated temperature, not for T-derivatives.
Baygi's data are the Avlund 2008 DIPPR correlations of Table 1 (Eqs. 9-10) on their 303.15-443.15 K fit
range, with AAD 0.62 % (p_sat) and 0.12 % (rho_L,sat) in Table 2. The packet conversion is the constant
kappa that best reproduces Baygi's d^3 model (Eq. 8 form) over the packet domain, 293.15-393.15 K.
"""

import copy
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path

import epcsaft
from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem, solve_equilibrium

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
PACKET = ROOT / "data/packets/mea-co2-h2o-nine-species-estimation/4/parameters/parameters.json"
CONVERTED = ROOT / "data/packets/mea-co2-h2o-nine-species-estimation/5/parameters/parameters.json"
MEA, WATER = "monoethanolamine", "water"
EDGE = "association/monoethanolamine/a/monoethanolamine/b/volume"
GRID = [293.15 + 5.0 * k for k in range(31)]
BAYGI = (303.15, 443.15)  # Baygi's fit range
DOMAIN = (293.15, 393.15)  # packet domain


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dippr_psat(T: float) -> float:  # Baygi Table 1 / Eq. 9, Pa
    return math.exp(92.624 - 10367.0 / T - 9.4699 * math.log(T) + 1.9e-18 * T**6)


def dippr_rho(T: float) -> float:  # Baygi Table 1 / Eq. 10, mol/L
    return 1.0011 / 0.22523 ** (1.0 + (1.0 - T / 678.2) ** 0.21515)


def ratio(epsilon: float, T: float) -> float:  # d/sigma, Gross-Sadowski 2001 Eq. 3 (Baygi Eq. 5)
    return 1.0 - 0.12 * math.exp(-3.0 * epsilon / T)


def water_sigma(document: dict, T: float) -> float:
    c = next(c for c in document["correlations"] if c["component_id"] == WATER and c["family"] == "segment_diameter")
    return c["constant"]["value"]["magnitude"] + sum(
        t["amplitude"]["value"]["magnitude"] * math.exp(t["exponent_coefficient"]["value"]["magnitude"] * T) for t in c["terms"])


def coefficient(document: dict, component: str, family: str) -> float:
    record = next(c for c in document["components"] if c["component_id"] == component)
    return next(k["value"]["magnitude"] for k in record["coefficients"] if k["family"] == family)


def with_kappa(document: dict, kappa: float) -> epcsaft.Pure:
    mapping = copy.deepcopy(document)
    next(e for e in mapping["topology"]["edges"] if e["volume"]["identity"] == EDGE)["volume"]["value"]["magnitude"] = kappa
    return epcsaft.Pure(epcsaft.Parameters.from_mapping(mapping, components=[MEA]))


def saturation(model: epcsaft.Pure, T: float) -> tuple[float, float]:
    result = solve_equilibrium(model, Problem(
        phases=[Phase("L", kind="liquid", composition_guess=[1.0]), Phase("V", kind="vapor", amount=Pinned(0.0), composition_guess=[1.0])],
        T=T, P=Free(dippr_psat(T)), feed=Amounts({MEA: 1.0}), row_tolerance=1.0e-10))
    if not result.success:
        raise RuntimeError(f"{T} K: {result.message}")
    return result.pressure, result.molar_densities[0] / 1000.0  # Pa, mol/L


def within(rows: list[dict], bounds: tuple[float, float]) -> list[dict]:
    return [r for r in rows if bounds[0] - 1e-9 <= r["T_K"] <= bounds[1] + 1e-9]


def aad(rows: list[dict], key: str, reference: str = "dippr") -> dict:
    p = [r[f"p_{key}_pa"] / r[f"p_{reference}_pa"] - 1.0 for r in rows]
    rho = [abs(r[f"rho_{key}_mol_l"] / r[f"rho_{reference}_mol_l"] - 1.0) for r in rows]
    return {"p_sat_aad_percent": 100 * sum(map(abs, p)) / len(p), "p_sat_signed_range_percent": [100 * min(p), 100 * max(p)],
            "rho_liq_aad_percent": 100 * sum(rho) / len(rho)}


def main() -> None:
    document = json.loads(PACKET.read_text())
    sigma, epsilon = coefficient(document, MEA, "segment_diameter"), coefficient(document, MEA, "dispersion_energy_over_k")
    eps_w = coefficient(document, WATER, "dispersion_energy_over_k")
    kappa = next(e for e in document["topology"]["edges"] if e["volume"]["identity"] == EDGE)["volume"]["value"]["magnitude"]
    assert (sigma, epsilon, kappa) == (3.0435, 277.174, 0.03747)  # Baygi Table 2, 2B

    diameters = []
    for T in [293.15 + 10.0 * k for k in range(11)] + [403.15 + 10.0 * k for k in range(5)]:
        s_w = water_sigma(document, T)
        d_ij = (sigma * ratio(epsilon, T) + s_w * ratio(eps_w, T)) / (sigma + s_w)
        diameters.append({"T_K": T, "mea_mea_d_over_sigma_cubed": ratio(epsilon, T) ** 3,
                          "water_sigma_A": s_w, "mea_water_d_over_sigma_cubed": d_ij**3})

    engine = with_kappa(document, kappa)
    rows = []
    for T in GRID:
        cube = ratio(epsilon, T) ** 3
        p_s, rho_s = saturation(engine, T)
        p_d, rho_d = saturation(with_kappa(document, kappa * cube), T)
        rows.append({"T_K": T, "d_over_sigma_cubed": cube, "p_dippr_pa": dippr_psat(T), "rho_dippr_mol_l": dippr_rho(T),
                     "p_sigma3_pa": p_s, "rho_sigma3_mol_l": rho_s, "p_d3_pa": p_d, "rho_d3_mol_l": rho_d})

    converted_document = json.loads(CONVERTED.read_text())
    packet_kappa = next(e for e in converted_document["topology"]["edges"] if e["volume"]["identity"] == EDGE)["volume"]["value"]["magnitude"]
    converted = epcsaft.Pure(epcsaft.Parameters.from_mapping(converted_document, components=[MEA]))
    for row in rows:
        row["p_converted_pa"], row["rho_converted_mol_l"] = saturation(converted, row["T_K"])

    for name, table in (("mea-saturation.csv", rows), ("association-diameter-ratio.csv", diameters)):
        with (HERE / "results" / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(table)
    wheel = Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))
    summary = {
        "baygi_table_2_reported": {"p_sat_aad_percent": 0.62, "rho_liq_aad_percent": 0.12},
        "vs_correlations_baygi_range_303_443_k": {key: aad(within(rows, BAYGI), key) for key in ("d3", "sigma3", "converted")},
        "vs_correlations_packet_domain_293_393_k": {key: aad(within(rows, DOMAIN), key) for key in ("d3", "sigma3", "converted")},
        "vs_d3_model_packet_domain_293_393_k": {key: aad(within(rows, DOMAIN), key, "d3") for key in ("sigma3", "converted")},
        "conversion": {"packet_v5": packet_kappa, "ratio_to_baygi": packet_kappa / kappa,
                       "d_over_sigma_cubed_packet_domain": [within(rows, DOMAIN)[0]["d_over_sigma_cubed"], within(rows, DOMAIN)[-1]["d_over_sigma_cubed"]]},
        "identity": {"wheel_sha256": sha(wheel), "packet_v4": sha(PACKET), "packet_v5": sha(CONVERTED)},
        "claim_limits": "Pure MEA, Baygi Table 2 2B set, 5 K grid 293.15-443.15 K against Baygi's Table 1 correlations "
                        "(Avlund 2008 DIPPR), not raw measurements; Baygi's own point set is unstated. The conversion is a "
                        "one-coordinate fit to the d^3 model with m, sigma, epsilon fixed.",
    }
    (HERE / "results/mea-saturation.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
