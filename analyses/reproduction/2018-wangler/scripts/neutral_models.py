"""Neutral subsets and one-at-a-time readings of Wangler Tables 4, 6–7."""

import copy
import math
import runpy
from pathlib import Path
import epcsaft

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
COMMON = runpy.run_path(str(LEAF / "figures/figure-14/scripts/generate_data.py"))
a = COMMON["COMMON"]


def coeff(mapping, cid, family, value):
    component = next(c for c in mapping["components"] if c["component_id"] == cid)
    parameter = next(p for p in component["coefficients"] if p["family"] == family)
    parameter["value"]["magnitude"] = value


def neutral_mapping(t, mode="baseline", selected=("water", "mdea")):
    m = COMMON["parameter_mapping"](t, 0.0, a["LANES"]["S0"])
    m["components"] = [c for c in m["components"] if c["component_id"] in selected]
    m["pairs"] = [
        p
        for p in m["pairs"]
        if p["component_id_a"] in selected and p["component_id_b"] in selected
    ]
    m["correlations"] = [
        c
        for c in m["correlations"]
        if c["family"] == "segment_diameter" and c["component_id"] in selected
    ]
    m["topology"]["sites"] = [
        s for s in m["topology"]["sites"] if s["component_id"] in selected
    ]
    m["topology"]["edges"] = [
        q
        for q in m["topology"]["edges"]
        if q["endpoint_a"]["component_id"] in selected
        and q["endpoint_b"]["component_id"] in selected
    ]
    for f in m["model_families"]:
        if f["kind"] in ("electrolyte", "permittivity"):
            f["choice"] = "none"
            f.pop("ionic_dispersion", None)
    mutate(m, t, mode)
    return m


def mutate(m, t, mode):
    dt = t - 298.15
    sw = a["_water_sigma"](t)
    sm = 3.563
    ratio = math.sqrt(sw * sm) / ((sw + sm) / 2)
    association_correction = a["_interaction"](t, (0.6, -0.00257, -0.000017))
    cross = [
        q
        for q in m["topology"]["edges"]
        if {q["endpoint_a"]["component_id"], q["endpoint_b"]["component_id"]}
        == {"water", "mdea"}
    ]
    if mode == "eq12-literal":
        for q in cross:
            q["volume"]["value"]["magnitude"] = math.sqrt(
                0.045 * 0.1238 * ratio**3 * (1 - association_correction)
            )
    elif mode == "volume-geometric":
        for q in cross:
            q["volume"]["value"]["magnitude"] = math.sqrt(0.045 * 0.1238) * (
                1 - association_correction
            )
    elif mode == "volume-l-under-root":
        for q in cross:
            q["volume"]["value"]["magnitude"] = (
                math.sqrt(0.045 * 0.1238 * (1 - association_correction)) * ratio**3
            )
    elif mode == "energy-geometric-mayer":
        en = t * math.log1p(math.sqrt(math.expm1(2425.6 / t) * math.expm1(2046.6 / t)))
        for q in cross:
            q["energy_over_k"]["value"]["magnitude"] = en
    elif mode in ("literal-k", "raw-T-k", "k-constant", "quadratic-k", "cubic-sign-k"):
        kv = {
            "literal-k": -0.184 + (-0.00121 + 0.0000027 + 0.000000153) * dt,
            "raw-T-k": -0.184 - 0.00121 * t + 0.0000027 * t * t + 0.000000153 * t**3,
            "k-constant": -0.184,
            "quadratic-k": -0.184 - 0.00121 * dt + 0.0000027 * dt * dt,
            "cubic-sign-k": -0.184
            - 0.00121 * dt
            + 0.0000027 * dt * dt
            - 0.000000153 * dt**3,
        }[mode]
        for p in m["pairs"]:
            if {p["component_id_a"], p["component_id_b"]} == {"water", "mdea"}:
                p["coefficients"][0]["value"]["magnitude"] = kv
    elif mode == "literal-l":
        ln = 0.6 + (-0.00257 - 0.000017) * dt
        for q in cross:
            q["volume"]["value"]["magnitude"] *= (1 - ln) / (1 - association_correction)
    elif mode in ("mdea-1:1", "mdea-4:2", "mdea-2:3"):
        na, nb = {"mdea-1:1": (1, 1), "mdea-4:2": (2, 4), "mdea-2:3": (2, 3)}[mode]
        for s in m["topology"]["sites"]:
            if s["component_id"] == "mdea":
                s["multiplicity"] = na if s["site_id"] == "a" else nb
    elif mode in ("induced-retain-cross", "induced-zero-self-energy"):
        m["topology"]["edges"] = [
            q
            for q in m["topology"]["edges"]
            if not (
                q["endpoint_a"]["component_id"]
                == q["endpoint_b"]["component_id"]
                == "mdea"
            )
        ]
        if mode == "induced-zero-self-energy":
            for q in cross:
                q["energy_over_k"]["value"]["magnitude"] = 1212.8
    elif mode == "water-sigma-constant":
        m["correlations"] = [
            c
            for c in m["correlations"]
            if not (c["component_id"] == "water" and c["family"] == "segment_diameter")
        ]
        w = next(c for c in m["components"] if c["component_id"] == "water")
        v = copy.deepcopy(w["coefficients"][0])
        v.update(
            identity="component/water/segment_diameter",
            family="segment_diameter",
            value={"magnitude": 2.7927, "unit": "angstrom"},
        )
        w["coefficients"].append(v)
        rn = math.sqrt(2.7927 * sm) / ((2.7927 + sm) / 2)
        for q in cross:
            q["volume"]["value"]["magnitude"] *= rn**3 / ratio**3
    elif mode in ("water-pahlav-2", "water-pahlav-4", "water-gross-2"):
        vals = {
            "water-pahlav-2": (1.9599, 2.362, 279.42, 0.1750, 2059.28, 1),
            "water-pahlav-4": (2.1945, 2.229, 141.66, 0.2039, 1804.17, 2),
            "water-gross-2": (1.0656, 3.0007, 366.51, 0.034868, 2500.7, 1),
        }[mode]
        mm, ss, ee, kk, aa, nn = vals
        m["correlations"] = [
            c
            for c in m["correlations"]
            if not (c["component_id"] == "water" and c["family"] == "segment_diameter")
        ]
        coeff(m, "water", "segment_count", mm)
        coeff(m, "water", "dispersion_energy_over_k", ee)
        w = next(c for c in m["components"] if c["component_id"] == "water")
        v = copy.deepcopy(w["coefficients"][0])
        v.update(
            identity="component/water/segment_diameter",
            family="segment_diameter",
            value={"magnitude": ss, "unit": "angstrom"},
        )
        w["coefficients"].append(v)
        for q in m["topology"]["edges"]:
            ends = [q["endpoint_a"]["component_id"], q["endpoint_b"]["component_id"]]
            if ends == ["water", "water"]:
                q["volume"]["value"]["magnitude"] = kk
                q["energy_over_k"]["value"]["magnitude"] = aa
        for st in m["topology"]["sites"]:
            if st["component_id"] == "water":
                st["multiplicity"] = nn
        for q in cross:
            rn = math.sqrt(ss * sm) / ((ss + sm) / 2)
            q["volume"]["value"]["magnitude"] = (
                math.sqrt(kk * 0.1238) * rn**3 * (1 - association_correction)
            )
            q["energy_over_k"]["value"]["magnitude"] = (aa + 2046.6) / 2
    elif mode in ("mdea-source-6site", "mdea-source-2site"):
        mm, ss, ee, kk, aa, nd, na = {
            "mdea-source-6site": (4.874, 3.211, 288.87, 0.016894, 513.2765, 2, 4),
            "mdea-source-2site": (3.3044, 3.5975, 237.44, 0.066454, 3709.9, 1, 1),
        }[mode]
        coeff(m, "mdea", "segment_count", mm)
        coeff(m, "mdea", "segment_diameter", ss)
        coeff(m, "mdea", "dispersion_energy_over_k", ee)
        for st in m["topology"]["sites"]:
            if st["component_id"] == "mdea":
                st["multiplicity"] = nd if st["site_id"] == "a" else na
        for q in m["topology"]["edges"]:
            ends = [q["endpoint_a"]["component_id"], q["endpoint_b"]["component_id"]]
            if ends == ["mdea", "mdea"]:
                q["volume"]["value"]["magnitude"] = kk
                q["energy_over_k"]["value"]["magnitude"] = aa
        for q in cross:
            rn = math.sqrt(sw * ss) / ((sw + ss) / 2)
            q["volume"]["value"]["magnitude"] = (
                math.sqrt(0.045 * kk) * rn**3 * (1 - association_correction)
            )
            q["energy_over_k"]["value"]["magnitude"] = (2425.6 + aa) / 2
    elif mode == "mdea-origin-digits":
        coeff(m, "mdea", "dispersion_energy_over_k", 228.711)
        for q in m["topology"]["edges"]:
            ends = [q["endpoint_a"]["component_id"], q["endpoint_b"]["component_id"]]
            if ends == ["mdea", "mdea"]:
                q["volume"]["value"]["magnitude"] = 0.123858
                q["energy_over_k"]["value"]["magnitude"] = 2046.624
        for q in cross:
            q["volume"]["value"]["magnitude"] *= math.sqrt(0.123858 / 0.1238)
            q["energy_over_k"]["value"]["magnitude"] = (2425.6 + 2046.624) / 2
    elif mode == "association-diameter-d":
        sig = {"water": sw, "mdea": sm}
        eps = {"water": 353.94, "mdea": 228.71}
        diam = {c: s * (1 - 0.12 * math.exp(-3 * eps[c] / t)) for c, s in sig.items()}
        for q in m["topology"]["edges"]:
            i, j = q["endpoint_a"]["component_id"], q["endpoint_b"]["component_id"]
            q["volume"]["value"]["magnitude"] *= (
                (diam[i] + diam[j]) / (sig[i] + sig[j])
            ) ** 3
    elif mode == "baseline":
        pass
    else:
        raise ValueError(mode)
    return m


def declarations():
    return [
        (epcsaft.Mixture(epcsaft.Parameters.from_mapping(neutral_mapping(353.0))), None)
    ]
