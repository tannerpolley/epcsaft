"""Bounded declarations; all state and equilibrium calculations use the installed native Engine."""

import copy
import csv
import hashlib
import importlib.util
import json
import math
import sys
import subprocess
from pathlib import Path
import epcsaft
from epcsaft import equilibrium as eq

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[1]
OUT = LEAF / "results/pressure-diagnosis"
OUT.mkdir(parents=True, exist_ok=True)
DOC = json.loads((LEAF / "inputs/parameters.json").read_text())
T = 323.15
SPECIES = ["water", "acetic-acid"]
spec = importlib.util.spec_from_file_location("acetic_case", LEAF / "scripts/run.py")
case = importlib.util.module_from_spec(spec)
spec.loader.exec_module(case)
rows = []


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save():
    (OUT / "probe-results.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False) + "\n"
    )


def state(st):
    return {
        k: list(getattr(st, k))
        if k in ("log_fugacity_coefficient", "site_fractions")
        else getattr(st, k)
        for k in (
            "pressure",
            "molar_density",
            "compressibility_factor",
            "residual_helmholtz",
            "association",
            "dispersion",
            "hard_chain",
            "log_fugacity_coefficient",
            "site_fractions",
        )
    }


def mixture(doc=DOC):
    return epcsaft.Mixture(epcsaft.Parameters.from_mapping(doc, components=SPECIES))


def record(m, q, label, xw=None):
    r = eq.solve_equilibrium(m, q)
    n = len(m.component_ids)
    d = dict(
        probe=label,
        x_water=xw,
        success=r.success,
        requested_tolerance_met=r.requested_tolerance_met,
        termination=str(r.solver_status),
        row_residual=r.convergence.residual,
        iterations=r.iterations,
        P_Pa=r.pressure,
        densities_mol_m3=list(r.molar_densities),
        mole_fractions=list(r.mole_fractions),
        amounts_mol=list(r.amounts),
    )
    rows.append(d)
    save()
    if not r.success or not r.requested_tolerance_met:
        raise RuntimeError(d)
    d["states"] = [
        state(
            m.state(
                T,
                rho=r.molar_densities[j],
                x=list(r.mole_fractions[j * n : (j + 1) * n]),
            )
        )
        for j in (0, 1)
    ]
    d["relative_pressure_closure"] = max(
        abs(s["pressure"] / r.pressure - 1) for s in d["states"]
    )
    d["log_fugacity_mismatch"] = [
        math.log(r.mole_fractions[i] / r.mole_fractions[n + i])
        + d["states"][0]["log_fugacity_coefficient"][i]
        - d["states"][1]["log_fugacity_coefficient"][i]
        for i in range(n)
        if r.mole_fractions[n + i] > 0
    ]
    assert d["relative_pressure_closure"] < 1e-7
    assert max(map(abs, d["log_fugacity_mismatch"])) < 1e-7
    print(label, xw, r.pressure, "Pa", r.convergence.residual, flush=True)
    save()
    return r, d


def pure(doc, name, label):
    m = epcsaft.Pure(epcsaft.Parameters.from_mapping(doc, components=[name]))
    return record(m, m.saturation(T), label)


def bubble(doc, label, xw):
    m = mixture(doc)
    return record(m, case.problem(xw, "no-dissociation"), label, xw)


def cross_doc(rule=None, correction=None):
    d = copy.deepcopy(DOC)
    for e in d["topology"]["edges"][2:]:
        if rule:
            e["source"]["rule_id"] = rule
        if correction is not None:
            e["volume_correction"] = {
                "identity": e["endpoint_a"]["component_id"]
                + "/"
                + e["endpoint_a"]["site_id"]
                + "/cross-volume-probe",
                "value": {"magnitude": correction, "unit": "dimensionless"},
            }
    return d


identity = {
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip(),
    "engine_source_tree": subprocess.check_output(
        ["git", "rev-parse", "HEAD:engine"], cwd=ROOT, text=True
    ).strip(),
    "python": sys.executable,
    "package": epcsaft.__file__,
    "extension": epcsaft._core.__file__,
    "input_sha256": {
        str(p.relative_to(LEAF)): digest(p)
        for p in (LEAF / "inputs").glob("*")
        if p.is_file()
    },
    "probe_sha256": digest(__file__),
    "temperature_K": T,
}
assert ROOT in Path(epcsaft._core.__file__).resolve().parents, (
    "Use this checkout's prepared .venv"
)

for name in SPECIES:
    pure(DOC, name, "pure-baseline-" + name)
raw = json.loads((LEAF / "inputs/bernatova-2006.json").read_text())
groups = {g["nPureOrMixtureDataNumber"]: g for g in raw["PureOrMixtureData"]}
by_y = {
    r["VariableValue"][0]["nVarValue"]: r["PropertyValue"][0]
    for r in groups[13]["NumValues"]
}
measurements = []
for i, r in enumerate(groups[14]["NumValues"], 1):
    xw = r["PropertyValue"][0]["nPropValue"]
    y = r["VariableValue"][0]["nVarValue"]
    p = by_y[y]
    measurements.append(
        dict(
            row=i,
            x_water=xw,
            y_water_reported=y,
            P_Pa=1000 * p["nPropValue"],
            U95_P_Pa=1000 * p["CombinedUncertainty"]["nCombExpandUncertValue"],
            U95_x=r["PropertyValue"][0]["CombinedUncertainty"][
                "nCombExpandUncertValue"
            ],
        )
    )
    bubble(DOC, "neutral-baseline", xw)
with (LEAF / "inputs/measurements.csv").open() as f:
    retained = list(csv.DictReader(f))
assert len(retained) == len(measurements)
mismatches = [
    {"row": m["row"], "quantity": key, "raw": m[key], "csv": float(r[key])}
    for m, r in zip(measurements, retained)
    for key in m
    if m[key] != float(r[key])
]
assert not mismatches, mismatches
(OUT / "measurement-check.json").write_text(
    json.dumps(
        {
            "temperature_K": T,
            "raw_sha256": digest(LEAF / "inputs/bernatova-2006.json"),
            "measurements": measurements,
            "retained_measurement_csv_check": {
                "rows": len(retained),
                "mismatches": mismatches,
                "gas_fraction_used_only_as_join_key": True,
                "pressure_unit": "kPa in group 13; multiplied by 1000 to Pa",
            },
            "pure_endpoints": [
                {
                    "reported_y_water": y,
                    "P_Pa": 1000 * by_y[y]["nPropValue"],
                    "U95_P_Pa": 1000
                    * by_y[y]["CombinedUncertainty"]["nCombExpandUncertValue"],
                }
                for y in (0, 1)
            ],
        },
        indent=2,
    )
    + "\n"
)

for xw in (0.5781, 0.6673, 0.7683):
    record(case.model("primary"), case.problem(xw, "primary"), "reactive-baseline", xw)
    for delta in (-0.01, 0.01):
        d = copy.deepcopy(DOC)
        d["pairs"][0]["coefficients"][0]["value"]["magnitude"] += delta
        bubble(d, f"kij-delta-{delta:+.2f}", xw)
    bubble(
        cross_doc(rule="arithmetic-energy-geometric-volume"), "cross-no-size-factor", xw
    )
    # Equation A7 as printed in Ascani 2023 has a square, compared with Gross 2002's cube.
    # This is parameter-input algebra at this one temperature, not a thermodynamic solve.
    sigw = 2.7927 + 10.11 * math.exp(-0.01775 * T) - 1.417 * math.exp(-0.01146 * T)
    siga = 3.8582
    size_ratio = math.sqrt(sigw * siga) / (0.5 * (sigw + siga))
    bubble(cross_doc(correction=1 - 1 / size_ratio), "cross-size-factor-square", xw)
    for correction in (-0.01, 0.01):
        bubble(
            cross_doc(correction=correction),
            f"cross-volume-correction-{correction:+.2f}",
            xw,
        )
    # Native same-state controls expose vapour association without replacing a phase model.
    base = next(
        r for r in rows if r["probe"] == "neutral-baseline" and r["x_water"] == xw
    )
    vap = base["states"][1]
    y = base["mole_fractions"][2:]
    for control in (
        "explicit-baseline",
        "acid-self-off",
        "cross-off",
        "all-association-off",
    ):
        d = copy.deepcopy(DOC)
        # Freeze the original cross strengths so a self-edge removal does not alter them.
        for edge in d["topology"]["edges"][2:]:
            edge["source"] = {
                "kind": "explicit",
                "provenance": [{"probe": "baseline cross strength at 323.15 K"}],
            }
            edge["energy_over_k"] = {
                "identity": edge["endpoint_a"]["site_id"] + "/cross-energy",
                "value": {"magnitude": 2735.05, "unit": "kelvin"},
            }
            edge["volume"] = {
                "identity": edge["endpoint_a"]["site_id"] + "/cross-volume",
                "value": {
                    "magnitude": math.sqrt(0.0451 * 0.07555) * size_ratio**3,
                    "unit": "dimensionless",
                },
            }
        if control == "acid-self-off":
            d["topology"]["edges"][1]["energy_over_k"]["value"]["magnitude"] = 0.0
        elif control == "cross-off":
            d["topology"]["edges"] = d["topology"]["edges"][:2]
            d["topology"]["no_cross_association"] = [
                {
                    "component_id_a": "water",
                    "component_id_b": "acetic-acid",
                    "provenance": {"probe": "remove-cross-edges-only"},
                }
            ]
        elif control == "all-association-off":
            for e in d["topology"]["edges"]:
                e["energy_over_k"]["value"]["magnitude"] = 0.0
        m = mixture(d)
        if control == "explicit-baseline":
            check = state(m.state(T, rho=vap["molar_density"], x=y))
            assert abs(check["pressure"] - base["P_Pa"]) < 1e-6
            assert (
                max(
                    abs(a - b)
                    for a, b in zip(
                        check["log_fugacity_coefficient"],
                        vap["log_fugacity_coefficient"],
                    )
                )
                < 1e-9
            )
        rows.append(
            dict(
                probe="vapor-" + control,
                x_water=xw,
                density_fixed=True,
                state=state(m.state(T, rho=vap["molar_density"], x=y)),
            )
        )
        save()

# Source-rounding and inconsistent source-table controls, not parameter fitting.
d = copy.deepcopy(DOC)
d["topology"]["edges"][0]["volume"]["value"]["magnitude"] = 0.04509
pure(d, "water", "pure-water-volume-0.04509")
for xw in (0.5781, 0.6673, 0.7683):
    bubble(d, "water-volume-0.04509", xw)
d = copy.deepcopy(DOC)
for c in d["components"]:
    if c["component_id"] == "acetic-acid":
        for p in c["coefficients"]:
            if p["family"] == "dispersion_energy_over_k":
                p["value"]["magnitude"] = 311.59
pure(d, "acetic-acid", "pure-acid-Ascani2023-printed-energy-311.59")
print("completed", len(rows), "retained results", flush=True)

paths = sorted(
    {
        line.split()[-1]
        for line in Path("/proc/self/maps").read_text().splitlines()
        if "epcsaft" in line and ".so" in line
    }
)
identity["loaded_binary_sha256"] = {p: digest(p) for p in paths}
(OUT / "run-identity.json").write_text(json.dumps(identity, indent=2) + "\n")
