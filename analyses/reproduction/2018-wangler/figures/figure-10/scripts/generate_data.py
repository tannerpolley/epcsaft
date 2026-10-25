"""Wangler Tables 8–9 (Figure 10): native fits of the water–MDEAH⁺ k_ij, one per isotherm (#330).

Wangler's OF3 (Eq. 33) on the ten osmotic coefficients at 273 and 313 K, with l_hb held at the
printed Table 7 values. Two lanes differ only in ion self-dispersion: S0 excludes it (Held 2014,
Uyan 2015; best-supported reading), S1 keeps it (literal l.146; all-ion-pairs is exact for 1:1).
Both use packing diameter 0.88σ (Eq. 7), Debye–Hückel diameter σ and MDEAH⁺ m = 3.675.
The #234 batch-5 comparison is retired; its batch-5-* files are frozen records.
"""

from __future__ import annotations

import csv
import json
import hashlib
from importlib import metadata
import math
from pathlib import Path
from urllib.parse import unquote, urlparse

import numpy as np
from epcsaft import Mixture, Parameters, regression

CAMPAIGN = Path(__file__).resolve().parents[3]
SOURCE = CAMPAIGN / "inputs/source" / "wangler-2018-tables-8-9.csv"
OUTPUT = Path(__file__).resolve().parents[1] / "output"
PROVENANCE = {
    "source_id": "wangler-2018",
    "locator": "Tables 4-9 and Eqs. 8-12, 33-35",
    "domain_id": "wangler-mdeah-chloride-osmotic",
}
UNITS = {
    "molar_mass": "kilogram / mole",
    "charge_number": "elementary-charge",
    "segment_count": "dimensionless",
    "segment_diameter": "angstrom",
    "dispersion_energy_over_k": "kelvin",
    "packing_diameter": "angstrom",
    "debye_huckel_diameter": "angstrom",
    "k_ij": "dimensionless",
}


def _value(magnitude: float, family: str) -> dict[str, object]:
    return {"magnitude": magnitude, "unit": UNITS[family]}


def _component(component_id: str, **families: float) -> dict[str, object]:
    fixed = {}
    coefficients = []
    for family, magnitude in families.items():
        record = {"value": _value(magnitude, family), "provenance": PROVENANCE}
        if family in {"molar_mass", "charge_number"}:
            fixed[family] = record
        else:
            coefficients.append(
                {
                    "identity": f"component/{component_id}/{family}",
                    "family": family,
                    **record,
                }
            )
    return {"component_id": component_id, "fixed": fixed, "coefficients": coefficients}


def _pair(component_a: str, component_b: str, kij: float) -> dict[str, object]:
    component_a, component_b = sorted((component_a, component_b))
    return {
        "component_id_a": component_a,
        "component_id_b": component_b,
        "coefficients": [
            {
                "identity": f"pair/{component_a}/{component_b}/k_ij",
                "family": "k_ij",
                "value": _value(kij, "k_ij"),
                "provenance": PROVENANCE,
            }
        ],
    }


def _water_sigma(temperature: float) -> float:
    return (
        2.7927
        + 10.11 * math.exp(-0.01775 * temperature)
        - 1.417 * math.exp(-0.01146 * temperature)
    )


def _water_permittivity(temperature: float) -> float:
    return (
        7.6555618295e-4 * temperature**2 - 8.1783881423e-1 * temperature + 254.19616803
    )


def _interaction(temperature: float, coefficients: tuple[float, ...]) -> float:
    delta = temperature - 298.15
    return math.fsum(value * delta**power for power, value in enumerate(coefficients))


def _edge(
    endpoint_a: tuple[str, str],
    endpoint_b: tuple[str, str],
    energy: float,
    volume: float,
) -> dict[str, object]:
    endpoint_a, endpoint_b = sorted((endpoint_a, endpoint_b))
    prefix = (
        f"association/{endpoint_a[0]}/{endpoint_a[1]}/{endpoint_b[0]}/{endpoint_b[1]}"
    )
    return {
        "endpoint_a": {"component_id": endpoint_a[0], "site_id": endpoint_a[1]},
        "endpoint_b": {"component_id": endpoint_b[0], "site_id": endpoint_b[1]},
        "energy_over_k": {
            "identity": f"{prefix}/energy_over_k",
            "value": {"magnitude": energy, "unit": "kelvin"},
        },
        "volume": {
            "identity": f"{prefix}/volume",
            "value": {"magnitude": volume, "unit": "dimensionless"},
        },
        "source": {"kind": "explicit", "provenance": [PROVENANCE]},
    }


def _mixing_factor(temperature: float) -> float:
    # Wangler Eqs. 33–35: freeze the published diameter mixing factor at each isotherm.
    sigma = _water_sigma(temperature)
    return math.sqrt(0.0450 * 0.1238) * (
        math.sqrt(sigma * 3.5630) / (0.5 * (sigma + 3.5630))
    ) ** 3


def _topology(temperature: float, l_hb: float) -> dict[str, object]:
    sites = [
        {
            "component_id": component,
            "site_id": site,
            "site_role": role,
            "multiplicity": multiplicity,
            "provenance": PROVENANCE,
        }
        for component, multiplicity in (("water", 1), ("mdeah-cation", 2))
        for site, role in (("a", "donor"), ("b", "acceptor"))
    ]
    water_energy = 2425.6
    water_volume = 0.0450
    cross_energy = 0.5 * water_energy
    cross_volume = _mixing_factor(temperature) * (1.0 - l_hb)
    return {
        "presets": [],
        "sites": sites,
        "edges": [
            _edge(("water", "a"), ("water", "b"), water_energy, water_volume),
            _edge(("water", "a"), ("mdeah-cation", "b"), cross_energy, cross_volume),
            _edge(("water", "b"), ("mdeah-cation", "a"), cross_energy, cross_volume),
        ],
    }


def _model(temperature: float, kij: float, l_hb: float, ionic_dispersion: str) -> Mixture:
    diameter_correlation = {
        "correlation_id": "component/water/segment_diameter/constant-plus-sum-of-exponentials",
        "component_id": "water",
        "family": "segment_diameter",
        "form": "constant-plus-sum-of-exponentials",
        "independent_variables": ["temperature"],
        "constant": {
            "identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/constant",
            "value": {"magnitude": 2.7927, "unit": "angstrom"},
        },
        "terms": [
            {
                "amplitude": {
                    "identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude",
                    "value": {"magnitude": 10.11, "unit": "angstrom"},
                },
                "exponent_coefficient": {
                    "identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient",
                    "value": {"magnitude": -0.01775, "unit": "1 / kelvin"},
                },
            },
            {
                "amplitude": {
                    "identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude",
                    "value": {"magnitude": -1.417, "unit": "angstrom"},
                },
                "exponent_coefficient": {
                    "identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient",
                    "value": {"magnitude": -0.01146, "unit": "1 / kelvin"},
                },
            },
        ],
        "provenance": PROVENANCE,
    }
    permittivity_correlation = {
        "correlation_id": "component/water/relative_permittivity/constant",
        "component_id": "water",
        "family": "relative_permittivity",
        "form": "constant",
        "independent_variables": [],
        "constant": {
            "identity": "component/water/relative_permittivity/constant/constant",
            "value": {
                "magnitude": _water_permittivity(temperature),
                "unit": "dimensionless",
            },
        },
        "provenance": PROVENANCE,
    }
    mapping = {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "document_id": "wangler-2018-mdeah-chloride",
        "document_version": 1,
        "purpose": "external-document",
        "sources": [
            {
                "source_id": "wangler-2018",
                "citation": "Wangler et al., Fluid Phase Equilibria 461 (2018) 15-27.",
                "use_basis": "Tables 4-9 and Eqs. 8-12, 33-35.",
            }
        ],
        "domains": [
            {
                "domain_id": "wangler-mdeah-chloride-osmotic",
                "kind": "reported-conditions",
                "temperature_min": {"magnitude": 273.0, "unit": "kelvin"},
                "temperature_max": {"magnitude": 313.0, "unit": "kelvin"},
                "pressure_min": {"magnitude": 1.0, "unit": "bar"},
                "pressure_max": {"magnitude": 1.0, "unit": "bar"},
            }
        ],
        "components": [
            _component(
                "water",
                molar_mass=0.01801528,
                charge_number=0,
                segment_count=1.2046,
                dispersion_energy_over_k=353.94,
            ),
            _component(
                "mdeah-cation",
                charge_number=1,
                # Wangler treats MDEAH+ as a nonspherical ion and transfers
                # the molecular MDEA segment count as well as sigma and u/k.
                segment_count=3.6750,
                segment_diameter=3.5630,
                dispersion_energy_over_k=228.71,
                packing_diameter=0.88 * 3.5630,
                debye_huckel_diameter=3.5630,
            ),
            _component(
                "chloride-anion",
                charge_number=-1,
                segment_count=1.0,
                segment_diameter=2.7560,
                dispersion_energy_over_k=170.0,
                packing_diameter=0.88 * 2.7560,
                debye_huckel_diameter=2.7560,
            ),
        ],
        "pairs": [
            _pair("water", "mdeah-cation", kij),
            _pair("water", "chloride-anion", -0.25),
            _pair("mdeah-cation", "chloride-anion", 0.0),
        ],
        "model_families": [
            {
                "family_id": "model/association",
                "kind": "association",
                "choice": "general-site",
                "provenance": PROVENANCE,
            },
            {
                "family_id": "model/base",
                "kind": "base",
                "choice": "pc-saft",
                "provenance": PROVENANCE,
            },
            {
                "family_id": "model/electrolyte",
                "kind": "electrolyte",
                "choice": "fully-dissociated-debye-huckel",
                "ionic_dispersion": ionic_dispersion,
                "provenance": PROVENANCE,
            },
            {
                "family_id": "model/relative_permittivity",
                "kind": "permittivity",
                "choice": "solvent-only",
                "provenance": PROVENANCE,
            },
        ],
        "model_coefficients": [],
        "correlations": [diameter_correlation, permittivity_correlation],
        "topology": _topology(temperature, l_hb),
    }
    return Mixture(Parameters.from_mapping(mapping))


def _read_rows() -> list[dict[str, float]]:
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]


def _observations(parameters: Parameters, rows: list[dict[str, float]]) -> list:
    return [regression.electrolyte(
        parameters, "osmotic_coefficient", row["temperature_k"], 1e5, {"water": 1.0},
        ions=("mdeah-cation", "chloride-anion"), stoichiometry=(1, 1),
        molality=row["molality_mol_per_kg"], observed=row["osmotic_coefficient"],
        form="difference", scale=row["osmotic_coefficient"],
    ) for row in rows]


def _prediction(model: Mixture, observation) -> float:
    prediction = regression.predict(model, observation)
    if prediction.status != regression.ObservationStatus.Available:
        raise RuntimeError(prediction.message)
    return prediction.value


LANES = {"S0": "exclude-same-sign-ion-pairs", "S1": "all-ion-pairs"}
PRINTED_K_IJ = (-0.184, -1.21e-3, 2.70e-6, 1.53e-7)  # Table 6 (a copy of water–MDEA)
L_HB = (0.372, 1.07e-3)  # Table 7, Eq. 35
U_PHI = 0.005  # Tables 8–9 standard uncertainty
ARD_LIMIT = 6.05  # F1: the paper's 5.05% plus one point


def _fit_isotherm(ionic_dispersion: str, temperature: float, rows: list[dict[str, float]]) -> dict:
    l_hb = _interaction(temperature, L_HB)
    parameters = _model(temperature, 0.0, l_hb, ionic_dispersion).parameters
    observations = _observations(parameters, rows)
    k_ij = regression.coordinate(parameters, [("k_ij", ("water", "mdeah-cation"), 1)],
                                 origin=0.0, scale=1.0, bounds=(-1.0, 1.0))
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=120.0,
                                      function_tolerance=1e-10, parameter_tolerance=1e-10, gradient_tolerance=1e-10)
    fit, other = regression.fit(parameters, [k_ij], observations, weights=[1.0] * len(observations),
                                controls=controls, holds=(), starts=((0.0,), (0.5,)),
                                shared_iterations=100, shared_seconds=120.0, singular_vectors=True)
    if not (fit.usable and other.usable):
        raise RuntimeError(f"{ionic_dispersion} {temperature} K: {fit.message} / {other.message}")
    # #245 endpoint judgment: |dk| <= sum_i |J_w^+|_i u(phi)/phi_obs,i, with r_i = phi_i/phi_obs,i - 1.
    jacobian = np.array(fit.optimizer_jacobian).reshape(len(observations), 1)
    allowance = float((np.abs(np.linalg.pinv(jacobian)) @ np.array([U_PHI / o.scale for o in observations]))[0])
    printed = _model(temperature, _interaction(temperature, PRINTED_K_IJ), l_hb, ionic_dispersion)
    cov = fit.covariance
    return {
        "l_hb_table_7": l_hb, "k_ij": fit.physical[0], "other_start_k_ij": other.physical[0],
        "start_difference": abs(fit.physical[0] - other.physical[0]), "u_phi_allowance": allowance,
        "status": fit.status, "message": fit.message, "iterations": fit.iterations,
        "residual_evaluations": fit.residual_evaluations, "relative_sse": 2.0 * fit.final_cost,
        "covariance": {"status": str(cov.status), "rank": cov.rank, "singular_values": list(cov.singular_values),
                       "active_bounds": list(cov.active_bounds), "physical": list(cov.physical)},
        "fitted_phi": list(fit.predictions),
        "printed_parameter_phi": [_prediction(printed, o) for o in observations],
    }


def main() -> None:
    rows = _read_rows()
    temperatures = sorted({row["temperature_k"] for row in rows})
    plotted, lanes = [], {}
    for lane, ionic_dispersion in LANES.items():
        isotherms = {}
        for temperature in temperatures:
            selected = [row for row in rows if row["temperature_k"] == temperature]
            result = _fit_isotherm(ionic_dispersion, temperature, selected)
            isotherms[f"{temperature:.0f}"] = result
            plotted += [[lane, temperature, row["molality_mol_per_kg"], row["osmotic_coefficient"], printed, fitted]
                        for row, printed, fitted in zip(selected, result.pop("printed_parameter_phi"),
                                                        result.pop("fitted_phi"), strict=True)]
            print(lane, temperature, json.dumps(result), flush=True)
        lane_rows = [r for r in plotted if r[0] == lane]
        ard = 100.0 * float(np.mean([abs(r[5] / r[3] - 1.0) for r in lane_rows]))
        low, high = isotherms["273"], isotherms["313"]
        # Eq. 34 truncated to two terms is linear in T: k(303 K) = 0.25 k(273 K) + 0.75 k(313 K).
        interior = all(i["covariance"]["status"] == "available" and not i["covariance"]["active_bounds"]
                       for i in (low, high))
        starts_agree = all(i["start_difference"] <= i["u_phi_allowance"] / 100.0 for i in (low, high))
        lanes[lane] = {
            "ionic_dispersion": ionic_dispersion, "isotherms": isotherms, "ard_percent": ard,
            "printed_parameter_ard_percent": 100.0 * float(np.mean([abs(r[4] / r[3] - 1.0) for r in lane_rows])),
            "k_ij_303_k": 0.25 * low["k_ij"] + 0.75 * high["k_ij"],
            "k_ij_303_k_allowance": 0.25 * low["u_phi_allowance"] + 0.75 * high["u_phi_allowance"],
            "F1": {"ard_limit_percent": ARD_LIMIT, "interior": interior, "starts_agree": starts_agree,
                   "passes": interior and starts_agree and ard <= ARD_LIMIT},
        }
    with (OUTPUT / "wangler-2018-mdeah-osmotic-plotted-data.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["lane", "temperature_k", "molality_mol_per_kg", "experimental_osmotic_coefficient",
                         "printed_parameter_model", "fitted_model"])
        writer.writerows(plotted)
    wheel_url = json.loads(metadata.distribution("epcsaft").read_text("direct_url.json"))
    if wheel_url.get("dir_info", {}).get("editable", False):
        raise RuntimeError("non-editable Engine wheel required")
    wheel = Path(unquote(urlparse(wheel_url["url"]).path))
    summary = {
        "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "paper_ard_percent": 5.05,
        "objective": "OF3 (Eq. 33): sum (phi/phi_obs - 1)^2, native Ceres Levenberg-Marquardt, exact Jacobian",
        "coordinate": "water-MDEAH+ k_ij per isotherm, bounds [-1, 1], declared starts 0.0 and 0.5",
        "held": "l_hb Table 7 (Eq. 35); MDEAH+ m 3.675, sigma 3.5630, u/k 228.71; packing 0.88 sigma; "
                "Debye-Hueckel diameter sigma; water permittivity; MDEAH+-Cl- k_ij 0; water-Cl- k_ij -0.25",
        "u_phi": U_PHI,
        "fit_controls": {"function_tolerance": 1e-10, "parameter_tolerance": 1e-10, "gradient_tolerance": 1e-10,
                         "maximum_iterations": 100, "shared_iterations": 100, "shared_seconds": 120.0},
        "lanes": lanes,
        "limits": "Calibration to ten Tables 8-9 osmotic coefficients at 1 bar, 273 and 313 K, molality <= 1.14 "
                  "mol/kg. Two isotherms determine one k_ij level each; no k_ij(T) form beyond linear and no "
                  "parameter adoption follows.",
    }
    (OUTPUT / "wangler-2018-mdeah-osmotic-summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def declarations():
    return [
        (
            _model(t, _interaction(t, PRINTED_K_IJ), _interaction(t, L_HB), d),
            _observations(
                _model(t, 0.0, _interaction(t, L_HB), d).parameters,
                [r for r in rows if r['temperature_k'] == t],
            ),
        )
        for rows in [_read_rows()]
        for d in LANES.values()
        for t in sorted({r['temperature_k'] for r in rows})
    ]


if __name__ == "__main__":
    main()
