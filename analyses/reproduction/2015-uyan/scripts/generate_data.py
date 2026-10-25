"""Uyan source-ledger reference inputs and native ordered pressure starts."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import hashlib
import math
import runpy
import time
from pathlib import Path
from typing import cast

import epcsaft
from epcsaft import equilibrium, regression

ROOT = Path(__file__).resolve().parents[4]
CAMPAIGN = Path(__file__).resolve().parents[1]
SOURCE = CAMPAIGN / "inputs/source" / "ermatchkov-2006-tables-1-3.csv"
OUTPUT = CAMPAIGN / "results" / "data" / "native" / "source-ledger"
COMMON = runpy.run_path(
    str(
        ROOT
        / "analyses/reproduction/2018-wangler/figures/figure-10/scripts/generate_data.py"
    )
)
component = COMMON["_component"]
pair = COMMON["_pair"]
edge = COMMON["_edge"]
water_sigma = COMMON["_water_sigma"]
archer_wang_water_permittivity = COMMON["_water_permittivity"]
interaction = COMMON["_interaction"]

PROVENANCE = {
    "source_id": "uyan-2015",
    "locator": "Uyan Tables 1-5 and Eqs. 5, 12-27; Ermatchkov Tables 1-3",
    "domain_id": "uyan-2015-ermatchkov-low-loading-mdea-co2",
}
for helper in (component, pair, edge):
    helper.__globals__["PROVENANCE"] = PROVENANCE

COMPONENTS = (
    "water",
    "mdea",
    "co2",
    "h-cation",
    "oh-anion",
    "hco3-anion",
    "carbonate-anion",
    "mdeah-cation",
)
CHARGES = (0, 0, 0, 1, -1, -1, -2, 1)
MOLAR_MASSES = (
    0.01801528,
    0.119163,
    0.0440095,
    0.00100794,
    0.01700734,
    0.06101684,
    0.0600089,
    0.12017094,
)
REACTIONS = (
    (-1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0),
    (-1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 0.0, 0.0),
    (0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
)
LN_K_COEFFICIENTS = (
    (132.899, -13445.9, -22.477, 0.0),
    (212.739, -11333.8, -33.844, -0.0018149),
    (287.444, -13648.7, -48.880, 0.030232),
    (-83.491, -819.7, 10.976, 0.0),
)


def temperature_polynomial_correlation(
    component_id: str,
    coefficients: tuple[float, ...],
    *,
    reference_temperature: float,
    provenance: dict[str, str],
) -> dict[str, object]:
    correlation_id = (
        f"component/{component_id}/relative_permittivity/"
        "constant-plus-polynomial-in-temperature"
    )
    return {
        "correlation_id": correlation_id,
        "component_id": component_id,
        "family": "relative_permittivity",
        "form": "constant-plus-polynomial-in-temperature",
        "independent_variables": ["temperature"],
        "constant": {
            "identity": f"{correlation_id}/constant",
            "value": {"magnitude": coefficients[0], "unit": "dimensionless"},
        },
        "reference_temperature": {
            "magnitude": reference_temperature,
            "unit": "kelvin",
        },
        "terms": [
            {
                "power": power,
                "coefficient": {
                    "identity": f"{correlation_id}/term-{power - 1}/coefficient",
                    "value": {
                        "magnitude": coefficient,
                        "unit": "1 / kelvin"
                        if power == 1
                        else f"1 / kelvin ** {power}",
                    },
                },
            }
            for power, coefficient in enumerate(coefficients[1:], start=1)
        ],
        "provenance": provenance,
    }


def mdea_permittivity(temperature: float) -> float:
    # Least-squares representation of the six pure-MDEA Hsieh et al. 2007
    # ThermoML rows at 303.15-328.15 K. Uyan also extrapolated that source
    # correlation above its measurement range.
    return 52.905597142857 - 0.102457142857 * temperature


def water_permittivity_transfer_approximation(temperature: float) -> float:
    """Return the explicit Archer--Wang approximation used in this replay.

    Uyan cites the pressure-dependent Floriano--Nascimento source without
    publishing the evaluated correlation or values used by the authors. The
    approximation is therefore named and provenance-labeled rather than
    silently presented as the cited source implementation.
    """
    return archer_wang_water_permittivity(temperature)


def cross_volume(
    left_volume: float,
    right_volume: float,
    left_sigma: float,
    right_sigma: float,
    correction: float,
) -> float:
    sigma_cross = 0.5 * (left_sigma + right_sigma)
    return (
        math.sqrt(left_volume * right_volume)
        * (math.sqrt(left_sigma * right_sigma) / sigma_cross) ** 3
        * (1.0 - correction)
    )


def association_topology(temperature: float) -> dict[str, object]:
    sites = [
        {
            "component_id": component_id,
            "site_id": site,
            "site_role": role,
            "multiplicity": multiplicity,
            "provenance": PROVENANCE,
        }
        for component_id, multiplicity in (
            ("water", 1),
            ("mdea", 2),
            ("mdeah-cation", 2),
        )
        for site, role in (("a", "donor"), ("b", "acceptor"))
    ]
    sigma_water = water_sigma(temperature)
    sigma_mdea = 3.563
    water_energy, water_volume = 2425.67, 0.04509
    mdea_energy, mdea_volume = 2046.624, 0.123858
    edges = [
        edge(("water", "a"), ("water", "b"), water_energy, water_volume),
        edge(("mdea", "a"), ("mdea", "b"), mdea_energy, mdea_volume),
    ]
    for left, right, energy, volume in (
        (
            "water",
            "mdea",
            0.5 * (water_energy + mdea_energy),
            cross_volume(water_volume, mdea_volume, sigma_water, sigma_mdea, 0.015),
        ),
        (
            "water",
            "mdeah-cation",
            0.5 * water_energy,
            cross_volume(water_volume, mdea_volume, sigma_water, sigma_mdea, 0.0),
        ),
        (
            "mdea",
            "mdeah-cation",
            0.5 * mdea_energy,
            cross_volume(mdea_volume, mdea_volume, sigma_mdea, sigma_mdea, 0.0),
        ),
    ):
        edges.extend(
            (
                edge((left, "a"), (right, "b"), energy, volume),
                edge((left, "b"), (right, "a"), energy, volume),
            )
        )
    return {"presets": [], "sites": sites, "edges": edges}


def build_parameters(temperature: float) -> epcsaft.Parameters:
    components = [
        component(
            "water",
            molar_mass=MOLAR_MASSES[0],
            charge_number=0,
            segment_count=1.2046,
            dispersion_energy_over_k=353.94,
        ),
        component(
            "mdea",
            molar_mass=MOLAR_MASSES[1],
            charge_number=0,
            segment_count=3.675,
            segment_diameter=3.563,
            dispersion_energy_over_k=228.711,
        ),
        component(
            "co2",
            molar_mass=MOLAR_MASSES[2],
            charge_number=0,
            segment_count=2.0729,
            segment_diameter=2.7852,
            dispersion_energy_over_k=169.21,
        ),
        component(
            "h-cation",
            molar_mass=MOLAR_MASSES[3],
            charge_number=1,
            segment_count=1.0,
            segment_diameter=3.4654,
            dispersion_energy_over_k=500.0,
            packing_diameter=0.88 * 3.4654,
            debye_huckel_diameter=0.88 * 3.4654,
        ),
        component(
            "oh-anion",
            molar_mass=MOLAR_MASSES[4],
            charge_number=-1,
            segment_count=1.0,
            segment_diameter=2.0177,
            dispersion_energy_over_k=650.0,
            packing_diameter=0.88 * 2.0177,
            debye_huckel_diameter=0.88 * 2.0177,
        ),
        component(
            "hco3-anion",
            molar_mass=MOLAR_MASSES[5],
            charge_number=-1,
            segment_count=1.0,
            segment_diameter=2.9296,
            dispersion_energy_over_k=70.0,
            packing_diameter=0.88 * 2.9296,
            debye_huckel_diameter=0.88 * 2.9296,
        ),
        component(
            "carbonate-anion",
            molar_mass=MOLAR_MASSES[6],
            charge_number=-2,
            segment_count=1.0,
            segment_diameter=2.4422,
            dispersion_energy_over_k=249.26,
            packing_diameter=0.88 * 2.4422,
            debye_huckel_diameter=0.88 * 2.4422,
        ),
        component(
            "mdeah-cation",
            molar_mass=MOLAR_MASSES[7],
            charge_number=1,
            segment_count=3.675,
            segment_diameter=3.563,
            dispersion_energy_over_k=228.711,
            packing_diameter=0.88 * 3.563,
            debye_huckel_diameter=0.88 * 3.563,
        ),
    ]
    water_mdea_k = interaction(temperature, (-0.015275, -5.24e-5, 7.9e-6))
    kij = {
        frozenset(("water", "mdea")): water_mdea_k,
        frozenset(("water", "co2")): interaction(
            temperature, (-0.022, 4.2e-4, -1.7e-6)
        ),
        frozenset(("water", "h-cation")): 0.25,
        frozenset(("water", "oh-anion")): -0.25,
        frozenset(("water", "hco3-anion")): 0.0,
        frozenset(("water", "carbonate-anion")): -0.25,
        # Uyan Section 4.3.2 and Table 5 inherit this interaction from
        # water/MDEA. The earlier Section 4.2.1 sentence saying zero is a
        # source inconsistency; the reported VLE calculation is decisive.
        frozenset(("water", "mdeah-cation")): water_mdea_k,
    }
    pairs = [
        pair(left, right, kij.get(frozenset((left, right)), 0.0))
        for left_index, left in enumerate(COMPONENTS)
        for right_index, right in enumerate(
            COMPONENTS[left_index + 1 :], left_index + 1
        )
        if CHARGES[left_index] * CHARGES[right_index] <= 0
    ]
    mapping = {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "document_id": "uyan-2015-explicit-ion-mdea-co2",
        "document_version": 1,
        "purpose": "external-document",
        "sources": [
            {
                "source_id": "uyan-2015",
                "citation": "Uyan et al., Fluid Phase Equilibria 393 (2015) 91-100.",
                "use_basis": "Tables 1-5 and Eqs. 5, 12-27.",
            },
            {
                "source_id": "hsieh-2007-thermoml",
                "citation": "Hsieh et al., J. Chem. Eng. Data 52 (2007) 619-623.",
                "use_basis": "six pure-MDEA relative-permittivity rows in NIST ThermoML.",
            },
            {
                "source_id": "uyan-2015-permittivity-transfer",
                "citation": "Uyan Eq. 5; Wangler campaign Archer-Wang approximation; Hsieh ThermoML values.",
                "use_basis": "explicitly labeled transfer approximation, not exact Floriano-Nascimento replay.",
            },
        ],
        "domains": [
            {
                "domain_id": PROVENANCE["domain_id"],
                "kind": "reported-conditions",
                "temperature_min": {"magnitude": 313.11, "unit": "kelvin"},
                "temperature_max": {"magnitude": 413.16, "unit": "kelvin"},
                "pressure_min": {"magnitude": 100.0, "unit": "pascal"},
                "pressure_max": {"magnitude": 1.0e7, "unit": "pascal"},
                "ion_mole_fraction_max": {"magnitude": 0.2, "unit": "dimensionless"},
            }
        ],
        "components": components,
        "pairs": pairs,
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
                "ionic_dispersion": "exclude-same-sign-ion-pairs",
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
        "correlations": [
            {
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
            },
            temperature_polynomial_correlation(
                "water",
                (
                    water_permittivity_transfer_approximation(298.15),
                    2.0 * 7.6555618295e-4 * 298.15 - 8.1783881423e-1,
                    7.6555618295e-4,
                ),
                reference_temperature=298.15,
                provenance={
                    "source_id": "uyan-2015-permittivity-transfer",
                    "locator": "Archer-Wang water approximation represented at Tref=298.15 K",
                    "domain_id": PROVENANCE["domain_id"],
                },
            ),
            temperature_polynomial_correlation(
                "mdea",
                (
                    mdea_permittivity(298.15),
                    -0.102457142857,
                ),
                reference_temperature=298.15,
                provenance={
                    "source_id": "hsieh-2007-thermoml",
                    "locator": "Least-squares linear fit of six NIST ThermoML rows",
                    "domain_id": PROVENANCE["domain_id"],
                },
            ),
        ],
        "topology": association_topology(temperature),
    }
    return epcsaft.Parameters.from_mapping(mapping)


def problem_for_row(
    row: dict[str, str], pressure_seed_pa: float,
    species_reference: dict[str, equilibrium.SpeciesReference] | None = None,
) -> equilibrium.Problem:
    temperature = float(row["temperature_k"])
    water = 1.0 / MOLAR_MASSES[0]
    mdea, co2 = (
        float(row[f"{name}_molality_mol_per_kg_water"]) for name in ("mdea", "co2")
    )
    hydrated = min(0.80 * co2, 0.80 * mdea)
    carbonate = min(0.01 * hydrated, 0.01 * co2)
    extents = (1.0e-8, hydrated, carbonate, -(hydrated + carbonate))
    feed = tuple(
        amount + math.fsum(REACTIONS[r][i] * extents[r] for r in range(4))
        for i, amount in enumerate((water, mdea, co2, 0.0, 0.0, 0.0, 0.0, 0.0))
    )
    total = math.fsum(feed)
    liquid = tuple(n / total for n in feed)
    neutral_total = math.fsum(liquid[:3])
    reactions = []
    for i, (nu, coefficients) in enumerate(
        zip(REACTIONS, LN_K_COEFFICIENTS, strict=True)
    ):
        a, b, c, d = coefficients
        # Polynomial uses ln(T/Tref); preserve the source ln(T) coefficients.
        a += c * math.log(313.15)
        reactions.append(
            equilibrium.Reaction(
                dict(zip(COMPONENTS, nu, strict=True)),
                name=f"uyan-2015-K{i + 1}",
                correlation=equilibrium.ReactionLogPolynomial(
                    a,
                    b,
                    c,
                    d,
                    313.15,
                    313.15,
                    393.15,
                    "uyan-raw-mole-fraction-infinite-dilution",
                ),
                reference=equilibrium.ReactionReference(
                    equilibrium.ReferenceSourceBasis.RawMoleFractionInfiniteDilution,
                    reference_pressure_pa=1.0e5,
                ),
            )
        )
    return equilibrium.Problem(
        T=temperature,
        P=equilibrium.Free(pressure_seed_pa),
        phases=[
            equilibrium.Phase(
                "liquid", kind="liquid", support=COMPONENTS, composition_guess=liquid
            ),
            equilibrium.Phase(
                "vapor",
                kind="vapor",
                support=COMPONENTS[:3],
                amount=equilibrium.Pinned(0.0),
                composition_guess=tuple(x / neutral_total for x in liquid[:3]),
            ),
        ],
        feed=equilibrium.Amounts(dict(zip(COMPONENTS, feed, strict=True))),
        reactions=reactions,
        row_tolerance=1e-10,
        neutral_reference=equilibrium.NeutralReference(
            solvent="water",
            solute_molality_ratios={
                "mdea": 0.125,
                "co2": 0.125,
                "h-cation": 0.5,
                "oh-anion": 0.25,
                "hco3-anion": 0.125,
                "carbonate-anion": 0.125,
                "mdeah-cation": 0.125,
            },
            standard_molality_mol_per_kg=1.0,
            species_reference=species_reference,
        ),
    )


def read_rows(row_ids: set[str] | None) -> list[dict[str, str]]:
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if row_ids is None:
        return rows
    selected = [row for row in rows if row["source_row"] in row_ids]
    missing = row_ids - {row["source_row"] for row in selected}
    if missing:
        raise ValueError(f"unknown source rows: {sorted(missing)}")
    return selected


def evaluate_row(payload: tuple[dict[str, str], float],
                 species_reference: dict[str, equilibrium.SpeciesReference] | None = None) -> dict[str, object]:
    row, pressure_seed_pa = payload
    row_started = time.perf_counter()
    attempted_seeds = tuple(dict.fromkeys((pressure_seed_pa, 1.0e5, 3.0e5, 5.0e5)))
    model = epcsaft.Mixture(build_parameters(float(row["temperature_k"])))
    problem = problem_for_row(row, pressure_seed_pa, species_reference)
    anchor = json.loads((CAMPAIGN / "inputs/starts/313-15-k-t1-r08.json").read_text())
    use_anchor = row["source_row"] in anchor["target_source_rows"]
    starts = []
    for value in attempted_seeds:  # physical input declaration only
        seed = equilibrium.EquilibriumStart()
        seed.pressure = max(
            value, float(row["co2_partial_pressure_kpa"]) * 1.0e3 + 1.0e5
        )
        if use_anchor:
            seed.compositions = [
                anchor["mole_fractions"][:8],
                anchor["mole_fractions"][8:11],
            ]
            seed.packing_fractions = anchor["packing_fractions"]
        starts.append(seed)
    limits = equilibrium.EquilibriumWorkLimits()
    solved = equilibrium.solve_equilibrium_starts(model, problem, starts, limits)
    result = solved.results[solved.selected] if solved.selected is not None else None
    selected_seed = (
        starts[solved.selected].pressure if solved.selected is not None else None
    )
    record: dict[str, object] = dict(row)
    mdea_mass = float(row["mdea_molality_mol_per_kg_water"]) * MOLAR_MASSES[1]
    record["mdea_mass_fraction_in_co2_free_solvent"] = mdea_mass / (1.0 + mdea_mass)
    record["grepe_status"] = "evaluated" if result is not None else "failed"
    record["grepe_solver_status"] = (
        result.solver_status if result else "; ".join(solved.causes)
    )
    record["reaction_constant_lane"] = "uyan-source-raw-mole-fraction-k"
    record["source_reference_lane"] = (
        "pure-water-co2-vapor-other-components-infinite-dilution"
        if species_reference and species_reference.get("co2") == equilibrium.SpeciesReference.PureGas
        else "pure-water-all-other-components-infinite-dilution"
    )
    record["fixed_composition_packing_anchor"] = (
        anchor["source_row"] if use_anchor else ""
    )
    record["wall_time_s"] = time.perf_counter() - row_started
    attempted_count = sum(cause != "unattempted" for cause in solved.causes)
    record["nominal_pressure_seed_attempts_kpa"] = ";".join(
        f"{v / 1e3:g}" for v in attempted_seeds[:attempted_count]
    )
    record["pressure_seed_attempts_kpa"] = ";".join(
        f"{s.pressure / 1e3:g}" for s in starts[:attempted_count]
    )
    record["selected_nominal_pressure_seed_kpa"] = (
        attempted_seeds[solved.selected] / 1e3 if solved.selected is not None else None
    )
    record["selected_pressure_seed_kpa"] = (
        selected_seed / 1e3 if selected_seed is not None else None
    )
    record["attempt_causes"] = "; ".join(solved.causes[:attempted_count])
    record["maximum_original_row"] = (
        max(abs(r) for r in result.residuals) if result else None
    )
    if result is not None:
        total_pressure_pa = result.pressure
        # Native PartialPressure observation projects P*y_CO2 through the shared driver.
        problem.P = equilibrium.Free(total_pressure_pa)
        for i, phase in enumerate(problem.phases):
            phase.composition_guess = tuple(
                result.mole_fractions[i * 8 : i * 8 + (8 if i == 0 else 3)]
            )
            phase.packing_guess = result.packing_fractions[i]
        observation = regression.observation(
            model.parameters,
            "partial_pressure",
            problem,
            observed=float(row["co2_partial_pressure_kpa"]) * 1.0e3,
            form="difference",
            scale=1.0,
            mass_basis=False,
            phase="vapor",
            species=("co2",),
        )
        prediction = regression.predict(model, observation)
        predicted_pa = prediction.value
        record["partial_pressure_native_message"] = prediction.message
        if total_pressure_pa is None or predicted_pa is None:
            raise RuntimeError("evaluated GREPE result omitted a requested output")
        observed_pa = float(row["co2_partial_pressure_kpa"]) * 1.0e3
        record.update(
            {
                "grepe_total_pressure_kpa": total_pressure_pa / 1.0e3,
                "grepe_co2_partial_pressure_kpa": predicted_pa / 1.0e3,
                "grepe_co2_partial_pressure_relative_error": predicted_pa / observed_pa
                - 1.0,
                "liquid_density_mol_m3": result.molar_densities[0],
                "vapor_density_mol_m3": result.molar_densities[1],
                "vapor_water_mole_fraction": result.mole_fractions[8],
                "vapor_mdea_mole_fraction": result.mole_fractions[9],
                "vapor_co2_mole_fraction": result.mole_fractions[10],
            }
        )
    else:
        record["grepe_failure"] = "; ".join(solved.causes)
    record["_attempts"] = [
        {
            "source_row": row["source_row"],
            "start_index": i,
            "pressure_seed_pa": attempted_seeds[i],
            "actual_pressure_seed_pa": starts[i].pressure,
            "selected": solved.selected == i,
            "cause": cause,
            "seconds": solved.attempt_seconds[i],
            "solver_status": attempt.solver_status if attempt else "",
            "iterations": attempt.iterations if attempt else 0,
            "maximum_original_row": max(map(abs, attempt.residuals), default=None)
            if attempt
            else None,
        }
        for i, (attempt, cause) in enumerate(
            zip(solved.results, solved.causes, strict=True)
        )
    ]
    record["_rows"] = [
        {
            "source_row": row["source_row"],
            "start_index": i,
            "row": j,
            "kind": str(label.kind),
            "equation": label.eqid,
            "residual": raw,
            "classified_residual": classified,
        }
        for i, attempt in enumerate(solved.results)
        if attempt is not None and attempt.residuals
        for j, (label, raw, classified) in enumerate(
            zip(
                attempt.rows,
                attempt.residuals,
                attempt.classified_residuals,
                strict=True,
            )
        )
    ]
    record["_references"] = [
        {
            "source_row": row["source_row"],
            "start_index": i,
            "reaction": diagnostic.reaction_index,
            "status": str(diagnostic.status),
            "reference_pressure_pa": diagnostic.reference_pressure_pa,
            "terminal_contraction": diagnostic.terminal_contraction,
            "observed_terminal_change": diagnostic.observed_terminal_change,
            "message": diagnostic.message,
        }
        for i, attempt in enumerate(solved.results)
        if attempt is not None
        for diagnostic in attempt.neutral_reference_diagnostics
    ]
    print(
        row["source_row"],
        record["grepe_status"],
        record["grepe_solver_status"],
        f"pCO2={record.get('grepe_co2_partial_pressure_kpa', '-')}",
        f"{record['wall_time_s']:.3f}s",
        flush=True,
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--row", action="append", dest="rows")
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    parser.add_argument("--maximum-mdea-molality", type=float)
    parser.add_argument("--minimum-mdea-mass-fraction", type=float)
    parser.add_argument("--maximum-mdea-mass-fraction", type=float)
    parser.add_argument("--pressure-seed-kpa", type=float)
    arguments = parser.parse_args()
    if arguments.workers < 1:
        raise ValueError("workers must be positive")
    rows = read_rows(None if arguments.rows is None else set(arguments.rows))
    if arguments.maximum_mdea_molality is not None:
        rows = [
            row
            for row in rows
            if float(row["mdea_molality_mol_per_kg_water"])
            <= arguments.maximum_mdea_molality
        ]
    if arguments.minimum_mdea_mass_fraction is not None:
        rows = [
            row
            for row in rows
            if (
                float(row["mdea_molality_mol_per_kg_water"])
                * MOLAR_MASSES[1]
                / (1.0 + float(row["mdea_molality_mol_per_kg_water"]) * MOLAR_MASSES[1])
                >= arguments.minimum_mdea_mass_fraction
            )
        ]
    if arguments.maximum_mdea_mass_fraction is not None:
        rows = [
            row
            for row in rows
            if (
                float(row["mdea_molality_mol_per_kg_water"])
                * MOLAR_MASSES[1]
                / (1.0 + float(row["mdea_molality_mol_per_kg_water"]) * MOLAR_MASSES[1])
                <= arguments.maximum_mdea_mass_fraction
            )
        ]
    started = time.perf_counter()
    payloads = [
        (
            row,
            (
                arguments.pressure_seed_kpa
                if arguments.pressure_seed_kpa is not None
                else (300.0 if float(row["temperature_k"]) >= 380.0 else 100.0)
            )
            * 1.0e3,
        )
        for row in rows
    ]
    if arguments.workers == 1:
        plotted_rows = [evaluate_row(payload) for payload in payloads]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=arguments.workers
        ) as executor:
            plotted_rows = list(executor.map(evaluate_row, payloads))
    if arguments.no_write:
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for field, filename in (
        ("_attempts", "attempts.csv"),
        ("_rows", "rows.csv"),
        ("_references", "references.csv"),
    ):
        diagnostics = [item for record in plotted_rows for item in record.pop(field)]
        if diagnostics:
            with (OUTPUT / filename).open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]))
                writer.writeheader()
                writer.writerows(diagnostics)
    table = OUTPUT / "uyan-2015-ermatchkov-2006-grepe-plotted-data.csv"
    fieldnames = list(dict.fromkeys(key for row in plotted_rows for key in row))
    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(plotted_rows)
    evaluated = [row for row in plotted_rows if row["grepe_status"] == "evaluated"]
    relative_errors = [
        abs(float(cast(float, row["grepe_co2_partial_pressure_relative_error"])))
        for row in evaluated
    ]
    temperatures = sorted(
        {cast(str, row["temperature_k"]) for row in evaluated}, key=float
    )
    ard_by_temperature = {
        temperature: 100.0
        * math.fsum(
            abs(float(cast(float, row["grepe_co2_partial_pressure_relative_error"])))
            for row in evaluated
            if row["temperature_k"] == temperature
        )
        / sum(row["temperature_k"] == temperature for row in evaluated)
        for temperature in temperatures
    }
    summary = {
        "wheel_sha256": hashlib.sha256(
            next((ROOT / "build/environment-wheel").glob("epcsaft-*.whl")).read_bytes()
        ).hexdigest(),
        "source_reference_convention": "water pure-component; all other components infinite dilution",
        "fixed_start_anchor": "inputs/starts/313-15-k-t1-r08.json",
        "fixed_start_anchor_rows": 4,
        "pressure_only_replay": "pressure-only/uyan-2015-ermatchkov-2006-grepe-summary.json",
        "pressure_seed_interpretation": (
            "Pressure seed fields report declared physical guesses "
            "max(nominal list entry, observed CO2 partial pressure + 100 kPa); "
            "nominal_pressure_seed_attempts_kpa and selected_nominal_pressure_seed_kpa "
            "retain the nominal list entries."
        ),
        "selected_pressure_guess_differs_from_nominal_count": sum(
            r["selected_pressure_seed_kpa"] != r["selected_nominal_pressure_seed_kpa"]
            for r in evaluated
        ),
        "reference_pressure_pa": 1.0e5,
        "row_tolerance": 1e-10,
        "maximum_original_row": max(
            (float(r["maximum_original_row"]) for r in evaluated), default=None
        ),
        "source_row_count": len(rows),
        "evaluated_row_count": len(evaluated),
        "failed_row_count": len(rows) - len(evaluated),
        "grepe_ard_percent": None
        if not relative_errors
        else 100.0 * math.fsum(relative_errors) / len(relative_errors),
        "grepe_ard_by_temperature_percent": ard_by_temperature,
        "liquid_denser_than_vapor_count": sum(
            float(cast(float, row["liquid_density_mol_m3"]))
            > float(cast(float, row["vapor_density_mol_m3"]))
            for row in evaluated
        ),
        "retried_row_count": sum(
            ";" in str(row["pressure_seed_attempts_kpa"]) for row in plotted_rows
        ),
        "wall_time_s": time.perf_counter() - started,
        "phase_topology": ["finite reacting liquid", "incipient neutral vapor"],
        "parameter_lane": "Uyan 2015 transfer reconstruction with explicit Archer-Wang water approximation and Hsieh ThermoML pure-MDEA fit",
        "interpretation": "Direct source-data GREPE reconstruction; not an exact author-code replay.",
    }
    (OUTPUT / "uyan-2015-ermatchkov-2006-grepe-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


def declarations():
    from epcsaft import Mixture

    return [
        (
            Mixture(build_parameters(t)),
            [problem_for_row(r, 1e5) for r in read_rows(None) if float(r['temperature_k']) == t],
        )
        for t in (313.15, 353.15, 393.15)
    ]


if __name__ == "__main__":
    main()
