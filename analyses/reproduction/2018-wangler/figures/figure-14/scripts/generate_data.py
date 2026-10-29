"""Wangler 2018 Table 10 (Figure 14): native reacting bubble pressures at 303 K and 60 wt% MDEA.

Mole-fraction reaction constants (Eq. 19, Table 2) referenced to pure water and to infinite
dilution in water for every other species (Eqs. 13-21), declared to the native neutral reference.
Every lane uses MDEAH+ m = 3.675, packing 0.88 sigma and Debye-Hueckel diameter sigma (#330).
The printed lane keeps Table 6's water-MDEAH+ k_ij; the fitted lanes predict Table 10 from the
Figure 10 k_ij(303 K) (Eq. 34, linear) and its +/- u(phi) allowance. S0 excludes ion
self-dispersion; S1 keeps it (all-ion-pairs with k_ij = 1 on same-sign cross pairs). No
parameter is fitted to Table 10; the per-row matching k_ij(303 K) is a diagnostic.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import runpy
from pathlib import Path

import epcsaft
from epcsaft import GAS_CONSTANT_J_PER_MOL_K, equilibrium, regression

ROOT = Path(__file__).resolve().parents[6]
CAMPAIGN = Path(__file__).resolve().parents[3]
SOURCE = CAMPAIGN / "inputs/source" / "wangler-2018-table-10.csv"
OUTPUT = Path(__file__).resolve().parents[1] / "output"
FIGURE_10 = CAMPAIGN / "figures/figure-10/output/wangler-2018-mdeah-osmotic-summary.json"
COMMON = runpy.run_path(
    str(ROOT / "analyses/reproduction/2018-wangler/figures/figure-10/scripts/generate_data.py")
)
component = COMMON["_component"]
pair = COMMON["_pair"]
edge = COMMON["_edge"]
water_sigma = COMMON["_water_sigma"]
water_permittivity = COMMON["_water_permittivity"]
interaction = COMMON["_interaction"]
provenance = {
    "source_id": "wangler-2018",
    "locator": "Tables 2, 4-7, 10 and Eqs. 3, 8, 13-21, 34-35",
    "domain_id": "wangler-2018-reactive-mdea-co2",
}
component.__globals__["PROVENANCE"] = provenance
pair.__globals__["PROVENANCE"] = provenance
edge.__globals__["PROVENANCE"] = provenance

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
# Eq. 19: ln K = c1 + c2/T + c3 ln T + c4 T, mole-fraction basis, unconverted.
LN_K_COEFFICIENTS = (
    (132.9, -13450.0, -22.48, 0.0),
    (212.7, -11330.0, -33.84, -1.815e-3),
    (287.4, -13650.0, -48.88, 3.023e-2),
    (-83.49, -819.7, 10.98, 0.0),
)


def correlation(component_id: str, family: str, value: float) -> dict[str, object]:
    return {
        "correlation_id": f"component/{component_id}/{family}/constant",
        "component_id": component_id,
        "family": family,
        "form": "constant",
        "independent_variables": [],
        "constant": {
            "identity": f"component/{component_id}/{family}/constant/constant",
            "value": {"magnitude": value, "unit": "dimensionless"},
        },
        "provenance": provenance,
    }


def association_topology(
    temperature: float, water_mdeah_l: float
) -> dict[str, object]:
    sites = [
        {
            "component_id": component_id,
            "site_id": site,
            "site_role": role,
            "multiplicity": multiplicity,
            "provenance": provenance,
        }
        for component_id, multiplicity in (("water", 1), ("mdea", 2), ("mdeah-cation", 2))
        for site, role in (("a", "donor"), ("b", "acceptor"))
    ]
    sigma_water = water_sigma(temperature)
    sigma_mdea = 3.5630
    water_energy, water_volume = 2425.6, 0.0450
    mdea_energy, mdea_volume = 2046.6, 0.1238

    def cross_volume(left_volume: float, right_volume: float, left_sigma: float,
                     right_sigma: float, correction: float) -> float:
        sigma_cross = 0.5 * (left_sigma + right_sigma)
        return (
            math.sqrt(left_volume * right_volume)
            * (math.sqrt(left_sigma * right_sigma) / sigma_cross) ** 3
            * (1.0 - correction)
        )

    water_mdea_l = interaction(temperature, (0.600, -2.57e-3, -1.70e-5))
    water_mdea_energy = 0.5 * (water_energy + mdea_energy)
    water_mdea_volume = cross_volume(
        water_volume, mdea_volume, sigma_water, sigma_mdea, water_mdea_l
    )
    water_mdeah_volume = cross_volume(
        water_volume, mdea_volume, sigma_water, sigma_mdea, water_mdeah_l
    )
    mdea_mdeah_volume = cross_volume(
        mdea_volume, mdea_volume, sigma_mdea, sigma_mdea, 0.0
    )
    edges = [
        edge(("water", "a"), ("water", "b"), water_energy, water_volume),
        edge(("mdea", "a"), ("mdea", "b"), mdea_energy, mdea_volume),
    ]
    for left, right, energy, volume in (
        ("water", "mdea", water_mdea_energy, water_mdea_volume),
        ("water", "mdeah-cation", 0.5 * water_energy, water_mdeah_volume),
        ("mdea", "mdeah-cation", 0.5 * mdea_energy, mdea_mdeah_volume),
    ):
        edges.extend(
            (
                edge((left, "a"), (right, "b"), energy, volume),
                edge((left, "b"), (right, "a"), energy, volume),
            )
        )
    return {"presets": [], "sites": sites, "edges": edges}


def parameter_mapping(temperature: float, water_mdeah_k: float, ionic_dispersion: str) -> dict:
    """Paper parameters (Tables 4-7) with a declared water-MDEAH+ k_ij."""
    eps_water = water_permittivity(temperature)
    eps_mdea = 22.4
    eps_blend = 0.4 * eps_water + 0.6 * eps_mdea
    water_mdeah_l = interaction(temperature, (0.372, 1.07e-3))
    components = [
        component(
            "water", molar_mass=MOLAR_MASSES[0], charge_number=0,
            segment_count=1.2046, dispersion_energy_over_k=353.94,
        ),
        component(
            "mdea", molar_mass=MOLAR_MASSES[1], charge_number=0,
            segment_count=3.6750, segment_diameter=3.5630,
            dispersion_energy_over_k=228.71,
        ),
        component(
            "co2", molar_mass=MOLAR_MASSES[2], charge_number=0,
            segment_count=2.0729, segment_diameter=2.7852,
            dispersion_energy_over_k=169.21,
        ),
        component(
            "h-cation", molar_mass=MOLAR_MASSES[3], charge_number=1,
            segment_count=1.0, segment_diameter=3.4654,
            dispersion_energy_over_k=500.0, packing_diameter=0.88 * 3.4654,
            debye_huckel_diameter=3.4654,
        ),
        component(
            "oh-anion", molar_mass=MOLAR_MASSES[4], charge_number=-1,
            segment_count=1.0, segment_diameter=2.0177,
            dispersion_energy_over_k=650.0, packing_diameter=0.88 * 2.0177,
            debye_huckel_diameter=2.0177,
        ),
        component(
            "hco3-anion", molar_mass=MOLAR_MASSES[5], charge_number=-1,
            segment_count=1.0, segment_diameter=2.9296,
            dispersion_energy_over_k=70.0, packing_diameter=0.88 * 2.9296,
            debye_huckel_diameter=2.9296,
        ),
        component(
            "carbonate-anion", molar_mass=MOLAR_MASSES[6], charge_number=-2,
            segment_count=1.0, segment_diameter=2.4422,
            dispersion_energy_over_k=249.26, packing_diameter=0.88 * 2.4422,
            debye_huckel_diameter=2.4422,
        ),
        component(
            "mdeah-cation", molar_mass=MOLAR_MASSES[7], charge_number=1,
            segment_count=3.6750, segment_diameter=3.5630,
            dispersion_energy_over_k=228.71, packing_diameter=0.88 * 3.5630,
            debye_huckel_diameter=3.5630,
        ),
    ]
    kij = {
        frozenset(("water", "mdea")): interaction(
            temperature, (-0.184, -1.21e-3, 2.70e-6, 1.53e-7)
        ),
        frozenset(("water", "co2")): interaction(
            temperature, (-2.15e-2, 4.20e-4, -1.70e-6)
        ),
        frozenset(("water", "h-cation")): 0.25,
        frozenset(("water", "oh-anion")): -0.25,
        frozenset(("water", "hco3-anion")): interaction(
            temperature, (-7.65e-2, -6.25e-4, -6.25e-5)
        ),
        frozenset(("water", "carbonate-anion")): -0.25,
        frozenset(("water", "mdeah-cation")): water_mdeah_k,
    }
    # S1 keeps ion self-dispersion only: same-sign cross pairs get k_ij = 1 (no dispersion).
    self_kept = ionic_dispersion == "all-ion-pairs"
    pairs = [
        pair(left, right, 1.0 if same_sign else kij.get(frozenset((left, right)), 0.0))
        for index, left in enumerate(COMPONENTS)
        for right in COMPONENTS[index + 1 :]
        if not (same_sign := CHARGES[index] * CHARGES[COMPONENTS.index(right)] > 0) or self_kept
    ]
    mapping = {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "document_id": "wangler-2018-reactive-mdea-co2-probe",
        "document_version": 1,
        "purpose": "external-document",
        "sources": [{
            "source_id": "wangler-2018",
            "citation": "Wangler et al., Fluid Phase Equilibria 461 (2018) 15-27.",
            "use_basis": "Tables 2, 4-7, 10 and Eqs. 3, 8, 13-21, 34-35.",
        }],
        "domains": [{
            "domain_id": "wangler-2018-reactive-mdea-co2",
            "kind": "reported-conditions",
            "temperature_min": {"magnitude": 293.0, "unit": "kelvin"},
            "temperature_max": {"magnitude": 413.0, "unit": "kelvin"},
            "pressure_min": {"magnitude": 1.0e4, "unit": "pascal"},
            "pressure_max": {"magnitude": 5.0e6, "unit": "pascal"},
            "ion_mole_fraction_max": {"magnitude": 0.60, "unit": "dimensionless"},
        }],
        "components": components,
        "pairs": pairs,
        "model_families": [
            {"family_id": "model/association", "kind": "association", "choice": "general-site", "provenance": provenance},
            {"family_id": "model/base", "kind": "base", "choice": "pc-saft", "provenance": provenance},
            {"family_id": "model/electrolyte", "kind": "electrolyte", "choice": "fully-dissociated-debye-huckel",
             "ionic_dispersion": ionic_dispersion, "provenance": provenance},
            {"family_id": "model/relative_permittivity", "kind": "permittivity", "choice": "solvent-only", "provenance": provenance},
        ],
        "model_coefficients": [],
        "correlations": [
            {
                "correlation_id": "component/water/segment_diameter/constant-plus-sum-of-exponentials",
                "component_id": "water", "family": "segment_diameter",
                "form": "constant-plus-sum-of-exponentials",
                "independent_variables": ["temperature"],
                "constant": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/constant", "value": {"magnitude": 2.7927, "unit": "angstrom"}},
                "terms": [
                    {"amplitude": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/amplitude", "value": {"magnitude": 10.11, "unit": "angstrom"}}, "exponent_coefficient": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-0/exponent_coefficient", "value": {"magnitude": -0.01775, "unit": "1 / kelvin"}}},
                    {"amplitude": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/amplitude", "value": {"magnitude": -1.417, "unit": "angstrom"}}, "exponent_coefficient": {"identity": "component/water/segment_diameter/constant-plus-sum-of-exponentials/term-1/exponent_coefficient", "value": {"magnitude": -0.01146, "unit": "1 / kelvin"}}},
                ],
                "provenance": provenance,
            },
            # Wangler/Uyan evaluate epsilon from the ion- and CO2-free
            # apparent 60 wt% MDEA solvent, not from the reacted true-species
            # composition. Assigning that fixed blend value to every neutral
            # carrier makes the Provider's generic solvent-only mass average
            # exactly reproduce that declared medium as speciation changes.
            correlation("water", "relative_permittivity", eps_blend),
            correlation("mdea", "relative_permittivity", eps_blend),
            correlation("co2", "relative_permittivity", eps_blend),
        ],
        "topology": association_topology(temperature, water_mdeah_l),
    }
    return mapping


def build_parameters(temperature: float, water_mdeah_k: float, ionic_dispersion: str) -> epcsaft.Parameters:
    return epcsaft.Parameters.from_mapping(parameter_mapping(temperature, water_mdeah_k, ionic_dispersion))


TEMPERATURE = 303.0
REFERENCE_TEMPERATURE = 298.15
REFERENCE = equilibrium.NeutralReference(
    solvent="water",
    # Charge-neutral dilution path to the infinite-dilution limit (the Uyan declaration).
    solute_molality_ratios={
        "mdea": 0.125, "co2": 0.125, "h-cation": 0.5, "oh-anion": 0.25,
        "hco3-anion": 0.125, "carbonate-anion": 0.125, "mdeah-cation": 0.125,
    },
    standard_molality_mol_per_kg=1.0,
)


def problem(alpha: float, pressure_seed_pa: float) -> equilibrium.Problem:
    mdea = 1.0
    water = (0.4 / 0.6) * MOLAR_MASSES[1] / MOLAR_MASSES[0]
    start_extents = (2.0e-6, 0.1, 1.0e-3, -9.9e-3)  # the historical declared start
    feed = tuple(
        amount + math.fsum(REACTIONS[r][i] * start_extents[r] for r in range(len(REACTIONS)))
        for i, amount in enumerate((water, mdea, alpha, 0.0, 0.0, 0.0, 0.0, 0.0))
    )
    liquid = tuple(n / math.fsum(feed) for n in feed)
    reactions = [
        equilibrium.Reaction(
            dict(zip(COMPONENTS, nu, strict=True)),
            name=f"wangler-2018-K{index}",
            correlation=equilibrium.ReactionLogPolynomial(
                # The native form uses ln(T/Tref); the source coefficient c1 absorbs c3 ln(Tref).
                c1 + c3 * math.log(REFERENCE_TEMPERATURE), c2, c3, c4,
                REFERENCE_TEMPERATURE, 293.0, 413.0, "wangler-raw-mole-fraction-infinite-dilution",
            ),
            reference=equilibrium.ReactionReference(
                equilibrium.ReferenceSourceBasis.RawMoleFractionInfiniteDilution,
                # Wangler p. 18, Eqs. 20–21: standards at the system pressure.
                reference_pressure_pa=None,
            ),
        )
        for index, (nu, (c1, c2, c3, c4)) in enumerate(zip(REACTIONS, LN_K_COEFFICIENTS, strict=True), 1)
    ]
    return equilibrium.Problem(
        T=TEMPERATURE,
        P=equilibrium.Free(pressure_seed_pa),
        phases=[
            equilibrium.Phase("liquid", kind="liquid", support=COMPONENTS, composition_guess=liquid),
            equilibrium.Phase(
                "vapor", kind="vapor", support=COMPONENTS[:3], amount=equilibrium.Pinned(0.0),
                # CO2-rich vapour start (#330): the liquid's water-rich neutral composition
                # stops every row at local infeasibility from seeds of 9 bar and above.
                composition_guess=(0.01, 1.0e-4, 0.9899),
            ),
        ],
        feed=equilibrium.Amounts(dict(zip(COMPONENTS, feed, strict=True))),
        reactions=reactions,
        row_tolerance=1e-10,
        neutral_reference=REFERENCE,
    )


def affinity(result) -> float:
    """Largest |sum_i nu_ri r_i| / sum|nu_ri| over the liquid chemical rows (moiety multipliers cancel)."""
    chemical = {row.index: value for row, value in zip(result.rows, result.residuals, strict=True)
                if row.kind.name == "Chemical" and row.phase == 0}
    return max(abs(math.fsum(v * chemical[i] for i, v in enumerate(nu) if v)) / math.fsum(map(abs, nu))
               for nu in REACTIONS)


def liquid_ln_phi(model, x, *, pressure=None, density=None) -> list[float]:
    state = model.state(TEMPERATURE, P=pressure, x=x, phase="liquid") if density is None else \
        model.state(TEMPERATURE, rho=density, x=x)
    pressure = state.pressure if pressure is None else pressure
    z = pressure / (state.molar_density * GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE)
    return [mu - math.log(z) for mu in state.residual_chemical_potential_over_rt]


def ln_gamma_star(model, result, index: int) -> float:
    """ln phi_i(T, P, x) - ln phi_i at 1e-12 mol/kg in pure water along the declared ratios (same T, P)."""
    x = list(result.mole_fractions[: len(COMPONENTS)])
    trace = [(1.0 / MOLAR_MASSES[0] if c == "water" else 0.0)
             + 1.0e-12 * REFERENCE.solute_molality_ratios.get(c, 0.0) for c in COMPONENTS]
    trace = [n / math.fsum(trace) for n in trace]
    return (liquid_ln_phi(model, x, density=result.molar_densities[0])[index]
            - liquid_ln_phi(model, trace, pressure=result.pressure)[index])


DISPERSION = {"S0": "exclude-same-sign-ion-pairs", "S1": "all-ion-pairs"}
PRINTED_SEED_BAR = 257.0  # printed lane: rows 2-3 stop at LOCAL_INFEASIBILITY from the observed-pressure seed


def lanes() -> list[tuple[str, str, float, float | None]]:
    """(lane, conventions, k_ij(303 K), declared seed in bar or None for the observed pressure)."""
    figure_10 = json.loads(FIGURE_10.read_text())["lanes"]
    result = [("printed-S0", "S0", interaction(TEMPERATURE, (-0.184, -1.21e-3, 2.70e-6, 1.53e-7)), PRINTED_SEED_BAR)]
    for name in DISPERSION:
        k, allowance = figure_10[name]["k_ij_303_k"], figure_10[name]["k_ij_303_k_allowance"]
        result += [(f"fitted-{name}", name, k, None), (f"fitted-{name}-minus", name, k - allowance, None),
                   (f"fitted-{name}-plus", name, k + allowance, None)]
    return result


def matching_k_ij(parameters, alpha: float, observed_pa: float, start: float) -> float:
    """Native one-coordinate fit of the water-MDEAH+ k_ij(303 K) to one Table 10 pressure (diagnostic)."""
    observation = regression.observation(parameters, "pressure", problem(alpha, observed_pa),
                                         observed=observed_pa, form="log_ratio", mass_basis=False)
    k_ij = regression.coordinate(parameters, [("k_ij", ("water", "mdeah-cation"), 1)],
                                 origin=0.0, scale=1.0, bounds=(-1.0, 1.0))
    controls = regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=120.0,
                                      function_tolerance=1e-10, parameter_tolerance=1e-10, gradient_tolerance=1e-10)
    fit = regression.fit(parameters, [k_ij], [observation], weights=[1.0], controls=controls, holds=(),
                         starts=((start,),))[0]
    if not fit.usable or abs(fit.residuals[0]) > 1e-8:
        raise RuntimeError(f"matching k_ij at alpha {alpha}: {fit.message}")
    return fit.physical[0]


def main() -> None:
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        source_rows = list(csv.DictReader(stream))
    rows = []
    for lane, conventions, k_ij, seed in lanes():
        rows += evaluate_lane(lane, conventions, k_ij, seed, source_rows)
    write_outputs(rows)


def evaluate_lane(lane, conventions, k_ij, seed, source_rows) -> list[dict]:
    parameters = build_parameters(TEMPERATURE, k_ij, DISPERSION[conventions])
    model = epcsaft.Mixture(parameters)
    rows = []
    for source in source_rows:
        alpha = float(source["co2_loading_mol_per_mol_mdea"])
        observed = float(source["experimental_total_pressure_bar"])
        uncertainty = float(source["pressure_standard_uncertainty_bar"])
        seed_bar = observed if seed is None else seed
        row = {"lane": lane, "conventions": conventions, "k_ij_303_k": k_ij, "temperature_k": TEMPERATURE,
               "co2_loading_mol_per_mol_mdea": alpha, "experimental_total_pressure_bar": observed,
               "pressure_standard_uncertainty_bar": uncertainty, "pressure_seed_bar": seed_bar}
        result = equilibrium.solve_equilibrium(model, problem(alpha, seed_bar * 1.0e5))
        row.update(success=result.success, solver_status=result.solver_status,
                   failure="" if result.success else result.message)
        if result.success:
            pressure = result.pressure / 1.0e5
            row.update({
                "native_pressure_bar": pressure,
                "relative_error": pressure / observed - 1.0,
                "deviation_over_uncertainty": (pressure - observed) / uncertainty,
                "liquid_density_mol_m3": result.molar_densities[0],
                "liquid_density_kg_m3": result.molar_densities[0] * math.fsum(
                    x * mw for x, mw in zip(result.mole_fractions[:8], MOLAR_MASSES)),
                "vapor_density_kg_m3": result.molar_densities[1] * math.fsum(
                    x * mw for x, mw in zip(result.mole_fractions[8:], MOLAR_MASSES[:3])),
                "mdeah_fraction_of_mdea": result.mole_fractions[7]
                / (result.mole_fractions[1] + result.mole_fractions[7]),
                "mdeah_ln_gamma_star": ln_gamma_star(model, result, 7),
                "vapor_water_mole_fraction": result.mole_fractions[8],
                "vapor_mdea_mole_fraction": result.mole_fractions[9],
                "vapor_co2_mole_fraction": result.mole_fractions[10],
                "maximum_original_row": max(map(abs, result.residuals)),
                "reaction_affinity_max": affinity(result),
                "reference_status": ";".join(sorted({d.status.name for d in result.neutral_reference_diagnostics})),
                "reference_admitted_change": ";".join(
                    f"K{d.reaction_index + 1}={d.observed_terminal_change}"
                    for d in result.neutral_reference_diagnostics),
                "iterations": result.iterations,
            })
        if lane in ("fitted-S0", "fitted-S1"):
            row["matching_k_ij_303_k"] = matching_k_ij(parameters, alpha, observed * 1.0e5, k_ij)
        rows.append(row)
        print(json.dumps(row), flush=True)
    return rows


def lane_summary(rows: list[dict]) -> dict[str, object]:
    solved = [r for r in rows if r["success"]]
    ard = 100.0 * math.fsum(abs(r["relative_error"]) for r in solved) / len(solved) if solved else None
    affinity_max = max((r["reaction_affinity_max"] for r in solved), default=None)
    admitted = all(r["reference_status"] == "Available" for r in solved)
    return {
        "conventions": rows[0]["conventions"], "k_ij_303_k": rows[0]["k_ij_303_k"],
        "pressure_seed_bar": [r["pressure_seed_bar"] for r in rows],
        "row_count": len(rows), "solved_row_count": len(solved),
        "reference_status": sorted({r["reference_status"] for r in solved}),
        "reference_admitted_change": sorted({r["reference_admitted_change"] for r in solved}),
        "maximum_original_row": max((r["maximum_original_row"] for r in solved), default=None),
        "reaction_affinity_max": affinity_max,
        "native_pressure_bar": [r.get("native_pressure_bar") for r in rows],
        "native_ard_percent": ard,
        # F2: 3/3 rows solve, affinity <= 1e-10, references admitted, ARD within 10.25 +/- 3 points (#308).
        "F2": {"all_rows_solve": len(solved) == len(rows),
               "affinity_closed": affinity_max is not None and affinity_max <= 1e-10,
               "references_admitted": admitted,
               "ard_within_tolerance": ard is not None and abs(ard - 10.25) <= 3.0},
        "mdeah_ln_gamma_star": [r.get("mdeah_ln_gamma_star") for r in rows],
        "mdeah_fraction_of_mdea": [r.get("mdeah_fraction_of_mdea") for r in rows],
        "liquid_density_kg_m3": [r.get("liquid_density_kg_m3") for r in rows],
        "matching_k_ij_303_k": [r.get("matching_k_ij_303_k") for r in rows],
    }


def write_outputs(rows: list[dict]) -> None:
    with (OUTPUT / "wangler-2018-table-10-native.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summaries = {lane: lane_summary([r for r in rows if r["lane"] == lane]) for lane in dict.fromkeys(r["lane"] for r in rows)}
    for name in DISPERSION:
        central = summaries[f"fitted-{name}"]
        central["F2"]["passes"] = all(central["F2"].values())
        # F3: Table 10 at k_ij(303 K) -/+ the Figure 10 u(phi) allowance.
        central["F3"] = {"k_ij_303_k_band": [summaries[f"fitted-{name}-{side}"]["k_ij_303_k"] for side in ("minus", "plus")],
                         "native_pressure_bar_band": [summaries[f"fitted-{name}-{side}"]["native_pressure_bar"]
                                                      for side in ("minus", "plus")],
                         "maximum_relative_pressure_change": max(
                             abs(b / c - 1.0) for side in ("minus", "plus")
                             for b, c in zip(summaries[f"fitted-{name}-{side}"]["native_pressure_bar"],
                                             central["native_pressure_bar"], strict=True))}
    wheel = next((ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
    summary = {
        "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "figure_10_summary_sha256": hashlib.sha256(FIGURE_10.read_bytes()).hexdigest(),
        "declaration": {"basis": "RawMoleFractionInfiniteDilution", "solvent": "water",
                        "species_reference": "water pure liquid; every other species infinite dilution in water",
                        "reference_pressure_pa": None, "reaction_constants": "Table 2, Eq. 19, unconverted",
                        "permittivity": "fixed 60 wt% blend 0.4 eps_water + 0.6 * 22.4 on every neutral species",
                        "ions": "MDEAH+ m 3.675; packing 0.88 sigma; Debye-Hueckel diameter sigma",
                        "l_hb": "Table 7 at 303 K",
                        "k_ij_303_k": "printed-S0: Table 6; fitted lanes: Figure 10 fits, Eq. 34 linear "
                                      "0.25 k(273 K) + 0.75 k(313 K), and -/+ its u(phi) allowance"},
        "paper_reported_model_ard_percent": 10.25,
        "working_tolerance_points": 3.0,
        "lanes": summaries,
    }
    (OUTPUT / "wangler-2018-table-10-native-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def declarations():
    from epcsaft import Mixture
    return [(Mixture(build_parameters(TEMPERATURE, 0.0, d)), problem(0.1, 1e5)) for d in DISPERSION.values()]


if __name__ == "__main__":
    main()
