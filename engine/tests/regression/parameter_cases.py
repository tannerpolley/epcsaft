from __future__ import annotations

from copy import deepcopy

_NEUTRAL_COMPONENTS = {
    "methane": (0.016043, 1.0, 3.7039, 150.03),
    "ethane": (0.03007, 1.6069, 3.5206, 191.42),
    "propane": (0.044096, 2.0020, 3.6184, 208.11),
}


def neutral_parameters(component_ids: tuple[str, ...]) -> dict[str, object]:
    values = tuple(_NEUTRAL_COMPONENTS[component_id] for component_id in component_ids)
    size = len(component_ids)
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": component_ids,
        "parameters": {
            "mw": tuple(value[0] for value in values),
            "m": tuple(value[1] for value in values),
            "s": tuple(value[2] for value in values),
            "e": tuple(value[3] for value in values),
            "k_ij": tuple(tuple(0.0 for _ in range(size)) for _ in range(size)),
        },
        "options": {"permittivity_model": "none"},
        "validity": {
            "kind": "reported-conditions",
            "temperature_min_k": 100.0,
            "temperature_max_k": 500.0,
            "pressure_min_pa": 1.0,
            "pressure_max_pa": 100_000_000.0,
        },
    }


def associating_parameters(
    *,
    component_id: str = "associating-fluid",
    sites: tuple[tuple[str, int], ...] = (("a", 1), ("b", 1)),
    pairs: tuple[tuple[str, str, float, float], ...] = (("a", "b", 2400.0, 0.03),),
    segment_count: float = 1.0,
    segment_diameter_angstrom: float = 3.0,
    dispersion_energy_over_k: float = 200.0,
) -> dict[str, object]:
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": (component_id,),
        "parameters": {
            "mw": (0.088,),
            "m": (segment_count,),
            "s": (segment_diameter_angstrom,),
            "e": (dispersion_energy_over_k,),
            "k_ij": ((0.0,),),
            "sites": tuple(
                {
                    "component_id": component_id,
                    "site_id": site_id,
                    "site_class": site_id,
                    "multiplicity": multiplicity,
                }
                for site_id, multiplicity in sites
            ),
            "association": tuple(
                {
                    "component_id_a": component_id,
                    "site_id_a": left,
                    "component_id_b": component_id,
                    "site_id_b": right,
                    "association_energy_over_k": energy,
                    "association_volume": volume,
                }
                for left, right, energy, volume in pairs
            ),
        },
        "options": {"relative_permittivity_formulation": "none"},
        "validity": {
            "kind": "reported-conditions",
            "temperature_min_k": 250.0,
            "temperature_max_k": 450.0,
            "pressure_min_pa": 1.0,
            "pressure_max_pa": 10_000_000.0,
        },
    }


_AQUEOUS_BASE = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("water", "sodium-cation", "chloride-anion"),
    "parameters": {
        "mw": (0.0180153, None, None),
        "m": (1.2047, 1.0, 1.0),
        "s": (None, 2.8232, 2.7560),
        "e": (353.95, 230.0, 170.0),
        "z": (0, 1, -1),
        "d_born": (None, 3.445, 4.100),
        "f_solv": (1.5, 1.0, 1.0),
        "k_ij": ((0.0, -0.3, -0.3), (-0.3, 0.0, 0.8), (-0.3, 0.8, 0.0)),
        "sites": (
            {
                "component_id": "water",
                "site_id": "a",
                "site_class": "donor",
                "multiplicity": 1,
            },
            {
                "component_id": "water",
                "site_id": "b",
                "site_class": "acceptor",
                "multiplicity": 1,
            },
        ),
        "association": (
            {
                "component_id_a": "water",
                "site_id_a": "a",
                "component_id_b": "water",
                "site_id_b": "b",
                "association_energy_over_k": 2425.7,
                "association_volume": 0.04509,
            },
        ),
        "correlations": (
            {
                "component_id": "water",
                "family": "segment_diameter",
                "form": "constant-plus-sum-of-exponentials",
                "constant": 2.7927,
                "terms": (
                    {"amplitude": 10.11, "exponent_coefficient_per_k": -0.01775},
                    {"amplitude": -1.417, "exponent_coefficient_per_k": -0.01146},
                ),
            },
            {
                "component_id": "water",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 78.09,
            },
        ),
    },
    "options": {
        "epsilon_r_ion": 8.0,
        "a_dh": 7.01,
        "permittivity_model": "ion-fraction-suppression",
    },
    "validity": {
        "kind": "reported-conditions",
        "temperature_min_k": 298.15,
        "temperature_max_k": 298.15,
        "pressure_min_pa": 100_000.0,
        "pressure_max_pa": 100_000.0,
        "ion_mole_fraction_max": 0.38,
    },
}


def aqueous_parameters(component_ids: tuple[str, str, str]) -> dict[str, object]:
    result = deepcopy(_AQUEOUS_BASE)
    result["components"] = component_ids
    solvent, cation, anion = component_ids
    parameters = result["parameters"]
    assert isinstance(parameters, dict)
    parameters["sites"] = tuple({**site, "component_id": solvent} for site in parameters["sites"])
    parameters["association"] = tuple(
        {
            **association,
            "component_id_a": solvent,
            "component_id_b": solvent,
        }
        for association in parameters["association"]
    )
    parameters["correlations"] = tuple(
        {**correlation, "component_id": solvent} for correlation in parameters["correlations"]
    )
    # The generic sentinels exercise structural role handling; source-specific
    # campaign values remain in ePC-SAFT-data and Validation.
    parameters["z"] = (0, 1, -1)
    result["components"] = (solvent, cation, anion)
    return result
