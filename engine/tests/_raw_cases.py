from __future__ import annotations

from copy import deepcopy


NEUTRAL_BINARY = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("methane", "ethane"),
    "parameters": {
        "mw": (0.016043, 0.03007),
        "m": (1.0, 1.6069),
        "s": (3.7039, 3.5206),
        "e": (150.03, 191.42),
        "k_ij": ((0.0, 0.0), (0.0, 0.0)),
    },
    "options": {"permittivity_model": "none"},
    "validity": {"kind": "unknown"},
}


ASSOCIATING_PURE = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("associating-fluid",),
    "parameters": {
        "m": (1.0,),
        "s": (3.0,),
        "e": (200.0,),
        "k_ij": ((0.0,),),
        "sites": (
            {
                "component_id": "associating-fluid",
                "site_id": "a",
                "site_class": "donor",
                "multiplicity": 1,
            },
            {
                "component_id": "associating-fluid",
                "site_id": "b",
                "site_class": "acceptor",
                "multiplicity": 1,
            },
        ),
        "association": (
            {
                "component_id_a": "associating-fluid",
                "site_id_a": "a",
                "component_id_b": "associating-fluid",
                "site_id_b": "b",
                "association_energy_over_k": 2400.0,
                "association_volume": 0.03,
            },
        ),
    },
    "options": {"relative_permittivity_formulation": "none"},
    "validity": {"kind": "unknown"},
}


AQUEOUS_ELECTROLYTE = {
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
        "k_ij": (
            (0.0, -0.3, -0.3),
            (-0.3, 0.0, 0.8),
            (-0.3, 0.8, 0.0),
        ),
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
        "pressure_min_pa": 100000.0,
        "pressure_max_pa": 100000.0,
        "ion_mole_fraction_max": 0.38,
    },
}


FULLY_DISSOCIATED_ELECTROLYTE = deepcopy(AQUEOUS_ELECTROLYTE)
for _field in ("d_born", "f_solv"):
    del FULLY_DISSOCIATED_ELECTROLYTE["parameters"][_field]
FULLY_DISSOCIATED_ELECTROLYTE["parameters"]["k_ij"] = (
    (0.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, 1.0, 0.0),
)
FULLY_DISSOCIATED_ELECTROLYTE["options"] = {
    "electrolyte": "fully-dissociated-debye-huckel",
    "permittivity_model": "solvent-only",
}


ORIGINAL_BORN_ELECTROLYTE = deepcopy(AQUEOUS_ELECTROLYTE)
for _field in ("d_born", "f_solv"):
    del ORIGINAL_BORN_ELECTROLYTE["parameters"][_field]
ORIGINAL_BORN_ELECTROLYTE["parameters"].update(
    epsilon_r=(None, 8.0, 8.0),
    l_ij=((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
)
ORIGINAL_BORN_ELECTROLYTE["options"] = {"permittivity_model": "component-permittivity-mixing"}


DIPOLAR_PERMITTIVITY_ELECTROLYTE = deepcopy(AQUEOUS_ELECTROLYTE)
DIPOLAR_PERMITTIVITY_ELECTROLYTE["parameters"].update(
    mu2=(1.0, 0.0, 0.0),
    alpha=(1.0, 0.0, 0.0),
    a2=(1.0, 1.0, 1.0),
)
DIPOLAR_PERMITTIVITY_ELECTROLYTE["options"] = {
    "epsilon_r_ion": 8.0,
    "permittivity_model": "dipolar-density-perturbation",
}


RAW_CASES = {
    "neutral-binary": NEUTRAL_BINARY,
    "associating-pure": ASSOCIATING_PURE,
    "corrected-born-electrolyte": AQUEOUS_ELECTROLYTE,
    "fully-dissociated-electrolyte": FULLY_DISSOCIATED_ELECTROLYTE,
    "original-born-electrolyte": ORIGINAL_BORN_ELECTROLYTE,
    "dipolar-permittivity-electrolyte": DIPOLAR_PERMITTIVITY_ELECTROLYTE,
}


_NEUTRAL_COMPONENTS = {
    "methane": (0.016043, 1.0, 3.7039, 150.03),
    "ethane": (0.03007, 1.6069, 3.5206, 191.42),
    "ethanol": (0.04607, 2.3827, 3.1771, 198.24),
    "propane": (0.044097, 2.0020, 3.6184, 208.11),
}


def neutral_pure(component_id: str) -> dict[str, object]:
    """Return one compact raw neutral case without retaining a model globally."""
    mw, segments, diameter, energy = _NEUTRAL_COMPONENTS.get(
        component_id, _NEUTRAL_COMPONENTS["methane"]
    )
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": (component_id,),
        "parameters": {
            "mw": (mw,),
            "m": (segments,),
            "s": (diameter,),
            "e": (energy,),
            "k_ij": ((0.0,),),
        },
        "options": {"permittivity_model": "none"},
        "validity": {"kind": "unknown"},
    }


def neutral_binary(component_ids: tuple[str, str]) -> dict[str, object]:
    values = [_NEUTRAL_COMPONENTS[component_id] for component_id in component_ids]
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": component_ids,
        "parameters": {
            "mw": tuple(value[0] for value in values),
            "m": tuple(value[1] for value in values),
            "s": tuple(value[2] for value in values),
            "e": tuple(value[3] for value in values),
            "k_ij": ((0.0, 0.0), (0.0, 0.0)),
        },
        "options": {"permittivity_model": "none"},
        "validity": {"kind": "unknown"},
    }


def electrolyte(component_ids: tuple[str, ...]) -> dict[str, object]:
    """Return the compact electrolyte case reordered and relabeled by charge role."""
    result = deepcopy(AQUEOUS_ELECTROLYTE)
    role = {
        component_id: (1 if "cation" in component_id else -1 if "anion" in component_id else 0)
        for component_id in component_ids
    }
    base_index = {0: 0, 1: 1, -1: 2}
    indices = tuple(base_index[role[component_id]] for component_id in component_ids)
    result["components"] = component_ids
    parameters = result["parameters"]
    for key in ("mw", "m", "s", "e", "z", "d_born", "f_solv"):
        source = parameters[key]
        parameters[key] = tuple(source[index] for index in indices)
    source_matrix = parameters["k_ij"]
    parameters["k_ij"] = tuple(
        tuple(source_matrix[left][right] for right in indices) for left in indices
    )
    solvent = next(component_id for component_id in component_ids if role[component_id] == 0)
    for site in parameters["sites"]:
        site["component_id"] = solvent
    for pair in parameters["association"]:
        pair["component_id_a"] = solvent
        pair["component_id_b"] = solvent
    for correlation in parameters["correlations"]:
        correlation["component_id"] = solvent
    return result


__all__ = (
    "AQUEOUS_ELECTROLYTE",
    "ASSOCIATING_PURE",
    "DIPOLAR_PERMITTIVITY_ELECTROLYTE",
    "FULLY_DISSOCIATED_ELECTROLYTE",
    "NEUTRAL_BINARY",
    "ORIGINAL_BORN_ELECTROLYTE",
    "RAW_CASES",
    "electrolyte",
    "neutral_binary",
    "neutral_pure",
)
