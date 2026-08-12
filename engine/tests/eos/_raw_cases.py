from __future__ import annotations

from copy import deepcopy

_SOURCE_ID = "engine-numerical-fixture"
_DOMAIN_ID = "engine-numerical-fixture-domain"
_PROVENANCE = {
    "source_id": _SOURCE_ID,
    "locator": "compact numerical fixture values",
    "domain_id": _DOMAIN_ID,
}
_UNITS = {
    "molar_mass": "kilogram / mole",
    "charge_number": "elementary-charge",
    "dipole_moment": "debye",
    "quadrupole_moment": "debye * angstrom",
    "segment_count": "dimensionless",
    "segment_diameter": "angstrom",
    "dispersion_energy_over_k": "kelvin",
    "relative_permittivity": "dimensionless",
    "packing_diameter": "angstrom",
    "debye_huckel_diameter": "angstrom",
    "born_diameter": "angstrom",
    "solvation_factor": "dimensionless",
    "dipole_scaling": "debye ** 2",
    "polarizability_scaling": "angstrom ** 3",
    "correlation_integral_parameter": "dimensionless",
    "k_ij": "dimensionless",
    "l_ij": "dimensionless",
    "ionic_region_relative_permittivity": "dimensionless",
    "ion_fraction_suppression_coefficient": "dimensionless",
}
_FIXED = frozenset({"molar_mass", "charge_number", "dipole_moment", "quadrupole_moment"})


def _value(magnitude: float, family: str) -> dict[str, object]:
    return {"magnitude": magnitude, "unit": _UNITS[family]}


def _component(component_id: str, **families: float) -> dict[str, object]:
    fixed = {}
    coefficients = []
    for family, magnitude in families.items():
        if family in _FIXED:
            fixed[family] = {
                "value": _value(magnitude, family),
                "provenance": deepcopy(_PROVENANCE),
            }
        else:
            coefficients.append(
                {
                    "identity": f"component/{component_id}/{family}",
                    "family": family,
                    "value": _value(magnitude, family),
                    "provenance": deepcopy(_PROVENANCE),
                }
            )
    return {
        "component_id": component_id,
        "fixed": fixed,
        "coefficients": sorted(coefficients, key=lambda item: str(item["family"])),
    }


def _pair(
    component_a: str,
    component_b: str,
    **families: float,
) -> dict[str, object]:
    component_a, component_b = sorted((component_a, component_b))
    return {
        "component_id_a": component_a,
        "component_id_b": component_b,
        "coefficients": [
            {
                "identity": f"pair/{component_a}/{component_b}/{family}",
                "family": family,
                "value": _value(magnitude, family),
                "provenance": deepcopy(_PROVENANCE),
            }
            for family, magnitude in sorted(families.items())
        ],
    }


def _model_families(
    *,
    association: str = "none",
    electrolyte: str = "none",
    permittivity: str = "none",
    polar: str = "none",
    born_corrections: tuple[float, float] | None = None,
) -> list[dict[str, object]]:
    choices = {
        "association": association,
        "base": "pc-saft",
        "electrolyte": electrolyte,
        "permittivity": permittivity,
        "polar": polar,
    }
    family_ids = {
        "association": "model/association",
        "base": "model/base",
        "electrolyte": "model/electrolyte",
        "permittivity": "model/relative_permittivity",
        "polar": "model/polar",
    }
    result = [
        {
            "family_id": family_ids[kind],
            "kind": kind,
            "choice": choice,
            "provenance": deepcopy(_PROVENANCE),
        }
        for kind, choice in choices.items()
    ]
    if electrolyte == "born":
        if born_corrections is None:
            raise ValueError("Born fixtures require explicit correction coefficients")
        electrolyte_family = next(item for item in result if item["kind"] == "electrolyte")
        electrolyte_family["c_shell"], electrolyte_family["c_dielectric"] = born_corrections
    elif born_corrections is not None:
        raise ValueError("Born correction coefficients require the Born model")
    return result


def _document(
    document_id: str,
    components: list[dict[str, object]],
    *,
    pairs: list[dict[str, object]] | None = None,
    model_families: list[dict[str, object]] | None = None,
    model_coefficients: list[dict[str, object]] | None = None,
    correlations: list[dict[str, object]] | None = None,
    topology: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "document_id": document_id,
        "document_version": 1,
        "purpose": "package-test-fixture",
        "sources": [
            {
                "source_id": _SOURCE_ID,
                "citation": "Project-authored compact numerical EOS fixture.",
                "use_basis": "Numerical invariants and independent equation oracles.",
            }
        ],
        "domains": [{"domain_id": _DOMAIN_ID, "kind": "unknown"}],
        "components": components,
        "pairs": [] if pairs is None else pairs,
        "model_families": _model_families() if model_families is None else model_families,
        "model_coefficients": [] if model_coefficients is None else model_coefficients,
        "correlations": [] if correlations is None else correlations,
        "topology": {"presets": [], "sites": [], "edges": []} if topology is None else topology,
    }


def _model_coefficient(family: str, magnitude: float) -> dict[str, object]:
    return {
        "identity": f"model/relative_permittivity/{family}",
        "family": family,
        "value": _value(magnitude, family),
        "provenance": deepcopy(_PROVENANCE),
    }


def _constant_correlation(component_id: str, family: str, magnitude: float) -> dict[str, object]:
    correlation_id = f"component/{component_id}/{family}/constant"
    return {
        "correlation_id": correlation_id,
        "component_id": component_id,
        "family": family,
        "form": "constant",
        "independent_variables": [],
        "constant": {
            "identity": f"{correlation_id}/constant",
            "value": _value(magnitude, family),
        },
        "provenance": deepcopy(_PROVENANCE),
    }


def _water_diameter_correlation() -> dict[str, object]:
    correlation_id = "component/water/segment_diameter/constant-plus-sum-of-exponentials"
    return {
        "correlation_id": correlation_id,
        "component_id": "water",
        "family": "segment_diameter",
        "form": "constant-plus-sum-of-exponentials",
        "independent_variables": ["temperature"],
        "constant": {
            "identity": f"{correlation_id}/constant",
            "value": _value(2.7927, "segment_diameter"),
        },
        "terms": [
            {
                "amplitude": {
                    "identity": f"{correlation_id}/term-0/amplitude",
                    "value": _value(10.11, "segment_diameter"),
                },
                "exponent_coefficient": {
                    "identity": f"{correlation_id}/term-0/exponent_coefficient",
                    "value": {"magnitude": -0.01775, "unit": "1 / kelvin"},
                },
            },
            {
                "amplitude": {
                    "identity": f"{correlation_id}/term-1/amplitude",
                    "value": _value(-1.417, "segment_diameter"),
                },
                "exponent_coefficient": {
                    "identity": f"{correlation_id}/term-1/exponent_coefficient",
                    "value": {"magnitude": -0.01146, "unit": "1 / kelvin"},
                },
            },
        ],
        "provenance": deepcopy(_PROVENANCE),
    }


def _association_topology(
    component_id: str,
    energy: float,
    volume: float,
) -> dict[str, object]:
    prefix = f"association/{component_id}/a/{component_id}/b"
    return {
        "presets": [],
        "sites": [
            {
                "component_id": component_id,
                "site_id": "a",
                "site_role": "donor",
                "multiplicity": 1,
                "provenance": deepcopy(_PROVENANCE),
            },
            {
                "component_id": component_id,
                "site_id": "b",
                "site_role": "acceptor",
                "multiplicity": 1,
                "provenance": deepcopy(_PROVENANCE),
            },
        ],
        "edges": [
            {
                "endpoint_a": {"component_id": component_id, "site_id": "a"},
                "endpoint_b": {"component_id": component_id, "site_id": "b"},
                "energy_over_k": {
                    "identity": f"{prefix}/energy_over_k",
                    "value": {"magnitude": energy, "unit": "kelvin"},
                },
                "volume": {
                    "identity": f"{prefix}/volume",
                    "value": {"magnitude": volume, "unit": "dimensionless"},
                },
                "source": {
                    "kind": "explicit",
                    "provenance": [deepcopy(_PROVENANCE)],
                },
            }
        ],
    }


NEUTRAL_BINARY = _document(
    "neutral-binary-fixture",
    [
        _component(
            "methane",
            molar_mass=0.016043,
            segment_count=1.0,
            segment_diameter=3.7039,
            dispersion_energy_over_k=150.03,
        ),
        _component(
            "ethane",
            molar_mass=0.03007,
            segment_count=1.6069,
            segment_diameter=3.5206,
            dispersion_energy_over_k=191.42,
        ),
    ],
    pairs=[_pair("methane", "ethane", k_ij=0.0)],
)

NEUTRAL_TERNARY_IMMISCIBLE = _document(
    "neutral-ternary-immiscible-fixture",
    [
        _component(
            component_id,
            molar_mass=0.030,
            segment_count=1.5,
            segment_diameter=segment_diameter,
            dispersion_energy_over_k=300.0,
        )
        for component_id, segment_diameter in zip(
            ("alpha", "beta", "gamma"),
            (3.5, 3.7, 3.9),
            strict=True,
        )
    ],
    pairs=[
        _pair(left, right, k_ij=0.20)
        for left, right in (
            ("alpha", "beta"),
            ("alpha", "gamma"),
            ("beta", "gamma"),
        )
    ],
)

CHARGED_TERNARY_IMMISCIBLE = _document(
    "charged-ternary-immiscible-fixture",
    [
        _component(
            component_id,
            molar_mass=0.030,
            charge_number=0,
            segment_count=1.5,
            segment_diameter=segment_diameter,
            dispersion_energy_over_k=300.0,
            solvation_factor=1.2,
        )
        for component_id, segment_diameter in zip(
            ("alpha", "beta", "gamma"),
            (3.5, 3.7, 3.9),
            strict=True,
        )
    ]
    + [
        _component(
            "cation",
            charge_number=1,
            segment_count=1.0,
            segment_diameter=2.8,
            dispersion_energy_over_k=230.0,
            packing_diameter=2.464,
            debye_huckel_diameter=2.464,
            born_diameter=3.4,
            solvation_factor=1.0,
        ),
        _component(
            "anion",
            charge_number=-1,
            segment_count=1.0,
            segment_diameter=2.8,
            dispersion_energy_over_k=170.0,
            packing_diameter=2.464,
            debye_huckel_diameter=2.464,
            born_diameter=4.1,
            solvation_factor=1.0,
        ),
    ],
    pairs=[
        _pair(
            left,
            right,
            k_ij=(
                0.20
                if left in {"alpha", "beta", "gamma"}
                and right in {"alpha", "beta", "gamma"}
                else 0.0
            ),
        )
        for index, left in enumerate(("alpha", "beta", "gamma", "cation", "anion"))
        for right in ("alpha", "beta", "gamma", "cation", "anion")[index + 1 :]
    ],
    model_families=_model_families(
        association="general-site",
        electrolyte="born",
        permittivity="ion-fraction-suppression",
        born_corrections=(1.0, 1.0),
    ),
    model_coefficients=[
        _model_coefficient("ionic_region_relative_permittivity", 8.0),
        _model_coefficient("ion_fraction_suppression_coefficient", 7.01),
    ],
    correlations=[
        _constant_correlation(component_id, "relative_permittivity", permittivity)
        for component_id, permittivity in (
            ("alpha", 78.0),
            ("beta", 55.0),
            ("gamma", 35.0),
        )
    ],
    topology=_association_topology("alpha", 2200.0, 0.03),
)

POLAR_BINARY = _document(
    "polar-binary-fixture",
    [
        _component(
            "dipolar",
            segment_count=1.5,
            segment_diameter=3.4,
            dispersion_energy_over_k=220.0,
            dipole_moment=1.7,
        ),
        _component(
            "quadrupolar",
            segment_count=2.25,
            segment_diameter=3.8,
            dispersion_energy_over_k=280.0,
            quadrupole_moment=4.1,
        ),
    ],
    pairs=[_pair("dipolar", "quadrupolar", k_ij=0.08, l_ij=0.17)],
    model_families=_model_families(polar="point-multipole"),
)

POLAR_PURE = _document(
    "polar-pure-fixture",
    [
        _component(
            "polar-fluid",
            segment_count=1.5,
            segment_diameter=3.4,
            dispersion_energy_over_k=220.0,
            dipole_moment=1.7,
            quadrupole_moment=4.1,
        )
    ],
    model_families=_model_families(polar="point-multipole"),
)

ASSOCIATING_PURE = _document(
    "associating-pure-fixture",
    [
        _component(
            "associating-fluid",
            segment_count=1.0,
            segment_diameter=3.0,
            dispersion_energy_over_k=200.0,
        )
    ],
    model_families=_model_families(association="general-site"),
    topology=_association_topology("associating-fluid", 2400.0, 0.03),
)

NEUTRAL_ASSOCIATING_MIXTURE = deepcopy(ASSOCIATING_PURE)
NEUTRAL_ASSOCIATING_MIXTURE["document_id"] = "neutral-associating-mixture-fixture"
_neutral_associating_components = NEUTRAL_ASSOCIATING_MIXTURE["components"]
_neutral_associating_pairs = NEUTRAL_ASSOCIATING_MIXTURE["pairs"]
assert isinstance(_neutral_associating_components, list)
assert isinstance(_neutral_associating_pairs, list)
_neutral_associating_components.append(
    _component(
        "inert-fluid",
        segment_count=1.4,
        segment_diameter=3.4,
        dispersion_energy_over_k=180.0,
    )
)
_neutral_associating_pairs.append(
    _pair("associating-fluid", "inert-fluid", k_ij=0.04)
)

AQUEOUS_ELECTROLYTE = _document(
    "aqueous-electrolyte-fixture",
    [
        _component(
            "water",
            molar_mass=0.0180153,
            charge_number=0,
            segment_count=1.2047,
            dispersion_energy_over_k=353.95,
            solvation_factor=1.5,
        ),
        _component(
            "sodium-cation",
            charge_number=1,
            segment_count=1.0,
            segment_diameter=2.8232,
            dispersion_energy_over_k=230.0,
            packing_diameter=2.484416,
            debye_huckel_diameter=2.484416,
            born_diameter=3.445,
            solvation_factor=1.0,
        ),
        _component(
            "chloride-anion",
            charge_number=-1,
            segment_count=1.0,
            segment_diameter=2.7560,
            dispersion_energy_over_k=170.0,
            packing_diameter=2.42528,
            debye_huckel_diameter=2.42528,
            born_diameter=4.100,
            solvation_factor=1.0,
        ),
    ],
    pairs=[
        _pair("water", "sodium-cation", k_ij=-0.3),
        _pair("water", "chloride-anion", k_ij=-0.3),
        _pair("sodium-cation", "chloride-anion", k_ij=0.8),
    ],
    model_families=_model_families(
        association="general-site",
        electrolyte="born",
        permittivity="ion-fraction-suppression",
        born_corrections=(1.0, 1.0),
    ),
    model_coefficients=[
        _model_coefficient("ionic_region_relative_permittivity", 8.0),
        _model_coefficient("ion_fraction_suppression_coefficient", 7.01),
    ],
    correlations=[
        _water_diameter_correlation(),
        _constant_correlation("water", "relative_permittivity", 78.09),
    ],
    topology=_association_topology("water", 2425.7, 0.04509),
)
AQUEOUS_ELECTROLYTE["domains"] = [
    {
        "domain_id": _DOMAIN_ID,
        "kind": "reported-conditions",
        "temperature_min": {"magnitude": 298.15, "unit": "kelvin"},
        "temperature_max": {"magnitude": 298.15, "unit": "kelvin"},
        "pressure_min": {"magnitude": 100000.0, "unit": "pascal"},
        "pressure_max": {"magnitude": 100000.0, "unit": "pascal"},
        "ion_mole_fraction_max": {"magnitude": 0.38, "unit": "dimensionless"},
    }
]

MULTINEUTRAL_REACTIVE = _document(
    "multineutral-reactive-fixture",
    [
        _component(
            "solvent",
            molar_mass=0.018,
            charge_number=0,
            segment_count=1.2,
            segment_diameter=2.8,
            dispersion_energy_over_k=350.0,
            solvation_factor=1.5,
        ),
        _component(
            "cosolvent",
            molar_mass=0.030,
            charge_number=0,
            segment_count=1.1,
            segment_diameter=3.0,
            dispersion_energy_over_k=280.0,
            solvation_factor=1.2,
        ),
        _component(
            "cation",
            charge_number=1,
            segment_count=1.0,
            segment_diameter=2.8,
            dispersion_energy_over_k=230.0,
            packing_diameter=2.464,
            debye_huckel_diameter=2.464,
            born_diameter=3.4,
            solvation_factor=1.0,
        ),
        _component(
            "anion",
            charge_number=-1,
            segment_count=1.0,
            segment_diameter=2.8,
            dispersion_energy_over_k=170.0,
            packing_diameter=2.464,
            debye_huckel_diameter=2.464,
            born_diameter=4.1,
            solvation_factor=1.0,
        ),
    ],
    pairs=[
        _pair(component_a, component_b, k_ij=0.0)
        for index, component_a in enumerate(("solvent", "cosolvent", "cation", "anion"))
        for component_b in ("solvent", "cosolvent", "cation", "anion")[index + 1 :]
    ],
    model_families=_model_families(
        association="general-site",
        electrolyte="born",
        permittivity="ion-fraction-suppression",
        born_corrections=(1.0, 1.0),
    ),
    model_coefficients=[
        _model_coefficient("ionic_region_relative_permittivity", 8.0),
        _model_coefficient("ion_fraction_suppression_coefficient", 7.01),
    ],
    correlations=[
        _constant_correlation("solvent", "relative_permittivity", 78.0),
        _constant_correlation("cosolvent", "relative_permittivity", 40.0),
    ],
    topology=_association_topology("solvent", 2425.7, 0.04509),
)
MULTINEUTRAL_REACTIVE["domains"] = [
    {
        "domain_id": _DOMAIN_ID,
        "kind": "reported-conditions",
        "temperature_min": {"magnitude": 298.0, "unit": "kelvin"},
        "temperature_max": {"magnitude": 350.0, "unit": "kelvin"},
        "pressure_min": {"magnitude": 1.0, "unit": "pascal"},
        "pressure_max": {"magnitude": 100_000_000.0, "unit": "pascal"},
    }
]

MULTISOLVENT_MULTIION = deepcopy(MULTINEUTRAL_REACTIVE)
MULTISOLVENT_MULTIION["document_id"] = "multisolvent-multiion-fixture"
next(
    family
    for family in MULTISOLVENT_MULTIION["model_families"]
    if family["kind"] == "permittivity"
)["choice"] = "component-permittivity-mixing"
MULTISOLVENT_MULTIION["model_coefficients"] = [
    coefficient
    for coefficient in MULTISOLVENT_MULTIION["model_coefficients"]
    if coefficient["family"] != "ion_fraction_suppression_coefficient"
]
MULTISOLVENT_MULTIION["correlations"] = [
    correlation
    for correlation in MULTISOLVENT_MULTIION["correlations"]
    if correlation["family"] != "relative_permittivity"
]
MULTISOLVENT_MULTIION["components"].extend(
    [
        _component(
            "divalent-cation",
            charge_number=2,
            segment_count=1.0,
            segment_diameter=3.0,
            dispersion_energy_over_k=210.0,
            packing_diameter=2.64,
            debye_huckel_diameter=2.64,
            born_diameter=3.2,
            solvation_factor=1.1,
        ),
        _component(
            "iodide-anion",
            charge_number=-1,
            segment_count=1.0,
            segment_diameter=3.2,
            dispersion_energy_over_k=160.0,
            packing_diameter=2.816,
            debye_huckel_diameter=2.816,
            born_diameter=4.4,
            solvation_factor=0.95,
        ),
    ]
)
for _component_id, _permittivity in (
    ("solvent", 78.0),
    ("cosolvent", 40.0),
    ("cation", 8.0),
    ("anion", 8.0),
    ("divalent-cation", 8.0),
    ("iodide-anion", 8.0),
):
    _record = next(
        component
        for component in MULTISOLVENT_MULTIION["components"]
        if component["component_id"] == _component_id
    )
    _record["coefficients"].extend(
        _component(_component_id, relative_permittivity=_permittivity)["coefficients"]
    )
    _record["coefficients"].sort(key=lambda item: str(item["family"]))
MULTISOLVENT_MULTIION["pairs"].extend(
    [
        _pair("divalent-cation", component_id, k_ij=0.0)
        for component_id in ("solvent", "cosolvent", "anion", "iodide-anion")
    ]
    + [
        _pair("iodide-anion", component_id, k_ij=0.0)
        for component_id in ("solvent", "cosolvent", "cation")
    ]
)

DIPOLAR_PERMITTIVITY_ELECTROLYTE = deepcopy(AQUEOUS_ELECTROLYTE)
DIPOLAR_PERMITTIVITY_ELECTROLYTE["document_id"] = "dipolar-permittivity-electrolyte-fixture"
DIPOLAR_PERMITTIVITY_ELECTROLYTE["correlations"] = [
    correlation
    for correlation in DIPOLAR_PERMITTIVITY_ELECTROLYTE["correlations"]
    if correlation["family"] != "relative_permittivity"
]
for _family, _magnitudes in (
    ("dipole_scaling", (1.0, 0.0, 0.0)),
    ("polarizability_scaling", (1.0, 0.0, 0.0)),
    ("correlation_integral_parameter", (1.0, 1.0, 1.0)),
):
    for _component_record, _magnitude in zip(
        DIPOLAR_PERMITTIVITY_ELECTROLYTE["components"], _magnitudes, strict=True
    ):
        _component_id = _component_record["component_id"]
        _component_record["coefficients"].append(
            {
                "identity": f"component/{_component_id}/{_family}",
                "family": _family,
                "value": _value(_magnitude, _family),
                "provenance": deepcopy(_PROVENANCE),
            }
        )
DIPOLAR_PERMITTIVITY_ELECTROLYTE["model_coefficients"] = [
    _model_coefficient("ionic_region_relative_permittivity", 8.0)
]
next(
    item
    for item in DIPOLAR_PERMITTIVITY_ELECTROLYTE["model_families"]
    if item["kind"] == "permittivity"
)["choice"] = "dipolar-density-perturbation"

_NEUTRAL_COMPONENTS = {
    "methane": (0.016043, 1.0, 3.7039, 150.03),
    "ethane": (0.03007, 1.6069, 3.5206, 191.42),
    "ethanol": (0.04607, 2.3827, 3.1771, 198.24),
    "propane": (0.044097, 2.0020, 3.6184, 208.11),
}


def neutral_pure(component_id: str) -> dict[str, object]:
    """Return one strict neutral Parameters fixture without retaining a model."""
    mw, segments, diameter, energy = _NEUTRAL_COMPONENTS.get(
        component_id, _NEUTRAL_COMPONENTS["methane"]
    )
    return _document(
        f"neutral-{component_id}-fixture",
        [
            _component(
                component_id,
                molar_mass=mw,
                segment_count=segments,
                segment_diameter=diameter,
                dispersion_energy_over_k=energy,
            )
        ],
    )


def neutral_binary(component_ids: tuple[str, str]) -> dict[str, object]:
    values = [_NEUTRAL_COMPONENTS[component_id] for component_id in component_ids]
    return _document(
        "neutral-" + "-".join(component_ids) + "-fixture",
        [
            _component(
                component_id,
                molar_mass=value[0],
                segment_count=value[1],
                segment_diameter=value[2],
                dispersion_energy_over_k=value[3],
            )
            for component_id, value in zip(component_ids, values, strict=True)
        ],
        pairs=[_pair(*component_ids, k_ij=0.0)],
    )


def component_values(mapping: dict[str, object], family: str) -> tuple[float | None, ...]:
    values = []
    for component in mapping["components"]:
        record = component["fixed"].get(family)
        if record is None:
            record = next(
                (item for item in component["coefficients"] if item["family"] == family),
                None,
            )
        values.append(None if record is None else float(record["value"]["magnitude"]))
    return tuple(values)


def set_component_values(
    mapping: dict[str, object], family: str, values: tuple[float | None, ...]
) -> None:
    for component, magnitude in zip(mapping["components"], values, strict=True):
        component["fixed"].pop(family, None)
        component["coefficients"] = [
            item for item in component["coefficients"] if item["family"] != family
        ]
        if magnitude is None:
            continue
        component_id = component["component_id"]
        if family in _FIXED:
            component["fixed"][family] = {
                "value": _value(magnitude, family),
                "provenance": deepcopy(_PROVENANCE),
            }
        else:
            component["coefficients"].append(
                {
                    "identity": f"component/{component_id}/{family}",
                    "family": family,
                    "value": _value(magnitude, family),
                    "provenance": deepcopy(_PROVENANCE),
                }
            )


def pair_matrix(mapping: dict[str, object], family: str) -> tuple[tuple[float, ...], ...]:
    component_ids = tuple(item["component_id"] for item in mapping["components"])
    indices = {component_id: index for index, component_id in enumerate(component_ids)}
    matrix = [[0.0 for _ in component_ids] for _ in component_ids]
    for pair in mapping["pairs"]:
        coefficient = next(
            (item for item in pair["coefficients"] if item["family"] == family), None
        )
        if coefficient is None:
            continue
        left = indices[pair["component_id_a"]]
        right = indices[pair["component_id_b"]]
        matrix[left][right] = matrix[right][left] = float(coefficient["value"]["magnitude"])
    return tuple(tuple(row) for row in matrix)


def set_binary_pair_value(mapping: dict[str, object], family: str, magnitude: float) -> None:
    pair = mapping["pairs"][0]
    pair["coefficients"] = [item for item in pair["coefficients"] if item["family"] != family]
    component_a = pair["component_id_a"]
    component_b = pair["component_id_b"]
    pair["coefficients"].append(
        {
            "identity": f"pair/{component_a}/{component_b}/{family}",
            "family": family,
            "value": _value(magnitude, family),
            "provenance": deepcopy(_PROVENANCE),
        }
    )


__all__ = (
    "AQUEOUS_ELECTROLYTE",
    "ASSOCIATING_PURE",
    "DIPOLAR_PERMITTIVITY_ELECTROLYTE",
    "MULTINEUTRAL_REACTIVE",
    "MULTISOLVENT_MULTIION",
    "NEUTRAL_ASSOCIATING_MIXTURE",
    "NEUTRAL_BINARY",
    "NEUTRAL_TERNARY_IMMISCIBLE",
    "POLAR_BINARY",
    "POLAR_PURE",
    "component_values",
    "neutral_binary",
    "neutral_pure",
    "pair_matrix",
    "set_binary_pair_value",
    "set_component_values",
)
