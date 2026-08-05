from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import combinations
from typing import Literal, SupportsFloat, cast

import pint

from .parameters import Parameters
from .records import (
    AssociationParameterRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    PolarFormulation,
    RelativePermittivityFormulation,
    SingleParameterRecord,
    SiteRecord,
)
from .units import ParameterError, _error

_ScalarRecord = SingleParameterRecord | ConstantCorrelation
_TemperatureRecord = (
    SingleParameterRecord | ConstantCorrelation | ConstantPlusSumOfExponentialsCorrelation
)


@dataclass(frozen=True, slots=True)
class _ResolvedAssociationPair:
    endpoint_a: tuple[str, str]
    endpoint_b: tuple[str, str]
    energy_over_k: float
    volume: float


@dataclass(frozen=True, slots=True)
class _ResolvedModel:
    topology: Literal[
        "neutral",
        "associating",
        "ionic-fully-dissociated-debye-huckel",
        "ionic-ssm-ds",
        "ionic-original-born",
    ]
    component_ids: tuple[str, ...]
    molar_mass_kg_per_mol: tuple[float | None, ...]
    segment_count: tuple[float, ...]
    segment_diameter: tuple[_TemperatureRecord, ...]
    dispersion_energy_over_k: tuple[float, ...]
    k_ij: tuple[float, ...]
    l_ij: tuple[float, ...]
    site_ids: tuple[tuple[str, str], ...]
    site_multiplicity: tuple[int, ...]
    association_pairs: tuple[_ResolvedAssociationPair, ...]
    charge_number: tuple[int, ...]
    relative_permittivity: tuple[_TemperatureRecord | None, ...]
    born_diameter_angstrom: tuple[float | None, ...]
    solvation_factor: tuple[_TemperatureRecord | None, ...]
    ionic_region_relative_permittivity: float | None
    ion_fraction_suppression_coefficient: float | None
    temperature_min_k: float | None
    temperature_max_k: float | None
    pressure_min_pa: float | None
    pressure_max_pa: float | None
    ion_mole_fraction_max: float | None
    relative_permittivity_formulation: Literal[
        "none",
        "ion-fraction-suppression",
        "component-permittivity-mixing",
        "density-volume-temperature",
        "ion-specific-suppression",
        "dipolar-density-perturbation",
        "solvent-only",
    ] = "none"
    dielectric_volume_m3_per_mol: tuple[float, ...] = ()
    dielectric_temperature_k: tuple[float, ...] = ()
    ion_specific_suppression_coefficient: tuple[float, ...] = ()
    dipole_scaling_debye2: tuple[float, ...] = ()
    polarizability_scaling_angstrom3: tuple[float, ...] = ()
    correlation_integral_parameter: tuple[float, ...] = ()
    polar_formulation: PolarFormulation | None = None
    dipole_moment_debye: tuple[float, ...] = ()
    quadrupole_moment_debye_angstrom: tuple[float, ...] = ()


def _missing(
    *, component_id: str | None = None, family: str, pair: tuple[str, str] | None = None
) -> ParameterError:
    if component_id is not None:
        return _error("missing-record", component_id=component_id, family=family)
    assert pair is not None
    return _error(
        "missing-record",
        component_id_a=pair[0],
        component_id_b=pair[1],
        family=family,
    )


def _magnitude(
    record: SingleParameterRecord
    | AssociationParameterRecord
    | ConstantCorrelation
    | ConstantPlusSumOfExponentialsCorrelation,
    unit: str,
) -> float:
    value = (
        record.constant
        if isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
        else record.value
    )
    if unit != "dimensionless":
        if not isinstance(value, pint.Quantity):
            raise TypeError(f"expected quantity convertible to {unit}")
        value = value.to(unit).magnitude
    return float(cast(SupportsFloat, value))


def _integer(record: SingleParameterRecord) -> int:
    if type(record.value) is not int:
        raise TypeError("expected integer parameter value")
    return record.value


def _required_matrix(values: tuple[float | None, ...], family: str) -> tuple[float, ...]:
    if any(value is None for value in values):
        raise _error("missing-record", family=family)
    return tuple(cast(float, value) for value in values)


def _domain_magnitude(value: object, unit: str) -> float:
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"expected domain quantity convertible to {unit}")
    return float(value.to(unit).magnitude)


def _component_record(
    component_id: str,
    family: str,
    singles: Mapping[tuple[str, str], SingleParameterRecord],
    correlations: Mapping[tuple[str, str], _TemperatureRecord],
    *,
    allow_correlation: bool,
) -> _TemperatureRecord:
    key = (component_id, family)
    candidates = tuple(
        record
        for record in (
            singles.get(key),
            correlations.get(key) if allow_correlation else None,
        )
        if record is not None
    )
    if not candidates:
        raise _missing(component_id=component_id, family=family)
    if len(candidates) != 1:
        raise _error(
            "invalid-value",
            component_id=component_id,
            family=family,
            reason="duplicate-owner",
        )
    return candidates[0]


def _pair_matrix(
    component_ids: tuple[str, ...],
    pairs: dict[tuple[str, str, str], PairParameterRecord],
    family: str,
    *,
    required: bool,
    not_applicable: frozenset[frozenset[str]] = frozenset(),
) -> tuple[float | None, ...]:
    size = len(component_ids)
    matrix: list[float | None] = [None] * (size * size)
    family_present = any(key[2] == family for key in pairs)
    if required or family_present:
        for index in range(size):
            matrix[index * size + index] = 0.0
    positions = {component_id: index for index, component_id in enumerate(component_ids)}
    for component_a, component_b in combinations(component_ids, 2):
        ordered = sorted((component_a, component_b))
        key = (ordered[0], ordered[1], family)
        record = pairs.get(key)
        if record is None:
            pair = frozenset((component_a, component_b))
            if pair in not_applicable:
                i = positions[component_a]
                j = positions[component_b]
                matrix[i * size + j] = matrix[j * size + i] = 0.0
                continue
            if required or family_present:
                raise _missing(pair=(key[0], key[1]), family=family)
            continue
        i = positions[component_a]
        j = positions[component_b]
        matrix[i * size + j] = matrix[j * size + i] = float(
            cast(SupportsFloat, record.value)
        )
    return tuple(matrix)


def _association(
    sites: tuple[SiteRecord, ...],
    records: tuple[AssociationParameterRecord, ...],
) -> tuple[tuple[tuple[str, str], ...], tuple[int, ...], tuple[_ResolvedAssociationPair, ...]]:
    if not sites:
        if records:
            raise _error("invalid-value", family="association", reason="missing-sites")
        return (), (), ()
    if not records:
        raise _error("missing-record", family="association", reason="missing-site-pairs")

    site_ids = tuple(sorted((record.component_id, record.site_id) for record in sites))
    multiplicity_by_site = {
        (record.component_id, record.site_id): record.multiplicity for record in sites
    }
    grouped: dict[
        tuple[tuple[str, str], tuple[str, str]],
        dict[str, AssociationParameterRecord],
    ] = {}
    for record in records:
        endpoints = (
            (record.component_id_a, record.site_id_a),
            (record.component_id_b, record.site_id_b),
        )
        grouped.setdefault(endpoints, {})[record.family] = record

    resolved = []
    participating: set[tuple[str, str]] = set()
    for endpoints, values in sorted(grouped.items()):
        for family in ("association_energy_over_k", "association_volume"):
            if family not in values:
                raise _error(
                    "missing-record",
                    site_a=endpoints[0],
                    site_b=endpoints[1],
                    family=family,
                )
        resolved.append(
            _ResolvedAssociationPair(
                endpoint_a=endpoints[0],
                endpoint_b=endpoints[1],
                energy_over_k=_magnitude(values["association_energy_over_k"], "kelvin"),
                volume=_magnitude(values["association_volume"], "dimensionless"),
            )
        )
        participating.update(endpoints)
    missing_site = next((site for site in site_ids if site not in participating), None)
    if missing_site is not None:
        raise _error("missing-record", site=missing_site, family="association")
    return (
        site_ids,
        tuple(multiplicity_by_site[site] for site in site_ids),
        tuple(resolved),
    )


def _domain_limits(
    parameters: Parameters,
) -> tuple[float | None, float | None, float | None, float | None, float | None]:
    domain_by_id = {domain.domain_id: domain for domain in parameters.domains}
    active_domains = tuple(
        domain_by_id[domain_id]
        for domain_id in {
            record.domain_id
            for record in parameters.records
            if record.domain_id in domain_by_id and domain_by_id[record.domain_id].kind != "unknown"
        }
    )
    lower_bounds = tuple(
        _domain_magnitude(domain.temperature_min, "kelvin")
        for domain in active_domains
        if domain.temperature_min is not None
    )
    upper_bounds = tuple(
        _domain_magnitude(domain.temperature_max, "kelvin")
        for domain in active_domains
        if domain.temperature_max is not None
    )
    pressure_lower_bounds = tuple(
        _domain_magnitude(domain.pressure_min, "pascal")
        for domain in active_domains
        if domain.pressure_min is not None
    )
    pressure_upper_bounds = tuple(
        _domain_magnitude(domain.pressure_max, "pascal")
        for domain in active_domains
        if domain.pressure_max is not None
    )
    ion_bounds = tuple(
        float(cast(SupportsFloat, domain.ion_mole_fraction_max))
        for domain in active_domains
        if domain.ion_mole_fraction_max is not None
    )
    temperature_min = max(lower_bounds) if lower_bounds else None
    temperature_max = min(upper_bounds) if upper_bounds else None
    pressure_min = max(pressure_lower_bounds) if pressure_lower_bounds else None
    pressure_max = min(pressure_upper_bounds) if pressure_upper_bounds else None
    if (
        temperature_min is not None
        and temperature_max is not None
        and temperature_min > temperature_max
    ):
        raise _error("invalid-domain", field="temperature", value="empty-intersection")
    if pressure_min is not None and pressure_max is not None and pressure_min > pressure_max:
        raise _error("invalid-domain", field="pressure", value="empty-intersection")
    return (
        temperature_min,
        temperature_max,
        pressure_min,
        pressure_max,
        min(ion_bounds) if ion_bounds else None,
    )


def _resolve_model(parameters: Parameters) -> _ResolvedModel:
    if not isinstance(parameters, Parameters):
        raise _error("invalid-value", field="parameters", expected="Parameters")

    singles: dict[tuple[str, str], SingleParameterRecord] = {
        (record.component_id, record.family): record
        for record in parameters.records
        if isinstance(record, SingleParameterRecord)
    }
    correlations: dict[tuple[str, str], _TemperatureRecord] = {
        (record.component_id, record.family): record
        for record in parameters.records
        if isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
    }
    pair_records: dict[tuple[str, str, str], PairParameterRecord] = {
        (record.component_id_a, record.component_id_b, record.family): record
        for record in parameters.records
        if isinstance(record, PairParameterRecord)
    }
    sites = tuple(record for record in parameters.records if isinstance(record, SiteRecord))
    association_records = tuple(
        record for record in parameters.records if isinstance(record, AssociationParameterRecord)
    )
    model_records = tuple(
        record for record in parameters.records if isinstance(record, ModelParameterRecord)
    )
    models: dict[str, ModelParameterRecord] = {
        record.family: record for record in model_records
    }
    polar_records = tuple(
        record for record in model_records if record.family == "polar_formulation"
    )
    if len(polar_records) > 1:
        raise _error("invalid-value", family="polar_formulation", reason="duplicate")
    polar_formulation: PolarFormulation | None = (
        None if not polar_records else polar_records[0].value  # type: ignore[assignment]
    )
    moment_families = {"dipole_moment", "quadrupole_moment"}
    if polar_formulation is None and any(family in moment_families for _, family in singles):
        raise _error(
            "invalid-value",
            family="polar_formulation",
            reason="missing-for-multipole",
        )
    formulation_records = tuple(
        record for record in model_records if record.family == "relative_permittivity_formulation"
    )
    if not formulation_records:
        raise _error("missing-record", family="relative_permittivity_formulation")
    if len(formulation_records) != 1:
        raise _error(
            "invalid-value",
            family="relative_permittivity_formulation",
            reason="duplicate",
        )
    relative_permittivity_formulation: RelativePermittivityFormulation = (
        formulation_records[0].value  # type: ignore[assignment]
    )
    formulation_model_families = {"ion_fraction_suppression_coefficient"}
    formulation_single_families = {
        "dielectric_volume",
        "dielectric_temperature",
        "ion_specific_suppression_coefficient",
        "dipole_scaling",
        "polarizability_scaling",
        "correlation_integral_parameter",
    }
    allowed_model_families = (
        {"ion_fraction_suppression_coefficient"}
        if relative_permittivity_formulation == "ion-fraction-suppression"
        else set()
    )
    allowed_single_families = {
        "density-volume-temperature": {"dielectric_volume", "dielectric_temperature"},
        "ion-specific-suppression": {"ion_specific_suppression_coefficient"},
        "dipolar-density-perturbation": {
            "dipole_scaling",
            "polarizability_scaling",
            "correlation_integral_parameter",
        },
    }.get(relative_permittivity_formulation, set())
    extraneous_model_family = next(
        (
            family
            for family in formulation_model_families
            if family in models and family not in allowed_model_families
        ),
        None,
    )
    extraneous_single_family = next(
        (
            family
            for _, family in singles
            if family in formulation_single_families and family not in allowed_single_families
        ),
        None,
    )
    if (extraneous_family := extraneous_model_family or extraneous_single_family) is not None:
        raise _error(
            "invalid-value",
            family=extraneous_family,
            reason="extraneous-formulation-family",
        )

    recognized = (
        SingleParameterRecord,
        ConstantCorrelation,
        ConstantPlusSumOfExponentialsCorrelation,
        PairParameterRecord,
        SiteRecord,
        AssociationParameterRecord,
        ModelParameterRecord,
    )
    unknown = next(
        (record for record in parameters.records if not isinstance(record, recognized)),
        None,
    )
    if unknown is not None:
        raise _error(
            "invalid-value",
            record_id=getattr(unknown, "record_id", type(unknown).__name__),
        )

    component_ids = parameters.component_ids
    molar_mass = tuple(
        None
        if (record := singles.get((component_id, "molar_mass"))) is None
        else _magnitude(record, "kilogram / mole")
        for component_id in component_ids
    )
    segment_count = tuple(
        _magnitude(
            _component_record(
                component_id,
                "segment_count",
                singles,
                correlations,
                allow_correlation=False,
            ),
            "dimensionless",
        )
        for component_id in component_ids
    )
    segment_diameter = tuple(
        _component_record(
            component_id,
            "segment_diameter",
            singles,
            correlations,
            allow_correlation=True,
        )
        for component_id in component_ids
    )
    dispersion_energy = tuple(
        _magnitude(
            _component_record(
                component_id,
                "dispersion_energy_over_k",
                singles,
                correlations,
                allow_correlation=False,
            ),
            "kelvin",
        )
        for component_id in component_ids
    )
    site_ids, site_multiplicity, association_pairs = _association(sites, association_records)
    temperature_min, temperature_max, pressure_min, pressure_max, ion_fraction_max = _domain_limits(
        parameters
    )

    selected_charges = tuple(
        singles.get((component_id, "charge_number")) for component_id in component_ids
    )
    dipole_moment = tuple(
        0.0
        if (record := singles.get((component_id, "dipole_moment"))) is None
        else _magnitude(record, "debye")
        for component_id in component_ids
    )
    quadrupole_moment = tuple(
        0.0
        if (record := singles.get((component_id, "quadrupole_moment"))) is None
        else _magnitude(record, "debye * angstrom")
        for component_id in component_ids
    )
    for index, component_id in enumerate(component_ids):
        if dipole_moment[index] < 0.0 or quadrupole_moment[index] < 0.0:
            raise _error(
                "invalid-value",
                component_id=component_id,
                family="multipole-moment",
                reason="negative",
            )
        charge_record = selected_charges[index]
        if (
            charge_record is not None
            and _integer(charge_record) != 0
            and (dipole_moment[index] != 0.0 or quadrupole_moment[index] != 0.0)
        ):
            raise _error(
                "invalid-value",
                component_id=component_id,
                family="multipole-moment",
                reason="charged-component",
            )
    ionic_topology = any(
        record is not None and _integer(record) != 0 for record in selected_charges
    ) or any(family == "born_diameter" for _, family in singles)
    if not ionic_topology:
        topology: Literal["neutral", "associating"] = "associating" if site_ids else "neutral"
        k_ij = _pair_matrix(component_ids, pair_records, "k_ij", required=True)
        l_ij = _pair_matrix(component_ids, pair_records, "l_ij", required=False)
        return _ResolvedModel(
            topology,
            component_ids,
            molar_mass,
            segment_count,
            segment_diameter,
            dispersion_energy,
            _required_matrix(k_ij, "k_ij"),
            tuple(0.0 if value is None else float(value) for value in l_ij),
            site_ids,
            site_multiplicity,
            association_pairs,
            (),
            (),
            (),
            (),
            None,
            None,
            temperature_min,
            temperature_max,
            pressure_min,
            pressure_max,
            None,
            "none",
            polar_formulation=polar_formulation,
            dipole_moment_debye=dipole_moment,
            quadrupole_moment_debye_angstrom=quadrupole_moment,
        )

    charge_records = tuple(
        _component_record(
            component_id,
            "charge_number",
            singles,
            correlations,
            allow_correlation=False,
        )
        for component_id in component_ids
    )
    charges = tuple(_integer(cast(SingleParameterRecord, record)) for record in charge_records)
    if not any(charge > 0 for charge in charges) or not any(charge < 0 for charge in charges):
        raise _error("invalid-value", family="charge_number", reason="missing-opposite-charges")

    same_sign_ion_pairs = frozenset(
        frozenset((component_ids[left], component_ids[right]))
        for left, right in combinations(range(len(component_ids)), 2)
        if charges[left] * charges[right] > 0
    )
    k_ij = _pair_matrix(
        component_ids,
        pair_records,
        "k_ij",
        required=True,
        not_applicable=same_sign_ion_pairs,
    )

    corrected_born_schema_present = any(
        family in {"born_diameter", "solvation_factor"} for _, family in singles
    ) or any(
        family
        in {
            "ionic_region_relative_permittivity",
            "ion_fraction_suppression_coefficient",
        }
        for family in models
    )
    l_ij_present = any(key[2] == "l_ij" for key in pair_records)
    formulation_record = models.get("electrolyte_formulation")
    formulation = None if formulation_record is None else formulation_record.value
    if corrected_born_schema_present and l_ij_present:
        raise _error(
            "invalid-value",
            family="electrolyte-model",
            reason="mixed-born-dielectric-schema",
        )
    if formulation is not None and (corrected_born_schema_present or l_ij_present):
        raise _error(
            "invalid-value",
            family="electrolyte-model",
            reason="mixed-formulation-schema",
        )
    if formulation is None and not corrected_born_schema_present and not l_ij_present:
        raise _error("missing-record", family="electrolyte_formulation")
    if formulation == "fully-dissociated-debye-huckel":
        if relative_permittivity_formulation != "solvent-only":
            raise _error(
                "invalid-value",
                family="relative_permittivity_formulation",
                reason="incompatible-electrolyte-topology",
            )
        neutral_components = tuple(
            component_id
            for component_id, charge in zip(component_ids, charges, strict=True)
            if charge == 0
        )
        if len(neutral_components) != 1:
            raise _error(
                "invalid-value",
                family="electrolyte_formulation",
                reason="exactly-one-neutral-solvent-required",
            )
        fully_dissociated_k_ij = tuple(value for value in k_ij if value is not None)
        if len(fully_dissociated_k_ij) != len(k_ij):
            raise _error("missing-record", family="k_ij")
        for left, right in combinations(range(len(component_ids)), 2):
            expected = 1.0 if charges[left] != 0 and charges[right] != 0 else 0.0
            if fully_dissociated_k_ij[left * len(component_ids) + right] != expected:
                raise _error(
                    "invalid-value",
                    family="k_ij",
                    reason="fully-dissociated-dispersion-rule",
                )
        solvent_relative_permittivity: list[_TemperatureRecord | None] = []
        for index, (component_id, charge) in enumerate(zip(component_ids, charges, strict=True)):
            if charge == 0:
                if molar_mass[index] is None:
                    raise _missing(component_id=component_id, family="molar_mass")
                solvent_relative_permittivity.append(
                    _component_record(
                        component_id,
                        "relative_permittivity",
                        singles,
                        correlations,
                        allow_correlation=True,
                    )
                )
            else:
                solvent_relative_permittivity.append(None)
        return _ResolvedModel(
            "ionic-fully-dissociated-debye-huckel",
            component_ids,
            molar_mass,
            segment_count,
            segment_diameter,
            dispersion_energy,
            fully_dissociated_k_ij,
            tuple(0.0 for _ in range(len(component_ids) ** 2)),
            site_ids,
            site_multiplicity,
            association_pairs,
            charges,
            tuple(solvent_relative_permittivity),
            (),
            (),
            None,
            None,
            temperature_min,
            temperature_max,
            pressure_min,
            pressure_max,
            ion_fraction_max,
            "solvent-only",
            polar_formulation=polar_formulation,
            dipole_moment_debye=dipole_moment,
            quadrupole_moment_debye_angstrom=quadrupole_moment,
        )
    if relative_permittivity_formulation == "component-permittivity-mixing":
        if corrected_born_schema_present:
            raise _error(
                "invalid-value",
                family="relative_permittivity_formulation",
                reason="mixed-formulation-schema",
            )
        l_ij = _pair_matrix(component_ids, pair_records, "l_ij", required=True)
        original_relative_permittivity = tuple(
            _component_record(
                component_id,
                "relative_permittivity",
                singles,
                correlations,
                allow_correlation=True,
            )
            for component_id in component_ids
        )
        return _ResolvedModel(
            "ionic-original-born",
            component_ids,
            molar_mass,
            segment_count,
            segment_diameter,
            dispersion_energy,
            _required_matrix(k_ij, "k_ij"),
            _required_matrix(l_ij, "l_ij"),
            site_ids,
            site_multiplicity,
            association_pairs,
            charges,
            original_relative_permittivity,
            (),
            (),
            None,
            None,
            temperature_min,
            temperature_max,
            pressure_min,
            pressure_max,
            ion_fraction_max,
            "component-permittivity-mixing",
            polar_formulation=polar_formulation,
            dipole_moment_debye=dipole_moment,
            quadrupole_moment_debye_angstrom=quadrupole_moment,
        )

    if relative_permittivity_formulation not in {
        "ion-fraction-suppression",
        "density-volume-temperature",
        "ion-specific-suppression",
        "dipolar-density-perturbation",
    }:
        raise _error(
            "invalid-value",
            family="relative_permittivity_formulation",
            reason="incompatible-electrolyte-topology",
        )

    relative_permittivity: list[_TemperatureRecord | None] = []
    born_diameter: list[float | None] = []
    solvation_factor: list[_TemperatureRecord | None] = []
    for index, (component_id, charge) in enumerate(zip(component_ids, charges, strict=True)):
        if charge == 0:
            if molar_mass[index] is None:
                raise _missing(component_id=component_id, family="molar_mass")
            if relative_permittivity_formulation in {
                "ion-fraction-suppression",
                "ion-specific-suppression",
            }:
                relative_permittivity.append(
                    _component_record(
                        component_id,
                        "relative_permittivity",
                        singles,
                        correlations,
                        allow_correlation=True,
                    )
                )
            else:
                relative_permittivity.append(None)
            born_diameter.append(None)
        else:
            relative_permittivity.append(None)
            born = _component_record(
                component_id,
                "born_diameter",
                singles,
                correlations,
                allow_correlation=False,
            )
            born_diameter.append(_magnitude(born, "angstrom"))
        solvation_factor.append(
            _component_record(
                component_id,
                "solvation_factor",
                singles,
                correlations,
                allow_correlation=True,
            )
        )

    required_models = ["ionic_region_relative_permittivity"]
    if relative_permittivity_formulation == "ion-fraction-suppression":
        required_models.append("ion_fraction_suppression_coefficient")
    model_values = []
    for family in required_models:
        model_record = models.get(family)
        if model_record is None:
            raise _error("missing-record", family=family)
        value = float(cast(SupportsFloat, model_record.value))
        if value <= 0.0:
            raise _error("invalid-value", family=family, value=value)
        model_values.append(value)

    dielectric_suppression = (
        model_values[1] if relative_permittivity_formulation == "ion-fraction-suppression" else None
    )

    def component_values(family: str, unit: str) -> tuple[float, ...]:
        values = []
        for component_id in component_ids:
            record = singles.get((component_id, family))
            values.append(0.0 if record is None else _magnitude(record, unit))
        return tuple(values)

    dielectric_volume = component_values("dielectric_volume", "meter ** 3 / mole")
    dielectric_temperature = component_values("dielectric_temperature", "kelvin")
    ion_specific_suppression = component_values(
        "ion_specific_suppression_coefficient", "dimensionless"
    )
    dipole_scaling = component_values("dipole_scaling", "debye ** 2")
    polarizability_scaling = component_values("polarizability_scaling", "angstrom ** 3")
    correlation_integral = component_values("correlation_integral_parameter", "dimensionless")

    for index, charge in enumerate(charges):
        if relative_permittivity_formulation == "density-volume-temperature" and charge == 0:
            if dielectric_volume[index] <= 0.0 or dielectric_temperature[index] <= 0.0:
                raise _error(
                    "missing-record",
                    component_id=component_ids[index],
                    family="density-permittivity-parameters",
                )
        if relative_permittivity_formulation == "ion-specific-suppression" and charge != 0:
            if ion_specific_suppression[index] <= 0.0:
                raise _error(
                    "missing-record",
                    component_id=component_ids[index],
                    family="ion_specific_suppression_coefficient",
                )
        if relative_permittivity_formulation == "dipolar-density-perturbation":
            if dipole_scaling[index] < 0.0 or polarizability_scaling[index] < 0.0:
                raise _error(
                    "invalid-value",
                    component_id=component_ids[index],
                    family="dipolar-permittivity-parameters",
                )
            if correlation_integral[index] <= 0.0:
                raise _error(
                    "missing-record",
                    component_id=component_ids[index],
                    family="correlation_integral_parameter",
                )
    if (
        relative_permittivity_formulation == "ion-specific-suppression"
        and sum(charge == 0 for charge in charges) != 1
    ):
        raise _error(
            "invalid-value",
            family="relative_permittivity_formulation",
            reason="mixed-solvent-volume-data-required",
        )

    return _ResolvedModel(
        "ionic-ssm-ds",
        component_ids,
        molar_mass,
        segment_count,
        segment_diameter,
        dispersion_energy,
        _required_matrix(k_ij, "k_ij"),
        tuple(0.0 for _ in range(len(component_ids) ** 2)),
        site_ids,
        site_multiplicity,
        association_pairs,
        charges,
        tuple(relative_permittivity),
        tuple(born_diameter),
        tuple(solvation_factor),
        model_values[0],
        dielectric_suppression,
        temperature_min,
        temperature_max,
        pressure_min,
        pressure_max,
        ion_fraction_max,
        relative_permittivity_formulation,
        dielectric_volume,
        dielectric_temperature,
        ion_specific_suppression,
        dipole_scaling,
        polarizability_scaling,
        correlation_integral,
        polar_formulation,
        dipole_moment,
        quadrupole_moment,
    )


__all__ = ()
