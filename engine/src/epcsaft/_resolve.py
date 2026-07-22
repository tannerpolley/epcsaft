from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import combinations
from typing import Literal

from .parameters import ParameterSet
from .records import (
    AssociationParameterRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
)
from .units import _error


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
        "ionic-held-2008-debye-huckel",
        "ionic-ssm-ds",
        "ionic-ascani-original-born",
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
    dielectric_ion_suppression_coefficient: float | None
    temperature_min_k: float | None
    temperature_max_k: float | None
    ion_mole_fraction_max: float | None


def _missing(*, component_id: str | None = None, family: str, pair: tuple[str, str] | None = None):
    if component_id is not None:
        return _error("missing-record", component_id=component_id, family=family)
    assert pair is not None
    return _error(
        "missing-record",
        component_id_a=pair[0],
        component_id_b=pair[1],
        family=family,
    )


def _magnitude(record: object, unit: str) -> float:
    value = record.constant if isinstance(record, ConstantCorrelation) else record.value  # type: ignore[attr-defined]
    if unit != "dimensionless":
        value = value.to(unit).magnitude
    return float(value)


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
        key = (*sorted((component_a, component_b)), family)
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
        matrix[i * size + j] = matrix[j * size + i] = float(record.value)
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


def _domain_limits(parameters: ParameterSet) -> tuple[float | None, float | None, float | None]:
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
        float(domain.temperature_min.to("kelvin").magnitude)
        for domain in active_domains
        if domain.temperature_min is not None
    )
    upper_bounds = tuple(
        float(domain.temperature_max.to("kelvin").magnitude)
        for domain in active_domains
        if domain.temperature_max is not None
    )
    ion_bounds = tuple(
        float(domain.ion_mole_fraction_max)
        for domain in active_domains
        if domain.ion_mole_fraction_max is not None
    )
    temperature_min = max(lower_bounds) if lower_bounds else None
    temperature_max = min(upper_bounds) if upper_bounds else None
    if (
        temperature_min is not None
        and temperature_max is not None
        and temperature_min > temperature_max
    ):
        raise _error("invalid-domain", field="temperature", value="empty-intersection")
    return temperature_min, temperature_max, min(ion_bounds) if ion_bounds else None


def _resolve_model(parameters: ParameterSet) -> _ResolvedModel:
    if not isinstance(parameters, ParameterSet):
        raise _error("invalid-value", field="parameters", expected="ParameterSet")

    singles = {
        (record.component_id, record.family): record
        for record in parameters.records
        if isinstance(record, SingleParameterRecord)
    }
    correlations = {
        (record.component_id, record.family): record
        for record in parameters.records
        if isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
    }
    pair_records = {
        (record.component_id_a, record.component_id_b, record.family): record
        for record in parameters.records
        if isinstance(record, PairParameterRecord)
    }
    sites = tuple(record for record in parameters.records if isinstance(record, SiteRecord))
    association_records = tuple(
        record for record in parameters.records if isinstance(record, AssociationParameterRecord)
    )
    models = {
        record.family: record
        for record in parameters.records
        if isinstance(record, ModelParameterRecord)
    }

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
    temperature_min, temperature_max, ion_fraction_max = _domain_limits(parameters)

    selected_charges = tuple(
        singles.get((component_id, "charge_number")) for component_id in component_ids
    )
    ionic_topology = any(
        record is not None and int(record.value) != 0 for record in selected_charges
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
            tuple(float(value) for value in k_ij),
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
            None,
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
    charges = tuple(int(record.value) for record in charge_records)  # type: ignore[union-attr]
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

    figiel_schema_present = any(
        family in {"born_diameter", "solvation_factor"} for _, family in singles
    ) or any(
        family
        in {
            "ionic_region_relative_permittivity",
            "dielectric_ion_suppression_coefficient",
        }
        for family in models
    )
    l_ij_present = any(key[2] == "l_ij" for key in pair_records)
    formulation_record = models.get("electrolyte_formulation")
    formulation = None if formulation_record is None else formulation_record.value
    if figiel_schema_present and l_ij_present:
        raise _error(
            "invalid-value",
            family="electrolyte-model",
            reason="mixed-born-dielectric-schema",
        )
    if formulation is not None and (figiel_schema_present or l_ij_present):
        raise _error(
            "invalid-value",
            family="electrolyte-model",
            reason="mixed-formulation-schema",
        )
    if formulation is None and not figiel_schema_present and not l_ij_present:
        raise _error("missing-record", family="electrolyte_formulation")
    if formulation == "held-2008-debye-huckel":
        held_charges = {
            "water": 0,
            "hydronium-cation": 1,
            "hydroxide-anion": -1,
        }
        if set(component_ids) != set(held_charges) or any(
            charge != held_charges[component_id]
            for component_id, charge in zip(component_ids, charges, strict=True)
        ):
            raise _error(
                "invalid-value",
                family="electrolyte_formulation",
                reason="component-charge-contract",
            )
        held_k_ij = tuple(value for value in k_ij if value is not None)
        if len(held_k_ij) != len(k_ij):
            raise _error("missing-record", family="k_ij")
        for left, right in combinations(range(len(component_ids)), 2):
            expected = 1.0 if charges[left] != 0 and charges[right] != 0 else 0.0
            if held_k_ij[left * len(component_ids) + right] != expected:
                raise _error(
                    "invalid-value",
                    family="k_ij",
                    reason="held-dispersion-rule",
                )
        held_relative_permittivity: list[_TemperatureRecord | None] = []
        for index, (component_id, charge) in enumerate(
            zip(component_ids, charges, strict=True)
        ):
            if charge == 0:
                if molar_mass[index] is None:
                    raise _missing(component_id=component_id, family="molar_mass")
                held_relative_permittivity.append(
                    _component_record(
                        component_id,
                        "relative_permittivity",
                        singles,
                        correlations,
                        allow_correlation=True,
                    )
                )
            else:
                held_relative_permittivity.append(None)
        return _ResolvedModel(
            "ionic-held-2008-debye-huckel",
            component_ids,
            molar_mass,
            segment_count,
            segment_diameter,
            dispersion_energy,
            held_k_ij,
            tuple(0.0 for _ in range(len(component_ids) ** 2)),
            site_ids,
            site_multiplicity,
            association_pairs,
            charges,
            tuple(held_relative_permittivity),
            (),
            (),
            None,
            None,
            temperature_min,
            temperature_max,
            ion_fraction_max,
        )
    if not figiel_schema_present:
        l_ij = _pair_matrix(component_ids, pair_records, "l_ij", required=True)
        relative_permittivity = tuple(
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
            "ionic-ascani-original-born",
            component_ids,
            molar_mass,
            segment_count,
            segment_diameter,
            dispersion_energy,
            tuple(float(value) for value in k_ij),
            tuple(float(value) for value in l_ij),
            site_ids,
            site_multiplicity,
            association_pairs,
            charges,
            relative_permittivity,
            (),
            (),
            None,
            None,
            temperature_min,
            temperature_max,
            ion_fraction_max,
        )

    relative_permittivity: list[_TemperatureRecord | None] = []
    born_diameter: list[float | None] = []
    solvation_factor: list[_TemperatureRecord | None] = []
    for index, (component_id, charge) in enumerate(zip(component_ids, charges, strict=True)):
        if charge == 0:
            if molar_mass[index] is None:
                raise _missing(component_id=component_id, family="molar_mass")
            relative_permittivity.append(
                _component_record(
                    component_id,
                    "relative_permittivity",
                    singles,
                    correlations,
                    allow_correlation=True,
                )
            )
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

    required_models = (
        "ionic_region_relative_permittivity",
        "dielectric_ion_suppression_coefficient",
    )
    model_values = []
    for family in required_models:
        record = models.get(family)
        if record is None:
            raise _error("missing-record", family=family)
        value = float(record.value)
        if value <= 0.0:
            raise _error("invalid-value", family=family, value=value)
        model_values.append(value)

    return _ResolvedModel(
        "ionic-ssm-ds",
        component_ids,
        molar_mass,
        segment_count,
        segment_diameter,
        dispersion_energy,
        tuple(float(value) for value in k_ij),
        tuple(0.0 for _ in range(len(component_ids) ** 2)),
        site_ids,
        site_multiplicity,
        association_pairs,
        charges,
        tuple(relative_permittivity),
        tuple(born_diameter),
        tuple(solvation_factor),
        model_values[0],
        model_values[1],
        temperature_min,
        temperature_max,
        ion_fraction_max,
    )


__all__ = ()
