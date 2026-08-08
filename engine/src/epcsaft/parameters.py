from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from typing import Literal, SupportsFloat, cast

from ._canonical import _document_fingerprint, _selected_document_fingerprint
from ._dictionary_io import _records_from_mapping
from ._document_io import _load_json
from .records import (
    AssociationTopology,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    CrossAssociationCorrectionRecord,
    ModelFamilySpec,
    ModelParameterRecord,
    ParameterRecord,
    ParameterSpec,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)
from .units import _error


@dataclass(frozen=True, slots=True)
class ActiveParameterDomain:
    """The intrinsic open/closed domain of one active coordinate value."""

    lower: float | None
    upper: float | None
    lower_closed: bool
    upper_closed: bool

    def contains(self, value: float) -> bool:
        return (
            self.lower is None or value > self.lower or self.lower_closed and value == self.lower
        ) and (
            self.upper is None or value < self.upper or self.upper_closed and value == self.upper
        )


@dataclass(frozen=True, slots=True)
class ActiveParameterApplicability:
    """Resolved physical application of one derivative-ready coordinate."""

    scope: Literal[
        "component",
        "unordered-component-pair",
        "association-edge",
        "model",
    ]
    component_ids: tuple[str, ...]
    component_indices: tuple[int, ...]
    model_domain: Literal["fixed-topology-mixture"]
    observation_family: Literal["fixed-state-helmholtz"]
    application: Literal[
        "reciprocal-cross-dispersion",
        "reciprocal-cross-size",
        "role-aware-cross-association-energy",
        "component-property",
        "typed-correlation-coefficient",
        "symmetric-association-property",
        "auxiliary-model-coefficient",
    ]
    association_role_pair: tuple[str, str] | None = None
    association_edges: tuple[tuple[int, int], ...] = ()
    association_energy_scales: tuple[float, ...] = ()


@dataclass(frozen=True, slots=True)
class ActiveParameter:
    """One exact continuous coordinate selected from ``Parameters``."""

    identity: str
    family: str
    value: float
    unit: str
    source_domain: ValidityDomain | None
    domain: ActiveParameterDomain | None
    applicability: ActiveParameterApplicability | None


@dataclass(frozen=True, slots=True, init=False)
class ActiveParameterSet:
    """An immutable ordered selection of canonical parameter identities."""

    parameters: Parameters
    coordinates: tuple[ActiveParameter, ...]
    fingerprint: str

    def __init__(
        self,
        parameters: Parameters,
        identities: Sequence[str],
        *,
        values: Sequence[object] | None = None,
    ) -> None:
        if not isinstance(parameters, Parameters):
            raise _error("invalid-value", field="parameters", expected="Parameters")
        if isinstance(identities, str):
            raise _error("invalid-selection", field="identities", value="string")
        requested = tuple(identities)
        if not requested or any(type(identity) is not str for identity in requested):
            raise _error("invalid-selection", field="identities", value=requested)
        if len(requested) != len(set(requested)):
            raise _error("duplicate-record", field="parameter_identity")
        available = {spec.identity: spec for spec in parameters.parameter_specs}
        selected = tuple(_active_spec(parameters, available, identity) for identity in requested)
        missing = tuple(
            identity for identity, spec in zip(requested, selected, strict=True) if spec is None
        )
        if missing:
            raise _error("missing-record", field="parameter_identity", value=missing)
        selected = cast(tuple[ParameterSpec, ...], selected)
        resolved_identities = tuple(spec.identity for spec in selected)
        if len(resolved_identities) != len(set(resolved_identities)):
            raise _error("duplicate-record", field="parameter_identity")
        discrete = tuple(spec.identity for spec in selected if not spec.continuous)
        if discrete:
            raise _error("invalid-selection", field="parameter_identity", value=discrete)
        if values is not None:
            if isinstance(values, str):
                raise _error("invalid-selection", field="parameter_values", value="string")
            supplied = tuple(values)
            if len(supplied) != len(selected):
                raise _error(
                    "invalid-selection",
                    field="parameter_values",
                    expected=len(selected),
                    value=len(supplied),
                )
        else:
            supplied = tuple(spec.value for spec in selected)
        coordinates = tuple(
            _active_coordinate(parameters, spec, value, explicit=values is not None)
            for spec, value in zip(selected, supplied, strict=True)
        )
        payload = json.dumps(
            {
                "schema": "epcsaft.active-parameter-set",
                "parameters": parameters.fingerprint,
                "identities": resolved_identities,
                "values": tuple(coordinate.value for coordinate in coordinates),
            },
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "fingerprint", f"sha256:{hashlib.sha256(payload).hexdigest()}")

    @property
    def identities(self) -> tuple[str, ...]:
        return tuple(coordinate.identity for coordinate in self.coordinates)

    @property
    def values(self) -> tuple[float, ...]:
        return tuple(coordinate.value for coordinate in self.coordinates)


def _record_components(record: ParameterRecord) -> tuple[str, ...]:
    if isinstance(record, ModelParameterRecord):
        return ()
    if isinstance(
        record,
        (
            SingleParameterRecord,
            SiteRecord,
            ConstantCorrelation,
            ConstantPlusSumOfExponentialsCorrelation,
        ),
    ):
        return (record.component_id,)
    return (record.component_id_a, record.component_id_b)


def _identity_components(identity: str) -> tuple[str, ...]:
    parts = identity.split("/")
    if parts[0] == "component" and len(parts) >= 2:
        return (parts[1],)
    if parts[0] == "pair" and len(parts) >= 3:
        return (parts[1], parts[2])
    if parts[0] == "association" and len(parts) >= 5:
        return (parts[1], parts[3])
    return ()


_K_IJ_DOMAIN = ActiveParameterDomain(None, 1.0, False, True)
_L_IJ_DOMAIN = ActiveParameterDomain(None, 1.0, False, False)
_K_HB_IJ_DOMAIN = ActiveParameterDomain(None, 1.0, False, False)
_POSITIVE_DOMAIN = ActiveParameterDomain(0.0, None, False, False)


def _active_spec(
    parameters: Parameters,
    available: Mapping[str, ParameterSpec],
    identity: str,
) -> ParameterSpec | None:
    if identity in available:
        return available[identity]
    parts = identity.split("/")
    if len(parts) != 5 or parts[0] != "pair" or parts[3:] != ["association", "k_hb_ij"]:
        return None
    if parts[1] >= parts[2]:
        raise _error(
            "invalid-selection",
            field="parameter_identity",
            value=identity,
            reason="noncanonical-pair-order",
        )
    component_pair = (parts[1], parts[2])
    site_roles = {
        (site.component_id, site.site_id): site.site_role for site in parameters.topology.sites
    }
    eligible_role_pairs: set[tuple[str, str]] = set()
    for edge in parameters.topology.edges:
        edge_components = tuple(sorted((edge.endpoint_a[0], edge.endpoint_b[0])))
        if edge.source.kind != "combining-rule" or edge_components != component_pair:
            continue
        roles = sorted((site_roles[edge.endpoint_a], site_roles[edge.endpoint_b]))
        eligible_role_pairs.add((roles[0], roles[1]))
    if len(eligible_role_pairs) > 1:
        raise _error(
            "invalid-selection",
            field="parameter_identity",
            value=identity,
            reason="ambiguous-association-role-pair",
        )
    if not eligible_role_pairs:
        return None
    role_a, role_b = next(iter(eligible_role_pairs))
    resolved = available.get(f"pair/{parts[1]}/{parts[2]}/association/{role_a}/{role_b}/k_hb_ij")
    return resolved if resolved is not None and resolved.family == "k_hb_ij" else None


def _active_value(spec: ParameterSpec, value: object) -> float:
    try:
        magnitude = (
            float(value.to(spec.unit).magnitude)
            if hasattr(value, "to")
            else float(cast(SupportsFloat, value))
        )
    except (AttributeError, TypeError, ValueError) as exc:
        raise _error(
            "invalid-active-coordinate",
            identity=spec.identity,
            reason="non-numeric",
        ) from exc
    if not math.isfinite(magnitude):
        raise _error(
            "invalid-active-coordinate",
            identity=spec.identity,
            reason="nonfinite",
        )
    return magnitude


def _active_coordinate(
    parameters: Parameters,
    spec: ParameterSpec,
    value: object,
    *,
    explicit: bool,
) -> ActiveParameter:
    magnitude = _active_value(spec, value)
    source_domain = (
        None
        if spec.domain_id is None
        else next(domain for domain in parameters.domains if domain.domain_id == spec.domain_id)
    )
    parts = spec.identity.split("/")
    if spec.family == "relative_permittivity":
        component_id = parts[1]
        charge = next(
            record.value
            for record in parameters.records
            if isinstance(record, SingleParameterRecord)
            and record.component_id == component_id
            and record.family == "charge_number"
        )
        formulation = next(
            family.choice for family in parameters.model_families if family.kind == "permittivity"
        )
        if formulation != "component-permittivity-mixing" and not (
            charge == 0
            and formulation
            in {
                "ion-fraction-suppression",
                "ion-specific-suppression",
                "solvent-only",
            }
        ):
            raise _error(
                "unsupported-active-coordinate",
                identity=spec.identity,
                reason="relative-permittivity-application-unresolved",
            )
    domain: ActiveParameterDomain | None = None
    applicability: ActiveParameterApplicability | None = None
    component_ids: tuple[str, ...]
    if (
        spec.family in {"k_ij", "l_ij"}
        and len(parts) == 4
        and parts[0] == "pair"
        and parts[3] == spec.family
    ):
        component_ids = (parts[1], parts[2])
        if component_ids[0] >= component_ids[1]:
            raise _error(
                "invalid-active-coordinate",
                identity=spec.identity,
                reason="noncanonical-pair-order",
            )
        component_indices = (
            parameters.component_ids.index(component_ids[0]),
            parameters.component_ids.index(component_ids[1]),
        )
        domain = _K_IJ_DOMAIN if spec.family == "k_ij" else _L_IJ_DOMAIN
        applicability = ActiveParameterApplicability(
            "unordered-component-pair",
            component_ids,
            component_indices,
            "fixed-topology-mixture",
            "fixed-state-helmholtz",
            "reciprocal-cross-dispersion" if spec.family == "k_ij" else "reciprocal-cross-size",
        )
    elif (
        spec.family == "k_hb_ij"
        and len(parts) == 7
        and parts[0] == "pair"
        and parts[3] == "association"
        and parts[6] == "k_hb_ij"
    ):
        component_ids = (parts[1], parts[2])
        role_pair = (parts[4], parts[5])
        if component_ids[0] >= component_ids[1] or role_pair[0] > role_pair[1]:
            raise _error(
                "invalid-active-coordinate",
                identity=spec.identity,
                reason="noncanonical-pair-order",
            )
        correction = next(
            (
                record
                for record in parameters.records
                if isinstance(record, CrossAssociationCorrectionRecord)
                and record.identity == spec.identity
            ),
            None,
        )
        if correction is None:
            raise _error(
                "invalid-active-coordinate",
                identity=spec.identity,
                reason="missing-cross-association-correction",
            )
        site_positions = {
            endpoint: index
            for index, endpoint in enumerate(
                sorted((site.component_id, site.site_id) for site in parameters.topology.sites)
            )
        }
        affected = tuple(
            edge
            for edge in parameters.topology.edges
            if spec.identity in edge.source.input_identities
        )
        if not affected:
            raise _error(
                "invalid-active-coordinate",
                identity=spec.identity,
                reason="missing-eligible-association-edge",
            )
        association_edges = tuple(
            tuple(sorted((site_positions[edge.endpoint_a], site_positions[edge.endpoint_b])))
            for edge in affected
        )
        association_energy_scales = tuple(
            -float(
                cast(
                    SupportsFloat,
                    next(
                        derivative.factor
                        for derivative in edge.source.derivative_sources
                        if derivative.output_identity == edge.energy_identity
                        and derivative.input_identity == spec.identity
                    ),
                )
            )
            for edge in affected
        )
        if any(scale <= 0.0 for scale in association_energy_scales):
            raise _error(
                "invalid-active-coordinate",
                identity=spec.identity,
                reason="invalid-combining-rule-derivative",
            )
        component_indices = (
            parameters.component_ids.index(component_ids[0]),
            parameters.component_ids.index(component_ids[1]),
        )
        domain = _K_HB_IJ_DOMAIN
        applicability = ActiveParameterApplicability(
            "unordered-component-pair",
            component_ids,
            component_indices,
            "fixed-topology-mixture",
            "fixed-state-helmholtz",
            "role-aware-cross-association-energy",
            role_pair,
            cast(tuple[tuple[int, int], ...], association_edges),
            association_energy_scales,
        )
    elif (
        len(parts) == 3
        and parts[0] == "component"
        and parts[2]
        in {
            "segment_count",
            "segment_diameter",
            "dispersion_energy_over_k",
            "relative_permittivity",
            "born_diameter",
            "solvation_factor",
        }
        and spec.family == parts[2]
    ):
        component_ids = (parts[1],)
        applicability = ActiveParameterApplicability(
            "component",
            component_ids,
            (parameters.component_ids.index(parts[1]),),
            "fixed-topology-mixture",
            "fixed-state-helmholtz",
            "component-property",
        )
        domain = _POSITIVE_DOMAIN
    elif (
        len(parts) in {5, 6}
        and parts[0] == "component"
        and (
            parts[2] == "segment_diameter"
            or parts[2] in {"relative_permittivity", "solvation_factor"}
            and len(parts) == 5
        )
    ):
        component_ids = (parts[1],)
        application_resolved = (
            len(parts) == 5
            and parts[4] == "constant"
            or len(parts) == 6
            and parts[4].startswith("term-")
            and parts[5] in {"amplitude", "exponent_coefficient"}
        )
        if not application_resolved:
            raise _error(
                "unsupported-active-coordinate",
                identity=spec.identity,
                reason="correlation-application-unresolved",
            )
        applicability = ActiveParameterApplicability(
            "component",
            component_ids,
            (parameters.component_ids.index(parts[1]),),
            "fixed-topology-mixture",
            "fixed-state-helmholtz",
            "typed-correlation-coefficient",
        )
        if len(parts) == 5:
            domain = _POSITIVE_DOMAIN
    elif (
        len(parts) == 6
        and parts[0] == "association"
        and parts[5] in {"energy_over_k", "volume"}
        and spec.family == f"association_{parts[5]}"
    ):
        endpoints = ((parts[1], parts[2]), (parts[3], parts[4]))
        site_positions = {
            endpoint: index
            for index, endpoint in enumerate(
                sorted((site.component_id, site.site_id) for site in parameters.topology.sites)
            )
        }
        if endpoints[0] not in site_positions or endpoints[1] not in site_positions:
            raise _error(
                "unsupported-active-coordinate",
                identity=spec.identity,
                reason="association-application-unresolved",
            )
        applicability = ActiveParameterApplicability(
            "association-edge",
            tuple(sorted({parts[1], parts[3]})),
            tuple(sorted((site_positions[endpoints[0]], site_positions[endpoints[1]]))),
            "fixed-topology-mixture",
            "fixed-state-helmholtz",
            "symmetric-association-property",
        )
        domain = _POSITIVE_DOMAIN
    elif (
        spec.identity == "model/relative_permittivity/ion_fraction_suppression_coefficient"
        and spec.family == "ion_fraction_suppression_coefficient"
    ) or (
        spec.identity == "model/relative_permittivity/ionic_region_relative_permittivity"
        and spec.family == "ionic_region_relative_permittivity"
    ):
        applicability = ActiveParameterApplicability(
            "model",
            (),
            (),
            "fixed-topology-mixture",
            "fixed-state-helmholtz",
            "auxiliary-model-coefficient",
        )
        domain = _POSITIVE_DOMAIN
    elif explicit:
        raise _error(
            "unsupported-active-coordinate",
            identity=spec.identity,
            reason="trial-value-application-unresolved",
        )
    if domain is not None and not domain.contains(magnitude):
        raise _error(
            "invalid-active-coordinate",
            identity=spec.identity,
            reason="out-of-domain",
            value=magnitude,
        )
    return ActiveParameter(
        spec.identity,
        spec.family,
        magnitude,
        spec.unit,
        source_domain,
        domain,
        applicability,
    )


@dataclass(frozen=True, slots=True, init=False)
class _ParameterDocument:
    """Private immutable owner for one canonical Parameters document."""

    document_id: str
    document_version: int
    purpose: str
    sources: tuple[SourceRecord, ...]
    domains: tuple[ValidityDomain, ...]
    components: tuple[ComponentRecord, ...]
    records: tuple[ParameterRecord, ...]
    model_families: tuple[ModelFamilySpec, ...]
    topology: AssociationTopology
    parameter_specs: tuple[ParameterSpec, ...]
    fingerprint: str
    _mapping_json: bytes

    def __init__(self) -> None:
        raise TypeError("use _ParameterDocument.from_mapping or _ParameterDocument.from_json")

    @classmethod
    def from_json(cls, path: str | PathLike[str]) -> _ParameterDocument:
        return cls.from_mapping(_load_json(path))

    @classmethod
    def from_mapping(cls, mapping: Mapping[str, object]) -> _ParameterDocument:
        values = _records_from_mapping(mapping)
        instance = object.__new__(cls)
        object.__setattr__(instance, "document_id", values.document_id)
        object.__setattr__(instance, "document_version", values.document_version)
        object.__setattr__(instance, "purpose", values.purpose)
        object.__setattr__(
            instance,
            "sources",
            tuple(sorted(values.sources, key=lambda item: item.source_id)),
        )
        object.__setattr__(
            instance,
            "domains",
            tuple(sorted(values.domains, key=lambda item: item.domain_id)),
        )
        object.__setattr__(instance, "components", values.components)
        records: tuple[ParameterRecord, ...] = (
            *values.singles,
            *values.pairs,
            *values.cross_association_corrections,
            *values.sites,
            *values.associations,
            *values.correlations,
            *values.models,
        )
        object.__setattr__(
            instance,
            "records",
            tuple(sorted(records, key=lambda item: item.record_id)),
        )
        object.__setattr__(instance, "model_families", values.model_families)
        object.__setattr__(instance, "topology", values.topology)
        object.__setattr__(instance, "parameter_specs", values.parameter_specs)
        object.__setattr__(instance, "fingerprint", _document_fingerprint(mapping))
        object.__setattr__(
            instance,
            "_mapping_json",
            json.dumps(mapping, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode(
                "utf-8"
            ),
        )
        return instance

    def select(self, component_ids: Sequence[str]) -> Parameters:
        if isinstance(component_ids, str):
            raise _error("invalid-selection", field="component_ids", value="string")
        requested = tuple(component_ids)
        if (
            not requested
            or any(type(component_id) is not str for component_id in requested)
            or len(requested) != len(set(requested))
        ):
            raise _error("invalid-selection", field="component_ids", value=requested)
        component_map = {record.component_id: record for record in self.components}
        if any(component_id not in component_map for component_id in requested):
            raise _error("invalid-selection", field="component_ids", value=requested)
        selected = set(requested)
        records = tuple(
            record
            for record in self.records
            if all(component_id in selected for component_id in _record_components(record))
        )
        specs = tuple(
            spec
            for spec in self.parameter_specs
            if all(component_id in selected for component_id in _identity_components(spec.identity))
        )
        topology = AssociationTopology(
            tuple(item for item in self.topology.presets if item[0] in selected),
            tuple(site for site in self.topology.sites if site.component_id in selected),
            tuple(
                edge
                for edge in self.topology.edges
                if edge.endpoint_a[0] in selected and edge.endpoint_b[0] in selected
            ),
        )
        return Parameters._create(
            document=self,
            component_ids=requested,
            components=tuple(component_map[component_id] for component_id in requested),
            records=records,
            parameter_specs=specs,
            topology=topology,
        )


@dataclass(frozen=True, slots=True, init=False)
class Parameters:
    component_ids: tuple[str, ...]
    components: tuple[ComponentRecord, ...]
    domains: tuple[ValidityDomain, ...]
    records: tuple[ParameterRecord, ...]
    model_families: tuple[ModelFamilySpec, ...]
    topology: AssociationTopology
    parameter_specs: tuple[ParameterSpec, ...]
    fingerprint: str
    _document: _ParameterDocument

    def __init__(self) -> None:
        raise TypeError("use Parameters.from_mapping or Parameters.from_json")

    @classmethod
    def from_json(
        cls,
        path: str | PathLike[str],
        *,
        components: Sequence[str] | None = None,
    ) -> Parameters:
        document = _ParameterDocument.from_json(path)
        selection: Sequence[str] = (
            tuple(component.component_id for component in document.components)
            if components is None
            else components
        )
        return document.select(selection)

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, object],
        *,
        components: Sequence[str] | None = None,
    ) -> Parameters:
        document = _ParameterDocument.from_mapping(mapping)
        selection: Sequence[str] = (
            tuple(component.component_id for component in document.components)
            if components is None
            else components
        )
        return document.select(selection)

    def select(self, component_ids: Sequence[str]) -> Parameters:
        return self._document.select(component_ids)

    def to_mapping(self) -> dict[str, object]:
        return cast(dict[str, object], json.loads(self._document._mapping_json))

    def with_values(
        self,
        values: Mapping[str, object],
        *,
        source_id: str,
        locator: str,
    ) -> Parameters:
        if not values or any(type(identity) is not str for identity in values):
            raise _error("invalid-selection", field="parameter_values")
        specs = {spec.identity: spec for spec in self.parameter_specs}
        if set(values) - set(specs):
            raise _error(
                "missing-record",
                field="parameter_identity",
                value=tuple(sorted(set(values) - set(specs))),
            )
        magnitudes: dict[str, float] = {}
        for identity, value in values.items():
            spec = specs[identity]
            try:
                magnitude = (
                    float(value.to(spec.unit).magnitude)
                    if hasattr(value, "to")
                    else float(cast(SupportsFloat, value))
                )
            except (AttributeError, TypeError, ValueError) as exc:
                raise _error("invalid-value", field=identity, value=value) from exc
            if not math.isfinite(magnitude):
                raise _error("invalid-value", field=identity, value=value)
            magnitudes[identity] = magnitude

        mapping = self.to_mapping()
        sources = mapping["sources"]
        domains = mapping["domains"]
        assert isinstance(sources, list) and isinstance(domains, list)
        domain_id = source_id
        provenance = {"source_id": source_id, "locator": locator, "domain_id": domain_id}
        sources.append(
            {
                "source_id": source_id,
                "citation": "Project-authored ePC-SAFT parameter fit result.",
                "use_basis": "Accepted fitted coefficient values and exact numerical evidence.",
            }
        )
        domains.append({"domain_id": domain_id, "kind": "unknown"})
        found: set[str] = set()

        def replace(item: object) -> bool:
            changed = False
            if isinstance(item, list):
                for value in item:
                    changed = replace(value) or changed
                return changed
            if not isinstance(item, dict):
                return False
            identity = item.get("identity")
            if identity in magnitudes:
                quantity = item.get("value")
                if not isinstance(quantity, dict) or "magnitude" not in quantity:
                    raise _error("invalid-schema", field=str(identity))
                quantity["magnitude"] = magnitudes[str(identity)]
                found.add(str(identity))
                changed = True
            for value in tuple(item.values()):
                changed = replace(value) or changed
            if changed and "provenance" in item:
                item["provenance"] = dict(provenance)
            if changed and isinstance(item.get("source"), dict):
                item["source"] = {"kind": "explicit", "provenance": [dict(provenance)]}
            return changed

        replace(mapping)
        if found != set(magnitudes):
            raise _error(
                "missing-record",
                field="parameter_identity",
                value=tuple(sorted(set(magnitudes) - found)),
            )
        return Parameters.from_mapping(mapping, components=self.component_ids)

    @classmethod
    def _create(
        cls,
        *,
        document: _ParameterDocument,
        component_ids: tuple[str, ...],
        components: tuple[ComponentRecord, ...],
        records: tuple[ParameterRecord, ...],
        parameter_specs: tuple[ParameterSpec, ...],
        topology: AssociationTopology,
    ) -> Parameters:
        instance = object.__new__(cls)
        object.__setattr__(instance, "component_ids", component_ids)
        object.__setattr__(instance, "components", components)
        object.__setattr__(instance, "domains", document.domains)
        object.__setattr__(instance, "records", records)
        object.__setattr__(instance, "model_families", document.model_families)
        object.__setattr__(instance, "topology", topology)
        object.__setattr__(instance, "parameter_specs", parameter_specs)
        object.__setattr__(
            instance,
            "fingerprint",
            _selected_document_fingerprint(document.fingerprint, component_ids),
        )
        object.__setattr__(instance, "_document", document)
        return instance


__all__ = (
    "ActiveParameter",
    "ActiveParameterApplicability",
    "ActiveParameterDomain",
    "ActiveParameterSet",
    "Parameters",
)
