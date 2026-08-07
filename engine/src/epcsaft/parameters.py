from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from typing import SupportsFloat, cast

from ._canonical import _document_fingerprint, _selected_document_fingerprint
from ._dictionary_io import _records_from_mapping
from ._document_io import _load_json
from .records import (
    AssociationTopology,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
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
class ActiveParameter:
    """One exact continuous coordinate selected from ``Parameters``."""

    identity: str
    family: str
    value: object
    unit: str


@dataclass(frozen=True, slots=True, init=False)
class ActiveParameterSet:
    """An immutable ordered selection of canonical parameter identities."""

    parameters: Parameters
    coordinates: tuple[ActiveParameter, ...]
    fingerprint: str

    def __init__(self, parameters: Parameters, identities: Sequence[str]) -> None:
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
        missing = tuple(identity for identity in requested if identity not in available)
        if missing:
            raise _error("missing-record", field="parameter_identity", value=missing)
        selected = tuple(available[identity] for identity in requested)
        discrete = tuple(spec.identity for spec in selected if not spec.continuous)
        if discrete:
            raise _error("invalid-selection", field="parameter_identity", value=discrete)
        coordinates = tuple(
            ActiveParameter(spec.identity, spec.family, spec.value, spec.unit) for spec in selected
        )
        payload = json.dumps(
            {
                "schema": "epcsaft.active-parameter-set",
                "parameters": parameters.fingerprint,
                "identities": requested,
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


__all__ = ("ActiveParameter", "ActiveParameterSet", "Parameters")
