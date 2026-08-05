from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from typing import Literal, TypeAlias, TypeVar, cast

from ._bundle_io import _create_template, _dump_path, _load_path
from ._canonical import _bundle_fingerprint, _parameter_fingerprint
from ._dictionary_io import _records_from_dictionary
from .records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    CorrelationRecord,
    ModelParameterRecord,
    PairParameterRecord,
    ParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
    _identifier,
)
from .units import _error

_BundlePurpose: TypeAlias = Literal[
    "external-bundle",
    "package-test-fixture",
    "regression-output",
    "user-provided",
]
_EXTERNAL_PURPOSES = frozenset(
    {
        "external-bundle",
        # Persisted Data v1 metadata; normalized below without lookup or path selection.
        "installed-catalog",
        "package-test-fixture",
        "regression-output",
        "user-provided",
    }
)


def _version(value: object) -> int:
    if type(value) is not int or value <= 0:
        raise _error("invalid-value", field="bundle_version", value=value)
    return value


def _purpose(value: object) -> _BundlePurpose:
    if not isinstance(value, str) or value not in _EXTERNAL_PURPOSES:
        raise _error("invalid-value", field="purpose", value=value)
    if value in {"package-test-fixture", "regression-output", "user-provided"}:
        return cast(_BundlePurpose, value)
    return "external-bundle"


_RecordT = TypeVar("_RecordT")


def _typed(
    values: Iterable[object], expected: type[_RecordT], name: str
) -> tuple[_RecordT, ...]:
    try:
        result = tuple(values)
    except TypeError as exc:
        raise _error("invalid-value", field=name, value="not-iterable") from exc
    if not all(isinstance(value, expected) for value in result):
        raise _error("invalid-value", field=name, value="wrong-record-type")
    return cast(tuple[_RecordT, ...], result)


def _unique(records: tuple[object, ...], attribute: str) -> None:
    seen: set[object] = set()
    for record in records:
        key = getattr(record, attribute)
        if key in seen:
            raise _error("duplicate-record", field=attribute, value=key)
        seen.add(key)


def _validate_records(
    sources: tuple[SourceRecord, ...],
    domains: tuple[ValidityDomain, ...],
    components: tuple[ComponentRecord, ...],
    records: tuple[ParameterRecord, ...],
) -> None:
    _unique(sources, "source_id")
    _unique(domains, "domain_id")
    _unique(components, "component_id")
    _unique(records, "record_id")

    source_ids = {record.source_id for record in sources}
    domain_ids = {record.domain_id for record in domains}
    component_ids = {record.component_id for record in components}
    sites = {
        (record.component_id, record.site_id)
        for record in records
        if isinstance(record, SiteRecord)
    }
    semantic_keys: set[tuple[object, ...]] = set()

    for record in records:
        if record.source_id not in source_ids:
            raise _error("invalid-source", record_id=record.record_id, source_id=record.source_id)
        if record.domain_id not in domain_ids:
            raise _error("invalid-domain", record_id=record.record_id, domain_id=record.domain_id)
        if isinstance(record, ModelParameterRecord):
            referenced_components: tuple[str, ...] = ()
        elif isinstance(
            record,
            (
                SingleParameterRecord,
                SiteRecord,
                ConstantCorrelation,
                ConstantPlusSumOfExponentialsCorrelation,
            ),
        ):
            referenced_components = (record.component_id,)
        else:
            referenced_components = (record.component_id_a, record.component_id_b)
        missing = next(
            (
                component_id
                for component_id in referenced_components
                if component_id not in component_ids
            ),
            None,
        )
        if missing is not None:
            raise _error("missing-record", field="component_id", value=missing)

        if isinstance(record, ModelParameterRecord):
            key: tuple[object, ...] = ("model", record.family)
        elif isinstance(record, SingleParameterRecord):
            key = ("single", record.component_id, record.family)
        elif isinstance(record, PairParameterRecord):
            key = ("pair", record.component_id_a, record.component_id_b, record.family)
        elif isinstance(record, SiteRecord):
            key = ("site", record.component_id, record.site_id)
        elif isinstance(record, AssociationParameterRecord):
            endpoints = (
                (record.component_id_a, record.site_id_a),
                (record.component_id_b, record.site_id_b),
            )
            missing_site = next((site for site in endpoints if site not in sites), None)
            if missing_site is not None:
                raise _error("missing-record", field="site", value=missing_site)
            key = ("association", *endpoints, record.family)
        else:
            key = ("correlation", record.component_id, record.family)
        if key in semantic_keys:
            raise _error("duplicate-record", field="semantic-key", value=key)
        semantic_keys.add(key)


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


@dataclass(frozen=True, slots=True, init=False)
class _ParameterBundle:
    bundle_id: str
    bundle_version: int
    purpose: _BundlePurpose
    sources: tuple[SourceRecord, ...]
    domains: tuple[ValidityDomain, ...]
    components: tuple[ComponentRecord, ...]
    records: tuple[ParameterRecord, ...]
    fingerprint: str

    def __init__(self) -> None:
        raise TypeError("use _ParameterBundle.from_path or _ParameterBundle.from_records")

    @classmethod
    def from_path(cls, path: str | PathLike[str]) -> _ParameterBundle:
        """Load and validate one bundle directory in the strict version-1 schema.

        The directory inventory must exactly match its manifest. Scientific records
        retain their declared units, provenance, domains, and canonical fingerprint.
        """
        return cls.from_records(**_load_path(path))

    @classmethod
    def from_records(
        cls,
        *,
        bundle_id: str,
        bundle_version: int,
        purpose: _BundlePurpose,
        sources: Iterable[SourceRecord],
        domains: Iterable[ValidityDomain],
        components: Iterable[ComponentRecord],
        singles: Iterable[SingleParameterRecord] = (),
        pairs: Iterable[PairParameterRecord] = (),
        sites: Iterable[SiteRecord] = (),
        associations: Iterable[AssociationParameterRecord] = (),
        correlations: Iterable[CorrelationRecord] = (),
        models: Iterable[ModelParameterRecord] = (),
    ) -> _ParameterBundle:
        """Create an immutable bundle from explicit typed parameter records.

        Use ``purpose="user-provided"`` for caller-owned inputs. Every scientific
        value retains its source and validity-domain reference; no physical value is
        inferred.
        """
        bundle_id = _identifier(bundle_id, field="bundle_id")
        bundle_version = _version(bundle_version)
        purpose = _purpose(purpose)
        source_records = _typed(sources, SourceRecord, "sources")
        domain_records = _typed(domains, ValidityDomain, "domains")
        component_records = _typed(components, ComponentRecord, "components")
        single_records = _typed(singles, SingleParameterRecord, "singles")
        pair_records = _typed(pairs, PairParameterRecord, "pairs")
        site_records = _typed(sites, SiteRecord, "sites")
        association_records = _typed(
            associations, AssociationParameterRecord, "associations"
        )
        correlation_records = tuple(correlations)
        if not all(
            isinstance(record, (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation))
            for record in correlation_records
        ):
            raise _error("invalid-value", field="correlations", value="wrong-record-type")
        model_records = _typed(models, ModelParameterRecord, "models")
        records: tuple[ParameterRecord, ...] = (
            *single_records,
            *pair_records,
            *site_records,
            *association_records,
            *correlation_records,
            *model_records,
        )
        _validate_records(source_records, domain_records, component_records, records)

        instance = object.__new__(cls)
        object.__setattr__(instance, "bundle_id", bundle_id)
        object.__setattr__(instance, "bundle_version", bundle_version)
        object.__setattr__(instance, "purpose", purpose)
        object.__setattr__(
            instance,
            "sources",
            tuple(sorted(source_records, key=lambda item: item.source_id)),
        )
        object.__setattr__(
            instance,
            "domains",
            tuple(sorted(domain_records, key=lambda item: item.domain_id)),
        )
        object.__setattr__(
            instance,
            "components",
            component_records,
        )
        object.__setattr__(
            instance, "records", tuple(sorted(records, key=lambda item: item.record_id))
        )
        object.__setattr__(instance, "fingerprint", _bundle_fingerprint(instance))
        return instance

    def to_path(self, path: str | PathLike[str]) -> None:
        """Persist this bundle in the strict version-1 directory schema.

        The target must not already exist. A failed write leaves no target directory.
        Reload the saved bundle with :meth:`Parameters.from_bundle`.
        """
        _dump_path(self, path)

    def select(self, component_ids: Sequence[str]) -> Parameters:
        """Select exact component IDs in the order used by every state and result.

        The returned set has an order-sensitive fingerprint. Unknown or duplicate IDs
        raise ``ParameterError``; model construction rejects incomplete source records.
        """
        if isinstance(component_ids, str):
            raise _error("invalid-selection", field="component_ids", value="string")
        try:
            requested = tuple(component_ids)
        except TypeError as exc:
            raise _error("invalid-selection", field="component_ids", value="not-sequence") from exc
        if not requested or any(type(component_id) is not str for component_id in requested):
            raise _error("invalid-selection", field="component_ids", value="invalid-element")
        if len(requested) != len(set(requested)):
            raise _error("invalid-selection", field="component_ids", value="empty-or-duplicate")
        component_map = {record.component_id: record for record in self.components}
        if any(component_id not in component_map for component_id in requested):
            raise _error("invalid-selection", field="component_ids", value=requested)
        selected_records = tuple(
            record
            for record in self.records
            if all(component_id in requested for component_id in _record_components(record))
        )
        return Parameters._create(
            bundle=self,
            component_ids=requested,
            components=tuple(component_map[component_id] for component_id in requested),
            records=selected_records,
        )


@dataclass(frozen=True, slots=True, init=False)
class Parameters:
    component_ids: tuple[str, ...]
    components: tuple[ComponentRecord, ...]
    domains: tuple[ValidityDomain, ...]
    records: tuple[ParameterRecord, ...]
    fingerprint: str
    _bundle: _ParameterBundle

    def __init__(self) -> None:
        raise TypeError("use Parameters.from_dictionary or Parameters.from_bundle")

    @classmethod
    def from_bundle(
        cls,
        path: str | PathLike[str],
        *,
        components: Sequence[str] | None = None,
    ) -> Parameters:
        """Load the strict manifest-plus-table bundle schema and select components."""
        bundle = _ParameterBundle.from_path(path)
        selection = (
            tuple(component.component_id for component in bundle.components)
            if components is None
            else components
        )
        return bundle.select(selection)

    @classmethod
    def from_dictionary(cls, mapping: Mapping[str, object]) -> Parameters:
        """Construct immutable parameters from the strict fixed-unit dictionary schema."""
        values = _records_from_dictionary(mapping)
        bundle = _ParameterBundle.from_records(
            bundle_id=values.bundle_id,
            bundle_version=values.bundle_version,
            purpose=values.purpose,
            sources=values.sources,
            domains=values.domains,
            components=values.components,
            singles=values.singles,
            pairs=values.pairs,
            sites=values.sites,
            associations=values.associations,
            correlations=values.correlations,
            models=values.models,
        )
        return bundle.select(values.component_ids)

    @staticmethod
    def create_bundle_template(
        path: str | PathLike[str],
        *,
        components: Sequence[str],
        profile: Literal[
            "neutral",
            "associating",
            "electrolyte",
            "associating-electrolyte",
        ],
        correlations: bool = False,
    ) -> None:
        """Create a strict non-overwriting parameter-bundle authoring skeleton."""
        if isinstance(components, str):
            raise _error("invalid-selection", field="components", value="string")
        try:
            component_ids = tuple(
                _identifier(component, field="component_id") for component in components
            )
        except TypeError as exc:
            raise _error("invalid-selection", field="components", value="not-sequence") from exc
        if type(correlations) is not bool:
            raise _error("invalid-value", field="correlations", value=correlations)
        _create_template(
            path,
            components=component_ids,
            profile=profile,
            correlations=correlations,
        )

    def to_bundle(self, path: str | PathLike[str]) -> None:
        """Persist the immutable source bundle using the strict directory schema."""
        self._bundle.to_path(path)

    @classmethod
    def _create(
        cls,
        *,
        bundle: _ParameterBundle,
        component_ids: tuple[str, ...],
        components: tuple[ComponentRecord, ...],
        records: tuple[ParameterRecord, ...],
    ) -> Parameters:
        instance = object.__new__(cls)
        object.__setattr__(instance, "component_ids", component_ids)
        object.__setattr__(instance, "components", components)
        object.__setattr__(instance, "domains", bundle.domains)
        object.__setattr__(instance, "records", records)
        object.__setattr__(
            instance,
            "fingerprint",
            _parameter_fingerprint(
                component_ids,
                records,
                bundle.domains,
            ),
        )
        object.__setattr__(instance, "_bundle", bundle)
        return instance


__all__ = ("Parameters",)
