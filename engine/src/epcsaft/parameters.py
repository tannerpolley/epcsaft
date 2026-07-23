from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Iterable, Literal, Sequence, TypeAlias, get_args

from ._bundle_io import _dump_path, _load_path
from ._canonical import _bundle_fingerprint, _provenance_fingerprint, _set_fingerprint
from .records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    CorrelationRecord,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
    _identifier,
)
from .units import _error


BundlePurpose: TypeAlias = Literal[
    "installed-catalog", "package-test-fixture", "user-provided"
]

_PURPOSES = frozenset(get_args(BundlePurpose))


def _version(value: object) -> int:
    if type(value) is not int or value <= 0:
        raise _error("invalid-value", field="bundle_version", value=value)
    return value


def _purpose(value: object) -> str:
    if not isinstance(value, str) or value not in _PURPOSES:
        raise _error("invalid-value", field="purpose", value=value)
    return value


def _typed(values: Iterable[object], expected: type | tuple[type, ...], name: str) -> tuple:
    try:
        result = tuple(values)
    except TypeError as exc:
        raise _error("invalid-value", field=name, value="not-iterable") from exc
    if not all(isinstance(value, expected) for value in result):
        raise _error("invalid-value", field=name, value="wrong-record-type")
    return result


def _unique(records: tuple, attribute: str) -> None:
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
    records: tuple,
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
    semantic_keys: set[tuple] = set()

    for record in records:
        if record.source_id not in source_ids:
            raise _error("invalid-source", record_id=record.record_id, source_id=record.source_id)
        if record.domain_id not in domain_ids:
            raise _error("invalid-domain", record_id=record.record_id, domain_id=record.domain_id)
        if isinstance(record, ModelParameterRecord):
            referenced_components = ()
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
            key = ("model", record.family)
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


def _record_components(record: object) -> tuple[str, ...]:
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
    return (record.component_id_a, record.component_id_b)  # type: ignore[attr-defined]


@dataclass(frozen=True, slots=True, init=False)
class ParameterBundle:
    bundle_id: str
    bundle_version: int
    purpose: BundlePurpose
    sources: tuple[SourceRecord, ...]
    domains: tuple[ValidityDomain, ...]
    components: tuple[ComponentRecord, ...]
    records: tuple
    fingerprint: str

    def __init__(self) -> None:
        raise TypeError("use ParameterBundle.from_catalog, from_path, or from_records")

    @classmethod
    def from_catalog(cls, bundle_id: str, *, version: int) -> ParameterBundle:
        """Load an immutable installed bundle by its exact ID and version.

        Missing identities and malformed installed records raise ``ParameterError``.
        Inspect ``components``, ``domains``, and ``sources`` before selecting a model.
        """
        bundle_id = _identifier(bundle_id, field="bundle_id")
        version = _version(version)
        path = Path(__file__).with_name("data") / "bundles" / bundle_id / str(version)
        if not path.is_dir():
            raise _error("missing-record", field="catalog", value=f"{bundle_id}/{version}")
        return cls._from_catalog_path(bundle_id, version, path)

    @classmethod
    def catalog(cls) -> tuple[ParameterBundle, ...]:
        """Return every installed catalog bundle in stable identity order.

        The returned tuple and each bundle are immutable. Catalog discovery performs
        the same schema, provenance, and fingerprint validation as ``from_catalog``.
        """
        root = Path(__file__).with_name("data") / "bundles"
        identities: list[tuple[str, int]] = []
        for bundle_path in root.iterdir():
            if not bundle_path.is_dir():
                continue
            for version_path in bundle_path.iterdir():
                if not version_path.is_dir():
                    continue
                try:
                    version = int(version_path.name)
                except ValueError as exc:
                    raise _error(
                        "invalid-schema",
                        field="catalog-version",
                        value=version_path.name,
                    ) from exc
                if version <= 0 or str(version) != version_path.name:
                    raise _error(
                        "invalid-schema",
                        field="catalog-version",
                        value=version_path.name,
                    )
                identities.append((bundle_path.name, version))
        return tuple(
            cls.from_catalog(bundle_id, version=version)
            for bundle_id, version in sorted(identities)
        )

    @classmethod
    def _from_catalog_path(
        cls,
        bundle_id: str,
        version: int,
        path: str | PathLike[str],
    ) -> ParameterBundle:
        bundle = cls.from_path(path)
        if bundle.bundle_id != bundle_id or bundle.bundle_version != version:
            raise _error(
                "invalid-schema",
                field="catalog-identity",
                value=f"{bundle.bundle_id}/{bundle.bundle_version}",
            )
        if bundle.purpose != "installed-catalog":
            raise _error(
                "invalid-schema",
                field="catalog-purpose",
                value=bundle.purpose,
            )
        return bundle

    @classmethod
    def from_path(cls, path: str | PathLike[str]) -> ParameterBundle:
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
        purpose: BundlePurpose,
        sources: Iterable[SourceRecord],
        domains: Iterable[ValidityDomain],
        components: Iterable[ComponentRecord],
        singles: Iterable[SingleParameterRecord] = (),
        pairs: Iterable[PairParameterRecord] = (),
        sites: Iterable[SiteRecord] = (),
        associations: Iterable[AssociationParameterRecord] = (),
        correlations: Iterable[CorrelationRecord] = (),
        models: Iterable[ModelParameterRecord] = (),
    ) -> ParameterBundle:
        """Create an immutable bundle from explicit typed parameter records.

        Use ``purpose="user-provided"`` for caller-owned inputs. Every scientific
        value retains its source and validity-domain reference; no physical value is
        inferred. Persist the result with :meth:`ParameterBundle.to_path`.
        """
        bundle_id = _identifier(bundle_id, field="bundle_id")
        bundle_version = _version(bundle_version)
        purpose = _purpose(purpose)
        source_records = _typed(sources, SourceRecord, "sources")
        domain_records = _typed(domains, ValidityDomain, "domains")
        component_records = _typed(components, ComponentRecord, "components")
        grouped_records = (
            _typed(singles, SingleParameterRecord, "singles"),
            _typed(pairs, PairParameterRecord, "pairs"),
            _typed(sites, SiteRecord, "sites"),
            _typed(associations, AssociationParameterRecord, "associations"),
            _typed(
                correlations,
                (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation),
                "correlations",
            ),
            _typed(models, ModelParameterRecord, "models"),
        )
        records = tuple(record for group in grouped_records for record in group)
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
            tuple(sorted(component_records, key=lambda item: item.component_id)),
        )
        object.__setattr__(
            instance, "records", tuple(sorted(records, key=lambda item: item.record_id))
        )
        object.__setattr__(instance, "fingerprint", _bundle_fingerprint(instance))
        return instance

    def to_path(self, path: str | PathLike[str]) -> None:
        """Persist this bundle in the strict version-1 directory schema.

        The target must not already exist. A failed write leaves no target directory.
        Reload the saved bundle with :meth:`ParameterBundle.from_path`.
        """
        _dump_path(self, path)

    def select(self, component_ids: Sequence[str]) -> ParameterSet:
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
        return ParameterSet._create(
            bundle=self,
            component_ids=requested,
            components=tuple(component_map[component_id] for component_id in requested),
            records=selected_records,
        )


@dataclass(frozen=True, slots=True, init=False)
class ParameterSet:
    bundle_id: str
    bundle_version: int
    bundle_purpose: BundlePurpose
    bundle_fingerprint: str
    component_ids: tuple[str, ...]
    components: tuple[ComponentRecord, ...]
    domains: tuple[ValidityDomain, ...]
    records: tuple
    fingerprint: str
    _source_domain_fingerprint: str

    def __init__(self) -> None:
        raise TypeError("use ParameterBundle.select")

    @classmethod
    def _create(
        cls,
        *,
        bundle: ParameterBundle,
        component_ids: tuple[str, ...],
        components: tuple[ComponentRecord, ...],
        records: tuple,
    ) -> ParameterSet:
        instance = object.__new__(cls)
        object.__setattr__(instance, "bundle_id", bundle.bundle_id)
        object.__setattr__(instance, "bundle_version", bundle.bundle_version)
        object.__setattr__(instance, "bundle_purpose", bundle.purpose)
        object.__setattr__(instance, "bundle_fingerprint", bundle.fingerprint)
        object.__setattr__(instance, "component_ids", component_ids)
        object.__setattr__(instance, "components", components)
        object.__setattr__(instance, "domains", bundle.domains)
        object.__setattr__(instance, "records", records)
        source_ids = {record.source_id for record in records}
        domain_ids = {record.domain_id for record in records}
        object.__setattr__(
            instance,
            "_source_domain_fingerprint",
            _provenance_fingerprint(
                tuple(source for source in bundle.sources if source.source_id in source_ids),
                tuple(domain for domain in bundle.domains if domain.domain_id in domain_ids),
            ),
        )
        object.__setattr__(
            instance,
            "fingerprint",
            _set_fingerprint(
                bundle.fingerprint,
                component_ids,
                tuple(record.record_id for record in records),
            ),
        )
        return instance


__all__ = ("BundlePurpose", "ParameterBundle", "ParameterSet")
