"""Read named, unit-bearing physical parameter declarations.

``parameter_document`` selects and combines retained packet records;
``Parameters.from_mapping`` maps declared names and units into the native
parameter owner (ALG-PARAMETERS). Source values, model choices and units
remain explicit. Selecting components does not fit or validate parameters.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, cast

from . import _core

ParameterError = _core.ParameterError
# Top-level keys: those read below, record metadata, and records the caller reads (reaction sets come from the caller).
_KEYS = {*"schema schema_version components correlations pairs model_families model_coefficients topology".split(),
         *"document_id document_version purpose sources domains materialization".split(),
         *"reaction_correlations empirical_density_correction cation_mixing".split()}


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise _core.ParameterError(f"{name} must be a mapping")
    return value


def _array(value: object, name: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise _core.ParameterError(f"{name} must be an array")
    return value


def _value(value: object, name: str) -> tuple[float, str]:
    item = _mapping(value, name)
    magnitude = item.get("magnitude")
    unit = item.get("unit")
    if not isinstance(unit, str) or magnitude is None:
        raise _core.ParameterError(f"{name} must contain magnitude and unit")
    return cast(float, magnitude), unit


def _entry(identity: object, field: str, index: Sequence[int], value: object, name: str) -> _core.ParameterEntry:
    if not isinstance(identity, str):
        raise _core.ParameterError(f"{name} has no parameter identity")
    magnitude, unit = _value(value, f"{name}.value")
    result = _core.ParameterEntry()
    result.identity = identity
    result.field = field
    result.index = list(index)
    result.magnitude = magnitude
    result.unit = unit
    return result


def _scalar(item: Mapping[str, object], field: str, index: Sequence[int], laws: dict[str, tuple[int, Mapping[str, object]]],
            entries: list[_core.ParameterEntry]) -> None:
    """One scalar record: exactly one of value or temperature_law; laws and volume corrections are numbered by identity."""
    identity = item.get("identity")
    if not isinstance(identity, str):
        raise _core.ParameterError(f"{field} record has no identity")
    if ("value" in item) == ("temperature_law" in item):
        raise _core.ParameterError(f"{identity} needs exactly one of value or temperature_law")
    if "value" in item and field != "association_volume_correction":
        entries.append(_entry(identity, field, index, item["value"], identity))
        return
    number, first = laws.setdefault(identity, (len(laws), item))
    if first != item:
        raise _core.ParameterError(f"conflicting declarations of {identity}")
    if "value" in item:
        entries.append(_entry(identity, field, (number, *index), item["value"], identity))
        return
    law = _mapping(item["temperature_law"], f"{identity}.temperature_law")
    if law.get("form") != "polynomial-plus-reciprocal-temperature":
        raise _core.ParameterError(f"unknown temperature law form {law.get('form')} for {identity}")
    entries.append(_entry(identity, f"{field}_temperature_law", (number, *index), law.get("reference_temperature"), identity))
    if first is not item:
        return
    for term_value in _array(law.get("terms", ()), f"{identity} terms"):
        term = _mapping(term_value, f"{identity} term")
        power, coefficient = term.get("power"), _mapping(term.get("coefficient"), f"{identity} coefficient")
        if not isinstance(power, int) or isinstance(power, bool) or not 0 <= power <= 3:
            raise _core.ParameterError(f"invalid temperature law power {power} for {identity}")
        entries.append(_entry(coefficient.get("identity"), "temperature_law_coefficient", (number, power), coefficient.get("value"), identity))
    if "reciprocal_temperature_slope" in law:
        slope = _mapping(law["reciprocal_temperature_slope"], f"{identity} reciprocal slope")
        entries.append(_entry(slope.get("identity"), "temperature_law_reciprocal_slope", (number,), slope.get("value"), identity))


def _choice(families: Sequence[object], kind: str, default: str | None = None) -> str:
    """The declared choice of one family kind; build_parameters owns the admissible choices."""
    choices = [record.get("choice") for record in (_mapping(item, "model_families") for item in families) if record.get("kind") == kind]
    if len(choices) > 1:
        raise _core.ParameterError(f"duplicate model choice {kind}")
    found = choices[0] if choices else default
    if not isinstance(found, str):
        raise _core.ParameterError(f"unknown {kind} model choice {found}")
    return found


# The model_families fields from_mapping reads.
_FAMILY_FIELDS = ("choice", "ionic_dispersion", "c_shell", "c_dielectric", "solvents", "normalization")


def _named(record: Mapping[str, Any]) -> set[str]:
    """The component ids one component, pair, correlation, site, edge or no-cross record names."""
    ends = [record[key] for key in ("endpoint_a", "endpoint_b") if key in record]
    return {item[key] for item in (record, *ends) for key in ("component_id", "component_id_a", "component_id_b") if key in item}


def _cation_mixing(document: Mapping[str, Any], selected: Sequence[str]) -> dict[str, Any] | None:
    if "cation_mixing" not in document:
        return None
    record = dict(_mapping(document["cation_mixing"], "cation_mixing"))
    names = [record.get("solvent"), *_array(record.get("cations"), "cations"), record.get("anion")]
    available = {c["component_id"] for c in document["components"]}
    if (
        set(record) != {"solvent", "cations", "anion", "coefficients"}
        or len(names) != 4
        or len(_array(record.get("coefficients"), "cation_mixing coefficients")) != 2
        or any(not isinstance(n, str) or n not in available for n in names)
        or len(set(names)) != 4
    ):
        raise ParameterError("cation_mixing requires four distinct known participants and two coefficients")
    record["cations"] = sorted(cast(list[str], names[1:3]))
    if not set(names[1:3]) <= set(selected):
        return None
    if not set(names) <= set(selected):
        raise ParameterError("cation_mixing selection omits its solvent or anion")
    return record


def parameter_document(root: str | Path, packets: Sequence[tuple[str, Sequence[str]]], *,
                       k_ij: Mapping[tuple[str, str], float] | None = None,
                       l_ij: Mapping[tuple[str, str], float] | None = None) -> dict[str, object]:
    """Combine selected packet components into one parameter document.

    Parameters
    ----------
    root : str or pathlib.Path
        Packet root containing ``<packet>/parameters/parameters.json``.
    packets : sequence of tuple
        (Packet identifier/version, component identifiers) in output order.
    k_ij : mapping of tuple of str to float, optional
        Dimensionless unlike dispersion-energy corrections. Each supplied
        value replaces that pair's constant and temperature-law declarations.
    l_ij : mapping of tuple of str to float, optional
        Dimensionless unlike diameter corrections, replacing existing values.

    Returns
    -------
    document : dict
        Version-1 ``epcsaft.parameters`` declaration with retained source/domain
        materialization records. Unit-bearing scalar values keep their units.

    Raises
    ------
    ParameterError
        A packet form is unsupported, a selected component is absent, model
        choices/coefficients conflict, cation mixing is incompletely selected,
        or an override does not name two distinct selected components.
    OSError
        A packet file cannot be read.
    json.JSONDecodeError
        A packet contains invalid JSON; missing required keys raise ``KeyError``.

    Notes
    -----
    Only records wholly within each packet's component selection are retained.
    A component with no association sites can join the general-site model;
    other conflicting model choices are refused. Duplicate components are
    refused later by ``Parameters.from_mapping``, not by this merge function.
    """
    merged: dict[str, Any] = {"schema": "epcsaft.parameters", "schema_version": 1, "components": [], "correlations": [],
                              "pairs": [], "model_coefficients": None, "materialization": [],
                              "topology": {"sites": [], "edges": [], "no_cross_association": []}}
    families: dict[str, tuple[str, Mapping[str, Any]]] = {}
    selected: set[str] = set()
    for packet, names in packets:
        document = json.loads((Path(root) / packet / "parameters" / "parameters.json").read_text())
        if document.get("schema") != "epcsaft.parameters" or document.get("schema_version") != 1 or set(document) - _KEYS:
            raise _core.ParameterError(f"{packet} is not a schema 1 parameter document the Engine reads")
        mixing = _cation_mixing(document, [name for _, ids in packets for name in ids])
        if mixing is not None:
            if not {mixing["solvent"], *mixing["cations"], mixing["anion"]} <= set(names):
                raise ParameterError("cation_mixing must come from one complete source selection")
            if merged.get("cation_mixing", mixing) != mixing:
                raise ParameterError("conflicting cation_mixing declarations")
            merged["cation_mixing"] = mixing
        chosen = set(names)
        if missing := sorted(chosen - {c["component_id"] for c in document["components"]}):
            raise _core.ParameterError(f"{packet} has no component {missing[0]}")
        selected |= chosen  # a component taken twice is refused by from_mapping as a duplicate component_id
        declared = {record["kind"]: record for record in document.get("model_families", ())}
        for kind in ("base", "association", "electrolyte", "permittivity", "polar"):
            record = declared.get(kind, {"kind": kind, "choice": "none"})
            first, kept = families.setdefault(kind, (packet, record))
            this, that = ({k: sorted(r[k]) if k == "solvents" else r[k] for k in _FAMILY_FIELDS if k in r} for r in (record, kept))
            if this == that or (kind == "association" and record["choice"] == "none"):
                continue
            # A component without sites joins general-site exactly; every other difference is refused.
            if kind != "association" or kept["choice"] != "none":
                raise _core.ParameterError(f"{kind} model {this} in {packet} conflicts with {that} in {first}")
            families[kind] = (packet, record)
        coefficients = document.get("model_coefficients", [])
        if merged["model_coefficients"] not in (None, coefficients):
            raise _core.ParameterError(f"model_coefficients of {packet} differ from those of {packets[0][0]}")
        merged["model_coefficients"] = coefficients
        merged["components"] += [c for name in names for c in document["components"] if c["component_id"] == name]
        topology = document.get("topology", {})
        for into, key, items in [(merged, "correlations", document), (merged, "pairs", document)] + [
                (merged["topology"], key, topology) for key in ("sites", "edges", "no_cross_association")]:
            into[key] += [r for r in items.get(key, ()) if _named(r) <= chosen]
        merged["materialization"].append({"packet": packet, "components": list(names)} | {
            key: document.get(key, []) for key in ("sources", "domains")})
    merged["model_families"] = [record for _, record in families.values() if "family_id" in record or record["choice"] != "none"]
    for family, values in (("k_ij", k_ij or {}), ("l_ij", l_ij or {})):
        for (a, b), magnitude in values.items():
            if a == b or not {a, b} <= selected:
                raise _core.ParameterError(f"{family} pair ({a}, {b}) needs two selected components")
            pair = next((p for p in merged["pairs"] if {p["component_id_a"], p["component_id_b"]} == {a, b}), None)
            if pair is None:
                pair = {"component_id_a": a, "component_id_b": b, "coefficients": []}
                merged["pairs"].append(pair)
            # A set k_ij also replaces the pair's k_ij temperature law and reciprocal-temperature slope.
            old = [c for c in pair["coefficients"] if c["family"].startswith(family)]
            pair["coefficients"] = [c for c in pair["coefficients"] if c not in old] + [{
                "identity": next((c["identity"] for c in old if c["family"] == family), f"pair/{a}/{b}/{family}"),
                "family": family, "value": {"magnitude": magnitude, "unit": "dimensionless"},
                "provenance": {"locator": "caller-declared"}}]
    return merged


class Parameters:
    """Keep ordered component names beside their native physical parameter record.

    Attributes
    ----------
    component_ids : tuple of str
        Component order for all native input and output arrays.
    site_ids : tuple of tuple of str
        Association sites in native order, each (component identifier, site name).
    molar_masses : tuple of float
        Component molar masses in kg/mol, in component_ids order, read from
        the native parameter record.

    Notes
    -----
    Prefer ``from_json``, ``from_packets`` or ``from_mapping`` for unit-bearing
    source records. The direct constructor only associates names with an
    existing native record; it does not read a source document.
    """

    __slots__ = ("component_ids", "site_ids", "_native", "_document")

    component_ids: tuple[str, ...]
    site_ids: tuple[tuple[str, str], ...]
    _native: _core.Parameters

    def __init__(
        self,
        native: _core.Parameters,
        component_ids: tuple[str, ...],
        site_ids: tuple[tuple[str, str], ...],
    ) -> None:
        """Associate an existing native record with ordered component and site names.

        Parameters
        ----------
        native : _core.Parameters
            Native physical parameters, already expressed in their native units.
        component_ids : tuple of str
            Component order matching ``native``.
        site_ids : tuple of tuple of str
            (Component, site) order matching native association coordinates.

        Notes
        -----
        This constructor assigns the inputs without validating their agreement.
        It retains no source document for later property/ASL serialization.
        """
        self._native = native
        self._document: _core.ParameterDocument | str | None = None
        self.component_ids = component_ids
        self.site_ids = site_ids

    @property
    def molar_masses(self) -> tuple[float, ...]:  # kg/mol, component_ids order, read from the native record
        """Read the native component molar masses.

        Returns
        -------
        masses : tuple of float
            Field documented with its unit and order in ``Parameters`` Attributes.
        """
        return tuple(self._native.molar_masses)

    @classmethod
    def from_json(cls, path: str | Path, *, components: Sequence[str] | None = None) -> Parameters:
        """Read a parameter JSON file and apply native parameter admission.

        Parameters
        ----------
        path : str or pathlib.Path
            Version-1 parameter document; each physical scalar declares its unit.
        components : sequence of str, optional
            Selected component identifiers in desired order; default is file order.

        Returns
        -------
        parameters : epcsaft.parameters.Parameters
            Ordered names, native physical values and retained parameter document.

        Raises
        ------
        ParameterError
            The document or selection is refused by ``from_mapping``.
        OSError
            The file cannot be read.
        json.JSONDecodeError
            The file contains invalid JSON.
        """
        return cls.from_mapping(cast(Mapping[str, object], json.loads(Path(path).read_text())), components=components)

    @classmethod
    def from_packets(cls, root: str | Path, packets: Sequence[tuple[str, Sequence[str]]], *,
                     k_ij: Mapping[tuple[str, str], float] | None = None,
                     l_ij: Mapping[tuple[str, str], float] | None = None) -> Parameters:
        """Load a native parameter set from selected source packets.

        Parameters
        ----------
        root : str or pathlib.Path
            Root used by ``parameter_document``.
        packets : sequence of tuple
            (Packet identifier/version, ordered component identifiers).
        k_ij : mapping of tuple of str to float, optional
            Dimensionless unlike dispersion corrections replacing source pair laws.
        l_ij : mapping of tuple of str to float, optional
            Dimensionless unlike diameter corrections replacing source pair laws.

        Returns
        -------
        parameters : epcsaft.parameters.Parameters
            Admitted native values with names in packet/component selection order.

        Raises
        ------
        ParameterError
            Packet combination or the resulting parameter declaration is refused.
        OSError
            A source packet is unreadable.
        json.JSONDecodeError
            A packet contains invalid JSON.

        Notes
        -----
        Packet combination follows ``parameter_document``; no fitting or adoption
        of undeclared unlike interactions is performed here.
        """
        return cls.from_mapping(parameter_document(root, packets, k_ij=k_ij, l_ij=l_ij))

    @classmethod
    def from_mapping(
        cls, mapping: Mapping[str, object], *, components: Sequence[str] | None = None
    ) -> Parameters:
        """Map a unit-bearing parameter document into native physical parameters.

        Parameters
        ----------
        mapping : mapping of str to object
            ``epcsaft.parameters`` version 1: components, model choices, correlations,
            pair coefficients and association topology. Scalars carry magnitude and
            unit; temperature intervals use K, molar masses kg/mol, and native
            segment/Born diameters Å. Unit conversion is owned by native admission.
        components : sequence of str, optional
            Unique selected identifiers in output order; default is document order.

        Returns
        -------
        parameters : epcsaft.parameters.Parameters
            Admitted native record, ordered component/site names and source document.

        Raises
        ------
        ParameterError
            Unsupported document form, unread fields, duplicate/unknown identifiers,
            invalid units or physical values, unsupported correlations, conflicting
            parameter identities, model choices, pairs or association topology.
        KeyError
            A required nested record field is missing.

        Notes
        -----
        Only the declared selection is materialized; complete participants are
        required for retained cation-mixing terms. Correlations and general-site
        edges retain named identities. Native ``build_parameters`` owns the
        physical domain, units and admitted model families (ALG-PARAMETERS).
        """
        document = _mapping(mapping, "parameter document")
        if document.get("schema") != "epcsaft.parameters" or document.get("schema_version") != 1:
            raise _core.ParameterError("unsupported parameter document schema")
        if unread := sorted(set(document) - _KEYS):
            raise _core.ParameterError(f"parameter document fields the Engine does not read: {', '.join(unread)}")
        raw_components = _array(document.get("components"), "components")
        component_records: dict[str, Mapping[str, object]] = {}
        for index, item in enumerate(raw_components):
            record = _mapping(item, "components")
            component_id = record.get("component_id")
            if not isinstance(component_id, str) or not component_id:
                raise _core.ParameterError(f"components[{index}].component_id must be a non-empty string")
            if component_id in component_records:
                raise _core.ParameterError(f"duplicate component_id {component_id}")
            component_records[component_id] = record
        available = list(component_records)
        selected = available if components is None else list(components)
        for index, component_id in enumerate(selected):
            if not isinstance(component_id, str) or not component_id:
                raise _core.ParameterError(f"component selection[{index}] must be a non-empty string")
        if len(set(selected)) != len(selected):
            raise _core.ParameterError("duplicate selected component")
        if any(item not in available for item in selected):
            missing = next(item for item in selected if item not in available)
            raise _core.ParameterError(f"unknown component_id {missing}")
        selected_ids = tuple(selected)
        positions = {item: index for index, item in enumerate(selected_ids)}
        entries: list[_core.ParameterEntry] = []
        laws: dict[str, tuple[int, Mapping[str, object]]] = {}
        if mixing := _cation_mixing(document, selected):
            indices = [positions[name] for name in (mixing["solvent"], *mixing["cations"], mixing["anion"])]
            for value in _array(mixing["coefficients"], "cation_mixing coefficients"):
                item = _mapping(value, "cation_mixing coefficient")
                if set(item) != {"identity", "family", "value"} or item.get("family") not in ("cation_pair_mixing", "cation_triple_mixing") or set(_mapping(item.get("value"), "cation_mixing value")) != {"magnitude", "unit"}:
                    raise ParameterError("cation_mixing needs pair and triple coefficients")
                entries.append(_entry(item["identity"], cast(str, item["family"]), indices, item["value"], "cation_mixing"))
        # Document families keep their native field name unless renamed here.
        renamed = {"segment_diameter": "sigma_constant", "ion_specific_suppression_coefficient": "ion_suppression_coefficient"}
        coefficient_families = {"segment_count", "dispersion_energy_over_k", "packing_diameter", "debye_huckel_diameter", "born_diameter",
                                "solvation_factor", "relative_permittivity", "dielectric_volume", "dielectric_temperature", "dipole_scaling",
                                "polarizability_scaling", "correlation_integral_parameter", "dipole_moment", "dipole_count",
                                "quadrupole_moment", "quadrupole_count", *renamed}
        for component_id in selected_ids:
            record = component_records[component_id]
            fixed = _mapping(record.get("fixed", {}), f"component/{component_id}/fixed")
            for family in ("molar_mass", "charge_number", "charged_molecule_ion"):
                if family in fixed:
                    item = _mapping(fixed[family], f"component/{component_id}/{family}")
                    entries.append(_entry(item.get("identity", f"component/{component_id}/{family}"), family, (positions[component_id],), item["value"], f"component/{component_id}/{family}"))
            coefficients = _array(record.get("coefficients", ()), f"component/{component_id}/coefficients")
            for item_value in coefficients:
                item = _mapping(item_value, f"component/{component_id}/coefficients")
                family_value = item.get("family")
                if not isinstance(family_value, str) or family_value not in coefficient_families:
                    raise _core.ParameterError(f"unknown component family {family_value} for {item.get('identity')}")
                identity = item.get("identity")
                if not isinstance(identity, str):
                    raise _core.ParameterError(f"component/{component_id}/coefficients has no identity")
                entries.append(_entry(identity, renamed.get(family_value, family_value), (positions[component_id],), item["value"], identity))

        raw_correlations = _array(document.get("correlations", ()), "correlations")
        sigma_term = 0
        permittivity_term = 0
        for correlation_index, item_value in enumerate(raw_correlations):
            item = _mapping(item_value, "correlation")
            component_id = item.get("component_id")
            if not isinstance(component_id, str) or not component_id:
                raise _core.ParameterError(f"correlation[{correlation_index}].component_id must be a non-empty string")
            if component_id not in positions:
                continue
            correlation_id = item.get("correlation_id")
            family_value = item.get("family")
            form_value = item.get("form")
            if not isinstance(family_value, str) or not isinstance(form_value, str):
                raise _core.ParameterError(f"unknown correlation {correlation_id}")
            family = family_value
            form = form_value
            correlation_position = (positions[component_id],)
            if family not in {"segment_diameter", "relative_permittivity"} or form not in {"constant", "constant-plus-sum-of-exponentials", "constant-plus-polynomial-in-temperature"}:
                raise _core.ParameterError(f"unknown correlation {correlation_id}")
            if form == "constant-plus-sum-of-exponentials" and family != "segment_diameter":
                raise _core.ParameterError(f"unknown correlation form {correlation_id}")
            if form == "constant-plus-polynomial-in-temperature" and family != "relative_permittivity":
                raise _core.ParameterError(f"unknown correlation form {correlation_id}")
            base = "sigma_constant" if family == "segment_diameter" else "relative_permittivity"
            constant = _mapping(item.get("constant"), f"{correlation_id}.constant")
            entries.append(_entry(constant.get("identity", correlation_id), base, correlation_position, constant["value"], cast(str, correlation_id)))
            if form == "constant-plus-sum-of-exponentials":
                terms = _array(item.get("terms", ()), f"{correlation_id} terms")
                for term in terms:
                    term_map = _mapping(term, cast(str, correlation_id))
                    amplitude = _mapping(term_map.get("amplitude"), "sigma amplitude")
                    exponent = _mapping(term_map.get("exponent_coefficient"), "sigma exponent")
                    entries.append(_entry(amplitude.get("identity"), "sigma_term_amplitude", (sigma_term, positions[component_id]), amplitude["value"], cast(str, correlation_id)))
                    entries.append(_entry(exponent.get("identity"), "sigma_term_exponent", (sigma_term, positions[component_id]), exponent["value"], cast(str, correlation_id)))
                    sigma_term += 1
            elif form == "constant-plus-polynomial-in-temperature":
                reference = _mapping(item.get("reference_temperature"), "permittivity reference")
                entries.append(_entry(f"{correlation_id}/reference_temperature", "permittivity_reference_temperature", correlation_position, reference, cast(str, correlation_id)))
                terms = _array(item.get("terms", ()), f"{correlation_id} terms")
                for term in terms:
                    term_map = _mapping(term, cast(str, correlation_id))
                    power = term_map.get("power")
                    if not isinstance(power, int):
                        raise _core.ParameterError(f"invalid polynomial power for {correlation_id}")
                    coefficient = _mapping(term_map.get("coefficient"), "permittivity coefficient")
                    entries.append(_entry(coefficient.get("identity"), "permittivity_term_coefficient", (permittivity_term, positions[component_id], power), coefficient["value"], cast(str, correlation_id)))
                    permittivity_term += 1

        raw_pairs = _array(document.get("pairs", ()), "pairs")
        pair_fields = {"k_ij", "l_ij", "k_ij_reciprocal_temperature_slope"}
        for index, pair_value in enumerate(raw_pairs):
            pair = _mapping(pair_value, "pair")
            a = pair.get("component_id_a")
            b = pair.get("component_id_b")
            if not isinstance(a, str) or not a:
                raise _core.ParameterError(f"pair[{index}].component_id_a must be a non-empty string")
            if not isinstance(b, str) or not b:
                raise _core.ParameterError(f"pair[{index}].component_id_b must be a non-empty string")
            if a not in positions or b not in positions:
                continue
            pair_entries = _array(pair.get("coefficients", ()), f"pair/{a}/{b} coefficients")
            for item_value in pair_entries:
                item = _mapping(item_value, f"pair/{a}/{b}")
                family_value = item.get("family")
                if not isinstance(family_value, str) or family_value not in pair_fields:
                    raise _core.ParameterError(f"unsupported pair coefficient {item.get('identity')}")
                family = family_value
                pair_index = (positions[a], positions[b])
                identity = item.get("identity")
                if not isinstance(identity, str):
                    raise _core.ParameterError(f"pair/{a}/{b} coefficient has no identity")
                if family != "k_ij_reciprocal_temperature_slope":
                    _scalar(item, family, pair_index, laws, entries)
                else:
                    entries.append(_entry(identity, family, pair_index, item["value"], identity))
                    if "reference_temperature" not in item:
                        raise _core.ParameterError(f"missing reference temperature for {item.get('identity')}")
                    reference = _mapping(item["reference_temperature"], "pair reference temperature")
                    entries.append(_entry(f"{identity}/reference_temperature", "k_ij_reference_temperature", pair_index, reference, identity))

        families = _array(document.get("model_families", ()), "model_families")
        for item in families:
            kind = _mapping(item, "model_families").get("kind")
            if kind not in ("base", "association", "electrolyte", "permittivity", "polar"):
                raise _core.ParameterError(f"unknown model family kind {kind}")
        electrolyte = _choice(families, "electrolyte")
        permittivity = _choice(families, "permittivity")
        if permittivity == "solvent-mole-weighting":
            record = next(_mapping(item, "permittivity") for item in families if _mapping(item, "permittivity").get("kind") == "permittivity")
            solvents = _array(record.get("solvents"), "permittivity solvents")
            mode = record.get("normalization")
            if mode not in ("unnormalized", "solvent-normalized"):
                raise ParameterError("explicit permittivity normalization required")
            if not solvents or any(not isinstance(s, str) or s not in component_records for s in solvents):
                raise ParameterError("permittivity solvents must name known components")
            if len(set(cast(Sequence[str], solvents))) != len(solvents):
                raise ParameterError("duplicate permittivity solvent")
            for solvent in cast(Sequence[str], solvents):
                fixed = _mapping(component_records[solvent].get("fixed", {}), "solvent fixed")
                charge = _mapping(fixed.get("charge_number"), "solvent charge")
                if _mapping(charge.get("value"), "solvent charge value").get("magnitude") != 0:
                    raise ParameterError("permittivity solvent must be neutral")
            for name in selected_ids:
                if name in solvents:
                    entries.append(_entry(f"permittivity/solvent/{name}", "permittivity_solvent", (positions[name],),
                                          {"magnitude": 1, "unit": "dimensionless"}, "permittivity"))
            entries.append(_entry("permittivity/normalization", "permittivity_normalization", (),
                                  {"magnitude": int(mode == "solvent-normalized"), "unit": "dimensionless"}, "permittivity"))
        base = _choice(families, "base")
        association = _choice(families, "association")
        polar = _choice(families, "polar", "none")  # optional family (ADR 0011)
        electrolyte_record = next(_mapping(item, "electrolyte family") for item in families if _mapping(item, "electrolyte family").get("kind") == "electrolyte")
        if electrolyte == "born":
            for name, field in (("c_shell", "born_c_shell"), ("c_dielectric", "born_c_dielectric")):
                if name not in electrolyte_record:
                    raise _core.ParameterError(f"missing electrolyte/{name}")
                entries.append(_entry(f"model/electrolyte/{name}", field, (), {"magnitude": electrolyte_record[name], "unit": "dimensionless"}, "model/electrolyte"))
        model_coefficients = _array(document.get("model_coefficients", ()), "model_coefficients")
        for item_value in model_coefficients:
            item = _mapping(item_value, "model coefficient")
            family_value = item.get("family")
            if family_value == "ionic_region_relative_permittivity":
                field = "ionic_region_relative_permittivity"
            elif family_value == "ion_fraction_suppression_coefficient":
                field = "ion_suppression_coefficient"
            else:
                raise _core.ParameterError(f"unknown model coefficient {item.get('identity')}")
            entries.append(_entry(item.get("identity"), field, (), item["value"], cast(str, item.get("identity"))))

        topology = _mapping(document.get("topology", {}), "topology")
        presets = topology.get("presets", ())
        if isinstance(presets, Sequence) and not isinstance(presets, (str, bytes)) and presets:
            raise _core.ParameterError("topology.presets is not supported")
        site_ids: list[tuple[str, str]] = []
        site_component: list[int] = []
        site_multiplicity: list[int] = []
        sites = _array(topology.get("sites", ()), "topology.sites")
        for index, item_value in enumerate(sites):
            item = _mapping(item_value, "topology site")
            component_id = item.get("component_id")
            site_id = item.get("site_id")
            if not isinstance(component_id, str) or not component_id:
                raise _core.ParameterError(f"topology site[{index}].component_id must be a non-empty string")
            if not isinstance(site_id, str) or not site_id:
                raise _core.ParameterError(f"topology site[{index}].site_id must be a non-empty string")
            if component_id not in positions:
                continue
            site_key = (component_id, site_id)
            if site_key in site_ids:
                raise _core.ParameterError(f"duplicate topology site {site_key}")
            site_ids.append(site_key)
            site_component.append(positions[component_id])
            site_multiplicity.append(cast(int, item.get("multiplicity")))
        site_positions = {pair: index for index, pair in enumerate(site_ids)}
        edges = _array(topology.get("edges", ()), "topology.edges")
        selected_edges: list[tuple[tuple[int, int], Mapping[str, object], Mapping[str, object]]] = []
        association_identities: dict[str, tuple[str, tuple[int, int]]] = {}
        for index, item_value in enumerate(edges):
            edge = _mapping(item_value, "topology edge")
            endpoint_a = _mapping(edge.get("endpoint_a"), "topology endpoint")
            endpoint_b = _mapping(edge.get("endpoint_b"), "topology endpoint")
            first_component_id = endpoint_a.get("component_id")
            first_site_id = endpoint_a.get("site_id")
            second_component_id = endpoint_b.get("component_id")
            second_site_id = endpoint_b.get("site_id")
            for value, name in ((first_component_id, "endpoint_a.component_id"),
                                (first_site_id, "endpoint_a.site_id"),
                                (second_component_id, "endpoint_b.component_id"),
                                (second_site_id, "endpoint_b.site_id")):
                if not isinstance(value, str) or not value:
                    raise _core.ParameterError(f"topology edge[{index}].{name} must be a non-empty string")
            first = (cast(str, first_component_id), cast(str, first_site_id))
            second = (cast(str, second_component_id), cast(str, second_site_id))
            if first[0] not in positions or second[0] not in positions:
                continue
            if first not in site_positions or second not in site_positions:
                raise _core.ParameterError(f"missing topology site {(first if first not in site_positions else second)[0]}/{(first if first not in site_positions else second)[1]}")
            source = _mapping(edge.get("source"), "topology edge source")
            edge_index = (site_positions[first], site_positions[second])
            selected_edges.append((edge_index, edge, source))
            if source.get("kind") != "explicit":
                continue
            if "volume_correction" in edge:
                raise _core.ParameterError(f"unconsumed volume_correction on explicit edge {first}/{second}")
            energy = _mapping(edge.get("energy_over_k"), "association energy")
            volume = _mapping(edge.get("volume"), "association volume")
            if not isinstance(energy.get("identity"), str) or not isinstance(volume.get("identity"), str):
                raise _core.ParameterError("association edge identities must be strings")
            for identity, family in ((cast(str, energy.get("identity")), "energy_over_k"), (cast(str, volume.get("identity")), "volume")):
                if identity in association_identities:
                    raise _core.ParameterError(f"duplicate association identity {identity}")
                association_identities[identity] = (family, edge_index)

        def resolve_inputs(values: Sequence[object], family: str) -> list[tuple[int, int]]:
            try:
                found = [association_identities[cast(str, value)] for value in values]
            except KeyError as error:
                raise _core.ParameterError(f"missing association input {error.args[0]}") from error
            if any(item[0] != family for item in found):
                raise _core.ParameterError(f"association {family} input has wrong family")
            return [item[1] for item in found]

        # Emit only after all explicit identities exist, including forward self-edge references.
        for edge_index, edge, source in selected_edges:
            kind = source.get("kind")
            if kind == "explicit":
                energy = _mapping(edge.get("energy_over_k"), "association energy")
                volume = _mapping(edge.get("volume"), "association volume")
                _scalar(energy, "association_energy_over_k", edge_index, laws, entries)
                _scalar(volume, "association_volume", edge_index, laws, entries)
            elif kind == "combining-rule":
                inputs = _mapping(source.get("inputs"), "combining rule inputs")
                energy_inputs = inputs.get("energy_over_k")
                volume_inputs = inputs.get("volume")
                if not isinstance(energy_inputs, Sequence) or not isinstance(volume_inputs, Sequence) or len(energy_inputs) != 2 or len(volume_inputs) != 2:
                    raise _core.ParameterError("combining rule requires two self edges")
                self_energy = resolve_inputs(energy_inputs, "energy_over_k")
                self_volume = resolve_inputs(volume_inputs, "volume")
                if self_energy != self_volume:
                    raise _core.ParameterError(f"association inputs disagree for {source.get('rule_id')}")
                entries.append(_entry(source.get("rule_id"), "association_combining_rule", [*edge_index, *self_energy[0], *self_energy[1]],
                                      {"magnitude": 0.0, "unit": "dimensionless"}, "association combining rule"))
                if "volume_correction" in edge:
                    _scalar(_mapping(edge["volume_correction"], "volume correction"), "association_volume_correction", edge_index, laws, entries)
            else:
                raise _core.ParameterError(f"unknown association edge source {kind}")

        # One unordered component pair whose unlike site edges are all inactive (preserved model, not proven absence).
        declared: set[frozenset[str]] = set()
        for index, value in enumerate(_array(topology.get("no_cross_association", ()), "no_cross_association")):
            item = _mapping(value, "no_cross_association")
            a, b = str(item.get("component_id_a")), str(item.get("component_id_b"))
            if {a, b} - set(component_records) or a == b or frozenset((a, b)) in declared or "provenance" not in item:
                raise ParameterError(f"no_cross_association[{index}]: two known components, once, with provenance")
            declared.add(frozenset((a, b)))
            if a in positions and b in positions:
                ab, zero = (positions[a], positions[b]), {"magnitude": 0.0, "unit": "dimensionless"}
                entries.append(_entry(f"no_cross_association/{a}/{b}", "no_cross_association", ab, zero, a))

        native_document = _core.ParameterDocument()
        native_document.component_ids = list(selected_ids)
        native_document.site_component = site_component
        native_document.site_multiplicity = site_multiplicity
        native_document.base_choice = base
        native_document.association_choice = association
        native_document.electrolyte_choice = electrolyte
        native_document.ionic_dispersion_choice = str(electrolyte_record.get("ionic_dispersion", ""))
        native_document.permittivity_choice = permittivity
        native_document.polar_choice = polar
        native_document.entries = entries
        result = cls(_core.build_parameters(native_document), selected_ids, tuple(site_ids))
        result._document = native_document
        return result
