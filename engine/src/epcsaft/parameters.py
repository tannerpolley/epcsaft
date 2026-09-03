from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import cast

from . import _core

ParameterError = _core.ParameterError


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise _core.ParameterError(f"{name} must be a mapping")
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


def _choice(families: Sequence[object], kind: str, allowed: set[str]) -> str:
    found: object | None = None
    for item in families:
        record = _mapping(item, "model_families")
        if record.get("kind") == kind:
            if found is not None:
                raise _core.ParameterError(f"duplicate model choice {kind}")
            found = record.get("choice")
    if not isinstance(found, str) or found not in allowed:
        raise _core.ParameterError(f"unknown {kind} model choice {found}")
    return found


class Parameters:
    __slots__ = ("component_ids", "site_ids", "_native")

    component_ids: tuple[str, ...]
    site_ids: tuple[tuple[str, str], ...]
    _native: _core.Parameters

    def __init__(
        self,
        native: _core.Parameters,
        component_ids: tuple[str, ...],
        site_ids: tuple[tuple[str, str], ...],
    ) -> None:
        self._native = native
        self.component_ids = component_ids
        self.site_ids = site_ids

    @classmethod
    def from_json(cls, path: str | Path, *, components: Sequence[str] | None = None) -> Parameters:
        return cls.from_mapping(cast(Mapping[str, object], json.loads(Path(path).read_text())), components=components)

    @classmethod
    def from_mapping(
        cls, mapping: Mapping[str, object], *, components: Sequence[str] | None = None
    ) -> Parameters:
        document = _mapping(mapping, "parameter document")
        if document.get("schema") != "epcsaft.parameters" or document.get("schema_version") != 1:
            raise _core.ParameterError("unsupported parameter document schema")
        raw_components = document.get("components")
        if not isinstance(raw_components, Sequence) or isinstance(raw_components, (str, bytes)):
            raise _core.ParameterError("components must be an array")
        records = [_mapping(item, "components") for item in raw_components]
        available = [cast(str, item.get("component_id")) for item in records]
        selected = available if components is None else list(components)
        if len(set(selected)) != len(selected):
            raise _core.ParameterError("duplicate selected component")
        if any(item not in available for item in selected):
            missing = next(item for item in selected if item not in available)
            raise _core.ParameterError(f"unknown component_id {missing}")
        selected_ids = tuple(selected)
        positions = {item: index for index, item in enumerate(selected_ids)}
        entries: list[_core.ParameterEntry] = []
        component_records = {cast(str, item.get("component_id")): item for item in records}
        component_fields = {
            "molar_mass": "molar_mass",
            "charge_number": "charge_number",
        }
        coefficient_fields = {
            "segment_count": "segment_count",
            "segment_diameter": "sigma_constant",
            "dispersion_energy_over_k": "dispersion_energy_over_k",
            "packing_diameter": "packing_diameter",
            "debye_huckel_diameter": "debye_huckel_diameter",
            "born_diameter": "born_diameter",
            "solvation_factor": "solvation_factor",
            "relative_permittivity": "relative_permittivity",
            "ion_specific_suppression_coefficient": "ion_suppression_coefficient",
        }
        for component_id in selected_ids:
            record = component_records[component_id]
            fixed = _mapping(record.get("fixed", {}), f"component/{component_id}/fixed")
            for family, field in component_fields.items():
                if family in fixed:
                    item = _mapping(fixed[family], f"component/{component_id}/{family}")
                    entries.append(_entry(item.get("identity", f"component/{component_id}/{family}"), field, (positions[component_id],), item["value"], f"component/{component_id}/{family}"))
                elif family == "charge_number":
                    entries.append(_entry(f"component/{component_id}/charge_number", field, (positions[component_id],), {"magnitude": 0, "unit": "elementary-charge"}, f"component/{component_id}/{family}"))
            coefficients = record.get("coefficients", ())
            if not isinstance(coefficients, Sequence) or isinstance(coefficients, (str, bytes)):
                raise _core.ParameterError(f"component/{component_id}/coefficients must be an array")
            for item_value in coefficients:
                item = _mapping(item_value, f"component/{component_id}/coefficients")
                family_value = item.get("family")
                if not isinstance(family_value, str) or family_value not in coefficient_fields:
                    raise _core.ParameterError(f"unknown component family {family_value} for {item.get('identity')}")
                identity = item.get("identity")
                if not isinstance(identity, str):
                    raise _core.ParameterError(f"component/{component_id}/coefficients has no identity")
                entries.append(_entry(identity, coefficient_fields[family_value], (positions[component_id],), item["value"], identity))

        raw_correlations = document.get("correlations", ())
        if not isinstance(raw_correlations, Sequence) or isinstance(raw_correlations, (str, bytes)):
            raise _core.ParameterError("correlations must be an array")
        sigma_term = 0
        permittivity_term = 0
        for item_value in raw_correlations:
            item = _mapping(item_value, "correlation")
            component_id = cast(str, item.get("component_id"))
            if component_id not in positions:
                continue
            correlation_id = item.get("correlation_id")
            family_value = item.get("family")
            form_value = item.get("form")
            if not isinstance(family_value, str) or not isinstance(form_value, str):
                raise _core.ParameterError(f"unknown correlation {correlation_id}")
            family = family_value
            form = form_value
            index = (positions[component_id],)
            if family not in {"segment_diameter", "relative_permittivity"} or form not in {"constant", "constant-plus-sum-of-exponentials", "constant-plus-polynomial-in-temperature"}:
                raise _core.ParameterError(f"unknown correlation {correlation_id}")
            if form == "constant-plus-sum-of-exponentials" and family != "segment_diameter":
                raise _core.ParameterError(f"unknown correlation form {correlation_id}")
            if form == "constant-plus-polynomial-in-temperature" and family != "relative_permittivity":
                raise _core.ParameterError(f"unknown correlation form {correlation_id}")
            base = "sigma_constant" if family == "segment_diameter" else "relative_permittivity"
            constant = _mapping(item.get("constant"), f"{correlation_id}.constant")
            entries.append(_entry(constant.get("identity", correlation_id), base, index, constant["value"], cast(str, correlation_id)))
            if form == "constant-plus-sum-of-exponentials":
                terms = item.get("terms", ())
                if not isinstance(terms, Sequence) or isinstance(terms, (str, bytes)):
                    raise _core.ParameterError(f"{correlation_id} terms must be an array")
                for term in terms:
                    term_map = _mapping(term, cast(str, correlation_id))
                    amplitude = _mapping(term.get("amplitude"), "sigma amplitude")
                    exponent = _mapping(term.get("exponent_coefficient"), "sigma exponent")
                    entries.append(_entry(amplitude.get("identity"), "sigma_term_amplitude", (sigma_term, positions[component_id]), amplitude["value"], cast(str, correlation_id)))
                    entries.append(_entry(exponent.get("identity"), "sigma_term_exponent", (sigma_term, positions[component_id]), exponent["value"], cast(str, correlation_id)))
                    sigma_term += 1
            elif form == "constant-plus-polynomial-in-temperature":
                reference = _mapping(item.get("reference_temperature"), "permittivity reference")
                entries.append(_entry(f"{correlation_id}/reference_temperature", "permittivity_reference_temperature", index, reference, cast(str, correlation_id)))
                terms = item.get("terms", ())
                if not isinstance(terms, Sequence) or isinstance(terms, (str, bytes)):
                    raise _core.ParameterError(f"{correlation_id} terms must be an array")
                for term in terms:
                    term_map = _mapping(term, cast(str, correlation_id))
                    power = term_map.get("power")
                    if not isinstance(power, int):
                        raise _core.ParameterError(f"invalid polynomial power for {correlation_id}")
                    coefficient = _mapping(term_map.get("coefficient"), "permittivity coefficient")
                    entries.append(_entry(coefficient.get("identity"), "permittivity_term_coefficient", (permittivity_term, positions[component_id], power), coefficient["value"], cast(str, correlation_id)))
                    permittivity_term += 1

        raw_pairs = document.get("pairs", ())
        if not isinstance(raw_pairs, Sequence) or isinstance(raw_pairs, (str, bytes)):
            raise _core.ParameterError("pairs must be an array")
        pair_fields = {"k_ij", "l_ij", "k_ij_reciprocal_temperature_slope"}
        for pair_value in raw_pairs:
            pair = _mapping(pair_value, "pair")
            a = cast(str, pair.get("component_id_a"))
            b = cast(str, pair.get("component_id_b"))
            if a not in positions or b not in positions:
                continue
            pair_entries = pair.get("coefficients", ())
            if not isinstance(pair_entries, Sequence) or isinstance(pair_entries, (str, bytes)):
                raise _core.ParameterError(f"pair/{a}/{b} coefficients must be an array")
            for item_value in pair_entries:
                item = _mapping(item_value, f"pair/{a}/{b}")
                family_value = item.get("family")
                if not isinstance(family_value, str) or family_value == "k_hb_ij" or family_value not in pair_fields:
                    raise _core.ParameterError(f"unsupported pair coefficient {item.get('identity')}")
                family = family_value
                pair_index = (positions[a], positions[b])
                identity = item.get("identity")
                if not isinstance(identity, str):
                    raise _core.ParameterError(f"pair/{a}/{b} coefficient has no identity")
                entries.append(_entry(identity, family, pair_index, item["value"], identity))
                if family == "k_ij_reciprocal_temperature_slope":
                    if "reference_temperature" not in item:
                        raise _core.ParameterError(f"missing reference temperature for {item.get('identity')}")
                    reference = _mapping(item["reference_temperature"], "pair reference temperature")
                    entries.append(_entry(f"{identity}/reference_temperature", "k_ij_reference_temperature", pair_index, reference, identity))

        families = document.get("model_families", ())
        if not isinstance(families, Sequence) or isinstance(families, (str, bytes)):
            raise _core.ParameterError("model_families must be an array")
        electrolyte = _choice(families, "electrolyte", {"none", "fully-dissociated-debye-huckel", "born"})
        permittivity = _choice(families, "permittivity", {"none", "solvent-only", "ion-fraction-suppression", "ion-specific-suppression", "component-permittivity-mixing"})
        if _choice(families, "base", {"pc-saft"}) != "pc-saft":
            raise _core.ParameterError("unknown base model choice")
        _choice(families, "association", {"none", "general-site"})
        electrolyte_record = next(_mapping(item, "electrolyte family") for item in families if _mapping(item, "electrolyte family").get("kind") == "electrolyte")
        if electrolyte == "born":
            for name, field in (("c_shell", "born_c_shell"), ("c_dielectric", "born_c_dielectric")):
                if name not in electrolyte_record:
                    raise _core.ParameterError(f"missing electrolyte/{name}")
                entries.append(_entry(f"model/electrolyte/{name}", field, (), {"magnitude": electrolyte_record[name], "unit": "dimensionless"}, "model/electrolyte"))
        model_coefficients = document.get("model_coefficients", ())
        if not isinstance(model_coefficients, Sequence) or isinstance(model_coefficients, (str, bytes)):
            raise _core.ParameterError("model_coefficients must be an array")
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
        sites = topology.get("sites", ())
        if not isinstance(sites, Sequence) or isinstance(sites, (str, bytes)):
            raise _core.ParameterError("topology.sites must be an array")
        for item_value in sites:
            item = _mapping(item_value, "topology site")
            component_id = cast(str, item.get("component_id"))
            site_id = cast(str, item.get("site_id"))
            if component_id not in positions:
                continue
            site_key = (component_id, site_id)
            if site_key in site_ids:
                raise _core.ParameterError(f"duplicate topology site {site_key}")
            site_ids.append(site_key)
            site_component.append(positions[component_id])
            site_multiplicity.append(cast(int, item.get("multiplicity")))
        site_positions = {pair: index for index, pair in enumerate(site_ids)}
        edges = topology.get("edges", ())
        if not isinstance(edges, Sequence) or isinstance(edges, (str, bytes)):
            raise _core.ParameterError("topology.edges must be an array")
        association_identities: dict[str, tuple[int, int]] = {}
        for item_value in edges:
            edge = _mapping(item_value, "topology edge")
            endpoint_a = _mapping(edge.get("endpoint_a"), "topology endpoint")
            endpoint_b = _mapping(edge.get("endpoint_b"), "topology endpoint")
            first = (cast(str, endpoint_a.get("component_id")), cast(str, endpoint_a.get("site_id")))
            second = (cast(str, endpoint_b.get("component_id")), cast(str, endpoint_b.get("site_id")))
            if first not in site_positions or second not in site_positions:
                continue
            source = _mapping(edge.get("source"), "topology edge source")
            kind = source.get("kind")
            edge_index = (site_positions[first], site_positions[second])
            if kind == "explicit":
                energy = _mapping(edge.get("energy_over_k"), "association energy")
                volume = _mapping(edge.get("volume"), "association volume")
                entries.append(_entry(energy.get("identity"), "association_energy_over_k", edge_index, energy["value"], "association edge"))
                entries.append(_entry(volume.get("identity"), "association_volume", edge_index, volume["value"], "association edge"))
                association_identities[cast(str, energy.get("identity"))] = edge_index
                association_identities[cast(str, volume.get("identity"))] = edge_index
            elif kind == "combining-rule":
                if source.get("rule_id") != "arithmetic-energy-geometric-volume":
                    raise _core.ParameterError(f"unknown association rule {source.get('rule_id')}")
                inputs = _mapping(source.get("inputs"), "combining rule inputs")
                energy_inputs = inputs.get("energy_over_k")
                volume_inputs = inputs.get("volume")
                if not isinstance(energy_inputs, Sequence) or not isinstance(volume_inputs, Sequence) or len(energy_inputs) != 2 or len(volume_inputs) != 2:
                    raise _core.ParameterError("combining rule requires two self edges")
                try:
                    self_energy = [association_identities[cast(str, value)] for value in energy_inputs]
                    self_volume = [association_identities[cast(str, value)] for value in volume_inputs]
                except KeyError as error:
                    raise _core.ParameterError(f"missing association input {error.args[0]}") from error
                if self_energy != self_volume:
                    raise _core.ParameterError(f"association inputs disagree for {source.get('rule_id')}")
                entries.append(_core.ParameterEntry())
                rule = entries[-1]
                rule.identity = cast(str, source.get("rule_id"))
                rule.field = "association_combining_rule"
                rule.index = [*edge_index, *self_energy[0], *self_energy[1]]
                rule.magnitude = 0.0
                rule.unit = "dimensionless"
            else:
                raise _core.ParameterError(f"unknown association edge source {kind}")

        native_document = _core.ParameterDocument()
        native_document.component_ids = list(selected_ids)
        native_document.site_component = site_component
        native_document.site_multiplicity = site_multiplicity
        native_document.electrolyte_choice = electrolyte
        native_document.ionic_dispersion_choice = str(electrolyte_record.get("ionic_dispersion", ""))
        native_document.permittivity_choice = permittivity
        native_document.entries = entries
        return cls(_core.build_parameters(native_document), selected_ids, tuple(site_ids))
