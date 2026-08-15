from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, cast

import epcsaft

from ..reference import SourceReferenceIonPair
from ._api import (
    AllComponents,
    ChemicalEquilibriumConstant,
    ChemicalEquilibriumProblem,
    ChemicalStandardState,
    DeclaredNeutralComponents,
    NeutralComponents,
    PhaseSupport,
)
from ._grepe import (
    EquilibriumOutput,
    FinitePhaseStart,
    Fixed,
    GeneralReactiveEquilibriumProblem,
    IdealGasModel,
    ProviderModel,
    ReactivePhase,
    Solved,
)
from ._observations import HomogeneousContinuationReference
from ._reactive_vle import ReactiveBubbleContinuationReference


def _mapping(value: object, fields: frozenset[str], label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise TypeError(f"{label} must be an object with string keys")
    if set(value) != fields:
        raise ValueError(
            f"{label} fields differ: expected {tuple(sorted(fields))}, "
            f"received {tuple(sorted(value))}"
        )
    return cast(Mapping[str, object], value)


def _sequence(value: object, label: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be a sequence")
    return value


def _number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{label} must be an integer")
    return value


def _text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise TypeError(f"{label} must be a nonempty string")
    return value


def _optional_text(value: object, label: str) -> str | None:
    return None if value is None else _text(value, label)


def _floats(value: object, label: str) -> tuple[float, ...]:
    return tuple(_number(item, label) for item in _sequence(value, label))


def _rows(value: object, label: str) -> tuple[tuple[float, ...], ...]:
    return tuple(
        _floats(row, f"{label}[{index}]") for index, row in enumerate(_sequence(value, label))
    )


def _reaction_system(value: object) -> ChemicalEquilibriumProblem:
    item = _mapping(
        value,
        frozenset(
            {
                "species_ids",
                "charges",
                "molar_masses_kg_per_mol",
                "balance_matrix",
                "conserved_totals",
                "reaction_matrix",
                "feed_amounts_mol",
                "equilibrium_constants",
                "strict_interior_amount_floor_mol",
                "source_standard_state",
            }
        ),
        "GREPE reaction system",
    )
    constants = []
    for raw in _sequence(item["equilibrium_constants"], "equilibrium constants"):
        record = _sequence(raw, "equilibrium constant")
        if len(record) != 6 or not isinstance(record[5], bool):
            raise ValueError("equilibrium constants require six typed fields")
        constants.append(
            ChemicalEquilibriumConstant(
                _number(record[0], "equilibrium constant"),
                _text(record[1], "equilibrium source"),
                _text(record[2], "equilibrium reference"),
                cast(Any, _text(record[3], "reaction orientation")),
                _text(record[4], "equilibrium conversion"),
                record[5],
            )
        )
    standard_value = item["source_standard_state"]
    standard = None
    if standard_value is not None:
        source = _mapping(
            standard_value,
            frozenset(
                {
                    "id",
                    "activity_scale_id",
                    "log_activity_scale_factors",
                    "reference_pressure_pa",
                    "component_ids",
                    "solvent_composition",
                    "ion_pairs",
                    "phase",
                    "reference_convention",
                    "activity_convention_id",
                    "standard_molality_mol_per_kg",
                }
            ),
            "GREPE source standard state",
        )
        ion_pairs = []
        for raw in _sequence(source["ion_pairs"], "source ion pairs"):
            record = _sequence(raw, "source ion pair")
            if len(record) != 2:
                raise ValueError("source ion pairs require components and stoichiometry")
            components = tuple(
                _text(item, "ion component") for item in _sequence(record[0], "ion components")
            )
            stoichiometry = tuple(
                _integer(item, "ion stoichiometry")
                for item in _sequence(record[1], "ion stoichiometry")
            )
            if len(components) != 2 or len(stoichiometry) != 2:
                raise ValueError("source ion pairs require two components")
            ion_pairs.append(SourceReferenceIonPair(components, stoichiometry))
        standard = ChemicalStandardState(
            _text(source["id"], "standard-state identity"),
            _text(source["activity_scale_id"], "activity-scale identity"),
            _floats(source["log_activity_scale_factors"], "activity-scale factors"),
            _number(source["reference_pressure_pa"], "reference pressure"),
            tuple(
                _text(item, "reference component")
                for item in _sequence(source["component_ids"], "reference components")
            ),
            _floats(source["solvent_composition"], "reference solvent composition"),
            tuple(ion_pairs),
            cast(Any, _text(source["phase"], "reference phase")),
            cast(Any, _text(source["reference_convention"], "reference convention")),
            _text(source["activity_convention_id"], "activity convention"),
            _number(source["standard_molality_mol_per_kg"], "standard molality"),
        )
    return ChemicalEquilibriumProblem(
        tuple(
            _text(item, "species identity")
            for item in _sequence(item["species_ids"], "species identities")
        ),
        tuple(_integer(item, "charge") for item in _sequence(item["charges"], "charges")),
        _floats(item["molar_masses_kg_per_mol"], "molar masses"),
        _rows(item["balance_matrix"], "balance matrix"),
        _floats(item["conserved_totals"], "conserved totals"),
        _rows(item["reaction_matrix"], "reaction matrix"),
        _floats(item["feed_amounts_mol"], "feed amounts"),
        tuple(constants),
        _number(item["strict_interior_amount_floor_mol"], "interior amount floor"),
        standard,
    )


def _homogeneous_reference(value: object) -> HomogeneousContinuationReference | None:
    if value is None:
        return None
    item = _mapping(
        value,
        frozenset(
            {
                "mole_fractions",
                "molar_volume_m3_per_mol",
                "maximum_log_composition_distance",
                "maximum_log_volume_distance",
                "phase_identity",
                "phase_role",
                "anchor_state_identity",
                "anchor_phase_role_receipt_fingerprint",
                "anchor_parameter_fingerprint",
                "anchor_reaction_fingerprint",
                "anchor_topology_fingerprint",
                "anchor_basin_states",
            }
        ),
        "GREPE homogeneous continuation reference",
    )
    basin_states = []
    for raw in _sequence(item["anchor_basin_states"], "anchor basin states"):
        record = _mapping(
            raw,
            frozenset({"mole_fractions", "molar_volume_m3_per_mol"}),
            "anchor basin state",
        )
        basin_states.append(
            (
                _floats(record["mole_fractions"], "anchor basin composition"),
                _number(record["molar_volume_m3_per_mol"], "anchor basin molar volume"),
            )
        )
    return HomogeneousContinuationReference(
        _floats(item["mole_fractions"], "reference composition"),
        _number(item["molar_volume_m3_per_mol"], "reference molar volume"),
        _number(item["maximum_log_composition_distance"], "composition distance"),
        _number(item["maximum_log_volume_distance"], "volume distance"),
        _text(item["phase_identity"], "reference phase"),
        cast(Any, _text(item["phase_role"], "reference phase role")),
        _text(item["anchor_state_identity"], "anchor state"),
        _text(item["anchor_phase_role_receipt_fingerprint"], "anchor phase-role receipt"),
        _text(item["anchor_parameter_fingerprint"], "anchor parameter fingerprint"),
        _text(item["anchor_reaction_fingerprint"], "anchor reaction fingerprint"),
        _text(item["anchor_topology_fingerprint"], "anchor topology fingerprint"),
        tuple(basin_states),
    )


def _bubble_reference(value: object) -> ReactiveBubbleContinuationReference | None:
    if value is None:
        return None
    item = _mapping(
        value,
        frozenset(
            {
                "pressure_pa",
                "vapor_mole_fractions",
                "liquid_molar_volume_m3_per_mol",
                "vapor_molar_volume_m3_per_mol",
                "maximum_log_pressure_distance",
                "maximum_log_composition_distance",
                "maximum_log_volume_distance",
                "anchor_state_identity",
                "anchor_problem_fingerprint",
                "anchor_topology_fingerprint",
                "anchor_parameter_fingerprint",
            }
        ),
        "GREPE coupled continuation reference",
    )
    return ReactiveBubbleContinuationReference(
        _number(item["pressure_pa"], "reference pressure"),
        _floats(item["vapor_mole_fractions"], "reference vapor composition"),
        _number(item["liquid_molar_volume_m3_per_mol"], "liquid molar volume"),
        _number(item["vapor_molar_volume_m3_per_mol"], "vapor molar volume"),
        _number(item["maximum_log_pressure_distance"], "pressure distance"),
        _number(item["maximum_log_composition_distance"], "composition distance"),
        _number(item["maximum_log_volume_distance"], "volume distance"),
        _text(item["anchor_state_identity"], "anchor state"),
        _text(item["anchor_problem_fingerprint"], "anchor problem fingerprint"),
        _text(item["anchor_topology_fingerprint"], "anchor topology fingerprint"),
        _text(item["anchor_parameter_fingerprint"], "anchor parameter fingerprint"),
    )


def _state_role(value: object, name: str) -> Fixed | Solved:
    item = cast(Mapping[str, object], value) if isinstance(value, Mapping) else {}
    role = item.get("role")
    if role == "fixed":
        fixed = _mapping(value, frozenset({"role", "value", "unit"}), f"GREPE {name}")
        return Fixed(
            epcsaft.unit_registry.Quantity(
                _number(fixed["value"], f"{name} value"),
                _text(fixed["unit"], f"{name} unit"),
            )
        )
    if role == "solved":
        solved = _mapping(
            value,
            frozenset({"role", "initial", "bounds", "starts", "unit"}),
            f"GREPE {name}",
        )
        unit = _text(solved["unit"], f"{name} unit")
        return Solved(
            epcsaft.unit_registry.Quantity(_number(solved["initial"], f"{name} initial"), unit),
            cast(
                tuple[Any, Any],
                tuple(
                    epcsaft.unit_registry.Quantity(item, unit)
                    for item in _floats(solved["bounds"], f"{name} bounds")
                ),
            ),
            tuple(
                epcsaft.unit_registry.Quantity(item, unit)
                for item in _floats(solved["starts"], f"{name} starts")
            ),
        )
    raise ValueError(f"GREPE {name} role must be fixed or solved")


def _phase(value: object) -> ReactivePhase:
    required = frozenset(
        {
            "identity",
            "fluid_role",
            "amount_role",
            "support",
            "model",
            "continuation",
        }
    )
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise TypeError("GREPE phase must be an object with string keys")
    if set(value) not in {required, required | {"start"}}:
        raise ValueError(
            "GREPE phase fields differ: expected required phase fields and optional start"
        )
    item = cast(Mapping[str, object], value)
    support_item = _mapping(
        item["support"], frozenset({"kind", "component_ids"}), "GREPE phase support"
    )
    support_kind = _text(support_item["kind"], "phase support kind")
    component_ids = tuple(
        _text(value, "support component")
        for value in _sequence(support_item["component_ids"], "support components")
    )
    if support_kind == "all_components" and not component_ids:
        support: PhaseSupport = AllComponents()
    elif support_kind == "neutral_components" and not component_ids:
        support = NeutralComponents()
    elif support_kind == "declared_neutral_components":
        support = DeclaredNeutralComponents(component_ids)
    else:
        raise ValueError("GREPE phase support declaration is inconsistent")
    model_item = _mapping(
        item["model"],
        frozenset({"kind", "reference_id", "admissible_packing_fraction_interval"}),
        "GREPE phase model",
    )
    model_kind = _text(model_item["kind"], "phase model kind")
    interval = model_item["admissible_packing_fraction_interval"]
    if model_kind == "provider":
        bounds = _floats(interval, "packing interval")
        if len(bounds) != 2:
            raise ValueError("provider packing interval requires two bounds")
        model: ProviderModel | IdealGasModel = ProviderModel(
            bounds, _text(model_item["reference_id"], "model reference")
        )
    elif model_kind == "ideal_gas" and interval is None:
        model = IdealGasModel(_text(model_item["reference_id"], "model reference"))
    else:
        raise ValueError("GREPE phase model declaration is inconsistent")
    continuation_value = item["continuation"]
    if continuation_value is None:
        continuation_identity = None
        branch_policy = None
        continuation_reference = None
    else:
        continuation = _mapping(
            continuation_value,
            frozenset({"identity", "branch_policy", "reference"}),
            "GREPE phase continuation",
        )
        continuation_identity = _text(continuation["identity"], "continuation identity")
        branch_policy = cast(Any, _text(continuation["branch_policy"], "branch policy"))
        continuation_reference = _homogeneous_reference(continuation["reference"])
    start_value = item.get("start")
    if start_value is None:
        start = None
    else:
        start_item = _mapping(
            start_value,
            frozenset({"amounts_mol", "molar_volume_m3_per_mol"}),
            "GREPE finite phase start",
        )
        start = FinitePhaseStart(
            _floats(start_item["amounts_mol"], "finite phase start amounts"),
            _number(start_item["molar_volume_m3_per_mol"], "finite phase start molar volume"),
        )
    return ReactivePhase(
        _text(item["identity"], "phase identity"),
        cast(Any, _text(item["fluid_role"], "fluid role")),
        cast(Any, _text(item["amount_role"], "amount role")),
        support,
        model,
        continuation_identity,
        branch_policy,
        continuation_reference,
        start,
    )


def _output(value: object) -> EquilibriumOutput:
    item = _mapping(
        value,
        frozenset(
            {
                "identity",
                "selector",
                "unit",
                "basis",
                "phase_identity",
                "coefficients",
                "solvent_mass_coefficients_kg_per_mol",
                "support",
                "censor_limit",
                "aggregate_identity",
                "covariance_identity",
            }
        ),
        "GREPE output",
    )
    return EquilibriumOutput(
        _text(item["identity"], "output identity"),
        cast(Any, _text(item["selector"], "output selector")),
        _text(item["unit"], "output unit"),
        _text(item["basis"], "output basis"),
        _optional_text(item["phase_identity"], "output phase"),
        _floats(item["coefficients"], "output coefficients"),
        _floats(item["solvent_mass_coefficients_kg_per_mol"], "solvent mass coefficients"),
        cast(Any, _text(item["support"], "output support")),
        None if item["censor_limit"] is None else _number(item["censor_limit"], "censor limit"),
        _optional_text(item["aggregate_identity"], "aggregate identity"),
        _optional_text(item["covariance_identity"], "covariance identity"),
    )


def general_reactive_equilibrium_problem_from_mapping(
    value: object,
) -> GeneralReactiveEquilibriumProblem:
    """Load one strict JSON-compatible, role-declared GREPE request."""

    item = _mapping(
        value,
        frozenset(
            {
                "identity",
                "temperature",
                "pressure",
                "phases",
                "reaction_system",
                "reaction_phase_ids",
                "outputs",
                "continuation",
            }
        ),
        "GREPE request",
    )
    continuation_value = item["continuation"]
    if continuation_value is None:
        continuation_identity = None
        continuation_reference = None
    else:
        continuation = _mapping(
            continuation_value,
            frozenset({"identity", "reference"}),
            "GREPE coupled continuation",
        )
        continuation_identity = _text(continuation["identity"], "coupled continuation identity")
        continuation_reference = _bubble_reference(continuation["reference"])
    return GeneralReactiveEquilibriumProblem(
        _text(item["identity"], "GREPE identity"),
        _state_role(item["temperature"], "temperature"),
        _state_role(item["pressure"], "pressure"),
        tuple(_phase(raw) for raw in _sequence(item["phases"], "GREPE phases")),
        _reaction_system(item["reaction_system"]),
        tuple(
            _text(raw, "reaction phase identity")
            for raw in _sequence(item["reaction_phase_ids"], "reaction phase identities")
        ),
        tuple(_output(raw) for raw in _sequence(item["outputs"], "GREPE outputs")),
        continuation_identity,
        continuation_reference,
    )


__all__ = ("general_reactive_equilibrium_problem_from_mapping",)
