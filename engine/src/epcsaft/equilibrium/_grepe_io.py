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
    DeclaredComponents,
    DeclaredNeutralComponents,
    NeutralComponents,
)
from ._feed import ResolvedFeed
from ._grepe import (
    EquilibriumOutput,
    FinitePhaseStart,
    Fixed,
    GeneralReactiveEquilibriumProblem,
    IdealGasModel,
    ProviderModel,
    ReactivePhase,
    ReactivePhaseSupport,
    Solved,
)
from ._observations import HomogeneousContinuationReference
from ._phase_reactions import (
    PhaseReactionRow,
    PhaseReactionStandardState,
    PhaseReactionSystem,
    PhaseSpeciesTerm,
)
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


def _mapping_with_optional(
    value: object,
    required: frozenset[str],
    optional: frozenset[str],
    label: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or any(type(key) is not str for key in value):
        raise TypeError(f"{label} must be an object with string keys")
    fields = set(value)
    if not required <= fields or fields - required - optional:
        raise ValueError(
            f"{label} fields differ: required {tuple(sorted(required))}, "
            f"optional {tuple(sorted(optional))}, received {tuple(sorted(fields))}"
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


def _phase_reaction_system(value: object) -> PhaseReactionSystem:
    item = _mapping(
        value,
        frozenset({"phase_ids", "species_ids", "charges", "balance_matrix", "reactions"}),
        "GREPE phase reactions",
    )
    reactions = []
    for row_index, raw_row in enumerate(_sequence(item["reactions"], "phase reaction rows")):
        row = _mapping(
            raw_row,
            frozenset(
                {
                    "identity",
                    "terms",
                    "ln_equilibrium_constant",
                    "standard_state",
                }
            ),
            f"phase reaction row {row_index}",
        )
        terms = tuple(
            PhaseSpeciesTerm(
                _text(term[0], "phase reaction term phase"),
                _text(term[1], "phase reaction term species"),
                _number(term[2], "phase reaction coefficient"),
            )
            for term in (
                _sequence(raw_term, "phase reaction term")
                for raw_term in _sequence(row["terms"], "phase reaction terms")
            )
            if len(term) == 3
        )
        if len(terms) != len(_sequence(row["terms"], "phase reaction terms")):
            raise ValueError("phase reaction terms require phase, species, and coefficient")
        standard_value = row["standard_state"]
        standard = None
        if standard_value is not None:
            standard_mapping = _mapping(
                standard_value,
                frozenset({"activity_scale_id", "reference_pressure_pa"}),
                "phase reaction standard state",
            )
            standard = PhaseReactionStandardState(
                _text(standard_mapping["activity_scale_id"], "activity scale identity"),
                (
                    None
                    if standard_mapping["reference_pressure_pa"] is None
                    else _number(
                        standard_mapping["reference_pressure_pa"],
                        "reference pressure",
                    )
                ),
            )
        reactions.append(
            PhaseReactionRow(
                _text(row["identity"], "phase reaction identity"),
                terms,
                (
                    None
                    if row["ln_equilibrium_constant"] is None
                    else _number(
                        row["ln_equilibrium_constant"],
                        "phase reaction equilibrium constant",
                    )
                ),
                standard,
            )
        )
    return PhaseReactionSystem(
        tuple(
            _text(raw, "phase identity") for raw in _sequence(item["phase_ids"], "phase identities")
        ),
        tuple(
            _text(raw, "species identity")
            for raw in _sequence(item["species_ids"], "species identities")
        ),
        tuple(reactions),
        _floats(item["charges"], "phase reaction charges"),
        _rows(item["balance_matrix"], "phase reaction balances"),
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
        support: ReactivePhaseSupport = AllComponents()
    elif support_kind == "neutral_components" and not component_ids:
        support = NeutralComponents()
    elif support_kind == "declared_neutral_components":
        support = DeclaredNeutralComponents(component_ids)
    elif support_kind == "declared_components":
        support = DeclaredComponents(component_ids)
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

    item = _mapping_with_optional(
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
        frozenset({"feed", "phase_reactions"}),
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
    feed_value = item.get("feed")
    feed = None
    if feed_value is not None:
        feed_mapping = _mapping(
            feed_value,
            frozenset(
                {
                    "component_ids",
                    "amounts_mol",
                    "coordinate_ids",
                    "amount_jacobian_mol_per_coordinate",
                }
            ),
            "GREPE resolved feed",
        )
        feed = ResolvedFeed(
            tuple(
                _text(raw, "feed component identity")
                for raw in _sequence(feed_mapping["component_ids"], "feed component identities")
            ),
            _floats(feed_mapping["amounts_mol"], "feed amounts"),
            tuple(
                _text(raw, "feed coordinate identity")
                for raw in _sequence(feed_mapping["coordinate_ids"], "feed coordinate identities")
            ),
            _rows(
                feed_mapping["amount_jacobian_mol_per_coordinate"],
                "feed amount Jacobian",
            ),
        )
    phase_reactions_value = item.get("phase_reactions")
    phase_reactions = (
        None if phase_reactions_value is None else _phase_reaction_system(phase_reactions_value)
    )
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
        feed,
        phase_reactions,
    )


def general_reactive_equilibrium_problem_to_mapping(
    problem: GeneralReactiveEquilibriumProblem,
) -> dict[str, object]:
    """Serialize one typed GREPE request for regression or JSON transport."""

    if not isinstance(problem, GeneralReactiveEquilibriumProblem):
        raise TypeError("problem must be a GeneralReactiveEquilibriumProblem")

    def state_role(role: Fixed | Solved, unit: str) -> dict[str, object]:
        if isinstance(role, Fixed):
            return {
                "role": "fixed",
                "value": float(role.value.to(unit).magnitude),
                "unit": unit,
            }
        return {
            "role": "solved",
            "initial": float(role.initial.to(unit).magnitude),
            "bounds": [float(value.to(unit).magnitude) for value in role.bounds],
            "starts": [float(value.to(unit).magnitude) for value in role.starts],
            "unit": unit,
        }

    def phase_mapping(phase: ReactivePhase) -> dict[str, object]:
        if isinstance(phase.support, AllComponents):
            support = {"kind": "all_components", "component_ids": []}
        elif isinstance(phase.support, NeutralComponents):
            support = {"kind": "neutral_components", "component_ids": []}
        elif isinstance(phase.support, DeclaredNeutralComponents):
            support = {
                "kind": "declared_neutral_components",
                "component_ids": list(phase.support.component_ids),
            }
        else:
            support = {
                "kind": "declared_components",
                "component_ids": list(phase.support.component_ids),
            }
        if isinstance(phase.model, ProviderModel):
            model: dict[str, object] = {
                "kind": "provider",
                "reference_id": phase.model.reference_id,
                "admissible_packing_fraction_interval": list(
                    phase.model.admissible_packing_fraction_interval
                ),
            }
        else:
            model = {
                "kind": "ideal_gas",
                "reference_id": phase.model.reference_id,
                "admissible_packing_fraction_interval": None,
            }
        continuation = (
            None
            if phase.continuation_identity is None
            else {
                "identity": phase.continuation_identity,
                "branch_policy": phase.branch_policy,
                "reference": (
                    None
                    if phase.continuation_reference is None
                    else phase.continuation_reference.to_mapping()
                ),
            }
        )
        start = (
            None
            if phase.start is None
            else {
                "amounts_mol": list(phase.start.amounts_mol),
                "molar_volume_m3_per_mol": phase.start.molar_volume_m3_per_mol,
            }
        )
        return {
            "identity": phase.identity,
            "fluid_role": phase.fluid_role,
            "amount_role": phase.amount_role,
            "support": support,
            "model": model,
            "continuation": continuation,
            "start": start,
        }

    reaction = problem.reaction_system
    standard = reaction.source_standard_state
    standard_mapping = (
        None
        if standard is None
        else {
            "id": standard.id,
            "activity_scale_id": standard.activity_scale_id,
            "log_activity_scale_factors": list(standard.log_activity_scale_factors),
            "reference_pressure_pa": standard.reference_pressure_pa,
            "component_ids": list(standard.source_reference_component_ids),
            "solvent_composition": list(standard.source_reference_solvent_composition),
            "ion_pairs": [
                [list(pair.component_ids), list(pair.stoichiometry)]
                for pair in standard.source_reference_ion_pairs
            ],
            "phase": standard.source_reference_phase,
            "reference_convention": standard.source_reference_convention,
            "activity_convention_id": (standard.source_reference_activity_convention_id),
            "standard_molality_mol_per_kg": (
                standard.source_reference_standard_molality_mol_per_kg
            ),
        }
    )
    return {
        "identity": problem.identity,
        "temperature": state_role(problem.temperature, "kelvin"),
        "pressure": state_role(problem.pressure, "pascal"),
        "phases": [phase_mapping(phase) for phase in problem.phases],
        "reaction_system": {
            "species_ids": list(reaction.species_ids),
            "charges": list(reaction.charges),
            "molar_masses_kg_per_mol": list(reaction.molar_masses_kg_per_mol),
            "balance_matrix": [list(row) for row in reaction.balance_matrix],
            "conserved_totals": list(reaction.conserved_totals),
            "reaction_matrix": [list(row) for row in reaction.reaction_matrix],
            "feed_amounts_mol": list(reaction.feed_amounts_mol),
            "equilibrium_constants": [
                [
                    record.ln_value,
                    record.source_id,
                    record.reference_id,
                    record.reaction_orientation,
                    record.conversion_id,
                    record.dimensionless,
                ]
                for record in reaction.equilibrium_constants
            ],
            "strict_interior_amount_floor_mol": (reaction.strict_interior_amount_floor_mol),
            "source_standard_state": standard_mapping,
        },
        "reaction_phase_ids": list(problem.reaction_phase_ids),
        "outputs": [
            {
                "identity": output.identity,
                "selector": output.selector,
                "unit": output.unit,
                "basis": output.basis,
                "phase_identity": output.phase_identity,
                "coefficients": list(output.coefficients),
                "solvent_mass_coefficients_kg_per_mol": list(
                    output.solvent_mass_coefficients_kg_per_mol
                ),
                "support": output.support,
                "censor_limit": output.censor_limit,
                "aggregate_identity": output.aggregate_identity,
                "covariance_identity": output.covariance_identity,
            }
            for output in problem.outputs
        ],
        "continuation": (
            None
            if problem.continuation_identity is None
            else {
                "identity": problem.continuation_identity,
                "reference": (
                    None
                    if problem.continuation_reference is None
                    else problem.continuation_reference.to_mapping()
                ),
            }
        ),
        "feed": (
            None
            if problem.feed is None
            else {
                "component_ids": list(problem.feed.component_ids),
                "amounts_mol": list(problem.feed.amounts_mol),
                "coordinate_ids": list(problem.feed.coordinate_ids),
                "amount_jacobian_mol_per_coordinate": [
                    list(row) for row in problem.feed.amount_jacobian_mol_per_coordinate
                ],
            }
        ),
        "phase_reactions": (
            None
            if problem.phase_reactions is None
            else {
                "phase_ids": list(problem.phase_reactions.phase_ids),
                "species_ids": list(problem.phase_reactions.species_ids),
                "charges": list(problem.phase_reactions.charges),
                "balance_matrix": [list(row) for row in problem.phase_reactions.balance_matrix],
                "reactions": [
                    {
                        "identity": reaction.identity,
                        "terms": [
                            [term.phase_id, term.species_id, term.coefficient]
                            for term in reaction.terms
                        ],
                        "ln_equilibrium_constant": reaction.ln_equilibrium_constant,
                        "standard_state": (
                            None
                            if reaction.standard_state is None
                            else {
                                "activity_scale_id": (reaction.standard_state.activity_scale_id),
                                "reference_pressure_pa": (
                                    reaction.standard_state.reference_pressure_pa
                                ),
                            }
                        ),
                    }
                    for reaction in problem.phase_reactions.reactions
                ],
            }
        ),
    }


__all__ = (
    "general_reactive_equilibrium_problem_from_mapping",
    "general_reactive_equilibrium_problem_to_mapping",
)
