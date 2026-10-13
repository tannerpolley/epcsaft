from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol, cast

from .. import _core


class _Species(Protocol):
    component_ids: tuple[str, ...]


class _Model(_Species, Protocol):
    _model: _core.HelmholtzModel


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Free:
    guess: float | None = None


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Pinned:
    value: float


class Phase:
    __slots__ = ("name", "support", "amount", "kind", "composition_guess", "packing_guess")
    name: str
    support: tuple[str, ...] | None
    amount: Free | Pinned
    kind: str | None
    composition_guess: tuple[float, ...] | dict[str, float]
    packing_guess: float | None

    def __init__(
        self,
        name: str,
        *,
        support: Sequence[str] | None = None,
        amount: Free | Pinned | None = None,
        kind: str | None = None,
        composition_guess: Sequence[float] | Mapping[str, float] | None = None,
        packing_guess: float | None = None,
    ) -> None:
        self.name = name
        self.support = None if support is None else tuple(support)
        self.amount = Free() if amount is None else amount
        self.kind = kind
        guess = () if composition_guess is None else composition_guess  # a mapping is ordered at compile time
        self.composition_guess = dict(guess) if isinstance(guess, Mapping) else tuple(guess)
        self.packing_guess = packing_guess


class Discover:
    __slots__ = ("trials", "declaration_budget", "trial_budget")

    def __init__(self, *, trials: Sequence[Phase] | None = None,
                 declaration_budget: int | None = None, trial_budget: int | None = None) -> None:
        self.trials = None if trials is None else tuple(trials)
        self.declaration_budget = declaration_budget
        self.trial_budget = trial_budget


class Feed:
    __slots__ = ("basis", "values", "solvent", "solvent_mass_kg")
    basis: _core.FeedKind
    values: Mapping[str, float]
    solvent: tuple[str, ...]
    solvent_mass_kg: float | None

    def __init__(
        self,
        basis: _core.FeedKind,
        values: Mapping[str, float],
        *,
        solvent_mass_kg: float | None,
        solvent: Sequence[str] | None = None,
    ) -> None:
        solvent = () if solvent is None else solvent
        self.basis = basis
        self.values = values
        self.solvent = tuple(solvent)
        self.solvent_mass_kg = solvent_mass_kg


class MoleFractions(Feed):
    def __init__(self, values: Mapping[str, float]) -> None:
        super().__init__(_core.FeedKind.MoleFractions, values, solvent_mass_kg=None)


class Amounts(Feed):
    def __init__(self, values: Mapping[str, float]) -> None:
        super().__init__(_core.FeedKind.Amounts, values, solvent_mass_kg=None)


class MassFractions(Feed):
    def __init__(self, values: Mapping[str, float]) -> None:
        super().__init__(_core.FeedKind.MassFractions, values, solvent_mass_kg=None)


class Molalities(Feed):
    def __init__(
        self,
        values: Mapping[str, float],
        *,
        solvent: Sequence[str],
        solvent_mass_kg: float,
    ) -> None:
        super().__init__(_core.FeedKind.Molalities, values, solvent=solvent, solvent_mass_kg=solvent_mass_kg)


class NeutralReference:
    """Pure neutral solvent and an explicit charge-neutral molality dilution path."""
    __slots__ = ("solvent", "solute_molality_ratios", "standard_molality_mol_per_kg", "species_reference")

    def __init__(self, *, solvent: str, solute_molality_ratios: Mapping[str, float],
                 standard_molality_mol_per_kg: float,
                 species_reference: Mapping[str, _core.SpeciesReference] | None = None) -> None:
        self.solvent = solvent
        self.solute_molality_ratios = dict(solute_molality_ratios)
        self.standard_molality_mol_per_kg = standard_molality_mol_per_kg
        self.species_reference = dict(species_reference or {})

    def _native(self, model: _Species) -> _core.NeutralReferenceDeclaration:
        item = _core.NeutralReferenceDeclaration()
        item.solvent_index = _index(model.component_ids, self.solvent)
        item.solute_molality_ratios = _vector(self.solute_molality_ratios, model.component_ids)
        item.standard_molality_mol_per_kg = self.standard_molality_mol_per_kg
        if self.species_reference:
            choices = [_core.SpeciesReference.InfiniteDilution] * len(model.component_ids)
            choices[item.solvent_index] = _core.SpeciesReference.PureLiquid
            for species, choice in self.species_reference.items():
                choices[_index(model.component_ids, species)] = choice
            item.species_reference = choices
        return item


class Reaction:
    __slots__ = ("name", "stoichiometry", "equilibrium_constant", "affinity_j_per_mol", "correlation", "reference")
    name: str
    stoichiometry: Mapping[str, float]
    equilibrium_constant: float | None
    affinity_j_per_mol: float | None
    correlation: _core.ReactionLogPolynomial | None
    reference: _core.ReactionReference | None

    def __init__(
        self,
        stoichiometry: Mapping[str, float],
        *,
        equilibrium_constant: float | None = None,
        name: str | None = None,
        affinity_j_per_mol: float | None = None,
        correlation: _core.ReactionLogPolynomial | None = None,
        reference: _core.ReactionReference | None = None,
    ) -> None:
        """Exactly one target; source-basis correlations also require an explicit reference."""
        if sum(item is not None for item in (equilibrium_constant, affinity_j_per_mol, correlation)) != 1:
            raise _core.ParameterError("reaction target must be one of equilibrium_constant or affinity_j_per_mol or correlation")
        name = "" if name is None else name
        self.name = name
        self.stoichiometry = stoichiometry
        self.equilibrium_constant = equilibrium_constant
        self.affinity_j_per_mol = affinity_j_per_mol
        self.correlation = correlation
        self.reference = reference


def _index(component_ids: tuple[str, ...], name: str) -> int:
    try:
        return component_ids.index(name)
    except ValueError as error:
        raise _core.ParameterError(f"unknown component_id {name}") from error


def _vector(values: Mapping[str, float], component_ids: tuple[str, ...]) -> list[float]:
    if unknown := [name for name in values if name not in component_ids]:
        raise _core.ParameterError(f"unknown species {unknown[0]}; expected one of {', '.join(component_ids)}")
    return [values.get(component_id, 0.0) for component_id in component_ids]


class Problem:
    __slots__ = (
        "phases", "T", "P", "feed", "reactions", "standard_gibbs_j_per_mol", "topology_declared",
        "row_tolerance", "neutral_reference", "enthalpy_j", "entropy_j_per_k", "critical", "envelope",
    )
    phases: tuple[Phase, ...] | Discover
    T: float | Free
    P: float | Free
    feed: Feed
    reactions: tuple[Reaction, ...]
    standard_gibbs_j_per_mol: tuple[float, ...]
    topology_declared: bool
    row_tolerance: float | None
    neutral_reference: NeutralReference | None
    enthalpy_j: float | None
    entropy_j_per_k: float | None
    critical: _core.CriticalDeclaration | None
    envelope: _core.EnvelopeCondition | None

    def __init__(
        self,
        *,
        phases: Sequence[Phase] | Discover | None = None,
        T: float | Free,
        P: float | Free,
        feed: Feed,
        reactions: Sequence[Reaction] | None = None,
        standard_gibbs_j_per_mol: Sequence[float] | None = None,
        topology_declared: bool | None = None,
        row_tolerance: float | None = None,
        enthalpy_j: float | None = None,
        entropy_j_per_k: float | None = None,
        neutral_reference: NeutralReference | None = None,
        critical: _core.CriticalDeclaration | None = None,
        envelope: _core.EnvelopeCondition | None = None,
    ) -> None:
        self.critical = critical
        self.envelope = envelope
        self.row_tolerance = row_tolerance
        self.enthalpy_j = enthalpy_j
        self.entropy_j_per_k = entropy_j_per_k
        self.neutral_reference = neutral_reference
        phases = () if phases is None else phases
        reactions = () if reactions is None else reactions
        standard_gibbs_j_per_mol = () if standard_gibbs_j_per_mol is None else standard_gibbs_j_per_mol
        topology_declared = True if topology_declared is None else topology_declared
        self.phases = phases if isinstance(phases, Discover) else tuple(phases)
        self.T = T
        self.P = P
        self.feed = feed
        self.reactions = tuple(reactions)
        self.standard_gibbs_j_per_mol = tuple(standard_gibbs_j_per_mol)
        self.topology_declared = topology_declared and not isinstance(phases, Discover)

    def _native(self, model: _Species) -> _core.EquilibriumProblem:
        native = _core.EquilibriumProblem()
        native.topology_declared = self.topology_declared
        native.feed.kind = self.feed.basis
        native.feed.values = _vector(self.feed.values, model.component_ids)
        native.feed.solvent = [_index(model.component_ids, name) for name in self.feed.solvent]
        if self.feed.solvent_mass_kg is not None:
            native.feed.solvent_mass_kg = self.feed.solvent_mass_kg
        native.standard_gibbs_j_per_mol = list(self.standard_gibbs_j_per_mol)
        native.enthalpy_j = self.enthalpy_j
        native.entropy_j_per_k = self.entropy_j_per_k
        native.critical = self.critical
        native.envelope = self.envelope
        native.temperature = _free_value(self.T)
        native.pressure = _free_value(self.P)
        native.neutral_reference = None if self.neutral_reference is None else self.neutral_reference._native(model)
        if self.row_tolerance is not None:
            native.row_tolerance = self.row_tolerance
        native_phases: list[_core.PhaseSpec] = []
        phases = (() if self.phases.trials is None else self.phases.trials) if isinstance(self.phases, Discover) else self.phases
        for phase in phases:
            item = _core.PhaseSpec()
            item.name = phase.name
            item.kind = _phase_kind(phase)
            item.support_declared = phase.support is not None
            names = () if phase.support is None else phase.support
            item.support = sorted(_index(model.component_ids, name) for name in names)  # native order
            item.amount = _amount_spec(phase.amount)
            guess = phase.composition_guess
            order = tuple(model.component_ids[i] for i in item.support) or model.component_ids
            item.composition_guess = _vector(guess, order) if isinstance(guess, dict) else list(guess)
            item.packing_guess = phase.packing_guess
            native_phases.append(item)
        if isinstance(self.phases, Discover):
            options = _core.StabilityOptions()
            options.trials = None if self.phases.trials is None else native_phases
            options.declaration_budget = self.phases.declaration_budget
            options.trial_budget = self.phases.trial_budget
            native.discovery = options
            native.phases = []
        else:
            native.phases = native_phases
        native_reactions: list[_core.ReactionSpec] = []
        for reaction in self.reactions:
            reaction_item = _core.ReactionSpec()
            reaction_item.name = reaction.name
            reaction_item.stoichiometry = _vector(reaction.stoichiometry, model.component_ids)
            reaction_item.equilibrium_constant = reaction.equilibrium_constant
            reaction_item.affinity_j_per_mol = reaction.affinity_j_per_mol
            reaction_item.correlation = reaction.correlation
            reaction_item.reference = reaction.reference
            native_reactions.append(reaction_item)
        native.reactions = native_reactions
        return native


def problem_from_mapping(mapping: Mapping[str, object]) -> Problem:
    allowed = {
        "phases", "T", "P", "feed", "reactions", "standard_gibbs_j_per_mol", "topology_declared",
        "row_tolerance", "neutral_reference", "enthalpy_j", "entropy_j_per_k",
    }
    if set(mapping) - allowed:
        raise _core.ParameterError("unknown Problem key")
    phases_value = mapping.get("phases")
    if isinstance(phases_value, Discover):
        phases: Sequence[Phase] | Discover = phases_value
    elif phases_value == "discover":
        phases = Discover()
    else:
        phases = tuple(_phase_from_mapping(item) for item in _records(phases_value, "phases"))
    feed = _feed_from_mapping(mapping.get("feed"))
    reactions = tuple(
        _reaction_from_mapping(item) for item in _records(mapping.get("reactions", ()), "reactions")
    )
    topology = mapping.get("topology_declared")
    if not isinstance(topology, bool | None):
        raise _core.ParameterError(f"topology_declared must be a bool, got {topology!r}")
    return Problem(
        phases=phases,
        T=_intensive_from_mapping(mapping.get("T"), "T"),
        P=_intensive_from_mapping(mapping.get("P"), "P"),
        feed=feed,
        reactions=reactions,
        standard_gibbs_j_per_mol=tuple(cast(Sequence[float], mapping.get("standard_gibbs_j_per_mol", ()))),
        topology_declared=topology,
        row_tolerance=cast(float | None, mapping.get("row_tolerance")),
        enthalpy_j=cast(float | None, mapping.get("enthalpy_j")),
        entropy_j_per_k=cast(float | None, mapping.get("entropy_j_per_k")),
        neutral_reference=_neutral_reference_from_mapping(mapping.get("neutral_reference")),
    )


def _records(value: object, name: str) -> Sequence[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise _core.ParameterError(f"{name} must be a sequence")
    return tuple(_record(item, name) for item in value)


def _record(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise _core.ParameterError(f"{name} entries must be mappings")
    return value


def _neutral_reference_from_mapping(value: object) -> NeutralReference | None:
    if value is None or isinstance(value, NeutralReference):
        return value
    record = _record(value, "neutral_reference")
    solvent, ratios = record.get("solvent"), record.get("solute_molality_ratios")
    if not isinstance(solvent, str) or not isinstance(ratios, Mapping) or "standard_molality_mol_per_kg" not in record:
        raise _core.ParameterError("neutral_reference requires solvent, solute_molality_ratios and standard_molality_mol_per_kg")
    return NeutralReference(solvent=solvent, solute_molality_ratios=ratios,
                            standard_molality_mol_per_kg=cast(float, record["standard_molality_mol_per_kg"]),
                            species_reference=cast(Mapping[str, _core.SpeciesReference] | None, record.get("species_reference")))


def _reaction_reference_from_mapping(value: object) -> _core.ReactionReference | None:
    if value is None or isinstance(value, _core.ReactionReference):
        return value
    record = _record(value, "reaction reference")
    basis = record.get("source_basis")
    if isinstance(basis, str):
        basis = getattr(_core.ReferenceSourceBasis, basis, None)
    if not isinstance(basis, _core.ReferenceSourceBasis):
        raise _core.ParameterError("reaction reference requires an explicit source_basis")
    return _core.ReactionReference(basis, cast(float | None, record.get("reference_pressure_pa")))


def _phase_from_mapping(value: Mapping[str, object]) -> Phase:
    amount_value = {"status": "free"} if value.get("amount") is None else value.get("amount")  # omitted: free
    if not isinstance(amount_value, Mapping):
        raise _core.ParameterError(f"phase {value.get('name')} amount must be a mapping, got {amount_value!r}")
    status, guess = amount_value.get("status"), amount_value.get("guess")
    if status == "pinned" and amount_value.get("value") is not None:
        amount: Free | Pinned = Pinned(float(amount_value["value"]))
    elif status == "free":
        amount = Free(None if guess is None else float(guess))
    else:
        raise _core.ParameterError(f"pinned amount of phase {value.get('name')} requires a value" if status == "pinned" else f"unknown amount status {status}")
    return Phase(
        str(value.get("name", "")),
        support=cast(Sequence[str] | None, value.get("support")),
        amount=amount,
        kind=cast(str | None, value.get("kind")),
        composition_guess=cast(Sequence[float], value.get("composition_guess", ())),
        packing_guess=cast(float | None, value.get("packing_guess")),
    )


def _intensive_from_mapping(value: object, name: str) -> float | Free:
    if isinstance(value, Mapping) and value.get("free") is True:
        guess = value.get("guess")
        return Free(None if guess is None else float(guess))
    if isinstance(value, Mapping):
        value = value.get("value")  # a fixed intensive value is never defaulted
    if isinstance(value, (int, float)):
        return float(value)
    raise _core.ParameterError(f"missing {name} value")


def _feed_from_mapping(value: object) -> Feed:
    record = _record(value, "feed")
    values = cast(Mapping[str, float], record.get("values"))
    basis = record.get("basis")
    if basis == "mole_fractions":
        return MoleFractions(values)
    if basis == "amounts":
        return Amounts(values)
    if basis == "mass_fractions":
        return MassFractions(values)
    if basis == "molalities":
        if "solvent_mass_kg" not in record:
            raise _core.ParameterError("missing solvent_mass_kg")
        return Molalities(
            values,
            solvent=cast(Sequence[str], record.get("solvent", ())),
            solvent_mass_kg=cast(float, record["solvent_mass_kg"]),
        )
    raise _core.ParameterError(f"unknown feed basis {basis}; expected mole_fractions/amounts/mass_fractions/molalities")


def _reaction_from_mapping(value: Mapping[str, object]) -> Reaction:
    return Reaction(
        cast(Mapping[str, float], value.get("stoichiometry")),
        name=str(value.get("name", "")),
        equilibrium_constant=cast(float | None, value.get("equilibrium_constant")),
        affinity_j_per_mol=cast(float | None, value.get("affinity_j_per_mol")),
        correlation=cast(_core.ReactionLogPolynomial | None, value.get("correlation")),
        reference=_reaction_reference_from_mapping(value.get("reference")),
    )


def _free_value(value: float | Free) -> _core.FreeValue:
    if not isinstance(value, Free):
        return _core.FreeValue.fixed(value)
    return _core.FreeValue.variable() if value.guess is None else _core.FreeValue.variable(value.guess)


def _amount_spec(value: Free | Pinned) -> _core.AmountSpec:
    return _core.AmountSpec.free() if isinstance(value, Free) else _core.AmountSpec.pinned(value.value)


def _phase_kind(phase: Phase) -> _core.PhaseKind:
    kind = phase.kind
    if kind is None:
        raise _core.ParameterError(f"phase {phase.name} requires an explicit kind (liquid or vapor)")
    if kind == "vapor":
        return _core.PhaseKind.Vapor
    if kind == "liquid":
        return _core.PhaseKind.Liquid
    raise _core.ParameterError(f"unknown phase kind {kind}; expected liquid or vapor")


def solve_equilibrium(
    model: _Model,
    problem: Problem | _core.EquilibriumProblem,
    start: _core.EquilibriumResult | None = None,
) -> _core.EquilibriumResult:
    native = problem._native(model) if isinstance(problem, Problem) else problem
    return _core.solve_equilibrium(model._model, native, start)


def compile_problem(model: _Model, problem: Problem | _core.EquilibriumProblem) -> _core.CompiledSystem:
    """Compile once for solved-state actions; the result keeps the native model alive."""
    native = problem._native(model) if isinstance(problem, Problem) else problem
    return _core.compile_problem(model._model, native)


solved_state_actions = _core.solved_state_actions


Continuation = _core.ContinuationOptions


def continue_equilibrium(
    model: _Model,
    start: Problem,
    target: Problem,
    options: Continuation,
    start_result: _core.EquilibriumResult | None = None,
    *,
    parameter_targets: Sequence[float] | None = None,
) -> _core.ContinuationResult:
    request = _core.ContinuationRequest()
    request.start_problem = start._native(model)
    request.target_problem = target._native(model)
    request.start_result = start_result
    request.parameter_targets = () if parameter_targets is None else parameter_targets
    request.options = options
    return _core.continue_equilibrium(model._model, request)


def bubble_point(*, liquid: Phase, T: float, feed: Feed, reactions: Sequence[Reaction] | None = None) -> Problem:
    reactions = () if reactions is None else reactions
    return Problem(
        phases=[liquid, Phase("V", kind="vapor", amount=Pinned(0.0))],
        T=T,
        P=Free(),
        feed=feed,
        reactions=reactions,
    )


def dew_point(*, vapor: Phase, T: float, feed: Feed, reactions: Sequence[Reaction] | None = None) -> Problem:
    reactions = () if reactions is None else reactions
    return Problem(
        phases=[Phase("L", kind="liquid", amount=Pinned(0.0)), vapor],
        T=T,
        P=Free(),
        feed=feed,
        reactions=reactions,
    )


def solve_equilibrium_starts(model: _Model, problem: Problem, starts: Sequence[_core.EquilibriumStart],
                             limits: _core.EquilibriumWorkLimits) -> _core.EquilibriumStartsResult:
    return _core.solve_equilibrium_starts(model._model, problem._native(model), list(starts), limits)
