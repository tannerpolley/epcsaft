"""Declare feed, phase and reaction inputs for the single native equilibrium driver.

Temperatures are K, pressure Pa, fluid amounts mol, and total energy
constraints J or J/K. Pinned phase amounts are fractions of the compiled
feed molar total, not absolute mol. Mapping keys are ordered into model
component order; unknown species are refused. Declarations do not solve
thermodynamic equations until compiled or evaluated by the native driver.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import KW_ONLY, dataclass, field
from typing import Literal, Protocol, cast

from .. import _core


class _Species(Protocol):
    component_ids: tuple[str, ...]


class _Model(_Species, Protocol):
    _model: _core.HelmholtzModel


class IdealMoleFraction:
    """Ideal liquid activities with explicit species order and dimensionless charges.

    Attributes
    ----------
    component_ids : tuple of str
        Species identifiers in native order.
    charge_number : tuple of int
        Dimensionless species charges.
    activity_basis : {"mole_fraction", "molarity"}
        Dimensionless ideal activity basis, defaulting to mole fractions.
    solution_volume_m3 : float or None
        Fixed positive solution volume required for molarity and absolute mol feeds.
    unit_activity_solvent : str or None
        Optional neutral species with unit activity on the molarity basis.
    """

    def __init__(
        self,
        component_ids: Sequence[str],
        charge_number: Sequence[int],
        *,
        solver_options: Mapping[str, float] | None = None,
        activity_basis: Literal["mole_fraction", "molarity"] | None = None,
        solution_volume_m3: float | None = None,
        unit_activity_solvent: str | None = None,
    ) -> None:
        self._model = _core.IdealMoleFraction(
            component_ids,
            charge_number,
            solver_options,
            activity_basis,
            solution_volume_m3,
            unit_activity_solvent,
        )
        self.activity_basis = self._model.activity_basis
        self.solution_volume_m3 = self._model.solution_volume_m3
        self.unit_activity_solvent = self._model.unit_activity_solvent
        self.component_ids = tuple(self._model.component_ids)
        self.charge_number = tuple(self._model.charge_number)


def reaction_activities(model: _Model | IdealMoleFraction, result: _core.EquilibriumResult) -> list[_core.ReactionActivities]:
    """Read native dimensionless species activities and ln Q at the returned state, on each declared basis."""
    return _core.reaction_activities(model._model, result)


def effective_water_amounts(
    model: IdealMoleFraction,
    result: _core.EquilibriumResult,
    *,
    water: str,
    ions: Sequence[str],
    molecular_species: Sequence[str],
) -> dict[str, float]:
    """Return native effective-water and preserved free molecular amounts in mol."""
    return _core.effective_water_amounts(model._model, result, water, ions, molecular_species)


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Free:
    """Declare a free temperature, pressure or phase amount.

    Attributes
    ----------
    guess : float or None
        Initial temperature in K or pressure in Pa when used for a free
        intensive. ``None`` leaves that guess absent. Phase-amount conversion
        refuses a non-``None`` guess; use ``Free()`` for phase amounts.

    Notes
    -----
    A free phase amount is solved in mol; a ``Free`` intensive uses the unit
    of its enclosing ``Problem`` field. Physical admissibility is checked by
    the native compiler, not this dataclass.

    Examples
    --------
    >>> from epcsaft.equilibrium import Free
    >>> Free(300.0).guess
    300.0
    >>> Free().guess is None
    True
    """

    guess: float | None = None


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Pinned:
    """Pin a phase amount to a fraction of the compiled feed molar total.

    Attributes
    ----------
    value : float
        Dimensionless fraction, admitted from zero to one for fluids.
        ``Pinned(0.0)`` declares an incipient fluid or an absent solid;
        solid declarations admit only the zero pin.

    Notes
    -----
    The native compiler forms the pinned molar amount as this fraction times
    the compiled feed total in mol. This dataclass performs no range check.

    Examples
    --------
    >>> from epcsaft.equilibrium import Pinned
    >>> Pinned(0.25).value
    0.25
    """

    value: float


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Phase:
    """Declare one liquid or vapor phase and its species support.

    Attributes
    ----------
    name : str
        Phase identifier used by observations and returned phase order.
    support : sequence of str or None
        Species with phase coordinates; omitted support is completed natively.
        Vapor completion excludes charged species. Explicit charged vapor requires
        DensityTemperature permittivity; an explicit empty support is refused.
    amount : Free or Pinned
        Free amount solved in mol, or pinned fraction of compiled feed mol.
        Default is ``Free()``; phase-amount guesses are not supported.
    kind : {'liquid', 'vapor'} or None
        Required physical branch when compiled; ``None`` is refused then.
    composition_guess : sequence or mapping
        Dimensionless mole-fraction guess in declared support's native component
        order; mappings are reordered and omitted supported species get zero.
    packing_guess : float or None
        Dimensionless initial packing fraction; ``None`` leaves it absent.

    Notes
    -----
    Sequences/mappings are copied at declaration. Native support, charge,
    branch and start admission follow ALG-EQUILIBRIUM. Declaring a branch does
    not establish its existence or global stability.

    Examples
    --------
    >>> from epcsaft.equilibrium import Phase
    >>> phase = Phase('L', kind='liquid', support=['water'])
    >>> phase.support
    ('water',)
    """

    name: str
    _: KW_ONLY
    support: Sequence[str] | None = None
    amount: Free | Pinned = field(default_factory=Free)
    kind: str | None = None
    composition_guess: Sequence[float] | Mapping[str, float] = ()  # a mapping is ordered at compile time
    packing_guess: float | None = None

    def __post_init__(self) -> None:  # None selects each default
        self.support = None if self.support is None else tuple(self.support)
        self.amount = Free() if self.amount is None else self.amount
        guess = () if self.composition_guess is None else self.composition_guess
        self.composition_guess = dict(guess) if isinstance(guess, Mapping) else tuple(guess)


@dataclass(slots=True, eq=False, repr=False, match_args=False, kw_only=True)
class Discover:
    """Request finite native phase discovery instead of a fixed phase list.

    Attributes
    ----------
    trials : sequence of Phase or None
        Trial phase declarations. ``None`` selects native defaults;
        a supplied sequence is copied and retains its support/branch choices.
    declaration_budget : int or None
        Maximum number of declaration-solve attempts, including recovery
        attempts; this is not a phase-count limit. None uses native policy.
    trial_budget : int or None
        Native limit on trial solves; ``None`` uses native policy.

    Notes
    -----
    Problem sets topology_declared false for discovery. Use solve_equilibrium
    at fixed positive T (K) and P (Pa), without solids; compile_problem requires
    declared phases. Free intensives and discovery with solids are refused.
    A finite search without a resolved descending split is not a global stability
    certificate; inspect diagnostics and exhausted budgets (ALG-STABILITY).
    """

    trials: Sequence[Phase] | None = None
    declaration_budget: int | None = None
    trial_budget: int | None = None

    def __post_init__(self) -> None:
        self.trials = None if self.trials is None else tuple(self.trials)


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Feed:
    """Declare feed values and their physical basis before native conversion.

    Attributes
    ----------
    basis : FeedKind
        Native amounts, mole fractions, mass fractions or molalities choice.
    values : mapping of str to float
        Species values: mol for amounts; dimensionless for mole or mass
        fractions; mol/kg solvent for solutes on the molality basis, with
        relative solvent mol amounts. Missing model species
        map to zero. Values are not copied by this dataclass.
    solvent_mass_kg : float or None
        Positive total solvent mass in kg for the molality basis; otherwise
        unused by the convenience constructors.
    solvent : sequence of str
        Solvent identifiers for molality conversion, copied to a tuple.

    Notes
    -----
    Prefer the four basis-specific subclasses. Native conversion checks finite
    nonnegative inputs and composition sums. Fractions are admitted within the
    native sum tolerance, without silent renormalization.
    """

    basis: _core.FeedKind
    values: Mapping[str, float]
    _: KW_ONLY
    solvent_mass_kg: float | None
    solvent: Sequence[str] = ()

    def __post_init__(self) -> None:
        self.solvent = () if self.solvent is None else tuple(self.solvent)


class MoleFractions(Feed):
    """Declare a feed by species mole fractions.

    Parameters
    ----------
    values : mapping of str to float
        Dimensionless species mole fractions; missing species are zero.

    Notes
    -----
    Native conversion uses a one-mol feed basis and requires the sum to be one
    within the composition tolerance; it does not renormalize.
    The inherited ``Feed`` attributes retain the supplied basis and values.
    Physical input validation occurs at native compilation.
    """

    def __init__(self, values: Mapping[str, float]) -> None:
        """Initialize the mole fractions feed declaration.

        Parameters
        ----------
        values : mapping of str to float
            Feed values on the basis documented by ``MoleFractions``.

        Notes
        -----
        No equilibrium calculation occurs here; native feed admission is deferred
        until compilation. See the class for units and normalization.
        """
        super().__init__(_core.FeedKind.MoleFractions, values, solvent_mass_kg=None)


class Amounts(Feed):
    """Declare an extensive feed in species mol.

    Parameters
    ----------
    values : mapping of str to float
        Nonnegative total species amounts in mol; missing species are zero.

    Notes
    -----
    The native feed total is the sum of the declared species mol. Reaction
    balances and phase amounts use this extensive basis.
    The inherited ``Feed`` attributes retain the supplied basis and values.
    Physical input validation occurs at native compilation.

    Examples
    --------
    >>> from epcsaft.equilibrium import Amounts
    >>> Amounts({'water': 2.0}).values['water']
    2.0
    """

    def __init__(self, values: Mapping[str, float]) -> None:
        """Initialize the amounts feed declaration.

        Parameters
        ----------
        values : mapping of str to float
            Feed values on the basis documented by ``Amounts``.

        Notes
        -----
        No equilibrium calculation occurs here; native feed admission is deferred
        until compilation. See the class for units and normalization.
        """
        super().__init__(_core.FeedKind.Amounts, values, solvent_mass_kg=None)


class MassFractions(Feed):
    """Declare a feed by species mass fractions.

    Parameters
    ----------
    values : mapping of str to float
        Dimensionless species mass fractions; missing species are zero.

    Notes
    -----
    Native conversion divides by component molar masses (kg/mol) on a one-kg feed basis;
    the resulting mol amounts are not normalized to one mol. Input mass fractions must sum to one
    within the composition tolerance.
    The inherited ``Feed`` attributes retain the supplied basis and values.
    Physical input validation occurs at native compilation.
    """

    def __init__(self, values: Mapping[str, float]) -> None:
        """Initialize the mass fractions feed declaration.

        Parameters
        ----------
        values : mapping of str to float
            Feed values on the basis documented by ``MassFractions``.

        Notes
        -----
        No equilibrium calculation occurs here; native feed admission is deferred
        until compilation. See the class for units and normalization.
        """
        super().__init__(_core.FeedKind.MassFractions, values, solvent_mass_kg=None)


class Molalities(Feed):
    """Declare solute molalities and a solvent mass basis.

    Parameters
    ----------
    values : mapping of str to float
        Solute molalities in mol/kg solvent; supplied solvent entries are
        relative solvent amounts in mol. Missing species are zero.
    solvent : sequence of str
        Solvent identifiers for native conversion.
    solvent_mass_kg : float
        Positive total solvent mass in kg.

    Notes
    -----
    Native conversion uses the solvent mass (kg) and component molar masses
    (kg/mol). Solute molality times total solvent kg gives solute mol; solvent
    amounts are rescaled to the declared solvent mass. A single omitted solvent
    amount is inferred from its molar mass; omitted mixed-solvent ratios are
    refused.
    The inherited ``Feed`` attributes retain the supplied basis and values.
    Physical input validation occurs at native compilation.
    """

    def __init__(
        self,
        values: Mapping[str, float],
        *,
        solvent: Sequence[str],
        solvent_mass_kg: float,
    ) -> None:
        """Initialize the molalities feed declaration.

        Parameters
        ----------
        values : mapping of str to float
            Feed values on the basis documented by ``Molalities``.
        solvent : sequence of str
            Solvent identifiers for native conversion.
        solvent_mass_kg : float
            Positive total solvent mass in kg.

        Notes
        -----
        No equilibrium calculation occurs here; native feed admission is deferred
        until compilation. See the class for units and normalization.
        """
        super().__init__(_core.FeedKind.Molalities, values, solvent=solvent, solvent_mass_kg=solvent_mass_kg)


class NeutralReference:
    """Declare a pure neutral solvent and a charge-neutral molality dilution path.

    Attributes
    ----------
    solvent : str
        Pure neutral solvent component identifier.
    solute_molality_ratios : mapping of str to float
        Dimensionless solute molality ratios defining the neutral dilution path,
        copied at construction. Native admission requires a finite positive ratio
        for every non-solvent model component and zero for the solvent. Omitted
        species map to zero, so omitting a non-solvent component is refused.
    standard_molality_mol_per_kg : float
        Positive standard molality in mol/kg solvent.
    species_reference : mapping of str to SpeciesReference
        Copied per-species standard choices. If nonempty, unspecified species
        default to infinite dilution and the solvent defaults to pure liquid.

    Notes
    -----
    Native reference admission checks charge neutrality, support, positive
    standard molality and branch realizability. This is a standard-state
    conversion declaration, not new fitted thermochemistry (ALG-REFERENCE).
    """
    __slots__ = ("solvent", "solute_molality_ratios", "standard_molality_mol_per_kg", "species_reference")

    def __init__(self, *, solvent: str, solute_molality_ratios: Mapping[str, float],
                 standard_molality_mol_per_kg: float,
                 species_reference: Mapping[str, _core.SpeciesReference] | None = None) -> None:
        """Copy a neutral-solvent reference declaration.

        Parameters
        ----------
        solvent : str
            Identifier of the pure neutral solvent.
        solute_molality_ratios : mapping of str to float
            Dimensionless neutral dilution-path ratios.
        standard_molality_mol_per_kg : float
            Standard molality in mol/kg solvent.
        species_reference : mapping of str to SpeciesReference, optional
            Explicit species standards; omitted mapping is empty.

        Notes
        -----
        Field meanings and native admission requirements are given by the class.
        Construction copies mappings without evaluating a dilution limit.
        """
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


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Reaction:
    """Declare one stoichiometric reaction with exactly one target.

    Attributes
    ----------
    stoichiometry : mapping of str to float
        Dimensionless signed coefficients: products positive, reactants
        negative. Missing species are zero, in model component order.
    equilibrium_constant : float or None
        Positive dimensionless K, not ln K or log10 K; absent when using another
        target. Without a source reference it uses the EOS density standard.
    name : str
        Reaction identifier for named reaction fitting; default is empty.
    affinity_j_per_mol : float or None
        Standard reaction affinity target in J/mol of reaction as written;
        native standard-energy constraint equals this value (``-RT ln K``).
    correlation : ReactionLogPolynomial or None
        Natural-log K(T) law; temperature in K, with its declared valid range
        and standard-state identifier.
    reference : ReactionReference or None
        Explicit source basis and reference pressure in Pa for source-basis K
        or correlations; no implicit molality/pressure conversion is inferred.

    Notes
    -----
    With no K, affinity or correlation, native compilation derives targets from
    participating species thermochemistry. More than one supplied target raises
    ``ParameterError``. Native compilation checks stoichiometry,
    reference consistency and source applicability. The EOS activity is
    ``(rho*x_i / (1 mol/m³))*exp(mu_res_i/(RT))``; K is the product of these
    activities raised to signed stoichiometric coefficients (ALG-EQUILIBRIUM).
    A fixed K alone does not supply arbitrary thermal derivatives.

    Examples
    --------
    >>> from epcsaft.equilibrium import Reaction
    >>> reaction = Reaction({'A': -1.0, 'B': 1.0}, equilibrium_constant=2.0)
    >>> reaction.equilibrium_constant
    2.0
    """
    stoichiometry: Mapping[str, float]
    _: KW_ONLY
    equilibrium_constant: float | None = None
    name: str = ""
    affinity_j_per_mol: float | None = None
    correlation: _core.ReactionLogPolynomial | None = None
    reference: _core.ReactionReference | None = None

    def __post_init__(self) -> None:
        if sum(item is not None for item in (self.equilibrium_constant, self.affinity_j_per_mol, self.correlation)) > 1:
            raise _core.ParameterError("reaction admits at most one equilibrium_constant, affinity_j_per_mol or correlation")
        self.name = "" if self.name is None else self.name


@dataclass(slots=True, eq=False, repr=False, match_args=False)
class Solid:
    """Declare a pure stoichiometric solid dissolving into fluid species.

    Attributes
    ----------
    name : str
        Solid identifier.
    formula : mapping of str to float
        Dimensionless dissolution coefficients for ``S -> sum_i nu_i A_i``;
        dissolved formula must be charge-neutral.
    reference : ReactionReference
        Source constant basis with reference pressure in Pa.
    molar_volume_m3_per_mol : float
        Positive pure-solid molar volume in m³/mol of formula units.
    equilibrium_constant : float or None
        Positive dimensionless source dissolution K, not its logarithm;
        supply exactly one of this and ``correlation``.
    correlation : ReactionLogPolynomial or None
        Source natural-log K(T) law with temperature in K and valid interval.
    amount : Free or Pinned
        Free solid amount solved in mol of formula units, or ``Pinned(0.0)``
        for the absent-solid condition; other pins are refused natively.

    Notes
    -----
    Native compilation owns source-reference conversion, molar-volume pressure
    correction and solid complementarity. No solid-solution or kinetic model is
    declared here (ALG-EQUILIBRIUM, pure stoichiometric solids).
    """
    name: str
    _: KW_ONLY
    formula: Mapping[str, float]
    reference: _core.ReactionReference
    molar_volume_m3_per_mol: float
    equilibrium_constant: float | None = None
    correlation: _core.ReactionLogPolynomial | None = None
    amount: Free | Pinned = field(default_factory=Free)

    def _native(self, model: _Species) -> _core.SolidSpec:
        item = _core.SolidSpec()
        item.name = self.name
        item.stoichiometry = _vector(self.formula, model.component_ids)
        item.equilibrium_constant = self.equilibrium_constant
        item.correlation = self.correlation
        item.reference = self.reference
        item.molar_volume_m3_per_mol = self.molar_volume_m3_per_mol
        item.amount = _amount_spec(self.amount, self.name)
        return item


def _index(component_ids: tuple[str, ...], name: str) -> int:
    try:
        return component_ids.index(name)
    except ValueError as error:
        raise _core.ParameterError(f"unknown component_id {name}") from error


def _vector(values: Mapping[str, float], component_ids: tuple[str, ...]) -> list[float]:
    if unknown := [name for name in values if name not in component_ids]:
        raise _core.ParameterError(f"unknown species {unknown[0]}; expected one of {', '.join(component_ids)}")
    return [values.get(component_id, 0.0) for component_id in component_ids]


@dataclass(slots=True, eq=False, repr=False, match_args=False, kw_only=True)
class Problem:
    """Declare the physical specification for the one native equilibrium driver.

    Attributes
    ----------
    phases : sequence of Phase or Discover
        Ordered fluid phases, or a finite phase-discovery request.
    T : float or Free
        Fixed temperature in K or a free temperature declaration with K guess.
    P : float or Free
        Fixed pressure in Pa or a free pressure declaration with Pa guess.
    feed : Feed
        Species inventory with explicit amounts/fractions/molality basis.
    reactions : sequence of Reaction
        Stoichiometric reaction declarations, empty for nonreacting fluids.
    standard_gibbs_j_per_mol : sequence of float
        Species standard Gibbs energies in J/mol, in component order; empty
        leaves native standard-reference completion to supplied inputs.
    topology_declared : bool or None
        Whether the phase topology is declared; ``None`` defaults true.
        ``Discover`` always sets this false.
    row_tolerance : float or None
        Dimensionless tolerance for native scaled equation rows; ``None`` keeps
        the native policy. It is not a physical uncertainty estimate.
    solver_options : mapping of str to float or None
        Native tolerances and work limits with the units of each named option.
    enthalpy_j : float or None
        Total extensive enthalpy constraint in J, on the compiled feed basis.
    entropy_j_per_k : float or None
        Total extensive entropy constraint in J/K, on the compiled feed basis.
    neutral_reference : NeutralReference or None
        Source standard-state conversion declaration.
    critical : CriticalDeclaration or None
        Critical-equation seed, including density in mol/m³; a solved necessary
        condition is not a stability certificate or unique critical point.
    envelope : EnvelopeCondition or None
        Local envelope equation using natural-log temperature, pressure or
        incipient/feed mole-fraction ratio with native dimensionless references.
    azeotrope : bool
        Request native azeotrope conditions; default false.
    solids : sequence of Solid
        Ordered pure stoichiometric solid declarations; default empty.

    Notes
    -----
    The constructor copies phase/reaction/solid sequences and standard energies.
    Native compilation checks support, charge, conserved inventory, number and
    rank of equations, reference consistency and admitted caloric constraints.
    All root finds, equilibrium iterations and continuation are native.
    Returned diagnostics, not an optimizer stop alone, establish numerical
    closure (ALG-EQUILIBRIUM and ALG-TERMINATION).

    Examples
    --------
    >>> from epcsaft.equilibrium import Amounts, Phase, Problem
    >>> problem = Problem(phases=[Phase('L', kind='liquid')], T=300.0,
    ...                   P=100000.0, feed=Amounts({'water': 1.0}))
    >>> problem.T, problem.P
    (300.0, 100000.0)
    """

    phases: Sequence[Phase] | Discover = ()
    T: float | Free
    P: float | Free
    feed: Feed
    reactions: Sequence[Reaction] = ()
    standard_gibbs_j_per_mol: Sequence[float] = ()
    topology_declared: bool | None = True
    row_tolerance: float | None = None
    solver_options: Mapping[str, float] | None = None
    enthalpy_j: float | None = None
    entropy_j_per_k: float | None = None
    neutral_reference: NeutralReference | None = None
    critical: _core.CriticalDeclaration | None = None
    envelope: _core.EnvelopeCondition | None = None
    azeotrope: bool = False
    solids: Sequence[Solid] = ()

    def __post_init__(self) -> None:  # None selects each default
        discover = isinstance(self.phases, Discover)
        if not isinstance(self.phases, Discover):
            self.phases = tuple(() if self.phases is None else self.phases)
        self.reactions = tuple(() if self.reactions is None else self.reactions)
        self.solids = tuple(() if self.solids is None else self.solids)
        self.standard_gibbs_j_per_mol = tuple(() if self.standard_gibbs_j_per_mol is None else self.standard_gibbs_j_per_mol)
        self.topology_declared = (True if self.topology_declared is None else self.topology_declared) and not discover

    def _native(self, model: _Species) -> _core.EquilibriumProblem:
        native = _core.EquilibriumProblem()
        native.topology_declared = bool(self.topology_declared)
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
        native.azeotrope = self.azeotrope
        native.temperature = _free_value(self.T)
        native.pressure = _free_value(self.P)
        native.neutral_reference = None if self.neutral_reference is None else self.neutral_reference._native(model)
        native.solver_options = self.solver_options
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
            item.amount = _amount_spec(phase.amount, phase.name)
            guess = phase.composition_guess
            order = tuple(model.component_ids[i] for i in item.support) or model.component_ids
            item.composition_guess = _vector(guess, order) if isinstance(guess, Mapping) else list(guess)
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
        native.solids = [solid._native(model) for solid in self.solids]
        return native


def problem_from_mapping(mapping: Mapping[str, object]) -> Problem:
    """Convert named problem fields into a Python equilibrium declaration.

    Parameters
    ----------
    mapping : mapping of str to object
        ``T`` in K, ``P`` in Pa, phase records, feed basis/values and optional
        reactions, standard Gibbs energies (J/mol), topology, row tolerance,
        solver options, neutral reference, total enthalpy (J) or entropy (J/K).
        Intensives can be scalars, ``{'value': ...}`` or free/guess mappings.
        Phase pins use dimensionless feed fractions.

    Returns
    -------
    problem : Problem
        Declaration only; no equilibrium solve has occurred.

    Raises
    ------
    ParameterError
        Unknown top-level key, missing intensive, invalid phase amount, feed
        basis, record structure, reaction target or reference basis.
    ValueError
        A requested numeric mapping value cannot be converted to float.

    Notes
    -----
    ``phases='discover'`` or a ``Discover`` object requests discovery. This
    mapping reader does not accept the direct ``Problem`` fields ``critical``,
    ``envelope``, ``azeotrope`` or ``solids``. Unknown species and physical
    constraints are checked when the declaration is compiled with a model.
    """
    allowed = {
        "phases", "T", "P", "feed", "reactions", "standard_gibbs_j_per_mol", "topology_declared",
        "row_tolerance", "solver_options", "neutral_reference", "enthalpy_j", "entropy_j_per_k",
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
        solver_options=cast(Mapping[str, float] | None, mapping.get("solver_options")),
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


def _amount_spec(value: Free | Pinned, name: str) -> _core.AmountSpec:
    if isinstance(value, Free) and value.guess is not None:
        raise _core.ParameterError(f"phase {name}: phase-amount guesses are not supported; use Free()")
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
    model: _Model | IdealMoleFraction,
    problem: Problem | _core.EquilibriumProblem,
    start: _core.EquilibriumResult | None = None,
) -> _core.EquilibriumResult:
    """Solve a declared phase/reaction problem in the native driver.

    Parameters
    ----------
    model : Mixture or Pure
        Fluid model whose component order names every declared species.
    problem : Problem or _core.EquilibriumProblem
        Physical specification, with T in K, P in Pa and declared feed basis.
    start : EquilibriumResult, optional
        Prior native result supplying a start; branch/layout admission remains
        native and is not bypassed by a supplied result.

    Returns
    -------
    result : EquilibriumResult
        Temperature in K, pressure in Pa, phase amounts in mol and
        dimensionless mole fractions in declared phase/component order, with
        achieved residuals, stop causes and availability diagnostics.

    Raises
    ------
    ParameterError
        Python species/phase mapping or native physical input admission fails.
    DeclarationError
        Native supports, structure, constraints or reference standards conflict.

    Notes
    -----
    Numerical nonconvergence retains native diagnostics. A fixed phase list is
    not a global stability certificate; discovery is a bounded native search
    (ALG-EQUILIBRIUM and ALG-STABILITY).
    """
    native = problem._native(model) if isinstance(problem, Problem) else problem
    return _core.solve_equilibrium(model._model, native, start)


def compile_problem(model: _Model | IdealMoleFraction, problem: Problem | _core.EquilibriumProblem) -> _core.CompiledSystem:
    """Compile a physical equilibrium declaration for reuse and exact actions.

    Parameters
    ----------
    model : Mixture or Pure
        Fluid model whose component order names every declared species.
    problem : Problem or _core.EquilibriumProblem
        Physical specification, with T in K, P in Pa and declared feed basis.

    Returns
    -------
    system : CompiledSystem
        Native balances, supports, coordinates and equation rows; keeps the
        native model alive for solved-state derivative actions.

    Raises
    ------
    ParameterError
        Python species/phase mapping or native physical input admission fails.
    DeclarationError
        Native supports, structure, constraints or reference standards conflict.

    Notes
    -----
    Numerical nonconvergence retains native diagnostics. A fixed phase list is
    not a global stability certificate; discovery is a bounded native search
    (ALG-EQUILIBRIUM and ALG-STABILITY).
    """
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
    """Continue between physical equilibrium declarations using native correctors.

    Parameters
    ----------
    model : Mixture or Pure
        Fluid model with fixed component/site order.
    start, target : Problem
        Endpoint declarations. Temperature uses K, pressure Pa, feed amounts
        mol, and total caloric constraints J or J/K where supplied.
    options : Continuation
        Declared path step controls, tolerances and work budgets.
    start_result : EquilibriumResult, optional
        Previously solved starting endpoint.
    parameter_targets : sequence of float, optional
        Target values in the model's active-parameter order and native family
        units; ``None`` supplies an empty parameter-target sequence.

    Returns
    -------
    result : ContinuationResult
        Accepted path steps and last solved state, with reached dimensionless
        path coordinate, stop cause and diagnostics. Incomplete paths remain
        incomplete even when the last accepted point is available.

    Raises
    ------
    ParameterError
        Species/phase mapping or physical path inputs are inadmissible.
    DeclarationError
        Endpoint supports, equations or reference declarations are incompatible.

    Notes
    -----
    Prediction, corrector solves, step rejection and work limits are native.
    Inspect ``success``, ``stop`` and ``reached_path`` before treating a returned
    state as the requested endpoint (ALG-EQUILIBRIUM).
    """
    request = _core.ContinuationRequest()
    request.start_problem = start._native(model)
    request.target_problem = target._native(model)
    request.start_result = start_result
    request.parameter_targets = () if parameter_targets is None else parameter_targets
    request.options = options
    return _core.continue_equilibrium(model._model, request)


def bubble_point(*, liquid: Phase, T: float, feed: Feed, reactions: Sequence[Reaction] | None = None,
                 vapor_support: Sequence[str] | None = None) -> Problem:
    """Declare an isothermal bubble point with free pressure.

    Parameters
    ----------
    liquid : Phase
        Existing liquid phase declaration, retaining its support and amount.
    T : float
        Fixed absolute temperature in K.
    feed : Feed
        Feed declaration with explicit amounts/fractions/molality basis.
    reactions : sequence of Reaction, optional
        Reaction declarations; ``None`` means no reactions.
    vapor_support : sequence of str, optional
        Explicit vapor species; ``None`` keeps neutral-only completion.

    Returns
    -------
    problem : Problem
        Original liquid phase plus incipient vapor pinned at zero mol,
        with pressure free in Pa. Incipient phase name is V.

    Notes
    -----
    No equilibrium calculation is performed. Native compilation/solve checks
    admitted phase kinds, supports, reactions and boundary existence.
    """
    reactions = () if reactions is None else reactions
    vapor = Phase("V", kind="vapor", support=vapor_support, amount=Pinned(0.0))
    return Problem(
        phases=[liquid, vapor],
        T=T,
        P=Free(),
        feed=feed,
        reactions=reactions,
    )


def dew_point(*, vapor: Phase, T: float, feed: Feed, reactions: Sequence[Reaction] | None = None) -> Problem:
    """Declare an isothermal dew point with free pressure.

    Parameters
    ----------
    vapor : Phase
        Existing vapor phase declaration, retaining its support and amount.
    T : float
        Fixed absolute temperature in K.
    feed : Feed
        Feed declaration with explicit amounts/fractions/molality basis.
    reactions : sequence of Reaction, optional
        Reaction declarations; ``None`` means no reactions.

    Returns
    -------
    problem : Problem
        Original vapor phase plus incipient liquid pinned at zero mol,
        with pressure free in Pa. Incipient phase name is L.

    Notes
    -----
    No equilibrium calculation is performed. Native compilation/solve checks
    admitted phase kinds, supports, reactions and boundary existence.
    """
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
    """Run declared starts under shared native equilibrium work limits.

    Parameters
    ----------
    model : Mixture or Pure
        Fluid model whose component order defines the inputs.
    problem : Problem
        Physical declaration: temperature K, pressure Pa, feed on its stated basis.
    starts : sequence of EquilibriumStart
        Ordered native starts, with temperature K, pressure Pa, physical amounts
        mol and dimensionless compositions/packing fractions where supplied.
    limits : EquilibriumWorkLimits
        Attempt/iteration limits and wall-time budgets in s.

    Returns
    -------
    result : EquilibriumStartsResult
        Per-attempt results, causes and consumed work, with an optional selected
        result. Selection does not establish global stability or uniqueness.

    Raises
    ------
    ParameterError
        Python mapping or native start/work-limit admission fails.
    DeclarationError
        The physical declaration is refused before start attempts.

    Notes
    -----
    The native owner enforces shared work limits and retains exhausted-budget
    causes rather than silently retrying without a bound (ALG-EQUILIBRIUM).
    """
    return _core.solve_equilibrium_starts(model._model, problem._native(model), list(starts), limits)
