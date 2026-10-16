from typing import Any, Mapping, Optional, Sequence, Union, overload

GAS_CONSTANT_J_PER_MOL_K: float
EOS_STANDARD_STATE_ID: str

class ParameterError(ValueError): ...
class EosError(RuntimeError): ...
class DensityStructureError(EosError): ...
class ConvergenceError(EosError): ...
class DeclarationError(ParameterError): ...

class Stop:
    converged: Stop
    converged_to_floor: Stop
    stalled: Stop
    budget_exhausted: Stop
    left_domain: Stop
    singular: Stop
    infeasible: Stop
    nonfinite: Stop

class Category:
    ThermodynamicConstraint: Category
    AlgorithmTolerance: Category

class Tolerance:
    name: str
    value: float
    category: Category
    source: str
    factor: float
    upstream: Optional[Tolerance]

def tolerances() -> list[Tolerance]: ...

class Convergence:
    residual: float
    criterion: Optional[Tolerance]
    floor: Optional[float]
    amplification: Optional[float]
    contraction: Optional[float]
    iterations: int
    evaluations: int
    stop: Stop
    @property
    def residual_over_criterion(self) -> Optional[float]: ...
    @property
    def residual_over_floor(self) -> Optional[float]: ...
    @property
    def amplification_over_declared(self) -> Optional[float]: ...

class Event:
    name: str
    span: int
    iteration: int
    values: list[float]
    reason: str

class Span:
    name: str
    parent: Optional[int]
    convergence: Optional[Convergence]
    wall_us: float

class Trace:
    spans: list[Span]
    events: list[Event]
    current: Optional[int]
    incomplete: bool

class Capture:
    def __init__(self, diagnostics: bool = False) -> None: ...
    trace: Optional[Trace]
    spans: list[Span]
    events: list[Event]
    def __enter__(self) -> Capture: ...
    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> bool: ...

class ParameterEntry:
    identity: str
    field: str
    index: list[int]
    magnitude: float
    unit: str
    def __init__(self) -> None: ...

class ParameterDocument:
    component_ids: list[str]
    site_component: list[int]
    site_multiplicity: list[int]
    base_choice: str
    association_choice: str
    electrolyte_choice: str
    ionic_dispersion_choice: str
    permittivity_choice: str
    polar_choice: str
    entries: list[ParameterEntry]
    def __init__(self) -> None: ...

class Parameters:
    @property
    def component_count(self) -> int: ...
    @property
    def site_count(self) -> int: ...
    @property
    def molar_masses(self) -> list[float]: ...  # kg/mol, component_ids order
    @property
    def temperature_law_coefficient_identities(self) -> list[str]: ...

def build_parameters(document: ParameterDocument) -> Parameters: ...

class IdealPolynomial:
    coefficients: list[float]
    reference_temperature: float
    reference_enthalpy: float
    reference_entropy: float
    def __init__(self, coefficients: Sequence[float], reference_temperature: float,
                 reference_enthalpy: float, reference_entropy: float) -> None: ...

class IdealShomate:
    coefficients: list[float]
    formation_enthalpy: float
    def __init__(self, coefficients: Sequence[float], formation_enthalpy: float) -> None: ...

class IdealNASA7:
    coefficients: list[float]
    def __init__(self, coefficients: Sequence[float]) -> None: ...

class IdealNASA9:
    coefficients: list[float]
    def __init__(self, coefficients: Sequence[float]) -> None: ...

class IdealDIPPR107:
    coefficients: list[float]
    reference_temperature: float
    reference_enthalpy: float
    reference_entropy: float
    def __init__(self, coefficients: Sequence[float], reference_temperature: float,
                 reference_enthalpy: float, reference_entropy: float) -> None: ...

IdealForm = Union[IdealPolynomial, IdealShomate, IdealNASA7, IdealNASA9, IdealDIPPR107]

class IdealCorrelation:
    form: IdealForm
    source_pressure: float
    def __init__(self, form: IdealForm, source_pressure: float) -> None: ...

class IdealInterval:
    lower: float
    upper: float
    lower_inclusive: bool
    upper_inclusive: bool
    correlation: IdealCorrelation
    def __init__(self, lower: float, upper: float, lower_inclusive: bool,
                 upper_inclusive: bool, correlation: IdealCorrelation) -> None: ...

class ThermochemistryRecord:
    standard_pressure: float
    component_intervals: list[list[IdealInterval]]
    def __init__(self, standard_pressure: float,
                 component_intervals: Sequence[Sequence[IdealInterval]]) -> None: ...

class ReactionLogPolynomial:
    """ln K(T) = a + b/T + c ln(T/reference_temperature) + d T, T in K, valid on [temperature_min, temperature_max].
    Without Reaction(reference=...), solve_equilibrium requires standard_state_id == EOS_STANDARD_STATE_ID; with
    a ReactionReference and the problem's neutral_reference, a source-basis correlation is converted natively."""
    a: float
    b: float
    c: float
    d: float
    reference_temperature: float
    temperature_min: float
    temperature_max: float
    standard_state_id: str
    def __init__(self, a: float, b: float, c: float, d: float, reference_temperature: float,
                 temperature_min: float, temperature_max: float, standard_state_id: str) -> None: ...

class ParameterFamily:
    BornCDielectric: ParameterFamily
    IonicRegionRelativePermittivity: ParameterFamily
    BornDiameter: ParameterFamily
    DebyeHuckelDiameter: ParameterFamily
    PackingDiameter: ParameterFamily
    IonSuppressionCoefficient: ParameterFamily
    SegmentCount: ParameterFamily
    SigmaConstant: ParameterFamily
    DispersionEnergyOverK: ParameterFamily
    Kij: ParameterFamily
    Lij: ParameterFamily
    AssociationEnergyOverK: ParameterFamily
    AssociationVolume: ParameterFamily
    KijReciprocalTemperatureSlope: ParameterFamily
    DielectricVolume: ParameterFamily
    DielectricTemperature: ParameterFamily
    TemperatureLawCoefficient: ParameterFamily

class ActiveParameter:
    family: ParameterFamily
    indices: list[int]
    def __init__(self, family: ParameterFamily, indices: Sequence[int]) -> None: ...

class DerivativeOrder:
    Value: DerivativeOrder
    Gradient: DerivativeOrder
    PressureSlope: DerivativeOrder
    Hessian: DerivativeOrder

class NeutralProperty:
    ResidualHelmholtz: NeutralProperty
    Pressure: NeutralProperty

class NeutralPropertyDerivatives:
    value: float
    gradient: list[float]
    hessian: list[float]

class PropertyObservable:
    IdealEnthalpy: PropertyObservable
    IdealInternalEnergy: PropertyObservable
    IdealEntropy: PropertyObservable
    IdealIsobaricHeatCapacity: PropertyObservable
    IdealIsochoricHeatCapacity: PropertyObservable
    ResidualEnthalpy: PropertyObservable
    ResidualInternalEnergy: PropertyObservable
    ResidualEntropy: PropertyObservable
    ResidualIsobaricHeatCapacity: PropertyObservable
    ResidualIsochoricHeatCapacity: PropertyObservable
    TotalEnthalpy: PropertyObservable
    TotalInternalEnergy: PropertyObservable
    TotalEntropy: PropertyObservable
    TotalIsobaricHeatCapacity: PropertyObservable
    TotalIsochoricHeatCapacity: PropertyObservable
    LogFugacityCoefficient: PropertyObservable

class PropertyActionStatus:
    Available: PropertyActionStatus
    MissingThermochemistry: PropertyActionStatus
    OutsideInterval: PropertyActionStatus
    ThermalBoundary: PropertyActionStatus
    SupportChange: PropertyActionStatus
    IonicInsertion: PropertyActionStatus
    AssociationSingular: PropertyActionStatus
    AssociationIllConditioned: PropertyActionStatus
    AssociationResidual: PropertyActionStatus
    DensityUnavailable: PropertyActionStatus
    PressureSingular: PropertyActionStatus
    CentralPropertyUnavailable: PropertyActionStatus
    Nonfinite: PropertyActionStatus

class ParameterActionDirection:
    family: ParameterFamily
    indices: list[int]
    value: float
    def __init__(self, family: ParameterFamily, indices: Sequence[int], value: float) -> None: ...

class DensityActionDirection:
    temperature: float
    molar_density: float
    composition: list[float]
    parameter: list[ParameterActionDirection]
    def __init__(self, temperature: float = 0.0, molar_density: float = 0.0,
                 composition: Sequence[float] = (),
                 parameter: Sequence[ParameterActionDirection] = ()) -> None: ...

class PressureActionDirection:
    temperature: float
    pressure: float
    composition: list[float]
    parameter: list[ParameterActionDirection]
    def __init__(self, temperature: float = 0.0, pressure: float = 0.0,
                 composition: Sequence[float] = (),
                 parameter: Sequence[ParameterActionDirection] = ()) -> None: ...

class PropertyActionResult:
    observable: PropertyObservable
    component: Optional[int]
    value: Optional[float]
    action: Optional[float]
    status: PropertyActionStatus
    energy_order: int

class DensityBranch:
    Unspecified: DensityBranch
    Vapor: DensityBranch
    Liquid: DensityBranch
    def __int__(self) -> int: ...

class ThermochemistryStatus:
    MissingRecord: ThermochemistryStatus
    Available: ThermochemistryStatus
    OutsideInterval: ThermochemistryStatus
    Nonfinite: ThermochemistryStatus

class PropertyActionDiagnostics:
    site_rank: int
    site_factorizations: int
    site_solves: int
    site_pivot_ratio: float
    association_residual: float
    pressure_density_derivative: float
    pressure_residual: float
    stability_determinant: float
    pressure_factorizations: int
    pressure_solves: int
    thermal_evaluations: int
    recordings: int
    directional_paths: int
    forward_sweeps: int
    branch: DensityBranch
    stable_root_count: Optional[int]

class PropertyActionBatch:
    solver_options: list[dict[str, Any]]
    results: list[PropertyActionResult]
    maximum_energy_order: int
    diagnostics: PropertyActionDiagnostics
    molar_density: Optional[float]
    density_action: list[float]
    density_mixed_action: Optional[float]
    density_status: PropertyActionStatus
    density_message: str
    density_convergence: Optional[Convergence]

class State:
    solver_options: list[dict[str, Any]]
    """One homogeneous state in SI molar units at the given T (K) and x (component_ids order).
    The contribution names through residual_helmholtz are A_res/(nRT), dimensionless."""
    hard_chain: float
    dispersion: float
    association: float
    debye_huckel: float
    born: float
    dipolar: float
    bulk_relative_permittivity: float
    residual_helmholtz: float
    compressibility_factor: float
    pressure: float  # Pa
    # Per component, component_ids order. mu_i^res(T, rho)/(RT), residual against the ideal gas at the same T and rho;
    # None for a charged component on an ion-free support; see component_property_status (IonicInsertion).
    residual_chemical_potential_over_rt: list[Optional[float]]
    # ln(fugacity coefficient) at the state's T and P: mu_i^res/(RT) - ln Z. The list is None unless Z is positive and
    # finite; an entry is None where residual_chemical_potential_over_rt is.
    log_fugacity_coefficient: Optional[list[Optional[float]]]
    component_property_status: list[PropertyActionStatus]
    log_fugacity_parameter_derivative: list[Optional[float]]
    # Residual (minus ideal gas at the same T and rho): J/mol for energies, J/(mol K) for entropy and heat capacities.
    residual_entropy: float
    residual_internal_energy: float
    residual_enthalpy: float
    residual_isochoric_heat_capacity: float
    residual_isobaric_heat_capacity: Optional[float]
    ideal_enthalpy: Optional[float]
    ideal_internal_energy: Optional[float]
    ideal_entropy: Optional[float]
    ideal_isobaric_heat_capacity: Optional[float]
    ideal_isochoric_heat_capacity: Optional[float]
    total_enthalpy: Optional[float]
    total_internal_energy: Optional[float]
    total_entropy: Optional[float]
    total_isobaric_heat_capacity: Optional[float]
    total_isochoric_heat_capacity: Optional[float]
    molar_density: float  # mol/m^3; times the mixture molar mass (kg/mol) for kg/m^3
    pressure_density_derivative: float  # (dP/drho)_T,x in Pa m^3/mol
    site_fractions: list[float]  # unbonded site fractions, Parameters.site_ids order
    density_convergence: Optional[Convergence]
    thermochemistry_status: ThermochemistryStatus
    density_branch: DensityBranch
    stable_root_count: Optional[int]
    # Requested by derived_properties=True, else None; fixed composition. None where the cause applies:
    # nonpositive (dP/drho)_T for all, missing total Cp/Cv (thermochemistry) for Cp-dependent values,
    # an unoccupied species without positive molar mass for sound.
    isothermal_compressibility: Optional[float]  # 1/Pa
    isobaric_expansivity: Optional[float]  # 1/K
    heat_capacity_ratio: Optional[float]  # total Cp/Cv
    isentropic_compressibility: Optional[float]  # 1/Pa
    speed_of_sound: Optional[float]  # m/s, thermodynamic composition-frozen limit
    joule_thomson_coefficient: Optional[float]  # K/Pa
    # Per component (component_ids order); entries None for IonicInsertion, enthalpy/entropy None without
    # Available thermochemistry, entropy None at zero mole fraction.
    partial_molar_volume: Optional[list[Optional[float]]]  # m^3/mol
    partial_molar_enthalpy: Optional[list[Optional[float]]]  # J/mol
    partial_molar_entropy: Optional[list[Optional[float]]]  # J/(mol K)

class DensityDiagnosticPoint:
    molar_density: float
    packing_fraction: float
    pressure: float
    pressure_density_derivative: float
    lower: float
    upper: float
    residual: float
    mechanical_sign: int
    linked: bool
    declared: bool
    cause: str
    stop: str

class DensityDiagnosticBracket:
    lower: float
    upper: float
    lower_residual: float
    upper_residual: float

class DensityDiagnostics:
    solver_options: list[dict[str, Any]]
    samples: list[DensityDiagnosticPoint]
    pressure_roots: list[DensityDiagnosticPoint]
    spinodals: list[DensityDiagnosticPoint]
    declared_brackets: list[DensityDiagnosticBracket]
    unresolved: list[DensityDiagnosticPoint]
    evaluation_order: DerivativeOrder
    finite_domain: Optional[DensityDiagnosticPoint]
    geometric_cap: float
    declared_crossings: int
    evaluations: int
    stop: str

class SweepCounters:
    forward: int
    reverse: int
    hessian_calls: int
    recordings: int

class HelmholtzModel:
    solver_options: list[dict[str, Any]]
    def __init__(self, parameters: Parameters, active_parameters: Sequence[ActiveParameter] = (),
                 thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> None: ...
    def state(self, temperature: float, molar_density: float, mole_fractions: Sequence[float],
              derived_properties: bool = False) -> State: ...
    def neutral_property(self, temperature: float, molar_density: float, property: NeutralProperty, order: DerivativeOrder) -> NeutralPropertyDerivatives: ...
    def density_actions(self, temperature: float, molar_density: float,
                        mole_fractions: Sequence[float], observables: Sequence[PropertyObservable],
                        directions: Sequence[DensityActionDirection]) -> PropertyActionBatch: ...
    def pressure_actions(self, temperature: float, pressure: float,
                         mole_fractions: Sequence[float], branch: DensityBranch,
                         observables: Sequence[PropertyObservable],
                         directions: Sequence[PressureActionDirection],
                         anchor_density: Optional[float] = None) -> PropertyActionBatch: ...
    def packing_fraction_coefficient(self, temperature: float, composition: list[float]) -> float: ...
    def state_at_pressure(self, temperature: float, pressure: float, mole_fractions: Sequence[float], branch: DensityBranch,
                          anchor_density: Optional[float] = None, derived_properties: bool = False) -> State: ...
    def molar_density_upper_bound(self, temperature: float, composition: Sequence[float]) -> float: ...
    def density_diagnostics(self, temperature: float, pressure: float, composition: Sequence[float],
                            nodes: Sequence[float], density_tolerance: float, log_tolerance: float,
                            max_evaluations: int, max_seconds: float, adaptive: bool,
                            finite_anchor: Optional[float], finite_upper: float,
                            evaluation_order: DerivativeOrder) -> DensityDiagnostics: ...
    def counters(self) -> SweepCounters: ...
    def reset_counters(self) -> None: ...

class PhaseKind:
    Liquid: PhaseKind
    Vapor: PhaseKind

class AmountKind:
    Free: AmountKind
    Pinned: AmountKind

class FeedKind:
    MoleFractions: FeedKind
    Amounts: FeedKind
    MassFractions: FeedKind
    Molalities: FeedKind

class AmountSpec:
    kind: AmountKind
    fraction: float
    def __init__(self) -> None: ...
    @staticmethod
    def free() -> AmountSpec: ...
    @staticmethod
    def pinned(value: float) -> AmountSpec: ...

class FreeValue:
    free: bool
    value: float
    def __init__(self) -> None: ...
    @staticmethod
    def fixed(value: float) -> FreeValue: ...
    @staticmethod
    def variable(guess: Optional[float] = None) -> FreeValue: ...

class PhaseSpec:
    name: str
    kind: PhaseKind
    support: list[int]
    support_declared: bool
    amount: AmountSpec
    composition_guess: list[float]
    packing_guess: Optional[float]
    def __init__(self) -> None: ...

class FeedSpec:
    kind: FeedKind
    values: list[float]
    solvent: list[int]
    solvent_mass_kg: float
    def __init__(self) -> None: ...

class SpeciesReference:
    InfiniteDilution: SpeciesReference
    PureLiquid: SpeciesReference
    PureGas: SpeciesReference

class NeutralReferenceDeclaration:
    solvent_index: Optional[int]
    solute_molality_ratios: list[float]
    standard_molality_mol_per_kg: float
    species_reference: list[SpeciesReference]
    def __init__(self) -> None: ...

class ReferenceSourceBasis:
    CommonMolalityInfiniteDilution: ReferenceSourceBasis
    RawMoleFractionInfiniteDilution: ReferenceSourceBasis

class NeutralReferenceStatus:
    Available: NeutralReferenceStatus
    InvalidDeclaration: NeutralReferenceStatus
    NonNeutralDeclaration: NeutralReferenceStatus
    UnsupportedSourceBasis: NeutralReferenceStatus
    FreeTemperature: NeutralReferenceStatus
    FreePressure: NeutralReferenceStatus
    DensityUnavailable: NeutralReferenceStatus
    BranchChange: NeutralReferenceStatus
    Nonfinite: NeutralReferenceStatus
    TerminalChangeExceeded: NeutralReferenceStatus

class NeutralReferenceReactionDiagnostic:
    solver_options: list[dict[str, Any]]
    reaction_index: int
    terminal_contraction: Optional[float]
    observed_terminal_change: Optional[float]
    terminal_contraction_pressure_derivative: Optional[float]
    terminal_contraction_pressure_second_derivative: Optional[float]
    reference_pressure_pa: float
    status: NeutralReferenceStatus
    message: str

class ReactionReference:
    source_basis: ReferenceSourceBasis
    reference_pressure_pa: Optional[float]
    def __init__(self, source_basis: ReferenceSourceBasis,
                 reference_pressure_pa: Optional[float] = None) -> None: ...

class ReactionSpec:
    name: str
    stoichiometry: list[float]
    equilibrium_constant: Optional[float]
    affinity_j_per_mol: Optional[float]
    correlation: Optional[ReactionLogPolynomial]
    reference: Optional[ReactionReference]
    def __init__(self) -> None: ...

class StabilityOptions:
    evaluate_supplied_trials: bool
    trials: Optional[list[PhaseSpec]]
    declaration_budget: Optional[int]
    trial_budget: Optional[int]
    standard_uncertainty_over_rt: list[float]
    def __init__(self) -> None: ...

class StabilityTrial:
    phase: PhaseSpec
    seed: str
    evidence: str
    message: str
    optimizer: EquilibriumResult
    raw_coordinates: list[float]
    raw_residuals: list[float]
    composition: list[float]
    split_coordinates: list[float]
    basis: list[list[float]]
    row_transformation: list[list[float]]
    quantities: dict[str, float]
    split_checks: list[dict[str, float]]
    donor: Optional[int]
    resolved_descent: bool

class StabilityReport:
    completed: bool
    message: str
    trials: list[StabilityTrial]
    candidate_coordinates: list[float]
    standard_over_rt: list[float]
    plane_over_rt: list[float]
    present_phases: list[dict[str, float]]
    work: dict[str, float]

class DiscoveryRecord:
    declaration: EquilibriumProblem
    candidate: EquilibriumResult
    rounds: list[StabilityReport]
    declarations: list[EquilibriumResult]
    declaration_roles: list[str]
    declaration_attempt_coordinates: list[list[float]]
    work: dict[str, float]


def cooling_onset(compiled: CompiledSystem) -> tuple[EquilibriumResult, list[StabilityReport]]: ...

def stability_test(compiled: CompiledSystem, central: EquilibriumResult,
                   options: StabilityOptions = ...) -> StabilityReport: ...

class CriticalDeclaration:
    density_guess_mol_per_m3: float
    mode_guess: list[float]
    def __init__(self) -> None: ...

class EnvelopeCoordinate:
    LogTemperature: EnvelopeCoordinate
    LogPressure: EnvelopeCoordinate
    LogK: EnvelopeCoordinate

class EnvelopeCondition:
    kind: EnvelopeCoordinate
    component: int
    target: float
    def __init__(self) -> None: ...

class CriticalRecord:
    eigenvalues: list[float]
    mode: list[float]
    cubic: float
    zero_band: float
    condition_estimate: float
    backward_error: float
    rank: int

class EquilibriumProblem:
    @property
    def solver_options(self) -> list[dict[str, Any]]: ...
    @solver_options.setter
    def solver_options(self, value: Mapping[str, float] | None) -> None: ...
    discovery: Optional[StabilityOptions]
    phases: list[PhaseSpec]
    temperature: FreeValue
    pressure: FreeValue
    feed: FeedSpec
    reactions: list[ReactionSpec]
    standard_gibbs_j_per_mol: list[float]
    topology_declared: bool
    row_tolerance: float
    enthalpy_j: float | None
    entropy_j_per_k: float | None
    critical: Optional[CriticalDeclaration]
    envelope: Optional[EnvelopeCondition]
    azeotrope: bool
    neutral_reference: Optional[NeutralReferenceDeclaration]
    def __init__(self) -> None: ...

class InventoryBound:
    species: int
    status: str
    termination: str
    maximum: Optional[float]
    feed_normalized_amount: Optional[float]
    log_fraction: Optional[float]
    balance_residual: float
    normalization_residual: float
    lower_bound_residual: float
    witness: list[float]
    amounts_mol: list[float]
    finite_witness: list[float]
    limiting: bool
    supported: bool
    floor: float
    lp_error: float
    finite_residual: float
    feasibility_termination: str

def inventory_bounds(compiled: CompiledSystem, species: Sequence[int], fraction: bool,
                     floor: float) -> list[InventoryBound]: ...

class RowKind:
    Chemical: RowKind
    Pressure: RowKind
    Moiety: RowKind
    Electro: RowKind
    Normalization: RowKind
    Enthalpy: RowKind
    Entropy: RowKind
    CriticalMode: RowKind
    CriticalNormalization: RowKind
    CriticalCubic: RowKind
    EnvelopeCondition: RowKind
    AzeotropeCondition: RowKind

class Row:
    kind: RowKind
    phase: int
    index: int
    eqid: str

class PointEvaluation:
    values: list[float]
    objective: float
    floor: float

class CompiledSystem:
    rows: list[Row]
    pinned_count: int
    free_intensive_count: int
    rank: int
    charge_row_count: int
    primal_size: int
    generic_nullity: int
    feed_total: float
    has_ions: bool
    has_reaction: bool
    conserves_total_moles: bool
    jacobian_rank: int
    removed_species: list[int]
    def generic_point(self) -> list[float]: ...
    def evaluate(self, point: Sequence[float]) -> PointEvaluation: ...
    def jacobian(self, point: Sequence[float]) -> list[float]: ...
    def lever_jacobian(self, point: Sequence[float]) -> list[float]: ...
    def objective_gradient(self, point: Sequence[float]) -> list[float]: ...
    def hessian(self, point: Sequence[float], constraint_multipliers: Sequence[float] = ...) -> list[float]: ...

class EquilibriumResult:
    solver_settings: dict[str, float]
    diagnostic_options_file: str
    solver_options: list[dict[str, Any]]
    """Solved phases in declared order (with Discover(): discovery.declaration.phases). Units: K, Pa, mol, mol/m^3.
    Per-phase-per-component lists are phase-major: entry p*n + i is phase p, component i of component_ids
    (n components). State lists are empty when success is False."""
    discovery: Optional[DiscoveryRecord]
    neutral_reference_diagnostics: list[NeutralReferenceReactionDiagnostic]
    success: bool
    message: str
    solver_status: str
    temperature: list[float]  # [T] in K
    pressure: list[float]  # [P] in Pa
    amounts: list[float]  # phase amounts in mol on the feed's basis (1 mol total for MoleFractions), per phase
    mole_fractions: list[float]  # phase-major, n per phase
    molar_densities: list[float]  # mol/m^3, per phase
    packing_fractions: list[float]  # per phase
    multipliers: list[float]
    constraint_multipliers: list[float]
    residuals: list[float]
    classified_residuals: list[float]
    coordinates: list[float]
    phase_count: int
    iterations: int
    evaluations: int
    hessian_calls: int
    objective: float
    reduced_hessian_minimum: float
    convergence: Convergence
    rows: list[Row]
    unknown_count: int
    jacobian_rank: int
    generic_nullity: int
    requested_tolerance_met: Optional[bool]
    topology_event: TopologyEvent
    topology_phase: Optional[int]
    topology_pair: Optional[tuple[int, int]]
    critical: Optional[CriticalRecord]
    def phase_mole_fractions(self, phase: int) -> list[float]: ...  # one declared phase, component_ids order

class TopologyEvent:
    None_: TopologyEvent
    PhaseDisappeared: TopologyEvent
    PhasesCoalesced: TopologyEvent

class ContinuationStop:
    reached_target: ContinuationStop
    initial_failure: ContinuationStop
    invalid_path: ContinuationStop
    singular_jacobian: ContinuationStop
    branch_lost: ContinuationStop
    phase_disappeared: ContinuationStop
    phases_coalesced: ContinuationStop
    left_domain: ContinuationStop
    infeasible: ContinuationStop
    nonfinite: ContinuationStop
    stalled: ContinuationStop
    step_too_small: ContinuationStop
    budget_exhausted: ContinuationStop
    critical_endpoint: ContinuationStop

class ContinuationOptions:
    """Ordinary path: steps are fractions of the start-to-target path; a rejected step halves and a step never
    grows, so initial_step is the largest step. Envelope mode (Problem.envelope set) ignores initial_step and
    the target problem: an ordinary step is at most envelope_max_log_step (log-coordinate units), the step after
    an accepted one is sized from the measured predictor density error, and a rejected ordinary step (branch-window
    miss, or a corrector that ends on coincident phases) halves; a rejected critical-point crossing, closing step or
    terminal failure ends the trace at once with the states accepted so far. It needs envelope_straddle > 0,
    2 * envelope_straddle <= envelope_max_log_step and T and P guesses through Free(guess). Both modes stop with
    budget_exhausted after more than max_retries_per_step consecutive rejections of retryable steps (reason in
    cause_stop and cause_message) and with step_too_small when the halved step falls below minimum_step."""
    initial_step: float
    minimum_step: float
    max_steps: int
    max_retries_per_step: int
    envelope_max_log_step: float
    envelope_straddle: float
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, *, initial_step: float, minimum_step: float, max_steps: int, max_retries_per_step: int,
                 envelope_max_log_step: float = 0.0, envelope_straddle: float = 0.0) -> None: ...

class ContinuationRequest:
    start_problem: EquilibriumProblem
    target_problem: EquilibriumProblem
    start_result: Optional[EquilibriumResult]
    options: ContinuationOptions
    parameter_targets: Sequence[float]
    def __init__(self) -> None: ...

class ContinuationResult:
    success: bool
    message: str
    stop: ContinuationStop
    path_values: list[float]  # per state: path fraction, or in envelope mode the conditions[k] target (ln T, ln P or ln K)
    parameter_values: list[list[float]]
    states: list[EquilibriumResult]
    reached_path: float
    attempted_path: Optional[float]
    event_phase: Optional[int]
    event_pair: Optional[tuple[int, int]]
    attempts: int
    accepted_steps: int
    retries: int
    last_result: Optional[EquilibriumResult]
    cause_stop: Optional[ContinuationStop]
    cause_convergence: Optional[Convergence]
    cause_message: Optional[str]
    conditions: list[EnvelopeCondition]
    critical_points: list[EquilibriumResult]
    critical_seeds: list[list[float]]

def continue_equilibrium(model: HelmholtzModel, request: ContinuationRequest) -> ContinuationResult: ...

class SolvedStateObservableKind:
    Temperature: SolvedStateObservableKind
    Pressure: SolvedStateObservableKind
    PhaseAmount: SolvedStateObservableKind
    PhaseComponentAmount: SolvedStateObservableKind
    PhaseMoleFraction: SolvedStateObservableKind
    PhaseVolume: SolvedStateObservableKind
    PhaseMolarDensity: SolvedStateObservableKind
    PhasePackingFraction: SolvedStateObservableKind
    PhasePressure: SolvedStateObservableKind
    PhaseResidualChemicalPotentialOverRt: SolvedStateObservableKind
    PhaseProperty: SolvedStateObservableKind
    PhaseLogActivity: SolvedStateObservableKind
    PhaseGibbsMolarVolume: SolvedStateObservableKind
    SystemEnthalpy: SolvedStateObservableKind
    EquilibriumIsobaricHeatCapacity: SolvedStateObservableKind

class SolvedStateObservable:
    kind: SolvedStateObservableKind
    phase: int
    component: int
    property: PropertyObservable
    def __init__(self, kind: SolvedStateObservableKind, phase: int = 0, component: int = 0,
                 property: PropertyObservable = ...) -> None: ...

class SolvedStateActionDirection:
    temperature_k: float
    pressure_pa: float
    feed_amounts_mol: list[float]
    active_parameters: list[float]
    reaction_coefficients: list[list[float]]
    enthalpy_j: float
    entropy_j_per_k: float
    def __init__(self, temperature_k: float, pressure_pa: float, feed_amounts_mol: Sequence[float],
                 active_parameters: Sequence[float], reaction_coefficients: Sequence[Sequence[float]] = ...,
                 enthalpy_j: float = ..., entropy_j_per_k: float = ...) -> None: ...

class SolvedStateActionRequest:
    observables: list[SolvedStateObservable]
    directions: list[SolvedStateActionDirection]
    def __init__(self, observables: Sequence[SolvedStateObservable],
                 directions: Sequence[SolvedStateActionDirection]) -> None: ...

class SolvedStateActionStatus:
    Available: SolvedStateActionStatus
    CentralSolveUnavailable: SolvedStateActionStatus
    DefiningEquationsNotSatisfied: SolvedStateActionStatus
    PhaseBoundaryNotDifferentiable: SolvedStateActionStatus
    NonuniqueCoalescedPhases: SolvedStateActionStatus
    RankDeficient: SolvedStateActionStatus
    NonfiniteJacobian: SolvedStateActionStatus
    LinearSolveResidual: SolvedStateActionStatus
    SupportChange: SolvedStateActionStatus
    IonicInsertion: SolvedStateActionStatus
    EnergyOrderUnavailable: SolvedStateActionStatus
    DensityUnavailable: SolvedStateActionStatus
    PressureSingular: SolvedStateActionStatus
    AssociationSingular: SolvedStateActionStatus
    AssociationIllConditioned: SolvedStateActionStatus
    AssociationResidual: SolvedStateActionStatus
    MissingThermochemistry: SolvedStateActionStatus
    ReferenceUnavailable: SolvedStateActionStatus
    OutsideInterval: SolvedStateActionStatus
    ThermalBoundary: SolvedStateActionStatus
    CentralPropertyUnavailable: SolvedStateActionStatus
    Nonfinite: SolvedStateActionStatus
    IllConditioned: SolvedStateActionStatus

class SolvedStateActionResult:
    observable: SolvedStateObservable
    unit: str
    order: int
    value: Optional[float]
    action: Optional[float]
    status: SolvedStateActionStatus
    property_diagnostics: Optional[PropertyActionDiagnostics]
    diagnostic_message: Optional[str]

class SolvedStateActionDiagnostics:
    equation_count: int
    state_dimension: int
    maximum_energy_order: int
    state_factorizations: int
    linear_solves: int
    rank: Optional[int]
    condition_estimate: Optional[float]
    maximum_linear_residual: Optional[float]
    central_residual: Optional[float]
    central_floor: Optional[float]

class SolvedStateActionBatch:
    results: list[SolvedStateActionResult]
    diagnostics: SolvedStateActionDiagnostics

def solved_state_actions(compiled: CompiledSystem, central: EquilibriumResult,
                         request: SolvedStateActionRequest) -> SolvedStateActionBatch: ...
def solved_state_first_actions(compiled: CompiledSystem, central: EquilibriumResult,
                               observables: Sequence[SolvedStateObservable],
                               directions: Sequence[SolvedStateActionDirection]) -> list[SolvedStateActionBatch]: ...
def compile_problem(model: HelmholtzModel, problem: EquilibriumProblem) -> CompiledSystem: ...
def equilibrium_asl_record(model: Any, problem: Any, continuation: ContinuationOptions, feed_offset: Sequence[float],
                           feed_matrix: Sequence[Sequence[float]], anchor_controls: Sequence[float],
                           identity_component: int, minimum_gap: float) -> str: ...
def property_asl_record(model: Any) -> str: ...
def sensitivity(compiled: CompiledSystem, point: Sequence[float], lever: Sequence[float]) -> list[float]: ...
def solve_equilibrium(
    model: HelmholtzModel,
    problem: EquilibriumProblem,
    start: Optional[EquilibriumResult] = ...,
) -> EquilibriumResult: ...

class InfiniteDilutionReference:
    solver_options: list[dict[str, Any]]
    """Ion reference fugacities at infinite dilution in the declared solvent (T in K, P in Pa)."""
    parameter_derivative: list[tuple[Optional[float], Optional[float]]]
    parameter_convergence_error: list[float]
    log_fugacity: tuple[float, float]
    molality: float
    convergence_error: float
    pure_solvent_log_fugacity: Optional[float]
    pure_solvent_state: Optional[State]
    solvent_status: PropertyActionStatus
    solvent_cause: str
    solvent_convergence: Optional[Convergence]

class SolventActivity:
    """Same-T/P pure-solvent activity and practical osmotic coefficient; None preserves unavailability."""
    solvent_activity: Optional[float]
    osmotic_coefficient: Optional[Union[float, Stop]]
    solvent_activity_parameter_derivative: list[Optional[float]]
    osmotic_parameter_derivative: list[Optional[float]]
    parameter_status: list[PropertyActionStatus]  # activity columns; osmotic columns also require empty osmotic_cause
    parameter_cause: list[str]
    status: PropertyActionStatus
    cause: str
    osmotic_cause: str
    state: State

class HenryConstant:
    """H_i^px = lim f_i/x_i = RT rho0 exp(raw_log_reference), rho0 = 1 mol/m^3; None unless status is Available."""
    component_index: int
    raw_log_reference: Optional[float]
    henry_constant_px: Optional[float]  # Pa
    observed_terminal_change: Optional[float]
    status: NeutralReferenceStatus
    message: str

def henry_constants(model: HelmholtzModel, temperature: float, pressure: float,
                    solvent_fractions: Sequence[float], solute_indices: Sequence[int]) -> list[HenryConstant]: ...

def solvent_reference(model: HelmholtzModel, temperature: float, pressure: float,
    solvent_fractions: Sequence[float]) -> InfiniteDilutionReference: ...
def solvent_activity(model: HelmholtzModel, temperature: float, pressure: float,
    composition: Sequence[float], reference: InfiniteDilutionReference) -> SolventActivity: ...

class ElectrolyteActivity:
    solver_options: list[dict[str, Any]]
    """Activity of one salt at molality m (mol per kg solvent), unsymmetric (infinite-dilution) convention."""
    single_ion: tuple[float, float]  # mole-fraction-scale coefficients of (cation, anion), the ions= order
    mean_ionic_mole_fraction_scale: float
    mean_ionic_log_molality_scale: float  # ln of mean_ionic_molality_scale
    mean_ionic_molality_scale: float  # gamma_pm on the molality scale
    osmotic_coefficient: Optional[Union[float, Stop]]  # practical (molality) osmotic coefficient

def infinite_dilution_reference(model: HelmholtzModel, temperature: float, pressure: float, solvent_fractions: Sequence[float], ion_indices: tuple[int, int], stoichiometry: tuple[int, int]) -> InfiniteDilutionReference: ...
def electrolyte_activity(model: HelmholtzModel, temperature: float, pressure: float, solvent_fractions: Sequence[float], ion_indices: tuple[int, int], stoichiometry: tuple[int, int], molality: float, reference: InfiniteDilutionReference) -> ElectrolyteActivity: ...

class ObservationKind:
    MeanIonicActivityCoefficient: ObservationKind
    OsmoticCoefficient: ObservationKind
    SolvationGibbsEnergy: ObservationKind
    TransferGibbsEnergy: ObservationKind
    BulkRelativePermittivity: ObservationKind
    LogEquilibriumConstant: ObservationKind
    IsobaricHeatCapacity: ObservationKind
    EquilibriumIsobaricHeatCapacity: ObservationKind
    Density: ObservationKind
    PhaseDensity: ObservationKind
    EquilibriumPressure: ObservationKind
    EquilibriumTemperature: ObservationKind
    PhaseMoleFraction: ObservationKind
    PhaseFraction: ObservationKind
    PartialPressure: ObservationKind
    PhaseAmountRatio: ObservationKind
    AdditionEnthalpy: ObservationKind

class ResidualForm:
    Difference: ResidualForm
    LogRatio: ResidualForm

class ObservationStatus:
    NotEvaluated: ObservationStatus
    ReferenceUnavailable: ObservationStatus
    ThermochemistryUnavailable: ObservationStatus
    Available: ObservationStatus
    DensityUnavailable: ObservationStatus
    CentralSolveUnavailable: ObservationStatus
    TopologyEvent: ObservationStatus
    MissingMolarMass: ObservationStatus
    NonpositiveLogArgument: ObservationStatus
    ActionUnavailable: ObservationStatus
    Nonfinite: ObservationStatus
    InadmissibleTrial: ObservationStatus

class Observation:
    ion_indices: tuple[int, int]
    stoichiometry: tuple[int, int]
    molality: float
    reference_composition: list[float]
    reaction: str
    kind: ObservationKind
    temperature: float
    pressure: float
    composition: list[float]
    branch: DensityBranch
    problem: EquilibriumProblem
    final_problem: EquilibriumProblem
    feed_amount: float
    phase: int
    species: list[int]
    denominator_phase: int
    denominator_species: list[int]
    mass_basis: bool
    observed: float
    form: ResidualForm
    scale: float

class Prediction:
    solver_options: list[dict[str, Any]]
    value: Optional[float]
    unit: str
    status: ObservationStatus
    message: str
    residual: Optional[float]
    convergence: Optional[Convergence]
    branch: Optional[DensityBranch]
    stable_root_count: Optional[int]
    requested_tolerance_met: Optional[bool]

class PredictionGradient:
    value: list[float]
    status: ObservationStatus
    message: str

class ReactionCoefficient:
    A: ReactionCoefficient
    B: ReactionCoefficient
    C: ReactionCoefficient
    D: ReactionCoefficient
    EnthalpyShift: ReactionCoefficient

class ReactionCoordinate:
    name: str
    coefficient: ReactionCoefficient
    pivot_temperature: float

class CoordinateTarget:
    family: ParameterFamily
    indices: list[int]
    coefficient: float

class Coordinate:
    targets: list[CoordinateTarget]
    reaction: Optional[ReactionCoordinate]
    origin: float
    scale: float
    lower: float
    upper: float

class CoordinateSet:
    coordinates: list[Coordinate]
    def active(self) -> list[ActiveParameter]: ...
    def tie(self) -> list[int]: ...
    def start(self, parameters: Parameters, observations: list[Observation] = ..., thermochemistry: Optional[ThermochemistryRecord] = ...) -> list[float]: ...

class FitControls:
    solver_options: list[dict[str, Any]]
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    maximum_iterations: int
    maximum_elapsed_time_seconds: float
    workers: int
    def __init__(self, *, maximum_iterations: int, maximum_elapsed_time_seconds: float, function_tolerance: float = 0.0,
                 gradient_tolerance: float = 0.0, parameter_tolerance: float = 0.0, workers: int = 0, solver_options: Mapping[str, float] | None = None) -> None: ...

class TrialFailure:
    trial: int
    observation: int
    optimizer: list[float]
    status: ObservationStatus
    message: str
    jacobian: bool

class FitCovariance:
    left_vectors: list[float]
    right_vectors: list[float]
    projected_gradient_norm: float
    status: str
    assumption: str
    rank: int
    singular_values: list[float]
    objective_rows: list[int]
    active_bounds: list[int]
    optimizer: list[float]
    physical: list[float]
    variance_factor: Optional[float]

class MinimaxRequest:
    targets: list[int]
    group_scales: list[float]
    q: float
    tolerance: float
    activity_distance: float
    margin: float
    def __init__(self) -> None: ...
class MinimaxRow:
    kind: str
    row: int
    group: int
    metric: CheckMetric
    sign: float
    threshold: float
    scale: float
    auxiliary: int | None
    members: list[int]
class MinimaxResult:
    availability: str
    application_status: str
    solver_status: str
    feasible: bool
    stationary: bool
    t: float
    target_maximum: float
    gap: float
    primal: float
    stationarity: float
    dual_sign: float
    complementarity: float
    point: list[float]
    lower: list[float]
    upper: list[float]
    values: list[float]
    physical_values: list[float]
    jacobian: list[float]
    multipliers: list[float]
    physical_multipliers: list[float]
    z_lower: list[float]
    z_upper: list[float]
    physical_z_lower: list[float]
    physical_z_upper: list[float]
    lower_distances: list[float]
    upper_distances: list[float]
    hold_multipliers: list[float]
    full_stationarity: list[float]
    root_errors: list[float]
    slack_groups: list[int]
    slack_rows: list[int]
    active_lower: list[int]
    active_upper: list[int]
    hold_causes: list[str]
    rows: list[MinimaxRow]
    frozen_checks: list[CheckResult]
class MeasuredInputKind:
    Temperature: MeasuredInputKind
    Pressure: MeasuredInputKind
    MoleFraction: MeasuredInputKind
    Molality: MeasuredInputKind
class MeasuredInput:
    def __init__(self) -> None: ...
    kind: MeasuredInputKind
    species: int
    sigma: float
    lower: float
    upper: float
    unit: str
class MeasurementGroup:
    def __init__(self) -> None: ...
    rows: list[int]
    inputs: list[MeasuredInput]
    dependent_species: list[int]
    feed_total_mol: float
class MeasurementUncertainty:
    def __init__(self) -> None: ...
    groups: list[MeasurementGroup]
    output_sigma: list[float]
    covariance: list[float] | None
    correlation: list[float] | None
class InputReconciliation:
    values: list[float]
    corrections: list[float]
    penalties: list[float]
    jacobian: list[float]
    row_groups: list[int]
    column_groups: list[int]
    column_inputs: list[int]
    output_cost: float
    input_cost: float
class FitResult:
    solver_options: list[dict[str, Any]]
    reconciliation: InputReconciliation | None
    minimax: MinimaxResult | None
    chart: EvaluatedCoordinates
    checks: list[CheckResult]
    acceptance: str
    optimizer: list[float]
    status: str
    message: str
    usable: bool
    # Native chart order: one value per target of each non-reaction coordinate (declaration order), then each held
    # Born diameter, then each reaction coordinate. Entry k < len(chart.active) is the value of chart.active[k].
    physical: list[float]
    parameters: Parameters
    predictions: list[float]  # per observation, in the observation's unit
    residuals: list[float]
    weighted_residuals: list[float]
    physical_jacobian: list[float]
    optimizer_jacobian: list[float]
    statuses: list[ObservationStatus]
    failures: list[TrialFailure]
    iterations: int
    residual_evaluations: int
    jacobian_evaluations: int
    trial_failures: int
    iteration_seconds: list[float]
    iteration_costs: list[float]
    initial_cost: float
    final_cost: float
    covariance: FitCovariance

def predict(model: HelmholtzModel, observation: Observation) -> Prediction: ...
def differentiate(model: HelmholtzModel, observation: Observation) -> PredictionGradient: ...
def fit(parameters: Parameters, thermochemistry: Optional[ThermochemistryRecord], coordinates: Sequence[Coordinate],
        observations: Sequence[Observation], weights: Sequence[float], controls: FitControls, holds: Sequence[HydrationHold] = ...,
        checks: Sequence[CheckGroup] = ..., baseline: Parameters | None = ..., singular_vectors: bool = ...,
        starts: Optional[Sequence[Sequence[float]]] = None, iterations: Optional[int] = None,
        seconds: Optional[float] = None, minimax: MinimaxRequest | None = ...,
        uncertainty: MeasurementUncertainty | None = ...) -> Union[FitResult, list[FitResult]]: ...

class CheckMetric:
    AbsoluteDifference: CheckMetric
    AbsoluteLog10Ratio: CheckMetric

class HydrationHold:
    def __init__(self) -> None: ...
    reference: Observation
    lower: float
    upper: float
    energy_tolerance: float
    diameter_tolerance: float
    maximum_iterations: int

class HoldResult:
    energy_error: float
    correction: float
    uncertainty: float
    derivative: float
    parameter_derivative: list[float]
    iterations: int

class EvaluatedCoordinates:
    reaction_laws: dict[str, ReactionLogPolynomial]
    parameters: Parameters
    active: list[ActiveParameter]
    physical: list[float]
    map: list[float]
    holds: list[HoldResult]

class CheckGroup:
    def __init__(self) -> None: ...
    name: str
    rows: list[int]
    metric: CheckMetric
    maximum: float | None
    mean: float | None
    baseline: bool
    report_only: bool

class CheckResult:
    name: str
    unit: str
    status: str
    rows: list[int]
    failed_rows: list[int]
    errors: list[float]
    slopes: list[float]
    maximum: float
    mean: float
    argmax: int
    maximum_threshold: float | None
    mean_threshold: float | None
    baseline_maximum: float
    baseline_mean: float
    report_only: bool

def evaluate_coordinates(parameters: Parameters, set: CoordinateSet, optimizer: Sequence[float], holds: Sequence[HydrationHold], thermochemistry: ThermochemistryRecord | None = ..., solver_options: Mapping[str, float] | None = None) -> EvaluatedCoordinates: ...
def evaluate_checks(stored: Parameters, observations: Sequence[Observation], checks: Sequence[CheckGroup], baseline: Parameters | None = ..., solver_options: Mapping[str, float] | None = None) -> list[CheckResult]: ...

def observations_at(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float], holds: Sequence[HydrationHold], observations: Sequence[Observation], errors: Sequence[float],
    thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> list[Observation]: ...

class ObservationEvaluation:
    solver_options: list[dict[str, Any]]
    reconciliation: InputReconciliation | None
    state_densities: list[list[float]]
    reference_values: list[list[float]]
    predictions: list[float]
    gradient: list[float]
    optimizer_gradient: list[float]
    residual_optimizer_gradient: list[float]
    chart: EvaluatedCoordinates
    statuses: list[ObservationStatus]
def evaluate_observations(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
    holds: Sequence[HydrationHold], observations: Sequence[Observation],
    thermochemistry: Optional[ThermochemistryRecord] = None,
    uncertainty: MeasurementUncertainty | None = ..., solver_options: Mapping[str, float] | None = None) -> ObservationEvaluation: ...
def evaluate_minimax(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
    holds: Sequence[HydrationHold], observations: Sequence[Observation], checks: Sequence[CheckGroup],
    baseline: Parameters | None, minimax: MinimaxRequest, thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> FitResult: ...
def formula_composition(solvent_fractions: Sequence[float], parameters: Parameters, ion_indices: tuple[int, int],
    stoichiometry: tuple[int, int], molality: float) -> list[float]: ...

class EquilibriumStart:
    coordinates: list[float]
    amounts: list[float]
    compositions: list[list[float]]
    packing_fractions: list[Optional[float]]
    temperature: Optional[float]
    pressure: Optional[float]
    previous: Optional[EquilibriumResult]
    def __init__(self, *, coordinates: Sequence[float] = ..., amounts: Sequence[float] = ...,
                 compositions: Sequence[Sequence[float]] = ..., packing_fractions: Sequence[Optional[float]] = ...,
                 temperature: Optional[float] = None, pressure: Optional[float] = None,
                 previous: Optional[EquilibriumResult] = None) -> None: ...

class EquilibriumWorkLimits:
    attempts: Optional[int]
    iterations: Optional[int]
    attempt_iterations: Optional[int]
    seconds: Optional[float]
    attempt_seconds: Optional[float]
    def __init__(self, *, attempts: Optional[int] = None, iterations: Optional[int] = None,
                 attempt_iterations: Optional[int] = None, seconds: Optional[float] = None,
                 attempt_seconds: Optional[float] = None) -> None: ...

class EquilibriumStartsResult:
    selected: Optional[int]
    results: list[Optional[EquilibriumResult]]
    causes: list[str]
    attempt_seconds: list[float]
    iterations: int
    seconds: float
    overrun_seconds: float
    budget_exhausted: bool

def solve_equilibrium_starts(model: HelmholtzModel, problem: EquilibriumProblem, starts: list[EquilibriumStart], limits: EquilibriumWorkLimits) -> EquilibriumStartsResult: ...

def reaction_log10_reciprocal(law: ReactionLogPolynomial) -> list[float]: ...
