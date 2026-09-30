from typing import Any, Optional, Sequence, Union

GAS_CONSTANT_J_PER_MOL_K: float
EOS_STANDARD_STATE_ID: str

class ParameterError(ValueError): ...
class EosError(RuntimeError): ...
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

class Span:
    name: str
    parent: Optional[int]
    convergence: Optional[Convergence]
    wall_us: float

class Trace:
    spans: list[Span]
    events: list[Event]
    current: Optional[int]

class Capture:
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
    SegmentCount: ParameterFamily
    SigmaConstant: ParameterFamily
    DispersionEnergyOverK: ParameterFamily
    Kij: ParameterFamily
    Lij: ParameterFamily
    AssociationEnergyOverK: ParameterFamily
    AssociationVolume: ParameterFamily
    KijReciprocalTemperatureSlope: ParameterFamily

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
    value: Optional[float]
    action: Optional[float]
    status: PropertyActionStatus
    energy_order: int

class DensityBranch:
    Unspecified: DensityBranch
    Vapor: DensityBranch
    Liquid: DensityBranch

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
    dense_hessian_calls: int
    branch: DensityBranch
    stable_root_count: Optional[int]

class PropertyActionBatch:
    results: list[PropertyActionResult]
    maximum_energy_order: int
    diagnostics: PropertyActionDiagnostics
    molar_density: Optional[float]
    density_action: list[float]
    density_status: PropertyActionStatus
    density_message: str
    density_convergence: Optional[Convergence]

class State:
    hard_chain: float
    dispersion: float
    association: float
    debye_huckel: float
    born: float
    bulk_relative_permittivity: float
    residual_helmholtz: float
    compressibility_factor: float
    pressure: float
    residual_chemical_potential_over_rt: list[float]
    log_fugacity_coefficient: Optional[list[float]]
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
    molar_density: float
    pressure_density_derivative: float
    site_fractions: list[float]
    density_convergence: Optional[Convergence]
    thermochemistry_status: ThermochemistryStatus
    density_branch: DensityBranch
    stable_root_count: Optional[int]

class SweepCounters:
    forward: int
    reverse: int
    hessian_calls: int
    recordings: int

class HelmholtzModel:
    def __init__(self, parameters: Parameters, active_parameters: Sequence[ActiveParameter] = (),
                 thermochemistry: Optional[ThermochemistryRecord] = None) -> None: ...
    def state(self, temperature: float, molar_density: float, mole_fractions: Sequence[float]) -> State: ...
    def neutral_property(self, temperature: float, molar_density: float, property: NeutralProperty, order: DerivativeOrder) -> NeutralPropertyDerivatives: ...
    def density_actions(self, temperature: float, molar_density: float,
                        mole_fractions: Sequence[float], observables: Sequence[PropertyObservable],
                        directions: Sequence[DensityActionDirection]) -> PropertyActionBatch: ...
    def pressure_actions(self, temperature: float, pressure: float,
                         mole_fractions: Sequence[float], branch: DensityBranch,
                         observables: Sequence[PropertyObservable],
                         directions: Sequence[PressureActionDirection],
                         anchor_density: Optional[float] = None) -> PropertyActionBatch: ...
    def state_at_pressure(self, temperature: float, pressure: float, mole_fractions: Sequence[float], branch: DensityBranch, anchor_density: Optional[float] = None) -> State: ...
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

class NeutralReferenceDeclaration:
    solvent_index: Optional[int]
    solute_molality_ratios: list[float]
    standard_molality_mol_per_kg: float
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

class EquilibriumProblem:
    phases: list[PhaseSpec]
    temperature: FreeValue
    pressure: FreeValue
    feed: FeedSpec
    reactions: list[ReactionSpec]
    standard_gibbs_j_per_mol: list[float]
    topology_declared: bool
    row_tolerance: float
    neutral_reference: Optional[NeutralReferenceDeclaration]
    def __init__(self) -> None: ...

class RowKind:
    Chemical: RowKind
    Pressure: RowKind
    Moiety: RowKind
    Electro: RowKind
    Normalization: RowKind

class Row:
    kind: RowKind
    phase: int
    index: int
    eqid: str

class PointEvaluation:
    values: list[float]
    jacobian: list[float]
    objective: float
    floor: float

class CompiledSystem:
    rows: list[Row]
    pinned_count: int
    free_intensive_count: int
    rank: int
    charge_row_count: int
    primal_size: int
    rows_minus_unknowns: int
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
    neutral_reference_diagnostics: list[NeutralReferenceReactionDiagnostic]
    success: bool
    message: str
    solver_status: str
    temperature: list[float]
    pressure: list[float]
    amounts: list[float]
    mole_fractions: list[float]
    molar_densities: list[float]
    packing_fractions: list[float]
    multipliers: list[float]
    constraint_multipliers: list[float]
    residuals: list[float]
    classified_residuals: list[float]
    sensitivity: list[float]
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

class ContinuationOptions:
    initial_step: float
    minimum_step: float
    max_steps: int
    max_retries_per_step: int
    def __init__(self) -> None: ...

class ContinuationRequest:
    start_problem: EquilibriumProblem
    target_problem: EquilibriumProblem
    start_result: Optional[EquilibriumResult]
    options: ContinuationOptions
    def __init__(self) -> None: ...

class ContinuationResult:
    success: bool
    message: str
    stop: ContinuationStop
    path_values: list[float]
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
    def __init__(self, temperature_k: float, pressure_pa: float, feed_amounts_mol: Sequence[float],
                 active_parameters: Sequence[float]) -> None: ...

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
def equilibrium_asl_record(model: Any, problem: Any, continuation: Any, feed_offset: Sequence[float],
                           feed_matrix: Sequence[Sequence[float]], anchor_controls: Sequence[float],
                           identity_component: int, minimum_gap: float) -> str: ...
def sensitivity(compiled: CompiledSystem, point: Sequence[float], lever: Sequence[float]) -> list[float]: ...
def solve_equilibrium(
    model: HelmholtzModel,
    problem: EquilibriumProblem,
    start: Optional[EquilibriumResult] = ...,
) -> EquilibriumResult: ...

class SensitivityInfo:
    equations: int
    parameters: int
    rank: int
    condition_estimate: float
    nonsingular: bool

class InfiniteDilutionReference:
    log_fugacity: tuple[float, float]
    molality: float
    convergence_error: float
    pure_solvent_log_fugacity: Optional[float]

class ElectrolyteActivity:
    single_ion: tuple[float, float]
    mean_ionic_mole_fraction_scale: float
    mean_ionic_log_molality_scale: float
    mean_ionic_molality_scale: float
    osmotic_coefficient: Optional[float]

def infinite_dilution_reference(model: HelmholtzModel, temperature: float, pressure: float, solvent_fractions: Sequence[float], ion_indices: tuple[int, int], stoichiometry: tuple[int, int]) -> InfiniteDilutionReference: ...
def electrolyte_activity(model: HelmholtzModel, temperature: float, pressure: float, solvent_fractions: Sequence[float], ion_indices: tuple[int, int], stoichiometry: tuple[int, int], molality: float, reference: InfiniteDilutionReference) -> ElectrolyteActivity: ...

class ObservationKind:
    IsobaricHeatCapacity: ObservationKind
    Density: ObservationKind
    PhaseDensity: ObservationKind
    EquilibriumPressure: ObservationKind
    EquilibriumTemperature: ObservationKind
    PhaseMoleFraction: ObservationKind
    PhaseFraction: ObservationKind
    PartialPressure: ObservationKind
    PhaseAmountRatio: ObservationKind

class ResidualForm:
    Difference: ResidualForm
    LogRatio: ResidualForm

class ObservationStatus:
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
    kind: ObservationKind
    temperature: float
    pressure: float
    composition: list[float]
    branch: DensityBranch
    problem: EquilibriumProblem
    phase: int
    species: list[int]
    denominator_phase: int
    denominator_species: list[int]
    mass_basis: bool
    observed: float
    form: ResidualForm
    scale: float

class Prediction:
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

class ReactionCoordinate:
    name: str
    coefficient: ReactionCoefficient

class Coordinate:
    family: ParameterFamily
    slots: list[list[int]]
    reaction: Optional[ReactionCoordinate]
    origin: float
    scale: float
    lower: float
    upper: float

class CoordinateSet:
    coordinates: list[Coordinate]
    def active(self) -> list[ActiveParameter]: ...
    def tie(self) -> list[int]: ...

class FitControls:
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    maximum_iterations: int
    maximum_elapsed_time_seconds: float

class TrialFailure:
    trial: int
    observation: int
    optimizer: list[float]
    status: ObservationStatus
    message: str
    jacobian: bool

class FitCovariance:
    status: str
    assumption: str
    rank: int
    singular_values: list[float]
    active_bounds: list[int]
    optimizer: list[float]
    physical: list[float]
    variance_factor: Optional[float]

class FitResult:
    status: str
    message: str
    usable: bool
    physical: list[float]
    parameters: Parameters
    predictions: list[float]
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
        observations: Sequence[Observation], weights: Sequence[float], controls: FitControls) -> FitResult: ...
