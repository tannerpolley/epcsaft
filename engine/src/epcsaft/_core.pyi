from typing import Optional, Sequence

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
    electrolyte_choice: str
    ionic_dispersion_choice: str
    permittivity_choice: str
    entries: list[ParameterEntry]
    def __init__(self) -> None: ...

class Parameters:
    @property
    def component_count(self) -> int: ...
    @property
    def site_count(self) -> int: ...

def build_parameters(document: ParameterDocument) -> Parameters: ...

class DensityBranch:
    Unspecified: DensityBranch
    Vapor: DensityBranch
    Liquid: DensityBranch

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
    molar_density: float
    pressure_density_derivative: float
    site_fractions: list[float]
    density_convergence: Optional[Convergence]

class SweepCounters:
    forward: int
    reverse: int
    hessian_calls: int
    recordings: int

class HelmholtzModel:
    def __init__(self, parameters: Parameters) -> None: ...
    def state(self, temperature: float, molar_density: float, mole_fractions: Sequence[float]) -> State: ...
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

class ReactionSpec:
    name: str
    stoichiometry: list[float]
    equilibrium_constant: float
    affinity_j_per_mol: Optional[float]
    def __init__(self) -> None: ...

class EquilibriumProblem:
    phases: list[PhaseSpec]
    temperature: FreeValue
    pressure: FreeValue
    feed: FeedSpec
    reactions: list[ReactionSpec]
    standard_gibbs_j_per_mol: list[float]
    topology_declared: bool
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
    jacobian_rank: int
    removed_species: list[int]
    def generic_point(self) -> list[float]: ...
    def evaluate(self, point: Sequence[float]) -> PointEvaluation: ...
    def jacobian(self, point: Sequence[float]) -> list[float]: ...
    def lever_jacobian(self, point: Sequence[float]) -> list[float]: ...
    def objective_gradient(self, point: Sequence[float]) -> list[float]: ...
    def hessian(self, point: Sequence[float]) -> list[float]: ...

class EquilibriumResult:
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
    residuals: list[float]
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

def compile_problem(model: HelmholtzModel, problem: EquilibriumProblem) -> CompiledSystem: ...
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
