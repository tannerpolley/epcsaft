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
