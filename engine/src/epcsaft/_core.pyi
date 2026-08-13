from collections.abc import Callable, Mapping, Sequence
from typing import Any

class NativeDensityDiagnostics:
    bracket_lower: float
    bracket_upper: float
    pressure_residual: float
    pressure_density_derivative: float
    iterations: int
    stable: bool
    branch: str
    stable_root_count: int
    selected_stable_root_index: int

class NativeResult:
    hard_chain: float
    dispersion: float
    association: float
    ares_polar: float
    debye_huckel: float
    born: float
    bulk_relative_permittivity: float
    residual_helmholtz: float
    compressibility_factor: float
    pressure: float
    residual_chemical_potential_over_rt: Sequence[float]
    log_fugacity_coefficient: Sequence[float] | None
    residual_entropy: float | None
    residual_internal_energy: float | None
    residual_enthalpy: float
    residual_isochoric_heat_capacity: float | None
    residual_isobaric_heat_capacity: float | None
    molar_density: float
    density_diagnostics: NativeDensityDiagnostics | None

class NativeReferenceSequence:
    log_fugacity: Sequence[float]
    pure_solvent_log_fugacity_coefficient: float | None
    reference_molality: float
    convergence_error: float

class NativeActiveCoordinate:
    identity: str
    family: str
    application: str
    indices: Sequence[int]
    association_edges: Sequence[tuple[int, int]]
    association_energy_scales: Sequence[float]
    value: float
    def __init__(
        self,
        identity: str,
        family: str,
        application: str,
        indices: Sequence[int],
        association_edges: Sequence[tuple[int, int]],
        association_energy_scales: Sequence[float],
        value: float,
    ) -> None: ...

class NativeMixturePhaseBlock:
    helmholtz_over_rt_reference_amount: float
    gradient: Sequence[float]
    hessian: Sequence[float]
    pressure_pa: float

class NativeModel:
    def __init__(self, *args: object) -> None: ...
    def evaluate(
        self, temperature_k: float, molar_density: float, composition: Sequence[float]
    ) -> NativeResult: ...
    def evaluate_pressure(
        self,
        temperature_k: float,
        pressure_pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> NativeResult: ...
    def with_parameter_identities(
        self, component_ids: Sequence[str], site_ids: Sequence[str]
    ) -> NativeModel: ...
    def with_active_coordinates(
        self, coordinates: Sequence[NativeActiveCoordinate]
    ) -> NativeModel: ...
    def evaluate_phase_derivatives(
        self,
        temperature: float,
        amounts: Sequence[float],
        volume: float,
        coordinates: Sequence[NativeActiveCoordinate],
        derivative_order: int,
    ) -> NativeMixturePhaseBlock: ...
    def molar_volume_bounds(
        self,
        temperature_k: float,
        composition: Sequence[float],
        packing_fraction_min: float,
        packing_fraction_max: float,
        ion_mole_fraction_max: float,
    ) -> Sequence[float]: ...
    def molar_density_upper_bound(
        self, temperature_k: float, composition: Sequence[float]
    ) -> float: ...
    def reference_composition(
        self,
        solvent_fractions: Sequence[float],
        ion_indices: tuple[int, int],
        stoichiometry: tuple[int, int],
        molality: float,
    ) -> Sequence[float]: ...
    def electrolyte_reference_sequence(
        self,
        temperature_k: float,
        pressure_pa: float,
        solvent_fractions: Sequence[float],
        ion_indices: tuple[int, int],
        stoichiometry: tuple[int, int],
        phase: str,
    ) -> NativeReferenceSequence: ...

def native_transport(*args: object, **kwargs: object) -> object: ...
def source_reference_transfer(*args: object) -> dict[str, Any]: ...
def _evaluate_saturation_nlp(
    model: NativeModel,
    temperature_k: float,
    variables: tuple[float, float, float],
    multipliers: tuple[float, float, float],
) -> dict[str, Any]: ...
def _solve_saturation(
    model: NativeModel,
    temperature_k: float,
    liquid_density_upper_mol_m3: float,
) -> dict[str, Any]: ...
def _solve_saturation_active(
    model: NativeModel,
    temperature_k: float,
    active_parameters: Sequence[NativeActiveCoordinate],
) -> dict[str, Any]: ...
def _fit_resolved(
    owner_evaluator: Callable[[int, int, Sequence[float], bool], Mapping[str, object]],
    active_parameter_count: int,
    identity_coordinates: Sequence[int],
    coordinates: Sequence[Mapping[str, object]],
    observations: Sequence[Mapping[str, object]],
    controls: Mapping[str, object],
    failure_codes: Sequence[str],
) -> dict[str, Any]: ...
