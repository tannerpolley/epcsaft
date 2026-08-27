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
    def with_provider_metadata(
        self,
        parameter_fingerprint: str,
        temperature_min_k: float,
        temperature_max_k: float,
        ion_mole_fraction_max: float,
        pressure_min_pa: float,
        pressure_max_pa: float,
        topology_fingerprint: str,
        component_order_fingerprint: str,
        domain_fingerprint: str,
        artifact_fingerprint: str,
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

def source_reference_transfer(*args: object) -> dict[str, Any]: ...
def evaluate_parameterized_phase(
    model: NativeModel,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    packing_fraction_min: float,
    packing_fraction_max: float,
) -> dict[str, object]: ...
def _audit_phase_state(
    model: NativeModel,
    role: str,
    temperature_k: float,
    mole_fractions: tuple[float, ...],
    molar_volume_m3_per_mol: float,
) -> dict[str, object]: ...
def evaluate_solved_state_sensitivity(
    model: NativeModel,
    temperature_k: float,
    pressure_pa: float,
    feed_amounts_mol: tuple[float, ...],
    phase_amounts_mol: tuple[tuple[float, ...], ...],
    phase_volumes_m3: tuple[float, ...],
) -> dict[str, object]: ...
def _solve_declared_fixed_tp(
    model: NativeModel,
    temperature_k: float,
    pressure_pa: float,
    feed_amounts_mol: tuple[float, ...],
    roles: tuple[str, ...],
    supports: tuple[str, ...],
    hint_indices: tuple[int, ...],
    hint_compositions: tuple[tuple[float, ...], ...],
    hint_densities_mol_m3: tuple[float, ...],
    initial_guess_indices: tuple[int, ...],
    initial_guess_compositions: tuple[tuple[float, ...], ...],
    initial_guess_densities_mol_m3: tuple[float, ...],
    initial_guess_phase_fractions: tuple[float, ...],
) -> dict[str, Any]: ...
def _solve_fixed_t_coexistence(
    reference_model: NativeModel,
    partner_model: NativeModel,
    temperature_k: float,
    roles: tuple[str, ...],
    supports: tuple[str, ...],
    reference_phase_index: int,
    reference_mole_fractions: tuple[float, ...],
    reference_direction_identities: tuple[str, ...],
    reference_composition_directions: tuple[tuple[float, ...], ...],
    hint_indices: tuple[int, ...],
    hint_compositions: tuple[tuple[float, ...], ...],
    hint_densities_mol_m3: tuple[float, ...],
) -> dict[str, Any]: ...
def _solve_held_discovery(
    model: NativeModel,
    temperature_k: float,
    pressure_pa: float,
    overall_mole_fractions: tuple[float, ...],
) -> dict[str, Any]: ...
def _held2_audit_farkas_certificate(
    matrix: list[list[float]],
    row_lower: list[float],
    row_upper: list[float],
    column_lower: list[float],
    column_upper: list[float],
    row_ray: list[float],
) -> dict[str, Any]: ...
def _held2_audit_step5_kkt(
    variables: list[float],
    lower: list[float],
    upper: list[float],
    objective_gradient: list[float],
    constraints: list[list[float]],
    constraint_upper: list[float],
    lower_bound_multipliers: list[float],
    upper_bound_multipliers: list[float],
    constraint_multipliers: list[float],
    pullback_residual: float,
    pullback_scale: float,
    pressure_residual: float,
    same_major_iteration: bool,
    step4_binding_valid: bool,
    pressure_branch_valid: bool,
) -> dict[str, Any]: ...
def _chemical_equilibrium(
    model: NativeModel | None,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_bounds: tuple[float, float] | None,
    trace_floor: float,
    sensitivities_requested: bool = ...,
    active_parameters_requested: bool = ...,
    max_iterations: int = ...,
) -> dict[str, Any]: ...
def _evaluate_provider_reaction_local(
    model: NativeModel,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_bounds: tuple[float, float],
    variables: Sequence[float],
    lower: Sequence[float],
    upper: Sequence[float],
    active_parameters_requested: bool = ...,
) -> dict[str, Any]: ...
def _reactive_finite_phases(
    model: NativeModel,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_intervals: tuple[tuple[float, float], ...],
    starts: tuple[dict[str, object], ...],
    trace_floor: float,
    active_parameters_requested: bool = ...,
    derivative_audits_requested: bool = ...,
    phase_direction_system: dict[str, object] | None = ...,
    phase_support_indices: tuple[tuple[int, ...], ...] = ...,
) -> dict[str, Any]: ...
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
    full_space_evaluator: Callable[[int, Sequence[float]], Mapping[str, object]] | None = ...,
    full_space_setups: Sequence[Mapping[str, object]] = ...,
) -> dict[str, Any]: ...
