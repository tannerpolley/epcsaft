from collections.abc import Mapping, Sequence
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
    def evaluate_fixed_pressure_composition_derivatives(
        self,
        temperature: float,
        mole_fractions: Sequence[float],
        molar_volume: float,
        ion_mole_fraction_max: float,
    ) -> Mapping[str, Any]: ...
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

class _ReactiveFinitePhaseResult:
    accepted: bool
    solver_status: str
    callback_error: str
    callback_failure_kind: str
    failure_reason: str
    numerical_status: str
    physical_status: str
    provider_domain_status: str
    resolved_policy: Sequence[tuple[str, float]]
    phase_amounts_mol: Sequence[Sequence[float]]
    phase_volumes_m3: Sequence[float]
    phase_pressures_pa: Sequence[float]
    phase_packing_fractions: Sequence[float]
    phase_chemical_potentials_over_rt: Sequence[Sequence[float]]
    phase_log_fugacities: Sequence[Sequence[float]]
    local_variables: Sequence[float]
    local_lower_bounds: Sequence[float]
    local_upper_bounds: Sequence[float]
    input_order: Sequence[str]
    kkt_input_jacobian: Sequence[float]
    input_amount_derivatives: Sequence[float]
    input_volume_derivatives: Sequence[float]
    input_pressure_derivatives: Sequence[float]
    input_chemical_potential_derivatives: Sequence[float]
    input_log_fugacity_derivatives: Sequence[float]
    kkt_root_values: Sequence[float]
    kkt_root_jacobian: Sequence[float]
    kkt_dimension: int
    kkt_rank: int
    condition_number_inf: float
    output_values: Sequence[float]
    output_input_derivatives: Sequence[float]
    output_local_derivatives: Sequence[float]
    output_parameter_partials: Sequence[float]
    objective_over_rt: float
    balance_inf_norm: float
    pressure_relative_inf_norm: float
    reaction_affinity_inf_norm: float
    transfer_chemical_potential_inf_norm: float
    kkt_stationarity_inf_norm: float
    phase_distinctness: float
    minimum_amount_mol: float
    minimum_bound_margin: float
    reduced_hessian_min_eigenvalue: float
    reduced_hessian_scale: float
    reduced_hessian_threshold: float
    declared_direction_rank: int
    admissible_nullity: int
    optimizer_iterations: int
    provider_evaluation_count: int
    derivative_audits_performed: bool
    derivative_audit_status: str
    objective_gradient_audit_max_relative_error: float
    constraint_jacobian_audit_max_relative_error: float
    lagrangian_hessian_audit_max_relative_error: float
    parameter_order: Sequence[str]
    sensitivity_status: str
    sensitivity_failure_reason: str
    amount_derivatives: Sequence[float]
    volume_derivatives: Sequence[float]
    pressure_derivatives: Sequence[float]
    chemical_potential_derivatives: Sequence[float]
    log_fugacity_derivatives: Sequence[float]

class _GrepeResult:
    formulation: str
    certificate_scope: str
    equilibrium: _ReactiveFinitePhaseResult
    failure_code: str
    temperature_k: float
    pressure_pa: float
    source_reference_representation_residual_inf_norm: float
    evidence: Sequence[tuple[str, float]]
    branch_kind: str
    branch_coordinates: Sequence[float]
    phase_mechanical_classes: Sequence[str]

def _solve_grepe(
    model: NativeModel,
    incipient_model: NativeModel,
    problem: object,
    active_parameters: object | None = ...,
    source_reference_fingerprint: str = ...,
    perform_derivative_audits: bool = ...,
    declared_starts_only: bool = ...,
) -> _GrepeResult: ...
def _initialize_grepe_local(
    model: NativeModel,
    problem: object,
    active_parameters: object | None,
    source_reference_fingerprint: str,
) -> _GrepeResult: ...
def _evaluate_grepe_local(
    model: NativeModel,
    problem: object,
    active_parameters: object | None,
    source_reference_fingerprint: str,
    local_variables: Sequence[float],
    local_lower_bounds: Sequence[float],
    local_upper_bounds: Sequence[float],
) -> _GrepeResult: ...
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
    active_parameter_count: int,
    identity_coordinates: Sequence[int],
    coordinates: Sequence[Mapping[str, object]],
    observations: Sequence[Mapping[str, object]],
    controls: Mapping[str, object],
    failure_codes: Sequence[str],
    experiments: Sequence[Mapping[str, object]],
) -> dict[str, Any]: ...
