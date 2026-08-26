from collections.abc import Sequence
from typing import Any

def evaluate_parameterized_phase(
    capsule: object,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    packing_fraction_min: float,
    packing_fraction_max: float,
) -> dict[str, object]: ...
def _audit_phase_state(
    capsule: object,
    role: str,
    temperature_k: float,
    mole_fractions: tuple[float, ...],
    molar_volume_m3_per_mol: float,
) -> dict[str, object]: ...
def evaluate_solved_state_sensitivity(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    feed_amounts_mol: tuple[float, ...],
    phase_amounts_mol: tuple[tuple[float, ...], ...],
    phase_volumes_m3: tuple[float, ...],
) -> dict[str, object]: ...
def _solve_declared_fixed_tp(
    capsule: object,
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
    reference_capsule: object,
    partner_capsule: object,
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
    capsule: object,
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
    capsule: object | None,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_bounds: tuple[float, float] | None,
    trace_floor: float,
    sensitivities_requested: bool = ...,
    active_parameters_requested: bool = ...,
) -> dict[str, Any]: ...
def _evaluate_provider_reaction_local(
    capsule: object,
    spec: dict[str, object],
    packing_fraction_bounds: tuple[float, float],
    trace_floor: float,
    variables: Sequence[float],
    lower: Sequence[float],
    upper: Sequence[float],
    active_parameters_requested: bool = ...,
) -> dict[str, Any]: ...
def _reactive_finite_phases(
    capsule: object,
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
