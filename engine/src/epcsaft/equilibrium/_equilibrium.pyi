from typing import Any

def transport_info(capsule: object) -> dict[str, object]: ...
def evaluate_phase(
    capsule: object,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
) -> dict[str, object]: ...
def evaluate_mixture_phase(
    capsule: object,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
) -> dict[str, object]: ...
def evaluate_parameterized_phase(
    capsule: object,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    packing_fraction_min: float,
    packing_fraction_max: float,
) -> dict[str, object]: ...
def evaluate_solved_state_sensitivity(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    feed_amounts_mol: tuple[float, ...],
    phase_amounts_mol: tuple[tuple[float, ...], ...],
    phase_volumes_m3: tuple[float, ...],
) -> dict[str, object]: ...
def _solve_tp_flash(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    overall_mole_fractions: tuple[float, ...],
    *,
    trace: bool = ...,
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
    active_parameters_requested: bool = ...,
) -> dict[str, Any]: ...
