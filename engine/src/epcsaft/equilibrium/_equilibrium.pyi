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
def _phase_equilibrium(
    capsule: object,
    temperature_k: float,
    pressure_pa: float,
    feed_amounts_mol: tuple[float, ...],
    mode: str,
    constrained_roles: tuple[str, ...],
    hint_roles: tuple[str, ...],
    hint_mole_fractions: tuple[tuple[float, ...], ...],
    hint_molar_densities_mol_m3: tuple[float, ...],
) -> dict[str, Any]: ...
def _chemical_equilibrium(
    capsule: object | None,
    spec: dict[str, object],
    source_standard_state: dict[str, object] | None,
    packing_fraction_bounds: tuple[float, float] | None,
    trace_floor: float,
    active_parameters_requested: bool = ...,
) -> dict[str, Any]: ...
