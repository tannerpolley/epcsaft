from collections.abc import Sequence
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
    born_ssm_ds: float
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
    def _polar_terms(self, *args: object) -> tuple[float, ...]: ...
    def _reduced_multipole_moments(self, *args: object) -> tuple[float, ...]: ...

def native_sdk_v1(*args: object, **kwargs: object) -> object: ...
def source_reference_transfer_v1(*args: object) -> dict[str, Any]: ...
def _polar_constants() -> dict[str, float]: ...
def _polar_pade(second: float, third: float, active: bool) -> float: ...
