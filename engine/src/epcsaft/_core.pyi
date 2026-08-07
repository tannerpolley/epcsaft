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

class NativeIonSolvationDerivativeBlock:
    solvation_gibbs_j_per_mol: Sequence[float]
    parameter_derivatives_j_per_mol: Sequence[float]
    reference_molality_mol_per_kg: float
    reference_convergence_error: float

class NativeAqueousMiacSolvationFactor:
    log_mean_ionic_activity_coefficient_molality: float
    derivative: float

class NativeIonFractionSuppression:
    bulk_relative_permittivity: float
    derivative: float

class NativeParameterRecord:
    identity: str
    value: float
    def __init__(self, identity: str, value: float) -> None: ...

class NativeFixedTopologyPhaseBlock:
    helmholtz_over_rt_reference_amount: float
    gradient: Sequence[float]
    hessian: Sequence[float]
    pressure_pa: float
    chemical_potential_over_rt: float

class NativeMixturePhaseBlock:
    helmholtz_over_rt_reference_amount: float
    gradient: Sequence[float]
    hessian: Sequence[float]
    pressure_pa: float

class NativeFugacityDerivativeBlock:
    log_fugacity_coefficient: Sequence[float]
    pressure_derivative_per_pa: Sequence[float]
    parameter_derivatives: Sequence[float]
    molar_density_mol_per_m3: float
    root_residual_pa: float
    pressure_density_derivative_pa_m3_per_mol: float
    density_condition_number: float
    branch: str
    bracket_lower_mol_per_m3: float
    bracket_upper_mol_per_m3: float
    stable_root_count: int
    selected_stable_root_index: int

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
    def evaluate_fixed_topology_association_active_parameters(
        self,
        temperature: float,
        amount: float,
        volume: float,
        active_parameters: Sequence[NativeParameterRecord],
    ) -> NativeFixedTopologyPhaseBlock: ...
    def evaluate_fixed_topology_associating_mixture_active_parameters(
        self,
        temperature: float,
        amounts: Sequence[float],
        volume: float,
        active_parameters: Sequence[NativeParameterRecord],
    ) -> NativeMixturePhaseBlock: ...
    def evaluate_mixture_phase_active_parameter(
        self,
        temperature: float,
        amounts: Sequence[float],
        volume: float,
        active_parameter: NativeParameterRecord,
    ) -> NativeMixturePhaseBlock: ...
    def evaluate_fugacity_derivatives_active(
        self,
        temperature: float,
        pressure: float,
        composition: Sequence[float],
        phase: str,
        active_parameters: Sequence[NativeParameterRecord],
    ) -> NativeFugacityDerivativeBlock: ...
    def evaluate_ion_solvation_derivatives(
        self,
        temperature: float,
        pressure: float,
        solvent_fractions: Sequence[float],
        ion_indices: Sequence[int],
        stoichiometry: Sequence[int],
        phase: str,
        active_parameters: Sequence[NativeParameterRecord],
    ) -> NativeIonSolvationDerivativeBlock: ...
    def evaluate_aqueous_miac_solvation_factor(
        self, temperature: float, pressure: float, molality: float, factor: float
    ) -> NativeAqueousMiacSolvationFactor: ...
    def evaluate_ion_fraction_suppression(
        self, ion_mole_fraction: float, coefficient: float
    ) -> NativeIonFractionSuppression: ...
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
    liquid_density_upper_mol_m3: float,
    active_parameters: Sequence[NativeParameterRecord],
) -> dict[str, Any]: ...
def _fit_resolved(
    model: NativeModel,
    identities: Sequence[str],
    identity_coordinates: Sequence[int],
    coordinates: Sequence[Mapping[str, object]],
    observations: Sequence[Mapping[str, object]],
    controls: Mapping[str, object],
) -> dict[str, Any]: ...
