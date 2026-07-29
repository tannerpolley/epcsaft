import ctypes
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import gc
import math
from pathlib import Path
import shutil
import subprocess
import time

import pytest

import epcsaft
from epcsaft import (
    Mixture,
    EosError,
    _ParameterBundle,
    export_native_model,
    native_sdk,
    unit_registry,
)
from epcsaft._core import NativeModel, native_sdk_v1 as _native_sdk_v1
from epcsaft.eos import _native_arguments
from epcsaft.records import (
    AssociationParameterRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    ValidityDomain,
)


CAPSULE_NAME = b"epcsaft.native_sdk.v1"
ABI_VERSION = 1
STATUS_OK = 0
STATUS_INVALID_ARGUMENT = 1
STATUS_UNSUPPORTED_MODEL = 2
STATUS_DOMAIN_ERROR = 3
STATUS_EVALUATION_ERROR = 4
STATUS_ABI_MISMATCH = 5
STATUS_CANCELLED = 6
STATUS_DEADLINE_EXCEEDED = 7


def test_resolved_model_artifact_is_native_loadable(tmp_path: Path) -> None:
    model = Mixture(
        _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
            ("methane", "ethane")
        )
    )
    artifact = tmp_path / "model.json"

    export_native_model(model, artifact)

    payload = artifact.read_text(encoding="utf-8")
    assert '"schema":"epcsaft-resolved-native-model-v1"' in payload
    assert model.parameter_fingerprint in payload


class PhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("gradient", ctypes.c_double * 2),
        ("hessian", ctypes.c_double * 4),
        ("third", ctypes.c_double * 8),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class ParameterizedPhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("gradient", ctypes.c_double * 5),
        ("hessian", ctypes.c_double * 25),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class MixturePhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("coordinate_count", ctypes.c_size_t),
        ("gradient_capacity", ctypes.c_size_t),
        ("hessian_capacity", ctypes.c_size_t),
        ("gradient", ctypes.POINTER(ctypes.c_double)),
        ("hessian", ctypes.POINTER(ctypes.c_double)),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AcceptedNativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
    )


class FixedMixtureNativeSdkTable(ctypes.Structure):
    _fields_ = AcceptedNativeSdkTable._fields_ + (
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
    )


class NativeSdkTable(ctypes.Structure):
    _fields_ = FixedMixtureNativeSdkTable._fields_ + (
        ("evaluate_mixture_phase_kij", ctypes.c_void_p),
    )


class ElectrolyteNativeSdkTable(ctypes.Structure):
    _fields_ = NativeSdkTable._fields_ + (
        ("component_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("component_charges", ctypes.POINTER(ctypes.c_int32)),
        ("evaluate_electrolyte_phase", ctypes.c_void_p),
    )


class VolumeBoundsNativeSdkTable(ctypes.Structure):
    _fields_ = ElectrolyteNativeSdkTable._fields_ + (
        ("evaluate_molar_volume_bounds", ctypes.c_void_p),
    )


class PackingFractionNativeSdkTable(ctypes.Structure):
    _fields_ = VolumeBoundsNativeSdkTable._fields_ + (
        ("evaluate_packing_fraction", ctypes.c_void_p),
    )


class ApplicabilityNativeSdkTable(ctypes.Structure):
    _fields_ = PackingFractionNativeSdkTable._fields_ + (
        ("source_temperature_min_k", ctypes.c_double),
        ("source_temperature_max_k", ctypes.c_double),
        ("total_ion_mole_fraction_max", ctypes.c_double),
    )


class IonSolvationBornResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("solvation_gibbs_j_per_mol", ctypes.c_double),
        ("derivative_j_per_mol_per_angstrom", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class ActiveBornNativeSdkTable(ctypes.Structure):
    _fields_ = ApplicabilityNativeSdkTable._fields_ + (
        ("ion_solvation_born_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_born", ctypes.c_void_p),
    )


class AqueousMiacKijResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("log_mean_ionic_activity_coefficient_molality", ctypes.c_double),
        ("derivative", ctypes.c_double * 3),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("reference_derivative_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AqueousMiacKijNativeSdkTable(ctypes.Structure):
    _fields_ = ActiveBornNativeSdkTable._fields_ + (
        ("aqueous_miac_kij_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij", ctypes.c_void_p),
    )


class NeutralReferenceResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("component_count", ctypes.c_size_t),
        ("neutral_basis_row_count", ctypes.c_size_t),
        ("neutral_basis_capacity", ctypes.c_size_t),
        ("contraction_capacity", ctypes.c_size_t),
        ("reference_composition_capacity", ctypes.c_size_t),
        ("neutral_basis", ctypes.POINTER(ctypes.c_double)),
        ("log_fugacity_contractions", ctypes.POINTER(ctypes.c_double)),
        ("reference_composition", ctypes.POINTER(ctypes.c_double)),
        ("derivative_availability", ctypes.c_uint32),
        ("temperature_k", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("solvent_molar_mass_kg_per_mol", ctypes.c_double),
        ("reference_amount_mol", ctypes.c_double),
        ("reference_number_density_mol_per_m3", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class NeutralReferenceNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacKijNativeSdkTable._fields_ + (
        ("neutral_reference_basis_row_count", ctypes.c_size_t),
        ("neutral_reference_result_size", ctypes.c_size_t),
        ("evaluate_neutral_reference", ctypes.c_void_p),
    )


class AqueousMiacSolvationFactorResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("log_mean_ionic_activity_coefficient_molality", ctypes.c_double),
        ("derivative", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("reference_derivative_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AqueousMiacSolvationFactorNativeSdkTable(ctypes.Structure):
    _fields_ = NeutralReferenceNativeSdkTable._fields_ + (
        ("aqueous_miac_solvation_factor_result_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_solvation_factor", ctypes.c_void_p),
    )


class AqueousMiacBatchNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacSolvationFactorNativeSdkTable._fields_ + (
        ("evaluate_aqueous_miac_kij_batch", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch", ctypes.c_void_p),
    )


class NativeEvaluationBudget(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("work_budget", ctypes.c_uint64),
        ("deadline_ns", ctypes.c_uint64),
    )


class AqueousMiacBoundedBatchNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacBatchNativeSdkTable._fields_ + (
        ("evaluation_budget_size", ctypes.c_size_t),
        ("evaluate_aqueous_miac_kij_batch_bounded", ctypes.c_void_p),
        ("evaluate_aqueous_miac_solvation_factor_batch_bounded", ctypes.c_void_p),
    )


class NativeCapabilityCoordinate(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("kind", ctypes.c_uint32),
        ("component_index", ctypes.c_int32),
        ("pair_component_index_a", ctypes.c_int32),
        ("pair_component_index_b", ctypes.c_int32),
        ("unit", ctypes.c_char_p),
    )


class NativeCapabilityDescriptor(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("schema_version", ctypes.c_uint32),
        ("capability", ctypes.c_uint32),
        ("parameter_family", ctypes.c_uint32),
        ("parameter_identity", ctypes.c_uint32),
        ("observation_contract", ctypes.c_uint32),
        ("model_domain", ctypes.c_uint32),
        ("tensor_layout", ctypes.c_uint32),
        ("derivative_order", ctypes.c_uint32),
        ("maturity", ctypes.c_uint32),
        ("authority_effect", ctypes.c_uint32),
        ("unsupported_status", ctypes.c_int32),
        ("domain_status", ctypes.c_int32),
        ("state_coordinate_count", ctypes.c_size_t),
        ("active_parameter_count", ctypes.c_size_t),
        ("coordinate_count", ctypes.c_size_t),
        ("coordinates", ctypes.POINTER(NativeCapabilityCoordinate)),
        ("component_count", ctypes.c_size_t),
        ("component_ids", ctypes.POINTER(ctypes.c_char_p)),
        ("temperature_min_k", ctypes.c_double),
        ("temperature_max_k", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("helmholtz_basis_id", ctypes.c_char * 72),
    )


class CapabilityNativeSdkTable(ctypes.Structure):
    _fields_ = AqueousMiacBoundedBatchNativeSdkTable._fields_ + (
        ("capability_count", ctypes.c_size_t),
        ("capabilities", ctypes.POINTER(NativeCapabilityDescriptor)),
    )


class LijNativeSdkTable(ctypes.Structure):
    _fields_ = CapabilityNativeSdkTable._fields_ + (
        ("evaluate_mixture_phase_lij", ctypes.c_void_p),
    )


class PureParameterNativeSdkTable(ctypes.Structure):
    _fields_ = LijNativeSdkTable._fields_ + (
        ("evaluate_pure_phase_parameter", ctypes.c_void_p),
    )


class IonFractionSuppressionResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("bulk_relative_permittivity", ctypes.c_double),
        ("derivative", ctypes.c_double),
        ("salt_free_relative_permittivity", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class IonFractionSuppressionNativeSdkTable(ctypes.Structure):
    _fields_ = PureParameterNativeSdkTable._fields_ + (
        ("ion_fraction_suppression_result_size", ctypes.c_size_t),
        ("evaluate_ion_fraction_suppression", ctypes.c_void_p),
    )


class IonSolvationKijResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("solvation_gibbs_j_per_mol", ctypes.c_double),
        ("derivative", ctypes.c_double * 3),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class IonSolvationKijNativeSdkTable(ctypes.Structure):
    _fields_ = IonFractionSuppressionNativeSdkTable._fields_ + (
        ("ion_solvation_kij_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_kij", ctypes.c_void_p),
    )


class IonSolvationIonicPermittivityResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("solvation_gibbs_j_per_mol", ctypes.c_double),
        ("derivative_j_per_mol", ctypes.c_double),
        ("reference_molality_mol_per_kg", ctypes.c_double),
        ("reference_convergence_error", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class IonSolvationIonicPermittivityNativeSdkTable(ctypes.Structure):
    _fields_ = IonSolvationKijNativeSdkTable._fields_ + (
        ("ion_solvation_ionic_permittivity_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_ionic_permittivity", ctypes.c_void_p),
    )


class IonSolvationSolventPermittivityResult(ctypes.Structure):
    _fields_ = IonSolvationIonicPermittivityResult._fields_


class IonSolvationSolventPermittivityNativeSdkTable(ctypes.Structure):
    _fields_ = IonSolvationIonicPermittivityNativeSdkTable._fields_ + (
        ("ion_solvation_solvent_permittivity_result_size", ctypes.c_size_t),
        ("evaluate_ion_solvation_solvent_permittivity", ctypes.c_void_p),
    )


class AssociatingParameterizedPhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("gradient", ctypes.c_double * 7),
        ("hessian", ctypes.c_double * 49),
        ("pressure_pa", ctypes.c_double),
        ("chemical_potential_over_rt", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("topology_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )


class AssociatingParameterizedNativeSdkTable(ctypes.Structure):
    _fields_ = IonSolvationSolventPermittivityNativeSdkTable._fields_ + (
        ("associating_parameterized_result_size", ctypes.c_size_t),
        (
            "evaluate_diameter_basis_associating_pure_phase_parameters",
            ctypes.c_void_p,
        ),
    )


class NativeSdkPrefix(ctypes.Structure):
    _fields_ = NativeSdkTable._fields_[:5]


FIXED_MIXTURE_PREFIX_SIZE = (
    FixedMixtureNativeSdkTable.evaluate_mixture_phase.offset + ctypes.sizeof(ctypes.c_void_p)
)
ACTIVE_KIJ_PREFIX_SIZE = NativeSdkTable.evaluate_mixture_phase_kij.offset + ctypes.sizeof(
    ctypes.c_void_p
)
ELECTROLYTE_PREFIX_SIZE = (
    ElectrolyteNativeSdkTable.evaluate_electrolyte_phase.offset + ctypes.sizeof(ctypes.c_void_p)
)
VOLUME_BOUNDS_PREFIX_SIZE = (
    VolumeBoundsNativeSdkTable.evaluate_molar_volume_bounds.offset + ctypes.sizeof(ctypes.c_void_p)
)
PACKING_FRACTION_PREFIX_SIZE = (
    PackingFractionNativeSdkTable.evaluate_packing_fraction.offset + ctypes.sizeof(ctypes.c_void_p)
)
APPLICABILITY_PREFIX_SIZE = (
    ApplicabilityNativeSdkTable.total_ion_mole_fraction_max.offset + ctypes.sizeof(ctypes.c_double)
)
ACTIVE_BORN_PREFIX_SIZE = (
    ActiveBornNativeSdkTable.evaluate_ion_solvation_born.offset + ctypes.sizeof(ctypes.c_void_p)
)
AQUEOUS_MIAC_KIJ_PREFIX_SIZE = (
    AqueousMiacKijNativeSdkTable.evaluate_aqueous_miac_kij.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
NEUTRAL_REFERENCE_PREFIX_SIZE = (
    NeutralReferenceNativeSdkTable.evaluate_neutral_reference.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
AQUEOUS_MIAC_SOLVATION_FACTOR_PREFIX_SIZE = (
    AqueousMiacSolvationFactorNativeSdkTable.evaluate_aqueous_miac_solvation_factor.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
AQUEOUS_MIAC_BATCH_PREFIX_SIZE = (
    AqueousMiacBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
AQUEOUS_MIAC_BOUNDED_BATCH_PREFIX_SIZE = (
    AqueousMiacBoundedBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch_bounded.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
CAPABILITY_PREFIX_SIZE = (
    CapabilityNativeSdkTable.capabilities.offset + ctypes.sizeof(ctypes.c_void_p)
)
ACTIVE_LIJ_PREFIX_SIZE = (
    LijNativeSdkTable.evaluate_mixture_phase_lij.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
ACTIVE_PURE_PARAMETER_PREFIX_SIZE = (
    PureParameterNativeSdkTable.evaluate_pure_phase_parameter.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
ION_FRACTION_SUPPRESSION_PREFIX_SIZE = (
    IonFractionSuppressionNativeSdkTable.evaluate_ion_fraction_suppression.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
ION_SOLVATION_KIJ_PREFIX_SIZE = (
    IonSolvationKijNativeSdkTable.evaluate_ion_solvation_kij.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
ION_SOLVATION_IONIC_PERMITTIVITY_PREFIX_SIZE = (
    IonSolvationIonicPermittivityNativeSdkTable.evaluate_ion_solvation_ionic_permittivity.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
ION_SOLVATION_SOLVENT_PERMITTIVITY_PREFIX_SIZE = (
    IonSolvationSolventPermittivityNativeSdkTable
    .evaluate_ion_solvation_solvent_permittivity.offset
    + ctypes.sizeof(ctypes.c_void_p)
)
ASSOCIATING_PARAMETERIZED_PREFIX_SIZE = (
    AssociatingParameterizedNativeSdkTable
    .evaluate_diameter_basis_associating_pure_phase_parameters.offset
    + ctypes.sizeof(ctypes.c_void_p)
)


EvaluatePurePhase = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(PhaseBlockResult),
)

EvaluateIonSolvationKij = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(IonSolvationKijResult),
)

EvaluatePurePhaseParameters = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ParameterizedPhaseBlockResult),
)

EvaluateAssociatingPurePhaseParameters = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(AssociatingParameterizedPhaseBlockResult),
)

EvaluateMixturePhase = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(MixturePhaseBlockResult),
)

EvaluateMixturePhaseKij = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(MixturePhaseBlockResult),
)
EvaluateMixturePhaseLij = EvaluateMixturePhaseKij
EvaluatePurePhaseParameter = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_uint32,
    ctypes.c_double,
    ctypes.POINTER(MixturePhaseBlockResult),
)

EvaluateMolarVolumeBounds = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)

EvaluatePackingFraction = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
)

EvaluateIonSolvationBorn = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(IonSolvationBornResult),
)

EvaluateIonSolvationIonicPermittivity = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(IonSolvationIonicPermittivityResult),
)

EvaluateIonSolvationSolventPermittivity = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(IonSolvationSolventPermittivityResult),
)

EvaluateAqueousMiacKij = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(AqueousMiacKijResult),
)

EvaluateNeutralReference = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(NeutralReferenceResult),
)

EvaluateAqueousMiacSolvationFactor = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(AqueousMiacSolvationFactorResult),
)

EvaluateIonFractionSuppression = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(IonFractionSuppressionResult),
)

EvaluateAqueousMiacKijBatch = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(AqueousMiacKijResult),
    ctypes.c_size_t,
)

EvaluateAqueousMiacSolvationFactorBatch = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(AqueousMiacSolvationFactorResult),
    ctypes.c_size_t,
)

EvaluateAqueousMiacKijBatchBounded = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.POINTER(AqueousMiacKijResult),
    ctypes.c_size_t,
    ctypes.POINTER(NativeEvaluationBudget),
)

EvaluateAqueousMiacSolvationFactorBatchBounded = ctypes.CFUNCTYPE(
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t,
    ctypes.c_double,
    ctypes.POINTER(AqueousMiacSolvationFactorResult),
    ctypes.c_size_t,
    ctypes.POINTER(NativeEvaluationBudget),
)


def _model(catalog: str, component: str) -> Mixture:
    return Mixture(_ParameterBundle.from_catalog(catalog, version=1).select((component,)))


def _mixture_model(component_ids: tuple[str, str] = ("methane", "ethane")) -> Mixture:
    parameters = _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    return Mixture(parameters.select(component_ids))


def _electrolyte_model(
    component_ids: tuple[str, str, str] = (
        "water",
        "sodium-cation",
        "chloride-anion",
    ),
) -> Mixture:
    parameters = _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    return Mixture(parameters.select(component_ids))


def _user_provided_copy(
    bundle: _ParameterBundle, *, unknown_applicability: bool = False
) -> _ParameterBundle:
    records = bundle.records
    return _ParameterBundle.from_records(
        bundle_id=f"user-{bundle.bundle_id}",
        bundle_version=bundle.bundle_version,
        purpose="user-provided",
        sources=bundle.sources,
        domains=(
            tuple(ValidityDomain(domain.domain_id, "unknown") for domain in bundle.domains)
            if unknown_applicability
            else bundle.domains
        ),
        components=bundle.components,
        singles=(record for record in records if isinstance(record, SingleParameterRecord)),
        pairs=(record for record in records if isinstance(record, PairParameterRecord)),
        sites=(record for record in records if isinstance(record, SiteRecord)),
        associations=(
            record for record in records if isinstance(record, AssociationParameterRecord)
        ),
        correlations=(
            record
            for record in records
            if isinstance(
                record,
                (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation),
            )
        ),
        models=(record for record in records if isinstance(record, ModelParameterRecord)),
    )


def _figiel_trial_model(
    component_ids: tuple[str, str, str],
    *,
    purpose: str = "user-provided",
    spoof_source_metadata: bool = False,
) -> Mixture:
    bundle = _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    records = tuple(
        replace(record, value=1.4)
        if record.record_id == "water-solvation-factor"
        else record
        for record in bundle.records
    )
    trial = _ParameterBundle.from_records(
        bundle_id="figiel-staged-trial",
        bundle_version=1,
        purpose=purpose,
        sources=(
            replace(source, citation="fabricated citation")
            if spoof_source_metadata and source.source_id == "figiel-yu-held-2025"
            else source
            for source in bundle.sources
        ),
        domains=bundle.domains,
        components=bundle.components,
        singles=(record for record in records if isinstance(record, SingleParameterRecord)),
        pairs=(record for record in records if isinstance(record, PairParameterRecord)),
        sites=(record for record in records if isinstance(record, SiteRecord)),
        associations=(
            record for record in records if isinstance(record, AssociationParameterRecord)
        ),
        correlations=(
            record
            for record in records
            if isinstance(
                record,
                (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation),
            )
        ),
        models=(record for record in records if isinstance(record, ModelParameterRecord)),
    )
    return Mixture(trial.select(component_ids))


def _held_water_ionization_model(
    component_ids: tuple[str, str, str] = (
        "water",
        "hydronium-cation",
        "hydroxide-anion",
    ),
) -> Mixture:
    parameters = _ParameterBundle.from_catalog("held-2008-water-self-ionization", version=1)
    return Mixture(parameters.select(component_ids))


def _khudaida_model(component_ids: tuple[str, ...]) -> Mixture:
    parameters = _ParameterBundle.from_catalog("khudaida-2026-figure-2-electrolyte-lle", version=1)
    return Mixture(parameters.select(component_ids))


def _ascani_model(component_ids: tuple[str, ...]) -> Mixture:
    parameters = _ParameterBundle.from_catalog("ascani-2022-case-study-2-electrolyte-lle", version=1)
    return Mixture(parameters.select(component_ids))


def _table(model: Mixture) -> tuple[object, NativeSdkTable]:
    capsule = native_sdk(model)
    get_name = ctypes.pythonapi.PyCapsule_GetName
    get_name.argtypes = (ctypes.py_object,)
    get_name.restype = ctypes.c_char_p
    assert get_name(capsule) == CAPSULE_NAME

    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    fixed_table = ctypes.cast(address, ctypes.POINTER(FixedMixtureNativeSdkTable)).contents
    assert fixed_table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
    return capsule, ctypes.cast(address, ctypes.POINTER(NativeSdkTable)).contents


def _capability_table(
    model: Mixture,
) -> tuple[object, CapabilityNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(CapabilityNativeSdkTable)).contents
    assert table.table_size >= CAPABILITY_PREFIX_SIZE
    return capsule, table


def _ion_fraction_suppression_table(
    model: Mixture,
) -> tuple[object, IonFractionSuppressionNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(IonFractionSuppressionNativeSdkTable)
    ).contents
    assert table.table_size >= ctypes.sizeof(
        IonFractionSuppressionNativeSdkTable
    )
    return capsule, table


def _ion_solvation_kij_table(
    model: Mixture,
) -> tuple[object, IonSolvationKijNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(IonSolvationKijNativeSdkTable)
    ).contents
    assert table.table_size >= ION_SOLVATION_KIJ_PREFIX_SIZE
    return capsule, table


def _ion_solvation_ionic_permittivity_table(
    model: Mixture,
) -> tuple[object, IonSolvationIonicPermittivityNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(IonSolvationIonicPermittivityNativeSdkTable)
    ).contents
    assert table.table_size >= ION_SOLVATION_IONIC_PERMITTIVITY_PREFIX_SIZE
    return capsule, table


def _ion_solvation_solvent_permittivity_table(
    model: Mixture,
) -> tuple[object, IonSolvationSolventPermittivityNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address,
        ctypes.POINTER(IonSolvationSolventPermittivityNativeSdkTable),
    ).contents
    assert table.table_size >= ION_SOLVATION_SOLVENT_PERMITTIVITY_PREFIX_SIZE
    return capsule, table


def _lij_table(model: Mixture) -> tuple[object, LijNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(LijNativeSdkTable)).contents
    assert table.table_size >= ACTIVE_LIJ_PREFIX_SIZE
    return capsule, table


def _pure_parameter_table(
    model: Mixture,
) -> tuple[object, PureParameterNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(PureParameterNativeSdkTable)
    ).contents
    assert table.table_size >= ACTIVE_PURE_PARAMETER_PREFIX_SIZE
    return capsule, table


def _associating_parameterized_table(
    model: Mixture,
) -> tuple[object, AssociatingParameterizedNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(AssociatingParameterizedNativeSdkTable)
    ).contents
    assert table.table_size >= ASSOCIATING_PARAMETERIZED_PREFIX_SIZE
    return capsule, table


def _electrolyte_table(
    model: Mixture,
) -> tuple[object, ElectrolyteNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(ElectrolyteNativeSdkTable)).contents
    assert table.table_size >= ELECTROLYTE_PREFIX_SIZE
    return capsule, table


def _volume_bounds_table(
    model: Mixture,
) -> tuple[object, VolumeBoundsNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(VolumeBoundsNativeSdkTable)).contents
    assert table.table_size >= VOLUME_BOUNDS_PREFIX_SIZE
    return capsule, table


def _packing_fraction_table(
    model: Mixture,
) -> tuple[object, PackingFractionNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(PackingFractionNativeSdkTable)).contents
    assert table.table_size >= PACKING_FRACTION_PREFIX_SIZE
    return capsule, table


def _applicability_table(
    model: Mixture,
) -> tuple[object, ApplicabilityNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(ApplicabilityNativeSdkTable)).contents
    assert table.table_size >= APPLICABILITY_PREFIX_SIZE
    return capsule, table


def _active_born_table(model: Mixture) -> tuple[object, ActiveBornNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(ActiveBornNativeSdkTable)).contents
    assert table.table_size >= ACTIVE_BORN_PREFIX_SIZE
    return capsule, table


def _aqueous_miac_kij_table(
    model: Mixture,
) -> tuple[object, AqueousMiacKijNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(AqueousMiacKijNativeSdkTable)).contents
    assert table.table_size >= AQUEOUS_MIAC_KIJ_PREFIX_SIZE
    return capsule, table


def _neutral_reference_table(
    model: Mixture,
) -> tuple[object, NeutralReferenceNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(address, ctypes.POINTER(NeutralReferenceNativeSdkTable)).contents
    assert table.table_size >= NEUTRAL_REFERENCE_PREFIX_SIZE
    return capsule, table


def _aqueous_miac_solvation_factor_table(
    model: Mixture,
) -> tuple[object, AqueousMiacSolvationFactorNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(AqueousMiacSolvationFactorNativeSdkTable)
    ).contents
    assert table.table_size >= AQUEOUS_MIAC_SOLVATION_FACTOR_PREFIX_SIZE
    return capsule, table


def _aqueous_miac_batch_table(
    model: Mixture,
) -> tuple[object, AqueousMiacBatchNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(AqueousMiacBatchNativeSdkTable)
    ).contents
    assert table.table_size >= AQUEOUS_MIAC_BATCH_PREFIX_SIZE
    return capsule, table


def _aqueous_miac_bounded_batch_table(
    model: Mixture,
) -> tuple[object, AqueousMiacBoundedBatchNativeSdkTable]:
    capsule = native_sdk(model)
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    assert address
    table = ctypes.cast(
        address, ctypes.POINTER(AqueousMiacBoundedBatchNativeSdkTable)
    ).contents
    assert table.table_size >= AQUEOUS_MIAC_BOUNDED_BATCH_PREFIX_SIZE
    return capsule, table


def _evaluate_aqueous_miac_solvation_factor(
    table: AqueousMiacSolvationFactorNativeSdkTable,
    fingerprint: str,
    molality: float,
    solvation_factor: float,
) -> AqueousMiacSolvationFactorResult:
    result = AqueousMiacSolvationFactorResult(
        struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
    )
    status = EvaluateAqueousMiacSolvationFactor(
        table.evaluate_aqueous_miac_solvation_factor
    )(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality,
        solvation_factor,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_aqueous_miac_solvation_factor_batch(
    table: AqueousMiacBatchNativeSdkTable,
    fingerprint: str,
    molalities: tuple[float, ...],
    solvation_factor: float,
) -> tuple[AqueousMiacSolvationFactorResult, ...]:
    molality_values = (ctypes.c_double * len(molalities))(*molalities)
    results = (AqueousMiacSolvationFactorResult * len(molalities))(
        *(
            AqueousMiacSolvationFactorResult(
                struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
            )
            for _ in molalities
        )
    )
    status = EvaluateAqueousMiacSolvationFactorBatch(
        table.evaluate_aqueous_miac_solvation_factor_batch
    )(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        solvation_factor,
        results,
        len(molalities),
    )
    assert all(result.status == status for result in results)
    return tuple(results)


def _evaluate_neutral_reference(
    table: NeutralReferenceNativeSdkTable,
    fingerprint: str,
    *,
    temperature_k: float = 298.15,
    pressure_pa: float = 100000.0,
) -> tuple[NeutralReferenceResult, tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    component_count = table.component_count
    row_count = table.neutral_reference_basis_row_count
    basis = (ctypes.c_double * (component_count * row_count))()
    contractions = (ctypes.c_double * row_count)()
    reference_composition = (ctypes.c_double * component_count)()
    result = NeutralReferenceResult(
        struct_size=ctypes.sizeof(NeutralReferenceResult),
        component_count=component_count,
        neutral_basis_row_count=row_count,
        neutral_basis_capacity=component_count * row_count,
        contraction_capacity=row_count,
        reference_composition_capacity=component_count,
        neutral_basis=basis,
        log_fugacity_contractions=contractions,
        reference_composition=reference_composition,
    )
    status = EvaluateNeutralReference(table.evaluate_neutral_reference)(
        table.model_context,
        fingerprint.encode(),
        temperature_k,
        pressure_pa,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(basis), tuple(contractions), tuple(reference_composition)


def _evaluate_aqueous_miac_kij(
    table: AqueousMiacKijNativeSdkTable,
    fingerprint: str,
    molality: float,
    k_ij: tuple[float, float, float],
) -> AqueousMiacKijResult:
    coordinates = (ctypes.c_double * 3)(*k_ij)
    result = AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult))
    status = EvaluateAqueousMiacKij(table.evaluate_aqueous_miac_kij)(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality,
        coordinates,
        3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_aqueous_miac_kij_batch(
    table: AqueousMiacBatchNativeSdkTable,
    fingerprint: str,
    molalities: tuple[float, ...],
    k_ij: tuple[float, float, float],
) -> tuple[AqueousMiacKijResult, ...]:
    molality_values = (ctypes.c_double * len(molalities))(*molalities)
    coordinates = (ctypes.c_double * 3)(*k_ij)
    results = (AqueousMiacKijResult * len(molalities))(
        *(
            AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult))
            for _ in molalities
        )
    )
    status = EvaluateAqueousMiacKijBatch(table.evaluate_aqueous_miac_kij_batch)(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        coordinates,
        3,
        results,
        len(molalities),
    )
    assert all(result.status == status for result in results)
    return tuple(results)


def _evaluate_active_born(
    table: ActiveBornNativeSdkTable,
    born_diameter_angstrom: float,
) -> IonSolvationBornResult:
    result = IonSolvationBornResult(struct_size=ctypes.sizeof(IonSolvationBornResult))
    status = EvaluateIonSolvationBorn(table.evaluate_ion_solvation_born)(
        table.model_context,
        298.15,
        100000.0,
        born_diameter_angstrom,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_ionic_permittivity(
    table: IonSolvationIonicPermittivityNativeSdkTable,
    fingerprint: str,
    ionic_region_relative_permittivity: float,
) -> IonSolvationIonicPermittivityResult:
    result = IonSolvationIonicPermittivityResult(
        struct_size=ctypes.sizeof(IonSolvationIonicPermittivityResult)
    )
    status = EvaluateIonSolvationIonicPermittivity(
        table.evaluate_ion_solvation_ionic_permittivity
    )(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        ionic_region_relative_permittivity,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_solvent_permittivity(
    table: IonSolvationSolventPermittivityNativeSdkTable,
    fingerprint: str,
    solvent_relative_permittivity: float,
) -> IonSolvationSolventPermittivityResult:
    result = IonSolvationSolventPermittivityResult(
        struct_size=ctypes.sizeof(IonSolvationSolventPermittivityResult)
    )
    status = EvaluateIonSolvationSolventPermittivity(
        table.evaluate_ion_solvation_solvent_permittivity
    )(
        table.model_context,
        fingerprint.encode(),
        298.15,
        100000.0,
        solvent_relative_permittivity,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate(
    table: NativeSdkTable,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
) -> PhaseBlockResult:
    result = PhaseBlockResult(struct_size=ctypes.sizeof(PhaseBlockResult))
    status = EvaluatePurePhase(table.evaluate_pure_phase)(
        table.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_parameters(
    table: NativeSdkTable,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
    segment_count: float,
    segment_diameter_angstrom: float,
    dispersion_energy_over_k: float,
) -> ParameterizedPhaseBlockResult:
    assert table.table_size >= ctypes.sizeof(NativeSdkTable)
    assert table.evaluate_pure_phase_parameters
    result = ParameterizedPhaseBlockResult(struct_size=ctypes.sizeof(ParameterizedPhaseBlockResult))
    status = EvaluatePurePhaseParameters(table.evaluate_pure_phase_parameters)(
        table.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _evaluate_associating_parameters(
    table: AssociatingParameterizedNativeSdkTable,
    temperature_k: float,
    values: tuple[float, float, float, float, float, float, float],
) -> AssociatingParameterizedPhaseBlockResult:
    result = AssociatingParameterizedPhaseBlockResult(
        struct_size=ctypes.sizeof(AssociatingParameterizedPhaseBlockResult)
    )
    status = EvaluateAssociatingPurePhaseParameters(
        table.evaluate_diameter_basis_associating_pure_phase_parameters
    )(
        table.model_context,
        temperature_k,
        *values,
        ctypes.byref(result),
    )
    assert status == result.status
    return result


def _mixture_result(coordinate_count: int) -> tuple[MixturePhaseBlockResult, object, object]:
    gradient = (ctypes.c_double * coordinate_count)()
    hessian = (ctypes.c_double * (coordinate_count * coordinate_count))()
    result = MixturePhaseBlockResult(
        struct_size=ctypes.sizeof(MixturePhaseBlockResult),
        coordinate_count=coordinate_count,
        gradient_capacity=coordinate_count,
        hessian_capacity=coordinate_count * coordinate_count,
        gradient=gradient,
        hessian=hessian,
    )
    return result, gradient, hessian


def _evaluate_mixture(
    table: NativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    assert table.table_size >= FIXED_MIXTURE_PREFIX_SIZE
    assert table.evaluate_mixture_phase
    coordinate_count = len(amounts_mol) + 1
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhase(table.evaluate_mixture_phase)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_mixture_kij(
    table: NativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    k_ij: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    assert table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
    assert table.evaluate_mixture_phase_kij
    coordinate_count = len(amounts_mol) + 2
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhaseKij(table.evaluate_mixture_phase_kij)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        k_ij,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_mixture_lij(
    table: LijNativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
    l_ij: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    assert table.evaluate_mixture_phase_lij
    coordinate_count = len(amounts_mol) + 2
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhaseLij(table.evaluate_mixture_phase_lij)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        l_ij,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_pure_parameter(
    table: PureParameterNativeSdkTable,
    temperature_k: float,
    amount_mol: float,
    volume_m3: float,
    parameter_family: int,
    parameter_value: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    result, gradient, hessian = _mixture_result(3)
    status = EvaluatePurePhaseParameter(table.evaluate_pure_phase_parameter)(
        table.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        parameter_family,
        parameter_value,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


def _evaluate_electrolyte(
    table: ElectrolyteNativeSdkTable,
    temperature_k: float,
    amounts_mol: tuple[float, ...],
    volume_m3: float,
) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
    coordinate_count = len(amounts_mol) + 1
    result, gradient, hessian = _mixture_result(coordinate_count)
    amounts = (ctypes.c_double * len(amounts_mol))(*amounts_mol)
    status = EvaluateMixturePhase(table.evaluate_electrolyte_phase)(
        table.model_context,
        temperature_k,
        amounts,
        len(amounts_mol),
        volume_m3,
        ctypes.byref(result),
    )
    assert status == result.status
    return result, tuple(gradient), tuple(hessian)


@pytest.mark.parametrize(
    "catalog,component,temperature_k,density_mol_m3",
    (
        ("gross-2001-methane-ethane", "methane", 190.0, 1000.0),
        ("gross-2001-methane-ethane", "ethane", 250.0, 1000.0),
        ("gross-2001-propane", "propane", 300.0, 1000.0),
    ),
)
def test_native_sdk_phase_block_matches_public_thermodynamics(
    catalog: str,
    component: str,
    temperature_k: float,
    density_mol_m3: float,
) -> None:
    model = _model(catalog, component)
    capsule, table = _table(model)
    assert capsule is not None
    assert table.abi_version == ABI_VERSION
    assert table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
    assert table.result_size == ctypes.sizeof(PhaseBlockResult)
    assert table.parameterized_result_size == ctypes.sizeof(ParameterizedPhaseBlockResult)
    assert table.model_context
    assert table.evaluate_pure_phase
    assert table.evaluate_pure_phase_parameters

    amount_mol = 1.0
    volume_m3 = amount_mol / density_mol_m3
    block = _evaluate(table, temperature_k, amount_mol, volume_m3)
    assert block.status == STATUS_OK
    assert block.error == b""
    assert block.parameter_fingerprint.decode() == model.parameter_fingerprint

    state = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density_mol_m3 * unit_registry.mole / unit_registry.meter**3,
        x=(1.0,),
    )
    expected_helmholtz = amount_mol * (math.log(density_mol_m3) - 1.0 + state._residual_helmholtz)
    expected_chemical_potential = (
        math.log(density_mol_m3) + state.residual_chemical_potential_over_rt[0]
    )
    assert block.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_helmholtz, rel=2.0e-12, abs=2.0e-12
    )
    assert block.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=2.0e-11, abs=1.0e-5
    )
    assert block.chemical_potential_over_rt == pytest.approx(
        expected_chemical_potential, rel=2.0e-11, abs=2.0e-11
    )


def test_native_sdk_derivative_tensors_match_directional_oracles() -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    capsule, table = _table(model)
    del model
    gc.collect()
    assert capsule is not None
    temperature_k = 190.0
    point = (1.2, 1.2e-3)
    direction = (0.05, -3.0e-5)

    def shifted(scale: float) -> PhaseBlockResult:
        result = _evaluate(
            table,
            temperature_k,
            point[0] + scale * direction[0],
            point[1] + scale * direction[1],
        )
        assert result.status == STATUS_OK
        return result

    center = shifted(0.0)
    step = 2.0e-4
    lower = shifted(-step)
    upper = shifted(step)

    gradient_contraction = sum(center.gradient[index] * direction[index] for index in range(2))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_contraction == pytest.approx(value_difference, rel=2.0e-8)

    for row in range(2):
        hessian_vector = sum(
            center.hessian[row * 2 + column] * direction[column] for column in range(2)
        )
        gradient_difference = (upper.gradient[row] - lower.gradient[row]) / (2.0 * step)
        assert hessian_vector == pytest.approx(gradient_difference, rel=2.0e-7)

    def hessian_contraction(result: PhaseBlockResult) -> float:
        return sum(
            direction[row] * result.hessian[row * 2 + column] * direction[column]
            for row in range(2)
            for column in range(2)
        )

    third_contraction = sum(
        direction[i] * direction[j] * direction[k] * center.third[(i * 2 + j) * 2 + k]
        for i in range(2)
        for j in range(2)
        for k in range(2)
    )
    hessian_difference = (hessian_contraction(upper) - hessian_contraction(lower)) / (2.0 * step)
    assert third_contraction == pytest.approx(hessian_difference, rel=2.0e-6)


def test_native_sdk_addition_preserves_the_original_table_prefix() -> None:
    capsule = native_sdk(_model("gross-2001-methane-ethane", "methane"))
    get_pointer = ctypes.pythonapi.PyCapsule_GetPointer
    get_pointer.argtypes = (ctypes.py_object, ctypes.c_char_p)
    get_pointer.restype = ctypes.c_void_p
    address = get_pointer(capsule, CAPSULE_NAME)
    prefix = ctypes.cast(address, ctypes.POINTER(NativeSdkPrefix)).contents
    assert prefix.abi_version == ABI_VERSION
    assert prefix.table_size >= ctypes.sizeof(NativeSdkPrefix)
    assert prefix.result_size == ctypes.sizeof(PhaseBlockResult)
    assert prefix.model_context
    assert prefix.evaluate_pure_phase
    for field_name, _ in AcceptedNativeSdkTable._fields_:
        assert (
            getattr(NativeSdkTable, field_name).offset
            == getattr(AcceptedNativeSdkTable, field_name).offset
        )
    for field_name, _ in FixedMixtureNativeSdkTable._fields_:
        assert (
            getattr(NativeSdkTable, field_name).offset
            == getattr(FixedMixtureNativeSdkTable, field_name).offset
        )
    assert FIXED_MIXTURE_PREFIX_SIZE == ctypes.sizeof(FixedMixtureNativeSdkTable)
    assert NativeSdkTable.evaluate_mixture_phase_kij.offset == FIXED_MIXTURE_PREFIX_SIZE
    for field_name, _ in ElectrolyteNativeSdkTable._fields_:
        assert (
            getattr(VolumeBoundsNativeSdkTable, field_name).offset
            == getattr(ElectrolyteNativeSdkTable, field_name).offset
        )
    assert VolumeBoundsNativeSdkTable.evaluate_molar_volume_bounds.offset == ELECTROLYTE_PREFIX_SIZE
    for field_name, _ in VolumeBoundsNativeSdkTable._fields_:
        assert (
            getattr(PackingFractionNativeSdkTable, field_name).offset
            == getattr(VolumeBoundsNativeSdkTable, field_name).offset
        )
    assert VOLUME_BOUNDS_PREFIX_SIZE == ctypes.sizeof(VolumeBoundsNativeSdkTable)
    assert (
        PackingFractionNativeSdkTable.evaluate_packing_fraction.offset == VOLUME_BOUNDS_PREFIX_SIZE
    )
    for field_name, _ in PackingFractionNativeSdkTable._fields_:
        assert (
            getattr(ApplicabilityNativeSdkTable, field_name).offset
            == getattr(PackingFractionNativeSdkTable, field_name).offset
        )
    assert PACKING_FRACTION_PREFIX_SIZE == ctypes.sizeof(PackingFractionNativeSdkTable) == 128
    assert ApplicabilityNativeSdkTable.source_temperature_min_k.offset == 128
    assert ApplicabilityNativeSdkTable.source_temperature_max_k.offset == 136
    assert ApplicabilityNativeSdkTable.total_ion_mole_fraction_max.offset == 144
    assert APPLICABILITY_PREFIX_SIZE == ctypes.sizeof(ApplicabilityNativeSdkTable) == 152
    for field_name, _ in AqueousMiacBoundedBatchNativeSdkTable._fields_:
        assert (
            getattr(CapabilityNativeSdkTable, field_name).offset
            == getattr(AqueousMiacBoundedBatchNativeSdkTable, field_name).offset
        )
    assert AQUEOUS_MIAC_BOUNDED_BATCH_PREFIX_SIZE == 264
    assert CapabilityNativeSdkTable.capability_count.offset == 264
    assert CapabilityNativeSdkTable.capabilities.offset == 272
    assert ctypes.sizeof(CapabilityNativeSdkTable) == 280
    for field_name, _ in CapabilityNativeSdkTable._fields_:
        assert (
            getattr(LijNativeSdkTable, field_name).offset
            == getattr(CapabilityNativeSdkTable, field_name).offset
        )
    assert LijNativeSdkTable.evaluate_mixture_phase_lij.offset == 280
    assert ctypes.sizeof(LijNativeSdkTable) == 288
    for field_name, _ in LijNativeSdkTable._fields_:
        assert (
            getattr(PureParameterNativeSdkTable, field_name).offset
            == getattr(LijNativeSdkTable, field_name).offset
        )
    assert PureParameterNativeSdkTable.evaluate_pure_phase_parameter.offset == 288
    assert ctypes.sizeof(PureParameterNativeSdkTable) == 296
    for field_name, _ in PureParameterNativeSdkTable._fields_:
        assert (
            getattr(IonFractionSuppressionNativeSdkTable, field_name).offset
            == getattr(PureParameterNativeSdkTable, field_name).offset
        )
    assert IonFractionSuppressionNativeSdkTable.ion_fraction_suppression_result_size.offset == 296
    assert IonFractionSuppressionNativeSdkTable.evaluate_ion_fraction_suppression.offset == 304
    assert ION_FRACTION_SUPPRESSION_PREFIX_SIZE == ctypes.sizeof(
        IonFractionSuppressionNativeSdkTable
    ) == 312


def test_native_sdk_describes_the_exact_neutral_binary_pair_derivative_contracts() -> None:
    model = _mixture_model()
    capsule, table = _capability_table(model)
    assert capsule is not None
    assert table.capability_count == 2
    assert table.capabilities

    descriptor = table.capabilities[0]
    assert descriptor.struct_size == ctypes.sizeof(NativeCapabilityDescriptor)
    assert (
        descriptor.schema_version,
        descriptor.capability,
        descriptor.parameter_family,
        descriptor.parameter_identity,
        descriptor.observation_contract,
        descriptor.model_domain,
        descriptor.tensor_layout,
        descriptor.derivative_order,
        descriptor.maturity,
        descriptor.authority_effect,
    ) == (1, 1, 1, 1, 1, 1, 1, 2, 1, 0)
    assert descriptor.unsupported_status == STATUS_UNSUPPORTED_MODEL
    assert descriptor.domain_status == STATUS_DOMAIN_ERROR
    assert (
        descriptor.state_coordinate_count,
        descriptor.active_parameter_count,
        descriptor.coordinate_count,
        descriptor.component_count,
    ) == (3, 1, 4, 2)
    assert tuple(
        descriptor.component_ids[index].decode()
        for index in range(descriptor.component_count)
    ) == ("methane", "ethane")
    assert (
        descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint
    )
    assert descriptor.topology_fingerprint.decode() == (
        "sha256:172e3c9dce280763cf785db5d6fc94deeb30a23359346b2222a836077930351a"
    )
    assert descriptor.topology_fingerprint != descriptor.parameter_fingerprint
    assert descriptor.helmholtz_basis_id.decode() == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert (
        descriptor.temperature_min_k,
        descriptor.temperature_max_k,
    ) == (
        model._model.temperature_min_k,
        model._model.temperature_max_k,
    )

    coordinates = tuple(
        descriptor.coordinates[index] for index in range(descriptor.coordinate_count)
    )
    assert all(
        coordinate.struct_size == ctypes.sizeof(NativeCapabilityCoordinate)
        for coordinate in coordinates
    )
    assert tuple(coordinate.kind for coordinate in coordinates) == (1, 1, 2, 3)
    assert tuple(coordinate.unit.decode() for coordinate in coordinates) == (
        "mol",
        "mol",
        "m3",
        "dimensionless",
    )
    assert tuple(coordinate.component_index for coordinate in coordinates) == (
        0,
        1,
        -1,
        -1,
    )
    assert tuple(
        (coordinate.pair_component_index_a, coordinate.pair_component_index_b)
        for coordinate in coordinates
    ) == ((-1, -1), (-1, -1), (-1, -1), (0, 1))

    _, reversed_table = _capability_table(_mixture_model(("ethane", "methane")))
    assert reversed_table.capability_count == 2
    assert (
        reversed_table.capabilities[0].topology_fingerprint
        != descriptor.topology_fingerprint
    )

    lij = table.capabilities[1]
    assert (
        lij.schema_version,
        lij.capability,
        lij.parameter_family,
        lij.parameter_identity,
        lij.observation_contract,
        lij.model_domain,
        lij.tensor_layout,
        lij.derivative_order,
        lij.maturity,
        lij.authority_effect,
    ) == (1, 2, 2, 1, 1, 1, 1, 2, 1, 0)
    assert lij.component_count == 2
    assert tuple(
        lij.component_ids[index].decode()
        for index in range(lij.component_count)
    ) == ("methane", "ethane")
    assert lij.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert lij.topology_fingerprint.decode() == descriptor.topology_fingerprint.decode()
    lij_coordinates = tuple(
        lij.coordinates[index] for index in range(lij.coordinate_count)
    )
    assert tuple(coordinate.kind for coordinate in lij_coordinates) == (
        1,
        1,
        2,
        4,
    )
    assert tuple(coordinate.unit.decode() for coordinate in lij_coordinates) == (
        "mol",
        "mol",
        "m3",
        "dimensionless",
    )


def test_native_sdk_describes_exact_scalar_pure_parameter_contracts() -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    _, table = _capability_table(model)
    assert table.capability_count == 3
    expected = (
        (3, 3, 5, "dimensionless"),
        (4, 4, 6, "angstrom"),
        (5, 5, 7, "kelvin"),
    )
    for descriptor, (
        capability_code,
        family_code,
        coordinate_kind,
        unit,
    ) in zip(table.capabilities[:3], expected, strict=True):
        assert (
            descriptor.schema_version,
            descriptor.capability,
            descriptor.parameter_family,
            descriptor.parameter_identity,
            descriptor.observation_contract,
            descriptor.model_domain,
            descriptor.tensor_layout,
            descriptor.derivative_order,
            descriptor.maturity,
            descriptor.authority_effect,
        ) == (1, capability_code, family_code, 2, 1, 2, 1, 2, 1, 0)
        assert (
            descriptor.state_coordinate_count,
            descriptor.active_parameter_count,
            descriptor.coordinate_count,
            descriptor.component_count,
        ) == (2, 1, 3, 1)
        coordinates = tuple(
            descriptor.coordinates[index]
            for index in range(descriptor.coordinate_count)
        )
        assert tuple(item.kind for item in coordinates) == (
            1,
            2,
            coordinate_kind,
        )
        assert tuple(item.unit.decode() for item in coordinates) == (
            "mol",
            "m3",
            unit,
        )
        assert coordinates[-1].component_index == 0
        assert (
            descriptor.parameter_fingerprint.decode()
            == model.parameter_fingerprint
        )


def test_native_sdk_describes_exact_associating_pure_parameter_contracts() -> None:
    model = _model("figiel-2025-reference-electrolytes", "ethanol")
    _, table = _capability_table(model)
    assert table.capability_count == 2
    expected = (
        (15, 9, 11, "kelvin"),
        (16, 10, 12, "dimensionless"),
    )
    for descriptor, (
        capability_code,
        family_code,
        coordinate_kind,
        unit,
    ) in zip(table.capabilities[:2], expected, strict=True):
        assert (
            descriptor.schema_version,
            descriptor.capability,
            descriptor.parameter_family,
            descriptor.parameter_identity,
            descriptor.observation_contract,
            descriptor.model_domain,
            descriptor.tensor_layout,
            descriptor.derivative_order,
            descriptor.maturity,
            descriptor.authority_effect,
        ) == (1, capability_code, family_code, 3, 1, 7, 1, 2, 1, 0)
        coordinates = tuple(
            descriptor.coordinates[index]
            for index in range(descriptor.coordinate_count)
        )
        assert tuple(item.kind for item in coordinates) == (
            1,
            2,
            coordinate_kind,
        )
        assert tuple(item.unit.decode() for item in coordinates) == (
            "mol",
            "m3",
            unit,
        )
        assert coordinates[-1].component_index == -1
        assert descriptor.topology_fingerprint


@pytest.mark.parametrize(
    ("component_ids", "expected"),
    (
        (
            ("water", "sodium-cation", "chloride-anion"),
            (6, 6, 2, 3, 8, 1),
        ),
        (
            ("water", "sodium-cation", "bromide-anion"),
            (7, 7, 3, 4, 9, 0),
        ),
    ),
)
def test_native_sdk_describes_exact_direct_observable_parameter_contracts(
    component_ids: tuple[str, ...],
    expected: tuple[int, int, int, int, int, int],
) -> None:
    model = _electrolyte_model(component_ids)
    _, table = _capability_table(model)
    (
        capability_code,
        family_code,
        observation_code,
        model_domain_code,
        coordinate_kind,
        active_component_index,
    ) = expected
    matches = tuple(
        descriptor
        for descriptor in table.capabilities[: table.capability_count]
        if descriptor.capability == capability_code
    )
    assert len(matches) == 1
    descriptor = matches[0]
    assert (
        descriptor.schema_version,
        descriptor.parameter_family,
        descriptor.parameter_identity,
        descriptor.observation_contract,
        descriptor.model_domain,
        descriptor.tensor_layout,
        descriptor.derivative_order,
        descriptor.maturity,
        descriptor.authority_effect,
    ) == (1, family_code, 2, observation_code, model_domain_code, 1, 1, 1, 0)
    assert (
        descriptor.state_coordinate_count,
        descriptor.active_parameter_count,
        descriptor.coordinate_count,
        descriptor.component_count,
    ) == (0, 1, 1, 3)
    assert tuple(
        descriptor.component_ids[index].decode()
        for index in range(descriptor.component_count)
    ) == component_ids
    coordinate = descriptor.coordinates[0]
    assert (
        coordinate.kind,
        coordinate.component_index,
        coordinate.pair_component_index_a,
        coordinate.pair_component_index_b,
    ) == (coordinate_kind, active_component_index, -1, -1)
    assert coordinate.unit.decode() == (
        "angstrom" if coordinate_kind == 8 else "dimensionless"
    )
    assert (
        descriptor.temperature_min_k,
        descriptor.temperature_max_k,
    ) == (298.15, 298.15)
    assert descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert descriptor.topology_fingerprint.decode().startswith("sha256:")


def test_native_sdk_describes_each_aqueous_kij_as_a_scalar_miac_coordinate() -> None:
    component_ids = ("water", "sodium-cation", "bromide-anion")
    model = _electrolyte_model(component_ids)
    _, table = _capability_table(model)

    descriptors = tuple(
        descriptor
        for descriptor in table.capabilities[: table.capability_count]
        if descriptor.capability in (8, 9, 10)
    )

    assert tuple(descriptor.capability for descriptor in descriptors) == (8, 9, 10)
    assert tuple(
        (
            descriptor.coordinates[0].pair_component_index_a,
            descriptor.coordinates[0].pair_component_index_b,
        )
        for descriptor in descriptors
    ) == ((0, 1), (0, 2), (1, 2))
    for descriptor in descriptors:
        assert (
            descriptor.schema_version,
            descriptor.parameter_family,
            descriptor.parameter_identity,
            descriptor.observation_contract,
            descriptor.model_domain,
            descriptor.tensor_layout,
            descriptor.derivative_order,
            descriptor.maturity,
            descriptor.authority_effect,
        ) == (1, 1, 1, 3, 4, 1, 1, 1, 0)
        assert (
            descriptor.state_coordinate_count,
            descriptor.active_parameter_count,
            descriptor.coordinate_count,
            descriptor.component_count,
        ) == (0, 1, 1, 3)
        assert tuple(
            descriptor.component_ids[index].decode()
            for index in range(descriptor.component_count)
        ) == component_ids
        coordinate = descriptor.coordinates[0]
        assert (
            coordinate.kind,
            coordinate.component_index,
            coordinate.unit.decode(),
        ) == (3, -1, "dimensionless")
        assert (
            descriptor.temperature_min_k,
            descriptor.temperature_max_k,
        ) == (298.15, 298.15)
        assert descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint
        assert descriptor.topology_fingerprint.decode().startswith("sha256:")


def test_native_sdk_describes_ion_fraction_suppression_coordinate() -> None:
    component_ids = ("water", "sodium-cation", "bromide-anion")
    model = _electrolyte_model(component_ids)
    _, table = _capability_table(model)

    matches = tuple(
        descriptor
        for descriptor in table.capabilities[: table.capability_count]
        if descriptor.capability == 11
    )

    assert len(matches) == 1
    descriptor = matches[0]
    assert (
        descriptor.schema_version,
        descriptor.parameter_family,
        descriptor.parameter_identity,
        descriptor.observation_contract,
        descriptor.model_domain,
        descriptor.tensor_layout,
        descriptor.derivative_order,
        descriptor.maturity,
        descriptor.authority_effect,
    ) == (1, 8, 3, 4, 5, 1, 1, 1, 0)
    assert (
        descriptor.state_coordinate_count,
        descriptor.active_parameter_count,
        descriptor.coordinate_count,
        descriptor.component_count,
    ) == (0, 1, 1, 3)
    coordinate = descriptor.coordinates[0]
    assert (
        coordinate.kind,
        coordinate.component_index,
        coordinate.pair_component_index_a,
        coordinate.pair_component_index_b,
        coordinate.unit.decode(),
    ) == (10, -1, -1, -1, "dimensionless")

    methanol = _electrolyte_model(
        ("methanol", "sodium-cation", "bromide-anion")
    )
    _, methanol_table = _capability_table(methanol)
    methanol_matches = tuple(
        descriptor
        for descriptor in methanol_table.capabilities[
            : methanol_table.capability_count
        ]
        if descriptor.capability == 11
    )
    assert len(methanol_matches) == 1
    assert tuple(
        methanol_matches[0].component_ids[index].decode()
        for index in range(methanol_matches[0].component_count)
    ) == ("methanol", "sodium-cation", "bromide-anion")


def test_native_sdk_evaluates_ion_fraction_suppression_exact_derivative() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _ion_fraction_suppression_table(model)
    assert table.ion_fraction_suppression_result_size == ctypes.sizeof(
        IonFractionSuppressionResult
    )
    callback = EvaluateIonFractionSuppression(table.evaluate_ion_fraction_suppression)

    def evaluate(
        coefficient: float,
        *,
        fingerprint: bytes | None = None,
        temperature: float = 298.15,
        pressure: float = 100000.0,
        mole_fraction: float = 0.1,
    ) -> IonFractionSuppressionResult:
        result = IonFractionSuppressionResult(
            struct_size=ctypes.sizeof(IonFractionSuppressionResult)
        )
        status = callback(
            table.model_context,
            fingerprint or model.parameter_fingerprint.encode(),
            temperature,
            pressure,
            mole_fraction,
            coefficient,
            ctypes.byref(result),
        )
        assert status == result.status
        return result

    coefficient = 7.01
    mole_fraction = 0.1
    salt_free = 78.09
    denominator = 1.0 + coefficient * mole_fraction
    result = evaluate(coefficient)
    assert result.status == STATUS_OK
    assert result.salt_free_relative_permittivity == salt_free
    assert result.bulk_relative_permittivity == pytest.approx(
        salt_free / denominator
    )
    assert result.derivative == pytest.approx(
        -salt_free * mole_fraction / denominator**2
    )
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert not result.error

    step = 1.0e-5
    finite_difference = (
        evaluate(coefficient + step).bulk_relative_permittivity
        - evaluate(coefficient - step).bulk_relative_permittivity
    ) / (2.0 * step)
    assert result.derivative == pytest.approx(finite_difference, rel=2.0e-10)

    for failure, expected in (
        (evaluate(coefficient, fingerprint=b"sha256:wrong"), STATUS_UNSUPPORTED_MODEL),
        (evaluate(coefficient, temperature=300.0), STATUS_DOMAIN_ERROR),
        (evaluate(coefficient, pressure=101325.0), STATUS_DOMAIN_ERROR),
        (evaluate(coefficient, mole_fraction=-0.1), STATUS_DOMAIN_ERROR),
        (evaluate(-1.0), STATUS_DOMAIN_ERROR),
    ):
        assert failure.status == expected
        assert math.isnan(failure.bulk_relative_permittivity)
        assert math.isnan(failure.derivative)
        assert failure.parameter_fingerprint == b""
        assert failure.error

    methanol = _electrolyte_model(
        ("methanol", "sodium-cation", "bromide-anion")
    )
    _, methanol_table = _ion_fraction_suppression_table(methanol)
    methanol_result = IonFractionSuppressionResult(
        struct_size=ctypes.sizeof(IonFractionSuppressionResult)
    )
    methanol_status = EvaluateIonFractionSuppression(
        methanol_table.evaluate_ion_fraction_suppression
    )(
        methanol_table.model_context,
        methanol.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        mole_fraction,
        coefficient,
        ctypes.byref(methanol_result),
    )
    assert methanol_status == methanol_result.status == STATUS_OK
    assert methanol_result.salt_free_relative_permittivity == 33.05
    assert (
        methanol_result.bulk_relative_permittivity
        / methanol_result.salt_free_relative_permittivity
    ) == pytest.approx(1.0 / denominator)


def test_native_sdk_evaluates_organic_ion_solvation_kij_derivatives() -> None:
    model = _electrolyte_model(
        ("methanol", "potassium-cation", "bromide-anion")
    )
    _, table = _ion_solvation_kij_table(model)
    assert table.ion_solvation_kij_result_size == ctypes.sizeof(
        IonSolvationKijResult
    )
    callback = EvaluateIonSolvationKij(table.evaluate_ion_solvation_kij)

    def evaluate(
        coordinates: tuple[float, float, float],
        *,
        active_ion: int = 1,
    ) -> IonSolvationKijResult:
        values = (ctypes.c_double * 3)(*coordinates)
        result = IonSolvationKijResult(
            struct_size=ctypes.sizeof(IonSolvationKijResult)
        )
        status = callback(
            table.model_context,
            model.parameter_fingerprint.encode(),
            298.15,
            100000.0,
            active_ion,
            values,
            3,
            ctypes.byref(result),
        )
        assert status == result.status
        return result

    coordinates = (0.32, 0.15, -0.35)
    result = evaluate(coordinates)
    assert result.status == STATUS_OK
    assert math.isfinite(result.solvation_gibbs_j_per_mol)
    assert all(math.isfinite(value) for value in result.derivative)
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint

    step = 1.0e-5
    for active in range(3):
        low = list(coordinates)
        high = list(coordinates)
        low[active] -= step
        high[active] += step
        finite_difference = (
            evaluate(tuple(high)).solvation_gibbs_j_per_mol
            - evaluate(tuple(low)).solvation_gibbs_j_per_mol
        ) / (2.0 * step)
        assert result.derivative[active] == pytest.approx(
            finite_difference, rel=2.0e-6, abs=1.0e-3
        )

    _, capability_table = _capability_table(model)
    matches = tuple(
        descriptor
        for descriptor in capability_table.capabilities[
            : capability_table.capability_count
        ]
        if descriptor.capability in (12, 13, 14)
    )
    assert len(matches) == 3
    assert all(
        descriptor.parameter_family == 1
        and descriptor.parameter_identity == 1
        and descriptor.observation_contract == 2
        and descriptor.model_domain == 6
        and descriptor.derivative_order == 1
        for descriptor in matches
    )
    assert {
        (
            descriptor.coordinates[0].pair_component_index_a,
            descriptor.coordinates[0].pair_component_index_b,
        )
        for descriptor in matches
    } == {(0, 1), (0, 2), (1, 2)}


def test_native_sdk_rejects_forced_organic_ion_solvation_kij_admission() -> None:
    unsupported = _electrolyte_model(
        ("water", "sodium-cation", "bromide-anion")
    )
    with pytest.raises(
        ValueError,
        match="unsupported active ion-solvation k_ij model identity",
    ):
        _native_sdk_v1(
            unsupported._native,
            unsupported.parameter_fingerprint,
            298.15,
            298.15,
            unsupported.component_ids,
            0.38,
            phase_callbacks_admitted=False,
            ion_solvation_kij_admitted=True,
            topology_fingerprint=(
                "sha256:a6879554a81191d2210ece2eed240f251a0f45327fa3877978228362bcf80443"
            ),
        )

    admitted = _electrolyte_model(
        ("methanol", "potassium-cation", "bromide-anion")
    )
    with pytest.raises(
        ValueError,
        match="unsupported active ion-solvation k_ij model identity",
    ):
        _native_sdk_v1(
            admitted._native,
            admitted.parameter_fingerprint,
            298.15,
            298.15,
            admitted.component_ids,
            0.38,
            phase_callbacks_admitted=False,
            ion_solvation_kij_admitted=True,
            topology_fingerprint="sha256:altered",
        )

    with pytest.raises(
        ValueError,
        match="unsupported active ion-solvation k_ij model identity",
    ):
        _native_sdk_v1(
            admitted._native,
            "sha256:altered",
            298.15,
            298.15,
            admitted.component_ids,
            0.38,
            phase_callbacks_admitted=False,
            ion_solvation_kij_admitted=True,
            topology_fingerprint=(
                "sha256:a6879554a81191d2210ece2eed240f251a0f45327fa3877978228362bcf80443"
            ),
        )


def test_native_sdk_ion_fraction_suppression_installed_header_layout() -> None:
    compiler = shutil.which("cc")
    assert compiler is not None
    include_root = Path(epcsaft.__file__).parent / "include"
    source = r"""
#include <stddef.h>
#include <epcsaft/native_sdk_v1.h>

_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_ION_FRACTION_SUPPRESSION_V1 == 11u,
    "dielectric capability"
);
_Static_assert(
    EPCSAFT_NATIVE_PARAMETER_FAMILY_ION_FRACTION_SUPPRESSION_V1 == 8u,
    "dielectric parameter family"
);
_Static_assert(
    EPCSAFT_NATIVE_OBSERVATION_RELATIVE_PERMITTIVITY_RATIO_V1 == 4u,
    "relative-permittivity-ratio observation"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, ion_fraction_suppression_result_size) == 296,
    "result-size offset"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, evaluate_ion_fraction_suppression) == 304,
    "callback offset"
);
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 312, "table prefix size");
_Static_assert(
    sizeof(epcsaft_ion_fraction_suppression_result_v1) == 264,
    "result size"
);
_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_CATION_KIJ_V1 == 12u,
    "ion-solvation solvent-cation capability"
);
_Static_assert(
    EPCSAFT_NATIVE_MODEL_DOMAIN_ORGANIC_SINGLE_ION_SOLVATION_V1 == 6u,
    "single-ion-solvation model domain"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, ion_solvation_kij_result_size) == 312,
    "ion-solvation result-size offset"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, evaluate_ion_solvation_kij) == 320,
    "ion-solvation callback offset"
);
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 328, "extended table size");
_Static_assert(
    sizeof(epcsaft_ion_solvation_kij_result_v1) == 288,
    "ion-solvation result size"
);
_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_IONIC_REGION_PERMITTIVITY_V1 == 17u,
    "ionic-region-permittivity capability"
);
_Static_assert(
    EPCSAFT_NATIVE_PARAMETER_FAMILY_IONIC_REGION_RELATIVE_PERMITTIVITY_V1 == 11u,
    "ionic-region-permittivity parameter family"
);
_Static_assert(
    offsetof(
        epcsaft_native_sdk_v1,
        ion_solvation_ionic_permittivity_result_size
    ) == 328,
    "ionic-region-permittivity result-size offset"
);
_Static_assert(
    offsetof(
        epcsaft_native_sdk_v1,
        evaluate_ion_solvation_ionic_permittivity
    ) == 336,
    "ionic-region-permittivity callback offset"
);
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 344, "final table size");
_Static_assert(
    sizeof(epcsaft_ion_solvation_ionic_permittivity_result_v1) == 272,
    "ionic-region-permittivity result size"
);
_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_PERMITTIVITY_V1 == 18u,
    "solvent-permittivity capability"
);
_Static_assert(
    EPCSAFT_NATIVE_PARAMETER_FAMILY_RELATIVE_PERMITTIVITY_V1 == 12u,
    "solvent-permittivity parameter family"
);
_Static_assert(
    offsetof(
        epcsaft_native_sdk_v1,
        ion_solvation_solvent_permittivity_result_size
    ) == 344,
    "solvent-permittivity result-size offset"
);
_Static_assert(
    offsetof(
        epcsaft_native_sdk_v1,
        evaluate_ion_solvation_solvent_permittivity
    ) == 352,
    "solvent-permittivity callback offset"
);
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 360, "final table size");
_Static_assert(
    sizeof(epcsaft_ion_solvation_solvent_permittivity_result_v1) == 272,
    "solvent-permittivity result size"
);

int main(void) {
    epcsaft_evaluate_ion_fraction_suppression_v1 callback = 0;
    epcsaft_evaluate_ion_solvation_kij_v1 ion_callback = 0;
    epcsaft_evaluate_ion_solvation_ionic_permittivity_v1
        ionic_permittivity_callback = 0;
    epcsaft_evaluate_ion_solvation_solvent_permittivity_v1
        solvent_permittivity_callback = 0;
    return callback != 0 || ion_callback != 0
        || ionic_permittivity_callback != 0
        || solvent_permittivity_callback != 0;
}
"""
    completed = subprocess.run(
        (
            compiler,
            "-std=c11",
            "-Werror",
            f"-I{include_root}",
            "-x",
            "c",
            "-fsyntax-only",
            "-",
        ),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_native_sdk_scalar_pure_parameter_installed_header_layout() -> None:
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler is required for the installed-header contract")
    include_root = Path(epcsaft.__file__).parent / "include"
    source = r"""
#include <stddef.h>
#include <epcsaft/native_sdk_v1.h>

_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1 == 3u,
    "segment-count capability"
);
_Static_assert(
    EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1 == 4u,
    "segment-diameter family"
);
_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1 == 7u,
    "dispersion-energy coordinate"
);
_Static_assert(
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1 == 6u,
    "Born direct-observable capability"
);
_Static_assert(
    EPCSAFT_NATIVE_OBSERVATION_AQUEOUS_MEAN_IONIC_ACTIVITY_V1 == 3u,
    "aqueous mean-ionic-activity observation"
);
_Static_assert(
    offsetof(epcsaft_native_sdk_v1, evaluate_pure_phase_parameter) == 288,
    "pure-parameter callback offset"
);
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 296, "table prefix size");

int main(void) {
    epcsaft_evaluate_pure_phase_parameter_v1 callback = 0;
    return callback != 0;
}
"""
    completed = subprocess.run(
        (
            compiler,
            "-std=c11",
            "-Werror",
            f"-I{include_root}",
            "-x",
            "c",
            "-fsyntax-only",
            "-",
        ),
        input=source,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    ("family_code", "parameter_value", "full_coordinate"),
    (
        (3, 1.0, 2),
        (4, 3.7039, 3),
        (5, 150.03, 4),
    ),
)
def test_native_sdk_scalar_pure_parameter_callback_matches_full_exact_block(
    family_code: int,
    parameter_value: float,
    full_coordinate: int,
) -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    capsule, table = _pure_parameter_table(model)
    assert capsule is not None
    temperature_k = 150.0
    point = (1.0, 2.0e-3, parameter_value)
    direction = (0.01, -2.0e-5, 0.02)

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(
            value + scale * delta
            for value, delta in zip(point, direction, strict=True)
        )
        return _evaluate_pure_parameter(
            table,
            temperature_k,
            values[0],
            values[1],
            family_code,
            values[2],
        )

    center, gradient, hessian = shifted(0.0)
    full = _evaluate_parameters(
        table,
        temperature_k,
        point[0],
        point[1],
        parameter_value if family_code == 3 else 1.0,
        parameter_value if family_code == 4 else 3.7039,
        parameter_value if family_code == 5 else 150.03,
    )
    selected = (0, 1, full_coordinate)
    assert center.status == STATUS_OK
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        full.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert center.pressure_pa == pytest.approx(full.pressure_pa, rel=2.0e-13)
    assert gradient == pytest.approx(tuple(full.gradient[index] for index in selected))
    assert hessian == pytest.approx(
        tuple(full.hessian[row * 5 + column] for row in selected for column in selected)
    )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    value_direction = sum(
        value * delta for value, delta in zip(gradient, direction, strict=True)
    )
    value_difference = (
        upper.helmholtz_over_rt_reference_amount
        - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert value_direction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(3):
        hessian_direction = sum(
            hessian[row * 3 + column] * direction[column]
            for column in range(3)
        )
        gradient_difference = (
            upper_gradient[row] - lower_gradient[row]
        ) / (2.0 * step)
        assert hessian_direction == pytest.approx(
            gradient_difference, rel=3.0e-7
        )


@pytest.mark.parametrize(
    ("family_code", "parameter_value", "parameter_step"),
    (
        (9, 2653.4, 2.0e-3),
        (10, 0.03238, 2.0e-8),
    ),
)
def test_native_sdk_associating_pure_parameter_callback_has_exact_hessian(
    family_code: int,
    parameter_value: float,
    parameter_step: float,
) -> None:
    model = _model("figiel-2025-reference-electrolytes", "ethanol")
    capsule, table = _pure_parameter_table(model)
    assert capsule is not None
    temperature_k = 298.15
    point = (1.0, 5.86e-5, parameter_value)
    direction = (1.0e-3, -2.0e-8, parameter_step)

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(
            value + scale * delta
            for value, delta in zip(point, direction, strict=True)
        )
        return _evaluate_pure_parameter(
            table,
            temperature_k,
            values[0],
            values[1],
            family_code,
            values[2],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert all(math.isfinite(value) for value in gradient)
    assert all(math.isfinite(value) for value in hessian)
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[1], rel=2.0e-13
    )
    provider_state = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=(point[0] / point[1])
        * unit_registry.mole
        / unit_registry.meter**3,
        x=(1.0,),
    )
    assert center.pressure_pa == pytest.approx(
        provider_state.pressure.to("pascal").magnitude,
        rel=3.0e-12,
    )

    difference_step = 1.0e-2
    lower, lower_gradient, _ = shifted(-difference_step)
    upper, upper_gradient, _ = shifted(difference_step)
    value_direction = sum(
        value * delta for value, delta in zip(gradient, direction, strict=True)
    )
    value_difference = (
        upper.helmholtz_over_rt_reference_amount
        - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * difference_step)
    assert value_direction == pytest.approx(value_difference, rel=2.0e-6)
    for row in range(3):
        hessian_direction = sum(
            hessian[row * 3 + column] * direction[column]
            for column in range(3)
        )
        gradient_difference = (
            upper_gradient[row] - lower_gradient[row]
        ) / (2.0 * difference_step)
        assert hessian_direction == pytest.approx(
            gradient_difference, rel=2.0e-5, abs=2.0e-8
        )


def test_native_sdk_joint_associating_parameters_have_exact_seven_coordinate_hessian(
) -> None:
    model = _model("figiel-2025-reference-electrolytes", "ethanol")
    capsule, table = _associating_parameterized_table(model)
    assert capsule is not None
    assert table.associating_parameterized_result_size == ctypes.sizeof(
        AssociatingParameterizedPhaseBlockResult
    )
    assert table.evaluate_diameter_basis_associating_pure_phase_parameters

    temperature_k = 298.15
    point = (1.0, 5.86e-5, 2.3827, 3.1771, 198.24, 2653.4, 0.03238)
    direction = (1.0e-3, -2.0e-8, 2.0e-4, -1.0e-4, 2.0e-2, 2.0e-1, 2.0e-6)

    def shifted(scale: float) -> AssociatingParameterizedPhaseBlockResult:
        return _evaluate_associating_parameters(
            table,
            temperature_k,
            tuple(
                value + scale * delta
                for value, delta in zip(point, direction, strict=True)
            ),
        )

    center = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert center.topology_fingerprint
    assert all(math.isfinite(value) for value in center.gradient)
    assert all(math.isfinite(value) for value in center.hessian)
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * center.gradient[1], rel=2.0e-13
    )

    state = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=(point[0] / point[1])
        * unit_registry.mole
        / unit_registry.meter**3,
        x=(1.0,),
    )
    density = point[0] / point[1]
    number_density = density * 6.02214076e-7
    diameter = point[3] * (
        1.0 - 0.12 * math.exp(-3.0 * point[4] / temperature_k)
    )
    zeta_2 = math.pi / 6.0 * number_density * point[2] * diameter**2
    zeta_3 = math.pi / 6.0 * number_density * point[2] * diameter**3
    contact = (
        1.0 / (1.0 - zeta_3)
        + 1.5 * diameter * zeta_2 / (1.0 - zeta_3) ** 2
        + 0.5 * diameter**2 * zeta_2**2 / (1.0 - zeta_3) ** 3
    )
    delta = (
        contact
        * math.expm1(point[5] / temperature_k)
        * diameter**3
        * point[6]
    )
    site_fraction = 2.0 / (
        1.0 + math.sqrt(1.0 + 4.0 * number_density * delta)
    )
    diameter_basis_association = 2.0 * (
        math.log(site_fraction) - 0.5 * site_fraction + 0.5
    )
    expected_phase = point[0] * (
        math.log(density)
        - 1.0
        + state.hard_chain
        + state.dispersion
        + diameter_basis_association
    )
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_phase, rel=2.0e-12
    )
    assert center.helmholtz_over_rt_reference_amount != pytest.approx(
        math.log(density)
        - 1.0
        + state._residual_helmholtz,
        rel=1.0e-5,
    )

    step = 1.0e-2
    lower = shifted(-step)
    upper = shifted(step)
    value_direction = sum(
        value * delta
        for value, delta in zip(center.gradient, direction, strict=True)
    )
    value_difference = (
        upper.helmholtz_over_rt_reference_amount
        - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert value_direction == pytest.approx(value_difference, rel=2.0e-6)
    for row in range(7):
        hessian_direction = sum(
            center.hessian[row * 7 + column] * direction[column]
            for column in range(7)
        )
        gradient_difference = (
            upper.gradient[row] - lower.gradient[row]
        ) / (2.0 * step)
        assert hessian_direction == pytest.approx(
            gradient_difference, rel=3.0e-5, abs=2.0e-8
        )


def test_native_sdk_joint_associating_parameters_reject_nonassociating_model(
) -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    _, table = _associating_parameterized_table(model)
    result = _evaluate_associating_parameters(
        table,
        150.0,
        (1.0, 2.0e-3, 1.0, 3.7039, 150.03, 2500.0, 0.03),
    )
    assert result.status == STATUS_UNSUPPORTED_MODEL
    assert b"2B model" in result.error


def test_native_sdk_joint_associating_parameters_require_one_of_each_2b_site(
) -> None:
    model = _model("figiel-2025-reference-electrolytes", "ethanol")
    malformed = NativeModel(
        *_native_arguments(replace(model._model, site_multiplicity=(2, 1)))
    )
    with pytest.raises(
        ValueError, match="native SDK v1 requires a supported phase-block model"
    ):
        _native_sdk_v1(
            malformed,
            model.parameter_fingerprint,
            298.15,
            298.15,
            model.component_ids,
            math.nan,
            phase_callbacks_admitted=True,
            topology_fingerprint=f"sha256:{'0' * 64}",
        )


@pytest.mark.parametrize(
    ("family_code", "parameter_value", "expected_status"),
    (
        (999, 1.0, STATUS_INVALID_ARGUMENT),
        (3, math.nan, STATUS_INVALID_ARGUMENT),
        (3, 0.0, STATUS_DOMAIN_ERROR),
    ),
)
def test_native_sdk_scalar_pure_parameter_callback_rejects_invalid_inputs(
    family_code: int,
    parameter_value: float,
    expected_status: int,
) -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    capsule, table = _pure_parameter_table(model)
    assert capsule is not None
    result, _, _ = _evaluate_pure_parameter(
        table, 150.0, 1.0, 2.0e-3, family_code, parameter_value
    )
    assert result.status == expected_status


def test_native_sdk_scalar_pure_parameter_callback_fails_closed_on_malformed_buffers(
) -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    capsule, table = _pure_parameter_table(model)
    assert capsule is not None
    callback = EvaluatePurePhaseParameter(table.evaluate_pure_phase_parameter)

    def call(
        result: MixturePhaseBlockResult,
        *,
        context: object = table.model_context,
        temperature_k: float = 150.0,
    ) -> int:
        return callback(
            context,
            temperature_k,
            1.0,
            2.0e-3,
            3,
            1.0,
            ctypes.byref(result),
        )

    result, _, _ = _mixture_result(3)
    result.gradient_capacity = 2
    assert call(result) == STATUS_INVALID_ARGUMENT
    assert result.status == STATUS_INVALID_ARGUMENT
    assert result.error

    result, _, _ = _mixture_result(3)
    assert call(result, context=None) == STATUS_INVALID_ARGUMENT
    assert result.status == STATUS_INVALID_ARGUMENT
    assert result.error

    result, _, _ = _mixture_result(3)
    assert call(result, temperature_k=0.0) == STATUS_INVALID_ARGUMENT
    assert result.status == STATUS_INVALID_ARGUMENT
    assert result.error

    result, _, _ = _mixture_result(3)
    result.struct_size -= 1
    assert call(result) == STATUS_ABI_MISMATCH


def test_native_sdk_held_neutral_reference_matches_the_public_owner() -> None:
    model = _held_water_ionization_model()
    capsule, table = _neutral_reference_table(model)
    assert capsule is not None
    observed, basis, contractions, composition = _evaluate_neutral_reference(
        table, model.parameter_fingerprint
    )
    reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=100000.0 * unit_registry.pascal,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    pure_water = model.state(
        T=reference.temperature,
        P=reference.pressure,
        x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    assert pure_water._log_fugacity_coefficient is not None

    assert observed.status == STATUS_OK
    assert observed.temperature_k == 298.15
    assert observed.pressure_pa == 100000.0
    assert basis == pytest.approx((1.0, 0.0, 0.0, 0.0, 1.0, 1.0))
    assert contractions == pytest.approx(
        (
            pure_water._log_fugacity_coefficient[0],
            math.fsum(reference.infinite_dilution_log_fugacity_coefficient),
        ),
        abs=2.0e-12,
    )
    assert composition == pytest.approx((1.0, 0.0, 0.0))
    assert observed.solvent_molar_mass_kg_per_mol == pytest.approx(0.01801528, abs=1.0e-12)
    assert observed.reference_molality_mol_per_kg == 1.0e-12
    assert 0.0 <= observed.reference_convergence_error <= 5.0e-5
    assert 1.0 / pure_water.molar_density.to("mole / meter**3").magnitude > 0.0
    assert observed.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert observed.helmholtz_basis_id.decode() == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert observed.error == b""


def test_native_sdk_held_neutral_reference_is_domain_bound_and_permutation_safe() -> None:
    model = _held_water_ionization_model()
    _, table = _neutral_reference_table(model)
    result, _, _, _ = _evaluate_neutral_reference(
        table,
        model.parameter_fingerprint,
        temperature_k=299.0,
    )
    assert result.status == STATUS_DOMAIN_ERROR
    assert result.error

    permuted = _held_water_ionization_model(
        ("hydroxide-anion", "water", "hydronium-cation")
    )
    _, permuted_table = _neutral_reference_table(permuted)
    permuted_result, permuted_basis, permuted_contractions, composition = (
        _evaluate_neutral_reference(
            permuted_table,
            permuted.parameter_fingerprint,
        )
    )
    assert permuted_result.status == STATUS_OK
    assert permuted_basis == pytest.approx((0.0, 1.0, 0.0, 1.0, 0.0, 1.0))
    canonical_reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=100000.0 * unit_registry.pascal,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    assert permuted_contractions[1] == pytest.approx(
        math.fsum(canonical_reference.infinite_dilution_log_fugacity_coefficient),
        abs=5.0e-12,
    )
    assert composition == pytest.approx((0.0, 1.0, 0.0))


@pytest.mark.parametrize(
    ("catalog", "component_ids", "domain_id", "expected_fingerprint", "expected"),
    (
        (
            "figiel-2025-reference-electrolytes",
            ("water", "sodium-cation", "chloride-anion"),
            "figiel-298-k-1-bar",
            "sha256:d01bc3e64f19293902c00d9cc52060fccc3945ce736f0cd3e2f756be226cf5bf",
            (298.15, 298.15, 0.38),
        ),
        (
            "khudaida-2026-figure-2-electrolyte-lle",
            ("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion"),
            "khudaida-293-15-k-0-1-mpa",
            "sha256:b43fac77754d9d5cca8b3db2cbe709892a786d97b756084b167ce126ab4c3007",
            (293.15, 293.15, math.nan),
        ),
    ),
)
def test_native_sdk_applicability_metadata_matches_source_and_resolved_model(
    catalog: str,
    component_ids: tuple[str, ...],
    domain_id: str,
    expected_fingerprint: str,
    expected: tuple[float, float, float],
) -> None:
    bundle = _ParameterBundle.from_catalog(catalog, version=1)
    source_domain = next(domain for domain in bundle.domains if domain.domain_id == domain_id)
    model = Mixture(bundle.select(component_ids))
    capsule, table = _applicability_table(model)
    assert capsule is not None
    assert table.table_size >= APPLICABILITY_PREFIX_SIZE
    assert model.parameter_fingerprint == expected_fingerprint

    source_minimum = float(source_domain.temperature_min.to("kelvin").magnitude)
    source_maximum = float(source_domain.temperature_max.to("kelvin").magnitude)
    source_ion_maximum = (
        math.nan
        if source_domain.ion_mole_fraction_max is None
        else float(source_domain.ion_mole_fraction_max)
    )
    observed = (
        table.source_temperature_min_k,
        table.source_temperature_max_k,
        table.total_ion_mole_fraction_max,
    )
    assert observed[:2] == expected[:2] == (source_minimum, source_maximum)
    assert observed[:2] == (
        model._model.temperature_min_k,
        model._model.temperature_max_k,
    )
    assert all(math.isfinite(value) and value > 0.0 for value in observed[:2])
    assert observed[0] <= observed[1]
    if math.isnan(expected[2]):
        assert math.isnan(observed[2])
        assert math.isnan(source_ion_maximum)
        assert model._model.ion_mole_fraction_max is None
    else:
        assert observed[2] == expected[2] == source_ion_maximum
        assert observed[2] == model._model.ion_mole_fraction_max
        assert 0.0 < observed[2] <= 1.0


@pytest.mark.parametrize(
    ("temperature_min_k", "temperature_max_k", "ion_mole_fraction_max"),
    (
        (math.nan, math.nan, 0.38),
        (math.nan, 298.15, 0.38),
        (298.15, math.nan, 0.38),
        (299.0, 298.15, 0.38),
        (298.15, 298.15, 0.0),
        (298.15, 298.15, 1.01),
        (298.15, 298.15, math.inf),
    ),
)
def test_native_sdk_applicability_metadata_rejects_missing_or_malformed_limits(
    temperature_min_k: float,
    temperature_max_k: float,
    ion_mole_fraction_max: float,
) -> None:
    model = _electrolyte_model(("water", "sodium-cation", "chloride-anion"))
    with pytest.raises(ValueError):
        _native_sdk_v1(
            model._native,
            model.parameter_fingerprint,
            temperature_min_k,
            temperature_max_k,
            model.component_ids,
            ion_mole_fraction_max,
        )


@pytest.mark.parametrize(
    ("pressure_min_pa", "pressure_max_pa"),
    (
        (0.0, math.nan),
        (math.nan, 0.0),
        (math.inf, math.nan),
        (math.nan, math.inf),
        (200000.0, 100000.0),
    ),
)
def test_native_sdk_rejects_malformed_pressure_limits(
    pressure_min_pa: float,
    pressure_max_pa: float,
) -> None:
    model = _electrolyte_model(("water", "sodium-cation", "chloride-anion"))
    with pytest.raises(
        ValueError, match="invalid neutral reference pressure applicability"
    ):
        _native_sdk_v1(
            model._native,
            model.parameter_fingerprint,
            298.15,
            298.15,
            model.component_ids,
            0.38,
            pressure_min_pa=pressure_min_pa,
            pressure_max_pa=pressure_max_pa,
            neutral_reference_admitted=True,
        )


def test_native_sdk_mixture_phase_block_matches_public_thermodynamics_and_derivatives() -> None:
    model = _mixture_model()
    capsule, table = _table(model)
    assert capsule is not None
    assert table.component_count == 2
    assert table.mixture_result_size == ctypes.sizeof(MixturePhaseBlockResult)
    assert table.evaluate_mixture_phase

    temperature_k = 220.0
    point = (0.4, 0.6, 1.25e-4)
    direction = (0.02, -0.01, -2.0e-6)

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        return _evaluate_mixture(
            table,
            temperature_k,
            (point[0] + scale * direction[0], point[1] + scale * direction[1]),
            point[2] + scale * direction[2],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 3
    assert center.gradient_capacity == 3
    assert center.hessian_capacity == 9
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint

    total = point[0] + point[1]
    density = total / point[2]
    composition = (point[0] / total, point[1] / total)
    state = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=composition,
    )
    expected = sum(amount * (math.log(amount / point[2]) - 1.0) for amount in point[:2])
    expected += total * state._residual_helmholtz
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected, rel=2.0e-12, abs=2.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=2.0e-11, abs=1.0e-5
    )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_contraction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_contraction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(3):
        hessian_direction = sum(
            hessian[row * 3 + column] * direction[column] for column in range(3)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


def test_native_sdk_mixture_phase_block_is_component_permutation_invariant() -> None:
    forward_capsule, forward_table = _table(_mixture_model())
    reverse_capsule, reverse_table = _table(_mixture_model(("ethane", "methane")))
    assert forward_capsule is not None
    assert reverse_capsule is not None
    forward, forward_gradient, forward_hessian = _evaluate_mixture(
        forward_table, 220.0, (0.4, 0.6), 1.25e-4
    )
    reverse, reverse_gradient, reverse_hessian = _evaluate_mixture(
        reverse_table, 220.0, (0.6, 0.4), 1.25e-4
    )
    permutation = (1, 0, 2)
    assert reverse.helmholtz_over_rt_reference_amount == pytest.approx(
        forward.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert reverse.pressure_pa == pytest.approx(forward.pressure_pa, rel=2.0e-13)
    for row in range(3):
        assert reverse_gradient[row] == pytest.approx(
            forward_gradient[permutation[row]], rel=2.0e-13, abs=2.0e-13
        )
        for column in range(3):
            assert reverse_hessian[row * 3 + column] == pytest.approx(
                forward_hessian[permutation[row] * 3 + permutation[column]],
                rel=2.0e-13,
                abs=2.0e-13,
            )


@pytest.mark.parametrize(
    ("temperature_k", "point", "direction"),
    (
        (220.0, (0.4, 0.6, 1.25e-4, 0.0), (0.02, -0.01, -2.0e-6, 0.01)),
        (280.0, (0.7, 0.3, 2.5e-4, -0.05), (-0.015, 0.02, 3.0e-6, -0.008)),
    ),
)
def test_native_sdk_active_kij_matches_fixed_and_directional_oracles(
    temperature_k: float,
    point: tuple[float, float, float, float],
    direction: tuple[float, float, float, float],
) -> None:
    model = _mixture_model()
    capsule, table = _table(model)
    assert capsule is not None
    assert table.evaluate_mixture_phase_kij

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(value + scale * delta for value, delta in zip(point, direction))
        return _evaluate_mixture_kij(
            table,
            temperature_k,
            (values[0], values[1]),
            values[2],
            values[3],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 4
    assert center.gradient_capacity == 4
    assert center.hessian_capacity == 16
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[2], rel=2.0e-13
    )

    if point[3] == 0.0:
        fixed, fixed_gradient, fixed_hessian = _evaluate_mixture(
            table, temperature_k, point[:2], point[2]
        )
        assert center.helmholtz_over_rt_reference_amount == pytest.approx(
            fixed.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
        )
        assert center.pressure_pa == pytest.approx(fixed.pressure_pa, rel=2.0e-13)
        for row in range(3):
            assert gradient[row] == pytest.approx(fixed_gradient[row], rel=2.0e-13, abs=2.0e-13)
            for column in range(3):
                assert hessian[row * 4 + column] == pytest.approx(
                    fixed_hessian[row * 3 + column], rel=2.0e-13, abs=2.0e-13
                )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_direction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_direction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(4):
        hessian_direction = sum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


@pytest.mark.parametrize(
    ("temperature_k", "point", "direction"),
    (
        (220.0, (0.4, 0.6, 1.25e-4, 0.0), (0.02, -0.01, -2.0e-6, 0.01)),
        (280.0, (0.7, 0.3, 2.5e-4, -0.02), (-0.015, 0.02, 3.0e-6, -0.008)),
    ),
)
def test_native_sdk_active_lij_matches_fixed_and_directional_oracles(
    temperature_k: float,
    point: tuple[float, float, float, float],
    direction: tuple[float, float, float, float],
) -> None:
    model = _mixture_model()
    capsule, table = _lij_table(model)
    assert capsule is not None

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(
            value + scale * delta
            for value, delta in zip(point, direction, strict=True)
        )
        return _evaluate_mixture_lij(
            table,
            temperature_k,
            (values[0], values[1]),
            values[2],
            values[3],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 4
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[2], rel=2.0e-13
    )
    if point[3] == 0.0:
        fixed, fixed_gradient, fixed_hessian = _evaluate_mixture(
            table, temperature_k, point[:2], point[2]
        )
        assert center.helmholtz_over_rt_reference_amount == pytest.approx(
            fixed.helmholtz_over_rt_reference_amount,
            rel=2.0e-13,
            abs=2.0e-13,
        )
        assert center.pressure_pa == pytest.approx(
            fixed.pressure_pa, rel=2.0e-13
        )
        for row in range(3):
            assert gradient[row] == pytest.approx(
                fixed_gradient[row], rel=2.0e-13, abs=2.0e-13
            )
            for column in range(3):
                assert hessian[row * 4 + column] == pytest.approx(
                    fixed_hessian[row * 3 + column],
                    rel=2.0e-13,
                    abs=2.0e-13,
                )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_direction = sum(
        value * delta for value, delta in zip(gradient, direction)
    )
    value_difference = (
        upper.helmholtz_over_rt_reference_amount
        - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_direction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(4):
        hessian_direction = sum(
            hessian[row * 4 + column] * direction[column]
            for column in range(4)
        )
        gradient_difference = (
            upper_gradient[row] - lower_gradient[row]
        ) / (2.0 * step)
        assert hessian_direction == pytest.approx(
            gradient_difference, rel=3.0e-7
        )


def test_native_sdk_active_kij_matches_independent_high_precision_oracle() -> None:
    temperature_k = 220.0
    point = (0.4, 0.6, 1.25e-4, 0.037)
    _, table = _table(_mixture_model())
    result, gradient, hessian = _evaluate_mixture_kij(
        table, temperature_k, point[:2], point[2], point[3]
    )
    # Frozen from an independent 70-digit Decimal transcription of Gross (2001),
    # Eqs. A.6-A.17, using fourth-order five-point derivatives with steps
    # (1e-5, 1e-5, 1e-9, 1e-5). The reference imports no provider implementation.
    expected_value = 6.001964399394894
    expected_gradient = (
        6.62799207335156,
        5.5565900739759755,
        134.50820534947894,
        0.9913299968451905,
    )
    expected_hessian = (
        1.98583616542871,
        -1.1078284190911964,
        -1037.0993177341288,
        2.4172796297654844,
        -1.1078284190911964,
        -0.9438004026648041,
        8075.292873882889,
        1.4899310731125321,
        -1037.0993177341288,
        8075.292873882889,
        -35442687.97788865,
        -6956.323991428179,
        2.4172796297654844,
        1.4899310731125321,
        -6956.323991428179,
        -0.20128500566962723,
    )

    assert result.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_value, rel=2.0e-12, abs=2.0e-12
    )
    for actual, expected in zip(gradient, expected_gradient):
        assert actual == pytest.approx(expected, rel=2.0e-11, abs=2.0e-11)
    for actual, expected in zip(hessian, expected_hessian):
        assert actual == pytest.approx(expected, rel=2.0e-11, abs=2.0e-10)


def test_native_sdk_parameterized_phase_block_has_exact_parameter_derivatives() -> None:
    model = _model("gross-2001-methane-ethane", "methane")
    _, table = _table(model)
    temperature_k = 150.0
    point = (1.0, 2.0e-3, 1.0, 3.7039, 150.03)
    direction = (0.01, -2.0e-5, 0.02, -0.03, 0.5)

    def shifted(scale: float) -> ParameterizedPhaseBlockResult:
        values = tuple(value + scale * delta for value, delta in zip(point, direction))
        result = _evaluate_parameters(table, temperature_k, *values)
        assert result.status == STATUS_OK
        return result

    center = shifted(0.0)
    fixed = _evaluate(table, temperature_k, point[0], point[1])
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        fixed.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert center.pressure_pa == pytest.approx(fixed.pressure_pa, rel=2.0e-13)
    assert center.chemical_potential_over_rt == pytest.approx(
        fixed.chemical_potential_over_rt, rel=2.0e-13
    )
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert (
        center.gradient[2].hex(),
        center.hessian[3].hex(),
        center.hessian[8].hex(),
        center.hessian[13].hex(),
    ) == (
        "-0x1.6fdef6ff64731p-3",
        "-0x1.22470c02f3168p-3",
        "0x1.17dc4a225752ep+5",
        "-0x1.247d050916fadp-3",
    )

    step = 2.0e-5
    lower = shifted(-step)
    upper = shifted(step)
    value_direction = sum(center.gradient[index] * direction[index] for index in range(5))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert value_direction == pytest.approx(value_difference, rel=3.0e-8)

    for row in range(5):
        hessian_direction = sum(
            center.hessian[row * 5 + column] * direction[column] for column in range(5)
        )
        gradient_difference = (upper.gradient[row] - lower.gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)

    gas_constant = 8.31446261815324
    assert center.pressure_pa == pytest.approx(
        -gas_constant * temperature_k * center.gradient[1], rel=2.0e-13
    )
    assert center.chemical_potential_over_rt == pytest.approx(center.gradient[0], rel=2.0e-13)


def test_native_sdk_rejects_unsupported_models_and_invalid_states() -> None:
    associating = _model("figiel-2025-reference-electrolytes", "water")
    with pytest.raises(EosError, match="unsupported-native-sdk"):
        native_sdk(associating)

    model = _model("gross-2001-methane-ethane", "methane")
    gross_bundle = _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    gross = gross_bundle.select(("methane",))
    synthetic = _ParameterBundle.from_records(
        bundle_id="synthetic-gross-copy",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=gross_bundle.sources,
        domains=gross.domains,
        components=gross.components,
        singles=tuple(
            record
            for record in gross.records
            if isinstance(record, SingleParameterRecord)
        ),
        models=tuple(
            record
            for record in gross.records
            if isinstance(record, ModelParameterRecord)
        ),
    )
    synthetic_model = Mixture(synthetic.select(("methane",)))
    assert synthetic_model.parameter_fingerprint != model.parameter_fingerprint
    with pytest.raises(EosError, match="unadmitted-parameter-fingerprint"):
        native_sdk(synthetic_model)

    class SpoofedEPCSAFT(Mixture):
        @property
        def parameter_fingerprint(self) -> str:
            return model.parameter_fingerprint

    spoofed_model = SpoofedEPCSAFT(synthetic.select(("methane",)))
    with pytest.raises(EosError, match="expected-exact-Mixture-model"):
        native_sdk(spoofed_model)

    _, table = _table(model)
    for temperature, amount, volume in (
        (96.0, 1.0, 1.0e-3),
        (190.0, 0.0, 1.0e-3),
        (190.0, 1.0, 0.0),
        (math.nan, 1.0, 1.0e-3),
    ):
        result = _evaluate(table, temperature, amount, volume)
        assert result.status != STATUS_OK
        assert result.error


def test_native_sdk_unknown_user_bundle_keeps_reference_domain_closed(
    tmp_path: Path,
) -> None:
    installed = _ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    )
    emitted = tmp_path / "user-electrolyte"
    _user_provided_copy(installed, unknown_applicability=True).to_path(emitted)
    loaded = _ParameterBundle.from_path(emitted)
    model = Mixture(
        loaded.select(("water", "sodium-cation", "chloride-anion"))
    )

    _, phase_table = _electrolyte_table(model)
    phase, gradient, hessian = _evaluate_electrolyte(
        phase_table, 298.15, (0.98, 0.01, 0.01), 1.0e-3
    )
    assert phase.status == STATUS_OK
    assert all(math.isfinite(value) for value in gradient)
    assert all(math.isfinite(value) for value in hessian)

    _, reference_table = _neutral_reference_table(model)
    assert reference_table.neutral_reference_basis_row_count == 2
    reference, _, _, _ = _evaluate_neutral_reference(
        reference_table, model.parameter_fingerprint
    )
    assert reference.status == STATUS_DOMAIN_ERROR
    assert reference.error


def test_native_sdk_mixture_callback_rejects_malformed_buffers_and_scope() -> None:
    capsule, table = _table(_mixture_model())
    assert capsule is not None
    callback = EvaluateMixturePhase(table.evaluate_mixture_phase)
    amounts = (ctypes.c_double * 2)(0.4, 0.6)

    def call(
        result: MixturePhaseBlockResult,
        *,
        amount_pointer: object = amounts,
        amount_count: int = 2,
        temperature_k: float = 220.0,
    ) -> int:
        return callback(
            table.model_context,
            temperature_k,
            amount_pointer,
            amount_count,
            1.25e-4,
            ctypes.byref(result),
        )

    malformed: list[tuple[MixturePhaseBlockResult, dict[str, object], int]] = []
    result, _, _ = _mixture_result(3)
    malformed.append((result, {"amount_pointer": None}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    malformed.append((result, {"amount_count": 1}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.gradient = ctypes.POINTER(ctypes.c_double)()
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.hessian = ctypes.POINTER(ctypes.c_double)()
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.gradient_capacity = 2
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    result.hessian_capacity = 8
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(3)
    malformed.append((result, {"temperature_k": 96.0}, STATUS_DOMAIN_ERROR))

    for result, arguments, expected_status in malformed:
        gradient_pointer = ctypes.cast(result.gradient, ctypes.c_void_p).value
        hessian_pointer = ctypes.cast(result.hessian, ctypes.c_void_p).value
        status = call(result, **arguments)
        assert status == expected_status
        assert result.status == status
        assert result.error
        assert ctypes.cast(result.gradient, ctypes.c_void_p).value == gradient_pointer
        assert ctypes.cast(result.hessian, ctypes.c_void_p).value == hessian_pointer

    result, _, _ = _mixture_result(3)
    result.struct_size -= 1
    assert call(result) == STATUS_ABI_MISMATCH

    result, _, _ = _mixture_result(3)
    nonpositive_amounts = (ctypes.c_double * 2)(0.0, 1.0)
    assert call(result, amount_pointer=nonpositive_amounts) == STATUS_INVALID_ARGUMENT
    assert result.status == STATUS_INVALID_ARGUMENT
    assert result.error

    pure_capsule, pure_table = _table(_model("gross-2001-methane-ethane", "methane"))
    assert pure_capsule is not None
    pure_mixture_result, _, _ = _mixture_result(2)
    pure_amount = (ctypes.c_double * 1)(1.0)
    status = EvaluateMixturePhase(pure_table.evaluate_mixture_phase)(
        pure_table.model_context,
        190.0,
        pure_amount,
        1,
        1.0e-3,
        ctypes.byref(pure_mixture_result),
    )
    assert status == STATUS_UNSUPPORTED_MODEL
    assert pure_mixture_result.status == status
    assert pure_mixture_result.error

    mixture_pure_result = _evaluate(table, 220.0, 1.0, 1.25e-4)
    assert mixture_pure_result.status == STATUS_UNSUPPORTED_MODEL
    assert mixture_pure_result.error

    for segment_count, segment_diameter, dispersion_energy in (
        (0.0, 3.7039, 150.03),
        (1.0, 0.0, 150.03),
        (1.0, 3.7039, 0.0),
        (1.0, 3.7039, math.nan),
    ):
        result = _evaluate_parameters(
            table,
            150.0,
            1.0,
            2.0e-3,
            segment_count,
            segment_diameter,
            dispersion_energy,
        )
        assert result.status != STATUS_OK
        assert result.error


@pytest.mark.parametrize(
    ("table_factory", "callback_type", "callback_field"),
    (
        (_table, EvaluateMixturePhaseKij, "evaluate_mixture_phase_kij"),
        (_lij_table, EvaluateMixturePhaseLij, "evaluate_mixture_phase_lij"),
    ),
)
def test_native_sdk_active_pair_rejects_malformed_buffers_models_and_nonfinite_inputs(
    table_factory,
    callback_type,
    callback_field: str,
) -> None:
    capsule, table = table_factory(_mixture_model())
    assert capsule is not None
    callback = callback_type(getattr(table, callback_field))
    amounts = (ctypes.c_double * 2)(0.4, 0.6)

    def call(
        result: MixturePhaseBlockResult,
        *,
        amount_pointer: object = amounts,
        amount_count: int = 2,
        temperature_k: float = 220.0,
        volume_m3: float = 1.25e-4,
        pair_parameter: float = 0.0,
    ) -> int:
        return callback(
            table.model_context,
            temperature_k,
            amount_pointer,
            amount_count,
            volume_m3,
            pair_parameter,
            ctypes.byref(result),
        )

    malformed: list[tuple[MixturePhaseBlockResult, dict[str, object], int]] = []
    result, _, _ = _mixture_result(3)
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    result.gradient = ctypes.POINTER(ctypes.c_double)()
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    result.hessian_capacity = 15
    malformed.append((result, {}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"amount_pointer": None}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"amount_count": 1}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"temperature_k": math.nan}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append((result, {"volume_m3": math.inf}, STATUS_INVALID_ARGUMENT))
    result, _, _ = _mixture_result(4)
    malformed.append(
        (result, {"pair_parameter": math.nan}, STATUS_INVALID_ARGUMENT)
    )
    result, _, _ = _mixture_result(4)
    malformed.append(
        (result, {"pair_parameter": math.inf}, STATUS_INVALID_ARGUMENT)
    )

    for result, arguments, expected_status in malformed:
        gradient_pointer = ctypes.cast(result.gradient, ctypes.c_void_p).value
        hessian_pointer = ctypes.cast(result.hessian, ctypes.c_void_p).value
        status = call(result, **arguments)
        assert status == expected_status
        assert result.status == status
        assert result.error
        assert ctypes.cast(result.gradient, ctypes.c_void_p).value == gradient_pointer
        assert ctypes.cast(result.hessian, ctypes.c_void_p).value == hessian_pointer

    foreign_size = ctypes.sizeof(ctypes.c_uint32)
    # Guard the four-byte foreign layout so any payload write is observable
    # without invoking undefined behavior in the test process.
    guarded_foreign_result = (ctypes.c_ubyte * ctypes.sizeof(MixturePhaseBlockResult))()
    ctypes.c_uint32.from_buffer(guarded_foreign_result).value = foreign_size
    for index in range(foreign_size, len(guarded_foreign_result)):
        guarded_foreign_result[index] = 0xA5
    canary = bytes(guarded_foreign_result)[foreign_size:]
    status = callback(
        table.model_context,
        220.0,
        amounts,
        2,
        1.25e-4,
        0.0,
        ctypes.cast(guarded_foreign_result, ctypes.POINTER(MixturePhaseBlockResult)),
    )
    assert status == STATUS_ABI_MISMATCH
    assert bytes(guarded_foreign_result)[foreign_size:] == canary

    finite_unbounded, _, _ = _mixture_result(4)
    assert call(finite_unbounded, pair_parameter=-0.2) == STATUS_OK
    assert finite_unbounded.status == STATUS_OK

    pure_capsule, pure_table = table_factory(
        _model("gross-2001-methane-ethane", "methane")
    )
    assert pure_capsule is not None
    pure_result, _, _ = _mixture_result(3)
    pure_amount = (ctypes.c_double * 1)(1.0)
    status = callback_type(getattr(pure_table, callback_field))(
        pure_table.model_context,
        190.0,
        pure_amount,
        1,
        1.0e-3,
        0.0,
        ctypes.byref(pure_result),
    )
    assert status == STATUS_UNSUPPORTED_MODEL
    assert pure_result.status == status
    assert pure_result.error


def test_native_sdk_electrolyte_phase_matches_public_state_and_directional_derivatives() -> None:
    model = _electrolyte_model()
    capsule, table = _electrolyte_table(model)
    assert capsule is not None
    assert table.table_size >= ctypes.sizeof(ElectrolyteNativeSdkTable) == 112
    assert table.component_count == 3
    assert tuple(table.component_ids[index].decode() for index in range(3)) == model.component_ids
    assert tuple(table.component_charges[index] for index in range(3)) == (0, 1, -1)
    assert table.evaluate_electrolyte_phase

    temperature_k = 298.15
    point = (0.98, 0.01, 0.01, 1.0e-3)
    center, gradient, hessian = _evaluate_electrolyte(table, temperature_k, point[:3], point[3])
    assert center.status == STATUS_OK
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    total_amount = math.fsum(point[:3])
    density = total_amount / point[3]
    fractions = tuple(amount / total_amount for amount in point[:3])
    public = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=fractions,
    )
    expected_value = (
        math.fsum(amount * (math.log(amount / point[3]) - 1.0) for amount in point[:3])
        + total_amount * public._residual_helmholtz
    )
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_value, rel=3.0e-12, abs=3.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        public.pressure.to("pascal").magnitude, rel=3.0e-11, abs=2.0e-5
    )
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[-1], rel=3.0e-13
    )

    direction = (0.02, -0.003, -0.003, 2.0e-5)
    step = 2.0e-4

    def shifted(scale: float):
        state = tuple(value + scale * delta for value, delta in zip(point, direction, strict=True))
        return _evaluate_electrolyte(table, temperature_k, state[:3], state[3])

    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    ) == pytest.approx(value_difference, rel=2.0e-7, abs=2.0e-8)
    for row in range(4):
        hessian_direction = math.fsum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=2.0e-6, abs=2.0e-6)

    permuted = _electrolyte_model(("water", "chloride-anion", "sodium-cation"))
    _, permuted_table = _electrolyte_table(permuted)
    permuted_result, permuted_gradient, permuted_hessian = _evaluate_electrolyte(
        permuted_table, temperature_k, (point[0], point[2], point[1]), point[3]
    )
    permutation = (0, 2, 1, 3)
    assert permuted_result.helmholtz_over_rt_reference_amount == pytest.approx(
        center.helmholtz_over_rt_reference_amount, rel=3.0e-13
    )
    assert permuted_result.pressure_pa == pytest.approx(center.pressure_pa, rel=3.0e-13)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * 4 + column] for row in permutation for column in permutation)
    )


def test_native_sdk_held_phase_matches_public_state_and_directional_derivatives() -> None:
    model = _held_water_ionization_model()
    capsule, table = _electrolyte_table(model)
    assert capsule is not None
    assert tuple(table.component_ids[index].decode() for index in range(3)) == model.component_ids
    assert tuple(table.component_charges[index] for index in range(3)) == (0, 1, -1)

    temperature_k = 298.15
    point = (0.98, 0.01, 0.01, 2.5e-5)
    center, gradient, hessian = _evaluate_electrolyte(table, temperature_k, point[:3], point[3])
    assert center.status == STATUS_OK
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    total_amount = math.fsum(point[:3])
    density = total_amount / point[3]
    fractions = tuple(amount / total_amount for amount in point[:3])
    public = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=fractions,
    )
    expected_value = (
        math.fsum(amount * (math.log(amount / point[3]) - 1.0) for amount in point[:3])
        + total_amount * public._residual_helmholtz
    )
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_value, rel=3.0e-12, abs=3.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        public.pressure.to("pascal").magnitude, rel=3.0e-11, abs=2.0e-5
    )
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[-1], rel=3.0e-13
    )

    direction = (0.02, -0.003, -0.003, 5.0e-7)
    step = 2.0e-4

    def shifted(scale: float):
        state = tuple(value + scale * delta for value, delta in zip(point, direction, strict=True))
        return _evaluate_electrolyte(table, temperature_k, state[:3], state[3])

    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    ) == pytest.approx(value_difference, rel=3.0e-7, abs=3.0e-8)
    for row in range(4):
        hessian_direction = math.fsum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-6, abs=3.0e-6)

    permuted = _held_water_ionization_model(
        ("hydroxide-anion", "water", "hydronium-cation")
    )
    _, permuted_table = _electrolyte_table(permuted)
    permutation = (2, 0, 1, 3)
    permuted_result, permuted_gradient, permuted_hessian = _evaluate_electrolyte(
        permuted_table,
        temperature_k,
        tuple(point[index] for index in permutation[:-1]),
        point[3],
    )
    assert permuted_result.helmholtz_over_rt_reference_amount == pytest.approx(
        center.helmholtz_over_rt_reference_amount, rel=3.0e-13
    )
    assert permuted_result.pressure_pa == pytest.approx(center.pressure_pa, rel=3.0e-13)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * 4 + column] for row in permutation for column in permutation)
    )


def test_native_sdk_electrolyte_phase_rejects_abi_domain_and_unsupported_models() -> None:
    capsule, table = _electrolyte_table(_electrolyte_model())
    assert capsule is not None
    callback = EvaluateMixturePhase(table.evaluate_electrolyte_phase)
    amounts = (ctypes.c_double * 3)(0.98, 0.01, 0.01)

    def call(result: MixturePhaseBlockResult, values: object = amounts) -> int:
        return callback(table.model_context, 298.15, values, 3, 1.0e-3, ctypes.byref(result))

    result, _, _ = _mixture_result(4)
    result.struct_size -= 1
    assert call(result) == STATUS_ABI_MISMATCH

    for invalid in (
        (ctypes.c_double * 3)(0.98, 0.01, 0.02),
        (ctypes.c_double * 3)(0.50, 0.25, 0.25),
        (ctypes.c_double * 3)(0.98, 0.0, 0.0),
        (ctypes.c_double * 3)(0.98, math.nan, 0.01),
    ):
        result, _, _ = _mixture_result(4)
        status = call(result, invalid)
        assert status == result.status != STATUS_OK
        assert result.error

    neutral_capsule, neutral_table = _electrolyte_table(_mixture_model())
    assert neutral_capsule is not None
    neutral_result, _, _ = _mixture_result(3)
    neutral_amounts = (ctypes.c_double * 2)(0.4, 0.6)
    status = EvaluateMixturePhase(neutral_table.evaluate_electrolyte_phase)(
        neutral_table.model_context,
        220.0,
        neutral_amounts,
        2,
        1.25e-4,
        ctypes.byref(neutral_result),
    )
    assert status == neutral_result.status == STATUS_UNSUPPORTED_MODEL
    assert neutral_result.error


def test_native_sdk_khudaida_five_component_phase_block_and_permutation() -> None:
    component_ids = (
        "water",
        "ethanol",
        "isobutanol",
        "sodium-cation",
        "chloride-anion",
    )
    model = _khudaida_model(component_ids)
    capsule, table = _electrolyte_table(model)
    assert capsule is not None
    assert table.component_count == 5
    assert tuple(table.component_ids[index].decode() for index in range(5)) == component_ids
    assert tuple(table.component_charges[index] for index in range(5)) == (0, 0, 0, 1, -1)

    temperature_k = 293.15
    point = (0.65905, 0.05095, 0.27015, 0.01985, 0.01985, 1.0e-3)
    center, gradient, hessian = _evaluate_electrolyte(table, temperature_k, point[:5], point[5])
    assert center.status == STATUS_OK
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    total_amount = math.fsum(point[:5])
    density = total_amount / point[5]
    fractions = tuple(amount / total_amount for amount in point[:5])
    public = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=fractions,
    )
    expected_value = (
        math.fsum(amount * (math.log(amount / point[5]) - 1.0) for amount in point[:5])
        + total_amount * public._residual_helmholtz
    )
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_value, rel=4.0e-12, abs=4.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        public.pressure.to("pascal").magnitude, rel=4.0e-11, abs=2.0e-5
    )
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[-1], rel=4.0e-13
    )

    direction = (0.01, -0.004, -0.003, -0.001, -0.001, 2.0e-5)
    step = 2.0e-4

    def shifted(scale: float):
        state = tuple(value + scale * delta for value, delta in zip(point, direction, strict=True))
        return _evaluate_electrolyte(table, temperature_k, state[:5], state[5])

    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    ) == pytest.approx(value_difference, rel=3.0e-7, abs=3.0e-8)
    for row in range(6):
        hessian_direction = math.fsum(
            hessian[row * 6 + column] * direction[column] for column in range(6)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=4.0e-6, abs=4.0e-6)

    permuted_ids = (
        "chloride-anion",
        "isobutanol",
        "water",
        "sodium-cation",
        "ethanol",
    )
    permuted_model = _khudaida_model(permuted_ids)
    permuted_capsule, permuted_table = _electrolyte_table(permuted_model)
    assert permuted_capsule is not None
    permutation = (4, 2, 0, 3, 1, 5)
    permuted_result, permuted_gradient, permuted_hessian = _evaluate_electrolyte(
        permuted_table,
        temperature_k,
        tuple(point[index] for index in permutation[:-1]),
        point[-1],
    )
    assert tuple(permuted_table.component_charges[index] for index in range(5)) == (
        -1,
        0,
        0,
        1,
        0,
    )
    assert permuted_result.helmholtz_over_rt_reference_amount == pytest.approx(
        center.helmholtz_over_rt_reference_amount, rel=5.0e-13
    )
    assert permuted_result.pressure_pa == pytest.approx(center.pressure_pa, rel=5.0e-13)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * 6 + column] for row in permutation for column in permutation)
    )
    with pytest.raises(EosError, match="unadmitted-parameter-fingerprint"):
        native_sdk(_khudaida_model(("water", "isobutanol", "sodium-cation", "chloride-anion")))


def test_native_sdk_ascani_phase_block_parity_derivatives_and_permutation() -> None:
    component_ids = (
        "water",
        "1-butanol",
        "sodium-cation",
        "potassium-cation",
        "chloride-anion",
    )
    model = _ascani_model(component_ids)
    capsule, table = _electrolyte_table(model)
    assert capsule is not None
    assert tuple(table.component_ids[index].decode() for index in range(5)) == component_ids
    assert tuple(table.component_charges[index] for index in range(5)) == (0, 0, 1, 1, -1)

    temperature_k = 298.15
    point = (0.70, 0.20, 0.02, 0.01, 0.03, 1.0e-3)
    center, gradient, hessian = _evaluate_electrolyte(table, temperature_k, point[:5], point[5])
    assert center.status == STATUS_OK
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    total_amount = math.fsum(point[:5])
    density = total_amount / point[5]
    fractions = tuple(amount / total_amount for amount in point[:5])
    public = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=fractions,
    )
    expected_value = (
        math.fsum(amount * (math.log(amount / point[5]) - 1.0) for amount in point[:5])
        + total_amount * public._residual_helmholtz
    )
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_value, rel=5.0e-12, abs=5.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        public.pressure.to("pascal").magnitude, rel=5.0e-11, abs=2.0e-5
    )
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[-1], rel=5.0e-13
    )

    direction = (0.007, -0.003, -0.001, -0.0005, -0.0015, 1.5e-5)
    step = 2.0e-4

    def shifted(scale: float):
        state = tuple(value + scale * delta for value, delta in zip(point, direction, strict=True))
        return _evaluate_electrolyte(table, temperature_k, state[:5], state[5])

    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    directional_gradient = math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    )
    centered_gradient = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert directional_gradient == pytest.approx(centered_gradient, rel=4.0e-7, abs=4.0e-8)
    for row in range(6):
        hessian_direction = math.fsum(
            hessian[row * 6 + column] * direction[column] for column in range(6)
        )
        centered_hessian = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(centered_hessian, rel=6.0e-6, abs=6.0e-6)

    permuted_ids = (
        "chloride-anion",
        "potassium-cation",
        "1-butanol",
        "water",
        "sodium-cation",
    )
    permutation = (4, 3, 1, 0, 2, 5)
    permuted_model = _ascani_model(permuted_ids)
    permuted_capsule, permuted_table = _electrolyte_table(permuted_model)
    assert permuted_capsule is not None
    permuted_result, permuted_gradient, permuted_hessian = _evaluate_electrolyte(
        permuted_table,
        temperature_k,
        tuple(point[index] for index in permutation[:-1]),
        point[-1],
    )
    assert permuted_result.helmholtz_over_rt_reference_amount == pytest.approx(
        center.helmholtz_over_rt_reference_amount, rel=7.0e-13
    )
    assert permuted_result.pressure_pa == pytest.approx(center.pressure_pa, rel=7.0e-13)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * 6 + column] for row in permutation for column in permutation)
    )

    with pytest.raises(EosError, match="unadmitted-parameter-fingerprint"):
        native_sdk(_ascani_model(("water", "1-butanol", "sodium-cation", "chloride-anion")))


def test_native_sdk_molar_volume_bounds_match_packing_oracles_and_permutation() -> None:
    eta_min = 1.0e-6
    eta_max = 0.74

    def oracle(
        temperature_k: float,
        composition: tuple[float, ...],
        segment_count: tuple[float, ...],
        sigma_angstrom: tuple[float, ...],
        energy_over_k: tuple[float, ...],
        charges: tuple[int, ...],
        *,
        charged_diameter_scale: float = 0.88,
    ) -> tuple[float, float]:
        diameters = tuple(
            charged_diameter_scale * sigma
            if charge
            else sigma * (1.0 - 0.12 * math.exp(-3.0 * energy / temperature_k))
            for sigma, energy, charge in zip(sigma_angstrom, energy_over_k, charges, strict=True)
        )
        coefficient = (
            math.pi
            / 6.0
            * 6.02214076e23
            / 1.0e30
            * math.fsum(
                fraction * segments * diameter**3
                for fraction, segments, diameter in zip(
                    composition, segment_count, diameters, strict=True
                )
            )
        )
        return coefficient / eta_max, coefficient / eta_min

    def evaluate(
        model: Mixture,
        temperature_k: float,
        composition: tuple[float, ...],
    ) -> tuple[object, VolumeBoundsNativeSdkTable, tuple[float, float]]:
        capsule, table = _volume_bounds_table(model)
        values = (ctypes.c_double * len(composition))(*composition)
        bounds = (ctypes.c_double * 2)()
        status = EvaluateMolarVolumeBounds(table.evaluate_molar_volume_bounds)(
            table.model_context,
            model.parameter_fingerprint.encode(),
            temperature_k,
            values,
            len(composition),
            eta_min,
            eta_max,
            bounds,
            2,
        )
        assert status == STATUS_OK
        return capsule, table, tuple(bounds)

    canonical_ids = (
        "water",
        "ethanol",
        "isobutanol",
        "sodium-cation",
        "chloride-anion",
    )
    amounts = (0.65905, 0.05095, 0.27015, 0.01985, 0.01985)
    total = math.fsum(amounts)
    composition = tuple(amount / total for amount in amounts)
    model = _khudaida_model(canonical_ids)
    capsule, table, bounds = evaluate(model, 293.15, composition)
    assert capsule is not None
    assert table.table_size >= ctypes.sizeof(VolumeBoundsNativeSdkTable) == 120
    water_sigma = 2.7927 + 10.11 * math.exp(-0.01775 * 293.15) - 1.417 * math.exp(-0.01146 * 293.15)
    expected = oracle(
        293.15,
        composition,
        (1.2047, 2.3827, 2.8440, 1.0, 1.0),
        (water_sigma, 3.1771, 3.5561, 2.8232, 2.7560),
        (353.95, 198.24, 257.13, 230.0, 170.0),
        (0, 0, 0, 1, -1),
    )
    assert bounds == pytest.approx(expected, rel=2.0e-15)
    assert tuple(table.component_ids[index].decode() for index in range(5)) == canonical_ids

    permutation = (4, 2, 0, 3, 1)
    permuted_ids = tuple(canonical_ids[index] for index in permutation)
    permuted_model = _khudaida_model(permuted_ids)
    permuted_capsule, _, permuted_bounds = evaluate(
        permuted_model,
        293.15,
        tuple(composition[index] for index in permutation),
    )
    assert permuted_capsule is not None
    assert permuted_bounds == pytest.approx(bounds, rel=2.0e-15)

    figiel = _electrolyte_model()
    figiel_capsule, _, figiel_bounds = evaluate(figiel, 298.15, (0.98, 0.01, 0.01))
    assert figiel_capsule is not None
    figiel_sigma = (
        2.7927 + 10.11 * math.exp(-0.01775 * 298.15) - 1.417 * math.exp(-0.01146 * 298.15)
    )
    assert figiel_bounds == pytest.approx(
        oracle(
            298.15,
            (0.98, 0.01, 0.01),
            (1.2047, 1.0, 1.0),
            (figiel_sigma, 2.8232, 2.7560),
            (353.95, 230.0, 170.0),
            (0, 1, -1),
        ),
        rel=2.0e-15,
    )

    ascani_ids = (
        "water",
        "1-butanol",
        "sodium-cation",
        "potassium-cation",
        "chloride-anion",
    )
    ascani_composition = (0.72, 0.20, 0.03, 0.01, 0.04)
    ascani = _ascani_model(ascani_ids)
    ascani_capsule, _, ascani_bounds = evaluate(ascani, 298.15, ascani_composition)
    assert ascani_capsule is not None
    ascani_water_sigma = (
        2.7927 + 10.11 * math.exp(-0.01775 * 298.15) - 1.417 * math.exp(-0.01146 * 298.15)
    )
    ascani_geometry = (
        298.15,
        ascani_composition,
        (1.2047, 2.7510, 1.0, 1.0, 1.0),
        (ascani_water_sigma, 3.6139, 2.8232, 3.3417, 2.7560),
        (353.95, 259.59, 230.0, 200.0, 170.0),
        (0, 0, 1, 1, -1),
    )
    assert ascani_bounds == pytest.approx(
        oracle(*ascani_geometry, charged_diameter_scale=1.0),
        rel=2.0e-15,
    )
    assert ascani_bounds != pytest.approx(
        oracle(*ascani_geometry, charged_diameter_scale=0.88),
        rel=1.0e-8,
    )


def test_native_sdk_molar_volume_bounds_reject_abi_fingerprint_and_domain_errors() -> None:
    model = _khudaida_model(("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion"))
    capsule, table = _volume_bounds_table(model)
    assert capsule is not None
    callback = EvaluateMolarVolumeBounds(table.evaluate_molar_volume_bounds)
    composition = (ctypes.c_double * 5)(0.64, 0.08, 0.24, 0.02, 0.02)

    def call(
        *,
        fingerprint: bytes | None = model.parameter_fingerprint.encode(),
        temperature_k: float = 293.15,
        values: object = composition,
        count: int = 5,
        eta_min: float = 1.0e-6,
        eta_max: float = 0.74,
        bounds: object,
        bound_count: int = 2,
    ) -> int:
        return callback(
            table.model_context,
            fingerprint,
            temperature_k,
            values,
            count,
            eta_min,
            eta_max,
            bounds,
            bound_count,
        )

    guarded = (ctypes.c_double * 2)(11.0, 22.0)
    assert call(bounds=guarded, bound_count=1) == STATUS_ABI_MISMATCH
    assert tuple(guarded) == (11.0, 22.0)
    assert call(bounds=guarded, fingerprint=b"sha256:wrong") == STATUS_UNSUPPORTED_MODEL
    assert call(bounds=guarded, temperature_k=math.nan) == STATUS_INVALID_ARGUMENT
    assert call(bounds=guarded, temperature_k=294.0) == STATUS_DOMAIN_ERROR
    assert (
        call(
            bounds=guarded,
            values=(ctypes.c_double * 5)(0.64, 0.08, 0.24, 0.03, 0.01),
        )
        == STATUS_DOMAIN_ERROR
    )
    assert call(bounds=guarded, eta_min=0.0) == STATUS_INVALID_ARGUMENT
    assert call(bounds=guarded, eta_min=0.74, eta_max=1.0e-6) == STATUS_INVALID_ARGUMENT
    assert call(bounds=guarded, eta_max=0.75) == STATUS_DOMAIN_ERROR
    assert tuple(guarded) == (11.0, 22.0)


def test_native_sdk_packing_fraction_matches_independent_derivatives_and_permutation() -> None:
    canonical_ids = (
        "water",
        "ethanol",
        "isobutanol",
        "sodium-cation",
        "chloride-anion",
    )
    temperature_k = 293.15
    amounts = (0.65905, 0.05095, 0.27015, 0.01985, 0.01985)
    volume_m3 = 1.0e-3
    model = _khudaida_model(canonical_ids)

    def evaluate(
        selected_model: Mixture,
        selected_amounts: tuple[float, ...],
    ) -> tuple[object, PackingFractionNativeSdkTable, float, tuple[float, ...], tuple[float, ...]]:
        capsule, table = _packing_fraction_table(selected_model)
        coordinates = len(selected_amounts) + 1
        amount_values = (ctypes.c_double * len(selected_amounts))(*selected_amounts)
        packing_fraction = ctypes.c_double()
        gradient = (ctypes.c_double * coordinates)()
        hessian = (ctypes.c_double * (coordinates * coordinates))()
        status = EvaluatePackingFraction(table.evaluate_packing_fraction)(
            table.model_context,
            selected_model.parameter_fingerprint.encode(),
            temperature_k,
            amount_values,
            len(selected_amounts),
            volume_m3,
            ctypes.byref(packing_fraction),
            gradient,
            coordinates,
            hessian,
            coordinates * coordinates,
        )
        assert status == STATUS_OK
        return capsule, table, packing_fraction.value, tuple(gradient), tuple(hessian)

    capsule, table, value, gradient, hessian = evaluate(model, amounts)
    assert capsule is not None
    assert table.table_size >= ctypes.sizeof(PackingFractionNativeSdkTable) == 128
    assert PackingFractionNativeSdkTable.evaluate_packing_fraction.offset == 120

    water_sigma = (
        2.7927
        + 10.11 * math.exp(-0.01775 * temperature_k)
        - 1.417 * math.exp(-0.01146 * temperature_k)
    )
    diameters = tuple(
        0.88 * sigma if charge else sigma * (1.0 - 0.12 * math.exp(-3.0 * energy / temperature_k))
        for sigma, energy, charge in zip(
            (water_sigma, 3.1771, 3.5561, 2.8232, 2.7560),
            (353.95, 198.24, 257.13, 230.0, 170.0),
            (0, 0, 0, 1, -1),
            strict=True,
        )
    )
    weights = tuple(
        math.pi / 6.0 * 6.02214076e23 / 1.0e30 * segments * diameter**3
        for segments, diameter in zip((1.2047, 2.3827, 2.8440, 1.0, 1.0), diameters, strict=True)
    )
    expected_value = (
        math.fsum(amount * weight for amount, weight in zip(amounts, weights, strict=True))
        / volume_m3
    )
    expected_gradient = tuple(weight / volume_m3 for weight in weights) + (
        -expected_value / volume_m3,
    )
    coordinates = 6
    expected_hessian = [0.0] * (coordinates * coordinates)
    for component, weight in enumerate(weights):
        expected_hessian[component * coordinates + 5] = -weight / volume_m3**2
        expected_hessian[5 * coordinates + component] = -weight / volume_m3**2
    expected_hessian[5 * coordinates + 5] = 2.0 * expected_value / volume_m3**2
    assert value == pytest.approx(expected_value, rel=2.0e-15)
    assert gradient == pytest.approx(expected_gradient, rel=2.0e-15, abs=1.0e-15)
    assert hessian == pytest.approx(expected_hessian, rel=3.0e-15, abs=1.0e-14)

    direction = (0.01, -0.002, 0.004, 0.0005, 0.0005, 2.0e-6)
    directional_gradient = math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    )
    oracle_directional_gradient = math.fsum(
        left * right for left, right in zip(expected_gradient, direction, strict=True)
    )
    assert directional_gradient == pytest.approx(oracle_directional_gradient, rel=2.0e-15)
    for row in range(coordinates):
        actual = math.fsum(
            hessian[row * coordinates + column] * direction[column] for column in range(coordinates)
        )
        expected = math.fsum(
            expected_hessian[row * coordinates + column] * direction[column]
            for column in range(coordinates)
        )
        assert actual == pytest.approx(expected, rel=3.0e-15, abs=1.0e-14)

    permutation = (4, 2, 0, 3, 1, 5)
    permuted_model = _khudaida_model(tuple(canonical_ids[index] for index in permutation[:-1]))
    permuted_capsule, _, permuted_value, permuted_gradient, permuted_hessian = evaluate(
        permuted_model, tuple(amounts[index] for index in permutation[:-1])
    )
    assert permuted_capsule is not None
    assert permuted_value == pytest.approx(value, rel=2.0e-15)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * coordinates + column] for row in permutation for column in permutation)
    )


def test_native_sdk_packing_fraction_fails_closed_for_abi_model_and_domain_errors() -> None:
    model = _khudaida_model(("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion"))
    capsule, table = _packing_fraction_table(model)
    assert capsule is not None
    callback = EvaluatePackingFraction(table.evaluate_packing_fraction)
    amounts = (ctypes.c_double * 5)(0.64, 0.08, 0.24, 0.02, 0.02)
    value = ctypes.c_double(11.0)
    gradient = (ctypes.c_double * 6)(*(12.0 for _ in range(6)))
    hessian = (ctypes.c_double * 36)(*(13.0 for _ in range(36)))

    def call(
        *,
        fingerprint: bytes | None = model.parameter_fingerprint.encode(),
        temperature_k: float = 293.15,
        amount_values: object = amounts,
        amount_count: int = 5,
        volume_m3: float = 1.0e-3,
        gradient_capacity: int = 6,
        hessian_capacity: int = 36,
    ) -> int:
        return callback(
            table.model_context,
            fingerprint,
            temperature_k,
            amount_values,
            amount_count,
            volume_m3,
            ctypes.byref(value),
            gradient,
            gradient_capacity,
            hessian,
            hessian_capacity,
        )

    def assert_untouched() -> None:
        assert value.value == 11.0
        assert tuple(gradient) == (12.0,) * 6
        assert tuple(hessian) == (13.0,) * 36

    assert call(gradient_capacity=5) == STATUS_ABI_MISMATCH
    assert_untouched()
    assert call(hessian_capacity=35) == STATUS_ABI_MISMATCH
    assert_untouched()
    assert call(fingerprint=b"sha256:wrong") == STATUS_UNSUPPORTED_MODEL
    assert_untouched()
    assert call(temperature_k=math.nan) == STATUS_INVALID_ARGUMENT
    assert_untouched()
    assert call(temperature_k=294.0) == STATUS_DOMAIN_ERROR
    assert_untouched()
    assert call(amount_count=4) == STATUS_INVALID_ARGUMENT
    assert_untouched()
    assert call(amount_values=(ctypes.c_double * 5)(0.64, 0.08, 0.24, 0.0, 0.0)) == (
        STATUS_INVALID_ARGUMENT
    )
    assert_untouched()
    assert call(amount_values=(ctypes.c_double * 5)(0.64, 0.08, 0.24, 0.03, 0.01)) == (
        STATUS_DOMAIN_ERROR
    )
    assert_untouched()
    assert call(volume_m3=math.nan) == STATUS_INVALID_ARGUMENT
    assert_untouched()
    assert call(volume_m3=1.0e-6) == STATUS_DOMAIN_ERROR
    assert_untouched()


def test_native_sdk_header_is_installed_package_data() -> None:
    header = (
        __import__("importlib.resources").resources.files(epcsaft)
        / "include"
        / "epcsaft"
        / "native_sdk_v1.h"
    )
    assert header.is_file()


@pytest.mark.parametrize(
    ("active_ion", "counter_ion", "source_diameter", "fingerprint", "source_target"),
    (
        (
            "lithium-cation",
            "chloride-anion",
            2.784,
            "sha256:dbed1c05a92d922369b9b631e89fe10add759e7c8d004ea53bc0a0c90d2d6f4a",
            -486200.0,
        ),
        (
            "sodium-cation",
            "chloride-anion",
            3.445,
            "sha256:d01bc3e64f19293902c00d9cc52060fccc3945ce736f0cd3e2f756be226cf5bf",
            -381100.0,
        ),
        (
            "potassium-cation",
            "chloride-anion",
            4.150,
            "sha256:05f2db21ad51a27e7b7b8fca267e8a38eb663145c980eae6bee7dce02700af5b",
            -309100.0,
        ),
        (
            "chloride-anion",
            "sodium-cation",
            4.100,
            "sha256:a8383fb8f70a1bb8c258f1662c792cd5507af1d23db71ad4575026b5dcab55ea",
            -314900.0,
        ),
        (
            "bromide-anion",
            "sodium-cation",
            4.480,
            "sha256:f313c30e9471daed0638ab70ff43eb033ed828724547c4b32857442ed35f1107",
            -290900.0,
        ),
    ),
)
def test_native_sdk_active_born_matches_public_reference_state(
    active_ion: str,
    counter_ion: str,
    source_diameter: float,
    fingerprint: str,
    source_target: float,
) -> None:
    model = _electrolyte_model(("water", active_ion, counter_ion))
    capsule, table = _active_born_table(model)
    assert capsule is not None
    assert model.parameter_fingerprint == fingerprint
    assert ACTIVE_BORN_PREFIX_SIZE == 168
    assert table.table_size >= AQUEOUS_MIAC_KIJ_PREFIX_SIZE == 184
    assert table.ion_solvation_born_result_size == ctypes.sizeof(IonSolvationBornResult) == 272
    assert table.evaluate_ion_solvation_born

    observed = _evaluate_active_born(table, source_diameter)
    reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    active_position = reference.ion_component_ids.index(active_ion)
    expected = reference.solvation_free_energy[active_position].to("joule / mole").magnitude

    assert observed.status == STATUS_OK
    assert observed.solvation_gibbs_j_per_mol == pytest.approx(expected, abs=2.0e-8)
    # Figiel et al. (2025), DOI 10.1021/acs.iecr.5c00475, SI Table S5.
    # These source targets are evidence only; they are not runtime constants.
    assert observed.solvation_gibbs_j_per_mol < 0.0
    assert abs(observed.solvation_gibbs_j_per_mol - source_target) <= 1.5e3
    assert math.isfinite(observed.derivative_j_per_mol_per_angstrom)
    assert observed.reference_molality_mol_per_kg == 1.0e-12
    assert 0.0 <= observed.reference_convergence_error <= 5.0e-5
    assert observed.parameter_fingerprint.decode() == fingerprint
    assert observed.error == b""


@pytest.mark.parametrize(
    ("active_ion", "counter_ion", "trial_diameter"),
    (
        ("sodium-cation", "chloride-anion", 3.2),
        ("bromide-anion", "sodium-cation", 4.8),
    ),
)
def test_native_sdk_active_born_derivative_matches_centered_difference(
    active_ion: str,
    counter_ion: str,
    trial_diameter: float,
) -> None:
    _, table = _active_born_table(_electrolyte_model(("water", active_ion, counter_ion)))
    center = _evaluate_active_born(table, trial_diameter)

    estimates = []
    for step in (2.0e-4, 1.0e-4):
        upper = _evaluate_active_born(table, trial_diameter + step)
        lower = _evaluate_active_born(table, trial_diameter - step)
        estimates.append(
            (upper.solvation_gibbs_j_per_mol - lower.solvation_gibbs_j_per_mol) / (2.0 * step)
        )
    tolerance = max(
        1.0e-3,
        20.0 * abs(estimates[0] - estimates[1]),
        2.0e-8 * abs(center.derivative_j_per_mol_per_angstrom),
    )
    assert abs(center.derivative_j_per_mol_per_angstrom - estimates[1]) <= tolerance


def test_native_sdk_active_born_partial_trends_to_zero_at_large_diameter() -> None:
    _, table = _active_born_table(_electrolyte_model(("water", "sodium-cation", "chloride-anion")))
    smaller = _evaluate_active_born(table, 16.0)
    larger = _evaluate_active_born(table, 32.0)

    assert abs(larger.derivative_j_per_mol_per_angstrom) < abs(
        smaller.derivative_j_per_mol_per_angstrom
    )


def test_native_sdk_active_born_preserves_prior_prefix_and_fails_closed() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "chloride-anion"))
    _, table = _active_born_table(model)
    for field_name, _ in ApplicabilityNativeSdkTable._fields_:
        assert (
            getattr(ActiveBornNativeSdkTable, field_name).offset
            == getattr(ApplicabilityNativeSdkTable, field_name).offset
        )
    assert ActiveBornNativeSdkTable.ion_solvation_born_result_size.offset == 152
    assert ActiveBornNativeSdkTable.evaluate_ion_solvation_born.offset == 160
    assert ctypes.sizeof(ActiveBornNativeSdkTable) == 168

    callback = EvaluateIonSolvationBorn(table.evaluate_ion_solvation_born)
    undersized = IonSolvationBornResult(struct_size=8, status=123)
    before = bytes(undersized)
    status = callback(
        table.model_context,
        298.15,
        100000.0,
        3.445,
        ctypes.byref(undersized),
    )
    assert status == STATUS_ABI_MISMATCH
    assert bytes(undersized) == before

    for temperature, pressure, diameter, expected_status in (
        (math.nan, 100000.0, 3.445, STATUS_INVALID_ARGUMENT),
        (298.15, 100000.0, 0.0, STATUS_INVALID_ARGUMENT),
        (299.0, 100000.0, 3.445, STATUS_DOMAIN_ERROR),
        (298.15, 101325.0, 3.445, STATUS_DOMAIN_ERROR),
    ):
        result = IonSolvationBornResult(struct_size=ctypes.sizeof(IonSolvationBornResult))
        status = callback(
            table.model_context,
            temperature,
            pressure,
            diameter,
            ctypes.byref(result),
        )
        assert status == result.status == expected_status
        assert math.isnan(result.solvation_gibbs_j_per_mol)
        assert math.isnan(result.derivative_j_per_mol_per_angstrom)
        assert math.isnan(result.reference_molality_mol_per_kg)
        assert math.isnan(result.reference_convergence_error)
        assert result.parameter_fingerprint == b""
        assert result.error

    null_result = IonSolvationBornResult(struct_size=ctypes.sizeof(IonSolvationBornResult))
    status = callback(None, 298.15, 100000.0, 3.445, ctypes.byref(null_result))
    assert status == null_result.status == STATUS_INVALID_ARGUMENT
    assert null_result.error

    neutral_capsule, neutral_table = _active_born_table(_mixture_model())
    assert neutral_capsule is not None
    neutral_result = IonSolvationBornResult(struct_size=ctypes.sizeof(IonSolvationBornResult))
    status = EvaluateIonSolvationBorn(neutral_table.evaluate_ion_solvation_born)(
        neutral_table.model_context,
        298.15,
        100000.0,
        3.445,
        ctypes.byref(neutral_result),
    )
    assert status == neutral_result.status == STATUS_UNSUPPORTED_MODEL
    assert neutral_result.error

    iodide_bundle = _ParameterBundle.from_catalog("figiel-2025-aqueous-iodides", version=1)
    iodide = Mixture(iodide_bundle.select(("water", "sodium-cation", "iodide-anion")))
    _, iodide_table = _active_born_table(iodide)
    iodide_result = _evaluate_active_born(iodide_table, 3.4)
    assert iodide_result.status == STATUS_UNSUPPORTED_MODEL

    lithium = _electrolyte_model(("water", "lithium-cation", "chloride-anion"))
    lithium_capsule, lithium_table = _electrolyte_table(lithium)
    assert lithium_capsule is not None
    phase_result, _, _ = _evaluate_electrolyte(
        lithium_table,
        298.15,
        (1.0, 1.0e-6, 1.0e-6),
        1.8e-5,
    )
    assert phase_result.status == STATUS_UNSUPPORTED_MODEL
    assert phase_result.error


@pytest.mark.parametrize(
    ("active_ion", "counter_ion", "source_diameter"),
    (
        ("sodium-cation", "chloride-anion", 3.445),
        ("bromide-anion", "sodium-cation", 4.480),
    ),
)
def test_native_sdk_ionic_region_permittivity_is_model_bound_and_exact(
    active_ion: str,
    counter_ion: str,
    source_diameter: float,
) -> None:
    model = _electrolyte_model(("water", active_ion, counter_ion))
    _, table = _ion_solvation_ionic_permittivity_table(model)
    descriptors = {
        descriptor.capability: descriptor
        for descriptor in table.capabilities[: table.capability_count]
    }
    descriptor = descriptors[17]
    assert descriptor.parameter_family == 11
    assert descriptor.parameter_identity == 3
    assert descriptor.observation_contract == 2
    assert descriptor.model_domain == 3
    assert descriptor.derivative_order == 1
    assert descriptor.coordinate_count == 1
    assert descriptor.coordinates[0].kind == 13
    assert descriptor.coordinates[0].unit == b"dimensionless"
    assert descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint

    center = _evaluate_ionic_permittivity(table, model.parameter_fingerprint, 8.0)
    source_born = _evaluate_active_born(
        table,
        source_diameter,
    )
    assert center.status == STATUS_OK
    assert center.solvation_gibbs_j_per_mol == pytest.approx(
        source_born.solvation_gibbs_j_per_mol, abs=2.0e-8
    )
    assert center.reference_molality_mol_per_kg == 1.0e-12
    assert 0.0 <= center.reference_convergence_error <= 5.0e-5
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint

    estimates = []
    for step in (2.0e-4, 1.0e-4):
        upper = _evaluate_ionic_permittivity(
            table, model.parameter_fingerprint, 8.0 + step
        )
        lower = _evaluate_ionic_permittivity(
            table, model.parameter_fingerprint, 8.0 - step
        )
        estimates.append(
            (upper.solvation_gibbs_j_per_mol - lower.solvation_gibbs_j_per_mol)
            / (2.0 * step)
        )
    tolerance = max(
        1.0e-3,
        20.0 * abs(estimates[0] - estimates[1]),
        2.0e-8 * abs(center.derivative_j_per_mol),
    )
    assert abs(center.derivative_j_per_mol - estimates[1]) <= tolerance

    wrong_fingerprint = _evaluate_ionic_permittivity(
        table, "sha256:" + "0" * 64, 8.0
    )
    assert wrong_fingerprint.status == STATUS_UNSUPPORTED_MODEL

    for invalid_value, expected_status in (
        (1.0, STATUS_DOMAIN_ERROR),
        (math.nan, STATUS_INVALID_ARGUMENT),
    ):
        invalid = _evaluate_ionic_permittivity(
            table, model.parameter_fingerprint, invalid_value
        )
        assert invalid.status == expected_status
        assert math.isnan(invalid.solvation_gibbs_j_per_mol)
        assert math.isnan(invalid.derivative_j_per_mol)
        assert invalid.error

    neutral = _mixture_model()
    _, neutral_table = _ion_solvation_ionic_permittivity_table(neutral)
    unsupported = _evaluate_ionic_permittivity(
        neutral_table, neutral.parameter_fingerprint, 8.0
    )
    assert unsupported.status == STATUS_UNSUPPORTED_MODEL
    assert unsupported.error


@pytest.mark.parametrize(
    ("active_ion", "counter_ion", "source_diameter"),
    (
        ("sodium-cation", "chloride-anion", 3.445),
        ("bromide-anion", "sodium-cation", 4.480),
    ),
)
def test_native_sdk_solvent_permittivity_is_component_bound_and_exact(
    active_ion: str,
    counter_ion: str,
    source_diameter: float,
) -> None:
    model = _electrolyte_model(("water", active_ion, counter_ion))
    _, table = _ion_solvation_solvent_permittivity_table(model)
    descriptors = {
        descriptor.capability: descriptor
        for descriptor in table.capabilities[: table.capability_count]
    }
    descriptor = descriptors[18]
    assert descriptor.parameter_family == 12
    assert descriptor.parameter_identity == 2
    assert descriptor.observation_contract == 2
    assert descriptor.model_domain == 3
    assert descriptor.derivative_order == 1
    assert descriptor.coordinate_count == 1
    assert descriptor.coordinates[0].kind == 14
    assert descriptor.coordinates[0].component_index == 0
    assert descriptor.coordinates[0].unit == b"dimensionless"
    assert descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint

    center = _evaluate_solvent_permittivity(
        table, model.parameter_fingerprint, 78.09
    )
    source_born = _evaluate_active_born(table, source_diameter)
    assert center.status == STATUS_OK
    assert center.solvation_gibbs_j_per_mol == pytest.approx(
        source_born.solvation_gibbs_j_per_mol, abs=2.0e-8
    )
    assert center.reference_molality_mol_per_kg == 1.0e-12
    assert 0.0 <= center.reference_convergence_error <= 5.0e-5

    estimates = []
    for step in (2.0e-3, 1.0e-3):
        upper = _evaluate_solvent_permittivity(
            table, model.parameter_fingerprint, 78.09 + step
        )
        lower = _evaluate_solvent_permittivity(
            table, model.parameter_fingerprint, 78.09 - step
        )
        estimates.append(
            (upper.solvation_gibbs_j_per_mol - lower.solvation_gibbs_j_per_mol)
            / (2.0 * step)
        )
    tolerance = max(
        1.0e-3,
        20.0 * abs(estimates[0] - estimates[1]),
        2.0e-8 * abs(center.derivative_j_per_mol),
    )
    assert abs(center.derivative_j_per_mol - estimates[1]) <= tolerance

    wrong_fingerprint = _evaluate_solvent_permittivity(
        table, "sha256:" + "0" * 64, 78.09
    )
    assert wrong_fingerprint.status == STATUS_UNSUPPORTED_MODEL

    for invalid_value, expected_status in (
        (1.0, STATUS_DOMAIN_ERROR),
        (math.nan, STATUS_INVALID_ARGUMENT),
    ):
        invalid = _evaluate_solvent_permittivity(
            table, model.parameter_fingerprint, invalid_value
        )
        assert invalid.status == expected_status
        assert math.isnan(invalid.solvation_gibbs_j_per_mol)
        assert math.isnan(invalid.derivative_j_per_mol)
        assert invalid.error

    callback = EvaluateIonSolvationSolventPermittivity(
        table.evaluate_ion_solvation_solvent_permittivity
    )
    assert (
        callback(
            table.model_context,
            model.parameter_fingerprint.encode(),
            298.15,
            100000.0,
            78.09,
            None,
        )
        == STATUS_INVALID_ARGUMENT
    )

    neutral = _mixture_model()
    _, neutral_table = _ion_solvation_solvent_permittivity_table(neutral)
    unsupported = _evaluate_solvent_permittivity(
        neutral_table, neutral.parameter_fingerprint, 78.09
    )
    assert unsupported.status == STATUS_UNSUPPORTED_MODEL
    assert unsupported.error


@pytest.mark.parametrize(
    ("cation", "anion", "molality"),
    (
        *(
            (cation, anion, molality)
            for cation, anion in (
                ("lithium-cation", "chloride-anion"),
                ("sodium-cation", "chloride-anion"),
                ("potassium-cation", "chloride-anion"),
                ("lithium-cation", "bromide-anion"),
                ("sodium-cation", "bromide-anion"),
                ("potassium-cation", "bromide-anion"),
            )
            for molality in (0.001, 1.0, 6.0)
        ),
    ),
)
def test_native_sdk_aqueous_miac_kij_matches_public_reference_state(
    cation: str,
    anion: str,
    molality: float,
) -> None:
    model = _electrolyte_model(("water", cation, anion))
    _, table = _aqueous_miac_kij_table(model)
    catalog = model._model.k_ij
    k_ij = (catalog[1], catalog[2], catalog[5])
    observed = _evaluate_aqueous_miac_kij(
        table, model.parameter_fingerprint, molality, k_ij
    )
    reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    ).state(molality=molality * unit_registry.mole / unit_registry.kilogram)

    assert table.table_size >= AQUEOUS_MIAC_KIJ_PREFIX_SIZE == 184
    assert table.aqueous_miac_kij_result_size == ctypes.sizeof(AqueousMiacKijResult) == 296
    assert observed.status == STATUS_OK
    # The derivative callback projects its value from the certified density
    # midpoint to the requested pressure; this bound covers only that correction.
    assert observed.log_mean_ionic_activity_coefficient_molality == pytest.approx(
        math.log(reference.activity.mean_molality), abs=8.0e-11
    )
    assert all(math.isfinite(value) for value in observed.derivative)
    assert observed.parameter_fingerprint.decode() == model.parameter_fingerprint


def _native_with_aqueous_kij(model: Mixture, k_ij: tuple[float, float, float]) -> NativeModel:
    matrix = list(model._model.k_ij)
    matrix[1] = matrix[3] = k_ij[0]
    matrix[2] = matrix[6] = k_ij[1]
    matrix[5] = matrix[7] = k_ij[2]
    return NativeModel(*_native_arguments(replace(model._model, k_ij=tuple(matrix))))


def _aqueous_composition(model: Mixture, molality: float) -> tuple[float, float, float]:
    solvent_moles = 1.0 / model._model.molar_mass_kg_per_mol[0]
    total_moles = solvent_moles + 2.0 * molality
    return (
        solvent_moles / total_moles,
        molality / total_moles,
        molality / total_moles,
    )


def test_native_sdk_aqueous_miac_kij_uses_canonical_dispersion_density_and_reference_owners() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "chloride-anion"))
    _, table = _aqueous_miac_kij_table(model)
    molality = 1.0
    center = (-0.15, 0.20, 0.35)
    step = 5.0e-5
    callback = _evaluate_aqueous_miac_kij(
        table, model.parameter_fingerprint, molality, center
    )
    composition = _aqueous_composition(model, molality)
    center_native = _native_with_aqueous_kij(model, center)
    center_state = center_native.evaluate_pressure(
        298.15, 100000.0, composition, "liquid"
    )

    pressure_derivative = []
    reference_derivative = []
    fixed_density_derivative = []
    for coordinate in range(3):
        upper = list(center)
        lower = list(center)
        upper[coordinate] += step
        lower[coordinate] -= step
        upper_native = _native_with_aqueous_kij(model, tuple(upper))
        lower_native = _native_with_aqueous_kij(model, tuple(lower))
        upper_pressure = upper_native.evaluate_pressure(
            298.15, 100000.0, composition, "liquid"
        )
        lower_pressure = lower_native.evaluate_pressure(
            298.15, 100000.0, composition, "liquid"
        )
        upper_fixed = upper_native.evaluate(298.15, center_state.molar_density, composition)
        lower_fixed = lower_native.evaluate(298.15, center_state.molar_density, composition)
        upper_reference = upper_native.electrolyte_reference_sequence(
            298.15, 100000.0, (1.0, 0.0, 0.0), (1, 2), (1, 1), "liquid"
        )
        lower_reference = lower_native.electrolyte_reference_sequence(
            298.15, 100000.0, (1.0, 0.0, 0.0), (1, 2), (1, 1), "liquid"
        )
        pressure_derivative.append(
            0.5
            * sum(
                upper_pressure.log_fugacity_coefficient[ion]
                - lower_pressure.log_fugacity_coefficient[ion]
                for ion in (1, 2)
            )
            / (2.0 * step)
        )
        fixed_density_derivative.append(
            0.5
            * sum(
                upper_fixed.log_fugacity_coefficient[ion]
                - lower_fixed.log_fugacity_coefficient[ion]
                for ion in (1, 2)
            )
            / (2.0 * step)
        )
        reference_derivative.append(
            0.5
            * sum(
                upper_reference.log_fugacity[ion] - lower_reference.log_fugacity[ion]
                for ion in (0, 1)
            )
            / (2.0 * step)
        )
        for unchanged in ("hard_chain", "association", "debye_huckel", "born_ssm_ds"):
            assert getattr(upper_fixed, unchanged) == pytest.approx(
                getattr(lower_fixed, unchanged), abs=2.0e-13
            )
        assert upper_fixed.dispersion != lower_fixed.dispersion

    total = [
        finite - reference
        for finite, reference in zip(pressure_derivative, reference_derivative, strict=True)
    ]
    coarse_total = []
    for coordinate in range(3):
        upper = list(center)
        lower = list(center)
        upper[coordinate] += 2.0 * step
        lower[coordinate] -= 2.0 * step
        upper_native = _native_with_aqueous_kij(model, tuple(upper))
        lower_native = _native_with_aqueous_kij(model, tuple(lower))
        upper_pressure = upper_native.evaluate_pressure(
            298.15, 100000.0, composition, "liquid"
        )
        lower_pressure = lower_native.evaluate_pressure(
            298.15, 100000.0, composition, "liquid"
        )
        upper_reference = upper_native.electrolyte_reference_sequence(
            298.15, 100000.0, (1.0, 0.0, 0.0), (1, 2), (1, 1), "liquid"
        )
        lower_reference = lower_native.electrolyte_reference_sequence(
            298.15, 100000.0, (1.0, 0.0, 0.0), (1, 2), (1, 1), "liquid"
        )
        finite = 0.5 * sum(
            upper_pressure.log_fugacity_coefficient[ion]
            - lower_pressure.log_fugacity_coefficient[ion]
            for ion in (1, 2)
        ) / (4.0 * step)
        reference = 0.5 * sum(
            upper_reference.log_fugacity[ion] - lower_reference.log_fugacity[ion]
            for ion in (0, 1)
        ) / (4.0 * step)
        coarse_total.append(finite - reference)
    for exact, fine, coarse in zip(
        callback.derivative, total, coarse_total, strict=True
    ):
        tolerance = max(1.0e-8, 20.0 * abs(coarse - fine), 2.0e-8 * abs(exact))
        assert abs(exact - fine) <= tolerance
    assert max(
        abs(fixed - pressure)
        for fixed, pressure in zip(fixed_density_derivative, pressure_derivative, strict=True)
    ) > 1.0e-7

    malformed_arguments = list(_native_arguments(model._model))
    malformed_arguments[11] = [1.0e308 if value else value for value in malformed_arguments[11]]
    with pytest.raises(RuntimeError, match="association"):
        NativeModel(*malformed_arguments).evaluate(
            298.15, center_state.molar_density, composition
        )


@pytest.mark.parametrize(
    ("cation", "anion", "molality"),
    (
        ("lithium-cation", "chloride-anion", 0.02),
        ("potassium-cation", "bromide-anion", 4.0),
    ),
)
def test_native_sdk_aqueous_miac_kij_exact_derivatives(
    cation: str,
    anion: str,
    molality: float,
) -> None:
    model = _electrolyte_model(("water", cation, anion))
    _, table = _aqueous_miac_kij_table(model)
    trial = (-0.15, 0.20, 0.35)
    center = _evaluate_aqueous_miac_kij(
        table, model.parameter_fingerprint, molality, trial
    )
    for coordinate in range(3):
        estimates = []
        for step in (1.0e-4, 5.0e-5):
            upper = list(trial)
            lower = list(trial)
            upper[coordinate] += step
            lower[coordinate] -= step
            high = _evaluate_aqueous_miac_kij(
                table, model.parameter_fingerprint, molality, tuple(upper)
            )
            low = _evaluate_aqueous_miac_kij(
                table, model.parameter_fingerprint, molality, tuple(lower)
            )
            estimates.append(
                (
                    high.log_mean_ionic_activity_coefficient_molality
                    - low.log_mean_ionic_activity_coefficient_molality
                )
                / (2.0 * step)
            )
        tolerance = max(
            1.0e-8,
            20.0 * abs(estimates[0] - estimates[1]),
            2.0e-8 * abs(center.derivative[coordinate]),
        )
        assert abs(center.derivative[coordinate] - estimates[1]) <= tolerance


@pytest.mark.parametrize(
    ("cation", "anion", "molality", "coordinate"),
    (
        ("sodium-cation", "bromide-anion", 0.6, 0),
        ("potassium-cation", "chloride-anion", 3.0, 1),
        ("potassium-cation", "bromide-anion", 0.01, 1),
        ("lithium-cation", "chloride-anion", 0.02, 2),
    ),
)
def test_native_sdk_aqueous_miac_kij_derivatives_cover_regression_discriminators(
    cation: str,
    anion: str,
    molality: float,
    coordinate: int,
) -> None:
    model = _electrolyte_model(("water", cation, anion))
    _, table = _aqueous_miac_kij_table(model)
    trial = (0.2, 0.2, 0.2)
    center = _evaluate_aqueous_miac_kij(
        table, model.parameter_fingerprint, molality, trial
    )
    estimates = []
    for step in (1.0e-4, 5.0e-5):
        upper = list(trial)
        lower = list(trial)
        upper[coordinate] += step
        lower[coordinate] -= step
        high = _evaluate_aqueous_miac_kij(
            table, model.parameter_fingerprint, molality, tuple(upper)
        )
        low = _evaluate_aqueous_miac_kij(
            table, model.parameter_fingerprint, molality, tuple(lower)
        )
        estimates.append(
            (
                high.log_mean_ionic_activity_coefficient_molality
                - low.log_mean_ionic_activity_coefficient_molality
            )
            / (2.0 * step)
        )
    tolerance = max(
        1.0e-8,
        20.0 * abs(estimates[0] - estimates[1]),
        2.0e-8 * abs(center.derivative[coordinate]),
    )
    assert abs(center.derivative[coordinate] - estimates[1]) <= tolerance


def test_native_sdk_aqueous_miac_kij_preserves_prefix_and_fails_closed() -> None:
    model = _electrolyte_model()
    capsule, table = _aqueous_miac_kij_table(model)
    assert capsule is not None
    for field_name, field_type in ActiveBornNativeSdkTable._fields_:
        assert field_type is not None
        assert (
            getattr(AqueousMiacKijNativeSdkTable, field_name).offset
            == getattr(ActiveBornNativeSdkTable, field_name).offset
        )
    assert AqueousMiacKijNativeSdkTable.aqueous_miac_kij_result_size.offset == 168
    assert AqueousMiacKijNativeSdkTable.evaluate_aqueous_miac_kij.offset == 176
    assert ctypes.sizeof(AqueousMiacKijNativeSdkTable) == 184

    callback = EvaluateAqueousMiacKij(table.evaluate_aqueous_miac_kij)
    coordinates = (ctypes.c_double * 3)(0.0, 0.0, 0.0)
    for fingerprint, temperature, molality, count, expected in (
        (model.parameter_fingerprint.encode(), math.nan, 1.0, 3, STATUS_INVALID_ARGUMENT),
        (model.parameter_fingerprint.encode(), 298.15, 0.0, 3, STATUS_DOMAIN_ERROR),
        (b"sha256:wrong", 298.15, 1.0, 3, STATUS_UNSUPPORTED_MODEL),
        (model.parameter_fingerprint.encode(), 298.15, 1.0, 2, STATUS_ABI_MISMATCH),
    ):
        result = AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult))
        status = callback(
            table.model_context,
            fingerprint,
            temperature,
            100000.0,
            molality,
            coordinates,
            count,
            ctypes.byref(result),
        )
        assert status == result.status == expected, result.error.decode()
        assert math.isnan(result.log_mean_ionic_activity_coefficient_molality)
        assert all(math.isnan(value) for value in result.derivative)
        assert result.parameter_fingerprint == b""
        assert result.error

    with pytest.raises(ValueError, match="unsupported aqueous active k_ij model identity"):
        _native_sdk_v1(
            model._native,
            model.parameter_fingerprint,
            298.15,
            298.15,
            ("water", "potassium-cation", "chloride-anion"),
            0.38,
            False,
            False,
            True,
        )

    lithium = _electrolyte_model(("water", "lithium-cation", "chloride-anion"))
    _, lithium_table = _aqueous_miac_kij_table(lithium)
    for trial, expected in (
        ((100.0, 100.0, 100.0), STATUS_DOMAIN_ERROR),
        ((-100.0, -100.0, -100.0), STATUS_EVALUATION_ERROR),
    ):
        result = _evaluate_aqueous_miac_kij(
            lithium_table, lithium.parameter_fingerprint, 1.0, trial
        )
        assert result.status == expected
        assert math.isnan(result.log_mean_ionic_activity_coefficient_molality)
        assert all(math.isnan(value) for value in result.derivative)
        assert result.parameter_fingerprint == b""
        assert result.error


def test_native_sdk_aqueous_miac_kij_installed_header_layout() -> None:
    compiler = shutil.which("cc")
    assert compiler is not None
    include_root = Path(epcsaft.__file__).parent / "include"
    source = r'''
#include <stddef.h>
#include "epcsaft/native_sdk_v1.h"

_Static_assert(offsetof(epcsaft_native_sdk_v1, aqueous_miac_kij_result_size) == 168, "result-size offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_kij) == 176, "callback offset");
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 184, "table prefix size");
_Static_assert(offsetof(epcsaft_aqueous_miac_kij_result_v1, log_mean_ionic_activity_coefficient_molality) == 8, "value offset");
_Static_assert(offsetof(epcsaft_aqueous_miac_kij_result_v1, derivative) == 16, "derivative offset");
_Static_assert(offsetof(epcsaft_aqueous_miac_kij_result_v1, parameter_fingerprint) == 64, "fingerprint offset");
_Static_assert(offsetof(epcsaft_aqueous_miac_kij_result_v1, error) == 136, "error offset");
_Static_assert(sizeof(epcsaft_aqueous_miac_kij_result_v1) == 296, "result size");
int main(void) { return 0; }
'''
    completed = subprocess.run(
        (compiler, "-std=c11", f"-I{include_root}", "-x", "c", "-fsyntax-only", "-"),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_native_sdk_aqueous_miac_solvation_factor_matches_public_value_and_exact_derivative() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _aqueous_miac_solvation_factor_table(model)
    assert table.aqueous_miac_solvation_factor_result_size == 280

    fixed = _evaluate_aqueous_miac_solvation_factor(
        table, model.parameter_fingerprint, 1.0, 1.5
    )
    reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    public = reference.state(molality=1.0 * unit_registry.mole / unit_registry.kilogram)
    assert fixed.status == STATUS_OK
    assert fixed.log_mean_ionic_activity_coefficient_molality == pytest.approx(
        math.log(public.activity.mean_molality), abs=8.0e-11
    )
    assert fixed.reference_molality_mol_per_kg == pytest.approx(
        reference.reference_molality.to("mole / kilogram").magnitude,
        rel=0.0,
        abs=0.0,
    )
    assert fixed.reference_convergence_error == pytest.approx(
        reference.convergence_error, rel=0.0, abs=0.0
    )
    assert math.isfinite(fixed.reference_derivative_convergence_error)
    assert fixed.parameter_fingerprint.decode() == model.parameter_fingerprint

    molalities = (
        0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1,
        0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4,
        3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
    )
    for solvation_factor in (1.35, 1.65):
        for molality in molalities:
            center = _evaluate_aqueous_miac_solvation_factor(
                table, model.parameter_fingerprint, molality, solvation_factor
            )
            estimates = []
            for step in (1.0e-4, 5.0e-5):
                high = _evaluate_aqueous_miac_solvation_factor(
                    table,
                    model.parameter_fingerprint,
                    molality,
                    solvation_factor + step,
                )
                low = _evaluate_aqueous_miac_solvation_factor(
                    table,
                    model.parameter_fingerprint,
                    molality,
                    solvation_factor - step,
                )
                estimates.append(
                    (
                        high.log_mean_ionic_activity_coefficient_molality
                        - low.log_mean_ionic_activity_coefficient_molality
                    )
                    / (2.0 * step)
                )
            tolerance = max(
                1.0e-8,
                20.0 * abs(estimates[0] - estimates[1]),
                2.0e-8 * abs(center.derivative),
            )
            assert abs(center.derivative - estimates[1]) <= tolerance


def test_native_sdk_aqueous_miac_batches_preserve_scalar_results() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _aqueous_miac_batch_table(model)

    assert AqueousMiacBatchNativeSdkTable.evaluate_aqueous_miac_kij_batch.offset == 224
    assert (
        AqueousMiacBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch.offset
        == 232
    )
    assert ctypes.sizeof(AqueousMiacBatchNativeSdkTable) == 240

    molalities = (0.001, 1.0, 6.0)
    solvation_factor = 1.35
    batch_solvation = _evaluate_aqueous_miac_solvation_factor_batch(
        table,
        model.parameter_fingerprint,
        molalities,
        solvation_factor,
    )
    scalar_solvation = tuple(
        _evaluate_aqueous_miac_solvation_factor(
            table,
            model.parameter_fingerprint,
            molality,
            solvation_factor,
        )
        for molality in molalities
    )
    for batch, scalar in zip(batch_solvation, scalar_solvation, strict=True):
        assert batch.status == scalar.status == STATUS_OK
        assert (
            batch.log_mean_ionic_activity_coefficient_molality
            == scalar.log_mean_ionic_activity_coefficient_molality
        )
        assert batch.derivative == scalar.derivative
        assert (
            batch.reference_molality_mol_per_kg
            == scalar.reference_molality_mol_per_kg
        )
        assert batch.reference_convergence_error == scalar.reference_convergence_error
        assert (
            batch.reference_derivative_convergence_error
            == scalar.reference_derivative_convergence_error
        )
        assert batch.parameter_fingerprint == scalar.parameter_fingerprint
        assert batch.error == scalar.error == b""


def test_native_sdk_aqueous_miac_bounded_batches_cancel_without_background_work() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _aqueous_miac_bounded_batch_table(model)
    assert table.evaluation_budget_size == ctypes.sizeof(NativeEvaluationBudget) == 24
    assert (
        AqueousMiacBoundedBatchNativeSdkTable.evaluation_budget_size.offset
        == 240
    )
    assert (
        AqueousMiacBoundedBatchNativeSdkTable.evaluate_aqueous_miac_kij_batch_bounded.offset
        == 248
    )
    assert (
        AqueousMiacBoundedBatchNativeSdkTable.evaluate_aqueous_miac_solvation_factor_batch_bounded.offset
        == 256
    )
    assert ctypes.sizeof(AqueousMiacBoundedBatchNativeSdkTable) == 264

    molalities = tuple(0.001 + (5.999 * index / 163.0) for index in range(164))
    molality_values = (ctypes.c_double * len(molalities))(*molalities)
    coordinates = (ctypes.c_double * 3)(-0.2, -0.3, 0.4)
    cancelled_results = (AqueousMiacKijResult * len(molalities))(
        *(
            AqueousMiacKijResult(
                struct_size=ctypes.sizeof(AqueousMiacKijResult)
            )
            for _ in molalities
        )
    )
    tiny_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=0,
        deadline_ns=0,
    )
    started = time.perf_counter()
    cancelled_status = EvaluateAqueousMiacKijBatchBounded(
        table.evaluate_aqueous_miac_kij_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        coordinates,
        3,
        cancelled_results,
        len(molalities),
        ctypes.byref(tiny_budget),
    )
    cancellation_elapsed = time.perf_counter() - started
    assert cancelled_status == STATUS_CANCELLED
    assert all(result.status == STATUS_CANCELLED for result in cancelled_results)
    assert all(result.error for result in cancelled_results)
    assert cancellation_elapsed < 5.0

    deadline_results = (AqueousMiacKijResult * 1)(
        AqueousMiacKijResult(struct_size=ctypes.sizeof(AqueousMiacKijResult))
    )
    expired_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=(1 << 64) - 1,
        deadline_ns=1,
    )
    deadline_status = EvaluateAqueousMiacKijBatchBounded(
        table.evaluate_aqueous_miac_kij_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        (ctypes.c_double * 1)(1.0),
        1,
        coordinates,
        3,
        deadline_results,
        1,
        ctypes.byref(expired_budget),
    )
    assert deadline_status == STATUS_DEADLINE_EXCEEDED
    assert deadline_results[0].status == STATUS_DEADLINE_EXCEEDED
    assert deadline_results[0].error

    solvation_cancelled_results = (AqueousMiacSolvationFactorResult * len(molalities))(
        *(
            AqueousMiacSolvationFactorResult(
                struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
            )
            for _ in molalities
        )
    )
    solvation_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=0,
        deadline_ns=0,
    )
    solvation_cancelled_status = EvaluateAqueousMiacSolvationFactorBatchBounded(
        table.evaluate_aqueous_miac_solvation_factor_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        molality_values,
        len(molalities),
        1.35,
        solvation_cancelled_results,
        len(molalities),
        ctypes.byref(solvation_budget),
    )
    assert solvation_cancelled_status == STATUS_CANCELLED
    assert all(
        result.status == STATUS_CANCELLED and result.error
        for result in solvation_cancelled_results
    )

    bounded_results = (AqueousMiacKijResult * 3)(
        *(
            AqueousMiacKijResult(
                struct_size=ctypes.sizeof(AqueousMiacKijResult)
            )
            for _ in range(3)
        )
    )
    sufficient_budget = NativeEvaluationBudget(
        struct_size=ctypes.sizeof(NativeEvaluationBudget),
        work_budget=(1 << 64) - 1,
        deadline_ns=0,
    )
    short_molalities = (ctypes.c_double * 3)(0.001, 1.0, 6.0)
    bounded_status = EvaluateAqueousMiacKijBatchBounded(
        table.evaluate_aqueous_miac_kij_batch_bounded
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        short_molalities,
        3,
        coordinates,
        3,
        bounded_results,
        3,
        ctypes.byref(sufficient_budget),
    )
    existing_results = _evaluate_aqueous_miac_kij_batch(
        ctypes.cast(
            ctypes.pointer(table), ctypes.POINTER(AqueousMiacBatchNativeSdkTable)
        ).contents,
        model.parameter_fingerprint,
        (0.001, 1.0, 6.0),
        (-0.2, -0.3, 0.4),
    )
    assert bounded_status == STATUS_OK
    for bounded, existing in zip(bounded_results, existing_results, strict=True):
        assert bytes(bounded) == bytes(existing)

    coordinates = (-0.2, -0.3, 0.4)
    batch_kij = _evaluate_aqueous_miac_kij_batch(
        table,
        model.parameter_fingerprint,
        molalities,
        coordinates,
    )
    scalar_kij = tuple(
        _evaluate_aqueous_miac_kij(
            table,
            model.parameter_fingerprint,
            molality,
            coordinates,
        )
        for molality in molalities
    )
    for batch, scalar in zip(batch_kij, scalar_kij, strict=True):
        assert batch.status == scalar.status == STATUS_OK
        assert (
            batch.log_mean_ionic_activity_coefficient_molality
            == scalar.log_mean_ionic_activity_coefficient_molality
        )
        assert tuple(batch.derivative) == tuple(scalar.derivative)
        assert (
            batch.reference_molality_mol_per_kg
            == scalar.reference_molality_mol_per_kg
        )
        assert batch.reference_convergence_error == scalar.reference_convergence_error
        assert (
            batch.reference_derivative_convergence_error
            == scalar.reference_derivative_convergence_error
        )
        assert batch.parameter_fingerprint == scalar.parameter_fingerprint
        assert batch.error == scalar.error == b""


def test_native_sdk_aqueous_miac_batches_fail_closed() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _aqueous_miac_batch_table(model)
    molalities = (ctypes.c_double * 2)(0.001, 6.0)
    malformed = (AqueousMiacSolvationFactorResult * 2)(
        AqueousMiacSolvationFactorResult(
            struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult) - 1,
            status=123,
        ),
        AqueousMiacSolvationFactorResult(
            struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult),
            status=456,
        ),
    )
    malformed_before = bytes(malformed)
    malformed_status = EvaluateAqueousMiacSolvationFactorBatch(
        table.evaluate_aqueous_miac_solvation_factor_batch
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        molalities,
        2,
        1.35,
        malformed,
        2,
    )
    assert malformed_status == STATUS_ABI_MISMATCH
    assert bytes(malformed) == malformed_before

    solvation_results = (AqueousMiacSolvationFactorResult * 2)(
        *(
            AqueousMiacSolvationFactorResult(
                struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
            )
            for _ in range(2)
        )
    )
    solvation_status = EvaluateAqueousMiacSolvationFactorBatch(
        table.evaluate_aqueous_miac_solvation_factor_batch
    )(
        table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        molalities,
        1,
        1.35,
        solvation_results,
        2,
    )
    assert solvation_status == STATUS_ABI_MISMATCH
    assert all(result.status == solvation_status for result in solvation_results)
    assert all(
        math.isnan(result.log_mean_ionic_activity_coefficient_molality)
        for result in solvation_results
    )
    assert all(result.error for result in solvation_results)

    coordinates = (ctypes.c_double * 3)(-0.2, -0.3, 0.4)
    kij_results = (AqueousMiacKijResult * 2)(
        *(
            AqueousMiacKijResult(
                struct_size=ctypes.sizeof(AqueousMiacKijResult)
            )
            for _ in range(2)
        )
    )
    kij_status = EvaluateAqueousMiacKijBatch(
        table.evaluate_aqueous_miac_kij_batch
    )(
        table.model_context,
        b"sha256:wrong",
        298.15,
        100000.0,
        molalities,
        2,
        coordinates,
        3,
        kij_results,
        2,
    )
    assert kij_status == STATUS_UNSUPPORTED_MODEL
    assert all(result.status == kij_status for result in kij_results)
    assert all(
        math.isnan(result.log_mean_ionic_activity_coefficient_molality)
        and all(math.isnan(value) for value in result.derivative)
        for result in kij_results
    )
    assert all(result.error for result in kij_results)


def test_native_sdk_aqueous_miac_batch_serializes_other_cppad_callbacks() -> None:
    electrolyte = _electrolyte_model(
        ("water", "sodium-cation", "bromide-anion")
    )
    batch_capsule, batch_table = _aqueous_miac_batch_table(electrolyte)
    pure_capsule, pure_table = _table(
        _model("gross-2001-methane-ethane", "methane")
    )
    assert batch_capsule is not None
    assert pure_capsule is not None

    with ThreadPoolExecutor(max_workers=2) as executor:
        batch_future = executor.submit(
            _evaluate_aqueous_miac_kij_batch,
            batch_table,
            electrolyte.parameter_fingerprint,
            (0.001, 1.0, 6.0),
            (-0.2, -0.3, 0.4),
        )
        phase_future = executor.submit(
            _evaluate,
            pure_table,
            150.0,
            1.0,
            1.0e-3,
        )
        batch = batch_future.result(timeout=60.0)
        phase = phase_future.result(timeout=60.0)

    assert all(result.status == STATUS_OK for result in batch)
    assert phase.status == STATUS_OK
    assert math.isfinite(phase.helmholtz_over_rt_reference_amount)


def test_native_sdk_staged_figiel_user_bundles_are_narrowly_admitted() -> None:
    nacl = _figiel_trial_model(("water", "sodium-cation", "chloride-anion"))
    _, nacl_table = _aqueous_miac_solvation_factor_table(nacl)
    assert nacl_table.evaluate_ion_solvation_born
    assert nacl_table.evaluate_aqueous_miac_kij
    assert nacl_table.evaluate_aqueous_miac_solvation_factor
    assert _evaluate_active_born(nacl_table, 3.4).status == STATUS_OK
    assert _evaluate_aqueous_miac_kij(
        nacl_table, nacl.parameter_fingerprint, 1.0, (-0.2, -0.2, 0.2)
    ).status == STATUS_OK

    nabr = _figiel_trial_model(("water", "sodium-cation", "bromide-anion"))
    _, nabr_table = _aqueous_miac_solvation_factor_table(nabr)
    assert _evaluate_aqueous_miac_solvation_factor(
        nabr_table, nabr.parameter_fingerprint, 1.0, 1.4
    ).status == STATUS_OK
    assert _evaluate_aqueous_miac_kij(
        nabr_table, nabr.parameter_fingerprint, 1.0, (-0.2, -0.2, 0.2)
    ).status == STATUS_OK

    for active_anion, counterion in (
        ("chloride-anion", "sodium-cation"),
        ("bromide-anion", "sodium-cation"),
    ):
        anion_model = _figiel_trial_model(("water", active_anion, counterion))
        _, anion_table = _aqueous_miac_solvation_factor_table(anion_model)
        assert _evaluate_active_born(anion_table, 4.2).status == STATUS_OK
        assert _evaluate_aqueous_miac_kij(
            anion_table,
            anion_model.parameter_fingerprint,
            1.0,
            (-0.2, -0.2, 0.2),
        ).status == STATUS_UNSUPPORTED_MODEL
        assert _evaluate_aqueous_miac_solvation_factor(
            anion_table,
            anion_model.parameter_fingerprint,
            1.0,
            1.4,
        ).status == STATUS_UNSUPPORTED_MODEL

    fixture = _figiel_trial_model(
        ("water", "sodium-cation", "bromide-anion"),
        purpose="package-test-fixture",
    )
    _, fixture_table = _capability_table(fixture)
    assert {
        descriptor.capability
        for descriptor in fixture_table.capabilities[
            : fixture_table.capability_count
        ]
    } == {11}
    _, fixture_active_table = _aqueous_miac_solvation_factor_table(fixture)
    assert (
        _evaluate_active_born(fixture_active_table, 3.4).status
        == STATUS_UNSUPPORTED_MODEL
    )
    assert (
        _evaluate_aqueous_miac_kij(
            fixture_active_table,
            fixture.parameter_fingerprint,
            1.0,
            (-0.2, -0.2, 0.2),
        )
        .status
        == STATUS_UNSUPPORTED_MODEL
    )
    provenance_renamed = _figiel_trial_model(
        ("water", "sodium-cation", "bromide-anion"),
        spoof_source_metadata=True,
    )
    _, provenance_renamed_table = _aqueous_miac_solvation_factor_table(
        provenance_renamed
    )
    assert _evaluate_aqueous_miac_solvation_factor(
        provenance_renamed_table,
        provenance_renamed.parameter_fingerprint,
        1.0,
        1.4,
    ).status == STATUS_OK
    reordered = _figiel_trial_model(("sodium-cation", "water", "bromide-anion"))
    _, reordered_table = _aqueous_miac_solvation_factor_table(reordered)
    assert _evaluate_aqueous_miac_solvation_factor(
        reordered_table, reordered.parameter_fingerprint, 1.0, 1.4
    ).status == STATUS_UNSUPPORTED_MODEL

    catalog = _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    unknown = Mixture(
        _user_provided_copy(catalog, unknown_applicability=True).select(
            ("water", "sodium-cation", "bromide-anion")
        )
    )
    _, unknown_table = _aqueous_miac_solvation_factor_table(unknown)
    assert _evaluate_aqueous_miac_solvation_factor(
        unknown_table, unknown.parameter_fingerprint, 1.0, 1.4
    ).status == STATUS_UNSUPPORTED_MODEL


def test_native_sdk_aqueous_miac_solvation_factor_preserves_abi_and_fails_closed() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    capsule, table = _aqueous_miac_solvation_factor_table(model)
    assert capsule is not None
    for field_name, field_type in NeutralReferenceNativeSdkTable._fields_:
        assert field_type is not None
        assert (
            getattr(AqueousMiacSolvationFactorNativeSdkTable, field_name).offset
            == getattr(NeutralReferenceNativeSdkTable, field_name).offset
        )
    assert AqueousMiacSolvationFactorNativeSdkTable.aqueous_miac_solvation_factor_result_size.offset == 208
    assert AqueousMiacSolvationFactorNativeSdkTable.evaluate_aqueous_miac_solvation_factor.offset == 216
    assert ctypes.sizeof(AqueousMiacSolvationFactorNativeSdkTable) == 224

    callback = EvaluateAqueousMiacSolvationFactor(
        table.evaluate_aqueous_miac_solvation_factor
    )
    for fingerprint, temperature, molality, solvation_factor, expected in (
        (model.parameter_fingerprint.encode(), math.nan, 1.0, 1.5, STATUS_INVALID_ARGUMENT),
        (model.parameter_fingerprint.encode(), 298.15, 0.0, 1.5, STATUS_DOMAIN_ERROR),
        (model.parameter_fingerprint.encode(), 298.15, 1.0, math.nan, STATUS_INVALID_ARGUMENT),
        (b"sha256:wrong", 298.15, 1.0, 1.5, STATUS_UNSUPPORTED_MODEL),
    ):
        result = AqueousMiacSolvationFactorResult(
            struct_size=ctypes.sizeof(AqueousMiacSolvationFactorResult)
        )
        status = callback(
            table.model_context,
            fingerprint,
            temperature,
            100000.0,
            molality,
            solvation_factor,
            ctypes.byref(result),
        )
        assert status == result.status == expected
        assert math.isnan(result.log_mean_ionic_activity_coefficient_molality)
        assert math.isnan(result.derivative)
        assert result.parameter_fingerprint == b""
        assert result.error


def test_native_sdk_aqueous_miac_solvation_factor_installed_header_layout() -> None:
    compiler = shutil.which("cc")
    assert compiler is not None
    include_root = Path(epcsaft.__file__).parent / "include"
    source = r'''
#include <stddef.h>
#include "epcsaft/native_sdk_v1.h"

_Static_assert(offsetof(epcsaft_native_sdk_v1, aqueous_miac_solvation_factor_result_size) == 208, "result-size offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_solvation_factor) == 216, "callback offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_kij_batch) == 224, "kij batch callback offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_solvation_factor_batch) == 232, "solvation-factor batch callback offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluation_budget_size) == 240, "budget-size offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_kij_batch_bounded) == 248, "bounded kij callback offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_aqueous_miac_solvation_factor_batch_bounded) == 256, "bounded solvation-factor callback offset");
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 264, "table prefix size");
_Static_assert(sizeof(epcsaft_native_evaluation_budget_v1) == 24, "budget size");
_Static_assert(offsetof(epcsaft_aqueous_miac_solvation_factor_result_v1, log_mean_ionic_activity_coefficient_molality) == 8, "value offset");
_Static_assert(offsetof(epcsaft_aqueous_miac_solvation_factor_result_v1, derivative) == 16, "derivative offset");
_Static_assert(offsetof(epcsaft_aqueous_miac_solvation_factor_result_v1, parameter_fingerprint) == 48, "fingerprint offset");
_Static_assert(offsetof(epcsaft_aqueous_miac_solvation_factor_result_v1, error) == 120, "error offset");
_Static_assert(sizeof(epcsaft_aqueous_miac_solvation_factor_result_v1) == 280, "result size");
int main(void) { return 0; }
'''
    completed = subprocess.run(
        (compiler, "-std=c11", f"-I{include_root}", "-x", "c", "-fsyntax-only", "-"),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_native_sdk_neutral_reference_is_complete_and_matches_public_owner() -> None:
    model = _electrolyte_model()
    capsule, table = _neutral_reference_table(model)
    assert capsule is not None
    assert table.table_size >= NEUTRAL_REFERENCE_PREFIX_SIZE == 208
    assert table.neutral_reference_basis_row_count == 2
    assert table.neutral_reference_result_size == ctypes.sizeof(NeutralReferenceResult) == 440

    result, basis, contractions, reference_composition = _evaluate_neutral_reference(
        table, model.parameter_fingerprint
    )
    public_reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    pure_water = model.state(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    assert pure_water._log_fugacity_coefficient is not None

    assert result.status == STATUS_OK
    assert tuple(table.component_ids[index].decode() for index in range(3)) == model.component_ids
    assert tuple(table.component_charges[index] for index in range(3)) == (0, 1, -1)
    assert basis == (1.0, 0.0, 0.0, 0.0, 1.0, 1.0)
    assert reference_composition == (1.0, 0.0, 0.0)
    assert contractions[0] == pytest.approx(pure_water._log_fugacity_coefficient[0], abs=2e-12)
    assert contractions[1] == pytest.approx(
        math.fsum(public_reference.infinite_dilution_log_fugacity_coefficient),
        abs=2e-12,
    )
    assert result.temperature_k == 298.15
    assert result.pressure_pa == 100000.0
    assert model._model.pressure_min_pa == result.pressure_pa
    assert model._model.pressure_max_pa == result.pressure_pa
    assert result.solvent_molar_mass_kg_per_mol == model._model.molar_mass_kg_per_mol[0]
    assert result.reference_amount_mol == 1.0
    assert result.reference_number_density_mol_per_m3 == 1.0
    assert result.reference_molality_mol_per_kg == 1.0e-12
    assert result.reference_convergence_error == public_reference.convergence_error
    assert result.derivative_availability == 0
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert result.helmholtz_basis_id.decode() == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert result.error == b""

    # Two independent charge-neutral rows span the same space under a basis change.
    reactions = ((-2.0, 1.0, 1.0), (1.0, 1.0, 1.0))
    basis_coefficients = ((-2.0, 1.0), (1.0, 1.0))
    reconstructed = tuple(
        tuple(
            math.fsum(coefficient * basis[row * 3 + component] for row, coefficient in enumerate(weights))
            for component in range(3)
        )
        for weights in basis_coefficients
    )
    assert reconstructed == reactions
    reaction_contractions = tuple(
        math.fsum(coefficient * contractions[row] for row, coefficient in enumerate(weights))
        for weights in basis_coefficients
    )
    transformed_rows = (
        tuple(left + right for left, right in zip(*reactions, strict=True)),
        reactions[1],
    )
    transformed_values = (reaction_contractions[0] + reaction_contractions[1], reaction_contractions[1])
    assert transformed_rows[0] == (-1.0, 2.0, 2.0)
    assert transformed_values[0] == pytest.approx(
        -contractions[0] + 2.0 * contractions[1], abs=2e-15
    )
    # For (water, H3O+, OH-), this is the accepted water self-ionization offset.
    assert reaction_contractions[0] == pytest.approx(contractions[1] - 2.0 * contractions[0])

    for cation in ("lithium-cation", "sodium-cation", "potassium-cation"):
        for anion in ("chloride-anion", "bromide-anion"):
            supported = _electrolyte_model(("water", cation, anion))
            _, supported_table = _neutral_reference_table(supported)
            assert supported_table.neutral_reference_basis_row_count == 2
            supported_result, *_ = _evaluate_neutral_reference(
                supported_table, supported.parameter_fingerprint
            )
            assert supported_result.status == STATUS_OK


def test_native_sdk_neutral_reference_spans_multicomponent_user_model() -> None:
    component_ids = (
        "water",
        "1-butanol",
        "sodium-cation",
        "potassium-cation",
        "chloride-anion",
    )
    installed = _ParameterBundle.from_catalog(
        "ascani-2022-case-study-2-electrolyte-lle", version=1
    )
    models = (
        Mixture(installed.select(component_ids)),
        Mixture(_user_provided_copy(installed).select(component_ids)),
    )
    expected_basis = (
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0, 1.0,
    )

    for model in models:
        _, table = _neutral_reference_table(model)
        assert table.neutral_reference_basis_row_count == 4
        result, basis, contractions, reference_composition = _evaluate_neutral_reference(
            table, model.parameter_fingerprint
        )
        pure_water = model._native.evaluate_pressure(
            298.15, 100000.0, (1.0, 0.0, 0.0, 0.0, 0.0), "liquid"
        )
        sodium_chloride = model._native.electrolyte_reference_sequence(
            298.15,
            100000.0,
            (1.0, 0.0, 0.0, 0.0, 0.0),
            (2, 4),
            (1, 1),
            "liquid",
        )
        potassium_chloride = model._native.electrolyte_reference_sequence(
            298.15,
            100000.0,
            (1.0, 0.0, 0.0, 0.0, 0.0),
            (3, 4),
            (1, 1),
            "liquid",
        )

        assert result.status == STATUS_OK
        assert tuple(table.component_charges[index] for index in range(5)) == (0, 0, 1, 1, -1)
        assert basis == expected_basis
        assert reference_composition == (1.0, 0.0, 0.0, 0.0, 0.0)
        assert contractions == pytest.approx(
            (
                pure_water.log_fugacity_coefficient[0],
                pure_water.log_fugacity_coefficient[1],
                math.fsum(sodium_chloride.log_fugacity),
                math.fsum(potassium_chloride.log_fugacity),
            ),
            abs=3.0e-12,
        )
        assert result.reference_convergence_error == pytest.approx(
            max(
                sodium_chloride.convergence_error,
                potassium_chloride.convergence_error,
            ),
            abs=3.0e-12,
        )


def test_native_sdk_neutral_reference_is_permutation_invariant_and_fails_closed() -> None:
    canonical = _electrolyte_model()
    _, canonical_table = _neutral_reference_table(canonical)
    canonical_result, _, canonical_contractions, _ = _evaluate_neutral_reference(
        canonical_table, canonical.parameter_fingerprint
    )

    permuted = _electrolyte_model(("chloride-anion", "water", "sodium-cation"))
    _, table = _neutral_reference_table(permuted)
    result, basis, contractions, reference_composition = _evaluate_neutral_reference(
        table, permuted.parameter_fingerprint
    )
    assert tuple(table.component_ids[index].decode() for index in range(3)) == permuted.component_ids
    assert tuple(table.component_charges[index] for index in range(3)) == (-1, 0, 1)
    assert basis == (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    assert reference_composition == (0.0, 1.0, 0.0)
    assert contractions == pytest.approx(canonical_contractions, abs=3e-12)
    assert result.reference_convergence_error == pytest.approx(
        canonical_result.reference_convergence_error, abs=5e-12
    )
    assert result.parameter_fingerprint.decode() == permuted.parameter_fingerprint

    callback = EvaluateNeutralReference(table.evaluate_neutral_reference)
    basis_buffer = (ctypes.c_double * 6)(*([123.0] * 6))
    contraction_buffer = (ctypes.c_double * 2)(123.0, 123.0)
    composition_buffer = (ctypes.c_double * 3)(123.0, 123.0, 123.0)

    def call(
        *, fingerprint: bytes = permuted.parameter_fingerprint.encode(), temperature: float = 298.15,
        pressure: float = 100000.0, basis_capacity: int = 6,
    ) -> NeutralReferenceResult:
        failure = NeutralReferenceResult(
            struct_size=ctypes.sizeof(NeutralReferenceResult),
            component_count=3,
            neutral_basis_row_count=2,
            neutral_basis_capacity=basis_capacity,
            contraction_capacity=2,
            reference_composition_capacity=3,
            neutral_basis=basis_buffer,
            log_fugacity_contractions=contraction_buffer,
            reference_composition=composition_buffer,
        )
        status = callback(
            table.model_context,
            fingerprint,
            temperature,
            pressure,
            ctypes.byref(failure),
        )
        assert status == failure.status
        return failure

    for failure, expected in (
        (call(fingerprint=b"sha256:wrong"), STATUS_UNSUPPORTED_MODEL),
        (call(temperature=299.0), STATUS_DOMAIN_ERROR),
        (call(pressure=101325.0), STATUS_DOMAIN_ERROR),
        (call(basis_capacity=5), STATUS_INVALID_ARGUMENT),
    ):
        assert failure.status == expected
        assert failure.parameter_fingerprint == b""
        assert failure.helmholtz_basis_id == b""
        assert failure.error
        assert tuple(basis_buffer) == (123.0,) * 6
        assert tuple(contraction_buffer) == (123.0,) * 2
        assert tuple(composition_buffer) == (123.0,) * 3

    undersized = NeutralReferenceResult(struct_size=8, status=123)
    before = bytes(undersized)
    status = callback(
        table.model_context,
        permuted.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        ctypes.byref(undersized),
    )
    assert status == STATUS_ABI_MISMATCH
    assert bytes(undersized) == before

    unsupported = _mixture_model()
    _, unsupported_table = _neutral_reference_table(unsupported)
    assert unsupported_table.neutral_reference_basis_row_count == 0
    unsupported_result = NeutralReferenceResult(
        struct_size=ctypes.sizeof(NeutralReferenceResult)
    )
    unsupported_status = EvaluateNeutralReference(
        unsupported_table.evaluate_neutral_reference
    )(
        unsupported_table.model_context,
        unsupported.parameter_fingerprint.encode(),
        220.0,
        100000.0,
        ctypes.byref(unsupported_result),
    )
    assert unsupported_status == unsupported_result.status == STATUS_UNSUPPORTED_MODEL
    assert unsupported_result.error


def test_native_sdk_neutral_reference_installed_header_layout() -> None:
    compiler = shutil.which("cc")
    assert compiler is not None
    include_root = Path(epcsaft.__file__).parent / "include"
    source = r'''
#include <stddef.h>
#include "epcsaft/native_sdk_v1.h"

_Static_assert(offsetof(epcsaft_native_sdk_v1, neutral_reference_basis_row_count) == 184, "row-count offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, neutral_reference_result_size) == 192, "result-size offset");
_Static_assert(offsetof(epcsaft_native_sdk_v1, evaluate_neutral_reference) == 200, "callback offset");
_Static_assert(sizeof(epcsaft_native_sdk_v1) >= 208, "table prefix size");
_Static_assert(offsetof(epcsaft_neutral_reference_result_v1, neutral_basis) == 48, "basis offset");
_Static_assert(offsetof(epcsaft_neutral_reference_result_v1, derivative_availability) == 72, "derivative offset");
_Static_assert(offsetof(epcsaft_neutral_reference_result_v1, parameter_fingerprint) == 136, "fingerprint offset");
_Static_assert(offsetof(epcsaft_neutral_reference_result_v1, error) == 280, "error offset");
_Static_assert(sizeof(epcsaft_neutral_reference_result_v1) == 440, "result size");
int main(void) { return 0; }
'''
    completed = subprocess.run(
        (compiler, "-std=c11", f"-I{include_root}", "-x", "c", "-fsyntax-only", "-"),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_native_sdk_internal_reference_and_density_failure_contracts(tmp_path: Path) -> None:
    compiler = shutil.which("c++")
    assert compiler is not None
    native_root = Path(__file__).parents[1] / "src" / "epcsaft" / "_native"
    source = r'''
#include "density.hpp"
#include "electrolyte_reference.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

struct State { double pressure; double pressure_density_derivative; };

int main() {
    bool rejected_nonfinite_reference = false;
    int call = 0;
    try {
        epcsaft::detail::electrolyte_reference_sequence(
            std::vector<double>{1.0, 0.0, 0.0},
            std::vector<double>{0.018, 0.007, 0.035},
            std::array<std::size_t, 2>{1, 2},
            std::array<int, 2>{1, 1},
            [&call](const std::vector<double>&) {
                ++call;
                std::vector<double> derivative(6, 1.0);
                if (call == 3) {
                    derivative[4] = std::numeric_limits<double>::quiet_NaN();
                }
                return epcsaft::detail::ReferenceObservation{{0.0, 0.0}, derivative};
            }
        );
    } catch (const std::runtime_error&) {
        rejected_nonfinite_reference = true;
    }
    if (!rejected_nonfinite_reference) return 1;

    bool domain_classified = false;
    try {
        epcsaft::detail::solve_density(
            1.0, 1.0, 10.0, "liquid",
            [](double density) { return State{density, 1.0}; }
        );
    } catch (const epcsaft::detail::DensityDomainError&) {
        domain_classified = true;
    }
    if (!domain_classified) return 2;

    bool evaluation_classified = false;
    try {
        epcsaft::detail::refine_density_root(
            1.0, 0.5, 1.5, -0.5,
            [](double) {
                return State{std::numeric_limits<double>::quiet_NaN(), 1.0};
            }
        );
    } catch (const epcsaft::detail::DensityDomainError&) {
        return 3;
    } catch (const std::runtime_error&) {
        evaluation_classified = true;
    }
    if (!evaluation_classified) return 4;

    const double adjacent_lower = 59456.796534681715;
    const double adjacent_upper = std::nextafter(
        adjacent_lower,
        std::numeric_limits<double>::infinity()
    );
    const auto adjacent_root = epcsaft::detail::refine_density_root(
        1.0, adjacent_lower, adjacent_upper, -2.0e-6,
        [adjacent_lower](double density) {
            return State{
                density <= adjacent_lower ? 1.0 - 2.0e-6 : 1.0 + 2.0e-6,
                1.0,
            };
        }
    );
    if (!(adjacent_root.molar_density >= adjacent_lower
          && adjacent_root.molar_density <= adjacent_upper)) return 5;
    return 0;
}
'''
    output = tmp_path / "owner-contract"
    compiled = subprocess.run(
        (compiler, "-std=c++17", f"-I{native_root}", "-x", "c++", "-", "-o", str(output)),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert compiled.returncode == 0, compiled.stderr
    executed = subprocess.run((str(output),), capture_output=True, text=True, check=False)
    assert executed.returncode == 0, executed.stderr
