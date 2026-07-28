#ifndef EPCSAFT_NATIVE_SDK_V1_H
#define EPCSAFT_NATIVE_SDK_V1_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME "epcsaft.native_sdk.v1"
#define EPCSAFT_NATIVE_SDK_V1_ABI_VERSION 1u
#define EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE 72u
#define EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE 160u
#define EPCSAFT_NATIVE_SDK_V1_HELMHOLTZ_BASIS_ID \
    "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"

enum epcsaft_native_status_v1 {
    EPCSAFT_NATIVE_STATUS_OK_V1 = 0,
    EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1 = 1,
    EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1 = 2,
    EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1 = 3,
    EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1 = 4,
    EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1 = 5,
    EPCSAFT_NATIVE_STATUS_CANCELLED_V1 = 6,
    EPCSAFT_NATIVE_STATUS_DEADLINE_EXCEEDED_V1 = 7
};

#define EPCSAFT_NATIVE_EVALUATION_BUDGET_UNLIMITED_V1 UINT64_MAX

typedef struct epcsaft_native_evaluation_budget_v1 {
    uint32_t struct_size;
    uint64_t work_budget;
    uint64_t deadline_ns;
} epcsaft_native_evaluation_budget_v1;

typedef struct epcsaft_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double helmholtz_over_rt_reference_amount;
    double gradient[2];
    double hessian[4];
    double third[8];
    double pressure_pa;
    double chemical_potential_over_rt;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_phase_block_result_v1;

enum epcsaft_parameterized_phase_coordinate_v1 {
    EPCSAFT_PHASE_AMOUNT_MOL_V1 = 0,
    EPCSAFT_PHASE_VOLUME_M3_V1 = 1,
    EPCSAFT_PHASE_SEGMENT_COUNT_V1 = 2,
    EPCSAFT_PHASE_SEGMENT_DIAMETER_ANGSTROM_V1 = 3,
    EPCSAFT_PHASE_DISPERSION_ENERGY_OVER_K_V1 = 4,
    EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1 = 5
};

typedef struct epcsaft_parameterized_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double helmholtz_over_rt_reference_amount;
    double gradient[EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1];
    double hessian[
        EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1 *
        EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1
    ];
    double pressure_pa;
    double chemical_potential_over_rt;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_parameterized_phase_block_result_v1;

typedef int (*epcsaft_evaluate_pure_phase_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    epcsaft_phase_block_result_v1* result
);

typedef int (*epcsaft_evaluate_pure_phase_parameters_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k,
    epcsaft_parameterized_phase_block_result_v1* result
);

enum epcsaft_associating_parameterized_phase_coordinate_v1 {
    EPCSAFT_ASSOCIATING_PHASE_AMOUNT_MOL_V1 = 0,
    EPCSAFT_ASSOCIATING_PHASE_VOLUME_M3_V1 = 1,
    EPCSAFT_ASSOCIATING_PHASE_SEGMENT_COUNT_V1 = 2,
    EPCSAFT_ASSOCIATING_PHASE_SEGMENT_DIAMETER_ANGSTROM_V1 = 3,
    EPCSAFT_ASSOCIATING_PHASE_DISPERSION_ENERGY_OVER_K_V1 = 4,
    EPCSAFT_ASSOCIATING_PHASE_ENERGY_OVER_K_V1 = 5,
    EPCSAFT_ASSOCIATING_PHASE_VOLUME_V1 = 6,
    EPCSAFT_ASSOCIATING_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1 = 7
};

typedef struct epcsaft_associating_parameterized_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double helmholtz_over_rt_reference_amount;
    double gradient[
        EPCSAFT_ASSOCIATING_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1
    ];
    double hessian[
        EPCSAFT_ASSOCIATING_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1 *
        EPCSAFT_ASSOCIATING_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1
    ];
    double pressure_pa;
    double chemical_potential_over_rt;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char topology_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_associating_parameterized_phase_block_result_v1;

/*
 * Source-model variant used by Baygi--Pahlavanzadeh 2015 Eq. (4): the
 * association-strength volume factor is d_ij^3 rather than sigma_ij^3.
 * The ordinary Provider EOS and scalar association callbacks retain their
 * source-owned sigma_ij^3 convention.
 */
typedef int (*epcsaft_evaluate_diameter_basis_associating_pure_phase_parameters_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k,
    double association_energy_over_k,
    double association_volume,
    epcsaft_associating_parameterized_phase_block_result_v1* result
);

/*
 * A mixture result struct-size mismatch is return-code authoritative. The
 * callbacks return ABI_MISMATCH without writing an untrusted result layout.
 * After exact-size acceptance, failures also set result status and error.
 */
typedef struct epcsaft_mixture_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    size_t coordinate_count;
    size_t gradient_capacity;
    size_t hessian_capacity;
    double* gradient;
    double* hessian;
    double helmholtz_over_rt_reference_amount;
    double pressure_pa;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_mixture_phase_block_result_v1;

typedef int (*epcsaft_evaluate_mixture_phase_v1)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    epcsaft_mixture_phase_block_result_v1* result
);

typedef int (*epcsaft_evaluate_mixture_phase_kij_v1)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    double k_ij,
    epcsaft_mixture_phase_block_result_v1* result
);

typedef int (*epcsaft_evaluate_mixture_phase_lij_v1)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    double l_ij,
    epcsaft_mixture_phase_block_result_v1* result
);

/*
 * Evaluates one active pure-component parameter while all other resolved
 * parameters remain fixed. parameter_family must be one of the scalar pure
 * families advertised by the bound model. Coordinates are
 * (amount_mol, volume_m3, active_parameter).
 */
typedef int (*epcsaft_evaluate_pure_phase_parameter_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    uint32_t parameter_family,
    double parameter_value,
    epcsaft_mixture_phase_block_result_v1* result
);

#define EPCSAFT_NATIVE_MOLAR_VOLUME_BOUND_COUNT_V1 2u

enum epcsaft_native_molar_volume_bound_v1 {
    EPCSAFT_NATIVE_MOLAR_VOLUME_LOWER_V1 = 0,
    EPCSAFT_NATIVE_MOLAR_VOLUME_UPPER_V1 = 1
};

/*
 * Converts caller-owned packing-fraction bounds to molar-volume bounds using
 * the immutable model's temperature-dependent PC-SAFT packing geometry.
 * The output order is [lower, upper] in m^3/mol. A bound-count mismatch is
 * return-code authoritative and leaves the caller buffer untouched.
 */
typedef int (*epcsaft_evaluate_molar_volume_bounds_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    const double* mole_fractions,
    size_t component_count,
    double packing_fraction_min,
    double packing_fraction_max,
    double* molar_volume_bounds_m3_per_mol,
    size_t bound_count
);

/*
 * Evaluates the model packing fraction and exact derivatives in ordered
 * coordinates (n_1, ..., n_C, V). Gradient and Hessian capacities must equal
 * C + 1 and (C + 1)^2. Any failure leaves all output buffers untouched.
 */
typedef int (*epcsaft_evaluate_packing_fraction_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    double* packing_fraction,
    double* gradient,
    size_t gradient_capacity,
    double* hessian,
    size_t hessian_capacity
);

/*
 * Model-bound Figiel water reference tracer with one active ion Born diameter.
 * The value is Delta_solv,x G_i^inf in J/mol and the derivative is its total
 * fixed-(T,P,reference-path) derivative in J/(mol angstrom). A result
 * struct-size mismatch is return-code authoritative and leaves the untrusted
 * result layout untouched. Exact-size failures also set status and error.
 */
typedef struct epcsaft_ion_solvation_born_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double solvation_gibbs_j_per_mol;
    double derivative_j_per_mol_per_angstrom;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_ion_solvation_born_result_v1;

typedef int (*epcsaft_evaluate_ion_solvation_born_v1)(
    const void* model_context,
    double temperature_k,
    double pressure_pa,
    double born_diameter_angstrom,
    epcsaft_ion_solvation_born_result_v1* result
);

#define EPCSAFT_AQUEOUS_ACTIVE_KIJ_COORDINATE_COUNT_V1 3u

enum epcsaft_aqueous_active_kij_coordinate_v1 {
    EPCSAFT_AQUEOUS_KIJ_WATER_CATION_V1 = 0,
    EPCSAFT_AQUEOUS_KIJ_WATER_ANION_V1 = 1,
    EPCSAFT_AQUEOUS_KIJ_CATION_ANION_V1 = 2
};

typedef struct epcsaft_aqueous_miac_kij_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double log_mean_ionic_activity_coefficient_molality;
    double derivative[EPCSAFT_AQUEOUS_ACTIVE_KIJ_COORDINATE_COUNT_V1];
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    double reference_derivative_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_aqueous_miac_kij_result_v1;

typedef int (*epcsaft_evaluate_aqueous_miac_kij_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    double formula_unit_molality_mol_per_kg,
    const double* k_ij,
    size_t k_ij_count,
    epcsaft_aqueous_miac_kij_result_v1* result
);

#define EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID_V1 \
    "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"

enum epcsaft_neutral_reference_derivative_v1 {
    EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1 = 0u,
    EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_TEMPERATURE_V1 = 1u << 0,
    EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PRESSURE_V1 = 1u << 1,
    EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_COMPOSITION_V1 = 1u << 2,
    EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PARAMETERS_V1 = 1u << 3
};

/*
 * Complete charge-neutral standard-reference basis in the model's component
 * order. neutral_basis is row-major [basis_row][component], and each
 * dimensionless contraction is sum_i basis[row,i] * ln(phi_i^reference).
 * reference_composition is the salt-free solvent composition used by the
 * canonical four-molality limit. No predictive single-ion reference is
 * exposed. An admitted model has exactly one water reference carrier and at
 * least one cation and anion. Its C-1 rows are every neutral component unit
 * row followed by a minimal cation/anion spanning tree. Neutral components use
 * their infinite-dilution fugacity in the reported pure-water composition;
 * ionic rows use separate canonical formula-unit limits. For
 * (water,H3O+,OH-), water self-ionization is row -2*b_0+b_1.
 *
 * All capacities and dimensions must exactly match the immutable table. A
 * struct-size mismatch is return-code authoritative and leaves the untrusted
 * result layout untouched. Other failures set status/error and leave the
 * caller-owned buffers untouched. derivative_availability is an explicit bit
 * mask; v1 returns NONE and no T, P, composition, or parameter derivatives.
 */
typedef struct epcsaft_neutral_reference_result_v1 {
    uint32_t struct_size;
    int32_t status;
    size_t component_count;
    size_t neutral_basis_row_count;
    size_t neutral_basis_capacity;
    size_t contraction_capacity;
    size_t reference_composition_capacity;
    double* neutral_basis;
    double* log_fugacity_contractions;
    double* reference_composition;
    uint32_t derivative_availability;
    double temperature_k;
    double pressure_pa;
    double solvent_molar_mass_kg_per_mol;
    double reference_amount_mol;
    double reference_number_density_mol_per_m3;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char helmholtz_basis_id[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_neutral_reference_result_v1;

typedef int (*epcsaft_evaluate_neutral_reference_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    epcsaft_neutral_reference_result_v1* result
);

/*
 * Model-bound NaBr molality-scale mean ionic activity with one active water
 * solvation factor. The derivative is total at fixed temperature and pressure
 * through the canonical density and electrolyte reference sequence.
 */
typedef struct epcsaft_aqueous_miac_solvation_factor_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double log_mean_ionic_activity_coefficient_molality;
    double derivative;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    double reference_derivative_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_aqueous_miac_solvation_factor_result_v1;

typedef int (*epcsaft_evaluate_aqueous_miac_solvation_factor_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    double formula_unit_molality_mol_per_kg,
    double solvation_factor,
    epcsaft_aqueous_miac_solvation_factor_result_v1* result
);

/*
 * Batch forms of the existing active aqueous MIAC callbacks. Molalities and
 * results are caller-owned arrays with identical nonzero counts. Provider
 * evaluates the active model and derivative reference sequence once, then
 * returns one ordinary scalar result record per molality in input order. The
 * exact catalog solvation-factor coordinate retains its existing value-parity
 * reference sequence once per batch.
 */
typedef int (*epcsaft_evaluate_aqueous_miac_kij_batch_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    const double* formula_unit_molalities_mol_per_kg,
    size_t molality_count,
    const double* k_ij,
    size_t k_ij_count,
    epcsaft_aqueous_miac_kij_result_v1* results,
    size_t result_count
);

typedef int (*epcsaft_evaluate_aqueous_miac_solvation_factor_batch_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    const double* formula_unit_molalities_mol_per_kg,
    size_t molality_count,
    double solvation_factor,
    epcsaft_aqueous_miac_solvation_factor_result_v1* results,
    size_t result_count
);

/*
 * Cooperative bounded forms of the two batch callbacks. The 240-byte batch
 * prefix remains unchanged; these callbacks are available only in the
 * appended tail. `work_budget` counts safe reference/row checkpoints and
 * UINT64_MAX disables that limit. `deadline_ns` is an absolute monotonic
 * clock deadline; zero disables it. Provider joins all row workers before
 * returning CANCELLED or DEADLINE_EXCEEDED.
 */
typedef int (*epcsaft_evaluate_aqueous_miac_kij_batch_bounded_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    const double* formula_unit_molalities_mol_per_kg,
    size_t molality_count,
    const double* k_ij,
    size_t k_ij_count,
    epcsaft_aqueous_miac_kij_result_v1* results,
    size_t result_count,
    const epcsaft_native_evaluation_budget_v1* budget
);

typedef int (*epcsaft_evaluate_aqueous_miac_solvation_factor_batch_bounded_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    const double* formula_unit_molalities_mol_per_kg,
    size_t molality_count,
    double solvation_factor,
    epcsaft_aqueous_miac_solvation_factor_result_v1* results,
    size_t result_count,
    const epcsaft_native_evaluation_budget_v1* budget
);

typedef struct epcsaft_figiel_permittivity_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double bulk_relative_permittivity;
    double derivative;
    double salt_free_relative_permittivity;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_figiel_permittivity_result_v1;

typedef int (*epcsaft_evaluate_figiel_permittivity_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    double total_ion_mole_fraction,
    double dielectric_ion_suppression_coefficient,
    epcsaft_figiel_permittivity_result_v1* result
);

/*
 * Infinite-dilution ion-solvation Gibbs energy and exact fixed-pressure
 * derivatives with respect to the three pair coordinates in exact component
 * order: solvent-cation, solvent-anion, cation-anion.
 */
typedef struct epcsaft_ion_solvation_kij_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double solvation_gibbs_j_per_mol;
    double derivative[3];
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_ion_solvation_kij_result_v1;

typedef int (*epcsaft_evaluate_ion_solvation_kij_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    size_t active_ion_component_index,
    const double* k_ij,
    size_t k_ij_count,
    epcsaft_ion_solvation_kij_result_v1* result
);

/*
 * Model-bound Figiel water reference tracer with one active SSM+DS ionic-
 * region relative permittivity. The value is Delta_solv,x G_i^inf in J/mol
 * and the derivative is its total fixed-(T,P,reference-path) derivative in
 * J/mol per dimensionless permittivity unit.
 */
typedef struct epcsaft_ion_solvation_ionic_permittivity_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double solvation_gibbs_j_per_mol;
    double derivative_j_per_mol;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_ion_solvation_ionic_permittivity_result_v1;

typedef int (*epcsaft_evaluate_ion_solvation_ionic_permittivity_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    double ionic_region_relative_permittivity,
    epcsaft_ion_solvation_ionic_permittivity_result_v1* result
);

/*
 * Model-bound water reference tracer with one active solvent relative
 * permittivity. The value is Delta_solv,x G_i^inf in J/mol and the derivative
 * is its total fixed-(T,P,reference-path) derivative in J/mol per
 * dimensionless permittivity unit.
 */
typedef struct epcsaft_ion_solvation_solvent_permittivity_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double solvation_gibbs_j_per_mol;
    double derivative_j_per_mol;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_ion_solvation_solvent_permittivity_result_v1;

typedef int (*epcsaft_evaluate_ion_solvation_solvent_permittivity_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    double solvent_relative_permittivity,
    epcsaft_ion_solvation_solvent_permittivity_result_v1* result
);

#define EPCSAFT_NATIVE_CAPABILITY_SCHEMA_VERSION_V1 1u

enum epcsaft_native_capability_v1 {
    EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_KIJ_HELMHOLTZ_V1 = 1u,
    EPCSAFT_NATIVE_CAPABILITY_NEUTRAL_BINARY_LIJ_HELMHOLTZ_V1 = 2u,
    EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_COUNT_HELMHOLTZ_V1 = 3u,
    EPCSAFT_NATIVE_CAPABILITY_PURE_SEGMENT_DIAMETER_HELMHOLTZ_V1 = 4u,
    EPCSAFT_NATIVE_CAPABILITY_PURE_DISPERSION_ENERGY_HELMHOLTZ_V1 = 5u,
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_BORN_V1 = 6u,
    EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_SOLVATION_FACTOR_MIAC_V1 = 7u,
    EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_WATER_CATION_KIJ_MIAC_V1 = 8u,
    EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_WATER_ANION_KIJ_MIAC_V1 = 9u,
    EPCSAFT_NATIVE_CAPABILITY_AQUEOUS_CATION_ANION_KIJ_MIAC_V1 = 10u,
    EPCSAFT_NATIVE_CAPABILITY_FIGIEL_DIELECTRIC_SUPPRESSION_V1 = 11u,
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_CATION_KIJ_V1 = 12u,
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_ANION_KIJ_V1 = 13u,
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_CATION_ANION_KIJ_V1 = 14u,
    EPCSAFT_NATIVE_CAPABILITY_PURE_ASSOCIATION_ENERGY_HELMHOLTZ_V1 = 15u,
    EPCSAFT_NATIVE_CAPABILITY_PURE_ASSOCIATION_VOLUME_HELMHOLTZ_V1 = 16u,
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_IONIC_REGION_PERMITTIVITY_V1 = 17u,
    EPCSAFT_NATIVE_CAPABILITY_ION_SOLVATION_SOLVENT_PERMITTIVITY_V1 = 18u
};

enum epcsaft_native_parameter_family_v1 {
    EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_KIJ_V1 = 1u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_BINARY_INTERACTION_LIJ_V1 = 2u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_COUNT_V1 = 3u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1 = 4u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1 = 5u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_BORN_DIAMETER_V1 = 6u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_SOLVATION_FACTOR_V1 = 7u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_DIELECTRIC_ION_SUPPRESSION_V1 = 8u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_ENERGY_OVER_K_V1 = 9u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_ASSOCIATION_VOLUME_V1 = 10u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_IONIC_REGION_RELATIVE_PERMITTIVITY_V1 = 11u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_RELATIVE_PERMITTIVITY_V1 = 12u
};

enum epcsaft_native_parameter_identity_v1 {
    EPCSAFT_NATIVE_PARAMETER_IDENTITY_UNORDERED_COMPONENT_PAIR_V1 = 1u,
    EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1 = 2u,
    EPCSAFT_NATIVE_PARAMETER_IDENTITY_MODEL_V1 = 3u
};

enum epcsaft_native_observation_contract_v1 {
    EPCSAFT_NATIVE_OBSERVATION_FIXED_COMPOSITION_HELMHOLTZ_PHASE_V1 = 1u,
    EPCSAFT_NATIVE_OBSERVATION_ION_SOLVATION_GIBBS_V1 = 2u,
    EPCSAFT_NATIVE_OBSERVATION_AQUEOUS_MEAN_IONIC_ACTIVITY_V1 = 3u,
    EPCSAFT_NATIVE_OBSERVATION_RELATIVE_PERMITTIVITY_RATIO_V1 = 4u
};

enum epcsaft_native_model_domain_v1 {
    EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_BINARY_V1 = 1u,
    EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_NONASSOCIATING_PURE_V1 = 2u,
    EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_WATER_SINGLE_ION_V1 = 3u,
    EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_AQUEOUS_NABR_V1 = 4u,
    EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_DIELECTRIC_V1 = 5u,
    EPCSAFT_NATIVE_MODEL_DOMAIN_FIGIEL_SINGLE_ION_SOLVATION_V1 = 6u,
    EPCSAFT_NATIVE_MODEL_DOMAIN_NEUTRAL_ASSOCIATING_PURE_V1 = 7u
};

enum epcsaft_native_tensor_layout_v1 {
    EPCSAFT_NATIVE_TENSOR_LAYOUT_ROW_MAJOR_V1 = 1u
};

enum epcsaft_native_capability_maturity_v1 {
    EPCSAFT_NATIVE_CAPABILITY_DERIVATIVE_READY_V1 = 1u
};

enum epcsaft_native_authority_effect_v1 {
    EPCSAFT_NATIVE_AUTHORITY_EFFECT_NONE_V1 = 0u
};

enum epcsaft_native_capability_coordinate_kind_v1 {
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_AMOUNT_V1 = 1u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_VOLUME_V1 = 2u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_KIJ_V1 = 3u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BINARY_INTERACTION_LIJ_V1 = 4u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_COUNT_V1 = 5u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SEGMENT_DIAMETER_V1 = 6u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DISPERSION_ENERGY_OVER_K_V1 = 7u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_BORN_DIAMETER_V1 = 8u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_SOLVATION_FACTOR_V1 = 9u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_DIELECTRIC_ION_SUPPRESSION_V1 = 10u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_ASSOCIATION_ENERGY_OVER_K_V1 = 11u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_ASSOCIATION_VOLUME_V1 = 12u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_IONIC_REGION_RELATIVE_PERMITTIVITY_V1 = 13u,
    EPCSAFT_NATIVE_CAPABILITY_COORDINATE_RELATIVE_PERMITTIVITY_V1 = 14u
};

typedef struct epcsaft_native_capability_coordinate_v1 {
    uint32_t struct_size;
    uint32_t kind;
    int32_t component_index;
    int32_t pair_component_index_a;
    int32_t pair_component_index_b;
    const char* unit;
} epcsaft_native_capability_coordinate_v1;

/*
 * Immutable model-bound metadata for an exact derivative callback already
 * present in this table. This is transport evidence only: maturity says which
 * derivative payload exists, while authority_effect remains NONE. Component
 * indices refer to the table's exact component order. The unordered pair
 * identity is therefore independent of that order even though amount
 * coordinates are ordered.
 */
typedef struct epcsaft_native_capability_descriptor_v1 {
    uint32_t struct_size;
    uint32_t schema_version;
    uint32_t capability;
    uint32_t parameter_family;
    uint32_t parameter_identity;
    uint32_t observation_contract;
    uint32_t model_domain;
    uint32_t tensor_layout;
    uint32_t derivative_order;
    uint32_t maturity;
    uint32_t authority_effect;
    int32_t unsupported_status;
    int32_t domain_status;
    size_t state_coordinate_count;
    size_t active_parameter_count;
    size_t coordinate_count;
    const epcsaft_native_capability_coordinate_v1* coordinates;
    size_t component_count;
    const char* const* component_ids;
    double temperature_min_k;
    double temperature_max_k;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char topology_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char helmholtz_basis_id[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
} epcsaft_native_capability_descriptor_v1;

typedef struct epcsaft_native_sdk_v1 {
    uint32_t abi_version;
    size_t table_size;
    size_t result_size;
    const void* model_context;
    epcsaft_evaluate_pure_phase_v1 evaluate_pure_phase;
    size_t parameterized_result_size;
    epcsaft_evaluate_pure_phase_parameters_v1 evaluate_pure_phase_parameters;
    size_t component_count;
    size_t mixture_result_size;
    epcsaft_evaluate_mixture_phase_v1 evaluate_mixture_phase;
    epcsaft_evaluate_mixture_phase_kij_v1 evaluate_mixture_phase_kij;
    const char* const* component_ids;
    const int32_t* component_charges;
    epcsaft_evaluate_mixture_phase_v1 evaluate_electrolyte_phase;
    epcsaft_evaluate_molar_volume_bounds_v1 evaluate_molar_volume_bounds;
    epcsaft_evaluate_packing_fraction_v1 evaluate_packing_fraction;
    /* Immutable limits resolved from the bound parameter-source domains. */
    double source_temperature_min_k;
    double source_temperature_max_k;
    /* NaN means that the bound source records declare no ionic-fraction cap. */
    double total_ion_mole_fraction_max;
    size_t ion_solvation_born_result_size;
    epcsaft_evaluate_ion_solvation_born_v1 evaluate_ion_solvation_born;
    size_t aqueous_miac_kij_result_size;
    epcsaft_evaluate_aqueous_miac_kij_v1 evaluate_aqueous_miac_kij;
    size_t neutral_reference_basis_row_count;
    size_t neutral_reference_result_size;
    epcsaft_evaluate_neutral_reference_v1 evaluate_neutral_reference;
    size_t aqueous_miac_solvation_factor_result_size;
    epcsaft_evaluate_aqueous_miac_solvation_factor_v1
        evaluate_aqueous_miac_solvation_factor;
    epcsaft_evaluate_aqueous_miac_kij_batch_v1
        evaluate_aqueous_miac_kij_batch;
    epcsaft_evaluate_aqueous_miac_solvation_factor_batch_v1
        evaluate_aqueous_miac_solvation_factor_batch;
    size_t evaluation_budget_size;
    epcsaft_evaluate_aqueous_miac_kij_batch_bounded_v1
        evaluate_aqueous_miac_kij_batch_bounded;
    epcsaft_evaluate_aqueous_miac_solvation_factor_batch_bounded_v1
        evaluate_aqueous_miac_solvation_factor_batch_bounded;
    size_t capability_count;
    const epcsaft_native_capability_descriptor_v1* capabilities;
    epcsaft_evaluate_mixture_phase_lij_v1 evaluate_mixture_phase_lij;
    epcsaft_evaluate_pure_phase_parameter_v1 evaluate_pure_phase_parameter;
    size_t figiel_permittivity_result_size;
    epcsaft_evaluate_figiel_permittivity_v1 evaluate_figiel_permittivity;
    size_t ion_solvation_kij_result_size;
    epcsaft_evaluate_ion_solvation_kij_v1 evaluate_ion_solvation_kij;
    size_t ion_solvation_ionic_permittivity_result_size;
    epcsaft_evaluate_ion_solvation_ionic_permittivity_v1
        evaluate_ion_solvation_ionic_permittivity;
    size_t ion_solvation_solvent_permittivity_result_size;
    epcsaft_evaluate_ion_solvation_solvent_permittivity_v1
        evaluate_ion_solvation_solvent_permittivity;
    size_t associating_parameterized_result_size;
    epcsaft_evaluate_diameter_basis_associating_pure_phase_parameters_v1
        evaluate_diameter_basis_associating_pure_phase_parameters;
} epcsaft_native_sdk_v1;

#ifdef __cplusplus
}
#endif

#endif
