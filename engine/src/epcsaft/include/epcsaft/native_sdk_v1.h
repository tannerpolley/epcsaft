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
    EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1 = 5
};

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
} epcsaft_native_sdk_v1;

#ifdef __cplusplus
}
#endif

#endif
