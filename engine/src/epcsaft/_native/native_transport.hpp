#pragma once

#include "eos.hpp"

#include <cstddef>
#include <cstdint>

inline constexpr std::size_t EPCSAFT_NATIVE_FINGERPRINT_SIZE = 72u;
inline constexpr std::size_t EPCSAFT_NATIVE_ERROR_SIZE = 160u;
inline constexpr const char* EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME =
    "epcsaft.native_transport";
inline constexpr const char* EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID =
    "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3";

enum epcsaft_native_status {
    EPCSAFT_NATIVE_STATUS_OK = 0,
    EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT = 1,
    EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL = 2,
    EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR = 3,
    EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR = 4,
    EPCSAFT_NATIVE_STATUS_BUFFER_ERROR = 5,
};

struct epcsaft_phase_block_result final {
    std::uint32_t struct_size;
    std::int32_t status;
    double helmholtz_over_rt_reference_amount;
    double gradient[2];
    double hessian[4];
    double third[8];
    double pressure_pa;
    double chemical_potential_over_rt;
    char parameter_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_ERROR_SIZE];
};

using epcsaft_evaluate_pure_phase = int (*)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    epcsaft_phase_block_result* result
);

struct epcsaft_mixture_phase_block_result final {
    std::uint32_t struct_size;
    std::int32_t status;
    std::size_t coordinate_count;
    std::size_t gradient_capacity;
    std::size_t hessian_capacity;
    double* gradient;
    double* hessian;
    double helmholtz_over_rt_reference_amount;
    double pressure_pa;
    char parameter_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_ERROR_SIZE];
};

using epcsaft_evaluate_mixture_phase = int (*)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    std::size_t amount_count,
    double volume_m3,
    epcsaft_mixture_phase_block_result* result
);

using epcsaft_evaluate_molar_volume_bounds = int (*)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    const double* mole_fractions,
    std::size_t component_count,
    double packing_fraction_min,
    double packing_fraction_max,
    double* molar_volume_bounds_m3_per_mol,
    std::size_t bound_count
);

using epcsaft_evaluate_packing_fraction = int (*)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    const double* amounts_mol,
    std::size_t amount_count,
    double volume_m3,
    double* packing_fraction,
    double* gradient,
    std::size_t gradient_capacity,
    double* hessian,
    std::size_t hessian_capacity
);

// CAP-07 transport-only surface; this is not the EOS derivative interface.
enum epcsaft_native_parameter_family {
    EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER = 4u,
    EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K = 5u,
};

enum epcsaft_native_parameter_identity {
    EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT = 2u,
};

struct epcsaft_active_parameter_request final {
    std::uint32_t parameter_family;
    std::uint32_t parameter_identity;
    std::int32_t component_index;
    double value;
    const char* unit;
};

struct epcsaft_reacting_phase_parameter_result final {
    std::uint32_t struct_size;
    std::int32_t status;
    std::size_t state_coordinate_count;
    std::size_t active_parameter_count;
    std::size_t state_gradient_capacity;
    std::size_t state_hessian_capacity;
    std::size_t parameter_gradient_capacity;
    std::size_t state_parameter_capacity;
    std::size_t packing_parameter_capacity;
    std::size_t molar_volume_bound_parameter_capacity;
    std::size_t pressure_parameter_capacity;
    std::size_t chemical_potential_parameter_capacity;
    double* state_gradient;
    double* state_hessian;
    double* parameter_gradient;
    double* state_parameter_derivatives;
    double* packing_parameter_derivatives;
    double* molar_volume_bound_parameter_derivatives_m3_per_mol;
    double* pressure_parameter_derivatives_pa;
    double* chemical_potential_parameter_derivatives_over_rt;
    double helmholtz_over_rt_reference_amount;
    double pressure_pa;
    double packing_fraction;
    double molar_volume_bounds_m3_per_mol[2];
    char parameter_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char topology_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char helmholtz_basis_id[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_ERROR_SIZE];
    std::size_t packing_state_gradient_capacity;
    std::size_t packing_state_hessian_capacity;
    double* packing_state_gradient;
    double* packing_state_hessian;
};

using epcsaft_evaluate_reacting_phase_parameters = int (*)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    const char* expected_topology_fingerprint,
    double temperature_k,
    const double* amounts_mol,
    std::size_t amount_count,
    double volume_m3,
    double packing_fraction_min,
    double packing_fraction_max,
    const epcsaft_active_parameter_request* active_parameters,
    std::size_t active_parameter_count,
    epcsaft_reacting_phase_parameter_result* result
);

struct epcsaft_inverse_packing_geometry_result final {
    std::uint32_t struct_size;
    std::int32_t status;
    std::size_t coordinate_count;
    std::size_t active_parameter_count;
    std::size_t gradient_capacity;
    std::size_t hessian_capacity;
    double* gradient;
    double* hessian;
    double volume_m3;
    char parameter_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char topology_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_ERROR_SIZE];
};

using epcsaft_evaluate_inverse_packing_geometry = int (*)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    const char* expected_topology_fingerprint,
    double temperature_k,
    const double* amounts_mol,
    std::size_t amount_count,
    double log_packing_fraction,
    const epcsaft_active_parameter_request* active_parameters,
    std::size_t active_parameter_count,
    epcsaft_inverse_packing_geometry_result* result
);

enum epcsaft_source_reference_status {
    EPCSAFT_SOURCE_REFERENCE_OK = 0,
    EPCSAFT_SOURCE_REFERENCE_UNSUPPORTED = 1,
    EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE = 2,
    EPCSAFT_SOURCE_REFERENCE_COMPONENT_MISMATCH = 3,
    EPCSAFT_SOURCE_REFERENCE_MISSING_MODEL = 4,
    EPCSAFT_SOURCE_REFERENCE_DOMAIN_REJECTION = 5,
    EPCSAFT_SOURCE_REFERENCE_LIMIT_FAILURE = 6,
    EPCSAFT_SOURCE_REFERENCE_DERIVATIVE_UNAVAILABLE = 7,
    EPCSAFT_SOURCE_REFERENCE_NONFINITE = 8,
    EPCSAFT_SOURCE_REFERENCE_IDENTITY_MISMATCH = 9,
    EPCSAFT_SOURCE_REFERENCE_EVALUATION_ERROR = 10,
    EPCSAFT_SOURCE_REFERENCE_BUFFER_ERROR = 11,
};

inline constexpr const char*
    EPCSAFT_SOURCE_REFERENCE_ACTIVITY_MOLALITY_INFINITE_DILUTION =
        "molality-infinite-dilution-v1";

struct epcsaft_source_reference_request final {
    std::uint32_t struct_size;
    std::size_t component_count;
    const char* const* component_ids;
    std::size_t source_composition_count;
    const double* source_solvent_mole_fractions;
    const char* reference_state_id;
    const char* activity_convention_id;
    const char* reference_state_fingerprint;
    double standard_molality_mol_per_kg;
    double source_reference_pressure_pa;
    std::uint32_t required_derivatives;
};

struct epcsaft_source_reference_transfer_result final {
    std::uint32_t struct_size;
    std::int32_t status;
    std::size_t component_count;
    std::size_t neutral_basis_row_count;
    std::size_t neutral_basis_capacity;
    std::size_t contraction_capacity;
    std::size_t source_composition_capacity;
    std::size_t pressure_derivative_capacity;
    double* neutral_basis;
    double* log_fugacity_contractions;
    double* activity_scale_log_contractions;
    double* transfer_log_contractions;
    double* source_composition;
    double* pressure_derivatives_per_pa;
    std::uint32_t derivative_availability;
    double temperature_k;
    double pressure_pa;
    double source_reference_pressure_pa;
    double standard_molality_mol_per_kg;
    double solvent_molar_mass_kg_per_mol;
    double reference_limit_molality_mol_per_kg;
    double reference_convergence_error;
    double source_temperature_min_k;
    double source_temperature_max_k;
    double source_pressure_min_pa;
    double source_pressure_max_pa;
    char parameter_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char topology_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char component_order_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char reference_state_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char domain_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char artifact_fingerprint[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char helmholtz_basis_id[EPCSAFT_NATIVE_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_ERROR_SIZE];
};

using epcsaft_evaluate_source_reference_transfer = int (*)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    const char* expected_topology_fingerprint,
    double temperature_k,
    double pressure_pa,
    const epcsaft_source_reference_request* request,
    epcsaft_source_reference_transfer_result* result
);

struct epcsaft_native_transport final {
    const void* model_context;
    epcsaft_evaluate_pure_phase evaluate_pure_phase;
    std::size_t component_count;
    epcsaft_evaluate_mixture_phase evaluate_mixture_phase;
    epcsaft_evaluate_mixture_phase evaluate_electrolyte_phase;
    epcsaft_evaluate_molar_volume_bounds evaluate_molar_volume_bounds;
    epcsaft_evaluate_packing_fraction evaluate_packing_fraction;
    epcsaft_evaluate_reacting_phase_parameters evaluate_reacting_phase_parameters;
    epcsaft_evaluate_inverse_packing_geometry evaluate_inverse_packing_geometry;
    epcsaft_evaluate_source_reference_transfer evaluate_source_reference_transfer;
    const char* const* component_ids;
    const std::int32_t* component_charges;
    const char* parameter_fingerprint;
    const char* topology_fingerprint;
    double source_temperature_min_k;
    double source_temperature_max_k;
    double total_ion_mole_fraction_max;
    std::size_t neutral_reference_basis_row_count;
};

namespace epcsaft {

const epcsaft_native_transport* create_native_transport(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k,
    const std::vector<std::string>& component_ids,
    double ion_mole_fraction_max,
    double pressure_min_pa,
    double pressure_max_pa,
    const std::string& topology_fingerprint,
    const std::string& component_order_fingerprint,
    const std::string& domain_fingerprint,
    const std::string& artifact_fingerprint
);

void destroy_native_transport(const epcsaft_native_transport* transport) noexcept;

}  // namespace epcsaft
