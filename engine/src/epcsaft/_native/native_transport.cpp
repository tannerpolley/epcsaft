#include "native_transport.hpp"
#include "density.hpp"
#include "electrolyte_reference.hpp"
#include "phase_block.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

// EqExclude: input and domain validation
// EqExclude: input and domain validation
struct NeutralReferenceBasis final {
    struct IonPair final {
        std::size_t cation;
        std::size_t anion;
        std::array<int, 2> stoichiometry;
    };

    std::vector<std::size_t> neutral_components;
    std::vector<IonPair> ion_pairs;
    std::vector<double> reference_composition;
    std::vector<double> rows;

    [[nodiscard]] std::size_t row_count() const {
        return reference_composition.empty()
            ? 0u
            : rows.size() / reference_composition.size();
    }
};

// EqExclude: input and domain validation
bool full_row_rank(
    const std::vector<double>& values,
    std::size_t rows,
    std::size_t columns
) {
    if (rows == 0 || columns == 0 || values.size() != rows * columns) {
        return false;
    }
    std::vector<double> matrix(values);
    std::size_t rank = 0;
    for (std::size_t column = 0; column < columns && rank < rows; ++column) {
        std::size_t pivot = rank;
        while (pivot < rows
               && std::abs(matrix[pivot * columns + column]) <= 1.0e-12) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        for (std::size_t item = column; item < columns; ++item) {
            std::swap(
                matrix[rank * columns + item],
                matrix[pivot * columns + item]
            );
        }
        const double divisor = matrix[rank * columns + column];
        for (std::size_t item = column; item < columns; ++item) {
            matrix[rank * columns + item] /= divisor;
        }
        for (std::size_t row = 0; row < rows; ++row) {
            if (row == rank) {
                continue;
            }
            const double factor = matrix[row * columns + column];
            for (std::size_t item = column; item < columns; ++item) {
                matrix[row * columns + item] -=
                    factor * matrix[rank * columns + item];
            }
        }
        ++rank;
    }
    return rank == rows;
}

// EqID: neutral_reference_stoichiometric_basis
NeutralReferenceBasis make_neutral_reference_basis(
    const std::vector<std::string>& component_ids,
    const std::vector<int>& charges,
    const std::vector<double>& reference_composition,
    const std::vector<NeutralReferenceBasis::IonPair>& ion_pairs
) {
    NeutralReferenceBasis result;
    const std::size_t components = component_ids.size();
    if (components < 3 || components != charges.size()
        || reference_composition.size() != components
        || std::any_of(
            component_ids.begin(), component_ids.end(), [](const auto& value) {
                return value.empty();
            }
        )
        || std::any_of(charges.begin(), charges.end(), [](int charge) {
            return charge == std::numeric_limits<int>::min();
        })) {
        return result;
    }
    std::vector<std::size_t> order(components);
    std::iota(order.begin(), order.end(), 0u);
    std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
        return component_ids[left] < component_ids[right];
    });
    for (std::size_t index = 1; index < order.size(); ++index) {
        if (component_ids[order[index - 1]] == component_ids[order[index]]) {
            return {};
        }
    }
    for (std::size_t component : order) {
        if (charges[component] == 0) {
            result.neutral_components.push_back(component);
        }
    }
    if (result.neutral_components.empty() || ion_pairs.empty()) {
        return {};
    }
    result.reference_composition = reference_composition;
    auto append_row = [&](std::size_t left, double left_value,
                          std::size_t right = 0u, double right_value = 0.0) {
        result.rows.resize(result.rows.size() + components, 0.0);
        const std::size_t offset = result.rows.size() - components;
        result.rows[offset + left] = left_value;
        if (right_value != 0.0) {
            result.rows[offset + right] = right_value;
        }
    };
    for (std::size_t component : result.neutral_components) {
        append_row(component, 1.0);
    }
    for (const auto& pair : ion_pairs) {
        if (pair.cation >= components || pair.anion >= components
            || charges[pair.cation] <= 0 || charges[pair.anion] >= 0
            || pair.stoichiometry[0] <= 0 || pair.stoichiometry[1] <= 0) {
            return {};
        }
        const std::int64_t charge_balance =
            static_cast<std::int64_t>(pair.stoichiometry[0])
                * static_cast<std::int64_t>(charges[pair.cation])
            + static_cast<std::int64_t>(pair.stoichiometry[1])
                * static_cast<std::int64_t>(charges[pair.anion]);
        if (charge_balance != 0) {
            return {};
        }
        result.ion_pairs.push_back(pair);
        append_row(
            pair.cation,
            static_cast<double>(pair.stoichiometry[0]),
            pair.anion,
            static_cast<double>(pair.stoichiometry[1])
        );
    }
    const std::size_t expected_rows = components - 1u;
    if (result.row_count() != expected_rows
        || !full_row_rank(result.rows, expected_rows, components)) {
        return {};
    }
    for (std::size_t row = 0; row < expected_rows; ++row) {
        double net_charge = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            net_charge += result.rows[row * components + component]
                * static_cast<double>(charges[component]);
        }
        if (net_charge != 0.0) {
            return {};
        }
    }
    return result;
}

// EqExclude: input and domain validation
bool sha256_fingerprint_valid(const std::string& fingerprint) {
    return fingerprint.size() == 71
        && fingerprint.rfind("sha256:", 0) == 0
        && std::all_of(
            fingerprint.begin() + 7,
            fingerprint.end(),
            [](unsigned char character) {
                return std::isxdigit(character) != 0;
            }
        );
}

struct NativeTransportEndpoint final {
    NativeModel model;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    double temperature_min_k;
    double temperature_max_k;
    double pressure_min_pa;
    double pressure_max_pa;
    double ion_mole_fraction_max;
    bool phase_callbacks_admitted;
    bool neutral_reference_admitted;
    std::vector<std::string> component_ids;
    std::vector<const char*> component_id_pointers;
    std::vector<std::int32_t> component_charges;
    std::vector<NativeActiveCoordinate> active_coordinates;
    std::vector<const char*> active_identity_pointers;
    epcsaft_native_transport table;

    NativeTransportEndpoint(
        const NativeModel& source_model,
        std::string fingerprint,
        double minimum_temperature,
        double maximum_temperature,
        std::vector<std::string> identifiers,
        double maximum_ion_fraction,
        double minimum_pressure,
        double maximum_pressure,
        std::string topology_identity,
        std::string component_order_identity,
        std::string domain_identity,
        std::string artifact_identity,
        std::vector<NativeActiveCoordinate> coordinates
    )
        : model(coordinates.empty()
              ? source_model
              : source_model.with_active_coordinates(coordinates)),
          parameter_fingerprint(std::move(fingerprint)),
          topology_fingerprint(std::move(topology_identity)),
          component_order_fingerprint(std::move(component_order_identity)),
          domain_fingerprint(std::move(domain_identity)),
          artifact_fingerprint(std::move(artifact_identity)),
          temperature_min_k(minimum_temperature),
          temperature_max_k(maximum_temperature),
          pressure_min_pa(minimum_pressure),
          pressure_max_pa(maximum_pressure),
          ion_mole_fraction_max(maximum_ion_fraction),
          phase_callbacks_admitted(
              source_model.supports_pure_transport()
              || source_model.supports_mixture_transport()
              || source_model.supports_fixed_topology_associating_mixture()
              || source_model.supports_electrolyte_transport()
          ),
          neutral_reference_admitted(
              source_model.supports_electrolyte_transport()
              && sha256_fingerprint_valid(topology_fingerprint)
              && sha256_fingerprint_valid(component_order_fingerprint)
              && sha256_fingerprint_valid(domain_fingerprint)
              && sha256_fingerprint_valid(artifact_fingerprint)
              && std::isfinite(minimum_temperature)
              && std::isfinite(maximum_temperature)
              && std::isfinite(minimum_pressure)
              && std::isfinite(maximum_pressure)
          ),
          component_ids(std::move(identifiers)),
          component_id_pointers([&] {
              std::vector<const char*> pointers;
              pointers.reserve(component_ids.size());
              for (const auto& identifier : component_ids) {
                  pointers.push_back(identifier.c_str());
              }
              return pointers;
          }()),
          component_charges(
              source_model.charge_numbers().begin(),
              source_model.charge_numbers().end()
          ),
          active_coordinates(std::move(coordinates)),
          active_identity_pointers([&] {
              std::vector<const char*> pointers;
              pointers.reserve(active_coordinates.size());
              for (const auto& coordinate : active_coordinates) {
                  pointers.push_back(coordinate.identity.c_str());
              }
              return pointers;
          }()),
          table{
              this,
              source_model.supports_pure_transport()
                  ? &NativeTransportEndpoint::evaluate
                  : nullptr,
              component_ids.size(),
              (source_model.supports_mixture_transport()
                  || source_model.supports_fixed_topology_associating_mixture())
                  ? &NativeTransportEndpoint::evaluate_mixture
                  : nullptr,
              source_model.supports_electrolyte_transport()
                  ? &NativeTransportEndpoint::evaluate_electrolyte
                  : nullptr,
              phase_callbacks_admitted
                  ? &NativeTransportEndpoint::evaluate_molar_volume_bounds
                  : nullptr,
              phase_callbacks_admitted
                  ? &NativeTransportEndpoint::evaluate_packing_fraction
                  : nullptr,
              phase_callbacks_admitted
                      && !topology_fingerprint.empty()
                  ? &NativeTransportEndpoint::evaluate_parameterized_phase
                  : nullptr,
              phase_callbacks_admitted
                  ? &NativeTransportEndpoint::evaluate_temperature_phase
                  : nullptr,
              phase_callbacks_admitted && !topology_fingerprint.empty()
                  ? &NativeTransportEndpoint::evaluate_inverse_packing_geometry
                  : nullptr,
              neutral_reference_admitted
                  ? &NativeTransportEndpoint::evaluate_source_reference_transfer
                  : nullptr,
              component_id_pointers.data(),
              component_charges.data(),
              active_coordinates.size(),
              active_identity_pointers.data(),
              parameter_fingerprint.c_str(),
              topology_fingerprint.c_str(),
              temperature_min_k,
              temperature_max_k,
              ion_mole_fraction_max,
          } {}

    static int evaluate_temperature_phase(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_temperature_phase_result* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if (result->struct_size != sizeof(epcsaft_temperature_phase_result)) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        result->status = EPCSAFT_NATIVE_STATUS_OK;
        result->pressure_temperature_derivative_pa_per_k =
            std::numeric_limits<double>::quiet_NaN();
        result->parameter_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || amounts_mol == nullptr) {
            return fail(EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                        "temperature phase required pointer is null");
        }
        const auto& endpoint =
            *static_cast<const NativeTransportEndpoint*>(context);
        const std::size_t support_count = static_cast<std::size_t>(std::count_if(
            amounts_mol, amounts_mol + amount_count,
            [](double amount) { return amount > 0.0; }
        ));
        if (amount_count != endpoint.component_ids.size()
            || result->supported_component_count != support_count
            || result->chemical_potential_temperature_capacity != support_count
            || result->chemical_potential_temperature_derivatives_over_rt_per_k
                == nullptr) {
            return fail(EPCSAFT_NATIVE_STATUS_BUFFER_ERROR,
                        "temperature phase buffers do not match immutable dimensions");
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(volume_m3) || volume_m3 <= 0.0
            || (std::isfinite(endpoint.temperature_min_k)
                && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                        "temperature phase state is outside the source domain");
        }
        try {
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            const auto phase = endpoint.model.evaluate_temperature_phase_derivatives(
                temperature_k, amounts, volume_m3,
                endpoint.ion_mole_fraction_max
            );
            const std::size_t state_count = support_count + 1u;
            const std::size_t width = state_count + 1u;
            const std::size_t volume_coordinate = support_count;
            const std::size_t temperature_coordinate = state_count;
            if (phase.gradient.size() != width
                || phase.hessian.size() != width * width) {
                throw std::runtime_error(
                    "temperature phase tensor dimensions changed"
                );
            }
            for (std::size_t component = 0; component < support_count; ++component) {
                result->chemical_potential_temperature_derivatives_over_rt_per_k[
                    component
                ] = phase.hessian[
                    component * width + temperature_coordinate
                ];
            }
            result->pressure_temperature_derivative_pa_per_k =
                -detail::GAS_CONSTANT * phase.gradient[volume_coordinate]
                -detail::GAS_CONSTANT * temperature_k * phase.hessian[
                    volume_coordinate * width + temperature_coordinate
                ];
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s", endpoint.parameter_fingerprint.c_str()
            );
            return result->status = EPCSAFT_NATIVE_STATUS_OK;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR, error.what());
        } catch (...) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR,
                        "unknown temperature phase evaluation failure");
        }
    }

    // EqID: phase_packing_fraction
    // EqExclude: transport and buffer mechanics
    static int evaluate_inverse_packing_geometry(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double log_packing_fraction,
        epcsaft_inverse_packing_geometry_result* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if (result->struct_size
            != sizeof(epcsaft_inverse_packing_geometry_result)) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        result->status = EPCSAFT_NATIVE_STATUS_OK;
        result->volume_m3 = std::numeric_limits<double>::quiet_NaN();
        result->parameter_fingerprint[0] = '\0';
        result->topology_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || amounts_mol == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "inverse packing required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeTransportEndpoint*>(context);
        const std::size_t request_count = endpoint.active_coordinates.size();
        if (endpoint.topology_fingerprint.empty()) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL,
                "inverse packing model or topology is not admitted"
            );
        }
        const std::size_t support_count = static_cast<std::size_t>(std::count_if(
            amounts_mol,
            amounts_mol + amount_count,
            [](double amount) { return amount > 0.0; }
        ));
        const std::size_t coordinate_count = support_count + 1 + request_count;
        if (amount_count != endpoint.table.component_count
            || result->coordinate_count != coordinate_count
            || result->active_parameter_count != request_count
            || result->gradient_capacity != coordinate_count
            || result->hessian_capacity
                != coordinate_count * coordinate_count
            || result->gradient == nullptr || result->hessian == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_BUFFER_ERROR,
                "inverse packing shape or capacity does not match"
            );
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::all_of(
                amounts_mol,
                amounts_mol + amount_count,
                [](double amount) {
                    return std::isfinite(amount) && amount >= 0.0;
                }
            )
            || !std::isfinite(log_packing_fraction)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "inverse packing state is invalid"
            );
        }
        const double packing_fraction = std::exp(log_packing_fraction);
        if (!std::isfinite(packing_fraction) || packing_fraction <= 0.0
            || packing_fraction >= detail::close_packing_fraction()) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                "inverse packing fraction is outside the model domain"
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            const auto block =
                endpoint.model.evaluate_inverse_packing_geometry(
                    temperature_k,
                    amounts,
                    log_packing_fraction,
                    endpoint.active_coordinates
                );
            std::copy(
                block.gradient.begin(),
                block.gradient.end(),
                result->gradient
            );
            std::copy(
                block.hessian.begin(),
                block.hessian.end(),
                result->hessian
            );
            result->volume_m3 = block.volume_m3;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            std::snprintf(
                result->topology_fingerprint,
                sizeof(result->topology_fingerprint),
                "%s",
                endpoint.topology_fingerprint.c_str()
            );
            return result->status = EPCSAFT_NATIVE_STATUS_OK;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR,
                "unknown inverse packing evaluation failure"
            );
        }
    }

    // EqID: neutral_reference_stoichiometric_basis
    static int evaluate_parameterized_phase(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double packing_fraction_min,
        double packing_fraction_max,
        epcsaft_parameterized_phase_result* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if (result->struct_size
            != sizeof(epcsaft_parameterized_phase_result)) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result->status = EPCSAFT_NATIVE_STATUS_OK;
        result->helmholtz_over_rt_reference_amount = nan;
        result->pressure_pa = nan;
        result->packing_fraction = nan;
        std::fill(
            std::begin(result->molar_volume_bounds_m3_per_mol),
            std::end(result->molar_volume_bounds_m3_per_mol),
            nan
        );
        result->parameter_fingerprint[0] = '\0';
        result->topology_fingerprint[0] = '\0';
        result->helmholtz_basis_id[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || amounts_mol == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "parameterized phase required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeTransportEndpoint*>(context);
        const std::size_t request_count = endpoint.active_coordinates.size();
        if (!endpoint.phase_callbacks_admitted
            || endpoint.topology_fingerprint.empty()) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL,
                "parameterized phase model identity is unsupported"
            );
        }
        if (!std::isfinite(temperature_k)
            || temperature_k < endpoint.temperature_min_k
            || temperature_k > endpoint.temperature_max_k) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                "parameterized phase temperature is outside the source domain"
            );
        }
        const std::size_t support_count = static_cast<std::size_t>(std::count_if(
            amounts_mol,
            amounts_mol + amount_count,
            [](double amount) { return amount > 0.0; }
        ));
        const std::size_t state_count = support_count + 1;
        if (amount_count != endpoint.component_ids.size()
            || result->state_coordinate_count != state_count
            || result->active_parameter_count != request_count
            || result->state_gradient_capacity != state_count
            || result->state_hessian_capacity != state_count * state_count
            || result->parameter_gradient_capacity != request_count
            || result->state_parameter_capacity != state_count * request_count
            || result->packing_parameter_capacity != request_count
            || result->packing_state_gradient_capacity != state_count
            || result->packing_state_hessian_capacity
                != state_count * state_count
            || result->molar_volume_bound_parameter_capacity
                != 2u * request_count
            || result->pressure_parameter_capacity != request_count
            || result->chemical_potential_parameter_capacity
                != support_count * request_count
            || result->state_gradient == nullptr
            || result->state_hessian == nullptr
            || (request_count != 0 && result->parameter_gradient == nullptr)
            || (request_count != 0 && result->state_parameter_derivatives == nullptr)
            || (request_count != 0 && result->packing_parameter_derivatives == nullptr)
            || result->packing_state_gradient == nullptr
            || result->packing_state_hessian == nullptr
            || (request_count != 0
                && result->molar_volume_bound_parameter_derivatives_m3_per_mol
                    == nullptr)
            || (request_count != 0
                && result->pressure_parameter_derivatives_pa == nullptr)
            || (request_count != 0
                && result->chemical_potential_parameter_derivatives_over_rt
                    == nullptr)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_BUFFER_ERROR,
                "parameterized phase buffers do not match immutable dimensions"
            );
        }
        try {
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            double total_amount = 0.0;
            for (const double amount : amounts) {
                total_amount += amount;
            }
            std::vector<double> composition(amount_count);
            for (std::size_t component = 0;
                 component < amount_count;
                 ++component) {
                composition[component] = amounts[component] / total_amount;
            }
            const auto packing =
                endpoint.model.evaluate_packing_fraction_parameters(
                    temperature_k,
                    amounts,
                    volume_m3,
                    endpoint.ion_mole_fraction_max,
                    endpoint.active_coordinates
                );
            const auto phase = endpoint.model.evaluate_parameterized_phase(
                temperature_k,
                amounts,
                volume_m3,
                endpoint.ion_mole_fraction_max,
                endpoint.active_coordinates
            );
            const auto volume_bounds =
                endpoint.model.molar_volume_bounds_parameters(
                    temperature_k,
                    composition,
                    packing_fraction_min,
                    packing_fraction_max,
                    endpoint.ion_mole_fraction_max,
                    endpoint.active_coordinates
                );
            const std::size_t joint_count = state_count + request_count;
            if (phase.gradient.size() != joint_count
                || phase.hessian.size() != joint_count * joint_count
                || packing.gradient.size() != joint_count
                || packing.hessian.size() != joint_count * joint_count) {
                throw std::runtime_error(
                    "parameterized phase tensor dimensions changed"
                );
            }
            std::vector<double> state_gradient(state_count);
            std::vector<double> state_hessian(state_count * state_count);
            std::vector<double> parameter_gradient(request_count);
            std::vector<double> mixed(state_count * request_count);
            std::vector<double> packing_parameter(request_count);
            std::vector<double> packing_state_gradient(state_count);
            std::vector<double> packing_state_hessian(
                state_count * state_count
            );
            std::vector<double> pressure_parameter(request_count);
            std::vector<double> chemical_parameter(support_count * request_count);
            for (std::size_t state = 0; state < state_count; ++state) {
                state_gradient[state] = phase.gradient[state];
                packing_state_gradient[state] = packing.gradient[state];
                for (std::size_t other = 0; other < state_count; ++other) {
                    state_hessian[state * state_count + other] =
                        phase.hessian[state * joint_count + other];
                    packing_state_hessian[state * state_count + other] =
                        packing.hessian[state * joint_count + other];
                }
                for (std::size_t parameter = 0;
                     parameter < request_count;
                     ++parameter) {
                    mixed[state * request_count + parameter] =
                        phase.hessian[
                            state * joint_count + state_count + parameter
                        ];
                }
            }
            for (std::size_t parameter = 0;
                 parameter < request_count;
                 ++parameter) {
                parameter_gradient[parameter] =
                    phase.gradient[state_count + parameter];
                packing_parameter[parameter] =
                    packing.gradient[state_count + parameter];
                pressure_parameter[parameter] =
                    -detail::GAS_CONSTANT * temperature_k
                    * mixed[
                        (state_count - 1) * request_count + parameter
                    ];
                for (std::size_t component = 0;
                     component < support_count;
                     ++component) {
                    chemical_parameter[
                        component * request_count + parameter
                    ] = mixed[component * request_count + parameter];
                }
            }
            const auto finite_values = [](const std::vector<double>& values) {
                return std::all_of(
                    values.begin(), values.end(), [](double value) {
                        return std::isfinite(value);
                    }
                );
            };
            if (!finite_values(state_gradient)
                || !finite_values(state_hessian)
                || !finite_values(parameter_gradient)
                || !finite_values(mixed)
                || !finite_values(packing_parameter)
                || !finite_values(packing_state_gradient)
                || !finite_values(packing_state_hessian)
                || !finite_values(
                    volume_bounds.parameter_derivatives_m3_per_mol
                )
                || !finite_values(pressure_parameter)
                || !finite_values(chemical_parameter)
                || !std::isfinite(phase.helmholtz_over_rt_reference_amount)
                || !std::isfinite(phase.pressure_pa)
                || !std::isfinite(packing.packing_fraction)) {
                throw std::runtime_error(
                    "parameterized phase result contains nonfinite values"
                );
            }
            std::copy(
                state_gradient.begin(),
                state_gradient.end(),
                result->state_gradient
            );
            std::copy(
                state_hessian.begin(),
                state_hessian.end(),
                result->state_hessian
            );
            std::copy(
                parameter_gradient.begin(),
                parameter_gradient.end(),
                result->parameter_gradient
            );
            std::copy(
                mixed.begin(),
                mixed.end(),
                result->state_parameter_derivatives
            );
            std::copy(
                packing_parameter.begin(),
                packing_parameter.end(),
                result->packing_parameter_derivatives
            );
            std::copy(
                packing_state_gradient.begin(),
                packing_state_gradient.end(),
                result->packing_state_gradient
            );
            std::copy(
                packing_state_hessian.begin(),
                packing_state_hessian.end(),
                result->packing_state_hessian
            );
            std::copy(
                volume_bounds.parameter_derivatives_m3_per_mol.begin(),
                volume_bounds.parameter_derivatives_m3_per_mol.end(),
                result
                    ->molar_volume_bound_parameter_derivatives_m3_per_mol
            );
            std::copy(
                pressure_parameter.begin(),
                pressure_parameter.end(),
                result->pressure_parameter_derivatives_pa
            );
            std::copy(
                chemical_parameter.begin(),
                chemical_parameter.end(),
                result->chemical_potential_parameter_derivatives_over_rt
            );
            result->helmholtz_over_rt_reference_amount =
                phase.helmholtz_over_rt_reference_amount;
            result->pressure_pa = phase.pressure_pa;
            result->packing_fraction = packing.packing_fraction;
            std::copy(
                volume_bounds.bounds_m3_per_mol.begin(),
                volume_bounds.bounds_m3_per_mol.end(),
                result->molar_volume_bounds_m3_per_mol
            );
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            std::snprintf(
                result->topology_fingerprint,
                sizeof(result->topology_fingerprint),
                "%s",
                endpoint.topology_fingerprint.c_str()
            );
            std::snprintf(
                result->helmholtz_basis_id,
                sizeof(result->helmholtz_basis_id),
                "%s",
                EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK;
            return result->status;
        } catch (const detail::DensityDomainError& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR, error.what());
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR, error.what());
        } catch (const std::logic_error& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL, error.what()
            );
        } catch (const std::exception& error) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR, error.what()
            );
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR,
                "unknown parameterized phase evaluation failure"
            );
        }
    }

    // EqID: declared_solvent_source_reference_transfer
    static int evaluate_source_reference_transfer(
        const void* context,
        double temperature_k,
        double pressure_pa,
        const epcsaft_source_reference_request* request,
        epcsaft_source_reference_transfer_result* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_SOURCE_REFERENCE_BUFFER_ERROR;
        }
        if (result->struct_size
            != sizeof(epcsaft_source_reference_transfer_result)) {
            return EPCSAFT_SOURCE_REFERENCE_BUFFER_ERROR;
        }
        result->status = EPCSAFT_SOURCE_REFERENCE_OK;
        result->derivative_availability =
            0u;
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr || request == nullptr) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_MISSING_MODEL,
                "source-reference model contract is incomplete"
            );
        }
        if (request->struct_size != sizeof(epcsaft_source_reference_request)) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_BUFFER_ERROR,
                "source-reference request shape is incompatible"
            );
        }
        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.neutral_reference_admitted
            || endpoint.topology_fingerprint.empty()) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_MISSING_MODEL,
                "source-reference model identity is incompatible"
            );
        }
        if (!std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(pressure_pa) || pressure_pa <= 0.0) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_NONFINITE,
                "source-reference temperature and pressure must be positive and finite"
            );
        }
        if (temperature_k < endpoint.temperature_min_k
            || temperature_k > endpoint.temperature_max_k
            || pressure_pa < endpoint.pressure_min_pa
            || pressure_pa > endpoint.pressure_max_pa) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_DOMAIN_REJECTION,
                "source-reference trial state is outside the Provider domain"
            );
        }
        const std::size_t component_count = endpoint.component_ids.size();
        const std::size_t row_count = component_count - 1u;
        if (request->component_count != component_count
            || request->source_composition_count != component_count
            || request->component_ids == nullptr
            || request->source_solvent_mole_fractions == nullptr
            || request->ion_pair_count == 0
            || request->ion_pair_cation_indices == nullptr
            || request->ion_pair_anion_indices == nullptr
            || request->ion_pair_cation_stoichiometry == nullptr
            || request->ion_pair_anion_stoichiometry == nullptr) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_COMPONENT_MISMATCH,
                "source-reference component dimensions do not match"
            );
        }
        for (std::size_t component = 0; component < component_count; ++component) {
            if (request->component_ids[component] == nullptr
                || endpoint.component_ids[component]
                    != request->component_ids[component]) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_COMPONENT_MISMATCH,
                    "source-reference component order does not match"
                );
            }
        }
        const bool molality_reference =
            request->activity_convention_id != nullptr
            && std::strcmp(
                request->activity_convention_id,
                EPCSAFT_SOURCE_REFERENCE_ACTIVITY_MOLALITY_INFINITE_DILUTION
            ) == 0
            && request->reference_convention != nullptr
            && std::strcmp(
                request->reference_convention,
                "molality-infinite-dilution"
            ) == 0;
        const bool mixed_mole_fraction_reference =
            request->activity_convention_id != nullptr
            && std::strcmp(
                request->activity_convention_id,
                EPCSAFT_SOURCE_REFERENCE_ACTIVITY_MOLE_FRACTION_MIXED_REFERENCE
            ) == 0
            && request->reference_convention != nullptr
            && std::strcmp(
                request->reference_convention,
                "mole-fraction-mixed-reference"
            ) == 0;
        if (request->reference_state_id == nullptr
            || request->reference_state_id[0] == '\0'
            || (!molality_reference && !mixed_mole_fraction_reference)
            || request->phase == nullptr
            || std::strcmp(request->phase, "liquid") != 0
            || request->reference_state_fingerprint == nullptr
            || !sha256_fingerprint_valid(request->reference_state_fingerprint)
            || !std::isfinite(request->standard_molality_mol_per_kg)
            || request->standard_molality_mol_per_kg <= 0.0
            || !std::isfinite(request->source_reference_pressure_pa)
            || request->source_reference_pressure_pa <= 0.0) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_UNSUPPORTED,
                "source-reference declaration or activity convention is unsupported"
            );
        }
        constexpr std::uint32_t pressure_derivative = 2u;
        constexpr std::uint32_t parameter_derivative = 8u;
        constexpr std::uint32_t supported_derivatives =
            pressure_derivative | parameter_derivative;
        const bool pressure_derivatives_requested =
            (request->required_derivatives & pressure_derivative) != 0u;
        const bool parameter_derivatives_requested =
            (request->required_derivatives & parameter_derivative) != 0u;
        if ((request->required_derivatives & ~supported_derivatives) != 0u
            || (parameter_derivatives_requested
                && endpoint.active_coordinates.empty())) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_DERIVATIVE_UNAVAILABLE,
                "required source-reference derivative is unavailable"
            );
        }
        const std::size_t active_parameter_count =
            endpoint.active_coordinates.size();
        if (result->component_count != component_count
            || result->neutral_basis_row_count != row_count
            || result->neutral_basis_capacity != row_count * component_count
            || result->contraction_capacity != row_count
            || result->source_composition_capacity != component_count
            || result->pressure_derivative_capacity
                != (pressure_derivatives_requested ? row_count : 0u)
            || result->parameter_derivative_capacity
                != (parameter_derivatives_requested
                        ? row_count * active_parameter_count
                        : 0u)
            || result->neutral_basis == nullptr
            || result->log_fugacity_contractions == nullptr
            || result->activity_scale_log_contractions == nullptr
            || result->transfer_log_contractions == nullptr
            || result->source_composition == nullptr
            || (pressure_derivatives_requested
                && result->pressure_derivatives_per_pa == nullptr)
            || (parameter_derivatives_requested
                && result->parameter_derivatives == nullptr)) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_BUFFER_ERROR,
                "source-reference result buffers do not match immutable dimensions"
            );
        }
        try {
            std::vector<double> source_composition(
                request->source_solvent_mole_fractions,
                request->source_solvent_mole_fractions + component_count
            );
            double sum = 0.0;
            for (std::size_t component = 0; component < component_count; ++component) {
                const double fraction = source_composition[component];
                if (!std::isfinite(fraction) || fraction < 0.0
                    || (endpoint.component_charges[component] != 0
                        && fraction != 0.0)) {
                    return fail(
                        EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                        "source composition must be finite, salt-free, and nonnegative"
                    );
                }
                sum += fraction;
            }
            if (std::abs(sum - 1.0) > 1.0e-12) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                    "source composition must sum to one"
                );
            }
            std::vector<bool> pure_component(component_count, false);
            if (request->pure_component_count > 0u
                && (request->pure_component_indices == nullptr
                    || request->pure_component_phases == nullptr)) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                    "pure-component reference identities or phases are unavailable"
                );
            }
            std::vector<std::string> pure_component_phase(component_count);
            for (std::size_t index = 0;
                 index < request->pure_component_count;
                 ++index) {
                const std::size_t component = request->pure_component_indices[index];
                const char* phase = request->pure_component_phases[index];
                if (component >= component_count
                    || endpoint.component_charges[component] != 0
                    || pure_component[component]
                    || phase == nullptr
                    || (std::strcmp(phase, "liquid") != 0
                        && std::strcmp(phase, "vapor") != 0)) {
                    return fail(
                        EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                        "pure-component references must be unique neutral components "
                        "with liquid or vapor phases"
                    );
                }
                pure_component[component] = true;
                pure_component_phase[component] = phase;
            }
            if (molality_reference && request->pure_component_count != 0u) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                    "molality reference does not admit pure-component overrides"
                );
            }
            std::vector<NeutralReferenceBasis::IonPair> ion_pairs;
            ion_pairs.reserve(request->ion_pair_count);
            for (std::size_t index = 0; index < request->ion_pair_count; ++index) {
                ion_pairs.push_back({
                    request->ion_pair_cation_indices[index],
                    request->ion_pair_anion_indices[index],
                    {
                        request->ion_pair_cation_stoichiometry[index],
                        request->ion_pair_anion_stoichiometry[index],
                    },
                });
            }
            const auto reference_basis = make_neutral_reference_basis(
                endpoint.component_ids,
                endpoint.model.charge_numbers(),
                source_composition,
                ion_pairs
            );
            if (reference_basis.row_count() != row_count) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                    "source-reference ion pairs do not define a full neutral basis"
                );
            }
            const auto& molar_masses = endpoint.model.molar_masses_kg_per_mol();
            const double solvent_molar_mass = std::inner_product(
                source_composition.begin(),
                source_composition.end(),
                molar_masses.begin(),
                0.0
            );
            if (!std::isfinite(solvent_molar_mass)
                || solvent_molar_mass <= 0.0) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE,
                    "source solvent molar mass is unavailable"
                );
            }
            std::vector<double> fugacity(row_count, 0.0);
            std::vector<double> activity(row_count, 0.0);
            std::vector<double> transfer(row_count, 0.0);
            std::vector<double> fugacity_pressure_derivatives(
                pressure_derivatives_requested ? row_count : 0u, 0.0
            );
            std::vector<double> fugacity_parameter_derivatives(
                parameter_derivatives_requested
                    ? row_count * active_parameter_count
                    : 0u,
                0.0
            );
            const std::size_t derivative_stride = 1u
                + (pressure_derivatives_requested ? 1u : 0u)
                + (parameter_derivatives_requested
                    ? active_parameter_count
                    : 0u);
            std::vector<std::size_t> neutral_solutes;
            for (std::size_t component : reference_basis.neutral_components) {
                if (source_composition[component] == 0.0
                    && !pure_component[component]) {
                    neutral_solutes.push_back(component);
                }
            }
            std::vector<double> neutral_trace_stoichiometry(
                component_count, 0.0
            );
            for (std::size_t component : neutral_solutes) {
                neutral_trace_stoichiometry[component] = 1.0;
            }
            const auto reference_composition = [&] (
                const std::vector<double>& trace_stoichiometry,
                double molality
            ) {
                const double trace_formula_moles = std::accumulate(
                    trace_stoichiometry.begin(),
                    trace_stoichiometry.end(),
                    0.0
                );
                const double solvent_moles = 1.0 / solvent_molar_mass;
                const double total_moles = solvent_moles
                    + molality * trace_formula_moles;
                std::vector<double> composition(component_count, 0.0);
                for (std::size_t component = 0;
                     component < component_count;
                     ++component) {
                    composition[component] =
                        source_composition[component] * solvent_moles
                        / total_moles
                        + molality * trace_stoichiometry[component]
                        / total_moles;
                }
                return composition;
            };
            std::optional<double> density_anchor;
            if (pressure_derivatives_requested
                || parameter_derivatives_requested) {
                const std::vector<double> density_anchor_composition =
                    reference_composition(
                        neutral_trace_stoichiometry,
                        detail::REFERENCE_MOLALITIES[2]
                    );
                density_anchor = endpoint.model.evaluate_pressure(
                    temperature_k,
                    pressure_pa,
                    density_anchor_composition,
                    "liquid"
                ).molar_density;
            }
            const auto phase_fugacity_derivatives =
                [&](const std::vector<double>& composition,
                    const char* phase_branch,
                    bool use_density_anchor) {
                    const auto state = density_anchor.has_value()
                            && use_density_anchor
                        ? endpoint.model.evaluate_pressure_continuation(
                            temperature_k,
                            pressure_pa,
                            composition,
                            *density_anchor
                        )
                        : endpoint.model.evaluate_pressure(
                            temperature_k,
                            pressure_pa,
                            composition,
                            phase_branch
                        );
                    if (!state.log_fugacity_coefficient.has_value()) {
                        throw std::runtime_error(
                            "source-reference fugacity is unavailable"
                        );
                    }
                    std::vector<double> values(
                        component_count * derivative_stride, 0.0
                    );
                    for (std::size_t component = 0;
                         component < component_count;
                         ++component) {
                        values[component * derivative_stride] =
                            (*state.log_fugacity_coefficient)[component];
                    }
                    if (!pressure_derivatives_requested
                        && !parameter_derivatives_requested) {
                        return values;
                    }
                    const auto absent_coordinate =
                        std::numeric_limits<std::size_t>::max();
                    std::vector<std::size_t> support_coordinates(
                        component_count, absent_coordinate
                    );
                    std::size_t support_count = 0u;
                    for (std::size_t component = 0;
                         component < component_count;
                         ++component) {
                        if (composition[component] > 0.0) {
                            support_coordinates[component] = support_count++;
                        }
                        std::fill(
                            values.begin()
                                + static_cast<std::ptrdiff_t>(
                                    component * derivative_stride + 1u
                                ),
                            values.begin()
                                + static_cast<std::ptrdiff_t>(
                                    (component + 1u) * derivative_stride
                                ),
                            std::numeric_limits<double>::quiet_NaN()
                        );
                    }
                    const double volume_m3 = 1.0 / state.molar_density;
                    const NativeMixturePhaseBlock phase =
                        endpoint.model.evaluate_supported_phase_derivatives(
                            temperature_k,
                            composition,
                            volume_m3,
                            endpoint.ion_mole_fraction_max,
                            endpoint.active_coordinates
                        );
                    const std::size_t state_count = support_count + 1u;
                    const std::size_t width =
                        state_count + active_parameter_count;
                    if (phase.hessian.size() != width * width) {
                        throw std::runtime_error(
                            "source-reference derivative tensor dimensions changed"
                        );
                    }
                    const std::size_t volume_coordinate = support_count;
                    const double pressure_volume_derivative =
                        -detail::GAS_CONSTANT * temperature_k
                        * phase.hessian[
                            volume_coordinate * width + volume_coordinate
                        ];
                    if (!std::isfinite(pressure_volume_derivative)
                        || pressure_volume_derivative == 0.0) {
                        throw std::runtime_error(
                            "source-reference pressure closure is singular"
                        );
                    }
                    for (std::size_t component = 0;
                         component < component_count;
                         ++component) {
                        const std::size_t component_coordinate =
                            support_coordinates[component];
                        if (component_coordinate == absent_coordinate) {
                            continue;
                        }
                        std::size_t offset = component * derivative_stride + 1u;
                        const double chemical_volume = phase.hessian[
                            component_coordinate * width + volume_coordinate
                        ];
                        if (pressure_derivatives_requested) {
                            values[offset++] = chemical_volume
                                / pressure_volume_derivative
                                - 1.0 / pressure_pa;
                        }
                        if (parameter_derivatives_requested) {
                            for (std::size_t parameter = 0;
                                 parameter < active_parameter_count;
                                 ++parameter) {
                                const std::size_t parameter_coordinate =
                                    state_count + parameter;
                                values[offset + parameter] = phase.hessian[
                                    component_coordinate * width
                                        + parameter_coordinate
                                ] - chemical_volume * phase.hessian[
                                    volume_coordinate * width
                                        + parameter_coordinate
                                ] / phase.hessian[
                                    volume_coordinate * width
                                        + volume_coordinate
                                ];
                            }
                        }
                    }
                    return values;
                };
            const auto selected_derivatives = [&] (
                const std::vector<double>& values,
                const std::vector<std::size_t>& components
            ) {
                std::vector<double> selected;
                selected.reserve(components.size() * derivative_stride);
                for (std::size_t component : components) {
                    const auto first = values.begin()
                        + static_cast<std::ptrdiff_t>(
                            component * derivative_stride
                        );
                    if (!std::all_of(
                            first,
                            first + static_cast<std::ptrdiff_t>(
                                derivative_stride
                            ),
                            [](double value) { return std::isfinite(value); }
                        )) {
                        throw std::runtime_error(
                            "source-reference derivative component is outside "
                            "the positive limiting support"
                        );
                    }
                    selected.insert(
                        selected.end(),
                        first,
                        first + static_cast<std::ptrdiff_t>(derivative_stride)
                    );
                }
                return selected;
            };
            const auto assign_derivatives = [&] (
                std::size_t target_row,
                const std::vector<double>& values,
                std::size_t selected_component
            ) {
                std::size_t offset =
                    selected_component * derivative_stride;
                fugacity[target_row] = values[offset++];
                if (pressure_derivatives_requested) {
                    fugacity_pressure_derivatives[target_row] =
                        values[offset++];
                }
                if (parameter_derivatives_requested) {
                    for (std::size_t parameter = 0;
                         parameter < active_parameter_count;
                         ++parameter) {
                        fugacity_parameter_derivatives[
                            target_row * active_parameter_count + parameter
                        ] = values[offset + parameter];
                    }
                }
            };
            double reference_molality =
                detail::REFERENCE_MOLALITIES.back();
            double convergence_error = 0.0;
            std::vector<double> neutral_values(
                reference_basis.neutral_components.size()
                    * derivative_stride,
                0.0
            );
            std::vector<std::size_t> non_pure_neutral_components;
            for (std::size_t component : reference_basis.neutral_components) {
                if (!pure_component[component]) {
                    non_pure_neutral_components.push_back(component);
                }
            }
            const auto copy_neutral_values = [&] (
                const std::vector<std::size_t>& selected_components,
                const std::vector<double>& values
            ) {
                for (std::size_t selected = 0;
                     selected < selected_components.size();
                     ++selected) {
                    const auto position = std::find(
                        reference_basis.neutral_components.begin(),
                        reference_basis.neutral_components.end(),
                        selected_components[selected]
                    );
                    if (position == reference_basis.neutral_components.end()) {
                        throw std::runtime_error(
                            "source-reference neutral basis changed"
                        );
                    }
                    const std::size_t target = static_cast<std::size_t>(
                        std::distance(
                            reference_basis.neutral_components.begin(),
                            position
                        )
                    );
                    std::copy_n(
                        values.begin()
                            + static_cast<std::ptrdiff_t>(
                                selected * derivative_stride
                            ),
                        derivative_stride,
                        neutral_values.begin()
                            + static_cast<std::ptrdiff_t>(
                                target * derivative_stride
                            )
                    );
                }
            };
            const auto reference_limit_derivatives = [&] (
                const std::vector<std::size_t>& selected_components,
                const std::vector<double>& trace_stoichiometry
            ) {
                const auto evaluate_limit =
                    [&](double molality) {
                        const std::vector<double> composition =
                            reference_composition(
                                trace_stoichiometry, molality
                            );
                        return selected_derivatives(
                            phase_fugacity_derivatives(
                                composition, "liquid", true
                            ),
                            selected_components
                        );
                    };
                return pressure_derivatives_requested
                        || parameter_derivatives_requested
                    ? detail::reference_limit_terminal_pair(evaluate_limit)
                    : detail::reference_limit_sequence(evaluate_limit);
            };
            const double neutral_trace_formula_moles = std::accumulate(
                neutral_trace_stoichiometry.begin(),
                neutral_trace_stoichiometry.end(),
                0.0
            );
            if (neutral_trace_formula_moles == 0.0) {
                copy_neutral_values(
                    non_pure_neutral_components,
                    selected_derivatives(
                        phase_fugacity_derivatives(
                            source_composition, "liquid", true
                        ),
                        non_pure_neutral_components
                    )
                );
            } else {
                const auto sequence = reference_limit_derivatives(
                    non_pure_neutral_components,
                    neutral_trace_stoichiometry
                );
                copy_neutral_values(
                    non_pure_neutral_components, sequence.values
                );
                reference_molality = sequence.molality;
                convergence_error = sequence.convergence_error;
            }
            for (std::size_t neutral = 0;
                 neutral < reference_basis.neutral_components.size();
                 ++neutral) {
                const std::size_t component =
                    reference_basis.neutral_components[neutral];
                if (!pure_component[component]) {
                    continue;
                }
                std::vector<double> composition(component_count, 0.0);
                composition[component] = 1.0;
                const std::vector<double> values = selected_derivatives(
                    phase_fugacity_derivatives(
                        composition,
                        pure_component_phase[component].c_str(),
                        false
                    ),
                    {component}
                );
                std::copy(
                    values.begin(),
                    values.end(),
                    neutral_values.begin()
                        + static_cast<std::ptrdiff_t>(neutral * derivative_stride)
                );
            }
            const double solute_scale = std::log(
                solvent_molar_mass
                * request->standard_molality_mol_per_kg
            );
            std::size_t row = 0;
            for (std::size_t component : reference_basis.neutral_components) {
                assign_derivatives(row, neutral_values, row);
                if (mixed_mole_fraction_reference) {
                    activity[row] = 0.0;
                } else if (source_composition[component] > 0.0) {
                    activity[row] = std::log(source_composition[component]);
                } else {
                    activity[row] = solute_scale;
                }
                transfer[row] = fugacity[row] + activity[row];
                ++row;
            }
            auto append_pair = [&](const NeutralReferenceBasis::IonPair& pair) {
                const std::size_t cation = pair.cation;
                const std::size_t anion = pair.anion;
                const auto& stoichiometry = pair.stoichiometry;
                auto pair_trace_stoichiometry = neutral_trace_stoichiometry;
                pair_trace_stoichiometry[cation] +=
                    static_cast<double>(stoichiometry[0]);
                pair_trace_stoichiometry[anion] +=
                    static_cast<double>(stoichiometry[1]);
                const auto sequence = reference_limit_derivatives(
                    {cation, anion}, pair_trace_stoichiometry
                );
                for (std::size_t derivative = 0;
                     derivative < derivative_stride;
                     ++derivative) {
                    const double contraction =
                        stoichiometry[0] * sequence.values[derivative]
                        + stoichiometry[1]
                            * sequence.values[derivative_stride + derivative];
                    if (derivative == 0u) {
                        fugacity[row] = contraction;
                    } else if (pressure_derivatives_requested
                               && derivative == 1u) {
                        fugacity_pressure_derivatives[row] = contraction;
                    } else {
                        const std::size_t parameter = derivative - 1u
                            - (pressure_derivatives_requested ? 1u : 0u);
                        fugacity_parameter_derivatives[
                            row * active_parameter_count + parameter
                        ] = contraction;
                    }
                }
                activity[row] = mixed_mole_fraction_reference
                    ? 0.0
                    : (static_cast<double>(stoichiometry[0])
                       + static_cast<double>(stoichiometry[1]))
                        * solute_scale;
                transfer[row] = fugacity[row] + activity[row];
                reference_molality = sequence.molality;
                convergence_error = std::max(
                    convergence_error, sequence.convergence_error
                );
                ++row;
            };
            for (const auto& pair : reference_basis.ion_pairs) {
                append_pair(pair);
            }
            const auto finite = [](const std::vector<double>& values) {
                return std::all_of(values.begin(), values.end(), [](double value) {
                    return std::isfinite(value);
                });
            };
            if (row != row_count || !finite(fugacity) || !finite(activity)
                || !finite(transfer)
                || !finite(fugacity_pressure_derivatives)
                || !finite(fugacity_parameter_derivatives)
                || !std::isfinite(convergence_error)) {
                return fail(
                    EPCSAFT_SOURCE_REFERENCE_NONFINITE,
                    "source-reference transfer is nonfinite"
                );
            }
            std::copy(
                reference_basis.rows.begin(),
                reference_basis.rows.end(),
                result->neutral_basis
            );
            std::copy(fugacity.begin(), fugacity.end(),
                      result->log_fugacity_contractions);
            std::copy(activity.begin(), activity.end(),
                      result->activity_scale_log_contractions);
            std::copy(transfer.begin(), transfer.end(),
                      result->transfer_log_contractions);
            std::copy(source_composition.begin(), source_composition.end(),
                      result->source_composition);
            if (pressure_derivatives_requested) {
                std::copy(
                    fugacity_pressure_derivatives.begin(),
                    fugacity_pressure_derivatives.end(),
                    result->pressure_derivatives_per_pa
                );
            }
            if (parameter_derivatives_requested) {
                std::copy(
                    fugacity_parameter_derivatives.begin(),
                    fugacity_parameter_derivatives.end(),
                    result->parameter_derivatives
                );
            }
            result->derivative_availability = request->required_derivatives;
            result->temperature_k = temperature_k;
            result->pressure_pa = pressure_pa;
            result->source_reference_pressure_pa =
                request->source_reference_pressure_pa;
            result->standard_molality_mol_per_kg =
                request->standard_molality_mol_per_kg;
            result->solvent_molar_mass_kg_per_mol = solvent_molar_mass;
            result->reference_limit_molality_mol_per_kg = reference_molality;
            result->reference_convergence_error = convergence_error;
            result->source_temperature_min_k = endpoint.temperature_min_k;
            result->source_temperature_max_k = endpoint.temperature_max_k;
            result->source_pressure_min_pa = endpoint.pressure_min_pa;
            result->source_pressure_max_pa = endpoint.pressure_max_pa;
            const auto copy_fingerprint = [](char* target, std::size_t capacity,
                                             const std::string& value) {
                std::snprintf(target, capacity, "%s", value.c_str());
            };
            copy_fingerprint(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                endpoint.parameter_fingerprint
            );
            copy_fingerprint(
                result->topology_fingerprint,
                sizeof(result->topology_fingerprint),
                endpoint.topology_fingerprint
            );
            copy_fingerprint(
                result->component_order_fingerprint,
                sizeof(result->component_order_fingerprint),
                endpoint.component_order_fingerprint
            );
            copy_fingerprint(
                result->reference_state_fingerprint,
                sizeof(result->reference_state_fingerprint),
                request->reference_state_fingerprint
            );
            copy_fingerprint(
                result->domain_fingerprint,
                sizeof(result->domain_fingerprint),
                endpoint.domain_fingerprint
            );
            copy_fingerprint(
                result->artifact_fingerprint,
                sizeof(result->artifact_fingerprint),
                endpoint.artifact_fingerprint
            );
            copy_fingerprint(
                result->helmholtz_basis_id,
                sizeof(result->helmholtz_basis_id),
                EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID
            );
            result->status = EPCSAFT_SOURCE_REFERENCE_OK;
            return result->status;
        } catch (const detail::DensityDomainError& error) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_DOMAIN_REJECTION, error.what()
            );
        } catch (const std::invalid_argument& error) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_INVALID_SOURCE, error.what()
            );
        } catch (const std::runtime_error& error) {
            const std::string message(error.what());
            const bool nonfinite = message.find("nonfinite") != std::string::npos;
            return fail(
                nonfinite
                    ? EPCSAFT_SOURCE_REFERENCE_NONFINITE
                    : (message.find("limit") == std::string::npos
                        ? EPCSAFT_SOURCE_REFERENCE_EVALUATION_ERROR
                        : EPCSAFT_SOURCE_REFERENCE_LIMIT_FAILURE),
                error.what()
            );
        } catch (const std::exception& error) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_EVALUATION_ERROR, error.what()
            );
        } catch (...) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_EVALUATION_ERROR,
                "unknown source-reference transfer failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        epcsaft_phase_block_result* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if (result->struct_size != sizeof(epcsaft_phase_block_result)) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;

        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "model context is null"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(amount_mol)
            || !std::isfinite(volume_m3) || temperature_k <= 0.0
            || amount_mol <= 0.0 || volume_m3 <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "temperature, amount, and volume must be positive and finite"
            );
        }

        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL,
                "model is not admitted for this phase callback"
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const auto block = endpoint.model.evaluate_pure_phase(
                temperature_k, amount_mol, volume_m3
            );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            std::copy(
                block.gradient.begin(), block.gradient.end(), result->gradient
            );
            std::copy(
                block.hessian.begin(), block.hessian.end(), result->hessian
            );
            std::copy(block.third.begin(), block.third.end(), result->third);
            result->pressure_pa = block.pressure_pa;
            result->chemical_potential_over_rt = block.chemical_potential_over_rt;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR,
                "unknown provider evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    template <typename Evaluator>
    // EqExclude: transport and buffer mechanics
    static int evaluate_mixture_common(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        std::size_t additional_coordinate_count,
        bool additional_inputs_valid,
        const char* invalid_state_message,
        epcsaft_mixture_phase_block_result* result,
        Evaluator evaluator
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if (result->struct_size != sizeof(epcsaft_mixture_phase_block_result)) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        if (context == nullptr) {
            result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "model context is null"
            );
            return result->status;
        }

        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            result->status = EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "model is not admitted for this phase callback"
            );
            return result->status;
        }
        const std::size_t coordinate_count =
            endpoint.table.component_count + additional_coordinate_count;
        if (result->coordinate_count != coordinate_count
            || result->gradient_capacity < coordinate_count
            || result->hessian_capacity < coordinate_count * coordinate_count
            || result->gradient == nullptr || result->hessian == nullptr) {
            result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "mixture result buffers do not match the model coordinate count"
            );
            return result->status;
        }

        result->status = EPCSAFT_NATIVE_STATUS_OK;
        result->helmholtz_over_rt_reference_amount = 0.0;
        result->pressure_pa = 0.0;
        result->parameter_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };

        if (amounts_mol == nullptr || amount_count != endpoint.table.component_count) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "amount buffer does not match the model component count"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(volume_m3)
            || temperature_k <= 0.0 || volume_m3 <= 0.0
            || !additional_inputs_valid
            || !std::all_of(amounts_mol, amounts_mol + amount_count, [](double amount) {
                return std::isfinite(amount) && amount > 0.0;
            })) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                invalid_state_message
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const std::vector<double> amounts(amounts_mol, amounts_mol + amount_count);
            const auto block = evaluator(endpoint.model, amounts);
            std::copy(block.gradient.begin(), block.gradient.end(), result->gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), result->hessian);
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            result->pressure_pa = block.pressure_pa;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR,
                "unknown provider evaluation failure"
            );
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_mixture(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_mixture_phase_block_result* result
    ) noexcept {
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            1,
            true,
            "temperature, amounts, and volume must be positive and finite",
            result,
            [temperature_k, volume_m3](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_mixture_phase(
                    temperature_k, amounts, volume_m3
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_electrolyte(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_mixture_phase_block_result* result
    ) noexcept {
        const auto* endpoint = static_cast<const NativeTransportEndpoint*>(context);
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            1,
            true,
            "temperature, amounts, and volume must be positive and finite",
            result,
            [temperature_k, volume_m3, endpoint](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_electrolyte_phase(
                    temperature_k,
                    amounts,
                    volume_m3,
                    endpoint->ion_mole_fraction_max
                );
            }
        );
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_molar_volume_bounds(
        const void* context,
        double temperature_k,
        const double* mole_fractions,
        size_t component_count,
        double packing_fraction_min,
        double packing_fraction_max,
        double* molar_volume_bounds_m3_per_mol,
        size_t bound_count
    ) noexcept {
        if (bound_count != 2u) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        if (context == nullptr || mole_fractions == nullptr
            || molar_volume_bounds_m3_per_mol == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
        }
        if (component_count != endpoint.table.component_count
            || !std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(packing_fraction_min)
            || !std::isfinite(packing_fraction_max)
            || packing_fraction_min <= 0.0
            || packing_fraction_min >= packing_fraction_max) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR;
        }

        try {
            const std::vector<double> composition(
                mole_fractions, mole_fractions + component_count
            );
            const auto bounds = endpoint.model.molar_volume_bounds(
                temperature_k,
                composition,
                packing_fraction_min,
                packing_fraction_max,
                endpoint.ion_mole_fraction_max
            );
            std::copy(
                bounds.begin(), bounds.end(), molar_volume_bounds_m3_per_mol
            );
            return EPCSAFT_NATIVE_STATUS_OK;
        } catch (const std::invalid_argument&) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR;
        } catch (...) {
            return EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR;
        }
    }

    // EqExclude: transport and buffer mechanics
    static int evaluate_packing_fraction(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double* packing_fraction,
        double* gradient,
        size_t gradient_capacity,
        double* hessian,
        size_t hessian_capacity
    ) noexcept {
        if (context == nullptr || amounts_mol == nullptr || packing_fraction == nullptr
            || gradient == nullptr || hessian == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
        }
        const std::size_t support_count = static_cast<std::size_t>(std::count_if(
            amounts_mol,
            amounts_mol + amount_count,
            [](double amount) { return amount > 0.0; }
        ));
        const std::size_t coordinate_count = support_count + 1;
        if (gradient_capacity != coordinate_count
            || hessian_capacity != coordinate_count * coordinate_count) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        if (amount_count != endpoint.table.component_count
            || !std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(volume_m3) || volume_m3 <= 0.0
            || !std::all_of(amounts_mol, amounts_mol + amount_count, [](double amount) {
                return std::isfinite(amount) && amount >= 0.0;
            })) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR;
        }

        try {
            const std::vector<double> amounts(amounts_mol, amounts_mol + amount_count);
            const auto block = endpoint.model.evaluate_packing_fraction_parameters(
                temperature_k,
                amounts,
                volume_m3,
                endpoint.ion_mole_fraction_max,
                {}
            );
            *packing_fraction = block.packing_fraction;
            std::copy(block.gradient.begin(), block.gradient.end(), gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), hessian);
            return EPCSAFT_NATIVE_STATUS_OK;
        } catch (const std::invalid_argument&) {
            return EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR;
        } catch (const std::logic_error&) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
        } catch (...) {
            return EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR;
        }
    }
};

}  // namespace

// EqExclude: transport and buffer mechanics
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
    const std::string& artifact_fingerprint,
    const std::vector<NativeActiveCoordinate>& active_parameters
) {
    if (!model.supports_pure_transport()
        && !model.supports_mixture_transport()
        && !model.supports_fixed_topology_associating_mixture()
        && !model.supports_electrolyte_transport()) {
        throw std::invalid_argument(
            "native transport requires a supported phase-block model"
        );
    }
    if (!sha256_fingerprint_valid(parameter_fingerprint)) {
        throw std::invalid_argument("invalid parameter fingerprint");
    }
    for (const auto& fingerprint : {
             std::cref(topology_fingerprint),
             std::cref(component_order_fingerprint),
             std::cref(domain_fingerprint),
             std::cref(artifact_fingerprint),
         }) {
        if (!fingerprint.get().empty()
            && !sha256_fingerprint_valid(fingerprint.get())) {
            throw std::invalid_argument("invalid transport fingerprint");
        }
    }
    if (component_ids.size() != model.component_count()
        || std::any_of(
            component_ids.begin(), component_ids.end(), [](const auto& identifier) {
                return identifier.empty();
            }
        )) {
        throw std::invalid_argument("invalid component identifiers");
    }
    for (std::size_t left = 0; left < component_ids.size(); ++left) {
        for (std::size_t right = left + 1; right < component_ids.size(); ++right) {
            if (component_ids[left] == component_ids[right]) {
                throw std::invalid_argument("duplicate component identifier");
            }
        }
    }
    const auto valid_optional_domain = [](double value) {
        return std::isnan(value) || (std::isfinite(value) && value > 0.0);
    };
    if (!valid_optional_domain(temperature_min_k)
        || !valid_optional_domain(temperature_max_k)
        || (std::isfinite(temperature_min_k)
            && std::isfinite(temperature_max_k)
            && temperature_min_k > temperature_max_k)
        || (!std::isnan(ion_mole_fraction_max)
            && (!std::isfinite(ion_mole_fraction_max)
                || ion_mole_fraction_max <= 0.0
                || ion_mole_fraction_max > 1.0))) {
        throw std::invalid_argument("invalid transport temperature domain");
    }
    if (!valid_optional_domain(pressure_min_pa)
        || !valid_optional_domain(pressure_max_pa)
        || (std::isfinite(pressure_min_pa) && std::isfinite(pressure_max_pa)
            && pressure_min_pa > pressure_max_pa)) {
        throw std::invalid_argument("invalid transport pressure domain");
    }
    if (!active_parameters.empty()) {
        static_cast<void>(model.with_active_coordinates(active_parameters));
        std::vector<std::string> identities;
        identities.reserve(active_parameters.size());
        for (const auto& coordinate : active_parameters) {
            if (coordinate.identity.empty()
                || std::find(
                       identities.begin(), identities.end(), coordinate.identity
                   ) != identities.end()) {
                throw std::invalid_argument(
                    "active parameter identities must be nonempty and unique"
                );
            }
            identities.push_back(coordinate.identity);
        }
    }
    auto* endpoint = new NativeTransportEndpoint(
        model,
        parameter_fingerprint,
        temperature_min_k,
        temperature_max_k,
        component_ids,
        ion_mole_fraction_max,
        pressure_min_pa,
        pressure_max_pa,
        topology_fingerprint,
        component_order_fingerprint,
        domain_fingerprint,
        artifact_fingerprint,
        active_parameters
    );
    return &endpoint->table;
}


// EqExclude: transport and buffer mechanics
void destroy_native_transport(const epcsaft_native_transport* transport) noexcept {
    if (transport != nullptr) {
        delete static_cast<const NativeTransportEndpoint*>(transport->model_context);
    }
}

}  // namespace epcsaft
