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
    const std::vector<int>& charges
) {
    NeutralReferenceBasis result;
    const std::size_t components = component_ids.size();
    if (components < 3 || components != charges.size()
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
    std::vector<std::size_t> cations;
    std::vector<std::size_t> anions;
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
        } else if (charges[component] > 0) {
            cations.push_back(component);
        } else {
            anions.push_back(component);
        }
    }
    if (result.neutral_components.empty() || cations.empty()
        || anions.empty()) {
        return {};
    }
    result.reference_composition.assign(components, 0.0);
    const double neutral_fraction =
        1.0 / static_cast<double>(result.neutral_components.size());
    for (std::size_t component : result.neutral_components) {
        result.reference_composition[component] = neutral_fraction;
    }
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
    auto append_pair = [&](std::size_t cation, std::size_t anion) {
        const int divisor = std::gcd(
            std::abs(charges[cation]), std::abs(charges[anion])
        );
        const std::array<int, 2> stoichiometry{
            std::abs(charges[anion]) / divisor,
            std::abs(charges[cation]) / divisor,
        };
        result.ion_pairs.push_back({cation, anion, stoichiometry});
        append_row(
            cation,
            static_cast<double>(stoichiometry[0]),
            anion,
            static_cast<double>(stoichiometry[1])
        );
    };
    for (std::size_t anion : anions) {
        append_pair(cations.front(), anion);
    }
    for (std::size_t index = 1; index < cations.size(); ++index) {
        append_pair(cations[index], anions.front());
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
    NeutralReferenceBasis neutral_reference_basis;
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
        std::string artifact_identity
    )
        : model(source_model),
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
              || source_model.supports_electrolyte_transport()
          ),
          neutral_reference_admitted(
              make_neutral_reference_basis(
                  identifiers, source_model.charge_numbers()
              ).row_count() == identifiers.size() - 1u
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
          neutral_reference_basis(make_neutral_reference_basis(
              component_ids, source_model.charge_numbers()
          )),
          table{
              this,
              source_model.supports_pure_transport()
                  ? &NativeTransportEndpoint::evaluate
                  : nullptr,
              component_ids.size(),
              source_model.supports_mixture_transport()
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
              source_model.supports_electrolyte_transport()
                      && !topology_fingerprint.empty()
                  ? &NativeTransportEndpoint::evaluate_reacting_phase_parameters
                  : nullptr,
              source_model.supports_electrolyte_transport()
                      && !topology_fingerprint.empty()
                  ? &NativeTransportEndpoint::evaluate_inverse_packing_geometry
                  : nullptr,
              neutral_reference_admitted
                  ? &NativeTransportEndpoint::evaluate_source_reference_transfer
                  : nullptr,
              component_id_pointers.data(),
              component_charges.data(),
              parameter_fingerprint.c_str(),
              topology_fingerprint.c_str(),
              temperature_min_k,
              temperature_max_k,
              ion_mole_fraction_max,
              neutral_reference_admitted
                  ? neutral_reference_basis.row_count()
                  : 0u,
          } {}

    // EqExclude: transport and buffer mechanics
    static std::vector<NativeActiveParameter> active_parameters(
        const NativeTransportEndpoint& endpoint,
        const epcsaft_active_parameter_request* requests,
        std::size_t request_count
    ) {
        if (request_count == 0) {
            return {};
        }
        if (requests == nullptr) {
            throw std::invalid_argument(
                "active parameter request pointer is null"
            );
        }
        std::vector<NativeActiveParameter> result;
        result.reserve(request_count);
        for (std::size_t index = 0; index < request_count; ++index) {
            const auto& request = requests[index];
            if (request.parameter_identity
                    != EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT
                || request.component_index < 0
                || static_cast<std::size_t>(request.component_index)
                    >= endpoint.component_ids.size()
                || request.unit == nullptr || !std::isfinite(request.value)
                || request.value <= 0.0) {
                throw std::invalid_argument(
                    "active parameter request is malformed"
                );
            }
            ActiveComponentParameter family;
            const char* expected_unit = nullptr;
            if (request.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER) {
                family = ActiveComponentParameter::SegmentDiameter;
                expected_unit = "angstrom";
            } else if (
                request.parameter_family
                == EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K
            ) {
                family = ActiveComponentParameter::DispersionEnergyOverK;
                expected_unit = "kelvin";
            } else {
                throw std::logic_error(
                    "active parameter family is unsupported"
                );
            }
            if (std::strcmp(request.unit, expected_unit) != 0) {
                throw std::invalid_argument(
                    "active parameter unit does not match its family"
                );
            }
            if (std::any_of(
                    result.begin(),
                    result.end(),
                    [&](const NativeActiveParameter& item) {
                        return item.family == family
                            && item.component
                                == static_cast<std::size_t>(
                                    request.component_index
                                );
                    }
                )) {
                throw std::invalid_argument(
                    "active parameter request is duplicated"
                );
            }
            result.push_back({
                family,
                static_cast<std::size_t>(request.component_index),
                request.value,
            });
        }
        return result;
    }

    // EqID: phase_packing_fraction
    // EqExclude: transport and buffer mechanics
    static int evaluate_inverse_packing_geometry(
        const void* context,
        const char* expected_fingerprint,
        const char* expected_topology_fingerprint,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double log_packing_fraction,
        const epcsaft_active_parameter_request* requests,
        size_t request_count,
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
        if (context == nullptr || expected_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr
            || amounts_mol == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "inverse packing required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeTransportEndpoint*>(context);
        if (endpoint.topology_fingerprint.empty()
            || std::strcmp(
                expected_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL,
                "inverse packing model or topology is not admitted"
            );
        }
        if (request_count > 2 * endpoint.table.component_count) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "inverse packing active parameter count exceeds model maximum"
            );
        }
        const std::size_t coordinate_count =
            endpoint.table.component_count + 1 + request_count;
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
                    return std::isfinite(amount) && amount > 0.0;
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
            const auto parameters =
                active_parameters(endpoint, requests, request_count);
            const std::vector<double> amounts(
                amounts_mol, amounts_mol + amount_count
            );
            const auto block =
                endpoint.model.evaluate_inverse_packing_geometry(
                    temperature_k,
                    amounts,
                    log_packing_fraction,
                    parameters
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
    // EqID: fixed_pressure_total_derivative
    static int evaluate_reacting_phase_parameters(
        const void* context,
        const char* expected_fingerprint,
        const char* expected_topology_fingerprint,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double packing_fraction_min,
        double packing_fraction_max,
        const epcsaft_active_parameter_request* requests,
        size_t request_count,
        epcsaft_reacting_phase_parameter_result* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        if (result->struct_size
            != sizeof(epcsaft_reacting_phase_parameter_result)) {
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
        if (context == nullptr || expected_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr
            || amounts_mol == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT,
                "reacting phase required pointer is null"
            );
        }
        const auto& endpoint =
            *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted
            || !endpoint.model.supports_electrolyte_transport()
            || endpoint.topology_fingerprint.empty()
            || std::strcmp(
                expected_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL,
                "reacting phase model identity is unsupported"
            );
        }
        if (!std::isfinite(temperature_k)
            || temperature_k < endpoint.temperature_min_k
            || temperature_k > endpoint.temperature_max_k) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR,
                "reacting phase temperature is outside the source domain"
            );
        }
        const std::size_t state_count = endpoint.component_ids.size() + 1;
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
                != endpoint.component_ids.size() * request_count
            || result->state_gradient == nullptr
            || result->state_hessian == nullptr
            || result->parameter_gradient == nullptr
            || result->state_parameter_derivatives == nullptr
            || result->packing_parameter_derivatives == nullptr
            || result->packing_state_gradient == nullptr
            || result->packing_state_hessian == nullptr
            || result->molar_volume_bound_parameter_derivatives_m3_per_mol
                == nullptr
            || result->pressure_parameter_derivatives_pa == nullptr
            || result->chemical_potential_parameter_derivatives_over_rt
                == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_BUFFER_ERROR,
                "reacting phase buffers do not match immutable dimensions"
            );
        }
        try {
            const auto parameters = active_parameters(
                endpoint, requests, request_count
            );
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
                    parameters
                );
            const auto phase = endpoint.model.evaluate_reacting_phase_parameters(
                temperature_k,
                amounts,
                volume_m3,
                endpoint.ion_mole_fraction_max,
                parameters
            );
            const auto volume_bounds =
                endpoint.model.molar_volume_bounds_parameters(
                    temperature_k,
                    composition,
                    packing_fraction_min,
                    packing_fraction_max,
                    endpoint.ion_mole_fraction_max,
                    parameters
                );
            const std::size_t joint_count = state_count + request_count;
            if (phase.gradient.size() != joint_count
                || phase.hessian.size() != joint_count * joint_count
                || packing.gradient.size() != joint_count
                || packing.hessian.size() != joint_count * joint_count) {
                throw std::runtime_error(
                    "reacting phase tensor dimensions changed"
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
            std::vector<double> chemical_parameter(
                amount_count * request_count
            );
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
                     component < amount_count;
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
                    "reacting phase result contains nonfinite values"
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
                "unknown reacting phase parameter failure"
            );
        }
    }

    // EqID: declared_solvent_source_reference_transfer
    static int evaluate_source_reference_transfer(
        const void* context,
        const char* expected_fingerprint,
        const char* expected_topology_fingerprint,
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
        if (context == nullptr || expected_fingerprint == nullptr
            || expected_topology_fingerprint == nullptr || request == nullptr) {
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
            || endpoint.topology_fingerprint.empty()
            || std::strcmp(
                expected_fingerprint, endpoint.parameter_fingerprint.c_str()
            ) != 0
            || std::strcmp(
                expected_topology_fingerprint,
                endpoint.topology_fingerprint.c_str()
            ) != 0) {
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
        const auto& reference_basis = endpoint.neutral_reference_basis;
        const std::size_t row_count = reference_basis.row_count();
        if (request->component_count != component_count
            || request->source_composition_count != component_count
            || request->component_ids == nullptr
            || request->source_solvent_mole_fractions == nullptr) {
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
        if (request->reference_state_id == nullptr
            || request->reference_state_id[0] == '\0'
            || request->activity_convention_id == nullptr
            || std::strcmp(
                request->activity_convention_id,
                EPCSAFT_SOURCE_REFERENCE_ACTIVITY_MOLALITY_INFINITE_DILUTION
            ) != 0
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
        if (request->required_derivatives
            != 0u) {
            return fail(
                EPCSAFT_SOURCE_REFERENCE_DERIVATIVE_UNAVAILABLE,
                "required source-reference derivative is unavailable"
            );
        }
        if (result->component_count != component_count
            || result->neutral_basis_row_count != row_count
            || result->neutral_basis_capacity != row_count * component_count
            || result->contraction_capacity != row_count
            || result->source_composition_capacity != component_count
            || result->pressure_derivative_capacity != 0
            || result->neutral_basis == nullptr
            || result->log_fugacity_contractions == nullptr
            || result->activity_scale_log_contractions == nullptr
            || result->transfer_log_contractions == nullptr
            || result->source_composition == nullptr) {
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
            double reference_molality =
                detail::REFERENCE_MOLALITIES.back();
            double convergence_error = 0.0;
            std::vector<std::size_t> neutral_solutes;
            for (std::size_t component : reference_basis.neutral_components) {
                if (source_composition[component] == 0.0) {
                    neutral_solutes.push_back(component);
                }
            }
            std::vector<double> neutral_log_fugacity(component_count, 0.0);
            if (neutral_solutes.empty()) {
                const auto state = endpoint.model.evaluate_pressure(
                    temperature_k, pressure_pa, source_composition, "liquid"
                );
                if (!state.log_fugacity_coefficient.has_value()) {
                    throw std::runtime_error(
                        "source-reference fugacity is unavailable"
                    );
                }
                neutral_log_fugacity = *state.log_fugacity_coefficient;
            } else {
                const auto sequence = detail::reference_limit_sequence(
                    [&](double molality) {
                        const double solvent_moles = 1.0 / solvent_molar_mass;
                        const double total_moles = solvent_moles
                            + molality
                                * static_cast<double>(neutral_solutes.size());
                        std::vector<double> composition(component_count, 0.0);
                        for (std::size_t component = 0;
                             component < component_count;
                             ++component) {
                            composition[component] =
                                source_composition[component] * solvent_moles
                                / total_moles;
                        }
                        for (std::size_t component : neutral_solutes) {
                            composition[component] = molality / total_moles;
                        }
                        const auto state = endpoint.model.evaluate_pressure(
                            temperature_k, pressure_pa, composition, "liquid"
                        );
                        if (!state.log_fugacity_coefficient.has_value()) {
                            throw std::runtime_error(
                                "neutral-solute reference fugacity is unavailable"
                            );
                        }
                        std::vector<double> values;
                        values.reserve(
                            reference_basis.neutral_components.size()
                        );
                        for (std::size_t component
                             : reference_basis.neutral_components) {
                            values.push_back(
                                (*state.log_fugacity_coefficient)[component]
                            );
                        }
                        return values;
                    }
                );
                for (std::size_t index = 0;
                     index < reference_basis.neutral_components.size();
                     ++index) {
                    neutral_log_fugacity[
                        reference_basis.neutral_components[index]
                    ] = sequence.values[index];
                }
                reference_molality = sequence.molality;
                convergence_error = sequence.convergence_error;
            }
            const double solute_scale = std::log(
                solvent_molar_mass
                * request->standard_molality_mol_per_kg
            );
            std::size_t row = 0;
            for (std::size_t component : reference_basis.neutral_components) {
                fugacity[row] = neutral_log_fugacity[component];
                if (source_composition[component] > 0.0) {
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
                const auto sequence = detail::reference_limit_sequence(
                    [&](double molality) {
                        const double solvent_moles = 1.0 / solvent_molar_mass;
                        const double added_moles = molality
                            * static_cast<double>(
                                stoichiometry[0] + stoichiometry[1]
                                + neutral_solutes.size()
                            );
                        const double total_moles = solvent_moles + added_moles;
                        std::vector<double> composition(component_count, 0.0);
                        for (std::size_t component = 0;
                             component < component_count;
                             ++component) {
                            composition[component] =
                                source_composition[component] * solvent_moles
                                / total_moles;
                        }
                        for (std::size_t component : neutral_solutes) {
                            composition[component] = molality / total_moles;
                        }
                        composition[cation] =
                            static_cast<double>(stoichiometry[0]) * molality
                            / total_moles;
                        composition[anion] =
                            static_cast<double>(stoichiometry[1]) * molality
                            / total_moles;
                        const auto state = endpoint.model.evaluate_pressure(
                            temperature_k, pressure_pa, composition, "liquid"
                        );
                        if (!state.log_fugacity_coefficient.has_value()) {
                            throw std::runtime_error(
                                "ionic source-reference fugacity is unavailable"
                            );
                        }
                        return std::vector<double>{
                            (*state.log_fugacity_coefficient)[cation],
                            (*state.log_fugacity_coefficient)[anion],
                        };
                    }
                );
                fugacity[row] =
                    stoichiometry[0] * sequence.values[0]
                    + stoichiometry[1] * sequence.values[1];
                activity[row] =
                    static_cast<double>(stoichiometry[0] + stoichiometry[1])
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
                || !finite(transfer) || !std::isfinite(convergence_error)) {
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
        const char* expected_parameter_fingerprint,
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
        if (context == nullptr || expected_parameter_fingerprint == nullptr
            || mole_fractions == nullptr
            || molar_volume_bounds_m3_per_mol == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
        }
        if (std::strcmp(
                expected_parameter_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0) {
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
    ) noexcept {
        if (context == nullptr || expected_parameter_fingerprint == nullptr
            || amounts_mol == nullptr || packing_fraction == nullptr
            || gradient == nullptr || hessian == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT;
        }
        const auto& endpoint = *static_cast<const NativeTransportEndpoint*>(context);
        if (!endpoint.phase_callbacks_admitted) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
        }
        const std::size_t coordinate_count = endpoint.table.component_count + 1;
        if (gradient_capacity != coordinate_count
            || hessian_capacity != coordinate_count * coordinate_count) {
            return EPCSAFT_NATIVE_STATUS_BUFFER_ERROR;
        }
        if (std::strcmp(
                expected_parameter_fingerprint,
                endpoint.parameter_fingerprint.c_str()
            ) != 0) {
            return EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL;
        }
        if (amount_count != endpoint.table.component_count
            || !std::isfinite(temperature_k) || temperature_k <= 0.0
            || !std::isfinite(volume_m3) || volume_m3 <= 0.0
            || !std::all_of(amounts_mol, amounts_mol + amount_count, [](double amount) {
                return std::isfinite(amount) && amount > 0.0;
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
            const auto block = endpoint.model.evaluate_packing_fraction(
                temperature_k,
                amounts,
                volume_m3,
                endpoint.ion_mole_fraction_max
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
    const std::string& artifact_fingerprint
) {
    if (!model.supports_pure_transport()
        && !model.supports_mixture_transport()
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
        artifact_fingerprint
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
