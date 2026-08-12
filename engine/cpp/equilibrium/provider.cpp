#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

std::string decode_provider_char_array(
    const char* value,
    std::size_t capacity,
    const char* field_name
) {
    const void* terminator = std::memchr(value, '\0', capacity);
    if (terminator == nullptr) {
        throw std::invalid_argument(std::string(field_name) + " is missing a NUL terminator");
    }
    const auto* end = static_cast<const char*>(terminator);
    return std::string(value, static_cast<std::size_t>(end - value));
}

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

MixturePhaseEvaluation evaluate_mixture_callback(
    const epcsaft_native_transport& transport,
    epcsaft_evaluate_mixture_phase callback,
    const std::string& fingerprint,
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) {
    if (transport.component_count == 0 || amounts_mol.size() != transport.component_count) {
        throw std::invalid_argument("provider mixture component count mismatch");
    }
    if (callback == nullptr) {
        throw std::invalid_argument("provider capsule is missing the mixture phase contract");
    }
    const std::size_t coordinate_count = transport.component_count + 1;
    std::vector<double> gradient(coordinate_count, 0.0);
    std::vector<double> hessian(coordinate_count * coordinate_count, 0.0);
    epcsaft_mixture_phase_block_result phase{};
    phase.struct_size = sizeof(phase);
    phase.coordinate_count = coordinate_count;
    phase.gradient_capacity = gradient.size();
    phase.hessian_capacity = hessian.size();
    phase.gradient = gradient.data();
    phase.hessian = hessian.data();
    const int status = callback(
        transport.model_context,
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        &phase
    );
    if (phase.struct_size != sizeof(epcsaft_mixture_phase_block_result)) {
        throw std::invalid_argument("provider mixture result struct size mismatch");
    }
    if (status != phase.status) {
        throw std::runtime_error("provider mixture evaluation returned inconsistent status values");
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK) {
        const std::string error = decode_provider_char_array(
            phase.error,
            sizeof(phase.error),
            "provider mixture error"
        );
        throw std::domain_error("provider mixture phase evaluation failed: " + error);
    }
    if (phase.coordinate_count != coordinate_count
        || phase.gradient_capacity != gradient.size()
        || phase.hessian_capacity != hessian.size()
        || phase.gradient != gradient.data()
        || phase.hessian != hessian.data()) {
        throw std::invalid_argument("provider mixture result buffer contract changed");
    }
    require_finite(phase.helmholtz_over_rt_reference_amount, "provider mixture value");
    require_finite(phase.pressure_pa, "provider mixture pressure");
    if (!std::all_of(gradient.begin(), gradient.end(), [](double value) {
            return std::isfinite(value);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("provider mixture tensors must be finite");
    }
    const std::string parameter_fingerprint = decode_provider_char_array(
        phase.parameter_fingerprint,
        sizeof(phase.parameter_fingerprint),
        "provider mixture parameter fingerprint"
    );
    if (parameter_fingerprint != fingerprint) {
        throw std::invalid_argument(
            "provider mixture result fingerprint does not match the requested model"
        );
    }
    return {
        phase.helmholtz_over_rt_reference_amount,
        std::move(gradient),
        std::move(hessian),
        phase.pressure_pa,
        std::move(parameter_fingerprint),
    };
}

}  // namespace

ProviderContext::ProviderContext(const epcsaft_native_transport& transport)
    : transport_(transport),
      fingerprint_(transport.parameter_fingerprint == nullptr
          ? ""
          : transport.parameter_fingerprint) {
    if (fingerprint_.empty() || transport_.model_context == nullptr) {
        throw std::invalid_argument("provider capsule model identity is incomplete");
    }
}

std::string ProviderContext::parameterized_topology_fingerprint() const {
    if (transport_.topology_fingerprint == nullptr
        || transport_.topology_fingerprint[0] == '\0') {
        throw std::invalid_argument(
            "Provider parameterized-phase topology metadata is incompatible"
        );
    }
    return transport_.topology_fingerprint;
}

PhaseEvaluation ProviderContext::evaluate(
    double temperature_k,
    double amount_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    epcsaft_phase_block_result phase{};
    phase.struct_size = sizeof(phase);
    const int status = transport_.evaluate_pure_phase(
        transport_.model_context,
        temperature_k,
        amount_mol,
        volume_m3,
        &phase
    );
    if (phase.struct_size != sizeof(epcsaft_phase_block_result)) {
        throw std::invalid_argument("provider phase result struct size mismatch");
    }
    if (status != phase.status) {
        throw std::runtime_error("provider phase evaluation returned inconsistent status values");
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK) {
        const std::string error = decode_provider_char_array(
            phase.error,
            sizeof(phase.error),
            "provider error"
        );
        throw std::domain_error("provider phase evaluation failed: " + error);
    }
    const std::string parameter_fingerprint = decode_provider_char_array(
        phase.parameter_fingerprint,
        sizeof(phase.parameter_fingerprint),
        "provider parameter fingerprint"
    );
    if (parameter_fingerprint != fingerprint_) {
        throw std::invalid_argument(
            "provider phase result fingerprint does not match the requested model"
        );
    }
    return {amount_mol, volume_m3, phase, std::move(parameter_fingerprint)};
}

MixturePhaseEvaluation ProviderContext::evaluate_mixture(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    return evaluate_mixture_callback(
        transport_, transport_.evaluate_mixture_phase, fingerprint_, temperature_k, amounts_mol, volume_m3
    );
}

MixturePhaseEvaluation ProviderContext::evaluate_electrolyte(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    return evaluate_mixture_callback(
        transport_, transport_.evaluate_electrolyte_phase, fingerprint_, temperature_k, amounts_mol, volume_m3
    );
}

std::array<double, 2> ProviderContext::evaluate_molar_volume_bounds(
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double packing_fraction_min,
    double packing_fraction_max
) const {
    ++evaluation_count_;
    if (transport_.evaluate_molar_volume_bounds == nullptr
        || mole_fractions.size() != transport_.component_count) {
        throw std::invalid_argument(
            "provider capsule is missing the molar-volume domain contract"
        );
    }
    require_finite(temperature_k, "temperature");
    require_finite(packing_fraction_min, "minimum packing fraction");
    require_finite(packing_fraction_max, "maximum packing fraction");
    if (temperature_k <= 0.0 || packing_fraction_min <= 0.0
        || packing_fraction_max <= packing_fraction_min) {
        throw std::invalid_argument(
            "temperature and packing-fraction bounds must be positive and ordered"
        );
    }
    if (!std::all_of(
            mole_fractions.begin(),
            mole_fractions.end(),
            [](double value) {
                return std::isfinite(value) && value >= 0.0;
            }
        )) {
        throw std::invalid_argument(
            "mole fractions must be finite and nonnegative"
        );
    }
    std::array<double, 2> bounds{};
    const int status = transport_.evaluate_molar_volume_bounds(
        transport_.model_context,
        temperature_k,
        mole_fractions.data(),
        mole_fractions.size(),
        packing_fraction_min,
        packing_fraction_max,
        bounds.data(),
        bounds.size()
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK) {
        throw std::domain_error(
            "provider molar-volume domain/fingerprint evaluation failed with status "
            + std::to_string(status)
        );
    }
    if (!std::isfinite(bounds[0]) || !std::isfinite(bounds[1])
        || bounds[0] <= 0.0 || bounds[1] <= bounds[0]) {
        throw std::invalid_argument(
            "provider molar-volume bounds must be finite, positive, and ordered"
        );
    }
    return bounds;
}

PackingFractionEvaluation ProviderContext::evaluate_packing_fraction(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    if (transport_.evaluate_packing_fraction == nullptr
        || amounts_mol.size() != transport_.component_count) {
        throw std::invalid_argument("provider capsule is missing the packing-fraction contract");
    }
    require_finite(temperature_k, "temperature");
    require_finite(volume_m3, "volume");
    if (temperature_k <= 0.0 || volume_m3 <= 0.0
        || !std::all_of(amounts_mol.begin(), amounts_mol.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument(
            "packing-fraction state must be finite, nonnegative, and positive-volume"
        );
    }

    const std::size_t coordinate_count = amounts_mol.size() + 1;
    PackingFractionEvaluation result;
    result.gradient.resize(coordinate_count);
    result.hessian.resize(coordinate_count * coordinate_count);
    const int status = transport_.evaluate_packing_fraction(
        transport_.model_context,
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        &result.value,
        result.gradient.data(),
        result.gradient.size(),
        result.hessian.data(),
        result.hessian.size()
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK) {
        throw std::domain_error(
            "provider packing-fraction evaluation failed with status "
            + std::to_string(status)
        );
    }
    require_finite(result.value, "provider packing fraction");
    if (!std::all_of(result.gradient.begin(), result.gradient.end(), [](double value) {
            return std::isfinite(value);
        })
        || !std::all_of(result.hessian.begin(), result.hessian.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument("provider packing-fraction tensors must be finite");
    }
    return result;
}

ProviderActiveParameterSet ProviderContext::active_parameters() const {
    ProviderActiveParameterSet result;
    if (transport_.active_parameter_count != 0
        && (transport_.active_parameter_identities == nullptr
            || transport_.evaluate_parameterized_phase == nullptr)) {
        throw std::invalid_argument(
            "Provider active-parameter transport is incomplete"
        );
    }
    result.identities.reserve(transport_.active_parameter_count);
    for (std::size_t index = 0; index < transport_.active_parameter_count; ++index) {
        const char* identity = transport_.active_parameter_identities[index];
        if (identity == nullptr || identity[0] == '\0'
            || std::find(result.identities.begin(), result.identities.end(), identity)
                != result.identities.end()) {
            throw std::invalid_argument(
                "Provider active-parameter identity metadata is invalid"
            );
        }
        result.identities.emplace_back(identity);
    }
    return result;
}

ParameterizedPhaseEvaluation ProviderContext::evaluate_parameterized_phase(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    double packing_fraction_min,
    double packing_fraction_max,
    const ProviderActiveParameterSet& active_parameters
) const {
    ++evaluation_count_;
    if (active_parameters.identities.empty()
        || transport_.evaluate_parameterized_phase == nullptr
        || amounts_mol.size() != transport_.component_count) {
        throw std::invalid_argument(
            "Provider atomic parameterized-phase contract is incomplete"
        );
    }
    const std::size_t state_count = amounts_mol.size() + 1;
    const std::size_t parameter_count = active_parameters.identities.size();
    ParameterizedPhaseEvaluation result;
    result.phase.gradient.resize(state_count);
    result.phase.hessian.resize(state_count * state_count);
    result.packing.gradient.resize(state_count);
    result.packing.hessian.resize(state_count * state_count);
    result.parameter_derivatives.resize(parameter_count);
    result.state_parameter_derivatives.resize(state_count * parameter_count);
    result.packing_parameter_derivatives.resize(parameter_count);
    result.molar_volume_bound_parameter_derivatives_m3_per_mol.resize(
        2 * parameter_count
    );
    result.pressure_parameter_derivatives_pa.resize(parameter_count);
    result.chemical_potential_parameter_derivatives_over_rt.resize(
        amounts_mol.size() * parameter_count
    );
    epcsaft_parameterized_phase_result native{};
    native.struct_size = sizeof(native);
    native.state_coordinate_count = state_count;
    native.active_parameter_count = parameter_count;
    native.state_gradient_capacity = result.phase.gradient.size();
    native.state_hessian_capacity = result.phase.hessian.size();
    native.parameter_gradient_capacity = result.parameter_derivatives.size();
    native.state_parameter_capacity = result.state_parameter_derivatives.size();
    native.packing_parameter_capacity = result.packing_parameter_derivatives.size();
    native.molar_volume_bound_parameter_capacity =
        result.molar_volume_bound_parameter_derivatives_m3_per_mol.size();
    native.pressure_parameter_capacity =
        result.pressure_parameter_derivatives_pa.size();
    native.chemical_potential_parameter_capacity =
        result.chemical_potential_parameter_derivatives_over_rt.size();
    native.state_gradient = result.phase.gradient.data();
    native.state_hessian = result.phase.hessian.data();
    native.parameter_gradient = result.parameter_derivatives.data();
    native.state_parameter_derivatives =
        result.state_parameter_derivatives.data();
    native.packing_parameter_derivatives =
        result.packing_parameter_derivatives.data();
    native.molar_volume_bound_parameter_derivatives_m3_per_mol =
        result.molar_volume_bound_parameter_derivatives_m3_per_mol.data();
    native.pressure_parameter_derivatives_pa =
        result.pressure_parameter_derivatives_pa.data();
    native.chemical_potential_parameter_derivatives_over_rt =
        result.chemical_potential_parameter_derivatives_over_rt.data();
    native.packing_state_gradient_capacity = result.packing.gradient.size();
    native.packing_state_hessian_capacity = result.packing.hessian.size();
    native.packing_state_gradient = result.packing.gradient.data();
    native.packing_state_hessian = result.packing.hessian.data();
    const int status = transport_.evaluate_parameterized_phase(
        transport_.model_context,
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        volume_m3,
        packing_fraction_min,
        packing_fraction_max,
        &native
    );
    if (status != native.status) {
        throw std::runtime_error(
            "Provider parameterized-phase evaluation returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK) {
        throw std::domain_error(
            "Provider parameterized-phase evaluation failed: "
            + decode_provider_char_array(
                native.error,
                sizeof(native.error),
                "Provider parameterized-phase error"
            )
        );
    }
    if (native.struct_size != sizeof(native)
        || native.state_coordinate_count != state_count
        || native.active_parameter_count != parameter_count
        || native.state_gradient_capacity != result.phase.gradient.size()
        || native.state_hessian_capacity != result.phase.hessian.size()
        || native.parameter_gradient_capacity != result.parameter_derivatives.size()
        || native.state_parameter_capacity
            != result.state_parameter_derivatives.size()
        || native.packing_parameter_capacity
            != result.packing_parameter_derivatives.size()
        || native.molar_volume_bound_parameter_capacity
            != result.molar_volume_bound_parameter_derivatives_m3_per_mol.size()
        || native.pressure_parameter_capacity
            != result.pressure_parameter_derivatives_pa.size()
        || native.chemical_potential_parameter_capacity
            != result.chemical_potential_parameter_derivatives_over_rt.size()
        || native.packing_state_gradient_capacity
            != result.packing.gradient.size()
        || native.packing_state_hessian_capacity
            != result.packing.hessian.size()
        || native.state_gradient != result.phase.gradient.data()
        || native.state_hessian != result.phase.hessian.data()
        || native.parameter_gradient != result.parameter_derivatives.data()
        || native.state_parameter_derivatives
            != result.state_parameter_derivatives.data()
        || native.packing_parameter_derivatives
            != result.packing_parameter_derivatives.data()
        || native.molar_volume_bound_parameter_derivatives_m3_per_mol
            != result.molar_volume_bound_parameter_derivatives_m3_per_mol.data()
        || native.pressure_parameter_derivatives_pa
            != result.pressure_parameter_derivatives_pa.data()
        || native.chemical_potential_parameter_derivatives_over_rt
            != result.chemical_potential_parameter_derivatives_over_rt.data()
        || native.packing_state_gradient != result.packing.gradient.data()
        || native.packing_state_hessian != result.packing.hessian.data()) {
        throw std::invalid_argument(
            "Provider parameterized-phase result buffer contract changed"
        );
    }
    const auto finite = [](const std::vector<double>& values) {
        return std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        });
    };
    if (!finite(result.phase.gradient) || !finite(result.phase.hessian)
        || !finite(result.packing.gradient) || !finite(result.packing.hessian)
        || !finite(result.parameter_derivatives)
        || !finite(result.state_parameter_derivatives)
        || !finite(result.packing_parameter_derivatives)
        || !finite(result.molar_volume_bound_parameter_derivatives_m3_per_mol)
        || !finite(result.pressure_parameter_derivatives_pa)
        || !finite(result.chemical_potential_parameter_derivatives_over_rt)
        || !std::isfinite(native.helmholtz_over_rt_reference_amount)
        || !std::isfinite(native.pressure_pa)
        || !std::isfinite(native.packing_fraction)
        || !std::isfinite(native.molar_volume_bounds_m3_per_mol[0])
        || !std::isfinite(native.molar_volume_bounds_m3_per_mol[1])
        || native.molar_volume_bounds_m3_per_mol[0] <= 0.0
        || native.molar_volume_bounds_m3_per_mol[1]
            <= native.molar_volume_bounds_m3_per_mol[0]) {
        throw std::invalid_argument(
            "Provider parameterized-phase tensors or bounds are invalid"
        );
    }
    result.phase.value = native.helmholtz_over_rt_reference_amount;
    result.phase.pressure_pa = native.pressure_pa;
    result.packing.value = native.packing_fraction;
    result.phase.parameter_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider parameterized-phase fingerprint"
    );
    result.topology_fingerprint = decode_provider_char_array(
        native.topology_fingerprint,
        sizeof(native.topology_fingerprint),
        "Provider parameterized-phase topology"
    );
    const std::string basis = decode_provider_char_array(
        native.helmholtz_basis_id,
        sizeof(native.helmholtz_basis_id),
        "Provider parameterized-phase basis"
    );
    if (result.phase.parameter_fingerprint != fingerprint_
        || result.topology_fingerprint != parameterized_topology_fingerprint()
        || basis != EPCSAFT_NATIVE_HELMHOLTZ_BASIS_ID) {
        throw std::invalid_argument(
            "Provider parameterized-phase identity changed"
        );
    }
    result.molar_volume_bounds_m3_per_mol = {
        native.molar_volume_bounds_m3_per_mol[0],
        native.molar_volume_bounds_m3_per_mol[1],
    };
    constexpr double kProjectionRelativeTolerance = 2.0e-11;
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        for (std::size_t species = 0; species < amounts_mol.size(); ++species) {
            const double mixed =
                result.state_parameter_derivatives[
                    species * parameter_count + parameter
                ];
            const double projected =
                result.chemical_potential_parameter_derivatives_over_rt[
                    species * parameter_count + parameter
                ];
            const double scale = std::max({1.0, std::abs(mixed), std::abs(projected)});
            if (std::abs(mixed - projected) > kProjectionRelativeTolerance * scale) {
                throw std::invalid_argument(
                    "Provider active chemical-potential projection is inconsistent"
                );
            }
        }
        const double projected_pressure =
            -8.31446261815324 * temperature_k
            * result.state_parameter_derivatives[
                amounts_mol.size() * parameter_count + parameter
            ];
        const double returned_pressure =
            result.pressure_parameter_derivatives_pa[parameter];
        const double scale = std::max(
            {1.0, std::abs(projected_pressure), std::abs(returned_pressure)}
        );
        if (std::abs(projected_pressure - returned_pressure)
            > kProjectionRelativeTolerance * scale) {
            throw std::invalid_argument(
                "Provider active pressure projection is inconsistent"
            );
        }
    }
    return result;
}

InversePackingGeometryEvaluation ProviderContext::evaluate_inverse_packing_geometry(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double log_packing_fraction,
    const ProviderActiveParameterSet& active_parameters
) const {
    ++evaluation_count_;
    if (transport_.evaluate_inverse_packing_geometry == nullptr
        || amounts_mol.size() != transport_.component_count) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry transport is incomplete"
        );
    }
    require_finite(temperature_k, "inverse-packing temperature");
    require_finite(log_packing_fraction, "log packing fraction");
    if (temperature_k <= 0.0
        || !std::all_of(amounts_mol.begin(), amounts_mol.end(), [](double value) {
            return std::isfinite(value) && value > 0.0;
        })) {
        throw std::invalid_argument(
            "inverse-packing geometry state must be finite and positive"
        );
    }
    const std::string topology = parameterized_topology_fingerprint();
    const std::size_t coordinate_count = amounts_mol.size() + 1
        + active_parameters.identities.size();
    InversePackingGeometryEvaluation result;
    result.gradient.resize(coordinate_count);
    result.hessian.resize(coordinate_count * coordinate_count);
    epcsaft_inverse_packing_geometry_result native{};
    native.struct_size = sizeof(native);
    native.coordinate_count = coordinate_count;
    native.active_parameter_count = active_parameters.identities.size();
    native.gradient_capacity = result.gradient.size();
    native.hessian_capacity = result.hessian.size();
    native.gradient = result.gradient.data();
    native.hessian = result.hessian.data();
    const int status = transport_.evaluate_inverse_packing_geometry(
        transport_.model_context,
        temperature_k,
        amounts_mol.data(),
        amounts_mol.size(),
        log_packing_fraction,
        &native
    );
    if (status != native.status) {
        throw std::runtime_error(
            "Provider inverse-packing geometry returned inconsistent status values"
        );
    }
    if (status != EPCSAFT_NATIVE_STATUS_OK) {
        throw std::domain_error(
            "Provider inverse-packing geometry failed: "
            + decode_provider_char_array(
                native.error,
                sizeof(native.error),
                "Provider inverse-packing geometry error"
            )
        );
    }
    if (native.struct_size != sizeof(native)
        || native.coordinate_count != coordinate_count
        || native.active_parameter_count != active_parameters.identities.size()
        || native.gradient_capacity != result.gradient.size()
        || native.hessian_capacity != result.hessian.size()
        || native.gradient != result.gradient.data()
        || native.hessian != result.hessian.data()) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry result buffer contract changed"
        );
    }
    const std::string returned_parameter_fingerprint = decode_provider_char_array(
        native.parameter_fingerprint,
        sizeof(native.parameter_fingerprint),
        "Provider inverse-packing geometry fingerprint"
    );
    const std::string returned_topology = decode_provider_char_array(
        native.topology_fingerprint,
        sizeof(native.topology_fingerprint),
        "Provider inverse-packing geometry topology"
    );
    if (returned_parameter_fingerprint != fingerprint_ || returned_topology != topology) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry identity changed"
        );
    }
    if (!std::isfinite(native.volume_m3) || native.volume_m3 <= 0.0
        || !std::all_of(result.gradient.begin(), result.gradient.end(),
            [](double value) { return std::isfinite(value); })
        || !std::all_of(result.hessian.begin(), result.hessian.end(),
            [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry tensors are invalid"
        );
    }
    result.volume_m3 = native.volume_m3;
    result.topology_fingerprint = returned_topology;
    return result;
}

const std::string& ProviderContext::fingerprint() const {
    return fingerprint_;
}

const epcsaft_native_transport& ProviderContext::transport() const noexcept {
    return transport_;
}

int ProviderContext::evaluation_count() const noexcept {
    return evaluation_count_;
}

}  // namespace epcsaft_engine::equilibrium
