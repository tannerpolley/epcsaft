#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

void require_temperature(
    const epcsaft::NativeProviderMetadata& metadata,
    double temperature_k
) {
    require_finite(temperature_k, "temperature");
    if (temperature_k <= 0.0
        || (std::isfinite(metadata.temperature_min_k)
            && temperature_k < metadata.temperature_min_k)
        || (std::isfinite(metadata.temperature_max_k)
            && temperature_k > metadata.temperature_max_k)) {
        throw ProviderDomainError(
            "temperature is outside the parameter source domain"
        );
    }
}

template <typename Evaluator>
auto provider_evaluation(Evaluator&& evaluator) {
    try {
        return evaluator();
    } catch (const std::invalid_argument& error) {
        throw ProviderDomainError(error.what());
    } catch (const std::logic_error& error) {
        throw ProviderDomainError(error.what());
    }
}

void require_amounts(
    const epcsaft::NativeModel& model,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    bool strictly_positive
) {
    require_finite(volume_m3, "volume");
    if (amounts_mol.size() != model.component_count() || volume_m3 <= 0.0
        || !std::all_of(
            amounts_mol.begin(), amounts_mol.end(),
            [strictly_positive](double amount) {
                return std::isfinite(amount)
                    && (strictly_positive ? amount > 0.0 : amount >= 0.0);
            }
        )) {
        throw std::invalid_argument(
            "provider amounts and volume do not match the native model"
        );
    }
}

std::vector<std::size_t> supported_components(
    const std::vector<double>& amounts_mol
) {
    std::vector<std::size_t> result;
    for (std::size_t component = 0; component < amounts_mol.size(); ++component) {
        if (amounts_mol[component] > 0.0) {
            result.push_back(component);
        }
    }
    return result;
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

}  // namespace

ProviderContext::ProviderContext(const epcsaft::NativeModel& model)
    : model_(model), fingerprint_(model.provider_metadata().parameter_fingerprint) {
    if (fingerprint_.empty() || model_.component_ids().size() != model_.component_count()) {
        throw std::invalid_argument("provider model identity is incomplete");
    }
}

std::string ProviderContext::parameterized_topology_fingerprint() const {
    const std::string& fingerprint = metadata().topology_fingerprint;
    if (fingerprint.empty()) {
        throw std::invalid_argument(
            "Provider parameterized-phase topology metadata is incompatible"
        );
    }
    return fingerprint;
}

PhaseEvaluation ProviderContext::evaluate(
    double temperature_k,
    double amount_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    require_temperature(metadata(), temperature_k);
    require_finite(amount_mol, "amount");
    require_finite(volume_m3, "volume");
    if (!model_.supports_pure_transport() || amount_mol <= 0.0 || volume_m3 <= 0.0) {
        throw std::invalid_argument(
            "provider pure phase requires positive amount and volume"
        );
    }
    auto phase = provider_evaluation([&] {
        return model_.evaluate_pure_phase(temperature_k, amount_mol, volume_m3);
    });
    return {amount_mol, volume_m3, std::move(phase), fingerprint_};
}

MixturePhaseEvaluation ProviderContext::evaluate_mixture(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    require_temperature(metadata(), temperature_k);
    require_amounts(model_, amounts_mol, volume_m3, true);
    if (!supports_mixture_phase()) {
        throw std::invalid_argument("provider model has no mixture phase evaluator");
    }
    auto phase = provider_evaluation([&] {
        return model_.evaluate_mixture_phase(temperature_k, amounts_mol, volume_m3);
    });
    return {
        phase.helmholtz_over_rt_reference_amount,
        std::move(phase.gradient),
        std::move(phase.hessian),
        phase.pressure_pa,
        fingerprint_,
    };
}

MixturePhaseEvaluation ProviderContext::evaluate_electrolyte(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    require_temperature(metadata(), temperature_k);
    require_amounts(model_, amounts_mol, volume_m3, true);
    if (!supports_electrolyte_phase()) {
        throw std::invalid_argument("provider model has no electrolyte phase evaluator");
    }
    auto phase = provider_evaluation([&] {
        return model_.evaluate_electrolyte_phase(
            temperature_k,
            amounts_mol,
            volume_m3,
            metadata().ion_mole_fraction_max
        );
    });
    return {
        phase.helmholtz_over_rt_reference_amount,
        std::move(phase.gradient),
        std::move(phase.hessian),
        phase.pressure_pa,
        fingerprint_,
    };
}

std::array<double, 2> ProviderContext::evaluate_molar_volume_bounds(
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double packing_fraction_min,
    double packing_fraction_max
) const {
    ++evaluation_count_;
    require_temperature(metadata(), temperature_k);
    require_finite(packing_fraction_min, "minimum packing fraction");
    require_finite(packing_fraction_max, "maximum packing fraction");
    if (mole_fractions.size() != model_.component_count()
        || packing_fraction_min <= 0.0
        || packing_fraction_max <= packing_fraction_min
        || !std::all_of(
            mole_fractions.begin(), mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value >= 0.0; }
        )) {
        throw std::invalid_argument("invalid provider molar-volume request");
    }
    const auto bounds = provider_evaluation([&] {
        return model_.molar_volume_bounds(
            temperature_k,
            mole_fractions,
            packing_fraction_min,
            packing_fraction_max,
            metadata().ion_mole_fraction_max
        );
    });
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
    require_temperature(metadata(), temperature_k);
    require_amounts(model_, amounts_mol, volume_m3, false);
    auto packing = provider_evaluation([&] {
        return model_.evaluate_packing_fraction_parameters(
            temperature_k,
            amounts_mol,
            volume_m3,
            metadata().ion_mole_fraction_max,
            {}
        );
    });
    if (!std::isfinite(packing.packing_fraction)
        || !finite(packing.gradient) || !finite(packing.hessian)) {
        throw std::invalid_argument("provider packing-fraction tensors must be finite");
    }
    return {
        packing.packing_fraction,
        std::move(packing.gradient),
        std::move(packing.hessian),
    };
}

ProviderActiveParameterSet ProviderContext::active_parameters() const {
    ProviderActiveParameterSet result;
    result.identities.reserve(model_.active_coordinates().size());
    for (const auto& coordinate : model_.active_coordinates()) {
        result.identities.push_back(coordinate.identity);
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
    require_temperature(metadata(), temperature_k);
    require_amounts(model_, amounts_mol, volume_m3, false);
    const auto supported = supported_components(amounts_mol);
    if (supported.empty()) {
        throw std::invalid_argument("parameterized phase support is empty");
    }
    if (active_parameters.identities != this->active_parameters().identities) {
        throw std::invalid_argument("Provider active-parameter identity changed");
    }
    const auto& coordinates = model_.active_coordinates();
    const std::size_t state_count = supported.size() + 1;
    const std::size_t parameter_count = coordinates.size();
    const std::size_t joint_count = state_count + parameter_count;
    const double total_amount = std::accumulate(
        amounts_mol.begin(), amounts_mol.end(), 0.0
    );
    if (!(total_amount > 0.0)) {
        throw std::invalid_argument("parameterized phase total amount must be positive");
    }
    std::vector<double> composition(amounts_mol.size());
    std::transform(
        amounts_mol.begin(), amounts_mol.end(), composition.begin(),
        [total_amount](double amount) { return amount / total_amount; }
    );
    auto phase = provider_evaluation([&] {
        return model_.evaluate_parameterized_phase(
            temperature_k,
            amounts_mol,
            volume_m3,
            metadata().ion_mole_fraction_max,
            coordinates
        );
    });
    auto packing = provider_evaluation([&] {
        return model_.evaluate_packing_fraction_parameters(
            temperature_k,
            amounts_mol,
            volume_m3,
            metadata().ion_mole_fraction_max,
            coordinates
        );
    });
    auto bounds = provider_evaluation([&] {
        return model_.molar_volume_bounds_parameters(
            temperature_k,
            composition,
            packing_fraction_min,
            packing_fraction_max,
            metadata().ion_mole_fraction_max,
            coordinates
        );
    });
    if (phase.gradient.size() != joint_count
        || phase.hessian.size() != joint_count * joint_count
        || packing.gradient.size() != joint_count
        || packing.hessian.size() != joint_count * joint_count) {
        throw std::runtime_error("parameterized phase tensor dimensions changed");
    }

    ParameterizedPhaseEvaluation result;
    result.phase.value = phase.helmholtz_over_rt_reference_amount;
    result.phase.pressure_pa = phase.pressure_pa;
    result.phase.parameter_fingerprint = fingerprint_;
    result.phase.gradient.resize(state_count);
    result.phase.hessian.resize(state_count * state_count);
    result.packing.value = packing.packing_fraction;
    result.packing.gradient.resize(state_count);
    result.packing.hessian.resize(state_count * state_count);
    result.parameter_derivatives.resize(parameter_count);
    result.state_parameter_derivatives.resize(state_count * parameter_count);
    result.packing_parameter_derivatives.resize(parameter_count);
    result.pressure_parameter_derivatives_pa.resize(parameter_count);
    result.chemical_potential_parameter_derivatives_over_rt.resize(
        supported.size() * parameter_count
    );
    result.supported_component_indices = supported;
    result.topology_fingerprint = parameterized_topology_fingerprint();
    result.molar_volume_bounds_m3_per_mol = bounds.bounds_m3_per_mol;
    result.molar_volume_bound_parameter_derivatives_m3_per_mol =
        std::move(bounds.parameter_derivatives_m3_per_mol);

    for (std::size_t state = 0; state < state_count; ++state) {
        result.phase.gradient[state] = phase.gradient[state];
        result.packing.gradient[state] = packing.gradient[state];
        for (std::size_t other = 0; other < state_count; ++other) {
            result.phase.hessian[state * state_count + other] =
                phase.hessian[state * joint_count + other];
            result.packing.hessian[state * state_count + other] =
                packing.hessian[state * joint_count + other];
        }
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            result.state_parameter_derivatives[state * parameter_count + parameter] =
                phase.hessian[state * joint_count + state_count + parameter];
        }
    }
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        result.parameter_derivatives[parameter] = phase.gradient[state_count + parameter];
        result.packing_parameter_derivatives[parameter] =
            packing.gradient[state_count + parameter];
        result.pressure_parameter_derivatives_pa[parameter] =
            -epcsaft::detail::GAS_CONSTANT * temperature_k
            * result.state_parameter_derivatives[
                (state_count - 1) * parameter_count + parameter
            ];
        for (std::size_t component = 0; component < supported.size(); ++component) {
            result.chemical_potential_parameter_derivatives_over_rt[
                component * parameter_count + parameter
            ] = result.state_parameter_derivatives[
                component * parameter_count + parameter
            ];
        }
    }
    if (!finite(result.phase.gradient) || !finite(result.phase.hessian)
        || !finite(result.packing.gradient) || !finite(result.packing.hessian)
        || !finite(result.parameter_derivatives)
        || !finite(result.state_parameter_derivatives)
        || !finite(result.packing_parameter_derivatives)
        || !finite(result.molar_volume_bound_parameter_derivatives_m3_per_mol)
        || !finite(result.pressure_parameter_derivatives_pa)
        || !finite(result.chemical_potential_parameter_derivatives_over_rt)
        || !std::isfinite(result.phase.value)
        || !std::isfinite(result.phase.pressure_pa)
        || !std::isfinite(result.packing.value)) {
        throw std::invalid_argument("Provider parameterized-phase tensors are invalid");
    }
    constexpr double projection_relative_tolerance = 2.0e-11;
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        for (std::size_t component = 0; component < supported.size(); ++component) {
            const double mixed = result.state_parameter_derivatives[
                component * parameter_count + parameter
            ];
            const double projected =
                result.chemical_potential_parameter_derivatives_over_rt[
                    component * parameter_count + parameter
                ];
            const double scale = std::max({1.0, std::abs(mixed), std::abs(projected)});
            if (std::abs(mixed - projected) > projection_relative_tolerance * scale) {
                throw std::invalid_argument(
                    "Provider active chemical-potential projection is inconsistent"
                );
            }
        }
        const double projected_pressure =
            -epcsaft::detail::GAS_CONSTANT * temperature_k
            * result.state_parameter_derivatives[
                (state_count - 1) * parameter_count + parameter
            ];
        const double returned_pressure =
            result.pressure_parameter_derivatives_pa[parameter];
        const double scale = std::max(
            {1.0, std::abs(projected_pressure), std::abs(returned_pressure)}
        );
        if (std::abs(projected_pressure - returned_pressure)
            > projection_relative_tolerance * scale) {
            throw std::invalid_argument(
                "Provider active pressure projection is inconsistent"
            );
        }
    }
    return result;
}

TemperaturePhaseEvaluation ProviderContext::evaluate_temperature_phase(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    ++evaluation_count_;
    require_temperature(metadata(), temperature_k);
    require_amounts(model_, amounts_mol, volume_m3, false);
    const std::size_t support_count = supported_components(amounts_mol).size();
    auto phase = provider_evaluation([&] {
        return model_.evaluate_temperature_phase_derivatives(
            temperature_k,
            amounts_mol,
            volume_m3,
            metadata().ion_mole_fraction_max
        );
    });
    const std::size_t state_count = support_count + 1;
    const std::size_t width = state_count + 1;
    if (phase.gradient.size() != width || phase.hessian.size() != width * width) {
        throw std::runtime_error("temperature phase tensor dimensions changed");
    }
    TemperaturePhaseEvaluation result;
    result.chemical_potential_temperature_derivatives_over_rt_per_k.resize(
        support_count
    );
    for (std::size_t component = 0; component < support_count; ++component) {
        result.chemical_potential_temperature_derivatives_over_rt_per_k[component] =
            phase.hessian[component * width + state_count];
    }
    result.pressure_temperature_derivative_pa_per_k =
        -epcsaft::detail::GAS_CONSTANT * phase.gradient[support_count]
        - epcsaft::detail::GAS_CONSTANT * temperature_k
            * phase.hessian[support_count * width + state_count];
    if (!finite(result.chemical_potential_temperature_derivatives_over_rt_per_k)
        || !std::isfinite(result.pressure_temperature_derivative_pa_per_k)) {
        throw std::invalid_argument("Provider temperature-phase derivatives are invalid");
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
    require_temperature(metadata(), temperature_k);
    require_finite(log_packing_fraction, "log packing fraction");
    if (amounts_mol.size() != model_.component_count()
        || !std::all_of(
            amounts_mol.begin(), amounts_mol.end(),
            [](double value) { return std::isfinite(value) && value >= 0.0; }
        ) || active_parameters.identities != this->active_parameters().identities) {
        throw std::invalid_argument("Provider inverse-packing request is invalid");
    }
    auto geometry = provider_evaluation([&] {
        return model_.evaluate_inverse_packing_geometry(
            temperature_k,
            amounts_mol,
            log_packing_fraction,
            model_.active_coordinates()
        );
    });
    if (!std::isfinite(geometry.volume_m3) || geometry.volume_m3 <= 0.0
        || !finite(geometry.gradient) || !finite(geometry.hessian)) {
        throw std::invalid_argument(
            "Provider inverse-packing geometry tensors are invalid"
        );
    }
    return {
        geometry.volume_m3,
        std::move(geometry.gradient),
        std::move(geometry.hessian),
        supported_components(amounts_mol),
        parameterized_topology_fingerprint(),
    };
}

const std::string& ProviderContext::fingerprint() const {
    return fingerprint_;
}

const epcsaft::NativeModel& ProviderContext::model() const noexcept {
    return model_;
}

const epcsaft::NativeProviderMetadata& ProviderContext::metadata() const {
    return model_.provider_metadata();
}

std::size_t ProviderContext::component_count() const noexcept {
    return model_.component_count();
}

const std::vector<std::string>& ProviderContext::component_ids() const noexcept {
    return model_.component_ids();
}

const std::vector<int>& ProviderContext::component_charges() const noexcept {
    return model_.charge_numbers();
}

bool ProviderContext::supports_pure_phase() const noexcept {
    return model_.supports_pure_transport();
}

bool ProviderContext::supports_mixture_phase() const noexcept {
    return model_.supports_mixture_transport()
        || model_.supports_fixed_topology_associating_mixture();
}

bool ProviderContext::supports_electrolyte_phase() const noexcept {
    return model_.supports_electrolyte_transport();
}

int ProviderContext::evaluation_count() const noexcept {
    return evaluation_count_;
}

}  // namespace epcsaft_engine::equilibrium
