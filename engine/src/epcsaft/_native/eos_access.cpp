#include "eos.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace epcsaft {
namespace {

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

template <typename Evaluator>
auto eos_evaluation(Evaluator&& evaluator) {
    try {
        return evaluator();
    } catch (const std::invalid_argument& error) {
        throw EosDomainError(error.what());
    } catch (const std::logic_error& error) {
        throw EosDomainError(error.what());
    }
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

bool same_active_coordinates(
    const std::vector<NativeActiveCoordinate>& left,
    const std::vector<NativeActiveCoordinate>& right
) {
    return left.size() == right.size()
        && std::equal(
            left.begin(), left.end(), right.begin(),
            [](const auto& first, const auto& second) {
                return first.identity == second.identity
                    && first.value == second.value;
            }
        );
}

}  // namespace

void NativeModel::validate_temperature(double temperature_k) const {
    require_finite(temperature_k, "temperature");
    const auto& metadata = eos_metadata();
    if (temperature_k <= 0.0
        || (std::isfinite(metadata.temperature_min_k)
            && temperature_k < metadata.temperature_min_k)
        || (std::isfinite(metadata.temperature_max_k)
            && temperature_k > metadata.temperature_max_k)) {
        throw EosDomainError(
            "temperature is outside the parameter source domain"
        );
    }
}

void NativeModel::validate_pressure(double pressure_pa) const {
    require_finite(pressure_pa, "pressure");
    const auto& metadata = eos_metadata();
    if (pressure_pa <= 0.0
        || (std::isfinite(metadata.pressure_min_pa)
            && pressure_pa < metadata.pressure_min_pa)
        || (std::isfinite(metadata.pressure_max_pa)
            && pressure_pa > metadata.pressure_max_pa)) {
        throw EosDomainError("pressure is outside the parameter source domain");
    }
}

void NativeModel::validate_composition(
    const std::vector<double>& mole_fractions
) const {
    if (mole_fractions.size() != component_count()
        || !std::all_of(
            mole_fractions.begin(), mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value >= 0.0; }
        )) {
        throw std::invalid_argument("EOS composition does not match the native model");
    }
    double total = 0.0;
    double charge = 0.0;
    double ion_fraction = 0.0;
    const auto& charges = charge_numbers();
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const double fraction = mole_fractions[component];
        total += fraction;
        charge += fraction * static_cast<double>(charges[component]);
        if (charges[component] != 0) {
            ion_fraction += fraction;
        }
    }
    if (std::abs(total - 1.0) > 1.0e-12
        || std::abs(charge) > 1.0e-12
        || (std::isfinite(eos_metadata().ion_mole_fraction_max)
            && ion_fraction > eos_metadata().ion_mole_fraction_max + 1.0e-12)) {
        throw EosDomainError("composition is outside the EOS domain");
    }
}

std::vector<std::size_t> NativeModel::validate_amounts(
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    require_finite(volume_m3, "volume");
    if (amounts_mol.size() != component_count() || volume_m3 <= 0.0
        || !std::all_of(
            amounts_mol.begin(), amounts_mol.end(),
            [](double amount) {
                return std::isfinite(amount) && amount >= 0.0;
            }
        )) {
        throw std::invalid_argument(
            "EOS amounts and volume do not match the native model"
        );
    }
    std::vector<std::size_t> supported;
    supported.reserve(amounts_mol.size());
    for (std::size_t component = 0; component < amounts_mol.size(); ++component) {
        if (amounts_mol[component] > 0.0) {
            supported.push_back(component);
        }
    }
    const double total = std::accumulate(
        amounts_mol.begin(), amounts_mol.end(), 0.0
    );
    if (!(total > 0.0)) {
        throw std::invalid_argument("EOS total amount must be positive");
    }
    std::vector<double> composition(amounts_mol.size());
    std::transform(
        amounts_mol.begin(), amounts_mol.end(), composition.begin(),
        [total](double amount) { return amount / total; }
    );
    validate_composition(composition);
    return supported;
}

std::vector<std::string> active_parameter_identities(
    const std::vector<NativeActiveCoordinate>& coordinates
) {
    std::vector<std::string> result;
    result.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) {
        result.push_back(coordinate.identity);
    }
    return result;
}

const std::string& NativeModel::parameter_fingerprint() const {
    const std::string& fingerprint = eos_metadata().parameter_fingerprint;
    if (fingerprint.empty() || component_ids().size() != component_count()) {
        throw std::invalid_argument("EOS model identity is incomplete");
    }
    return fingerprint;
}

std::vector<std::string> NativeModel::active_parameter_identities() const {
    return epcsaft::active_parameter_identities(active_coordinates());
}

std::string NativeModel::parameterized_topology_fingerprint() const {
    const std::string& fingerprint = eos_metadata().topology_fingerprint;
    if (fingerprint.empty()) {
        throw std::invalid_argument(
            "EOS parameterized-phase topology metadata is incompatible"
        );
    }
    return fingerprint;
}

NativePhaseState NativeModel::evaluate_phase_state(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    validate_temperature(temperature_k);
    const auto supported = validate_amounts(amounts_mol, volume_m3);
    if (supported.empty()) {
        throw std::invalid_argument("EOS phase support is empty");
    }
    NativeMixturePhaseBlock phase;
    if (amounts_mol.size() == 1 && supports_pure_transport()) {
        const auto pure = eos_evaluation([&] {
            return evaluate_pure_phase(temperature_k, amounts_mol[0], volume_m3);
        });
        phase = {
            pure.helmholtz_over_rt_reference_amount,
            {pure.gradient.begin(), pure.gradient.end()},
            {pure.hessian.begin(), pure.hessian.end()},
            pure.pressure_pa,
        };
    } else if (supported.size() != amounts_mol.size()) {
        phase = eos_evaluation([&] {
            return evaluate_supported_phase_derivatives(
                temperature_k,
                amounts_mol,
                volume_m3,
                eos_metadata().ion_mole_fraction_max,
                {},
                supported
            );
        });
    } else if (supports_electrolyte_transport()) {
        phase = eos_evaluation([&] {
            return evaluate_electrolyte_phase(
                temperature_k, amounts_mol, volume_m3,
                eos_metadata().ion_mole_fraction_max
            );
        });
    } else if (supports_mixture_transport()
               || supports_fixed_topology_associating_mixture()) {
        phase = eos_evaluation([&] {
            return evaluate_mixture_phase(temperature_k, amounts_mol, volume_m3);
        });
    } else {
        throw std::invalid_argument("EOS model has no phase-state evaluator");
    }
    const double total_amount = std::accumulate(
        amounts_mol.begin(), amounts_mol.end(), 0.0
    );
    std::vector<double> composition(amounts_mol.size());
    std::transform(
        amounts_mol.begin(), amounts_mol.end(), composition.begin(),
        [total_amount](double amount) { return amount / total_amount; }
    );
    const double packing_fraction = detail::packing_fraction_coefficient(
        temperature_k, composition, parameters_
    ) * total_amount / volume_m3;
    const std::size_t volume_coordinate = phase.gradient.size() - 1;
    const double mechanical_derivative = -detail::GAS_CONSTANT * temperature_k
        * volume_m3 * phase.hessian[
            volume_coordinate * phase.gradient.size() + volume_coordinate
        ];
    const double scale = std::max(1.0, std::abs(phase.pressure_pa));
    const std::string mechanical_class =
        mechanical_derivative < -1.0e-6 * scale ? "strict_stable"
        : mechanical_derivative > 1.0e-6 * scale ? "unstable"
        : "marginal";
    return {
        phase.helmholtz_over_rt_reference_amount,
        std::move(phase.gradient),
        std::move(phase.hessian),
        phase.pressure_pa,
        packing_fraction,
        mechanical_class,
        supported,
        eos_metadata().topology_fingerprint,
        {1, 0, 0},
    };
}

NativeMolarVolumeBoundsParameterBlock NativeModel::evaluate_molar_volume_bounds(
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double packing_fraction_min,
    double packing_fraction_max
) const {
    validate_temperature(temperature_k);
    require_finite(packing_fraction_min, "minimum packing fraction");
    require_finite(packing_fraction_max, "maximum packing fraction");
    validate_composition(mole_fractions);
    if (packing_fraction_min <= 0.0
        || packing_fraction_max <= packing_fraction_min
    ) {
        throw std::invalid_argument("invalid EOS molar-volume request");
    }
    auto bounds = eos_evaluation([&] {
        return molar_volume_bounds_parameters(
            temperature_k,
            mole_fractions,
            packing_fraction_min,
            packing_fraction_max,
            active_coordinates()
        );
    });
    if (!std::isfinite(bounds.bounds_m3_per_mol[0])
        || !std::isfinite(bounds.bounds_m3_per_mol[1])
        || bounds.bounds_m3_per_mol[0] <= 0.0
        || bounds.bounds_m3_per_mol[1] <= bounds.bounds_m3_per_mol[0]
        || !finite(bounds.parameter_derivatives_m3_per_mol)) {
        throw std::invalid_argument(
            "EOS molar-volume bounds must be finite, positive, and ordered"
        );
    }
    bounds.eos_evaluation_count = 1;
    return bounds;
}

NativePackingFractionBlock NativeModel::evaluate_state_packing_fraction(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    validate_temperature(temperature_k);
    const auto supported = validate_amounts(amounts_mol, volume_m3);
    auto packing = eos_evaluation([&] {
        return evaluate_packing_fraction_parameters(
            temperature_k,
            amounts_mol,
            volume_m3,
            {},
            supported
        );
    });
    if (!std::isfinite(packing.packing_fraction)
        || !finite(packing.gradient) || !finite(packing.hessian)) {
        throw std::invalid_argument("EOS packing-fraction tensors must be finite");
    }
    return packing;
}

NativeParameterizedPhaseState NativeModel::evaluate_parameterized_phase_state(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    const std::vector<NativeActiveCoordinate>& active_parameters
) const {
    validate_temperature(temperature_k);
    const auto supported = validate_amounts(amounts_mol, volume_m3);
    if (supported.empty()) {
        throw std::invalid_argument("parameterized phase support is empty");
    }
    if (!same_active_coordinates(active_parameters, active_coordinates())) {
        throw std::invalid_argument("EOS active-parameter identity changed");
    }
    const auto& coordinates = active_coordinates();
    const std::size_t state_count = supported.size() + 1;
    const std::size_t parameter_count = coordinates.size();
    const std::size_t joint_count = state_count + parameter_count;
    if (coordinates.empty()) {
        NativeParameterizedPhaseState result;
        const auto phase = evaluate_phase_state(
            temperature_k, amounts_mol, volume_m3
        );
        result.phase = {
            phase.value,
            phase.gradient,
            phase.hessian,
            phase.pressure_pa,
        };
        result.supported_component_indices = phase.supported_component_indices;
        result.topology_fingerprint = phase.topology_fingerprint;
        result.evaluations = phase.evaluations;
        return result;
    }
    auto phase = eos_evaluation([&] {
        return evaluate_supported_phase_derivatives(
            temperature_k,
            amounts_mol,
            volume_m3,
            eos_metadata().ion_mole_fraction_max,
            coordinates,
            supported
        );
    });
    if (phase.gradient.size() != joint_count
        || phase.hessian.size() != joint_count * joint_count) {
        throw std::runtime_error("parameterized phase tensor dimensions changed");
    }

    NativeParameterizedPhaseState result;
    result.phase.helmholtz_over_rt_reference_amount =
        phase.helmholtz_over_rt_reference_amount;
    result.phase.pressure_pa = phase.pressure_pa;
    result.supported_component_indices = supported;
    result.topology_fingerprint = parameterized_topology_fingerprint();
    result.phase.gradient.resize(state_count);
    result.phase.hessian.resize(state_count * state_count);
    result.parameter_derivatives.resize(parameter_count);
    result.state_parameter_derivatives.resize(state_count * parameter_count);
    result.pressure_parameter_derivatives_pa.resize(parameter_count);
    result.chemical_potential_parameter_derivatives_over_rt.resize(
        supported.size() * parameter_count
    );
    for (std::size_t state = 0; state < state_count; ++state) {
        result.phase.gradient[state] = phase.gradient[state];
        for (std::size_t other = 0; other < state_count; ++other) {
            result.phase.hessian[state * state_count + other] =
                phase.hessian[state * joint_count + other];
        }
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            result.state_parameter_derivatives[state * parameter_count + parameter] =
                phase.hessian[state * joint_count + state_count + parameter];
        }
    }
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        result.parameter_derivatives[parameter] = phase.gradient[state_count + parameter];
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
        || !finite(result.parameter_derivatives)
        || !finite(result.state_parameter_derivatives)
        || !finite(result.pressure_parameter_derivatives_pa)
        || !finite(result.chemical_potential_parameter_derivatives_over_rt)
        || !std::isfinite(result.phase.helmholtz_over_rt_reference_amount)
        || !std::isfinite(result.phase.pressure_pa)) {
        throw std::invalid_argument("EOS parameterized-phase tensors are invalid");
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
                    "EOS active chemical-potential projection is inconsistent"
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
                "EOS active pressure projection is inconsistent"
            );
        }
    }
    result.evaluations.state = 1;
    return result;
}

NativeParameterizedPackingFraction
NativeModel::evaluate_parameterized_packing_fraction(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    const std::vector<NativeActiveCoordinate>& active_parameters
) const {
    validate_temperature(temperature_k);
    const auto supported = validate_amounts(amounts_mol, volume_m3);
    const std::size_t support_count = supported.size();
    if (support_count == 0
        || !same_active_coordinates(active_parameters, active_coordinates())) {
        throw std::invalid_argument("EOS parameterized packing request is invalid");
    }
    if (active_parameters.empty()) {
        NativeParameterizedPackingFraction result;
        result.packing = evaluate_state_packing_fraction(
            temperature_k, amounts_mol, volume_m3
        );
        result.evaluations.packing = 1;
        return result;
    }
    auto packing = eos_evaluation([&] {
        return evaluate_packing_fraction_parameters(
            temperature_k,
            amounts_mol,
            volume_m3,
            active_parameters,
            supported
        );
    });
    const std::size_t state_count = support_count + 1;
    const std::size_t parameter_count = active_parameters.size();
    const std::size_t joint_count = state_count + parameter_count;
    if (packing.gradient.size() != joint_count
        || packing.hessian.size() != joint_count * joint_count) {
        throw std::runtime_error("parameterized packing tensor dimensions changed");
    }
    NativeParameterizedPackingFraction result;
    result.packing.packing_fraction = packing.packing_fraction;
    result.packing.gradient.resize(state_count);
    result.packing.hessian.resize(state_count * state_count);
    result.parameter_derivatives.resize(parameter_count);
    for (std::size_t state = 0; state < state_count; ++state) {
        result.packing.gradient[state] = packing.gradient[state];
        for (std::size_t other = 0; other < state_count; ++other) {
            result.packing.hessian[state * state_count + other] =
                packing.hessian[state * joint_count + other];
        }
    }
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        result.parameter_derivatives[parameter] =
            packing.gradient[state_count + parameter];
    }
    if (!std::isfinite(result.packing.packing_fraction)
        || !finite(result.packing.gradient) || !finite(result.packing.hessian)
        || !finite(result.parameter_derivatives)) {
        throw std::invalid_argument("EOS parameterized packing tensors are invalid");
    }
    result.evaluations.packing = 1;
    return result;
}

NativeTemperaturePhaseDerivatives NativeModel::evaluate_temperature_phase(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3
) const {
    validate_temperature(temperature_k);
    const auto supported = validate_amounts(amounts_mol, volume_m3);
    const std::size_t support_count = supported.size();
    auto phase = eos_evaluation([&] {
        return evaluate_temperature_phase_derivatives(
            temperature_k,
            amounts_mol,
            volume_m3,
            eos_metadata().ion_mole_fraction_max,
            supported
        );
    });
    const std::size_t state_count = support_count + 1;
    const std::size_t width = state_count + 1;
    if (phase.gradient.size() != width || phase.hessian.size() != width * width) {
        throw std::runtime_error("temperature phase tensor dimensions changed");
    }
    NativeTemperaturePhaseDerivatives result;
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
        throw std::invalid_argument("EOS temperature-phase derivatives are invalid");
    }
    result.evaluations.state = 1;
    return result;
}

NativeInversePackingGeometryResult NativeModel::evaluate_state_inverse_packing_geometry(
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double log_packing_fraction,
    const std::vector<NativeActiveCoordinate>& active_parameters
) const {
    validate_temperature(temperature_k);
    require_finite(log_packing_fraction, "log packing fraction");
    if (amounts_mol.size() != component_count()
        || !std::all_of(
            amounts_mol.begin(), amounts_mol.end(),
            [](double value) { return std::isfinite(value) && value >= 0.0; }
        ) || !same_active_coordinates(active_parameters, active_coordinates())) {
        throw std::invalid_argument("EOS inverse-packing request is invalid");
    }
    static_cast<void>(validate_amounts(amounts_mol, 1.0));
    auto geometry = eos_evaluation([&] {
        return evaluate_inverse_packing_geometry(
            temperature_k,
            amounts_mol,
            log_packing_fraction,
            active_coordinates()
        );
    });
    if (!std::isfinite(geometry.volume_m3) || geometry.volume_m3 <= 0.0
        || !finite(geometry.gradient) || !finite(geometry.hessian)) {
        throw std::invalid_argument(
            "EOS inverse-packing geometry tensors are invalid"
        );
    }
    return {
        geometry.volume_m3,
        std::move(geometry.gradient),
        std::move(geometry.hessian),
        {0, 0, 1},
    };
}

}  // namespace epcsaft
