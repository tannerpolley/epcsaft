#include "held.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
}

void require_finite(const std::array<double, 2>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::overflow_error(std::string(name) + " must be finite");
    }
}

void require_finite(const std::array<double, 4>& values, const char* name) {
    if (!std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::overflow_error(std::string(name) + " must be finite");
    }
}

}  // namespace

double held_tunnel_minimum_finite_exponential_distance() {
    // exp(p*/distance) is representable only below log(max double).
    return kHeldTunnelPoleStrength / std::log(std::numeric_limits<double>::max());
}

HeldStateEvaluation evaluate_held_state(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    double x_component_zero,
    double log_volume
) {
    require_finite(temperature_k, "temperature");
    require_finite(pressure_pa, "pressure");
    require_finite(x_component_zero, "component-zero mole fraction");
    require_finite(log_volume, "log volume");
    if (temperature_k <= 0.0 || pressure_pa <= 0.0) {
        throw std::invalid_argument("temperature and pressure must be positive");
    }
    if (x_component_zero <= 0.0 || x_component_zero >= 1.0) {
        throw std::invalid_argument("component-zero mole fraction must be strictly between zero and one");
    }
    const double volume_m3 = std::exp(log_volume);
    require_finite(volume_m3, "molar volume");
    if (volume_m3 <= 0.0) {
        throw std::invalid_argument("molar volume must be positive");
    }

    const std::array<double, 2> amounts{x_component_zero, 1.0 - x_component_zero};
    epcsaft::NativePhaseState phase = eos.evaluate_phase_state(
        temperature_k,
        {amounts[0], amounts[1]},
        volume_m3
    );
    const double pressure_over_rt = pressure_pa / (kGasConstantJPerMolK * temperature_k);
    const auto& eos_gradient = phase.gradient;
    const auto& eos_hessian = phase.hessian;

    HeldStateEvaluation result;
    result.x_component_zero = x_component_zero;
    result.log_volume = log_volume;
    result.volume_m3 = volume_m3;
    result.amounts_mol = amounts;
    result.g_bar = phase.value + pressure_over_rt * volume_m3;
    result.gradient = {
        eos_gradient[0] - eos_gradient[1],
        volume_m3 * (eos_gradient[2] + pressure_over_rt),
    };
    result.hessian = {
        eos_hessian[0] - eos_hessian[1] - eos_hessian[3]
            + eos_hessian[4],
        volume_m3 * (eos_hessian[2] - eos_hessian[5]),
        volume_m3 * (eos_hessian[6] - eos_hessian[7]),
        volume_m3 * volume_m3 * eos_hessian[8]
            + volume_m3 * (eos_gradient[2] + pressure_over_rt),
    };
    result.pressure_stationarity_relative = (phase.pressure_pa - pressure_pa) / pressure_pa;
    result.eos = std::move(phase);
    require_finite(result.g_bar, "HELD g_bar");
    require_finite(result.gradient, "HELD gradient");
    require_finite(result.hessian, "HELD Hessian");
    require_finite(result.pressure_stationarity_relative, "pressure stationarity");
    return result;
}

HeldTpdEvaluation evaluate_held_tpd(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const HeldStateEvaluation& reference,
    double x_component_zero,
    double log_volume
) {
    HeldStateEvaluation state = evaluate_held_state(
        eos,
        temperature_k,
        pressure_pa,
        x_component_zero,
        log_volume
    );
    HeldTpdEvaluation result;
    result.d_bar = state.g_bar - reference.g_bar
        - reference.gradient[0] * (x_component_zero - reference.x_component_zero);
    result.gradient = {
        state.gradient[0] - reference.gradient[0],
        state.gradient[1],
    };
    result.hessian = state.hessian;
    result.state = std::move(state);
    require_finite(result.d_bar, "HELD tangent-plane distance");
    require_finite(result.gradient, "HELD tangent-plane gradient");
    require_finite(result.hessian, "HELD tangent-plane Hessian");
    return result;
}

HeldTunnelingEvaluation evaluate_held_tunneling(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const HeldStateEvaluation& reference,
    double minimum_x_component_zero,
    double minimum_tpd,
    double x_component_zero,
    double log_volume
) {
    require_finite(minimum_x_component_zero, "tunneling minimum composition");
    require_finite(minimum_tpd, "tunneling minimum TPD");
    const double composition_delta = x_component_zero - minimum_x_component_zero;
    const double composition_distance = std::abs(composition_delta);
    if (composition_distance <= held_tunnel_minimum_finite_exponential_distance()) {
        throw std::invalid_argument(
            "tunneling point lies in the singular composition neighborhood"
        );
    }
    HeldTpdEvaluation tpd = evaluate_held_tpd(
        eos,
        temperature_k,
        pressure_pa,
        reference,
        x_component_zero,
        log_volume
    );
    const double exponent = kHeldTunnelPoleStrength / composition_distance;
    const double exponential = std::exp(exponent);
    const double difference = tpd.d_bar - minimum_tpd;
    const double inverse_distance_cubed = 1.0
        / (composition_distance * composition_distance * composition_distance);
    const double log_gradient_x =
        -kHeldTunnelPoleStrength * composition_delta * inverse_distance_cubed;
    const double log_hessian_xx =
        2.0 * kHeldTunnelPoleStrength * inverse_distance_cubed;

    HeldTunnelingEvaluation result;
    result.objective = difference * exponential;
    result.gradient = {
        exponential * (tpd.gradient[0] + difference * log_gradient_x),
        exponential * tpd.gradient[1],
    };
    result.hessian = {
        exponential
            * (tpd.hessian[0] + 2.0 * tpd.gradient[0] * log_gradient_x
               + difference
                   * (log_hessian_xx + log_gradient_x * log_gradient_x)),
        exponential * (tpd.hessian[1] + log_gradient_x * tpd.gradient[1]),
        exponential * (tpd.hessian[2] + tpd.gradient[1] * log_gradient_x),
        exponential * tpd.hessian[3],
    };
    result.tpd = std::move(tpd);
    require_finite(result.objective, "HELD tunneling objective");
    require_finite(result.gradient, "HELD tunneling gradient");
    require_finite(result.hessian, "HELD tunneling Hessian");
    return result;
}

}  // namespace epcsaft_engine::equilibrium
