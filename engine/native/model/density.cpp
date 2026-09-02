#include "model/helmholtz_model.hpp"

#include "physics/numerical_policy.hpp"
#include "physics/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

struct PressurePoint final {
    double residual;
    double pressure;
    double slope;
};

PressurePoint at_density(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    double density
) {
    const auto derivatives = model.evaluate(
        temperature, 1.0 / density, composition, DerivativeOrder::PressureSlope
    );
    const double phi_v = derivatives.gradient[1];
    const double calculated_pressure = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (density - phi_v);
    const double z = calculated_pressure / (density * GAS_CONSTANT_J_PER_MOL_K * temperature);
    return {
        z - pressure / (density * GAS_CONSTANT_J_PER_MOL_K * temperature),
        calculated_pressure,
        GAS_CONSTANT_J_PER_MOL_K * temperature
            * (1.0 + derivatives.volume_second_derivative / (density * density)),
    };
}

bool association_failure(const std::runtime_error& error) {
    const std::string message(error.what());
    return message == "association solve did not converge"
        || message == "association solve left physical domain";
}

DensityRoot refine(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    double lower,
    double upper,
    double lower_residual
) {
    const std::size_t maximum = static_cast<std::size_t>(
        policy::DENSITY_MAX_ITERATIONS.value
    );
    double density = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
    PressurePoint point = at_density(model, temperature, pressure, composition, density);
    std::size_t iterations = 0;
    for (; iterations < maximum; ++iterations) {
        if (std::abs(point.residual) <= policy::DENSITY_ROOT_RESIDUAL.value) {
            break;
        }
        double candidate = density;
        const double derivative = (
            point.slope * density - (point.pressure - pressure)
        ) / (density * density * GAS_CONSTANT_J_PER_MOL_K * temperature);
        if (std::isfinite(derivative) && derivative != 0.0) {
            candidate = density - point.residual / derivative;
        }
        if (!(candidate > lower && candidate < upper) || !std::isfinite(candidate)) {
            candidate = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
        }
        const PressurePoint next = at_density(
            model, temperature, pressure, composition, candidate
        );
        if (lower_residual * next.residual <= 0.0) {
            upper = candidate;
        } else {
            lower = candidate;
            lower_residual = next.residual;
        }
        density = candidate;
        point = next;
    }
    if (iterations == maximum) {
        throw std::runtime_error("density root refinement did not converge");
    }
    return {
        density,
        point.residual,
        point.slope,
        iterations + 1,
        0,
        DensityBranch::Unspecified,
    };
}

}  // namespace

DensityRoot solve_density(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& mole_fractions,
    DensityBranch branch,
    std::optional<double> anchor_density
) {
    // Source: density.hpp on main; one logarithmic bracket scan with a
    // safeguarded Newton refinement is retained for the native model.
    if (!std::isfinite(temperature) || temperature <= 0.0) {
        throw std::invalid_argument("invalid temperature");
    }
    if (!std::isfinite(pressure) || pressure <= 0.0) {
        throw std::invalid_argument("invalid pressure");
    }
    model.validate_composition(mole_fractions);
    const double ideal_density = pressure
        / (GAS_CONSTANT_J_PER_MOL_K * temperature);
    const double lower = ideal_density * policy::DENSITY_LOWER_FRACTION.value;
    const double upper = model.molar_density_upper_bound(temperature, mole_fractions)
        * policy::DENSITY_UPPER_FRACTION.value;
    if (!(lower > 0.0) || !(lower < upper)) {
        throw std::runtime_error("density domain is empty");
    }

    if (anchor_density.has_value() && std::isfinite(*anchor_density)
        && *anchor_density > lower && *anchor_density < upper) {
        double density = *anchor_density;
        for (std::size_t iteration = 0;
             iteration < static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
             ++iteration) {
            const auto point = at_density(
                model, temperature, pressure, mole_fractions, density
            );
            if (std::abs(point.residual) <= policy::DENSITY_ROOT_RESIDUAL.value
                && point.slope > 0.0) {
                return {density, point.residual, point.slope, iteration + 1, 1, branch};
            }
            const double derivative = (
                point.slope * density - (point.pressure - pressure)
            ) / (density * density * GAS_CONSTANT_J_PER_MOL_K * temperature);
            if (!std::isfinite(derivative) || derivative == 0.0) break;
            const double candidate = density - point.residual / derivative;
            if (!(candidate > lower && candidate < upper)
                || !std::isfinite(candidate)) break;
            density = candidate;
        }
    }

    std::vector<DensityRoot> roots;
    bool previous_available = false;
    double previous_density = 0.0;
    double previous_residual = 0.0;
    const std::size_t intervals = static_cast<std::size_t>(
        policy::DENSITY_BRACKET_INTERVALS.value
    );
    const double log_lower = std::log(lower);
    const double log_upper = std::log(upper);
    for (std::size_t index = 0; index <= intervals; ++index) {
        const double fraction = static_cast<double>(index) / intervals;
        const double density = std::exp(log_lower + fraction * (log_upper - log_lower));
        PressurePoint point;
        try {
            point = at_density(model, temperature, pressure, mole_fractions, density);
        } catch (const std::runtime_error& error) {
            if (association_failure(error)) {
                previous_available = false;
                continue;
            }
            throw;
        }
        if (previous_available && previous_residual * point.residual <= 0.0) {
            roots.push_back(refine(
                model, temperature, pressure, mole_fractions,
                previous_density, density, previous_residual
            ));
        }
        previous_available = true;
        previous_density = density;
        previous_residual = point.residual;
    }
    std::vector<DensityRoot> stable;
    for (auto& root : roots) {
        if (std::isfinite(root.pressure_density_derivative)
            && root.pressure_density_derivative > 0.0) {
            stable.push_back(root);
        }
    }
    if (stable.empty()) {
        throw std::runtime_error("no stable density root (count 0)");
    }
    if (stable.size() == 1) {
        stable.front().stable_root_count = 1;
        stable.front().branch = branch;
        return stable.front();
    }
    if (branch == DensityBranch::Unspecified) {
        throw std::runtime_error(
            "density root is ambiguous (stable count " + std::to_string(stable.size()) + ")"
        );
    }
    std::sort(stable.begin(), stable.end(), [](const auto& left, const auto& right) {
        return left.molar_density < right.molar_density;
    });
    DensityRoot result = branch == DensityBranch::Vapor ? stable.front() : stable.back();
    result.stable_root_count = stable.size();
    result.branch = branch;
    return result;
}

}  // namespace epcsaft
