#include "properties/helmholtz_model.hpp"

#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

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
    double gradient_floor;
};

PressurePoint at_density(const HelmholtzModel& model, double temperature, double pressure,
                         const std::vector<double>& composition, double density,
                         std::size_t& evaluations, DerivativeOrder order) {
    ++evaluations;
    const auto derivatives = model.evaluate(temperature, 1.0 / density, composition, order);
    const double phi_v = derivatives.gradient[1];
    const double calculated_pressure = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (density - phi_v);
    const double z = calculated_pressure / (density * GAS_CONSTANT_J_PER_MOL_K * temperature);
    return {
        z - pressure / (density * GAS_CONSTANT_J_PER_MOL_K * temperature),
        calculated_pressure,
        order == DerivativeOrder::PressureSlope
            ? GAS_CONSTANT_J_PER_MOL_K * temperature
                * (1.0 + derivatives.volume_second_derivative / (density * density))
            : 0.0,
        derivatives.gradient_floor,
    };
}

double residual_derivative(
    const PressurePoint& point, double density, double pressure, double temperature
) {
    return (point.slope * density - (point.pressure - pressure))
        / (density * density * GAS_CONSTANT_J_PER_MOL_K * temperature);
}

bool association_failure(const std::runtime_error& error) {
    const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error);
    // Who threw, not how it stopped. Enumerating stops silently narrowed this when
    // `classify` gained `stalled`: the scan propagated instead of skipping the point.
    return failure != nullptr
        && failure->record.criterion == &policy::ASSOCIATION_RESIDUAL;
}

bool far_side_has_no_crossing(const HelmholtzModel& model, double temperature, double pressure,
                              const std::vector<double>& composition, double density,
                              DensityBranch branch, double log_lower, double log_upper,
                              std::size_t intervals, core::Scope& anchor,
                              std::size_t& evaluations) {
    for (std::size_t index = 0; index <= intervals; ++index) {
        const double grid_density = std::exp(log_lower + static_cast<double>(index)
            / intervals * (log_upper - log_lower));
        if ((branch == DensityBranch::Liquid && grid_density <= density)
            || (branch == DensityBranch::Vapor && grid_density >= density)) {
            continue;
        }
        try {
            const auto point = at_density(model, temperature, pressure, composition, grid_density,
                                          evaluations, DerivativeOrder::Gradient);
            anchor.probe("far_side_residual", {grid_density, point.residual}, static_cast<std::uint32_t>(index));
            if ((branch == DensityBranch::Liquid && !(point.residual > 0.0))
                || (branch == DensityBranch::Vapor && !(point.residual < 0.0))) return false;
        } catch (const std::runtime_error& error) {
            if (association_failure(error)) return false;
            throw;
        }
    }
    return true;
}

DensityRoot refine(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    double lower,
    double upper,
    double lower_residual,
    core::Scope& parent_scope,
    std::size_t& evaluations
) {
    const std::size_t maximum = static_cast<std::size_t>(
        policy::DENSITY_MAX_ITERATIONS.value
    );
    double density = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
    PressurePoint point = at_density(model, temperature, pressure, composition, density,
                                     evaluations, DerivativeOrder::PressureSlope);
    core::ResidualHistory history;
    history.push(std::abs(point.residual));
    std::size_t iterations = 0;
    for (; iterations < maximum; ++iterations) {
        if (std::abs(point.residual) <= policy::DENSITY_ROOT_RESIDUAL.value) break;
        const double derivative = residual_derivative(point, density, pressure, temperature);
        double candidate = density;
        if (std::isfinite(derivative) && derivative != 0.0) {
            candidate = density - point.residual / derivative;
        }
        if (!(candidate > lower && candidate < upper) || !std::isfinite(candidate)) {
            candidate = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
        }
        parent_scope.probe(
            "density_newton_step", {density, point.residual, derivative, candidate},
            static_cast<std::uint32_t>(iterations)
        );
        const PressurePoint next = at_density(model, temperature, pressure, composition, candidate,
                                              evaluations, DerivativeOrder::PressureSlope);
        history.push(std::abs(next.residual));
        if (lower_residual * next.residual <= 0.0) {
            upper = candidate;
        } else {
            lower = candidate;
            lower_residual = next.residual;
        }
        density = candidate;
        point = next;
    }
    const double derivative = residual_derivative(
        point, density, pressure, temperature
    );
    core::Convergence record{
        std::abs(point.residual), &policy::DENSITY_ROOT_RESIDUAL,
        point.gradient_floor / density,
        std::isfinite(derivative) && derivative != 0.0
            ? std::optional<double>(1.0 / std::abs(derivative)) : std::nullopt,
        history.contraction(), iterations + 1, evaluations, core::Stop::stalled
    };
    record.stop = core::classify(record.residual, *record.criterion, record.floor,
                                 record.contraction, iterations == maximum,
                                 record.amplification);
    const auto root = DensityRoot::publish(
        density, point.residual, point.slope, iterations + 1, std::nullopt,
        DensityBranch::Unspecified, record, iterations == maximum
    );
    parent_scope.probe(
        "density_root",
        {density, point.residual, point.slope, point.slope > 0.0 ? 1.0 : 0.0},
        static_cast<std::uint32_t>(iterations)
    );
    return root;
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
    core::Scope scope("solve_density");
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

    const std::size_t intervals = static_cast<std::size_t>(
        policy::DENSITY_BRACKET_INTERVALS.value
    );
    const double log_lower = std::log(lower);
    const double log_upper = std::log(upper);
    std::size_t evaluations = 0;
    DensityRoot result;
    bool have_result = false;
    if (branch != DensityBranch::Unspecified && anchor_density.has_value()
        && std::isfinite(*anchor_density) && *anchor_density > lower && *anchor_density < upper) {
        core::Scope anchor("density_anchor");
        double density = *anchor_density;
        core::ResidualHistory history;
        core::Convergence anchor_record;
        std::size_t iteration = 0;
        for (; iteration < static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
             ++iteration) {
            const auto point = at_density(model, temperature, pressure, mole_fractions, density,
                                          evaluations, DerivativeOrder::PressureSlope);
            history.push(std::abs(point.residual));
            const double derivative = residual_derivative(point, density, pressure, temperature);
            const double candidate = std::isfinite(derivative) && derivative != 0.0
                ? density - point.residual / derivative
                : density;
            anchor.probe("density_newton_step", {density, point.residual, derivative, candidate},
                         static_cast<std::uint32_t>(iteration));
            anchor_record = {
                std::abs(point.residual), &policy::DENSITY_ROOT_RESIDUAL,
                point.gradient_floor / density,
                std::isfinite(derivative) && derivative != 0.0
                    ? std::optional<double>(1.0 / std::abs(derivative)) : std::nullopt,
                history.contraction(), iteration + 1, evaluations, core::Stop::stalled
            };
            if (anchor_record.residual <= anchor_record.criterion->value
                && point.slope > 0.0) {
                const bool far_side_clear = far_side_has_no_crossing(
                    model, temperature, pressure, mole_fractions, density, branch,
                    log_lower, log_upper, intervals, anchor, evaluations
                );
                anchor_record.evaluations = evaluations;
                if (far_side_clear) {
                    anchor_record.stop = core::classify(
                        anchor_record.residual, *anchor_record.criterion,
                        anchor_record.floor, anchor_record.contraction, false,
                        anchor_record.amplification
                    );
                    try {
                        result = DensityRoot::publish(
                            density, point.residual, point.slope, iteration + 1,
                            std::nullopt, branch, anchor_record
                        );
                    } catch (const core::ConvergenceFailure&) {
                        anchor.finish(anchor_record);
                        throw;
                    }
                    anchor.finish(anchor_record);
                    have_result = true;
                    break;
                }
                anchor_record.stop = core::classify(
                    anchor_record.residual, *anchor_record.criterion,
                    anchor_record.floor, anchor_record.contraction, false,
                    anchor_record.amplification
                );
                anchor.finish(anchor_record);
                break;
            }
            if (!std::isfinite(derivative) || derivative == 0.0) {
                anchor_record.stop = !std::isfinite(anchor_record.residual)
                    || (anchor_record.floor.has_value()
                        && !std::isfinite(*anchor_record.floor))
                    || (anchor_record.contraction.has_value()
                        && !std::isfinite(*anchor_record.contraction))
                    || (anchor_record.amplification.has_value()
                        && !std::isfinite(*anchor_record.amplification))
                    ? core::Stop::nonfinite
                    : core::Stop::singular;
                anchor.finish(anchor_record);
                break;
            }
            if (!(candidate > lower && candidate < upper)
                || !std::isfinite(candidate)) {
                anchor_record.stop = core::Stop::left_domain;
                anchor.finish(anchor_record);
                break;
            }
            density = candidate;
        }
        if (iteration == static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value)) {
            anchor_record.stop = core::classify(anchor_record.residual, *anchor_record.criterion,
                                                anchor_record.floor, anchor_record.contraction, true,
                                                anchor_record.amplification);
            anchor.finish(anchor_record);
        }
    }
    if (!have_result) {
        std::vector<DensityRoot> roots;
        bool previous_available = false;
        double previous_density = 0.0;
        double previous_residual = 0.0;
        for (std::size_t index = 0; index <= intervals; ++index) {
            const double density = std::exp(log_lower + static_cast<double>(index) / intervals
                                            * (log_upper - log_lower));
            PressurePoint point;
            try {
                point = at_density(model, temperature, pressure, mole_fractions, density,
                                   evaluations, DerivativeOrder::Gradient);
            } catch (const std::runtime_error& error) {
                if (association_failure(error)) {
                    previous_available = false;
                    continue;
                }
                throw;
            }
            scope.probe(
                "density_scan_point",
                {density, point.residual},
                static_cast<std::uint32_t>(index)
            );
            if (previous_available && previous_residual * point.residual <= 0.0) {
                scope.probe(
                    "density_bracket",
                    {previous_density, density, previous_residual, point.residual},
                    static_cast<std::uint32_t>(index)
                );
                roots.push_back(refine(
                    model, temperature, pressure, mole_fractions,
                    previous_density, density, previous_residual, scope, evaluations
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
        scope.probe("stable_root_count", static_cast<double>(stable.size()));
        if (stable.size() == 1) {
            result = stable.front();
            result.stable_root_count = 1;
            result.branch = branch;
        } else {
            if (branch == DensityBranch::Unspecified) {
                throw std::runtime_error(
                    "density root is ambiguous (stable count "
                    + std::to_string(stable.size()) + ")"
                );
            }
            std::sort(stable.begin(), stable.end(), [](const auto& left, const auto& right) {
                return left.molar_density < right.molar_density;
            });
            result = branch == DensityBranch::Vapor ? stable.front() : stable.back();
            result.stable_root_count = stable.size();
            result.branch = branch;
        }
    }
    scope.probe("density_root", result.molar_density);
    scope.finish(result.convergence);
    return result;
}

}  // namespace epcsaft
