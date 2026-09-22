#include "properties/helmholtz_model.hpp"

#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

// 2 - phi: the golden-section interior fraction, a mathematical constant.
constexpr double GOLDEN_SECTION = 0.3819660112501051;

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

DensityRoot publish_root(
    double density, double pressure_residual, double pressure_slope, std::size_t iteration_count,
    std::optional<std::size_t> stable_count, DensityBranch root_branch,
    const core::Convergence& record, bool at_iteration_cap = false
) {
    if (record.stop == core::Stop::nonfinite
        || (at_iteration_cap && record.stop != core::Stop::converged)) {
        throw core::ConvergenceFailure(record);
    }
    return {
        density, pressure_residual, pressure_slope, iteration_count,
        stable_count, root_branch, record
    };
}

// Safeguarded Newton from a known iterate inside [lower, upper]; lower_residual carries
// the sign at the lower end so the bracket keeps shrinking around one crossing.
DensityRoot refine(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, double lower, double upper, double lower_residual,
    double density, PressurePoint point, core::Scope& parent_scope, std::size_t& evaluations
) {
    const std::size_t maximum = static_cast<std::size_t>(
        policy::DENSITY_MAX_ITERATIONS.value
    );
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
    const auto root = publish_root(
        density, point.pressure - pressure, point.slope, iterations + 1, std::nullopt,
        DensityBranch::Unspecified, record, iterations == maximum
    );
    parent_scope.probe(
        "density_root",
        {density, point.residual, point.slope, point.slope > 0.0 ? 1.0 : 0.0},
        static_cast<std::uint32_t>(iterations)
    );
    return root;
}

// One scan cell. A sign change is refined from both ends: safeguarded Newton from the vapor
// end lands on the smallest crossing (the reduced residual is concave there) and from the
// liquid end on the largest (convex there away from the spinodal), so a loop narrower than
// the cell still yields both stable branches. A same-sign cell whose slope changes sign
// hides an extremum that is bisected and refined on both halves. One loop per cell.
void scan_cell(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, double lower, const PressurePoint& low,
    double upper, const PressurePoint& high, core::Scope& scope, std::size_t& evaluations,
    std::vector<DensityRoot>& roots
) {
    if (low.residual * high.residual <= 0.0) {
        scope.probe("density_bracket", {lower, upper, low.residual, high.residual});
        const auto left = refine(model, temperature, pressure, composition, lower, upper,
                                 low.residual, lower, low, scope, evaluations);
        const auto right = refine(model, temperature, pressure, composition, lower, upper,
                                  low.residual, upper, high, scope, evaluations);
        // Two converged refinements name one root when they agree within the sum of
        // their first-order radii tol * |d rho / d residual| (merge-biased by design).
        const double coincidence = policy::DENSITY_ROOT_RESIDUAL.value
            * (left.convergence.amplification.value_or(0.0)
               + right.convergence.amplification.value_or(0.0));
        if (std::abs(right.molar_density - left.molar_density) > coincidence) {
            roots.push_back(left);
            roots.push_back(right);
        } else {
            roots.push_back(right.convergence.residual < left.convergence.residual ? right : left);
        }
        return;
    }
    if (!(low.slope * high.slope < 0.0)) return;
    // Bisect the spinodal's slope sign until the residual there changes sign (a hidden pair) or
    // no crossing is possible: across one spinodal |r| cannot fall by more than max|dr/drho|(b-a).
    double a = lower, b = upper, density = lower;
    PressurePoint pa = low, pb = high, point = low;
    for (std::size_t i = 0; i < static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value)
             && !(low.residual * point.residual < 0.0); ++i) {
        const double reach = (b - a) * std::max(std::abs(residual_derivative(pa, a, pressure, temperature)),
                                                std::abs(residual_derivative(pb, b, pressure, temperature)));
        if (std::min(std::abs(pa.residual), std::abs(pb.residual)) > reach) return;
        density = std::sqrt(a * b);
        point = at_density(model, temperature, pressure, composition, density, evaluations,
                           DerivativeOrder::PressureSlope);
        scope.probe("density_spinodal", {density, point.residual, point.slope});
        (pa.slope * point.slope <= 0.0 ? b : a) = density;
        (pa.slope * point.slope <= 0.0 ? pb : pa) = point;
    }
    if (!(low.residual * point.residual < 0.0)) return;
    scan_cell(model, temperature, pressure, composition, lower, low, density, point, scope,
              evaluations, roots);
    scan_cell(model, temperature, pressure, composition, density, point, upper, high, scope,
              evaluations, roots);
}

// linked: the cell ending at this point has no skipped or exactly-zero predecessor.
struct ScanPoint final { double density; PressurePoint point; bool linked; };

// Deterministic log-density grid, then the cells. A loop with positive slopes at both ends
// of its cell is invisible to the cells; with one mechanical loop P_rho has one interior
// minimum, so a positive grid minimum of the slope is bracketed and searched by golden
// section in log density until a negative slope exposes the loop or the bracket closes.
std::vector<DensityRoot> scan_roots(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, double log_lower, double log_upper,
    std::size_t intervals, core::Scope& scope, std::size_t& evaluations
) {
    std::vector<ScanPoint> points;
    bool linked = false;
    for (std::size_t index = 0; index <= intervals; ++index) {
        const double density = std::exp(log_lower + static_cast<double>(index) / intervals
                                        * (log_upper - log_lower));
        try {
            points.push_back({density, at_density(model, temperature, pressure, composition,
                density, evaluations, DerivativeOrder::PressureSlope), linked});
        } catch (const std::runtime_error& error) {
            if (!association_failure(error)) throw;
            linked = false;
            continue;
        }
        scope.probe("density_scan_point", {density, points.back().point.residual,
                    points.back().point.slope}, static_cast<std::uint32_t>(index));
        linked = points.back().point.residual != 0.0;
    }
    std::size_t lowest = 0;
    for (std::size_t i = 1; i + 1 < points.size(); ++i) {
        if (points[i].linked && points[i + 1].linked
            && (lowest == 0 || points[i].point.slope < points[lowest].point.slope)) lowest = i;
    }
    if (lowest != 0 && points[lowest].point.slope > 0.0) {
        double a = std::log(points[lowest - 1].density), b = std::log(points[lowest + 1].density);
        double x = std::log(points[lowest].density), slope_x = points[lowest].point.slope;
        while (b - a > policy::DENSITY_ANCHOR_LOG_WINDOW.value && slope_x > 0.0) {
            const double u = x - a > b - x ? x - GOLDEN_SECTION * (x - a) : x + GOLDEN_SECTION * (b - x);
            const ScanPoint trial{std::exp(u), at_density(model, temperature, pressure, composition,
                std::exp(u), evaluations, DerivativeOrder::PressureSlope), true};
            scope.probe("density_slope_search", {trial.density, trial.point.residual, trial.point.slope});
            points.insert(std::upper_bound(points.begin(), points.end(), trial,
                [](const ScanPoint& l, const ScanPoint& r) { return l.density < r.density; }), trial);
            if (trial.point.slope < slope_x) {
                (u < x ? b : a) = x;
                x = u;
                slope_x = trial.point.slope;
            } else {
                (u < x ? a : b) = u;
            }
        }
    }
    std::vector<DensityRoot> roots;
    for (std::size_t i = 1; i < points.size(); ++i) {
        if (points[i].linked) {
            scan_cell(model, temperature, pressure, composition, points[i - 1].density,
                      points[i - 1].point, points[i].density, points[i].point, scope, evaluations,
                      roots);
        }
    }
    return roots;
}

DensityRoot anchor(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, DensityBranch branch, double lower, double upper,
    double density, std::size_t& evaluations
) {
    core::Scope anchor_scope("density_anchor");
    const double window = policy::DENSITY_ANCHOR_LOG_WINDOW.value, log_anchor = std::log(density);
    const double local_lower = std::exp(std::max(std::log(lower), log_anchor - window));
    const double local_upper = std::exp(std::min(std::log(upper), log_anchor + window));
    const auto maximum = static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
    core::ResidualHistory history;
    core::Convergence anchor_record;
    std::size_t iteration = 0;
    for (; iteration < maximum; ++iteration) {
        const auto point = at_density(model, temperature, pressure, composition, density,
                                      evaluations, DerivativeOrder::PressureSlope);
        history.push(std::abs(point.residual));
        const double derivative = residual_derivative(point, density, pressure, temperature);
        const double candidate = std::isfinite(derivative) && derivative != 0.0
            ? density - point.residual / derivative
            : density;
        anchor_scope.probe("density_newton_step", {density, point.residual, derivative, candidate},
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
            std::size_t local_stable_roots = 0;
            for (const auto& root : scan_roots(model, temperature, pressure, composition,
                     std::log(local_lower), std::log(local_upper),
                     static_cast<std::size_t>(policy::DENSITY_ANCHOR_SCAN_INTERVALS.value),
                     anchor_scope, evaluations)) {
                if (root.pressure_density_derivative > 0.0) ++local_stable_roots;
            }
            anchor_record.evaluations = evaluations;
            anchor_record.stop = core::classify(anchor_record.residual, *anchor_record.criterion,
                anchor_record.floor, anchor_record.contraction, false, anchor_record.amplification);
            if (local_stable_roots == 1) {
                const auto result = publish_root(density, point.pressure - pressure, point.slope,
                                                 iteration + 1, std::nullopt, branch, anchor_record);
                anchor_scope.finish(anchor_record);
                return result;
            }
            anchor_record.stop = core::Stop::infeasible;
            break;
        }
        if (!std::isfinite(derivative) || derivative == 0.0) {
            // classify reports nonfinite for any nonfinite residual, floor, contraction or
            // amplification; a finite record with a dead derivative is singular.
            anchor_record.stop = core::classify(anchor_record.residual, *anchor_record.criterion,
                anchor_record.floor, anchor_record.contraction, false, anchor_record.amplification)
                == core::Stop::nonfinite ? core::Stop::nonfinite : core::Stop::singular;
            break;
        }
        if (!(candidate > local_lower && candidate < local_upper)
            || !std::isfinite(candidate)) {
            anchor_record.stop = core::Stop::left_domain;
            break;
        }
        density = candidate;
    }
    if (iteration == maximum) {
        anchor_record.stop = core::classify(anchor_record.residual, *anchor_record.criterion,
            anchor_record.floor, anchor_record.contraction, true, anchor_record.amplification);
    }
    anchor_scope.finish(anchor_record);
    throw core::ConvergenceFailure(anchor_record);
}

}  // namespace

DensityRoot solve_density(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& mole_fractions, DensityBranch branch,
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
    std::size_t evaluations = 0;
    DensityRoot result;
    if (anchor_density.has_value()) {
        if (branch == DensityBranch::Unspecified) {
            throw std::invalid_argument("density anchor requires explicit branch");
        }
        if (!std::isfinite(*anchor_density) || *anchor_density <= 0.0) {
            throw std::invalid_argument("invalid density anchor");
        }
        if (!(*anchor_density > lower && *anchor_density < upper)) {
            core::Convergence record{0.0, &policy::DENSITY_ROOT_RESIDUAL, std::nullopt,
                std::nullopt, std::nullopt, 0, evaluations, core::Stop::left_domain};
            throw core::ConvergenceFailure(record);
        }
        result = anchor(model, temperature, pressure, mole_fractions, branch, lower, upper,
                        *anchor_density, evaluations);
    } else {
        const auto roots = scan_roots(model, temperature, pressure, mole_fractions,
            std::log(lower), std::log(upper),
            static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value), scope, evaluations);
        std::vector<DensityRoot> stable;
        std::copy_if(roots.begin(), roots.end(), std::back_inserter(stable), [](const auto& root) {
            return std::isfinite(root.pressure_density_derivative) && root.pressure_density_derivative > 0.0;
        });
        if (stable.empty()) {
            throw std::runtime_error("no stable density root (count 0)");
        }
        std::sort(stable.begin(), stable.end(), [](const auto& left, const auto& right) {
            return left.molar_density < right.molar_density;
        });
        scope.probe("stable_root_count", static_cast<double>(stable.size()));
        if (stable.size() == 1) {
            result = stable.front();
            result.stable_root_count = 1;
            result.branch = DensityBranch::Unspecified;
        } else {
            if (branch == DensityBranch::Unspecified) {
                throw std::runtime_error("density root is ambiguous (stable count "
                                         + std::to_string(stable.size()) + ")");
            }
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
