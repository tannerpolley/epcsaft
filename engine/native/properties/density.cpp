#include "properties/helmholtz_model.hpp"

#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <iterator>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft {
namespace {

// The diagnostic scan's evaluation or time budget ran out; it ends the scan instead of
// becoming one sample's cause.
struct DiagnosticBudgetExhausted final : std::runtime_error {
    DiagnosticBudgetExhausted() : std::runtime_error("density diagnostic budget exhausted") {}
};

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
    if (!std::all_of(derivatives.gradient.begin(), derivatives.gradient.end(), [](double v) { return std::isfinite(v); })
        || !std::all_of(derivatives.hessian.begin(), derivatives.hessian.end(), [](double v) { return std::isfinite(v); }))
        throw std::runtime_error("nonfinite density derivatives");
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature, calculated = rt * (density - derivatives.gradient[1]);
    const double scale = density * GAS_CONSTANT_J_PER_MOL_K * temperature;  // z - P/(rho R T), as recorded
    return {calculated / scale - pressure / scale, calculated, order == DerivativeOrder::PressureSlope
            ? rt * (1.0 + derivatives.volume_second_derivative / (density * density)) : 0.0, derivatives.gradient_floor};
}

double residual_derivative(
    const PressurePoint& point, double density, double pressure, double temperature
) {
    return (point.slope * density - (point.pressure - pressure))
        / (density * density * GAS_CONSTANT_J_PER_MOL_K * temperature);
}

bool association_failure(const std::exception& error) {
    const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error);
    // Who threw, not how it stopped. Enumerating stops silently narrowed this when
    // `classify` gained `stalled`: the scan propagated instead of skipping the point.
    return failure != nullptr
        && failure->record.criterion == &policy::ASSOCIATION_RESIDUAL;
}

// Safeguarded Newton from the bracket endpoint with smaller absolute pressure error.
// Both endpoint evaluations enclose one crossing in [lower, upper]. The caller
// isolates a mechanically stable interval, or supplies an exact scan-node root.
DensityRoot refine(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, double lower, PressurePoint low,
    double upper, PressurePoint high, core::Scope& parent_scope, std::size_t& evaluations
) {
    const auto maximum = static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value);
    const double tolerance = policy::DENSITY_ROOT_RESIDUAL.value;
    const bool from_lower = std::abs(low.pressure - pressure) < std::abs(high.pressure - pressure);
    double density = from_lower ? lower : upper;
    PressurePoint point = from_lower ? low : high;
    core::ResidualHistory history;
    history.push(std::abs(point.residual));
    std::size_t iterations = 0;
    for (; iterations < maximum; ++iterations) {
        const double derivative = residual_derivative(point, density, pressure, temperature);
        const bool residual_met = std::abs(point.residual) <= tolerance;
        if (point.residual == 0.0 || (residual_met && std::isfinite(derivative) && derivative != 0.0
            && (upper - lower) * std::abs(derivative) <= tolerance)) break;
        // Once the residual passes, the bracket must still witness the tolerance/amplification
        // radius: step a quarter radius past the Newton root toward the far end (bisection fallback).
        double candidate = policy::DENSITY_BISECTION_FRACTION.value * (lower + upper);
        if (std::isfinite(derivative) && derivative != 0.0) {
            const double far = upper - density > density - lower ? upper : lower;
            candidate = density - point.residual / derivative
                + (residual_met ? std::copysign(0.25 * tolerance / std::abs(derivative), far - density) : 0.0);
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
        if (low.residual * next.residual <= 0.0) {
            upper = candidate;
            high = next;
        } else {
            lower = candidate;
            low = next;
        }
        const bool from_lower = std::abs(low.pressure - pressure) < std::abs(high.pressure - pressure);
        density = from_lower ? lower : upper;
        point = from_lower ? low : high;
    }
    const double derivative = residual_derivative(point, density, pressure, temperature);
    parent_scope.probe("density_root_bracket", {lower, upper, derivative, point.residual});
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
    if (iterations == maximum && record.stop == core::Stop::converged)
        throw std::runtime_error("density root bracket unresolved within the search budget");
    if (record.stop == core::Stop::nonfinite
        || (iterations == maximum && record.stop != core::Stop::converged)) {
        throw core::ConvergenceFailure(record);
    }
    const DensityRoot root{density, point.pressure - pressure, point.slope, iterations + 1,
                           std::nullopt, DensityBranch::Unspecified, record};
    parent_scope.probe("density_root", {density, point.residual, point.slope, point.slope > 0.0 ? 1.0 : 0.0},
                       static_cast<std::uint32_t>(iterations));
    return root;
}

// Once the single mechanical loop is exposed, each cell contains at most one fold.
// Isolate its positive-slope interval before refining one strict pressure crossing.
// Exact scan-node roots belong to scan_roots, not to either adjacent cell.
void scan_cell(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, double lower, PressurePoint low,
    double upper, PressurePoint high, core::Scope& scope, std::size_t& evaluations,
    std::vector<DensityRoot>& roots
) {
    if (low.slope <= 0.0 && high.slope <= 0.0) return;
    if ((low.slope > 0.0) != (high.slope > 0.0)) {
        const bool stable_left = low.slope > 0.0;
        double a = lower, b = upper;
        PressurePoint pa = low, pb = high;
        // The stable side's pressure runs monotonically from its cell end toward the fold, so a
        // cell end already past the target leaves that side rootless without locating the fold.
        const double margin = policy::DENSITY_ROOT_RESIDUAL.value * GAS_CONSTANT_J_PER_MOL_K * temperature * upper;
        bool bracketed = stable_left ? low.pressure > pressure + margin : high.pressure < pressure - margin;
        for (std::size_t i = 0; i < static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value)
                 && !bracketed && pa.slope != 0.0 && pb.slope != 0.0; ++i) {
            const double density = a + policy::DENSITY_BISECTION_FRACTION.value * (b - a);
            if (!(density > a && density < b)) break;
            const auto point = at_density(model, temperature, pressure, composition, density,
                                          evaluations, DerivativeOrder::PressureSlope);
            scope.probe("density_spinodal", {density, point.residual, point.slope});
            if ((pa.slope > 0.0) == (point.slope > 0.0)) { a = density; pa = point; }
            else { b = density; pb = point; }
            // P rises across the stable side: a stable point past the target already
            // brackets the unique root; only a rootless side needs the spinodal itself.
            bracketed = point.slope > 0.0 && std::abs(point.residual) > policy::DENSITY_ROOT_RESIDUAL.value
                && (stable_left ? low.residual < 0.0 : high.residual > 0.0) && (point.residual > 0.0) == stable_left;
        }
        if (!bracketed && pa.slope != 0.0 && pb.slope != 0.0 && std::nextafter(a, b) < b)
            throw std::runtime_error("density spinodal unresolved within the search budget");
        const bool use_a = pa.slope == 0.0 || (pb.slope != 0.0 && stable_left);
        const double density = use_a ? a : b;
        const auto point = use_a ? pa : pb;
        if (std::abs(point.residual) <= policy::DENSITY_ROOT_RESIDUAL.value)
            throw std::runtime_error("density branch unresolved at spinodal pressure within residual tolerance");
        if (stable_left) { upper = density; high = point; }
        else { lower = density; low = point; }
    }
    if (low.residual < 0.0 && high.residual > 0.0) {
        scope.probe("density_bracket", {lower, upper, low.residual, high.residual});
        roots.push_back(refine(model, temperature, pressure, composition, lower, low,
            upper, high, scope, evaluations));
    }
}

// linked: the cell ending at this point has no skipped predecessor.
struct ScanPoint final { double density; PressurePoint point; bool linked; };

void search_slope(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, std::vector<ScanPoint>& points,
    core::Scope& scope, std::size_t& evaluations,
    const std::function<PressurePoint(double)>& evaluate = {},
    const std::function<void(double, double, const std::string&)>& unresolved = {}
) {
    std::size_t lowest = 0;
    for (std::size_t i = 1; i + 1 < points.size(); ++i) {
        if (points[i].linked && points[i + 1].linked
            && points[i].point.slope <= std::min(points[i - 1].point.slope, points[i + 1].point.slope)
            && (lowest == 0 || points[i].point.slope < points[lowest].point.slope)) lowest = i;
    }
    if (lowest != 0 && points[lowest].point.slope > 0.0) {
        double a = std::log(points[lowest - 1].density), b = std::log(points[lowest + 1].density);
        double x = std::log(points[lowest].density), slope_x = points[lowest].point.slope;
        double slope_a = points[lowest - 1].point.slope, slope_b = points[lowest + 1].point.slope;
        // A positive lower bound for the interpolating quadratic stops this model-based search.
        // It is not a bound on unsampled EOS curvature, even near criticality in log density.
        try {
            for (std::size_t searched = 0;; ++searched) {
                if (!(b - a > std::numeric_limits<double>::epsilon() * (1.0 + std::abs(x))))
                    throw std::runtime_error("density slope search bracket collapsed before resolution");
                const double curvature = 2.0 * (slope_a / ((a - x) * (a - b)) + slope_x / ((x - a) * (x - b))
                                                + slope_b / ((b - a) * (b - x)));
                if (slope_x - 0.5 * std::max(curvature, 0.0) * (b - a) * (b - a) > 0.0) break;
                if (searched == static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value))
                    throw std::runtime_error("density slope minimum unresolved within the search budget");
                const double u = x - a > b - x ? x - GOLDEN_SECTION * (x - a) : x + GOLDEN_SECTION * (b - x);
                const ScanPoint trial{std::exp(u), evaluate ? evaluate(std::exp(u)) : at_density(model, temperature, pressure, composition,
                    std::exp(u), evaluations, DerivativeOrder::PressureSlope), true};
                scope.probe("density_slope_search", {trial.density, trial.point.residual, trial.point.slope});
                if (trial.point.slope < 0.0) {
                    points.insert(std::upper_bound(points.begin(), points.end(), trial,
                        [](const ScanPoint& l, const ScanPoint& r) { return l.density < r.density; }), trial);
                    break;
                }
                if (trial.point.slope < slope_x) {
                    (u < x ? slope_b : slope_a) = slope_x;
                    (u < x ? b : a) = x;
                    x = u;
                    slope_x = trial.point.slope;
                } else {
                    (u < x ? slope_a : slope_b) = trial.point.slope;
                    (u < x ? a : b) = u;
                }
            }
        } catch (const std::runtime_error& error) {
            if (unresolved) unresolved(std::exp(a), std::exp(b), error.what());
            throw;
        }
    }
}

// Deterministic log-density grid, then the cells. A loop with positive slopes at both ends
// of its cell is invisible to the cells; with one mechanical loop P_rho has one interior
// minimum, so a positive grid minimum of the slope is bracketed and searched by golden
// section in log density until a negative slope exposes the loop or the quadratic slope model
// estimates its minimum positive; this finite search is not a completeness certificate.
std::vector<DensityRoot> scan_roots(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, double log_lower, double log_upper,
    std::size_t intervals, core::Scope& scope, std::size_t& evaluations
) {
    std::vector<ScanPoint> points;
    std::string interruptions;
    bool linked = false;
    if (scope.diagnostics_on()) scope.probe("density_scan_input",
        {temperature, pressure, std::exp(log_lower), std::exp(log_upper), static_cast<double>(intervals)});
    for (std::size_t index = 0; index <= intervals; ++index) {
        const double density = std::exp(log_lower + static_cast<double>(index) / intervals
                                        * (log_upper - log_lower));
        try {
            points.push_back({density, at_density(model, temperature, pressure, composition,
                density, evaluations, DerivativeOrder::PressureSlope), linked});
        } catch (const std::exception& error) {
            const bool association = association_failure(error);
            if (scope.diagnostics_on()) scope.probe(association ? "density_scan_skip" : "density_scan_failure",
                {density, static_cast<double>(index), 0.0}, static_cast<std::uint32_t>(index), error.what());
            if (!association && interruptions.empty()) throw;
            interruptions += "node " + std::to_string(index) + " rho "
                + core::optional_number(density) + ": " + error.what() + "\n";
            if (!association)
                throw std::runtime_error("density scan interrupted by association: " + interruptions);
            linked = false;
            continue;
        }
        scope.probe("density_scan_point", {density, points.back().point.residual,
            points.back().point.slope, linked ? 1.0 : 0.0}, static_cast<std::uint32_t>(index));
        linked = true;
    }
    if (!interruptions.empty())
        throw std::runtime_error("density scan interrupted by association: " + interruptions);
    search_slope(model, temperature, pressure, composition, points, scope, evaluations);
    std::vector<DensityRoot> roots;
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& node = points[i];
        if (node.point.residual == 0.0 && node.point.slope > 0.0)
            roots.push_back(refine(model, temperature, pressure, composition, node.density,
                node.point, node.density, node.point, scope, evaluations));
        if (i != 0 && points[i].linked) {
            scan_cell(model, temperature, pressure, composition, points[i - 1].density,
                      points[i - 1].point, points[i].density, points[i].point, scope, evaluations,
                      roots);
        }
    }
    const auto witness = std::find_if(points.begin(), points.end(), [](const auto& node) {
        return std::isfinite(node.point.slope) && node.point.slope < 0.0;
    });
    if (witness != points.end())
        for (auto& root : roots) root.negative_slope_density = witness->density;
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
            auto roots = scan_roots(model, temperature, pressure, composition,
                std::log(local_lower), std::log(local_upper),
                static_cast<std::size_t>(policy::DENSITY_ANCHOR_SCAN_INTERVALS.value),
                anchor_scope, evaluations);
            roots.erase(std::remove_if(roots.begin(), roots.end(),
                [](const auto& root) { return !(root.pressure_density_derivative > 0.0); }), roots.end());
            anchor_record.evaluations = evaluations;
            if (roots.size() == 1) {
                auto result = roots.front();
                result.branch = branch;
                result.convergence.iterations += iteration + 1;
                result.convergence.evaluations = evaluations;
                result.iterations = result.convergence.iterations;
                anchor_scope.finish(result.convergence);
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

DensityDiagnostics HelmholtzModel::density_diagnostics(
    double temperature, double pressure, const std::vector<double>& composition,
    const std::vector<double>& nodes, double density_tolerance, double log_tolerance,
    std::size_t max_evaluations, double max_seconds, bool adaptive,
    std::optional<double> finite_anchor, double finite_upper, DerivativeOrder evaluation_order
) const {
    core::Scope scope("density_diagnostics");
    validate_composition(composition);
    if (!std::isfinite(temperature) || temperature <= 0.0 || !std::isfinite(pressure) || pressure <= 0.0
        || !std::isfinite(density_tolerance) || density_tolerance <= 0.0
        || !std::isfinite(log_tolerance) || log_tolerance <= 0.0
        || !std::isfinite(max_seconds) || max_seconds <= 0.0 || max_evaluations == 0)
        throw std::invalid_argument("invalid density diagnostic inputs or budget");
    DensityDiagnostics result;
    result.evaluation_order = evaluation_order;
    result.geometric_cap = molar_density_upper_bound(temperature, composition);
    if (nodes.size() < 2) throw std::invalid_argument("density diagnostics require at least two declared nodes");
    for (std::size_t i = 0; i < nodes.size(); ++i)
        if (!std::isfinite(nodes[i]) || nodes[i] <= 0.0 || nodes[i] >= result.geometric_cap
            || (i && nodes[i] <= nodes[i - 1])) throw std::invalid_argument("invalid declared density grid");
    if (finite_anchor && (!std::isfinite(*finite_anchor) || *finite_anchor <= 0.0
        || *finite_anchor >= result.geometric_cap || !std::isfinite(finite_upper) || finite_upper <= *finite_anchor))
        throw std::invalid_argument("invalid finite-domain interval");
    if (evaluation_order != DerivativeOrder::Gradient && evaluation_order != DerivativeOrder::PressureSlope
        && evaluation_order != DerivativeOrder::Hessian)
        throw std::invalid_argument("finite-domain evaluation requires Gradient, PressureSlope or Hessian");
    const auto started = std::chrono::steady_clock::now();
    auto evaluate = [&](double rho, DerivativeOrder order = DerivativeOrder::PressureSlope) {
        if (result.evaluations >= max_evaluations
            || std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() >= max_seconds)
            throw DiagnosticBudgetExhausted();
        if (rho >= result.geometric_cap) throw std::runtime_error("geometric domain exit");
        auto point = at_density(*this, temperature, pressure, composition, rho, result.evaluations, order);
        if (!std::isfinite(point.pressure) || !std::isfinite(point.slope) || !std::isfinite(point.gradient_floor))
            throw std::runtime_error("nonfinite density evaluation");
        return point;
    };
    auto row = [&](double rho, const PressurePoint& point, double lo, double hi, bool spinodal) {
        DensityDiagnosticPoint value;
        value.molar_density = rho;
        value.packing_fraction = rho * packing_fraction_coefficient(temperature, composition);
        value.pressure = point.pressure;
        value.pressure_density_derivative = point.slope;
        value.mechanical_sign = spinodal ? 0 : (point.slope > 0.0) - (point.slope < 0.0);
        value.lower = lo;
        value.upper = hi;
        value.residual = spinodal ? point.slope / (GAS_CONSTANT_J_PER_MOL_K * temperature) : point.pressure - pressure;
        return value;
    };
    auto refine_crossing = [&](double lo, PressurePoint low, double hi, PressurePoint high, bool spinodal) {
        auto value = [&](const PressurePoint& p) { return spinodal ? p.slope : p.pressure - pressure; };
        std::string stop = "bracket_tolerance";
        try {
            while (hi - lo > density_tolerance || std::log(hi / lo) > log_tolerance) {
                if (value(low) == 0.0 || value(high) == 0.0) break;
                double rho = lo + 0.5 * (hi - lo);
                if (!(rho > lo && rho < hi))
                    throw std::runtime_error("density diagnostic bracket unresolved at adjacent densities");
                const auto point = evaluate(rho);
                if ((value(point) > 0.0) == (value(low) > 0.0)) { lo = rho; low = point; }
                else { hi = rho; high = point; }
            }
        } catch (const std::runtime_error& error) {
            auto record = row(lo, low, lo, hi, spinodal);
            record.stop = "unresolved_bracket";
            record.cause = error.what();
            result.unresolved.push_back(record);
            throw;
        }
        const bool use_low = std::abs(value(low)) <= std::abs(value(high));
        if (value(use_low ? low : high) == 0.0) stop = "exact_numerical_zero";
        auto record = row(use_low ? lo : hi, use_low ? low : high, lo, hi, spinodal);
        record.stop = stop;
        return record;
    };
    std::vector<ScanPoint> points;
    bool linked = false;
    try {
        for (double rho : nodes) {
            try {
                const auto point = evaluate(rho);
                auto sample = row(rho, point, rho, rho, false);
                sample.linked = linked;
                result.samples.push_back(sample);
                if (point.pressure == pressure) {
                    sample.stop = "exact_numerical_zero";
                    result.pressure_roots.push_back(sample);
                }
                if (point.slope == 0.0) {
                    auto fold = row(rho, point, rho, rho, true);
                    fold.stop = "exact_numerical_zero";
                    result.spinodals.push_back(fold);
                }
                if (linked && (points.back().point.pressure - pressure) * (point.pressure - pressure) < 0.0)
                    result.declared_brackets.push_back({points.back().density, rho,
                        points.back().point.pressure - pressure, point.pressure - pressure});
                result.declared_crossings = result.declared_brackets.size();
                points.push_back({rho, point, linked});
                linked = true;
            } catch (const DiagnosticBudgetExhausted&) { throw;
            } catch (const std::runtime_error& error) {
                DensityDiagnosticPoint sample;
                sample.molar_density = sample.lower = sample.upper = rho;
                sample.cause = error.what();
                result.samples.push_back(sample);
                linked = false;
            }
        }
        auto crossings = [&](const std::vector<ScanPoint>& cells) {
            for (std::size_t i = 0; i < cells.size(); ++i) {
                const auto& node = cells[i];
                if (node.point.pressure == pressure && std::none_of(result.pressure_roots.begin(),
                    result.pressure_roots.end(), [&](const auto& root) { return root.molar_density == node.density; })) {
                    auto record = row(node.density, node.point, node.density, node.density, false);
                    record.stop = "exact_numerical_zero";
                    result.pressure_roots.push_back(record);
                }
                if (i && node.linked && (cells[i - 1].point.pressure - pressure) * (node.point.pressure - pressure) < 0.0) {
                    const auto& low = cells[i - 1];
                    // A retained sign bracket or exact numerical zero owns its crossing;
                    // enclosure in this connected partition avoids re-refinement, not distance merging.
                    const bool owned = std::any_of(result.pressure_roots.begin(), result.pressure_roots.end(), [&](const auto& root) {
                        return (root.lower >= low.density && root.upper <= node.density)
                            || (root.stop == "exact_numerical_zero" && root.molar_density > low.density && root.molar_density < node.density);
                    });
                    if (!owned) result.pressure_roots.push_back(refine_crossing(low.density, low.point, node.density, node.point, false));
                }
            }
        };
        crossings(points);
        if (adaptive) {
            search_slope(*this, temperature, pressure, composition, points, scope, result.evaluations,
                         [&](double rho) {
                auto point = evaluate(rho);
                auto sample = row(rho, point, rho, rho, false);
                sample.declared = false;
                sample.linked = true;
                result.samples.push_back(sample);
                return point;
            }, [&](double lo, double hi, const std::string& cause) {
                DensityDiagnosticPoint bracket;
                bracket.lower = lo; bracket.upper = hi;
                bracket.stop = "unresolved_slope_search"; bracket.cause = cause;
                result.unresolved.push_back(bracket);
            });
        }
        std::vector<ScanPoint> partitions;
        for (std::size_t i = 0; i < points.size(); ++i) {
            const auto& node = points[i];
            if (i && node.linked && points[i - 1].point.slope * node.point.slope < 0.0) {
                const auto& low = points[i - 1];
                const auto fold = refine_crossing(low.density, low.point, node.density, node.point, true);
                result.spinodals.push_back(fold);
                if (adaptive) partitions.push_back({fold.molar_density, evaluate(fold.molar_density), true});
            }
            partitions.push_back(node);
        }
        if (adaptive) crossings(partitions);
        result.stop = result.pressure_roots.empty() ? "no_sampled_crossing" : "resolved_crossings";
        for (const auto& fold : result.spinodals)
            if (std::abs(fold.pressure - pressure) <= policy::DENSITY_ROOT_RESIDUAL.value
                * fold.molar_density * GAS_CONSTANT_J_PER_MOL_K * temperature)
                result.stop = "unresolved_tangent_pressure_root";
        if (finite_anchor) {
            double lo = *finite_anchor, hi = lo;
            evaluate(lo, evaluation_order);
            result.finite_domain = DensityDiagnosticPoint{};
            result.finite_domain->lower = result.finite_domain->upper = lo;
            const double cap = std::min(finite_upper, result.geometric_cap);
            bool unavailable = false;
            while (hi < cap) {
                hi = std::min(2.0 * lo, cap);
                result.finite_domain->upper = hi;
                try { evaluate(hi, evaluation_order); lo = hi; }
                catch (const DiagnosticBudgetExhausted&) { throw; }
                catch (const std::runtime_error& error) {
                    result.finite_domain->cause = error.what(); unavailable = true; break;
                }
                result.finite_domain->lower = lo;
            }
            while (unavailable && std::nextafter(lo, hi) < hi) {
                const double rho = lo + 0.5 * (hi - lo);
                try { evaluate(rho, evaluation_order); lo = rho; }
                catch (const DiagnosticBudgetExhausted&) { throw; }
                catch (const std::runtime_error& error) {
                    hi = rho; result.finite_domain->cause = error.what();
                }
                result.finite_domain->lower = lo; result.finite_domain->upper = hi;
            }
            result.finite_domain->molar_density = lo;
            result.finite_domain->stop = unavailable ? "finite_to_unavailable" : "declared_interval_finite";
        }
    } catch (const std::runtime_error& error) {
        result.stop = error.what();
        if (result.finite_domain) result.finite_domain->stop = "unresolved_bracket";
    }
    for (auto* roots : {&result.pressure_roots, &result.spinodals})
        std::sort(roots->begin(), roots->end(), [](const auto& a, const auto& b) { return a.molar_density < b.molar_density; });
    return result;
}

bool properties::regular_density(double density, double slope, const std::optional<core::Convergence>& record) {
    if (!(density > 0.0) || !std::isfinite(density) || !(slope > 0.0) || !std::isfinite(slope)
        || !record || !record->floor || !record->amplification) return false;
    const auto& c = *record;
    return std::isfinite(c.residual) && c.residual <= policy::DENSITY_ROOT_RESIDUAL.value
        && *c.floor >= 0.0 && std::isfinite(*c.floor)
        && *c.amplification > 0.0 && std::isfinite(*c.amplification)
        && (policy::DENSITY_ROOT_RESIDUAL.value + *c.floor) * *c.amplification / density
            <= policy::DENSITY_ANCHOR_REGULARITY.value;
}


DensityRoot solve_density(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& mole_fractions, DensityBranch branch,
    std::optional<double> anchor_density, DensityRoot* sibling_root
) {
    core::Scope scope("solve_density");
    if (anchor_density.has_value() && sibling_root != nullptr)
        throw std::invalid_argument("density sibling requires an unanchored root search");
    if (!std::isfinite(temperature) || temperature <= 0.0) throw std::invalid_argument("invalid temperature");
    if (!std::isfinite(pressure) || pressure <= 0.0) throw std::invalid_argument("invalid pressure");
    model.validate_composition(mole_fractions);
    const double ideal_density = pressure
        / (GAS_CONSTANT_J_PER_MOL_K * temperature);
    const double lower = ideal_density * policy::DENSITY_LOWER_FRACTION.value;
    const double upper = model.molar_density_upper_bound(temperature, mole_fractions)
        * policy::DENSITY_UPPER_FRACTION.value;
    if (!(lower > 0.0) || !(lower < upper)) throw std::runtime_error("density domain is empty");
    std::size_t evaluations = 0;
    DensityRoot result;
    if (anchor_density.has_value()) {
        try {
            if (branch == DensityBranch::Unspecified) throw std::invalid_argument("density anchor requires explicit branch");
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
        } catch (const std::exception& error) {
            if (scope.diagnostics_on()) scope.probe("anchor_failure", {*anchor_density}, 0, error.what());
            throw;
        }
    } else {
        const auto roots = scan_roots(model, temperature, pressure, mole_fractions,
            std::log(lower), std::log(upper),
            static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value), scope, evaluations);
        std::vector<DensityRoot> stable;
        std::copy_if(roots.begin(), roots.end(), std::back_inserter(stable), [](const auto& root) {
            return std::isfinite(root.pressure_density_derivative) && root.pressure_density_derivative > 0.0;
        });
        if (stable.empty()) throw std::runtime_error("no stable density root (count 0)");
        std::sort(stable.begin(), stable.end(), [](const auto& left, const auto& right) {
            return left.molar_density < right.molar_density;
        });
        scope.probe("stable_root_count", static_cast<double>(stable.size()));
        if (stable.size() > 1 && branch == DensityBranch::Unspecified)
            throw std::runtime_error("density root is ambiguous (stable count "
                                     + std::to_string(stable.size()) + ")");
        result = branch == DensityBranch::Vapor ? stable.front() : stable.back();
        result.stable_root_count = stable.size();
        result.branch = stable.size() == 1 ? DensityBranch::Unspecified : branch;
        if (sibling_root != nullptr) {
            *sibling_root = branch == DensityBranch::Vapor ? stable.back() : stable.front();
            sibling_root->stable_root_count = stable.size();
            sibling_root->branch = stable.size() == 1 ? DensityBranch::Unspecified
                : branch == DensityBranch::Vapor ? DensityBranch::Liquid : DensityBranch::Vapor;
        }
    }
    scope.probe("density_root", result.molar_density);
    scope.finish(result.convergence);
    return result;
}

}  // namespace epcsaft
