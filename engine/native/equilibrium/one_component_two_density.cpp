#include "equilibrium/gibbs_minimizer.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

#include <Eigen/Dense>

namespace epcsaft::equilibrium {

namespace {
template <typename Base>
struct Reduction final {
    std::vector<Base> rows, jacobian;
    Base floor;
};

template <typename Base>
Reduction<Base> reduction(const std::vector<PureIsothermalValues<Base>>& phase,
    const std::vector<Base>& volume, const std::vector<Base>& u, std::size_t liquid, Base temperature) {
    const auto vapor = 1 - liquid;
    const auto& l = phase[liquid]; const auto& v = phase[vapor];
    const auto& dl = l.derivative; const auto& dv = v.derivative;
    const Base rt = Base(GAS_CONSTANT_J_PER_MOL_K) * temperature, ratio = std::exp(u[vapor] - u[liquid]);
    const Base fp = (l.pressure - v.pressure) * volume[liquid] / rt;
    using Sum = std::conditional_t<std::is_same_v<Base, double>, long double, Base>;
    const Base fm = static_cast<Base>((Sum(u[liquid]) - u[vapor]) + (Sum(dl.value) - dv.value)
        - (Sum(volume[liquid]) * dl.gradient[1] - Sum(volume[vapor]) * dv.gradient[1]));
    const Base kl = l.pressure_density_derivative / rt, kv = v.pressure_density_derivative / rt;
    const Base vl = volume[liquid], vv = volume[vapor];
    const Base scale = std::max(1 + std::abs(vl * dl.gradient[1]) + ratio * (1 + std::abs(vv * dv.gradient[1])),
        std::abs(u[liquid] - u[vapor]) + std::abs(dl.value) + std::abs(dv.value)
            + std::abs(vl * dl.gradient[1]) + std::abs(vv * dv.gradient[1]));
    Reduction<Base> out{{fp, fm}, {kl - fp, -ratio * kv, kl, -kv},
        std::numeric_limits<Base>::epsilon() * scale + std::max(vl * dl.gradient_floor + ratio * vv * dv.gradient_floor,
            vl * dl.gradient_floor + vv * dv.gradient_floor)};
    if (liquid != 0) { std::swap(out.jacobian[0], out.jacobian[1]); std::swap(out.jacobian[2], out.jacobian[3]); }
    return out;
}

template <typename Base>
struct ReducedCorrection final {
    using Vector = Eigen::Matrix<Base, 2, 1>;
    Vector delta;
    Base inverse, condition;
    explicit ReducedCorrection(const std::vector<Base>& jacobian, const std::vector<Base>& values) {
        using Matrix = Eigen::Matrix<Base, 2, 2>;
        const Matrix j = Eigen::Map<const Eigen::Matrix<Base, 2, 2, Eigen::RowMajor>>(jacobian.data());
        const Vector f(values[0], values[1]);
        if (!j.allFinite() || !f.allFinite()) throw std::runtime_error("one-component two-density nonfinite rows/Jacobian");
        const Vector columns = j.cwiseAbs().colwise().maxCoeff().transpose().unaryExpr(
            [](Base v) { return v > 0.0 ? 1.0 / v : 1.0; });
        const Matrix column_scaled = j * columns.asDiagonal();
        const Vector rows = column_scaled.cwiseAbs().rowwise().maxCoeff().unaryExpr(
            [](Base v) { return v > 0.0 ? 1.0 / v : 1.0; });
        const Eigen::FullPivLU<Matrix> lu(rows.asDiagonal() * column_scaled);
        if (lu.rank() != 2) throw std::runtime_error("one-component two-density reduced-Jacobian rank");
        delta = columns.asDiagonal() * lu.solve(-(rows.asDiagonal() * f));
        if (!delta.allFinite()) throw std::runtime_error("one-component two-density nonfinite direction");
        Vector sums = Vector::Zero();
        for (Eigen::Index i = 0; i < 2; ++i) {
            const Vector rhs = Vector::Unit(i);
            sums += (columns.asDiagonal() * lu.solve(rows.asDiagonal() * rhs)).cwiseAbs();
        }
        inverse = sums.maxCoeff(); condition = j.cwiseAbs().rowwise().sum().maxCoeff() * inverse;
    }
};

template <typename Base>
LinearSolveDiagnostics diagnostics(const Reduction<Base>& at, const ReducedCorrection<Base>& step, Base displacement = 0) {
    const Base norm = std::max(std::abs(at.rows[0]), std::abs(at.rows[1]));
    return {"one_component_two_density_reduction", "", double(norm), double(at.floor), double(step.inverse), double(step.condition),
        double(step.inverse * (norm + at.floor) + displacement), std::numeric_limits<Base>::digits, 0};
}
}  // namespace

bool one_component_two_density_eligible(const CompiledSystem& c) {
    if (c.ideal || !c.model || c.model->parameters().ideal_gas || !newton_eligible(c)
        || c.model->parameters().component_count() != 1 || c.problem.source_reference
        || !c.problem.solids.empty() || c.problem.temperature.free || !c.problem.temperature.is_positive()
        || !c.problem.pressure.free) return false;
    return c.problem.phases[0].support == std::vector<std::size_t>{0}
        && c.problem.phases[1].support == std::vector<std::size_t>{0};
}

// EqID one_component_two_density: eliminate normalization, feed, P and the common chemical multiplier.
OneComponentTwoDensityPoint one_component_two_density_point(const CompiledSystem& c, const std::vector<double>& u) {
    if (u.size() != 2 || !std::all_of(u.begin(), u.end(), [](double v) { return std::isfinite(v); }))
        throw std::invalid_argument("nonfinite pure log densities");
    const auto offsets = c.phase_offsets();
    auto point = c.generic_point();
    const double coefficient = c.model->packing_fraction_coefficient(c.problem.temperature.value, {1.0});
    std::size_t liquid = 0, vapor = 1;
    for (std::size_t p = 0; p < offsets.size(); ++p) {
        point[offsets[p]] = 0.0;
        if (c.problem.phases[p].is_free()) point[offsets[p] + 1] = c.feed_total;
        const auto eta = offsets[p] + 1 + (c.problem.phases[p].is_free() ? 1 : 0);
        point[eta] = std::log(coefficient) + u[p];
        (c.problem.phases[p].kind == PhaseKind::Liquid ? liquid : vapor) = p;
    }
    if (!(u[liquid] > u[vapor])) throw std::domain_error("one-component two-density density ordering");
    OneComponentTwoDensityPoint out;
    out.phase = c.prepare(point, true);
    auto& at = out.phase;
    const auto& v = at.pure[vapor];
    const auto reduced = reduction(at.pure, at.volume, u, liquid, at.temperature);
    out.rows = reduced.rows; out.jacobian = reduced.jacobian;
    at.pressure = at.point[c.primal_size - 1] = v.pressure;
    at.point[c.primal_size] = v.chemical_potential_over_rt + at.standard.derivative[0][0];
    return out;
}

EquilibriumResult one_component_two_density_attempt(const CompiledSystem& c, const std::vector<double>& initial,
    const EquilibriumResult& structural, core::Scope& scope, const DriverLimits* limits) {
    auto result = structural;
    result.message = "one-component two-density initialization";
    result.linear_solve_diagnostics.emplace();
    try {
        const auto offsets = c.phase_offsets();
        const double log_c = std::log(c.model->packing_fraction_coefficient(c.problem.temperature.value, {1.0}));
        std::vector<double> u;
        for (std::size_t p = 0; p < offsets.size(); ++p)
            u.push_back(initial[offsets[p] + 1 + (c.problem.phases[p].is_free() ? 1 : 0)] - log_c);
        GibbsMinimizer bounds(c, initial, &scope, limits);
        std::vector<double> lower(initial.size()), upper(initial.size());
        if (!bounds.physical_bounds(lower, upper)) {
            if (limits) limits->check();
            throw std::runtime_error("one-component two-density bounds unavailable");
        }
        const auto inside = [&](const PhasePoint& at) {
            for (std::size_t i = 0; i < lower.size(); ++i)
                if (!std::isfinite(at.point[i]) || at.point[i] < lower[i] || at.point[i] > upper[i]) return false;
            return true;
        };
        const auto evaluate = [&](const std::vector<double>& logs) {
            if (limits) {
                limits->check();
                if (!limits->spend(limits->evaluations, limits->evaluation_cap)) throw std::runtime_error("one-component two-density evaluation budget");
            }
            ++result.evaluations;
            return one_component_two_density_point(c, logs);
        };
        auto at = evaluate(u);
        const double tolerance = c.problem.solver_options.value(policy::row);
        bool polishing = false;
        double polish_residual = std::numeric_limits<double>::infinity();
        for (std::size_t step = 0; step < static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value); ++step) {
            if (limits && (limits->expired() || limits->consumed >= limits->iterations))
                throw std::runtime_error("one-component two-density iteration/elapsed budget");
            const auto full = c.evaluate_at(at.phase);
            polishing = polishing || maximum_abs(c.classified_residuals(at.phase, full.values)) <= tolerance;
            const bool final = polishing;
            const ReducedCorrection<double> correction(at.jacobian, at.rows);
            const auto& delta = correction.delta;
            const double forward = polish_forward_bound(at.jacobian, at.rows, {delta[0], delta[1]}, {1.0, 1.0});
            scope.probe("one_component_two_density_correction", {delta.lpNorm<Eigen::Infinity>(), forward,
                policy::EQUILIBRIUM_POLISH_CORRECTION.value}, step);
            ++result.iterations;
            if (limits) ++limits->consumed;
            bool accepted = false;
            std::string cause = "one-component two-density contraction";
            for (const double damping : {1.0, 0.5, 0.25}) {
                if (final && damping != 1.0) break;
                const std::vector<double> trial{u[0] + damping * delta[0], u[1] + damping * delta[1]};
                try {
                    auto candidate = evaluate(trial);
                    if (!inside(candidate.phase)) { cause = "one-component two-density bounds"; continue; }
                    if (!final && !(maximum_abs(candidate.rows) < maximum_abs(at.rows))) continue;
                    u = trial; at = std::move(candidate); accepted = true; break;
                } catch (const std::exception& error) { cause = error.what(); }
            }
            if (!accepted) throw std::runtime_error(cause);
            if (!final) continue;
            std::size_t liquid = c.problem.phases[0].kind == PhaseKind::Liquid ? 0 : 1;
            const auto reduced = reduction(at.phase.pure, at.phase.volume, u, liquid, at.phase.temperature);
            const ReducedCorrection<double> estimate(at.jacobian, at.rows);
            result.linear_solve_diagnostics = diagnostics(reduced, estimate);
            scope.probe("one_component_two_density_error", {estimate.inverse, result.linear_solve_diagnostics->forward_error_estimate, reduced.floor});
            const double achieved = maximum_abs(at.rows);
            const bool contracting = achieved < polish_residual;
            polish_residual = achieved;
            if (estimate.delta.lpNorm<Eigen::Infinity>() > policy::EQUILIBRIUM_POLISH_CORRECTION.value
                && achieved > reduced.floor && contracting) continue;
            if (result.linear_solve_diagnostics->forward_error_estimate > policy::TWO_DENSITY_FORWARD_ERROR.value
                || estimate.delta.lpNorm<Eigen::Infinity>() > policy::TWO_DENSITY_FORWARD_ERROR.value) {
                if constexpr (std::numeric_limits<long double>::digits <= std::numeric_limits<double>::digits) {
                    result.linear_solve_diagnostics->reason = "long double has no additional significand bits";
                } else {
                    using Wide = long double;
                    std::vector<Wide> logs{Wide(u[0]), Wide(u[1])};
                    const Wide temperature = c.problem.temperature.value;
                    const auto wide_evaluate = [&](const std::vector<Wide>& values) {
                        core::Scope extended("one_component_two_density_extended_evaluate");
                        if (!(values[liquid] > values[1 - liquid])) throw std::runtime_error("one-component two-density density ordering");
                        if (limits) { limits->check(); if (!limits->spend(limits->evaluations, limits->evaluation_cap))
                            throw std::runtime_error("one-component two-density evaluation budget"); }
                        ++result.evaluations;
                        std::vector<Wide> volume; std::vector<PureIsothermalValues<Wide>> phases;
                        for (Wide value : values) {
                            volume.push_back(1 / std::exp(value));
                            phases.push_back(c.model->pure_isothermal(temperature, volume.back()));
                        }
                        return std::tuple{reduction(phases, volume, values, liquid, temperature), phases[1 - liquid]};
                    };
                    const auto lift = [&](const std::vector<Wide>& values, const PureIsothermalValues<Wide>& vapor) {
                        PhasePoint point; point.point = at.phase.point;
                        for (std::size_t p = 0; p < 2; ++p)
                            point.point[offsets[p] + 1 + (c.problem.phases[p].is_free() ? 1 : 0)] = double(Wide(log_c) + values[p]);
                        point.point[c.primal_size - 1] = double(vapor.pressure);
                        point.point[c.primal_size] = double(vapor.chemical_potential_over_rt + Wide(at.phase.standard.derivative[0][0]));
                        return point;
                    };
                    auto [wide, vapor] = wide_evaluate(logs);
                    std::size_t refinements = 0;
                    for (;;) {
                        if (result.iterations >= static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)
                            || (limits && (limits->expired() || limits->consumed >= limits->iterations)))
                            throw std::runtime_error("one-component two-density extended iteration/elapsed budget");
                        const ReducedCorrection<Wide> correction(wide.jacobian, wide.rows);
                        const Wide norm = std::max(std::abs(wide.rows[0]), std::abs(wide.rows[1]));
                        const bool final_wide = norm <= Wide(tolerance);
                        bool accepted_wide = false;
                        ++result.iterations; ++refinements; if (limits) ++limits->consumed;
                        for (Wide damping : {Wide(1), Wide(0.5), Wide(0.25)}) {
                            if (final_wide && damping != 1) break;
                            const std::vector<Wide> trial{logs[0] + damping * correction.delta[0], logs[1] + damping * correction.delta[1]};
                            auto [next, next_vapor] = wide_evaluate(trial);
                            if (!inside(lift(trial, next_vapor))) continue;
                            if (!final_wide && !(std::max(std::abs(next.rows[0]), std::abs(next.rows[1])) < norm)) continue;
                            logs = trial; wide = std::move(next); vapor = std::move(next_vapor); accepted_wide = true; break;
                        }
                        if (!accepted_wide) throw std::runtime_error("one-component two-density extended contraction/bounds");
                        const ReducedCorrection<Wide> next(wide.jacobian, wide.rows);
                        auto record = diagnostics(wide, next); record.refinement_iterations = refinements;
                        result.linear_solve_diagnostics = record;
                        const Wide achieved = std::max(std::abs(wide.rows[0]), std::abs(wide.rows[1]));
                        const Wide error = next.inverse * (achieved + wide.floor), correction_norm = next.delta.template lpNorm<Eigen::Infinity>();
                        if (!final_wide || achieved > Wide(tolerance)
                            || !((correction_norm <= policy::EQUILIBRIUM_POLISH_CORRECTION.value && error <= policy::EQUILIBRIUM_POLISH_CORRECTION.value)
                                || (achieved <= wide.floor && correction_norm <= next.inverse * wide.floor))) continue;
                        if (!(logs[liquid] - logs[1 - liquid] > 2 * error)) throw std::runtime_error("one-component two-density unresolved density intervals");
                        if (limits && !limits->spend(limits->evaluations, limits->evaluation_cap)) throw std::runtime_error("one-component two-density evaluation budget");
                        ++result.evaluations;
                        at.phase = c.prepare(lift(logs, vapor).point, true);
                        std::vector<Wide> exported; Wide displacement = 0;
                        for (std::size_t p = 0; p < 2; ++p) {
                            exported.push_back(std::log(Wide(1.0 / at.phase.volume[p])));
                            displacement = std::max(displacement, std::abs(exported.back() - logs[p]));
                        }
                        const auto [output, output_vapor] = wide_evaluate(exported);
                        const ReducedCorrection<Wide> output_step(output.jacobian, output.rows);
                        result.linear_solve_diagnostics = diagnostics(output, output_step, displacement);
                        result.linear_solve_diagnostics->refinement_iterations = refinements;
                        break;
                    }
                }
            }
            const auto endpoint = c.evaluate_at(at.phase);
            result.residuals = endpoint.values;
            result.classified_residuals = c.classified_residuals(at.phase, endpoint.values);
            result.objective = endpoint.objective;
            core::Convergence& convergence = result.convergence;
            convergence = {};
            convergence.residual = maximum_abs(result.classified_residuals);
            core::stopping_criterion(convergence, c.problem.solver_options, policy::row);
            if (endpoint.floor > 0.0) convergence.floor = endpoint.floor;
            convergence.iterations = result.iterations; convergence.evaluations = result.evaluations;
            convergence.stop = core::classify(convergence.residual, *convergence.criterion, convergence.floor,
                convergence.contraction, false, convergence.amplification);
            result.requested_tolerance_met = convergence.residual <= tolerance;
            result.success = convergence.stop == core::Stop::converged && inside(at.phase);
            result.message = core::describe(convergence);
            result.solver_status = "NEWTON_CONVERGED";
            result.constraint_multipliers.assign(c.rows.size(), 0.0);
            scope.finish(convergence);
            result = publish_candidate(c, at.phase, std::move(result), scope, limits);
            if (result.success) {
                const auto branch = check_density_branches(*c.model, c, result, {}, true, scope, limits, &at.phase);
                if (!branch.accepted) { result.success = false; result.message = branch.message; clear_state(result); }
            }
            return result;
        }
        result.message = "one-component two-density iteration budget";
    } catch (const std::exception& error) { result.message = error.what(); }
    scope.probe("one_component_two_density_rejected", {static_cast<double>(result.iterations), static_cast<double>(result.evaluations)}, 0, result.message.c_str());
    if (result.linear_solve_diagnostics) result.linear_solve_diagnostics->reason = "rejected attempt: " + result.message;
    result.success = false; clear_state(result);
    return result;
}
}  // namespace epcsaft::equilibrium
