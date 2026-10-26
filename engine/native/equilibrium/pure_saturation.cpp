#include "equilibrium/gibbs_minimizer.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

#include <Eigen/Dense>

namespace epcsaft::equilibrium {

bool pure_saturation_eligible(const CompiledSystem& c) {
    if (c.ideal || !c.model || c.model->parameters().ideal_gas || !newton_eligible(c)
        || c.model->parameters().component_count() != 1 || c.problem.source_reference
        || !c.problem.solids.empty() || c.problem.temperature.free || !c.problem.temperature.is_positive()
        || !c.problem.pressure.free) return false;
    return c.problem.phases[0].support == std::vector<std::size_t>{0}
        && c.problem.phases[1].support == std::vector<std::size_t>{0};
}

// EqID pure_saturation_reduction: eliminate normalization, feed, P and the common chemical multiplier.
PureSaturationPoint pure_saturation_point(const CompiledSystem& c, const std::vector<double>& u) {
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
    if (!(u[liquid] > u[vapor])) throw std::domain_error("pure saturation density ordering");
    PureSaturationPoint out;
    out.phase = c.prepare(point, true);
    auto& at = out.phase;
    const auto& l = at.pure[liquid]; const auto& v = at.pure[vapor];
    const double rt = GAS_CONSTANT_J_PER_MOL_K * at.temperature;
    const double ratio = std::exp(u[vapor] - u[liquid]);
    const double fp = (l.pressure - v.pressure) * at.volume[liquid] / rt;
    const auto& dl = l.derivative; const auto& dv = v.derivative;
    const double fm = static_cast<double>((static_cast<long double>(u[liquid]) - u[vapor])
        + (static_cast<long double>(dl.value) - dv.value)
        - (static_cast<long double>(at.volume[liquid]) * dl.gradient[1]
            - static_cast<long double>(at.volume[vapor]) * dv.gradient[1]));
    const double kl = l.pressure_density_derivative / rt, kv = v.pressure_density_derivative / rt;
    out.rows = {fp, fm};
    out.jacobian = {kl - fp, -ratio * kv, kl, -kv};
    if (liquid != 0) { std::swap(out.jacobian[0], out.jacobian[1]); std::swap(out.jacobian[2], out.jacobian[3]); }
    at.pressure = at.point[c.primal_size - 1] = v.pressure;
    at.point[c.primal_size] = v.chemical_potential_over_rt + at.standard.derivative[0][0];
    return out;
}

EquilibriumResult pure_saturation_attempt(const CompiledSystem& c, const std::vector<double>& initial,
    const EquilibriumResult& structural, core::Scope& scope, const DriverLimits* limits) {
    auto result = structural;
    result.message = "pure saturation initialization";
    try {
        const auto offsets = c.phase_offsets();
        const double log_c = std::log(c.model->packing_fraction_coefficient(c.problem.temperature.value, {1.0}));
        std::vector<double> u;
        for (std::size_t p = 0; p < offsets.size(); ++p)
            u.push_back(initial[offsets[p] + 1 + (c.problem.phases[p].is_free() ? 1 : 0)] - log_c);
        GibbsMinimizer bounds(c, initial, &scope, limits);
        std::vector<double> lower(initial.size()), upper(initial.size());
        if (!bounds.physical_bounds(lower, upper)) throw std::runtime_error("pure saturation bounds unavailable");
        const auto inside = [&](const PhasePoint& at) {
            for (std::size_t i = 0; i < lower.size(); ++i)
                if (!std::isfinite(at.point[i]) || at.point[i] < lower[i] || at.point[i] > upper[i]) return false;
            return true;
        };
        const auto evaluate = [&](const std::vector<double>& logs) {
            if (limits) {
                limits->check();
                if (!limits->spend(limits->evaluations, limits->evaluation_cap)) throw std::runtime_error("pure saturation evaluation budget");
            }
            ++result.evaluations;
            return pure_saturation_point(c, logs);
        };
        auto at = evaluate(u);
        const double tolerance = c.problem.solver_options.value(policy::row);
        for (std::size_t step = 0; step < static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value); ++step) {
            if (limits && (limits->expired() || limits->consumed >= limits->iterations))
                throw std::runtime_error("pure saturation iteration/elapsed budget");
            const auto full = c.evaluate_at(at.phase);
            const bool final = maximum_abs(c.classified_residuals(at.phase, full.values)) <= tolerance;
            using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            const Eigen::MatrixXd j = Eigen::Map<const Matrix>(at.jacobian.data(), 2, 2);
            const Eigen::Vector2d f(at.rows[0], at.rows[1]);
            if (!j.allFinite() || !f.allFinite()) throw std::runtime_error("pure saturation nonfinite rows/Jacobian");
            const Eigen::Vector2d columns = j.cwiseAbs().colwise().maxCoeff().transpose().unaryExpr(
                [](double v) { return v > 0.0 ? 1.0 / v : 1.0; });
            const Eigen::Matrix2d column_scaled = j * columns.asDiagonal();
            const Eigen::Vector2d rows = column_scaled.cwiseAbs().rowwise().maxCoeff().unaryExpr(
                [](double v) { return v > 0.0 ? 1.0 / v : 1.0; });
            const Eigen::FullPivLU<Eigen::Matrix2d> lu(rows.asDiagonal() * column_scaled);
            scope.probe("pure_saturation_condition", {static_cast<double>(lu.rank()), lu.rcond(),
                std::abs(u[0] - u[1]), j.row(1).cwiseAbs().minCoeff()}, step);
            if (lu.rank() != 2) throw std::runtime_error("pure saturation reduced-Jacobian rank");
            const Eigen::Vector2d delta = columns.asDiagonal() * lu.solve(-(rows.asDiagonal() * f));
            if (!delta.allFinite()) throw std::runtime_error("pure saturation nonfinite direction");
            const double forward = polish_forward_bound(at.jacobian, at.rows, {delta[0], delta[1]}, {1.0, 1.0});
            scope.probe("pure_saturation_correction", {delta.lpNorm<Eigen::Infinity>(), forward,
                policy::EQUILIBRIUM_POLISH_CORRECTION.value}, step);
            ++result.iterations;
            if (limits) ++limits->consumed;
            bool accepted = false;
            std::string cause = "pure saturation contraction";
            for (const double damping : {1.0, 0.5, 0.25}) {
                if (final && damping != 1.0) break;
                const std::vector<double> trial{u[0] + damping * delta[0], u[1] + damping * delta[1]};
                try {
                    auto candidate = evaluate(trial);
                    if (!inside(candidate.phase)) { cause = "pure saturation bounds"; continue; }
                    if (!final && !(maximum_abs(candidate.rows) < maximum_abs(at.rows))) continue;
                    u = trial; at = std::move(candidate); accepted = true; break;
                } catch (const std::exception& error) { cause = error.what(); }
            }
            if (!accepted) throw std::runtime_error(cause);
            if (!final) continue;
            const auto endpoint = c.evaluate_at(at.phase);
            result.residuals = endpoint.values;
            result.classified_residuals = c.classified_residuals(at.phase, endpoint.values);
            result.objective = endpoint.objective;
            auto& convergence = result.convergence;
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
        result.message = "pure saturation iteration budget";
    } catch (const std::exception& error) { result.message = error.what(); }
    scope.probe("pure_saturation_rejected", {static_cast<double>(result.iterations), static_cast<double>(result.evaluations)}, 0, result.message.c_str());
    result.success = false; clear_state(result);
    return result;
}
}  // namespace epcsaft::equilibrium
