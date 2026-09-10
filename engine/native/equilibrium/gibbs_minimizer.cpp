#include "equilibrium/gibbs_minimizer.hpp"

#include "core/trace.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <IpIpoptApplication.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>

namespace epcsaft::equilibrium {
namespace {

std::size_t phase_width(const CompiledSystem& compiled) {
    std::size_t result = 0;
    for (const auto& phase : compiled.problem.phases) {
        result += phase.support.size() + 1;
        if (phase.amount.kind == AmountKind::Free) ++result;
    }
    return result;
}

std::vector<std::size_t> phase_offsets(const CompiledSystem& compiled) {
    std::vector<std::size_t> result;
    std::size_t offset = 0;
    for (const auto& phase : compiled.problem.phases) {
        result.push_back(offset);
        offset += phase.support.size() + 1;
        if (phase.amount.kind == AmountKind::Free) ++offset;
    }
    return result;
}

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) result = std::max(result, std::abs(value));
    return result;
}

double specified_amount(const PhaseSpec& phase, const std::vector<double>& feed) {
    if (phase.amount.kind == AmountKind::Free) return 1.0;
    double total = 0.0;
    for (const double value : feed) total += value;
    return phase.amount.fraction * total;
}

double reduced_hessian_minimum_at(
    const CompiledSystem& compiled, const std::vector<double>& point,
    const std::vector<double>& multipliers
) {
    if (compiled.pinned_count != 0 || point.empty()) return 0.0;
    const auto rows = [&compiled] {
        std::vector<std::size_t> result;
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            if (compiled.rows[row].kind != RowKind::Chemical
                && compiled.rows[row].kind != RowKind::Pressure) result.push_back(row);
        }
        return result;
    }();
    const auto jacobian = compiled.jacobian(point);
    Eigen::MatrixXd constraints(static_cast<Eigen::Index>(rows.size()),
                                static_cast<Eigen::Index>(point.size()));
    for (std::size_t row = 0; row < rows.size(); ++row) {
        for (std::size_t column = 0; column < point.size(); ++column) {
            constraints(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                jacobian[rows[row] * point.size() + column];
        }
    }
    const auto svd = constraints.jacobiSvd(Eigen::ComputeFullV);
    const Eigen::Index rank = svd.rank();
    if (rank >= static_cast<Eigen::Index>(point.size())) return 1.0;
    const Eigen::MatrixXd nullspace = svd.matrixV().rightCols(
        static_cast<Eigen::Index>(point.size()) - rank);
    const auto exact = compiled.hessian(point, multipliers);
    Eigen::MatrixXd hessian(static_cast<Eigen::Index>(point.size()),
                            static_cast<Eigen::Index>(point.size()));
    for (std::size_t row = 0; row < point.size(); ++row)
        for (std::size_t column = 0; column < point.size(); ++column)
            hessian(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                exact[row * point.size() + column];
    const Eigen::MatrixXd reduced = nullspace.transpose() * hessian * nullspace;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(reduced);
    if (eigen.info() != Eigen::Success || eigen.eigenvalues().size() == 0) return 0.0;
    return eigen.eigenvalues()(0);
}

double largest_difference(
    const std::vector<double>& left, const std::vector<double>& right
) {
    if (left.size() != right.size()) return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result = std::max(result, std::abs(left[index] - right[index]));
    }
    return result;
}

bool numeric_state(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

double michelsen_deficit(
    const CompiledSystem& compiled,
    const std::vector<double>& composition,
    std::size_t phase,
    const std::vector<double>& residuals
) {
    double sum = 0.0;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        const auto& item = compiled.rows[row];
        if (item.kind != RowKind::Chemical || item.phase != phase) continue;
        if (item.index >= composition.size() || row >= residuals.size()) continue;
        sum += composition[item.index] * std::exp(-residuals[row]);
    }
    return 1.0 - sum;
}

void clear_state(EquilibriumResult& result) {
    result.temperature.clear();
    result.pressure.clear();
    result.amounts.clear();
    result.mole_fractions.clear();
    result.molar_densities.clear();
    result.packing_fractions.clear();
    result.multipliers.clear();
    result.constraint_multipliers.clear();
    result.coordinates.clear();
}

void require_solution_shape(
    const EquilibriumResult& result, std::size_t phases, std::size_t components
) {
    if (!numeric_state(result.residuals) || result.rows.empty()) {
        throw std::runtime_error("equilibrium result has invalid residuals");
    }
    if (result.amounts.size() != phases || result.molar_densities.size() != phases
        || result.mole_fractions.size() != phases * components) {
        throw std::runtime_error("equilibrium result has invalid topology");
    }
    for (const double amount : result.amounts) {
        if (!std::isfinite(amount) || amount < 0.0) {
            throw std::runtime_error("equilibrium result has invalid amount");
        }
    }
    for (const double density : result.molar_densities) {
        if (!std::isfinite(density) || density <= 0.0) {
            throw std::runtime_error("equilibrium result has invalid density");
        }
    }
    for (std::size_t phase = 0; phase < phases; ++phase) {
        double total = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            const double value = result.mole_fractions[phase * components + component];
            if (!std::isfinite(value) || value < 0.0) {
                throw std::runtime_error("equilibrium result has invalid composition");
            }
            total += value;
        }
        if (std::abs(total - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value) {
            throw std::runtime_error("equilibrium result is not normalized");
        }
    }
}

std::vector<double> add_multipliers(
    const CompiledSystem& compiled, const std::vector<double>& point
) {
    if (point.size() != compiled.primal_size) return point;
    const std::size_t width = compiled.rank + compiled.charge_row_count;
    std::vector<double> result = point;
    result.resize(point.size() + width, 0.0);
    if (width == 0) return result;
    const auto values = compiled.evaluate(result).values;
    std::vector<std::size_t> chemical;
    std::vector<std::size_t> phi_phase;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.rows[row].kind == RowKind::Chemical) chemical.push_back(row);
        if (compiled.rows[row].kind == RowKind::Electro) phi_phase.push_back(compiled.rows[row].phase);
    }
    Eigen::MatrixXd matrix(
        static_cast<Eigen::Index>(chemical.size()), static_cast<Eigen::Index>(width)
    );
    Eigen::VectorXd rhs(static_cast<Eigen::Index>(chemical.size()));
    const auto z = species_z(compiled.model->parameters());
    for (std::size_t row = 0; row < chemical.size(); ++row) {
        const auto& item = compiled.rows[chemical[row]];
        rhs(static_cast<Eigen::Index>(row)) = values[chemical[row]];
        for (std::size_t column = 0; column < compiled.rank; ++column) {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                compiled.B[column][item.index];
        }
        for (std::size_t column = 0; column < phi_phase.size(); ++column) {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(compiled.rank + column)) =
                z[item.index] * (item.phase == phi_phase[column] ? 1.0 : 0.0);
        }
    }
    const Eigen::VectorXd solved = matrix.colPivHouseholderQr().solve(rhs);
    for (std::size_t column = 0; column < width; ++column) {
        result[point.size() + column] = solved(static_cast<Eigen::Index>(column));
    }
    return result;
}

bool valid_previous_phase(
    const EquilibriumResult& start, std::size_t phase, std::size_t components
);

void copy_previous_phase(
    const CompiledSystem& compiled, const EquilibriumResult& start,
    std::size_t phase, std::vector<double>& point,
    const std::vector<std::size_t>& offsets
);

std::vector<double> start_point(
    const CompiledSystem& compiled, const EquilibriumResult& start,
    std::vector<double> point
) {
    const std::size_t components = compiled.model->parameters().component_count();
    if (start.phase_count != compiled.problem.phases.size()) return point;
    const auto offsets = phase_offsets(compiled);
    for (std::size_t phase = 0; phase < start.phase_count; ++phase) {
        if (valid_previous_phase(start, phase, components)) {
            copy_previous_phase(compiled, start, phase, point, offsets);
        }
    }
    std::size_t intensive = phase_width(compiled);
    if (compiled.problem.temperature.free && !start.temperature.empty()) point[intensive++] = start.temperature[0];
    if (compiled.problem.pressure.free && !start.pressure.empty()) point[intensive++] = start.pressure[0];
    if (compiled.pinned_count != 0 && start.multipliers.size() == compiled.rank + compiled.charge_row_count) {
        std::copy(start.multipliers.begin(), start.multipliers.end(), point.begin() + compiled.primal_size);
    }
    return point;
}

bool valid_previous_phase(
    const EquilibriumResult& start, std::size_t phase, std::size_t components
) {
    if (!start.success || phase >= start.phase_count
        || start.mole_fractions.size() != start.phase_count * components
        || start.molar_densities.size() != start.phase_count) return false;
    if (!(std::isfinite(start.molar_densities[phase])
          && start.molar_densities[phase] > 0.0)) return false;
    for (std::size_t species = 0; species < components; ++species) {
        const double value = start.mole_fractions[phase * components + species];
        if (!std::isfinite(value) || value < 0.0) return false;
    }
    return true;
}

void copy_previous_phase(
    const CompiledSystem& compiled, const EquilibriumResult& start,
    std::size_t phase, std::vector<double>& point,
    const std::vector<std::size_t>& offsets
) {
    const auto& spec = compiled.problem.phases[phase];
    for (std::size_t i = 0; i < spec.support.size(); ++i) {
        point[offsets[phase] + i] = std::log(start.mole_fractions[
            phase * compiled.model->parameters().component_count() + spec.support[i]]);
    }
    std::size_t cursor = offsets[phase] + spec.support.size();
    if (spec.amount.kind == AmountKind::Free && start.amounts.size() == start.phase_count) {
        if (std::isfinite(start.amounts[phase]) && start.amounts[phase] > 0.0) {
            point[cursor] = start.amounts[phase];
        }
        ++cursor;
    }
    std::vector<double> composition(compiled.model->parameters().component_count(), 0.0);
    for (std::size_t species = 0; species < composition.size(); ++species)
        composition[species] = start.mole_fractions[phase * composition.size() + species];
    double density = start.molar_densities[phase];
    if (spec.amount.kind == AmountKind::Pinned && !start.temperature.empty()
        && !start.pressure.empty()) {
        const auto branch = spec.kind == PhaseKind::Vapor
            ? DensityBranch::Vapor : DensityBranch::Liquid;
        density = solve_density(*compiled.model, start.temperature[0], start.pressure[0],
                                composition, branch, std::nullopt).molar_density;
    }
    point[cursor] = std::log(compiled.model->packing_fraction_coefficient(
        start.temperature.empty() ? compiled.problem.temperature.value : start.temperature[0],
        composition) * density);
}

double active_constraint_residual(
    const CompiledSystem& compiled, const std::vector<double>& point
) {
    const auto values = compiled.evaluate(point).values;
    double result = 0.0;
    for (std::size_t row = 0; row < values.size(); ++row) {
        if (compiled.pinned_count == 0 && compiled.rows[row].kind == RowKind::Chemical) continue;
        result = std::max(result, std::abs(values[row]));
    }
    return result;
}

struct ResidualSummary final {
    double maximum = 0.0;
    double intensive = 0.0;
    double balance = 0.0;
    double normalization = 0.0;
    double charge = 0.0;
};

ResidualSummary summarize_residuals(
    const std::vector<Row>& rows, const std::vector<double>& values
) {
    ResidualSummary result;
    const std::size_t count = std::min(rows.size(), values.size());
    for (std::size_t index = 0; index < count; ++index) {
        const double magnitude = std::abs(values[index]);
        result.maximum = std::max(result.maximum, magnitude);
        if (rows[index].kind == RowKind::Chemical || rows[index].kind == RowKind::Pressure) {
            result.intensive = std::max(result.intensive, magnitude);
        } else if (rows[index].kind == RowKind::Moiety) {
            result.balance = std::max(result.balance, magnitude);
        } else if (rows[index].kind == RowKind::Normalization) {
            result.normalization = std::max(result.normalization, magnitude);
        } else if (rows[index].kind == RowKind::Electro) {
            result.charge = std::max(result.charge, magnitude);
        }
    }
    return result;
}

}  // namespace

GibbsMinimizer::GibbsMinimizer(
    const CompiledSystem& compiled, std::vector<double> start, core::Scope* scope
)
    : compiled_(compiled), start_(std::move(start)), scope_(scope) {}

bool GibbsMinimizer::get_nlp_info(
    Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style
) {
    n = static_cast<Ipopt::Index>(compiled_.pinned_count == 0
        ? compiled_.primal_size : compiled_.primal_size + compiled_.rank + compiled_.charge_row_count);
    m = static_cast<Ipopt::Index>(constraint_rows().size());
    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;
    index_style = C_STYLE;
    return true;
}

bool GibbsMinimizer::get_bounds_info(
    Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
    Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u
) {
    const double infinity = 1.0e19;
    std::fill(x_l, x_l + n, -infinity);
    std::fill(x_u, x_u + n, infinity);
    std::size_t offset = 0;
    for (const auto& phase : compiled_.problem.phases) {
        for (std::size_t i = 0; i < phase.support.size(); ++i) {
            x_l[offset] = std::log(policy::STABILITY_COMPOSITION_FLOOR.value);
            x_u[offset++] = 0.0;
        }
        if (phase.amount.kind == AmountKind::Free) x_l[offset++] = 0.0;
        x_l[offset] = std::log(1.0e-12);
        x_u[offset++] = std::log(policy::PACKING_FRACTION_CEILING.value);
    }
    if (compiled_.problem.temperature.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    if (compiled_.problem.pressure.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    const auto rows = constraint_rows();
    std::fill(g_l, g_l + m, 0.0);
    std::fill(g_u, g_u + m, 0.0);
    (void)rows;
    return true;
}

bool GibbsMinimizer::get_starting_point(
    Ipopt::Index n, bool init_x, Ipopt::Number* x,
    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda
) {
    if (!init_x) return false;
    const std::vector<double> point = compiled_.pinned_count == 0
        ? std::vector<double>(start_.begin(), start_.begin() + compiled_.primal_size)
        : start_;
    for (Ipopt::Index i = 0; i < n; ++i) x[i] = i < static_cast<Ipopt::Index>(point.size()) ? point[i] : 0.0;
    if (init_z) {
        std::fill(z_L, z_L + n, 0.0);
        std::fill(z_U, z_U + n, 0.0);
    }
    if (init_lambda) std::fill(lambda, lambda + m, 0.0);
    return true;
}

std::vector<double> GibbsMinimizer::point(Ipopt::Index n, const Ipopt::Number* x) const {
    return std::vector<double>(x, x + n);
}

std::vector<std::size_t> GibbsMinimizer::constraint_rows() const {
    std::vector<std::size_t> result;
    for (std::size_t row = 0; row < compiled_.rows.size(); ++row) {
        if (compiled_.pinned_count != 0
            || (compiled_.rows[row].kind != RowKind::Chemical
                && compiled_.rows[row].kind != RowKind::Pressure)) {
            result.push_back(row);
        }
    }
    return result;
}

bool GibbsMinimizer::eval_f(
    Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& value
) {
    (void)new_x;
    ++evaluations_;
    value = compiled_.pinned_count == 0 ? compiled_.evaluate(point(n, x)).objective : 0.0;
    return std::isfinite(value);
}

std::vector<double> GibbsMinimizer::objective_gradient(const std::vector<double>& x) const {
    return compiled_.objective_gradient(x);
}

bool GibbsMinimizer::eval_grad_f(
    Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f
) {
    (void)new_x;
    const auto gradient = objective_gradient(point(n, x));
    std::copy(gradient.begin(), gradient.end(), grad_f);
    return true;
}

bool GibbsMinimizer::eval_g(
    Ipopt::Index n, const Ipopt::Number* x, bool new_x,
    Ipopt::Index m, Ipopt::Number* g
) {
    (void)new_x;
    const auto values = compiled_.evaluate(point(n, x)).values;
    const auto rows = constraint_rows();
    for (Ipopt::Index row = 0; row < m; ++row) g[row] = values[rows[row]];
    return true;
}

bool GibbsMinimizer::eval_jac_g(
    Ipopt::Index n, const Ipopt::Number* x, bool new_x,
    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    (void)new_x;
    const auto rows = constraint_rows();
    if (values == nullptr) {
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < m; ++row) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                iRow[entry] = row;
                jCol[entry++] = column;
            }
        }
        return entry == nele_jac;
    }
    const auto jacobian = compiled_.jacobian(point(n, x));
    Ipopt::Index entry = 0;
    for (const std::size_t row : rows) {
        for (Ipopt::Index column = 0; column < n; ++column) {
            values[entry++] = jacobian[row * static_cast<std::size_t>(n) + static_cast<std::size_t>(column)];
        }
    }
    return entry == nele_jac;
}

bool GibbsMinimizer::eval_h(
    Ipopt::Index n, const Ipopt::Number* x, bool new_x,
    Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
    bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    (void)new_x; (void)new_lambda;
    (void)m; (void)lambda;
    ++hessian_calls_;
    if (values == nullptr) {
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < n; ++row) {
            for (Ipopt::Index column = 0; column <= row; ++column) {
                iRow[entry] = row;
                jCol[entry++] = column;
            }
        }
        return entry == nele_hess;
    }
    std::fill(values, values + nele_hess, 0.0);
    const auto base = point(n, x);
    const auto objective_hessian = compiled_.hessian(base);
    const auto lagrangian_hessian = compiled_.hessian(
        base, std::vector<double>(lambda, lambda + m));
    Ipopt::Index entry = 0;
    for (Ipopt::Index row = 0; row < n; ++row) {
        for (Ipopt::Index column = 0; column <= row; ++column) {
            const std::size_t index = static_cast<std::size_t>(row) * n
                + static_cast<std::size_t>(column);
            values[entry++] = obj_factor * objective_hessian[index]
                + lagrangian_hessian[index] - objective_hessian[index];
        }
    }
    return entry == nele_hess;
}

bool GibbsMinimizer::intermediate_callback(
    Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
    Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
    Ipopt::Number d_norm, Ipopt::Number regularization_size,
    Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
    const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq
) {
    (void)mode; (void)obj_value; (void)ip_data; (void)ip_cq;
    iterations_ = static_cast<std::size_t>(iter);
    regularization_size_ = regularization_size;
    if (scope_ != nullptr) scope_->probe("ipopt_iteration",
        std::vector<double>{inf_pr, inf_du, mu, d_norm, regularization_size,
                            alpha_pr, alpha_du, static_cast<double>(ls_trials)},
        static_cast<std::uint32_t>(iter));
    return true;
}

void GibbsMinimizer::finalize_solution(
    Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
    const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
    const Ipopt::Number* g, const Ipopt::Number* lambda,
    Ipopt::Number objective,
    const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq
) {
    (void)z_L; (void)z_U; (void)g;
    (void)objective; (void)ip_data; (void)ip_cq;
    status_ = status;
    solution_.assign(x, x + n);
    constraint_multipliers_.assign(lambda, lambda + m);
}

EquilibriumResult solve_equilibrium(
    const HelmholtzModel& model,
    const EquilibriumProblem& request,
    const std::optional<EquilibriumResult>& start
) {
    core::Scope scope("solve_equilibrium");
    CompiledSystem compiled = compile_problem(model, request);
    std::vector<double> initial = boundary_start(compiled);
    if (start.has_value()) initial = start_point(compiled, *start, std::move(initial));
    Ipopt::SmartPtr<GibbsMinimizer> nlp = new GibbsMinimizer(compiled, std::move(initial), &scope);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = ::IpoptApplicationFactory();
    const double objective_scale = 1.0 / compiled.feed_total;
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("obj_scaling_factor", objective_scale);
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetNumericValue("tol", policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    application->Options()->SetNumericValue("acceptable_tol", policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue("dual_inf_tol", policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    application->Options()->SetNumericValue("constr_viol_tol", policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    application->Options()->SetIntegerValue("max_iter",
        static_cast<int>(policy::EQUILIBRIUM_MAX_ITERATIONS.value));
    application->Initialize();
    application->OptimizeTNLP(nlp);
    EquilibriumResult result;
    result.rows = compiled.rows;
    result.phase_count = compiled.problem.phases.size();
    std::vector<double> point = nlp->solution();
    if (point.empty()) point = boundary_start(compiled);
    std::vector<double> full = compiled.pinned_count == 0 ? add_multipliers(compiled, point) : point;
    const auto evaluation = compiled.evaluate(full);
    result.residuals = evaluation.values;
    result.objective = evaluation.objective;
    const auto residuals = summarize_residuals(result.rows, result.residuals);
    scope.probe("maximum_residual", residuals.maximum);
    scope.probe("intensive_residual", residuals.intensive);
    scope.probe("balance_residual", residuals.balance);
    scope.probe("normalization_residual", residuals.normalization);
    scope.probe("charge_residual", residuals.charge);
    scope.probe("row_count", static_cast<double>(result.rows.size()));
    scope.probe("primal_size", static_cast<double>(compiled.primal_size));
    scope.probe("balance_rank", static_cast<double>(compiled.rank));
    scope.probe("pinned_count", static_cast<double>(compiled.pinned_count));
    scope.probe("active_residual", active_constraint_residual(compiled, full));
    scope.probe("objective", result.objective);
    scope.probe("hessian_calls", static_cast<double>(nlp->hessian_calls()));
    core::Convergence& convergence = result.convergence;
    convergence = {};
    convergence.residual = maximum_abs(result.residuals);
    convergence.criterion = &policy::EQUILIBRIUM_ROW_RESIDUAL;
    if (evaluation.floor > 0.0) convergence.floor = evaluation.floor;
    convergence.amplification = nlp->regularization_size();
    result.iterations = nlp->iterations();
    result.evaluations = nlp->evaluations();
    result.hessian_calls = nlp->hessian_calls();
    convergence.iterations = result.iterations;
    convergence.evaluations = result.evaluations;
    if (nlp->status() == Ipopt::LOCAL_INFEASIBILITY) {
        convergence.stop = core::Stop::infeasible;
    } else {
        convergence.stop = core::classify(
            convergence.residual, *convergence.criterion,
            convergence.floor, convergence.contraction,
            nlp->status() != Ipopt::SUCCESS
        );
    }
    result.message = core::describe(convergence);
    scope.finish(convergence);
    const auto offsets = phase_offsets(compiled);
    std::vector<std::vector<double>> compositions;
    std::vector<double> amounts;
    std::vector<double> volumes;
    std::size_t intensive_offset = phase_width(compiled);
    const double candidate_temperature = compiled.problem.temperature.free
        ? full[intensive_offset] : compiled.problem.temperature.value;
    bool topology_ok = true;
    std::string topology_failure;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> composition(model.parameters().component_count(), 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const double value = std::exp(full[offsets[phase] + i]);
            composition[spec.support[i]] = value;
            total += value;
        }
        if (!(total > 0.0)) topology_ok = false;
        for (double& value : composition) value /= std::max(total, std::numeric_limits<double>::min());
        compositions.push_back(composition);
        std::size_t cursor = offsets[phase] + spec.support.size();
        const double amount = spec.amount.kind == AmountKind::Free
            ? full[cursor++] : specified_amount(spec, compiled.n_feed);
        amounts.push_back(amount);
        volumes.push_back(model.packing_fraction_coefficient(
            candidate_temperature, composition) / std::exp(full[cursor]));
        const double deficit = michelsen_deficit(
            compiled, compositions[phase], phase, result.residuals);
        scope.probe("michelsen_deficit", deficit, static_cast<std::uint32_t>(phase));
        if (spec.amount.kind == AmountKind::Free
            && amount / compiled.feed_total <= policy::STABILITY_PHASE_AMOUNT_FLOOR.value) {
            topology_ok = false;
            topology_failure = "phase " + spec.name + " vanished; Michelsen deficit "
                + std::to_string(deficit);
        }
    }
    if (compiled.problem.phases.size() > 1) {
        for (std::size_t left = 0; left < compositions.size(); ++left) {
            for (std::size_t right = left + 1; right < compositions.size(); ++right) {
                const auto left_eta = offsets[left] + compiled.problem.phases[left].support.size()
                    + (compiled.problem.phases[left].amount.kind == AmountKind::Free ? 1 : 0);
                const auto right_eta = offsets[right] + compiled.problem.phases[right].support.size()
                    + (compiled.problem.phases[right].amount.kind == AmountKind::Free ? 1 : 0);
                const double distance = std::max(
                    std::abs(full[left_eta] - full[right_eta]),
                    largest_difference(compositions[left], compositions[right]));
                scope.probe("phase_distance", distance);
                if (distance <= policy::PHASE_COINCIDENCE_DISTANCE.value) {
                    topology_ok = false;
                    topology_failure = "phases_coincide: phase distance between "
                        + compiled.problem.phases[left].name + " and "
                        + compiled.problem.phases[right].name + " is too small: "
                        + std::to_string(distance);
                }
            }
        }
    }
    const bool solver_ok = nlp->status() == Ipopt::SUCCESS
        || (nlp->status() == Ipopt::STOP_AT_TINY_STEP
            && result.convergence.stop == core::Stop::converged);
    const bool accepted = solver_ok
        && (result.convergence.stop == core::Stop::converged
            || result.convergence.stop == core::Stop::converged_to_floor)
        && topology_ok;
    if (!accepted) {
        result.success = false;
        if (!topology_ok) result.message = topology_failure;
        clear_state(result);
        return result;
    }
    const double temperature = compiled.problem.temperature.free
        ? full[intensive_offset++] : compiled.problem.temperature.value;
    const double pressure = compiled.problem.pressure.free
        ? full[intensive_offset++] : compiled.problem.pressure.value;
    result.success = true;
    result.temperature = {temperature};
    result.pressure = {pressure};
    result.amounts = amounts;
    result.molar_densities.reserve(volumes.size());
    result.packing_fractions.reserve(volumes.size());
    for (std::size_t phase = 0; phase < volumes.size(); ++phase) {
        result.molar_densities.push_back(1.0 / volumes[phase]);
        result.packing_fractions.push_back(
            model.packing_fraction_coefficient(temperature, compositions[phase]) / volumes[phase]);
        result.mole_fractions.insert(result.mole_fractions.end(),
                                     compositions[phase].begin(), compositions[phase].end());
    }
    result.multipliers.assign(full.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size), full.end());
    result.constraint_multipliers = nlp->constraint_multipliers();
    result.coordinates = full;
    std::vector<double> primal_point(
        full.begin(), full.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size));
    result.reduced_hessian_minimum = reduced_hessian_minimum_at(
        compiled, primal_point, result.constraint_multipliers);
    if (compiled.pinned_count == 0 && result.reduced_hessian_minimum <= 0.0) {
        result.success = false;
        result.message = "reduced_hessian_stability is not positive definite";
        clear_state(result);
        return result;
    }
    result.rows = compiled.rows;
    require_solution_shape(result, compiled.problem.phases.size(),
                           model.parameters().component_count());
    scope.probe("reduced_hessian_minimum", result.reduced_hessian_minimum);
    scope.probe("phase_count", static_cast<double>(result.phase_count));
    scope.probe("residual", result.convergence.residual);
    scope.probe("iterations", static_cast<double>(result.iterations));
    scope.probe("evaluations", static_cast<double>(result.evaluations));
    scope.probe("rows", static_cast<double>(result.rows.size()));
    scope.probe("multipliers", static_cast<double>(result.multipliers.size()));
    scope.probe("amount_total", std::accumulate(result.amounts.begin(), result.amounts.end(), 0.0));
    scope.probe("topology_ok", topology_ok ? 1.0 : 0.0);
    scope.probe("feasible", compiled.feasible ? 1.0 : 0.0);
    scope.probe("square", compiled.is_square() ? 1.0 : 0.0);
    return result;
}

}  // namespace epcsaft::equilibrium
