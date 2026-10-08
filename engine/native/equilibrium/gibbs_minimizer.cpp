#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "equilibrium/ipopt.hpp"
#include "contributions/helmholtz.hpp"

#include "core/trace.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <IpIpoptApplication.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <string>

namespace epcsaft::equilibrium {

// Infinity-norm forward-error estimate in the scaled physical correction chart.
// Eigen estimates the 1-norm condition of the transpose, hence the infinity norm of J D.
double polish_forward_bound(const std::vector<double>& jacobian, const std::vector<double>& rows,
    const std::vector<double>& delta, const std::vector<double>& scale) {
    const auto n = static_cast<Eigen::Index>(delta.size());
    const double unavailable = std::numeric_limits<double>::infinity();
    if (n == 0 || rows.size() != delta.size() || scale.size() != delta.size()
        || jacobian.size() != delta.size() * delta.size()) return unavailable;
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> J(jacobian.data(), n, n);
    const Eigen::Map<const Eigen::VectorXd> d(delta.data(), n), r(rows.data(), n), D(scale.data(), n);
    if (!J.allFinite() || !d.allFinite() || !r.allFinite() || !D.allFinite() || !(D.minCoeff() > 0.0)) return unavailable;
    const Eigen::MatrixXd scaled = J * D.asDiagonal();
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(scaled.transpose());
    const double denominator = lu.rcond() * scaled.cwiseAbs().rowwise().sum().maxCoeff();
    if (lu.rank() != n || !(denominator > 0.0) || !std::isfinite(denominator)) return unavailable;
    return d.cwiseQuotient(D).lpNorm<Eigen::Infinity>() + (J * d + r).lpNorm<Eigen::Infinity>() / denominator;
}
namespace {

// EqID reduced_hessian_stability on the balance nullspace. `minimum` is the smallest eigenvalue in the solved
// chart. The decision reads the congruent equilibration D = diag(S)^-1/2, S = |objective| + |weighted constraint|
// curvature, whose inertia is the same, against the rounding band eps (n ||DSD|| + 2 kappa(AD) ||DHD||)
// (ALG-EQUILIBRIUM, #159): trace log-chart directions have curvature far below eps ||H||.
struct Curvature final { double minimum = 0.0, equilibrated = 1.0, band = 0.0; };

Curvature reduced_curvature_at(
    const CompiledSystem& compiled, PhasePoint prepared,
    const std::vector<double>& multipliers
) {
    if (compiled.feasibility_cell() || prepared.point.empty()) return {};
    prepared.point.resize(compiled.primal_size);
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto n = static_cast<Eigen::Index>(prepared.point.size());
    const auto jacobian = compiled.jacobian_at(prepared);
    Eigen::MatrixXd constraints(0, n);
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.rows[row].kind == RowKind::Chemical || compiled.rows[row].kind == RowKind::Pressure) continue;
        constraints.conservativeResize(constraints.rows() + 1, n);
        constraints.bottomRows(1) = Eigen::Map<const Dense>(jacobian.data() + row * prepared.point.size(), 1, n);
    }
    const Eigen::MatrixXd objective = Eigen::Map<const Dense>(compiled.hessian_at(prepared, {}).data(), n, n);
    const Eigen::MatrixXd weighted = Eigen::Map<const Dense>(
        compiled.constraint_hessian_at(prepared, multipliers).data(), n, n);
    const Eigen::MatrixXd hessian = objective + weighted;
    const auto reduce = [n](const Eigen::MatrixXd& h, const Eigen::MatrixXd& a) -> std::pair<double, double> {
        const auto svd = a.jacobiSvd(Eigen::ComputeFullV);
        const Eigen::Index rank = svd.rank();
        if (rank >= n) return {1.0, 1.0};
        const Eigen::MatrixXd z = svd.matrixV().rightCols(n - rank);
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(z.transpose() * h * z, Eigen::EigenvaluesOnly);
        return {eigen.info() == Eigen::Success ? eigen.eigenvalues()(0) : std::numeric_limits<double>::quiet_NaN(),
                rank > 0 ? svd.singularValues()(0) / svd.singularValues()(rank - 1) : 1.0};
    };
    const Eigen::MatrixXd magnitude = objective.cwiseAbs() + weighted.cwiseAbs();
    const Eigen::VectorXd d = magnitude.diagonal().unaryExpr([](double s) { return s > 0.0 ? 1.0 / std::sqrt(s) : 1.0; });
    const Eigen::MatrixXd scaled = d.asDiagonal() * hessian * d.asDiagonal();
    const auto [equilibrated, kappa] = reduce(scaled, constraints * d.asDiagonal());
    const double rounding = static_cast<double>(n) * (d.asDiagonal() * magnitude * d.asDiagonal()).norm();
    return {reduce(hessian, constraints).first, equilibrated,
            std::numeric_limits<double>::epsilon() * (rounding + 2.0 * kappa * scaled.norm())};
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

void require_solution_shape(
    const EquilibriumResult& result, std::size_t phases, std::size_t components
) {
    if (!all_numeric(result.residuals) || result.rows.empty()) {
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

std::vector<double> start_point(
    const CompiledSystem& compiled, const EquilibriumResult& start,
    std::vector<double> point
) {
    const std::size_t components = compiled.model->parameters().component_count();
    if (start.phase_count != compiled.problem.phases.size()) return point;
    const auto offsets = compiled.phase_offsets();
    for (std::size_t phase = 0; phase < start.phase_count; ++phase) {
        if (valid_previous_phase(start, phase, components)) {
            copy_previous_phase(compiled, start, phase, point, offsets);
        }
    }
    std::size_t intensive = compiled.primal_size - compiled.free_intensive_count;
    if (compiled.problem.temperature.free && !start.temperature.empty()) point[intensive++] = start.temperature[0];
    if (compiled.problem.pressure.free && !start.pressure.empty()) point[intensive++] = start.pressure[0];
    if (compiled.feasibility_cell() && start.multipliers.size() == compiled.rank + compiled.charge_row_count) {
        std::copy(start.multipliers.begin(), start.multipliers.end(), point.begin() + compiled.primal_size);
    }
    return point;
}

bool usable_start(const CompiledSystem& compiled, const EquilibriumResult& start) {
    if (!supports_split_seed(compiled) || !start.success
        || start.phase_count != compiled.problem.phases.size()
        || start.amounts.size() != start.phase_count) return false;
    const std::size_t components = compiled.model->parameters().component_count();
    for (std::size_t phase = 0; phase < start.phase_count; ++phase) {
        if (!valid_previous_phase(start, phase, components)
            || !(std::isfinite(start.amounts[phase]) && start.amounts[phase] > 0.0)) return false;
        std::vector<double> composition(components);
        for (std::size_t i = 0; i < components; ++i) {
            composition[i] = start.mole_fractions[phase * components + i];
            if (!(composition[i] > 0.0)) return false;
        }
        const double packing = compiled.model->packing_fraction_coefficient(
            start.temperature.empty() ? compiled.problem.temperature.value : start.temperature[0],
            composition) * start.molar_densities[phase];
        if (!admitted_packing_fraction(packing)) return false;
    }
    return true;
}

double active_constraint_residual(
    const CompiledSystem& compiled, const std::vector<double>& values
) {
    double result = 0.0;
    for (std::size_t row = 0; row < values.size(); ++row) {
        if (!compiled.feasibility_cell() && compiled.rows[row].kind == RowKind::Chemical) continue;
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
    const auto keep = [](double& slot, double value) {  // a nonfinite row stays infinite
        slot = std::isfinite(value) ? std::max(slot, std::abs(value)) : std::numeric_limits<double>::infinity();
    };
    for (std::size_t index = 0; index < std::min(rows.size(), values.size()); ++index) {
        const RowKind kind = rows[index].kind;
        keep(result.maximum, values[index]);
        if (kind == RowKind::Chemical || kind == RowKind::Pressure) keep(result.intensive, values[index]);
        else if (kind == RowKind::Moiety) keep(result.balance, values[index]);
        else if (kind == RowKind::Enthalpy || kind == RowKind::Entropy) keep(result.balance, values[index]);
        else if (kind == RowKind::Normalization) keep(result.normalization, values[index]);
        else if (kind == RowKind::Electro) keep(result.charge, values[index]);
    }
    return result;
}

}  // namespace

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

std::vector<double> add_multipliers(
    const CompiledSystem& compiled, const std::vector<double>& point
) {
    if (point.size() != compiled.primal_size) return point;
    const std::size_t width = compiled.rank + compiled.charge_row_count;
    std::vector<double> result = point;
    result.resize(point.size() + width, 0.0);
    if (width == 0) return result;
    // Recover the balance and charge multipliers by the KKT-weighted least squares of the
    // chemical rows: the free-amount minimizer's stationarity conditions carry the weights
    // N_p x_i / N_feed, so the weighted fit is the one whose classified residual is minimal.
    const auto prepared = compiled.prepare(result, true);
    const auto values = compiled.evaluate_at(prepared).values;
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
    const auto& z = compiled.model->parameters().charge_number;
    for (std::size_t row = 0; row < chemical.size(); ++row) {
        const auto& item = compiled.rows[chemical[row]];
        const double weight = prepared.amount[item.phase] * prepared.composition[item.phase][item.index]
            / compiled.feed_total;
        rhs(static_cast<Eigen::Index>(row)) = weight * values[chemical[row]];
        for (std::size_t column = 0; column < compiled.rank; ++column) {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                weight * compiled.B[column][item.index];
        }
        for (std::size_t column = 0; column < phi_phase.size(); ++column) {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(compiled.rank + column)) =
                weight * z[item.index] * (item.phase == phi_phase[column] ? 1.0 : 0.0);
        }
    }
    const Eigen::VectorXd solved = matrix.colPivHouseholderQr().solve(rhs);
    for (std::size_t column = 0; column < width; ++column) {
        result[point.size() + column] = solved(static_cast<Eigen::Index>(column));
    }
    return result;
}

GibbsMinimizer::GibbsMinimizer(
    const CompiledSystem& compiled, std::vector<double> start, core::Scope* scope, const DriverLimits* limits
)
    : compiled_(compiled), start_(std::move(start)), scope_(scope), limits_(limits) {}

bool GibbsMinimizer::get_nlp_info(
    Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style
) {
    n = static_cast<Ipopt::Index>(!compiled_.feasibility_cell()
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
    if (limits_ && limits_->expired()) return false;
    const double infinity = 1.0e19;
    std::fill(x_l, x_l + n, -infinity);
    std::fill(x_u, x_u + n, infinity);
    // Feasibility cells impose chemical equality, so their box is the representability limit. Other
    // cells cap the box at a fraction of each species' attainable inventory u_i/N_feed, the LP
    // max y_i s.t. B y = b/N_feed, y >= 0 on the supports (B has mixed signs, so no row ratio
    // bounds y_i): the phases are never forced to hold a trace species that does not exist (#172).
    const std::size_t species = compiled_.n_feed.size();
    const bool feasibility = compiled_.feasibility_cell();
    std::vector<double> box(species, feasibility ? std::numeric_limits<double>::min() : policy::EQUILIBRIUM_BOX_FLOOR.value);
    if (!feasibility) {
        if (limits_ && limits_->expired()) return false;
        for (std::size_t i = 0; i < species; ++i)
            box[i] = std::max(box[i], std::min(policy::EQUILIBRIUM_COMPOSITION_BOX.value,
                policy::EQUILIBRIUM_BOX_INVENTORY_FRACTION.value * compiled_.inventory_upper_bounds[i]));
    }
    if (limits_ && limits_->expired()) return false;
    std::size_t offset = 0;
    for (const auto& phase : compiled_.problem.phases) {
        for (const std::size_t i : phase.support) {
            x_l[offset] = std::log(box[i]);
            x_u[offset++] = 0.0;
        }
        if (phase.amount.kind == AmountKind::Free) x_l[offset++] = 0.0;
        if (!compiled_.log_packing_unbounded()) {
            x_l[offset] = std::log(policy::PACKING_FRACTION_FLOOR.value);
            x_u[offset] = std::log(policy::PACKING_FRACTION_CEILING.value);
        }
        ++offset;
    }
    offset += compiled_.mode_count;  // the critical mode u is unbounded
    if (compiled_.problem.temperature.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    if (compiled_.problem.pressure.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    std::fill(g_l, g_l + m, 0.0);
    std::fill(g_u, g_u + m, 0.0);
    return true;
}

bool GibbsMinimizer::get_starting_point(
    Ipopt::Index n, bool init_x, Ipopt::Number* x,
    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda
) {
    if (!init_x) return false;
    const std::vector<double> point = !compiled_.feasibility_cell()
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

std::vector<std::size_t> GibbsMinimizer::constraint_rows() const {
    std::vector<std::size_t> result;
    for (std::size_t row = 0; row < compiled_.rows.size(); ++row) {
        if (compiled_.feasibility_cell()
            || (compiled_.rows[row].kind != RowKind::Chemical
                && compiled_.rows[row].kind != RowKind::Pressure)) {
            result.push_back(row);
        }
    }
    return result;
}

const PhasePoint& GibbsMinimizer::prepared(Ipopt::Index n, const Ipopt::Number* x) {
    if (limits_ && limits_->expired()) throw std::runtime_error("elapsed work budget exhausted");
    std::vector<double> candidate(x, x + n);
    if (!prepared_.has_value() || prepared_->point != candidate) {
        rows_.reset();
        jacobian_.reset();
        gradient_.reset();
        objective_hessian_.reset();
        constraint_hessian_.reset();
        try {
            prepared_ = compiled_.prepare(candidate, true);
        } catch (const std::exception& error) {
            prepared_.reset();
            preparation_error_ = error.what();
            throw;
        }
        ++work_.prepared_points;
    }
    if (limits_ && limits_->expired()) throw std::runtime_error("elapsed work budget exhausted");
    return *prepared_;
}

const PointEvaluation& GibbsMinimizer::rows_at(Ipopt::Index n, const Ipopt::Number* x) {
    const auto& point = prepared(n, x);
    if (!rows_.has_value()) {
        rows_ = compiled_.evaluate_at(point);
        ++work_.row_evaluations;
    }
    return *rows_;
}

bool GibbsMinimizer::eval_f(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& value
) {
    ++evaluations_;
    try {
        value = !compiled_.feasibility_cell() ? rows_at(n, x).objective : 0.0;
        return std::isfinite(value);
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_grad_f(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f
) {
    try {
        const auto& point = prepared(n, x);
        if (!gradient_.has_value()) {
            gradient_ = compiled_.objective_gradient_at(point);
            ++work_.gradients;
        }
        std::copy(gradient_->begin(), gradient_->end(), grad_f);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_g(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Index m, Ipopt::Number* g
) {
    try {
        const auto& values = rows_at(n, x).values;
        const auto rows = constraint_rows();
        for (Ipopt::Index row = 0; row < m; ++row) g[row] = values[rows[row]];
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_jac_g(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
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
    try {
        const auto& point = prepared(n, x);
        if (!jacobian_.has_value()) {
            jacobian_ = compiled_.jacobian_at(point);
            ++work_.jacobians;
        }
        Ipopt::Index entry = 0;
        for (const std::size_t row : rows) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                values[entry++] = (*jacobian_)[row * static_cast<std::size_t>(n) + static_cast<std::size_t>(column)];
            }
        }
        return entry == nele_jac;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_h(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
    bool /*new_lambda*/, Ipopt::Index nele_hess, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
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
    try {
        const auto& point = prepared(n, x);
        // Objective curvature is multiplier independent; the weighted constraint curvature is
        // recombined whenever the optimizer multipliers change at the same point.
        if (!objective_hessian_.has_value()) {
            objective_hessian_ = compiled_.hessian_at(point, {});
            ++work_.objective_hessians;
        }
        const std::vector<double> multipliers(lambda, lambda + m);
        if (!constraint_hessian_.has_value() || constraint_hessian_multipliers_ != multipliers) {
            constraint_hessian_ = compiled_.constraint_hessian_at(point, multipliers);
            constraint_hessian_multipliers_ = multipliers;
            ++work_.constraint_hessians;
        }
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < n; ++row) {
            for (Ipopt::Index column = 0; column <= row; ++column) {
                const std::size_t index = static_cast<std::size_t>(row) * n
                    + static_cast<std::size_t>(column);
                values[entry++] = obj_factor * (*objective_hessian_)[index]
                    + (*constraint_hessian_)[index];
            }
        }
        return entry == nele_hess;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::intermediate_callback(
    Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, Ipopt::Number /*obj_value*/,
    Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
    Ipopt::Number d_norm, Ipopt::Number regularization_size,
    Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
    const Ipopt::IpoptData* /*ip_data*/, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
) {
    iterations_ = static_cast<std::size_t>(iter);
    if (limits_) limits_->consumed = std::max(limits_->consumed, consumed_before_ + iterations_);
    regularization_size_ = regularization_size;
    if (scope_ != nullptr) scope_->probe("ipopt_iteration",
        {inf_pr, inf_du, mu, d_norm, regularization_size,
                            alpha_pr, alpha_du, static_cast<double>(ls_trials)},
        static_cast<std::uint32_t>(iter));
    return !limits_ || !limits_->expired();  // max_iter owns iteration termination
}

void GibbsMinimizer::finalize_solution(
    Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
    const Ipopt::Number* /*z_L*/, const Ipopt::Number* /*z_U*/, Ipopt::Index m,
    const Ipopt::Number* /*g*/, const Ipopt::Number* lambda,
    Ipopt::Number /*objective*/,
    const Ipopt::IpoptData* /*ip_data*/, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
) {
    status_ = status;
    solution_.assign(x, x + n);
    constraint_multipliers_.assign(lambda, lambda + m);
}

namespace {

void row_criterion(core::Convergence& record, double tolerance) {
    if (tolerance == policy::EQUILIBRIUM_ROW_RESIDUAL.value) record.criterion = &policy::EQUILIBRIUM_ROW_RESIDUAL;
    else {
        record.criterion_owner = std::make_shared<const policy::Tolerance>(policy::Tolerance{
            "EQUILIBRIUM_REQUESTED_ROW_RESIDUAL", tolerance,
            policy::Category::AlgorithmTolerance, "caller-supplied equilibrium row tolerance"});
        record.criterion = record.criterion_owner.get();
    }
}

void reject_elapsed(EquilibriumResult& result, core::Scope& scope) {
    auto& record = result.convergence;
    if (!record.criterion) row_criterion(record, result.declaration->row_tolerance);
    record.stop = core::classify(record.residual, *record.criterion, record.floor,
                                record.contraction, true, record.amplification);
    result.success = false; clear_state(result); scope.finish(record);
}

}  // namespace

EquilibriumResult structural_result(const CompiledSystem& compiled, core::Scope& scope) {
    EquilibriumResult result;
    result.model_identity = compiled.model;
    result.declaration = std::make_shared<const EquilibriumProblem>(compiled.problem);
    result.rows = compiled.rows;
    result.reference_diagnostics = compiled.reference_diagnostics;
    result.phase_count = compiled.problem.phases.size();
    result.unknown_count = compiled.unknown_count();
    result.jacobian_rank = compiled.jacobian_rank;
    result.generic_nullity = compiled.generic_nullity;
    if (!compiled.feasible) {
        // Structural rank loss at the generic point: no start, no optimizer, no candidate.
        result.solver_status = "NOT_RUN";
        result.message = "rank-deficient declared equations; nonlinear solve not run";
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        result.convergence.stop = core::Stop::infeasible;
        scope.probe("generic_nullity", static_cast<double>(compiled.generic_nullity));
        scope.finish(result.convergence);
    }
    return result;
}

// The one Ipopt lifecycle shared by the public seed route and the private coordinate route.
EquilibriumResult run_driver(
    const HelmholtzModel& model, const CompiledSystem& compiled,
    std::vector<double> initial, EquilibriumResult result, core::Scope& scope,
    std::vector<double>* attempt, EquilibriumResult* audit, const DriverLimits* limits
) {
    const double row_tolerance = compiled.problem.row_tolerance;
    // The optimizer works one policy decade below the caller criterion, which alone classifies.
    // A feasibility cell has a constant objective and is Newton-polished below, so it stops there.
    const double inner_tolerance = std::min(row_tolerance, policy::EQUILIBRIUM_ROW_RESIDUAL.value)
        * (compiled.feasibility_cell() ? 1.0 : policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value);
    const auto expired = [&] { return limits && limits->expired(); };
    if (expired()) { result.message = "elapsed work budget exhausted before driver";
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        reject_elapsed(result, scope); return result; }
    std::optional<EquilibriumResult> limit_audit;
    if (limits && !audit) audit = &limit_audit.emplace(result);
    auto fallback = limits ? initial : std::vector<double>{};
    Ipopt::SmartPtr<GibbsMinimizer> nlp = new GibbsMinimizer(compiled, std::move(initial), &scope, limits);
    if (audit) audit->convergence.residual = std::numeric_limits<double>::quiet_NaN();
    const auto capture_driver = [&] {
        if (!audit) return;
        audit->solver_status = solver_status_name(nlp->status());
        audit->convergence.iterations = audit->iterations = nlp->iterations();
        audit->convergence.evaluations = audit->evaluations = nlp->evaluations();
        audit->hessian_calls = nlp->hessian_calls();
        if (nlp->last_rows()) {
            audit->residuals = nlp->last_rows()->values;
            audit->classified_residuals = nlp->last_classified();
            if (nlp->last_rows()->floor > 0.0) audit->convergence.floor = nlp->last_rows()->floor;
        }
        audit->convergence.residual = audit->classified_residuals.empty()
            ? std::numeric_limits<double>::quiet_NaN() : maximum_abs(audit->classified_residuals);
    };
    const auto [initialized, returned] = [&] {
        try { return run_ipopt(nlp, inner_tolerance,
            static_cast<int>(limits ? std::min(limits->iterations - limits->consumed,
                static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)) : policy::EQUILIBRIUM_MAX_ITERATIONS.value),
            limits ? std::max(0.0, std::chrono::duration<double>(limits->deadline - std::chrono::steady_clock::now()).count()) : 0.0,
            true, 1.0 / compiled.feed_total);
        } catch (...) { capture_driver(); throw; }
    }();
    capture_driver();
    const Ipopt::SolverReturn status = nlp->status();
    // ALG-TERMINATION: a failed initialization (status stays UNASSIGNED) or fatal return leaves no candidate.
    const bool fatal = initialized != Ipopt::Solve_Succeeded
        || (returned <= Ipopt::Invalid_Problem_Definition && returned != Ipopt::Invalid_Number_Detected);
    // The one accepted-outcome rule: only a non-error stop polishes, classifies its rows and publishes.
    const bool non_error = !expired() && !fatal && (status == Ipopt::SUCCESS
        || status == Ipopt::STOP_AT_ACCEPTABLE_POINT || status == Ipopt::STOP_AT_TINY_STEP
        || status == Ipopt::FEASIBLE_POINT_FOUND || status == Ipopt::RESTORATION_FAILURE
        || status == Ipopt::MAXITER_EXCEEDED || status == Ipopt::CPUTIME_EXCEEDED);
    result.solver_status = solver_status_name(status);  // both exact names tell a solver failure from infeasibility
    if (fatal) result.solver_status += (initialized != Ipopt::Solve_Succeeded ? "; Initialize returned "
        : "; OptimizeTNLP returned ") + application_return_name(returned);
    const std::string ipopt_failure = fatal ? "; Ipopt " + result.solver_status : "";
    if (audit) audit->solver_status = result.solver_status;
    const auto limited_result = [&] {
        auto failed = *audit;
        failed.message = "elapsed work budget exhausted";
        reject_elapsed(failed, scope); return failed;
    };
    if (expired()) return limited_result();
    std::vector<double> point = nlp->solution();
    if (point.empty()) point = limits ? std::move(fallback) : boundary_start(compiled);
    std::vector<double> full = !compiled.feasibility_cell() ? add_multipliers(compiled, point) : point;
    // The optimizer's amount-weighted rows underresolve intensive trace defects.
    // Exact Newton steps contract raw rows at any positive supported composition.
    // Curvature may raise classified rows temporarily; return the best classified
    // point visited, including the optimizer start, rather than the last step.
    std::size_t polish_evaluations = 0;
    std::string rejected_preparation = nlp->preparation_error().value_or("");
    // A budget stop polishes too: an unattainable request can exhaust Ipopt after passing a
    // better point, and the polish never returns a worse classified point than its start.
    if (non_error && !expired()) {
        auto at = compiled.prepare(full, true);
        ++polish_evaluations;
        if (audit) audit->convergence.evaluations = ++audit->evaluations;
        if (expired()) return limited_result();
        auto rows_before = compiled.evaluate_at(at).values;
        std::vector<double> best = full;
        double best_classified = maximum_abs(compiled.classified_residuals(at, rows_before));
        // Past the step budget only raw rows above the caller criterion continue, under a separate polish cap,
        // and only while the raw contraction ratio keeps falling, as it does in Newton's quadratic region.
        double best_raw = scope.diagnostics_on() ? maximum_abs(rows_before) : 0.0;
        double ratio = 1.0;
        const bool detailed = scope.diagnostics_on();
        std::size_t final_step = 0;
        bool stopped = false;
        const auto stop = [&](const char* reason) {
            if (stopped) return;
            stopped = true;
            if (detailed) scope.observe([&] { scope.probe("polish_stop",
                {maximum_abs(rows_before), maximum_abs(compiled.classified_residuals(at, rows_before)),
                 static_cast<double>(final_step)}, static_cast<std::uint32_t>(final_step), reason); });
        };
        for (std::size_t step = 0; step < static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value); ++step) {
            final_step = step;
            if (expired()) { stop("elapsed budget"); break; }
            const double defect = maximum_abs(rows_before);
            const bool extension = step >= static_cast<std::size_t>(policy::EQUILIBRIUM_POLISH_STEPS.value);
            if (extension && defect <= row_tolerance) { stop("raw criterion"); break; }
            const double classified_before = maximum_abs(compiled.classified_residuals(at, rows_before));
            if (!(defect > 0.0) || !std::isfinite(defect)) { stop("nonfinite/zero defect"); break; }
            const auto jacobian = compiled.jacobian_at(at);
            if (expired()) { stop("elapsed budget"); return limited_result(); }
            const auto n = static_cast<Eigen::Index>(full.size());
            if (jacobian.size() != full.size() * full.size() || rows_before.size() != full.size()) { stop("Jacobian shape"); break; }
            const auto lu = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor>>(jacobian.data(), n, n).fullPivLu();
            const Eigen::VectorXd delta = lu.solve(-Eigen::Map<const Eigen::VectorXd>(rows_before.data(), n));
            Eigen::VectorXd scale = Eigen::VectorXd::Ones(n);
            bool positive = true;
            for (std::size_t p = 0; p < at.shape.size(); ++p) if (at.shape[p].amount) {
                const auto i = static_cast<Eigen::Index>(*at.shape[p].amount);
                positive = positive && at.amount[p] > 0.0;
                if (at.amount[p] > 0.0) scale(i) = at.amount[p];
            }
            auto intensive = static_cast<Eigen::Index>(compiled.primal_size - compiled.free_intensive_count);
            if (compiled.problem.temperature.free) scale(intensive++) = at.temperature;
            if (compiled.problem.pressure.free) scale(intensive) = at.pressure;
            Eigen::Index controlling = 0;
            const double scaled_correction = delta.allFinite() ? delta.cwiseQuotient(scale).cwiseAbs().maxCoeff(&controlling)
                : std::numeric_limits<double>::infinity();
            double bound = std::numeric_limits<double>::infinity();
            const bool size_admitted = positive && full == best && classified_before <= row_tolerance
                && scaled_correction <= policy::EQUILIBRIUM_POLISH_CORRECTION.value;
            if (size_admitted && lu.rank() == n)
                bound = polish_forward_bound(jacobian, rows_before, {delta.data(), delta.data() + n},
                    {scale.data(), scale.data() + n});
            const bool correction_stop = size_admitted && bound <= policy::EQUILIBRIUM_POLISH_CORRECTION.value;
            if (detailed) scope.observe([&] {
                double feed_relative = 0.0;
                for (const auto& shape : at.shape) if (shape.amount)
                    feed_relative = std::max(feed_relative, std::abs(delta(*shape.amount)) / compiled.feed_total);
                scope.probe("polish_correction", {scaled_correction, static_cast<double>(controlling), feed_relative, bound,
                classified_before, full == best ? 1.0 : 0.0, positive ? 1.0 : 0.0,
                size_admitted ? static_cast<double>(lu.rank()) : -1.0, correction_stop ? 1.0 : 0.0},
                static_cast<std::uint32_t>(step)); });
            double correction = std::numeric_limits<double>::quiet_NaN(), control = correction;
            if (detailed) scope.observe([&] {
                Eigen::Index index = 0;
                correction = delta.cwiseAbs().maxCoeff(&index); control = static_cast<double>(index);
                scope.probe("polish_direction", {correction, control, defect, classified_before,
                    extension ? 1.0 : 0.0, ratio}, step, delta.allFinite() ? "" : "nonfinite direction");
                double phase = std::numeric_limits<double>::quiet_NaN(), component = phase;
                const char* quantity = "multiplier in RT units";
                const auto offsets = compiled.phase_offsets();
                for (std::size_t p = 0; p < offsets.size(); ++p) {
                    const auto& spec = compiled.problem.phases[p];
                    const auto offset = offsets[p], eta = offset + spec.support.size() + (spec.is_free() ? 1 : 0);
                    if (index < static_cast<Eigen::Index>(offset) || index > static_cast<Eigen::Index>(eta)) continue;
                    phase = static_cast<double>(p);
                    if (index < static_cast<Eigen::Index>(offset + spec.support.size())) {
                        component = static_cast<double>(spec.support[index - offset]); quantity = "log mole fraction";
                    } else quantity = index == static_cast<Eigen::Index>(eta) ? "log packing fraction" : "phase amount in mol";
                }
                auto intensive = compiled.primal_size - compiled.free_intensive_count;
                if (compiled.problem.temperature.free && index == static_cast<Eigen::Index>(intensive++)) quantity = "T in K";
                if (compiled.problem.pressure.free && index == static_cast<Eigen::Index>(intensive)) quantity = "P in Pa";
                scope.probe("polish_coordinate", {control, phase, component, full[index]}, step, quantity);
            });
            if (correction_stop) { stop("scaled correction"); break; }

            bool accepted = false;
            // Bounded damping retains raw contraction; preparation admits only positive finite x.
            // Best-point selection below protects the returned classified residual.
            std::uint32_t trial = 0;
            for (const double alpha : {1.0, 0.5, 0.25}) {
                const auto trial_index = trial++;
                if (expired()) { stop("elapsed budget"); break; }
                std::vector<double> candidate = full;
                for (Eigen::Index i = 0; i < n; ++i)
                    candidate[static_cast<std::size_t>(i)] += alpha * delta(i);
                PhasePoint at_candidate;
                std::vector<double> rows_after;
                try {
                    at_candidate = compiled.prepare(candidate, true);
                    ++polish_evaluations;
                    if (audit) audit->convergence.evaluations = ++audit->evaluations;
                    if (expired()) { stop("elapsed budget"); return limited_result(); }
                    rows_after = compiled.evaluate_at(at_candidate).values;
                } catch (const std::exception& error) {
                    rejected_preparation = error.what();
                    if (detailed) scope.probe("polish_trial", {correction, control, defect, classified_before,
                        NAN, NAN, alpha, 0.0, extension ? 1.0 : 0.0, ratio,
                        static_cast<double>(trial_index)}, step, error.what());
                    continue;
                }
                const double after = maximum_abs(rows_after);
                const double classified_after = maximum_abs(compiled.classified_residuals(at_candidate, rows_after));
                const bool admit = after < (extension ? ratio : 1.0) * defect;
                if (detailed) scope.probe("polish_trial", {correction, control, defect, classified_before,
                    after, classified_after, alpha, admit ? 1.0 : 0.0, extension ? 1.0 : 0.0, ratio,
                    static_cast<double>(trial_index)}, step,
                    admit ? "raw_contraction" : extension ? "extension ratio" : "raw noncontraction");
                else scope.probe("polish_step", {defect, after, classified_before, classified_after, alpha},
                    static_cast<std::uint32_t>(step));
                if (!admit) continue;
                ratio = after / defect;
                full = std::move(candidate);
                at = std::move(at_candidate);
                rows_before = std::move(rows_after);
                if (classified_after <= best_classified) {
                    best = full;
                    best_classified = classified_after;
                    if (detailed) best_raw = after;
                }
                accepted = true;
                break;
            }
            if (!accepted) { stop("no accepted trial"); break; }
            if (step + 1 == static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)) stop("iteration cap");
        }
        if (detailed && best != full) scope.probe("polish_stop", {best_raw, best_classified,
            static_cast<double>(final_step)}, final_step, "best classified point returned");
        full = std::move(best);
    }
    if (expired()) return limited_result();
    PhasePoint prepared;
    try {
        prepared = compiled.prepare(full);
    } catch (const std::invalid_argument& error) {
        if (std::string(error.what()).rfind("reference_unavailable:", 0) != 0) throw;
        result.success = false;
        result.message = error.what() + ipopt_failure;
        result.convergence.stop = core::Stop::left_domain;
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        result.iterations = result.convergence.iterations = nlp->iterations();
        result.evaluations = result.convergence.evaluations = nlp->evaluations() + polish_evaluations;
        clear_state(result);
        scope.finish(result.convergence);
        return result;
    }
    if (expired()) return limited_result();
    if (attempt && !fatal && !nlp->solution().empty()) *attempt = full;
    const auto evaluation = compiled.evaluate_at(prepared);
    if (expired()) return limited_result();
    result.residuals = evaluation.values;
    result.classified_residuals = compiled.classified_residuals(prepared, evaluation.values);
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
    scope.probe("active_residual", active_constraint_residual(compiled, evaluation.values));
    scope.probe("objective", result.objective);
    scope.probe("hessian_calls", static_cast<double>(nlp->hessian_calls()));
    scope.probe("prepared_points", static_cast<double>(nlp->work().prepared_points));
    core::Convergence& convergence = result.convergence;
    convergence = {};
    convergence.residual = maximum_abs(result.classified_residuals);
    row_criterion(convergence, row_tolerance);
    if (evaluation.floor > 0.0) convergence.floor = evaluation.floor;
    convergence.amplification = nlp->regularization_size();
    result.iterations = convergence.iterations = nlp->iterations();
    result.evaluations = convergence.evaluations = nlp->evaluations() + polish_evaluations;
    result.hessian_calls = nlp->hessian_calls();
    switch (status) {
    case Ipopt::INVALID_NUMBER_DETECTED: convergence.stop = core::Stop::nonfinite; break;
    case Ipopt::DIVERGING_ITERATES: convergence.stop = core::Stop::left_domain; break;
    case Ipopt::ERROR_IN_STEP_COMPUTATION: convergence.stop = core::Stop::singular; break;
    case Ipopt::LOCAL_INFEASIBILITY:
    case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM: convergence.stop = core::Stop::infeasible; break;
    default:
        // No candidate after a fatal return: infeasible like an unrun structural solve, whatever the start's rows.
        convergence.stop = !non_error ? core::Stop::infeasible : core::classify(convergence.residual,
            *convergence.criterion, convergence.floor, convergence.contraction,
            status == Ipopt::MAXITER_EXCEEDED || status == Ipopt::CPUTIME_EXCEEDED, convergence.amplification);
    }
    if (std::isfinite(convergence.residual)) {
        result.requested_tolerance_met = convergence.residual <= row_tolerance;
    }
    result.message = core::describe(convergence);
    scope.finish(convergence);
    std::vector<double> volumes;
    std::size_t intensive_offset = compiled.primal_size - compiled.free_intensive_count;
    const double candidate_temperature = compiled.problem.temperature.free
        ? full[intensive_offset] : compiled.problem.temperature.value;
    bool topology_ok = true;
    std::string topology_failure;
    // Preparation above admitted every composition with a positive total.
    const auto topology = check_topology(compiled, full);
    const auto& compositions = topology.compositions;
    const auto& amounts = topology.amounts;
    result.topology_event = topology.event;
    result.topology_phase = topology.phase;
    result.topology_pair = topology.pair;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        volumes.push_back(model.packing_fraction_coefficient(
            candidate_temperature, compositions[phase]) / std::exp(topology.log_packing[phase]));
        const double deficit = michelsen_deficit(
            compiled, compositions[phase], phase, result.residuals);
        scope.probe("michelsen_deficit", deficit, static_cast<std::uint32_t>(phase));
        if (topology.event == TopologyEvent::PhaseDisappeared && topology.phase == phase) {
            topology_ok = false;
            topology_failure = "phase " + spec.name + " vanished; Michelsen deficit "
                + std::to_string(deficit);
        }
    }
    for (const double distance : topology.distances) scope.probe("phase_distance", distance);
    if (topology.event == TopologyEvent::PhasesCoalesced) {
        topology_ok = false;
        topology_failure = "phases_coincide: phase distance between "
            + compiled.problem.phases[topology.pair->first].name + " and "
            + compiled.problem.phases[topology.pair->second].name + " is too small: "
            + std::to_string(topology.distance);
    }
    // The recomputed rows classify: a criterion-meeting point is accepted after any non-error
    // stop (error stops never classify as converged); a floor-limited point needs success.
    const bool accepted = !expired() && non_error && topology_ok && (result.convergence.stop == core::Stop::converged
        || (result.convergence.stop == core::Stop::converged_to_floor && status == Ipopt::SUCCESS));
    if (!accepted) {
        result.success = false;
        if (!topology_ok) result.message = topology_failure;
        result.message += ipopt_failure;
        if (!rejected_preparation.empty()) result.message += "; trial preparation failed: " + rejected_preparation;
        clear_state(result);
        return result;
    }
    if (result.requested_tolerance_met.has_value() && !*result.requested_tolerance_met) {
        result.message += "; requested row tolerance " + std::to_string(row_tolerance)
            + " not met by an available floor-limited state";
    }
    const double temperature = compiled.problem.temperature.free
        ? full[intensive_offset++] : compiled.problem.temperature.value;
    const double pressure = compiled.problem.pressure.free
        ? full[intensive_offset++] : compiled.problem.pressure.value;
    result.success = true;
    result.temperature = {temperature};
    result.pressure = {pressure};
    if (!prepared.standard.reference_diagnostics.empty()) {
        result.reference_diagnostics = prepared.standard.reference_diagnostics;
    }
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
    // Extensive-objective constraint multipliers at the returned (polished) point, from the row
    // multipliers: a material row (B n - b)/N_feed carries -N_feed lambda_j, a charge row sum z x
    // carries -n_p phi_p, and a free-amount phase's normalization row is zero because the
    // objective and material rows obey joint extensive homogeneity in amount, composition
    // and volume, together with stationarity in the free amount. The raw
    // optimizer values stay a trace probe.
    scope.probe("optimizer_constraint_multipliers", nlp->constraint_multipliers());
    if (!compiled.feasibility_cell()) {
        result.constraint_multipliers.clear();
        std::size_t electro = 0;
        for (const auto& row : compiled.rows) {
            if (row.kind == RowKind::Moiety) {
                result.constraint_multipliers.push_back(-compiled.feed_total * full[compiled.primal_size + row.index]);
            } else if (row.kind == RowKind::Electro) {
                result.constraint_multipliers.push_back(
                    -prepared.amount[row.phase] * full[compiled.primal_size + compiled.rank + electro++]);
            } else if (row.kind != RowKind::Chemical && row.kind != RowKind::Pressure) {
                result.constraint_multipliers.push_back(0.0);
            }
        }
    } else {
        result.constraint_multipliers = nlp->constraint_multipliers();
    }
    result.coordinates = full;
    if (expired()) return limited_result();
    Curvature curvature = reduced_curvature_at(compiled, prepared, result.constraint_multipliers);
    if (compiled.has_caloric_condition()) {
        const bool pure = model.parameters().component_count() == 1 && compiled.problem.phases.size() == 2;
        if (pure) {
            if (!(amounts[0] > 0.0 && amounts[1] > 0.0) || volumes[0] == volumes[1]) {
                result.success = false;
                result.message = "caloric coexistence requires an interior distinct split";
                clear_state(result);
                return result;
            }
            result.message += "; coexistence equality root; fixed-TP Gibbs split is analytically flat";
        } else {
            auto diagnostic = compiled.problem;
            diagnostic.enthalpy_j.reset();
            diagnostic.entropy_j_per_k.reset();
            diagnostic.temperature = FreeValue::fixed(temperature);
            auto tp = compile_problem(model, diagnostic);
            auto restricted = tp.generic_point();
            const auto layout = tp.prepare(restricted).shape;
            for (std::size_t p = 0; p < layout.size(); ++p) {
                for (std::size_t i = 0; i < layout[p].x.size(); ++i)
                    restricted[layout[p].x[i]] = full[prepared.shape[p].x[i]];
                restricted[layout[p].volume] = full[prepared.shape[p].volume];
                restricted[*layout[p].amount] = amounts[p];
            }
            restricted.resize(tp.primal_size);
            restricted = add_multipliers(tp, restricted);
            std::vector<double> weights;
            for (const auto& row : tp.rows) {
                if (row.kind == RowKind::Chemical || row.kind == RowKind::Pressure) continue;
                weights.push_back(row.kind == RowKind::Moiety ? -tp.feed_total * restricted[tp.primal_size + row.index] : 0.0);
            }
            curvature = reduced_curvature_at(tp, tp.prepare(restricted), weights);
        }
    }
    result.reduced_hessian_minimum = curvature.minimum;
    scope.probe("reduced_hessian_equilibrated", curvature.equilibrated);
    scope.probe("reduced_hessian_band", curvature.band);
    if ((!compiled.feasibility_cell() || compiled.has_caloric_condition()) && !(curvature.equilibrated > curvature.band)) {
        // Inside the band the inertia is unresolved at working precision: a typed singular outcome, not a noise sign.
        const bool singular = curvature.equilibrated >= -curvature.band;
        if (singular || std::isnan(curvature.equilibrated))
            result.convergence.stop = singular ? core::Stop::singular : core::Stop::nonfinite;
        result.success = false;
        result.message = std::string("reduced_hessian_stability ") + (singular ? "is indeterminate" : "is not positive definite")
            + ": equilibrated minimum " + core::optional_number(curvature.equilibrated)
            + ", rounding band " + core::optional_number(curvature.band);
        if (!rejected_preparation.empty()) {
            result.message += "; trial preparation failed: " + rejected_preparation;
        }
        clear_state(result);
        return result;
    }
    if (expired()) { result.message = "elapsed work budget exhausted";
        reject_elapsed(result, scope); return result; }
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
    return result;
}

namespace {

// Direct admission of a critical-equation candidate (EqID critical_equations) in place of the density-root
// rule, which refuses P_rho = 0: a regular augmented Jacobian and one lowest zero eigenvalue of M, read from the
// mode rows' u columns, within |lambda| <= ||M u||/||u|| + C eps ||M||, above no resolved negative one. No
// stability certificate; the record is retained either way.
BranchCheck admit_critical(const CompiledSystem& compiled, EquilibriumResult& result) {
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto n = static_cast<Eigen::Index>(result.coordinates.size()), c = static_cast<Eigen::Index>(compiled.mode_count);
    const auto jacobian = compiled.jacobian(result.coordinates);
    const Eigen::Map<const Dense> J(jacobian.data(), n, n);
    const Eigen::MatrixXd M = J.block(static_cast<Eigen::Index>(compiled.rows.size()) - c - 2,
                                      static_cast<Eigen::Index>(compiled.mode_offset()), c, c);
    const Eigen::Map<const Eigen::VectorXd> u(result.coordinates.data() + compiled.mode_offset(), c);
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(M, Eigen::EigenvaluesOnly);
    const auto lu = J.fullPivLu();
    const Eigen::VectorXd right = Eigen::VectorXd::Ones(n), solution = lu.solve(right);
    const auto singular = J.jacobiSvd().singularValues();
    CriticalRecord record;
    record.eigenvalues.assign(spectrum.eigenvalues().data(), spectrum.eigenvalues().data() + c);
    record.mode.assign(u.data(), u.data() + c);
    record.cubic = result.residuals.back();
    record.zero_band = (M * u).norm() / u.norm() + static_cast<double>(c) * std::numeric_limits<double>::epsilon() * M.norm();
    record.rank = static_cast<std::size_t>(lu.rank());
    record.condition_estimate = singular(0) / singular(n - 1);
    record.backward_error = (J * solution - right).cwiseAbs().maxCoeff()
        / (J.cwiseAbs().rowwise().sum().maxCoeff() * solution.cwiseAbs().maxCoeff() + 1.0);
    result.critical = record;
    BranchCheck outcome;
    if (record.rank != static_cast<std::size_t>(n) || spectrum.info() != Eigen::Success)
        outcome.message = "augmented critical Jacobian is singular (rank " + std::to_string(record.rank) + ")";
    else if (!(std::abs(record.eigenvalues[0]) <= record.zero_band))
        outcome.message = "lowest eigenvalue of M is not a resolved zero";
    else if (c > 1 && !(record.eigenvalues[1] > record.zero_band))
        outcome.message = "second eigenvalue of M is not resolved positive";
    else {
        outcome.accepted = true;
        result.message += "; critical-equation candidate, no stability certificate";
    }
    return outcome;
}

bool zero_amount(const PhaseSpec& phase) {
    return phase.amount.kind == AmountKind::Pinned && phase.amount.fraction == 0.0;
}

// A cold reactive boundary at fixed T and free P: the phases present (nonzero amount) alone, in reaction
// equilibrium at the pressure seed, are the solution's present phases up to their pressure dependence. Their
// compositions replace the caller's guesses before the ordinary boundary start, and an unguessed zero-amount
// phase starts from the first present phase on its support (#152); none if any of it fails. It spends `limits`.
std::optional<std::vector<double>> speciated_start(const CompiledSystem& compiled, const DriverLimits* limits, StageWork& stages) {
    const auto& problem = compiled.problem;
    EquilibriumProblem present = problem;
    present.phases.erase(std::remove_if(present.phases.begin(), present.phases.end(), zero_amount), present.phases.end());
    if (!compiled.has_reaction || !problem.pressure.free || problem.temperature.free || present.phases.empty()
        || present.phases.size() == problem.phases.size()) return std::nullopt;
    present.pressure = FreeValue::fixed(compiled.seed_pressure());
    try {
        const auto solved = solve_compiled(compile_problem(*compiled.model, present), std::nullopt, nullptr, limits);
        stages.add(solved);
        if (!solved.success) return std::nullopt;
        const std::size_t components = compiled.model->parameters().component_count();
        CompiledSystem seeded = compiled;
        std::size_t index = 0;
        for (auto& phase : seeded.problem.phases) {
            const bool guessed = phase.composition_guess.size() == phase.support.size();
            if (zero_amount(phase) && guessed) continue;
            const std::size_t source = zero_amount(phase) ? 0 : index++;
            phase.composition_guess.clear();
            for (const std::size_t i : phase.support) phase.composition_guess.push_back(solved.mole_fractions[source * components + i]);
        }
        return boundary_start(seeded, true, limits, &stages);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

}  // namespace

TopologyCheck check_topology(const CompiledSystem& compiled, const std::vector<double>& full) {
    TopologyCheck result;
    if (full.size() < compiled.primal_size) return result;
    const std::size_t components = compiled.model->parameters().component_count();
    const auto offsets = compiled.phase_offsets();
    auto& compositions = result.compositions;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> composition(components, 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            composition[spec.support[i]] = std::exp(full[offsets[phase] + i]);
            total += composition[spec.support[i]];
        }
        for (double& value : composition) value /= std::max(total, std::numeric_limits<double>::min());
        std::size_t cursor = offsets[phase] + spec.support.size();
        const double amount = spec.amount.kind == AmountKind::Free
            ? full[cursor++] : compiled.phase_amount(spec);
        result.amounts.push_back(amount);
        result.log_packing.push_back(full[cursor]);
        compositions.push_back(std::move(composition));
        if (spec.amount.kind == AmountKind::Free
            && amount / compiled.feed_total <= policy::PHASE_DISAPPEARANCE_AMOUNT.value
            && result.event == TopologyEvent::None) {
            result.event = TopologyEvent::PhaseDisappeared;
            result.phase = phase;
        }
    }
    for (std::size_t left = 0; left + 1 < compositions.size(); ++left) {
        for (std::size_t right = left + 1; right < compositions.size(); ++right) {
            const double distance = std::max(
                std::abs(result.log_packing[left] - result.log_packing[right]),
                largest_difference(compositions[left], compositions[right]));
            result.distances.push_back(distance);
            if (distance <= policy::PHASE_COINCIDENCE_DISTANCE.value
                && result.event == TopologyEvent::None) {
                result.event = TopologyEvent::PhasesCoalesced;
                result.pair = std::make_pair(left, right);
                result.distance = distance;
            }
        }
    }
    return result;
}

// Observe only frozen returned states; isolated work never enters search accounting.
static void observe_endpoint(const HelmholtzModel& model, const CompiledSystem* compiled,
    const EquilibriumResult& result, core::Scope& scope, std::uint32_t attempt = 0) noexcept {
    if (!scope.diagnostics_on() || !result.has_state()) return;
    scope.observe([&] {
        SensitivityInfo info;
        try {
            core::Capture isolated(false);
            info = inspect_sensitivity(compiled ? *compiled : compile_problem(model, *result.declaration), result.coordinates, {});
        } catch (const std::exception& error) {
            scope.probe("solved_endpoint_jacobian", {NAN, static_cast<double>(result.coordinates.size()), NAN, NAN}, attempt, error.what());
            return;
        }
        scope.probe("solved_endpoint_jacobian", {static_cast<double>(info.equations),
            static_cast<double>(result.coordinates.size()), static_cast<double>(info.rank), info.condition_estimate},
            attempt, "raw defining Jacobian");
    });
}

EquilibriumResult solve_equilibrium(const HelmholtzModel& model, const EquilibriumProblem& request,
                                    const std::optional<EquilibriumResult>& start) {
    core::Scope scope("solve_equilibrium");  // compilation belongs to the solve's span
    std::optional<CompiledSystem> compiled;
    if (!request.discovery) compiled = compile_problem(model, request);
    auto result = request.discovery ? discover_phases(model, request) : solve_compiled(*compiled, start, &scope);
    if (scope.diagnostics_on() && !scope.nested("solve_equilibrium"))
        observe_endpoint(model, compiled ? &*compiled : nullptr, result, scope);
    return result;
}

EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start,
                                 core::Scope* outer, const DriverLimits* limits) {
    std::optional<core::Scope> own;
    auto& scope = outer != nullptr ? *outer : own.emplace("solve_equilibrium");
    const auto& model = *compiled.model;
    const EquilibriumResult structural = structural_result(compiled, scope);
    if (!compiled.feasible) return structural;
    const bool use_start = start.has_value() && usable_start(compiled, *start);
    // A previous result that start_point copies replaces the start's composition; skip its speciation.
    const bool copied = start.has_value() && start->phase_count == compiled.problem.phases.size()
        && valid_previous_phase(*start, 0, model.parameters().component_count());
    StageWork stages;
    auto seeded = start.has_value() ? std::nullopt : speciated_start(compiled, limits, stages);
    EquilibriumResult result;
    // The speciated seed when there is one; the caller's start alone or after the seed's failure.
    for (int attempt = seeded ? 0 : 1; attempt < 2 && !result.success; ++attempt) {
        stages.add(result);  // a failed seeded attempt
        std::vector<double> initial;
        if (attempt == 0) {
            initial = std::move(*seeded);
        } else {
            initial = use_start ? compiled.generic_point() : boundary_start(compiled, !copied, limits, &stages);
            if (start.has_value()) initial = start_point(compiled, *start, std::move(initial));
        }
        try {
            result = run_driver(model, compiled, std::move(initial), structural, scope, nullptr, nullptr, limits);
        } catch (const std::exception&) {
            if (attempt == 0) continue;  // the seeded attempt's exception falls back to the declared start
            throw;
        }
        if (!result.success) continue;
        BranchCheck branches;
        try {
            branches = compiled.problem.critical ? admit_critical(compiled, result)
                : check_density_branches(model, compiled, result, {}, true, scope, limits);
        } catch (const std::exception& error) {
            branches.message = std::string("density branch admission failed: ") + error.what();
        }
        if (!branches.accepted) {
            result.success = false;
            result.message = branches.message;
            clear_state(result);
        }
    }
    stages.report(result);
    return result;
}

EquilibriumStartsResult solve_equilibrium_starts(const HelmholtzModel& model,
    const EquilibriumProblem& problem, const std::vector<EquilibriumStart>& starts,
    const EquilibriumWorkLimits& limits) {
    core::Scope scope("solve_equilibrium_starts");
    using Clock = std::chrono::steady_clock;
    const auto begun = Clock::now();
    const auto deadline_at = [](Clock::time_point now, const std::optional<double>& seconds) {
        if (seconds && (!std::isfinite(*seconds) || *seconds < 0.0))
            throw std::invalid_argument("work seconds must be finite and nonnegative");
        if (!seconds || *seconds >= std::chrono::duration<double>(Clock::time_point::max() - now).count())
            return Clock::time_point::max();
        return now + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(*seconds));
    };
    const auto deadline = deadline_at(begun, limits.seconds);
    deadline_at(begun, limits.attempt_seconds);  // validate even an unattempted set
    EquilibriumStartsResult out;
    out.results.resize(starts.size()); out.causes.assign(starts.size(), "unattempted");
    out.attempt_seconds.resize(starts.size(), 0.0);
    const auto exhausted = [&] { return Clock::now() >= deadline
        || (limits.iterations && out.iterations >= *limits.iterations); };
    std::optional<CompiledSystem> retained;
    if (!exhausted()) {
        auto& compiled = retained.emplace(compile_problem(model, problem));
        for (const auto& seed : starts) if (seed.previous) {
            const auto& old = *seed.previous;
            if (!old.declaration || old.model_identity != &model)
                throw std::invalid_argument("previous start model/parameters or declaration differs");
            if (const auto difference = declared_state_difference(*old.declaration, compiled.problem); !difference.empty())
                throw std::invalid_argument("previous start: " + difference);
        }
        for (std::size_t i = 0; i < starts.size(); ++i) {
            if (exhausted() || (limits.attempts && i >= *limits.attempts)) { out.budget_exhausted = true; break; }
            core::Scope scope("equilibrium_start_attempt");
            scope.probe("declared_start_index", static_cast<double>(i));
            const auto& seed = starts[i];
            const auto attempt_begun = Clock::now();
            DriverLimits work{std::min(limits.attempt_iterations.value_or(static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)),
                limits.iterations ? *limits.iterations - out.iterations : std::numeric_limits<std::size_t>::max()),
                std::min(deadline, deadline_at(Clock::now(), limits.attempt_seconds))};
            if (work.iterations == 0 || work.expired()) { out.causes[i] = "attempt budget exhausted"; continue; }
            std::optional<EquilibriumResult> audit;
            StageWork stages;
            try {
                auto initial = seed.coordinates;
                const auto check = [](bool ok, const std::string& field = "invalid declared seed") {
                    if (!ok) throw std::invalid_argument(field);
                };
                check((initial.empty() && !seed.previous) || (seed.compositions.empty() && seed.packing_fractions.empty()
                    && seed.amounts.empty() && !seed.temperature && !seed.pressure && !(seed.previous && !initial.empty())));
                if (seed.previous) { check(seed.previous->has_state()); initial = seed.previous->coordinates; }
                else if (initial.empty()) {
                    auto seeded = compiled;
                    check(seed.compositions.empty() || seed.compositions.size() == compiled.problem.phases.size());
                    check(seed.packing_fractions.empty() || seed.packing_fractions.size() == compiled.problem.phases.size());
                    check(seed.amounts.empty() || seed.amounts.size() == compiled.problem.phases.size());
                    // Complete: the seed supplies every coordinate the plain solve's cold start would generate.
                    bool uncomposed = true, complete = (!seed.amounts.empty() || compiled.problem.phases.size() == compiled.pinned_count)
                        && (seed.temperature || !problem.temperature.free) && (seed.pressure || !problem.pressure.free);
                    const auto composed = [&](std::size_t p) { return !seed.compositions.empty() && !seed.compositions[p].empty(); };
                    const auto packed = [&](std::size_t p) { return !seed.packing_fractions.empty() && seed.packing_fractions[p]; };
                    for (std::size_t p = 0; p < seeded.problem.phases.size(); ++p) {
                        auto& phase = seeded.problem.phases[p];
                        if (composed(p)) {
                            check(seed.compositions[p].size() == compiled.declared_supports[p].size()
                                  && std::all_of(seed.compositions[p].begin(), seed.compositions[p].end(), [](double x) { return std::isfinite(x) && x >= 0.0 && x <= 1.0; }));
                            phase.composition_guess = seed.compositions[p];
                            restrict_composition_guess(phase, compiled.declared_supports[p]);
                        }
                        if (packed(p)) phase.packing_guess = seed.packing_fractions[p];
                        check(phase.composition_guess.empty() || phase.composition_guess.size() == phase.support.size());
                        for (const double x : phase.composition_guess)
                            check(std::isfinite(x) && x > 0.0 && x <= 1.0,
                                  "composition for phase " + phase.name + " is outside the admitted range");
                        check(!phase.packing_guess || admitted_packing_fraction(*phase.packing_guess),
                              "packing for phase " + phase.name + " is outside the admitted range");
                        if (!seed.amounts.empty() && (!std::isfinite(seed.amounts[p]) || seed.amounts[p] < 0.0
                            || (phase.is_pinned() && seed.amounts[p] != compiled.phase_amount(phase))))
                            throw std::invalid_argument("amount for phase " + phase.name + " is nonfinite, negative or conflicts with its pin");
                        complete = complete && composed(p) && packed(p);
                        uncomposed = uncomposed && !composed(p);
                    }
                    if (seed.temperature) { check(problem.temperature.free); seeded.problem.temperature = FreeValue::variable(*seed.temperature); }
                    if (seed.pressure) { check(problem.pressure.free); seeded.problem.pressure = FreeValue::variable(*seed.pressure); }
                    if (work.expired()) throw std::runtime_error("elapsed work budget exhausted in initialization");
                    auto speciated = complete || !uncomposed ? std::nullopt : speciated_start(seeded, &work, stages);
                    initial = complete ? seeded.generic_point() : speciated ? std::move(*speciated) : boundary_start(seeded, uncomposed, &work, &stages);
                    const auto offsets = compiled.phase_offsets();
                    for (std::size_t p = 0; p < seeded.problem.phases.size(); ++p) {
                        const auto& phase = seeded.problem.phases[p];
                        if (composed(p)) for (std::size_t j = 0; j < phase.composition_guess.size(); ++j)
                            initial[offsets[p] + j] = std::log(phase.composition_guess[j]);
                        auto cursor = offsets[p] + phase.support.size();
                        if (phase.is_free()) { if (!seed.amounts.empty()) initial[cursor] = seed.amounts[p]; ++cursor; }
                        if (packed(p)) initial[cursor] = std::log(*phase.packing_guess);
                    }
                    auto cursor = compiled.primal_size - compiled.free_intensive_count;
                    if (problem.temperature.free) { if (seed.temperature) initial[cursor] = *seed.temperature; ++cursor; }
                    if (problem.pressure.free && seed.pressure) initial[cursor] = *seed.pressure;
                }
                if (work.consumed >= work.iterations) throw std::runtime_error("iteration budget exhausted by the seed solves");
                check(initial.size() == compiled.unknown_count() && all_numeric(initial));
                const auto offsets = compiled.phase_offsets();
                for (std::size_t p = 0; p < offsets.size(); ++p) {
                    const auto& phase = compiled.problem.phases[p]; auto cursor = offsets[p];
                    for (std::size_t j = 0; j < phase.support.size(); ++j) {
                        const double x = std::exp(initial[cursor++]);
                        check(x > 0.0);
                    }
                    if (phase.is_free()) check(initial[cursor++] >= 0.0);
                    check(admitted_packing_fraction(std::exp(initial[cursor])));
                }
                for (auto j = compiled.primal_size - compiled.free_intensive_count; j < compiled.primal_size; ++j) check(initial[j] > 0.0);
                audit = structural_result(compiled, scope);
                audit->convergence.residual = std::numeric_limits<double>::quiet_NaN();
                auto result = compiled.feasible ? run_driver(model, compiled, std::move(initial), *audit, scope, nullptr, &*audit, &work)
                                                : *audit;
                stages.report(result);
                if (work.expired()) { result.message += "; elapsed work budget exhausted";
                    reject_elapsed(result, scope); }
                if (result.success) {
                    BranchCheck branch;
                    try { branch = compiled.problem.critical ? admit_critical(compiled, result)
                              : check_density_branches(model, compiled, result, {}, true, scope, &work); }
                    catch (const std::exception& error) { branch.message = error.what(); }
                    if (!branch.accepted || work.expired()) {
                        result.success = false;
                        result.message = work.expired() ? "elapsed work budget exhausted in admission" : branch.message;
                        if (work.expired()) reject_elapsed(result, scope); else clear_state(result);
                    }
                }
                scope.finish(result.convergence);
                out.causes[i] = result.message; out.results[i] = std::move(result);
                if (out.results[i]->success && out.results[i]->requested_tolerance_met.value_or(false)) out.selected = i;
            } catch (const std::exception& error) {
                out.causes[i] = std::string(audit ? "driver/post-driver: " : "initialization: ") + error.what();
                if (work.expired()) out.causes[i] += "; elapsed deadline exhausted";
                if (audit) { audit->message = out.causes[i]; stages.report(*audit);
                    if (work.expired()) reject_elapsed(*audit, scope); else clear_state(*audit);
                    out.results[i] = std::move(audit); }
            }
            out.iterations += work.consumed;
            out.attempt_seconds[i] = std::chrono::duration<double>(Clock::now() - attempt_begun).count();
            if (out.selected) break;
        }
    } else out.budget_exhausted = true;
    out.budget_exhausted = out.budget_exhausted || (!out.selected && exhausted());
    out.seconds = std::chrono::duration<double>(Clock::now() - begun).count();
    out.overrun_seconds = limits.seconds ? std::max(0.0, out.seconds - *limits.seconds) : 0.0;
    if (limits.attempt_seconds) for (const double t : out.attempt_seconds)
        out.overrun_seconds = std::max(out.overrun_seconds, t - *limits.attempt_seconds);
    if (scope.diagnostics_on() && retained) for (std::size_t i = 0; i < out.results.size(); ++i)
        if (out.results[i]) observe_endpoint(model, &*retained, *out.results[i], scope, i);
    return out;
}

EquilibriumResult solve_equilibrium_from_coordinates(
    const HelmholtzModel& model,
    const EquilibriumProblem& request,
    const std::vector<double>& coordinates, std::vector<double>* attempt, const DriverLimits* limits
) {
    core::Scope scope("solve_equilibrium");
    CompiledSystem compiled = compile_problem(model, request);
    EquilibriumResult result = structural_result(compiled, scope);
    if (!compiled.feasible) return result;
    if (coordinates.size() != compiled.unknown_count() || !all_numeric(coordinates)) {
        throw std::invalid_argument("continuation coordinates do not match the compiled layout");
    }
    return run_driver(model, compiled, coordinates, std::move(result), scope, attempt, nullptr, limits);
}

}  // namespace epcsaft::equilibrium
