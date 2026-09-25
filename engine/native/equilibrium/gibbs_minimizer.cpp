#include "equilibrium/gibbs_minimizer.hpp"

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
namespace {

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) {
        const double magnitude = std::abs(value);
        if (!std::isfinite(magnitude)) return std::numeric_limits<double>::infinity();
        result = std::max(result, magnitude);
    }
    return result;
}

double specified_amount(const PhaseSpec& phase, const std::vector<double>& feed) {
    if (phase.amount.kind == AmountKind::Free) return 1.0;
    double total = 0.0;
    for (const double value : feed) total += value;
    return phase.amount.fraction * total;
}

double reduced_hessian_minimum_at(
    const CompiledSystem& compiled, PhasePoint prepared,
    const std::vector<double>& multipliers
) {
    if (compiled.pinned_count != 0 || prepared.point.empty()) return 0.0;
    prepared.point.resize(compiled.primal_size);
    const auto& point = prepared.point;
    const auto rows = [&compiled] {
        std::vector<std::size_t> result;
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            if (compiled.rows[row].kind != RowKind::Chemical
                && compiled.rows[row].kind != RowKind::Pressure) result.push_back(row);
        }
        return result;
    }();
    const auto jacobian = compiled.jacobian_at(prepared);
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
    const auto exact = compiled.hessian_at(prepared, multipliers);
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

const char* solver_status_name(Ipopt::SolverReturn status) noexcept {
    switch (status) {
    case Ipopt::SUCCESS: return "SUCCESS";
    case Ipopt::MAXITER_EXCEEDED: return "MAXITER_EXCEEDED";
    case Ipopt::CPUTIME_EXCEEDED: return "CPUTIME_EXCEEDED";
    case Ipopt::STOP_AT_TINY_STEP: return "STOP_AT_TINY_STEP";
    case Ipopt::STOP_AT_ACCEPTABLE_POINT: return "STOP_AT_ACCEPTABLE_POINT";
    case Ipopt::LOCAL_INFEASIBILITY: return "LOCAL_INFEASIBILITY";
    case Ipopt::USER_REQUESTED_STOP: return "USER_REQUESTED_STOP";
    case Ipopt::FEASIBLE_POINT_FOUND: return "FEASIBLE_POINT_FOUND";
    case Ipopt::DIVERGING_ITERATES: return "DIVERGING_ITERATES";
    case Ipopt::RESTORATION_FAILURE: return "RESTORATION_FAILURE";
    case Ipopt::ERROR_IN_STEP_COMPUTATION: return "ERROR_IN_STEP_COMPUTATION";
    case Ipopt::INVALID_NUMBER_DETECTED: return "INVALID_NUMBER_DETECTED";
    case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM: return "TOO_FEW_DEGREES_OF_FREEDOM";
    case Ipopt::INVALID_OPTION: return "INVALID_OPTION";
    case Ipopt::OUT_OF_MEMORY: return "OUT_OF_MEMORY";
    case Ipopt::INTERNAL_ERROR: return "INTERNAL_ERROR";
    case Ipopt::UNASSIGNED: return "UNASSIGNED";
    }
    return "UNKNOWN";
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
    const auto z = species_z(compiled.model->parameters());
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
    if (compiled.pinned_count != 0 && start.multipliers.size() == compiled.rank + compiled.charge_row_count) {
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
        if (!(std::isfinite(packing) && packing >= policy::PACKING_FRACTION_FLOOR.value
              && packing <= policy::PACKING_FRACTION_CEILING.value)) return false;
    }
    return true;
}

double active_constraint_residual(
    const CompiledSystem& compiled, const std::vector<double>& values
) {
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
    const auto keep = [](double& slot, double value) {  // a nonfinite row stays infinite
        slot = std::isfinite(value) ? std::max(slot, std::abs(value)) : std::numeric_limits<double>::infinity();
    };
    for (std::size_t index = 0; index < std::min(rows.size(), values.size()); ++index) {
        const RowKind kind = rows[index].kind;
        keep(result.maximum, values[index]);
        if (kind == RowKind::Chemical || kind == RowKind::Pressure) keep(result.intensive, values[index]);
        else if (kind == RowKind::Moiety) keep(result.balance, values[index]);
        else if (kind == RowKind::Normalization) keep(result.normalization, values[index]);
        else if (kind == RowKind::Electro) keep(result.charge, values[index]);
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
    const double composition_floor = compiled_.pinned_count != 0 && compiled_.is_square()
        ? policy::STABILITY_COMPOSITION_FLOOR.value : policy::EQUILIBRIUM_COMPOSITION_BOX.value;
    std::size_t offset = 0;
    for (const auto& phase : compiled_.problem.phases) {
        // Pinned square cells impose chemical equality and need the full physical support.
        // All-free cells retain the optimizer box and recover trace states by polishing.
        for (std::size_t i = 0; i < phase.support.size(); ++i) {
            x_l[offset] = std::log(composition_floor);
            x_u[offset++] = 0.0;
        }
        if (phase.amount.kind == AmountKind::Free) x_l[offset++] = 0.0;
        x_l[offset] = std::log(policy::PACKING_FRACTION_FLOOR.value);
        x_u[offset++] = std::log(policy::PACKING_FRACTION_CEILING.value);
    }
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

const PhasePoint& GibbsMinimizer::prepared(Ipopt::Index n, const Ipopt::Number* x) {
    std::vector<double> candidate = point(n, x);
    if (!prepared_.has_value() || prepared_->point != candidate) {
        try {
            prepared_ = compiled_.prepare(candidate, true);
        } catch (const std::exception& error) {
            prepared_.reset();
            rows_.reset();
            jacobian_.reset();
            gradient_.reset();
            objective_hessian_.reset();
            constraint_hessian_.reset();
            preparation_error_ = error.what();
            throw;
        }
        rows_.reset();
        jacobian_.reset();
        gradient_.reset();
        objective_hessian_.reset();
        constraint_hessian_.reset();
        ++work_.prepared_points;
    }
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
    Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& value
) {
    (void)new_x;
    ++evaluations_;
    try {
        value = compiled_.pinned_count == 0 ? rows_at(n, x).objective : 0.0;
        return std::isfinite(value);
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_grad_f(
    Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f
) {
    (void)new_x;
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
    Ipopt::Index n, const Ipopt::Number* x, bool new_x,
    Ipopt::Index m, Ipopt::Number* g
) {
    (void)new_x;
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
    Ipopt::Index n, const Ipopt::Number* x, bool new_x,
    Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
    bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    (void)new_x; (void)new_lambda;
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

namespace {

EquilibriumResult structural_result(const CompiledSystem& compiled, core::Scope& scope) {
    EquilibriumResult result;
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
    std::vector<double> initial, EquilibriumResult result, core::Scope& scope
) {
    const double row_tolerance = compiled.problem.row_tolerance;
    // The optimizer works one policy decade below the caller criterion, which alone classifies.
    // A pinned square system has a constant objective and is Newton-polished below, so it stops there.
    const double inner_tolerance = std::min(row_tolerance, policy::EQUILIBRIUM_ROW_RESIDUAL.value)
        * (compiled.pinned_count != 0 && compiled.is_square() ? 1.0 : policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value);
    Ipopt::SmartPtr<GibbsMinimizer> nlp = new GibbsMinimizer(compiled, std::move(initial), &scope);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> application = ::IpoptApplicationFactory();
    const double objective_scale = 1.0 / compiled.feed_total;
    application->Options()->SetIntegerValue("print_level", 0);
    application->Options()->SetStringValue("nlp_scaling_method", "none");
    application->Options()->SetNumericValue("obj_scaling_factor", objective_scale);
    application->Options()->SetStringValue("hessian_approximation", "exact");
    application->Options()->SetNumericValue("tol", inner_tolerance);
    application->Options()->SetNumericValue("acceptable_tol", inner_tolerance);
    application->Options()->SetIntegerValue("acceptable_iter", 0);
    application->Options()->SetNumericValue("dual_inf_tol", inner_tolerance);
    application->Options()->SetNumericValue("constr_viol_tol", inner_tolerance);
    application->Options()->SetIntegerValue("max_iter",
        static_cast<int>(policy::EQUILIBRIUM_MAX_ITERATIONS.value));
    application->Initialize();
    application->OptimizeTNLP(nlp);
    result.solver_status = solver_status_name(nlp->status());
    std::vector<double> point = nlp->solution();
    if (point.empty()) point = boundary_start(compiled);
    std::vector<double> full = compiled.pinned_count == 0 ? add_multipliers(compiled, point) : point;
    // The optimizer's amount-weighted rows underresolve intensive trace defects.
    // Exact Newton steps contract raw rows inside the physical composition cutoff.
    // Curvature may raise classified rows temporarily; return the best classified
    // point visited, including the optimizer start, rather than the last step.
    std::size_t polish_evaluations = 0;
    std::string rejected_preparation = nlp->preparation_error().value_or("");
    // A budget stop polishes too: an unattainable request can exhaust Ipopt after passing a
    // better point, and the polish never returns a worse classified point than its start.
    if ((nlp->status() == Ipopt::SUCCESS || nlp->status() == Ipopt::STOP_AT_ACCEPTABLE_POINT
         || nlp->status() == Ipopt::STOP_AT_TINY_STEP || nlp->status() == Ipopt::MAXITER_EXCEEDED
         || nlp->status() == Ipopt::CPUTIME_EXCEEDED || nlp->status() == Ipopt::RESTORATION_FAILURE)
        && compiled.is_square()) {
        auto at = compiled.prepare(full, true);
        auto rows_before = compiled.evaluate_at(at).values;
        std::vector<double> best = full;
        double best_classified = maximum_abs(compiled.classified_residuals(at, rows_before));
        ++polish_evaluations;
        for (std::size_t step = 0; step < static_cast<std::size_t>(policy::EQUILIBRIUM_POLISH_STEPS.value); ++step) {
            const double defect = maximum_abs(rows_before);
            const double classified_before = maximum_abs(compiled.classified_residuals(at, rows_before));
            if (!(defect > 0.0) || !std::isfinite(defect)) break;
            const auto jacobian = compiled.jacobian_at(at);
            const auto n = static_cast<Eigen::Index>(full.size());
            if (jacobian.size() != full.size() * full.size() || rows_before.size() != full.size()) break;
            const Eigen::VectorXd delta = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor>>(jacobian.data(), n, n).fullPivLu().solve(
                    -Eigen::Map<const Eigen::VectorXd>(rows_before.data(), n));
            bool accepted = false;
            // Bounded damping retains raw contraction and the physical cutoff.
            // Best-point selection below protects the returned classified residual.
            for (const double alpha : {1.0, 0.5, 0.25}) {
                std::vector<double> candidate = full;
                for (Eigen::Index i = 0; i < n; ++i)
                    candidate[static_cast<std::size_t>(i)] += alpha * delta(i);
                PhasePoint at_candidate;
                std::vector<double> rows_after;
                try {
                    at_candidate = compiled.prepare(candidate, true);
                    rows_after = compiled.evaluate_at(at_candidate).values;
                    ++polish_evaluations;
                } catch (const std::exception& error) {
                    rejected_preparation = error.what();
                    continue;
                }
                bool inside = true;
                for (const auto& phase : at_candidate.composition)
                    for (const double x : phase)
                        inside = inside && (x == 0.0 || x >= policy::STABILITY_COMPOSITION_FLOOR.value);
                const double after = maximum_abs(rows_after);
                const double classified_after = maximum_abs(compiled.classified_residuals(at_candidate, rows_after));
                scope.probe("polish_step", {defect, after, classified_before, classified_after, alpha},
                    static_cast<std::uint32_t>(step));
                if (!inside || !(after < defect)) continue;
                full = std::move(candidate);
                at = std::move(at_candidate);
                rows_before = std::move(rows_after);
                if (classified_after <= best_classified) {
                    best = full;
                    best_classified = classified_after;
                }
                accepted = true;
                break;
            }
            if (!accepted) break;
        }
        full = std::move(best);
    }
    PhasePoint prepared;
    try {
        prepared = compiled.prepare(full);
    } catch (const std::invalid_argument& error) {
        if (std::string(error.what()).rfind("reference_unavailable:", 0) != 0) throw;
        result.success = false;
        result.message = error.what();
        result.convergence.stop = core::Stop::left_domain;
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        result.iterations = result.convergence.iterations = nlp->iterations();
        result.evaluations = result.convergence.evaluations = nlp->evaluations() + polish_evaluations;
        clear_state(result);
        scope.finish(result.convergence);
        return result;
    }
    const auto evaluation = compiled.evaluate_at(prepared);
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
    if (row_tolerance == policy::EQUILIBRIUM_ROW_RESIDUAL.value) {
        convergence.criterion = &policy::EQUILIBRIUM_ROW_RESIDUAL;
    } else {
        // The caller-valued criterion lives with the record so copies outlive the request.
        convergence.criterion_owner = std::make_shared<const policy::Tolerance>(policy::Tolerance{
            "EQUILIBRIUM_REQUESTED_ROW_RESIDUAL", row_tolerance,
            policy::Category::AlgorithmTolerance, "caller-supplied equilibrium row tolerance"});
        convergence.criterion = convergence.criterion_owner.get();
    }
    if (evaluation.floor > 0.0) convergence.floor = evaluation.floor;
    convergence.amplification = nlp->regularization_size();
    result.iterations = nlp->iterations();
    result.evaluations = nlp->evaluations() + polish_evaluations;
    result.hessian_calls = nlp->hessian_calls();
    convergence.iterations = result.iterations;
    convergence.evaluations = result.evaluations;
    switch (nlp->status()) {
    case Ipopt::INVALID_NUMBER_DETECTED: convergence.stop = core::Stop::nonfinite; break;
    case Ipopt::DIVERGING_ITERATES: convergence.stop = core::Stop::left_domain; break;
    case Ipopt::ERROR_IN_STEP_COMPUTATION: convergence.stop = core::Stop::singular; break;
    case Ipopt::LOCAL_INFEASIBILITY:
    case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM: convergence.stop = core::Stop::infeasible; break;
    default:
        convergence.stop = core::classify(
            convergence.residual, *convergence.criterion,
            convergence.floor, convergence.contraction,
            (nlp->status() == Ipopt::MAXITER_EXCEEDED ||
             nlp->status() == Ipopt::CPUTIME_EXCEEDED), convergence.amplification
        );
    }
    if (std::isfinite(convergence.residual)) {
        result.requested_tolerance_met = convergence.residual <= row_tolerance;
    }
    result.message = core::describe(convergence);
    scope.finish(convergence);
    const auto offsets = compiled.phase_offsets();
    std::vector<std::vector<double>> compositions;
    std::vector<double> amounts;
    std::vector<double> volumes;
    std::size_t intensive_offset = compiled.primal_size - compiled.free_intensive_count;
    const double candidate_temperature = compiled.problem.temperature.free
        ? full[intensive_offset] : compiled.problem.temperature.value;
    bool topology_ok = true;
    std::string topology_failure;
    const auto topology = check_topology(compiled, full);
    result.topology_event = topology.event;
    result.topology_phase = topology.phase;
    result.topology_pair = topology.pair;
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
    const bool solver_ok = nlp->status() == Ipopt::SUCCESS
        || result.convergence.stop == core::Stop::converged;
    const bool accepted = solver_ok
        && (result.convergence.stop == core::Stop::converged
            || result.convergence.stop == core::Stop::converged_to_floor)
        && topology_ok;
    if (!accepted) {
        result.success = false;
        if (!topology_ok) result.message = topology_failure;
        if (!rejected_preparation.empty()) {
            result.message += "; trial preparation failed: " + rejected_preparation;
        }
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
    if (compiled.pinned_count == 0) {
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
    result.reduced_hessian_minimum = reduced_hessian_minimum_at(
        compiled, std::move(prepared), result.constraint_multipliers);
    if (compiled.pinned_count == 0 && result.reduced_hessian_minimum <= 0.0) {
        result.success = false;
        result.message = "reduced_hessian_stability is not positive definite";
        if (!rejected_preparation.empty()) {
            result.message += "; trial preparation failed: " + rejected_preparation;
        }
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

bool zero_amount(const PhaseSpec& phase) {
    return phase.amount.kind == AmountKind::Pinned && phase.amount.fraction == 0.0;
}

// A cold reactive boundary at fixed T and free P: the phases present (nonzero amount) alone, in reaction
// equilibrium at the pressure seed, are the solution's present phases up to their pressure dependence. Their
// compositions replace the caller's guesses before the ordinary boundary start, and an unguessed zero-amount
// phase starts from the first present phase on its support (#152); none if any of it fails.
std::optional<std::vector<double>> speciated_start(const CompiledSystem& compiled) {
    const auto& problem = compiled.problem;
    EquilibriumProblem present = problem;
    present.phases.erase(std::remove_if(present.phases.begin(), present.phases.end(), zero_amount), present.phases.end());
    if (!compiled.has_reaction || !problem.pressure.free || problem.temperature.free || present.phases.empty()
        || present.phases.size() == problem.phases.size()) return std::nullopt;
    present.pressure = FreeValue::fixed(problem.pressure.seed_or(1.0e5));
    try {
        const auto solved = solve_equilibrium(*compiled.model, present);
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
        return boundary_start(seeded);
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
    std::vector<std::vector<double>> compositions;
    std::vector<std::size_t> eta_index;
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
            ? full[cursor++] : specified_amount(spec, compiled.n_feed);
        eta_index.push_back(cursor);
        compositions.push_back(std::move(composition));
        if (spec.amount.kind == AmountKind::Free
            && amount / compiled.feed_total <= policy::STABILITY_PHASE_AMOUNT_FLOOR.value
            && result.event == TopologyEvent::None) {
            result.event = TopologyEvent::PhaseDisappeared;
            result.phase = phase;
        }
    }
    for (std::size_t left = 0; left + 1 < compositions.size(); ++left) {
        for (std::size_t right = left + 1; right < compositions.size(); ++right) {
            const double distance = std::max(
                std::abs(full[eta_index[left]] - full[eta_index[right]]),
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

EquilibriumResult solve_equilibrium(const HelmholtzModel& model, const EquilibriumProblem& request,
                                    const std::optional<EquilibriumResult>& start) {
    core::Scope scope("solve_equilibrium");  // compilation belongs to the solve's span
    return solve_compiled(compile_problem(model, request), start, &scope);
}

EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start,
                                 core::Scope* outer) {
    std::optional<core::Scope> own;
    auto& scope = outer != nullptr ? *outer : own.emplace("solve_equilibrium");
    const auto& model = *compiled.model;
    const EquilibriumResult structural = structural_result(compiled, scope);
    if (!compiled.feasible) return structural;
    const bool use_start = start.has_value() && usable_start(compiled, *start);
    auto seeded = start.has_value() ? std::nullopt : speciated_start(compiled);
    EquilibriumResult result;
    // The speciated seed when there is one; the caller's start alone or after the seed's failure.
    for (int attempt = seeded ? 0 : 1; attempt < 2 && !result.success; ++attempt) {
        std::vector<double> initial;
        if (attempt == 0) {
            initial = std::move(*seeded);
        } else {
            initial = use_start ? compiled.generic_point() : boundary_start(compiled, !(start && start->success));
            if (start.has_value()) initial = start_point(compiled, *start, std::move(initial));
        }
        try {
            result = run_driver(model, compiled, std::move(initial), structural, scope);
        } catch (const std::exception&) {
            if (attempt == 0) continue;  // the seeded attempt's exception falls back to the declared start
            throw;
        }
        if (!result.success) continue;
        BranchCheck branches;
        try {
            branches = check_density_branches(model, compiled, result, {}, true, scope);
        } catch (const std::exception& error) {
            branches.message = std::string("density branch admission failed: ") + error.what();
        }
        if (!branches.accepted) {
            result.success = false;
            result.message = branches.message;
            clear_state(result);
        }
    }
    return result;
}

EquilibriumResult solve_equilibrium_from_coordinates(
    const HelmholtzModel& model,
    const EquilibriumProblem& request,
    const std::vector<double>& coordinates
) {
    core::Scope scope("solve_equilibrium");
    CompiledSystem compiled = compile_problem(model, request);
    EquilibriumResult result = structural_result(compiled, scope);
    if (!compiled.feasible) return result;
    if (coordinates.size() != compiled.unknown_count() || !numeric_state(coordinates)) {
        throw std::invalid_argument("continuation coordinates do not match the compiled layout");
    }
    return run_driver(model, compiled, coordinates, std::move(result), scope);
}

}  // namespace epcsaft::equilibrium
