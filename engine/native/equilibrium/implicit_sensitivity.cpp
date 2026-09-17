#include "equilibrium/gibbs_minimizer.hpp"

#include "core/trace.hpp"
#include "core/numerical_policy.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

Eigen::MatrixXd dense_matrix(
    const std::vector<double>& values, std::size_t rows, std::size_t columns
) {
    if (values.size() != rows * columns) throw std::invalid_argument("matrix size");
    Eigen::MatrixXd result(static_cast<Eigen::Index>(rows),
                           static_cast<Eigen::Index>(columns));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < columns; ++column) {
            result(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                values[row * columns + column];
        }
    }
    return result;
}

double condition_estimate(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() == 0 || matrix.cols() == 0) return 0.0;
    const auto singular = matrix.jacobiSvd().singularValues();
    if (singular.size() == 0 || singular(singular.size() - 1) <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return singular(0) / singular(singular.size() - 1);
}

std::size_t parameter_count(
    const std::vector<double>& parameter_rows, std::size_t rows
) {
    if (rows == 0 || parameter_rows.empty() || parameter_rows.size() % rows != 0) return 0;
    return parameter_rows.size() / rows;
}

void validate_parameter_rows(
    const std::vector<double>& parameter_rows, std::size_t rows
) {
    if (!std::all_of(parameter_rows.begin(), parameter_rows.end(),
                     [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument("invalid parameter derivative");
    }
    if (!parameter_rows.empty() && (rows == 0 || parameter_rows.size() % rows != 0)) {
        throw std::invalid_argument("parameter derivative shape");
    }
}

std::vector<double> flatten(const Eigen::MatrixXd& matrix) {
    std::vector<double> result(
        static_cast<std::size_t>(matrix.rows() * matrix.cols()), 0.0);
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            result[static_cast<std::size_t>(row * matrix.cols() + column)] =
                matrix(row, column);
        }
    }
    return result;
}

bool valid_point(const std::vector<double>& point) {
    return !point.empty() && std::all_of(
        point.begin(), point.end(), [](double value) { return std::isfinite(value); });
}

double maximum_abs(const Eigen::MatrixXd& matrix) {
    double result = 0.0;
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            result = std::max(result, std::abs(matrix(row, column)));
        }
    }
    return result;
}

bool residual_is_small(const Eigen::MatrixXd& matrix,
                       const Eigen::MatrixXd& solution,
                       const Eigen::MatrixXd& right) {
    const double scale = std::max(1.0, maximum_abs(right));
    return maximum_abs(matrix * solution - right)
        <= policy::EQUILIBRIUM_ROW_RESIDUAL.value * scale;
}

bool square_and_numeric(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() != matrix.cols()) return false;
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            if (!std::isfinite(matrix(row, column))) return false;
        }
    }
    return true;
}

double relative_residual(const Eigen::MatrixXd& matrix,
                         const Eigen::MatrixXd& solution,
                         const Eigen::MatrixXd& right) {
    const double scale = std::max(1.0, maximum_abs(right));
    return maximum_abs(matrix * solution - right) / scale;
}

}  // namespace

SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
) {
    core::Scope scope("inspect_implicit_sensitivity");
    const std::size_t rows = compiled.rows.size();
    const std::size_t dimension = point.size();
    if (!valid_point(point)) throw std::invalid_argument("invalid sensitivity point");
    validate_parameter_rows(parameter_rows, rows);
    const auto jacobian = compiled.jacobian(point);
    if (jacobian.size() != rows * dimension) {
        throw std::invalid_argument("equilibrium Jacobian shape");
    }
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    const auto lu = matrix.fullPivLu();
    SensitivityInfo result;
    result.equations = rows;
    result.parameters = parameter_count(parameter_rows, rows);
    result.rank = static_cast<std::size_t>(lu.rank());
    result.condition_estimate = condition_estimate(matrix);
    result.nonsingular = square_and_numeric(matrix) && result.rank == dimension;
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    scope.probe("condition", result.condition_estimate);
    scope.probe("jacobian_entries", static_cast<double>(jacobian.size()));
    scope.probe("parameter_entries", static_cast<double>(parameter_rows.size()));
    scope.probe("parameters", static_cast<double>(result.parameters));
    scope.probe("nonsingular", result.nonsingular ? 1.0 : 0.0);
    scope.probe("rank", static_cast<double>(result.rank));
    scope.probe("condition", result.condition_estimate);
    return result;
}

std::vector<double> solve_implicit(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
) {
    core::Scope scope("implicit_sensitivity");
    const std::size_t rows = compiled.rows.size();
    const std::size_t dimension = point.size();
    if (!valid_point(point)) throw std::invalid_argument("invalid sensitivity point");
    validate_parameter_rows(parameter_rows, rows);
    const std::size_t parameters = parameter_count(parameter_rows, rows);
    if (parameters == 0) return {};
    const auto jacobian = compiled.jacobian(point);
    if (jacobian.size() != rows * dimension || rows != dimension) return {};
    const auto info = inspect_sensitivity(compiled, point, parameter_rows);
    if (!info.nonsingular) return {};
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    Eigen::MatrixXd right(static_cast<Eigen::Index>(rows),
                          static_cast<Eigen::Index>(parameters));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
            right(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(parameter)) =
                -parameter_rows[row * parameters + parameter];
        }
    }
    const Eigen::MatrixXd result = matrix.fullPivLu().solve(right);
    if (!result.allFinite() || !residual_is_small(matrix, result, right)) return {};
    scope.probe("parameters", static_cast<double>(parameters));
    scope.probe("max_solution", maximum_abs(result));
    scope.probe("relative_residual", relative_residual(matrix, result, right));
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    return flatten(result);
}

std::vector<double> sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& lever
) {
    core::Scope scope("sensitivity");
    const auto F_theta = compiled.lever_jacobian(point);
    if (compiled.rows.empty() || F_theta.size() % compiled.rows.size() != 0) {
        throw std::invalid_argument("lever Jacobian shape");
    }
    const std::size_t levers = F_theta.size() / compiled.rows.size();
    if (lever.size() != levers)
        throw std::invalid_argument("lever must contain T, P, feed, and active parameter weights");
    if (F_theta.size() != compiled.rows.size() * levers) {
        throw std::invalid_argument("lever Jacobian shape");
    }
    std::vector<double> selected(compiled.rows.size(), 0.0);
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        double F_theta_value = 0.0;
        for (std::size_t column = 0; column < levers; ++column)
            F_theta_value += F_theta[levers * row + column] * lever[column];
        selected[row] = -(-F_theta_value);
    }
    return solve_implicit(compiled, point, selected);
}

// ------------------------------------------------------------- solved-state actions (E5)

namespace {

using Series = std::vector<double>;

Series product(const Series& left, const Series& right) {
    Series result(left.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i)
        for (std::size_t j = 0; j <= i; ++j) result[i] += left[j] * right[i - j];
    return result;
}

Series inverse(const Series& value) {
    if (value.empty() || !std::isfinite(value.front()) || value.front() == 0.0)
        throw std::domain_error("series inverse has a zero or nonfinite center");
    Series result(value.size(), 0.0);
    result.front() = 1.0 / value.front();
    for (std::size_t i = 1; i < value.size(); ++i) {
        for (std::size_t j = 1; j <= i; ++j) result[i] -= value[j] * result[i - j];
        result[i] /= value.front();
    }
    return result;
}

Series combine(const Series& left, const Series& right, double sign) {
    Series result(left.size());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = left[i] + sign * right[i];
    return result;
}

// Logarithm of a series with a positive centre: k s_0 L_k = k s_k - sum_{j<k} j L_j s_{k-j}.
Series logarithm(const Series& value) {
    if (value.empty() || !std::isfinite(value.front()) || !(value.front() > 0.0))
        throw std::domain_error("series logarithm has a nonpositive or nonfinite center");
    Series result(value.size(), 0.0);
    result.front() = std::log(value.front());
    for (std::size_t k = 1; k < value.size(); ++k) {
        double sum = static_cast<double>(k) * value[k];
        for (std::size_t j = 1; j < k; ++j) sum -= static_cast<double>(j) * result[j] * value[k - j];
        result[k] = sum / (static_cast<double>(k) * value.front());
    }
    return result;
}

bool is_phase_kind(SolvedStateObservableKind kind) {
    return kind != SolvedStateObservableKind::Temperature && kind != SolvedStateObservableKind::Pressure;
}

bool has_component(SolvedStateObservableKind kind) {
    return kind == SolvedStateObservableKind::PhaseComponentAmount
        || kind == SolvedStateObservableKind::PhaseMoleFraction
        || kind == SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt
        || kind == SolvedStateObservableKind::PhaseLogActivity;
}

bool same_observable(const SolvedStateObservable& a, const SolvedStateObservable& b) {
    if (a.kind != b.kind) return false;
    if (is_phase_kind(a.kind) && a.phase != b.phase) return false;
    if (has_component(a.kind) && a.component != b.component) return false;
    return a.kind != SolvedStateObservableKind::PhaseProperty || a.property == b.property;
}

const char* observable_unit(const SolvedStateObservable& observable) {
    switch (observable.kind) {
    case SolvedStateObservableKind::Temperature: return "K";
    case SolvedStateObservableKind::Pressure: case SolvedStateObservableKind::PhasePressure: return "Pa";
    case SolvedStateObservableKind::PhaseAmount: case SolvedStateObservableKind::PhaseComponentAmount: return "mol";
    case SolvedStateObservableKind::PhaseVolume: return "m^3";
    case SolvedStateObservableKind::PhaseMolarDensity: return "mol/m^3";
    case SolvedStateObservableKind::PhaseProperty:
        return static_cast<std::size_t>(observable.property) % 5 < 2 ? "J/mol" : "J/(mol K)";
    default: return "1";
    }
}

const char* property_name(PropertyObservable observable) {
    static const char* const names[] = {
        "IdealEnthalpy", "IdealInternalEnergy", "IdealEntropy", "IdealIsobaricHeatCapacity",
        "IdealIsochoricHeatCapacity", "ResidualEnthalpy", "ResidualInternalEnergy", "ResidualEntropy",
        "ResidualIsobaricHeatCapacity", "ResidualIsochoricHeatCapacity", "TotalEnthalpy",
        "TotalInternalEnergy", "TotalEntropy", "TotalIsobaricHeatCapacity", "TotalIsochoricHeatCapacity"};
    const auto index = static_cast<std::size_t>(observable);
    return index < 15 ? names[index] : "unknown";
}

// Every observable except a fixed T or P input depends on the solved state.
bool state_dependent(const SolvedStateObservable& observable, const EquilibriumProblem& problem) {
    if (observable.kind == SolvedStateObservableKind::Temperature) return problem.temperature.free;
    if (observable.kind == SolvedStateObservableKind::Pressure) return problem.pressure.free;
    return true;
}

void validate_request(const CompiledSystem& compiled, const EquilibriumResult& central,
                      const SolvedStateActionRequest& request) {
    const auto& problem = compiled.problem;
    const std::size_t components = compiled.model->parameters().component_count();
    const std::size_t actives = compiled.model->active_parameters().size();
    if (request.observables.empty()) throw std::invalid_argument("no solved-state observables");
    for (std::size_t i = 0; i < request.observables.size(); ++i) {
        const auto& item = request.observables[i];
        if (is_phase_kind(item.kind) && item.phase >= problem.phases.size())
            throw std::invalid_argument("solved-state observable phase ordinal");
        if (has_component(item.kind) && item.component >= components)
            throw std::invalid_argument("solved-state observable component ordinal");
        if (item.kind == SolvedStateObservableKind::PhaseProperty && static_cast<std::size_t>(item.property) >= 15)
            throw std::invalid_argument("solved-state property observable");
        const bool log_activity = item.kind == SolvedStateObservableKind::PhaseLogActivity;
        if (log_activity || item.kind == SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt) {
            const auto& support = problem.phases[item.phase].support;
            if (std::find(support.begin(), support.end(), item.component) == support.end())
                throw std::invalid_argument(std::string(log_activity ? "log activity" : "residual chemical potential")
                                            + " outside the phase support");
        }
        for (std::size_t j = 0; j < i; ++j)
            if (same_observable(item, request.observables[j]))
                throw std::invalid_argument("duplicate solved-state observable");
    }
    if (request.directions.empty() || request.directions.size() > 3)
        throw std::invalid_argument("solved-state action order");
    for (const auto& direction : request.directions) {
        if (!std::isfinite(direction.temperature_k) || !std::isfinite(direction.pressure_pa)
            || direction.feed_amounts_mol.size() != components || direction.active_parameters.size() != actives
            || !std::all_of(direction.feed_amounts_mol.begin(), direction.feed_amounts_mol.end(),
                            [](double v) { return std::isfinite(v); })
            || !std::all_of(direction.active_parameters.begin(), direction.active_parameters.end(),
                            [](double v) { return std::isfinite(v); }))
            throw std::invalid_argument("solved-state action direction");
        if (problem.temperature.free && direction.temperature_k != 0.0)
            throw std::invalid_argument("free temperature is a solved output, not an external control");
        if (problem.pressure.free && direction.pressure_pa != 0.0)
            throw std::invalid_argument("free pressure is a solved output, not an external control");
    }
    if (central.rows.size() != compiled.rows.size() || central.phase_count != problem.phases.size())
        throw std::invalid_argument("central result does not match the compiled system");
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        const auto& a = central.rows[row];
        const auto& b = compiled.rows[row];
        if (a.kind != b.kind || a.phase != b.phase || a.index != b.index)
            throw std::invalid_argument("central rows do not match the compiled system");
    }
    if (central.has_state() && central.coordinates.size() != compiled.unknown_count())
        throw std::invalid_argument("central coordinates do not match the compiled layout");
}

// A constant supplied Gibbs vector cannot carry a temperature-dependent reaction target.
bool thermal_reference_missing(const CompiledSystem& compiled, double temperature) {
    const auto& problem = compiled.problem;
    if (problem.standard_gibbs_j_per_mol.empty() || problem.reactions.empty()) return false;
    const auto targets = reaction_target_series(temperature, problem.reactions, 3);
    for (std::size_t row = 0; row < problem.reactions.size(); ++row) {
        const double scale = std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(targets[0][row]));
        double power = 1.0;
        for (std::size_t k = 1; k <= 3; ++k) {
            power *= temperature;
            if (std::abs(targets[k][row]) * power > scale) return true;
        }
    }
    return false;
}

struct Center final {
    PhasePoint prepared;
    Eigen::MatrixXd jacobian;
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
};

struct Path final {
    RowSeries series;
    SolvedStateActionStatus status = SolvedStateActionStatus::Available;
};

SolvedStateActionStatus energy_status(PropertyActionStatus status) {
    return status == PropertyActionStatus::Nonfinite ? SolvedStateActionStatus::NonfiniteJacobian
        : solved_state_status(status);
}

// The coefficient recurrence: J y_k = -[t^k] F(y(t), q(t)) with y_k unknown set to zero,
// re-evaluated after each solve so the retained residual is the owned row residual.
Path follow(const CompiledSystem& compiled, const Center& center, const std::vector<double>& control,
            std::size_t order, SolvedStateActionDiagnostics& diagnostics) {
    Path path;
    const std::size_t n = center.prepared.point.size();
    std::vector<std::vector<double>> state(order + 1, std::vector<double>(n, 0.0));
    state[0] = center.prepared.point;
    std::vector<double> scale(order + 1, 1.0);
    for (std::size_t k = 1; k <= order; ++k) {
        const std::vector<std::vector<double>> prefix(state.begin(), state.begin() + static_cast<std::ptrdiff_t>(k + 1));
        const auto series = compiled.row_series(center.prepared, prefix, control, k);
        diagnostics.maximum_energy_order = std::max(diagnostics.maximum_energy_order, series.maximum_energy_order);
        if (series.energy != PropertyActionStatus::Available) { path.status = energy_status(series.energy); return path; }
        Eigen::VectorXd right(static_cast<Eigen::Index>(n));
        for (std::size_t row = 0; row < n; ++row) right(static_cast<Eigen::Index>(row)) = -series.rows[k][row];
        if (!right.allFinite()) { path.status = SolvedStateActionStatus::NonfiniteJacobian; return path; }
        scale[k] = std::max(1.0, right.cwiseAbs().maxCoeff());
        const Eigen::VectorXd solution = center.lu.solve(right);
        ++diagnostics.linear_solves;
        if (!solution.allFinite()) { path.status = SolvedStateActionStatus::NonfiniteJacobian; return path; }
        for (std::size_t row = 0; row < n; ++row) state[k][row] = solution(static_cast<Eigen::Index>(row));
    }
    path.series = compiled.row_series(center.prepared, state, control, order);
    diagnostics.maximum_energy_order = std::max(diagnostics.maximum_energy_order, path.series.maximum_energy_order);
    if (path.series.energy != PropertyActionStatus::Available) { path.status = energy_status(path.series.energy); return path; }
    for (std::size_t k = 1; k <= order; ++k) {
        double residual = 0.0;
        for (const double value : path.series.rows[k]) residual = std::max(residual, std::abs(value));
        residual /= scale[k];
        if (!std::isfinite(residual)) { path.status = SolvedStateActionStatus::NonfiniteJacobian; return path; }
        diagnostics.maximum_linear_residual = std::max(diagnostics.maximum_linear_residual.value_or(0.0), residual);
        if (residual > policy::EQUILIBRIUM_ROW_RESIDUAL.value) path.status = SolvedStateActionStatus::LinearSolveResidual;
    }
    return path;
}

// Physical map of one phase along a path: normalised composition, density and pressure.
struct PhaseMap final {
    Series amount, volume, density, packing, pressure;
    std::vector<Series> mole_fraction, potential;
};

PhaseMap phase_map(const RowSeries& series, std::size_t phase, std::size_t components) {
    const auto& item = series.phase[phase];
    const std::size_t order = item.point.size();
    PhaseMap result;
    result.amount = item.amount;
    result.packing = item.packing;
    Series volume(order), total(order, 0.0), potential_v(order);
    for (std::size_t k = 0; k < order; ++k) {
        volume[k] = item.point[k][1];
        potential_v[k] = item.gradient[k][1];
        for (std::size_t i = 0; i < components; ++i) total[k] += item.point[k][2 + i];
    }
    result.volume = volume;
    result.density = inverse(volume);
    const auto inverse_total = inverse(total);
    for (std::size_t i = 0; i < components; ++i) {
        Series x(order), mu(order);
        for (std::size_t k = 0; k < order; ++k) { x[k] = item.point[k][2 + i]; mu[k] = item.gradient[k][2 + i]; }
        result.mole_fraction.push_back(product(x, inverse_total));
        result.potential.push_back(mu);
    }
    Series pressure = product(series.temperature, combine(result.density, potential_v, -1.0));
    for (double& value : pressure) value *= GAS_CONSTANT_J_PER_MOL_K;
    result.pressure = pressure;
    return result;
}

Series observable_series(const SolvedStateObservable& observable, const RowSeries& series,
                         const std::vector<PhaseMap>& maps) {
    switch (observable.kind) {
    case SolvedStateObservableKind::Temperature: return series.temperature;
    case SolvedStateObservableKind::Pressure: return series.pressure;
    case SolvedStateObservableKind::PhaseAmount: return maps[observable.phase].amount;
    case SolvedStateObservableKind::PhaseComponentAmount:
        return product(maps[observable.phase].amount, maps[observable.phase].mole_fraction[observable.component]);
    case SolvedStateObservableKind::PhaseMoleFraction: return maps[observable.phase].mole_fraction[observable.component];
    case SolvedStateObservableKind::PhaseVolume: return product(maps[observable.phase].amount, maps[observable.phase].volume);
    case SolvedStateObservableKind::PhaseMolarDensity: return maps[observable.phase].density;
    case SolvedStateObservableKind::PhasePackingFraction: return maps[observable.phase].packing;
    case SolvedStateObservableKind::PhasePressure: return maps[observable.phase].pressure;
    case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:
        return maps[observable.phase].potential[observable.component];
    case SolvedStateObservableKind::PhaseLogActivity:
        return combine(logarithm(product(maps[observable.phase].mole_fraction[observable.component],
                                         maps[observable.phase].density)),
                       maps[observable.phase].potential[observable.component], 1.0);
    default: return {};
    }
}

// The P4 chain of one phase along a curved path: [t^1] g = Dg[z1]; [t^2] g = Dg[z2] +
// D2g[z1,z1]/2; [t^3] g = Dg[z3] + D2g[z1,z2] + D3g[z1,z1,z1]/6 in the (T, rho, x, theta) chart.
struct PropertyChain final {
    std::vector<std::optional<double>> value, coefficient;
    std::vector<PropertyActionStatus> status;
    PropertyActionDiagnostics diagnostics;
    std::size_t maximum_energy_order = 0;
};

PropertyChain property_chain(const HelmholtzModel& model, const RowSeries& series, const PhaseMap& map,
                             std::size_t phase, const std::vector<PropertyObservable>& observables,
                             std::size_t order) {
    const auto& item = series.phase[phase];
    const std::size_t components = model.parameters().component_count();
    // Coefficient k as a chart direction; order zero is the zero direction (values only).
    const auto direction = [&](std::size_t k) {
        DensityActionDirection result;
        result.temperature = k == 0 ? 0.0 : series.temperature[k];
        result.molar_density = k == 0 ? 0.0 : map.density[k];
        for (std::size_t i = 0; i < components; ++i) result.composition.push_back(k == 0 ? 0.0 : map.mole_fraction[i][k]);
        // The solved coefficients satisfy the normalisation and charge rows to the linear
        // residual policy; the P4 trust boundary requires them exactly, so the roundoff
        // component along (1, z) is removed before the call.
        const auto charges = species_z(model.parameters());
        double sum = 0.0, charge = 0.0, zz = 0.0, z1 = 0.0;
        for (std::size_t i = 0; i < components; ++i) {
            sum += result.composition[i]; charge += charges[i] * result.composition[i];
            zz += charges[i] * charges[i]; z1 += charges[i];
        }
        const double c = static_cast<double>(components), det = c * zz - z1 * z1;
        const double alpha = det != 0.0 ? (sum * zz - charge * z1) / det : sum / c;
        const double beta = det != 0.0 ? (c * charge - z1 * sum) / det : 0.0;
        for (std::size_t i = 0; i < components; ++i) result.composition[i] -= alpha + beta * charges[i];
        for (std::size_t a = 0; a < model.active_parameters().size(); ++a) {
            const auto& active = model.active_parameters()[a];
            result.parameter.push_back({active.family, active.indices, k == 0 ? 0.0 : item.point[k][2 + components + a]});
        }
        return result;
    };
    std::vector<double> central_composition;
    for (std::size_t i = 0; i < components; ++i) central_composition.push_back(map.mole_fraction[i][0]);
    struct Term final { std::vector<std::size_t> orders; double weight; };
    std::vector<Term> terms;
    if (order == 0) terms = {{{0}, 1.0}};
    else if (order == 1) terms = {{{1}, 1.0}};
    else if (order == 2) terms = {{{2}, 1.0}, {{1, 1}, 0.5}};
    else terms = {{{3}, 1.0}, {{1, 2}, 1.0}, {{1, 1, 1}, 1.0 / 6.0}};
    PropertyChain result;
    result.value.assign(observables.size(), std::nullopt);
    result.coefficient.assign(observables.size(), 0.0);
    result.status.assign(observables.size(), PropertyActionStatus::Available);
    for (const auto& term : terms) {
        std::vector<DensityActionDirection> directions;
        for (const std::size_t k : term.orders) directions.push_back(direction(k));
        const auto batch = model.density_actions(series.temperature[0], map.density[0], central_composition,
                                                 observables, directions);
        result.diagnostics = batch.diagnostics;
        result.maximum_energy_order = std::max(result.maximum_energy_order, batch.maximum_energy_order);
        for (std::size_t i = 0; i < observables.size(); ++i) {
            const auto& entry = batch.results[i];
            if (entry.value.has_value()) result.value[i] = entry.value;
            if (entry.status != PropertyActionStatus::Available && result.status[i] == PropertyActionStatus::Available)
                result.status[i] = entry.status;
            if (result.coefficient[i].has_value() && entry.action.has_value())
                *result.coefficient[i] += term.weight * *entry.action;
            else result.coefficient[i].reset();
        }
    }
    return result;
}

}  // namespace

std::string solved_state_observable_id(const SolvedStateObservable& observable) {
    const std::string phase = "phase[" + std::to_string(observable.phase) + "]";
    const std::string component = "[" + std::to_string(observable.component) + "]";
    switch (observable.kind) {
    case SolvedStateObservableKind::Temperature: return "temperature";
    case SolvedStateObservableKind::Pressure: return "pressure";
    case SolvedStateObservableKind::PhaseAmount: return phase + ".amount";
    case SolvedStateObservableKind::PhaseComponentAmount: return phase + ".amount" + component;
    case SolvedStateObservableKind::PhaseMoleFraction: return phase + ".mole_fraction" + component;
    case SolvedStateObservableKind::PhaseVolume: return phase + ".volume";
    case SolvedStateObservableKind::PhaseMolarDensity: return phase + ".molar_density";
    case SolvedStateObservableKind::PhasePackingFraction: return phase + ".packing_fraction";
    case SolvedStateObservableKind::PhasePressure: return phase + ".pressure";
    case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:
        return phase + ".residual_chemical_potential_over_rt" + component;
    case SolvedStateObservableKind::PhaseLogActivity: return phase + ".log_activity" + component;
    case SolvedStateObservableKind::PhaseProperty: return phase + ".property/" + property_name(observable.property);
    }
    return "unknown";
}

SolvedStateActionBatch solved_state_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central, const SolvedStateActionRequest& request
) {
    core::Scope scope("solved_state_actions");
    validate_request(compiled, central, request);
    const auto& problem = compiled.problem;
    const auto& model = *compiled.model;
    const std::size_t components = model.parameters().component_count();
    const std::size_t actives = model.active_parameters().size();
    const std::size_t width = 2 + components + actives;
    const std::size_t order = request.directions.size();
    SolvedStateActionBatch batch;
    batch.diagnostics.equation_count = compiled.rows.size();
    batch.diagnostics.state_dimension = compiled.unknown_count();
    for (const auto& observable : request.observables) {
        SolvedStateActionResult item;
        item.observable = observable;
        item.unit = observable_unit(observable);
        item.order = order;
        batch.results.push_back(std::move(item));
    }
    const auto finish = [&](SolvedStateActionStatus status) {
        for (auto& item : batch.results)
            if (item.status == SolvedStateActionStatus::Available) item.status = status;
        return batch;
    };
    if (!central.success || !central.has_state()
        || !std::all_of(central.coordinates.begin(), central.coordinates.end(), [](double v) { return std::isfinite(v); })) {
        for (auto& item : batch.results) item.status = SolvedStateActionStatus::CentralSolveUnavailable;
        return batch;
    }
    for (auto& item : batch.results) item.status = SolvedStateActionStatus::Available;
    Center center;
    center.prepared = compiled.prepare(central.coordinates);
    // The gate quantity is the classified residual the driver certifies (owner decisions 2 and
    // 6, 2026-09-16); the unscaled rows stay in the central result's diagnostics.
    const auto evaluation = compiled.evaluate_at(center.prepared);
    double classified_residual = 0.0;
    for (const double value : compiled.classified_residuals(center.prepared, evaluation.values))
        classified_residual = std::isfinite(value) ? std::max(classified_residual, std::abs(value)) : std::numeric_limits<double>::infinity();
    if (std::isfinite(classified_residual)) batch.diagnostics.central_residual = classified_residual;
    if (evaluation.floor > 0.0) batch.diagnostics.central_floor = evaluation.floor;
    scope.probe("central_residual", classified_residual);
    // Controls, their scales and the direction magnitudes.
    std::vector<double> scale(width, 1.0);
    scale[0] = std::max(center.prepared.temperature, 1.0);
    scale[1] = std::max(center.prepared.pressure, 1.0);
    for (std::size_t i = 0; i < components; ++i) scale[2 + i] = std::max(std::abs(compiled.n_feed[i]), 1.0);
    std::vector<std::vector<double>> controls;
    std::vector<double> magnitude;
    for (const auto& direction : request.directions) {
        std::vector<double> control = {direction.temperature_k, direction.pressure_pa};
        control.insert(control.end(), direction.feed_amounts_mol.begin(), direction.feed_amounts_mol.end());
        control.insert(control.end(), direction.active_parameters.begin(), direction.active_parameters.end());
        double largest = 0.0;
        for (std::size_t i = 0; i < width; ++i) largest = std::max(largest, std::abs(control[i]) / scale[i]);
        controls.push_back(std::move(control));
        magnitude.push_back(largest);
    }
    // Central values through the physical map; fixed inputs keep their identity actions.
    const std::vector<double> no_control(width, 0.0);
    const auto base = compiled.row_series(center.prepared, {center.prepared.point}, no_control, 0);
    std::vector<PhaseMap> base_maps;
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) base_maps.push_back(phase_map(base, phase, components));
    for (auto& item : batch.results) {
        if (item.observable.kind == SolvedStateObservableKind::PhaseProperty) continue;
        const auto series = observable_series(item.observable, base, base_maps);
        if (!series.empty() && std::isfinite(series.front())) item.value = series.front();
        if (!state_dependent(item.observable, problem)) {
            const std::size_t index = item.observable.kind == SolvedStateObservableKind::Temperature ? 0 : 1;
            item.action = order == 1 ? controls[0][index] : 0.0;
        }
    }
    // Property central values and dependency statuses through the P4 owner at a zero direction.
    std::vector<std::vector<std::size_t>> property_rows(problem.phases.size());
    for (std::size_t i = 0; i < batch.results.size(); ++i)
        if (batch.results[i].observable.kind == SolvedStateObservableKind::PhaseProperty)
            property_rows[batch.results[i].observable.phase].push_back(i);
    std::vector<std::vector<PropertyObservable>> property_lists(problem.phases.size());
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        for (const std::size_t row : property_rows[phase]) property_lists[phase].push_back(batch.results[row].observable.property);
        if (property_rows[phase].empty()) continue;
        const auto chain = property_chain(model, base, base_maps[phase], phase, property_lists[phase], 0);
        batch.diagnostics.maximum_energy_order = std::max(batch.diagnostics.maximum_energy_order, chain.maximum_energy_order);
        for (std::size_t i = 0; i < property_rows[phase].size(); ++i) {
            auto& item = batch.results[property_rows[phase][i]];
            item.value = chain.value[i];
            item.property_diagnostics = chain.diagnostics;
            if (chain.status[i] != PropertyActionStatus::Available) {
                item.status = solved_state_status(chain.status[i]);
                item.diagnostic_message = std::string("phase_property dependency: ") + solved_state_status_name(item.status);
            }
        }
    }
    // Admission of the central state, in the owner's precedence.
    SolvedStateActionStatus gate = SolvedStateActionStatus::Available;
    if (!(classified_residual <= std::min(problem.row_tolerance, policy::EQUILIBRIUM_ROW_RESIDUAL.value))) {
        gate = SolvedStateActionStatus::DefiningEquationsNotSatisfied;
    } else {
        const auto topology = check_topology(compiled, central.coordinates);
        if (topology.event == TopologyEvent::PhaseDisappeared) gate = SolvedStateActionStatus::PhaseBoundaryNotDifferentiable;
        else if (topology.event == TopologyEvent::PhasesCoalesced) gate = SolvedStateActionStatus::NonuniqueCoalescedPhases;
        for (std::size_t phase = 0; phase < problem.phases.size() && gate == SolvedStateActionStatus::Available; ++phase)
            for (const std::size_t species : problem.phases[phase].support)
                if (center.prepared.composition[phase][species] <= policy::STABILITY_COMPOSITION_FLOOR.value)
                    gate = SolvedStateActionStatus::PhaseBoundaryNotDifferentiable;
    }
    if (gate == SolvedStateActionStatus::Available) {
        const auto charges = species_z(model.parameters());
        bool support = false, nonzero_z = false;
        for (const std::size_t species : compiled.removed_species)
            for (const auto& control : controls)
                if (control[2 + species] != 0.0) (charges[species] == 0 ? support : nonzero_z) = true;
        if (nonzero_z) gate = SolvedStateActionStatus::IonicInsertion;
        else if (support) gate = SolvedStateActionStatus::SupportChange;
    }
    if (gate == SolvedStateActionStatus::Available) {
        const std::size_t n = central.coordinates.size();
        const auto jacobian = compiled.jacobian_at(center.prepared);
        center.jacobian = dense_matrix(jacobian, compiled.rows.size(), n);
        if (!square_and_numeric(center.jacobian)) {
            gate = SolvedStateActionStatus::NonfiniteJacobian;
        } else {
            center.lu = center.jacobian.fullPivLu();
            batch.diagnostics.state_factorizations = 1;
            batch.diagnostics.rank = static_cast<std::size_t>(center.lu.rank());
            const double condition = condition_estimate(center.jacobian);
            if (std::isfinite(condition)) batch.diagnostics.condition_estimate = condition;
            if (*batch.diagnostics.rank < n) gate = SolvedStateActionStatus::RankDeficient;
        }
    }
    const bool zero = std::any_of(magnitude.begin(), magnitude.end(), [](double m) { return m == 0.0; });
    const bool thermal = std::any_of(controls.begin(), controls.end(), [](const auto& c) { return c[0] != 0.0; });
    if (gate == SolvedStateActionStatus::Available && !zero && thermal && thermal_reference_missing(compiled, center.prepared.temperature))
        gate = SolvedStateActionStatus::MissingThermochemistry;
    // A central-state gate precedes the selected dependency: it overrides every
    // state-dependent row's status, keeping the mapped value and property diagnostics.
    const auto apply_state_gate = [&](SolvedStateActionStatus status) {
        for (auto& item : batch.results) {
            if (!state_dependent(item.observable, problem)) continue;
            item.status = status;
            item.action.reset();
            item.diagnostic_message.reset();
        }
        return batch;
    };
    if (gate != SolvedStateActionStatus::Available) return apply_state_gate(gate);
    if (zero) {
        for (auto& item : batch.results)
            if (item.status == SolvedStateActionStatus::Available && item.value.has_value()) item.action = 0.0;
        return finish(SolvedStateActionStatus::Available);
    }
    // Polarised diagonal paths: half the sign patterns with the first sign fixed positive.
    std::vector<double> accumulated(batch.results.size(), 0.0);
    std::vector<bool> usable(batch.results.size(), true);
    const std::size_t paths = std::size_t(1) << (order - 1);
    for (std::size_t mask = 0; mask < paths; ++mask) {
        std::vector<double> control(width, 0.0);
        double weight = 1.0;
        for (std::size_t j = 0; j < order; ++j) {
            const double sign = j != 0 && (mask & (std::size_t(1) << (j - 1))) ? -1.0 : 1.0;
            weight *= sign;
            for (std::size_t i = 0; i < width; ++i) control[i] += sign * controls[j][i] / magnitude[j];
        }
        core::Scope path_scope("solved_state_action_path");
        const auto path = follow(compiled, center, control, order, batch.diagnostics);
        if (path.status != SolvedStateActionStatus::Available) return apply_state_gate(path.status);
        std::vector<PhaseMap> maps;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) maps.push_back(phase_map(path.series, phase, components));
        for (std::size_t i = 0; i < batch.results.size(); ++i) {
            auto& item = batch.results[i];
            if (item.observable.kind == SolvedStateObservableKind::PhaseProperty || !state_dependent(item.observable, problem)) continue;
            const auto series = observable_series(item.observable, path.series, maps);
            accumulated[i] += weight * series[order];
        }
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            if (property_rows[phase].empty()) continue;
            const auto chain = property_chain(model, path.series, maps[phase], phase, property_lists[phase], order);
            batch.diagnostics.maximum_energy_order = std::max(batch.diagnostics.maximum_energy_order, chain.maximum_energy_order);
            for (std::size_t i = 0; i < property_rows[phase].size(); ++i) {
                auto& item = batch.results[property_rows[phase][i]];
                item.property_diagnostics = chain.diagnostics;
                if (chain.status[i] != PropertyActionStatus::Available && item.status == SolvedStateActionStatus::Available) {
                    item.status = solved_state_status(chain.status[i]);
                    item.diagnostic_message = std::string("phase_property dependency: ") + solved_state_status_name(item.status);
                }
                if (chain.coefficient[i].has_value()) accumulated[property_rows[phase][i]] += weight * *chain.coefficient[i];
                else usable[property_rows[phase][i]] = false;
            }
        }
    }
    double product_of_magnitudes = 1.0;
    for (const double m : magnitude) product_of_magnitudes *= m;
    for (std::size_t i = 0; i < batch.results.size(); ++i) {
        auto& item = batch.results[i];
        if (!state_dependent(item.observable, problem) || item.status != SolvedStateActionStatus::Available) continue;
        const double action = accumulated[i] / static_cast<double>(paths) * product_of_magnitudes;
        if (!usable[i] || !std::isfinite(action) || !item.value.has_value()) {
            item.status = SolvedStateActionStatus::Nonfinite;
            continue;
        }
        item.action = action;
    }
    scope.probe("results", static_cast<double>(batch.results.size()));
    return batch;
}

}  // namespace epcsaft::equilibrium
