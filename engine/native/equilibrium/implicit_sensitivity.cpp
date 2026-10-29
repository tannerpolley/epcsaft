#include "equilibrium/gibbs_minimizer.hpp"

#include "core/trace.hpp"
#include "core/linear_error.hpp"
#include "core/numerical_policy.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <numeric>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

using RowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Eigen::MatrixXd dense_matrix(const std::vector<double>& values, std::size_t rows, std::size_t columns) {
    if (values.size() != rows * columns) throw std::invalid_argument("matrix size");
    return Eigen::Map<const RowMajor>(values.data(), static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(columns));
}

double condition_estimate(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() == 0 || matrix.cols() == 0) return 0.0;
    const auto singular = matrix.jacobiSvd().singularValues();
    if (singular.size() == 0 || singular(singular.size() - 1) <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return singular(0) / singular(singular.size() - 1);
}

bool valid_point(const std::vector<double>& point) {
    return !point.empty() && std::all_of(
        point.begin(), point.end(), [](double value) { return std::isfinite(value); });
}

using core::backward_error;
using core::compensated_defect;

bool square_and_numeric(const Eigen::MatrixXd& matrix) { return matrix.rows() == matrix.cols() && matrix.allFinite(); }

}  // namespace

SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point
) {
    core::Scope scope("inspect_implicit_sensitivity");
    const std::size_t rows = compiled.rows.size();
    const std::size_t dimension = point.size();
    if (!valid_point(point)) throw std::invalid_argument("invalid sensitivity point");
    const auto jacobian = compiled.jacobian(point);
    if (jacobian.size() != rows * dimension) {
        throw std::invalid_argument("equilibrium Jacobian shape");
    }
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    const auto lu = matrix.fullPivLu();
    SensitivityInfo result;
    result.equations = rows;
    result.rank = static_cast<std::size_t>(lu.rank());
    result.threshold = lu.threshold();
    result.condition_estimate = condition_estimate(matrix);
    result.nonsingular = square_and_numeric(matrix) && result.rank == dimension;
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    scope.probe("condition", result.condition_estimate);
    scope.probe("jacobian_entries", static_cast<double>(jacobian.size()));
    scope.probe("nonsingular", result.nonsingular ? 1.0 : 0.0);
    scope.probe("rank", static_cast<double>(result.rank));
    return result;
}

std::vector<double> solve_implicit(
    const CompiledSystem& compiled,
    const PhasePoint& prepared,
    const std::vector<double>& column, std::optional<core::Convergence>* failure
) {
    core::Scope scope("implicit_sensitivity");
    const std::size_t rows = compiled.rows.size();
    const std::size_t dimension = prepared.point.size();
    if (!valid_point(prepared.point)) throw std::invalid_argument("invalid sensitivity point");
    if (!all_numeric(column)) throw std::invalid_argument("invalid parameter derivative");
    const auto jacobian = compiled.jacobian_at(prepared);
    const auto refuse = [&](const char* cause, core::Stop stop, const auto* factor, std::optional<double> backward = std::nullopt) {
        if (failure) { core::Convergence r; r.cause = cause;
            if (factor) r.rank = core::numerical_rank(*factor, "implicit defining rows");
            r.residual = backward.value_or(NAN); r.backward_error = backward; r.stop = stop;
            if (backward) r.criterion = &policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR;
            core::retain_first(*failure, r);
        }
        return std::vector<double>{};
    };
    if (jacobian.size() != rows * dimension || rows != dimension) return refuse("implicit Jacobian shape", core::Stop::infeasible, static_cast<const Eigen::FullPivLU<Eigen::MatrixXd>*>(nullptr));
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    const auto lu = matrix.fullPivLu();
    if (!square_and_numeric(matrix) || static_cast<std::size_t>(lu.rank()) != dimension)
        return refuse(matrix.allFinite() ? "implicit Jacobian rank" : "implicit nonfinite Jacobian", matrix.allFinite() ? core::Stop::singular : core::Stop::nonfinite, &lu);
    // One column, kept a matrix so that Eigen takes the same solve path.
    const Eigen::MatrixXd right = -Eigen::Map<const Eigen::VectorXd>(column.data(), static_cast<Eigen::Index>(column.size()));
    const Eigen::MatrixXd result = lu.solve(right);
    const double error = std::max(0.0, backward_error((matrix * result.col(0) - right.col(0)).cwiseAbs().maxCoeff(),
        matrix.cwiseAbs().rowwise().sum().maxCoeff(), result.cwiseAbs().maxCoeff(), right.cwiseAbs().maxCoeff()));
    if (!result.allFinite() || !(error <= policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR.value))
        return refuse(result.allFinite() ? "implicit backward error" : "implicit nonfinite response", result.allFinite() ? core::Stop::stalled : core::Stop::nonfinite, &lu,
            result.allFinite() ? std::optional<double>(error) : std::nullopt);
    scope.probe("max_solution", result.cwiseAbs().maxCoeff());
    scope.probe("backward_error", error);
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    return {result.data(), result.data() + rows};
}

std::vector<double> sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& lever, std::optional<core::Convergence>* failure
) {
    core::Scope scope("sensitivity");
    if (lever.size() != compiled.control_width() && lever.size() != compiled.extended_control_width())
        throw std::invalid_argument("lever must contain T, P, feed, and active parameter weights");
    // F_theta * lever: the first row coefficient along the lever at a fixed state, the actions' directional
    // series, which carries source-reference T/P chains (the dense lever Jacobian does not).
    const auto prepared = compiled.prepare(point);
    const std::vector<std::vector<double>> state = {prepared.point, std::vector<double>(point.size(), 0.0)};
    const auto series = compiled.row_series(prepared, state, lever, 1);
    return solve_implicit(compiled, prepared, series.rows[1], failure);
}

// ------------------------------------------------------------- solved-state actions

namespace {

using Series = std::vector<double>;

template <typename Base>
std::vector<Base> product(const std::vector<Base>& left, const std::vector<Base>& right) {
    std::vector<Base> result(left.size(), 0.0);
    for (std::size_t i = 0; i < result.size(); ++i)
        for (std::size_t j = 0; j <= i; ++j) result[i] += left[j] * right[i - j];
    return result;
}

template <typename Base>
std::vector<Base> inverse(const std::vector<Base>& value) {
    if (value.empty() || !std::isfinite(value.front()) || value.front() == 0.0)
        throw std::domain_error("series inverse has a zero or nonfinite center");
    std::vector<Base> result(value.size(), 0.0);
    result.front() = 1.0 / value.front();
    for (std::size_t i = 1; i < value.size(); ++i) {
        for (std::size_t j = 1; j <= i; ++j) result[i] -= value[j] * result[i - j];
        result[i] /= value.front();
    }
    return result;
}

template <typename Base>
std::vector<Base> combine(const std::vector<Base>& left, const std::vector<Base>& right, double sign) {
    std::vector<Base> result(left.size());
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

bool is_solid_kind(SolvedStateObservableKind kind) {
    return kind == SolvedStateObservableKind::SolidAmount || kind == SolvedStateObservableKind::SolidSaturationIndex;
}

bool is_phase_kind(SolvedStateObservableKind kind) {
    return kind != SolvedStateObservableKind::Temperature && kind != SolvedStateObservableKind::Pressure
        && kind != SolvedStateObservableKind::Loading && kind != SolvedStateObservableKind::SystemEnthalpy
        && kind != SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity && !is_solid_kind(kind);
}

bool has_component(SolvedStateObservableKind kind) {
    return kind == SolvedStateObservableKind::PhaseComponentAmount
        || kind == SolvedStateObservableKind::PhaseMoleFraction
        || kind == SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt
        || kind == SolvedStateObservableKind::PhaseLogActivity;
}

bool same_observable(const SolvedStateObservable& a, const SolvedStateObservable& b) {
    if (a.kind != b.kind) return false;
    if ((is_phase_kind(a.kind) || is_solid_kind(a.kind)) && a.phase != b.phase) return false;
    if (has_component(a.kind) && a.component != b.component) return false;
    return a.kind != SolvedStateObservableKind::PhaseProperty || a.property == b.property;
}

const char* observable_unit(const SolvedStateObservable& observable) {
    switch (observable.kind) {
    case SolvedStateObservableKind::Loading: return "mol/mol";
    case SolvedStateObservableKind::SystemEnthalpy: return "J";
    case SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity: return "J/(mol-feed K)";
    case SolvedStateObservableKind::Temperature: return "K";
    case SolvedStateObservableKind::Pressure: case SolvedStateObservableKind::PhasePressure: return "Pa";
    case SolvedStateObservableKind::PhaseAmount: case SolvedStateObservableKind::PhaseComponentAmount:
    case SolvedStateObservableKind::SolidAmount: return "mol";
    case SolvedStateObservableKind::PhaseVolume: return "m^3";
    case SolvedStateObservableKind::PhaseGibbsMolarVolume:
        return "m^3/mol";
    case SolvedStateObservableKind::PhaseMolarDensity: return "mol/m^3";
    case SolvedStateObservableKind::PhaseProperty:
        return property_kind(observable.property) < 2 ? "J/mol" : "J/(mol K)";
    default: return "1";
    }
}

std::string property_name(PropertyObservable observable) {
    static const char* const family[] = {"Ideal", "Residual", "Total"};
    static const char* const kind[] = {"Enthalpy", "InternalEnergy", "Entropy", "IsobaricHeatCapacity", "IsochoricHeatCapacity"};
    if (static_cast<std::size_t>(observable) >= CALORIC_PROPERTY_COUNT) return "unknown";
    return std::string(family[property_family(observable)]) + kind[property_kind(observable)];
}

// Every observable except a fixed T or P input depends on the solved state.
bool state_dependent(const SolvedStateObservable& observable, const EquilibriumProblem& problem) {
    if (observable.kind == SolvedStateObservableKind::Loading) return problem.loading.has_value();
    if (observable.kind == SolvedStateObservableKind::Temperature) return problem.temperature.free;
    if (observable.kind == SolvedStateObservableKind::Pressure) return problem.pressure.free;
    return true;
}

void validate_request(const CompiledSystem& compiled, const EquilibriumResult& central,
                      const SolvedStateActionRequest& request) {
    if (compiled.ideal) throw std::invalid_argument("ideal solved-state actions unavailable");
    const auto& problem = compiled.problem;
    if (central.has_state() && (!central.declaration ||
        !(central.declaration->solver_options == problem.solver_options)))
        throw std::invalid_argument("central solver options do not match the compiled system");
    const std::size_t components = compiled.model->parameters().component_count();
    const std::size_t actives = compiled.model->active_parameters().size();
    if (request.observables.empty()) throw std::invalid_argument("no solved-state observables");
    for (std::size_t i = 0; i < request.observables.size(); ++i) {
        const auto& item = request.observables[i];
        if (is_phase_kind(item.kind) && item.phase >= problem.phases.size())
            throw std::invalid_argument("solved-state observable phase ordinal");
        if (is_solid_kind(item.kind) && item.phase >= problem.solids.size())
            throw std::invalid_argument("solved-state observable solid ordinal");
        if (has_component(item.kind) && item.component >= components)
            throw std::invalid_argument("solved-state observable component ordinal");
        if (item.kind == SolvedStateObservableKind::PhaseProperty && static_cast<std::size_t>(item.property) >= CALORIC_PROPERTY_COUNT)
            throw std::invalid_argument("solved-state property observable");
        if (item.kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity
            && (problem.temperature.free || problem.pressure.free || problem.temperature.absent
                || problem.pressure.absent || request.directions.size() != 1))
            throw std::invalid_argument("apparent Cp requires driven T/P and first actions only");
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
    for (const auto& item : request.observables)
        if (item.kind == SolvedStateObservableKind::Loading && !problem.loading)
            throw std::invalid_argument("loading action requires a declared affine coordinate");
    for (const auto& direction : request.directions) {
        if (!std::isfinite(direction.temperature_k) || !std::isfinite(direction.pressure_pa)
            || direction.feed_amounts_mol.size() != components || direction.active_parameters.size() != actives
            || !all_numeric(direction.feed_amounts_mol) || !all_numeric(direction.active_parameters))
            throw std::invalid_argument("solved-state action direction");
        if (!direction.reaction_coefficients.empty() &&
            direction.reaction_coefficients.size() != problem.reactions.size())
            throw std::invalid_argument("reaction coefficient direction shape");
        for (const auto &coefficients : direction.reaction_coefficients)
            for (double v : coefficients)
                if (!std::isfinite(v))
                    throw std::invalid_argument("nonfinite reaction coefficient direction");
                else if (v != 0.0 && !problem.standard_gibbs_j_per_mol.empty())
                    throw std::invalid_argument("reaction direction conflicts with supplied standard Gibbs values");
        if (!std::isfinite(direction.enthalpy_j) || !std::isfinite(direction.entropy_j_per_k)
            || (!problem.enthalpy_j && direction.enthalpy_j != 0.0)
            || (!problem.entropy_j_per_k && direction.entropy_j_per_k != 0.0))
            throw std::invalid_argument("caloric control without its active condition or nonfinite");
        if (problem.temperature.free && direction.temperature_k != 0.0)
            throw std::invalid_argument("free temperature is a solved output, not an external control");
        if (problem.pressure.free && direction.pressure_pa != 0.0)
            throw std::invalid_argument("free pressure is a solved output, not an external control");
        if (!problem.solids.empty() && request.directions.size() > 1 && (direction.temperature_k != 0.0 || direction.pressure_pa != 0.0))
            throw std::invalid_argument("solids admit T and P actions through order one");
    }
    // A referenced problem: compilation resolves a diagnostic for every source-referenced reaction.
    const bool referenced = !compiled.reference_diagnostics.empty();
    std::size_t parameter_directions = 0;
    for (const auto &d : request.directions) {
        bool parameter = false;
        for (double v : d.active_parameters) parameter |= v != 0.0;
        for (const auto &field : d.reaction_coefficients)
            for (double v : field)
                parameter |= v != 0.0;
        parameter_directions += std::size_t(parameter);
        if (referenced && request.directions.size() == 3 &&
            (parameter || d.temperature_k != 0.0 || d.pressure_pa != 0.0 || problem.temperature.free ||
             problem.pressure.free))
            throw std::invalid_argument("a referenced third-order action admits only feed directions at fixed T and P");
    }
    if (referenced && request.directions.size() == 2 && parameter_directions > 1)
        throw std::invalid_argument("a referenced second-order action admits at most one parameter or reaction-law direction");
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

// Ideal and total caloric properties carry the standard-state (record) terms; residual ones do not.
bool caloric_property(PropertyObservable property) {
    return static_cast<std::size_t>(property) < CALORIC_PROPERTY_COUNT && property_family(property) != 1;
}

struct Center final {
    PhasePoint prepared;
    Eigen::MatrixXd jacobian;  // equilibrated R J C with power-of-two diagonals row, column
    Eigen::VectorXd row, column;
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
};

struct Path final {
    RowSeries series;
    SolvedStateActionStatus status = SolvedStateActionStatus::Available;
    std::optional<std::string> diagnostic_message;
};

bool is_reference_unavailable(const std::exception& error) {
    const std::string message = error.what();
    return message.rfind("reference_unavailable:", 0) == 0;
}

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
    std::vector<double> right_norm(order + 1, 0.0), solution_norm(order + 1, 0.0);
    // The series through order k, or false with the path refused when a reference is unavailable.
    const auto rows_through = [&](const std::vector<std::vector<double>>& prefix, std::size_t k, RowSeries& series) {
        try {
            series = compiled.row_series(center.prepared, prefix, control, k);
            return true;
        } catch (const std::invalid_argument& error) {
            if (!is_reference_unavailable(error)) throw;
            path.status = SolvedStateActionStatus::ReferenceUnavailable;
            path.diagnostic_message = error.what();
            return false;
        }
    };
    for (std::size_t k = 1; k <= order; ++k) {
        RowSeries series;
        if (!rows_through({state.begin(), state.begin() + static_cast<std::ptrdiff_t>(k + 1)}, k, series)) return path;
        diagnostics.maximum_energy_order = std::max(diagnostics.maximum_energy_order, series.maximum_energy_order);
        if (series.energy != PropertyActionStatus::Available) { path.status = energy_status(series.energy); return path; }
        Eigen::VectorXd right(static_cast<Eigen::Index>(n));
        for (std::size_t row = 0; row < n; ++row) right(static_cast<Eigen::Index>(row)) = -center.row(static_cast<Eigen::Index>(row)) * series.rows[k][row];
        if (!right.allFinite()) { path.status = SolvedStateActionStatus::NonfiniteJacobian; return path; }
        right_norm[k] = right.cwiseAbs().maxCoeff();
        // Solve R J C z = R b, then refine once with the double-double defect.
        Eigen::VectorXd scaled = center.lu.solve(right);
        scaled += center.lu.solve(compensated_defect(center.jacobian, scaled, right));
        solution_norm[k] = scaled.cwiseAbs().maxCoeff();
        const Eigen::VectorXd solution = center.column.cwiseProduct(scaled);
        ++diagnostics.linear_solves;
        if (!solution.allFinite()) { path.status = SolvedStateActionStatus::NonfiniteJacobian; return path; }
        for (std::size_t row = 0; row < n; ++row) state[k][row] = solution(static_cast<Eigen::Index>(row));
    }
    if (!rows_through(state, order, path.series)) return path;
    if (compiled.has_caloric_condition() && std::any_of(path.series.temperature.begin() + 1,
            path.series.temperature.end(), [](double v) { return v != 0.0; })) {
        const auto calorics = reference_calorics(compiled, center.prepared.temperature, center.prepared.pressure, false, false, false, 1);
        if (calorics.thermal_boundary) {
            path.status = SolvedStateActionStatus::ThermalBoundary;
            path.diagnostic_message = "caloric condition derivative crosses an interval endpoint";
            return path;
        }
    }
    diagnostics.maximum_energy_order = std::max(diagnostics.maximum_energy_order, path.series.maximum_energy_order);
    if (path.series.energy != PropertyActionStatus::Available) { path.status = energy_status(path.series.energy); return path; }
    for (std::size_t k = 1; k <= order; ++k) {
        double residual = 0.0;
        for (std::size_t row = 0; row < n; ++row) residual = std::max(residual, std::abs(center.row(static_cast<Eigen::Index>(row)) * path.series.rows[k][row]));
        residual = backward_error(residual, center.jacobian.cwiseAbs().rowwise().sum().maxCoeff(), solution_norm[k], right_norm[k]);
        if (!std::isfinite(residual)) { path.status = SolvedStateActionStatus::NonfiniteJacobian; return path; }
        diagnostics.maximum_linear_residual = std::max(diagnostics.maximum_linear_residual.value_or(0.0), residual);
        // The first refused coefficient names the cause; 2 cond eta is the first-order forward bound.
        if (path.status != SolvedStateActionStatus::Available) continue;
        if (residual > policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR.value) path.status = SolvedStateActionStatus::LinearSolveResidual;
        else if (!(2.0 * *diagnostics.condition_estimate * residual <= policy::EQUILIBRIUM_ACTION_FORWARD_ERROR.value))
            path.status = SolvedStateActionStatus::IllConditioned;
    }
    return path;
}

// Physical map of one phase along a path: normalised composition, density and pressure.
template <typename Base>
struct PhaseMapValues final {
    std::vector<Base> amount, volume, density, packing, pressure;
    std::vector<std::vector<Base>> mole_fraction, potential;
};

using PhaseMap = PhaseMapValues<double>;
template <typename Base>
PhaseMapValues<Base> phase_map(const HelmholtzModel& model, const PhaseSeriesValues<Base>& item,
    const std::vector<Base>& temperature, std::size_t components) {
    const std::size_t order = item.point.size();
    PhaseMapValues<Base> result;
    result.amount = item.amount;
    result.packing = item.packing;
    std::vector<Base> volume(order), total(order, 0.0), potential_v(order);
    for (std::size_t k = 0; k < order; ++k) {
        volume[k] = item.point[k][1];
        potential_v[k] = item.gradient[k][1];
        for (std::size_t i = 0; i < components; ++i) total[k] += item.point[k][2 + i];
    }
    result.volume = volume;
    result.density = model.parameters().ideal_gas ? product(total, inverse(volume)) : inverse(volume);
    const auto inverse_total = inverse(total);
    for (std::size_t i = 0; i < components; ++i) {
        std::vector<Base> x(order), mu(order);
        for (std::size_t k = 0; k < order; ++k) { x[k] = item.point[k][2 + i]; mu[k] = item.gradient[k][2 + i]; }
        result.mole_fraction.push_back(product(x, inverse_total));
        result.potential.push_back(mu);
    }
    std::vector<Base> pressure = product(temperature, combine(result.density, potential_v, -1.0));
    for (Base& value : pressure) value *= GAS_CONSTANT_J_PER_MOL_K;
    result.pressure = pressure;
    return result;
}

Series system_enthalpy_series(const HelmholtzModel &, const RowSeries &, const std::vector<PhaseMap> &,
                              const ReferenceCalorics *, SolvedStateActionResult *, SolvedStateActionDiagnostics *);

Series standard_series(const ReferenceCalorics &, const RowSeries &, std::size_t, std::size_t, std::size_t,
                       std::size_t, std::size_t);

Series observable_series(const HelmholtzModel &model, const SolvedStateObservable &observable, const RowSeries &series,
                         const std::vector<PhaseMap> &maps, const ReferenceCalorics *calorics,
                         SolvedStateActionResult *result, SolvedStateActionDiagnostics *diagnostics) {
    switch (observable.kind) {
    case SolvedStateObservableKind::SystemEnthalpy:
        return system_enthalpy_series(model, series, maps, calorics, result, diagnostics);
    case SolvedStateObservableKind::Temperature: return series.temperature;
    case SolvedStateObservableKind::Loading: return series.loading;
    case SolvedStateObservableKind::Pressure: return series.pressure;
    case SolvedStateObservableKind::PhaseAmount: return maps[observable.phase].amount;
    case SolvedStateObservableKind::PhaseComponentAmount:
        return product(maps[observable.phase].amount, maps[observable.phase].mole_fraction[observable.component]);
    case SolvedStateObservableKind::PhaseMoleFraction: return maps[observable.phase].mole_fraction[observable.component];
    case SolvedStateObservableKind::PhaseVolume: return product(maps[observable.phase].amount, maps[observable.phase].volume);
    case SolvedStateObservableKind::PhaseGibbsMolarVolume: {
        auto volume = maps[observable.phase].volume;
        if (calorics)
            for (std::size_t i = 0; i < maps[observable.phase].mole_fraction.size(); ++i)
                volume = combine(
                    volume,
                    product(maps[observable.phase].mole_fraction[i],
                            standard_series(*calorics, series, observable.phase, i, 0, 1, volume.size() - 1)),
                    1.0);
        return volume;
    }
    case SolvedStateObservableKind::PhaseMolarDensity: return maps[observable.phase].density;
    case SolvedStateObservableKind::PhasePackingFraction: return maps[observable.phase].packing;
    case SolvedStateObservableKind::PhasePressure: return maps[observable.phase].pressure;
    case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:
        return maps[observable.phase].potential[observable.component];
    case SolvedStateObservableKind::PhaseLogActivity:
        return combine(logarithm(product(maps[observable.phase].mole_fraction[observable.component],
                                         maps[observable.phase].density)),
                       maps[observable.phase].potential[observable.component], 1.0);
    case SolvedStateObservableKind::SolidAmount: return series.solid_amount[observable.phase];
    case SolvedStateObservableKind::SolidSaturationIndex: {  // EqID solid_saturation_index
        auto index = series.solid_affinity[observable.phase];
        for (double& value : index) value /= std::log(10.0);
        return index;
    }
    default: return {};
    }
}

// The P4 chain of one phase along a curved path: [t^1] g = Dg[z1]; [t^2] g = Dg[z2] +
// D2g[z1,z1]/2; [t^3] g = Dg[z3] + D2g[z1,z2] + D3g[z1,z1,z1]/6 in the (T, rho, x, theta) chart.
struct PropertyChain final {
    std::vector<std::optional<double>> value, coefficient;
    std::vector<Series> series;
    std::vector<SolvedStateActionStatus> status;
    std::vector<std::string> message;
    PropertyActionDiagnostics diagnostics;
    std::size_t maximum_energy_order = 0;
};
struct CaloricUnavailable final : std::runtime_error {
    SolvedStateActionStatus status;
    CaloricUnavailable(SolvedStateActionStatus s, const std::string& message)
        : std::runtime_error(message), status(s) {}
};

// Compose a raw anchored Gibbs derivative on the solved T/P/parameter path.
Series standard_series(const ReferenceCalorics &calorics, const RowSeries &path, std::size_t phase,
                       std::size_t species, std::size_t a, std::size_t b, std::size_t order) {
    if (calorics.status != SolvedStateActionStatus::Available)
        throw CaloricUnavailable(calorics.status, calorics.message);
    if (calorics.thermal_boundary && std::any_of(path.temperature.begin() + 1, path.temperature.end(),
                                               [](double v) { return v != 0.0; }))
        throw CaloricUnavailable(SolvedStateActionStatus::ThermalBoundary, "record derivative crosses an interval endpoint");
    Series result(order + 1, 0.0), dt(order + 1, 0.0), dp(order + 1, 0.0), one(order + 1, 0.0);
    one[0] = 1.0;
    for (std::size_t k = 1; k <= order; ++k) {
        dt[k] = path.temperature[k];
        dp[k] = path.pressure[k];
    }
    for (const auto &entry : calorics.gibbs) {
        if (entry.temperature_order < a || entry.pressure_order < b)
            continue;
        const auto ta = entry.temperature_order - a, pb = entry.pressure_order - b;
        if (ta + pb + std::size_t(entry.parameter.has_value()) > order)
            continue;
        Series weight = one;
        for (std::size_t k = 1; k <= ta; ++k) {
            weight = product(weight, dt);
            for (auto &v : weight)
                v /= double(k);
        }
        for (std::size_t k = 1; k <= pb; ++k) {
            weight = product(weight, dp);
            for (auto &v : weight)
                v /= double(k);
        }
        if (entry.parameter) {
            Series theta(order + 1, 0.0);
            const auto index = 2 + entry.column.values.size() + *entry.parameter;
            for (std::size_t k = 1; k <= order; ++k)
                theta[k] = path.phase[phase].point[k][index];
            weight = product(weight, theta);
        }
        if (std::all_of(weight.begin(), weight.end(), [](double v) { return v == 0.0; }))
            continue;
        if (entry.column.status != SolvedStateActionStatus::Available)
            throw CaloricUnavailable(entry.column.status, entry.column.message);
        for (std::size_t k = 0; k <= order; ++k)
            result[k] += entry.column.values[species] * weight[k];
    }
    return result;
}

// P4 supplies the EOS functions; the anchored standard jets supply only their caloric additions.
PropertyChain property_chain(const HelmholtzModel &model, const RowSeries &series, const PhaseMap &map,
                             std::size_t phase, const std::vector<PropertyObservable> &observables,
                             std::size_t order, const ReferenceCalorics *calorics = nullptr) {
    const auto& item = series.phase[phase];
    const std::size_t components = model.parameters().component_count();
    const auto direction = [&](std::size_t k) {
        DensityActionDirection result;
        result.temperature = k == 0 ? 0.0 : series.temperature[k];
        result.molar_density = k == 0 ? 0.0 : map.density[k];
        for (std::size_t i = 0; i < components; ++i) result.composition.push_back(k == 0 ? 0.0 : map.mole_fraction[i][k]);
        // The solved coefficients satisfy the normalisation and charge rows to the linear
        // residual policy; the P4 trust boundary requires them exactly, so the roundoff
        // component along (1, z) is removed within the positive support before the call.
        const auto& charges = model.parameters().charge_number;
        double sum = 0.0, charge = 0.0, zz = 0.0, z1 = 0.0, c = 0.0;
        for (std::size_t i = 0; i < components; ++i) {
            if (!(map.mole_fraction[i][0] > 0.0)) continue;
            c += 1.0;
            sum += result.composition[i]; charge += charges[i] * result.composition[i];
            zz += charges[i] * charges[i]; z1 += charges[i];
        }
        const double det = c * zz - z1 * z1;
        const double alpha = det != 0.0 ? (sum * zz - charge * z1) / det : sum / c;
        const double beta = det != 0.0 ? (c * charge - z1 * sum) / det : 0.0;
        for (std::size_t i = 0; i < components; ++i)
            if (map.mole_fraction[i][0] > 0.0) result.composition[i] -= alpha + beta * charges[i];
        for (std::size_t a = 0; a < model.active_parameters().size(); ++a) {
            const auto &active = model.active_parameters()[a];
            result.parameter.push_back(
                {active.family, active.indices, k == 0 ? 0.0 : item.point[k][2 + components + a]});
        }
        return result;
    };
    std::vector<double> x;
    for (const auto &value : map.mole_fraction)
        x.push_back(value[0]);
    struct Term final { std::vector<std::size_t> orders; double weight; };
    std::vector<Term> terms;
    if (order == 0)
        terms = {{{0}, 1.0}};
    else if (order == 1)
        terms = {{{1}, 1.0}};
    else if (order == 2)
        terms = {{{2}, 1.0}, {{1, 1}, 0.5}};
    else
        terms = {{{3}, 1.0}, {{1, 2}, 1.0}, {{1, 1, 1}, 1.0 / 6.0}};
    std::vector<PropertyObservable> requested;
    for (auto observable : observables) {
        if (calorics && caloric_property(observable))
            observable = property_observable(1, property_kind(observable));
        if (std::find(requested.begin(), requested.end(), observable) == requested.end())
            requested.push_back(observable);
    }
    const bool cv = calorics && std::any_of(observables.begin(), observables.end(), [](auto o) {
                        return caloric_property(o) && property_kind(o) == 4;
                    });
    if (cv) {
        requested.push_back(PropertyObservable::MolarVolumeTemperatureDerivative);
        requested.push_back(PropertyObservable::MolarVolumePressureDerivative);
    }
    std::vector<Series> eos(requested.size(), Series(order + 1, 0.0));
    PropertyChain result;
    result.value.assign(observables.size(), std::nullopt);
    result.coefficient.assign(observables.size(), std::nullopt);
    result.status.assign(observables.size(), SolvedStateActionStatus::Available);
    result.message.resize(observables.size());
    result.series.resize(observables.size());
    std::vector<SolvedStateActionStatus> statuses(requested.size(), SolvedStateActionStatus::Available);
    for (const auto &term : terms) {
        std::vector<DensityActionDirection> directions;
        for (auto k : term.orders)
            directions.push_back(direction(k));
        const auto computed =
            model.density_actions(series.temperature[0], map.density[0], x, requested, directions);
        result.diagnostics = computed.diagnostics;
        result.maximum_energy_order = std::max(result.maximum_energy_order, computed.maximum_energy_order);
        for (std::size_t i = 0; i < requested.size(); ++i) {
            const auto &entry = computed.results[i];
            if (entry.status != PropertyActionStatus::Available)
                statuses[i] = solved_state_status(entry.status);
            if (entry.value)
                eos[i][0] = *entry.value;
            if (order && entry.action)
                eos[i][order] += term.weight * *entry.action;
        }
    }
    // Lower coefficients enter nonlinear products and the frozen-volume quotient.
    for (std::size_t k = 1; k < order; ++k) {
        RowSeries prefix = series;
        prefix.temperature.resize(k + 1);
        prefix.pressure.resize(k + 1);
        for (auto &p : prefix.phase) {
            p.point.resize(k + 1);
            p.gradient.resize(k + 1);
            p.amount.resize(k + 1);
            p.packing.resize(k + 1);
        }
        auto local = map;
        for (auto *v : {&local.amount, &local.volume, &local.density, &local.packing, &local.pressure})
            v->resize(k + 1);
        for (auto *family : {&local.mole_fraction, &local.potential})
            for (auto &v : *family)
                v.resize(k + 1);
        const auto lower = property_chain(model, prefix, local, phase, requested, k);
        for (std::size_t i = 0; i < requested.size(); ++i) {
            if (lower.coefficient[i])
                eos[i][k] = *lower.coefficient[i];
            if (lower.status[i] != SolvedStateActionStatus::Available)
                statuses[i] = lower.status[i];
        }
    }
    const auto eos_series = [&](PropertyObservable o) -> const Series & {
        return eos[std::find(requested.begin(), requested.end(), o) - requested.begin()];
    };
    const double R = GAS_CONSTANT_J_PER_MOL_K;
    for (std::size_t i = 0; i < observables.size(); ++i) {
        const auto o = observables[i];
        const auto kind = property_kind(o);
        const bool referenced = calorics && caloric_property(o);
        const auto residual = referenced ? property_observable(1, kind) : o;
        const auto slot = std::find(requested.begin(), requested.end(), residual) - requested.begin();
        result.status[i] = statuses[slot];
        Series value = eos[slot];
        try {
        if (referenced) {
            Series ideal(order + 1, 0.0), gp(order + 1, 0.0), gtp(order + 1, 0.0), gpp(order + 1, 0.0);
            for (std::size_t species = 0; species < components; ++species) {
                const auto xg = [&](std::size_t a, std::size_t b) {
                    return product(map.mole_fraction[species],
                        standard_series(*calorics, series, phase, species, a, b, order));
                };
                if (kind < 2)
                    ideal = combine(ideal, combine(xg(0, 0), product(series.temperature, xg(1, 0)), -1.0), 1.0);
                else if (kind == 2) {
                    ideal = combine(ideal, xg(1, 0), -1.0);
                    if (x[species] > 0.0)
                        ideal = combine(ideal, product(map.mole_fraction[species],
                            logarithm(product(map.mole_fraction[species], map.density))), -R);
                } else ideal = combine(ideal, product(series.temperature, xg(2, 0)), -1.0);
                if (kind == 1 || kind == 4) gp = combine(gp, xg(0, 1), 1.0);
                if (kind == 4) {
                    gtp = combine(gtp, xg(1, 1), 1.0);
                    gpp = combine(gpp, xg(0, 2), 1.0);
                }
            }
            if (!series.reaction_shift[0].empty()) {
                // Admit the combined source RHS, never each unit column separately: compatible
                // correlated motions may cancel its component outside the anchored inverse's range.
                for (std::size_t j = 0; j <= (kind >= 3 ? 2U : 1U); ++j)
                    for (std::size_t row = 0; row < calorics->reaction_defect.size(); ++row)
                        for (std::size_t k = 0; k <= order; ++k) {
                            double defect = 0.0, scale = 0.0;
                            for (std::size_t r = 0; r < calorics->reaction_map.size(); ++r) {
                                const double rhs = R * series.temperature[0] * series.reaction_shift[j][r][k];
                                defect += calorics->reaction_defect[row][r] * rhs;
                                scale += std::abs(rhs);
                            }
                            if (std::abs(defect) > std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, scale))
                                throw CaloricUnavailable(SolvedStateActionStatus::ReferenceUnavailable,
                                    "source direction inconsistent with fixed record anchors: reaction " + std::to_string(row));
                        }
                for (std::size_t r = 0; r < calorics->reaction_map.size(); ++r) {
                    Series amount(order + 1, 0.0);
                    for (std::size_t species = 0; species < components; ++species)
                        amount = combine(amount, map.mole_fraction[species], calorics->reaction_map[r][species]);
                    const auto& L = series.reaction_shift[0][r];
                    const auto thermal = product(series.temperature, series.reaction_shift[1][r]);
                    const auto shift = kind < 2 ? product(series.temperature, thermal)
                        : kind == 2 ? combine(L, thermal, 1.0)
                        : combine(thermal, product(product(series.temperature, series.temperature), series.reaction_shift[2][r]), 0.5);
                    ideal = combine(ideal, product(amount, shift), kind >= 3 ? 2.0 * R : R);
                }
            }
            if (kind == 0) ideal = combine(ideal, series.temperature, R);
            else if (kind == 2 || kind == 3) ideal[0] += R;
            if (kind == 1) ideal = combine(ideal, product(series.pressure, gp), -1.0);
            if (kind == 4) {
                const auto &vt = eos_series(PropertyObservable::MolarVolumeTemperatureDerivative);
                const auto &vp = eos_series(PropertyObservable::MolarVolumePressureDerivative);
                const auto vG = combine(map.volume, gp, 1.0), vGT = combine(vt, gtp, 1.0),
                           vGP = combine(vp, gpp, 1.0);
                const auto tslot = std::find(requested.begin(), requested.end(),
                                             PropertyObservable::MolarVolumeTemperatureDerivative) -
                                   requested.begin();
                if (statuses[tslot] != SolvedStateActionStatus::Available ||
                    statuses[tslot + 1] != SolvedStateActionStatus::Available)
                    result.status[i] = statuses[tslot] != SolvedStateActionStatus::Available
                                           ? statuses[tslot]
                                           : statuses[tslot + 1];
                else if (!(std::isfinite(vG[0]) && std::isfinite(vGP[0]) && vGP[0] != 0.0))
                    result.status[i] = SolvedStateActionStatus::PressureSingular;
                else
                    ideal = combine(
                        ideal,
                        product(series.temperature, combine(product(product(vGT, vGT), inverse(vGP)),
                                                            product(product(vt, vt), inverse(vp)), -1.0)),
                        1.0);
            }
            value = property_family(o) == 0 ? ideal : combine(value, ideal, 1.0);
        }
        if (result.status[i] == SolvedStateActionStatus::Available) {
            result.series[i] = value;
            result.value[i] = value[0];
            result.coefficient[i] = order ? value[order] : 0.0;
        }
        } catch (const CaloricUnavailable& error) {
            result.status[i] = error.status;
            result.message[i] = error.what();
        }
    }
    return result;
}

// Extensive enthalpy uses the same complete phase jets and Taylor product as all actions.
Series system_enthalpy_series(const HelmholtzModel &model, const RowSeries &series,
                              const std::vector<PhaseMap> &maps, const ReferenceCalorics *calorics,
                              SolvedStateActionResult *result, SolvedStateActionDiagnostics *diagnostics) {
    Series enthalpy(series.temperature.size(), 0.0);
    for (std::size_t phase = 0; phase < maps.size(); ++phase) {
        const auto &amount = maps[phase].amount;
        if (std::all_of(amount.begin(), amount.end(), [](double v) { return v == 0.0; })) continue;
        const auto chain = property_chain(model, series, maps[phase], phase,
                                          {PropertyObservable::TotalEnthalpy}, enthalpy.size() - 1, calorics);
        if (result) result->property_diagnostics = chain.diagnostics;
        if (diagnostics) diagnostics->maximum_energy_order =
            std::max(diagnostics->maximum_energy_order, chain.maximum_energy_order);
        if (chain.status[0] != SolvedStateActionStatus::Available)
            throw CaloricUnavailable(chain.status[0], "phase[" + std::to_string(phase) + "]: " + chain.message[0]);
        enthalpy = combine(enthalpy, product(amount, chain.series[0]), 1.0);
    }
    return enthalpy;
}

// Reference calorics an action reads, at least `floor`: {Gibbs degree order + 1, + 2 for Cp/Cv (xg(2,0)), 2 if nothing
// moves; parameter columns; pressure columns if P moves or V, U, Cv read one}. Apparent Cp is the mixed (T, d) action.
using CaloricNeed = std::array<std::size_t, 3>;
CaloricNeed caloric_need(const EquilibriumProblem &problem, const SolvedStateActionRequest &request, bool zero,
                         CaloricNeed floor) {
    std::size_t extra = 1, apparent = 0, parameters = 0, pressure = 0;
    for (const auto &o : request.observables) {
        apparent |= o.kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity;
        const bool property = o.kind == SolvedStateObservableKind::PhaseProperty && caloric_property(o.property);
        const auto kind = property ? property_kind(o.property) : 0;
        pressure |= o.kind == SolvedStateObservableKind::PhaseGibbsMolarVolume || kind == 1 || kind == 4;
        if (kind >= 3) extra = 2;
    }
    bool moving = !zero && (apparent || problem.temperature.free || problem.pressure.free);
    for (const auto &d : request.directions) {
        for (double v : d.active_parameters) parameters |= !zero && v != 0.0;
        moving |= parameters || (!zero && (d.temperature_k != 0.0 || d.pressure_pa != 0.0));
        pressure |= !zero && (problem.pressure.free || d.pressure_pa != 0.0);
    }
    const std::size_t degree = moving ? std::min<std::size_t>(request.directions.size() + apparent + extra, 4) : 2;
    return {std::max(floor[0], degree), std::max(floor[1], parameters), std::max(floor[2], pressure)};
}
struct ActionCenter final {
    Center center;
    std::vector<SolvedStateActionResult> central_results;
    std::vector<std::vector<std::size_t>> property_rows;
    std::vector<std::vector<PropertyObservable>> property_lists;
    SolvedStateActionDiagnostics diagnostics;
    std::optional<ReferenceCalorics> calorics;
    std::optional<SolvedStateActionBatch> thermal_enthalpy;
    CaloricNeed caloric{};
    SolvedStateActionStatus central_gate = SolvedStateActionStatus::Available;
    SolvedStateActionStatus matrix_gate = SolvedStateActionStatus::Available;
    double classified_residual = 0.0;
    bool topology_checked = false;
    bool matrix_checked = false;
};

}  // namespace

std::string solved_state_observable_id(const SolvedStateObservable& observable) {
    const std::string phase = "phase[" + std::to_string(observable.phase) + "]";
    const std::string component = "[" + std::to_string(observable.component) + "]";
    switch (observable.kind) {
    case SolvedStateObservableKind::SystemEnthalpy: return "system.enthalpy";
    case SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity: return "system.isobaric_heat_capacity";
    case SolvedStateObservableKind::Temperature: return "temperature";
    case SolvedStateObservableKind::Loading: return "loading";
    case SolvedStateObservableKind::Pressure: return "pressure";
    case SolvedStateObservableKind::PhaseAmount: return phase + ".amount";
    case SolvedStateObservableKind::PhaseComponentAmount: return phase + ".amount" + component;
    case SolvedStateObservableKind::PhaseMoleFraction: return phase + ".mole_fraction" + component;
    case SolvedStateObservableKind::PhaseVolume: return phase + ".volume";
    case SolvedStateObservableKind::PhaseGibbsMolarVolume: return phase + ".gibbs_molar_volume";
    case SolvedStateObservableKind::PhaseMolarDensity: return phase + ".molar_density";
    case SolvedStateObservableKind::PhasePackingFraction: return phase + ".packing_fraction";
    case SolvedStateObservableKind::PhasePressure: return phase + ".pressure";
    case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:
        return phase + ".residual_chemical_potential_over_rt" + component;
    case SolvedStateObservableKind::PhaseLogActivity: return phase + ".log_activity" + component;
    case SolvedStateObservableKind::PhaseProperty: return phase + ".property/" + property_name(observable.property);
    case SolvedStateObservableKind::SolidAmount: return "solid[" + std::to_string(observable.phase) + "].amount";
    case SolvedStateObservableKind::SolidSaturationIndex: return "solid[" + std::to_string(observable.phase) + "].saturation_index";
    }
    return "unknown";
}

static SolvedStateActionBatch evaluate_action(
    const CompiledSystem& compiled, const EquilibriumResult& central,
    const SolvedStateActionRequest& request, std::optional<ActionCenter>& shared, core::Scope& scope,
    const CaloricNeed& covering = {}
) {
    const auto& problem = compiled.problem;
    const auto& model = *compiled.model;
    const auto species_target = std::find_if(problem.reactions.begin(), problem.reactions.end(), [](const auto& r) { return r.uses_species_thermochemistry(); });
    if ((model.parameters().ideal_gas || species_target != problem.reactions.end()) && std::any_of(request.observables.begin(), request.observables.end(), [](const auto& o) {
        return o.kind == SolvedStateObservableKind::SystemEnthalpy || o.kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity
            || o.kind == SolvedStateObservableKind::PhaseGibbsMolarVolume
            || (o.kind == SolvedStateObservableKind::PhaseProperty && caloric_property(o.property));
    })) {
        SolvedStateActionBatch refused;
        for (const auto& o : request.observables) {
            SolvedStateActionResult item;
            item.observable = o; item.unit = observable_unit(o); item.order = request.directions.size();
            item.status = model.parameters().ideal_gas ? SolvedStateActionStatus::EnergyOrderUnavailable : SolvedStateActionStatus::ReferenceUnavailable;
            item.diagnostic_message = model.parameters().ideal_gas ? "ideal-gas solved-state calorics are deferred" : "species-derived calorics deferred for reaction " + reaction_name(*species_target, species_target - problem.reactions.begin());
            refused.results.push_back(std::move(item));
        }
        return refused;
    }
    const bool apparent = std::any_of(request.observables.begin(), request.observables.end(), [](const auto &o) {
        return o.kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity;
    });
    if (apparent) {
        auto enthalpy = request;
        for (auto &o : enthalpy.observables)
            if (o.kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity)
                o.kind = SolvedStateObservableKind::SystemEnthalpy;
        SolvedStateActionDirection thermal{1.0, 0.0, std::vector<double>(compiled.n_feed.size()),
                                           std::vector<double>(model.active_parameters().size())};
        const bool new_thermal = !shared || !shared->thermal_enthalpy;
        if (new_thermal) {
            auto first = evaluate_action(compiled, central, {enthalpy.observables, {thermal}}, shared, scope, covering);
            if (!shared) {
                for (std::size_t i = 0; i < first.results.size(); ++i) {
                    first.results[i].observable = request.observables[i];
                    first.results[i].unit = observable_unit(request.observables[i]);
                }
                return first;
            }
            shared->thermal_enthalpy = std::move(first);
        }
        auto ordinary = evaluate_action(compiled, central, enthalpy, shared, scope, covering);
        enthalpy.directions.insert(enthalpy.directions.begin(), thermal);
        const auto mixed = evaluate_action(compiled, central, enthalpy, shared, scope, covering);
        const double feed = std::accumulate(compiled.n_feed.begin(), compiled.n_feed.end(), 0.0);
        const auto &direction = request.directions.front();
        const double dfeed = std::accumulate(direction.feed_amounts_mol.begin(), direction.feed_amounts_mol.end(), 0.0);
        for (std::size_t i = 0; i < ordinary.results.size(); ++i) {
            if (request.observables[i].kind != SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity) continue;
            auto &item = ordinary.results[i];
            const auto &first = shared->thermal_enthalpy->results[i];
            const auto &second = mixed.results[i];
            item.observable = request.observables[i];
            item.unit = observable_unit(item.observable);
            item.value.reset();
            item.action.reset();
            item.status = first.status;
            item.diagnostic_message = first.diagnostic_message;
            if (first.status != SolvedStateActionStatus::Available || !first.action) continue;
            item.value = *first.action / (feed * thermal.temperature_k);
            item.status = second.status;
            item.diagnostic_message = second.diagnostic_message;
            item.property_diagnostics = second.property_diagnostics;
            if (second.status == SolvedStateActionStatus::Available && second.action)
                item.action = *second.action / (feed * thermal.temperature_k) - *item.value * dfeed / feed;
            if (!std::isfinite(*item.value) || (item.action && !std::isfinite(*item.action))) {
                item.status = SolvedStateActionStatus::Nonfinite;
                item.action.reset();
            }
        }
        auto &diagnostics = ordinary.diagnostics;
        for (const auto *part : {&mixed, new_thermal ? static_cast<const SolvedStateActionBatch *>(&*shared->thermal_enthalpy) : nullptr}) {
            if (!part) continue;
            diagnostics.state_factorizations += part->diagnostics.state_factorizations;
            diagnostics.linear_solves += part->diagnostics.linear_solves;
            diagnostics.maximum_energy_order = std::max(diagnostics.maximum_energy_order, part->diagnostics.maximum_energy_order);
            if (part->diagnostics.maximum_linear_residual)
                diagnostics.maximum_linear_residual = std::max(diagnostics.maximum_linear_residual.value_or(0.0),
                                                               *part->diagnostics.maximum_linear_residual);
        }
        return ordinary;
    }
    const std::size_t components = model.parameters().component_count();
    const std::size_t source_offset = compiled.control_width();
    const std::size_t width = compiled.extended_control_width();
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
    // Every result starts CentralSolveUnavailable (the struct default).
    if (!central.success || !central.has_state() || !all_numeric(central.coordinates)) return batch;
    for (auto& item : batch.results) item.status = SolvedStateActionStatus::Available;
    if (!shared) {
        shared.emplace();
        auto& common = *shared;
        common.caloric = covering;
        common.center.prepared = compiled.prepare(central.coordinates);
        // The gate quantity is the classified residual the driver certifies; the unscaled
        // rows stay in the central result's diagnostics.
        const auto evaluation = compiled.evaluate_at(common.center.prepared);
        common.classified_residual = maximum_abs(compiled.classified_residuals(common.center.prepared, evaluation.values));
        common.diagnostics.equation_count = compiled.rows.size();
        common.diagnostics.state_dimension = compiled.unknown_count();
        if (std::isfinite(common.classified_residual)) common.diagnostics.central_residual = common.classified_residual;
        if (evaluation.floor > 0.0) common.diagnostics.central_floor = evaluation.floor;
    }
    auto& center = shared->center;
    batch.diagnostics = shared->diagnostics;
    scope.probe("central_residual", shared->classified_residual);
    std::vector<double> scale(width, 1.0);
    scale[0] = std::max(center.prepared.temperature, 1.0);
    scale[1] = std::max(center.prepared.pressure, 1.0);
    for (std::size_t i = 0; i < components; ++i) scale[2 + i] = std::max(std::abs(compiled.n_feed[i]), 1.0);
    for (std::size_t r = 0; r < problem.reactions.size(); ++r) {
        scale[source_offset + 4 * r + 1] = center.prepared.temperature;
        scale[source_offset + 4 * r + 3] = 1.0 / center.prepared.temperature;
    }
    scale[width - 2] = compiled.feed_total * GAS_CONSTANT_J_PER_MOL_K * center.prepared.temperature;
    scale[width - 1] = compiled.feed_total * GAS_CONSTANT_J_PER_MOL_K;
    std::vector<std::vector<double>> controls;
    std::vector<double> magnitude;
    for (const auto& direction : request.directions) {
        std::vector<double> control = {direction.temperature_k, direction.pressure_pa};
        control.insert(control.end(), direction.feed_amounts_mol.begin(), direction.feed_amounts_mol.end());
        control.insert(control.end(), direction.active_parameters.begin(), direction.active_parameters.end());
        control.resize(width, 0.0);
        for (std::size_t r = 0; r < direction.reaction_coefficients.size(); ++r)
            std::copy(direction.reaction_coefficients[r].begin(), direction.reaction_coefficients[r].end(),
                      control.begin() +
                          static_cast<std::ptrdiff_t>(source_offset + 4 * r));
        control[width - 2] = direction.enthalpy_j;
        control[width - 1] = direction.entropy_j_per_k;
        double largest = 0.0;
        for (std::size_t i = 0; i < width; ++i) largest = std::max(largest, std::abs(control[i]) / scale[i]);
        controls.push_back(std::move(control));
        magnitude.push_back(largest);
    }
    const bool reference_problem = !compiled.reference_diagnostics.empty();
    const ReferenceCalorics *calorics = nullptr;
    const auto other_caloric = [](const auto &o) {
        return o.kind == SolvedStateObservableKind::PhaseGibbsMolarVolume ||
               (o.kind == SolvedStateObservableKind::PhaseProperty && caloric_property(o.property));
    };
    const auto caloric = [&](const auto &o) {
        return o.kind == SolvedStateObservableKind::SystemEnthalpy || other_caloric(o);
    };
    const bool needs_calorics = std::any_of(request.observables.begin(), request.observables.end(), caloric);
    const bool zero = std::any_of(magnitude.begin(), magnitude.end(), [](double m) { return m == 0.0; });
    if (reference_problem && needs_calorics) {
        const auto need = caloric_need(problem, request, zero, shared->caloric);
        if (!shared->calorics || need != shared->caloric) {
            shared->caloric = need;
            shared->calorics = reference_calorics(compiled, center.prepared.temperature, center.prepared.pressure,
                                                  true, need[2], need[1], need[0]);
        }
        calorics = &*shared->calorics;
    }
    const auto refuse = [](SolvedStateActionResult &item, SolvedStateActionStatus status, std::string message,
                           bool value) {
        item.status = status;
        item.action.reset();
        if (value)
            item.value.reset();
        item.diagnostic_message = std::move(message);
    };
    auto& property_rows = shared->property_rows;
    auto& property_lists = shared->property_lists;
    if (shared->central_results.empty()) {
        // Central values and property dependencies are identical for every independent first-order direction.
        const std::vector<double> no_control(width, 0.0);
        const auto base = compiled.row_series(center.prepared, {center.prepared.point}, no_control, 0);
        std::vector<PhaseMap> base_maps;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase)
            base_maps.push_back(phase_map(model, base.phase[phase], base.temperature, components));
        for (auto& item : batch.results) {
            if (item.observable.kind == SolvedStateObservableKind::PhaseProperty) continue;
            try {
                const auto series = observable_series(model, item.observable, base, base_maps, calorics, &item, &batch.diagnostics);
                if (!series.empty() && std::isfinite(series.front())) item.value = series.front();
            } catch (const CaloricUnavailable& error) {
                refuse(item, error.status, error.what(), true);
            }
        }
        property_rows.resize(problem.phases.size());
        property_lists.resize(problem.phases.size());
        for (std::size_t i = 0; i < batch.results.size(); ++i)
            if (batch.results[i].observable.kind == SolvedStateObservableKind::PhaseProperty)
                property_rows[batch.results[i].observable.phase].push_back(i);
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            for (const std::size_t row : property_rows[phase])
                property_lists[phase].push_back(batch.results[row].observable.property);
            if (property_rows[phase].empty())
                continue;
            const auto chain =
                property_chain(model, base, base_maps[phase], phase, property_lists[phase], 0, calorics);
            batch.diagnostics.maximum_energy_order = std::max(batch.diagnostics.maximum_energy_order, chain.maximum_energy_order);
            for (std::size_t i = 0; i < property_rows[phase].size(); ++i) {
                auto& item = batch.results[property_rows[phase][i]];
                item.value = chain.value[i];
                item.property_diagnostics = chain.diagnostics;
                // A one-species critical candidate has P_rho = 0: its isobaric heat capacities are pressure-closed.
                const auto property = item.observable.property;
                if (problem.critical && problem.phases[0].support.size() == 1 && (property == PropertyObservable::ResidualIsobaricHeatCapacity
                                         || property == PropertyObservable::TotalIsobaricHeatCapacity)) {
                    item.status = SolvedStateActionStatus::PressureSingular;
                    item.value.reset();
                } else if (chain.status[i] != SolvedStateActionStatus::Available) {
                    item.status = chain.status[i];
                    item.diagnostic_message = chain.message[i].empty()
                        ? std::string("phase_property dependency: ") + solved_state_status_name(item.status) : chain.message[i];
                }
            }
        }
        if (calorics && calorics->status != SolvedStateActionStatus::Available)
            for (auto &item : batch.results)
                if (caloric(item.observable))
                    refuse(item, calorics->status, calorics->message, true);
        shared->central_results = batch.results;
        shared->diagnostics = batch.diagnostics;
    } else {
        batch.results = shared->central_results;
    }
    for (auto& item : batch.results) {
        item.order = order;
        if (!state_dependent(item.observable, problem)) {
            const std::size_t index = item.observable.kind == SolvedStateObservableKind::Temperature ? 0 : 1;
            item.action = order == 1 ? controls[0][index] : 0.0;
        }
    }
    // Admission of the central state, in the owner's precedence.
    if (!shared->topology_checked) {
        if (!(shared->classified_residual <= std::min(problem.solver_options.value(policy::row), policy::EQUILIBRIUM_ROW_RESIDUAL.value))) {
            shared->central_gate = SolvedStateActionStatus::DefiningEquationsNotSatisfied;
        } else {
            const auto topology = check_topology(compiled, central.coordinates);
            if (topology.event == TopologyEvent::PhaseDisappeared)
                shared->central_gate = SolvedStateActionStatus::PhaseBoundaryNotDifferentiable;
            else if (topology.event == TopologyEvent::PhasesCoalesced)
                shared->central_gate = SolvedStateActionStatus::NonuniqueCoalescedPhases;
        }
        shared->topology_checked = true;
    }
    auto gate = shared->central_gate;
    batch.diagnostics.rank.reset();
    batch.diagnostics.condition_estimate.reset();
    if (gate == SolvedStateActionStatus::Available) {
        const auto& charges = model.parameters().charge_number;
        bool support = false, nonzero_z = false;
        bool changes_inventory = compiled.inventory_restrictions.empty();
        for (const auto& control : controls) {
            const std::vector<double> feed(control.begin() + 2, control.begin() + 2 + compiled.n_feed.size());
            for (const auto& restriction : compiled.inventory_restrictions)
                changes_inventory = changes_inventory || !exact_zero_dot(restriction, feed);
        }
        for (const std::size_t species : compiled.removed_species)
            for (const auto& control : controls)
                if (changes_inventory && control[2 + species] != 0.0) (charges[species] == 0 ? support : nonzero_z) = true;
        if (nonzero_z) gate = SolvedStateActionStatus::IonicInsertion;
        else if (support) gate = SolvedStateActionStatus::SupportChange;
    }
    if (gate == SolvedStateActionStatus::Available) {
        if (!shared->matrix_checked) {
            const std::size_t n = central.coordinates.size();
            const auto jacobian = compiled.jacobian_at(center.prepared);
            center.jacobian = dense_matrix(jacobian, compiled.rows.size(), n);
            if (!square_and_numeric(center.jacobian)) {
                shared->matrix_gate = SolvedStateActionStatus::NonfiniteJacobian;
            } else {
                // Exact power-of-two column then row equilibration; a zero line keeps scale one.
                const auto power = [](double value) { return value > 0.0 ? std::exp2(-std::round(std::log2(value))) : 1.0; };
                center.column = center.jacobian.cwiseAbs().colwise().maxCoeff().transpose().unaryExpr(power);
                center.row = (center.jacobian * center.column.asDiagonal()).cwiseAbs().rowwise().maxCoeff().unaryExpr(power);
                center.jacobian = center.row.asDiagonal() * center.jacobian * center.column.asDiagonal();
                center.lu = center.jacobian.fullPivLu();
                batch.diagnostics.state_factorizations = 1;
                shared->diagnostics.rank = static_cast<std::size_t>(center.lu.rank());
                // cond_inf(R J C), the norm of the backward error it multiplies in the forward bound.
                if (*shared->diagnostics.rank < n) shared->matrix_gate = SolvedStateActionStatus::RankDeficient;
                else shared->diagnostics.condition_estimate = center.jacobian.cwiseAbs().rowwise().sum().maxCoeff()
                    * center.lu.inverse().cwiseAbs().rowwise().sum().maxCoeff();
            }
            shared->matrix_checked = true;
        }
        batch.diagnostics.rank = shared->diagnostics.rank;
        batch.diagnostics.condition_estimate = shared->diagnostics.condition_estimate;
        gate = shared->matrix_gate;
    }
    const bool thermal = problem.temperature.free ||
        std::any_of(controls.begin(), controls.end(), [](const auto &c) { return c[0] != 0.0; });
    if (gate == SolvedStateActionStatus::Available && !zero && thermal
        && thermal_reference_missing(center.prepared.temperature, problem.reactions, problem.standard_gibbs_j_per_mol, model.thermochemistry() ? &*model.thermochemistry() : nullptr))
        gate = SolvedStateActionStatus::MissingThermochemistry;
    // A critical cell's row series reaches Taylor order two (five-direction energy contractions);
    // a third-order action needs order three of the critical rows, which needs six-direction ones.
    if (gate == SolvedStateActionStatus::Available && problem.critical && order > 2)
        gate = SolvedStateActionStatus::EnergyOrderUnavailable;
    // A central-state gate precedes the selected dependency: it overrides every
    // state-dependent row's status, keeping the mapped value and property diagnostics.
    const auto apply_state_gate = [&](SolvedStateActionStatus status,
                                      std::optional<std::string> diagnostic_message = std::nullopt) {
        for (auto& item : batch.results) {
            if (!state_dependent(item.observable, problem)) continue;
            item.status = status;
            item.action.reset();
            item.diagnostic_message = diagnostic_message;
        }
        return batch;
    };
    if (gate != SolvedStateActionStatus::Available)
        return apply_state_gate(gate);
    if (zero) {
        for (auto& item : batch.results)
            if (item.status == SolvedStateActionStatus::Available && item.value.has_value()) item.action = 0.0;
        return batch;
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
        if (path.status != SolvedStateActionStatus::Available)
            return apply_state_gate(path.status, path.diagnostic_message);
        std::vector<PhaseMap> maps;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) maps.push_back(phase_map(model, path.series.phase[phase], path.series.temperature, components));
        for (std::size_t i = 0; i < batch.results.size(); ++i) {
            auto& item = batch.results[i];
            if (item.observable.kind == SolvedStateObservableKind::PhaseProperty || !state_dependent(item.observable, problem)) continue;
            try {
                const auto series = observable_series(model, item.observable, path.series, maps, calorics, &item, &batch.diagnostics);
                accumulated[i] += weight * series[order];
            } catch (const CaloricUnavailable& error) {
                refuse(item, error.status, error.what(), false);
            }
        }
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            if (property_rows[phase].empty()) continue;
            const auto chain = property_chain(model, path.series, maps[phase], phase, property_lists[phase],
                                              order, calorics);
            batch.diagnostics.maximum_energy_order = std::max(batch.diagnostics.maximum_energy_order, chain.maximum_energy_order);
            for (std::size_t i = 0; i < property_rows[phase].size(); ++i) {
                auto& item = batch.results[property_rows[phase][i]];
                item.property_diagnostics = chain.diagnostics;
                if (chain.status[i] != SolvedStateActionStatus::Available && item.status == SolvedStateActionStatus::Available) {
                    item.status = chain.status[i];
                    item.diagnostic_message = chain.message[i].empty()
                        ? std::string("phase_property dependency: ") + solved_state_status_name(item.status) : chain.message[i];
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

SolvedStateActionBatch solved_state_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central, const SolvedStateActionRequest& request
) {
    core::Scope scope("solved_state_actions");
    validate_request(compiled, central, request);
    std::optional<ActionCenter> shared;
    if (compiled.solids.empty()) return evaluate_action(compiled, central, request, shared, scope);
    return evaluate_action(classify_solids(compiled, central.coordinates), central, request, shared, scope);
}

std::vector<SolvedStateActionBatch> solved_state_first_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central,
    const std::vector<SolvedStateObservable>& observables,
    const std::vector<SolvedStateActionDirection>& directions
) {
    core::Scope scope("solved_state_first_actions");
    if (directions.empty()) throw std::invalid_argument("solved-state action order");
    std::vector<SolvedStateActionRequest> requests;
    requests.reserve(directions.size());
    CaloricNeed cover{};  // the center's first reference request covers every direction, a zero one as moving
    for (const auto& direction : directions) {
        requests.push_back({observables, {direction}});
        validate_request(compiled, central, requests.back());
        cover = caloric_need(compiled.problem, requests.back(), false, cover);
    }
    std::optional<ActionCenter> shared;
    std::vector<SolvedStateActionBatch> batches;
    batches.reserve(requests.size());
    const auto classified = compiled.solids.empty() ? std::nullopt : std::optional(classify_solids(compiled, central.coordinates));
    for (const auto& request : requests)
        batches.push_back(evaluate_action(classified ? *classified : compiled, central, request, shared, scope, cover));
    return batches;
}

template <typename Base>
PhysicalProjection<Base> physical_projection(const CompiledSystem& c, const PhasePointValues<Base>& at) {
    const auto n = at.point.size(), components = c.n_feed.size(), per_phase = 3 + components;
    const auto outputs = 2 + at.amount.size() * per_phase + c.problem.loading.has_value();
    PhysicalProjection<Base> out;
    out.gradient = decltype(out.gradient)::Zero(outputs, n);
    out.values = decltype(out.values)::Zero(outputs); out.scale = decltype(out.scale)::Ones(n);
    const auto pressure = c.primal_size - 1, temperature = c.primal_size - c.free_intensive_count;
    out.values[0] = at.temperature; out.values[1] = at.pressure;
    if (c.problem.loading) { out.values[outputs - 1] = at.point[c.loading_column()]; out.gradient(outputs - 1, c.loading_column()) = 1; }
    if (c.problem.pressure.free) { out.gradient(1, pressure) = 1; out.scale[pressure] = std::max(std::abs(at.pressure), Base(1)); }
    if (c.problem.temperature.free) { out.gradient(0, temperature) = 1; out.scale[temperature] = std::max(std::abs(at.temperature), Base(1)); }
    for (const auto& solid : c.solids) if (solid.column) out.scale[*solid.column] = c.feed_total;
    for (std::size_t p = 0; p < at.shape.size(); ++p) {
        const auto& shape = at.shape[p]; const auto o = 2 + p * per_phase;
        if (shape.amount) out.scale[*shape.amount] = c.feed_total;
        for (std::size_t j = 0; j < n; ++j) {
            PhaseSeriesValues<Base> path;
            path.point.assign(2, std::vector<Base>(2 + components, 0)); path.gradient = path.point;
            path.amount = {at.amount[p], Base(shape.amount && j == *shape.amount)};
            path.packing = {at.eta[p], j == shape.volume ? at.eta[p] : Base(0)};
            path.point[0][1] = at.volume[p];
            path.gradient[0] = at.derivative[p].gradient;
            for (std::size_t i = 0; i < components; ++i) path.point[0][2 + i] = at.composition[p][i];
            for (std::size_t i = 0; i < shape.x.size(); ++i) if (j == shape.x[i]) {
                const auto species = c.problem.phases[p].support[i];
                path.point[1][2 + species] = at.composition[p][species];
                path.point[1][1] = at.packing[p].component[species] / at.eta[p];
            }
            if (j == shape.volume) path.point[1][1] = -at.volume[p];
            if (c.problem.temperature.free && j == temperature) path.point[1][1] = at.packing[p].temperature_derivative / at.eta[p];
            const auto map = phase_map(*c.model, path, std::vector<Base>{at.temperature, Base(c.problem.temperature.free && j == temperature)}, components);
            out.gradient(o, j) = map.density[1] / map.density[0];
            out.gradient(o + 1, j) = c.model->parameters().ideal_gas ? map.packing[1] : map.packing[1] / map.packing[0];
            out.gradient(o + 2, j) = map.amount[1];
            out.values[o] = map.density[0]; out.values[o + 1] = map.packing[0]; out.values[o + 2] = map.amount[0];
            for (std::size_t i = 0; i < components; ++i) {
                out.gradient(o + 3 + i, j) = map.mole_fraction[i][1]; out.values[o + 3 + i] = map.mole_fraction[i][0];
            }
        }
    }
    return out;
}
template PhysicalProjection<double> physical_projection(const CompiledSystem&, const PhasePoint&);
template PhysicalProjection<long double> physical_projection(const CompiledSystem&, const PhasePointValues<long double>&);

template <typename Base>
EquilibriumNumericalEstimates numerical_estimates(const CompiledSystem& c, const EquilibriumResult& result, const std::vector<Base>& internal, const PhasePointValues<Base>* prepared = nullptr, bool replay = true, const DriverLimits* limits = nullptr) {
    core::Capture isolated(false);
    EquilibriumNumericalEstimates out;
    out.dimension = result.coordinates.size(); out.precision_bits = std::numeric_limits<Base>::digits;
    using Vector = Eigen::Matrix<Base, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic>;
    if (!result.has_state() || !result.declaration) { out.cause = "original_state_unavailable"; return out; }
    try {
        const auto prepare = [&](const std::vector<Base>& point, bool state_only) {
            if (limits && (!limits->spend(limits->evaluations, limits->evaluation_cap) || limits->expired()))
                throw std::runtime_error("work budget exhausted during returned-point estimate");
            return c.prepare_values(point, state_only);
        };
        const auto at = prepared ? *prepared : prepare(std::vector<Base>(result.coordinates.begin(), result.coordinates.end()), false);
        const auto evaluation = c.evaluate_at(at);
        const auto uncertainty = [&] {
            if constexpr (std::is_same_v<Base, double>) if (replay) return c.row_uncertainty(at);
            std::vector<std::vector<Base>> errors;
            for (const auto& d : at.derivative) errors.push_back(d.gradient_error.empty()
                ? std::vector<Base>(d.gradient.size(), d.gradient_floor) : d.gradient_error);
            return c.row_uncertainty_at(at, errors);
        }();
        out.residuals.assign(evaluation.values.begin(), evaluation.values.end());
        out.row_terms.assign(uncertainty.terms.begin(), uncertainty.terms.end()); out.row_inner.assign(uncertainty.inner.begin(), uncertainty.inner.end());
        if (!c.reference_diagnostics.empty() || !c.solids.empty()) out.source_reference_method =
            "local higher-precision re-evaluation plus first-order pressure/site propagation; "
            "can underestimate shared errors and internal cancellation; fixed-input arithmetic excludes supplied operand uncertainty; "
            "reference replay precision: " + std::to_string(std::numeric_limits<long double>::digits) + " bits";
        if (!uncertainty.cause.empty()) {
            out.temperature.cause = out.pressure.cause = out.cause = uncertainty.cause; return out;
        }
        const auto n = out.dimension, phases = at.amount.size(), components = c.n_feed.size();
        const auto per_phase = 3 + components, outputs = 2 + phases * per_phase + c.problem.loading.has_value();
        Matrix g = Matrix::Zero(outputs, n);
        Vector values = Vector::Zero(outputs), h = Vector::Ones(n);
        values[0] = *result.temperature; values[1] = *result.pressure;
        if (c.problem.loading) {
            values[outputs - 1] = at.point[c.loading_column()];
            g(outputs - 1, c.loading_column()) = 1;
        }
        for (std::size_t p = 0; p < phases; ++p) {
            const auto o = 2 + p * per_phase;
            if (!c.ideal) { values[o] = result.molar_densities[p]; values[o + 1] = result.packing_fractions[p]; }
            values[o + 2] = result.amounts[p];
            for (std::size_t j = 0; j < components; ++j) values[o + 3 + j] = result.mole_fractions[p * components + j];
        }
        // Reuse the actions' exact physical map; no chart is reconstructed in a reader.
        std::vector<Base> mapped(outputs, 0.0);
        if (c.ideal) {
            mapped[0] = at.temperature; mapped[1] = at.pressure; mapped[4] = at.amount[0];
            for (std::size_t i = 0; i < components; ++i) mapped[5 + i] = at.composition[0][i];
            for (std::size_t j = 0; j < n; ++j) {
                const auto species = c.problem.phases[0].support[j];
                const double x = at.composition[0][species];
                g(4, j) = at.amount[0] * x;
                for (std::size_t i = 0; i < components; ++i)
                    g(5 + i, j) = at.composition[0][i] * (double(i == species) - x);
            }
        } else {
            const auto projection = physical_projection(c, at);
            g = projection.gradient; h = projection.scale;
            mapped.assign(projection.values.data(), projection.values.data() + projection.values.size());
        }

        const auto entries = c.jacobian_at(at);
        const Matrix j = Eigen::Map<const Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(entries.data(), n, n);
        const Vector f = Eigen::Map<const Vector>(evaluation.values.data(), n);
        const Vector terms = Eigen::Map<const Vector>(uncertainty.terms.data(), n);
        const Vector inner = Eigen::Map<const Vector>(uncertainty.inner.data(), n);
        const Vector eta = core::canonical_row_floor(j, h, terms, inner);
        out.physical_scale = {h.data(), h.data() + n}; out.row_floor = {eta.data(), eta.data() + n};
        const auto e = core::linear_error(j, f, g, eta, Vector::Ones(n), h);
        out.rank = e.rank; out.threshold = e.threshold; out.condition = e.condition;
        out.inverse_norm = e.inverse_norm; out.backward_error = e.backward;
        if (!e.cause.empty()) { out.cause = e.cause; return out; }
        if (c.problem.loading) mapped.back() = values[outputs - 1];
        const auto internal_map = internal.empty() ? std::vector<Base>{} : [&] {
            const auto projection = physical_projection(c, prepare(internal, true));
            return std::vector<Base>(projection.values.data(), projection.values.data() + projection.values.size());
        }();
        Vector radius = e.radius;
        const auto project = [&](std::size_t i, bool positive) {
            const Base displacement = internal.empty() ? Base(0) : positive
                ? std::abs(std::log1p((mapped[i] - internal_map[i]) / internal_map[i])) : std::abs(mapped[i] - internal_map[i]);
            const auto estimate = core::project_numerical_estimate(double(values[i]), mapped[i], e.radius[i] + displacement, positive, g.row(i).isZero(0.0));
            radius[i] = positive ? estimate.chart : estimate.absolute;
            return estimate;
        };
        out.temperature = project(0, false); out.pressure = project(1, false);
        if (c.problem.loading) out.loading = project(outputs - 1, false);
        bool domain = true;
        for (std::size_t p = 0; p < phases; ++p) {
            const auto o = 2 + p * per_phase;
            if (!c.ideal) { out.densities.push_back(project(o, true)); out.packing.push_back(project(o + 1, !c.model->parameters().ideal_gas)); }
            out.amounts.push_back(project(o + 2, false));
            domain = domain && values[o + 2] - out.amounts.back().absolute >= 0;
            domain = domain && (c.ideal || c.model->parameters().ideal_gas || (values[o + 1] * std::exp(-out.packing.back().chart) >= policy::PACKING_FRACTION_FLOOR.value
                && values[o + 1] * std::exp(out.packing.back().chart) <= policy::PACKING_FRACTION_CEILING.value));
            for (std::size_t i = 0; i < components; ++i) {
                out.composition.push_back(project(o + 3 + i, false));
                domain = domain && values[o + 3 + i] - out.composition.back().absolute >= 0
                    && values[o + 3 + i] + out.composition.back().absolute <= 1;
            }
        }
        for (std::size_t p = 0; p < phases; ++p) for (std::size_t q = p + 1; q < phases; ++q)
            if (c.problem.phases[p].kind != c.problem.phases[q].kind)
                domain = domain && std::abs(result.molar_densities[p] - result.molar_densities[q])
                    > out.densities[p].absolute + out.densities[q].absolute;
        if (limits) limits->check();
        out.physical_radius = physical_radius(c, at, radius);
        out.cause = domain ? "" : "unresolved_domain";
        if (!domain) {
            out.loading.resolved = out.temperature.resolved = out.pressure.resolved = false;
            for (auto* entries : {&out.densities, &out.packing, &out.composition, &out.amounts})
                for (auto& estimate : *entries) { estimate.resolved = false; estimate.cause = out.cause; }
        }
    } catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) { out.cause = error.what(); }
    return out;
}

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem& c, const EquilibriumResult& result, const PhasePoint* prepared) {
    if (result.numerical_estimates && result.numerical_estimates->precision_bits > std::numeric_limits<double>::digits) return *result.numerical_estimates;
    return numerical_estimates<double>(c, result, {}, prepared, !prepared);
}
EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem& c, const EquilibriumResult& result, const std::vector<long double>& internal, const DriverLimits* limits) {
    return numerical_estimates(c, result, internal, static_cast<const PhasePointValues<long double>*>(nullptr), false, limits);
}

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const HelmholtzModel& model, const EquilibriumResult& result) {
    if (result.model_identity != &model || !result.declaration) throw std::invalid_argument("numerical estimate model/result mismatch");
    return equilibrium_numerical_estimates(compile_problem(model, *result.declaration), result);
}

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const IdealMoleFraction& model, const EquilibriumResult& result) {
    if (result.ideal_identity != &model || !result.declaration) throw std::invalid_argument("numerical estimate model/result mismatch");
    return equilibrium_numerical_estimates(compile_problem(model, *result.declaration), result);
}

}  // namespace epcsaft::equilibrium
