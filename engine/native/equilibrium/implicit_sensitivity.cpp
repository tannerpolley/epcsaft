#include "equilibrium/gibbs_minimizer.hpp"

#include "contributions/electrolyte.hpp"
#include "core/bitwise_equal.hpp"
#include "core/trace.hpp"
#include "core/linear_error.hpp"
#include "core/numerical_policy.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

// Instantiated once in equilibrium/singular_value_decompositions.cpp.
extern template class Eigen::JacobiSVD<Eigen::MatrixXd>;

namespace epcsaft::equilibrium {
namespace {

using RowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Eigen::MatrixXd dense_matrix(const std::vector<double>& values, std::size_t rows, std::size_t columns) {
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

// Exact power-of-two column then row equilibration of a; a zero line keeps scale one. Returns {row, column}.
template <typename Matrix>
std::pair<Eigen::Matrix<typename Matrix::Scalar, Eigen::Dynamic, 1>, Eigen::Matrix<typename Matrix::Scalar, Eigen::Dynamic, 1>>
power_of_two_equilibration(const Matrix& a) {
    using Real = typename Matrix::Scalar;
    const auto power = [](Real value) { using std::exp2, std::log2, std::round; return value > 0 ? Real(exp2(-round(log2(value)))) : Real(1); };
    const Eigen::Matrix<Real, Eigen::Dynamic, 1> column = a.cwiseAbs().colwise().maxCoeff().transpose().unaryExpr(power);
    return {(a * column.asDiagonal()).cwiseAbs().rowwise().maxCoeff().unaryExpr(power), column};
}

struct EquilibratedRows {
    Eigen::MatrixXd jacobian;
    Eigen::VectorXd row, column;
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
    void equilibrate() {
        std::tie(row, column) = power_of_two_equilibration(jacobian);
        jacobian = row.asDiagonal() * jacobian * column.asDiagonal();
    }

    SolvedStateActionStatus factor() {
        lu = jacobian.fullPivLu();
        return lu.rank() == jacobian.rows() ? SolvedStateActionStatus::Available : SolvedStateActionStatus::RankDeficient;
    }

    double norm() const { return jacobian.cwiseAbs().rowwise().sum().maxCoeff(); }
    double condition() const { return norm() * lu.inverse().cwiseAbs().rowwise().sum().maxCoeff(); }
};

SolvedStateActionStatus linear_solve_status(double error, double condition) {
    if (!std::isfinite(error)) return SolvedStateActionStatus::NonfiniteJacobian;
    if (error > policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR.value) return SolvedStateActionStatus::LinearSolveResidual;
    return 2.0 * condition * error <= policy::EQUILIBRIUM_ACTION_FORWARD_ERROR.value
        ? SolvedStateActionStatus::Available : SolvedStateActionStatus::IllConditioned;
}

}  // namespace

SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    bool estimate_condition
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
    const auto [row_scale, column_scale] = power_of_two_equilibration(matrix);
    const Eigen::MatrixXd equilibrated = row_scale.asDiagonal() * matrix * column_scale.asDiagonal();
    const auto lu = equilibrated.fullPivLu();
    SensitivityInfo result;
    result.rank = static_cast<std::size_t>(lu.rank());
    result.threshold = lu.threshold();
    if (estimate_condition || scope.on()) result.condition_estimate = condition_estimate(matrix);
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
    const std::vector<double>& column, std::optional<core::Convergence>* failure, std::vector<double> jacobian
) {
    core::Scope scope("implicit_sensitivity");
    const std::size_t rows = compiled.rows.size(), dimension = prepared.point.size();
    if (!valid_point(prepared.point)) throw std::invalid_argument("invalid sensitivity point");
    if (!all_numeric(column)) throw std::invalid_argument("invalid parameter derivative");
    if (jacobian.empty()) jacobian = compiled.jacobian_at(prepared);
    const auto refuse = [&](const char* cause, core::Stop stop, const auto* factor, std::optional<double> backward = std::nullopt) {
        if (failure) { core::Convergence r; r.cause = cause;
            if (factor) r.rank = std::make_shared<core::NumericalRank>(core::numerical_rank(*factor, "implicit defining rows"));
            r.residual = backward.value_or(NAN); r.backward_error = backward; r.stop = stop;
            if (backward) r.criterion = {std::shared_ptr<const policy::Tolerance>{}, &policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR};
            core::retain_first(*failure, r);
        }
        return std::vector<double>{};
    };
    if (jacobian.size() != rows * dimension || rows != dimension) return refuse("implicit Jacobian shape", core::Stop::infeasible, static_cast<const Eigen::FullPivLU<Eigen::MatrixXd>*>(nullptr));
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    const auto [row_scale, column_scale] = power_of_two_equilibration(matrix);
    const Eigen::MatrixXd equilibrated = row_scale.asDiagonal() * matrix * column_scale.asDiagonal();
    const auto lu = equilibrated.fullPivLu();
    if (!matrix.allFinite() || static_cast<std::size_t>(lu.rank()) != dimension)
        return refuse(matrix.allFinite() ? "implicit Jacobian rank" : "implicit nonfinite Jacobian", matrix.allFinite() ? core::Stop::singular : core::Stop::nonfinite, &lu);
    // One column, kept a matrix so that Eigen takes the same solve path.
    const Eigen::MatrixXd right = -Eigen::Map<const Eigen::VectorXd>(column.data(), static_cast<Eigen::Index>(column.size()));
    const Eigen::MatrixXd result = column_scale.asDiagonal() * lu.solve(row_scale.asDiagonal() * right);
    const double error = std::max(0.0, backward_error((matrix * result.col(0) - right.col(0)).cwiseAbs().maxCoeff(),
        matrix.cwiseAbs().rowwise().sum().maxCoeff(), result.cwiseAbs().maxCoeff(), right.cwiseAbs().maxCoeff()));
    if (!result.allFinite() || !(error <= policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR.value))
        return refuse(result.allFinite() ? "implicit backward error" : "implicit nonfinite response", result.allFinite() ? core::Stop::stalled : core::Stop::nonfinite, &lu,
            result.allFinite() ? std::optional<double>(error) : std::nullopt);
    scope.observe([&] { scope.probe("max_solution", result.cwiseAbs().maxCoeff()); });
    scope.probe("backward_error", error);
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    return {result.data(), result.data() + rows};
}

bool check_held_responses(const CompiledSystem& compiled, const std::function<Eigen::MatrixXd()>& jacobian,
    const Eigen::VectorXd& scale, EquilibriumResult& result) {
    using Status = SolvedStateActionStatus;
    if (const auto& parameters = compiled.model->parameters(); uses_solvent_ratio(parameters)
        && std::any_of(parameters.charge_number.begin(), parameters.charge_number.end(), [](int z) { return z != 0; })) {
        for (std::size_t p = 0; p < compiled.shape.size(); ++p)
            if (std::count_if(result.held_species.begin(), result.held_species.end(), [&](const auto& cell) {
                return cell.phase == p && dielectric_solvent_component(cell.species, parameters);
            }) >= 2)
                return held_refusal(result, core::Stop::infeasible, "hold_dependent", "phase " + compiled.problem.phases[p].name + " uses a ratio of held solvents");
    }
    EquilibratedRows factor;
    factor.jacobian = jacobian();
    if (factor.jacobian.rows() != factor.jacobian.cols())
        return held_refusal(result, core::Stop::infeasible, "held_response_unavailable", "classified Jacobian is not square");
    if (!factor.jacobian.allFinite())
        return held_refusal(result, core::Stop::nonfinite, "held_response_unavailable", "nonfinite_jacobian");
    factor.equilibrate();
    if (!factor.jacobian.allFinite())
        return held_refusal(result, core::Stop::nonfinite, "held_response_unavailable", "scaled Jacobian is nonfinite");
    if (factor.factor() == Status::RankDeficient)
        return held_refusal(result, core::Stop::singular, "held_response_unavailable", "rank_deficient");
    const double condition = factor.condition();
    auto record = result.held_species.begin();
    for (std::size_t r = 0; r < compiled.rows.size(); ++r) {
        const auto& row = compiled.rows[r];
        if (!row.held_column) continue;
        const Eigen::VectorXd right = factor.row[r] * Eigen::VectorXd::Unit(scale.size(), r);
        const Eigen::VectorXd scaled = factor.lu.solve(right), response = factor.column.cwiseProduct(scaled).cwiseQuotient(scale);
        const double error = backward_error((factor.jacobian * scaled - right).lpNorm<Eigen::Infinity>(),
            factor.norm(), scaled.lpNorm<Eigen::Infinity>(), right.lpNorm<Eigen::Infinity>());
        const auto status = response.allFinite() ? linear_solve_status(error, condition) : Status::NonfiniteJacobian;
        const auto cell = "species " + std::to_string(row.index) + " in phase " + compiled.problem.phases[row.phase].name;
        if (status != Status::Available)
            return held_refusal(result, status == Status::NonfiniteJacobian ? core::Stop::nonfinite : core::Stop::singular,
                "held_response_unavailable", cell + ": " + solved_state_status_name(status));
        double norm = 0;
        for (std::size_t k = 0; k < compiled.primal_size; ++k)
            if (k != *row.held_column) norm = std::max(norm, std::abs(response[k]));
        (record++)->response = norm;
        const auto stop = core::classify(norm, policy::EQUILIBRIUM_POLISH_CORRECTION, {}, {}, false);
        if (stop != core::Stop::converged)
            return held_refusal(result, stop, "hold_dependent", "state depends on the held fraction of " + cell + ": " + core::optional_number(norm));
    }
    return true;
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
    const auto prepared = compiled.prepare(point, false, nullptr, {}, {}, true);
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
    if (value.empty() || !wide_math::isfinite(value.front()) || value.front() == 0.0)
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
    if (value.empty() || !wide_math::isfinite(value.front()) || !(value.front() > 0.0))
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
        && kind != SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity && kind != SolvedStateObservableKind::SpeciesStandardMinusPlane && !is_solid_kind(kind);
}

bool has_component(SolvedStateObservableKind kind) {
    return kind == SolvedStateObservableKind::SpeciesStandardMinusPlane || kind == SolvedStateObservableKind::PhaseComponentAmount
        || kind == SolvedStateObservableKind::PhaseMoleFraction
        || kind == SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt
        || kind == SolvedStateObservableKind::PhaseLogActivity;
}

bool same_observable(const SolvedStateObservable& a, const SolvedStateObservable& b) {
    if (a.kind != b.kind) return false;
    if ((is_phase_kind(a.kind) || is_solid_kind(a.kind)) && a.phase != b.phase) return false;
    if (has_component(a.kind) && a.component != b.component) return false;
    if (a.kind == SolvedStateObservableKind::AssociationSiteFraction) return a.site == b.site;
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

bool moves_eos_parameters(const SolvedStateActionDirection& d) {
    return std::any_of(d.active_parameters.begin(), d.active_parameters.end(), [](double v) { return v != 0.0; });
}

void validate_request(const CompiledSystem& compiled, const EquilibriumResult& central,
                      const SolvedStateActionRequest& request) {
    if (compiled.problem.prescribed_liquid_composition) throw std::invalid_argument("solved-state actions unavailable for prescribed liquid composition");
    if (compiled.ideal) throw std::invalid_argument("ideal solved-state actions unavailable");
    const auto& problem = compiled.problem;
    if (central.has_state() && !same_result_origin(compiled, central))
        throw std::invalid_argument("central result does not match the compiled system");
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
            && (!compiled.admits_apparent_cp() || request.directions.size() != 1))
            throw std::invalid_argument("apparent Cp requires driven T/P or free T with a fixed-P enthalpy target, and first actions only");
        if (item.kind == SolvedStateObservableKind::PhaseBulkRelativePermittivity && request.directions.size() > 1)
            throw std::invalid_argument("bulk relative permittivity supports first actions only");
        if (item.kind == SolvedStateObservableKind::AssociationSiteFraction) {
            validate_phase_site(compiled.model->parameters(), problem.phases[item.phase], item.site);
            if (request.directions.size() != 1)
                throw std::invalid_argument("association site fraction supports first actions only");
        }
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
        if ((!direction.solid_log_k.empty() && direction.solid_log_k.size() != problem.solids.size()) || !all_numeric(direction.solid_log_k))
            throw std::invalid_argument("solid ln K direction shape or nonfinite");
        if (!std::isfinite(direction.internal_energy_j) || (!problem.internal_energy_j && direction.internal_energy_j != 0.0))
            throw std::invalid_argument("internal energy control without its active condition or nonfinite");
        if (!std::isfinite(direction.volume_m3) || (!problem.volume_m3 && direction.volume_m3 != 0.0))
            throw std::invalid_argument("invalid volume target direction");
        if (!std::isfinite(direction.enthalpy_j) || !std::isfinite(direction.entropy_j_per_k)
            || (!problem.enthalpy_j && direction.enthalpy_j != 0.0)
            || (!problem.entropy_j_per_k && direction.entropy_j_per_k != 0.0))
            throw std::invalid_argument("caloric control without its active condition or nonfinite");
        if (problem.temperature.free && direction.temperature_k != 0.0)
            throw std::invalid_argument("free temperature is a solved output, not an external control");
        if (problem.pressure.free && direction.pressure_pa != 0.0)
            throw std::invalid_argument("free pressure is a solved output, not an external control");
    }
    // A referenced problem: compilation resolves a diagnostic for every source-referenced reaction.
    const bool referenced = !compiled.reference_diagnostics.empty();
    std::size_t eos_parameter_directions = 0;
    bool source_directions = false;
    for (const auto &d : request.directions) {
        bool parameter = moves_eos_parameters(d);
        eos_parameter_directions += std::size_t(parameter);
        // ln K^rho_s is compiled to second order in (T, theta); V_s (P - P0)/RT is exact.
        if (!problem.solids.empty() && request.directions.size() == 3 && (parameter || d.temperature_k != 0.0 || d.pressure_pa != 0.0))
            throw std::invalid_argument("solids admit T, P and parameter directions through order two");
        for (const auto &field : d.reaction_coefficients)
            for (double v : field)
                parameter |= v != 0.0, source_directions |= v != 0.0;
        for (double v : d.solid_log_k) parameter |= v != 0.0;
        if (referenced && request.directions.size() == 3 &&
            (parameter || d.temperature_k != 0.0 || d.pressure_pa != 0.0 || problem.temperature.free ||
             problem.pressure.free))
            throw std::invalid_argument("a referenced third-order action admits only feed directions at fixed T and P");
    }
    // Every observable of an anchored PH/PS cell follows its caloric row, whose g carries one EOS-parameter
    // coordinate; the row has no reaction-law source shift.
    if (compiled.anchored_calorics() && (eos_parameter_directions > 1 || source_directions))
        throw std::invalid_argument("anchored PH/PS admits at most one EOS-parameter direction and no reaction-law direction");
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

struct Center final : EquilibratedRows {
    PhasePoint prepared;
};

struct FirstCoefficient final {
    std::vector<double> state;
    double right_norm, solution_norm;
    std::size_t energy_order;
};

struct Path final {
    std::optional<double> refused_backward_error;
    RowSeries series;
    SolvedStateActionStatus status = SolvedStateActionStatus::Available;
    std::optional<std::string> diagnostic_message;
};

SolvedStateActionStatus energy_status(PropertyActionStatus status) {
    return status == PropertyActionStatus::Nonfinite ? SolvedStateActionStatus::NonfiniteJacobian
        : solved_state_status(status);
}

// The coefficient recurrence: J y_k = -[t^k] F(y(t), q(t)) with y_k unknown set to zero,
// re-evaluated after each solve so the retained residual is the owned row residual.
Path follow(const CompiledSystem& compiled, const Center& center, const std::vector<double>& control,
            std::size_t order, SolvedStateActionDiagnostics& diagnostics, bool plane, bool parameter_pairs,
            const std::vector<std::size_t>& parameter_union, const FirstCoefficient* first = nullptr,
            FirstCoefficient* retain = nullptr) {
    Path path;
    const std::size_t n = center.prepared.point.size();
    std::vector<std::vector<double>> state(order + 1, std::vector<double>(n, 0.0));
    state[0] = center.prepared.point;
    std::vector<double> right_norm(order + 1, 0.0), solution_norm(order + 1, 0.0);
    // The series through order k, or false with the path refused when a reference is unavailable.
    const auto rows_through = [&](const std::vector<std::vector<double>>& prefix, std::size_t k, RowSeries& series) {
        try {
            series = compiled.row_series(center.prepared, prefix, control, k, plane, parameter_pairs, parameter_union);
            return true;
        } catch (const ReferenceUnavailable& error) {
            path.status = SolvedStateActionStatus::ReferenceUnavailable;
            path.diagnostic_message = error.what();
            return false;
        }
    };
    if (first) {
        state[1] = first->state;
        right_norm[1] = first->right_norm;
        solution_norm[1] = first->solution_norm;
        diagnostics.maximum_energy_order = std::max(diagnostics.maximum_energy_order, first->energy_order);
    }
    for (std::size_t k = first ? 2 : 1; k <= order; ++k) {
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
        if (retain && k == 1) *retain = {state[k], right_norm[k], solution_norm[k], series.maximum_energy_order};
    }
    if (!rows_through(state, order, path.series)) return path;
    if (compiled.has_caloric_condition() && std::any_of(path.series.temperature.begin() + 1,
            path.series.temperature.end(), [](double v) { return v != 0.0; })) {
        const auto calorics = reference_calorics(compiled, center.prepared.temperature, center.prepared.pressure, {});
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
        residual = backward_error(residual, center.norm(), solution_norm[k], right_norm[k]);
        const auto status = linear_solve_status(residual, *diagnostics.condition_estimate);
        if (status == SolvedStateActionStatus::NonfiniteJacobian) { path.status = status; return path; }
        diagnostics.maximum_linear_residual = std::max(diagnostics.maximum_linear_residual.value_or(0.0), residual);
        // The first refused coefficient names the cause; 2 cond eta is the first-order forward bound.
        if (path.status != SolvedStateActionStatus::Available || status == SolvedStateActionStatus::Available) continue;
        path.status = status;
        path.refused_backward_error = residual;
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
    case SolvedStateObservableKind::SpeciesStandardMinusPlane: return series.standard_minus_plane.at(observable.component);
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
    case SolvedStateObservableKind::PhaseBulkRelativePermittivity: {
        const auto& map = maps[observable.phase];
        const auto& point = series.phase[observable.phase].point;
        return model.bulk_relative_permittivity_series(series.temperature, map.density, map.mole_fraction,
            point.size() == 1 ? Series{} : Series(point[1].end() - model.active_parameters().size(), point[1].end()));
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
        if (ta + pb + std::size_t(entry.parameter.has_value()) + std::size_t(entry.second_parameter.has_value()) > order)
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
        for (const auto parameter : {entry.parameter, entry.second_parameter}) {
            if (!parameter)
                continue;
            Series theta(order + 1, 0.0);
            const auto index = 2 + entry.column.values.size() + *parameter;
            for (std::size_t k = 1; k <= order; ++k)
                theta[k] = path.phase[phase].point[k][index];
            weight = product(weight, theta);
        }
        if (entry.second_parameter && entry.second_parameter == entry.parameter)
            for (auto &v : weight)
                v *= 0.5;
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
    std::vector<std::size_t> pairs;
    SolvedStateActionStatus central_gate = SolvedStateActionStatus::Available;
    SolvedStateActionStatus matrix_gate = SolvedStateActionStatus::Available;
    double classified_residual = 0.0;
    bool preparation_failed = false;
    bool topology_checked = false;
    bool matrix_checked = false;
};

using FirstCoefficients = std::map<std::string, FirstCoefficient>;

}  // namespace

void validate_phase_site(const Parameters& parameters, const PhaseSpec& phase, std::optional<std::size_t> site) {
    if (!site) throw std::invalid_argument("association site selector is required");
    if (*site >= parameters.site_component.size()) throw std::invalid_argument("association site ordinal outside the model");
    if (std::find(phase.support.begin(), phase.support.end(), parameters.site_component[*site]) == phase.support.end())
        throw std::invalid_argument("association site owner outside the phase support");
}

HelmholtzDerivatives phase_site_values(const HelmholtzModel& model, const EquilibriumResult& central, std::size_t phase) {
    const auto count = model.parameters().component_count();
    const auto start = central.mole_fractions.begin() + phase * count;
    return model.evaluate(*central.temperature, 1.0 / central.molar_densities[phase],
                          std::vector<double>(start, start + count), DerivativeOrder::Value);
}

std::string solved_state_observable_id(const SolvedStateObservable& observable) {
    const std::string phase = "phase[" + std::to_string(observable.phase) + "]";
    const std::string component = "[" + std::to_string(observable.component) + "]";
    switch (observable.kind) {
    case SolvedStateObservableKind::SystemEnthalpy: return "system.enthalpy";
    case SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity: return "system.isobaric_heat_capacity";
    case SolvedStateObservableKind::SpeciesStandardMinusPlane: return "standard_minus_plane" + component;
    case SolvedStateObservableKind::Temperature: return "temperature";
    case SolvedStateObservableKind::Loading: return "loading";
    case SolvedStateObservableKind::Pressure: return "pressure";
    case SolvedStateObservableKind::PhaseAmount: return phase + ".amount";
    case SolvedStateObservableKind::PhaseComponentAmount: return phase + ".amount" + component;
    case SolvedStateObservableKind::PhaseMoleFraction: return phase + ".mole_fraction" + component;
    case SolvedStateObservableKind::PhaseVolume: return phase + ".volume";
    case SolvedStateObservableKind::PhaseGibbsMolarVolume: return phase + ".gibbs_molar_volume";
    case SolvedStateObservableKind::PhaseMolarDensity: return phase + ".molar_density";
    case SolvedStateObservableKind::PhaseBulkRelativePermittivity: return phase + ".bulk_relative_permittivity";
    case SolvedStateObservableKind::AssociationSiteFraction:
        return phase + ".site_fraction[" + (observable.site ? std::to_string(*observable.site) : "unset") + "]";
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
    const CaloricNeed& covering, std::optional<PhasePoint>* prepared, const std::vector<std::size_t>& parameter_union = {},
    FirstCoefficients* first_coefficients = nullptr, bool batch_parameters = false
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
        // Driven T/P: dH along T. Free T with an enthalpy target: DT along H, so c = delta/(N_f DT[w_H]).
        const bool enthalpy_target = problem.temperature.free;
        const double feed = std::accumulate(compiled.n_feed.begin(), compiled.n_feed.end(), 0.0);
        const auto kind = enthalpy_target ? SolvedStateObservableKind::Temperature : SolvedStateObservableKind::SystemEnthalpy;
        auto enthalpy = request;
        for (auto &o : enthalpy.observables)
            if (o.kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity) o.kind = kind;
        // Once the thermal column is retained, a lone Cp row consumes only the mixed derivative.
        const bool mixed_only = request.observables.size() == 1 && shared && !shared->preparation_failed && shared->thermal_enthalpy.has_value();
        auto ordinary = mixed_only ? SolvedStateActionBatch{} :
            evaluate_action(compiled, central, enthalpy, shared, scope, covering, prepared, parameter_union, nullptr, batch_parameters);
        SolvedStateActionDirection thermal{enthalpy_target ? 0.0 : 1.0, 0.0, std::vector<double>(compiled.n_feed.size()),
                                           std::vector<double>(model.active_parameters().size())};
        if (enthalpy_target) thermal.enthalpy_j = feed * GAS_CONSTANT_J_PER_MOL_K * central.temperature.value_or(1.0);
        const double delta = enthalpy_target ? thermal.enthalpy_j : thermal.temperature_k;
        const bool center_available = shared && !shared->preparation_failed;
        const bool new_thermal = center_available && !shared->thermal_enthalpy;
        const auto &direction = request.directions.front();
        if (new_thermal)
            shared->thermal_enthalpy = evaluate_action(compiled, central, {{{kind}}, {thermal}}, shared, scope, covering, prepared, parameter_union);
        const auto mixed = center_available ? evaluate_action(compiled, central, {{{kind}}, {thermal, direction}}, shared, scope, covering, prepared)
                                  : SolvedStateActionBatch{};
        if (mixed_only) {
            ordinary = mixed;
            ordinary.results.front().order = request.directions.size();
        }
        const double dfeed = std::accumulate(direction.feed_amounts_mol.begin(), direction.feed_amounts_mol.end(), 0.0);
        for (std::size_t i = 0; i < ordinary.results.size(); ++i) {
            if (request.observables[i].kind != SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity) continue;
            auto &item = ordinary.results[i];
            item.observable = request.observables[i];
            item.unit = observable_unit(item.observable);
            if (!center_available) continue;
            const auto &first = shared->thermal_enthalpy->results.front();
            const auto &second = mixed.results.front();
            item.value.reset();
            item.action.reset();
            item.status = first.status;
            item.diagnostic_message = first.diagnostic_message;
            item.refusal = first.refusal;
            if (first.status != SolvedStateActionStatus::Available || !first.action) continue;
            item.value = enthalpy_target ? delta / (feed * *first.action) : *first.action / (feed * delta);
            item.status = second.status;
            item.diagnostic_message = second.diagnostic_message;
            item.property_diagnostics = second.property_diagnostics;
            item.refusal = second.refusal;
            if (second.status == SolvedStateActionStatus::Available && second.action)
                item.action = (enthalpy_target ? -delta * *second.action / (feed * *first.action * *first.action)
                                               : *second.action / (feed * delta)) - *item.value * dfeed / feed;
            if (!std::isfinite(*item.value) || (item.action && !std::isfinite(*item.action))) {
                item.status = SolvedStateActionStatus::Nonfinite;
                item.action.reset();
            }
        }
        if (!shared) return ordinary;
        auto &diagnostics = ordinary.diagnostics;
        for (const auto *part : {&mixed, new_thermal ? static_cast<const SolvedStateActionBatch *>(&*shared->thermal_enthalpy) : nullptr}) {
            if (!part || (mixed_only && part == &mixed)) continue;
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
    const bool plane_uncertainty = std::any_of(request.observables.begin(), request.observables.end(),
        [](const auto& o) { return o.kind == SolvedStateObservableKind::SpeciesStandardMinusPlane && o.plane_uncertainty; });
    const auto parameter_directions = std::count_if(request.directions.begin(), request.directions.end(), moves_eos_parameters);
    for (auto& item : batch.results) item.status = SolvedStateActionStatus::Available;
    if (shared && shared->preparation_failed) {
        batch.results = shared->central_results;
        return batch;
    }
    if (!shared) {
        shared.emplace();
        auto& common = *shared;
        common.caloric = covering;
        // A caller's point is the fresh preparation of this same system and state; each center works on its own copy.
        std::vector<bool> site_phases;
        for (const auto& item : request.observables) if (item.kind == SolvedStateObservableKind::AssociationSiteFraction) {
            site_phases.resize(problem.phases.size());
            site_phases[item.phase] = true;
        }
        const auto refuse_preparation = [&](SolvedStateActionStatus status, const std::exception& error) {
            std::map<std::size_t, HelmholtzDerivatives> values;
            for (auto& item : batch.results) {
                item.status = status;
                item.diagnostic_message = error.what();
                auto& diag = item.property_diagnostics.emplace();
                if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error))
                    diag.first_numerical_failure = failure->record;
                else core::retain_failure(diag.first_numerical_failure, error, "solved-state action center");
                if (item.observable.kind != SolvedStateObservableKind::AssociationSiteFraction) continue;
                const auto phase = item.observable.phase;
                auto [it, inserted] = values.try_emplace(phase);
                if (inserted) try {
                    it->second = phase_site_values(model, central, phase);
                } catch (const core::ConvergenceFailure&) {
                    // The center's first failure already identifies the refused action.
                }
                if (!it->second.site_fractions.empty()) item.value = it->second.site_fractions.at(*item.observable.site);
            }
            common.central_results = batch.results;
            common.preparation_failed = true;
            return batch;
        };
        try {
            if (order == 0 && central.prepared_endpoint && central.prepared_endpoint->point == central.coordinates)
                common.center.prepared = *central.prepared_endpoint;
            else if (order == 0)
                common.center.prepared = compiled.prepare_values(central.coordinates, true, nullptr, {}, {}, false, true);
            else if (!site_phases.empty()) common.center.prepared = compiled.prepare(central.coordinates, false, nullptr, {}, site_phases, true);
            else {
                if (prepared && !*prepared) *prepared = compiled.prepare(central.coordinates, false, nullptr, {}, {}, true);
                common.center.prepared = prepared ? **prepared : compiled.prepare(central.coordinates, false, nullptr, {}, {},
                    !plane_uncertainty || batch_parameters || parameter_directions != 0);
            }
        } catch (const AssociationJacobianFailure& error) {
            if (site_phases.empty()) throw;
            return refuse_preparation(solved_state_status(error.status), error);
        } catch (const core::ConvergenceFailure& error) {
            if (site_phases.empty() || !error.record.criterion
                || std::string_view(error.record.criterion->name) != policy::ASSOCIATION_RESIDUAL.name) throw;
            return refuse_preparation(SolvedStateActionStatus::AssociationResidual, error);
        }
        // The gate quantity is the classified residual the driver certifies; the unscaled
        // rows stay in the central result's diagnostics.
        if (order != 0) {
            const auto evaluation = compiled.evaluate_at(common.center.prepared);
            common.classified_residual = maximum_abs(compiled.classified_residuals(common.center.prepared, evaluation.values));
            if (evaluation.floor > 0.0) common.diagnostics.central_floor = evaluation.floor;
        } else common.classified_residual = maximum_abs(central.classified_residuals);
        common.diagnostics.equation_count = compiled.rows.size();
        common.diagnostics.state_dimension = compiled.unknown_count();
        if (std::isfinite(common.classified_residual)) common.diagnostics.central_residual = common.classified_residual;
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
    const auto target = width - 2 - problem.volume_m3.has_value() - problem.internal_energy_j.has_value();
    if (problem.volume_m3) scale[target + 2] = *problem.volume_m3;
    if (problem.internal_energy_j) scale.back() = compiled.feed_total * GAS_CONSTANT_J_PER_MOL_K * center.prepared.temperature;
    scale[target] = compiled.feed_total * GAS_CONSTANT_J_PER_MOL_K * center.prepared.temperature;
    scale[target + 1] = compiled.feed_total * GAS_CONSTANT_J_PER_MOL_K;
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
        std::copy(direction.solid_log_k.begin(), direction.solid_log_k.end(),
                  control.begin() + static_cast<std::ptrdiff_t>(source_offset + 4 * problem.reactions.size()));
        control[target] = direction.enthalpy_j;
        control[target + 1] = direction.entropy_j_per_k;
        if (problem.volume_m3) control[target + 2] = direction.volume_m3;
        if (problem.internal_energy_j) control.back() = direction.internal_energy_j;
        double largest = 0.0;
        for (std::size_t i = 0; i < width; ++i) largest = std::max(largest, std::abs(control[i]) / scale[i]);
        controls.push_back(std::move(control));
        magnitude.push_back(largest);
    }
    const bool plane = std::any_of(request.observables.begin(), request.observables.end(), [](const auto& o) {
        return o.kind == SolvedStateObservableKind::SpeciesStandardMinusPlane; });
    const bool reference_problem = !compiled.reference_diagnostics.empty() || compiled.anchored_calorics();
    const ReferenceCalorics *calorics = nullptr;
    const auto other_caloric = [](const auto &o) {
        return o.kind == SolvedStateObservableKind::PhaseGibbsMolarVolume ||
               (o.kind == SolvedStateObservableKind::PhaseProperty && caloric_property(o.property));
    };
    const auto caloric = [&](const auto &o) {
        return o.kind == SolvedStateObservableKind::SystemEnthalpy || other_caloric(o);
    };
    const bool needs_calorics = std::any_of(request.observables.begin(), request.observables.end(), caloric);
    const bool zero = order == 0 || std::any_of(magnitude.begin(), magnitude.end(), [](double m) { return m == 0.0; });
    if (reference_problem && needs_calorics) {
        const auto need = order == 0 ? CaloricNeed{1, 0, 0} : caloric_need(problem, request, zero, shared->caloric);
        std::vector<std::size_t> pairs;  // the coordinates the directions move, as row_series
        for (std::size_t k = 0; parameter_directions > 1 && !zero && k < model.active_parameters().size(); ++k)
            if (std::any_of(request.directions.begin(), request.directions.end(), [k](const auto& d) { return d.active_parameters[k] != 0.0; }))
                pairs.push_back(k);
        if (!shared->calorics || need != shared->caloric || pairs != shared->pairs) {
            shared->caloric = need;
            shared->pairs = pairs;
            ReferenceJets jets;
            jets.temperature_order = need[0];
            jets.pressure_order = need[2] ? jets.temperature_order : 0;
            if (need[1]) jets.request_every_parameter(model);
            jets.pairs = pairs;
            shared->calorics = reference_calorics(compiled, center.prepared.temperature, center.prepared.pressure, jets);
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
    // A present solid has no caloric record: the system enthalpy, and apparent Cp through it, is refused, never fluid-only.
    const auto observe = [&](SolvedStateActionResult& item, const RowSeries& series, const std::vector<PhaseMap>& maps) {
        for (std::size_t s = 0; item.observable.kind == SolvedStateObservableKind::SystemEnthalpy && s < compiled.solids.size(); ++s)
            if (compiled.solids[s].present)
                throw CaloricUnavailable(SolvedStateActionStatus::MissingThermochemistry, "solid " + problem.solids[s].name + " has no enthalpy record");
        if (item.observable.kind == SolvedStateObservableKind::AssociationSiteFraction) {
            const auto phase = item.observable.phase;
            const auto site = *item.observable.site;
            const auto& derivative = center.prepared.derivative[phase];
            Series projected{derivative.site_fractions[site]};
            const auto& point = series.phase[phase].point;
            if (point.size() > 1)
                projected.push_back(derivative.site_sensitivity.row(site).dot(
                    Eigen::Map<const Eigen::VectorXd>(point[1].data(), point[1].size())));
            return projected;
        }
        return observable_series(model, item.observable, series, maps, calorics, &item, &batch.diagnostics);
    };
    auto& property_rows = shared->property_rows;
    auto& property_lists = shared->property_lists;
    // Apparent Cp's auxiliary requests select only enthalpy or temperature from the ordinary cache.
    const bool thermal_subset = batch.results.size() < shared->central_results.size();
    if (shared->central_results.empty()) {
        // Central values and property dependencies are identical for every independent first-order direction.
        const std::vector<double> no_control(width, 0.0);
        const auto base = compiled.row_series(center.prepared, {center.prepared.point}, no_control, 0, plane);
        std::vector<PhaseMap> base_maps;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase)
            base_maps.push_back(phase_map(model, base.phase[phase], base.temperature, components));
        for (auto& item : batch.results) {
            if (item.observable.kind == SolvedStateObservableKind::PhaseProperty) continue;
            try {
                const auto series = observe(item, base, base_maps);
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
    } else if (thermal_subset) {
        batch.results.front() = *std::find_if(shared->central_results.begin(), shared->central_results.end(),
            [&](const auto& item) { return same_observable(item.observable, request.observables.front()); });
    } else {
        batch.results = shared->central_results;
    }
    if (order == 0) return batch;
    for (auto& item : batch.results) {
        item.order = order;
        if (!state_dependent(item.observable, problem)) {
            const std::size_t index = item.observable.kind == SolvedStateObservableKind::Temperature ? 0 : 1;
            item.action = order == 1 ? controls[0][index] : 0.0;
        }
    }
    // Admission of the central state, in the owner's precedence.
    const auto* row = problem.solver_options.criterion(policy::row);
    if (!shared->topology_checked) {
        if (!(shared->classified_residual <= std::min(row->value, policy::EQUILIBRIUM_ROW_RESIDUAL.value))) {
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
            for (const auto& restriction : compiled.inventory_restrictions) {
                if (exact_zero_dot(restriction, feed)) continue;
                changes_inventory = true;
                // A discarded charge balance constrains feed directions even with every species retained.
                nonzero_z = nonzero_z || std::equal(restriction.begin(), restriction.end(), charges.begin(), charges.end());
            }
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
            shared->matrix_gate = SolvedStateActionStatus::NonfiniteJacobian;
            if (square_and_numeric(center.jacobian)) {
                center.equilibrate();
                shared->matrix_gate = center.factor();
                batch.diagnostics.state_factorizations = 1;
                shared->diagnostics.rank = static_cast<std::size_t>(center.lu.rank());
                if (shared->matrix_gate == SolvedStateActionStatus::Available) shared->diagnostics.condition_estimate = center.condition();
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
    if (gate == SolvedStateActionStatus::Available && (problem.critical || problem.spinodal) && order > 2)
        gate = SolvedStateActionStatus::EnergyOrderUnavailable;
    // A central-state gate precedes the selected dependency: it overrides every
    // state-dependent row's status, keeping the mapped value and property diagnostics.
    const auto apply_state_gate = [&](SolvedStateActionStatus status,
                                      std::optional<std::string> diagnostic_message = std::nullopt,
                                      std::optional<double> backward = {}) {
        auto refusal = core::refusal(solved_state_status_name(status));
        if (status == SolvedStateActionStatus::DefiningEquationsNotSatisfied) {
            refusal.residual = shared->classified_residual;
            refusal.criterion = {problem.solver_options.owners[policy::row],
                row->value < policy::EQUILIBRIUM_ROW_RESIDUAL.value ? row : &policy::EQUILIBRIUM_ROW_RESIDUAL};
        }
        if (status == SolvedStateActionStatus::NonfiniteJacobian) refusal.stop = core::Stop::nonfinite;
        if (status == SolvedStateActionStatus::RankDeficient) refusal.stop = core::Stop::singular;
        if (status == SolvedStateActionStatus::RankDeficient || backward) {
            refusal.rank = std::make_shared<core::NumericalRank>(core::NumericalRank{
                "RJC", central.coordinates.size(), *batch.diagnostics.rank, center.lu.threshold()});
            refusal.rank->observation.condition = batch.diagnostics.condition_estimate.value_or(INFINITY);
        }
        if (backward) {
            const bool forward = status == SolvedStateActionStatus::IllConditioned;
            if (forward) refusal.stop = core::Stop::ill_conditioned;
            refusal.backward_error = backward;
            refusal.residual = forward ? 2.0 * *batch.diagnostics.condition_estimate * *backward : *backward;
            refusal.criterion = {std::shared_ptr<const policy::Tolerance>{}, forward
                ? &policy::EQUILIBRIUM_ACTION_FORWARD_ERROR : &policy::EQUILIBRIUM_ACTION_BACKWARD_ERROR};
        }
        for (auto& item : batch.results) {
            if (!state_dependent(item.observable, problem)) continue;
            item.status = status;
            item.action.reset();
            item.diagnostic_message = diagnostic_message;
            item.refusal = refusal;
        }
        return batch;
    };
    if (gate != SolvedStateActionStatus::Available)
        return apply_state_gate(gate);
    if (plane_uncertainty) {
        const auto n = central.coordinates.size();
        Eigen::MatrixXd projection = Eigen::MatrixXd::Zero(components, n);
        std::vector<double> unit(n), control(width);
        for (std::size_t j = 0; j < n; ++j) {
            unit[j] = 1.0;
            const auto series = compiled.row_series(center.prepared, {central.coordinates, unit}, control, 1, true);
            for (std::size_t i = 0; i < components; ++i) projection(i, j) = series.standard_minus_plane[i][1] * center.column[j];
            unit[j] = 0.0;
        }
        const auto estimates = equilibrium_numerical_estimates(compiled, central, &center.prepared, true, true);
        if (estimates.row_floor.size() != n)
            return apply_state_gate(SolvedStateActionStatus::IllConditioned,
                estimates.cause.empty() ? "plane row estimates unavailable" : estimates.cause);
        const Eigen::VectorXd f = center.row.cwiseProduct(Eigen::Map<const Eigen::VectorXd>(central.residuals.data(), n));
        const Eigen::VectorXd floor = center.row.cwiseAbs().cwiseProduct(Eigen::Map<const Eigen::VectorXd>(estimates.row_floor.data(), n));
        const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);
        const auto error = core::linear_error(center.jacobian, f, projection, floor, ones, ones, &center.lu);
        if (!error.cause.empty()) return apply_state_gate(SolvedStateActionStatus::IllConditioned, error.cause);
        for (auto& item : batch.results) if (item.observable.kind == SolvedStateObservableKind::SpeciesStandardMinusPlane && item.observable.plane_uncertainty)
            item.numerical_error = error.radius[item.observable.component];
    }
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
        std::string key;
        const FirstCoefficient* first = nullptr;
        if (first_coefficients) {
            key.assign(reinterpret_cast<const char*>(control.data()), control.size() * sizeof(double));
            const auto found = first_coefficients->find(key);
            if (found != first_coefficients->end()) first = &found->second;
        }
        FirstCoefficient retained{};
        const bool retain = first_coefficients && !first && (order == 1 || key.find_first_not_of('\0') == std::string::npos);
        const auto path = follow(compiled, center, control, order, batch.diagnostics, plane, parameter_directions > 1, parameter_union,
            first, retain ? &retained : nullptr);
        if (path.status != SolvedStateActionStatus::Available)
            return apply_state_gate(path.status, path.diagnostic_message,
                path.status == SolvedStateActionStatus::NonfiniteJacobian ? std::nullopt : path.refused_backward_error);
        if (retain) first_coefficients->emplace(std::move(key), std::move(retained));
        std::vector<PhaseMap> maps;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) maps.push_back(phase_map(model, path.series.phase[phase], path.series.temperature, components));
        for (std::size_t i = 0; i < batch.results.size(); ++i) {
            auto& item = batch.results[i];
            if (item.observable.kind == SolvedStateObservableKind::PhaseProperty || !state_dependent(item.observable, problem)) continue;
            try {
                const auto series = observe(item, path.series, maps);
                accumulated[i] += weight * series[order];
            } catch (const CaloricUnavailable& error) {
                refuse(item, error.status, error.what(), false);
            }
        }
        for (std::size_t phase = 0; !thermal_subset && phase < problem.phases.size(); ++phase) {
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

static SolvedStateActionBatch actions_at(const CompiledSystem& compiled, const EquilibriumResult& central,
    const SolvedStateActionRequest& request, std::optional<PhasePoint>* prepared, core::Scope& scope) {
    validate_request(compiled, central, request);
    std::optional<ActionCenter> shared;
    if (compiled.solids.empty() && central.held_species.empty()) return evaluate_action(compiled, central, request, shared, scope, {}, prepared);
    return evaluate_action(classify_result(compiled, central), central, request, shared, scope, {}, prepared);
}

SolvedStateActionBatch solved_state_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central, const SolvedStateActionRequest& request
) {
    core::Scope scope("solved_state_actions");
    return actions_at(compiled, central, request, nullptr, scope);
}

struct SolvedStateActionWorkspace final {
    std::shared_ptr<const char> model;
    std::shared_ptr<const EquilibriumProblem> declaration;
    std::vector<double> point;
    std::optional<PhasePoint> prepared;
    std::optional<ActionCenter> amount_center;
    FirstCoefficients first_coefficients;
};

SolvedStateActionResult solved_state_enthalpy_value(const CompiledSystem& compiled, const EquilibriumResult& central) {
    core::Scope scope("solved_state_enthalpy_value");
    std::optional<ActionCenter> shared;
    const auto value = [&](const auto& system) {
        return evaluate_action(system, central, {{{SolvedStateObservableKind::SystemEnthalpy}}, {}}, shared, scope, {}, nullptr).results.front();
    };
    return compiled.solids.empty() && central.held_species.empty() ? value(compiled) : value(classify_result(compiled, central));
}

SolvedStateActionBatch solved_state_actions_prepared(const CompiledSystem& compiled, const EquilibriumResult& central,
    const SolvedStateActionRequest& request, std::shared_ptr<SolvedStateActionWorkspace>& workspace) {
    core::Scope scope("solved_state_actions");
    validate_request(compiled, central, request);
    if (!workspace || !central.success || !central.has_state() || workspace->model != central.model_identity
        || workspace->declaration != central.declaration || !core::bitwise_equal(workspace->point, central.coordinates)) {
        workspace = std::make_shared<SolvedStateActionWorkspace>();
        workspace->model = central.model_identity;
        workspace->declaration = central.declaration;
        workspace->point = central.coordinates;
    }
    const auto& problem = compiled.problem;
    const bool fixed_feed_actions = compiled.solids.empty() && central.held_species.empty() && compiled.model->active_parameters().empty()
        && !problem.temperature.free && !problem.pressure.free && !compiled.has_caloric_condition()
        && request.directions.size() <= 2
        && std::all_of(request.directions.begin(), request.directions.end(), [](const auto& d) {
            return d.temperature_k == 0 && d.pressure_pa == 0 && d.reaction_coefficients.empty();
        })
        && std::all_of(request.observables.begin(), request.observables.end(), [](const auto& o) {
            return o.kind == SolvedStateObservableKind::PhaseComponentAmount;
        });
    if (!fixed_feed_actions) return actions_at(compiled, central, request, &workspace->prepared, scope);
    auto batch = evaluate_action(compiled, central, request, workspace->amount_center, scope, {}, &workspace->prepared, {},
        &workspace->first_coefficients);
    if (workspace->amount_center) workspace->amount_center->center.prepared.energy_actions.clear();
    return batch;
}

struct SolvedStateCaloricCenter final {
    SolvedStateObservableKind kind;
    std::shared_ptr<const HelmholtzModel> model;
    std::shared_ptr<const EquilibriumProblem> declaration;
    ActionCenter value;
};

std::vector<SolvedStateActionBatch> solved_state_first_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central,
    const std::vector<SolvedStateObservable>& observables,
    const std::vector<SolvedStateActionDirection>& directions, std::shared_ptr<const SolvedStateCaloricCenter>* caloric
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
    const bool retain = caloric && central.declaration && central.model_identity == compiled.model->identity() && observables.size() == 1
        && (observables.front().kind == SolvedStateObservableKind::SystemEnthalpy
            || observables.front().kind == SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity);
    const bool reuse = retain && central.success && central.has_state() && *caloric && (*caloric)->model->identity() == compiled.model->identity()
        && (*caloric)->declaration == central.declaration && (*caloric)->kind == observables.front().kind
        && core::bitwise_equal((*caloric)->value.center.prepared.point, central.coordinates);
    if (reuse) shared = (*caloric)->value;
    else if (caloric) caloric->reset();

    std::vector<SolvedStateActionBatch> batches;
    batches.reserve(requests.size());
    const auto classified = compiled.solids.empty() && central.held_species.empty() ? std::nullopt : std::optional(classify_result(compiled, central));
    // Solids share one first-order source reference. Reacting directions share only when their T/P requests agree.
    const auto group = [&](const SolvedStateActionDirection& d) {
        return compiled.reference_diagnostics.empty() ? 0 : (d.temperature_k != 0.0) + 2 * (d.pressure_pa != 0.0);
    };
    std::array<std::vector<std::size_t>, 4> parameters;
    if (!compiled.reference_diagnostics.empty() || !compiled.solids.empty())
        for (std::size_t a = 0; a < compiled.model->active_parameters().size(); ++a)
            for (const auto& d : directions)
                if (auto& positions = parameters[group(d)]; d.active_parameters[a] != 0.0 && (positions.empty() || positions.back() != a))
                    positions.push_back(a);
    const bool batch_parameters = std::any_of(directions.begin(), directions.end(), moves_eos_parameters);
    for (const auto& request : requests)
        batches.push_back(evaluate_action(classified ? *classified : compiled, central, request, shared, scope,
            cover, nullptr, parameters[group(request.directions.front())], nullptr, batch_parameters));
    if (retain && !reuse && shared && !shared->preparation_failed)
        if (auto owner = compiled.model->weak_from_this().lock())
            *caloric = std::make_shared<const SolvedStateCaloricCenter>(SolvedStateCaloricCenter{
                observables.front().kind, std::move(owner), central.declaration, std::move(*shared)});
    return batches;
}

template <typename Base>
static auto physical_scale(const CompiledSystem& c, const PhasePointValues<Base>& at) {
    using Vector = typename PhysicalProjection<Base>::Vector;
    Vector scale = Vector::Ones(at.point.size());
    if (c.problem.pressure.free) scale[c.primal_size - 1] = pressure_scale(at.pressure);
    if (c.problem.temperature.free) scale[c.primal_size - c.free_intensive_count] = std::max(wide_math::abs(at.temperature), Base(1));
    for (const auto& solid : c.solids) if (solid.column) scale[*solid.column] = c.feed_total;
    for (const auto& shape : at.shape) if (shape.amount) scale[*shape.amount] = c.feed_total;
    return scale;
}

template <typename Base>
PhysicalProjection<Base> physical_projection(const CompiledSystem& c, const PhasePointValues<Base>& at) {
    const auto n = at.point.size(), components = c.component_count(), per_phase = 3 + components;
    const auto outputs = 2 + at.amount.size() * per_phase + c.problem.loading.has_value();
    PhysicalProjection<Base> out;
    out.gradient = decltype(out.gradient)::Zero(outputs, n);
    out.values = decltype(out.values)::Zero(outputs); out.scale = physical_scale(c, at);
    const auto pressure = c.primal_size - 1, temperature = c.primal_size - c.free_intensive_count;
    out.values[0] = at.temperature; out.values[1] = at.pressure;
    if (c.problem.loading) { out.values[outputs - 1] = at.point[c.loading_column()]; out.gradient(outputs - 1, c.loading_column()) = 1; }
    if (c.problem.pressure.free) out.gradient(1, pressure) = 1;
    if (c.problem.temperature.free) out.gradient(0, temperature) = 1;
    for (std::size_t p = 0; p < at.shape.size(); ++p) {
        const auto& shape = at.shape[p]; const auto o = 2 + p * per_phase;
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
template PhysicalProjection<WiderScalar> physical_projection(const CompiledSystem&, const PhasePointValues<WiderScalar>&);

template <typename Base>
EquilibriumNumericalEstimates numerical_estimates(const CompiledSystem& declared, const EquilibriumResult& result, const std::vector<Base>& internal, const PhasePointValues<Base>* prepared = nullptr, bool replay = true, const DriverLimits* limits = nullptr, bool counted = false, bool plane_rows = false,
    const std::vector<double>& prescribed_errors = {}) {
    core::Capture isolated(false);
    core::CalculationWork excluded(counted);
    EquilibriumNumericalEstimates out;
    out.dimension = result.coordinates.size(); out.precision_bits = std::numeric_limits<Base>::digits;
    using Vector = Eigen::Matrix<Base, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic>;
    if (!result.has_state() || !result.declaration) { out.cause = "original_state_unavailable"; return out; }
    try {
        const auto classified = declared.solids.empty() && result.held_species.empty() ? std::nullopt : std::optional(classify_result(declared, result));
        const auto& c = classified ? *classified : declared;
        const auto prepare = [&](const std::vector<Base>& point, bool state_only) {
            if (limits && (!limits->spend(limits->evaluations, limits->evaluation_cap) || limits->expired()))
                throw std::runtime_error("work budget exhausted during returned-point estimate");
            return c.prepare_values(point, state_only);
        };
        const auto at = prepared ? *prepared : prepare(std::vector<Base>(result.coordinates.begin(), result.coordinates.end()), false);
        const auto evaluation = plane_rows ? PointEvaluationValues<Base>{} : c.evaluate_at(at);
        const auto uncertainty = [&] {
            if constexpr (std::is_same_v<Base, double>) if (replay) return c.row_uncertainty(at);
            std::vector<std::vector<Base>> errors;
            for (const auto& d : at.derivative) errors.push_back(d.gradient_error.empty()
                ? std::vector<Base>(d.gradient.size(), d.gradient_floor) : d.gradient_error);
            return c.row_uncertainty_at(at, errors);
        }();
        out.residuals = narrow_values(evaluation.values.begin(), evaluation.values.end());
        out.row_terms = narrow_values(uncertainty.terms.begin(), uncertainty.terms.end()); out.row_inner = narrow_values(uncertainty.inner.begin(), uncertainty.inner.end());
        if (!c.reference_diagnostics.empty() || !c.solids.empty()) out.source_reference_method =
            "local higher-precision re-evaluation plus first-order pressure/site propagation; "
            "can underestimate shared errors and internal cancellation; fixed-input arithmetic excludes supplied operand uncertainty; "
            "reference replay precision: " + std::to_string(std::numeric_limits<WiderScalar>::digits) + " bits";
        if (!uncertainty.cause.empty()) {
            out.temperature.cause = out.pressure.cause = out.cause = uncertainty.cause; return out;
        }
        const auto n = out.dimension, phases = at.amount.size(), components = c.component_count();
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
            values[o + 2] = c.problem.prescribed_liquid_composition ? at.amount[p] : Base(result.amounts[p]);
            for (std::size_t j = 0; j < components; ++j) values[o + 3 + j] = result.mole_fractions[p * components + j];
        }
        // Reuse the actions' exact physical map; no chart is reconstructed in a reader.
        std::vector<Base> mapped(outputs, 0.0);
        if (plane_rows) h = physical_scale(c, at);
        else if (c.ideal) {
            mapped[0] = at.temperature; mapped[1] = at.pressure; mapped[4] = at.amount[0];
            for (std::size_t i = 0; i < components; ++i) mapped[5 + i] = at.composition[0][i];
            for (std::size_t j = 0; j < n; ++j) {
                const auto species = c.problem.phases[0].support[j];
                const double x = static_cast<double>(at.composition[0][species]);
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

        const Vector terms = Eigen::Map<const Vector>(uncertainty.terms.data(), n);
        const Vector inner = Eigen::Map<const Vector>(uncertainty.inner.data(), n);
        Vector eta = core::canonical_row_floor(j, h, terms, inner);
        if (!prescribed_errors.empty()) {
            if (!c.problem.prescribed_liquid_composition || prescribed_errors.size() != components)
                throw std::invalid_argument("prescribed composition uncertainty shape");
            for (std::size_t i = 0; i < components; ++i)
                if ((*c.problem.prescribed_liquid_composition)[i] == 0.0 && prescribed_errors[i] != 0.0)
                    throw std::invalid_argument("prescribed composition uncertainty changes support");
            for (std::size_t i = 0; i < c.rows.size(); ++i) {
                if (c.rows[i].kind != RowKind::PrescribedComposition) continue;
                const auto species = c.problem.phases[0].support[c.rows[i].index];
                const double x = (*c.problem.prescribed_liquid_composition)[species];
                const double u = prescribed_errors[species];
                if (!(u >= 0 && u < x)) throw std::invalid_argument("prescribed composition uncertainty reaches support boundary");
                eta[i] -= std::log1p(-u / x);
            }
        }
        out.physical_scale = narrow_values(h.data(), h.data() + n); out.row_floor = narrow_values(eta.data(), eta.data() + n);
        if (plane_rows) return out;
        const Vector f = Eigen::Map<const Vector>(evaluation.values.data(), n);
        // Equilibrate J diag(h) before the rank decision; eta stays from the physical h.
        const auto [row, column] = power_of_two_equilibration(Matrix(j * h.asDiagonal()));
        const auto e = core::linear_error(j, f, g, eta, row, h.cwiseProduct(column));
        out.rank = e.rank; out.threshold = static_cast<double>(e.threshold); out.condition = static_cast<double>(e.condition);
        out.inverse_norm = static_cast<double>(e.inverse_norm); out.backward_error = static_cast<double>(e.backward);
        if (!e.cause.empty()) { out.cause = e.cause; return out; }
        if (c.problem.loading) mapped.back() = values[outputs - 1];
        const auto internal_map = internal.empty() ? std::vector<Base>{} : [&] {
            const auto projection = physical_projection(c, prepare(internal, true));
            return std::vector<Base>(projection.values.data(), projection.values.data() + projection.values.size());
        }();
        Vector radius = e.radius;
        const auto project = [&](std::size_t i, bool positive) {
            const Base displacement = internal.empty() ? Base(0) : positive
                ? wide_math::abs(wide_math::log1p((mapped[i] - internal_map[i]) / internal_map[i])) : wide_math::abs(mapped[i] - internal_map[i]);
            const auto estimate = equilibrium::project_numerical_estimate(double(values[i]), mapped[i], e.radius[i] + displacement, positive, g.row(i).isZero(0.0));
            radius[i] = positive ? estimate.chart : estimate.absolute;
            return estimate;
        };
        out.temperature = project(0, false); out.pressure = project(1, false);
        if (c.problem.loading) out.loading = project(outputs - 1, false);
        bool domain = true;
        for (std::size_t p = 0; p < phases; ++p) {
            const auto o = 2 + p * per_phase;
            if (!c.ideal) { out.densities.push_back(project(o, true)); out.packing.push_back(project(o + 1, !c.model->parameters().ideal_gas)); }
            if (!c.problem.prescribed_liquid_composition) {
                out.amounts.push_back(project(o + 2, false));
                domain = domain && values[o + 2] - out.amounts.back().absolute >= 0;
            }
            domain = domain && (c.ideal || c.model->parameters().ideal_gas || (values[o + 1] * wide_math::exp(-out.packing.back().chart) >= policy::PACKING_FRACTION_FLOOR.value
                && values[o + 1] * wide_math::exp(out.packing.back().chart) <= policy::PACKING_FRACTION_CEILING.value));
            for (std::size_t i = 0; i < components; ++i) {
                out.composition.push_back(project(o + 3 + i, false));
                domain = domain && values[o + 3 + i] - out.composition.back().absolute >= 0
                    && values[o + 3 + i] + out.composition.back().absolute <= 1;
            }
        }
        for (std::size_t p = 0; p < phases; ++p) for (std::size_t q = p + 1; q < phases; ++q)
            if (c.problem.phases[p].kind != c.problem.phases[q].kind)
                domain = domain && wide_math::abs(result.molar_densities[p] - result.molar_densities[q])
                    > out.densities[p].absolute + out.densities[q].absolute;
        if (limits) limits->check();
        out.physical_radius = static_cast<double>(physical_radius(c, at, radius));
        out.cause = domain ? "" : "unresolved_domain";
        if (!domain) {
            out.loading.resolved = out.temperature.resolved = out.pressure.resolved = false;
            for (auto* entries : {&out.densities, &out.packing, &out.composition, &out.amounts})
                for (auto& estimate : *entries) { estimate.resolved = false; estimate.cause = out.cause; }
        }
    } catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) { out.cause = error.what(); }
    return out;
}

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem& c, const EquilibriumResult& result, const PhasePoint* prepared, bool counted, bool plane_rows,
    const std::vector<double>& prescribed_errors, bool replay) {
    if (result.has_state() && !same_result_origin(c, result))
        throw std::invalid_argument("numerical estimate result does not match the compiled system");
    const auto* endpoint = result.prepared_endpoint.get();
    if (!endpoint || !core::bitwise_equal(endpoint->point, result.coordinates)) endpoint = nullptr;
    if (prepared != endpoint) prepared = nullptr;
    if (!replay && result.success && result.has_state() && result.declaration && prescribed_errors.empty() && !plane_rows
        && endpoint && result.numerical_estimates && result.numerical_estimates->precision_bits > std::numeric_limits<double>::digits)
        return *result.numerical_estimates;
    return numerical_estimates<double>(c, result, {}, prepared, replay || plane_rows || !prepared, nullptr, counted, plane_rows, prescribed_errors);
}
EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem& c, const EquilibriumResult& result, const std::vector<WiderScalar>& internal, const DriverLimits* limits) {
    if (result.has_state() && !same_result_origin(c, result))
        throw std::invalid_argument("numerical estimate result does not match the compiled system");
    return numerical_estimates(c, result, internal, static_cast<const PhasePointValues<WiderScalar>*>(nullptr), false, limits, core::CalculationWork::active_limit());
}

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const HelmholtzModel& model, const EquilibriumResult& result) {
    if (result.model_identity != model.identity() || !result.declaration) throw std::invalid_argument("numerical estimate model/result mismatch");
    return equilibrium_numerical_estimates(compile_problem(model, *result.declaration), result);
}

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const IdealMoleFraction& model, const EquilibriumResult& result) {
    if (result.ideal_identity != model.identity() || !result.declaration) throw std::invalid_argument("numerical estimate model/result mismatch");
    return equilibrium_numerical_estimates(compile_problem(model, *result.declaration), result);
}

// ALG-ACTIONS verification: each neighbor re-solves from `central`; a failed, coalesced or branch-changed neighbor
// is a status of every record, never a value.
std::vector<core::DifferenceCheck> solved_state_difference(const CompiledSystem& compiled, const EquilibriumResult& central,
    const std::vector<SolvedStateObservable>& observables, const SolvedStateActionDirection& d, double step, double tau, double kappa) {
    core::Scope scope("solved_state_difference");
    core::admit_difference_parameters(step, tau, kappa);
    const auto& p = compiled.problem;
    const auto zero = [](double value) { return value == 0.0; };
    const bool feed_zero = std::all_of(d.feed_amounts_mol.begin(), d.feed_amounts_mol.end(), zero);
    const bool supported = compiled.model && std::all_of(d.active_parameters.begin(), d.active_parameters.end(), zero)
        && std::all_of(d.reaction_coefficients.begin(), d.reaction_coefficients.end(),
            [&](const auto& c) { return std::all_of(c.begin(), c.end(), zero); })
        && std::all_of(d.solid_log_k.begin(), d.solid_log_k.end(), zero)
        && (zero(d.temperature_k) || !p.temperature.free) && (zero(d.pressure_pa) || !p.pressure.free)
        && (zero(d.internal_energy_j) || p.internal_energy_j) && (zero(d.volume_m3) || p.volume_m3) && (zero(d.enthalpy_j) || p.enthalpy_j) && (zero(d.entropy_j_per_k) || p.entropy_j_per_k)
        && (feed_zero || d.feed_amounts_mol.size() == compiled.n_feed.size());
    std::vector<core::DifferenceCheck> out;
    for (const auto& o : observables) out.push_back({solved_state_observable_id(o), "unsupported_direction"});
    if (!supported) return out;
    const SolvedStateActionRequest request{observables, {d}};
    const auto exact = solved_state_actions(compiled, central, request);
    std::array<std::vector<double>, 2> values;
    for (std::size_t side = 0; side < 2; ++side) {
        const double h = side == 0 ? step : -step;
        auto q = p;
        q.temperature.value += h * d.temperature_k;
        q.pressure.value += h * d.pressure_pa;
        if (q.internal_energy_j) *q.internal_energy_j += h * d.internal_energy_j;
        if (q.volume_m3) *q.volume_m3 += h * d.volume_m3;
        if (q.enthalpy_j) *q.enthalpy_j += h * d.enthalpy_j;
        if (q.entropy_j_per_k) *q.entropy_j_per_k += h * d.entropy_j_per_k;
        if (!feed_zero) {  // the resolved molar feed, whatever the declared basis (design r4)
            q.feed = FeedSpec{FeedKind::Amounts, compiled.n_feed, {}, 1.0};
            for (std::size_t i = 0; i < d.feed_amounts_mol.size(); ++i) q.feed.values[i] += h * d.feed_amounts_mol[i];
        }
        const auto neighbor_system = compile_problem(*compiled.model, q);
        const auto neighbor = solve_compiled(neighbor_system, central);
        // Availability first: a failed result keeps phase_count but clears its state arrays.
        const std::string cause = !neighbor.success ? (neighbor.convergence.cause.empty()
            ? core::stop_name(neighbor.convergence.stop) : neighbor.convergence.cause)
            : neighbor.topology_event != TopologyEvent::None ? "topology_event" : "";
        bool same = cause.empty() && neighbor.phase_count == central.phase_count
            && neighbor.molar_densities.size() == central.molar_densities.size()
            && neighbor.solid_amounts.size() == central.solid_amounts.size();
        for (std::size_t i = 0; same && i < central.molar_densities.size(); ++i)
            same = std::abs(std::log(neighbor.molar_densities[i] / central.molar_densities[i])) <= policy::DENSITY_ANCHOR_LOG_WINDOW.value;
        for (std::size_t i = 0; same && i < central.solid_amounts.size(); ++i)
            same = (neighbor.solid_amounts[i] > 0.0) == (central.solid_amounts[i] > 0.0);
        if (!same) {
            for (auto& c : out) { c.status = cause.empty() ? "branch_changed" : "neighbor_unavailable"; c.cause = cause; }
            return out;
        }
        for (const auto& r : solved_state_actions(neighbor_system, neighbor, request).results) values[side].push_back(r.value.value_or(NAN));
    }
    const double eps = std::numeric_limits<double>::epsilon();  // one result per observable at each point
    for (std::size_t k = 0; k < out.size(); ++k)
        out[k] = core::difference_check(out[k].label, exact.results[k].action.value_or(NAN), values[0][k], values[1][k],
            eps * std::abs(values[0][k]), eps * std::abs(values[1][k]), step, tau, kappa);
    return out;
}

}  // namespace epcsaft::equilibrium
