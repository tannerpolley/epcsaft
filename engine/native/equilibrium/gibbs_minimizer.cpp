#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "contributions/helmholtz.hpp"

#include "core/trace.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <IpCompoundVector.hpp>
#include <IpIpoptApplication.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <IpIpoptData.hpp>
#include <IpOrigIpoptNLP.hpp>
#include <IpTNLPAdapter.hpp>
#include <IpoptConfig.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <string>

// Instantiated once in equilibrium/singular_value_decompositions.cpp.
extern template class Eigen::JacobiSVD<Eigen::MatrixXd>;
extern template class Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

namespace epcsaft::equilibrium {
namespace { constexpr double bound_infinity = 1.0e19; }

HelmholtzSpectrum helmholtz_spectrum(const HelmholtzDerivatives& derivative,
    const std::vector<double>& x, const std::vector<std::size_t>& support, const std::vector<double>& normalization) {
    const auto& z = normalization.empty() ? x : normalization;
    const auto c = static_cast<Eigen::Index>(support.size());
    const auto width = derivative.gradient.size();
    HelmholtzSpectrum out;
    out.matrix.resize(c, c);
    for (Eigen::Index a = 0; a < c; ++a)
        for (Eigen::Index b = 0; b < c; ++b)
            out.matrix(a, b) = std::sqrt(z[support[a]] * z[support[b]])
                * derivative.hessian[(2 + support[a]) * width + 2 + support[b]] + (a == b ? z[support[a]] / x[support[a]] : 0.0);
    if (!out.matrix.allFinite()) throw std::domain_error("nonfinite Helmholtz amount Hessian");
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(out.matrix);
    if (eigen.info() != Eigen::Success) throw std::domain_error("Helmholtz eigensolve failed");
    out.values = eigen.eigenvalues(); out.mode = eigen.eigenvectors().col(0);
    out.band = (out.matrix * out.mode - out.values(0) * out.mode).norm()
        + static_cast<double>(c) * std::numeric_limits<double>::epsilon() * out.matrix.norm();
    return out;
}

namespace {

// EqID reduced_hessian_stability on the balance nullspace. `minimum` is the smallest eigenvalue in the solved
// chart. The decision reads the congruent equilibration D = diag(S)^-1/2, S = |objective| + |weighted constraint|
// curvature, whose inertia is the same, against the rounding band eps (n ||DSD|| + 2 kappa(AD) ||DHD||)
// (ALG-EQUILIBRIUM, #159): trace log-chart directions have curvature far below eps ||H||. `jacobian` holds the
// physical rows at `at`, assembled here when empty; `frozen` replaces D to recertify a corrected point (#1064).
struct Curvature final { double minimum = 0.0, equilibrated = 1.0, band = 0.0, magnitude = 0.0; Eigen::VectorXd d; };

Curvature reduced_curvature_at(const CompiledSystem& compiled, const PhasePoint& at, const std::vector<double>& multipliers,
    std::vector<double>& jacobian, const Eigen::VectorXd* frozen = nullptr) {
    if ((compiled.feasibility_cell() && !compiled.problem.volume_m3) || at.point.empty()) return {};
    if (jacobian.empty()) jacobian = compiled.jacobian_at(at);
    auto prepared = at;
    prepared.point.resize(compiled.primal_size);
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto width = static_cast<Eigen::Index>(prepared.point.size());
    const auto n = width - compiled.problem.volume_m3.has_value();
    const auto stride = static_cast<Eigen::Index>(jacobian.size() / compiled.rows.size());
    Eigen::MatrixXd constraints(0, n);
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        const auto& item = compiled.rows[row];  // an absent solid's row n_s = 0 constrains
        if (stationarity_row(item.kind) && !(item.kind == RowKind::SolidAffinity && !compiled.solids[item.index].present)) continue;
        constraints.conservativeResize(constraints.rows() + 1, n);
        constraints.bottomRows(1) = Eigen::Map<const Dense>(jacobian.data() + row * stride, 1, stride).leftCols(n);
    }
    std::vector<double> entropy_magnitude;
    const Eigen::MatrixXd objective = Eigen::Map<const Dense>((compiled.problem.internal_energy_j ? compiled.entropy_hessian_at(prepared, &entropy_magnitude) : compiled.hessian_at(prepared, {}, compiled.problem.volume_m3.has_value())).data(), width, width).topLeftCorner(n, n);
    Eigen::MatrixXd weighted = Eigen::MatrixXd::Zero(n, n);
    if (!compiled.problem.internal_energy_j) weighted = Eigen::Map<const Dense>(
        compiled.constraint_hessian_at(prepared, multipliers, compiled.problem.volume_m3.has_value()).data(), width, width).topLeftCorner(n, n);
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
    Eigen::MatrixXd magnitude = objective.cwiseAbs() + weighted.cwiseAbs();
    if (!entropy_magnitude.empty()) magnitude = Eigen::Map<const Dense>(entropy_magnitude.data(), width, width).topLeftCorner(n,n);
    Eigen::VectorXd d = magnitude.diagonal().unaryExpr([](double s) { return s > 0.0 ? 1.0 / std::sqrt(s) : 1.0; });
    if (frozen) d = *frozen;
    const Eigen::MatrixXd scaled = d.asDiagonal() * hessian * d.asDiagonal();
    const auto [equilibrated, kappa] = reduce(scaled, constraints * d.asDiagonal());
    const double dsd = (d.asDiagonal() * magnitude * d.asDiagonal()).norm(), rounding = static_cast<double>(n) * dsd;
    return {reduce(hessian, constraints).first, equilibrated,
            std::numeric_limits<double>::epsilon() * (rounding + 2.0 * kappa * scaled.norm()), dsd, d};
}

// Extensive-objective constraint weights from the row multipliers of `point`: a material row (B n - b)/N_feed
// carries -N_feed lambda_j, a charge row sum z x carries -n_p phi_p (zero when `moiety_only`), and a free-amount
// phase's normalization row is zero because the objective and material rows obey joint extensive homogeneity in
// amount, composition and volume, together with stationarity in the free amount.
std::vector<double> constraint_weights(const CompiledSystem& compiled, const PhasePoint& at, const std::vector<double>& point,
    bool moiety_only) {
    std::vector<double> out;
    std::size_t electro = 0;
    const auto multiplier = point.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size);
    for (const auto& row : compiled.rows) {
        if (stationarity_row(row.kind)) continue;
        out.push_back(row.kind == RowKind::Moiety ? -compiled.feed_total * multiplier[row.index] : 0.0);
        if (row.kind != RowKind::Electro || moiety_only) continue;
        const auto& free = at.shape[row.phase].amount;
        out.back() = -(free ? point[*free] : at.amount[row.phase]) * multiplier[compiled.rank + electro++];
    }
    return out;
}

// The polish's coordinate scale: a free phase amount, a free temperature and a free pressure by their own values.
template <typename Base>
auto coordinate_scale(const CompiledSystem& compiled, const PhasePointValues<Base>& at, Eigen::Index width) {
    auto scale = Eigen::Matrix<Base, Eigen::Dynamic, 1>::Ones(width).eval();
    for (std::size_t p = 0; p < at.shape.size(); ++p)
        if (at.shape[p].amount && at.amount[p] > 0.0) scale(static_cast<Eigen::Index>(*at.shape[p].amount)) = at.amount[p];
    auto intensive = static_cast<Eigen::Index>(compiled.primal_size - compiled.free_intensive_count);
    if (compiled.problem.temperature.free) scale(intensive++) = at.temperature;
    if (compiled.problem.pressure.free) scale(intensive) = at.pressure;
    return scale;
}

// Local phase stability (ALG-EQUILIBRIUM, #352) of a neutral nonreacting liquid/vapor pair with an incipient
// zero-amount phase: every phase's M = D Q D, Q = F_nn at unit amount as in EqID critical_equations, D = diag
// sqrt(x) on its support, may hold no resolved negative eigenvalue, lambda < -(||M v - lambda v|| + C eps ||M||_F).
// Returns the refusal, empty when admitted or outside that domain.
std::string local_phase_instability(const CompiledSystem& compiled, const PhasePoint& prepared, Candidate& disposition) {
    const auto& phases = compiled.problem.phases;
    if (compiled.has_ions || compiled.has_reaction || (!compiled.problem.volume_m3 && (phases.size() != 2 || phases[0].kind == phases[1].kind
        || std::none_of(phases.begin(), phases.end(), [](const PhaseSpec& phase) { return phase.is_zero_pinned(); })))) return {};
    for (std::size_t phase = 0; phase < phases.size(); ++phase) {
        const auto& support = phases[phase].support;
        const auto& x = prepared.composition[phase];
        HelmholtzSpectrum spectrum;
        try { spectrum = helmholtz_spectrum(prepared.derivative[phase], x, support); }
        catch (const std::bad_alloc&) { throw; }
        catch (const std::exception& error) { return "phase " + phases[phase].name + " spectrum unavailable: " + error.what(); }
        const double lowest = spectrum.values(0), band = spectrum.band;
        if (lowest < -band) { disposition = Candidate::Unstable; return "local phase instability: phase " + phases[phase].name + " has lowest eigenvalue " + core::optional_number(lowest)
            + " below the rounding band -" + core::optional_number(band); }
    }
    return {};
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
        sum += composition[item.index] * std::exp(-residuals[row]);
    }
    return 1.0 - sum;
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
    if (spec.amount.kind == AmountKind::Pinned && start.temperature && start.pressure) {
        const auto branch = branch_of(spec.kind);
        density = solve_density(*compiled.model, *start.temperature, *start.pressure,
                                composition, branch, std::nullopt).molar_density;
    }
    point[cursor] = std::log((compiled.model->parameters().ideal_gas ? 1.0 : compiled.model->packing_fraction_coefficient(
        start.temperature.value_or(compiled.problem.temperature.value),
        composition)) * density);
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
    if (compiled.problem.temperature.free && start.temperature) point[intensive++] = *start.temperature;
    if (compiled.problem.pressure.free && start.pressure) point[intensive++] = *start.pressure;
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
            start.temperature.value_or(compiled.problem.temperature.value),
            composition) * start.molar_densities[phase];
        if (!admitted_packing_fraction(packing)) return false;
    }
    return true;
}

}  // namespace

void require_solution_shape(
    const EquilibriumResult& result, std::size_t phases, std::size_t components
) {
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

static void clear_state(EquilibriumResult& result) {
    result.phase_pressures.clear();
    result.optimizer_state.reset();
    result.prepared_endpoint.reset();
    result.temperature.reset();
    result.pressure.reset();
    result.amounts.clear();
    result.mole_fractions.clear();
    result.molar_densities.clear();
    result.packing_fractions.clear();
    result.multipliers.clear();
    result.constraint_multipliers.clear();
    result.coordinates.clear();
    result.solid_amounts.clear();
    result.solid_affinity_decades.clear();
}

void refuse(EquilibriumResult& result, Candidate as) {
    result.candidate = as;
    result.success = false;
    result.failure = nullptr;
    clear_state(result);
}

std::vector<double> add_multipliers(
    const CompiledSystem& compiled, const std::vector<double>& point, EquilibriumResult* observations,
    const PhasePoint* cached_prepared, const PointEvaluation* cached_rows
) {
    if (point.size() != compiled.primal_size) return point;
    const std::size_t width = compiled.rank + compiled.charge_row_count;
    std::vector<double> result = point;
    result.resize(point.size() + width, 0.0);
    if (width == 0) return result;
    // Recover the balance and charge multipliers by the KKT-weighted least squares of the
    // chemical rows: the free-amount minimizer's stationarity conditions carry the weights
    // N_p x_i / N_feed, so the weighted fit is the one whose classified residual is minimal.
    const auto recomputed = cached_prepared ? PhasePoint{} : compiled.prepare(result, true, observations ? &observations->first_numerical_failure : nullptr);
    const auto& prepared = cached_prepared ? *cached_prepared : recomputed;
    const auto recomputed_rows = cached_rows ? PointEvaluation{} : compiled.evaluate_at(prepared);
    const auto& values = cached_rows ? cached_rows->values : recomputed_rows.values;
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
    const auto factor = matrix.colPivHouseholderQr();
    if (observations) core::retain_rank(observations->rank_records, factor, "weighted_multiplier_lift");
    const Eigen::VectorXd solved = factor.solve(rhs);
    for (std::size_t column = 0; column < width; ++column) {
        result[point.size() + column] = solved(static_cast<Eigen::Index>(column));
    }
    return result;
}

// TNLP endpoint before polish; duals use the normalized objective and feed-relative amounts.
struct OptimizerState final {
    EquilibriumProblem declaration;
    std::vector<std::vector<double>> basis;
    std::vector<double> identity, point, lambda, lower, upper;
};

namespace {
EquilibriumProblem optimizer_declaration(EquilibriumProblem p) {
    p.temperature.value = p.pressure.value = 1.0;
    p.temperature.absent = p.pressure.absent = false;
    p.feed.solvent_mass_kg = 1.0;
    for (auto& phase : p.phases) if (phase.is_pinned()) phase.amount.fraction = phase.is_zero_pinned() ? 0.0 : 1.0;
    std::fill(p.standard_gibbs_j_per_mol.begin(), p.standard_gibbs_j_per_mol.end(), 0.0);
    for (auto& r : p.reactions) {
        if (r.equilibrium_constant) r.equilibrium_constant = 1.0;
        if (r.affinity_j_per_mol) r.affinity_j_per_mol = 0.0;
        if (r.correlation) r.correlation->a = r.correlation->b = r.correlation->c = r.correlation->d = 0.0;
    }
    return p;
}

}

GibbsMinimizer::~GibbsMinimizer() = default;

bool GibbsMinimizer::initialize_warm(const OptimizerState& state, double cap, bool predicted) {
    const auto refuse = [&](const char* cause) { if (scope_) scope_->probe("warm_ineligible", {}, 0, cause); return false; };
    if (compiled_.ideal) return refuse("ideal previous-result starts unsupported");
    if (!(cap > 0.0)) return refuse("explicit options file or nonpositive effective multiplier cap");
    if (state.identity != compiled_.optimizer_identity() || state.basis != compiled_.B
        || !declaration_difference(state.declaration, optimizer_declaration(compiled_.problem)).empty()) return refuse("layout changed");
    const auto n = chart_.size(), m = constraint_rows_.size();
    if (state.point.size() != n || state.lower.size() != n || state.upper.size() != n || state.lambda.size() != m)
        return refuse("tuple shape");
    auto chart = chart_;
    auto start = predicted ? std::vector<double>(start_.begin(), start_.begin() + n) : state.point;
    if (!predicted && !compiled_.feasibility_cell()) for (auto i : amount_columns_)
        chart[i] = std::max(start[i], policy::PHASE_DISAPPEARANCE_AMOUNT.value * compiled_.feed_total);
    std::vector<double> lower(n), upper(n), gl(m), gu(m);
    if (!get_bounds_info(n, lower.data(), upper.data(), m, gl.data(), gu.data())) return refuse("bounds unavailable");
    const double a = 1.0 / compiled_.feed_total, push = policy::EQUILIBRIUM_WARM_PUSH.value;
    double minimum = INFINITY, maximum = 0.0, sum = 0.0, products = 0.0;
    warm_lower_.assign(n, 0.0); warm_upper_.assign(n, 0.0); warm_lambda_.resize(m);
    for (std::size_t j = 0; j < m; ++j) {
        warm_lambda_[j] = state.lambda[j] / a;
        if (!std::isfinite(warm_lambda_[j]) || std::abs(state.lambda[j]) > cap) return refuse("equality multiplier cap or nonfinite");
    }
    for (std::size_t i = 0; i < n; ++i) {
        double y = start[i] / chart[i];
        const bool lo = lower[i] > -bound_infinity, hi = upper[i] < bound_infinity;
        const double width = hi && lo ? upper[i] - lower[i] : INFINITY;
        if (!(width > 0.0)) return refuse("fixed bound");
        if (lo) y = std::max(y, lower[i] + std::min(push * std::max(1.0, std::abs(lower[i])), push * width));
        if (hi) y = std::min(y, upper[i] - std::min(push * std::max(1.0, std::abs(upper[i])), push * width));
        if (!std::isfinite(y) || (lo && !(y > lower[i])) || (hi && !(y < upper[i]))) return refuse("unrepresentable interior");
        const double unit = std::find(amount_columns_.begin(), amount_columns_.end(), i) != amount_columns_.end() ? compiled_.feed_total : 1.0;
        for (int side = 0; side < 2; ++side) {
            if (!(side == 0 ? lo : hi)) continue;
            const double slack = side == 0 ? y - lower[i] : upper[i] - y;
            double z = (side == 0 ? state.lower[i] : state.upper[i]) * chart[i] / unit;
            if (!(z >= 0.0) || !std::isfinite(z)) return refuse("invalid bound multiplier");
            z = std::max(z, push);
            const double product = z * slack;
            if (!std::isfinite(product) || z > cap) return refuse("bound multiplier cap or nonfinite product");
            (side == 0 ? warm_lower_[i] : warm_upper_[i]) = z / a;
            minimum = std::min(minimum, product); maximum = std::max(maximum, product); sum += product; ++products;
        }
        start[i] = y * chart[i];
    }
    std::copy(start.begin(), start.end(), start_.begin());
    chart_ = std::move(chart);
    warm_ = true;
    if (scope_) scope_->probe("warm_initialized", {minimum, products ? sum / products : 0.0, maximum});
    return true;
}

GibbsMinimizer::GibbsMinimizer(
    const CompiledSystem& compiled, std::vector<double> start, core::Scope* scope, const DriverLimits* limits, bool physical_amounts
)
    : compiled_(compiled), start_(std::move(start)), scope_(scope), limits_(limits) {
    for (std::size_t row = 0; row < compiled_.rows.size(); ++row)
        if (compiled_.feasibility_cell() || !stationarity_row(compiled_.rows[row].kind)) constraint_rows_.push_back(row);
    const auto offsets = compiled_.phase_offsets();
    for (std::size_t p = 0; p < offsets.size(); ++p)
        if (!compiled_.ideal && compiled_.problem.phases[p].is_free()) amount_columns_.push_back(offsets[p] + compiled_.problem.phases[p].support.size());
    // ALG-EQUILIBRIUM: Ipopt measures each free amount in its start amount, floored at the disappearance cutoff.
    const bool cell = compiled_.feasibility_cell();
    chart_.assign(cell ? compiled_.unknown_count() : compiled_.primal_size, 1.0);
    if (!cell && !physical_amounts) for (const auto i : amount_columns_)
        chart_[i] = std::max(start_[i], policy::PHASE_DISAPPEARANCE_AMOUNT.value * compiled_.feed_total);
    if (scope_) scope_->probe("amount_chart", {}, 0, physical_amounts ? "physical_mol" : cell ? "physical_cell" : "start_amount");
}

bool GibbsMinimizer::get_nlp_info(
    Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style
) {
    n = static_cast<Ipopt::Index>(chart_.size());
    m = static_cast<Ipopt::Index>(constraint_rows_.size());
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
    const double infinity = bound_infinity;
    std::fill(x_l, x_l + n, -infinity);
    std::fill(x_u, x_u + n, infinity);
    if (compiled_.problem.loading) {
        const auto i = compiled_.loading_column();
        x_l[i] = compiled_.problem.loading->lower;
        x_u[i] = compiled_.problem.loading->upper;
    }
    lower_.assign(static_cast<std::size_t>(n), -infinity);  // the amount step rule's bounds, read at each callback
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
            if (!(compiled_.reaction_relaxed_split && compiled_.model->parameters().charge_number[i] != 0.0))
                box[i] = std::max(box[i], std::min(policy::EQUILIBRIUM_COMPOSITION_BOX.value,
                    policy::EQUILIBRIUM_BOX_INVENTORY_FRACTION.value * compiled_.inventory_upper_bounds[i]));
    }
    if (limits_ && limits_->expired()) return false;
    std::size_t offset = 0;
    for (const auto& phase : compiled_.problem.phases) {
        // ALG-EQUILIBRIUM: a one-species phase's normalization row fixes ln x = 0, so an upper bound there is always
        // active; with unrelaxed bounds it would hold the solution on the barrier.
        const double top = phase.support.size() == 1 ? infinity : 0.0;
        for (const std::size_t i : phase.support) {
            x_l[offset] = std::log(box[i]);
            x_u[offset++] = compiled_.ideal ? std::log(compiled_.inventory_upper_bounds[i]) : top;
        }
        if (compiled_.ideal) continue;
        if (phase.amount.kind == AmountKind::Free) x_l[offset++] = 0.0;
        if (!compiled_.log_packing_unbounded()) {
            x_l[offset] = std::log(policy::PACKING_FRACTION_BOX_FLOOR.value);
            x_u[offset] = std::log(policy::PACKING_FRACTION_CEILING.value);
        }
        if (compiled_.problem.spinodal) {
            const double coefficient = compiled_.model->packing_fraction_coefficient(compiled_.seed_temperature(), compiled_.feed_composition);
            const auto& scan = spinodal_search(compiled_, limits_);
            x_l[offset] = std::max(x_l[offset], std::log(coefficient * scan.lower));
            x_u[offset] = std::min(x_u[offset], std::log(coefficient * scan.upper));
        }
        ++offset;
    }
    for (const auto& solid : compiled_.solids)
        if (solid.column) x_l[offset++] = 0.0;  // EqID solid_amount_complementarity: n_s >= 0, no upper bound
    offset += compiled_.problem.loading.has_value();
    offset += compiled_.mode_count;  // the critical mode u is unbounded
    if (compiled_.problem.temperature.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    if (compiled_.problem.pressure.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    std::fill(g_l, g_l + m, 0.0);
    std::fill(g_u, g_u + m, 0.0);
    return true;
}

bool GibbsMinimizer::physical_bounds(std::vector<double>& lower, std::vector<double>& upper) {
    const auto n = static_cast<Ipopt::Index>(chart_.size());
    const auto m = static_cast<Ipopt::Index>(constraint_rows_.size());
    std::vector<double> gl(m), gu(m);
    if (!get_bounds_info(n, lower.data(), upper.data(), m, gl.data(), gu.data())) return false;
    for (std::size_t i = 0; i < chart_.size(); ++i) {
        lower[i] = lower[i] <= -bound_infinity ? -INFINITY : lower[i] * chart_[i];
        upper[i] = upper[i] >= bound_infinity ? INFINITY : upper[i] * chart_[i];
    }
    return true;
}

bool GibbsMinimizer::get_starting_point(
    Ipopt::Index n, bool init_x, Ipopt::Number* x,
    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda
) {
    if (!init_x && !warm_) return false;
    const std::size_t size = !compiled_.feasibility_cell() ? compiled_.primal_size : start_.size();
    for (Ipopt::Index i = 0; init_x && i < n; ++i) x[i] = (static_cast<std::size_t>(i) < size ? start_[i] : 0.0) / chart_[i];
    if (init_z) {
        if (warm_) std::copy(warm_lower_.begin(), warm_lower_.end(), z_L); else std::fill(z_L, z_L + n, 0.0);
        if (warm_) std::copy(warm_upper_.begin(), warm_upper_.end(), z_U); else std::fill(z_U, z_U + n, 0.0);
    }
    if (init_lambda) { if (warm_) std::copy(warm_lambda_.begin(), warm_lambda_.end(), lambda); else std::fill(lambda, lambda + m, 0.0); }
    return true;
}

static_assert(IPOPT_VERSION_MAJOR == 3 && IPOPT_VERSION_MINOR == 11 && IPOPT_VERSION_RELEASE == 9,
    "the amount step rule reads Ipopt 3.11.9 internals (ALG-EQUILIBRIUM): requalify it on an Ipopt upgrade");

// Ipopt's current iterate in TNLP order; in restoration, the original block of the restoration iterate.
// Cached by Ipopt's vector tag, so a callback at an unchanged base costs one tag comparison.
const std::vector<double>& GibbsMinimizer::search_base() {
    const Ipopt::Vector* current = Ipopt::GetRawPtr(search_->curr()->x());
    if (const auto* compound = dynamic_cast<const Ipopt::CompoundVector*>(current)) current = Ipopt::GetRawPtr(compound->GetComp(0));
    if (base_.empty() || current->GetTag() != base_tag_) {
        base_.resize(lower_.size());
        adapter_->ResortX(*current, base_.data());
        base_tag_ = current->GetTag();
    }
    return base_;
}

bool GibbsMinimizer::step_admitted(const Ipopt::Number* x) {
    if (limits_) limits_->check();
    // ALG-EQUILIBRIUM: a trial point keeps at least kappa_N of each free amount's distance to its relaxed lower bound.
    if (Ipopt::IsValid(search_)) {
        const auto& base = search_base();
        bool refused = false;
        for (const std::size_t i : amount_columns_)
            refused = refused || (base[i] > lower_[i]
                && x[i] - lower_[i] < policy::EQUILIBRIUM_AMOUNT_STEP_RETENTION.value * (base[i] - lower_[i]));
        if (refused) return false;
    }
    return true;
}

const PhasePoint& GibbsMinimizer::prepared(Ipopt::Index n, const Ipopt::Number* x) {
    if (!prepared_.has_value() || !std::equal(x, x + n, prepared_->point.begin(), prepared_->point.end(),
        [&, i = std::size_t{0}](double value, double stored) mutable { return value * chart_[i++] == stored; })) {
        std::vector<double> candidate(x, x + n);
        for (Ipopt::Index i = 0; i < n; ++i) candidate[i] *= chart_[i];
        prepared_ = compiled_.prepare(candidate, true, &observations_.first_numerical_failure);
        rows_.reset();
        jacobian_.reset();
        gradient_.reset();
        objective_hessian_.reset();
        constraint_hessian_.reset();
        ++prepared_points_;
        if (limits_ && limits_->admit) if (const auto refused = limits_->admit(*prepared_)) {
            prepared_.reset();
            throw std::domain_error(*refused);
        }
    }
    if (limits_) limits_->check();
    return *prepared_;
}

const PointEvaluation& GibbsMinimizer::rows_at(Ipopt::Index n, const Ipopt::Number* x) {
    const auto& point = prepared(n, x);
    if (!rows_.has_value()) rows_ = compiled_.evaluate_at(point);
    return *rows_;
}

bool GibbsMinimizer::eval_f(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& value
) {
    core::Scope callback("ipopt_eval_f");
    if (limits_ && !limits_->spend(limits_->evaluations, limits_->evaluation_cap)) return false;
    ++evaluations_;
    try {
        if (!compiled_.feasibility_cell() && !step_admitted(x)) return false;
        value = !compiled_.feasibility_cell() ? rows_at(n, x).objective : 0.0;
        return std::isfinite(value);
    } catch (const std::exception& error) {
        retain_callback_failure(error, "eval_f");
        return false;
    }
}

bool GibbsMinimizer::eval_grad_f(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f
) {
    core::Scope callback("ipopt_eval_grad_f");
    try {
        if (!step_admitted(x)) return false;
        const auto& point = prepared(n, x);
        if (!gradient_.has_value()) gradient_ = compiled_.objective_gradient_at(point);
        for (Ipopt::Index i = 0; i < n; ++i) grad_f[i] = (*gradient_)[i] * chart_[i];
        return true;
    } catch (const std::exception& error) {
        retain_callback_failure(error, "eval_grad_f");
        return false;
    }
}

bool GibbsMinimizer::eval_g(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Index m, Ipopt::Number* g
) {
    core::Scope callback("ipopt_eval_g");
    try {
        if (!step_admitted(x)) return false;
        const auto& values = rows_at(n, x).values;
        const auto& rows = constraint_rows_;
        for (Ipopt::Index row = 0; row < m; ++row) g[row] = values[rows[row]];
        return true;
    } catch (const std::exception& error) {
        retain_callback_failure(error, "eval_g");
        return false;
    }
}

bool GibbsMinimizer::eval_jac_g(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    core::Scope callback("ipopt_eval_jac_g");
    const auto& rows = constraint_rows_;
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
        if (!step_admitted(x)) return false;
        const auto& point = prepared(n, x);
        if (!jacobian_.has_value()) jacobian_ = compiled_.jacobian_at(point);
        Ipopt::Index entry = 0;
        for (const std::size_t row : rows) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                values[entry++] = (*jacobian_)[row * static_cast<std::size_t>(n) + static_cast<std::size_t>(column)] * chart_[column];
            }
        }
        return entry == nele_jac;
    } catch (const std::exception& error) {
        retain_callback_failure(error, "eval_jac_g");
        return false;
    }
}

bool GibbsMinimizer::eval_h(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
    bool /*new_lambda*/, Ipopt::Index nele_hess, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    core::Scope callback("ipopt_eval_h");
    if (limits_ && !limits_->spend(limits_->hessian_calls, limits_->hessian_cap)) return false;
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
        if (!step_admitted(x)) return false;
        const auto& point = prepared(n, x);
        // Objective curvature is multiplier independent; the weighted constraint curvature is
        // recombined whenever the optimizer multipliers change at the same point.
        if (!objective_hessian_.has_value()) objective_hessian_ = compiled_.hessian_at(point, {});
        const std::vector<double> multipliers(lambda, lambda + m);
        if (!constraint_hessian_.has_value() || constraint_hessian_multipliers_ != multipliers) {
            constraint_hessian_ = compiled_.constraint_hessian_at(point, multipliers);
            constraint_hessian_multipliers_ = multipliers;
        }
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < n; ++row) {
            for (Ipopt::Index column = 0; column <= row; ++column) {
                const std::size_t index = static_cast<std::size_t>(row) * n
                    + static_cast<std::size_t>(column);
                values[entry++] = (obj_factor * (*objective_hessian_)[index]
                    + (*constraint_hessian_)[index]) * chart_[row] * chart_[column];
            }
        }
        return entry == nele_hess;
    } catch (const std::exception& error) {
        retain_callback_failure(error, "eval_h");
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
    core::Scope callback("ipopt_intermediate");
    observe_ipopt(*observations_.ipopt_diagnostics, scope_, {int(mode), iter, ls_trials, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du, alpha_pr});
    restored_ = restored_ || mode == Ipopt::RestorationPhaseMode;
    // Ipopt adjusts bounds only when it accepts a point, before this callback, so the active problem's bound read
    // here holds for every search until the next callback; in restoration that is restoration's own bound block, whose
    // original-variable part restoration does not update when it adjusts the original problem's bounds.
    if (!amount_columns_.empty()) {  // without free amounts the rule keeps nothing: no search base, no bound reads
        Ipopt::IpoptNLP& active = *ip_cq->GetIpoptNLP();
        if (adapter_ == nullptr) {  // the first callback is regular mode, before any search
            auto* original = dynamic_cast<Ipopt::OrigIpoptNLP*>(&active);
            adapter_ = original ? dynamic_cast<Ipopt::TNLPAdapter*>(Ipopt::GetRawPtr(original->nlp())) : nullptr;
            if (adapter_ == nullptr) step_rule_error_ = "amount step rule: Ipopt's TNLP adapter is not reachable";
        }
        if (!step_rule_error_ && active.NLP_scaling()->have_x_scaling())
            step_rule_error_ = "amount step rule: Ipopt scales the variables, so its iterate and bounds are not the callback coordinates";
        if (step_rule_error_) return false;
        Ipopt::SmartPtr<const Ipopt::Vector> lower = active.x_L();
        if (const auto* compound = dynamic_cast<const Ipopt::CompoundVector*>(Ipopt::GetRawPtr(lower))) lower = compound->GetComp(0);
        if (lower->GetTag() != lower_tag_) {  // Ipopt relaxed or adjusted the bounds since the last callback
            adapter_->ResortBnds(*lower, lower_.data(), *lower, nullptr);
            lower_tag_ = lower->GetTag();
        }
        search_ = ip_data;
    }
    iterations_ = static_cast<std::size_t>(iter);
    if (limits_) limits_->consumed = std::max(limits_->consumed, consumed_before_ + iterations_);
    regularization_size_ = regularization_size;
    if (warm_ && restored_) core::retain_failure(observations_.first_numerical_failure,
        "warm attempt stopped on restoration", "Ipopt callback");
    return !(warm_ && restored_) && (!limits_ || !limits_->expired());  // max_iter owns iteration termination
}

void GibbsMinimizer::finalize_solution(
    Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
    const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
    const Ipopt::Number* /*g*/, const Ipopt::Number* lambda,
    Ipopt::Number /*objective*/,
    const Ipopt::IpoptData* /*ip_data*/, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
) {
    if (status != Ipopt::SUCCESS && observations_.ipopt_diagnostics) {
        auto& d = *observations_.ipopt_diagnostics;
        IpoptPoint point;
        point.tnlp_coordinates.assign(x, x + n);
        point.tnlp_lower_multipliers.assign(z_L, z_L + n);
        point.tnlp_upper_multipliers.assign(z_U, z_U + n);
        point.tnlp_constraint_multipliers.assign(lambda, lambda + m);
        point.coordinate_scale = chart_; point.constraint_rows = constraint_rows_;
        point.declaration = std::make_shared<const EquilibriumProblem>(compiled_.problem);
        point.rows = compiled_.rows; point.objective_factor = 1.0 / compiled_.feed_total;
        d.failed_point = std::move(point);
    }
    status_ = status;
    solution_.assign(x, x + n);
    for (Ipopt::Index i = 0; i < n; ++i) solution_[i] *= chart_[i];
    constraint_multipliers_.assign(lambda, lambda + m);
    lower_multipliers_.assign(z_L, z_L + n);
    for (Ipopt::Index i = 0; i < n; ++i) lower_multipliers_[i] /= chart_[i];
    if (compiled_.ideal || status != Ipopt::SUCCESS || restored_) return;
    auto endpoint = std::make_shared<OptimizerState>();
    endpoint->declaration = optimizer_declaration(compiled_.problem);
    endpoint->basis = compiled_.B; endpoint->identity = compiled_.optimizer_identity();
    endpoint->point = solution_; endpoint->lambda = constraint_multipliers_;
    const double a = 1.0 / compiled_.feed_total;
    for (auto& v : endpoint->lambda) v *= a;
    for (Ipopt::Index i = 0; i < n; ++i) {
        const double unit = std::find(amount_columns_.begin(), amount_columns_.end(), i) != amount_columns_.end() ? compiled_.feed_total : 1.0;
        endpoint->lower.push_back(a * unit / chart_[i] * z_L[i]);
        endpoint->upper.push_back(a * unit / chart_[i] * z_U[i]);
    }
    if (all_numeric(endpoint->point) && all_numeric(endpoint->lambda) && all_numeric(endpoint->lower) && all_numeric(endpoint->upper)
        && std::all_of(endpoint->lower.begin(), endpoint->lower.end(), [](double v) { return v >= 0.0; })
        && std::all_of(endpoint->upper.begin(), endpoint->upper.end(), [](double v) { return v >= 0.0; })) endpoint_ = std::move(endpoint);

}

namespace {

void reject_elapsed(EquilibriumResult& result, core::Scope& scope) {
    auto& record = result.convergence;
    if (!record.criterion) core::stopping_criterion(record, result.declaration->solver_options, policy::row);
    record.stop = core::classify(record.residual, *record.criterion, record.floor,
                                record.contraction, true, record.amplification);
    refuse(result); scope.finish(record);
}

}  // namespace

EquilibriumResult structural_result(const CompiledSystem& compiled, core::Scope& scope) {
    EquilibriumResult result;
    result.model_identity = compiled.model;
    result.ideal_identity = compiled.ideal;
    result.declaration = std::make_shared<const EquilibriumProblem>(compiled.problem);
    result.rows = compiled.rows;
    result.reference_diagnostics = compiled.reference_diagnostics;
    result.phase_count = compiled.problem.phases.size();
    result.unknown_count = compiled.unknown_count();
    result.jacobian_rank = compiled.jacobian_rank;
    result.last_rank = compiled.last_rank;
    result.rank_records = compiled.rank_records;
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

static void capture_driver(const GibbsMinimizer& nlp, const EquilibriumResult& result, EquilibriumResult& record, const core::AssociationObservations& observations) {
    observations.merge(record.rank_records);
    nlp.publish_observations(record);
    record.attempt = result.attempt;
    record.retry = result.retry;
    record.solver_settings = result.solver_settings;
    record.diagnostic_options_file = result.diagnostic_options_file;
    record.linear_solver = result.linear_solver;
    record.solver_status = solver_status_name(nlp.status());
    record.convergence.iterations = record.iterations = nlp.iterations();
    record.convergence.evaluations = record.evaluations = nlp.evaluations();
    record.hessian_calls = nlp.hessian_calls();
    if (nlp.last_rows()) {
        record.residuals = nlp.last_rows()->values;
        record.classified_residuals = nlp.last_classified();
        if (nlp.last_rows()->floor > 0.0) record.convergence.floor = nlp.last_rows()->floor;
    }
    record.convergence.residual = record.classified_residuals.empty()
        ? std::numeric_limits<double>::quiet_NaN() : maximum_abs(record.classified_residuals);
}

static bool mixture_precision_eligible(const CompiledSystem& c) {
    return !c.ideal && newton_eligible(c) && !c.problem.temperature.free && c.problem.pressure.free
        && !c.problem.source_reference && c.solids.empty() && c.model->parameters().component_count() > 1;
}

std::vector<IpoptIterationRecord>*& ipopt_iteration_log() noexcept {
    static thread_local std::vector<IpoptIterationRecord>* log = nullptr;
    return log;
}

// ALG-TERMINATION: the leading sentence of a topology rejection states what check_topology observed, at the
// declared intensives; it claims neither global stability nor a critical point (#976).
static std::string topology_sentence(const CompiledSystem& compiled, const TopologyCheck& topology) {
    const auto& p = compiled.problem;
    const auto number = [](double value) { char text[32]; std::snprintf(text, sizeof text, "%.6g", value); return std::string(text); };
    if (topology.event == TopologyEvent::PhaseDisappeared) return "phase " + p.phases[*topology.phase].name
        + " vanished: its amount reached at most PHASE_DISAPPEARANCE_AMOUNT of the feed; the declared split was rejected";
    const auto& a = p.phases[topology.pair->first];
    const auto& b = p.phases[topology.pair->second];
    const std::string pair = "declared phases " + a.name + " and " + b.name + " coalesced",
        distance = " (phase distance " + number(topology.distance) + " at most PHASE_COINCIDENCE_DISTANCE)";
    if (p.phases.size() == 2 && a.kind != b.kind && a.is_zero_pinned() != b.is_zero_pinned() && p.temperature.free != p.pressure.free)
        return std::string("no distinct ") + ((a.is_zero_pinned() ? a : b).kind == PhaseKind::Vapor ? "bubble" : "dew")
            + " point was found at " + (p.pressure.free ? "T = " + number(p.temperature.value) + " K" : "P = " + number(p.pressure.value) + " Pa")
            + ": " + pair + distance + "; a critical limit near these conditions is one possible, untested explanation";
    std::string at;
    if (!p.temperature.free) at += " at T = " + number(p.temperature.value) + " K";
    if (!p.pressure.free) at += (at.empty() ? " at P = " : ", P = ") + number(p.pressure.value) + " Pa";
    return pair + at + distance + "; the declared split was rejected; Discover() searches other phase sets";
}

EquilibriumResult publish_candidate(const CompiledSystem& compiled, const PhasePoint& prepared,
    EquilibriumResult result, core::Scope& scope, const DriverLimits* limits) {
    const auto* model = compiled.model;
    const auto& full = prepared.point;
    const double row_tolerance = compiled.problem.solver_options.value(policy::row);
    const auto expired = [&] { return limits && limits->expired(); };
    std::vector<double> volumes;
    const double temperature = prepared.temperature;
    bool topology_ok = true;
    std::string topology_failure, cause;
    // One Engine verdict: a rejected point that classified as converged stops as infeasible; other stops stand.
    const auto infeasible = [&result](const char* reason) {
        auto& stop = result.convergence.stop;
        if (stop != core::Stop::converged && stop != core::Stop::converged_to_floor) return;
        stop = core::Stop::infeasible;
        result.convergence.cause = reason;
    };
    // Preparation above admitted every composition with a positive total.
    const auto topology = check_topology(compiled, full);
    const auto& compositions = topology.compositions;
    const auto& amounts = topology.amounts;
    result.topology_event = topology.event;
    result.topology_phase = topology.phase;
    result.topology_pair = topology.pair;
    for (std::size_t phase = 0; !compiled.ideal && phase < compiled.problem.phases.size(); ++phase) {
        volumes.push_back(model->parameters().ideal_gas ? prepared.volume[phase] / sum(prepared.composition[phase])
            : model->packing_fraction_coefficient(temperature, compositions[phase]) / std::exp(topology.log_packing[phase]));
        const double deficit = michelsen_deficit(
            compiled, compositions[phase], phase, result.residuals);
        scope.probe("michelsen_deficit", deficit, static_cast<std::uint32_t>(phase));
        // The deficit is a probe only: a vanished phase's coordinates are undetermined.
        if (topology.event == TopologyEvent::PhaseDisappeared && topology.phase == phase) {
            topology_ok = false;
            cause = "phase_disappeared";
        }
    }
    for (const double distance : topology.distances) scope.probe("phase_distance", distance);
    if (topology.event == TopologyEvent::PhasesCoalesced) {
        topology_ok = false;
        cause = "phases_coalesced";
        topology_failure = "phases_coincide: phase distance between "
            + compiled.problem.phases[topology.pair->first].name + " and "
            + compiled.problem.phases[topology.pair->second].name + " is too small: "
            + core::optional_number(topology.distance);
    }
    // The recomputed rows classify: a criterion-meeting point is accepted after any non-error
    // stop (error stops never classify as converged); a floor-limited point needs success.
    std::string present;
    for (std::size_t s = 0; s < compiled.solids.size(); ++s)
        if (compiled.solids[s].present) present += (present.empty() ? "" : ", ") + compiled.problem.solids[s].name;
    std::string violations;
    for (std::size_t i = 0; i < compiled.solids.size(); ++i) {
        const auto& solid = compiled.solids[i];
        if (solid.column && !(solid.present ? full[*solid.column] > 0.0 : compiled.solid_affinity(prepared, i) <= row_tolerance))
            violations += (violations.empty() ? "" : ", ") + compiled.problem.solids[i].name;
    }
    if (!violations.empty()) {
        topology_ok = false;
        cause = "solid_classification";
        topology_failure = "solid classification violated after one flip round: " + violations;
    } else if (!present.empty()) {
        const auto info = inspect_sensitivity(compiled, full);
        result.last_rank = core::NumericalRank{"canonical", full.size(), info.rank, info.threshold};
        if (!info.nonsingular) {
            topology_ok = false;
            result.convergence.stop = core::Stop::singular;
            cause = "present_solid_rank";
            topology_failure = "rank-deficient classified equations with present solids: " + present;
        }
    }
    // One full-width physical row Jacobian at the returned point and one LU of its numerical rows serve the
    // charge rank, the Newton correction and the certificate.
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto n = static_cast<Eigen::Index>(full.size());
    std::vector<double> physical;
    Eigen::FullPivLU<Matrix> numerical;
    const auto factor = [&] { if (physical.empty()) physical = compiled.jacobian_at(prepared); numerical.compute(Eigen::Map<const Matrix>(compiled.numerical_rows(prepared, physical, full.size()).data(), n, n)); };
    if (!compiled.log_charge_rows.empty()) {
        factor();
        const auto raw = Eigen::Map<const Matrix>(physical.data(), n, n).fullPivLu();
        const auto physical_rank = raw.rank();
        core::retain_rank(result.rank_records, core::conditioning(raw, "raw_charge_endpoint", 0, "endpoint"), full);
        result.last_rank = core::conditioning(numerical, "log_charge", 0, "endpoint");
        const auto numerical_rank = static_cast<Eigen::Index>(result.last_rank->rank);
        scope.probe("charge_endpoint_ranks", {static_cast<double>(physical_rank), static_cast<double>(numerical_rank)});
        if (numerical_rank != n) { topology_ok = false; cause = "log_charge_rank"; topology_failure = "rank-deficient log-charge equations"; }
    }
    bool accepted = !expired() && result.candidate == Candidate::Admitted && topology_ok;
    Candidate disposition = result.candidate == Candidate::Admitted ? Candidate::Refused : result.candidate;
    if (accepted) {
        topology_failure = local_phase_instability(compiled, prepared, disposition);
        topology_ok = accepted = topology_failure.empty();
        cause = disposition == Candidate::Unstable ? "local_phase_instability" : "phase_spectrum_unavailable";
    }
    if (!accepted) {
        if (!topology_ok) {
            infeasible(cause.c_str());
            result.message = cause == "phases_coalesced" ? topology_sentence(compiled, topology) + "; " + topology_failure
                : cause == "phase_disappeared" ? topology_sentence(compiled, topology) : topology_failure;
        }
        refuse(result, disposition);
        return result;
    }
    if (result.requested_tolerance_met.has_value() && !*result.requested_tolerance_met) {
        result.message += "; requested row tolerance " + core::optional_number(row_tolerance)
            + " not met by an available floor-limited state";
    }
    result.candidate = Candidate::Published; result.success = true;
    result.prepared_endpoint = std::make_shared<const PhasePoint>(prepared);
    result.temperature = temperature;
    result.pressure = prepared.pressure;
    if (!prepared.standard.reference_diagnostics.empty()) {
        result.reference_diagnostics = prepared.standard.reference_diagnostics;
    }
    result.amounts = amounts;
    result.molar_densities.reserve(volumes.size());
    result.packing_fractions.reserve(volumes.size());
    for (std::size_t phase = 0; phase < volumes.size(); ++phase) {
        result.molar_densities.push_back(1.0 / volumes[phase]);
        result.packing_fractions.push_back(
            model->packing_fraction_coefficient(temperature, compositions[phase]) / volumes[phase]);
        result.mole_fractions.insert(result.mole_fractions.end(),
                                     compositions[phase].begin(), compositions[phase].end());
    }
    if (compiled.ideal) result.mole_fractions = compositions[0];
    result.multipliers.assign(full.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size), full.end());
    // Constraint multipliers at the returned (polished) point; the raw optimizer values stay a trace probe.
    const bool moiety_only = compiled.problem.volume_m3.has_value();
    std::vector<double> physical_weights;
    if (!compiled.feasibility_cell() || moiety_only) physical_weights = constraint_weights(compiled, prepared, full, moiety_only);
    if (!compiled.feasibility_cell()) result.constraint_multipliers = physical_weights;
    result.coordinates = full;
    if (expired()) { result.message = "elapsed work budget exhausted"; reject_elapsed(result, scope); return result; }
    if (compiled.ideal) {
        const auto info = inspect_sensitivity(compiled, full);
        result.jacobian_rank = info.rank;
        result.last_rank = core::NumericalRank{"canonical", full.size(), info.rank, info.threshold};
        if (!info.nonsingular) {
             result.message = "rank-deficient ideal endpoint";
            result.convergence.stop = core::Stop::singular; refuse(result); return result;
        }
    }
    if (compiled.problem.internal_energy_j) for (std::size_t p = 0; p < amounts.size(); ++p) {
        const auto cv = model->state(temperature, 1.0 / volumes[p], compositions[p]).total_isochoric_heat_capacity;
        if (!cv || !(*cv > 0.0)) throw std::domain_error("UV phase has no positive isochoric heat capacity");
    }
    const bool caloric = compiled.has_caloric_condition() && !compiled.problem.internal_energy_j;
    const bool restricted = caloric && !(model->parameters().component_count() == 1 && compiled.problem.phases.size() == 2);
    if (caloric && !restricted) {
        if (!(amounts[0] > 0.0 && amounts[1] > 0.0) || volumes[0] == volumes[1]) {

            infeasible("caloric_split");
            result.message = "caloric coexistence requires an interior distinct split";
            refuse(result);
            return result;
        }
        result.message += "; coexistence equality root; fixed-TP Gibbs split is analytically flat";
    }
    // The certificate at the returned point, or at a moved point under frozen D. A caloric mixture certifies its
    // fixed-T,P Gibbs problem, rebuilt by one construction at either point.
    const auto certify = [&](const PhasePoint& at, const std::vector<double>& point, EquilibriumResult* observed,
        const Eigen::VectorXd* frozen) {
        std::vector<double> jacobian;
        if (!restricted) return reduced_curvature_at(compiled, at, frozen ? constraint_weights(compiled, at, point, moiety_only) : physical_weights, frozen ? jacobian : physical, frozen);
        auto diagnostic = compiled.problem;
        diagnostic.enthalpy_j.reset();
        diagnostic.entropy_j_per_k.reset();
        diagnostic.temperature = FreeValue::fixed(at.temperature);
        const auto tp = compile_problem(*model, diagnostic);
        auto* failure = observed ? &observed->first_numerical_failure : nullptr;
        auto point_tp = tp.generic_point();
        const auto layout = tp.prepare(point_tp, false, failure).shape;
        const auto phase_amounts = frozen ? check_topology(compiled, point).amounts : amounts;
        for (std::size_t p = 0; p < layout.size(); ++p) {
            for (std::size_t i = 0; i < layout[p].x.size(); ++i) point_tp[layout[p].x[i]] = point[at.shape[p].x[i]];
            point_tp[layout[p].volume] = point[at.shape[p].volume];
            point_tp[*layout[p].amount] = phase_amounts[p];
        }
        point_tp.resize(tp.primal_size);
        point_tp = add_multipliers(tp, point_tp, observed);
        const auto at_tp = tp.prepare(point_tp, false, failure);
        return reduced_curvature_at(tp, at_tp, constraint_weights(tp, at_tp, point_tp, true), jacobian, frozen);
    };
    const Curvature curvature = certify(prepared, full, &result, nullptr);
    // Solve error (ALG-EQUILIBRIUM): the next Newton correction d = -J^-1 F of the returned point moves the
    // certificate by e = |lambda~(x + d; D(x)) - lambda~(x)|, recertified only where |lambda~| <= tau + G.
    const bool decides = !compiled.feasibility_cell() || compiled.has_caloric_condition() || compiled.problem.volume_m3;
    double change = NAN, estimate = 0.0, step = 0.0;
    std::string refused;
    if (decides && curvature.d.size() != 0 && !std::isnan(curvature.equilibrated)) {
        if (numerical.rows() != n) factor();  // the row system is square: rows match the KKT coordinates
        const auto rows = compiled.numerical_rows(prepared, result.residuals);
        const Eigen::VectorXd delta = numerical.solve(-Eigen::Map<const Eigen::VectorXd>(rows.data(), n));
        const auto primal = static_cast<Eigen::Index>(compiled.primal_size);
        const Eigen::VectorXd scaled = delta.head(primal).cwiseQuotient(coordinate_scale(compiled, prepared, n).head(primal));
        step = delta.allFinite() ? scaled.lpNorm<Eigen::Infinity>() : INFINITY;
        auto moved = full;
        for (Eigen::Index i = 0; i < n; ++i) moved[static_cast<std::size_t>(i)] += delta(i);
        const auto weights = constraint_weights(compiled, prepared, moved, moiety_only);
        const auto size = static_cast<Eigen::Index>(physical_weights.size());  // relative weight step; none for caloric or (U,V)
        const Eigen::Map<const Eigen::VectorXd> before(physical_weights.data(), size), after(weights.data(), size);
        const double norm = before.lpNorm<Eigen::Infinity>(), shift = (after - before).lpNorm<Eigen::Infinity>();
        estimate = curvature.magnitude * (step + (norm > 0.0 && !restricted && !compiled.problem.internal_energy_j ? shift / norm : 0.0));
        if (!(std::abs(curvature.equilibrated) > curvature.band + estimate)) {
            change = std::isinf(step) ? INFINITY : moved == full ? 0.0 : NAN;
            if (std::isnan(change)) try {
                change = std::abs(curvature.equilibrated - certify(compiled.prepare(moved), moved, nullptr, &curvature.d).equilibrated);
            } catch (const std::bad_alloc&) { throw; }
            catch (const std::exception& error) { refused = std::string("; corrected point refused: ") + error.what(); }
            if (std::isnan(change)) change = INFINITY;
        }
        scope.probe("reduced_hessian_correction", {change, estimate, step, std::isnan(change) ? 0.0 : 1.0});
    }
    result.reduced_hessian_minimum = curvature.minimum;
    scope.probe("reduced_hessian_equilibrated", curvature.equilibrated);
    scope.probe("reduced_hessian_band", curvature.band);
    const double margin = curvature.band + (std::isnan(change) ? 0.0 : change);
    if (decides && !(curvature.equilibrated > margin)) {
        // Inside the band the inertia is unresolved at working precision: a typed singular outcome, not a noise sign.
        const bool singular = curvature.equilibrated >= -margin;
        if (singular || std::isnan(curvature.equilibrated))
            result.convergence.stop = singular ? core::Stop::singular : core::Stop::nonfinite;
        infeasible("reduced_hessian_not_positive_definite");

        result.message = std::string("reduced_hessian_stability ") + (singular ? "is indeterminate" : "is not positive definite")
            + ": equilibrated minimum " + core::optional_number(curvature.equilibrated)
            + ", rounding band " + core::optional_number(curvature.band) + ", solve-error change " + (std::isnan(change)
                ? "not evaluated (estimate " + core::optional_number(estimate) : core::optional_number(change) + " (Newton correction " + core::optional_number(step)) + ")" + refused;
        refuse(result);
        return result;
    }
    if (expired()) { result.message = "elapsed work budget exhausted";
        reject_elapsed(result, scope); return result; }
    scope.probe("returned_temperature", *result.temperature);
    scope.probe("returned_pressure", *result.pressure);
    scope.probe("returned_amounts", result.amounts);
    scope.probe("returned_mole_fractions", result.mole_fractions);
    scope.probe("returned_densities", result.molar_densities);
    for (std::size_t s = 0; s < compiled.solids.size(); ++s) {  // EqID solid_saturation_index
        result.solid_amounts.push_back(compiled.solids[s].column ? full[*compiled.solids[s].column] : 0.0);
        result.solid_affinity_decades.push_back(compiled.solid_affinity(prepared, s) / std::log(10.0));
    }
    if (!all_numeric(result.residuals) || result.rows.empty()) {
        throw std::runtime_error("equilibrium result has invalid residuals");
    }
    require_solution_shape(result, compiled.problem.phases.size(),
                           compiled.n_feed.size());
    return result;
}

template <typename Base>
static bool inside_correction_box(const CompiledSystem& compiled, std::vector<Base>& candidate,
    const std::vector<double>& lower, const std::vector<double>& upper, bool free_amount_newton) {
        // The same four-ulp x=1 snap serves both Newton domains.
        const auto offsets = compiled.phase_offsets();
        for (std::size_t p = 0; p < offsets.size(); ++p) for (std::size_t j = 0; j < compiled.problem.phases[p].support.size(); ++j) {
            const auto i = offsets[p] + j;
            double xmax = 1.0;
            for (int ulp = 0; ulp < 4; ++ulp) xmax = std::nextafter(xmax, INFINITY);
            if (candidate[i] > upper[i] && candidate[i] <= wide_math::log(xmax)) candidate[i] = upper[i];
        }
        bool inside = true;
        for (std::size_t i = 0; i < candidate.size(); ++i)
            inside = inside && (!free_amount_newton || wide_math::isfinite(candidate[i])) && candidate[i] >= lower[i] && candidate[i] <= upper[i];
        return inside;
}

template <typename Base, typename Box, typename Branch, typename Solids>
static bool polish_correction(const CompiledSystem& compiled, std::vector<Base>& full,
    EquilibriumResult& result, core::Scope& scope, const DriverLimits* limits,
    bool newton, bool refine_trial, bool non_error, const Box& inside_box, const Branch& trial_branch,
    const Solids& solid_violations, std::size_t& polish_evaluations, bool& corrected, bool& capped,
    core::Stop& initial_branch_stop) {
    const auto* model = compiled.model;
    const auto expired = [&] { return limits && limits->expired(); };
    const auto row_tolerance = compiled.problem.solver_options.value(policy::row);
    const bool free_amount_newton = newton && !compiled.feasibility_cell();
    const bool guarded_newton = free_amount_newton || !std::is_same_v<Base, double>;
    const bool preserve_trial_branch = refine_trial && !compiled.log_charge_rows.empty();
    const auto numerical = [&](const auto& at, auto rows, std::size_t width = 1, const PhasePointValues<Base>* sizes = nullptr) {
        if constexpr (std::is_same_v<Base, double>) return compiled.numerical_rows(at, std::move(rows), width, sizes);
        else return compiled.trace_moiety_relative(at, std::move(rows), width, sizes);
    };
    const auto classified_rows = [&](const auto& at, const auto& rows) {
        if constexpr (std::is_same_v<typename std::decay_t<decltype(rows)>::value_type, double>) return compiled.classified_residuals(at, rows);
        else return compiled.trace_moiety_relative(at, rows);
    };
    using Vector = Eigen::Matrix<Base, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic>;
    std::size_t polish_round = 0;
    while (non_error && !expired()) {
        if (limits && !limits->spend(limits->evaluations, limits->evaluation_cap)) return false;
        auto at = compiled.prepare_values(full, true, &result.first_numerical_failure);
        ++polish_evaluations;
        if (expired()) return false;
        auto rows_before = compiled.evaluate_at(at).values;
        if constexpr (std::is_same_v<Base, double>) if (preserve_trial_branch && !trial_branch(at, full)) {
            if (expired()) return false;
            initial_branch_stop = result.convergence.stop;
            break;
        }
        std::vector<Base> best = full;
        Base best_classified = maximum_abs(classified_rows(at, rows_before));
        // Past the step budget only raw rows above the caller criterion continue, under a separate polish cap,
        // and only while the raw contraction ratio keeps falling, as it does in Newton's quadratic region.
        Base best_raw = scope.diagnostics_on() ? maximum_abs(numerical(at, rows_before)) : 0.0;
        Base ratio = 1.0;
        bool trial_recovery = false;
        Base trial_scale = 1.0;
        core::ResidualHistory damped;  // raw row maxima over consecutive damped recovery corrections
        const bool detailed = std::is_same_v<Base, double> && scope.diagnostics_on();
        std::size_t final_step = 0;
        bool stopped = false;
        const auto stop = [&](const char* reason) {
            if (stopped) return;
            stopped = true;
            if (std::string_view(reason) == "Jacobian rank") core::retain_first(result.first_numerical_failure, [&] {
                core::Convergence refusal; refusal.stop = core::Stop::singular; if (result.last_rank) refusal.rank = std::make_shared<core::NumericalRank>(*result.last_rank);
                refusal.residual = static_cast<double>(maximum_abs(numerical(at, rows_before)));
                refusal.cause = reason; return refusal;
            });
            if constexpr (!std::is_same_v<Base, double>) { result.message = reason; core::retain_failure(result.first_numerical_failure, reason, "wide correction"); }
            if (detailed) scope.observe([&] { scope.probe("polish_stop",
                {double(maximum_abs(numerical(at, rows_before))),
                 double(maximum_abs(classified_rows(at, rows_before))),
                 static_cast<double>(final_step)}, static_cast<std::uint32_t>(final_step), reason); });
        };
        for (std::size_t step = 0; step < static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value); ++step) {
            final_step = step;
            if constexpr (!std::is_same_v<Base, double>) {
                if (limits && (limits->iterations_spent() || !limits->spend(limits->consumed, limits->iterations))) { stop("iteration budget"); return false; }
                ++result.iterations;
            }
            if (expired()) { stop("elapsed budget"); break; }
            const auto rows_scaled = numerical(at, rows_before);
            if constexpr (!std::is_same_v<Base, double>) if (!std::all_of(rows_scaled.begin(), rows_scaled.end(),
                [](Base value) { return wide_math::isfinite(value); })) { stop("nonfinite_or_unavailable_rows"); break; }
            const Base defect = maximum_abs(rows_scaled);
            const bool extension = step >= static_cast<std::size_t>(policy::EQUILIBRIUM_POLISH_STEPS.value);
            const bool final_correction = (newton || !std::is_same_v<Base, double>) && defect <= row_tolerance;
            if (extension && defect <= row_tolerance && !final_correction) { stop("raw criterion"); break; }
            const Base classified_before = maximum_abs(classified_rows(at, rows_before));
            if ((!final_correction && !(defect > 0.0)) || !wide_math::isfinite(defect)) { stop("nonfinite/zero defect"); break; }
            // Recovery stops when its damped contraction, carried over the corrections left under the cap, misses the criterion.
            if (trial_recovery) {
                damped.push(static_cast<double>(defect));
                const auto rate = damped.contraction();
                if (rate && static_cast<double>(defect) * std::pow(*rate, policy::EQUILIBRIUM_MAX_ITERATIONS.value - double(step))
                        > static_cast<double>(row_tolerance)) {
                    stop("recovery stalled"); break;
                }
            }
            const auto jacobian = numerical(at, compiled.jacobian_at(at), full.size());
            if (expired()) { stop("elapsed budget"); return false; }
            const auto n = static_cast<Eigen::Index>(full.size());
            if (jacobian.size() != full.size() * full.size() || rows_before.size() != full.size()) { stop("Jacobian shape"); break; }
            if (guarded_newton && !std::all_of(jacobian.begin(), jacobian.end(),
                [](Base value) { return wide_math::isfinite(value); })) { stop("nonfinite Jacobian"); break; }
            const Eigen::Map<const Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> j(jacobian.data(), n, n);
            Vector column = Vector::Ones(n);
            if constexpr (!std::is_same_v<Base, double>) column = physical_projection(compiled, at).scale;
            const auto a = [&]() -> decltype(auto) { if constexpr (std::is_same_v<Base, double>) return (j); else return Matrix(j * column.asDiagonal()); }();
            const auto lu = a.fullPivLu();
            result.last_rank = core::numerical_rank(lu, "numerical");
            core::ConditioningObservation observation{"Newton", "", "", "LU rank",
                result.last_rank->rank == result.last_rank->dimension ? "endpoint_factor_unavailable" : "rank_deficient"};
            observation.iteration = step; observation.precision_bits = std::numeric_limits<Base>::digits;
            observation.controls = {double(at.temperature), double(at.pressure)};
            core::observe_conditioning(*result.last_rank, std::move(observation));
            for (Eigen::Index i = 0; i < n; ++i) result.last_rank->observation.column_scale.push_back(double(column[i]));
            if (auto& warning = result.last_rank->first_conditioning_warning) warning->column_scale = result.last_rank->observation.column_scale;
            if (result.last_rank->first_conditioning_warning)
                core::retain_rank(result.rank_records, *result.last_rank, narrow_values(full.begin(), full.end()));
            else core::retain_rank(result.rank_records, *result.last_rank);
            if (scope.diagnostics_on()) scope.probe("numerical_rank", {double(result.last_rank->rank), double(result.last_rank->dimension), result.last_rank->threshold}, 0, "numerical");
            // An optimizer polish keeps a rank-deficient full-pivot step: a free trace ion split under the absolute charge row (#1071).
            if ((guarded_newton || !compiled.log_charge_rows.empty()) && lu.rank() != n) { stop("Jacobian rank"); break; }
            const auto retain_error = [&](Base error, const policy::Tolerance& criterion, const char* cause) {
                auto observation = result.last_rank->observation;
                observation.error = double(error); observation.criterion = criterion.name;
                observation.criterion_value = criterion.value; observation.cause = cause;
                core::observe_conditioning(*result.last_rank, std::move(observation));
                core::retain_rank(result.rank_records, *result.last_rank,
                    result.last_rank->first_conditioning_warning ? narrow_values(full.begin(), full.end()) : std::vector<double>{});
            };
            const Vector right = -Eigen::Map<const Vector>(rows_scaled.data(), n), d = lu.solve(right);
            const Vector delta = [&]() -> Vector { if constexpr (std::is_same_v<Base, double>) return d; else return column.asDiagonal() * d; }();
            if ((guarded_newton || refine_trial) && !delta.allFinite()) { stop("nonfinite direction"); break; }
            if constexpr (!std::is_same_v<Base, double>) {
                const Vector f = -right;
                const auto e = physical_error(compiled, at, j, f, lu, d, column);
                if (!e.cause.empty()) { stop(e.cause.c_str()); break; }
                const Base radius = physical_radius(compiled, at, e.radius);
                result.linear_solve_diagnostics = LinearSolveDiagnostics{
                    "canonical_defining_rows", "", double(defect), double(e.floor.maxCoeff()), double(e.inverse_norm),
                    double(e.condition), double(radius), std::numeric_limits<Base>::digits, step};
                auto& observation = result.last_rank->observation;
                observation.condition = double(e.condition); observation.available = std::isfinite(observation.condition);
                observation.norm = "infinity"; observation.method = "inverse_action";
                retain_error(radius, policy::MIXTURE_FORWARD_ERROR, "");
                const bool floor = e.cause.empty() && (right.cwiseAbs().array() <= e.row_floor.array()).all();
                if (defect <= row_tolerance && radius <= policy::MIXTURE_FORWARD_ERROR.value
                    && (delta.cwiseQuotient(column).cwiseAbs().maxCoeff() <= policy::EQUILIBRIUM_POLISH_CORRECTION.value
                        || (floor && e.correction.cwiseAbs().maxCoeff() <= e.floor.maxCoeff()))) {
                    best = full; corrected = true; break;
                }
            }
            Base trial_limit = 1.0;
            if (refine_trial) {
                Base log_step = 0.0;
                for (const auto& shape : at.shape) {
                    for (auto i : shape.x) log_step = std::max(log_step, wide_math::abs(delta(i)));
                    log_step = std::max(log_step, wide_math::abs(delta(shape.volume)));
                }
                if (log_step > policy::EQUILIBRIUM_TRIAL_LOG_STEP.value)
                    trial_limit = policy::EQUILIBRIUM_TRIAL_LOG_STEP.value / log_step;
            }
            const Vector scale = coordinate_scale(compiled, at, n);
            bool positive = true;
            for (std::size_t p = 0; p < at.shape.size(); ++p) if (at.shape[p].amount) positive = positive && at.amount[p] > 0.0;
            Eigen::Index controlling = 0;
            const Base scaled_correction = delta.allFinite() ? delta.cwiseQuotient(scale).cwiseAbs().maxCoeff(&controlling)
                : std::numeric_limits<double>::infinity();
            Base bound = std::numeric_limits<double>::infinity();
            const bool size_admitted = positive && full == best && classified_before <= row_tolerance
                && scaled_correction <= policy::EQUILIBRIUM_POLISH_CORRECTION.value;
            if constexpr (std::is_same_v<Base, double>) if (size_admitted && lu.rank() == n && !final_correction) {
                bound = core::polish_forward_bound(j, Eigen::Map<const Vector>(rows_scaled.data(), n), delta, scale);
                retain_error(bound, policy::EQUILIBRIUM_POLISH_CORRECTION, "inner_floor_unavailable");
            }
            const bool correction_stop = size_admitted && bound <= policy::EQUILIBRIUM_POLISH_CORRECTION.value;
            if constexpr (std::is_same_v<Base, double>) if (corrected || (extension && defect <= row_tolerance && !final_correction)) { stop("raw criterion"); break; }
            if (detailed) scope.observe([&] {
                Base feed_relative = 0.0;
                for (const auto& shape : at.shape) if (shape.amount)
                    feed_relative = std::max(feed_relative, wide_math::abs(delta(*shape.amount)) / compiled.feed_total);
                scope.probe("polish_correction", {double(scaled_correction), static_cast<double>(controlling), double(feed_relative), double(bound),
                double(classified_before), full == best ? 1.0 : 0.0, positive ? 1.0 : 0.0,
                size_admitted ? static_cast<double>(lu.rank()) : -1.0, correction_stop ? 1.0 : 0.0},
                static_cast<std::uint32_t>(step)); });
            double correction = std::numeric_limits<double>::quiet_NaN(), control = correction;
            if (detailed) scope.observe([&] {
                Eigen::Index index = 0;
                correction = static_cast<double>(delta.cwiseAbs().maxCoeff(&index)); control = static_cast<double>(index);
                scope.probe("polish_direction", {correction, control, double(defect), double(classified_before),
                    extension ? 1.0 : 0.0, double(ratio)}, step, delta.allFinite() ? "" : "nonfinite direction");
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
                    } else quantity = index == static_cast<Eigen::Index>(eta) ? (model && model->parameters().ideal_gas ? "log inverse unit-amount volume" : "log packing fraction") : "phase amount in mol";
                }
                auto intensive = compiled.primal_size - compiled.free_intensive_count;
                if (compiled.problem.temperature.free && index == static_cast<Eigen::Index>(intensive++)) quantity = "T in K";
                if (compiled.problem.pressure.free && index == static_cast<Eigen::Index>(intensive)) quantity = "P in Pa";
                scope.probe("polish_coordinate", {control, phase, component, double(full[index])}, step, quantity);
            });
            // The endpoint estimate reuses this step's factor: the returned point is `full == best` with the row criterion met.
            const auto publish_endpoint = [&] {
                if constexpr (std::is_same_v<Base, double>) if (!compiled.problem.solver_options.nested_solve)
                    result.last_rank = core::conditioning(lu, "numerical", step, "endpoint", {double(at.temperature), double(at.pressure)});
            };
            if (correction_stop && !final_correction) {
                if constexpr (std::is_same_v<Base, double>) {
                    publish_endpoint();
                    retain_error(bound, policy::EQUILIBRIUM_POLISH_CORRECTION, "inner_floor_unavailable");
                }
                stop("scaled correction"); break;
            }

            bool accepted = false, representable = false, trial_in_box = false;
            // Bounded damping retains raw contraction; preparation admits only positive finite x.
            // Best-point selection below protects the returned classified residual.
            std::uint32_t trial = 0;
            // Newton attempts and wide correction also try the fourth fraction.
            for (const Base damping : {1.0, 0.5, 0.25, 0.125}) {
                if constexpr (std::is_same_v<Base, double>) if (trial == 3 && !newton) break;
                const Base alpha = std::min((trial_recovery ? trial_scale : Base(1)) * damping, trial_limit);
                if (final_correction && alpha != 1.0) break;
                const auto trial_index = trial++;
                if (expired()) { stop("elapsed budget"); break; }
                std::vector<Base> candidate = full;
                for (Eigen::Index i = 0; i < n; ++i)
                    candidate[static_cast<std::size_t>(i)] += alpha * delta(i);
                representable = representable || candidate != full;
                if constexpr (!std::is_same_v<Base, double>) if (candidate == full) continue;
                if (trial_recovery && candidate == full) continue;
                if (guarded_newton && !inside_box(candidate)) {
                    if constexpr (std::is_same_v<Base, double>) if (detailed) scope.probe("polish_trial_bounds", candidate, static_cast<std::uint32_t>(step));
                    continue;
                }
                trial_in_box = true;
                PhasePointValues<Base> at_candidate;
                std::vector<Base> rows_after;
                if (limits && !limits->spend(limits->evaluations, limits->evaluation_cap)) return false;
                try {
                    at_candidate = compiled.prepare_values(candidate, true, &result.first_numerical_failure);
                    ++polish_evaluations;
                    if (expired()) { stop("elapsed budget"); return false; }
                    rows_after = compiled.evaluate_at(at_candidate).values;
                } catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) {
                    core::retain_failure(result.first_numerical_failure, error, "polish trial");
                    if (detailed) scope.probe("polish_trial", {double(correction), double(control), double(defect), double(classified_before),
                        NAN, NAN, double(alpha), 0.0, extension ? 1.0 : 0.0, double(ratio),
                        static_cast<double>(trial_index)}, step, error.what());
                    continue;
                }
                // Judged with the step's moiety sizes: a size crossing the trace threshold rescales a row by 1/s_k (#1094).
                const Base after = maximum_abs(numerical(at_candidate, rows_after, 1, &at));
                const Base classified_after = maximum_abs(classified_rows(at_candidate, rows_after));
                const bool admit = final_correction || (free_amount_newton && after <= row_tolerance)
                    || (after < defect && after <= row_tolerance)
                    || after < (extension && !trial_recovery ? ratio : 1.0) * defect;
                if (detailed) scope.probe("polish_trial", {double(correction), double(control), double(defect), double(classified_before),
                    double(after), double(classified_after), double(alpha), admit ? 1.0 : 0.0, extension ? 1.0 : 0.0, double(ratio),
                    static_cast<double>(trial_index)}, step,
                    final_correction ? "final correction" : admit ? "raw_contraction" : extension && !trial_recovery ? "extension ratio" : "raw noncontraction");
                else scope.probe("polish_step", {double(defect), double(after), double(classified_before), double(classified_after), double(alpha)},
                    static_cast<std::uint32_t>(step));
                if (!admit) continue;
                if constexpr (std::is_same_v<Base, double>) if (limits && limits->admit) if (const auto refused = limits->admit(at_candidate)) {
                    core::retain_failure(result.first_numerical_failure, *refused, "guarded polish");
                    continue;
                }
                if constexpr (std::is_same_v<Base, double>) if (preserve_trial_branch && !trial_branch(at_candidate, candidate)) {
                    if (expired()) return false;
                    continue;
                }
                if (trial_recovery) trial_scale = std::min(Base(1), 2.0 * alpha);
                if (alpha == trial_limit) damped = {};
                if (!final_correction) ratio = after / defect;
                full = std::move(candidate);
                at = std::move(at_candidate);
                rows_before = std::move(rows_after);
                if (final_correction || classified_after <= best_classified) {
                    best = full;
                    best_classified = classified_after;
                    if (detailed) best_raw = after;
                }
                accepted = true;
                break;
            }
            if (!accepted) {
                if (expired()) { stop("elapsed budget"); break; }
                if ((refine_trial || !std::is_same_v<Base, double>) && !representable) { stop("no representable progress"); break; }
                if (!trial_in_box) { stop("trial bounds"); break; }
                if (trial_recovery) trial_scale *= 0.125;
                else if (refine_trial && best_classified > row_tolerance) {
                    trial_recovery = true;
                    damped.push(static_cast<double>(defect));
                    scope.probe("trial_recovery_begin", {double(defect), double(classified_before)}, step);
                } else {
                    // Another linear solver can leave the correction just above its criterion; the endpoint conditioning is still the factor's.
                    if (full == best && classified_before <= row_tolerance && lu.rank() == n) {
                        publish_endpoint();
                        if constexpr (std::is_same_v<Base, double>) core::retain_rank(result.rank_records, *result.last_rank);
                    }
                    stop("no accepted trial"); break;
                }
            }
            if (step + 1 == static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)) { capped = true; stop("iteration cap"); }
            if constexpr (std::is_same_v<Base, double>) if (final_correction) { corrected = true; stop("final correction"); break; }
        }
        if (detailed && best != full) scope.probe("polish_stop", {double(best_raw), double(best_classified),
            static_cast<double>(final_step)}, final_step, "best classified point returned");
        full = std::move(best);
        if constexpr (!std::is_same_v<Base, double>) break;
        else if (polish_round++ || solid_violations(true).empty()) break;  // one flip round, then one more polish
    }
    return true;
}

EquilibriumResult spinodal_search_failure(const CompiledSystem& compiled, EquilibriumResult result, const std::string& cause) {
    refuse(result);
    result.message = cause; result.convergence.stop = core::Stop::left_domain;
    result.spinodal.emplace();
    if (compiled.spinodal_search) result.spinodal->search = *compiled.spinodal_search;
    core::retain_failure(result.first_numerical_failure, cause, "spinodal search");
    return result;
}

static EquilibriumResult association_result(EquilibriumResult result, const core::AssociationObservations& observations) {
    observations.merge(result.rank_records);
    return result;
}

EquilibriumResult refine_selected_endpoint(const CompiledSystem& compiled, EquilibriumResult result,
    core::Scope& scope, const DriverLimits* limits) {
    core::AssociationObservations observations;
    if (!result.has_state() || !mixture_precision_eligible(compiled) || !result.prepared_endpoint) return association_result(std::move(result), observations);
    DriverLimits defaults{static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value),
        std::chrono::steady_clock::time_point::max(), result.iterations};
    if (!limits) limits = &defaults;
    const auto evaluated_before = limits->evaluations;
    core::Scope span("mixture_precision");
    result.numerical_estimates = equilibrium_numerical_estimates(compiled, result, result.prepared_endpoint.get());
    const auto& error = *result.numerical_estimates;
    const double radius = error.physical_radius;
    scope.probe("mixture_precision_decision", {*result.temperature, radius}, 0, error.cause.c_str());
    if (!error.cause.empty() || !std::isfinite(radius) || radius <= policy::MIXTURE_FORWARD_ERROR.value) return association_result(std::move(result), observations);
    result.linear_solve_diagnostics = LinearSolveDiagnostics{"canonical_defining_rows"};
    auto& diagnostic = *result.linear_solve_diagnostics;
    diagnostic.forward_error_estimate = radius;
    if constexpr (std::numeric_limits<WiderScalar>::digits <= std::numeric_limits<double>::digits) {
        diagnostic.reason = "extra precision unavailable"; return association_result(std::move(result), observations);
    }
    scope.probe("mixture_precision_activate", {*result.temperature, *result.pressure, radius});
    auto refined = result;
    refined.message.clear();
    bool accepted = false, corrected = false, capped = false;
    std::size_t evaluations = 0;
    std::vector<WiderScalar> point(result.coordinates.begin(), result.coordinates.end());
    std::vector<double> lower(point.size(), -INFINITY), upper(point.size(), INFINITY);
    try {
        GibbsMinimizer bounds(compiled, result.coordinates, nullptr, limits, false);
        const auto box = [&](auto& y) { return inside_correction_box(compiled, y, lower, upper, true); };
        accepted = bounds.physical_bounds(lower, upper) && box(point) && polish_correction(compiled, point, refined, scope, limits,
            true, false, true, box, [](const auto&, const auto&) { return true; }, [](bool) { return std::string{}; },
            evaluations, corrected, capped, refined.convergence.stop) && corrected && !capped;
        if (accepted) {
            std::vector<double> exported;
            for (const auto& value : point) exported.push_back(static_cast<double>(value));
            if (!limits->spend(limits->evaluations, limits->evaluation_cap) || limits->expired())
                limits->check("work budget exhausted before export admission");
            const auto at = compiled.prepare(exported, false, &refined.first_numerical_failure);
            const auto rows = compiled.evaluate_at(at);
            refined.residuals = rows.values; refined.classified_residuals = compiled.classified_residuals(at, rows.values);
            refined.convergence.residual = maximum_abs(refined.classified_residuals);
            refined.requested_tolerance_met = refined.convergence.residual <= compiled.problem.solver_options.value(policy::row);
            accepted = *refined.requested_tolerance_met && maximum_abs(compiled.numerical_rows(at, rows.values))
                <= compiled.problem.solver_options.value(policy::row) && box(exported);
            clear_state(refined);
            refined.candidate = accepted ? Candidate::Admitted : Candidate::Refused;
            refined = publish_candidate(compiled, at, std::move(refined), scope, limits);
            accepted = refined.success;
            if (accepted) {
                refined.numerical_estimates = equilibrium_numerical_estimates(compiled, refined, point, limits);
                const auto& returned = *refined.numerical_estimates;
                const auto admission = check_density_branches(compiled, refined, {}, true, scope, limits, &at);
                accepted = admission.accepted && returned.cause.empty() && returned.physical_radius <= policy::MIXTURE_FORWARD_ERROR.value;
                if (!accepted) refined.message = admission.accepted ? "returned-point error target not met" : admission.message;
            }
        } else if (refined.message.empty()) refined.message = limits && limits->expired() ? "work budget exhausted" : "wide correction refused";
    } catch (const std::bad_alloc&) { throw; } catch (const std::exception& failure) {
        accepted = false; refined.message = failure.what();
        core::retain_failure(refined.first_numerical_failure, failure, "mixture refinement");
    }
    refined.evaluations += limits->evaluations - evaluated_before;
    if (!accepted) {
        core::retain_failure(refined.first_numerical_failure, refined.message, "mixture refinement");
        StageWork::publish_observations(refined, result);
        result.iterations = refined.iterations; result.evaluations = refined.evaluations;
        diagnostic.reason = "refinement refused: " + refined.message; return association_result(std::move(result), observations);
    }
    refined.message = result.message;
    return association_result(std::move(refined), observations);
}

Retry retry_after(Attempt kind, bool reaction, bool non_error, bool fatal, Ipopt::SolverReturn status,
    TopologyEvent topology, core::Stop stop, Candidate candidate) {
    if (kind == Attempt::Newton || kind == Attempt::PhysicalAmounts) return Retry::None;
    if (reaction) return non_error && (topology == TopologyEvent::PhaseDisappeared || stop == core::Stop::budget_exhausted)
        ? Retry::PhysicalAmounts : Retry::None;
    return !fatal && topology != TopologyEvent::PhaseDisappeared && (status == Ipopt::LOCAL_INFEASIBILITY || status == Ipopt::MAXITER_EXCEEDED
        || ((status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT)
            && (topology == TopologyEvent::PhasesCoalesced || candidate == Candidate::Unstable))) ? Retry::Source : Retry::None;
}

// The one Ipopt lifecycle shared by the public seed route and the private coordinate route.
static EquilibriumResult run_driver_once(
    const CompiledSystem& declared,
    std::vector<double> initial, EquilibriumResult result, core::Scope& scope,
    std::vector<double>* attempt, const DriverLimits* limits, Attempt kind,
    std::shared_ptr<const OptimizerState> warm, bool predicted, bool refine_trial
) {
    core::AssociationObservations observations;
    result.attempt = kind == Attempt::Warm ? Attempt::Cold : kind;
    result.retry = Retry::None;
    // The work a failure returns: the driver's counts and observations, never a candidate.
    auto record = result;
    std::size_t polish_evaluations = 0;
    const auto failed = [&](std::string message) {
        if (polish_evaluations) record.convergence.evaluations = record.evaluations += polish_evaluations;
        record.message = std::move(message);
        return record;
    };
    bool publishing = false;
    try {
    const bool newton = kind == Attempt::Newton, physical_amounts = kind == Attempt::PhysicalAmounts;
    const double row_tolerance = declared.problem.solver_options.value(policy::row);
    // The optimizer works one policy decade below the caller criterion, which alone classifies.
    // A feasibility cell has a constant objective and is Newton-polished below, so it stops there.
    const double inner_tolerance = row_tolerance
        * (declared.feasibility_cell() ? 1.0 : policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value);
    if (limits && limits->attempt_cap && limits->attempts >= limits->attempt_cap) limits->exhausted = true;
    const auto expired = [&] { return limits && limits->expired(); };
    if (expired()) { result.message = "elapsed work budget exhausted before driver";
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        reject_elapsed(result, scope); return association_result(std::move(result), observations); }
    if (limits) ++limits->attempts;
    if (declared.problem.spinodal) {
        try { spinodal_search(declared, limits); }
        catch (const std::bad_alloc&) { throw; }
        catch (const std::exception& error) {
            result = spinodal_search_failure(declared, std::move(result), error.what());
            scope.finish(result.convergence); return result;
        }
    }
    auto fallback = limits || newton ? initial : std::vector<double>{};
    Ipopt::SmartPtr<GibbsMinimizer> nlp = new GibbsMinimizer(declared, std::move(initial), &scope, limits, physical_amounts);
    record.convergence.residual = std::numeric_limits<double>::quiet_NaN();
    const std::function<bool(double)> warm_start = warm ? std::function<bool(double)>([&](double cap) {
        const bool enabled = nlp->initialize_warm(*warm, cap, predicted);
        if (enabled) result.attempt = Attempt::Warm;
        return enabled;
    }) : std::function<bool(double)>{};
    const auto [initialized, returned] = [&] {
        if (newton) return std::make_pair(Ipopt::Solve_Succeeded, Ipopt::Solve_Succeeded);
        if (const char* file = std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE")) result.diagnostic_options_file = file;
        try { nlp->begin_ipopt(); return run_ipopt(nlp, inner_tolerance,
            static_cast<int>(limits ? std::min(limits->iterations - limits->consumed,
                static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)) : policy::EQUILIBRIUM_MAX_ITERATIONS.value),
            limits ? std::max(0.0, std::chrono::duration<double>(limits->deadline - std::chrono::steady_clock::now()).count()) : 0.0,
            true, 1.0 / declared.feed_total, &result.solver_settings, row_tolerance,
            warm_start, &result.linear_solver);
        } catch (...) { nlp->publish_observations(result); capture_driver(*nlp, result, record, observations); throw; }
    }();
    if (!newton) {
        nlp->finish_ipopt(int(initialized), int(returned), predicted, physical_amounts);
    }
    nlp->publish_observations(result);
    capture_driver(*nlp, result, record, observations);
    if (nlp->step_rule_error()) throw std::logic_error(*nlp->step_rule_error());
    const Ipopt::SolverReturn status = newton ? Ipopt::SUCCESS : nlp->status();
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
    record.solver_status = result.solver_status;
    const auto limited_result = [&] {
        observations.merge(result.rank_records); StageWork::publish_observations(result, record);
        auto limited = failed("elapsed work budget exhausted");
        reject_elapsed(limited, scope); return limited;
    };
    if (expired()) return limited_result();
    std::vector<double> point = nlp->solution();
    if (point.empty()) point = limits || newton ? std::move(fallback) : declared.ideal ? declared.generic_point() : boundary_start(declared);
    const auto& cached = nlp->last_prepared();
    const auto& rows = nlp->last_rows();
    const bool reuse = !nlp->solution().empty() && point.size() == declared.primal_size && cached && rows
        && cached->point.size() == point.size()
        && std::memcmp(cached->point.data(), point.data(), point.size() * sizeof(double)) == 0;
    std::vector<double> full = !declared.feasibility_cell()
        ? add_multipliers(declared, point, &result, reuse ? &*cached : nullptr, reuse ? &*rows : nullptr) : point;
    // ALG-EQUILIBRIUM: each free solid is classified once at the Ipopt endpoint by complementarity: present (saturation
    // row) when n_s/N_feed exceeds its bound multiplier z_s = -A_s/RT, absent (n_s = 0) otherwise. The classified
    // declaration is square for the polish; after it a present solid keeps n_s > 0 and an absent one A_s/RT <= criterion.
    std::optional<CompiledSystem> classified;
    if (!declared.solids.empty()) {
        classified = declared;
        const auto& bound = nlp->lower_multipliers();
        for (std::size_t s = 0; s < declared.solids.size(); ++s) {
            auto& solid = classified->solids[s];
            solid.present = solid.column
                && full[*solid.column] / declared.feed_total > (bound.empty() ? 0.0 : bound[*solid.column]);
            if (solid.column && !solid.present) full[*solid.column] = 0.0;
        }
    }
    const CompiledSystem& compiled = classified ? *classified : declared;
    // Names the solids violating their classification and, with `flip`, flips them. An absent amount is exactly zero.
    const auto solid_violations = [&](bool flip) {
        std::string names;
        if (!classified) return names;
        const auto at = compiled.prepare(full, true, &result.first_numerical_failure);
        for (std::size_t s = 0; s < classified->solids.size(); ++s) {
            auto& solid = classified->solids[s];
            if (!solid.column) continue;
            if (!solid.present) full[*solid.column] = 0.0;
            if (solid.present ? full[*solid.column] > 0.0 : compiled.solid_affinity(at, s) <= row_tolerance) continue;
            names += (names.empty() ? "" : ", ") + declared.problem.solids[s].name;
            if (!flip) continue;
            solid.present = !solid.present;
            if (!solid.present) full[*solid.column] = 0.0;
        }
        return names;
    };
    const bool free_amount_newton = newton && !compiled.feasibility_cell();
    std::vector<double> lower(newton ? full.size() : 0, -INFINITY), upper(lower.size(), INFINITY);
    const auto inside_box = [&](auto& candidate) { return inside_correction_box(compiled, candidate, lower, upper, free_amount_newton); };
    if (free_amount_newton && (!nlp->physical_bounds(lower, upper) || !inside_box(full)
        || check_topology(compiled, full).event != TopologyEvent::None)) {
        refuse(result); result.message = "Newton seed bounds/topology"; return association_result(std::move(result), observations);
    }
    // The optimizer's amount-weighted rows underresolve intensive trace defects.
    // Exact Newton steps contract raw rows at any positive supported composition.
    // Curvature may raise classified rows temporarily; return the best classified
    // point visited, including the optimizer start, rather than the last step.
    bool capped = false, in_bounds = true, corrected = false;

    core::Stop initial_branch_stop = core::Stop::converged;
    const auto trial_branch = [&](const PhasePoint& prepared, const std::vector<double>& coordinates) {
        result.convergence.stop = core::Stop::left_domain;
        if (scope.diagnostics_on()) scope.probe("trial_branch_coordinates", coordinates);
        try {
            EquilibriumResult view; view.coordinates = coordinates;
            for (std::size_t p = 0; p < prepared.shape.size(); ++p) {
                const double total = sum(prepared.composition[p]);
                if (!std::isfinite(total)) result.convergence.stop = core::Stop::nonfinite;
                if (!(total > 0.0) || !std::isfinite(total)) throw std::domain_error("invalid trial normalization");
                for (const auto i : prepared.shape[p].x) view.coordinates[i] -= std::log(total);
            }
            const auto branch = check_density_branches(compiled, view, {}, true, scope, limits);
            if (branch.accepted && !expired()) return true;
            core::retain_failure(result.first_numerical_failure, branch.message, "trial branch");
        } catch (const core::ConvergenceFailure& error) {
            result.convergence = error.record; core::retain_failure(result.first_numerical_failure, error, "trial branch");
        } catch (const std::exception& error) { core::retain_failure(result.first_numerical_failure, error, "trial branch"); }
        if (expired()) return false;
        scope.probe("trial_branch_cause", {static_cast<double>(result.convergence.stop)}, 0, result.first_numerical_failure ? result.first_numerical_failure->cause.c_str() : "");
        return false;
    };
    // A budget stop polishes too: an unattainable request can exhaust Ipopt after passing a
    // better point, and the polish never returns a worse classified point than its start.
    const auto polish = [&](auto& coordinates) {
        return polish_correction(compiled, coordinates, result, scope, limits, newton, refine_trial,
            non_error, inside_box, trial_branch, solid_violations, polish_evaluations, corrected, capped, initial_branch_stop);
    };
    const auto end_point = refine_trial && !newton && status == Ipopt::SUCCESS ? full : std::vector<double>{};
    if (!polish(full)) return limited_result();
    if (expired()) return limited_result();
    if (newton) {
        if (!free_amount_newton) {
            std::vector<double> row_lower(compiled.rows.size()), row_upper(compiled.rows.size());
            in_bounds = nlp->get_bounds_info(static_cast<Ipopt::Index>(full.size()), lower.data(), upper.data(),
                static_cast<Ipopt::Index>(compiled.rows.size()), row_lower.data(), row_upper.data());
        }
        in_bounds = inside_box(full) && in_bounds;
    }
    PhasePoint prepared;
    try {
        prepared = compiled.prepare(full, false, &result.first_numerical_failure);
    } catch (const ReferenceUnavailable& error) {

        result.message = error.what() + ipopt_failure;
        result.convergence.stop = core::Stop::left_domain;
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        result.iterations = result.convergence.iterations = nlp->iterations();
        result.evaluations = result.convergence.evaluations = nlp->evaluations() + polish_evaluations;
        refuse(result);
        result.retry = retry_after(result.attempt, declared.has_reaction, non_error, fatal, status, result.topology_event, result.convergence.stop, result.candidate);
        scope.finish(result.convergence);
        return association_result(std::move(result), observations);
    }
    if (expired()) return limited_result();
    auto evaluation = compiled.evaluate_at(prepared);
    // A negative distance at the optimal end point is a witness; a polish that leaves it for a distinct point
    // of no lower distance (a stationary point of another basin) returns the end point instead (#49).
    if (!end_point.empty()) {
        const auto e = check_topology(compiled, end_point), p = check_topology(compiled, full);
        const double apart = std::max(largest_difference(e.log_packing, p.log_packing),
            largest_difference(e.compositions[0], p.compositions[0]));  // a trial has one phase
        if (apart > policy::PHASE_COINCIDENCE_DISTANCE.value) {
            auto at = compiled.prepare(end_point, false, &result.first_numerical_failure);
            auto kept = compiled.evaluate_at(at);
            ++polish_evaluations;
            if (kept.objective < 0.0 && !(evaluation.objective < kept.objective)) {
                scope.probe("trial_end_point_kept", {kept.objective, evaluation.objective, apart});
                std::tie(full, prepared, evaluation) = std::tuple(end_point, std::move(at), std::move(kept));
            }
        }
    }
    if (attempt && !fatal && !nlp->solution().empty()) *attempt = full;
    result.residuals = evaluation.values;
    result.classified_residuals = compiled.classified_residuals(prepared, evaluation.values);
    result.objective = evaluation.objective;
    if (expired()) return limited_result();
    scope.probe("hessian_calls", static_cast<double>(nlp->hessian_calls()));
    scope.probe("prepared_points", static_cast<double>(nlp->prepared_points()));
    core::Convergence& convergence = result.convergence;
    convergence = {};
    convergence.residual = maximum_abs(result.classified_residuals);
    core::stopping_criterion(convergence, compiled.problem.solver_options, policy::row);
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
    if (initial_branch_stop != core::Stop::converged) convergence.stop = initial_branch_stop;
    if (std::isfinite(convergence.residual)) {
        result.requested_tolerance_met = convergence.residual <= row_tolerance;
    }
    result.message = core::describe(convergence);
    scope.finish(convergence);
    const bool admitted = non_error && (!newton || (!capped && in_bounds && corrected))
        && (convergence.stop == core::Stop::converged
            || (!newton && convergence.stop == core::Stop::converged_to_floor && status == Ipopt::SUCCESS));
    result.candidate = admitted ? Candidate::Admitted : fatal ? Candidate::Absent : Candidate::Refused;
    if (newton && admitted && !free_amount_newton) result.solver_status = "NEWTON_CONVERGED";

    if (!newton) scope.probe("optimizer_constraint_multipliers", nlp->constraint_multipliers());
    if (compiled.feasibility_cell()) result.constraint_multipliers = newton
        ? std::vector<double>(compiled.rows.size(), 0.0) : nlp->constraint_multipliers();
    observations.merge(result.rank_records); StageWork::publish_observations(result, record);
    publishing = true;
    result = publish_candidate(compiled, prepared, std::move(result), scope, limits);
    if (!result.success) {
        result.message += ipopt_failure;
        if (result.first_numerical_failure) result.message += "; trial preparation failed: " + core::describe(*result.first_numerical_failure);
        result.retry = retry_after(result.attempt, declared.has_reaction, non_error, fatal, status, result.topology_event, result.convergence.stop, result.candidate);
    } else if (!newton && result.requested_tolerance_met.value_or(false)) result.optimizer_state = nlp->optimizer_state();
    return association_result(std::move(result), observations);
    } catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) {
        auto& observed = publishing ? record : result;
        observations.merge(observed.rank_records);
        core::retain_failure(observed.first_numerical_failure, error, publishing ? "equilibrium publication" : "equilibrium preparation");
        if (!publishing) StageWork::publish_observations(result, record);
        auto failure = failed(error.what());
        failure.failure = std::current_exception();
        return failure;
    }
}

EquilibriumResult recover_physical_amounts(
    const CompiledSystem& compiled, EquilibriumResult result, bool eligible, core::Scope& scope,
    const std::vector<double>* restart, const DriverLimits* limits, std::vector<double>* attempt,
    bool admit_density, std::size_t* declarations, std::size_t declaration_budget,
    bool* invoked, Work counts
) {
    const auto exhausted = [&] {
        return limits && (limits->expired() || limits->iterations_spent()
            || (limits->recovery && limits->caps_spent()));
    };
    if (result.success || !eligible || !compiled.has_reaction || compiled.feasibility_cell()
        || !compiled.solids.empty() || exhausted()
        || (declarations && *declarations >= declaration_budget)) return result;
    if (declarations) ++*declarations;
    if (invoked) *invoked = false;
    core::Scope retry_scope("reacting_physical_amount_retry");
    StageWork work;
    if (counts == Work::Accumulated) work.add(result); else work.observe(result);
    auto retried = structural_result(compiled, scope);
    retried.convergence.residual = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> raw;
    const char* stage = "initialization";
    std::optional<std::string> cause;  // the failed stage's cause
    try {
        auto initial = restart ? *restart : boundary_start(compiled, true, limits, &work);
        if (exhausted()) cause = "work budget exhausted before physical-mol retry";
        else {
            stage = "driver";
            if (invoked) *invoked = true;
            auto driven = run_driver_once(compiled, std::move(initial), retried, retry_scope,
                &raw, limits, Attempt::PhysicalAmounts, {}, false, false);
            if (driven.failure) { cause = std::move(driven.message); driven.message = retried.message; }
            retried = std::move(driven);
        }
        if (!cause && retried.success && admit_density) {
            stage = "density branch admission";
            auto branches = check_density_branches(compiled, retried, {}, true, retry_scope, limits);
            if (!branches.accepted) cause = std::move(branches.message);
        }
        if (!cause && limits && limits->expired()) cause = "elapsed work budget exhausted";
    } catch (const std::exception& error) {
        core::retain_failure(retried.first_numerical_failure, error, "physical amount retry");
        cause = error.what();
    }
    if (cause) {
        core::retain_failure(retried.first_numerical_failure, *cause, "physical amount retry");
        retried.message = std::string(stage) + " failed: " + *cause + "; driver: " + retried.message;
        refuse(retried);
    }
    if (attempt && !raw.empty()) *attempt = std::move(raw);
    else if (attempt) {
        retried.topology_event = result.topology_event;
        retried.topology_phase = result.topology_phase;
        retried.topology_pair = result.topology_pair;
    }
    retried.message += "; amount chart physical_mol; retry after start_amount: " + result.message;
    std::move(work).report(retried);
    return retried;
}

EquilibriumResult run_driver(
    const CompiledSystem& compiled,
    std::vector<double> initial, EquilibriumResult result, core::Scope& scope,
    std::vector<double>* attempt, const DriverLimits* limits, Attempt kind,
    std::shared_ptr<const OptimizerState> warm, bool predicted, bool refine_trial, Retry runs
) {
    const auto restart = initial;
    DriverLimits shared = limits ? *limits : DriverLimits{static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value),
        std::chrono::steady_clock::time_point::max()};
    const auto* work = limits ? limits : warm ? &shared : nullptr;
    const auto primary = [&]() -> EquilibriumResult {
    if (!warm || kind == Attempt::Newton) return run_driver_once(compiled, std::move(initial), std::move(result), scope, attempt, limits, kind, {}, predicted, refine_trial);
    auto tried = [&] {
        core::Scope warm_scope("equilibrium_warm_attempt");
        return run_driver_once(compiled, initial, result, warm_scope, attempt, work, Attempt::Warm, warm, predicted, refine_trial);
    }();
    if (tried.failure) scope.probe("warm_failed", {}, 0, tried.message.c_str());
    else if (tried.attempt != Attempt::Warm || tried.success) return tried;
    scope.probe("warm_fallback", {static_cast<double>(tried.iterations), static_cast<double>(tried.evaluations)}, 0, tried.message.c_str());
    core::Scope cold_scope("equilibrium_cold_fallback");
    StageWork stages; stages.add(tried);
    auto cold = run_driver_once(compiled, std::move(initial), result, cold_scope, attempt, work, Attempt::Cold, {}, predicted, refine_trial);
    std::move(stages).report(cold);
    return cold;
    };
    auto selected = primary();
    const bool retry = selected.retry == Retry::PhysicalAmounts;
    if (runs == Retry::PhysicalAmounts && kind != Attempt::Newton) selected = recover_physical_amounts(compiled, std::move(selected), retry,
        scope, &restart, work, attempt);
    return selected;
}


// Direct admission of a critical-equation candidate (EqID critical_equations) in place of the density-root
// rule, which refuses P_rho = 0: a regular augmented Jacobian and one lowest zero eigenvalue of M, read from the
// mode rows' u columns, within |lambda| <= ||M u||/||u|| + C eps ||M||, above no resolved negative one. No
// stability certificate; the record is retained either way.
BranchCheck admit_null_mode(const CompiledSystem& compiled, EquilibriumResult& result) {
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto n = static_cast<Eigen::Index>(result.coordinates.size()), c = static_cast<Eigen::Index>(compiled.mode_count);
    const auto jacobian = compiled.jacobian(result.coordinates);
    const Eigen::Map<const Dense> J(jacobian.data(), n, n);
    const auto prepared = compiled.prepare(result.coordinates);
    const auto spectrum = helmholtz_spectrum(prepared.derivative[0], prepared.composition[0], compiled.problem.phases[0].support, compiled.feed_composition);
    const auto& M = spectrum.matrix;
    const Eigen::Map<const Eigen::VectorXd> u(result.coordinates.data() + compiled.mode_offset(), c);
    const auto lu = J.fullPivLu();
    const Eigen::VectorXd right = Eigen::VectorXd::Ones(n), solution = lu.solve(right);
    const auto singular = J.jacobiSvd().singularValues();
    HelmholtzModeRecord record;
    record.eigenvalues.assign(spectrum.values.data(), spectrum.values.data() + c);
    record.mode.assign(u.data(), u.data() + c);
    record.zero_band = (M * u).norm() / u.norm() + static_cast<double>(c) * std::numeric_limits<double>::epsilon() * M.norm();
    record.rank = static_cast<std::size_t>(lu.rank());
    result.last_rank = core::numerical_rank(lu, compiled.problem.critical ? "critical" : "diffusive_spinodal");
    record.condition_estimate = singular(0) / singular(n - 1);
    record.backward_error = core::backward_error((J * solution - right).cwiseAbs().maxCoeff(),
        J.cwiseAbs().rowwise().sum().maxCoeff(), solution.cwiseAbs().maxCoeff(), 1.0);
    if (compiled.problem.critical) {
        result.critical.emplace();
        static_cast<HelmholtzModeRecord&>(*result.critical) = record;
        result.critical->cubic = result.residuals.back();
    } else {
        const auto& selection = spinodal_search(compiled, nullptr, result.molar_densities[0], result.spinodal ? &result.spinodal->search : nullptr);
        result.spinodal.emplace();
        static_cast<HelmholtzModeRecord&>(*result.spinodal) = record;
        result.spinodal->search = selection;
        result.spinodal->normalized_null_residual = (M * u).norm() / M.norm();
    }
    BranchCheck outcome;
    if (record.rank != static_cast<std::size_t>(n))
        outcome.message = "augmented null-mode Jacobian is singular (rank " + std::to_string(record.rank) + ")";
    else if (compiled.problem.spinodal && (!((M * u).norm() / M.norm() <= policy::EQUILIBRIUM_ROW_RESIDUAL.value)
        || std::abs(u.squaredNorm() - 1) > policy::EQUILIBRIUM_ROW_RESIDUAL.value
        || result.molar_densities[0] < result.spinodal->search.lower || result.molar_densities[0] > result.spinodal->search.upper))
        outcome.message = "spinodal null mode or density bracket not satisfied";
    else if (!(std::abs(record.eigenvalues[0]) <= record.zero_band))
        outcome.message = "lowest eigenvalue of M is not a resolved zero";
    else if (c > 1 && !(record.eigenvalues[1] > record.zero_band))
        outcome.message = "second eigenvalue of M is not resolved positive";
    else {
        outcome.accepted = true;
        result.message += compiled.problem.critical ? "; critical-equation candidate, no stability certificate" : "; diffusive spinodal, no global stability certificate";
    }
    return outcome;
}

namespace {

// A cold reactive boundary at fixed T and free P: the phases present (nonzero amount) alone, in reaction
// equilibrium at the pressure seed, are the solution's present phases up to their pressure dependence. Their
// compositions replace the caller's guesses before the ordinary boundary start, and an unguessed zero-amount
// phase starts from the first present phase on its support (#152); none if any of it fails. It spends `limits`.
std::optional<std::vector<double>> speciated_start(const CompiledSystem& compiled, const DriverLimits* limits, StageWork& stages, bool correct_dense = false) {
    const auto& problem = compiled.problem;
    EquilibriumProblem present = problem;
    present.phases.erase(std::remove_if(present.phases.begin(), present.phases.end(),
        [](const PhaseSpec& phase) { return phase.is_zero_pinned(); }), present.phases.end());
    if (!compiled.has_reaction || !problem.pressure.free || problem.temperature.free || present.phases.empty()
        || present.phases.size() == problem.phases.size()) return std::nullopt;
    present.pressure = FreeValue::fixed(compiled.seed_pressure());
    try {
        const auto solved = try_solve_compiled(compile_problem(*compiled.model, present), std::nullopt, nullptr, limits);
        if (solved.failure) {
            stages.observe(solved);
            with_failure(solved, [&](const std::exception& error) { core::retain_failure(stages.first_numerical_failure, error, "speciated start"); });
            return std::nullopt;
        }
        stages.add(solved);
        if (!solved.success) return std::nullopt;
        const std::size_t components = compiled.model->parameters().component_count();
        CompiledSystem seeded = compiled;
        std::size_t index = 0;
        for (auto& phase : seeded.problem.phases) {
            const bool guessed = phase.composition_guess.size() == phase.support.size();
            if (phase.is_zero_pinned() && guessed) continue;
            const std::size_t source = phase.is_zero_pinned() ? 0 : index++;
            phase.composition_guess.clear();
            phase.declared_composition_guess = false;
            for (const std::size_t i : phase.support) phase.composition_guess.push_back(solved.mole_fractions[source * components + i]);
        }
        return boundary_start(seeded, true, limits, &stages, {}, correct_dense);
    } catch (const std::exception& error) {
        core::retain_failure(stages.first_numerical_failure, error, "speciated start");
        return std::nullopt;
    }
}

}  // namespace

TopologyCheck check_topology(const CompiledSystem& compiled, const std::vector<double>& full) {
    TopologyCheck result;
    if (full.size() < compiled.primal_size) return result;
    if (compiled.ideal) {
        const auto at = compiled.prepare(full, true);
        result.compositions = at.composition; result.amounts = at.amount;
        return result;
    }
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
template <typename Model>
static void observe_endpoint(const Model& model, const CompiledSystem* compiled,
    EquilibriumResult& result, core::Scope& scope, std::uint32_t attempt = 0) noexcept {
    if (!scope.diagnostics_on() || !result.has_state() || result.numerical_estimates) return;
    scope.observe([&] {
        result.numerical_estimates = compiled ? equilibrium_numerical_estimates(*compiled, result)
            : equilibrium_numerical_estimates(model, result);
        const auto& e = *result.numerical_estimates;
        scope.probe("solved_endpoint_jacobian", {static_cast<double>(e.dimension), static_cast<double>(e.dimension),
            static_cast<double>(e.rank), e.condition}, attempt, e.cause.c_str());
    });
}

EquilibriumResult solve_equilibrium(const HelmholtzModel& model, const EquilibriumProblem& request,
                                    const std::optional<EquilibriumResult>& start) {
    core::Scope scope("solve_equilibrium");  // compilation belongs to the solve's span
    std::optional<CompiledSystem> compiled;
    if (!request.discovery) compiled = compile_problem(model, request);
    auto result = request.discovery ? discover_phases(model, request) : solve_compiled(*compiled, start, &scope, nullptr,
        start && start->success && start->requested_tolerance_met.value_or(false) ? start->optimizer_state : nullptr, !start, true);
    if (scope.diagnostics_on() && !scope.nested("solve_equilibrium"))
        observe_endpoint(model, compiled ? &*compiled : nullptr, result, scope);
    return result;
}

EquilibriumResult solve_equilibrium(const IdealMoleFraction& model, const EquilibriumProblem& request,
                                    const std::optional<EquilibriumResult>& start) {
    if (start) throw std::invalid_argument("ideal previous-result starts unsupported before calculation");
    core::Scope scope("solve_equilibrium");
    const auto compiled = compile_problem(model, request);
    auto result = solve_compiled(compiled, std::nullopt, &scope);
    observe_endpoint(model, &compiled, result, scope);
    return result;
}

bool newton_eligible(const CompiledSystem& compiled) {
    const auto& p = compiled.problem;
    if (!compiled.feasibility_cell() || compiled.unknown_count() != compiled.rows.size()
        || compiled.has_caloric_condition() || p.critical || p.spinodal || p.envelope || p.azeotrope || p.discovery
        || compiled.has_reaction || compiled.has_ions || p.phases.size() != 2 || compiled.free_intensive_count != 1) return false;
    const auto& a = p.phases[0]; const auto& b = p.phases[1];
    return a.kind != b.kind && ((a.is_zero_pinned() && b.is_free()) || (b.is_zero_pinned() && a.is_free()));
}

static bool free_amount_newton_eligible(const CompiledSystem& compiled, bool discovery = false) {
    const auto& p = compiled.problem;
    if (compiled.feasibility_cell() || compiled.unknown_count() != compiled.rows.size()
        || compiled.free_intensive_count != 0 || !p.temperature.is_positive() || !p.pressure.is_positive()
        || compiled.has_caloric_condition() || p.critical || p.spinodal || p.envelope || p.discovery || p.source_reference
        || compiled.has_reaction || compiled.has_ions
        || (p.phases.size() != 2 && !(discovery && p.phases.size() == 1))) return false;
    if (std::any_of(compiled.model->parameters().charge_number.begin(), compiled.model->parameters().charge_number.end(),
        [](double z) { return z != 0.0; })) return false;
    const auto& a = p.phases[0]; const auto& b = p.phases.back();
    if (!a.is_free() || !b.is_free() || a.support != b.support || (p.phases.size() == 2 && (a.support.size() < 2
        || (a.kind == PhaseKind::Vapor && b.kind == PhaseKind::Vapor)))) return false;
    for (std::size_t i = 0; i < compiled.n_feed.size(); ++i)
        if ((compiled.n_feed[i] > 0.0) != (std::find(a.support.begin(), a.support.end(), i) != a.support.end())) return false;
    return true;
}

bool discovery_newton_eligible(const CompiledSystem& compiled) {
    return !compiled.ideal && !compiled.problem.azeotrope && compiled.problem.solids.empty()
        && !std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE") && free_amount_newton_eligible(compiled, true);
}

EquilibriumResult newton_attempt(const CompiledSystem& compiled, const std::vector<double>& initial,
    const EquilibriumResult& structural, const DriverLimits* limits, std::vector<double>* attempt) {
    core::Scope scope("newton_fast_path");
    auto fast = structural;
    std::vector<double> coordinates;
    try {
        fast = run_driver(compiled, initial, structural, scope, &coordinates, limits, Attempt::Newton);
        const auto branch = fast.success ? check_density_branches(compiled, fast, {}, true, scope, limits) : BranchCheck{};
        const bool accepted = fast.success && branch.accepted;
        if (fast.success && !accepted) core::retain_failure(fast.first_numerical_failure, branch.message, "Newton admission");
        scope.probe("newton_fast_path_outcome", {static_cast<double>(fast.evaluations)}, 0,
            accepted ? "accepted" : fast.success ? "admission" : fast.message.c_str());
        if (accepted) { if (attempt) *attempt = std::move(coordinates); return fast; }
    } catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) {
        core::retain_failure(fast.first_numerical_failure, error, "Newton attempt");
        scope.probe("newton_fast_path_outcome", {static_cast<double>(fast.evaluations)}, 0, error.what());
    }
    refuse(fast);
    return fast;
}

// The #565 recovery domain: a neutral nonreacting non-associating L+V boundary or TP flash on the full positive feed,
// never an azeotrope cell.
bool recovery_eligible(const CompiledSystem& compiled) {
    if (compiled.ideal) return false;
    const auto& p = compiled.problem;
    const auto components = compiled.model->parameters().component_count();
    return (newton_eligible(compiled) || free_amount_newton_eligible(compiled)) && p.phases[0].kind != p.phases[1].kind
        && !p.source_reference && !p.azeotrope && compiled.model->parameters().site_component.empty()
        && p.phases[0].support.size() == components && p.phases[1].support.size() == components
        && std::all_of(compiled.n_feed.begin(), compiled.n_feed.end(), [](double n) { return n > 0.0; });
}

EquilibriumResult continue_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start,
    core::Scope& scope, const DriverLimits* limits, std::shared_ptr<const OptimizerState> warm,
    bool recover, const EquilibriumResult& structural, StageWork stages,
    std::optional<std::vector<double>> original, std::shared_ptr<const EquilibriumResult> rejected, bool select_endpoint) {
    const auto& model = *compiled.model;
    const bool pure_source = start && start->success && limits && limits->recovery && one_component_two_density_eligible(compiled);
    const bool use_start = pure_source || (start.has_value() && usable_start(compiled, *start));
    // A previous result that start_point copies replaces the start's composition; skip its speciation.
    const bool copied = start.has_value() && start->phase_count == compiled.problem.phases.size()
        && valid_previous_phase(*start, 0, model.parameters().component_count());
    if (compiled.problem.spinodal && start && start->has_state() && start->spinodal) {
        try { spinodal_search(compiled, limits, start->molar_densities[0], &start->spinodal->search); }
        catch (const std::bad_alloc&) { throw; }
        catch (const std::exception& error) { return spinodal_search_failure(compiled, structural, error.what()); }
    }
    auto cold = start && !use_start ? std::make_optional(compiled) : std::nullopt;
    if (cold && start->phase_count == cold->problem.phases.size()) {
        for (std::size_t phase = 0; phase < start->phase_count; ++phase) {
            if (valid_previous_phase(*start, phase, model.parameters().component_count()))
                cold->problem.phases[phase].declared_composition_guess = false;
        }
    }
    const auto selected = [&](EquilibriumResult endpoint, const DriverLimits* remaining) {
        if ((!select_endpoint && !recover) || scope.nested("solve_equilibrium")) return endpoint;
        return refine_selected_endpoint(compiled, std::move(endpoint), scope, remaining);
    };
    auto seeded = start.has_value() || original ? std::nullopt : speciated_start(compiled, limits, stages, true);
    EquilibriumResult result;
    bool admission_refused = false;
    // The speciated seed when there is one; the caller's start alone or after the seed's failure.
    for (int attempt = seeded ? 0 : 1; attempt < 2 && !result.success; ++attempt) {
        stages.add(result);  // a failed seeded attempt
        std::vector<double> initial;
        if (original) initial = *original;
        else if (attempt == 0) {
            initial = std::move(*seeded);
        } else {
            try { initial = use_start ? compiled.generic_point() : boundary_start(cold ? *cold : compiled, !copied, limits, &stages); }
            catch (const std::bad_alloc&) { throw; }
            catch (const std::exception& error) {
                core::retain_failure(stages.first_numerical_failure, error, "solve initialization");
                if (compiled.problem.spinodal && compiled.spinodal_search) {
                    result = spinodal_search_failure(compiled, structural, error.what());
                    std::move(stages).report(result); return result;
                }
                result = structural;
                stages.report_observations(result);
                result.message = error.what();
                result.failure = std::current_exception();
                return result;
            }
            if (start.has_value()) initial = start_point(compiled, *start, std::move(initial));
            if (pure_source) {
                const double shift = std::log(model.packing_fraction_coefficient(compiled.seed_temperature(), {1.0})
                    / model.packing_fraction_coefficient(*start->temperature, {1.0}));
                const auto offsets = compiled.phase_offsets();
                for (std::size_t p = 0; p < offsets.size(); ++p)
                    initial[offsets[p] + 1 + (compiled.problem.phases[p].is_free() ? 1 : 0)] += shift;
            }
        }
        if (!original && (newton_eligible(compiled) || (!warm && !std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE") && free_amount_newton_eligible(compiled)))) {
            auto fast = one_component_two_density_eligible(compiled)
                ? one_component_two_density_attempt(compiled, initial, structural, scope, limits)
                : newton_attempt(compiled, initial, structural, limits);
            if (one_component_two_density_eligible(compiled) && !fast.success)
                rejected = std::make_shared<const EquilibriumResult>(fast);
            if (fast.success) { fast.rejected_reduction = rejected; std::move(stages).report(fast); return selected(std::move(fast), limits); }
            stages.add(fast);
        }
        result = run_driver(compiled, std::move(initial), structural, scope, nullptr, limits, Attempt::Cold, warm, false, false, Retry::None);
        if (result.failure && attempt == 0) continue;  // the seeded attempt's failure falls back to the declared start
        if (result.failure) { stages.report_observations(result); return result; }
        if (!result.success) continue;
        BranchCheck branches;
        try {
            branches = compiled.mode_count ? admit_null_mode(compiled, result)
                : check_density_branches(compiled, result, {}, true, scope, limits);
        } catch (const std::bad_alloc&) { throw; }
        catch (const std::exception& error) {
            branches.message = std::string("density branch admission failed: ") + error.what();
        }
        if (!branches.accepted) {
            admission_refused = true;

            result.message = branches.message;
            refuse(result);
        }
    }
    std::move(stages).report(result);
    const bool retry = result.retry == Retry::PhysicalAmounts && !admission_refused;
    result = recover_physical_amounts(compiled, std::move(result), retry,
        scope, nullptr, limits, nullptr, true);
    result.rejected_reduction = rejected;
    // ALG-EQUILIBRIUM recovery work: its own caps, within any caller remainder.
    const auto allowance = [&] {
        const auto cap = static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value);
        const auto left = limits ? limits->iterations - std::min(limits->consumed, limits->iterations) : cap;
        const auto deadline = limits ? limits->deadline : std::chrono::steady_clock::time_point::max();
        return DriverLimits::recovery_allowance(deadline, std::min(cap, left));
    };
    // ALG-STABILITY: a vanished declared phase at fixed T and P is named from the surviving declaration under that
    // allowance; the declared split stays refused.
    if (!result.success && recover && result.topology_event == TopologyEvent::PhaseDisappeared && result.topology_phase
        && compiled.problem.phases.size() > 1 && !compiled.problem.temperature.free && !compiled.problem.pressure.free
        && !compiled.ideal && !model.parameters().ideal_gas && compiled.problem.solids.empty()) {
        auto probe = allowance();
        auto report = vanished_phase_stability(compiled, *result.topology_phase, probe);
        if (limits) limits->consumed += probe.consumed;
        // Driver plus tangent-plane work, as the report keeps them; the original solve's work is the remainder.
        auto& w = report.work;
        result.iterations += static_cast<std::size_t>(w["iterations"] + w["trial_tangent_plane_iterations"]);
        result.evaluations += static_cast<std::size_t>(w["evaluations"] + w["trial_tangent_plane_evaluations"]);
        result.hessian_calls += static_cast<std::size_t>(w["hessian_calls"] + w["trial_tangent_plane_hessian_evaluations"]);
        result.message += "; " + report.message;
        result.vanished_phase_stability = std::make_shared<const StabilityReport>(std::move(report));
    }
    if (result.success || !recover || one_component_two_density_eligible(compiled) || admission_refused || (limits && (limits->expired() || limits->iterations_spent()))
        || !recovery_eligible(compiled) || result.retry != Retry::Source) return selected(std::move(result), limits);
    // Neutral recovery: both causes and all work are kept.
    auto work = allowance();
    // ALG-STABILITY (#987): after coalescence at fixed T and P, the finite trial search under the same allowance
    // construction; no descent in it ends the solve without recovery. It is not a global stability certificate.
    if (result.topology_event == TopologyEvent::PhasesCoalesced && !compiled.problem.temperature.free
        && !compiled.problem.pressure.free && !model.parameters().ideal_gas && compiled.problem.solids.empty()) {
        auto probe = work;  // the same allowance construction; recovery's remainder shrinks by its use below
        result.recovery_stability = std::make_shared<const StabilityReport>(recovery_stability(compiled, probe));
        if (limits) limits->consumed += probe.consumed;
        if (limits) work.iterations = std::min(work.iterations, limits->iterations - limits->consumed);
        result.iterations += probe.consumed;
        result.evaluations += probe.evaluations;
        result.hessian_calls += probe.hessian_calls;
        result.message += "; " + result.recovery_stability->message;
        if (result.recovery_stability->no_descent_in_trial_set) {
            result.message += "; not a global stability certificate: a split outside this finite trial set is not excluded; recovery not run";
            return selected(std::move(result), limits);
        }
    }
    auto recovered = selected(recover_from_source(compiled, work), &work);
    StageWork observations; observations.observe(result); observations.observe(recovered);
    std::move(observations).publish_observations(recovered.success ? recovered : result);
    if (limits) limits->consumed += work.consumed;
    result.iterations += work.consumed;
    result.evaluations += work.evaluations;
    result.hessian_calls += work.hessian_calls;
    if (!recovered.success) { result.message += "; recovery: " + recovered.message; return result; }
    recovered.rejected_reduction = rejected;
    recovered.recovery_stability = result.recovery_stability;
    recovered.message += "; recovered after: " + result.message;
    recovered.iterations = result.iterations;
    recovered.evaluations = result.evaluations;
    recovered.hessian_calls = result.hessian_calls;
    return recovered;
}

// Reduced temperature and packing starts; coefficients are initialization heuristics (ALG-EQUILIBRIUM).
static std::vector<std::array<double, 2>> critical_seed_scales(const HelmholtzModel& model) {
    const auto& p = model.parameters();
    const double theta = p.dispersion_energy_over_k[0], base = theta * (1 + 0.35 * p.segment_count[0]);
    double energy = 0;
    for (std::size_t i = 0; i < p.association_active.size(); ++i)
        if (p.association_active[i] && p.association_volume[i] > 0) energy = std::max(energy, p.association_energy_over_k[i]);
    std::vector<std::array<double, 2>> seeds{{4 * theta, 0.025}, {3 * theta, 0.04}, {base, 0.12}, {1.5 * base, 0.072}};
    if (energy > 0) {
        const double associated = std::max(base, energy / 4);
        seeds.push_back({associated, 0.025});
        seeds.push_back({associated, 0.12});
    }
    return seeds;
}

// One driver, finite physical starts, and no selection between distinguishable critical candidates.
static EquilibriumResult recover_critical(const CompiledSystem& compiled, EquilibriumResult original, core::Scope& scope) {
    auto search = std::make_shared<CriticalSearch>();
    search->seeds.push_back({compiled.seed_temperature(), compiled.problem.critical->density_guess_mol_per_m3});
    search->results.push_back(original);
    search->added_work.push_back({0, 0, 0});
    const auto seeds = critical_seed_scales(*compiled.model);
    search->scheduled = seeds.size();
    const auto cap = static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value);
    auto work = DriverLimits::recovery_allowance(std::chrono::steady_clock::time_point::max(), cap);
    StageWork stages;
    stages.add(original);
    bool comparable = true;
    for (const auto& seed : seeds) {
        if (work.expired() || work.caps_spent()) break;
        auto declaration = compiled.problem;
        declaration.temperature = FreeValue::variable(seed[0]);
        declaration.critical->mode_guess.clear();
        EquilibriumResult result;
        double rho = NAN;
        try {
            rho = seed[1] / compiled.model->packing_fraction_coefficient(seed[0], compiled.feed_composition);
            declaration.critical->density_guess_mol_per_m3 = rho;
            const auto cell = compile_problem(*compiled.model, declaration);
            result = try_solve_compiled(cell, {}, &scope, &work);
            if (result.failure) with_failure(result, [&](const std::exception& error) {
                core::retain_failure(result.first_numerical_failure, error, "critical start");
                result.message = std::string("; critical start: ") + error.what(); });
            if (result.success && work.spend(work.evaluations, work.evaluation_cap) && !work.expired()) {
                result.numerical_estimates = equilibrium_numerical_estimates(cell, result);
                work.check();
            }
        } catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) {
            core::retain_failure(result.first_numerical_failure, error, "critical start");
            result.message += std::string("; critical start: ") + error.what();
        }
        stages.observe(result);
        search->seeds.push_back({seed[0], rho});
        search->results.push_back(std::move(result));
        search->added_work.push_back({work.consumed, work.evaluations, work.hessian_calls});
        const auto index = search->results.size() - 1;
        const auto& current = search->results.back();
        if (!current.has_state()) continue;
        const auto& estimate = current.numerical_estimates;
        if (!estimate || !estimate->cause.empty() || !estimate->temperature.resolved || !estimate->pressure.resolved
            || estimate->densities.empty() || !estimate->densities[0].resolved) comparable = false;
        if (!search->selected) search->selected = index;
        if (!comparable) continue;
        const auto& b = *estimate;
        const auto agrees = [](equilibrium::Quantity quantity, double x, double y, const equilibrium::NumericalEstimate& ex, const equilibrium::NumericalEstimate& ey) {
            return equilibrium::compare_equilibrium_quantity(quantity, x, y, 0, ex, ey, {}, equilibrium::ComparisonMode::EstimatedError).judgment != "failed";
        };
        for (std::size_t i = 1; comparable && i < index; ++i) {
            const auto& previous = search->results[i];
            if (!previous.has_state()) continue;
            const auto& a = *previous.numerical_estimates;
            comparable = agrees(equilibrium::Quantity::Linear, *previous.temperature, *current.temperature, a.temperature, b.temperature)
                && agrees(equilibrium::Quantity::Linear, *previous.pressure, *current.pressure, a.pressure, b.pressure)
                && agrees(equilibrium::Quantity::Positive, previous.molar_densities[0], current.molar_densities[0], a.densities[0], b.densities[0]);
        }
    }
    search->budget_exhausted = work.expired() || work.caps_spent();
    if (!comparable) search->selected.reset();
    auto result = search->selected ? search->results[*search->selected] : std::move(original);
    if (!search->selected) { refuse(result); }
    result.message += search->selected ? "; finite critical start search" : "; critical starts unresolved or no admitted candidate";
    if (search->budget_exhausted) result.message += "; critical start allowance exhausted";
    result.iterations = work.consumed;
    result.evaluations = work.evaluations;
    result.hessian_calls = work.hessian_calls;
    std::move(stages).report(result);
    result.critical_search = std::move(search);
    return result;
}

EquilibriumResult try_solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start,
                                 core::Scope* outer, const DriverLimits* limits, std::shared_ptr<const OptimizerState> warm,
                                 bool recover, bool select_endpoint) {
    std::optional<core::Scope> own;
    auto& scope = outer != nullptr ? *outer : own.emplace("solve_equilibrium");
    if (compiled.ideal && start) throw std::invalid_argument("ideal previous-result starts unsupported before calculation");
    const EquilibriumResult structural = structural_result(compiled, scope);
    if (!compiled.feasible) return structural;
    if (compiled.ideal) return run_driver_once(compiled, compiled.generic_point(), structural, scope,
        nullptr, limits, Attempt::Cold, {}, false, false);
    const bool pure_recovery = recover && !limits && !start && !warm && !std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE") && one_component_two_density_eligible(compiled);
    const bool critical = recover && !start && !limits && !warm && compiled.problem.critical
        && compiled.n_feed.size() == 1 && !std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE");
    auto result = [&] {
        try { return continue_compiled(compiled, start, scope, limits, std::move(warm), recover, structural, {}, std::nullopt, {}, select_endpoint); }
        catch (const std::bad_alloc&) { throw; } catch (const std::exception& error) {
            auto failed = structural;
            failed.message = error.what();
            failed.failure = std::current_exception();
            return failed;
        }
    }();
    if (pure_recovery) {
        if (result.failure) with_failure(result, [&](const std::exception& error) {
            core::retain_failure(result.first_numerical_failure, error, "saturation initialization");
            if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error)) result.terminal_numerical_failure = failure->record;
            refuse(result);
            if (!result.last_full_evaluation) result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        });
        return recover_from_source(compiled, std::move(result));
    }
    if (critical && result.failure) with_failure(result, [&](const std::exception& error) {
        core::retain_failure(result.first_numerical_failure, error, "critical initialization"); });
    return critical && !result.success ? recover_critical(compiled, std::move(result), scope) : result;
}

EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start,
                                 core::Scope* outer, const DriverLimits* limits, std::shared_ptr<const OptimizerState> warm,
                                 bool recover, bool select_endpoint) {
    auto result = try_solve_compiled(compiled, start, outer, limits, std::move(warm), recover, select_endpoint);
    if (result.failure) std::rethrow_exception(result.failure);
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
                const auto offsets = compiled.phase_offsets();
                const auto check = [](bool ok, const std::string& text) { if (!ok) throw std::invalid_argument(text); };
                const auto field = [](const std::string& text) { return "invalid declared seed: " + text; };
                check((initial.empty() && !seed.previous) || (seed.compositions.empty() && seed.packing_fractions.empty()
                    && seed.amounts.empty() && !seed.temperature && !seed.pressure && !(seed.previous && !initial.empty())),
                    field("coordinates or a previous state exclude every other seed field"));
                if (seed.previous) { check(seed.previous->has_state(), field("previous state has no solved state")); initial = seed.previous->coordinates; }
                else if (initial.empty()) {
                    auto seeded = compiled;
                    check(seed.compositions.empty() || seed.compositions.size() == compiled.problem.phases.size(),
                          field("compositions: one entry per declared phase"));
                    check(seed.packing_fractions.empty() || seed.packing_fractions.size() == compiled.problem.phases.size(),
                          field("packing_fractions: one entry per declared phase"));
                    check(seed.amounts.empty() || seed.amounts.size() == compiled.problem.phases.size(),
                          field("amounts: one entry per declared phase"));
                    // Complete: the seed supplies every coordinate the plain solve's cold start would generate.
                    bool uncomposed = true, complete = (!seed.amounts.empty() || compiled.problem.phases.size() == compiled.pinned_count)
                        && (seed.temperature || !problem.temperature.free) && (seed.pressure || !problem.pressure.free);
                    const auto composed = [&](std::size_t p) { return !seed.compositions.empty() && !seed.compositions[p].empty(); };
                    const auto packed = [&](std::size_t p) { return !seed.packing_fractions.empty() && seed.packing_fractions[p]; };
                    for (std::size_t p = 0; p < seeded.problem.phases.size(); ++p) {
                        auto& phase = seeded.problem.phases[p];
                        if (composed(p)) {
                            check(seed.compositions[p].size() == compiled.declared_supports[p].size()
                                  && std::all_of(seed.compositions[p].begin(), seed.compositions[p].end(), [](double x) { return std::isfinite(x) && x >= 0.0 && x <= 1.0; }),
                                  field("composition for phase " + phase.name + ": one fraction in [0, 1] per declared species"));
                            phase.composition_guess = seed.compositions[p];
                            phase.declared_composition_guess = false;
                            restrict_composition_guess(phase, compiled.declared_supports[p]);
                        }
                        check(!model.parameters().ideal_gas || !packed(p), "ideal-gas packing seeds are undefined");
                        if (packed(p)) phase.packing_guess = seed.packing_fractions[p];
                        check(phase.composition_guess.empty() || phase.composition_guess.size() == phase.support.size(),
                              field("composition for phase " + phase.name + " does not match its species support"));
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
                    if (seed.temperature) { check(problem.temperature.free, field("temperature seed on a fixed temperature")); seeded.problem.temperature = FreeValue::variable(*seed.temperature); }
                    if (seed.pressure) { check(problem.pressure.free, field("pressure seed on a fixed pressure")); seeded.problem.pressure = FreeValue::variable(*seed.pressure); }
                    work.check("elapsed work budget exhausted in initialization");
                    auto speciated = complete || !uncomposed ? std::nullopt : speciated_start(seeded, &work, stages);
                    std::vector<bool> packings(seeded.problem.phases.size());
                    for (std::size_t p = 0; p < packings.size(); ++p) packings[p] = packed(p);
                    initial = complete ? seeded.generic_point() : speciated ? std::move(*speciated)
                        : boundary_start(seeded, uncomposed, &work, &stages, packings);
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
                check(initial.size() == compiled.unknown_count() && all_numeric(initial), field("coordinates: one finite value per unknown"));
                for (std::size_t p = 0; p < offsets.size(); ++p) {
                    const auto& phase = compiled.problem.phases[p]; auto cursor = offsets[p];
                    for (std::size_t j = 0; j < phase.support.size(); ++j) {
                        const double x = std::exp(initial[cursor++]);
                        check(x > 0.0, field("composition coordinate for phase " + phase.name + " is not positive"));
                    }
                    if (phase.is_free()) check(initial[cursor++] >= 0.0, field("amount for phase " + phase.name + " is negative"));
                    check(model.parameters().ideal_gas ? std::isfinite(std::exp(initial[cursor])) && std::exp(initial[cursor]) > 0.0
                          : admitted_packing_fraction(std::exp(initial[cursor])),
                          field("packing for phase " + phase.name + " is outside the admitted range"));
                }
                for (auto j = compiled.primal_size - compiled.free_intensive_count; j < compiled.primal_size; ++j)
                    check(initial[j] > 0.0, field("free temperature or pressure is not positive"));
                audit = structural_result(compiled, scope);
                audit->convergence.residual = std::numeric_limits<double>::quiet_NaN();
                if (compiled.feasible) audit = run_driver(compiled, std::move(initial), *audit, scope, nullptr, &work);
                if (audit->failure) std::rethrow_exception(audit->failure);  // the attempt's one failure handler below
                auto result = *audit;
                StageWork(stages).report(result);
                if (work.expired()) { result.message += "; elapsed work budget exhausted";
                    reject_elapsed(result, scope); }
                if (result.success) {
                    BranchCheck branch;
                    try {
                        branch = compiled.mode_count ? admit_null_mode(compiled, result)
                              : check_density_branches(compiled, result, {}, true, scope, &work);
                        if (!branch.accepted || work.expired()) {
                            result.message = (work.expired() ? "elapsed work budget exhausted in admission" : branch.message) + "; driver: " + result.message;
                            if (work.expired()) reject_elapsed(result, scope); else refuse(result);
                        }
                    } catch (const std::exception& error) {
                        const bool budget = dynamic_cast<const BudgetExhausted*>(&error);
                        result.message = std::string(budget ? "elapsed work budget exhausted in admission" : error.what()) + "; driver: " + result.message;
                        if (budget) reject_elapsed(result, scope); else refuse(result);
                    }
                }
                scope.finish(result.convergence);
                out.causes[i] = result.message; out.results[i] = std::move(result);
                if (out.results[i]->success && out.results[i]->requested_tolerance_met.value_or(false)) out.selected = i;
            } catch (const std::exception& error) {
                const bool budget = dynamic_cast<const BudgetExhausted*>(&error);
                out.causes[i] = std::string(audit ? "driver/post-driver: " : "initialization: ") + error.what();
                if (budget) out.causes[i] += "; elapsed deadline exhausted";
                if (audit) { audit->message = out.causes[i]; std::move(stages).report(*audit);
                    if (budget) reject_elapsed(*audit, scope); else refuse(*audit);
                    out.results[i] = std::move(audit); }
            }
            if (out.results[i]) for (auto& rank : out.results[i]->rank_records)
                if (rank.first_conditioning_warning && !rank.first_conditioning_warning->attempt) (rank.first_conditioning_warning = std::make_shared<core::ConditioningObservation>(*rank.first_conditioning_warning))->attempt = i;
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
    const HelmholtzModel& model, const EquilibriumProblem& request, const std::vector<double>& coordinates,
    std::vector<double>* attempt, const DriverLimits* limits, std::shared_ptr<const OptimizerState> warm
) {
    const auto compiled = compile_problem(model, request);
    return solve_equilibrium_from_coordinates(compiled, coordinates, attempt, limits, std::move(warm));
}

EquilibriumResult solve_equilibrium_from_coordinates(
    const CompiledSystem& compiled, const std::vector<double>& coordinates, std::vector<double>* attempt,
    const DriverLimits* limits, std::shared_ptr<const OptimizerState> warm
) {
    core::Scope scope("solve_equilibrium");
    EquilibriumResult result = structural_result(compiled, scope);
    if (!compiled.feasible) return result;
    if (coordinates.size() != compiled.unknown_count() || !all_numeric(coordinates))
        throw std::invalid_argument("continuation coordinates do not match the compiled layout");
    return run_driver(compiled, coordinates, std::move(result), scope, attempt, limits, Attempt::Cold, std::move(warm), true);
}

}  // namespace epcsaft::equilibrium
