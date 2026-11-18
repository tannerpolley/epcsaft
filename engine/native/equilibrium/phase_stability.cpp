#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "contributions/helmholtz.hpp"

#include "core/trace.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cstring>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <variant>

// ALG-STABILITY: the explicit tangent-plane test and fixed-T/P phase discovery.
namespace epcsaft::equilibrium {

double present_phase_residual_distance(double chemical_contraction, double pressure_row) {
    return chemical_contraction - pressure_row;
}
bool descent_resolved(double tpd, double slope, double error_estimate) {
    return std::isfinite(tpd) && std::isfinite(slope) && std::isfinite(error_estimate)
        && error_estimate >= 0.0 && tpd < 0.0 && slope < -error_estimate;
}
namespace {
using Quantities = std::map<std::string, double>;

void report_work(Quantities& work, const core::CalculationWork& interval) {
    const auto counts = interval.counts();
    work["root_calls"] = counts.root_calls;
    work["root_evaluations"] = counts.root_evaluations;
    work["energy_evaluations"] = counts.energy_evaluations;
}

void count_driver(Quantities& work, const EquilibriumResult& result, const char* role) {
    ++work["driver_invocations"]; ++work[role];
    work["iterations"] += result.iterations;
    work["evaluations"] += result.evaluations;
    work["hessian_calls"] += result.hessian_calls;
}

DensityRoot stability_root(const HelmholtzModel& model, double T, double P, const std::vector<double>& w, PhaseKind kind) {
    return solve_density(model, T, P, w, kind == PhaseKind::Liquid ? DensityBranch::Liquid : DensityBranch::Vapor, std::nullopt);
}

void packing_admission(double eta) {
    if (!admitted_packing_fraction(eta))
        throw std::invalid_argument("stability point outside packing domain");
}

// ALG-STABILITY, EqID michelsen_modified_tangent_plane (all-neutral, non-reacting systems); empty: the driver fallback.
// `limits` (#987) is charged before every iteration, point and Hessian; a refusal ends the trial with its work so far.
using TrialPoint = std::pair<PhasePoint, double>;  // composition and packing fraction; dual -ln sum W
std::optional<TrialPoint> tangent_plane_trial(const HelmholtzModel& model, double T, double P, const PhaseSpec& phase,
    const std::vector<double>& seed, double seed_density, const std::vector<double>& shift, StabilityTrial& trial, Quantities& work,
    const DriverLimits* limits = nullptr) {
    // DriverLimits::spend reads no clock: the deadline is checked first.
    const auto charge = [limits](std::size_t& used, std::size_t cap) {
        if (limits->check(); !limits->spend(used, cap)) limits->check();
    };
    const auto& S = phase.support;
    const auto n = static_cast<Eigen::Index>(S.size());
    struct Iterate { Eigen::VectorXd W, g; std::vector<double> w; double density = 0.0, objective = 0.0; };
    auto& o = trial.optimizer;
    bool first = true;
    const auto at = [&](Eigen::VectorXd W) {
        Iterate p{std::move(W), Eigen::VectorXd(n), std::vector<double>(seed.size(), 0.0)};
        for (Eigen::Index a = 0; a < n; ++a) p.w[S[a]] = p.W(a) / p.W.sum();
        if (limits) charge(limits->evaluations, limits->evaluation_cap);
        ++o.evaluations;
        if (first && p.w == seed) p.density = seed_density;
        else { const auto root = stability_root(model, T, P, p.w, phase.kind);
            core::retain_first(o.first_numerical_failure, root.first_numerical_failure); p.density = root.molar_density; }
        first = false;
        packing_admission(model.packing_fraction_coefficient(T, p.w) * p.density);
        const auto d = model.evaluate(T, 1.0 / p.density, p.w, DerivativeOrder::Gradient);
        core::retain_first(o.first_numerical_failure, d.first_numerical_failure);
        for (Eigen::Index a = 0; a < n; ++a)  // g_i = ln W_i + ln phi_i + ln(P/RT) + standard_i - plane_i
            p.g(a) = std::log(p.W(a)) + d.gradient[2 + S[a]] + std::log(p.density) + shift[S[a]];
        if (!p.g.allFinite()) throw std::domain_error("non-finite or non-positive Michelsen amounts");  // via ln W
        p.objective = 1.0 + p.W.dot(p.g) - p.W.sum();  // g* (eqs. 13-17)
        return p;
    };
    const auto& criterion = policy::EQUILIBRIUM_ROW_RESIDUAL;
    std::optional<TrialPoint> result;
    try {
        auto p = at(Eigen::Map<const Eigen::VectorXd>(seed.data(), static_cast<Eigen::Index>(seed.size()))(S));
        while (p.g.lpNorm<Eigen::Infinity>() > criterion.value && o.iterations < policy::EQUILIBRIUM_MAX_ITERATIONS.value) {
            std::optional<Iterate> q;
            if (limits) charge(limits->consumed, limits->iterations);
            if (o.iterations++ >= policy::TANGENT_PLANE_SUBSTITUTIONS.value) {
                if (limits) charge(limits->hessian_calls, limits->hessian_cap);
                ++work["trial_tangent_plane_hessian_evaluations"];
                const auto d = model.evaluate(T, 1.0 / p.density, p.w, DerivativeOrder::IsothermalPhaseHessian);
                core::retain_first(o.first_numerical_failure, d.first_numerical_failure);
                const Eigen::MatrixXd J = Eigen::Map<const Eigen::MatrixXd>(
                    tangent_plane_jacobian(d, 1.0 / p.density, S).data(), n, n);
                const Eigen::VectorXd root = p.W.array().sqrt();
                Eigen::MatrixXd H = J.cwiseProduct(root * root.transpose()) / p.W.sum();
                H.diagonal().array() += 1.0 + 0.5 * p.g.array();
                const Eigen::LLT<Eigen::MatrixXd> factor(H);
                if (factor.info() != Eigen::Success) core::retain_first(o.first_numerical_failure, [&] {
                    core::Convergence r; r.cause = "tangent-plane matrix not positive definite";
                    r.residual = p.g.lpNorm<Eigen::Infinity>(); r.iterations = o.iterations;
                    return r;
                });
                if (factor.info() == Eigen::Success) {
                    const Eigen::VectorXd alpha = 2.0 * root - factor.solve(root.cwiseProduct(p.g));
                    if ((alpha.array() > 0.0).all()) q = at((0.5 * alpha).array().square());
                }
            }
            const bool accepted = q && q->objective <= p.objective;
            o.hessian_calls += accepted;
            p = accepted ? std::move(*q) : at(p.W.array() * (-p.g.array()).exp());
        }
        const double residual = p.g.lpNorm<Eigen::Infinity>();
        const auto stop = core::classify(residual, criterion, std::nullopt, std::nullopt, true);  // converged or the cap
        o.convergence = core::Convergence{residual, {std::shared_ptr<const policy::Tolerance>{}, &criterion}, std::nullopt, std::nullopt, std::nullopt, o.iterations,
            o.evaluations, stop};
        if (o.convergence.stop == core::Stop::converged) {
            o.success = true; o.message = trial.message = "Michelsen tangent-plane trial " + core::describe(o.convergence);
            o.packing_fractions = {model.packing_fraction_coefficient(T, p.w) * p.density};
            for (auto i : S) trial.raw_coordinates.push_back(std::log(p.w[i]));
            trial.raw_coordinates.insert(trial.raw_coordinates.end(), {1.0, std::log(o.packing_fractions[0])});
            PhasePoint point; point.composition = {p.w}; point.eta = o.packing_fractions;
            result.emplace(std::move(point), -std::log(p.W.sum()));
        }
    } catch (const std::exception& error) { retain_density_failure(o.first_numerical_failure, error, "tangent-plane trial"); }
    ++work["trial_tangent_plane_invocations"];
    work["trial_tangent_plane_iterations"] += o.iterations; work["trial_tangent_plane_newton_steps"] += o.hessian_calls;
    work["trial_tangent_plane_evaluations"] += o.evaluations;
    return result;
}

bool same_stationary_trial(const StabilityTrial& a, const StabilityTrial& b, bool same_kind = true) {
    return a.phase.support == b.phase.support && (!same_kind || a.phase.kind == b.phase.kind)
        && largest_difference(a.composition, b.composition) <= policy::PHASE_COINCIDENCE_DISTANCE.value
        && std::abs(std::log(a.quantities.at("volume_m3_per_mol") / b.quantities.at("volume_m3_per_mol")))
            <= policy::PHASE_COINCIDENCE_DISTANCE.value;
}

// The density owner's stable root of w within the anchor window of v, else none (#423).
std::optional<DensityRoot> anchored_root(const HelmholtzModel& model, double T, double P, const std::vector<double>& w,
    double v, PhaseKind kind, Quantities& q, const std::string& prefix) {
    DensityRoot sibling;
    auto root = solve_density(model, T, P, w, kind == PhaseKind::Liquid ? DensityBranch::Liquid : DensityBranch::Vapor,
        std::nullopt, &sibling);
    if (std::abs(std::log(root.molar_density * v)) > policy::DENSITY_ANCHOR_LOG_WINDOW.value) root = sibling;
    const double displacement = q[prefix + "root_log_displacement"] = std::abs(std::log(root.molar_density * v));
    q[prefix + "root_residual"] = root.convergence.residual;
    if (!(displacement <= policy::DENSITY_ANCHOR_LOG_WINDOW.value && root.convergence.residual <= root.convergence.criterion->value))
        return std::nullopt;
    return root;
}

constexpr auto density_refusal = "stability density branch unresolved";
std::optional<DensityRoot> retain_density_side(const HelmholtzModel& model, double T, double P, const std::vector<double>& w,
    double v, PhaseKind kind, Quantities& q) {
    const auto root = anchored_root(model, T, P, w, v, kind, q, "identity_");
    if (!root) return root;
    const auto side = properties::density_side(root->branch, root->molar_density, root->negative_slope_density);
    q["root_density_side"] = static_cast<double>(side);
    q["density_side"] = static_cast<double>(side == DensityBranch::Unspecified ? (kind == PhaseKind::Liquid ? DensityBranch::Liquid : DensityBranch::Vapor) : side);
    return root;
}

// Reduced energy and Euler/association estimates use unit-mol amount columns and the SI volume column.
struct PhaseEnergyPoint {
    HelmholtzDerivatives derivative;
    Quantities diagnostics;
};
PhaseEnergyPoint phase_energy_point(const HelmholtzModel& model, double T, double v,
    const std::vector<double>& w, Quantities& work) {
    packing_admission(model.packing_fraction_coefficient(T, w) / v);
    PhaseEnergyPoint out;
    out.derivative = model.evaluate(T, v, w, DerivativeOrder::IsothermalPhaseHessian);
    if (!model.parameters().site_component.empty()) {
        ++work["association_diagnostic_evaluations"];
        const auto diagnostics = model.density_actions(T, 1.0 / v, w,
            {PropertyObservable::ResidualEntropy}, {DensityActionDirection{0.0, 0.0,
                std::vector<double>(w.size(), 0.0), {}}}).diagnostics;
        out.diagnostics["association_residual"] = diagnostics.association_residual;
        out.diagnostics["association_site_pivot_ratio"] = diagnostics.site_pivot_ratio;
        out.diagnostics["association_site_rank"] = static_cast<double>(diagnostics.site_rank);
    }
    return out;
}
Quantities phase_energy(const HelmholtzModel& model, double T, double P, double v,
    const std::vector<double>& w, const std::vector<double>& standard, Quantities& work,
    HelmholtzDerivatives* retained = nullptr, const PhaseEnergyPoint* retained_point = nullptr) {
    const auto own = retained_point ? PhaseEnergyPoint{} : phase_energy_point(model, T, v, w, work);
    const auto& point = retained_point ? *retained_point : own;
    const auto& d = point.derivative;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const double ideal = ideal_mechanical_helmholtz(w, v, 1.0);
    double euler = v * d.gradient[1], magnitude = std::abs(d.value);
    for (std::size_t i = 0; i < w.size(); ++i) if (w[i] > 0.0) {
        euler += w[i] * d.gradient[2 + i];
        magnitude += std::abs(w[i] * (std::log(w[i] / v) - 1.0 + standard[i]));
    }
    const double pv = P * v / rt;
    const double value = d.value + ideal + pv
        + std::inner_product(w.begin(), w.end(), standard.begin(), 0.0);
    // gamma covers final summation only, not the EOS. Euler defect and inner error are separate.
    const double m = 2.0 * w.size() + 4.0;
    const double gamma = m * std::numeric_limits<double>::epsilon()
        / (1.0 - m * std::numeric_limits<double>::epsilon());
    Quantities q{{"g_over_rt", value}, {"volume_m3_per_mol", v},
        {"pressure_row", v * (1.0 / v - d.gradient[1] - P / rt)},
        {"gradient_floor_max_si_columns", d.gradient_floor},
        {"association_energy_estimate_over_rt", (1.0 + std::abs(v)) * d.gradient_floor},
        {"euler_defect_over_rt", std::abs(d.value - euler)},
        {"summation_estimate_over_rt", gamma * (magnitude + std::abs(pv))}};
    q.insert(point.diagnostics.begin(), point.diagnostics.end());
    if (retained) *retained = d;
    return q;
}

double energy_estimate(const Quantities& q) {
    return q.at("association_energy_estimate_over_rt") + q.at("euler_defect_over_rt")
        + q.at("summation_estimate_over_rt");
}

std::vector<PhaseSpec> stability_supports(const CompiledSystem& c, const StabilityOptions& options) {
    const auto& z = c.model->parameters().charge_number;
    std::vector<std::size_t> united, feasible, zero_charge;
    const auto participant = reaction_participants(c.n_feed.size(), c.problem.reactions);
    for (std::size_t i = 0; i < c.n_feed.size(); ++i)
        if (c.n_feed[i] > 0.0 || participant[i]) {
            united.push_back(i);
            if (std::find(c.removed_species.begin(), c.removed_species.end(), i) != c.removed_species.end()) continue;
            feasible.push_back(i); if (z[i] == 0.0) zero_charge.push_back(i);
        }
    std::vector<PhaseSpec> phases;
    if (options.trials) phases = *options.trials;
    else {
        PhaseSpec liquid; liquid.name = "trial-liquid"; liquid.support = feasible;
        phases.push_back(liquid);
        if (!zero_charge.empty()) {
            liquid.support = zero_charge; phases.push_back(liquid);
            liquid.name = "trial-vapor"; liquid.kind = PhaseKind::Vapor; phases.push_back(liquid);
        }
    }
    if (phases.empty()) throw std::invalid_argument("empty stability trial family");
    std::vector<PhaseSpec> admitted;
    for (auto phase : phases) {
        phase.support = resolve_support(phase, feasible, z, &c.model->parameters());
        bool positive = false, negative = false;
        for (auto i : phase.support) {
            if (std::find(united.begin(), united.end(), i) == united.end())
                throw std::invalid_argument("stability support outside fed/reaction union");
            positive |= z[i] > 0; negative |= z[i] < 0;
        }
        if (positive != negative)
            throw std::invalid_argument("inadmissible ionic stability support/label");
        phase.support_declared = true; phase.amount = AmountSpec::free();
        if (options.evaluate_supplied_trials || std::none_of(admitted.begin(), admitted.end(), [&](const auto& p) {
                return p.kind == phase.kind && p.support == phase.support; })) admitted.push_back(std::move(phase));
    }
    return admitted;
}

struct Candidate final {
    PhasePoint physical;
    std::vector<Quantities> quantities;
    double correction_norm = 0.0, feasibility_estimate = 0.0;
};

// Equality repair acts on supported amounts. Redundant charge/balance rows are solved together;
// the achieved defect, positivity and structural support decide admission, not a row-count guess.
Candidate stability_candidate(const CompiledSystem& c, const EquilibriumResult& result,
    StabilityReport& report) {
    if (!result.has_state() || result.model_identity != c.model || !result.declaration)
        throw std::invalid_argument("stability test requires this model's accepted declared result");
    if (const auto difference = declared_state_difference(*result.declaration, c.problem); !difference.empty())
        throw std::invalid_argument("stability result/declaration mismatch: " + difference);
    Candidate candidate;
    auto& at = candidate.physical;
    at = c.prepare(result.coordinates);
    const auto raw = c.evaluate_at(at);
    if (raw.values != result.residuals) throw std::invalid_argument("stability result/input identity mismatch");
    report.standard_over_rt = at.standard.derivative[0];
    report.plane_over_rt.assign(c.n_feed.size(), 0.0);
    for (std::size_t k = 0; k < c.rank; ++k)
        for (std::size_t i = 0; i < c.n_feed.size(); ++i)
            report.plane_over_rt[i] += c.B[k][i] * result.multipliers[k];
    std::vector<std::pair<std::size_t, std::size_t>> cells;
    for (std::size_t p = 0; p < at.amount.size(); ++p) if (at.amount[p] > 0.0)
        for (auto i : c.problem.phases[p].support) cells.emplace_back(p, i);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(c.rank + at.amount.size(), cells.size());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(J.rows());
    Eigen::VectorXd amount(cells.size()), mu(cells.size());
    for (std::size_t k = 0; k < c.rank; ++k) rhs(k) = c.balance[k];
    const auto& z = c.model->parameters().charge_number;
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        for (std::size_t k = 0; k < c.rank; ++k) J(k, j) = c.B[k][i];
        J(c.rank + p, j) = z[i];
        amount(j) = at.amount[p] * at.composition[p][i];
        mu(j) = std::log(at.composition[p][i] / at.volume[p])
            + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i];
        rhs -= J.col(j) * at.amount[p] * at.composition[p][i];
    }
    Eigen::MatrixXd weighted = J * amount.asDiagonal();
    Eigen::VectorXd scaled_rhs = rhs;
    for (Eigen::Index k = 0; k < weighted.rows(); ++k) {
        const double scale = weighted.row(k).cwiseAbs().maxCoeff();
        if (scale > 0.0) { weighted.row(k) /= scale; scaled_rhs(k) /= scale; }
    }
    const auto factor = weighted.completeOrthogonalDecomposition();
    report.last_rank = core::numerical_rank(factor, "weighted_amount_correction");
    const auto refuse = [&](const std::string& cause, double defect) {
        core::Convergence r; r.cause = cause; r.residual = defect; r.stop = core::Stop::infeasible; if (report.last_rank) r.rank = std::make_shared<core::NumericalRank>(*report.last_rank);
        core::retain_first(report.first_numerical_failure, r);
        throw core::DeclarationError(cause, r);
    };
    const Eigen::VectorXd intended = amount.cwiseProduct(factor.solve(scaled_rhs));
    if (!intended.allFinite())
        refuse("stability amount correction nonfinite response", NAN);
    const auto lost = [&](std::size_t p, std::size_t i) {
        return "stability correction destroys positive support: species " + std::to_string(i) + " in " + c.problem.phases[p].name;
    };
    Eigen::VectorXd applied(cells.size());
    std::vector<std::vector<double>> inventory(at.amount.size(), std::vector<double>(c.n_feed.size(), 0.0));
    std::optional<std::size_t> bounded;
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        inventory[p][i] = amount(j) + intended(j);
        // #7: a trace cell whose linear step is not positive keeps n exp(dn/n); the applied balance decides.
        if (!(inventory[p][i] > 0.0)) { inventory[p][i] = amount(j) * std::exp(intended(j) / amount(j)); if (!bounded) bounded = j; }
        if (!(inventory[p][i] > 0.0) || !std::isfinite(inventory[p][i])) refuse(lost(p, i), std::abs(inventory[p][i]));
        applied(j) = inventory[p][i] - amount(j);
        candidate.feasibility_estimate += (std::abs(applied(j)) + std::abs(intended(j) - applied(j))) * std::abs(mu(j));
        candidate.correction_norm += std::abs(applied(j));
    }
    if (const double defect = (J * applied - rhs).lpNorm<Eigen::Infinity>(); defect > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
        refuse(bounded ? lost(cells[*bounded].first, cells[*bounded].second) : "stability amount correction applied balance defect", defect);
    at.energy_actions.clear();
    report.candidate_coordinates.clear();
    for (std::size_t p = 0; p < at.amount.size(); ++p) {
        if (at.amount[p] > 0.0) {
            at.amount[p] = std::accumulate(inventory[p].begin(), inventory[p].end(), 0.0);
            for (auto i : c.problem.phases[p].support) at.composition[p][i] = inventory[p][i] / at.amount[p];
        }
        auto q = phase_energy(*c.model, at.temperature, at.pressure, at.volume[p], at.composition[p],
            report.standard_over_rt, report.work, &at.derivative[p]);
        double normalization = -1.0, charge = 0.0, potential = 0.0, residual = 0.0, log_normalization = -1.0, log_charge = 0.0;
        for (auto i : c.problem.phases[p].support) {
            normalization += at.composition[p][i]; charge += z[i] * at.composition[p][i];
            potential += at.composition[p][i] * report.plane_over_rt[i];
            residual += at.composition[p][i] * (std::log(at.composition[p][i] / at.volume[p])
                + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i] - report.plane_over_rt[i]);
            report.candidate_coordinates.push_back(std::log(at.composition[p][i]));
            const double logged = std::exp(report.candidate_coordinates.back());
            log_normalization += logged; log_charge += z[i] * logged;
        }
        if (std::max({std::abs(normalization), std::abs(charge), std::abs(log_normalization), std::abs(log_charge)}) > policy::COMPOSITION_SUM_TOLERANCE.value)
            throw std::invalid_argument("stability candidate arithmetic feasibility failure");
        report.candidate_coordinates.push_back(at.amount[p]);
        report.candidate_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature,
            at.composition[p]) / at.volume[p]));
        q["tpd_over_rt"] = q.at("g_over_rt") - potential;
        q["chemical_plus_electric_contraction"] = residual;
        q["present_identity_defect_over_rt"] = q["tpd_over_rt"] - present_phase_residual_distance(residual, q.at("pressure_row"));
        q["normalization_residual"] = normalization; q["charge_residual"] = charge;
        candidate.quantities.push_back(std::move(q));
    }
    Eigen::VectorXd represented_rhs = Eigen::VectorXd::Zero(J.rows());
    for (std::size_t k = 0; k < c.rank; ++k) represented_rhs(k) = c.balance[k];
    Eigen::VectorXd logged_rhs = represented_rhs;
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        const double represented = at.amount[p] * at.composition[p][i];
        const double logged = at.amount[p] * std::exp(std::log(at.composition[p][i]));
        if (!(represented > 0.0) || !(logged > 0.0)) refuse(lost(p, i), std::abs(inventory[p][i]));
        represented_rhs -= J.col(j) * represented; logged_rhs -= J.col(j) * logged;
        const double displacement = std::abs(represented - inventory[p][i]) + std::abs(logged - represented);
        candidate.correction_norm += displacement;
        candidate.feasibility_estimate += displacement * std::abs(mu(j));
    }
    if (!represented_rhs.allFinite() || !logged_rhs.allFinite()
        || std::max(represented_rhs.lpNorm<Eigen::Infinity>(), logged_rhs.lpNorm<Eigen::Infinity>()) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
        throw std::invalid_argument("stability candidate arithmetic feasibility failure");
    for (auto& q : candidate.quantities) {
        q["correction_l1_mol"] = candidate.correction_norm;
        q["feasibility_estimate_over_rt"] = candidate.feasibility_estimate;
    }
    report.present_phases = candidate.quantities;
    return candidate;
}

double repair_effect(const CompiledSystem &c, const PhasePoint &at, const PhasePoint &corrected,
                     const std::vector<double> &plane, Quantities &diagnostics, const std::string &prefix) {
    const auto &m = *c.model;
    long double allowance = 0;
    size_t ns = c.n_feed.size();
    double T = at.temperature, P = at.pressure, rt = GAS_CONSTANT_J_PER_MOL_K * T;
    std::vector<double> phi(at.amount.size(), 0);
    size_t electro = 0;
    for (auto &row : c.rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = at.point[c.primal_size + c.rank + electro++];
    const long double limit = std::sqrt(std::numeric_limits<double>::epsilon());
    for (size_t p = 0; p < at.amount.size(); ++p) {
        if (!(at.amount[p] > 0)) continue;
        std::vector<long double> n(ns), actual(ns);
        long double N = 0, V = (long double)at.amount[p] * at.volume[p],
                    dV = (long double)corrected.amount[p] * corrected.volume[p] - V, rho = std::abs(dV) / V,
                    plane_movement = 0;
        std::vector<double> dn(ns), x(ns);
        for (size_t i = 0; i < ns; ++i) {
            n[i] = (long double)at.amount[p] * at.composition[p][i];
            actual[i] = (long double)corrected.amount[p] * corrected.composition[p][i] - n[i];
            plane_movement += plane[i] * actual[i];
            dn[i] = actual[i];
            N += n[i];

            if (n[i] > 0)
                rho = std::max(rho, std::abs(actual[i]) / n[i]);
            else if (actual[i] != 0)
                throw std::domain_error("repair changes structural support");
        }
        for (size_t i = 0; i < ns; ++i) x[i] = n[i] / N;
        double v = V / N;
        ReferenceTaylorError e;
        auto action = m.reference_taylor(T, v, x, dn, PressureActionDirection{}, 0, false, &e);
        if (action.status != PropertyActionStatus::Available) throw std::domain_error("repair action unavailable");
        long double projected = action.contraction.at(0) + (P - e.pressure[0]) * dV / rt;
        for (size_t i = 0; i < ns; ++i)
            if (x[i] > 0)
                projected +=
                    (long double)dn[i] * (std::log((long double)x[i]) + at.standard.derivative[0][i] - plane[i]);

        diagnostics[prefix + std::to_string(p) + "_projected_over_rt_mol"] = projected;
        diagnostics[prefix + std::to_string(p) + "_plane_over_rt_mol"] = plane_movement;
        allowance += std::abs(projected);
        allowance += e.fixed[0] + std::abs(dV / rt) * e.pressure_floor[0];

        size_t width = 1 + ns;
        std::vector<long double> delta = actual, scale = n;
        delta.insert(delta.begin(), dV);
        scale.insert(scale.begin(), V);
        using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;
        Matrix first, envelope = Matrix::Zero(width, width), spread = envelope;
        Eigen::Map<const Eigen::Matrix<long double, Eigen::Dynamic, 1>> scaling(scale.data(), width),
            movement(delta.data(), width);
        const auto scale_matrix = scaling.asDiagonal();
        long double norm = 0, change = 0;
        for (long double t : {0.L, .5L, 1.L}) {
            std::vector<double> amount(ns);
            for (size_t i = 0; i < ns; ++i) amount[i] = n[i] + t * actual[i];
            auto d = m.evaluate(T, V + t * dV, amount, DerivativeOrder::IsothermalPhaseHessian);
            size_t dw = d.gradient.size();
            Matrix h(width, width);
            long double total = std::accumulate(amount.begin(), amount.end(), 0.L);
            long double vol = V + t * dV;
            for (size_t i = 0; i < width; ++i)
                for (size_t j = 0; j < width; ++j) {
                    long double a = d.hessian[(1+i)*dw+1+j];
                    if (i == 0 && j == 0)
                        a += total / (vol * vol);
                    else if (i == 0 || j == 0)
                        a -= 1 / vol;
                    else if (i == j && amount[i - 1] > 0)
                        a += 1 / amount[i - 1];
                    h(i, j) = a;
                }
            if (!h.allFinite()) throw std::domain_error("repair Hessian unavailable");
            if (t == 0) first = h;
            const Matrix variation = h - first;
            norm = std::max(norm, (scale_matrix * h * scale_matrix).cwiseAbs().rowwise().sum().maxCoeff());
            change = std::max(change, (scale_matrix * variation * scale_matrix).cwiseAbs().rowwise().sum().maxCoeff());
            envelope = envelope.cwiseMax(h.cwiseAbs());
            spread = spread.cwiseMax(variation.cwiseAbs());
        }
        const long double remainder =
            .5L * (movement.cwiseAbs().transpose() * (envelope + spread) * movement.cwiseAbs())(0, 0);
        long double hchange = norm ? change / norm : (change ? INFINITY : 0);

        if (!(std::isfinite(rho) && rho <= limit && std::isfinite(hchange) && hchange <= limit))
            throw std::domain_error("Gibbs repair outside local small-repair domain");
        allowance += remainder;
        long double charge = 0, rounding = 0;
        auto d = m.evaluate(T, v, x, DerivativeOrder::IsothermalPhaseHessian);
        for (size_t i = 0; i < ns; ++i) {
            if (x[i] > 0)
                rounding += std::abs(actual[i] - dn[i]) *
                            (std::abs(d.gradient[2 + i]) + std::abs(std::log((long double)x[i] / v)) +
                             std::abs(at.standard.derivative[0][i]));
            charge += (long double)corrected.amount[p] * corrected.composition[p][i] * m.parameters().charge_number[i];
        }

        allowance += rounding + std::abs(phi[p] * charge);
        long double magnitude = std::abs(projected) + std::abs(action.contraction[0]) + std::abs(dV * P / rt);
        for (size_t i = 0; i < ns; ++i)
            magnitude += std::abs((long double)dn[i]) * (std::abs(plane[i]) + std::abs(at.standard.derivative[0][i]));
        allowance += (8 * ns + 24) * std::numeric_limits<long double>::epsilon() * magnitude;
    }
    if (!(std::isfinite(allowance) && allowance >= 0)) throw std::domain_error("repair estimate unavailable");
    return allowance;
}

struct CondensedRefusal final : std::invalid_argument { using std::invalid_argument::invalid_argument; };
// Without a unit vacuum limit, a charged trial counts only on its own liquid side.
void condensed_admission(const std::optional<DensityRoot>& root, bool quotient, Quantities* q = nullptr) {
    if (!root || properties::density_side(root->branch, root->molar_density, root->negative_slope_density) != DensityBranch::Liquid)
        throw CondensedRefusal(quotient
            ? "charged trial outside the condensed domain of permittivity without a unit vacuum limit"
            : "charged trial outside the condensed domain of density-independent permittivity");
    if (q && root->negative_slope_density) (*q)["negative_slope_density_mol_m3"] = *root->negative_slope_density;
}

std::vector<double> trial_coordinates(const CompiledSystem& c, const std::vector<double>& w, double eta) {
    auto point = c.generic_point();
    const auto offsets = c.phase_offsets();
    for (std::size_t i = 0; i < c.problem.phases[0].support.size(); ++i)
        point[offsets[0] + i] = std::log(w[c.problem.phases[0].support[i]]);
    point[offsets[0] + c.problem.phases[0].support.size()] = 1.0;
    point[offsets[0] + c.problem.phases[0].support.size() + 1] = std::log(eta);
    return point;
}
} // namespace

// ALG-STABILITY: cache the original solved-state row uncertainty and supported
// repair once.
GibbsDifferenceReference conserved_split_reference(const CompiledSystem &c, const EquilibriumResult &result) {
    GibbsDifferenceReference out;
    try {

        out.physical = stability_candidate(c, result, out.report).physical;
        std::vector<SolvedStateObservable> o;
        for (size_t i = 0; i < c.n_feed.size(); ++i)
            o.push_back({SolvedStateObservableKind::SpeciesStandardMinusPlane, 0, i});
        SolvedStateActionDirection zero;
        zero.feed_amounts_mol.assign(c.n_feed.size(), 0);
        zero.active_parameters.assign(c.model->active_parameters().size(), 0);
        auto actions = solved_state_first_actions(c, result, o, {zero}).front();

        for (auto &a : actions.results) {
            if (!a.numerical_error) throw std::domain_error("parent row uncertainty unavailable");
            out.row_errors.push_back(*a.numerical_error);
        }

        out.repair_error = repair_effect(c, c.prepare(result.coordinates), out.physical, out.report.plane_over_rt,
                                         out.repair_quantities, "parent_repair_");
    } catch (const std::exception &error) {
        out.cause = error.what();
    }
    return out;
}
GibbsDifference conserved_split_difference(const CompiledSystem &c, const GibbsDifferenceReference &parent,
                                           const CompiledSystem &split_c, const EquilibriumResult &result,
                                           std::size_t donor, const StabilityOptions &options) {
    GibbsDifference out;
    try {
        if (!parent.cause.empty()) {
            out.cause = parent.cause;
            return out;
        }
        StabilityReport split_report;
        const auto split = stability_candidate(split_c, result, split_report);
        out.coordinates = split_report.candidate_coordinates;
        out.quantities = parent.repair_quantities;
        const double repair =
            parent.repair_error + repair_effect(split_c, split_c.prepare(result.coordinates), split.physical,
                                                parent.report.plane_over_rt, out.quantities, "split_repair_");
        // Affine donor path in extensive species amounts and volume; exact native
        // gradients.
        const auto &m = *c.model;
        const auto &r = parent.physical;
        const auto &er = split.physical;
        const auto &s = parent.report;
        const size_t ns = c.n_feed.size();
        double rt = GAS_CONSTANT_J_PER_MOL_K * r.temperature;
        std::vector<long double> n0(ns), dn(ns);
        long double v0 = r.amount[donor] * r.volume[donor], dv = er.amount[0] * er.volume[0] - v0;
        long double plane_defect = 0, trial_plane = 0, parent_error = 0, standard_error = 0, plane_mag = 0;
        for (size_t i = 0; i < ns; ++i) {
            n0[i] = (long double)r.amount[donor] * r.composition[donor][i];
            dn[i] = (long double)er.amount[0] * er.composition[0][i] - n0[i];
            trial_plane += (long double)er.composition.back()[i] * s.plane_over_rt[i];
            plane_defect += (dn[i] + (long double)er.amount.back() * er.composition.back()[i]) * s.plane_over_rt[i];
            const long double weight =
                std::abs(dn[i]) + (long double)er.amount.back() * std::abs(er.composition.back()[i]);
            parent_error += weight * parent.row_errors.at(i);
            if (!options.standard_uncertainty_over_rt.empty())
                standard_error += weight * options.standard_uncertainty_over_rt[i];
            plane_mag += weight * std::abs(s.plane_over_rt[i]);
        }
        long double last = NAN, last_floor = NAN, last_quadrature = NAN;
        bool previous_refined = false;
        // Simpson ladder: first eligible rung from 32 panels with two consecutive refinements, capped at 128.
        for (int steps = 4; steps <= 128; steps *= 2) {
            long double integral = 0, owner_error = 0, arithmetic = 0, volumeerror = 0, charterror = 0;
            for (int j = 0; j <= steps; ++j) {
                long double t = (long double)j / steps, N = 0, V = v0 + t * dv;
                std::vector<double> x(ns), nu(ns);
                for (size_t i = 0; i < ns; ++i) N += n0[i] + t * dn[i];
                for (size_t i = 0; i < ns; ++i) {
                    x[i] = (n0[i] + t * dn[i]) / N;
                    nu[i] = dn[i];
                }
                if (!(V > 0 && N > 0 && std::isfinite(V) && std::isfinite(N)))
                    throw std::domain_error("nonpositive Gibbs path");
                for (size_t i = 0; i < ns; ++i)
                    if ((n0[i] > 0 && !(x[i] > 0)) || (n0[i] == 0 && dn[i] != 0))
                        throw std::domain_error("Gibbs path changes positive support");
                double v = V / N;
                ReferenceTaylorError error;
                auto action = m.reference_taylor(r.temperature, v, x, nu, PressureActionDirection{}, 0, false, &error);
                if (action.contraction.empty() || action.status != PropertyActionStatus::Available)
                    throw std::domain_error("reference contraction unavailable");
                if (!(error.fixed.at(0) >= 0 && error.pressure_floor.at(0) >= 0))
                    throw std::domain_error("directional Gibbs error unavailable");
                long double val = action.contraction[0] + (r.pressure - error.pressure[0]) * dv / rt,
                            mag =
                                std::abs(action.contraction[0]) + std::abs((r.pressure - error.pressure[0]) * dv / rt),
                            logerror = 0;
                for (size_t i = 0; i < ns; ++i)
                    if (x[i] > 0) {
                        long double term = std::log((long double)x[i]) + s.standard_over_rt[i] - s.plane_over_rt[i];
                        val += dn[i] * term;
                        mag += std::abs(dn[i]) * (std::abs(std::log((long double)x[i])) +
                                                  std::abs(s.standard_over_rt[i]) + std::abs(s.plane_over_rt[i]));
                        logerror += std::abs(dn[i]) * (std::abs(std::log((long double)x[i]) - std::log(x[i])) +
                                                       std::numeric_limits<long double>::epsilon() *
                                                           std::abs(std::log((long double)x[i])));
                    }
                auto d = m.evaluate(r.temperature, v, x, DerivativeOrder::IsothermalPhaseHessian);
                long double chart = 0;
                size_t width = 2 + ns;
                for (size_t col = 1; col < width; ++col) {
                    long double movement = col == 1 ? (long double)v - V / N
                                                    : (long double)x[col - 2] - (n0[col - 2] + t * dn[col - 2]) / N;
                    long double deriv = dv * (d.hessian[width + col] + (col == 1 ? 1 / (v * v) : -1 / v));
                    for (size_t i = 0; i < ns; ++i)
                        if (x[i] > 0)
                            deriv += dn[i] * (d.hessian[(2 + i) * width + col] + (col == 1       ? -1 / v
                                                                                  : col == 2 + i ? 1 / x[i]
                                                                                                 : 0));
                    chart += std::abs(movement * deriv);
                }
                for (size_t i = 0; i < ns; ++i)
                    if (x[i] > 0) chart += std::abs(dn[i] - nu[i]) * std::abs(d.gradient[2 + i] - std::log(v));
                long double gamma = (8 * ns + 24) * std::numeric_limits<long double>::epsilon();
                int weight = (j == 0 || j == steps) ? 1 : (j % 2 ? 4 : 2);
                integral += weight * val;
                const long double simpson_weight = (long double)weight / (3 * steps);
                charterror += simpson_weight * chart;
                owner_error += simpson_weight * error.fixed[0];
                volumeerror += simpson_weight * std::abs(dv / rt) * error.pressure_floor[0];
                arithmetic += simpson_weight * (logerror + gamma / (1 - gamma) * mag);
            }
            long double trial = (long double)er.amount.back() * (split.quantities.back().at("g_over_rt") - trial_plane);
            long double current_integral = integral / (3 * steps), delta = current_integral + trial + plane_defect;
            double trialerror = er.amount.back() * energy_estimate(split.quantities.back());
            long double quadrature = std::abs(current_integral - last) / 15;
            long double round = (8 * ns + 2 * steps + 40) * std::numeric_limits<long double>::epsilon() *
                                (plane_mag + std::abs(current_integral) + std::abs(trial));

            long double integration_floor = owner_error + volumeerror + arithmetic + round + charterror;
            bool refined = std::isfinite(quadrature) && quadrature <= integration_floor + last_floor;

            long double complete = integration_floor + trialerror + quadrature + std::abs(plane_defect) + parent_error +
                                   standard_error + repair;
            const bool qualified = previous_refined && refined;
            if (steps >= 32 && (qualified || steps == 128)) {
                out.quantities.merge(Quantities{{"gibbs_change_over_rt_mol", (double)delta},
                                                {"quadrature_panels", (double)steps},
                                                {"previous_quadrature_error_over_rt_mol", (double)last_quadrature},
                                                {"integration_floor_over_rt_mol", (double)integration_floor},
                                                {"previous_integration_floor_over_rt_mol", (double)last_floor},
                                                {"complete_error_over_rt_mol", (double)complete},
                                                {"gradient_error_over_rt_mol", (double)(owner_error)},
                                                {"pressure_error_over_rt_mol", (double)(volumeerror)},
                                                {"parent_error_over_rt_mol", (double)parent_error},
                                                {"standard_error_over_rt_mol", (double)standard_error},
                                                {"trial_error_over_rt_mol", trialerror},
                                                {"repair_error_over_rt_mol", repair},
                                                {"inventory_plane_over_rt_mol", (double)plane_defect},
                                                {"quadrature_error_over_rt_mol", (double)quadrature},
                                                {"ideal_error_over_rt_mol", (double)(arithmetic)},
                                                {"summation_error_over_rt_mol", (double)round},
                                                {"chart_error_over_rt_mol", (double)charterror}});
                if (!qualified) {
                    out.cause = "incomplete Gibbs integral refinement";
                    return out;
                }
                for (const auto &[name, value] : out.quantities)
                    if (!std::isfinite(value) || (name.find("error_") != std::string::npos && value < 0)) {
                        out.cause = "unavailable Gibbs difference estimate: " + name;
                        return out;
                    }

                out.status = delta + complete < 0   ? GibbsDifferenceStatus::Decrease
                             : delta - complete > 0 ? GibbsDifferenceStatus::NoDecrease
                                                    : GibbsDifferenceStatus::Overlap;
                return out;
            }
            last = current_integral;
            last_floor = integration_floor;
            last_quadrature = quadrature;
            previous_refined = refined;
        }
    } catch (const std::exception &error) {
        out.cause = error.what();
    }
    return out;
}

// EqID michelsen_modified_tangent_plane: J_ij = F_ij + 1 + p_i p_j / p_V, p_i = 1/v - F_iV, p_V = -F_VV - 1/v^2.
std::vector<double> tangent_plane_jacobian(const HelmholtzDerivatives& d, double v, const std::vector<std::size_t>& support) {
    const std::size_t width = d.gradient.size(), n = support.size();
    const double pv = -d.hessian[width + 1] - 1.0 / (v * v);
    std::vector<double> p(n), J(n * n);
    for (std::size_t a = 0; a < n; ++a) p[a] = 1.0 / v - d.hessian[(2 + support[a]) * width + 1];
    for (std::size_t a = 0; a < n; ++a) for (std::size_t b = 0; b < n; ++b)
        J[a * n + b] = d.hessian[(2 + support[a]) * width + 2 + support[b]] + 1.0 + p[a] * p[b] / pv;
    return J;
}

// Charge repair increases the opposite-sign pivot before common normalization (ties by index).
std::vector<double> reconstruct_trial_composition(std::vector<double> w,
    const std::vector<std::size_t>& support, const std::vector<int>& z, PhaseKind kind) {
    if (support.empty() || z.size() != w.size()) throw std::invalid_argument("invalid reconstruction support");
    double positive = 0.0, negative = 0.0;
    std::optional<std::size_t> plus, minus;
    for (auto i : support) {
        if (i >= w.size() || !(std::isfinite(w[i]) && w[i] > 0.0))
            throw std::invalid_argument("nonpositive/nonfinite trial composition");
        if (z[i] > 0.0) { positive += z[i] * w[i]; if (!plus || w[i] > w[*plus]) plus = i; }
        if (z[i] < 0.0) { negative -= z[i] * w[i]; if (!minus || w[i] > w[*minus]) minus = i; }
    }
    if (plus.has_value() != minus.has_value()) throw std::invalid_argument("one-sign trial support");
    if (plus) {
        if (kind == PhaseKind::Vapor) {
            const double shift = 0.5 * (std::log(negative) - std::log(positive));
            for (auto i : support) if (z[i]) w[i] = std::exp(std::log(w[i]) + (z[i] > 0 ? shift : -shift));
        } else {
            const double ratio = std::sqrt(negative / positive);
            for (auto i : support) { if (z[i] > 0) w[i] *= ratio; if (z[i] < 0) w[i] /= ratio; }
        }
    }
    double total = 0.0;
    for (auto i : support) total += w[i];
    if (!(std::isfinite(total) && total > 0.0)) throw std::invalid_argument("trial normalization failed");
    for (auto i : support) w[i] /= total;
    if (plus && kind != PhaseKind::Vapor) {
        double charge = 0.0;
        for (auto i : support) charge += z[i] * w[i];
        if (charge > 0.0) w[*minus] += charge / (-z[*minus]);
        if (charge < 0.0) w[*plus] += (-charge) / z[*plus];
        total = 0.0;
        for (auto i : support) total += w[i];
        if (!(std::isfinite(total) && total > 0.0)) throw std::invalid_argument("trial normalization failed");
        for (auto i : support) w[i] /= total;
    } else if (!plus) {
        const auto pivot = *std::max_element(support.begin(), support.end(), [&](auto a, auto b) { return w[a] < w[b]; });
        double sum = 0.0; for (auto i : support) sum += w[i]; w[pivot] += 1.0 - sum;
    }
    double sum = 0.0, charge = 0.0;
    for (auto i : support) {
        if (!(std::isfinite(w[i]) && w[i] > 0.0)) throw std::invalid_argument("reconstruction destroys support");
        sum += w[i]; charge += z[i] * w[i];
    }
    if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value)
        throw std::invalid_argument("trial reconstruction arithmetic residual");
    return w;
}
namespace {
// #21: the second pass spends one evaluation credit before each density or finite-split energy call.
void spend_credit(const DriverLimits& limits, Quantities& work) {
    if (!limits.spend(limits.evaluations, limits.evaluation_cap)) limits.check();
    ++work["recovery_density_credits"];
}
// For a relaxed trial the frozen slope can classify no descent only within plane arithmetic error.
void donor_descent(const CompiledSystem& c, const Candidate& candidate, StabilityTrial& trial,
    StabilityReport& report, const StabilityOptions& options, const DriverLimits* limits, bool relaxed, bool credited = false) {
    const auto& at = candidate.physical;
    const auto& w = trial.composition;
    auto& q = trial.quantities;
    double capacity = 0.0;
    for (std::size_t p = 0; p < at.amount.size(); ++p) {
        if (!(at.amount[p] > 0.0)) continue;
        double beta = std::numeric_limits<double>::infinity();
        for (auto i : trial.phase.support) {
            if (at.composition[p][i] <= 0.0) { beta = 0.0; break; }
            beta = std::min(beta, at.amount[p] * at.composition[p][i] / w[i]);
        }
        if (std::isfinite(beta) && beta > capacity) { capacity = beta; trial.donor = p; }
    }
    if (!trial.donor) {
        trial.evidence = options.evaluate_supplied_trials ? "evaluated_without_positive_donor" : "negative_plane_without_positive_donor";
        if (options.evaluate_supplied_trials) trial.message = "no_positive_full_support_donor";
        return;
    }
    const auto p = *trial.donor;
    const auto& donor = candidate.quantities[p];
    double chemical = 0.0, magnitude = std::abs(q.at("g_over_rt"));
    double standard_error = 0.0;
    for (auto i : trial.phase.support) {
        const double term = w[i] * (std::log(at.composition[p][i] / at.volume[p])
            + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i]);
        chemical += term; magnitude += std::abs(term);
        if (!options.standard_uncertainty_over_rt.empty())
            standard_error += 2.0 * w[i] * options.standard_uncertainty_over_rt[i];
    }
    const double slope = q.at("g_over_rt") - chemical + donor.at("pressure_row");
    const double gamma = (w.size() + 3.0) * std::numeric_limits<double>::epsilon();
    const double arithmetic = gamma * (magnitude + std::abs(donor.at("pressure_row"))) / (1.0 - gamma);
    const double uncertainty = energy_estimate(q) + (1.0 + std::abs(at.volume[p]))
        * at.derivative[p].gradient_floor + candidate.feasibility_estimate + arithmetic + standard_error;
    q["donor_capacity_mol"] = capacity;
    q["donor_chemical_contraction_over_rt"] = chemical;
    q["donor_pressure_row"] = donor.at("pressure_row");
    q["donor_correction_over_rt"] = slope - q.at("tpd_over_rt");
    q["donor_slope_over_rt"] = slope;
    q["descent_error_estimate_over_rt"] = uncertainty;
    q["standard_error_estimate_over_rt"] = standard_error;
    q["feasibility_error_estimate_over_rt"] = candidate.feasibility_estimate;
    q["arithmetic_error_estimate_over_rt"] = arithmetic;
    if (options.evaluate_supplied_trials) { trial.evidence = "evaluated"; return; }
    trial.evidence = "indeterminate";
    if (std::isfinite(slope) && std::isfinite(uncertainty) && uncertainty >= 0.0 && slope >= -uncertainty) {
        if (relaxed) {
            double plane_magnitude = std::abs(q.at("g_over_rt"));
            for (auto i : trial.phase.support) plane_magnitude += std::abs(w[i] * report.plane_over_rt[i]);
            const double roundoff = (2.0 * trial.phase.support.size() + 1.0) * std::numeric_limits<double>::epsilon();
            const double error = q.at("summation_estimate_over_rt") + roundoff * plane_magnitude / (1.0 - roundoff);
            if (!(std::isfinite(error) && error >= 0.0 && q.at("tpd_over_rt") >= -error)) return;
            q["tpd_arithmetic_error_estimate_over_rt"] = error;
        }
        trial.evidence = "no_resolved_descent";
        return;
    }
    if (relaxed || !descent_resolved(q.at("tpd_over_rt"), slope, uncertainty)) return;
    const auto& z = c.model->parameters().charge_number;
    const double limit = policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total);
    // The donor after it gives `moved` mol of composition `y`: its amount and composition, with the split's conservation checks.
    const auto give = [&](double moved, const std::vector<double>& y, Quantities& split) {
        const double remaining = at.amount[p] - moved;
        std::vector<double> x(w.size(), 0.0), error(w.size(), 0.0);
        double sum = 0.0, charge = 0.0, inventory_error = 0.0;
        for (auto i : c.problem.phases[p].support) {
            const double n = at.amount[p] * at.composition[p][i] - moved * y[i];
            if (!(n > 0.0 && remaining > 0.0)) throw std::invalid_argument("split destroys donor support");
            x[i] = n / remaining; sum += x[i]; charge += z[i] * x[i];
            error[i] = remaining * x[i] + moved * y[i] - at.amount[p] * at.composition[p][i];
            inventory_error = std::max(inventory_error, std::abs(error[i]));
        }
        const auto moiety_error = matvec(c.B, error);
        split["component_error_mol"] = inventory_error;
        split["moiety_error_mol"] = maximum_abs(moiety_error);
        split["normalization_residual"] = sum - 1.0; split["charge_residual"] = charge;
        if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value
            || inventory_error > limit || maximum_abs(moiety_error) > limit)
            throw std::invalid_argument("finite split arithmetic conservation failed");
        return std::pair(remaining, std::move(x));
    };
    // The split's optimizer coordinates appended to `out`: every parent phase with the donor as given, then the added phase.
    const auto write = [&](std::vector<double>& out, double remaining, const std::vector<double>& x, double donor_volume,
                           const std::vector<std::size_t>& support, double moved, const std::vector<double>& y) {
        for (std::size_t phase = 0; phase < at.amount.size(); ++phase) {
            const auto& composition = phase == p ? x : at.composition[phase];
            for (auto i : c.problem.phases[phase].support) out.push_back(std::log(composition[i]));
            out.push_back(phase == p ? remaining : at.amount[phase]);
            out.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature,
                composition) / (phase == p ? donor_volume : at.volume[phase])));
        }
        for (auto i : support) out.push_back(std::log(y[i]));
        out.push_back(moved);
        out.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature, y) / q.at("volume_m3_per_mol")));
    };
    double beta = capacity / 2.0;
    for (std::size_t check = 0; check <= 8; ++check, beta /= 2.0) {
        Quantities split{{"beta_mol", beta}, {"admitted", 0.0}};
        if (limits) limits->check();
        try {
            const auto [remaining, x] = give(beta, w, split);
            ++report.work["finite_split_evaluations"];
            if (credited) spend_credit(*limits, report.work);
            const double donor_volume = 1.0 / declared_root(*c.model, at.temperature, at.pressure, x, c.problem.phases[p].kind).molar_density;
            if (credited) spend_credit(*limits, report.work);
            const auto changed = phase_energy(*c.model, at.temperature, at.pressure, donor_volume, x,
                report.standard_over_rt, report.work);
            const double delta = remaining * changed.at("g_over_rt") + beta * q.at("g_over_rt")
                - at.amount[p] * donor.at("g_over_rt");
            const double error_estimate = remaining * energy_estimate(changed) + beta * energy_estimate(q)
                + at.amount[p] * energy_estimate(donor) + candidate.feasibility_estimate
                + gamma * (std::abs(remaining * changed.at("g_over_rt")) + std::abs(beta * q.at("g_over_rt"))
                    + std::abs(at.amount[p] * donor.at("g_over_rt"))) / (1.0 - gamma) + beta * standard_error;
            split["delta_g_over_rt"] = delta; split["error_estimate_over_rt"] = error_estimate;
            split["admitted"] = 1.0;
            if (descent_resolved(q.at("tpd_over_rt"), delta, error_estimate)) {
                trial.resolved_descent = true;
                trial.evidence = "numerically_resolved_descending_split_with_estimated_uncertainty";
                q["split_beta_mol"] = beta; q["split_delta_g_over_rt"] = delta;
                q["split_error_estimate_over_rt"] = error_estimate;
                write(trial.split_coordinates, remaining, x, donor_volume, trial.phase.support, beta, w);
                // ALG-STABILITY: the same split with the neutral `bound_species` left in the donor: transfers beta w_i on
                // the retained species, none on the omitted ones; same volumes, no density or energy evaluation.
                if (!trial.optimizer.bound_species.empty()) try {
                    auto y = w;
                    auto support = trial.phase.support;
                    for (auto j : trial.optimizer.bound_species) {
                        if (z[j] != 0.0) throw std::invalid_argument("an ion ends at its composition bound");
                        y[j] = 0.0;
                        support.erase(std::find(support.begin(), support.end(), j));
                    }
                    for (auto j : support) y[j] *= beta;
                    const double moved = std::accumulate(y.begin(), y.end(), 0.0);
                    if (!(moved > 0.0)) throw std::invalid_argument("reduced split transfers no inventory");
                    for (auto j : support) y[j] /= moved;
                    Quantities unused;
                    const auto [rest, kept] = give(moved, y, unused);
                    std::vector<double> seed, defect(w.size(), 0.0);
                    write(seed, rest, kept, donor_volume, support, moved, y);
                    // Check the inventory actually represented by the stored log fractions.
                    std::size_t column = 0;
                    for (std::size_t phase = 0; phase < p; ++phase) column += c.problem.phases[phase].support.size() + 2;
                    for (auto i : c.problem.phases[p].support)
                        defect[i] = rest * std::exp(seed[column++]);
                    column = seed.size() - support.size() - 2;
                    for (auto i : support) defect[i] += moved * std::exp(seed[column++]);
                    for (auto i : c.problem.phases[p].support) defect[i] -= at.amount[p] * at.composition[p][i];
                    if (!std::isfinite(maximum_abs(seed)) || !(maximum_abs(defect) <= limit && maximum_abs(matvec(c.B, defect)) <= limit))
                        throw std::invalid_argument("represented amounts lose inventory");
                    trial.reduced_split_coordinates = std::move(seed);
                } catch (const std::exception& error) { trial.message += "; reduced split: " + std::string(error.what()); }
                trial.split_checks.push_back(std::move(split));
                return;
            }
        } catch (const std::exception& error) { trial.message += "; finite split: " + std::string(error.what()); }
        trial.split_checks.push_back(std::move(split));
    }
}

// ALG-STABILITY: the donor and the trial as one conserved fixed-T/P Gibbs cell of the original reactions and
// standards, one cold shared-driver solve under one allowance, then the supplied-split Gibbs difference.
GibbsDifference relaxed_split(const CompiledSystem& c, const PhasePoint& at, const StabilityTrial& trial, std::size_t p,
    const GibbsDifferenceReference& parent, const StabilityOptions& options, Quantities& work) {
    GibbsDifference out;
    StageWork stages;
    try {
        EquilibriumProblem local;
        local.temperature = FreeValue::fixed(at.temperature);
        local.pressure = FreeValue::fixed(at.pressure);
        local.reactions = c.problem.reactions;
        local.standard_gibbs_j_per_mol = c.problem.standard_gibbs_j_per_mol;
        local.source_reference = c.problem.source_reference;
        local.solver_options = c.problem.solver_options;
        local.phases = {c.problem.phases[p], trial.phase};
        local.phases[1].kind = PhaseKind::Vapor;
        local.phases[1].name = "reaction-relaxed-trial";
        const std::vector<double>* guess[] = {&at.composition[p], &trial.composition};
        const double volume[] = {at.volume[p], trial.quantities.at("volume_m3_per_mol")};
        for (std::size_t j = 0; j < 2; ++j) {
            local.phases[j].amount = AmountSpec{};
            local.phases[j].composition_guess.clear();
            for (const auto i : local.phases[j].support)  // compilation reads guesses in support order
                local.phases[j].composition_guess.push_back((*guess[j])[i]);
            local.phases[j].packing_guess = c.model->packing_fraction_coefficient(at.temperature, *guess[j]) / volume[j];
        }
        for (const double x : at.composition[p]) local.feed.values.push_back(at.amount[p] * x);
        auto compiled = compile_problem(*c.model, local);
        compiled.reaction_relaxed_split = true;
        const auto limits = DriverLimits::recovery_allowance(std::chrono::steady_clock::time_point::max(),
            static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value));
        const auto start = [&] {
            try { return boundary_start(compiled, true, &limits, &stages); }
            catch (...) {  // a refused initialization keeps its sub-solve work, with no driver invocation
                work["iterations"] += stages.iterations;
                work["evaluations"] += stages.evaluations;
                work["hessian_calls"] += stages.hessian_calls;
                throw;
            }
        }();
        core::Scope scope("reaction_relaxed_split");
        auto split = run_driver(compiled, start, structural_result(compiled, scope), scope, nullptr,
            &limits, Attempt::Cold, {}, false, false, Retry::None);
        std::move(stages).report(split);
        count_driver(work, split, "split_driver_invocations");
        if (!split.success) throw std::runtime_error(split.message);
        limits.check();
        out = conserved_split_difference(c, parent, compiled, split, p, options);
    } catch (const std::exception& error) { out.cause = error.what(); }
    return out;
}

// The one mapping of a split comparison to trial evidence: decrease, no decrease, or unavailable with its cause.
void classify_split(const CompiledSystem& c, const StabilityReport& report, StabilityTrial& trial, const GibbsDifference& d) {
    auto check = d.quantities;
    check["admitted"] = d.status != GibbsDifferenceStatus::Unavailable;
    trial.evidence = d.status == GibbsDifferenceStatus::Decrease ? "numerically_resolved_descending_split_with_estimated_uncertainty"
        : d.status == GibbsDifferenceStatus::NoDecrease ? "no_resolved_descent" : "indeterminate";
    trial.resolved_descent = d.status == GibbsDifferenceStatus::Decrease;
    if (d.status == GibbsDifferenceStatus::Overlap) trial.message += "; reaction-relaxed split: split Gibbs difference overlaps zero";
    if (d.status == GibbsDifferenceStatus::Unavailable) trial.message += "; reaction-relaxed split: " + d.cause;
    else {
        // Parent phases in their order with the repaired local donor substituted, then the repaired trial.
        std::size_t offset = 0, width = c.problem.phases[*trial.donor].support.size() + 2;
        for (std::size_t p = 0; p < *trial.donor; ++p) offset += c.problem.phases[p].support.size() + 2;
        trial.split_coordinates = report.candidate_coordinates;
        std::copy(d.coordinates.begin(), d.coordinates.begin() + width, trial.split_coordinates.begin() + offset);
        trial.split_coordinates.insert(trial.split_coordinates.end(), d.coordinates.begin() + width, d.coordinates.end());
        auto& q = trial.quantities;
        q["split_beta_mol"] = check["beta_mol"] = d.coordinates.end()[-2];
        q["split_delta_g_over_rt"] = d.quantities.at("gibbs_change_over_rt_mol");
        q["split_error_estimate_over_rt"] = d.quantities.at("complete_error_over_rt_mol");
    }
    trial.split_checks.push_back(std::move(check));
}

enum class SuppliedStage { Seed, Inventory, Coordinates, Raw, Composition, State, Energy, Complete };
struct SuppliedPoint {
    SuppliedStage next = SuppliedStage::Seed;
    std::string failure;
    std::vector<double> seed, coordinates, composition;
    PhasePoint raw;
    double volume = 0.0;
    State state;
    PhaseEnergyPoint energy;
};
void advance_supplied(SuppliedPoint& point, SuppliedStage through, const HelmholtzModel& model,
    double T, const PhaseSpec& phase, CompiledSystem* c, double eta, Quantities& work) {
    const auto& z = model.parameters().charge_number;
    while (point.next <= through) {
        if (!point.failure.empty()) throw std::invalid_argument(point.failure);
        try {
            switch (point.next) {
            case SuppliedStage::Seed:
                point.seed = phase.composition_guess;
                if (point.seed.size() != z.size()) throw std::invalid_argument("supplied trial composition size");
                for (std::size_t i = 0; i < point.seed.size(); ++i)
                    if (std::find(phase.support.begin(), phase.support.end(), i) == phase.support.end() && point.seed[i] != 0.0)
                        throw std::invalid_argument("supplied trial composition outside support");
                point.seed = reconstruct_trial_composition(std::move(point.seed), phase.support, z, phase.kind);
                ++work["reconstruction_evaluations"]; break;
            case SuppliedStage::Inventory:
                for (const auto& bound : inventory_bounds(*c, phase.support, false, 0.0)) admit_inventory_bound(bound);
                break;
            case SuppliedStage::Coordinates:
                point.coordinates = trial_coordinates(*c, point.seed, eta);
                break;
            case SuppliedStage::Raw:
                point.raw = c->prepare(point.coordinates, true);
                for (auto& d : point.raw.derivative) std::vector<double>().swap(d.hessian);
                break;
            case SuppliedStage::Composition:
                point.composition = reconstruct_trial_composition(point.raw.composition[0], phase.support, z, phase.kind);
                ++work["reconstruction_evaluations"];
                point.volume = model.packing_fraction_coefficient(T, point.composition) / point.raw.eta[0]; break;
            case SuppliedStage::State:
                if (phase.kind == PhaseKind::Vapor && std::any_of(phase.support.begin(), phase.support.end(), [&](auto i) { return z[i] != 0; }))
                    point.state = model.state(T, 1.0 / point.volume, point.composition);
                break;
            case SuppliedStage::Energy:
                point.energy = phase_energy_point(model, T, point.volume, point.composition, work);
                std::vector<double>().swap(point.energy.derivative.hessian); break;
            case SuppliedStage::Complete: break;
            }
            point.next = static_cast<SuppliedStage>(static_cast<int>(point.next) + 1);
        } catch (const std::exception& error) { point.failure = error.what(); throw; }
    }
}
using TrialSystems = std::map<std::pair<PhaseKind, std::vector<std::size_t>>, std::pair<CompiledSystem, std::vector<std::vector<double>>>>;
CompiledSystem& supplied_system(TrialSystems& systems, const HelmholtzModel& model, double T, double P,
    const PhaseSpec& phase, SuppliedPoint& point, const std::vector<double>& shift,
    std::vector<std::vector<double>>& transformation, Quantities& work) {
    auto key = std::make_pair(phase.kind, phase.support);
    auto found = systems.find(key);
    if (found == systems.end()) {
        auto c = compile_stability_trial(model, T, P, phase, point.seed, shift, transformation);
        found = systems.emplace(key, std::make_pair(std::move(c), transformation)).first;
        if (point.next == SuppliedStage::Inventory) point.next = SuppliedStage::Coordinates;
    } else {
        transformation = found->second.second;
        auto& c = found->second.first;
        c.problem.feed.values = point.seed;
        c.n_feed = resolve_feed(model, c.problem.feed);
        admit_trial_feed(c); normalize_trial_feed(c);
        advance_supplied(point, SuppliedStage::Inventory, model, T, phase, &c, 0.0, work);
    }
    return found->second.first;
}
} // namespace
struct PreparedSuppliedPoints { std::vector<SuppliedPoint> entries; };
std::shared_ptr<const PreparedSuppliedTrials> prepare_supplied_trials(const HelmholtzModel& model, double T,
    const std::vector<PhaseSpec>& trials) {
    core::Scope scope("prepare_supplied_trials");
    auto set = std::make_shared<PreparedSuppliedTrials>();
    set->model = &model; set->temperature = T; set->trials = trials;
    auto points = std::make_shared<PreparedSuppliedPoints>(); set->points = points;
    TrialSystems systems; core::CalculationWork work;
    for (auto& phase : set->trials) {
        if (phase.support.empty()) throw std::invalid_argument("prepared trial requires declared support");
        phase.support = resolve_support(phase, phase.support, model.parameters().charge_number, &model.parameters());
        auto& point = points->entries.emplace_back();
        try {
            advance_supplied(point, SuppliedStage::Seed, model, T, phase, nullptr, 0.0, set->work);
            if (!phase.packing_guess) continue;
            packing_admission(*phase.packing_guess);
            std::vector<std::vector<double>> transformation;
            CompiledSystem* c;
            try { c = &supplied_system(systems, model, T, 1e5, phase, point,
                std::vector<double>(model.parameters().component_count(), 0.0), transformation, set->work); }
            catch (const std::exception&) { if (!point.failure.empty()) throw; point.next = SuppliedStage::Inventory; continue; }
            advance_supplied(point, SuppliedStage::Energy, model, T, phase, c, *phase.packing_guess, set->work);
        } catch (const std::exception&) {}
    }
    report_work(set->work, work); return set;
}
namespace {
StabilityReport test_snapshot(const CompiledSystem& input, const EquilibriumResult& result,
    const StabilityOptions& requested, const std::vector<std::pair<PhaseSpec, std::vector<double>>>& removed,
    std::vector<std::string>* unsearched = nullptr, const DriverLimits* limits = nullptr, bool diagnostics = false,
    std::optional<GibbsDifferenceReference>* reference = nullptr, DriverLimits* recovery = nullptr) {
    auto options = requested;
    if (options.prepared_trials) {
        const auto& set = *options.prepared_trials;
        if (options.trials) throw std::invalid_argument("prepared trials given together with trials");
        if (set.model != input.model) throw std::invalid_argument("prepared trials model mismatch");
        if (!result.temperature || std::memcmp(&set.temperature, &*result.temperature, sizeof(double)) != 0)
            throw std::invalid_argument("prepared trials temperature mismatch");
        options.trials = set.trials; options.evaluate_supplied_trials = true;
    }
    TrialSystems systems;
    std::size_t entry = 0;
    if (diagnostics && limits) throw std::invalid_argument("stability diagnostics require an unlimited explicit test");
    auto c = input; c.phase_calculation = c.problem.phase_calculation = PhaseCalculation::Discovery;
    StabilityReport report;
    const bool fast = options.mode == DiscoveryMode::Fast;
    bool stopped = false;
    core::CalculationWork work;
    if (options.evaluate_supplied_trials && !options.trials)
        throw std::invalid_argument("supplied trial evaluation requires explicit trials");
    const auto supports = stability_supports(c, options);
    const auto& z = c.model->parameters().charge_number;
    if (!options.standard_uncertainty_over_rt.empty()
        && (options.standard_uncertainty_over_rt.size() != z.size()
            || !std::all_of(options.standard_uncertainty_over_rt.begin(), options.standard_uncertainty_over_rt.end(),
                [](double x) { return std::isfinite(x) && x >= 0.0; })))
        throw std::invalid_argument("invalid standard-term uncertainty estimates");
    const auto candidate = stability_candidate(c, result, report);
    const auto& at = candidate.physical;
    const auto& form = c.model->parameters().permittivity_form;
    std::vector<double> shift = report.standard_over_rt;
    for (std::size_t i = 0; i < shift.size(); ++i) shift[i] -= report.plane_over_rt[i];
    std::optional<GibbsDifferenceReference> own;
    auto& parent = reference ? *reference : own;
    std::vector<std::pair<StabilityTrial, GibbsDifference>> splits;  // one comparison per coincident gas endpoint
    // #21: after the completed schedule, one guarded pass over its condensed-domain refusals, funded by `recovery`.
    // A slot refused before optimization restarts from the parent restriction, a refused end point from its own seed.
    enum class Start { ParentRestriction, FirstAdmittedSeed, OwnSeed };
    struct Slot { std::size_t trial; Start start; std::vector<double> seed; };
    std::vector<Slot> slots;
    std::vector<std::vector<std::pair<std::vector<double>, double>>> admitted(supports.size());
    std::map<std::vector<double>, std::variant<std::string, DensityRoot>> liquid_roots;  // root or cause, by composition
    std::map<std::vector<double>, std::optional<std::string>> decisions;  // refusal cause, by prepared point
    // A proposal keeps the admitted liquid root of its own composition: denser than its negative-slope sample,
    // mechanically stable and within the anchor window of that root.
    const auto guard = [&](const PhasePoint& point) {
        const auto [decision, fresh] = decisions.try_emplace(point.point);
        if (!fresh) { ++report.work["recovery_guard_cache_hits"]; return decision->second; }
        auto& refusal = decision->second;
        try {
            auto w = point.composition[0];  // the unnormalized chart: its volume holds the chart total
            const double total = std::accumulate(w.begin(), w.end(), 0.0), rho = total / point.volume[0];
            for (double& x : w) x /= total;
            auto [root, missing] = liquid_roots.try_emplace(w);
            if (!missing) ++report.work["recovery_root_cache_hits"];
            else try {
                spend_credit(*recovery, report.work);
                auto found = stability_root(*c.model, at.temperature, at.pressure, w, PhaseKind::Liquid);
                condensed_admission(found, false);
                if (found.convergence.residual <= found.convergence.criterion->value) root->second = std::move(found);
                else root->second = "residual above its criterion";
            } catch (const std::exception& error) { root->second = error.what(); }
            const auto* liquid = std::get_if<DensityRoot>(&root->second);
            if (!liquid) refusal = "liquid root unavailable: " + std::get<std::string>(root->second);
            else if (!liquid->negative_slope_density) refusal = "liquid root without a negative-slope sample";
            else if (const double sample = *liquid->negative_slope_density; !(liquid->molar_density > sample && rho > sample))
                refusal = "point or root not denser than the negative-slope sample";
            else if (const double slope = c.model->state(at.temperature, rho, w).pressure_density_derivative; !(std::isfinite(slope) && slope > 0.0))
                refusal = "pressure-density slope not positive";
            else if (std::abs(std::log(rho / liquid->molar_density)) > policy::DENSITY_ANCHOR_LOG_WINDOW.value)
                refusal = "outside the anchor window of the liquid root";
        } catch (const std::exception& error) { refusal = error.what(); }
        if (refusal) ++report.work["recovery_guard_refusals"];
        return refusal;
    };
    for (std::size_t visit = 0, k = 0; visit < 2 * supports.size(); k = ++visit % supports.size()) {
        const bool retry = visit >= supports.size();
        const auto& phase = supports[k];
        if (retry && (!recovery || slots.empty())) break;
        if (retry) limits = recovery, recovery->admit = guard;
        std::string identity = phase.kind == PhaseKind::Liquid ? "Liquid/" : "Vapor/";
        for (auto i : phase.support) identity += std::to_string(i) + ",";
        const bool quotient = form.is(PermittivityForm::DensityTemperatureIonSuppression);
        const bool condensed = (quotient || !density_temperature_permittivity(form)) && !form.is(PermittivityForm::DipolarPerturbation)
            && std::any_of(phase.support.begin(), phase.support.end(), [&](auto i) { return z[i] != 0; });
        std::vector<std::pair<std::string, std::vector<double>>> seeds;
        std::vector<const Slot*> sources;
        if (options.evaluate_supplied_trials) seeds.emplace_back("supplied-composition", phase.composition_guess);
        else {
            std::vector<double> base(z.size(), 0.0), remainder;
            for (auto i : phase.support) base[i] = 1.0;
            base = reconstruct_trial_composition(base, phase.support, z, phase.kind);
            remainder = base;
            for (std::size_t p = 0; p < at.amount.size(); ++p) {
                if (std::all_of(phase.support.begin(), phase.support.end(), [&](auto i) { return at.composition[p][i] > 0.0; })) {
                    std::vector<double> restriction(z.size(), 0.0);
                    for (auto i : phase.support) restriction[i] = at.composition[p][i];
                    remainder = reconstruct_trial_composition(std::move(restriction), phase.support, z, phase.kind); break;
                }
            }
            const bool has_charge = std::any_of(phase.support.begin(), phase.support.end(), [&](auto i) { return z[i] != 0; });
            if (!fast || has_charge) seeds.emplace_back("positive-neutral-base", base);
            else if (phase.support.size() > 1) unsearched->push_back(identity + "/positive-neutral-base");
            if (fast && phase.kind == PhaseKind::Liquid && phase.support.size() == 1)
                seeds.emplace_back("near-pure-" + std::to_string(phase.support.front()), base);
            if (phase.support.size() > 1) for (auto j : phase.support) if (z[j] == 0.0) {
                if (fast && phase.kind == PhaseKind::Vapor) {
                    unsearched->push_back(identity + "/near-pure-" + std::to_string(j)); continue;
                }
                auto seed = remainder;
                double sum = 0.0; for (auto i : phase.support) if (i != j) sum += seed[i];
                for (auto i : phase.support) seed[i] = i == j ? 0.99 : 0.01 * seed[i] / sum;
                seeds.emplace_back("near-pure-" + std::to_string(j), std::move(seed));
            }
            if (phase.kind == PhaseKind::Vapor) {
                auto seed = base;
                double maximum = -std::numeric_limits<double>::infinity();
                for (auto i : phase.support) maximum = std::max(maximum, -shift[i]);
                for (auto i : phase.support) seed[i] = std::exp(-shift[i] - maximum);
                seeds.emplace_back("ideal-gas", std::move(seed));
            }
            for (const auto& [previous, composition] : removed)
                if (previous.kind == phase.kind && previous.support == phase.support)
                    seeds.emplace_back("removed-phase", composition);
            if (retry) seeds.clear();
            for (const auto& slot : slots) if (const auto& original = report.trials[slot.trial]; retry
                    && original.phase.kind == phase.kind && original.phase.support == phase.support) {
                seeds.emplace_back(original.seed, slot.start == Start::OwnSeed ? slot.seed : remainder);
                sources.push_back(&slot);
            }
        }
        auto& rooted = admitted[k];
        std::size_t visited = 0;
        for (auto& [seed_name, seed] : seeds) {
            const Slot* slot = retry ? sources[visited++] : nullptr;
            if (stopped) { unsearched->push_back(identity + "/" + seed_name); continue; }
            StabilityTrial trial; trial.phase = phase; trial.seed = seed_name; trial.evidence = "unresolved";
            auto start = slot ? slot->start : Start::ParentRestriction;  // in the first pass: the retry a refusal here gets
            bool searched = false;  // the start coincides with an already-searched original seed
            SuppliedPoint point = options.prepared_trials ? options.prepared_trials->points->entries[entry++] : SuppliedPoint{};
            const auto advance = [&](SuppliedStage stage, CompiledSystem* c = nullptr, double eta = 0.0) {
                advance_supplied(point, stage, *input.model, at.temperature, phase, c, eta, report.work);
            };
            try {
                if (options.evaluate_supplied_trials) { advance(SuppliedStage::Seed); seed = point.seed; }
                else if (!slot) { seed = reconstruct_trial_composition(std::move(seed), phase.support, z, phase.kind);
                    ++report.work["reconstruction_evaluations"]; }  // a second-pass start is reconstructed already
                double eta, seed_density = 0.0;
                if (options.evaluate_supplied_trials && phase.packing_guess) eta = *phase.packing_guess;
                else {
                    std::optional<DensityRoot> root;
                    while (!root) {
                        if (slot) spend_credit(*recovery, report.work);
                        root = stability_root(*c.model, at.temperature, at.pressure, seed, phase.kind);
                        core::retain_first(trial.optimizer.first_numerical_failure, root->first_numerical_failure);
                        if (condensed) try { condensed_admission(root, quotient); }
                        catch (const CondensedRefusal&) {  // a refused parent restriction: the first admitted seed
                            if (!slot || start != Start::ParentRestriction || rooted.empty()) throw;
                            start = Start::FirstAdmittedSeed;
                            root.reset(), seed = rooted.front().first;
                        }
                    }
                    seed_density = root->molar_density;
                    eta = c.model->packing_fraction_coefficient(at.temperature, seed) * root->molar_density;
                }
                packing_admission(eta);
                searched = std::any_of(rooted.begin(), rooted.end(), [&](const auto& prior) {
                    return prior.first == seed && std::abs(std::log(prior.second / eta)) <= policy::DENSITY_ANCHOR_LOG_WINDOW.value; });
                if (!slot && searched) continue;
                if (!slot) rooted.emplace_back(seed, eta);
                const bool spent = limits && limits->expired();  // #987: the trials so far stand
                if (spent || (options.trial_budget && report.work["trial_invocations"] >= *options.trial_budget)) {
                    if (slot) throw BudgetExhausted(spent ? "recovery allowance exhausted" : "trial budget exhausted");
                    report.message = std::string("search-incomplete: ") + (spent ? "work" : "trial") + " budget exhausted";
                    if (!fast) { report_work(report.work, work); return report; }
                    stopped = true; unsearched->push_back(identity + "/" + seed_name); continue;
                }
                const bool native = !c.has_ions && !c.has_reaction && !options.evaluate_supplied_trials;
                CompiledSystem compiled;  // the shared driver's trial; none for a converged Michelsen trial (rank 0)
                const auto compile = [&] {
                    if (options.evaluate_supplied_trials) compiled = supplied_system(systems, *c.model,
                        at.temperature, at.pressure, phase, point, shift, trial.row_transformation, report.work);
                    else compiled = compile_stability_trial(*c.model, at.temperature, at.pressure, phase, seed, shift,
                        trial.row_transformation);
                    trial.basis = compiled.B;
                };
                if (!native) compile();
                core::Scope scope("tangent_plane_trial");
                ++report.work["trial_invocations"];
                const auto found = native ? tangent_plane_trial(*c.model, at.temperature, at.pressure, phase, seed, seed_density, shift,
                    trial, report.work, limits) : std::nullopt;
                if (limits) limits->check();  // an interrupted trial starts no driver fallback
                if (native && !found) compile();
                if (options.evaluate_supplied_trials) {
                    advance(SuppliedStage::Coordinates, &compiled, eta);
                    trial.raw_coordinates = point.coordinates;
                }
                else if (!found) {
                    const auto initial = trial_coordinates(compiled, seed, eta);
                    StageWork prior; prior.observe(trial.optimizer);
                    trial.optimizer = run_driver(compiled, initial, structural_result(compiled, scope),
                        scope, &trial.raw_coordinates, limits, Attempt::Cold, {}, false, true);
                    prior.report_observations(trial.optimizer);
                    count_driver(report.work, trial.optimizer, "trial_driver_invocations");
                    trial.message = trial.optimizer.message;
                }
                if (!slot) start = Start::OwnSeed;
                if (trial.raw_coordinates.empty()) throw std::invalid_argument("no admitted nonfatal optimizer point");
                if (options.evaluate_supplied_trials) advance(SuppliedStage::Raw, &compiled, eta);
                auto raw = options.evaluate_supplied_trials ? point.raw : found ? found->first : compiled.prepare(trial.raw_coordinates, true);
                if (options.evaluate_supplied_trials) { raw.pressure = at.pressure; raw.standard = *compiled.fixed_standard; }
                if (!found) trial.raw_residuals = compiled.evaluate_at(raw).values;
                if (!options.evaluate_supplied_trials) { point.raw = raw; point.next = SuppliedStage::Composition; }
                advance(SuppliedStage::Composition);
                trial.composition = point.composition;
                const double volume = point.volume;
                if (phase.kind == PhaseKind::Vapor && std::any_of(phase.support.begin(), phase.support.end(),
                        [&](auto i) { return z[i] != 0; })) {
                    auto& q = trial.quantities;
                    q["packing_fraction"] = raw.eta[0]; q["molar_density_mol_m3"] = 1.0 / volume;
                    advance(SuppliedStage::State);
                    const auto& state = point.state;
                    q["pressure_pa"] = state.pressure;
                    q["pressure_density_derivative"] = state.pressure_density_derivative;
                    const auto root = declared_root(*c.model, at.temperature, at.pressure, trial.composition, PhaseKind::Vapor);
                    q["vapor_root_log_displacement"] = std::abs(std::log(root.molar_density * volume));
                    q["vapor_root_residual"] = root.convergence.residual;
                    if (!(q["vapor_root_log_displacement"] <= policy::DENSITY_ANCHOR_LOG_WINDOW.value
                            && root.convergence.residual <= root.convergence.criterion->value
                            && std::isfinite(state.pressure_density_derivative) && state.pressure_density_derivative > 0.0))
                        throw std::invalid_argument("charged Vapor point outside its own vapor pressure-root window");
                }
                if (condensed && slot) spend_credit(*recovery, report.work);
                if (condensed) condensed_admission(anchored_root(*c.model, at.temperature, at.pressure, trial.composition,
                    volume, PhaseKind::Liquid, trial.quantities, "condensed_"), quotient, &trial.quantities);
                advance(SuppliedStage::Energy);
                trial.quantities.merge(phase_energy(*c.model, at.temperature, at.pressure, volume, trial.composition,
                    report.standard_over_rt, report.work, nullptr, &point.energy));
                auto& q = trial.quantities;
                double charge = 0.0, sum = -1.0, plane = 0.0, displacement = 0.0;
                for (auto i : phase.support) {
                    sum += trial.composition[i]; charge += z[i] * trial.composition[i];
                    plane += trial.composition[i] * report.plane_over_rt[i];
                    displacement += std::abs(trial.composition[i] - raw.composition[0][i]);
                }
                q["normalization_residual"] = sum; q["charge_residual"] = charge;
                q["reconstruction_l1"] = displacement;
                q["tpd_over_rt"] = q.at("g_over_rt") - plane;
                if (options.evaluate_supplied_trials) {
                    q["packing_fraction"] = raw.eta[0]; q["molar_density_mol_m3"] = 1.0 / volume;
                    double displacement = 0.0;
                    for (std::size_t i = 0; i < z.size(); ++i)
                        displacement += std::abs(trial.composition[i] - phase.composition_guess[i]);
                    q["input_reconstruction_l1"] = displacement;
                } else {
                    double dual = found ? found->second : 0.0;
                    for (std::size_t k = 0; k < compiled.rank; ++k)
                        dual += trial.raw_coordinates[compiled.primal_size + k] * compiled.balance[k];
                    q["trial_dual_over_rt"] = dual; q["trial_identity_defect_over_rt"] = q["tpd_over_rt"] - dual;
                    if (!found) q["raw_objective_over_rt"] = trial.optimizer.objective;
                }
                trial.evidence = "no_resolved_descent";
                const bool negative = q["tpd_over_rt"] < 0.0, explicit_negative = diagnostics && negative && !options.evaluate_supplied_trials;
                bool gas = false, coincident = false;
                if (explicit_negative) {
                    for (std::size_t p = 0; p < at.amount.size(); ++p) {
                        StabilityPhaseComparison comparison; comparison.phase = p;
                        const auto& support = c.problem.phases[p].support;
                        std::set_difference(phase.support.begin(), phase.support.end(), support.begin(), support.end(),
                            std::back_inserter(comparison.added_species));
                        std::set_difference(support.begin(), support.end(), phase.support.begin(), phase.support.end(),
                            std::back_inserter(comparison.missing_species));
                        comparison.composition_max_distance = largest_difference(trial.composition, at.composition[p]);
                        comparison.volume_log_distance = std::abs(std::log(volume / at.volume[p]));
                        coincident |= support == phase.support && c.problem.phases[p].kind == PhaseKind::Vapor
                            && comparison.composition_max_distance <= policy::PHASE_COINCIDENCE_DISTANCE.value
                            && comparison.volume_log_distance <= policy::PHASE_COINCIDENCE_DISTANCE.value;
                        trial.present_phase_comparisons.push_back(std::move(comparison));
                    }
                    if (slot) spend_credit(*recovery, report.work);
                    try {
                        if (!retain_density_side(*c.model, at.temperature, at.pressure, trial.composition, volume, phase.kind, q))
                            trial.diagnostic_cause = density_refusal;
                        else gas = q.at("root_density_side") == static_cast<double>(DensityBranch::Vapor);
                    } catch (const std::exception& error) { trial.diagnostic_cause = error.what(); }
                }
                // A distinct reacting gas-side trial uses the relaxed split if the frozen check leaves it indeterminate.
                const bool reacting = explicit_negative && c.has_reaction && trial.optimizer.success;
                if (options.evaluate_supplied_trials || negative)
                    donor_descent(c, candidate, trial, report, options, limits, reacting && gas && !coincident, slot);
                if (reacting && gas && !coincident && trial.evidence == "indeterminate") {
                    auto shared = std::find_if(splits.begin(), splits.end(), [&](const auto& s) {
                        return s.first.donor == trial.donor && same_stationary_trial(s.first, trial, false);
                    });
                    if (shared == splits.end()) {
                        if (!parent) {  // the original plane and repair, their numerical estimates counted here
                            auto counted = result;
                            const auto prepared = input.prepare(result.coordinates);
                            if (!counted.numerical_estimates) counted.numerical_estimates = equilibrium_numerical_estimates(input, result, &prepared, true);
                            parent = conserved_split_reference(input, counted);
                        }
                        shared = splits.emplace(splits.end(), trial, relaxed_split(c, at, trial, *trial.donor, *parent, options, report.work));
                    }
                    const auto& d = shared->second;
                    // A rank-refused two-phase cell has rank one: each phase holds its own reaction equilibrium,
                    // the donor potentials are the plane and the split change is at most beta * tpd (equal when the
                    // reactions conserve mole number), so a frozen descent is a witness and the frozen slope and splits decide.
                    if (d.status == GibbsDifferenceStatus::Unavailable && d.cause.rfind("phase count exceeds balance rank", 0) == 0) {
                        donor_descent(c, candidate, trial, report, options, limits, false, slot);
                        trial.message += "; reaction-relaxed split: " + d.cause;
                    } else classify_split(c, report, trial, d);
                } else if (reacting && !trial.diagnostic_cause.empty() && trial.evidence == "no_resolved_descent") {
                    trial.evidence = "indeterminate";  // without an admitted side the frozen slope decides nothing
                    trial.message += "; " + trial.diagnostic_cause;
                }
            } catch (const std::exception& error) {
                trial.message += "; " + std::string(error.what());
                trial.evidence = "unresolved";
                if (recovery && !retry && !quotient && dynamic_cast<const CondensedRefusal*>(&error))
                    slots.push_back({report.trials.size(), start, seed});
            }
            if (slot) {  // only a determined optimized retry replaces the original slot; its message keeps the refusal
                auto& original = report.trials[slot->trial];
                // A start at an already-searched seed or a returned present phase is no independent coverage.
                bool present = false;
                for (std::size_t i = 0; i < at.amount.size() && !trial.composition.empty(); ++i) present |= c.problem.phases[i].kind
                    == phase.kind && largest_difference(trial.composition, at.composition[i]) <= policy::PHASE_COINCIDENCE_DISTANCE.value;
                const char* from[] = {"the parent-liquid restriction", "the first admitted seed", "its own seed"};
                trial.message = original.message + "; second pass from " + from[static_cast<int>(start)]
                    + (searched ? ", an already-searched seed" : "") + (present ? ", returned a present phase" : "") + ": " + trial.message;
                if (trial.optimizer.success && (trial.resolved_descent || trial.evidence == "no_resolved_descent")) original = std::move(trial);
                else original.message = std::move(trial.message);
                continue;
            }
            report.trials.push_back(std::move(trial));
            if (fast && report.trials.back().resolved_descent) stopped = true;
        }
    }
    report.completed = !stopped;
    // The one trial classification (#987): a trial without a resolved descent supplies no-descent evidence only from a
    // converged point; an unconverged driver point, a failure or an unresolved descent is unresolved.
    const auto unresolved = std::count_if(report.trials.begin(), report.trials.end(), [](const auto& t) {
        return t.evidence == "unresolved" || t.evidence == "indeterminate" || t.evidence == "negative_plane_without_positive_donor"
            || (t.evidence == "no_resolved_descent" && !t.optimizer.success); });
    const auto descending = std::count_if(report.trials.begin(), report.trials.end(), [](const auto& t) { return t.resolved_descent; });
    report.no_descent_in_trial_set = report.completed && !report.trials.empty() && !descending && !unresolved;
    if (report.message.empty()) report.message = descending ? "numerically resolved descending split with estimated uncertainty"
        : "no witness in the recorded trial set";
    report.message += "; unresolved trials " + std::to_string(unresolved);
    if (options.evaluate_supplied_trials) report.message = "supplied trial evaluation; unavailable trials " + std::to_string(unresolved);
    if (!options.trials && std::all_of(supports.front().support.begin(), supports.front().support.end(),
        [&](auto j) { return z[j] != 0; })) report.message += "; no neutral species";
    report_work(report.work, work);
    return report;
}

// `reference` keeps the original plane and repair a reaction-relaxed split computed, for the measure.
StabilityReport explicit_test(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options,
    std::optional<GibbsDifferenceReference>* reference, DriverLimits* recovery = nullptr) {
    if (options.mode != DiscoveryMode::Full) throw std::invalid_argument("discovery mode is only admitted in Discover");
    if (c.ideal) throw std::invalid_argument("ideal phase stability unavailable");
    if (c.model->parameters().ideal_gas) throw std::invalid_argument("ideal-gas stability testing is deferred");
    if (!c.problem.solids.empty()) throw std::invalid_argument("stability trials are not admitted with solids");
    return test_snapshot(c, result, options, {}, nullptr, nullptr, true, reference,
        options.evaluate_supplied_trials || options.prepared_trials ? nullptr : recovery);
}
} // namespace

StabilityReport recovered_stability_test(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options,
    DriverLimits recovery) {
    return explicit_test(c, result, options, nullptr, &recovery);
}
StabilityReport stability_test(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options) {
    core::Scope scope("stability_test");
    auto recovery = DriverLimits::recovery_allowance(std::chrono::steady_clock::time_point::max(),
        static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value));
    return explicit_test(c, result, options, nullptr, options.retry_condensed_refusals ? &recovery : nullptr);
}

namespace {
StabilityMeasure measure_value(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options,
    std::optional<GibbsDifferenceReference>& parent) {
    StabilityMeasure out;
    if (options.mode != DiscoveryMode::Full || options.evaluate_supplied_trials || options.prepared_trials)
        throw std::invalid_argument("stability requirements need Full optimized trials");
    {
        core::Scope test("stability_test");
        out.report = explicit_test(c, result, options, &parent);
    }
    if (!out.report.completed || out.report.trials.empty()) { out.cause = "incomplete trial schedule"; return out; }
    for (std::size_t i = 0; i < out.report.trials.size(); ++i) {
        const auto& t = out.report.trials[i];
        const auto at = t.quantities.find("tpd_over_rt");
        if (!t.optimizer.success || at == t.quantities.end() || !std::isfinite(at->second)
            || (t.evidence != "no_resolved_descent" && !t.resolved_descent)) {
            out.cause = t.message.empty() ? "failed final trial optimizer or undetermined evidence: " + std::to_string(i)
                : "trial " + std::to_string(i) + " (" + t.evidence + "): " + t.message;
        }
        if (at == t.quantities.end() || !std::isfinite(at->second)) continue;
        if (!out.raw_minimum || at->second < *out.raw_minimum) {
            out.raw_minimum = at->second; out.minimizing_trial = i;
        }
    }
    if (!out.cause.empty()) return out;
    try {
        const auto at = c.prepare(result.coordinates, true);
        for (std::size_t p = 0; p < out.report.present_phases.size(); ++p)
            if (!retain_density_side(*c.model, at.temperature, at.pressure, at.composition[p], at.volume[p],
                c.problem.phases[p].kind, out.report.present_phases[p])) { out.cause = density_refusal; return out; }
        for (auto& t : out.report.trials)
            if (!retain_density_side(*c.model, at.temperature, at.pressure, t.composition,
                t.quantities.at("volume_m3_per_mol"), t.phase.kind, t.quantities)) { out.cause = density_refusal; return out; }
    } catch (const std::exception& error) { out.cause = error.what(); return out; }
    auto rows = parent ? parent->row_errors : std::vector<double>{};  // the split's original plane errors, when made
    if (rows.size() != c.n_feed.size()) {
        rows.clear();
        std::vector<SolvedStateObservable> plane;
        for (std::size_t i = 0; i < c.n_feed.size(); ++i) plane.push_back({SolvedStateObservableKind::SpeciesStandardMinusPlane, 0, i});
        SolvedStateActionDirection zero; zero.feed_amounts_mol.assign(c.n_feed.size(), 0.0);
        zero.active_parameters.assign(c.model->active_parameters().size(), 0.0);
        const auto actions = solved_state_first_actions(c, result, plane, {zero}).front();
        for (const auto& action : actions.results) {
            if (!action.numerical_error) { out.cause = "parent plane uncertainty unavailable"; return out; }
            rows.push_back(*action.numerical_error);
        }
    }
    std::vector<double> error(out.report.trials.size());
    for (std::size_t j = 0; j < out.report.trials.size(); ++j) {
        auto& t = out.report.trials[j];
        error[j] = energy_estimate(t.quantities);
        for (std::size_t i = 0; i < t.composition.size(); ++i)
            error[j] += t.composition[i] * (rows[i]
                + (options.standard_uncertainty_over_rt.empty() ? 0.0 : options.standard_uncertainty_over_rt[i])
                + std::numeric_limits<double>::epsilon() * std::abs(out.report.plane_over_rt[i]));
        for (const auto* name : {"donor_correction_over_rt", "feasibility_error_estimate_over_rt"})
            if (t.quantities.count(name)) error[j] += std::abs(t.quantities.at(name));
        t.quantities["plane_value_error_estimate_over_rt"] = error[j];
    }
    const auto best = *out.minimizing_trial; out.uncertainty = error[best];
    const auto& minimum = out.report.trials[best];
    for (std::size_t j = 0; j < out.report.trials.size(); ++j) {
        const auto& t = out.report.trials[j];
        if (t.quantities.at("tpd_over_rt") - *out.raw_minimum > error[j] + error[best]) continue;
        out.minima.push_back(j);
        const bool same = same_stationary_trial(t, minimum, false);
        if (!same) out.event = "distinct stationary minima have overlapping value intervals";
    }
    if (!out.event.empty()) { out.cause = out.event; return out; }
    out.value = out.report.no_descent_in_trial_set ? 0.0 : *out.raw_minimum;
    out.available = true;
    return out;
}
} // namespace

// Work totals include the search and every evaluation after it.
StabilityMeasure stability_measure(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options) {
    core::Scope scope("stability_measure");
    core::CalculationWork work;
    std::optional<GibbsDifferenceReference> parent;
    auto out = measure_value(c, result, options, parent);
    report_work(out.report.work, work);
    return out;
}

// ALG-STABILITY before neutral recovery (#987): the Full trial set on the lower-objective single-phase L1/V1
// declaration of the feed, every solve and native trial charged to `limits`. Every outcome returns its report;
// nothing from it feeds recovery. A failed parent keeps its cause and work.
StabilityReport recovery_stability(const CompiledSystem& input, const DriverLimits& limits) {
    auto c = input; c.phase_calculation = c.problem.phase_calculation = PhaseCalculation::Discovery;
    StabilityReport parents, report;
    std::optional<std::pair<CompiledSystem, EquilibriumResult>> best;
    for (const auto kind : {PhaseKind::Liquid, PhaseKind::Vapor}) {
        auto single = c.problem;
        single.phases = {PhaseSpec{kind == PhaseKind::Liquid ? "L1" : "V1", kind}};
        auto compiled = compile_problem(*c.model, single);
        auto parent = try_solve_compiled(compiled, std::nullopt, nullptr, &limits);  // a failed parent keeps its work
        count_driver(parents.work, parent, "parent_declarations");
        core::retain_first(parents.first_numerical_failure, parent.first_numerical_failure);
        if (!parent.success) parents.message += "; " + single.phases[0].name + ": " + parent.message;
        if (parent.success && (!best || parent.objective < best->second.objective)) best.emplace(std::move(compiled), std::move(parent));
    }
    core::CalculationWork work;
    try {
        if (!best) throw std::runtime_error("no single-phase candidate");
        report = test_snapshot(best->first, best->second, {}, {}, nullptr, &limits);
        report.message = " from single-phase " + best->first.problem.phases[0].name + ", "
            + std::to_string(report.trials.size()) + " trials: " + report.message;
    } catch (const std::bad_alloc&) { throw; }
    catch (const std::exception& error) {  // only before the first trial: the trials' own catches keep theirs
        report.message = std::string(": ") + error.what();
        core::retain_failure(report.first_numerical_failure, error, "recovery stability");
        report_work(report.work, work);
    }
    for (const auto& [name, value] : parents.work) report.work[name] += value;
    core::retain_first(report.first_numerical_failure, parents.first_numerical_failure, true);
    report.message = "tangent-plane search before recovery" + report.message + parents.message;
    return report;
}

// ALG-STABILITY: can the declared vanished phase appear? The parent is the declaration without it, accepted only on
// success; one explicit trial of its kind and declared support, its existing seeds, every solve charged to `limits`.
StabilityReport vanished_phase_stability(const CompiledSystem& input, std::size_t phase, const DriverLimits& limits) {
    auto declaration = input.problem;
    const auto vanished = declaration.phases[phase];
    declaration.phases.erase(declaration.phases.begin() + static_cast<std::ptrdiff_t>(phase));
    const std::string kind = vanished.kind == PhaseKind::Vapor ? "vapor" : "liquid";
    StabilityReport report;
    EquilibriumResult parent;  // also a failed parent's work
    try {
        const auto compiled = compile_problem(*input.model, declaration);
        parent = try_solve_compiled(compiled, std::nullopt, nullptr, &limits);
        if (!parent.success) throw std::runtime_error(parent.message);
        StabilityOptions options;
        options.trials = std::vector<PhaseSpec>{vanished};
        report = test_snapshot(compiled, parent, options, {}, nullptr, &limits);
        std::optional<double> raw;  // descending trials, or all when none descends or is unresolved; not clipped
        std::string support;
        for (const auto& t : report.trials) if (t.resolved_descent || report.no_descent_in_trial_set)
            raw = std::min(raw.value_or(t.quantities.at("tpd_over_rt")), t.quantities.at("tpd_over_rt"));
        for (const auto i : report.trials.empty() ? vanished.support : report.trials.front().phase.support) support += (support.empty() ? "" : ", ") + std::to_string(i);
        report.message = report.no_descent_in_trial_set ? "no descent against " + kind + " trials of support {" + support
                + "} in this finite search: minimum TPD/RT = " + core::optional_number(raw) + "; not a stability certificate"
            : raw ? "a " + kind + " trial descends from the surviving declaration: TPD/RT = " + core::optional_number(raw)
            : kind + " trial search unresolved: " + report.message;
    } catch (const std::bad_alloc&) { throw; }
    catch (const std::exception& error) {  // before the first trial: the trials' own catches keep theirs
        report.message = (parent.success ? kind + " trial search unresolved: " : "surviving declaration unresolved: ") + error.what();
    }
    count_driver(report.work, parent, "parent_declarations");
    core::retain_first(report.first_numerical_failure, parent.first_numerical_failure, true);
    return report;
}

std::pair<EquilibriumResult, std::vector<StabilityReport>> cooling_onset(const CompiledSystem& input) {
    core::Scope scope("cooling_onset");
    auto c = input; c.phase_calculation = c.problem.phase_calculation = PhaseCalculation::Discovery;
    if (c.ideal) throw std::invalid_argument("ideal cooling onset unavailable");
    if (c.model->parameters().ideal_gas) throw std::invalid_argument("ideal-gas cooling onset is deferred");
    const auto& request = c.problem;
    if (c.has_reaction || c.has_ions || request.discovery || request.critical || request.spinodal || request.envelope
        || request.enthalpy_j || request.entropy_j_per_k || !request.temperature.free || request.pressure.free
        || request.phases.size() != 2 || c.feed_composition.size() < 2)
        throw std::invalid_argument("cooling onset requires a neutral fixed-pressure dew declaration");
    const auto& phases = request.phases;
    const auto liquid = phases[0].kind == PhaseKind::Liquid ? 0u : 1u, vapor = 1u - liquid;
    if (phases[liquid].kind != PhaseKind::Liquid || phases[vapor].kind != PhaseKind::Vapor
        || !phases[liquid].is_zero_pinned()
        || phases[vapor].amount.kind != AmountKind::Free
        || !std::all_of(phases.begin(), phases.end(), [&](const auto& phase) {
            return phase.composition_guess.empty() && !phase.packing_guess && phase.support.size() == c.feed_composition.size();
        }) || !std::all_of(c.feed_composition.begin(), c.feed_composition.end(), [](double z) { return z > 0.0; }))
        throw std::invalid_argument("cooling onset requires full supports and generated phase starts");
    auto failed = solve_compiled(c, std::nullopt, nullptr, nullptr, {}, true);
    // Search only after a resolved local instability or the optimizer's own infeasibility at a point with no
    // topology event; a rejected publication (coalesced, vanished), an unavailable phase spectrum and a
    // structurally rank-deficient declaration (no solve run) return unchanged.
    if (!failed.success && failed.candidate != equilibrium::Candidate::Unstable
        && !(c.feasible && failed.convergence.stop == core::Stop::infeasible && failed.convergence.cause.empty()
             && failed.topology_event == TopologyEvent::None))
        return {std::move(failed), {}};
    std::vector<StabilityReport> reports; StageWork work; std::vector<double> composition;
    // The existing starts owner retains work when a reference or corrector throws after entering the driver.
    const auto solve = [&](const EquilibriumProblem& problem) {
        auto attempt = solve_equilibrium_starts(*c.model, problem, {EquilibriumStart{}}, {});
        if (!attempt.results[0]) throw std::runtime_error(attempt.causes[0]);
        return std::move(*attempt.results[0]);
    };
    auto reference = c.problem; reference.phases = {c.problem.phases[vapor]};
    StabilityOptions options; options.trials = std::vector<PhaseSpec>{c.problem.phases[liquid]};
    options.trial_budget = c.model->parameters().component_count() + 1;
    const auto sample = [&](double temperature) {
        reference.temperature = FreeValue::fixed(temperature);
        // Sufficient, not necessary, vapor-side evidence: more than one stable root, or a sole root below a negative-slope sample (#423).
        const auto root = [&] {  // the whole added scan counts, work before a throw included
            core::CalculationWork scan;
            const auto count = [&] { work.evaluations += scan.counts().root_evaluations; };
            try { auto found = stability_root(*c.model, temperature, request.pressure.value, c.feed_composition, PhaseKind::Vapor); count(); return found; }
            catch (...) { count(); throw; }
        }();
        const auto central = solve(reference);
        work.add(central);
        if (!central.success || (root.stable_root_count.value_or(0) < 2
            && !(root.negative_slope_density && central.molar_densities[0] < *root.negative_slope_density)))
            throw std::runtime_error("vapor reference unavailable: " + (central.success ? "sole density root not below a negative-slope sample" : central.message));
        reports.push_back(stability_test(compile_problem(*c.model, reference), central, options));
        auto& report = reports.back();
        work.iterations += report.work["iterations"]; work.evaluations += report.work["evaluations"]; work.hessian_calls += report.work["hessian_calls"];
        if (!report.completed || report.trials.empty()) throw std::runtime_error("incomplete liquid trial family: " + report.message);
        double best = std::numeric_limits<double>::infinity();
        for (const auto& trial : report.trials) {
            const auto distance = trial.quantities.find("tpd_over_rt");
            if (!trial.optimizer.success || distance == trial.quantities.end() || !std::isfinite(distance->second))
                throw std::runtime_error("liquid trial unavailable: " + trial.message);
            if (distance->second < best) { best = distance->second; composition = trial.composition; }
        }
        return best;
    };
    const auto competing = [&](double temperature) {
        sample(temperature);
        return std::any_of(reports.back().trials.begin(), reports.back().trials.end(), [](const auto& t) { return t.resolved_descent; });
    };
    EquilibriumResult recovered;
    try {
        // An ordinary dew stands only when no competing liquid descends at its temperature.
        if (failed.success && !competing(*failed.temperature)) {
            StageWork(work).report(failed);
            return {std::move(failed), std::move(reports)};
        }
        double upper = c.seed_temperature();
        for (int step = 1; sample(upper) <= 0.0; ++step) {
            if (step == 16) throw std::runtime_error("no vapor-side upper bound in 16 heating samples");
            upper *= 1.01;
        }
        const double increment = 0.005 * upper;
        double lower = upper - increment;
        for (int step = 1; sample(lower) > 0.0; ++step) {
            if (step == 64) throw std::runtime_error("no cooling onset in 64 downward samples");
            upper = lower; lower = upper - increment;
        }
        auto boundary = c.problem;
        boundary.temperature = FreeValue::variable(0.5 * (upper + lower));
        boundary.phases[liquid].composition_guess = composition;
        boundary.phases[liquid].declared_composition_guess = false;
        boundary.phases[vapor].composition_guess = c.feed_composition;
        boundary.phases[vapor].declared_composition_guess = false;
        recovered = solve(boundary);
        if (!recovered.success) throw std::runtime_error("boundary corrector unavailable: " + recovered.message);
        if (*recovered.temperature < lower - 1e-10 * upper || *recovered.temperature > upper + 1e-10 * upper)
            throw std::runtime_error("boundary corrector left cooling-onset bracket");
        if (competing(*recovered.temperature)) throw std::runtime_error("competing liquid descent at corrected dew");
        work.add(failed); StageWork(work).report(recovered);
        return {std::move(recovered), std::move(reports)};
    } catch (const std::exception& error) {
        work.add(recovered);
        if (failed.success) {  failed.message = "ordinary dew unconfirmed"; refuse(failed); }
        failed.message += "; cooling onset refused: " + std::string(error.what());
        std::move(work).report(failed);
        return {std::move(failed), std::move(reports)};
    }
}

// The event belongs solely to check_topology. Merge finite inventories and extensive volumes,
// then reclose the recipient on its declared density branch before the reduced solve.
std::pair<EquilibriumProblem, std::vector<double>> remove_discovery_phase(
    const CompiledSystem& c, const std::vector<double>& coordinates, Quantities* retained_work) {
    Quantities local_work;
    auto& work = retained_work ? *retained_work : local_work;
    const auto event = check_topology(c, coordinates);
    if (event.event == TopologyEvent::None) throw std::invalid_argument("no topology removal event");
    const auto at = c.prepare(coordinates, true);
    std::size_t dropped = event.phase.value_or(event.pair ? event.pair->second : 0);
    std::optional<std::size_t> recipient;
    if (event.pair) recipient = event.pair->first;
    const auto contains = [&](std::size_t p) {
        const auto& support = c.problem.phases[p].support;
        return std::all_of(c.problem.phases[dropped].support.begin(), c.problem.phases[dropped].support.end(),
            [&](auto i) { return std::find(support.begin(), support.end(), i) != support.end(); });
    };
    if (!recipient) for (std::size_t p = 0; p < at.amount.size(); ++p)
        if (p != dropped && contains(p)) { recipient = p; break; }
    if (!recipient || !contains(*recipient)) throw std::invalid_argument("removal has no inventory recipient");
    if (!(at.amount[dropped] >= 0.0 && at.amount[*recipient] > 0.0))
        throw std::invalid_argument("removal attempt has invalid amounts");
    std::vector<std::vector<double>> n(at.amount.size(), std::vector<double>(c.n_feed.size(), 0.0));
    std::vector<double> before(c.n_feed.size(), 0.0), after(before);
    for (std::size_t p = 0; p < n.size(); ++p) for (std::size_t i = 0; i < before.size(); ++i) {
        n[p][i] = at.amount[p] * at.composition[p][i]; before[i] += n[p][i];
    }
    for (std::size_t i = 0; i < before.size(); ++i) n[*recipient][i] += n[dropped][i];
    const double merged_volume = at.amount[*recipient] * at.volume[*recipient] + at.amount[dropped] * at.volume[dropped];
    EquilibriumProblem reduced = c.problem;
    reduced.phases.erase(reduced.phases.begin() + static_cast<std::ptrdiff_t>(dropped));
    std::vector<double> seed;
    const auto& z = c.model->parameters().charge_number;
    for (std::size_t p = 0; p < n.size(); ++p) if (p != dropped) {
        const double amount = std::accumulate(n[p].begin(), n[p].end(), 0.0);
        std::vector<double> w(n[p].size(), 0.0);
        double charge = 0.0, sum = 0.0;
        for (auto i : c.problem.phases[p].support) {
            w[i] = n[p][i] / amount;
            if (!(w[i] > 0.0 && std::isfinite(w[i]))) throw std::invalid_argument("removal destroys support");
            sum += w[i]; charge += z[i] * w[i]; after[i] += amount * w[i]; seed.push_back(std::log(w[i]));
        }
        if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value)
            throw std::invalid_argument("removal arithmetic feasibility failed");
        seed.push_back(amount);
        double volume = p == *recipient ? merged_volume / amount : at.volume[p];
        packing_admission(c.model->packing_fraction_coefficient(at.temperature, w) / volume);
        if (p == *recipient) {
            const auto root = declared_root(*c.model, at.temperature, at.pressure, w, c.problem.phases[p].kind);
            volume = 1.0 / root.molar_density;
        }
        seed.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature, w) / volume));
    }
    std::vector<double> error(before.size());
    for (std::size_t i = 0; i < before.size(); ++i) error[i] = after[i] - before[i];
    work["removal_component_error_mol"] = maximum_abs(error);
    work["removal_moiety_error_mol"] = maximum_abs(matvec(c.B, error));
    if (maximum_abs(error) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total)
        || maximum_abs(matvec(c.B, error)) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
        throw std::invalid_argument("removal loses inventory");
    return {std::move(reduced), std::move(seed)};
}

EquilibriumResult discover_phases(const HelmholtzModel& model, const EquilibriumProblem& input) {
    auto request = input; request.phase_calculation = PhaseCalculation::Discovery;
    if (model.parameters().ideal_gas) throw std::invalid_argument("ideal-gas phase discovery is deferred");
    if (!request.discovery || request.topology_declared || !request.phases.empty()
        || !request.temperature.is_positive() || !request.pressure.is_positive()
        || request.temperature.free || request.pressure.free)
        throw std::invalid_argument("discovery requires explicit unknown topology at fixed positive T/P");
    if (!request.solids.empty()) throw std::invalid_argument("Discover() is not admitted with solids");
    core::CalculationWork work;
    const double row_tolerance = request.solver_options.resolved(model.solver_options()).value(policy::row);
    const auto& options = *request.discovery;
    const bool fast = options.mode != DiscoveryMode::Full;
    if (fast && (!request.reactions.empty() || options.trials))
        throw std::invalid_argument("fast discovery requires no reactions and default trial supports");
    if (options.prepared_trials) throw std::invalid_argument("prepared trials are not phase discovery");
    if (options.evaluate_supplied_trials) throw std::invalid_argument("supplied trial evaluation is not phase discovery");
    const auto feed = resolve_feed(model, request.feed);
    const auto& z = model.parameters().charge_number;
    std::vector<std::size_t> united;
    bool has_charge = false;
    const auto participant = reaction_participants(feed.size(), request.reactions);
    for (std::size_t i = 0; i < feed.size(); ++i)
        if (feed[i] > 0.0 || participant[i]) { united.push_back(i); has_charge |= z[i] != 0.0; }
    auto record = std::make_shared<DiscoveryRecord>();
    record->mode = options.mode;
    if (fast) record->unsearched = {"removed-phase starts and retries", "alternative descending additions", "post-flash stability",
        "three-or-more-phase topologies"};
    EquilibriumProblem declaration = request;
    declaration.discovery.reset(); declaration.topology_declared = true;
    PhaseSpec initial; initial.name = "L1"; initial.support = united; initial.support_declared = true;
    if (!request.reactions.empty()) {
        const auto base = reconstruct_trial_composition(std::vector<double>(feed.size(), 1.0), united, z);
        for (auto i : united) initial.composition_guess.push_back(base[i]);
        initial.declared_composition_guess = false;
    }
    declaration.phases = {initial};
    auto compiled = compile_problem(model, declaration);
    // Admit the requested trial physics before any declaration solve.
    stability_supports(compiled, options);
    const std::size_t rank = compiled.rank;
    const std::size_t budget = options.declaration_budget.value_or(2 * rank);
    std::size_t calls = 0;
    std::vector<std::pair<PhaseSpec, std::vector<double>>> removed;
    auto invoke = [&](const CompiledSystem& c, std::vector<double> seed, const std::string& role,
                      std::vector<double>& attempt) {
        EquilibriumResult result;
        if (calls >= budget) { result.message = "declaration budget exhausted"; return result; }
        ++calls;
        core::Scope scope("phase_discovery_declaration");
        EquilibriumResult audit = structural_result(c, scope);
        bool invoked = false;
        StageWork newton_work;
        bool newton_accepted = false;
        const auto before = calls;
        std::vector<double> retry_attempt;
        try {
            if (seed.empty()) seed = boundary_start(c, true, nullptr, &newton_work);
            else seed = add_multipliers(c, seed, &audit);
            newton_work.report_observations(audit);
            if (discovery_newton_eligible(c)) {
                ++record->work["declaration_newton_attempts"];
                result = newton_attempt(c, seed, audit, nullptr, &attempt);
                newton_work.add(result);
                newton_accepted = result.success;
            }
            if (!newton_accepted) {
                invoked = true;
                audit = result = run_driver(c, seed, audit, scope, &attempt, nullptr, Attempt::Cold, {}, false, false, Retry::None);
                auto retried = fast ? result : recover_physical_amounts(c, result, result.retry == Retry::PhysicalAmounts, scope, &seed, nullptr,
                    &retry_attempt, true, &calls, budget, &invoked, Work::Invocation);
                if (calls != before) {
                    count_driver(record->work, result, "declaration_invocations");
                    StageWork(newton_work).report(result);
                    record->declarations.push_back(result); record->declaration_roles.push_back(role);
                    record->declaration_attempt_coordinates.push_back(attempt);
                    result = std::move(retried);
                    if (!retry_attempt.empty()) attempt = retry_attempt;
                }
                if (result.success && calls == before) {
                    const auto branch = check_density_branches(c, result, {}, true, scope);
                    if (!branch.accepted) {  result.message = branch.message; refuse(result); }
                }
            }
        } catch (const std::exception& error) { result = audit; result.message = error.what(); }
        record->work["iterations"] += newton_work.iterations;
        record->work["evaluations"] += newton_work.evaluations;
        record->work["hessian_calls"] += newton_work.hessian_calls;
        if (invoked) count_driver(record->work, result, "declaration_invocations");
        else if (!newton_accepted) ++record->work["declaration_initialization_failures"];
        if (!newton_accepted && calls == before) std::move(newton_work).report(result);
        record->declarations.push_back(result); record->declaration_roles.push_back(role + (calls != before ? "-physical-mol" : ""));
        record->declaration_attempt_coordinates.push_back(calls != before ? retry_attempt : attempt);
        return result;
    };
    EquilibriumResult best;
    auto best_compiled = compiled;
    for (const auto kind : {PhaseKind::Liquid, PhaseKind::Vapor}) {
        if (has_charge && kind == PhaseKind::Vapor) continue;
        declaration.phases[0].kind = kind; declaration.phases[0].name = kind == PhaseKind::Liquid ? "L1" : "V1";
        if (kind == PhaseKind::Vapor) compiled = compile_problem(model, declaration);
        std::vector<double> attempt;
        auto result = invoke(compiled, {}, "initial", attempt);
        if (result.success && (!best.success || result.objective < best.objective)) {
            best = std::move(result); best_compiled = compiled;
        }
    }
    const auto finish = [&](bool success, const std::string& message) {
        record->declaration = best_compiled.problem; record->candidate = best;
        // The accepted state also supplies the fixed declaration's ordinary cold
        // initializer, including same-kind and three-phase declarations.
        if (best.has_state()) for (std::size_t p = 0; p < best.phase_count; ++p) {
            auto& phase = record->declaration.phases[p];
            phase.composition_guess.clear();
            phase.declared_composition_guess = false;
            for (auto j : phase.support)
                phase.composition_guess.push_back(best.mole_fractions[p * feed.size() + j]);
            phase.packing_guess = best.packing_fractions[p];
        }
        report_work(record->work, work);
        EquilibriumResult result = best; result.discovery = record;
        if (fast) {
            record->termination = success ? (best.phase_count == 2 ? "admitted_two_phase" : "no_admitted_split") : "incomplete";
            result.message = success ? (best.phase_count == 2 ? "locally admitted two-phase state; post-flash stability unsearched"
                : "no admitted split in the fast trial set") : "no admitted split: " + message;
        }
        if (!success) {  if (!fast) result.message = message; refuse(result); }
        return result;
    };
    if (!best.success) {
        if (fast) record->unsearched.push_back("entire trial search (positive-neutral-base, near-pure Liquid/Vapor, ideal-gas Vapor starts)");
        return finish(false, fast && budget == 0 ? "declaration budget exhausted" : "discovery has no accepted initial candidate");
    }
    bool checking = false;  // fast-tested: one more round on the admitted endpoint with Full's trial set
    for (std::size_t round = 0; round < 2 * rank + 1; ++round) {
        auto trial_options = options;
        if (fast) trial_options.mode = checking ? DiscoveryMode::Full : DiscoveryMode::Fast;
        if (options.trial_budget) {
            const auto used = static_cast<std::size_t>(record->work["trial_invocations"]);
            trial_options.trial_budget = used < *options.trial_budget ? *options.trial_budget - used : 0;
        }
        auto report = test_snapshot(best_compiled, best, trial_options, removed, fast ? &record->unsearched : nullptr);
        for (const auto& [name, value] : report.work) record->work[name] += value;
        record->rounds.push_back(std::move(report));
        auto& current = record->rounds.back();
        if (!fast && !current.completed) return finish(false, current.message);
        std::vector<std::size_t> descending;
        for (std::size_t i = 0; i < current.trials.size(); ++i) if (current.trials[i].resolved_descent) {
            const auto& trial = current.trials[i];
            const bool duplicate = std::any_of(descending.begin(), descending.end(), [&](auto j) {
                const auto& prior = current.trials[j];
                return same_stationary_trial(prior, trial);
            });
            if (!duplicate) descending.push_back(i);
        }
        const bool unresolved = (!fast || checking) && std::any_of(current.trials.begin(), current.trials.end(), [](const auto& t) {
            return t.quantities.count("donor_slope_over_rt") && descent_resolved(t.quantities.at("tpd_over_rt"),
                t.quantities.at("donor_slope_over_rt"), t.quantities.at("descent_error_estimate_over_rt")); });
        if (checking) {  // flags only; the returned split stays the fast split
            const bool flagged = !descending.empty() || unresolved;
            auto result = finish(!flagged, "");
            if (flagged) {
                record->termination = "descent_after_split";
                result.message = descending.empty() ? "negative tangent plane with unresolved finite descent"
                    : "resolved descent at the fast admitted split; mode='full' can add the phase";
            } else if (current.completed) {
                record->unsearched.erase(std::find(record->unsearched.begin(), record->unsearched.end(), "post-flash stability"));
                result.message = "locally admitted two-phase state; no witness in the Full endpoint trial set";
            } else result.message += "; post-flash check incomplete";
            return result;
        }
        if (descending.empty()) {
            if (fast) {
                const bool available = current.completed && std::all_of(current.trials.begin(), current.trials.end(),
                    [](const auto& t) { return t.evidence == "no_resolved_descent"; });
                return finish(available, current.message);
            }
            return finish(!unresolved, unresolved ? "negative tangent plane with unresolved finite descent" : current.message);
        }
        if (best.phase_count >= rank) return finish(false, "resolved descent at balance-rank phase limit; replacement deferred");
        if (!fast) std::stable_sort(descending.begin(), descending.end(), [&](auto i, auto j) {
            return current.trials[i].quantities.at("donor_slope_over_rt") < current.trials[j].quantities.at("donor_slope_over_rt"); });
        bool added = false;
        for (auto i : descending) {
            if (calls >= budget) break;
            auto next = best_compiled.problem;
            auto phase = current.trials[i].phase;
            auto& trial = current.trials[i];
            // The optimizer may cross its starting root's label. Reuse the density owner's two admitted end
            // branches; a unique root takes the side its negative-slope sample establishes, else its source label (#423).
            try {
                const auto root = anchored_root(model, *best.temperature, *best.pressure, trial.composition,
                    trial.quantities.at("volume_m3_per_mol"), phase.kind, trial.quantities, "addition_");
                if (!root) throw std::runtime_error("descending trial has no admitted density branch");
                if (const auto side = properties::density_side(root->branch, root->molar_density, root->negative_slope_density);
                    side != DensityBranch::Unspecified)
                    phase.kind = side == DensityBranch::Liquid ? PhaseKind::Liquid : PhaseKind::Vapor;
                if (phase.kind == PhaseKind::Vapor && std::any_of(phase.support.begin(), phase.support.end(),
                    [&](auto j) { return z[j] != 0; })) throw std::runtime_error("ionic descending trial left liquid branch");
            } catch (const std::exception& error) { trial.message += "; addition: " + std::string(error.what()); continue; }
            phase.name = (phase.kind == PhaseKind::Liquid ? "L" : "V") + std::to_string(next.phases.size() + 1);
            next.phases.push_back(phase);
            // ALG-STABILITY: a trial with a reduced split adds its phase without the species its attempt left at their
            // composition bounds.
            const auto& omitted = trial.optimizer.bound_species;
            std::optional<CompiledSystem> reduced_system;
            if (!trial.reduced_split_coordinates.empty()) try {
                auto without = next;
                auto& last = without.phases.back();
                last.support.erase(std::remove_if(last.support.begin(), last.support.end(),
                    [&](auto j) { return std::count(omitted.begin(), omitted.end(), j) != 0; }), last.support.end());
                last.composition_guess.clear();
                reduced_system = compile_problem(model, without);
                if (!reduced_system->feasible) throw std::invalid_argument("rank-deficient reduced declaration");
            } catch (const std::exception& error) {
                trial.message += "; addition keeps the trial support: " + std::string(error.what());
                reduced_system.reset();
            }
            const bool omitting = reduced_system.has_value();
            auto system = omitting ? std::move(*reduced_system) : compile_problem(model, next);
            std::vector<double> attempt;
            auto result = invoke(system, omitting ? trial.reduced_split_coordinates : trial.split_coordinates, "addition", attempt);
            if (omitting && result.topology_event != TopologyEvent::None) {
                trial.message += "; addition: reduced declaration encountered a topology event";
                continue;
            }
            while (!fast && !result.success && result.topology_event != TopologyEvent::None && !attempt.empty() && calls < budget) {
                const auto before = system.prepare(attempt, true);
                const auto event = check_topology(system, attempt);
                const auto dropped = event.phase.value_or(event.pair ? event.pair->second : 0);
                removed.emplace_back(system.problem.phases[dropped], before.composition[dropped]);
                try {
                    auto [reduced, seed] = remove_discovery_phase(system, attempt, &record->work);
                    system = compile_problem(model, reduced); attempt.clear();
                    result = invoke(system, std::move(seed), "removal", attempt);
                } catch (const std::exception& error) {
                    trial.message += "; removal of " + removed.back().first.name + ": " + error.what(); break; }
            }
            if (!result.success || result.convergence.residual > row_tolerance) continue;
            // ALG-STABILITY: each omitted species' zero-row log fraction at the converged state must be finite and below
            // the composition box, the bound an unreduced addition applies at default settings; else no admission.
            if (omitting) {
                const auto box = composition_box(system);
                std::string cause;
                for (auto j : omitted) {
                    const double estimate = -system.chemical_row(*result.prepared_endpoint, system.problem.phases.size() - 1, j, 1.0);
                    const double lower = std::log(box[j]);
                    const auto species = std::to_string(j);
                    trial.quantities["addition_omitted_log_fraction_" + species] = estimate;
                    trial.quantities["addition_omitted_log_bound_" + species] = lower;
                    if (!std::isfinite(estimate)) cause = "omitted species " + species + " has no finite zero-row estimate";
                    else if (!(estimate < lower)) cause = "omitted species " + species + " is not below its composition bound";
                }
                if (!cause.empty()) { trial.message += "; addition: " + cause; continue; }
            }
            double error_estimate = 0.0;
            for (std::size_t p = 0; p < current.present_phases.size(); ++p)
                error_estimate += best.amounts[p] * energy_estimate(current.present_phases[p]);
            trial.quantities["addition_gibbs_estimate_over_rt_mol"] = error_estimate;
            if (!(result.objective < best.objective - error_estimate)) continue;
            best = std::move(result); best_compiled = std::move(system); added = true; break;
        }
        if (fast && added && (checking = options.mode == DiscoveryMode::FastTested)) continue;
        if (fast && added) return finish(true, "");
        if (!added) return finish(false, calls >= budget ? "resolved descent with exhausted declaration budget"
            : "resolved descent without an admitted lower-G declaration");
    }
    return finish(false, "search-incomplete: trial-round budget exhausted");
}

}  // namespace epcsaft::equilibrium
