#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "contributions/helmholtz.hpp"

#include "core/trace.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>

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

// Explicit searches alone capture the shared owners' work, including initialization,
// reference roots, failed roots, branch admission and compiler evaluations.
void count_stability_work(Quantities& work, const core::Trace& trace, std::size_t begin = 0) {
    std::size_t roots = 0, density_evaluations = 0, energy_evaluations = 0;
    for (std::size_t i = begin; i < trace.spans.size(); ++i) {
        const auto& span = trace.spans[i];
        if (std::string(span.name) == "solve_density") ++roots;
        if (std::string(span.name) != "evaluate") continue;
        ++energy_evaluations;
        auto parent = span.parent;
        while (parent) {
            if (std::string(trace.spans[*parent].name) == "solve_density") { ++density_evaluations; break; }
            parent = trace.spans[*parent].parent;
        }
    }
    work["root_calls"] = roots; work["root_evaluations"] = density_evaluations;
    work["energy_evaluations"] = energy_evaluations;
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

// Reduced energy and Euler/association estimates use unit-mol amount columns and the SI volume column.
Quantities phase_energy(const HelmholtzModel& model, double T, double P, double v,
    const std::vector<double>& w, const std::vector<double>& standard, Quantities& work,
    HelmholtzDerivatives* retained = nullptr) {
    packing_admission(model.packing_fraction_coefficient(T, w) / v);
    auto d = model.evaluate(T, v, w, DerivativeOrder::IsothermalPhaseHessian);
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
    if (!model.parameters().site_component.empty()) {
        ++work["association_diagnostic_evaluations"];
        const auto diagnostics = model.density_actions(T, 1.0 / v, w,
            {PropertyObservable::ResidualEntropy}, {DensityActionDirection{0.0, 0.0,
                std::vector<double>(w.size(), 0.0), {}}}).diagnostics;
        q["association_residual"] = diagnostics.association_residual;
        q["association_site_pivot_ratio"] = diagnostics.site_pivot_ratio;
        q["association_site_rank"] = static_cast<double>(diagnostics.site_rank);
    }
    if (retained) *retained = std::move(d);
    return q;
}

double energy_estimate(const Quantities& q) {
    return q.at("association_energy_estimate_over_rt") + q.at("euler_defect_over_rt")
        + q.at("summation_estimate_over_rt");
}

std::vector<PhaseSpec> stability_supports(const CompiledSystem& c, const StabilityOptions& options) {
    const auto& z = c.model->parameters().charge_number;
    std::vector<std::size_t> united, zero_charge;
    const auto participant = reaction_participants(c.n_feed.size(), c.problem.reactions);
    for (std::size_t i = 0; i < c.n_feed.size(); ++i)
        if (c.n_feed[i] > 0.0 || participant[i]) { united.push_back(i); if (z[i] == 0.0) zero_charge.push_back(i); }
    std::vector<PhaseSpec> phases;
    if (options.trials) phases = *options.trials;
    else {
        PhaseSpec liquid; liquid.name = "trial-liquid"; liquid.support = united;
        phases.push_back(liquid);
        if (!zero_charge.empty()) {
            liquid.support = zero_charge; phases.push_back(liquid);
            liquid.name = "trial-vapor"; liquid.kind = PhaseKind::Vapor; phases.push_back(liquid);
        }
    }
    if (phases.empty()) throw std::invalid_argument("empty stability trial family");
    std::vector<PhaseSpec> admitted;
    for (auto phase : phases) {
        phase.support = resolve_support(phase, united, z);
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
    for (std::size_t k = 0; k < c.rank; ++k) rhs(k) = c.balance[k];
    const auto& z = c.model->parameters().charge_number;
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        for (std::size_t k = 0; k < c.rank; ++k) J(k, j) = c.B[k][i];
        J(c.rank + p, j) = z[i];
        rhs -= J.col(j) * at.amount[p] * at.composition[p][i];
    }
    const Eigen::VectorXd delta = J.completeOrthogonalDecomposition().solve(rhs);
    if (!delta.allFinite() || (J * delta - rhs).lpNorm<Eigen::Infinity>() > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
        throw std::invalid_argument("stability amount correction rank/consistency failure");
    std::vector<std::vector<double>> inventory(at.amount.size(), std::vector<double>(c.n_feed.size(), 0.0));
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        inventory[p][i] = at.amount[p] * at.composition[p][i] + delta(j);
        if (!(inventory[p][i] > 0.0)) throw std::invalid_argument("stability correction destroys positive support");
        const double mu = std::log(at.composition[p][i] / at.volume[p])
            + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i];
        candidate.feasibility_estimate += std::abs(delta(j) * mu);
        candidate.correction_norm += std::abs(delta(j));
    }
    report.candidate_coordinates.clear();
    for (std::size_t p = 0; p < at.amount.size(); ++p) {
        if (at.amount[p] > 0.0) {
            at.amount[p] = std::accumulate(inventory[p].begin(), inventory[p].end(), 0.0);
            for (auto i : c.problem.phases[p].support) at.composition[p][i] = inventory[p][i] / at.amount[p];
        }
        auto q = phase_energy(*c.model, at.temperature, at.pressure, at.volume[p], at.composition[p],
            report.standard_over_rt, report.work, &at.derivative[p]);
        double normalization = -1.0, charge = 0.0, potential = 0.0, residual = 0.0;
        for (auto i : c.problem.phases[p].support) {
            normalization += at.composition[p][i]; charge += z[i] * at.composition[p][i];
            potential += at.composition[p][i] * report.plane_over_rt[i];
            residual += at.composition[p][i] * (std::log(at.composition[p][i] / at.volume[p])
                + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i] - report.plane_over_rt[i]);
            report.candidate_coordinates.push_back(std::log(at.composition[p][i]));
        }
        if (std::abs(normalization) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value)
            throw std::invalid_argument("stability candidate arithmetic feasibility failure");
        report.candidate_coordinates.push_back(at.amount[p]);
        report.candidate_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature,
            at.composition[p]) / at.volume[p]));
        q["tpd_over_rt"] = q.at("g_over_rt") - potential;
        q["chemical_plus_electric_contraction"] = residual;
        q["present_identity_defect_over_rt"] = q["tpd_over_rt"] - present_phase_residual_distance(residual, q.at("pressure_row"));
        q["normalization_residual"] = normalization; q["charge_residual"] = charge;
        q["correction_l1_mol"] = candidate.correction_norm;
        q["feasibility_estimate_over_rt"] = candidate.feasibility_estimate;
        candidate.quantities.push_back(std::move(q));
    }
    report.present_phases = candidate.quantities;
    return candidate;
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

// Density-independent permittivity is a condensed-liquid correlation: a charged trial counts only on its own liquid side.
void condensed_admission(const std::optional<DensityRoot>& root, Quantities* q = nullptr) {
    if (!root || properties::density_side(root->branch, root->molar_density, root->negative_slope_density) != DensityBranch::Liquid)
        throw std::invalid_argument("charged trial outside the condensed domain of density-independent permittivity");
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

// Charge rescaling followed by a two-pivot arithmetic correction (largest entry, ties by index).
std::vector<double> reconstruct_trial_composition(std::vector<double> w,
    const std::vector<std::size_t>& support, const std::vector<int>& z) {
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
        const double ratio = std::sqrt(negative / positive);
        for (auto i : support) { if (z[i] > 0) w[i] *= ratio; if (z[i] < 0) w[i] /= ratio; }
    }
    double total = 0.0;
    for (auto i : support) total += w[i];
    if (!(std::isfinite(total) && total > 0.0)) throw std::invalid_argument("trial normalization failed");
    for (auto i : support) w[i] /= total;
    if (plus) {
        double sum = 0.0, charge = 0.0;
        for (auto i : support) { sum += w[i]; charge += z[i] * w[i]; }
        const double dp = (-charge + z[*minus] * (sum - 1.0)) / (z[*plus] - z[*minus]);
        w[*plus] += dp; w[*minus] += 1.0 - sum - dp;
    } else {
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
void donor_descent(const CompiledSystem& c, const Candidate& candidate, StabilityTrial& trial,
    StabilityReport& report, const StabilityOptions& options) {
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
    if (!descent_resolved(q.at("tpd_over_rt"), slope, uncertainty)) return;
    const auto& z = c.model->parameters().charge_number;
    double beta = capacity / 2.0;
    for (std::size_t check = 0; check <= 8; ++check, beta /= 2.0) {
        Quantities split{{"beta_mol", beta}, {"admitted", 0.0}};
        try {
            const double remaining = at.amount[p] - beta;
            std::vector<double> x(w.size(), 0.0), error(w.size(), 0.0);
            double sum = 0.0, charge = 0.0, inventory_error = 0.0;
            for (auto i : c.problem.phases[p].support) {
                const double n = at.amount[p] * at.composition[p][i] - beta * w[i];
                if (!(n > 0.0 && remaining > 0.0)) throw std::invalid_argument("split destroys donor support");
                x[i] = n / remaining; sum += x[i]; charge += z[i] * x[i];
                error[i] = remaining * x[i] + beta * w[i] - at.amount[p] * at.composition[p][i];
                inventory_error = std::max(inventory_error, std::abs(error[i]));
            }
            const auto moiety_error = matvec(c.B, error);
            split["component_error_mol"] = inventory_error;
            split["moiety_error_mol"] = maximum_abs(moiety_error);
            split["normalization_residual"] = sum - 1.0; split["charge_residual"] = charge;
            if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value
                || inventory_error > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total)
                || maximum_abs(moiety_error) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
                throw std::invalid_argument("finite split arithmetic conservation failed");
            ++report.work["finite_split_evaluations"];
            const double donor_volume = 1.0 / declared_root(*c.model, at.temperature, at.pressure, x, c.problem.phases[p].kind).molar_density;
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
                for (std::size_t phase = 0; phase < at.amount.size(); ++phase) {
                    const auto& composition = phase == p ? x : at.composition[phase];
                    for (auto i : c.problem.phases[phase].support) trial.split_coordinates.push_back(std::log(composition[i]));
                    trial.split_coordinates.push_back(phase == p ? remaining : at.amount[phase]);
                    trial.split_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature,
                        composition) / (phase == p ? donor_volume : at.volume[phase])));
                }
                for (auto i : trial.phase.support) trial.split_coordinates.push_back(std::log(w[i]));
                trial.split_coordinates.push_back(beta);
                trial.split_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature, w)
                    / q.at("volume_m3_per_mol")));
                trial.split_checks.push_back(std::move(split));
                return;
            }
        } catch (const std::exception& error) { trial.message += "; finite split: " + std::string(error.what()); }
        trial.split_checks.push_back(std::move(split));
    }
}

StabilityReport test_snapshot(const CompiledSystem& c, const EquilibriumResult& result,
    const StabilityOptions& options, const std::vector<std::pair<PhaseSpec, std::vector<double>>>& removed,
    const core::Trace& trace) {
    StabilityReport report;
    const auto begin = trace.spans.size();
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
    for (const auto& phase : supports) {
        const bool condensed = !form.is(PermittivityForm::DensityTemperature) && !form.is(PermittivityForm::DipolarPerturbation)
            && std::any_of(phase.support.begin(), phase.support.end(), [&](auto i) { return z[i] != 0; });
        std::vector<std::pair<std::string, std::vector<double>>> seeds;
        if (options.evaluate_supplied_trials) seeds.emplace_back("supplied-composition", phase.composition_guess);
        else {
            std::vector<double> base(z.size(), 0.0), remainder;
            for (auto i : phase.support) base[i] = 1.0;
            base = reconstruct_trial_composition(base, phase.support, z);
            remainder = base;
            for (std::size_t p = 0; p < at.amount.size(); ++p) {
                if (std::all_of(phase.support.begin(), phase.support.end(), [&](auto i) { return at.composition[p][i] > 0.0; })) {
                    std::vector<double> restriction(z.size(), 0.0);
                    for (auto i : phase.support) restriction[i] = at.composition[p][i];
                    remainder = reconstruct_trial_composition(std::move(restriction), phase.support, z); break;
                }
            }
            seeds.emplace_back("positive-neutral-base", base);
            if (phase.support.size() > 1) for (auto j : phase.support) if (z[j] == 0.0) {
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
        }
        std::vector<std::pair<std::vector<double>, double>> rooted;
        for (auto& [seed_name, seed] : seeds) {
            StabilityTrial trial; trial.phase = phase; trial.seed = seed_name; trial.evidence = "unresolved";
            try {
                if (options.evaluate_supplied_trials) {
                    if (seed.size() != z.size()) throw std::invalid_argument("supplied trial composition size");
                    for (std::size_t i = 0; i < seed.size(); ++i)
                        if (std::find(phase.support.begin(), phase.support.end(), i) == phase.support.end() && seed[i] != 0.0)
                            throw std::invalid_argument("supplied trial composition outside support");
                }
                seed = reconstruct_trial_composition(std::move(seed), phase.support, z);
                ++report.work["reconstruction_evaluations"];
                double eta;
                if (options.evaluate_supplied_trials && phase.packing_guess) eta = *phase.packing_guess;
                else {
                    const auto root = stability_root(*c.model, at.temperature, at.pressure, seed, phase.kind);
                    if (condensed) condensed_admission(root);
                    eta = c.model->packing_fraction_coefficient(at.temperature, seed) * root.molar_density;
                }
                packing_admission(eta);
                if (std::any_of(rooted.begin(), rooted.end(), [&](const auto& prior) {
                        return prior.first == seed && std::abs(std::log(prior.second / eta))
                            <= policy::DENSITY_ANCHOR_LOG_WINDOW.value; })) continue;
                rooted.emplace_back(seed, eta);
                if (options.trial_budget && report.work["trial_invocations"] >= *options.trial_budget) {
                    report.message = "search-incomplete: trial budget exhausted";
                    count_stability_work(report.work, trace, begin);
                    return report;
                }
                auto compiled = compile_stability_trial(*c.model, at.temperature, at.pressure, phase, seed, shift,
                    trial.row_transformation);
                trial.basis = compiled.B;
                core::Scope scope("tangent_plane_trial");
                const auto initial = trial_coordinates(compiled, seed, eta);
                ++report.work["trial_invocations"];
                if (options.evaluate_supplied_trials) trial.raw_coordinates = initial;
                else {
                    EquilibriumResult audit;
                    try {
                        trial.optimizer = run_driver(*c.model, compiled, initial, structural_result(compiled, scope),
                            scope, &trial.raw_coordinates, &audit);
                    } catch (const std::exception& error) { trial.optimizer = audit; trial.optimizer.message = error.what(); }
                    count_driver(report.work, trial.optimizer, "trial_driver_invocations");
                    trial.message = trial.optimizer.message;
                }
                if (trial.raw_coordinates.empty()) throw std::invalid_argument("no admitted nonfatal optimizer point");
                const auto raw = compiled.prepare(trial.raw_coordinates, true);
                trial.raw_residuals = compiled.evaluate_at(raw).values;
                trial.composition = reconstruct_trial_composition(raw.composition[0], phase.support, z);
                ++report.work["reconstruction_evaluations"];
                const double volume = c.model->packing_fraction_coefficient(at.temperature, trial.composition) / raw.eta[0];
                if (condensed) condensed_admission(anchored_root(*c.model, at.temperature, at.pressure, trial.composition,
                    volume, PhaseKind::Liquid, trial.quantities, "condensed_"), &trial.quantities);
                trial.quantities.merge(phase_energy(*c.model, at.temperature, at.pressure, volume, trial.composition,
                    report.standard_over_rt, report.work));
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
                    double dual = 0.0;
                    for (std::size_t k = 0; k < compiled.rank; ++k)
                        dual += trial.raw_coordinates[compiled.primal_size + k] * compiled.balance[k];
                    q["trial_dual_over_rt"] = dual; q["trial_identity_defect_over_rt"] = q["tpd_over_rt"] - dual;
                    q["raw_objective_over_rt"] = trial.optimizer.objective;
                }
                trial.evidence = "no_resolved_descent";
                if (options.evaluate_supplied_trials || q["tpd_over_rt"] < 0.0) donor_descent(c, candidate, trial, report, options);
            } catch (const std::exception& error) { trial.message += "; " + std::string(error.what()); }
            report.trials.push_back(std::move(trial));
        }
    }
    report.completed = true;
    const auto unresolved = std::count_if(report.trials.begin(), report.trials.end(), [](const auto& t) {
        return t.evidence == "unresolved" || t.evidence == "indeterminate" || t.evidence == "negative_plane_without_positive_donor"; });
    const auto descending = std::count_if(report.trials.begin(), report.trials.end(), [](const auto& t) { return t.resolved_descent; });
    report.message = descending ? "numerically resolved descending split with estimated uncertainty"
        : "no witness in the recorded trial set";
    report.message += "; unresolved trials " + std::to_string(unresolved);
    if (options.evaluate_supplied_trials) report.message = "supplied trial evaluation; unavailable trials " + std::to_string(unresolved);
    if (!options.trials && std::all_of(supports.front().support.begin(), supports.front().support.end(),
        [&](auto j) { return z[j] != 0; })) report.message += "; no neutral species";
    count_stability_work(report.work, trace, begin);
    return report;
}
} // namespace

StabilityReport stability_test(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options) {
    core::Scope scope("stability_test");
    if (!c.problem.solids.empty()) throw std::invalid_argument("stability trials are not admitted with solids");
    core::Capture capture;
    return test_snapshot(c, result, options, {}, capture.trace());
}

std::pair<EquilibriumResult, std::vector<StabilityReport>> cooling_onset(const CompiledSystem& c) {
    core::Scope scope("cooling_onset");
    const auto& request = c.problem;
    if (c.has_reaction || c.has_ions || request.discovery || request.critical || request.envelope
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
    if (!failed.success && failed.message.rfind("local phase instability:", 0) != 0 && failed.message.rfind("infeasible:", 0) != 0)
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
            core::Capture scan; Quantities scanned;
            const auto count = [&] { count_stability_work(scanned, scan.trace()); work.evaluations += scanned["root_evaluations"]; };
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
        if (failed.success && !competing(failed.temperature[0])) {
            work.report(failed);
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
        boundary.phases[vapor].composition_guess = c.feed_composition;
        recovered = solve(boundary);
        if (!recovered.success) throw std::runtime_error("boundary corrector unavailable: " + recovered.message);
        if (recovered.temperature[0] < lower - 1e-10 * upper || recovered.temperature[0] > upper + 1e-10 * upper)
            throw std::runtime_error("boundary corrector left cooling-onset bracket");
        if (competing(recovered.temperature[0])) throw std::runtime_error("competing liquid descent at corrected dew");
        work.add(failed); work.report(recovered);
        return {std::move(recovered), std::move(reports)};
    } catch (const std::exception& error) {
        work.add(recovered);
        if (failed.success) { failed.success = false; failed.message = "ordinary dew unconfirmed"; clear_state(failed); }
        failed.message += "; cooling onset refused: " + std::string(error.what());
        work.report(failed);
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
    work["energy_evaluations"] += at.amount.size();
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
            ++work["root_calls"];
            work["root_evaluations"] += root.convergence.evaluations;
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

EquilibriumResult discover_phases(const HelmholtzModel& model, const EquilibriumProblem& request) {
    if (!request.discovery || request.topology_declared || !request.phases.empty()
        || !request.temperature.is_positive() || !request.pressure.is_positive()
        || request.temperature.free || request.pressure.free)
        throw std::invalid_argument("discovery requires explicit unknown topology at fixed positive T/P");
    if (!request.solids.empty()) throw std::invalid_argument("Discover() is not admitted with solids");
    core::Capture capture;
    const double row_tolerance = request.solver_options.resolved(model.solver_options()).value(policy::row);
    const auto& options = *request.discovery;
    if (options.evaluate_supplied_trials) throw std::invalid_argument("supplied trial evaluation is not phase discovery");
    const auto feed = resolve_feed(model, request.feed);
    const auto& z = model.parameters().charge_number;
    std::vector<std::size_t> united;
    bool has_charge = false;
    const auto participant = reaction_participants(feed.size(), request.reactions);
    for (std::size_t i = 0; i < feed.size(); ++i)
        if (feed[i] > 0.0 || participant[i]) { united.push_back(i); has_charge |= z[i] != 0.0; }
    auto record = std::make_shared<DiscoveryRecord>();
    EquilibriumProblem declaration = request;
    declaration.discovery.reset(); declaration.topology_declared = true;
    PhaseSpec initial; initial.name = "L1"; initial.support = united; initial.support_declared = true;
    if (!request.reactions.empty()) {
        const auto base = reconstruct_trial_composition(std::vector<double>(feed.size(), 1.0), united, z);
        for (auto i : united) initial.composition_guess.push_back(base[i]);
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
        try {
            if (seed.empty()) seed = boundary_start(c);
            else seed = add_multipliers(c, seed);
            invoked = true;
            result = run_driver(model, c, std::move(seed), audit, scope, &attempt, &audit);
            if (result.success) {
                const auto branch = check_density_branches(model, c, result, {}, true, scope);
                if (!branch.accepted) { result.success = false; result.message = branch.message; clear_state(result); }
            }
        } catch (const std::exception& error) { result = audit; result.message = error.what(); }
        record->declarations.push_back(result); record->declaration_roles.push_back(role);
        record->declaration_attempt_coordinates.push_back(attempt);
        if (invoked) count_driver(record->work, result, "declaration_invocations");
        else ++record->work["declaration_initialization_failures"];
        return result;
    };
    EquilibriumResult best;
    auto best_compiled = compiled;
    for (const auto kind : {PhaseKind::Liquid, PhaseKind::Vapor}) {
        if (has_charge && kind == PhaseKind::Vapor) continue;
        declaration.phases[0].kind = kind; declaration.phases[0].name = kind == PhaseKind::Liquid ? "L1" : "V1";
        compiled = compile_problem(model, declaration);
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
            for (auto j : phase.support)
                phase.composition_guess.push_back(best.mole_fractions[p * feed.size() + j]);
            phase.packing_guess = best.packing_fractions[p];
        }
        count_stability_work(record->work, capture.trace());
        EquilibriumResult result = best; result.discovery = record;
        if (!success) { result.success = false; result.message = message; clear_state(result); }
        return result;
    };
    if (!best.success) return finish(false, "discovery has no accepted initial candidate");
    for (std::size_t round = 0; round < 2 * rank + 1; ++round) {
        auto trial_options = options;
        if (options.trial_budget) {
            const auto used = static_cast<std::size_t>(record->work["trial_invocations"]);
            trial_options.trial_budget = used < *options.trial_budget ? *options.trial_budget - used : 0;
        }
        auto report = test_snapshot(best_compiled, best, trial_options, removed, capture.trace());
        for (const auto& [name, value] : report.work) record->work[name] += value;
        record->rounds.push_back(std::move(report));
        auto& current = record->rounds.back();
        if (!current.completed) return finish(false, current.message);
        std::vector<std::size_t> descending;
        for (std::size_t i = 0; i < current.trials.size(); ++i) if (current.trials[i].resolved_descent) {
            const auto& trial = current.trials[i];
            const bool duplicate = std::any_of(descending.begin(), descending.end(), [&](auto j) {
                const auto& prior = current.trials[j];
                return prior.phase.support == trial.phase.support && prior.phase.kind == trial.phase.kind
                    && largest_difference(prior.composition, trial.composition) <= policy::PHASE_COINCIDENCE_DISTANCE.value
                    && std::abs(std::log(prior.quantities.at("volume_m3_per_mol") / trial.quantities.at("volume_m3_per_mol")))
                        <= policy::PHASE_COINCIDENCE_DISTANCE.value;
            });
            if (!duplicate) descending.push_back(i);
        }
        if (descending.empty()) {
            const bool unresolved = std::any_of(current.trials.begin(), current.trials.end(), [](const auto& t) {
                return t.quantities.count("donor_slope_over_rt") && descent_resolved(t.quantities.at("tpd_over_rt"),
                    t.quantities.at("donor_slope_over_rt"), t.quantities.at("descent_error_estimate_over_rt")); });
            return finish(!unresolved, unresolved ? "negative tangent plane with unresolved finite descent" : current.message);
        }
        if (best.phase_count >= rank) return finish(false, "resolved descent at balance-rank phase limit; replacement deferred");
        std::stable_sort(descending.begin(), descending.end(), [&](auto i, auto j) {
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
                const auto root = anchored_root(model, best.temperature[0], best.pressure[0], trial.composition,
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
            auto system = compile_problem(model, next);
            std::vector<double> attempt;
            auto result = invoke(system, current.trials[i].split_coordinates, "addition", attempt);
            while (!result.success && result.topology_event != TopologyEvent::None && !attempt.empty() && calls < budget) {
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
            double error_estimate = 0.0;
            for (std::size_t p = 0; p < current.present_phases.size(); ++p)
                error_estimate += best.amounts[p] * energy_estimate(current.present_phases[p]);
            trial.quantities["addition_gibbs_estimate_over_rt_mol"] = error_estimate;
            if (!(result.objective < best.objective - error_estimate)) continue;
            best = std::move(result); best_compiled = std::move(system); added = true; break;
        }
        if (!added) return finish(false, calls >= budget ? "resolved descent with exhausted declaration budget"
            : "resolved descent without an admitted lower-G declaration");
    }
    return finish(false, "search-incomplete: trial-round budget exhausted");
}

}  // namespace epcsaft::equilibrium
