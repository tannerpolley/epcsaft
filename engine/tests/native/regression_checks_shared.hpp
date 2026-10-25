#pragma once
// Includes and helpers shared by the regression_checks translation units.
// Numerical falsifiers of the regression layer: observations, fits and their exact columns
// (#46 R1, #86 density actions, #47 R2, #48 R3 and later issues), and the native rows the
// public R4 transport reproduces (#49).
#include "quad_base.hpp"
#include "regression/fit.hpp"

#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_support.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"
#include "surface/parameter_builder.hpp"
#include "properties/helmholtz_tape.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "properties/thermochemistry.hpp"
#include "binary128_support.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
namespace regression_checks {

using namespace epcsaft;
using namespace epcsaft::regression;
using equilibrium::AmountSpec;
using equilibrium::EquilibriumProblem;
using equilibrium::EquilibriumResult;
using equilibrium::FreeValue;
using equilibrium::PhaseKind;
using recording_support::field;
using recording_support::fields;
using recording_support::number;

constexpr double EPS = std::numeric_limits<double>::epsilon();
// The work budget every fit in this program runs under; FitControls has no default budget.
inline FitControls budgeted_controls() {
    FitControls controls;
    controls.maximum_iterations = 100; controls.maximum_elapsed_time_seconds = 600.0;
    return controls;
}
// #572: source-law ln K = A + d T rows (d = 0.0045/K) with fitted A (origin 0, scale 1, box [-1, 1]), output
// sigmas and the frozen Ceres controls; the caller declares input groups and any measurement matrix.
inline FitRequest intercept_request(const std::vector<double>& temperatures, const std::vector<double>& observed,
                                    const std::vector<double>& sigma) {
    FitRequest request; request.parameters = properties_support::neutral_ternary(); request.controls = budgeted_controls();
    auto& options = request.controls.solver_options.values;
    options[policy::function] = 1e-16; options[policy::gradient] = 1e-14; options[policy::parameter] = 1e-12;
    regression::Coordinate a; a.reaction = regression::ReactionCoordinate{"reaction_a", regression::ReactionCoefficient::A};
    a.origin = 0; a.scale = 1; a.lower = -1; a.upper = 1;
    request.coordinates.coordinates = {a}; request.uncertainty.emplace(); request.uncertainty->output_sigma = sigma;
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    problem.reactions[0].equilibrium_constant.reset();
    problem.reactions[0].correlation = ReactionLogPolynomial{0, 0, 0, 0.0045, 300, 250, 450, equilibrium::EOS_STANDARD_STATE_ID};
    for (std::size_t i = 0; i < temperatures.size(); ++i) {
        problem.temperature = FreeValue::fixed(temperatures[i]);
        Observation row; row.kind = ObservationKind::LogEquilibriumConstant; row.problem = problem; row.reaction = "reaction_a";
        row.observed = observed[i]; request.observations.push_back(row);
    }
    return request;
}
constexpr double FIRST_ORDER = 5.0e-6;  // P4's frozen first-order action criterion

// |model - reference| <= max(atol, rtol max(|model|, |reference|)); records and returns the verdict.
inline bool agree(const std::string& id, const char* quantity, const char* unit, double model, double reference,
           double atol, double rtol, const std::string& inputs, const char* claim, const char* group) {
    const bool verdict = std::isfinite(model) && std::isfinite(reference)
        && std::abs(model - reference) <= std::max(atol, rtol * std::max(std::abs(model), std::abs(reference)));
    recording_support::compare(id.c_str(), quantity, unit, model, reference, inputs, "Relationship", "regression",
                               claim, "symmetric_close", atol, rtol, verdict, group, id.c_str());
    return verdict;
}

inline bool reference_scale_agrees(const std::string& id, const char* quantity, const char* unit, double model,
                            double reference, double atol, double rtol, const std::string& inputs,
                            const char* claim, const char* group) {
    const bool verdict = std::isfinite(model) && std::isfinite(reference)
        && std::abs(model - reference) <= atol + rtol * std::abs(reference);
    recording_support::compare(id.c_str(), quantity, unit, model, reference, inputs, "Identity", "direct_reference",
        claim, "reference_scale", atol, rtol, verdict, group, id.c_str());
    return verdict;
}

inline bool truth(const std::string& id, const char* claim, bool ok, const std::string& inputs, const char* group) {
    recording_support::compare(id.c_str(), "verdict", "bool", ok ? 1.0 : 0.0, 1.0, inputs, "Relationship",
                               "regression", claim, "absolute", 0.0, 0.0, ok, group, id.c_str());
    return ok;
}

inline std::string rejection(const std::function<void()>& call) {
    try { call(); } catch (const std::invalid_argument& error) { return error.what(); }
    return "";
}

// ---- declared specifications

inline EquilibriumProblem pure_saturation(double temperature, double guess, bool dew = false) {
    EquilibriumProblem problem;
    problem.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, dew ? AmountSpec::pinned(0.0) : AmountSpec::free()),
                      equilibrium_support::phase("V", PhaseKind::Vapor, dew ? AmountSpec::free() : AmountSpec::pinned(0.0))};
    for (auto& phase : problem.phases) phase.composition_guess = {1.0};
    problem.temperature = FreeValue::fixed(temperature);
    problem.pressure = FreeValue::variable(guess);
    problem.feed.kind = equilibrium::FeedKind::Amounts;
    problem.feed.values = {1.0};
    problem.solver_options.values[epcsaft::policy::row] = 1.0e-10;
    return problem;
}

inline Observation solved(ObservationKind kind, const EquilibriumProblem& problem, std::size_t phase = 0,
                   std::vector<std::size_t> species = {}, bool mass = false) {
    Observation result;
    result.kind = kind;
    result.problem = problem;
    result.phase = phase;
    result.species = std::move(species);
    result.mass_basis = mass;
    result.observed = 1.0;
    result.form = ResidualForm::LogRatio;
    return result;
}

inline Observation density(double temperature, double pressure, std::vector<double> composition, DensityBranch branch,
                    bool mass = false) {
    Observation result;
    result.kind = ObservationKind::Density;
    result.temperature = temperature;
    result.pressure = pressure;
    result.composition = std::move(composition);
    result.branch = branch;
    result.mass_basis = mass;
    result.observed = 1.0;
    result.form = ResidualForm::LogRatio;
    return result;
}

// ---- R1 F1: units and basis on the retained NIST SRD 69 ethane rows (packet v1,
// source b01333e827933c0a…, parameters 7dd1899712aacfbc…, E6 probe inputs).
struct NistRow final { double T, p_source, rho_source, p_model, rho_model, rho_molar; };
const NistRow NIST[10] = {
    {100.0, 11.080787, 640.94852, 11.123831887357404, 645.7410638010867, 21474.59473897861},
    {120.0, 352.30167, 618.94997, 351.9423497001842, 622.0818292190638, 20687.78946521662},
    {140.0, 3813.5564, 596.58156, 3814.3604750250333, 598.8924184547142, 19916.60852859043},
    {160.0, 21405.224, 573.55144, 21443.595352485714, 575.8127513131142, 19149.077196977527},
    {180.0, 78638.137, 549.50874, 78838.3060997493, 552.2582967935026, 18365.75646137355},
    {200.0, 217232.94, 523.97698, 217800.6949469248, 527.4654263154658, 17541.25129083691},
    {220.0, 492046.39, 496.27145, 493283.9144780184, 500.41506993891124, 16641.671763848062},
    {240.0, 966787.79, 465.30887, 969152.4807171348, 469.5944555701788, 15616.70953010239},
    {260.0, 1711835.4, 429.07617, 1715342.640248023, 432.36903897862464, 14378.750880566167},
    {280.0, 2806735.8, 382.72712, 2807103.3101879596, 382.94892442644885, 12735.248567557328},
};

struct Ladder final { double coarse = 0.0, fine = 0.0; };
template <class Value>
Ladder ladder(const Value& value, double step) {
    Ladder result;
    result.coarse = (value(step) - value(-step)) / (2.0 * step);
    result.fine = (value(0.5 * step) - value(-0.5 * step)) / step;
    return result;
}

inline bool ladder_agrees(const std::string& id, const char* quantity, double action, const Ladder& rungs, const std::string& inputs,
                   const char* claim, const char* group, bool magnitude_scaled = false) {
    // The refined rung meets the criterion; the coarse rung may carry up to the fourfold
    // truncation of a second-order central difference.
    const double scale = magnitude_scaled ? std::abs(rungs.fine) : std::max(1.0, std::abs(rungs.fine));
    const bool fine = agree(id + "/fine", quantity, "native", action, rungs.fine, FIRST_ORDER * scale, 0.0, inputs, claim, group);
    const bool coarse = agree(id + "/coarse", quantity, "native", action, rungs.coarse, 4.0 * FIRST_ORDER * scale, 0.0, inputs, claim, group);
    return fine && coarse;
}

inline ObservationEvaluation evaluated_columns(const FitRequest& request, const std::vector<double>& v) {
    std::vector<std::size_t> rows(request.observations.size()); std::iota(rows.begin(), rows.end(), 0);
    FitResult failures; ObservationEvaluator evaluator(request, failures, rows); ObservationEvaluation out;
    CHECK(evaluator.evaluate(v, true, out));
    return out;
}

// (J^T J)_{ab} of a result's n x 2 optimizer Jacobian.
inline double gram_entry(const FitResult& result, std::size_t a, std::size_t b) {
    double gram = 0.0;
    const std::size_t rows = result.optimizer_jacobian.size() / 2;
    for (std::size_t i = 0; i < rows; ++i) gram += result.optimizer_jacobian[i * 2 + a] * result.optimizer_jacobian[i * 2 + b];
    return gram;
}

struct Truth final { Parameters parameters; std::vector<Observation> observations; };
inline Truth manufactured() {
    Truth result;
    result.parameters = properties_support::neutral_binary();
    result.parameters.dispersion_energy_over_k[0] = 152.0;
    result.parameters.k_ij[1] = result.parameters.k_ij[2] = -0.02;
    const HelmholtzModel model(result.parameters);
    for (const double T : {170.0, 180.0, 190.0}) result.observations.push_back(solved(ObservationKind::EquilibriumPressure, equilibrium_support::bubble_at(T)));
    auto y = solved(ObservationKind::PhaseMoleFraction, equilibrium_support::bubble_at(180.0), 1, {0});
    y.form = ResidualForm::Difference;
    y.scale = 0.01;
    result.observations.push_back(y);
    const double P_bubble = *equilibrium::solve_equilibrium(model, equilibrium_support::bubble_at(180.0)).pressure;
    result.observations.push_back(density(180.0, 2.0 * P_bubble, {0.5, 0.5}, DensityBranch::Liquid, true));
    for (auto& observation : result.observations) observation.observed = *predict(model, observation).value;
    return result;
}

inline std::string summary(const FitResult& result) {
    return fields({field("status", result.status), field("termination", result.termination), field("iterations", result.iterations),
                   field("residual_evaluations", result.residual_evaluations), field("jacobian_evaluations", result.jacobian_evaluations),
                   field("trial_failures", result.trial_failures), field("final_cost", result.final_cost),
                   field("covariance", result.covariance.status)});
}

}  // namespace regression_checks
