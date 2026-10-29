// Reaction constants, analytic minimax, ln K, addition enthalpies and the reactive Cp fit.
#include "regression_checks_shared.hpp"
#include <cstring>

namespace regression_checks {

// #187: source-basis reaction constants through the existing solved-state first actions.
Coordinate reaction_coordinate(const char* name, ReactionCoefficient coefficient = ReactionCoefficient::A,
                               double scale = 1.0) {
    Coordinate result;
    result.reaction = ReactionCoordinate{name, coefficient};
    result.scale = scale;
    return result;
}

Observation amount_ratio(const EquilibriumProblem& problem, std::size_t numerator, std::size_t denominator,
                         std::vector<std::size_t> species, std::vector<std::size_t> denominator_species) {
    auto result = solved(ObservationKind::PhaseAmountRatio, problem, numerator, species);
    result.denominator_phase = denominator;
    result.denominator_species = std::move(denominator_species);
    return result;
}

// A state is used only after the same owner's full-rank, residual and topology gates admit a nonzero action.
void regular_reaction_state(const HelmholtzModel& model, const Observation& observation, const std::string& id) {
    const auto prediction = predict(model, observation);
    CHECK(prediction.status == ObservationStatus::Available);
    equilibrium::SolvedStateActionDirection direction;
    direction.feed_amounts_mol.assign(model.parameters().component_count(), 0.0);
    direction.active_parameters.assign(model.active_parameters().size(), 0.0);
    direction.reaction_coefficients.assign(observation.problem.reactions.size(), {});
    direction.reaction_coefficients[0][0] = 1.0;
    equilibrium::SolvedStateObservable observable;
    observable.kind = equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
    observable.component = observation.problem.phases[0].support_declared ? observation.problem.phases[0].support.front() : 0;
    const auto batch = equilibrium::solved_state_first_actions(*prediction.compiled, *prediction.solved, {observable}, {direction}).front();
    const auto& result = batch.results.front();
    truth(id + "/regular", "the central state admits a finite nonzero first reaction action",
        result.status == equilibrium::SolvedStateActionStatus::Available && result.action && std::abs(*result.action) > 0.0
        && prediction.solved->topology_event == equilibrium::TopologyEvent::None,
        fields({field("central_residual", batch.diagnostics.central_residual.value_or(std::numeric_limits<double>::quiet_NaN())),
                field("dx_dlnK", result.action.value_or(std::numeric_limits<double>::quiet_NaN())),
                field("condition", batch.diagnostics.condition_estimate.value_or(std::numeric_limits<double>::quiet_NaN()))}), "reaction-fit");
}

void reaction_data(FitRequest& request) {
    const HelmholtzModel model(request.parameters, {}, request.thermochemistry);
    for (auto& row : request.observations) {
        const auto prediction = predict(model, row);
        CHECK(prediction.status == ObservationStatus::Available);
        row.observed = *prediction.value;
    }
}

FitResult reaction_columns(const std::string& id, FitRequest request) {
    reaction_data(request);  // zero residuals keep the returned Jacobian at this exact centre
    const auto result = fit(request);
    CHECK(result.usable && result.status == "converged");
    const auto u0 = request.coordinates.start(request.parameters, request.observations);
    const auto columns = u0.size();
    for (std::size_t row = 0; row < request.observations.size(); ++row) {
        for (std::size_t c = 0; c < columns; ++c) {
            const auto value = [&](double shift) {
                auto u = u0;
                u[c] += shift;
                auto observations = request.observations;
                request.coordinates.apply(observations, u);
                const HelmholtzModel rebuilt(request.coordinates.apply(request.parameters, u));
                const auto prediction = predict(rebuilt, observations[row]);
                CHECK(prediction.status == ObservationStatus::Available);
                return std::log(*prediction.value);
            };
            const auto rungs = ladder(value, 0.005);
            const auto key = id + "/" + std::to_string(row) + "/" + std::to_string(c);
            const auto inputs = fields({field("fd_signal", std::abs(rungs.fine)), field("u_step", 0.005),
                                      field("prediction", result.predictions[row])});
            truth(key + "/signal", "the dimensionless finite-difference column has signal at least 1e-2",
                  std::abs(rungs.fine) >= 1.0e-2, inputs, "reaction-fit");
            ladder_agrees(key, "dlny_du", result.optimizer_jacobian[row * columns + c], rungs, inputs,
                         "the exact log-residual column agrees with independently rebuilt central differences", "reaction-fit");
        }
    }
    return result;
}

EquilibriumProblem reacting_bubble(const HelmholtzModel& model) {
    auto problem = equilibrium_support::bubble_at(180.0);
    const auto center = equilibrium::solve_equilibrium(model, problem);
    CHECK(center.success && center.topology_event == equilibrium::TopologyEvent::None);
    const auto compiled = equilibrium::compile_problem(model, problem);
    equilibrium::SolvedStateObservable first, second;
    first.kind = second.kind = equilibrium::SolvedStateObservableKind::PhaseLogActivity;
    second.component = 1;
    const auto batch = equilibrium::solved_state_actions(compiled, center,
        {{first, second}, {{0.0, 0.0, {0.0, 0.0}, {}}}});
    CHECK(batch.results[0].value && batch.results[1].value);
    const double L = *batch.results[1].value - *batch.results[0].value;
    problem.reactions = {{"iso", {-1.0, 1.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{L - 300.0 / 180.0 - 0.5 * std::log(180.0 / 200.0) - 0.001 * 180.0,
                            300.0, 0.5, 0.001, 200.0, 150.0, 220.0, equilibrium::EOS_STANDARD_STATE_ID}}};
    const auto reacting = equilibrium::solve_equilibrium(model, problem);
    CHECK(reacting.success && reacting.topology_event == equilibrium::TopologyEvent::None);
    agree("reaction-a3/pinned_pressure", "P", "Pa", *reacting.pressure, *center.pressure, 0.0, 1.0e-10,
          "T=180 K", "the reacting manufactured bubble retains the original state", "reaction-fit");
    for (std::size_t i = 0; i < center.mole_fractions.size(); ++i)
        agree("reaction-a3/pinned_x/" + std::to_string(i), "x", "1", reacting.mole_fractions[i], center.mole_fractions[i],
              1.0e-10, 0.0, "T=180 K", "the reacting manufactured bubble retains the original compositions", "reaction-fit");
    return problem;
}

EquilibriumProblem referenced_formula(std::optional<double> source_pressure) {
    auto problem = equilibrium_support::flash(3);
    problem.temperature = FreeValue::fixed(298.15);
    problem.feed.values = {1.0, 1.0e-6, 1.0e-6};
    problem.source_reference = NeutralReferenceDeclaration{0, {0.0, 1.0, 1.0}, 1.0};
    problem.reactions = {{"H2O=Na++Cl-", {-1.0, 1.0, 1.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{std::log(1.0e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0, "analytic-source"},
        ReactionReference{source_pressure, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
    return problem;
}

struct AnalyticMinimaxProblem final { FitRequest request; std::vector<CheckResult> frozen; };

HydrationHold analytic_minimax_hold(double upper) {
    HydrationHold hold;
    hold.reference.kind = ObservationKind::SolvationGibbsEnergy;
    hold.reference.temperature = 298.15; hold.reference.pressure = 1.0e5;
    hold.reference.composition = {1.0}; hold.reference.ion_indices = {0, 1};
    hold.reference.stoichiometry = {1, 1}; hold.reference.branch = DensityBranch::Liquid;
    hold.reference.species = {0}; hold.reference.form = ResidualForm::Difference;
    hold.lower = 0.0; hold.upper = upper;
    hold.energy_tolerance = 1.0e-6; hold.diameter_tolerance = 1.0e-10; hold.maximum_iterations = 100;
    return hold;
}

AnalyticMinimaxProblem analytic_minimax_problem(bool b1, double q, double sigma) {
    AnalyticMinimaxProblem problem;
    auto& request = problem.request;
    request.coordinates.coordinates.push_back({{{ParameterFamily::Kij, {1, 2}}}, 0.0, 1.0, 0.0, b1 ? 2.0 : 3.0});
    if (b1) request.coordinates.coordinates.push_back({{{ParameterFamily::SegmentCount, {0}}}, 0.0, 1.0, 0.0, 2.0});
    request.observations.resize(5);
    for (auto& row : request.observations) {
        row.kind = ObservationKind::MeanIonicActivityCoefficient;
        row.form = ResidualForm::LogRatio; row.observed = 1.0; row.scale = 1.0;
    }
    request.weights.assign(5, 1.0);
    request.holds = {analytic_minimax_hold(b1 ? 4.0 : 6.0)};
    request.controls.maximum_iterations = 100; request.controls.maximum_elapsed_time_seconds = 5.0;
    request.minimax = MinimaxRequest{};
    request.minimax->targets = b1 ? std::vector<std::size_t>{0, 1} : std::vector<std::size_t>{0, 1, 2, 3, 4};
    request.minimax->group_scales = b1 ? std::vector<double>{sigma, sigma} : std::vector<double>{};
    request.minimax->q = q; request.minimax->tolerance = 1.0e-10;
    request.minimax->activity_distance = 1.0e-6; request.minimax->margin = 0.1;
    if (b1) {
        CheckGroup maximum; maximum.name = "maximum"; maximum.rows = {2};
        maximum.metric = CheckMetric::AbsoluteLog10Ratio; maximum.maximum = 0.0075;
        request.checks.push_back(maximum);
        CheckResult frozen_maximum; frozen_maximum.name = maximum.name; frozen_maximum.unit = "decade";
        frozen_maximum.status = "pass"; frozen_maximum.rows = maximum.rows;
        frozen_maximum.maximum_threshold = maximum.maximum; problem.frozen.push_back(frozen_maximum);
        CheckGroup mean; mean.name = "mean"; mean.rows = {3, 4};
        mean.metric = CheckMetric::AbsoluteLog10Ratio; mean.mean = 0.0075;
        request.checks.push_back(mean);
        CheckResult frozen_mean; frozen_mean.name = mean.name; frozen_mean.unit = "decade";
        frozen_mean.status = "pass"; frozen_mean.rows = mean.rows;
        frozen_mean.mean_threshold = mean.mean; problem.frozen.push_back(frozen_mean);
    }
    return problem;
}

ObservationFunction analytic_minimax_evaluator(const FitRequest& request) {
    return [&request](const std::vector<double>& u, bool, ObservationEvaluation& out) {
        constexpr double a = 0.01;
        const double ln10 = std::log(10.0);
        const bool b1 = u.size() == 2;
        const std::size_t physical = b1 ? 3 : 2, n = request.observations.size();
        const double x = u[0], z = b1 ? u[1] : 0.0, d = 2.0 * x;
        out = ObservationEvaluation{};
        out.predictions.resize(n); out.residuals.resize(n); out.weighted.resize(n);
        out.statuses.assign(n, ObservationStatus::Available);
        out.gradient.assign(n * physical, 0.0); out.optimizer_gradient.assign(n * u.size(), 0.0);
        out.chart.parameters = request.parameters; out.chart.active = request.coordinates.active_parameters();
        out.chart.active.push_back({ParameterFamily::BornDiameter, {0}});
        out.chart.physical = b1 ? std::vector<double>{x, z, d} : std::vector<double>{x, d};
        out.chart.map = b1 ? std::vector<double>{1.0, 0.0, 0.0, 1.0, 2.0, 0.0}
                           : std::vector<double>{1.0, 2.0};
        HoldResult hold; hold.derivative = 1.0;
        hold.parameter_derivative = b1 ? std::vector<double>{-2.0, 0.0, 1.0} : std::vector<double>{-2.0, 1.0};
        // A bracket-limited hold: its achieved uncertainty is below its last Newton correction.
        hold.correction = 2.0e-3; hold.uncertainty = 1.0e-3;
        hold.iterations = 1; out.chart.holds = {hold};
        for (std::size_t i = 0; i < n; ++i) {
            double e = 0.0; std::vector<double> de(physical, 0.0);
            if (b1 && i == 0) { e = a * (x + z); de = {a, a, 0.0}; }
            else if (b1 && i == 1) { e = a * (3.0 - d / 2.0 - z); de = {0.0, -a, -a / 2.0}; }
            else if (b1 && i == 2) { e = a * x; de = {a, 0.0, 0.0}; }
            else if (b1 && i == 3) { e = a * z; de = {0.0, a, 0.0}; }
            else if (b1) { e = a * (z + 0.5); de = {0.0, a, 0.0}; }
            else if (i < 4) { e = a * x; de = {a, 0.0}; }
            else { e = a * (d / 2.0 - 3.0); de = {0.0, a / 2.0}; }
            const double y = std::pow(10.0, e), weight = request.weights.empty() ? 1.0 : request.weights[i];
            out.predictions[i] = y; out.residuals[i] = std::log(y / request.observations[i].observed);
            out.weighted[i] = std::sqrt(weight) * out.residuals[i]; out.cost += 0.5 * out.weighted[i] * out.weighted[i];
            for (std::size_t p = 0; p < physical; ++p) {
                out.gradient[i * physical + p] = ln10 * y * de[p];
                for (std::size_t c = 0; c < u.size(); ++c)
                    out.optimizer_gradient[i * u.size() + c] += out.gradient[i * physical + p] * out.chart.map[p * u.size() + c];
            }
        }
        return true;
    };
}

void check_minimax_kkt(const std::string& id, const MinimaxResult& m, const std::string& inputs) {
    agree(id + "/primal", "scaled_primal_norm", "1", m.primal, 0.0, 1.0e-8, 0.0,
        inputs, "scaled primal feasibility is within tolerance", "r2-minimax");
    agree(id + "/stationarity", "scaled_stationarity_norm", "1", m.stationarity, 0.0, 1.0e-8, 0.0,
        inputs, "scaled stationarity is within tolerance", "r2-minimax");
    agree(id + "/dual_sign", "scaled_dual_sign_norm", "1", m.dual_sign, 0.0, 1.0e-8, 0.0,
        inputs, "scaled dual signs are feasible", "r2-minimax");
    agree(id + "/complementarity", "scaled_complementarity_norm", "1", m.complementarity, 0.0,
        1.0e-8, 0.0, inputs, "scaled complementarity is within tolerance", "r2-minimax");
}

void check_b1_minimax(const std::string& id, const FitRequest& request, const FitResult& result,
                      const std::vector<double>& start) {
    const auto inputs = fields({field("q", request.minimax->q), field("sigma", request.minimax->group_scales[0]),
        field("start_x", start[0]), field("start_z", start[1]), field("a_per_decade", 0.01)});
    CHECK(result.minimax.has_value());
    const auto& m = *result.minimax;
    truth(id + "/available", "the analytic constrained candidate and KKT diagnostics are available",
  result.usable && m.availability == "available" && m.feasible && m.stationary, inputs, "r2-minimax");
    const std::vector<double> coordinates = {0.75, 0.5, 1.5}, errors = {0.0125, 0.0175, 0.0075, 0.005, 0.010};
    for (std::size_t i = 0; i < coordinates.size(); ++i)
        agree(id + "/coordinate/" + std::to_string(i), "physical_coordinate", "1", result.chart.physical[i],
      coordinates[i], 1.0e-7, 0.0, inputs, "the unique guarded minimax point is recovered", "r2-minimax");
    agree(id + "/H", "H", "1", result.chart.physical[2] - 2.0 * result.chart.physical[0], 0.0,
  1.0e-12, 0.0, inputs, "the held equality is eliminated", "r2-minimax");
    for (std::size_t i = 0; i < errors.size(); ++i)
        agree(id + "/error/" + std::to_string(i), "log10_prediction_ratio", "decade",
      check_error(CheckMetric::AbsoluteLog10Ratio, request.observations[i], result.predictions[i]).first,
      errors[i], 1.0e-9, 0.0, inputs, "the analytic observation error is recovered", "r2-minimax");
    agree(id + "/t", "minimax_t", "decade", m.t, 0.0175, 1.0e-9, 0.0,
  inputs, "the active target sets the epigraph", "r2-minimax");
    agree(id + "/target_maximum", "target_maximum", "decade", m.target_maximum, 0.0175, 1.0e-9, 0.0,
  inputs, "the reported maximum matches the active target", "r2-minimax");
    agree(id + "/gap", "epigraph_gap", "decade", m.gap, 0.0, 1.0e-9, 0.0,
  inputs, "the epigraph is tight", "r2-minimax");
    CHECK(m.rows.size() == m.physical_values.size() && m.rows.size() == m.physical_multipliers.size());
    double mu_e2 = 0.0, mu_max = 0.0, mu_e4 = 0.0, mu_e5 = 0.0;
    for (std::size_t i = 0; i < m.rows.size(); ++i) {
        const auto& r = m.rows[i]; double value = 0.0, multiplier = 0.0;
        if (r.kind == "target") {
            value = r.sign * errors[r.row] - m.t;
            if (r.row == 1 && r.sign == 1.0) { multiplier = 1.0; mu_e2 = m.physical_multipliers[i]; }
        } else if (r.kind == "maximum") {
            value = r.sign * errors[r.row] - 0.0075;
            if (r.group == 0 && r.sign == 1.0) { multiplier = 1.0; mu_max = m.physical_multipliers[i]; }
        } else if (r.kind == "link") {
            const double slack = errors[r.row]; value = r.sign * errors[r.row] - slack;
            if (r.sign == 1.0) {
                multiplier = 0.5;
                if (r.row == 3) mu_e4 = m.physical_multipliers[i]; else mu_e5 = m.physical_multipliers[i];
                CHECK(r.auxiliary.has_value());
                agree(id + "/slack/" + std::to_string(r.row), "mean_slack", "decade",
              m.point.at(*r.auxiliary), slack, 1.0e-9, 0.0, inputs,
              "each mean slack equals its absolute row error", "r2-minimax");
            }
        } else { value = 0.0; multiplier = 1.0; }
        agree(id + "/row/" + std::to_string(i), "physical_constraint", "decade", m.physical_values[i],
      value, 1.0e-9, 0.0, inputs, "the physical target and guard rows match the analytic values", "r2-minimax");
        agree(id + "/multiplier/" + std::to_string(i), "physical_row_multiplier", "1", m.physical_multipliers[i],
      multiplier, 1.0e-6, 0.0, inputs, "physical row multipliers match the analytic active set", "r2-minimax");
        // Only e2 depends on d (de2/dd = -a/2): its scaled root error is (a/2)/scale times the hold uncertainty.
        agree(id + "/root_error/" + std::to_string(i), "scaled_root_error", "1", m.root_errors.at(i),
      r.members.empty() && r.row == 1 ? 0.005 / r.scale * 1.0e-3 : 0.0, 0.0, 1.0e-12, inputs,
      "the root error counts the achieved hold uncertainty once", "r2-minimax");
    }
    for (const double z : m.physical_z_lower)
        agree(id + "/lower_multiplier", "physical_bound_multiplier", "1", z, 0.0, 1.0e-6, 0.0,
      inputs, "the interior point has no lower-bound multiplier", "r2-minimax");
    for (const double z : m.physical_z_upper)
        agree(id + "/upper_multiplier", "physical_bound_multiplier", "1", z, 0.0, 1.0e-6, 0.0,
      inputs, "the interior point has no upper-bound multiplier", "r2-minimax");
    agree(id + "/nu", "equality_hold_multiplier", "decade_per_H", m.hold_multipliers.at(0), 0.005,
  1.0e-6, 0.0, inputs, "the eliminated equality reconstructs nu=a/2", "r2-minimax");
    CHECK(m.full_stationarity.size() == 3);
    for (std::size_t i = 0; i < m.full_stationarity.size(); ++i)
        agree(id + "/full_stationarity/" + std::to_string(i), "full_space_stationarity", "decade_per_coordinate",
      m.full_stationarity[i], 0.0, 1.0e-8, 0.0, inputs, "full-space stationarity vanishes", "r2-minimax");
    agree(id + "/stationarity_x", "stationarity_x", "decade_per_x", 0.01 * mu_max - 2.0 * m.hold_multipliers[0],
  0.0, 1.0e-8, 0.0, inputs, "a*mu_max-2*nu=0", "r2-minimax");
    agree(id + "/stationarity_z", "stationarity_z", "decade_per_z", -0.01 * mu_e2 + 0.01 * (mu_e4 + mu_e5),
  0.0, 1.0e-8, 0.0, inputs, "-a*mu_e2+a*(mu_e4+mu_e5)=0", "r2-minimax");
    agree(id + "/stationarity_d", "stationarity_d", "decade_per_d", -0.01 * mu_e2 / 2.0 + m.hold_multipliers[0],
  0.0, 1.0e-8, 0.0, inputs, "-a*mu_e2/2+nu=0", "r2-minimax");
    CHECK(result.chart.map == std::vector<double>({1.0, 0.0, 0.0, 1.0, 2.0, 0.0}));
    CHECK(result.chart.holds[0].parameter_derivative == std::vector<double>({-2.0, 0.0, 1.0}));
    agree(id + "/H_d", "H_d", "1", result.chart.holds[0].derivative, 1.0, 0.0, 0.0,
  inputs, "the analytic hold has H_d=1", "r2-minimax");
    agree(id + "/hold_residual", "H", "1", result.chart.holds[0].energy_error, 0.0, 1.0e-12, 0.0,
  inputs, "the analytic equality has zero residual", "r2-minimax");
    check_minimax_kkt(id, m, inputs);
}

void check_d7_minimax() {
    for (bool typed : {false, true}) {
        const auto p = analytic_minimax_problem(false, 0.02, 0.02);
        const auto real = analytic_minimax_evaluator(p.request);
        int calls = 0;
        core::Convergence record; record.stop = core::Stop::left_domain; record.residual = 0.25;
        record.criterion = &policy::DENSITY_ROOT_RESIDUAL; record.cause = "R2-D7 observation action refused";
        const ObservationFunction refused = [&](const auto& x, bool columns, auto& out) {
            if (calls++ == 0) return real(x, columns, out);
            if (typed) throw core::ConvergenceFailure(record);
            return false;
        };
        const auto result = solve_minimax(p.request, {0.5}, p.frozen, refused);
        CHECK(result.first_numerical_failure && result.minimax && result.minimax->ipopt_diagnostics);
        CHECK(result.minimax->ipopt_diagnostics->initialization_return == int(Ipopt::Solve_Succeeded));
        CHECK(result.minimax->ipopt_diagnostics->application_return.has_value());
        if (typed) {
            CHECK(result.first_numerical_failure->stop == record.stop && result.first_numerical_failure->residual == record.residual);
            CHECK(result.first_numerical_failure->criterion == record.criterion && result.first_numerical_failure->cause == record.cause);
        } else CHECK(result.first_numerical_failure->cause.find("observation callback returned false") != std::string::npos);
    }
    std::optional<equilibrium::IpoptDiagnostics> observed(std::in_place);
    for (int mode : {0, 1, 0}) equilibrium::observe_ipopt(*observed, nullptr, {mode, 7, 3, 1, 2, 3, 4, 5, 6, 7, 8});
    CHECK(observed->callbacks == 3 && observed->restoration_entries == 1 && observed->first_restoration_iteration == 7);
    CHECK(observed->last->mode == 0 && observed->last->alpha_pr == 8 && observed->last->line_search_trials == 3);
    for (const auto& start : std::vector<std::vector<double>>{{0.25, 0.25}, {1.0, 1.0}}) {
        const auto problem = analytic_minimax_problem(true, 0.02, 0.02);
        const auto evaluate = analytic_minimax_evaluator(problem.request);
        const auto result = solve_minimax(problem.request, start, problem.frozen, evaluate);
        check_b1_minimax("r2-d7-b1/start-" + number(start[0]), problem.request, result, start);
    }
    {
        const std::vector<double> start = {0.25, 0.25};
        const auto problem = analytic_minimax_problem(true, 0.01, 0.04);
        const auto evaluate = analytic_minimax_evaluator(problem.request);
        const auto result = solve_minimax(problem.request, start, problem.frozen, evaluate);
        check_b1_minimax("r2-d7-b1/rescaled", problem.request, result, start);
    }
    const auto problem = analytic_minimax_problem(false, 0.02, 0.02);
    const auto evaluate = analytic_minimax_evaluator(problem.request);
    ObservationEvaluation ls;
    const std::vector<double> ls_point = {3.0 / 5.0};
    CHECK(evaluate(ls_point, true, ls));
    double ls_max = 0.0, ls_gradient = 0.0;
    for (std::size_t i = 0; i < problem.request.observations.size(); ++i) {
        const auto [error, slope] = check_error(CheckMetric::AbsoluteLog10Ratio, problem.request.observations[i], ls.predictions[i]);
        ls_max = std::max(ls_max, std::abs(error));
        ls_gradient += error * slope * ls.optimizer_gradient[i];
    }
    const auto ls_inputs = fields({field("a_per_decade", 0.01), field("equal_LogRatio_weights", 1.0)});
    agree("r2-d7-b2/least_squares_x", "x", "1", ls.chart.physical[0], 0.6, 1.0e-12, 0.0,
        ls_inputs, "four ax rows and one a(x-3) row have their analytic least-squares minimum", "r2-minimax");
    agree("r2-d7-b2/least_squares_d", "d", "1", ls.chart.physical[1], 1.2, 1.0e-12, 0.0,
        ls_inputs, "the least-squares point obeys d=2x", "r2-minimax");
    agree("r2-d7-b2/least_squares_gradient", "least_squares_gradient", "decade_squared_per_x", ls_gradient, 0.0,
        1.0e-14, 0.0, ls_inputs, "4x+(x-3)=0 is the equal-weight least-squares condition", "r2-minimax");
    agree("r2-d7-b2/least_squares_maximum", "maximum_target_error", "decade", ls_max, 0.024, 1.0e-12, 0.0,
        ls_inputs, "the least-squares maximum is 0.024 decade", "r2-minimax");
    truth("r2-d7-b2/least_squares_guard", "the least-squares point misses the 0.02 decade maximum guard",
        ls_max > 0.02, fields({ls_inputs, field("guard_decade", 0.02)}), "r2-minimax");
    for (const double x : {0.5, 2.5}) {
        const std::vector<double> start = {x};
        const auto result = solve_minimax(problem.request, start, problem.frozen, evaluate);
        const std::string id = "r2-d7-b2/start-" + number(x);
        const auto inputs = fields({field("q", problem.request.minimax->q), field("start_x", x), field("a_per_decade", 0.01)});
        CHECK(result.minimax.has_value());
        const auto& m = *result.minimax;
        truth(id + "/available", "the analytic minimax candidate and KKT diagnostics are available",
            result.usable && m.availability == "available" && m.feasible && m.stationary, inputs, "r2-minimax");
        agree(id + "/x", "x", "1", result.chart.physical[0], 1.5, 1.0e-7, 0.0,
            inputs, "the opposing maximum errors balance at x=1.5", "r2-minimax");
        agree(id + "/d", "d", "1", result.chart.physical[1], 3.0, 1.0e-7, 0.0,
            inputs, "the eliminated equality gives d=2x", "r2-minimax");
        agree(id + "/H", "H", "1", result.chart.physical[1] - 2.0 * result.chart.physical[0], 0.0,
            1.0e-12, 0.0, inputs, "the held equality is eliminated", "r2-minimax");
        agree(id + "/t", "minimax_t", "decade", m.t, 0.015, 1.0e-9, 0.0,
            inputs, "the epigraph is the maximum target error", "r2-minimax");
        agree(id + "/target_maximum", "target_maximum", "decade", m.target_maximum, 0.015, 1.0e-9, 0.0,
            inputs, "the reported maximum matches the analytic target", "r2-minimax");
        agree(id + "/gap", "epigraph_gap", "decade", m.gap, 0.0, 1.0e-9, 0.0,
            inputs, "the epigraph is tight", "r2-minimax");
        CHECK(m.rows.size() == m.physical_values.size());
        for (std::size_t i = 0; i < m.rows.size(); ++i) {
            const auto& row = m.rows[i];
            const double e = check_error(row.metric, problem.request.observations[row.row], result.predictions[row.row]).first;
            agree(id + "/row/" + std::to_string(i), "physical_constraint", "decade", m.physical_values[i],
                row.sign * e - m.t, 1.0e-9, 0.0, inputs,
                "the physical max rows retain their analytic errors", "r2-minimax");
        }
        check_minimax_kkt(id, m, inputs);
    }
}

void check_reaction_derivatives() {
    FitRequest a1;
    a1.controls = budgeted_controls();
    a1.parameters = equilibrium_support::synthetic_reactive_ternary();
    auto liquids = equilibrium_support::synthetic_reactive_problem();
    liquids.standard_gibbs_j_per_mol.clear();
    a1.observations = {solved(ObservationKind::PhaseMoleFraction, liquids, 0, {2}),
                       solved(ObservationKind::PhaseFraction, liquids, 0), amount_ratio(liquids, 0, 1, {2}, {2})};
    a1.coordinates.coordinates = {reaction_coordinate("A+B=C")};
    regular_reaction_state(HelmholtzModel(a1.parameters), a1.observations[0], "reaction-a1");
    reaction_columns("reaction-a1/scalar", a1);
    a1.coordinates.coordinates = {{{{ParameterFamily::Kij, {2, 6}}}, 0.045, 0.05}, reaction_coordinate("A+B=C"),
                                 {{{ParameterFamily::SigmaConstant, {1}}}, 3.0, 0.5}};
    const auto interleaved = reaction_columns("reaction-a1/joint", a1);
    const auto& j = interleaved.physical_jacobian;
    for (std::size_t row = 0; row < a1.observations.size(); ++row) {
        const double slope = 1.0 / interleaved.predictions[row];
        for (std::size_t c = 0; c < 3; ++c) {
            const std::size_t physical = c == 1 ? 2 : c == 2 ? 1 : 0;
            agree("reaction-a1/order/" + std::to_string(row) + "/" + std::to_string(c), "dlny_du", "1",
                interleaved.optimizer_jacobian[row * 3 + c], slope * j[row * 3 + physical] * a1.coordinates.coordinates[c].scale, 0.0, 1.0e-12, "",
                "optimizer declaration order and EOS-first physical-row order carry the same derivatives", "reaction-fit");
        }
    }
    for (const auto source_pressure : {std::optional<double>{}, std::optional<double>{2.0e5}}) {
        FitRequest a2;
        a2.controls = budgeted_controls();
        a2.parameters = properties_support::aqueous_electrolyte();
        a2.observations = {solved(ObservationKind::PhaseMoleFraction, referenced_formula(source_pressure), 0, {1})};
        a2.coordinates.coordinates = {reaction_coordinate("H2O=Na++Cl-"),
                                     reaction_coordinate("H2O=Na++Cl-", ReactionCoefficient::B, 298.15)};
        const auto direct = predict(HelmholtzModel(a2.parameters), a2.observations.front());
        regular_reaction_state(HelmholtzModel(a2.parameters), a2.observations.front(), source_pressure ? "reaction-a2/fixed" : "reaction-a2/system");
        const auto result = reaction_columns(source_pressure ? "reaction-a2/fixed" : "reaction-a2/system", a2);
        agree(source_pressure ? "reaction-a2/fixed_basis" : "reaction-a2/system_basis", "x_Na", "1",
              result.predictions[0], *direct.value, 0.0, 1.0e-12, "",
              "rewriting the native coefficient retains the declared source reference", "reaction-fit");
        a2.coordinates.coordinates[0].reaction->coefficient = ReactionCoefficient::EnthalpyShift;
        a2.coordinates.coordinates[0].reaction->pivot_temperature = 313.15;
        a2.coordinates.coordinates[0].scale = 2500.0;
        reaction_columns(source_pressure ? "reaction-pair/fixed" : "reaction-pair/system", a2);
    }
    FitRequest a3;
    a3.controls = budgeted_controls();
    a3.parameters = properties_support::neutral_binary();
    const auto bubble = reacting_bubble(HelmholtzModel(a3.parameters));
    a3.observations = {solved(ObservationKind::EquilibriumPressure, bubble), solved(ObservationKind::PartialPressure, bubble, 1, {0})};
    a3.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::B, 180.0),
        reaction_coordinate("iso", ReactionCoefficient::C), reaction_coordinate("iso", ReactionCoefficient::D, 1.0 / 180.0)};
    regular_reaction_state(HelmholtzModel(a3.parameters), a3.observations.front(), "reaction-a3");
    const auto result = reaction_columns("reaction-a3/free_pressure", a3);
    auto pair = a3;
    pair.coordinates.coordinates.resize(2);
    pair.coordinates.coordinates[1].reaction = ReactionCoordinate{"iso", ReactionCoefficient::EnthalpyShift, 200.0};
    pair.coordinates.coordinates[1].scale = 2500.0;
    reaction_columns("reaction-pair/free_pressure", pair);
    const double basis[] = {1.0, 1.0 / 180.0, std::log(180.0 / 200.0), 180.0};
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t c = 1; c < 4; ++c)
            agree("reaction-a3/basis/" + std::to_string(row) + "/" + std::to_string(c), "dy_dtheta", "native",
                  result.physical_jacobian[row * 4 + c], result.physical_jacobian[row * 4] * basis[c], 0.0, 1.0e-12,
                  "T=180 K; Tref=200 K", "physical columns obey the source ln K basis vector", "reaction-fit");
    auto free_t = bubble;
    free_t.pressure = FreeValue::fixed(result.predictions.front());
    free_t.temperature = FreeValue::variable(180.0);
    const HelmholtzModel model(a3.parameters);
    const auto j1 = equilibrium::solve_equilibrium(model, bubble);
    for (std::size_t q = 0; q < free_t.phases.size(); ++q) {
        free_t.phases[q].composition_guess = {j1.mole_fractions[q * 2], j1.mole_fractions[q * 2 + 1]};
        free_t.phases[q].packing_guess = j1.packing_fractions[q];
    }
    const auto recovered = equilibrium::solve_equilibrium(model, free_t);
    CHECK(recovered.success && recovered.topology_event == equilibrium::TopologyEvent::None);
    agree("reaction-a4/pinned_temperature", "T", "K", *recovered.temperature, 180.0, 0.0, 1.0e-10,
          "declared J1 compositions and packing fractions", "the free-temperature solve returns the specified J1 state", "reaction-fit");
    agree("reaction-a4/pinned_pressure", "P", "Pa", *recovered.pressure, *j1.pressure, 0.0, 1.0e-10,
          "", "the free-temperature solve retains the J1 pressure", "reaction-fit");
    for (std::size_t i = 0; i < j1.mole_fractions.size(); ++i)
        agree("reaction-a4/pinned_x/" + std::to_string(i), "x", "1", recovered.mole_fractions[i], j1.mole_fractions[i],
              1.0e-10, 0.0, "", "the free-temperature solve retains the J1 compositions", "reaction-fit");
    for (std::size_t q = 0; q < j1.packing_fractions.size(); ++q)
        agree("reaction-a4/pinned_packing/" + std::to_string(q), "eta", "1", recovered.packing_fractions[q], j1.packing_fractions[q],
              0.0, 1.0e-10, "", "the free-temperature solve retains the J1 packing fractions", "reaction-fit");
    a3.observations = {solved(ObservationKind::EquilibriumTemperature, free_t)};
    a3.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::D, 1.0 / 180.0)};
    regular_reaction_state(HelmholtzModel(a3.parameters), a3.observations.front(), "reaction-a4");
    reaction_columns("reaction-a4/free_temperature", a3);
}

void recover_reactions(const std::string& id, const FitRequest& request, const std::vector<double>& expected) {
    const auto result = fit(request);
    truth(id + "/converged", "the manufactured fit converges without failed trials",
          result.status == "converged" && result.usable, summary(result), "reaction-fit");
    for (std::size_t row = 0; row < expected.size(); ++row) {
        const bool b = request.coordinates.coordinates[request.coordinates.tie()[row]].reaction
            && request.coordinates.coordinates[request.coordinates.tie()[row]].reaction->coefficient == ReactionCoefficient::B;
        agree(id + "/value/" + std::to_string(row), "theta", b ? "K" : "1", result.physical[row], expected[row],
              b ? 0.0 : 1.0e-6, b ? 1.0e-6 : 0.0, summary(result), "the known physical coefficient is recovered", "reaction-fit");
    }
    agree(id + "/cost", "cost", "1", result.final_cost, 0.0, 1.0e-14, 0.0, summary(result),
          "the noise-free final objective vanishes", "reaction-fit");
    truth(id + "/covariance", "the interior fit has full rank and an available conditional covariance",
          result.covariance.status == "available" && result.covariance.rank == expected.size(), summary(result), "reaction-fit");
}

void check_reaction_recovery() {
    FitRequest b1;
    b1.controls = budgeted_controls();
    b1.parameters = equilibrium_support::synthetic_reactive_ternary();
    auto problem = equilibrium_support::synthetic_reactive_problem();
    problem.standard_gibbs_j_per_mol.clear();
    problem.reactions[0].equilibrium_constant = std::exp(0.1);
    b1.parameters.k_ij[2] = b1.parameters.k_ij[6] = 0.05;
    b1.observations = {solved(ObservationKind::PhaseMoleFraction, problem, 0, {2}),
        solved(ObservationKind::PhaseMoleFraction, problem, 1, {0}), solved(ObservationKind::PhaseFraction, problem, 0),
        amount_ratio(problem, 0, 1, {2}, {2})};
    b1.coordinates.coordinates = {{{{ParameterFamily::Kij, {2, 6}}}, 0.045, 0.1}, reaction_coordinate("A+B=C")};
    reaction_data(b1);
    b1.parameters.k_ij[2] = b1.parameters.k_ij[6] = 0.045;
    for (auto& row : b1.observations) row.problem.reactions[0].equilibrium_constant = 1.0;
    recover_reactions("reaction-b1", b1, {0.05, 0.1});
    FitRequest b2;
    b2.controls = budgeted_controls();
    b2.parameters = properties_support::neutral_binary();
    const auto bubble = reacting_bubble(HelmholtzModel(b2.parameters));
    const double a_truth = bubble.reactions[0].correlation->a;
    for (const double T : {170.0, 180.0, 190.0}) {
        auto cell = bubble;
        cell.temperature = FreeValue::fixed(T);
        b2.observations.push_back(solved(ObservationKind::EquilibriumPressure, cell));
        b2.observations.push_back(solved(ObservationKind::PhaseMoleFraction, cell, 1, {0}));
    }
    b2.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::B, 180.0)};
    reaction_data(b2);
    for (auto& row : b2.observations) { row.problem.reactions[0].correlation->a += 0.1; row.problem.reactions[0].correlation->b -= 20.0; }
    recover_reactions("reaction-b2", b2, {a_truth, 300.0});
    FitRequest b3;
    b3.controls = budgeted_controls();
    b3.parameters = properties_support::neutral_ternary();
    auto vapor = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    b3.coordinates.coordinates = {reaction_coordinate("reaction_a"), reaction_coordinate("reaction_b")};
    for (const auto& feed : {std::vector<double>{0.6, 0.3, 0.1}, std::vector<double>{0.2, 0.3, 0.5}}) {
        vapor.feed.values = feed;
        if (!b3.observations.empty()) std::reverse(vapor.reactions.begin(), vapor.reactions.end());
        for (auto& r : vapor.reactions) r.equilibrium_constant = 1.0;
        regular_reaction_state(HelmholtzModel(b3.parameters), solved(ObservationKind::PhaseMoleFraction, vapor, 0, {0}), "reaction-b3/feed_" + std::to_string(b3.observations.size() / 2));
        for (auto& r : vapor.reactions) r.equilibrium_constant = std::exp(r.name == "reaction_a" ? 0.2 : -0.1);
        b3.observations.push_back(solved(ObservationKind::PhaseMoleFraction, vapor, 0, {0}));
        b3.observations.push_back(solved(ObservationKind::PhaseMoleFraction, vapor, 0, {2}));
    }
    reaction_data(b3);
    for (auto& row : b3.observations) for (auto& r : row.problem.reactions) r.equilibrium_constant = 1.0;
    recover_reactions("reaction-b3/name_tying", b3, {0.2, -0.1});
    // A large affine offset loses the supplied 0.1 ln K during resolution. Validation's
    // compiled problem is therefore stale even though its EOS active-value key is empty.
    FitRequest preparation;
    preparation.controls = budgeted_controls();
    preparation.parameters = equilibrium_support::synthetic_reactive_ternary();
    auto source = equilibrium_support::synthetic_reactive_problem();
    source.standard_gibbs_j_per_mol.clear();
    preparation.observations = {solved(ObservationKind::PhaseMoleFraction, source, 0, {2})};
    preparation.coordinates.coordinates = {reaction_coordinate("A+B=C")};
    preparation.coordinates.coordinates[0].origin = std::ldexp(1.0, 54);
    reaction_data(preparation);  // data at the independently solved resolved value ln K = 0
    preparation.observations[0].problem.reactions[0].equilibrium_constant = std::exp(0.1);
    const auto u = preparation.coordinates.start(preparation.parameters, preparation.observations);
    const double resolved = preparation.coordinates.physical(u).front();
    CHECK(resolved == 0.0);
    const auto supplied = predict(HelmholtzModel(preparation.parameters), preparation.observations[0]);
    CHECK(supplied.status == ObservationStatus::Available);
    const double signal = std::abs(std::log(*supplied.value / preparation.observations[0].observed));
    CHECK(signal >= 1.0e-2);
    const auto rebuilt = fit(preparation);
    const auto inputs = fields({field("supplied_lnK", 0.1), field("resolved_lnK", resolved), field("log_prediction_change", signal)});
    agree("reaction-preparation/prediction", "x_C", "1", rebuilt.predictions.front(), preparation.observations[0].observed,
          0.0, 1.0e-12, inputs, "a changed resolved reaction constant invalidates validation preparation", "reaction-fit");
    agree("reaction-preparation/initial_cost", "cost", "1", rebuilt.initial_cost, 0.0, 1.0e-14, 0.0, inputs,
          "initial prediction uses the resolved constant rather than the supplied compiled problem", "reaction-fit");
}

void check_reaction_rejections() {
    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = properties_support::neutral_ternary();
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    request.coordinates.coordinates = {reaction_coordinate("reaction_a")};
    request.observations = {solved(ObservationKind::PhaseMoleFraction, problem, 0, {0})};
    const auto reject = [&](FitRequest candidate, const std::string& cause) {
        const auto result = fit(candidate);
        truth("reaction-rejection/" + cause, "an overdetermined or fixed reaction input is refused by name",
            result.status == "invalid_setup" && result.message.find(cause) != std::string::npos
            && result.message.find("reaction_a") != std::string::npos,
            fields({field("message", result.message)}), "reaction-fit");
    };
    auto affinity = request;
    affinity.observations[0].problem.reactions[0].equilibrium_constant.reset();
    affinity.observations[0].problem.reactions[0].affinity_j_per_mol = 0.0;
    reject(affinity, "affinity target");
    auto gibbs = request;
    gibbs.observations[0].problem.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    reject(gibbs, "supplied standard Gibbs");
    auto dependent = request;
    dependent.observations[0].problem.reactions[1].stoichiometry = problem.reactions[0].stoichiometry;
    reject(dependent, "dependent reaction");
    auto fixed = request;
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(3, {{200.0, 400.0, true, true, {IdealPolynomial{{30.0}, 300.0, 0.0, 0.0}, 1.0e5}}});
    fixed.thermochemistry = record;
    reject(fixed, "K fixed by thermochemistry");
    // Two named source constants with different declarations are physically different coordinates.
    auto conflicting = request;
    conflicting.observations.push_back(request.observations[0]);
    conflicting.observations[1].problem.reactions[0].equilibrium_constant = std::exp(0.01);
    reject(conflicting, "tied declarations or starts disagree");
    auto scalar_b = request;
    scalar_b.coordinates.coordinates[0].reaction->coefficient = ReactionCoefficient::B;
    reject(scalar_b, "scalar K admits only A");
    auto reference = request;
    reference.parameters = properties_support::aqueous_electrolyte();
    auto source = referenced_formula(std::nullopt);
    source.reactions[0].name = "reaction_a";
    source.reactions[0].reference->source_basis = ReferenceSourceBasis::CommonMolalityInfiniteDilution;
    reference.observations = {solved(ObservationKind::PhaseMoleFraction, source, 0, {1}),
                              solved(ObservationKind::PhaseMoleFraction, source, 0, {1})};
    reference.observations[1].problem.source_reference->standard_molality_mol_per_kg = 2.0;
    reject(reference, "tied source bases disagree");
    reference.observations[1].problem.source_reference = source.source_reference;
    reference.observations[1].problem.reactions[0].correlation->reference_temperature = 300.0;
    reject(reference, "tied declarations or starts disagree");
    // #448: a tie needs one convention per reacting species; the explicit pure-liquid solvent is the default.
    // Raw mole fractions admit mixed references, so only the tie refuses the pure-gas row.
    reference.observations[1].problem.reactions[0].correlation->reference_temperature = 298.15;
    for (auto& row : reference.observations)
        row.problem.reactions[0].reference->source_basis = ReferenceSourceBasis::RawMoleFractionInfiniteDilution;
    reference.observations[1].problem.source_reference->species_reference = {SpeciesReference::PureGas,
        SpeciesReference::InfiniteDilution, SpeciesReference::InfiniteDilution};
    reject(reference, "tied species references disagree: component 0 pure liquid vs pure gas");
    reference.observations[1].problem.source_reference->species_reference[0] = SpeciesReference::PureLiquid;
    const auto admitted = [&] {
        try { return reference.coordinates.start(reference.parameters, reference.observations).size() == 1; }
        catch (const std::invalid_argument&) { return false; }
    };
    truth("reaction-rejection/matching-species-reference", "a tie whose reacting species share one effective convention is admitted", admitted(),
          fields({field("water_reference", "explicit pure liquid equal to the default")}), "issue-448");
    // #433: rows that differ only in species_reference are different states; each keeps its own solve
    // (fixed ln K and an EOS coordinate, so no tie across conventions).
    FitRequest shared = reference;
    shared.observations = {solved(ObservationKind::PhaseMoleFraction, referenced_formula(std::nullopt), 0, {1})};
    shared.observations.push_back(shared.observations[0]);
    shared.observations[1].problem.source_reference->species_reference = {SpeciesReference::PureGas,
        SpeciesReference::InfiniteDilution, SpeciesReference::InfiniteDilution};
    shared.coordinates.coordinates = {{{{ParameterFamily::SegmentCount, {0}}}, 0.0, 1.0}};
    FitResult failures; ObservationEvaluation rows;
    CHECK(ObservationEvaluator(shared, failures, {0, 1}).evaluate(shared.coordinates.start(shared.parameters, shared.observations), false, rows));
    agree("reaction-rejections/species-reference", "x_Na", "1", rows.predictions[1],
          *predict(HelmholtzModel(shared.parameters), shared.observations[1]).value, 0.0, 1.0e-10,
          "T=298.15 K, water pure-gas reference", "a row differing only in species_reference keeps its own solve", "issue-433");
}

// #211: the source constant is measured directly; no EOS reference transfer enters it.
Observation log_constant(const EquilibriumProblem& problem, const char* name) {
    Observation row;
    row.kind = ObservationKind::LogEquilibriumConstant;
    row.problem = problem;
    row.reaction = name;
    return row;
}

void check_log_constant() {
    const auto parameters = properties_support::neutral_ternary();
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    problem.reactions[0].equilibrium_constant.reset();
    problem.reactions[0].correlation = ReactionLogPolynomial{0.4, 80.0, 0.5, 0.001, 300.0, 150.0, 400.0,
                                                           equilibrium::EOS_STANDARD_STATE_ID};
    CoordinateSet coordinates{{reaction_coordinate("reaction_a", ReactionCoefficient::B, 300.0),
        {{{ParameterFamily::SegmentCount, {0}}}, 0.0, 1.0}, reaction_coordinate("reaction_a"),
        reaction_coordinate("reaction_b"), reaction_coordinate("reaction_a", ReactionCoefficient::C),
        reaction_coordinate("reaction_a", ReactionCoefficient::D, 1.0 / 300.0)}};
    const HelmholtzModel model(parameters, coordinates.active_parameters());
    for (const double T : {200.0, 250.0, 300.0, 350.0}) {
        problem.temperature = FreeValue::fixed(T);
        std::reverse(problem.reactions.begin(), problem.reactions.end());
        auto row = log_constant(problem, "reaction_a");
        validate(row, model);
        const auto prediction = predict(model, row);
        const double value = 0.4 + 80.0 / T + 0.5 * std::log(T / 300.0) + 0.001 * T;
        const auto id = "ln-k/" + std::to_string(T);
        CHECK(prediction.status == ObservationStatus::Available);
        agree(id + "/value", "ln_K_src", "1", *prediction.value, value, 0.0, 0.0,
              fields({field("T_K", T), field("Tref_K", 300.0)}), "declared source polynomial", "ln-k");
        const auto gradient = differentiate_group(model, {&row}, prediction, coordinates).front();
        CHECK(gradient.status == ObservationStatus::Available);
        // Physical rows: EOS, b, a, other reaction a, c, d; optimizer order is interleaved.
        const std::vector<double> expected{0.0, 1.0 / T, 1.0, 0.0, std::log(T / 300.0), T};
        for (std::size_t k = 0; k < expected.size(); ++k)
            agree(id + "/column/" + std::to_string(k), "dln_K_src_dtheta", "native",
                gradient.value[k], expected[k], 0.0, 0.0, fields({field("T_K", T), field("Tref_K", 300.0)}),
                "exact coefficient basis; EOS and other reaction columns vanish", "ln-k");
    }
    // A referenced source value remains the source value even though its density-basis target differs.
    auto reference = log_constant(referenced_formula(std::nullopt), "H2O=Na++Cl-");
    const HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    validate(reference, aqueous);
    agree("ln-k/source", "ln_K_src", "1", *predict(aqueous, reference).value, std::log(1.0e-8),
          0.0, 0.0, "T=298.15 K; raw mole-fraction source", "no EOS reference transfer in direct data", "ln-k");
    auto scalar = log_constant(problem, "reaction_b");
    agree("ln-k/scalar", "ln_K_src", "1", *predict(model, scalar).value,
          std::log(*problem.reactions[1].equilibrium_constant), 0.0, 0.0,
          "T=350 K", "scalar K is observed as ln K", "ln-k");

    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = properties_support::neutral_binary();
    const auto bubble = reacting_bubble(HelmholtzModel(request.parameters));
    const auto law = *bubble.reactions[0].correlation;
    request.coordinates.coordinates = {reaction_coordinate("iso"), reaction_coordinate("iso", ReactionCoefficient::B, 180.0)};
    for (const double T : {170.0, 180.0, 190.0}) {
        auto cell = bubble;
        cell.temperature = FreeValue::fixed(T);
        auto row = log_constant(cell, "iso");
        row.observed = law.a + law.b / T + law.c * std::log(T / law.reference_temperature) + law.d * T;
        request.observations.push_back(row);
    }
    auto equilibrium = solved(ObservationKind::EquilibriumPressure, bubble);
    equilibrium.observed = *predict(HelmholtzModel(request.parameters), equilibrium).value;
    request.observations.push_back(equilibrium);
    for (auto& row : request.observations) { row.problem.reactions[0].correlation->a += 0.1; row.problem.reactions[0].correlation->b -= 20.0; }
    recover_reactions("ln-k/recovery", request, {law.a, law.b});

    // Observation admission must reject fixed targets even when only an EOS coordinate is fitted.
    FitRequest rejected;
    rejected.controls = budgeted_controls();
    rejected.parameters = parameters;
    rejected.coordinates.coordinates = {{{{ParameterFamily::SegmentCount, {0}}}, 0.0, 1.0}};
    rejected.observations = {log_constant(problem, "reaction_a")};
    const auto reject = [&](FitRequest candidate, const std::string& cause) {
        const auto result = fit(candidate);
        truth("ln-k/rejection/" + cause, "the direct source target is refused by reaction name",
            result.status == "invalid_setup" && result.message.find(cause) != std::string::npos
            && result.message.find("reaction_a") != std::string::npos,
            fields({field("message", result.message)}), "ln-k");
    };
    // At the last reversal reaction_a is first.
    auto affinity = rejected;
    affinity.observations[0].problem.reactions[0].correlation.reset();
    affinity.observations[0].problem.reactions[0].affinity_j_per_mol = 0.0;
    reject(affinity, "affinity target");
    auto dependent = rejected;
    dependent.observations[0].problem.reactions[1].stoichiometry = problem.reactions[0].stoichiometry;
    reject(dependent, "dependent reaction");
    auto gibbs = rejected;
    gibbs.observations[0].problem.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    reject(gibbs, "supplied standard Gibbs");
    auto fixed = rejected;
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(3, {{200.0, 400.0, true, true, {IdealPolynomial{{30.0}, 300.0, 0.0, 0.0}, 1.0e5}}});
    fixed.thermochemistry = record;
    reject(fixed, "K fixed by thermochemistry");
}

// #210: finite enthalpy increments on a referenced, reaction-completed ethane pair.
Parameters addition_parameters() {
    auto result = properties_support::parameters(2);
    const auto ethane = properties_support::ethane();
    for (std::size_t i = 0; i < 2; ++i) {
        result.segment_count[i] = ethane.segment_count[0];
        result.sigma_constant[i] = ethane.sigma_constant[0];
        result.dispersion_energy_over_k[i] = ethane.dispersion_energy_over_k[0];
        result.molar_mass_kg_per_mol[i] = ethane.molar_mass_kg_per_mol[0];
    }
    return result;
}

ThermochemistryRecord addition_record() {
    ThermochemistryRecord result;
    result.standard_pressure = 1.0e5;
    result.component_intervals = {
        {{200.0, 400.0, true, true, {IdealPolynomial{{52.0, 0.12}, 300.0, 0.0, 0.0}, 1.0e5}}},
        {}};
    return result;
}

ReactionLogPolynomial addition_law(double reference_temperature = 300.0) {
    const double a = std::log(0.3) + 3.0 * std::log(reference_temperature / 300.0);
    return {a, 200.0, 3.0, -0.004, reference_temperature, 250.0, 350.0, "analytic-source"};
}

EquilibriumProblem addition_problem(double temperature, double pressure, std::vector<double> feed,
                                    std::optional<double> reference_pressure = {}, bool free_pressure = false,
                                    std::optional<double> vapor_fraction = {}) {
    EquilibriumProblem result;
    result.phases = {equilibrium_support::phase("L", PhaseKind::Liquid)};
    if (vapor_fraction) {
        auto vapor = equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(*vapor_fraction));
        result.phases.push_back(std::move(vapor));
    } else if (free_pressure) {
        result.phases.push_back(equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0)));
    }
    result.temperature = FreeValue::fixed(temperature);
    result.pressure = free_pressure ? FreeValue::variable(pressure) : FreeValue::fixed(pressure);
    result.feed.kind = equilibrium::FeedKind::Amounts;
    result.feed.values = std::move(feed);
    result.source_reference = NeutralReferenceDeclaration{0, {0.0, 1.0}, 1.0};
    result.reactions = {{"2A=B", {-2.0, 1.0}, std::nullopt, std::nullopt, addition_law(),
        ReactionReference{reference_pressure, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
    result.solver_options.values[epcsaft::policy::row] = 1.0e-10;
    return result;
}

Observation addition_observation(const EquilibriumProblem& initial, const EquilibriumProblem& final,
                                 double amount, double feed_temperature = 260.0,
                                 double feed_pressure = 1.0e5, std::vector<double> feed = {1.0, 0.0}) {
    Observation result;
    result.kind = ObservationKind::AdditionEnthalpy;
    result.problem = initial;
    result.final_problem = final;
    result.feed_amount = amount;
    result.temperature = feed_temperature;
    result.pressure = feed_pressure;
    result.composition = std::move(feed);
    result.branch = DensityBranch::Vapor;
    result.form = ResidualForm::Difference;
    result.scale = 1.0;
    return result;
}

CoordinateSet addition_coordinates() {
    CoordinateSet result;
    result.coordinates = {
        {{{ParameterFamily::Kij, {1, 2}}}, 0.0, 1.0, -0.2, 0.2},
        {{{ParameterFamily::DispersionEnergyOverK, {0}}, {ParameterFamily::DispersionEnergyOverK, {1}}}, 191.42, 1.0, 150.0, 250.0},
        reaction_coordinate("2A=B"),
        reaction_coordinate("2A=B", ReactionCoefficient::B),
        reaction_coordinate("2A=B", ReactionCoefficient::C),
        reaction_coordinate("2A=B", ReactionCoefficient::D)};
    return result;
}

struct AdditionEndpoint final {
    double enthalpy = 0.0;
    std::vector<double> phase_amount, phase_enthalpy;
};

AdditionEndpoint endpoint_enthalpy(const HelmholtzModel& model, const EquilibriumProblem& problem,
                                   const equilibrium::CompiledSystem& compiled,
                                   const EquilibriumResult& solved_state) {
    AdditionEndpoint result;
    result.phase_amount.assign(problem.phases.size(), 0.0);
    result.phase_enthalpy.assign(problem.phases.size(), 0.0);
    std::vector<equilibrium::SolvedStateObservable> observables;
    std::vector<std::size_t> phases;
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        if (problem.phases[phase].is_pinned() && problem.phases[phase].amount.fraction == 0.0) continue;
        equilibrium::SolvedStateObservable amount, enthalpy;
        amount.kind = equilibrium::SolvedStateObservableKind::PhaseAmount;
        amount.phase = phase;
        enthalpy.kind = equilibrium::SolvedStateObservableKind::PhaseProperty;
        enthalpy.phase = phase;
        enthalpy.property = PropertyObservable::TotalEnthalpy;
        observables.push_back(amount);
        observables.push_back(enthalpy);
        phases.push_back(phase);
    }
    equilibrium::SolvedStateActionDirection zero;
    zero.feed_amounts_mol.assign(model.parameters().component_count(), 0.0);
    zero.active_parameters.assign(model.active_parameters().size(), 0.0);
    const auto batch = equilibrium::solved_state_actions(compiled, solved_state, {observables, {zero}});
    CHECK(batch.results.size() == 2 * phases.size());
    for (std::size_t i = 0; i < phases.size(); ++i) {
        const auto& amount = batch.results[2 * i];
        const auto& enthalpy = batch.results[2 * i + 1];
        CHECK(amount.status == equilibrium::SolvedStateActionStatus::Available && amount.value);
        CHECK(enthalpy.status == equilibrium::SolvedStateActionStatus::Available && enthalpy.value);
        result.phase_amount[phases[i]] = *amount.value;
        result.phase_enthalpy[phases[i]] = *enthalpy.value;
        result.enthalpy += *amount.value * *enthalpy.value;
    }
    return result;
}

double addition_oracle_phase_enthalpy(const HelmholtzModel& pure, const EquilibriumResult& solved_state,
                                     const ReactionLogPolynomial& law, std::size_t phase) {
    const double T = *solved_state.temperature;
    const double P_ref = *solved_state.pressure;
    const double dt = T - 300.0;
    const double h_a = 52.0 * dt + 0.06 * dt * dt;
    const double ell = -law.b / (T * T) + law.c / T + law.d;
    const double h_b = 2.0 * h_a + GAS_CONSTANT_J_PER_MOL_K * T * T * ell
        + pure.state_at_pressure(T, P_ref, {1.0}, DensityBranch::Liquid).residual_enthalpy;
    const double x_a = solved_state.mole_fractions[phase * 2];
    const double x_b = solved_state.mole_fractions[phase * 2 + 1];
    const double residual = pure.state(T, solved_state.molar_densities[phase], {1.0}).residual_enthalpy;
    return x_a * h_a + x_b * h_b + residual;
}

double addition_oracle_endpoint(const HelmholtzModel& pure, const EquilibriumResult& solved_state,
                                const ReactionLogPolynomial& law) {
    double result = 0.0;
    for (std::size_t phase = 0; phase < solved_state.phase_count; ++phase)
        result += solved_state.amounts[phase] * addition_oracle_phase_enthalpy(pure, solved_state, law, phase);
    return result;
}

double total_b_amount(const EquilibriumResult& solved_state) {
    double result = 0.0;
    for (std::size_t phase = 0; phase < solved_state.phase_count; ++phase)
        result += solved_state.amounts[phase] * solved_state.mole_fractions[phase * 2 + 1];
    return result;
}

struct AdditionVariant final {
    std::string name;
    Observation observation;
    Prediction prediction;
    std::vector<double> gradient;
};

AdditionVariant addition_variant(const std::string& name, double amount, std::optional<double> reference_pressure,
                                 bool free_pressure = false, std::optional<double> vapor_fraction = {},
                                 bool unequal_temperature = false) {
    const double pressure = free_pressure ? 1.5e6 : 3.0e6;  // original Case N regularity input
    auto initial = addition_problem(260.0, pressure, {2.0, 0.2}, reference_pressure, free_pressure, vapor_fraction);
    auto final = addition_problem(unequal_temperature ? 265.0 : 260.0, pressure, {2.0 + amount, 0.2},
                                  reference_pressure, free_pressure, vapor_fraction);
    if (unequal_temperature) {
        initial.reactions[0].correlation = addition_law(260.0);
        final.reactions[0].correlation = addition_law(260.0);
    }
    auto row = addition_observation(initial, final, amount);
    return {name, std::move(row), {}, {}};
}

double addition_coordinate_gradient(const std::vector<double>& gradient, std::size_t column) {
    return column == 1 ? gradient[1] + gradient[2] : gradient[column < 2 ? column : column + 1];
}

// Case I of #210, composed into the actual two-endpoint fitting observation.
void check_ionic_addition_enthalpy() {
    const auto parameters = properties_support::aqueous_electrolyte();
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals = {
        {{200.0, 500.0, true, true, {IdealPolynomial{{33.59}, 298.15, -241826.0, 188.83}, 1.0e5}}}, {}, {}};
    constexpr double T = 298.15, P = 1.0e5, amount = 0.2;
    auto initial = equilibrium_support::flash(3);
    initial.temperature = FreeValue::fixed(T);
    initial.pressure = FreeValue::fixed(P);
    initial.feed.values = {1.0, 1.0e-6, 1.0e-6};
    initial.source_reference = NeutralReferenceDeclaration{0, {0.0, 1.0, 1.0}, 1.0};
    initial.reactions = {{"H2O=Na++Cl-", {-1.0, 1.0, 1.0}, std::nullopt, std::nullopt,
        ReactionLogPolynomial{std::log(1.0e-8) + 2000.0 / T, -2000.0, 0.0, 0.0, T, 250.0, 350.0, "analytic-source"},
        ReactionReference{std::nullopt, ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
    auto final = initial;
    final.feed.values[0] += amount;
    auto row = addition_observation(initial, final, amount, T, P, {1.0, 0.0, 0.0});
    row.branch = DensityBranch::Liquid;
    CoordinateSet coordinates;
    coordinates.coordinates = {
        // Solvent segment count changes packing/density, rather than a Born
        // share that nearly cancels between the ionic infinite-dilution standards.
        {{{ParameterFamily::SegmentCount, {0}}}, parameters.segment_count[0], 1.0, 1.0, 1.5},
        reaction_coordinate("H2O=Na++Cl-", ReactionCoefficient::B, T)};
    coordinates.coordinates[1].origin = -2000.0;
    coordinates.coordinates[1].lower = -3000.0;
    coordinates.coordinates[1].upper = -1000.0;
    auto shift = reaction_coordinate("H2O=Na++Cl-", ReactionCoefficient::EnthalpyShift, 2500.0);
    shift.reaction->pivot_temperature = 313.15;
    coordinates.coordinates.push_back(shift);
    const HelmholtzModel model(parameters, coordinates.active_parameters(), record);
    const auto prediction = predict(model, row);
    CHECK(prediction.status == ObservationStatus::Available && prediction.value);
    const auto exact = differentiate_group(model, {&row}, prediction, coordinates).front();
    CHECK(exact.status == ObservationStatus::Available && exact.value.size() == 3);
    const auto endpoint = [&](const HelmholtzModel& rebuilt, const EquilibriumProblem& problem, double gauge = 0.0, bool centered = false) {
        const auto solved = equilibrium::solve_equilibrium(rebuilt, problem);
        CHECK(solved.success);
        const auto compiled = equilibrium::compile_problem(rebuilt, problem);
        const auto source = equilibrium::reference_calorics(compiled, T, P, false, false, false, 1);
        CHECK(source.status == equilibrium::SolvedStateActionStatus::Available);
        const auto standard = [&](std::size_t order) -> const std::vector<double>& {
            const auto found = std::find_if(source.gibbs.begin(), source.gibbs.end(), [&](const auto& column) {
                return column.temperature_order == order && column.pressure_order == 0 && !column.parameter;
            });
            CHECK(found != source.gibbs.end() && found->column.status == equilibrium::SolvedStateActionStatus::Available
                && found->column.values.size() == 3);
            return found->column.values;
        };
        const auto& g = standard(0);
        const auto& gT = standard(1);
        const auto state = rebuilt.state(T, solved.molar_densities[0], solved.mole_fractions);
        // Accumulate the independently composed oracle in extended precision:
        // the trace ionic increment subtracts ~3e5 J endpoints.
        long double h = state.residual_enthalpy;
        // Subtract the feed enthalpy on the conserved water-equivalent basis
        // (1, 1/2, 1/2), not total moles: H2O -> Na+ + Cl- changes N.
        // Then (H'_final-H'_initial)/n_add is exactly the same increment,
        // without amplifying endpoint amount/normalization rounding by 3e5 J/mol.
        const double feed_residual = centered
            ? rebuilt.state_at_pressure(T, P, {1.0, 0.0, 0.0}, DensityBranch::Liquid).residual_enthalpy : 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            const double basis = i == 0 ? 1.0 : 0.5;
            // g_rho includes RT ln(RT/P0); its physical standard enthalpy restores RT.
            const long double source_h = static_cast<long double>(g[i]) - T * static_cast<long double>(gT[i])
                + GAS_CONSTANT_J_PER_MOL_K * T + gauge * parameters.charge_number[i]
                + (centered ? basis * 241826.0 : 0.0);
            h += static_cast<long double>(solved.mole_fractions[i]) * (source_h - basis * feed_residual);
        }
        return solved.amounts[0] * h;
    };
    const auto feed_h = [&](const HelmholtzModel& rebuilt) {
        return rebuilt.state_at_pressure(T, P, {1.0, 0.0, 0.0}, DensityBranch::Liquid).residual_enthalpy - 241826.0;
    };
    const double H_initial = endpoint(model, initial), H_final = endpoint(model, final);
    const double oracle = (H_final - H_initial) / amount - feed_h(model);
    const double allowance = 1.0e-9 * (std::abs(H_initial) + std::abs(H_final)) / amount + 1.0e-6;
    agree("addition-ionic/value", "delta_h_add", "J/mol", *prediction.value, oracle, allowance, 0.0,
        fields({field("H_initial_J", H_initial), field("H_final_J", H_final), field("n_add_mol", amount)}),
        "the ionic observation equals independently rebuilt endpoint and feed enthalpies", "addition-ionic");
    const double centered_oracle = (endpoint(model, final, 0.0, true) - endpoint(model, initial, 0.0, true)) / amount;
    agree("addition-ionic/conserved-basis-value", "delta_h_add", "J/mol", centered_oracle, oracle, allowance, 0.0,
        field("conserved_basis", std::vector<double>{1.0, 0.5, 0.5}),
        "subtracting the feed on the conserved basis preserves the two-endpoint increment", "addition-ionic");
    const double shifted = (endpoint(model, final, 1.0e6) - endpoint(model, initial, 1.0e6)) / amount - feed_h(model);
    agree("addition-ionic/charge-gauge", "delta_h_add", "J/mol", shifted, oracle, allowance, 0.0,
        field("charge_gauge_J_per_mol", 1.0e6), "a charge-neutral increment is invariant to ion enthalpy gauge", "addition-ionic");
    const auto residual_endpoint = [&](const EquilibriumProblem& problem) {
        const auto solved = equilibrium::solve_equilibrium(model, problem);
        CHECK(solved.success);
        return solved.amounts[0] * model.state(T, solved.molar_densities[0], solved.mole_fractions).residual_enthalpy;
    };
    const double reference_share = (H_final - H_initial - residual_endpoint(final) + residual_endpoint(initial)) / amount + 241826.0;
    truth("addition-ionic/reference-signal", "the ionic reference contribution exceeds the subtraction allowance",
        std::abs(reference_share) > allowance, fields({field("reference_share_J_per_mol", reference_share),
            field("allowance_J_per_mol", allowance)}), "addition-ionic");
    const auto start = coordinates.start(parameters, {row}, record);
    using binary128_support::Wide;
    const std::vector<binary128_support::Quantity> enthalpy{{"H", "J", 0.0, {}}};
    for (std::size_t column = 0; column < coordinates.coordinates.size(); ++column) {
        // Delta h_add composed in binary128: WideIonic re-solves each endpoint's density and
        // dissociation and the added water's state from the EOS kernels; the double solves only
        // seed and check it, so no double endpoint rounding (~7e-10 J/mol) enters the differences.
        const auto composed = [&](double shift) {
            auto u = start;
            u[column] += shift;
            auto changed = std::vector<Observation>{row};
            coordinates.apply(changed, u);
            const auto shifted = coordinates.apply(parameters, u);
            const HelmholtzModel rebuilt(shifted, coordinates.active_parameters(), record);
            Wide h[2], feed = 0;
            for (std::size_t k = 0; k < 2; ++k) {
                const auto& problem = k ? changed[0].final_problem : changed[0].problem;
                const auto solved = equilibrium::solve_equilibrium(rebuilt, problem);
                CHECK(solved.success);
                const double rho = solved.molar_densities[0];
                binary128_support::WideIonic oracle(shifted, record, problem, solved.mole_fractions[1]);
                h[k] = oracle.values(enthalpy, {}, 1.0, 0, solved, rho)[0];
                if (k) feed = oracle.solvent_enthalpy(rho);
            }
            return (h[1] - h[0]) / Wide(amount) - feed;
        };
        const auto central = [&](double h) { return (composed(h) - composed(-h)) / Wide(2.0 * h); };
        const double step = 0.02;
        const Wide d1 = central(step), d2 = central(step / 2), d4 = central(step / 4);
        const double scale = coordinates.coordinates[column].scale;
        const double R1 = double((4 * d2 - d1) / 3) / scale, reference = double((4 * d4 - d2) / 3) / scale;
        const double action = exact.value[column];
        const double tolerance = 1.0e-6 * std::abs(reference);
        const double uncertainty = std::abs(reference - R1);
        truth("addition-ionic/column-refinement/" + std::to_string(column), "the binary128 endpoint/feed derivative is resolved to 1e-7",
            uncertainty <= 1.0e-7 * std::abs(reference),
            fields({field("coarse", double(d1) / scale), field("middle", double(d2) / scale), field("fine", double(d4) / scale),
                field("uncertainty", uncertainty), field("relative_uncertainty", uncertainty / std::abs(reference))}), "addition-ionic");
        agree("addition-ionic/column/" + std::to_string(column), "d_delta_h_add_dtheta", "native", action, reference,
            tolerance, 0.0, fields({field("column", column), field("step", step), field("uncertainty", uncertainty)}),
            "one EOS and one source-law column agree with binary128 composed endpoint/feed differences", "addition-ionic");
    }
}

// Issue #225 F6: the same fitting consumer on two referenced free-T endpoints.
void check_free_temperature_addition() {
    const auto parameters=addition_parameters();const auto record=addition_record();
    const auto all=addition_coordinates();CoordinateSet coordinates;
    coordinates.coordinates={all.coordinates[1],all.coordinates[3]};
    const HelmholtzModel model(parameters,coordinates.active_parameters(),record);
    for(bool fixed:{false,true}) {
        auto initial=addition_problem(270,1.952823e6,{2,0.2},fixed?std::optional<double>(5e6):std::nullopt);
        initial.temperature=FreeValue::variable(270);
        initial.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
        auto final=initial;final.feed.values={2.01,0.2};
        const auto observation=addition_observation(initial,final,0.01,300,1.5e6);
        const auto center=predict(model,observation);
        if(center.status!=ObservationStatus::Available||!center.value)
            throw std::runtime_error(std::string("F6/free-T-addition/")+(fixed?"fix: ":"sys: ")+center.message);
        CHECK(center.solved&&center.final_solved&&center.compiled&&center.final_compiled);
        const auto exact=differentiate_group(model,{&observation},center,coordinates).front();
        if(exact.status!=ObservationStatus::Available||exact.value.size()!=3)
            throw std::runtime_error("F6/free-T-addition action: "+exact.message);
        const auto start=coordinates.start(parameters,{observation},record);
        const std::string id=std::string("F6/free-T-addition/")+(fixed?"fix":"sys");
        const auto inputs=fields({field("n_add_mol",0.01),field("initial_feed",initial.feed.values),
            field("final_feed",final.feed.values),field("initial_T_K",*center.solved->temperature),
            field("final_T_K",*center.final_solved->temperature),field("system_P_Pa",1.952823e6),
            field("source_P_Pa",fixed?5e6:1.952823e6),field("feed_T_K",300),field("feed_P_Pa",1.5e6),
            field("initial_roots",center.solved->molar_densities),field("final_roots",center.final_solved->molar_densities),
            field("initial_residual",center.solved->residual_norm()),field("final_residual",center.final_solved->residual_norm()),
            field("initial_iterations",center.solved->iterations),field("final_iterations",center.final_solved->iterations)});
        const auto first=endpoint_enthalpy(model,initial,*center.compiled,*center.solved);
        const auto second=endpoint_enthalpy(model,final,*center.final_compiled,*center.final_solved);
        auto feed_record=record;feed_record.component_intervals.resize(1);
        const HelmholtzModel feed(properties_support::ethane(),{},feed_record);
        const auto incoming=feed.state_at_pressure(300,1.5e6,{1},DensityBranch::Vapor);
        CHECK(incoming.total_enthalpy);
        const double value=(second.enthalpy-first.enthalpy)/0.01-*incoming.total_enthalpy;
        agree(id+"/value","delta_h_add","J/mol",*center.value,value,1e-5,1e-6,inputs,
            "two solved endpoint extensive enthalpies and declared feed enthalpy compose the increment","reactive-caloric");
        for(std::size_t column=0;column<2;++column) {
            const double action=column==0?exact.value[0]+exact.value[1]:exact.value[2];
            const double sigma=column==0?std::abs(exact.value[0])+std::abs(exact.value[1]):std::abs(action);
            const auto value_at=[&](double shift) {
                auto u=start;u[column]+=shift;std::vector<Observation> rows={observation};
                coordinates.apply(rows,u);
                const HelmholtzModel rebuilt(coordinates.apply(parameters,u),coordinates.active_parameters(),record);
                const auto result=predict(rebuilt,rows.front());
                if(result.status!=ObservationStatus::Available||!result.value)
                    throw std::runtime_error(id+" rebuilt: "+result.message);
                CHECK(result.solved&&result.final_solved&&result.solved->topology_event==equilibrium::TopologyEvent::None
                    &&result.final_solved->topology_event==equilibrium::TopologyEvent::None);
                return *result.value;
            };
            double step=column==0?0.02:2;
            const auto difference=[&](double h){return (value_at(h)-value_at(-h))/(2*h);};
            double coarse=difference(step);bool passed=false;
            for(std::size_t rung=0;rung<6;++rung) {
                step/=2;const double fine=difference(step),rich=(4*fine-coarse)/3,refinement=std::abs(fine-coarse)/3;
                const double target=1e-5+1e-6*std::max(std::abs(rich),sigma);
                const bool resolved=refinement<=target/4;
                passed=resolved&&std::isfinite(action)&&std::abs(action-rich)<=target;
                const auto row=fields({inputs,field("coordinate",column==0?"tied-epsilon":"reaction-b"),
                    field("step",step),field("D_h",coarse),field("D_h2",fine),field("Richardson",rich),
                    field("refinement",refinement),field("sigma",sigma),field("target",target),
                    field("error",std::abs(action-rich)),field("error_over_target",std::abs(action-rich)/target)});
                const auto comparison=id+"/column-"+std::to_string(column)+"/rung-"+std::to_string(rung);
                if(!resolved&&rung<5) recording_support::observe(comparison.c_str(),"d_delta_h_add","J/mol",action,rich,row,
                    "Refinement","native-rebuilt-prediction","Unresolved coarse step; retained and refined","reactive-caloric");
                else recording_support::compare_custom(comparison.c_str(),"d_delta_h_add","J/mol",action,rich,row,
                    "Relationship","native-rebuilt-prediction","Free-T endpoints and feed use exact native coefficient motion",
                    "abs(exact-R)<=atol+1e-6*max(abs(R),sigma) && refinement<=target/4",passed,"reactive-caloric");
                if(resolved) break;coarse=fine;
            }
        }
    }
}

void check_addition_enthalpy() {
    const auto parameters = addition_parameters();
    const auto record = addition_record();
    const auto coordinates = addition_coordinates();
    const HelmholtzModel model(parameters, coordinates.active_parameters(), record);
    const double amount = 0.2;
    std::vector<AdditionVariant> variants;
    variants.push_back(addition_variant("fixed_system", amount, {}));
    variants.push_back(addition_variant("fixed_source", amount, 5.0e6));
    variants.push_back(addition_variant("free_pressure", amount, {}, true));
    variants.push_back(addition_variant("two_phase", amount, {}, true, 0.3));
    variants.push_back(addition_variant("unequal_temperature", amount, {}, false, {}, true));

    const double steps[] = {1.0e-4, 0.05, 5.0e-3, 1.0, 5.0e-3, 2.0e-5};
    const char* names[] = {"kij", "epsilon_tied", "reaction_a", "reaction_b", "reaction_c", "reaction_d"};
    for (auto& variant : variants) {
        auto center = predict(model, variant.observation);
        CHECK(center.status == ObservationStatus::Available && center.value && center.solved
            && center.final_solved && center.compiled && center.final_compiled);
        CHECK(center.solved->success && center.final_solved->success
            && center.solved->topology_event == equilibrium::TopologyEvent::None
            && center.final_solved->topology_event == equilibrium::TopologyEvent::None
            && center.solved->requested_tolerance_met.value_or(false)
            && center.final_solved->requested_tolerance_met.value_or(false)
            && center.solved->jacobian_rank == center.solved->unknown_count
            && center.final_solved->jacobian_rank == center.final_solved->unknown_count);
        variant.prediction = center;
        const auto start = coordinates.start(parameters, {variant.observation}, record);
        const auto exact = differentiate_group(model, {&variant.observation}, center, coordinates).front();
        CHECK(exact.status == ObservationStatus::Available && exact.value.size() == 7);
        variant.gradient = exact.value;

        for (std::size_t column = 0; column < 6; ++column) {
            const auto value = [&](double shift) {
                auto u = start;
                u[column] += shift;
                auto changed_observations = std::vector<Observation>{variant.observation};
                coordinates.apply(changed_observations, u);
                const HelmholtzModel rebuilt(coordinates.apply(parameters, u), coordinates.active_parameters(), record);
                const auto rebuilt_prediction = predict(rebuilt, changed_observations.front());
                if (rebuilt_prediction.status != ObservationStatus::Available || !rebuilt_prediction.value)
                    throw std::runtime_error("addition-f3 " + variant.name + " " + names[column] + ": "
                        + observation_status_name(rebuilt_prediction.status) + " " + rebuilt_prediction.message);
                return *rebuilt_prediction.value;
            };
            auto rungs = ladder(value, steps[column]);
            const double action = addition_coordinate_gradient(exact.value, column);
            const auto inputs = fields({field("variant", variant.name), field("T_initial_K", 260.0),
                field("T_final_K", variant.observation.final_problem.temperature.value), field("step", steps[column]),
                field("gradient", exact.value)});
            const std::string id = "addition-f3/" + variant.name + "/" + names[column];
            ladder_agrees(id, "d_delta_h_add_dtheta", action, rungs, inputs,
                "the exact finite-addition gradient agrees with rebuilt-model prediction re-solves", "addition-f3");
            if (column == 1 || column == 3) {
                const double fine_tolerance = FIRST_ORDER * std::max(1.0, std::abs(rungs.fine));
                double omitted = 0.0;
                if (column == 1) {
                    const HelmholtzModel pure(properties_support::ethane(),
                        {{ParameterFamily::DispersionEnergyOverK, {0}}});
                    const auto pure_action = [&](double T, double P) {
                        const PressureActionDirection direction{0.0, 0.0, {0.0},
                            {{ParameterFamily::DispersionEnergyOverK, {0}, 1.0}}};
                        const auto batch = pure.pressure_actions(T, P, {1.0}, DensityBranch::Liquid,
                            {PropertyObservable::ResidualEnthalpy}, {direction});
                        CHECK(batch.density_status == PropertyActionStatus::Available && batch.results.size() == 1
                            && batch.results[0].status == PropertyActionStatus::Available && batch.results[0].action);
                        return *batch.results[0].action;
                    };
                    const double P1 = variant.observation.problem.reactions[0].reference->reference_pressure_pa
                        .value_or(*center.solved->pressure);
                    const double P2 = variant.observation.final_problem.reactions[0].reference->reference_pressure_pa
                        .value_or(*center.final_solved->pressure);
                    omitted = (total_b_amount(*center.final_solved) * pure_action(*center.final_solved->temperature, P2)
                        - total_b_amount(*center.solved) * pure_action(*center.solved->temperature, P1)) / amount;
                } else {
                    omitted = -GAS_CONSTANT_J_PER_MOL_K
                        * (total_b_amount(*center.final_solved) - total_b_amount(*center.solved)) / amount;
                }
                const std::string direct_name = column == 1 ? "C_Ttheta_borne_term" : "reaction_b_endpoint_term";
                truth(id + "/omission_discrimination", "the independently calculated omitted endpoint term exceeds 1e3 times the fine criterion",
                    std::isfinite(omitted) && std::abs(omitted) > 1.0e3 * fine_tolerance,
                    fields({field("direct_term", omitted), field("fine_tolerance", fine_tolerance), field("term", direct_name)}),
                    "addition-f3");
            }
        }
        if (variant.observation.problem.temperature.value == variant.observation.final_problem.temperature.value) {
            const double temperature = variant.observation.problem.temperature.value;
            const double tref = variant.observation.problem.reactions[0].correlation->reference_temperature;
            const double G = exact.value[3];
            // g_b = G/T - D/T^2, so recover the source-slope coefficient D from a and b.
            const double D = temperature * G - temperature * temperature * exact.value[4];
            agree("addition-f3/" + variant.name + "/basis_c", "d_delta_h_add_dc", "J/mol",
                exact.value[5], G * std::log(temperature / tref) + D / temperature, 0.0, 1.0e-10,
                fields({field("T_K", temperature), field("Tref_K", tref), field("G", G), field("D", D)}),
                "the c column follows the source reaction basis identity before coordinate scaling", "addition-f3");
            agree("addition-f3/" + variant.name + "/basis_d", "d_delta_h_add_dd", "J/mol",
                exact.value[6], G * temperature + D, 0.0, 1.0e-10,
                fields({field("T_K", temperature), field("G", G), field("D", D)}),
                "the d column follows the source reaction basis identity before coordinate scaling", "addition-f3");
        } else {
            const double initial_b = total_b_amount(*center.solved);
            const double initial_c_term = -GAS_CONSTANT_J_PER_MOL_K * 260.0 * initial_b / amount;
            truth("addition-f3/unequal_temperature/initial_c_term", "c has a nonzero initial-endpoint caloric term at Tref with no state motion",
                std::isfinite(initial_c_term) && initial_c_term != 0.0,
                fields({field("phi_c", 0.0), field("phi_c_prime_per_K", 1.0 / 260.0),
                    field("initial_caloric_term_J_per_mol", initial_c_term), field("initial_n_B_mol", initial_b)}),
                "addition-f3");
        }
    }

    auto& reference = variants.front();
    const auto first_h = endpoint_enthalpy(model, reference.observation.problem,
        *reference.prediction.compiled, *reference.prediction.solved);
    const auto second_h = endpoint_enthalpy(model, reference.observation.final_problem,
        *reference.prediction.final_compiled, *reference.prediction.final_solved);
    auto shifted_record = record;
    std::get<IdealPolynomial>(shifted_record.component_intervals[0][0].correlation.form).reference_enthalpy = 1.0e4;
    const HelmholtzModel shifted_model(parameters, coordinates.active_parameters(), shifted_record);
    const auto shifted_prediction = predict(shifted_model, reference.observation);
    CHECK(shifted_prediction.status == ObservationStatus::Available && shifted_prediction.value);
    agree("addition-f4/value", "delta_h_add", "J/mol", *shifted_prediction.value, *reference.prediction.value,
        0.0, 1.0e-9, fields({field("anchor_h_J_per_mol", 1.0e4)}),
        "a common record enthalpy shift cancels from the finite addition", "addition-f4");
    const auto base_gradient = differentiate_group(model, {&reference.observation}, reference.prediction, coordinates).front();
    const auto shifted_gradient = differentiate_group(shifted_model, {&reference.observation}, shifted_prediction, coordinates).front();
    CHECK(base_gradient.status == ObservationStatus::Available && shifted_gradient.status == ObservationStatus::Available
        && base_gradient.value.size() == shifted_gradient.value.size());
    for (std::size_t column = 0; column < coordinates.coordinates.size(); ++column)
        agree("addition-f4/gradient/" + std::to_string(column), "d_delta_h_add_dtheta", "native",
            addition_coordinate_gradient(shifted_gradient.value, column),
            addition_coordinate_gradient(base_gradient.value, column), 0.0, 1.0e-9,
            fields({field("anchor_h_J_per_mol", 1.0e4), field("column", column)}),
            "a common record enthalpy shift preserves every gradient column", "addition-f4");
    const double condition = (std::abs(first_h.enthalpy) + std::abs(second_h.enthalpy))
        / (amount * std::abs(*reference.prediction.value));
    recording_support::observe("addition-f4/conditioning", "enthalpy_cancellation_condition", "1", condition,
        std::nullopt, fields({field("H_initial_J", first_h.enthalpy), field("H_final_J", second_h.enthalpy),
            field("n_add_mol", amount), field("delta_h_add_J_per_mol", *reference.prediction.value)}),
        "Observation", "endpoint_enthalpy_ratio", "the reference cancellation conditioning ratio", "addition-f4");
    truth("addition-f4/conditioning_limit", "the endpoint enthalpy cancellation ratio is at most 1e4",
        std::isfinite(condition) && condition <= 1.0e4, field("condition", condition), "addition-f4");

    for (const double increment : {0.2, 2.0e-3}) {
        auto initial = addition_problem(260.0, 3.0e6, {2.0, 0.2});
        auto final = addition_problem(260.0, 3.0e6, {2.0 + increment, 0.2});
        auto row = addition_observation(initial, final, increment);
        const auto prediction = predict(model, row);
        CHECK(prediction.status == ObservationStatus::Available && prediction.value && prediction.solved
            && prediction.final_solved && prediction.compiled && prediction.final_compiled);
        equilibrium::SolvedStateObservable activity;
        activity.kind = equilibrium::SolvedStateObservableKind::PhaseLogActivity;
        activity.phase = 0;
        activity.component = 0;
        equilibrium::SolvedStateActionDirection thermal;
        thermal.temperature_k = 1.0;
        thermal.feed_amounts_mol.assign(2, 0.0);
        thermal.active_parameters.assign(model.active_parameters().size(), 0.0);
        const auto action = equilibrium::solved_state_first_actions(*prediction.compiled, *prediction.solved,
            {activity}, {thermal}).front();
        CHECK(action.results.size() == 1 && action.results[0].status == equilibrium::SolvedStateActionStatus::Available
            && action.results[0].action);
        const HelmholtzModel pure(properties_support::ethane());
        const double h_feed = pure.state_at_pressure(260.0, 1.0e5, {1.0}, DensityBranch::Vapor).residual_enthalpy;
        const double q0 = -GAS_CONSTANT_J_PER_MOL_K * 260.0 * 260.0 * *action.results[0].action
            - GAS_CONSTANT_J_PER_MOL_K * 260.0 - h_feed;
        const auto h_initial = endpoint_enthalpy(model, row.problem, *prediction.compiled, *prediction.solved);
        const auto h_final = endpoint_enthalpy(model, row.final_problem, *prediction.final_compiled, *prediction.final_solved);
        const double roundoff = 64.0 * EPS * (std::abs(h_initial.enthalpy) + std::abs(h_final.enthalpy)) / increment;
        reference_scale_agrees("addition-f5/" + number(increment), "finite_addition_limit", "J/mol", *prediction.value, q0,
            roundoff, 1.0e-8, fields({field("n_add_mol", increment), field("D_n", *prediction.value),
                field("q0", q0), field("H_initial_J", h_initial.enthalpy), field("H_final_J", h_final.enthalpy),
                field("dln_activity_dT_per_K", *action.results[0].action), field("h_feed_residual_J_per_mol", h_feed)}),
            "each finite quotient equals the Gibbs-Helmholtz limit for the linear case", "addition-f5");
    }

    FitRequest recovery;
    recovery.controls = budgeted_controls();
    auto truth_parameters = parameters;
    truth_parameters.dispersion_energy_over_k = {191.42, 191.42};
    const HelmholtzModel truth_model(truth_parameters, {}, record);
    recovery.parameters = parameters;
    recovery.parameters.dispersion_energy_over_k = {195.2484, 195.2484};
    recovery.thermochemistry = record;
    recovery.coordinates.coordinates = {
        {{{ParameterFamily::DispersionEnergyOverK, {0}}, {ParameterFamily::DispersionEnergyOverK, {1}}}, 191.42, 10.0, 150.0, 250.0},
        reaction_coordinate("2A=B", ReactionCoefficient::B, 260.0)};
    recovery.coordinates.coordinates[1].origin = 200.0;
    recovery.coordinates.coordinates[1].lower = -2000.0;
    recovery.coordinates.coordinates[1].upper = 2000.0;
    for (const double T : {255.0, 260.0, 265.0, 270.0}) {
        auto initial = addition_problem(T, 3.0e6, {2.0, 0.2});
        auto final = addition_problem(T, 3.0e6, {2.2, 0.2});
        auto row = addition_observation(initial, final, 0.2, T);
        row.scale = 100.0;
        const auto truth_prediction = predict(truth_model, row);
        CHECK(truth_prediction.status == ObservationStatus::Available && truth_prediction.value);
        row.observed = *truth_prediction.value;
        row.problem.reactions[0].correlation->b = 180.0;
        row.final_problem.reactions[0].correlation->b = 180.0;
        recovery.observations.push_back(std::move(row));
    }
    const HelmholtzModel start_model(recovery.parameters, recovery.coordinates.active_parameters(), record);
    for (std::size_t row_index = 0; row_index < recovery.observations.size(); ++row_index) {
        const auto& row = recovery.observations[row_index];
        const auto prediction = predict(start_model, row);
        CHECK(prediction.status == ObservationStatus::Available);
        const auto gradient = differentiate_group(start_model, {&row}, prediction, recovery.coordinates).front();
        CHECK(gradient.status == ObservationStatus::Available && gradient.value.size() == 3);
    }
    const auto recovered = fit(recovery);
    truth("addition-f6/converged", "the noise-free known-parameter fit converges",
        recovered.status == "converged" && recovered.usable,
        fields({summary(recovered), field("trial_failures", recovered.trial_failures)}), "addition-f6");
    for (std::size_t row = 0; row < recovery.observations.size(); ++row)
        recording_support::observe(("addition-f6/data/" + std::to_string(row)).c_str(), "delta_h_add", "J/mol",
            recovered.predictions[row], recovery.observations[row].observed,
            fields({field("T_K", recovery.observations[row].temperature), field("P_Pa", 3.0e6),
                field("scale_J_per_mol", recovery.observations[row].scale)}),
            "Observation", "known_parameter_recovery", "synthetic datum and fitted prediction", "addition-f6");
    CHECK(recovered.physical.size() == 3);
    agree("addition-f6/epsilon_s0", "epsilon_over_k", "K", recovered.physical[0], 191.42, 0.0, 1.0e-6,
        summary(recovered), "the tied epsilon coordinate recovers the known value on s0", "addition-f6");
    agree("addition-f6/epsilon_s1", "epsilon_over_k", "K", recovered.physical[1], 191.42, 0.0, 1.0e-6,
        summary(recovered), "the tied epsilon coordinate recovers the known value on s1", "addition-f6");
    agree("addition-f6/reaction_b", "reaction_b", "K", recovered.physical[2], 200.0, 0.0, 1.0e-6,
        summary(recovered), "the reaction b coordinate recovers the known value", "addition-f6");
    agree("addition-f6/cost", "cost", "1", recovered.final_cost, 0.0, 1.0e-14, 0.0,
        summary(recovered), "the noise-free residual cost is at most 1e-14", "addition-f6");
    truth("addition-f6/covariance", "the fit has an available full-rank covariance",
        recovered.covariance.status == "available" && recovered.covariance.rank == 2,
        fields({field("covariance", recovered.covariance.status), field("rank", recovered.covariance.rank),
            field("singular_values", recovered.covariance.singular_values)}), "addition-f6");

    auto two_phase = std::find_if(variants.begin(), variants.end(), [](const AdditionVariant& item) {
        return item.name == "two_phase";
    });
    CHECK(two_phase != variants.end());
    const auto pure = HelmholtzModel(properties_support::ethane());
    const double oracle_initial = addition_oracle_endpoint(pure, *two_phase->prediction.solved,
        *two_phase->observation.problem.reactions[0].correlation);
    const double oracle_final = addition_oracle_endpoint(pure, *two_phase->prediction.final_solved,
        *two_phase->observation.final_problem.reactions[0].correlation);
    const auto public_initial = endpoint_enthalpy(model, two_phase->observation.problem,
        *two_phase->prediction.compiled, *two_phase->prediction.solved);
    const auto public_final = endpoint_enthalpy(model, two_phase->observation.final_problem,
        *two_phase->prediction.final_compiled, *two_phase->prediction.final_solved);
    for (const auto& [label, oracle, public_value] : {
            std::tuple<const char*, double, double>{"initial", oracle_initial, public_initial.enthalpy},
            {"final", oracle_final, public_final.enthalpy}}) {
        reference_scale_agrees(std::string("addition-f8/") + label + "/endpoint_enthalpy", "H", "J", public_value, oracle,
            1.0e-6, 1.0e-9, fields({field("oracle_H_J", oracle), field("public_H_J", public_value)}),
            "each endpoint total enthalpy agrees with the independent pure-fluid oracle", "addition-f8");
    }
    const auto& solved_initial = *two_phase->prediction.solved;
    const auto& solved_final = *two_phase->prediction.final_solved;
    const double hfeed = pure.state_at_pressure(260.0, 1.0e5, {1.0}, DensityBranch::Vapor).residual_enthalpy
        + 52.0 * (260.0 - 300.0) + 0.06 * (260.0 - 300.0) * (260.0 - 300.0);
    const double oracle_delta = (oracle_final - oracle_initial) / amount - hfeed;
    const double value_tolerance = 1.0e-9 * (std::abs(oracle_initial) + std::abs(oracle_final)) / amount + 1.0e-6;
    agree("addition-f8/value", "delta_h_add", "J/mol", *two_phase->prediction.value, oracle_delta,
        value_tolerance, 0.0, fields({field("H_initial_J", oracle_initial), field("H_final_J", oracle_final),
            field("h_feed_J_per_mol", hfeed), field("delta_h_add_J_per_mol", *two_phase->prediction.value)}),
        "the total-system finite addition agrees with the independent endpoint oracle", "addition-f8");
    const double h_initial_vapor = addition_oracle_phase_enthalpy(pure, solved_initial,
        *two_phase->observation.problem.reactions[0].correlation, 1);
    const double h_final_vapor = addition_oracle_phase_enthalpy(pure, solved_final,
        *two_phase->observation.final_problem.reactions[0].correlation, 1);
    const double vapor_contribution = std::abs(solved_final.amounts[1] * h_final_vapor
        - solved_initial.amounts[1] * h_initial_vapor) / amount;
    truth("addition-f8/vapor_discrimination", "the changed vapor contribution exceeds 1e3 times the value tolerance",
        std::isfinite(vapor_contribution) && vapor_contribution > 1.0e3 * value_tolerance,
        fields({field("vapor_contribution_J_per_mol", vapor_contribution), field("delta_h_tolerance_J_per_mol", value_tolerance),
            field("nV_initial_mol", solved_initial.amounts[1]), field("nV_final_mol", solved_final.amounts[1]),
            field("hV_initial_J_per_mol", h_initial_vapor), field("hV_final_J_per_mol", h_final_vapor)}), "addition-f8");

    auto failure = variants.front().observation;
    failure.temperature = 200.0;
    failure.pressure = 1.0e10;
    failure.branch = DensityBranch::Liquid;
    auto failed_feed = predict(model, failure);
    truth("addition-f7/feed_density", "the declared feed root failure is typed and names feed",
        failed_feed.status == ObservationStatus::DensityUnavailable && !failed_feed.value
            && failed_feed.message.find("feed:") != std::string::npos,
        fields({field("status", observation_status_name(failed_feed.status)), field("message", failed_feed.message)}),
        "addition-f7");
    auto missing_feed = variants.front().observation;
    missing_feed.composition = {0.0, 1.0};
    missing_feed.final_problem.feed.values = {2.0, 0.4};
    const auto missing = predict(model, missing_feed);
    truth("addition-f7/feed_record", "the balanced feed species without an ideal record is typed and names feed",
        missing.status == ObservationStatus::ThermochemistryUnavailable && !missing.value
            && missing.message.find("feed:") != std::string::npos,
        fields({field("status", observation_status_name(missing.status)), field("message", missing.message),
            field("feed_B_mol", 1.0), field("n_add_mol", amount)}), "addition-f7");
    ThermochemistryRecord no_records;
    no_records.standard_pressure = 1.0e5;
    no_records.component_intervals = {{}, {}};
    const HelmholtzModel undetermined_model(parameters, coordinates.active_parameters(), no_records);
    const auto undetermined = predict(undetermined_model, variants.front().observation);
    truth("addition-f7/initial_record", "undetermined endpoint calorics are typed and name initial",
        undetermined.status == ObservationStatus::ThermochemistryUnavailable && !undetermined.value
            && undetermined.message.find("initial:") != std::string::npos,
        fields({field("status", observation_status_name(undetermined.status)), field("message", undetermined.message)}),
        "addition-f7");
    auto mismatch = variants.front().observation;
    mismatch.final_problem.feed.values[0] += 0.01;
    const auto mismatch_message = rejection([&] { validate(mismatch, model); });
    truth("addition-f7/material_balance", "a mismatched component balance is rejected by component name",
        mismatch_message.find("component") != std::string::npos,
        fields({field("message", mismatch_message), field("delta_A_mol", 0.01)}), "addition-f7");
}

#include "reactive_cp_regression_checks.hpp"

}  // namespace regression_checks

namespace regression_checks {
// #484: van 't Hoff coefficient motion, supplied-record preservation and decadic replay.
void check_reaction_enthalpy_shift() {
    FitRequest request;
    request.parameters = properties_support::neutral_ternary();
    request.controls = budgeted_controls();
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    problem.reactions[0].equilibrium_constant.reset();
    problem.reactions[1].equilibrium_constant.reset();
    problem.reactions[0].correlation = ReactionLogPolynomial{0.4, 80.0, 0.5, 0.001, 300.0, 250.0, 400.0,
        equilibrium::EOS_STANDARD_STATE_ID};
    // A supplied reciprocal coefficient already shifted from a manufactured earlier record.
    const double log10 = std::log(10.0);
    problem.reactions[1].correlation = ReactionLogPolynomial{-9.004 * log10, -(2500.0 + 117.0) * log10,
        0.0, 0.009 * log10, 300.0, 250.0, 400.0, equilibrium::EOS_STANDARD_STATE_ID};
    for (const auto& reaction : problem.reactions) {
        Coordinate c;
        c.reaction = ReactionCoordinate{reaction.name, ReactionCoefficient::EnthalpyShift, 313.15};
        c.origin = 100.0; c.scale = 1000.0; c.lower = -5000.0; c.upper = 5000.0;
        request.coordinates.coordinates.push_back(c);
        for (const double t : {293.15, 313.15, 353.15, 393.15}) {
            auto p = problem; p.temperature = FreeValue::fixed(t);
            request.observations.push_back(log_constant(p, reaction.name.c_str()));
        }
    }
    const auto start = request.coordinates.start(request.parameters, request.observations);
    CHECK(start == std::vector<double>({-0.1, -0.1}));
    for (const double h : {0.0, -25.0, 25.0, -2500.0, 2500.0, -5000.0, 5000.0}) {
        auto rows = request.observations;
        const auto laws = request.coordinates.apply(rows, {(h-100.0)/1000.0, (h-100.0)/1000.0});
        for (const auto& reaction : problem.reactions) {
            const auto& base = *reaction.correlation;
            const auto& law = laws.at(reaction.name);
            const auto key = "reaction-enthalpy/" + reaction.name + "/" + std::to_string(h);
            agree(key + "/pivot", "ln_K", "1", reaction_log_value(law, 313.15), reaction_log_value(base, 313.15),
                2e-12, 0.0, "pivot=313.15 K", "the supplied correlation is preserved at the declared pivot", "issue-484");
            CHECK(law.c == base.c && law.d == base.d && law.reference_temperature == base.reference_temperature
                && law.temperature_min == base.temperature_min && law.temperature_max == base.temperature_max
                && law.standard_state_id == base.standard_state_id);
            if (h == 0.0) {
                for (const auto member:{&ReactionLogPolynomial::a,&ReactionLogPolynomial::b,&ReactionLogPolynomial::c,&ReactionLogPolynomial::d}) {
                    const double x=law.*member,y=base.*member;
                    CHECK(std::memcmp(&x,&y,sizeof(double))==0);
                }
            }
            const auto direction = reaction_direction(*request.coordinates.coordinates[0].reaction);
            const ReactionLogPolynomial action{direction[0], direction[1], direction[2], direction[3], 300.0, 250.0, 400.0, {}};
            for (const double t : {293.15, 313.15, 353.15, 393.15}) {
                agree(key + "/enthalpy/" + std::to_string(t), "d_h_reaction/d_h", "1",
                    GAS_CONSTANT_J_PER_MOL_K * t * t * reaction_log_value(action, t, 1), 1.0,
                    1e-12, 0.0, "constant source basis", "the analytic mixed thermal action changes reaction enthalpy by h", "issue-484");
                const double a = t*t*reaction_log_value(law, t, 1), b = t*t*reaction_log_value(base, t, 1);
                const double tolerance = 32*std::numeric_limits<double>::epsilon()*std::max({std::abs(a),std::abs(b),std::abs(h/GAS_CONSTANT_J_PER_MOL_K)});
                agree(key + "/slope/" + std::to_string(t), "delta_dlnK_dinverseT", "K", b-a, -h/GAS_CONSTANT_J_PER_MOL_K,
                    tolerance, 0.0, "propagated subtraction round-off", "the inverse-temperature slope changes by -h/R", "issue-484");
                if (law.c == 0.0) {
                    const auto source = reaction_log10_reciprocal(law);
                    agree(key + "/replay/" + std::to_string(t), "ln_K", "1", -log10*(source[0]/t+source[1]+source[2]*t),
                        reaction_log_value(law,t), 2e-12, 0.0, "negative base-10 reciprocal polynomial",
                        "both fitted coefficients replay the supplied correlation", "issue-484");
                }
            }
        }
    }
    // Manufacture values at nonzero shifts, keeping the original law as the fit's offset.
    auto shifted = request.observations;
    request.coordinates.apply(shifted, {1.1, -0.8});  // 1200 and -700 J/mol
    for (std::size_t i=0; i<shifted.size(); ++i)
        request.observations[i].observed = *predict(HelmholtzModel(request.parameters), shifted[i]).value;
    const auto result = fit(request);
    CHECK(result.usable && result.covariance.rank == 2);
    for (std::size_t i=0; i<2; ++i)
        agree("reaction-enthalpy/recovery/"+std::to_string(i), "h", "J/mol", result.physical[i], i==0 ? 1200.0 : -700.0,
            1e-8, 1e-10, "noise-free source observations at four temperatures", "one physical shift is recovered per reaction", "issue-484");
    CHECK(result.chart.reaction_laws.size() == 2);
    auto joint = request;
    Coordinate eos{{{ParameterFamily::SegmentCount, {0}}}, request.parameters.segment_count[0], 0.1};
    eos.lower = eos.origin-0.5; eos.upper = eos.origin+0.5;
    joint.coordinates.coordinates.insert(joint.coordinates.coordinates.begin()+1, eos);
    joint.observations.push_back(density(300.0, 1e5, {1.0,0.0,0.0}, DensityBranch::Vapor));
    const auto truth_point = std::vector<double>{1.1, 0.2, -0.8};
    auto truth_rows = joint.observations;
    joint.coordinates.apply(truth_rows,truth_point);
    const HelmholtzModel truth_model(joint.coordinates.apply(joint.parameters,truth_point));
    for (std::size_t i=0;i<truth_rows.size();++i)
        joint.observations[i].observed = *predict(truth_model,truth_rows[i]).value;
    const auto joint_fit = fit(joint);
    CHECK(joint_fit.usable && joint_fit.covariance.rank == 3);
    for (std::size_t i=0;i<3;++i)
        agree("reaction-enthalpy/joint/"+std::to_string(i), "optimizer_coordinate", "1", joint_fit.optimizer[i], truth_point[i],
            1e-8, 1e-10, "interleaved EOS and reaction coordinates", "the affine chain applies each coordinate scale once", "issue-484");
    bool refused=false;
    try { (void)reaction_log10_reciprocal(*problem.reactions[0].correlation); }
    catch (const std::invalid_argument&) { refused=true; }
    CHECK(refused);
    for (const auto member : {&ReactionLogPolynomial::e, &ReactionLogPolynomial::f}) {
        auto changed = problem.reactions[1];
        changed.correlation.value().*member = 1.0;
        CHECK(!same_reaction(changed, problem.reactions[1]));
        refused = false;
        try { (void)reaction_log10_reciprocal(*changed.correlation); }
        catch (const std::invalid_argument&) { refused = true; }
        CHECK(refused);
    }
    refused=false;
    try { (void)fit_sequence(request, {start}, 10, 10.0); } catch (const std::invalid_argument& e) { refused=std::string(e.what()).find("fit_sequence") != std::string::npos; }
    CHECK(refused);
    request.minimax = MinimaxRequest{};
    CHECK(fit(request).message.find("not supported in minimax") != std::string::npos);
    CHECK(fit_minimax(request, start).message.find("not supported in minimax") != std::string::npos);
    request.minimax.reset();
    auto wrong = request;
    wrong.coordinates.coordinates[0].reaction->pivot_temperature = 0.0;
    CHECK(fit(wrong).status == "invalid_setup");
    for (const auto coefficient : {ReactionCoefficient::A, ReactionCoefficient::B}) {
        auto pair = request;
        pair.coordinates.coordinates.resize(1);
        auto plain = reaction_coordinate("reaction_a", coefficient, coefficient == ReactionCoefficient::A ? 1.0 : 80.0);
        plain.origin = coefficient == ReactionCoefficient::A ? problem.reactions[0].correlation->a : problem.reactions[0].correlation->b;
        pair.coordinates.coordinates.push_back(plain);
        for (const bool reverse : {false, true}) {
            if (reverse) std::reverse(pair.coordinates.coordinates.begin(), pair.coordinates.coordinates.end());
            auto u = pair.coordinates.start(pair.parameters, pair.observations);
            const std::size_t h = reverse ? 1 : 0, q = 1-h;
            u[q] += 0.4;
            auto unshifted = pair.observations;
            const auto base = pair.coordinates.apply(unshifted, u).at("reaction_a");
            u[h] = 2.4; // 2500 J/mol with the inherited affine origin.
            auto rows = pair.observations;
            const auto law = pair.coordinates.apply(rows, u).at("reaction_a");
            agree("reaction-pair/pivot", "ln_K", "1", reaction_log_value(law, 313.15), reaction_log_value(base, 313.15),
                2e-12, 0.0, "same plain value", "the shift preserves the pivot after the plain assignment", "issue-557");
            for (const double t : {293.15, 313.15, 353.15, 393.15})
                agree("reaction-pair/enthalpy", "delta_H/h", "1",
                    GAS_CONSTANT_J_PER_MOL_K*t*t*(reaction_log_value(law,t,1)-reaction_log_value(base,t,1))/2500.0,
                    1.0, 1e-12, 0.0, "applied native laws", "the completed law adds the physical enthalpy increment", "issue-557");
            const HelmholtzModel model(pair.parameters);
            for (std::size_t i=0; i<rows.size(); ++i) {
                const auto prediction = predict(model, rows[i]);
                CHECK(prediction.status == ObservationStatus::Available);
                pair.observations[i].observed = *prediction.value;
            }
            const auto recovered = fit(pair);
            CHECK(recovered.usable && recovered.covariance.rank == 2);
            for (std::size_t c=0; c<2; ++c)
                agree("reaction-pair/recovery", "u", "1", recovered.optimizer[c], u[c], 1e-8, 1e-10,
                    "noise-free log K at four temperatures", "the pair is recovered jointly", "issue-557");
            auto replay = pair.observations;
            for (auto& row : replay) for (auto& reaction : row.problem.reactions)
                if (reaction.name == "reaction_a") reaction.correlation = recovered.chart.reaction_laws.at(reaction.name);
            for (std::size_t i=0; i<replay.size(); ++i)
                CHECK(*predict(model, replay[i]).value == recovered.predictions[i]);
            auto duplicate = pair;
            duplicate.coordinates.coordinates.push_back(plain);
            CHECK(fit(duplicate).message.find("duplicate target") != std::string::npos);
            duplicate = pair;
            duplicate.coordinates.coordinates.push_back(request.coordinates.coordinates[0]);
            duplicate.coordinates.coordinates.back().reaction->pivot_temperature = 300.0;
            CHECK(fit(duplicate).message.find("duplicate target") != std::string::npos);
        }
    }
    auto preapplied = request;
    const std::vector<double> offset{-8201.9, 8418.5};
    for (auto& coordinate : preapplied.coordinates.coordinates) { coordinate.lower = -10000.0; coordinate.upper = 10000.0; }
    preapplied.coordinates.apply(preapplied.observations, {(offset[0]-100.0)/1000.0, (offset[1]-100.0)/1000.0});
    for (std::size_t c=0; c<2; ++c) {
        preapplied.coordinates.coordinates[c].lower -= offset[c];
        preapplied.coordinates.coordinates[c].upper -= offset[c];
    }
    const auto zero = preapplied.coordinates.start(preapplied.parameters, preapplied.observations);
    CHECK(preapplied.coordinates.physical(zero) == std::vector<double>({0.0, 0.0}));
    auto shifted_rows = preapplied.observations;
    const auto total_laws = preapplied.coordinates.apply(shifted_rows, {0.15, -0.3});
    for (std::size_t c=0; c<2; ++c)
        agree("reaction-enthalpy/preapplied", "total_h", "J/mol",
            GAS_CONSTANT_J_PER_MOL_K*313.15*313.15*(reaction_log_value(total_laws.at(problem.reactions[c].name),313.15,1)
                -reaction_log_value(*problem.reactions[c].correlation,313.15,1)), offset[c]+(c==0 ? 250.0 : -200.0),
            1e-8, 0.0, fields({field("pre_applied_J_mol", offset[c]), field("increment_J_mol", c==0 ? 250.0 : -200.0),
                field("law_a", total_laws.at(problem.reactions[c].name).a), field("law_b_K", total_laws.at(problem.reactions[c].name).b),
                field("lower_increment_J_mol", preapplied.coordinates.coordinates[c].lower), field("upper_increment_J_mol", preapplied.coordinates.coordinates[c].upper)}),
            "translated bounds retain the total physical shift", "issue-557");
    // #448's effective species-reference refusal reaches point readers and heat final endpoints.
    wrong = request; wrong.parameters = properties_support::aqueous_electrolyte();
    auto source = referenced_formula(std::nullopt);
    source.reactions[0].name = "reaction_a";
    source.reactions[0].reference->source_basis = ReferenceSourceBasis::RawMoleFractionInfiniteDilution;
    wrong.coordinates.coordinates.resize(1);
    wrong.observations = {solved(ObservationKind::PhaseMoleFraction, source, 0, {1})};
    wrong.observations.push_back(wrong.observations[0]);
    wrong.observations[1].kind = ObservationKind::AdditionEnthalpy;
    wrong.observations[1].final_problem = source;
    wrong.observations[1].final_problem.source_reference->species_reference = {SpeciesReference::PureGas,
        SpeciesReference::InfiniteDilution, SpeciesReference::InfiniteDilution};
    CHECK(fit(wrong).message.find("tied species references disagree") != std::string::npos);
    refused=false;
    try { FitResult errors; ObservationEvaluator reader(wrong, errors, {0,1}); }
    catch (const std::invalid_argument& e) { refused=std::string(e.what()).find("tied species references disagree") != std::string::npos; }
    CHECK(refused);
    wrong.observations[1].final_problem.source_reference->species_reference[0] = SpeciesReference::PureLiquid;
    CHECK(wrong.coordinates.start(wrong.parameters,wrong.observations).size() == 1);
}
}  // namespace regression_checks
