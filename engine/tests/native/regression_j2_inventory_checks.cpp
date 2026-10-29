// R2 J2, R3 and Cp recovery; hydration failures, minimax consistency, inventory boundaries and electrostatics.
#include "regression_checks_shared.hpp"

namespace regression_checks {

// ---- R2 J2 and R3: manufactured recovery, bounds, failures, weights, rank and covariance.
FitRequest recovery_request(const Truth& truth_case, double epsilon, double kij, double upper = 170.0) {
    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = properties_support::neutral_binary();
    request.parameters.dispersion_energy_over_k[0] = epsilon;
    request.parameters.k_ij[1] = request.parameters.k_ij[2] = kij;
    request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 150.0, 10.0, 140.0, upper},
                                       {{{ParameterFamily::Kij, {1, 2}}}, 0.0, 0.1, -0.2, 0.2}};
    request.observations = truth_case.observations;
    return request;
}

// Held, Cameretti and Sadowski (2008), Table 1 / Eq. (7), packet v3 neutral water subset.
Parameters water_parameters() {
    auto p = properties_support::parameters(1);
    p.segment_count = {1.204659};
    p.sigma_constant = {2.7927};
    p.dispersion_energy_over_k = {353.9449};
    p.molar_mass_kg_per_mol = {0.01801528};
    p.sigma_term_component = {0, 0};
    p.sigma_term_amplitude = {10.11, -1.417};
    p.sigma_term_exponent = {-0.01775, -0.01146};
    p.site_component = {0, 0};
    p.site_multiplicity = {1, 1};
    p.association_active = {0, 1, 1, 0};
    p.association_energy_over_k = {0.0, 2425.6714, 2425.6714, 0.0};
    p.association_volume = {0.0, 0.0450989, 0.0450989, 0.0};
    return p;
}

ThermochemistryRecord constant_cp(std::size_t components) {
    ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    record.component_intervals.assign(components, {{150.0, 500.0, true, true,
        {IdealPolynomial{{4.0 * GAS_CONSTANT_J_PER_MOL_K}, 300.0, 0.0, 0.0}, 1.0e5}}});
    return record;
}

CoordinateSet water_coordinates() {
    CoordinateSet set;
    set.coordinates = {
        {{{ParameterFamily::SegmentCount, {0}}}, 1.204659, 0.1, 0.8, 2.0},
        {{{ParameterFamily::SigmaConstant, {0}}}, 2.7927, 0.05, 2.4, 3.2},
        {{{ParameterFamily::DispersionEnergyOverK, {0}}}, 353.9449, 10.0, 200.0, 500.0},
        {{{ParameterFamily::AssociationEnergyOverK, {1, 2}}}, 2425.6714, 100.0, 1500.0, 3500.0},
        {{{ParameterFamily::AssociationVolume, {1, 2}}}, 0.0450989, 0.005, 0.005, 0.2}};
    return set;
}

void check_cp_derivatives() {
    struct Case { const char* name; Parameters parameters; double T, P; std::vector<double> x;
                  std::vector<ActiveParameter> active; std::vector<double> steps; };
    const auto water_active = water_coordinates().active_parameters();
    const std::vector<Case> cases = {
        {"a1", properties_support::ethane(), 200.0, 217800.6949469248, {1.0},
         {{ParameterFamily::SegmentCount, {0}}, {ParameterFamily::SigmaConstant, {0}},
          {ParameterFamily::DispersionEnergyOverK, {0}}}, {5e-4, 5e-4, 0.05}},
        {"a2-298", water_parameters(), 298.15, 1e5, {1.0}, water_active, {5e-4, 5e-4, 0.05, 1.0, 1e-5}},
        {"a2-400", water_parameters(), 400.0, 1e6, {1.0}, water_active, {5e-4, 5e-4, 0.05, 1.0, 1e-5}},
        {"a3", properties_support::associating_binary(), 313.15, 1e5, {0.25, 0.75},
         {{ParameterFamily::Kij, {1, 2}}, {ParameterFamily::SegmentCount, {0}},
          {ParameterFamily::AssociationEnergyOverK, {3, 12}},
          {ParameterFamily::AssociationVolume, {11, 14}}}, {1e-4, 5e-4, 1.0, 1e-5}}};
    for (const auto& c : cases) {
        const auto thermal = constant_cp(c.x.size());
        const HelmholtzModel model(c.parameters, c.active, thermal);
        auto item = density(c.T, c.P, c.x, DensityBranch::Liquid);
        item.kind = ObservationKind::IsobaricHeatCapacity;
        const auto value = predict(model, item);
        const auto gradient = differentiate(model, item, value);
        CHECK(value.status == ObservationStatus::Available && gradient.status == ObservationStatus::Available);
        const auto state = model.state_at_pressure(c.T, c.P, c.x, item.branch);
        const std::string id = std::string("cp-a/") + c.name;
        const auto inputs = fields({field("T_K", c.T), field("P_Pa", c.P), field("x", c.x)});
        agree(id + "/ideal-value", "ideal_cp", "J/(mol K)",
            *value.value - *state.residual_isobaric_heat_capacity, 4.0 * GAS_CONSTANT_J_PER_MOL_K,
            1e-12 * std::max(1.0, *value.value), 0.0, inputs,
            "the total prediction contains the known constant ideal contribution", "cp-a");
        double lift_ratio = 0.0;
        for (std::size_t k = 0; k < c.active.size(); ++k) {
            const auto& active = c.active[k];
            const auto closure = [&](double shift) {
                auto p = c.parameters;
                for (const auto index : active.indices) detail::family_values(p, active.family).first[index] += shift;
                return *predict(HelmholtzModel(p, {}, thermal), item).value;
            };
            const auto rungs = ladder(closure, c.steps[k]);
            const double richardson = (4.0 * rungs.fine - rungs.coarse) / 3.0;
            const auto column_inputs = fields({inputs, field("step", c.steps[k]), field("fd_signal", std::abs(rungs.fine)),
                                               field("richardson", richardson), field("fine_minus_coarse", rungs.fine - rungs.coarse),
                                               field("richardson_correction", richardson - rungs.fine),
                                               field("richardson_error", std::abs(gradient.value[k] - richardson))});
            const auto column = id + "/" + std::to_string(k);
            // Owner amendment: these two physical sensitivities lie below the native-unit signal floor.
            const bool weak = k == 2 && (std::string(c.name) == "a2-298" || std::string(c.name) == "a3");
            if (weak) {
                agree(column + "/fine", "dcp_dtheta", "native", gradient.value[k], rungs.fine,
                      FIRST_ORDER * std::abs(rungs.fine), 0.0, column_inputs,
                      "the weak physical action agrees within the stricter relative criterion", "cp-a");
                agree(column + "/coarse", "dcp_dtheta", "native", gradient.value[k], rungs.coarse,
                      4.0 * FIRST_ORDER * std::abs(rungs.fine), 0.0, column_inputs,
                      "the coarse relative rung carries at most four times the criterion", "cp-a");
                truth(column + "/resolved", "finite-difference refinement resolves the relative criterion",
                    std::isfinite(rungs.fine) && rungs.fine != 0.0
                    && std::abs(rungs.fine - rungs.coarse) < FIRST_ORDER * std::abs(rungs.fine)
                    && std::abs(richardson - rungs.fine) < FIRST_ORDER * std::abs(rungs.fine), column_inputs, "cp-a");
            } else {
                ladder_agrees(column, "dcp_dtheta", gradient.value[k], rungs, column_inputs,
                             "the fixed-pressure action agrees with rebuilt-model central differences", "cp-a");
            }
            agree(column + "/richardson", "dcp_dtheta", "native", gradient.value[k], richardson,
                  FIRST_ORDER * (weak ? std::abs(richardson) : std::max(1.0, std::abs(richardson))), 0.0, column_inputs,
                  "the action agrees with Richardson extrapolation", "cp-a");
            truth(column + "/signal", "the difference carries a nonzero parameter signal",
                  weak ? std::isfinite(rungs.fine) && rungs.fine != 0.0 : std::abs(rungs.fine) >= 1e-2,
                  column_inputs, "cp-a");
            const ParameterActionDirection parameter{active.family, active.indices, 1.0};
            const auto pressure = model.pressure_actions(c.T, c.P, c.x, item.branch,
                {PropertyObservable::TotalIsobaricHeatCapacity, PropertyObservable::ResidualIsobaricHeatCapacity},
                {{0.0, 0.0, std::vector<double>(c.x.size(), 0.0), {parameter}}});
            const auto fixed = model.density_actions(c.T, state.molar_density, c.x,
                {PropertyObservable::TotalIsobaricHeatCapacity},
                {{0.0, 0.0, std::vector<double>(c.x.size(), 0.0), {parameter}}});
            CHECK(pressure.results[0].action && pressure.results[1].action && fixed.results[0].action);
            agree(column + "/ideal-action", "dcp_dtheta", "native", *pressure.results[0].action,
                  *pressure.results[1].action, 0.0, 1e-12, column_inputs,
                  "fixed ideal thermochemistry carries no parameter action", "cp-a");
            lift_ratio = std::max(lift_ratio, std::abs(gradient.value[k] - *fixed.results[0].action)
                / (FIRST_ORDER * std::max(1.0, std::abs(rungs.fine))));
        }
        truth(id + "/density-lift", "at least one omitted density lift exceeds 1000 criteria",
              lift_ratio > 1e3, fields({inputs, field("lift_criterion_ratio", lift_ratio)}), "cp-a");
    }
    auto item = density(200.0, 217800.6949469248, {1.0}, DensityBranch::Liquid);
    item.kind = ObservationKind::IsobaricHeatCapacity;
    const auto missing = predict(HelmholtzModel(properties_support::ethane()), item);
    truth("cp-a/missing", "missing thermochemistry is typed and never a value",
        missing.status == ObservationStatus::ThermochemistryUnavailable && !missing.value,
        field("message", missing.message), "cp-a");
    item.pressure = 1e10;
    const auto absent = predict(HelmholtzModel(properties_support::ethane(), {}, constant_cp(1)), item);
    truth("cp-a/root", "an absent density root retains its typed cause",
        absent.status == ObservationStatus::DensityUnavailable && !absent.value && !absent.message.empty(),
        field("message", absent.message), "cp-a");
}

void check_cp_recovery() {
    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = water_parameters();
    request.coordinates = water_coordinates();
    request.thermochemistry = constant_cp(1);
    const HelmholtzModel known(request.parameters, {}, request.thermochemistry);
    for (double T : {300.0, 340.0, 380.0, 420.0}) {
        auto rho = density(T, 1e6, {1.0}, DensityBranch::Liquid);
        auto cp = rho;
        cp.kind = ObservationKind::IsobaricHeatCapacity;
        auto sat = solved(ObservationKind::EquilibriumPressure, pure_saturation(T, 1e5));
        for (auto item : {rho, cp, sat}) {
            const auto prediction = predict(known, item);
            CHECK(prediction.status == ObservationStatus::Available);
            item.observed = *prediction.value;
            request.observations.push_back(item);
        }
    }
    for (const auto& active : request.coordinates.active_parameters())
        for (const auto index : active.indices)
            detail::family_values(request.parameters, active.family).first[index] *=
                active.family == ParameterFamily::AssociationVolume ? 1.05 : 1.02;
    const auto result = fit(request);
    const auto inputs = fields({summary(result), field("singular_values", result.covariance.singular_values),
                               field("rank", result.covariance.rank)});
    truth("cp-b/converged", "five-coordinate recovery converges",
        result.status == "converged" && result.usable, inputs, "cp-b");
    for (std::size_t k = 0; k < 5; ++k)
        agree("cp-b/value/" + std::to_string(k), "parameter", "native", result.physical[k],
            request.coordinates.coordinates[k].origin, 0.0, 1e-6, inputs,
            "the known physical coordinate is recovered", "cp-b");
    agree("cp-b/cost", "cost", "1", result.final_cost, 0.0, 1e-14, 0.0, inputs,
          "the manufactured residual cost vanishes", "cp-b");
    truth("cp-b/rank", "the conditional covariance is available with rank five",
          result.covariance.status == "available" && result.covariance.rank == 5, inputs, "cp-b");
}

void check_j2() {
    const auto truth_case = manufactured();
    const auto result = fit(recovery_request(truth_case, 148.0, 0.0));
    truth("r2-j2/converged", "the noise-free recovery converges with no trial failure",
          result.status == "converged" && result.usable, summary(result), "r2-j2");
    agree("r2-j2/epsilon", "dispersion_energy_over_k", "K", result.physical[0], 152.0, 0.0, 1.0e-6, summary(result),
          "the manufactured epsilon is recovered", "r2-j2");
    agree("r2-j2/kij", "k_ij", "1", result.physical[1], -0.02, 1.0e-6, 0.0, summary(result), "the manufactured k_ij is recovered", "r2-j2");
    agree("r2-j2/final_cost", "cost", "1", result.final_cost, 0.0, 1.0e-14, 0.0, summary(result), "the final cost vanishes", "r2-j2");
    truth("r2-j2/identifiable", "the interior full-rank fit carries a covariance",
          result.covariance.status == "available" && result.covariance.rank == 2, summary(result), "r2-j2");
    // The bounded-stationarity measure max_c |u_c - clamp(u_c - g_c)|, with g = J^T r from the published columns.
    const auto projected_gradient = [](const FitRequest& request, const FitResult& result) {
        const auto [lower, upper] = request.coordinates.optimizer_bounds();
        const std::size_t width = result.optimizer.size();
        double projected = 0.0;
        for (std::size_t c = 0; c < width; ++c) {
            double gradient = 0.0;
            for (std::size_t i = 0; i < result.weighted_residuals.size(); ++i)
                gradient += result.optimizer_jacobian[i * width + c] * result.weighted_residuals[i];
            projected = std::max(projected, std::abs(result.optimizer[c] - std::clamp(result.optimizer[c] - gradient, lower[c], upper[c])));
        }
        return projected;
    };
    auto unfinished = recovery_request(truth_case, 148.0, 0.0);
    unfinished.controls.maximum_iterations = 2;
    const auto iterate = fit(unfinished);
    truth("r2-j2/not_converged", "an iteration-limited fit keeps its candidate and rank but no covariance",
          iterate.status == "limit_reached" && iterate.usable && iterate.covariance.status == "not_converged"
              && iterate.covariance.rank == 2 && iterate.covariance.optimizer.empty() && iterate.covariance.physical.empty()
              && !iterate.covariance.variance_factor.has_value() && iterate.final_cost > 1.0e-10, summary(iterate), "r2-j2");
    agree("r2-j2/not_converged/projected_gradient", "projected_gradient_max", "1",
          iterate.covariance.projected_gradient_norm, projected_gradient(unfinished, iterate), 0.0, 1.0e-10, summary(iterate),
          "the projected gradient is reported at the returned point", "r2-j2");
    CHECK(iterate.rank_records.size() == 1 && iterate.rank_records.front().rank == 2);
    CHECK(iterate.rank_records.front().dimension == 2);
    const auto other = fit(recovery_request(truth_case, 156.0, 0.05));
    agree("r2-j2/second_start/epsilon", "dispersion_energy_over_k", "K", other.physical[0], 152.0, 0.0, 1.0e-6, summary(other),
          "a second start recovers the same parameters", "r2-j2");
    agree("r2-j2/second_start/kij", "k_ij", "1", other.physical[1], -0.02, 1.0e-6, 0.0, summary(other),
          "a second start recovers the same parameters", "r2-j2");
    const auto bounded = fit(recovery_request(truth_case, 148.0, 0.0, 151.0));
    truth("r2-j2/active_bound", "a bound below the truth ends at the bound with covariance withheld",
          bounded.usable && bounded.physical[0] == 151.0 && bounded.covariance.status == "active_bound"
              && bounded.covariance.active_bounds == std::vector<std::size_t>{0}, summary(bounded), "r2-j2");
    // The same manufactured targets, starting on the upper epsilon bound with descent initially
    // inward. Coupling to k_ij turns that gradient outward without a new bound contact (#193).
    auto bound_start = recovery_request(truth_case, 151.0, -0.1, 151.0);
    bound_start.controls.maximum_iterations = 12;  // shared work limit; no iteration-count pass criterion
    const auto initial_epsilon_gradient = [](const FitRequest& start) {
        const HelmholtzModel model(start.parameters, start.coordinates.active_parameters());
        double gradient = 0.0;
        for (std::size_t row = 0; row < start.observations.size(); ++row) {
            const auto& observation = start.observations[row];
            const auto prediction = predict(model, observation);
            const auto action = differentiate(model, observation, prediction);
            if (!(prediction.value && prediction.residual && action.status == ObservationStatus::Available))
                throw std::runtime_error("initial epsilon gradient row " + std::to_string(row) + ": prediction "
                    + observation_status_name(prediction.status) + " " + prediction.message + "; action "
                    + observation_status_name(action.status) + " " + action.message);
            gradient += *prediction.residual * residual_slope(observation, *prediction.value)
                        * action.value[0] * start.coordinates.coordinates[0].scale;
        }
        return gradient;
    };
    CHECK(initial_epsilon_gradient(bound_start) > 0.0);  // descent from the upper bound is initially inward
    // Removing epsilon from the fit gives the independent scalar minimization at epsilon = 151 K.
    auto reduced = bound_start;
    reduced.coordinates.coordinates.erase(reduced.coordinates.coordinates.begin());
    const auto reference = fit(reduced);
    const auto refreshed = fit(bound_start);
    CHECK(reference.usable && refreshed.usable && refreshed.optimizer.size() == 2);
    agree("r2-j2/bound_start/cost", "cost", "1", refreshed.final_cost, reference.final_cost,
          1.0e-12, 1.0e-8, summary(refreshed),
          "the coupled fit reaches the independently reduced bounded minimum within the work limit", "r2-j2");
    agree("r2-j2/bound_start/kij", "k_ij", "1", refreshed.physical[1], reference.physical[0],
          1.0e-7, 0.0, summary(refreshed), "the free coordinate agrees with the scalar minimum", "r2-j2");
    const double projected = projected_gradient(bound_start, refreshed);
    truth("r2-j2/bound_start/feasible", "the accepted optimum is feasible with epsilon at its upper bound",
          bound_start.coordinates.within_bounds(refreshed.optimizer) && refreshed.physical[0] == 151.0,
          summary(refreshed), "r2-j2");
    agree("r2-j2/bound_start/projected_gradient", "projected_gradient_max", "1", projected, 0.0,
          1.0e-6, 0.0, summary(refreshed), "the full weighted gradient satisfies bounded stationarity", "r2-j2");
    // Initially held epsilon must be released as k_ij moves and the gradient turns inward.
    auto interior_truth = truth_case;
    interior_truth.parameters.dispersion_energy_over_k[0] = 150.0;
    const HelmholtzModel interior_model(interior_truth.parameters);
    for (auto& observation : interior_truth.observations) {
        const auto prediction = predict(interior_model, observation);
        if (!prediction.value)
            throw std::runtime_error("interior truth row " + std::to_string(&observation - interior_truth.observations.data())
                + ": " + observation_status_name(prediction.status) + " " + prediction.message);
        observation.observed = *prediction.value;
        // The default cold seed coalesces at the release start (+0.1); these measured two-phase
        // guesses admit the same T/feed/bounds/parameter start, without a topology-event bypass.
        if (prediction.solved) {
            const auto& state = *prediction.solved;
            for (std::size_t q = 0; q < observation.problem.phases.size(); ++q) {
                observation.problem.phases[q].composition_guess.assign(state.mole_fractions.begin() + q * 2,
                                                                        state.mole_fractions.begin() + (q + 1) * 2);
                observation.problem.phases[q].packing_guess = state.packing_fractions[q];
            }
            observation.problem.pressure = FreeValue::variable(*state.pressure);
        }
    }
    auto inward = recovery_request(interior_truth, 151.0, 0.1, 151.0);
    CHECK(initial_epsilon_gradient(inward) < 0.0);  // outward descent at the upper bound: initially held
    const auto released = fit(inward);
    CHECK(released.usable);
    agree("r2-j2/bound_start/inward", "dispersion_energy_over_k", "K", released.physical[0], 150.0,
          0.0, 1.0e-6, summary(released), "inward descent releases the held coordinate and recovers the interior truth", "r2-j2");
    agree("r2-j2/bound_start/inward_cost", "cost", "1", released.final_cost, 0.0,
          1.0e-14, 0.0, summary(released), "the released interior fit recovers the manufactured targets", "r2-j2");
    auto failing = recovery_request(truth_case, 148.0, 0.0);
    failing.observations.push_back(density(200.0, 1.0e10, {0.5, 0.5}, DensityBranch::Liquid));
    const auto failed = fit(failing);
    const std::vector<ObservationStatus> expected_statuses = {ObservationStatus::Available, ObservationStatus::Available,
        ObservationStatus::Available, ObservationStatus::Available, ObservationStatus::Available, ObservationStatus::DensityUnavailable};
    truth("r2-j2/initial_failure", "a failing observation at the start is retained with its cause and index",
          failed.status == "initial_evaluation_failed" && !failed.usable && failed.failures.size() == 1
              && failed.failures[0].observation == 5 && failed.failures[0].status == ObservationStatus::DensityUnavailable
              && failed.predictions.empty() && failed.statuses == expected_statuses, summary(failed), "r2-j2");
    // R4 rows: the public transport reproduces these native values for the same explicit inputs.
    const auto row = [&](const char* id, const char* unit, double value) {
        recording_support::observe(id, "native", unit, value, std::nullopt, summary(result), "Relationship", "regression",
                                   "native fit value for the public R4 comparison", "r4");
    };
    row("r4/recovery/epsilon", "K", result.physical[0]);
    row("r4/recovery/kij", "1", result.physical[1]);
    row("r4/recovery/final_cost", "1", result.final_cost);
    row("r4/recovery/initial_cost", "1", result.initial_cost);
    for (std::size_t i = 0; i < result.predictions.size(); ++i)
        row((std::string("r4/recovery/prediction/") + std::to_string(i)).c_str(), "native", result.predictions[i]);
    for (std::size_t k = 0; k < 4; ++k) {
        row((std::string("r4/recovery/covariance/") + std::to_string(k)).c_str(), "1", result.covariance.optimizer[k]);
        row((std::string("r4/recovery/physical_covariance/") + std::to_string(k)).c_str(), "1", result.covariance.physical[k]);
    }
    row("r4/recovery/variance_factor", "1", result.covariance.variance_factor.value_or(-1.0));
    for (std::size_t k = 0; k < result.optimizer_jacobian.size(); ++k)
        row((std::string("r4/recovery/optimizer_jacobian/") + std::to_string(k)).c_str(), "1", result.optimizer_jacobian[k]);
}

void check_r3() {
    const auto truth_case = manufactured();
    auto request = recovery_request(truth_case, 148.0, 0.0);
    const auto base = fit(request);
    request.weights.assign(request.observations.size(), 2.0);
    const auto doubled = fit(request);
    CHECK(base.usable && doubled.usable && base.covariance.status == "available" && doubled.covariance.status == "available");
    for (std::size_t c = 0; c < 2; ++c)
        agree("r3/weight_scale/parameter/" + std::to_string(c), "theta", "native", doubled.physical[c], base.physical[c], 0.0, 1.0e-9,
              summary(doubled), "doubling every weight leaves the fitted parameters unchanged", "r3");
    for (std::size_t k = 0; k < 4; ++k)
        agree("r3/weight_scale/covariance/" + std::to_string(k), "covariance", "1", doubled.covariance.optimizer[k],
              0.5 * base.covariance.optimizer[k], 0.0, 1.0e-6, summary(doubled), "doubling every weight halves the covariance", "r3");
    agree("r3/weight_scale/cost", "cost", "1", doubled.initial_cost, 2.0 * base.initial_cost, 0.0, 1.0e-12, summary(doubled),
          "doubling every weight doubles the cost", "r3");
    // A zero weight removes a row from the objective but not from the report.
    auto zero = recovery_request(truth_case, 148.0, 0.0);
    auto inconsistent = solved(ObservationKind::EquilibriumPressure, equilibrium_support::bubble_at(175.0));
    inconsistent.observed = 1.5 * *predict(HelmholtzModel(truth_case.parameters), inconsistent).value;
    zero.observations.push_back(inconsistent);
    zero.weights.assign(zero.observations.size(), 1.0);
    zero.weights.back() = 0.0;
    const auto ignored = fit(zero);
    CHECK(ignored.usable && ignored.residuals.size() == 6 && ignored.statuses.size() == 6);
    agree("r3/zero_weight/epsilon", "dispersion_energy_over_k", "K", ignored.physical[0], base.physical[0], 0.0, 1.0e-9, summary(ignored),
          "a zero-weighted inconsistent row does not move the fit", "r3");
    agree("r3/zero_weight/reported_residual", "residual", "1", ignored.residuals[5], std::log(1.0 / 1.5), 0.0, 1.0e-6, summary(ignored),
          "the zero-weighted row keeps its residual in the report", "r3");
    agree("r3/zero_weight/optimizer_residual", "weighted_residual", "1", ignored.weighted_residuals[5], 0.0, 0.0, 0.0, summary(ignored),
          "the zero-weighted row contributes nothing to the objective", "r3");
    double ignored_sum = 0.0;
    for (const double r : ignored.weighted_residuals) ignored_sum += r * r;
    agree("r3/zero_weight/variance_factor", "variance_factor", "1", *ignored.covariance.variance_factor, ignored_sum / (5.0 - 2.0), 0.0, 1.0e-12,
          summary(ignored), "the variance factor divides by the five weighted rows minus the rank, not six", "r3");
    // Rank deficiency: one observation cannot identify two coordinates.
    auto single = recovery_request(truth_case, 148.0, 0.0);
    single.observations = {truth_case.observations[0]};
    const auto deficient = fit(single);
    truth("r3/rank_deficient", "one row and two coordinates withhold covariance as rank deficient",
          deficient.usable && deficient.covariance.status == "rank_deficient" && deficient.covariance.rank == 1
              && deficient.covariance.optimizer.empty() && !deficient.covariance.variance_factor.has_value(), summary(deficient), "r3");
    CHECK(deficient.rank_records.size() == 1);
    const auto& rank = deficient.rank_records.front();
    CHECK(rank.rank == 1 && rank.dimension == 2);
    agree("r3/rank_cutoff", "absolute_singular_value_cutoff", "optimizer_J", rank.threshold,
        5.47456335865895e-15, 1e-29, 0.0, summary(deficient), "the existing singular spectrum determines the absolute cutoff", "r3");
    CHECK(deficient.first_numerical_failure && deficient.first_numerical_failure->stop == core::Stop::singular);
    CHECK(deficient.first_numerical_failure->rank->threshold == rank.threshold);
    CHECK(ignored.rank_records.front().rank == ignored.covariance.rank);
    // Covariance transformation through a tie and scales, and the variance factor.
    auto shared = recovery_request(truth_case, 148.0, 0.0);
    shared.parameters.sigma_constant = {3.6, 3.6};
    shared.coordinates.coordinates = {{{{ParameterFamily::SigmaConstant, {0}}, {ParameterFamily::SigmaConstant, {1}}}, 3.6, 0.1, 3.0, 4.0},
                                      {{{ParameterFamily::DispersionEnergyOverK, {0}}}, 150.0, 10.0, 130.0, 170.0}};
    const auto tied = fit(shared);  // converged: an unfinished iterate carries no covariance (r2-j2/not_converged)
    CHECK(tied.usable && tied.covariance.status == "available" && tied.covariance.physical.size() == 9);
    const auto& cu = tied.covariance.optimizer;
    const double hand[9] = {0.01 * cu[0], 0.01 * cu[0], 1.0 * cu[1], 0.01 * cu[0], 0.01 * cu[0], 1.0 * cu[1],
                            1.0 * cu[2], 1.0 * cu[2], 100.0 * cu[3]};
    for (std::size_t k = 0; k < 9; ++k)
        agree("r3/physical_covariance/" + std::to_string(k), "covariance", "native", tied.covariance.physical[k], hand[k], 0.0, 1.0e-12,
              summary(tied), "Cov_theta = D Cov_u D^T with D = S diag(scale)", "r3");
    // Cov_u is the inverse of J^T J of the independently qualified optimizer Jacobian.
    double worst = 0.0;
    for (std::size_t a = 0; a < 2; ++a)
        for (std::size_t b = 0; b < 2; ++b)
            worst = std::max(worst, std::abs(cu[a * 2] * gram_entry(tied, 0, b) + cu[a * 2 + 1] * gram_entry(tied, 1, b) - (a == b ? 1.0 : 0.0)));
    agree("r3/covariance_inverse", "max_abs(Cov_u J^T J - I)", "1", worst, 0.0, 1.0e-10, 0.0, summary(tied),
          "the optimizer covariance inverts the Gram matrix of the qualified Jacobian", "r3");
    double weighted = 0.0;
    for (const double r : tied.weighted_residuals) weighted += r * r;
    agree("r3/variance_factor", "variance_factor", "1", *tied.covariance.variance_factor, weighted / (5.0 - 2.0), 0.0, 1.0e-12,
          summary(tied), "the variance factor is sum(w r^2) / (n - rank)", "r3");
    // A declared bound the model itself refuses (m = 0) is a typed trial failure, not an exception: an
    // unattainable target (pressures fifty times the truth) drives the optimizer onto that bound.
    FitRequest inadmissible;
    inadmissible.controls = budgeted_controls();
    inadmissible.parameters = properties_support::neutral_binary();
    inadmissible.parameters.dispersion_energy_over_k[0] = 148.0;
    inadmissible.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 150.0, 10.0, 140.0, 170.0},
                                            {{{ParameterFamily::SegmentCount, {0}}}, 1.0, 2.0, 0.0, 5.0}};
    for (std::size_t i = 0; i < 3; ++i) {
        inadmissible.observations.push_back(truth_case.observations[i]);
        inadmissible.observations.back().observed *= 50.0;
    }
    const auto refused = fit(inadmissible);
    const bool typed = std::any_of(refused.failures.begin(), refused.failures.end(), [&](const TrialFailure& item) {
        return item.status == ObservationStatus::InadmissibleTrial && item.observation == 3
            && item.message.find("segment_count") != std::string::npos; });
    // The optimizer backs off the refused bound and converges there; the refusal is retained by name.
    truth("r3/inadmissible_trial", "a trial the model refuses by name is retained as inadmissible_trial and the fit converges at the bound",
          typed && refused.trial_failures > 0 && refused.status == "converged" && refused.usable
              && refused.covariance.status == "active_bound", summary(refused), "r3");
}

// Native rejected holds: a cached point followed by no root, and a Born-free singular relation.
void check_hydration_failures() {
    auto parameters = properties_support::aqueous_electrolyte();
    Observation energy;
    energy.kind = ObservationKind::SolvationGibbsEnergy; energy.temperature = 298.15; energy.pressure = 1e5;
    energy.composition = {1, 0, 0}; energy.ion_indices = {1, 2}; energy.stoichiometry = {1, 1};
    energy.species = {1}; energy.branch = DensityBranch::Liquid; energy.observed = -381930;
    FitRequest request; request.parameters = parameters;
    request.coordinates.coordinates = {{{{ParameterFamily::Kij, {1, 3}}}, 0, .2, -1, 1}};
    request.observations = {energy}; request.holds = {{energy, 1, 8, 1e-6, 1e-10, 50}};
    const auto start = request.coordinates.start(parameters);
    FitResult failures;
    ObservationEvaluator good(request, failures, {0});
    ObservationEvaluation values;
    CHECK(good.evaluate(start, true, values));
    const auto saved = values.chart.physical;
    CHECK(values.chart.holds.front().rank.rank == 1);
    auto exhausted = request; exhausted.controls = budgeted_controls(); exhausted.holds[0].maximum_iterations = 1;
    const auto limited = fit(exhausted);
    CHECK(limited.status == "initial_evaluation_failed" && limited.trial_failures == 1 && limited.jacobian_evaluations == 1);
    CHECK(limited.first_numerical_failure && limited.first_numerical_failure->stop == core::Stop::budget_exhausted);
    CHECK(limited.first_numerical_failure->iterations == 1 && std::isfinite(limited.first_numerical_failure->residual));
    CHECK(!limited.first_numerical_failure->criterion && limited.first_numerical_failure->rank->rank == 1);
    for (const bool singular : {false, true}) {
        auto bad = request;
        bad.holds[0].lower = saved.back() - 1e-6;
        bad.holds[0].upper = saved.back() + 1e-6;
        if (singular) {
            bad.parameters.electrolyte_form = ElectrolyteForm::DebyeHuckel; // manufactured H independent of d_Born
            HelmholtzModel constant(bad.parameters);
            bad.holds[0].reference.observed = *predict(constant, energy).value;
        }
        ObservationEvaluator evaluator(bad, failures, {0});
        auto rejected = start;
        if (!singular) { CHECK(evaluator.evaluate(start, true, values)); rejected[0] += 1.0; }
        const auto admitted = values.chart.physical;
        const bool accepted = evaluator.evaluate(rejected, true, values);
        const auto& cause = failures.failures.back().message;
        const auto& record = failures.failures.back().first_numerical_failure;
        CHECK(record && record->stop == (singular ? core::Stop::singular : core::Stop::left_domain));
        if (singular) {
            CHECK(record->rank && record->rank->rank == 0 && record->rank->dimension == 1);
            CHECK(record->rank->threshold == std::numeric_limits<double>::epsilon()
                * std::max(1.0, std::abs(bad.holds[0].reference.observed)) / std::max(1.0, std::abs(std::clamp(parameters.born_diameter_angstrom[1], bad.holds[0].lower, bad.holds[0].upper))));
        }
        truth(singular ? "hydration-failure/singular" : "hydration-failure/no-root",
            "A failed native held-coordinate trial publishes no diameter or Jacobian", !accepted && values.chart.physical.empty()
            && values.chart.map.empty() && values.gradient.empty() && cause.find(singular ? "singular" : "outside bracket") != std::string::npos,
            fields({field("target_J_per_mol", number(bad.holds[0].reference.observed)), field("cause", cause), field("optimizer", rejected)}), "issue-219");
        CHECK((singular ? good : evaluator).evaluate(start, true, values));
        truth(singular ? "hydration-failure/reuse-singular" : "hydration-failure/reuse-no-root",
            "A rejected trial preserves the admitted chart", values.chart.physical == admitted,
            "{}", "issue-219");
        CHECK(good.evaluate(start, true, values));
    }
}

// B3 selects the actual native electrolyte evaluator and complete nonunit tied chart.
void check_minimax_consistency() {
    // EqIDs: density_residual, pressure_density_derivative; #245 independent cold closures.
    const auto compare_state = [&](const HelmholtzModel& model, const State& seeded, const State& cold,
                                   double temperature, const std::vector<double>& x, const std::string& id) {
        const auto inputs = fields({field("T_K", temperature), field("x", x),
            field("seed_closure", core::describe(*seeded.density_convergence)), field("cold_closure", core::describe(*cold.density_convergence))});
        const auto compare = [&](const std::string& quantity, const char* unit, double value, double reference) {
            agree(id + "/" + quantity, quantity.c_str(), unit, value, reference,
          1e-12, 1e-10, inputs, "seeding agrees within solver accuracy", "issue-245");
        };
        compare("density", "mol/m3", seeded.molar_density, cold.molar_density);
        compare("slope", "Pa m3/mol", seeded.pressure_density_derivative, cold.pressure_density_derivative);
        for (std::size_t i = 0; i < cold.site_fractions.size(); ++i)
            compare("site/" + std::to_string(i), "1", seeded.site_fractions[i], cold.site_fractions[i]);
        const double h = cold.molar_density * 1e-5;
        const double slope = (model.state(temperature, cold.molar_density + h, x).pressure
            - model.state(temperature, cold.molar_density - h, x).pressure) / (2 * h);
        agree(id + "/slope-derivative", "P_rho", "Pa m3/mol", cold.pressure_density_derivative,
      slope, 1e-6 * std::max(1.0, std::abs(slope)), 0.0, inputs, "independent central pressure derivative", "issue-245");
    };
    auto ethane = properties_support::ethane();
    Observation density; density.kind = ObservationKind::Density;
    density.temperature = 200; density.pressure = 217800.69494692178;
    density.composition = {1}; density.branch = DensityBranch::Liquid;
    FitRequest density_request; density_request.parameters = ethane; density_request.observations = {density};
    density_request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 0, 1, 100, 250}};
    const auto u0 = density_request.coordinates.start(ethane);
    auto u1 = u0; u1[0] = 191.43;
    FitResult density_failures; ObservationEvaluator sequence(density_request, density_failures, {0});
    ObservationEvaluation first, density_seeded, density_promoted, density_cold;
    CHECK(sequence.evaluate(u0, true, first) && sequence.evaluate(u1, false, density_seeded)
        && sequence.evaluate(u1, true, density_promoted));
    FitResult fresh_failures; ObservationEvaluator density_fresh(density_request, fresh_failures, {0});
    CHECK(density_fresh.evaluate(u1, true, density_cold));
    agree("density-seed/column", "density_column", "native", density_promoted.gradient[0], density_cold.gradient[0],
        1e-12, 1e-10, "u0 to u1 then promotion", "independent cold density column", "issue-245");
    HelmholtzModel old_model(first.chart.parameters), target(density_promoted.chart.parameters);
    const auto previous = predict(old_model, density);
    const auto prediction = predict(target, density, std::nullopt, std::nullopt, nullptr, &previous, &old_model, false);
    compare_state(target, *prediction.state, *predict(target, density).state, 200, {1}, "density-seed/ethane");
    for (const std::string mode : {"missing-witness", "wrong-side", "irregular", "outside-window", "empty-domain"}) {
        auto offered = previous; auto row = density;
        if (mode == "missing-witness") offered.state->negative_slope_density.reset();
        if (mode == "wrong-side") row.branch = DensityBranch::Vapor;
        if (mode == "irregular") offered.state->density_convergence->amplification = 1e20;
        if (mode == "outside-window") offered.state->molar_density *= 1.02;
        if (mode == "empty-domain") row.pressure = 1e18;
        const auto candidate = predict(target, row, {}, {}, nullptr, &offered, &old_model,
            mode == "outside-window" || mode == "empty-domain"), reference = predict(target, row);
        truth("density-seed/" + mode, "an offered density seed returns the cold route's value or typed cause",
            candidate.status == reference.status && (candidate.value ? candidate.value == reference.value : candidate.message == reference.message),
            fields({field("offered_density", offered.state->molar_density), field("cold_cause", reference.message)}), "issue-245");
    }
    auto refused = u1; refused[0] = -1;
    CHECK(!sequence.evaluate(refused, true, density_seeded) && sequence.evaluate(u1, true, density_cold));
    CHECK(density_cold.predictions == density_promoted.predictions && density_cold.gradient == density_promoted.gradient
        && density_cold.chart.physical == density_promoted.chart.physical && density_cold.optimizer_gradient == density_promoted.optimizer_gradient);
    // The inherited same-liquid 130 -> 140 K candidate must still expose the new upper branch.
    auto lower_parameters = ethane; lower_parameters.dispersion_energy_over_k[0] = 130;
    auto upper_parameters = ethane; upper_parameters.dispersion_energy_over_k[0] = 140;
    HelmholtzModel lower(lower_parameters), upper(upper_parameters);
    const auto lower_result = predict(lower, density), upper_cold = predict(upper, density);
    CHECK(lower_result.stable_root_count == 1 && upper_cold.stable_root_count && *upper_cold.stable_root_count > 1);
    CHECK(predict(upper, density, {}, {}, nullptr, &lower_result, &lower, false).value == upper_cold.value);
    auto alcohol = properties_support::parameters(2);
    alcohol.segment_count = {2.3827, 2.3316}; alcohol.sigma_constant = {3.1771, 3.7086};
    alcohol.dispersion_energy_over_k = {198.24, 222.88}; alcohol.molar_mass_kg_per_mol = {.04606844, .0581222};
    alcohol.k_ij = {0, .028, .028, 0}; alcohol.site_component = {0, 0}; alcohol.site_multiplicity = {1, 1};
    alcohol.association_active = {0, 1, 1, 0}; alcohol.association_energy_over_k = {0, 2653.4, 2653.4, 0};
    alcohol.association_volume = {0, .032384, .032384, 0};
    HelmholtzModel alcohol_model(alcohol);
    const auto alcohol_cold = alcohol_model.state_at_pressure(45, 1e-4, {.5, .5}, DensityBranch::Vapor);
    compare_state(alcohol_model, alcohol_model.state_at_pressure(45, 1e-4, {.5, .5}, DensityBranch::Vapor,
        alcohol_cold.molar_density), alcohol_cold, 45, {.5, .5}, "density-seed/natural-275");
    auto parameters = properties_support::aqueous_electrolyte();
    Observation energy;
    energy.kind = ObservationKind::SolvationGibbsEnergy; energy.temperature = 298.15; energy.pressure = 1e5;
    energy.composition = {1, 0, 0}; energy.ion_indices = {1, 2}; energy.stoichiometry = {1, 1};
    energy.species = {1}; energy.branch = DensityBranch::Liquid; energy.observed = -381930;
    FitRequest request; request.parameters = parameters;
    request.controls = budgeted_controls();
    request.coordinates.coordinates = {{{{ParameterFamily::Kij, {1, 3}}}, 0, .2, -1, 1},
        {{{ParameterFamily::IonSuppressionCoefficient, {1}}, {ParameterFamily::IonSuppressionCoefficient, {2}}}, 0, .5, 0, 20}};
    request.holds = {{energy, 1, 8, 1e-6, 1e-10, 50}};
    auto gamma = energy; gamma.kind = ObservationKind::MeanIonicActivityCoefficient;
    gamma.species.clear(); gamma.molality = .1; gamma.form = ResidualForm::LogRatio; gamma.observed = 1;
    auto phi = gamma; phi.kind = ObservationKind::OsmoticCoefficient;
    auto epsilon = gamma; epsilon.kind = ObservationKind::BulkRelativePermittivity;
    epsilon.form = ResidualForm::Difference;
    epsilon.composition = formula_composition({1, 0, 0}, parameters, {1, 2}, {1, 1}, .1);
    request.observations = {gamma, phi, epsilon, energy};
    {
        HelmholtzModel columns(parameters, {{ParameterFamily::DispersionEnergyOverK, {0}}, {ParameterFamily::Kij, {1, 3}}});
        auto prediction = predict(columns, energy);
        const auto intact = *prediction.reference;
        auto& failed = prediction.reference->parameter_derivative[0][0];
        failed.action.reset(); failed.status = PropertyActionStatus::AssociationResidual;
        failed.diagnostics.association_residual = 2e-8;
        const auto check = [&](const PredictionGradient& gradient) {
            CHECK(gradient.status == ObservationStatus::ActionUnavailable && gradient.value.size() == 2);
            CHECK(!gradient.parameter_actions[0].action && gradient.parameter_actions[0].status == PropertyActionStatus::AssociationResidual);
            CHECK(gradient.parameter_actions[0].diagnostics.association_residual == 2e-8);
            CHECK(gradient.parameter_actions[1].action && std::isfinite(gradient.value[1]));
        };
        check(differentiate(columns, energy, prediction));
        auto transfer = energy; transfer.kind = ObservationKind::TransferGibbsEnergy;
        prediction.second_reference = prediction.reference; prediction.reference = intact;
        check(differentiate(columns, transfer, prediction));
        prediction.reference = prediction.second_reference;
        prediction.activity = electrolyte_activity(columns, gamma.temperature, gamma.pressure, {1, 0, 0}, {1, 2}, {1, 1}, .1, *prediction.reference);
        CHECK(std::isfinite(prediction.activity->mean_ionic_molality_scale));
        CHECK(prediction.reference->log_fugacity == intact.log_fugacity);
        CHECK(prediction.reference->parameter_derivative[0][1].action == intact.parameter_derivative[0][1].action);
        check(differentiate(columns, gamma, prediction));
        bool refused = false;
        try {
            HelmholtzModel nine(properties_support::nine_species());
            (void)infinite_dilution_reference(nine, 298.15, 1e5, {.7, 0, 0, 0, 0, 0, .1, .2, 0}, {3, 4}, {1, 1});
        } catch (const IonicReferenceBackground&) { refused = true; }
        CHECK(refused);
    }
    {
        FitRequest mixtures; mixtures.parameters = properties_support::nine_species();
        mixtures.parameters.charge_number[7] = 2; // HCl–MgCl2 inventory: H=3, Cl=4, Mg=7.
        mixtures.coordinates.coordinates = {{{{ParameterFamily::Kij, {3, 27}}}, 0, 1, -5, 5},
            {{{ParameterFamily::Kij, {31, 39}}}, 0, 1, -5, 5}};
        auto o = gamma; o.composition = {1,0,0,0,0,0,0,0,0}; o.ion_indices = {3,4};
        for (const auto [acid, background] : {std::pair{.2,.4}, {.7,1.1}, {.2,1.1}, {.2,0.}}) {
            o.molality = acid; o.finite_solute_molalities = std::vector<double>{0,0,0,acid,acid+2*background,0,0,background,0};
            mixtures.observations.push_back(o);
        }
        FitResult failures; ObservationEvaluation values;
        ObservationEvaluator grouped(mixtures, failures, {0,1,2,3});
        CHECK(grouped.evaluate(mixtures.coordinates.start(mixtures.parameters), true, values));
        HelmholtzModel mixed(values.chart.parameters, values.chart.active);
        const auto ref = infinite_dilution_reference(mixed, o.temperature, o.pressure, o.composition, o.ion_indices, o.stoichiometry);
        for (std::size_t i = 0; i < mixtures.observations.size(); ++i) {
            const auto& row = mixtures.observations[i];
            const auto direct = electrolyte_activity(mixed, row.temperature, row.pressure, row.composition, row.ion_indices,
                row.stoichiometry, row.molality, ref, {}, row.finite_solute_molalities);
            CHECK(values.statuses[i] == ObservationStatus::Available);
            agree("HCl-MgCl2/"+std::to_string(i), "log_activity", "1", std::log(values.predictions[i]),
                direct.mean_ionic_log_molality_scale, 1e-10, 0, "298.15 K,100000 Pa", "Properties mapping", "ALG-FIT");
            for (std::size_t c = 0; c < 2; ++c)
                agree("HCl-MgCl2/"+std::to_string(i)+"/"+std::to_string(c), "parameter_column", "1",
                    values.gradient[2*i+c], *direct.mean_ionic_parameter_derivative[c].action, 1e-12, 1e-10,
                    "H-water,H-Cl", "Properties mapping", "ALG-FIT");
        }
        auto binary = o; binary.finite_solute_molalities.reset();
        CHECK(std::abs(std::log(*predict(mixed, binary).value / values.predictions.back())) < 1e-10);
        (*o.finite_solute_molalities)[4] += .1;
        bool refused = false; try { validate(o, mixed); } catch (const std::invalid_argument&) { refused = true; }
        CHECK(refused);
        HelmholtzModel dilute_model(mixtures.parameters);
        const auto dilute_ref = infinite_dilution_reference(dilute_model, o.temperature, o.pressure, o.composition, o.ion_indices, o.stoichiometry);
        for (double acid : {1e-15, 1.0000005e-8/2}) {
            const auto dilute = electrolyte_activity(dilute_model, o.temperature, o.pressure, o.composition, o.ion_indices,
                o.stoichiometry, acid, dilute_ref, {}, std::vector<double>{0,0,0,acid,acid,0,0,0,0});
            CHECK((*dilute.state.log_fugacity_coefficient)[3].has_value());
            CHECK(std::holds_alternative<double>(*dilute.osmotic_coefficient) == (acid > 1e-9));
        }
    }
    const auto u = request.coordinates.start(parameters);
    request.observations = observations_at(request, u, {0, 0, .1, 0});
    request.weights = {1, 1, 0, 0};
    request.checks = {{"decades", {0, 1}, CheckMetric::AbsoluteLog10Ratio, .02, {}, false, false},
        {"permittivity", {2}, CheckMetric::AbsoluteDifference, .2, .2, false, false},
        {"energy", {3}, CheckMetric::AbsoluteDifference, 10, {}, false, false}};
    request.minimax = MinimaxRequest{{0, 1}, {.02, 1, 1000}, .02, 1e-10, 1e-6, .1};
    FitResult failures; ObservationEvaluator evaluator(request, failures, {0, 1, 2, 3}); ObservationEvaluation a;
    CHECK(evaluator.evaluate(u, true, a));
    auto shifted = u; shifted[0] += .001;
    ObservationEvaluation seeded, promoted, cold;
    FitResult independent_failures;
    ObservationEvaluator fresh(request, independent_failures, {0, 1, 2, 3});
    CHECK(evaluator.evaluate(shifted, false, seeded) && evaluator.evaluate(shifted, true, promoted) && fresh.evaluate(shifted, true, cold));
    // A column entry that cancels an observation's sensitivities (the held energy's optimizer column cancels
    // terms near 1.5e4) carries their rounding: its absolute floor is four ulps of the row's largest physical
    // sensitivity (#478); values keep the 1e-12 floor.
    const std::size_t sensitivities = cold.gradient.size() / cold.predictions.size();
    const auto floor = [&](std::size_t row) {
        double scale = 0.0;
        for (std::size_t k = 0; k < sensitivities; ++k) scale = std::max(scale, std::abs(cold.gradient[row * sensitivities + k]));
        return std::max(1e-12, 4.0 * std::numeric_limits<double>::epsilon() * scale);
    };
    for (const auto& pair : {std::pair{promoted.predictions, cold.predictions}, std::pair{promoted.gradient, cold.gradient}, std::pair{promoted.optimizer_gradient, cold.optimizer_gradient}})
        for (std::size_t i = 0, width = pair.first.size() / cold.predictions.size(); i < pair.first.size(); ++i)
            agree("density-seed/electrolyte/" + std::to_string(i), "value_or_column", "native", pair.first[i], pair.second[i],
                width == 1 ? 1e-12 : floor(i / width), 1e-10, "held aqueous changed point and promotion", "fresh-point identity", "issue-245");
    HelmholtzModel source_model(a.chart.parameters), shifted_model(promoted.chart.parameters);
    const auto prior = predict(source_model, gamma);
    const auto warm = predict(shifted_model, gamma, {}, {}, nullptr, &prior, &source_model, false);
    const auto independent = predict(shifted_model, gamma);
    compare_state(shifted_model, *warm.reference->pure_solvent_state, *independent.reference->pure_solvent_state,
        gamma.temperature, gamma.composition, "density-seed/pure-solvent");
    compare_state(shifted_model, warm.activity->state, independent.activity->state, gamma.temperature,
        formula_composition(gamma.composition, promoted.chart.parameters, gamma.ion_indices, gamma.stoichiometry, gamma.molality), "density-seed/activity");
    const auto b = fit_minimax(request, u);
    CHECK(b.minimax.has_value());
    std::size_t quantity = 0;
    for (const auto& pair : {std::pair{b.predictions, a.predictions},
            std::pair{b.physical_jacobian, a.gradient}, std::pair{b.optimizer_jacobian, a.optimizer_gradient}}) {
        CHECK(pair.first.size() == pair.second.size());
        const auto columns = quantity == 0 ? 1 : quantity == 1 ? a.chart.physical.size() : u.size();
        for (std::size_t i = 0; i < pair.first.size(); ++i) {
            const bool energy_row = request.observations[i / columns].kind == ObservationKind::SolvationGibbsEnergy;
            const bool diameter_column = quantity == 1 && i % columns == columns - 1;
            const char* unit = energy_row ? (diameter_column ? "J/(mol Angstrom)" : "J/mol")
                : (diameter_column ? "1/Angstrom" : "1");
            const double tolerance = pair.second[i] == 0.0 ? 1e-12 : 1e-12 * std::abs(pair.second[i]);
            agree("minimax-consistency/" + std::to_string(quantity) + "/" + std::to_string(i),
          quantity == 0 ? "prediction" : quantity == 1 ? "physical_column" : "optimizer_column",
          unit, pair.first[i], pair.second[i], tolerance, 0.0,
          "T=298.15 K,P=100000 Pa,m=.1 mol/kg", "B preserves A at the same held point", "r2-minimax");
        }
        ++quantity;
    }
    const auto checks = evaluate_checks(request, b.parameters);
    for (std::size_t g = 0; g < checks.size(); ++g) {
        agree("minimax-consistency/check/" + std::to_string(g), "maximum", checks[g].unit.c_str(),
            b.checks[g].maximum, checks[g].maximum, checks[g].maximum == 0 ? 1e-12 : 1e-12 * std::abs(checks[g].maximum), 0.0, "same stored diameters",
            "B final checks use the A check owner", "r2-minimax");
        agree("minimax-consistency/threshold/" + std::to_string(g), "maximum_threshold", checks[g].unit.c_str(),
            *b.minimax->frozen_checks[g].maximum_threshold, *checks[g].maximum_threshold,
            *checks[g].maximum_threshold == 0 ? 1e-12 : 1e-12 * std::abs(*checks[g].maximum_threshold), 0.0,
            "same explicit thresholds", "B freezes declared unweighted checks", "r2-minimax");
    }
}

void check_inventory_boundary(const HelmholtzModel& model, const EquilibriumProblem& problem,
                              const EquilibriumResult& result, const std::string& id, double material_limit, bool rank_required = true) {
    if (!result.success) throw std::runtime_error(result.message);
    CHECK(result.requested_tolerance_met.value_or(false));
    const auto compiled = equilibrium::compile_problem(model, problem);
    const auto point = compiled.prepare(result.coordinates);
    const auto values = compiled.jacobian(result.coordinates);
    const Eigen::Index size = static_cast<Eigen::Index>(result.coordinates.size());
    CHECK(values.size() == static_cast<std::size_t>(size * size));
    const Eigen::MatrixXd jacobian = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(values.data(), size, size);
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian);
    const double condition = svd.singularValues()(0) / svd.singularValues()(size - 1);
    double raw = 0, classified = 0, normalization = 0, charge = 0, material = 0;
    std::vector<double> amounts(compiled.n_feed.size(), 0), slopes;
    bool exact = true;
    for (const double value : result.residuals) raw = std::max(raw, std::abs(value));
    for (const double value : result.classified_residuals) classified = std::max(classified, std::abs(value));
    for (std::size_t p = 0; p < point.composition.size(); ++p) {
        normalization = std::max(normalization, std::abs(std::accumulate(point.composition[p].begin(), point.composition[p].end(), 0.0) - 1));
        double phase_charge = 0;
        for (std::size_t i = 0; i < amounts.size(); ++i) {
            amounts[i] += result.amounts[p] * point.composition[p][i];
            phase_charge += result.amounts[p] * point.composition[p][i] * model.parameters().charge_number[i];
        }
        charge = std::max(charge, std::abs(phase_charge));
        for (const auto i : compiled.removed_species) exact = exact && point.composition[p][i] == 0 && result.amounts[p] * point.composition[p][i] == 0;
        slopes.push_back(model.state(point.temperature, result.molar_densities[p], point.composition[p]).pressure_density_derivative);
    }
    for (const auto& row : equilibrium::balance_basis(amounts.size(), problem.reactions)) {
        double error = 0;
        for (std::size_t i = 0; i < amounts.size(); ++i) error += row[i] * (amounts[i] - compiled.n_feed[i]);
        material = std::max(material, std::abs(error));
    }
    const auto inputs = fields({field("T_K", point.temperature), field("P_Pa", point.pressure), field("phase_amounts_mol", result.amounts),
        field("fractions", result.mole_fractions), field("densities_mol_m3", result.molar_densities), field("slopes_Pa_m3_mol", slopes),
        field("endpoint_rank", jacobian.fullPivLu().rank()), field("endpoint_size", size), field("condition_2", condition), field("balance_rank", compiled.rank)});
    agree(id + "/raw", "maximum_raw_row", "1", raw, 0, 1e-10, 0, inputs, "raw stationarity is resolved", "issue-269");
    agree(id + "/classified", "maximum_classified_row", "1", classified, 0, 1e-10, 0, inputs, "classified rows meet the request", "issue-269");
    agree(id + "/material", "original_balance_error", "mol", material, 0, material_limit, 0, inputs, "original conserved inventories are preserved", "issue-269");
    agree(id + "/charge", "phase_charge_error", "charge-equivalent mol", charge, 0, material_limit, 0, inputs, "phase charge is conserved", "issue-269");
    agree(id + "/normalization", "normalization_error", "1", normalization, 0, 1e-12, 0, inputs, "phase fractions sum to one", "issue-269");
    if (rank_required) truth(id + "/endpoint-rank", "the reduced Jacobian at the solved endpoint is full rank",
  jacobian.fullPivLu().rank() == size && std::isfinite(condition), inputs, "issue-269");
    truth(id + "/regular", "distinct mechanical branches have finite positive slopes and exact absent zeros",
  std::isfinite(point.pressure) && point.pressure > 0 && exact
  && std::isfinite(slopes[0]) && std::isfinite(slopes[1]) && slopes[0] > 0 && slopes[1] > 0
  && result.molar_densities[0] > result.molar_densities[1], inputs, "issue-269");
}

void check_propane_inventory_boundary() {
    auto parameters = properties_support::parameters(3);
    parameters.segment_count.assign(3, 2.0020); parameters.sigma_constant.assign(3, 3.6184);
    parameters.dispersion_energy_over_k.assign(3, 208.11);
    parameters.molar_mass_kg_per_mol.assign(3, 0.044096);
    const HelmholtzModel model(parameters);
    auto problem = pure_saturation(300, 1e6);
    problem.feed.values = {1, 0, 0}; problem.reactions = {{"D=E", {0, -1, 1}, 1.0}};
    for (auto& phase : problem.phases) { phase.support = {0, 1, 2}; phase.support_declared = true; phase.composition_guess = {.8, .1, .1}; }
    for (const bool dew : {false, true}) {
        problem.phases[0].amount = dew ? AmountSpec::pinned(0) : AmountSpec::free();
        problem.phases[1].amount = dew ? AmountSpec::free() : AmountSpec::pinned(0);
        const std::string id = dew ? "incipient/propane-dew" : "incipient/propane-bubble";
        const auto result = equilibrium::solve_equilibrium(model, problem);
        check_inventory_boundary(model, problem, result, id, 1e-10);
        truth(id + "/absent-pair", "both manufactured carriers have literal zero composition in both phases",
            result.mole_fractions[1] == 0 && result.mole_fractions[2] == 0 && result.mole_fractions[4] == 0 && result.mole_fractions[5] == 0,
            fields({field("fractions", result.mole_fractions)}), "issue-269");
        agree(id + "/pressure", "boundary_pressure", "Pa", *result.pressure, 998660.8957398997, .01, 0, "teqp same-model reference", "independent propane coexistence pressure", "issue-269");
        for (std::size_t p = 0; p < 2; ++p)
            agree(id + "/density/" + std::to_string(p), "molar_density", "mol/m^3", result.molar_densities[p], p == 0 ? 11100.25123320962 : 482.5121268150332, 0, 1e-7, "teqp same-model reference", "independent mechanical branches", "issue-269");
        auto control = problem;
        for (auto& phase : control.phases) { phase.support = {0}; phase.composition_guess = {1}; }
        const auto reduced = equilibrium::solve_equilibrium(model, control);
        equilibrium::EquilibriumStart start; start.compositions = {{.7, .2, .1}, {.6, .1, .3}};
        const auto started = equilibrium::solve_equilibrium_starts(model, problem, {start}, {});
        CHECK(started.selected && started.results[*started.selected]);
        for (const auto& comparison : {std::pair{"explicit", &reduced}, std::pair{"declared-start", &*started.results[*started.selected]}}) {
            const auto* other = comparison.second;
            const std::string prefix = id + "/" + comparison.first;
            check_inventory_boundary(model, problem, *other, prefix, 1e-10);
            reference_scale_agrees(prefix + "/pressure", "boundary_pressure", "Pa", *other->pressure, *result.pressure, 1e-4, 1e-8, "same wheel", "full-width guesses and explicit face agree", "issue-269");
            for (std::size_t i = 0; i < result.mole_fractions.size(); ++i)
                reference_scale_agrees(prefix + "/fraction/" + std::to_string(i), "mole_fraction", "1", other->mole_fractions[i], result.mole_fractions[i], 1e-12, 1e-8, "same wheel", "supported fractions and literal zeros agree", "issue-269");
            for (std::size_t p = 0; p < 2; ++p)
                reference_scale_agrees(prefix + "/amount/" + std::to_string(p), "phase_amount", "mol", other->amounts[p], result.amounts[p], 1e-12, 1e-8, "same wheel", "phase amounts agree", "issue-269");
        }
    }
    const auto refuse = [&](const HelmholtzModel& m, const EquilibriumProblem& request, const std::string& id, const std::string& reason) {
        const auto cause = rejection([&] { equilibrium::compile_problem(m, request); });
        truth("incipient/refusal/" + id, "excluded incipient domain returns no fabricated boundary", cause.find("zero_inventory_incipient_support") != std::string::npos && cause.find(reason) != std::string::npos,
            fields({field("cause", cause)}), "issue-269");
    };
    auto excluded = problem; excluded.phases[1].kind = PhaseKind::Liquid;
    refuse(model, excluded, "same-kind", "opposite liquid/vapor");
    excluded = problem; excluded.temperature = FreeValue::variable(300); excluded.pressure = FreeValue::fixed(1e6);
    refuse(model, excluded, "free-temperature", "fixed positive T and free P");
    const HelmholtzModel ionic(properties_support::nine_species());
    // The incipient phase is the liquid: a vapor support holds no ion (ALG-EQUILIBRIUM).
    excluded = equilibrium_support::boundary(); excluded.feed.values = {1, 0, 0, .01, .01, 0, 0, 0, 0};
    excluded.phases[0].amount = AmountSpec::pinned(0); excluded.phases[1].amount = AmountSpec::free();
    refuse(ionic, excluded, "charged-incipient", "uncharged species");
    excluded.phases[0].support = {0, 3}; excluded.phases[0].support_declared = true; excluded.phases[1].kind = PhaseKind::Liquid;
    const auto cause = rejection([&] { equilibrium::compile_problem(ionic, excluded); });
    truth("incipient/refusal/charge-infeasible", "one charge sign cannot form an electroneutral phase", cause.find("support z sign") != std::string::npos,
        fields({field("cause", cause)}), "issue-269");
}

void check_inventory_face() {
    truth("inventory-face/exact-dot", "error-free conservation arithmetic distinguishes cancellation from a lost positive unit",
        !equilibrium::exact_zero_dot({1e16, 1, -1e16}, {1, 1, 1})
        && equilibrium::exact_zero_dot({1e16, 1, -1e16, -1}, {1, 1, 1, 1}),
        fields({field("nonzero_left", std::vector<double>{1e16, 1, -1e16}), field("nonzero_right", std::vector<double>{1, 1, 1}),
            field("zero_left", std::vector<double>{1e16, 1, -1e16, -1}), field("zero_right", std::vector<double>{1, 1, 1, 1})}), "issue-231");
    HelmholtzModel model(properties_support::neutral_ternary());
    EquilibriumProblem problem;
    problem.temperature = FreeValue::fixed(300);
    problem.pressure = FreeValue::fixed(1e5);
    problem.feed.kind = equilibrium::FeedKind::Amounts;
    problem.feed.values = {1, 0, 0};
    problem.phases = {equilibrium_support::phase("L", PhaseKind::Liquid)};
    problem.reactions = {{"A+B=C", {-1, -1, 1}, 1.0}};
    const auto declarations = [&] {
        return fields({field("T_K", problem.temperature.value), field("P_Pa", problem.pressure.value),
            field("feed_mol", problem.feed.values), field("stoichiometry", problem.reactions[0].stoichiometry),
            field("L_support", problem.phases[0].support_declared ? "A" : "A,B,C")});
    };
    const auto compiled = equilibrium::compile_problem(model, problem);
    const std::vector<double> formal{0, 1, -1};
    CHECK(compiled.B.size() == 1);
    agree("inventory-face/formal-feed", "d_conserved_A_inventory", "mol",
        std::inner_product(compiled.B[0].begin(), compiled.B[0].end(), formal.begin(), 0.0), -1, 1e-12, 0,
        declarations() + ";" + fields({field("delta_feed_mol", formal), field("interpretation", "algebraic balance derivative only")}),
        "the full reduced feed map retains the A+C inventory derivative", "issue-231");
    truth("inventory-face/vanished-feed", "the same signed feed direction preserves the vanished B+C inventory",
        std::all_of(compiled.inventory_restrictions.begin(), compiled.inventory_restrictions.end(),
            [&](const auto& row) { return equilibrium::exact_zero_dot(row, formal); }),
        declarations() + ";" + field("delta_feed_mol", formal), "issue-231");
    const auto value = equilibrium::solve_equilibrium(model, problem);
    CHECK(value.success);
    agree("inventory-face/absent-B", "n_B", "mol", value.amounts[0] * value.mole_fractions[1], 0, 0, 0,
        declarations(), "nonnegative B+C conservation forces both carriers absent", "issue-231");
    problem.phases[0].support = {0};
    problem.phases[0].support_declared = true;
    CHECK(equilibrium::compile_problem(model, problem).feasible);
    auto uncovered = problem;
    uncovered.reactions[0].stoichiometry = {-1, 1, 1};
    truth("inventory-face/uncovered", "support-conditional mixed-sign certificates cannot erase attainable products",
        rejection([&] { equilibrium::compile_problem(model, uncovered); }).find("has no support") != std::string::npos,
        fields({field("stoichiometry", uncovered.reactions[0].stoichiometry), field("feed_mol", uncovered.feed.values),
            field("L_support", "A")}), "issue-231");
    problem.phases[0].support.clear();
    problem.phases[0].support_declared = false;
    problem.reactions[0].stoichiometry = {0, -1, 1};
    for (double inventory : {0.0, 1e-10, 1e-12, 1e-14}) {
        problem.feed.values[1] = inventory;
        const auto at = equilibrium::compile_problem(model, problem);
        truth("inventory-face/trace/" + number(inventory), "positive conserved inventories remain positive-species charts",
            inventory == 0 ? at.removed_species == std::vector<std::size_t>{1, 2} : at.removed_species.empty(),
            declarations() + ";" + field("B_plus_C_mol", inventory), "issue-231");
    }
    auto padded = properties_support::parameters(5);
    const auto ternary = equilibrium_support::synthetic_reactive_ternary();
    const auto ethane = properties_support::ethane();
    padded.segment_count = ternary.segment_count; padded.segment_count.resize(5, ethane.segment_count[0]);
    padded.sigma_constant = ternary.sigma_constant; padded.sigma_constant.resize(5, ethane.sigma_constant[0]);
    padded.dispersion_energy_over_k = ternary.dispersion_energy_over_k;
    padded.dispersion_energy_over_k.resize(5, ethane.dispersion_energy_over_k[0]);
    padded.molar_mass_kg_per_mol = ternary.molar_mass_kg_per_mol;
    padded.molar_mass_kg_per_mol.resize(5, ethane.molar_mass_kg_per_mol[0]);
    for (std::size_t i = 0; i < 3; ++i) for (std::size_t j = 0; j < 3; ++j)
        padded.k_ij[i * 5 + j] = ternary.k_ij[i * 3 + j];
    padded.site_component = ternary.site_component; padded.site_multiplicity = ternary.site_multiplicity;
    padded.association_active = ternary.association_active;
    padded.association_energy_over_k = ternary.association_energy_over_k;
    padded.association_volume = ternary.association_volume;
    HelmholtzModel padded_model(padded);
    auto full_guess = equilibrium_support::synthetic_reactive_problem();
    full_guess.feed.values.resize(5, 0); full_guess.standard_gibbs_j_per_mol.resize(5, 0);
    full_guess.reactions[0].stoichiometry.resize(5, 0);
    full_guess.reactions.push_back({"D=E", {0, 0, 0, -1, 1}, 1.0});
    auto reduced_guess = full_guess;
    for (std::size_t q = 0; q < 2; ++q) {
        full_guess.phases[q].composition_guess.resize(5, 0);
        reduced_guess.phases[q].support = {0, 1, 2}; reduced_guess.phases[q].support_declared = true;
    }
    const auto full_value = equilibrium::solve_equilibrium(padded_model, full_guess);
    const auto reduced_value = equilibrium::solve_equilibrium(padded_model, reduced_guess);
    CHECK(full_value.success && reduced_value.success);
    double scaled_state_error = 0;
    for (std::size_t i = 0; i < full_value.mole_fractions.size(); ++i)
        scaled_state_error = std::max(scaled_state_error,
            std::abs(full_value.mole_fractions[i] - reduced_value.mole_fractions[i])
                / (1e-12 + 1e-8 * std::max(std::abs(full_value.mole_fractions[i]), std::abs(reduced_value.mole_fractions[i]))));
    for (std::size_t q = 0; q < full_value.amounts.size(); ++q)
        scaled_state_error = std::max(scaled_state_error,
            std::abs(full_value.amounts[q] - reduced_value.amounts[q])
                / (1e-12 + 1e-8 * std::max(std::abs(full_value.amounts[q]), std::abs(reduced_value.amounts[q]))));
    agree("inventory-face/supplied-guess", "maximum_scaled_physical_state_error", "1", scaled_state_error, 0, 1, 0,
        fields({field("fixture", "synthetic_reactive_ternary plus two ethane records"),
            field("T_K", full_guess.temperature.value), field("P_Pa", full_guess.pressure.value), field("feed_mol", full_guess.feed.values),
            field("A_plus_B_to_C", full_guess.reactions[0].stoichiometry), field("D_to_E", full_guess.reactions[1].stoichiometry),
            field("K_each", 1), field("full_support", "0,1,2,3,4"), field("reduced_support", "0,1,2"),
            field("full_L1_guess", full_guess.phases[0].composition_guess), field("full_L2_guess", full_guess.phases[1].composition_guess),
            field("reduced_L1_guess", reduced_guess.phases[0].composition_guess), field("reduced_L2_guess", reduced_guess.phases[1].composition_guess),
            field("full_phase_fractions", full_value.mole_fractions), field("reduced_phase_fractions", reduced_value.mole_fractions),
            field("full_amounts_mol", full_value.amounts), field("reduced_amounts_mol", reduced_value.amounts),
            field("full_classified_residual", full_value.convergence.residual), field("reduced_classified_residual", reduced_value.convergence.residual),
            field("full_iterations", full_value.convergence.iterations), field("reduced_iterations", reduced_value.convergence.iterations)}),
        "restricting a supplied full-support guess preserves explicit-face phase fractions and amounts", "issue-231");
    problem.feed.values[1] = 0;
    problem.pressure = FreeValue::variable(1e5);
    auto incipient = equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0));
    incipient.support = {2}; incipient.support_declared = true;
    problem.phases.push_back(incipient);
    const auto cause = rejection([&] { equilibrium::compile_problem(model, problem); });
    truth("inventory-face/incipient", "certified absence leaves the incipient effective support empty and is refused before solving",
        cause.find("zero_inventory_incipient_support") != std::string::npos && cause.find("component 2") != std::string::npos
        && cause.find("phase V") != std::string::npos && cause.find("empty effective support") != std::string::npos,
        declarations() + ";" + fields({field("V_support", "C"), field("V_amount_fraction", 0.0), field("cause", cause)}), "issue-231");
}

void check_neutral_electrostatics() {
    auto p = properties_support::aqueous_electrolyte();
    p.born_c_shell = p.born_c_dielectric = 1.0;
    const std::vector<ActiveParameter> active{{ParameterFamily::BornDiameter, {1}}};
    HelmholtzModel model(p, active);
    const auto water_inputs = fields({field("T_K", 298.15), field("P_Pa", 1e5), field("composition", std::vector<double>{1, 0, 0}),
        field("components", "water,Na+,Cl-"), field("Born_c_shell", 1.0), field("Born_c_dielectric", 1.0)});
    const auto neutral = model.state_at_pressure(298.15, 1e5, {1, 0, 0}, DensityBranch::Liquid);
    auto selected = properties_support::parameters(1);
    selected.segment_count[0] = p.segment_count[0]; selected.sigma_constant[0] = p.sigma_constant[0];
    selected.dispersion_energy_over_k[0] = p.dispersion_energy_over_k[0]; selected.molar_mass_kg_per_mol[0] = p.molar_mass_kg_per_mol[0];
    selected.sigma_term_component = p.sigma_term_component; selected.sigma_term_amplitude = p.sigma_term_amplitude;
    selected.sigma_term_exponent = p.sigma_term_exponent; selected.site_component = p.site_component;
    selected.site_multiplicity = p.site_multiplicity; selected.association_active = p.association_active;
    selected.association_energy_over_k = p.association_energy_over_k; selected.association_volume = p.association_volume;
    const auto control = HelmholtzModel(selected).state_at_pressure(298.15, 1e5, {1}, DensityBranch::Liquid);
    agree("neutral-electrostatics/density", "rho", "mol/m3", neutral.molar_density, control.molar_density,
        0, 1e-9, water_inputs, "full ionic declaration agrees with neutral selection", "issue-231");
    reference_scale_agrees("neutral-electrostatics/enthalpy", "h_res", "J/mol", neutral.residual_enthalpy,
        control.residual_enthalpy, 1e-6, 1e-9, water_inputs, "neutral support retains the same residual calorics", "issue-231");
    auto boundary = equilibrium_support::flash(3);
    boundary.temperature = FreeValue::fixed(298.15); boundary.pressure = FreeValue::variable(1e5);
    boundary.feed.values = {1, 0, 0};
    boundary.phases = {equilibrium_support::phase("L", PhaseKind::Liquid),
        equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0))};
    for (auto& phase : boundary.phases) { phase.support = {0}; phase.support_declared = true; }
    auto selected_boundary = boundary; selected_boundary.feed.values = {1};
    const HelmholtzModel water_model(selected);
    const auto full_cell = equilibrium::compile_problem(model, boundary);
    const auto selected_cell = equilibrium::compile_problem(water_model, selected_boundary);
    const auto seed = full_cell.prepare(equilibrium::boundary_start(full_cell));
    const auto selected_seed = selected_cell.prepare(equilibrium::boundary_start(selected_cell));
    const auto seed_inputs = fields({field("T_K", seed.temperature), field("P_declared_seed_Pa", 1e5),
        field("P_returned_full_Pa", seed.pressure), field("P_returned_selected_Pa", selected_seed.pressure),
        field("feed_mol", boundary.feed.values), field("L_support", "water"), field("V_support", "water"),
        field("V_fraction", 0), field("composition_guesses", "empty"), field("fixture", "aqueous_electrolyte versus selected water")});
    for (const auto pressure : {std::pair{"full", seed.pressure}, std::pair{"selected", selected_seed.pressure}})
        recording_support::observe((std::string("neutral-electrostatics/cold-seed/pressure/") + pressure.first).c_str(),
            "P_returned_seed", "Pa", pressure.second, std::nullopt, seed_inputs, "Observation", "boundary_start",
            "existing full-declaration and pure-selected starts retain their pressure policies", "issue-231");
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto id = "neutral-electrostatics/cold-seed/" + std::to_string(phase);
        const double density = 1 / seed.volume[phase];
        const auto selected_state = water_model.state_at_pressure(seed.temperature, seed.pressure, {1},
            phase == 0 ? DensityBranch::Liquid : DensityBranch::Vapor);
        agree(id + "/density", "rho_seed", "mol/m3", density, selected_state.molar_density,
            0, 1e-9, seed_inputs, "decoded full cold density agrees with selected water at the same T/P/branch", "issue-231");
        reference_scale_agrees(id + "/enthalpy", "h_res_seed", "J/mol",
            model.state(seed.temperature, density, seed.composition[phase]).residual_enthalpy,
            selected_state.residual_enthalpy,
            1e-6, 1e-9, seed_inputs, "cold seeds retain the same supported residual calorics", "issue-231");
        agree(id + "/amount", "n_phase_seed", "mol", seed.amount[phase], selected_seed.amount[phase],
            0, 0, seed_inputs, "cold seed phase amounts retain their declared basis", "issue-231");
        agree(id + "/water", "x_water_seed", "1", seed.composition[phase][0], selected_seed.composition[phase][0],
            0, 0, seed_inputs, "supported cold composition agrees with selected water", "issue-231");
        for (std::size_t species = 1; species < 3; ++species)
            agree(id + "/absent/" + std::to_string(species), "x_absent_seed", "1", seed.composition[phase][species],
                0, 0, 0, seed_inputs, "globally removed cold seed entries remain literal zero", "issue-231");
    }
    for (const bool temperature : {true, false}) {
        const auto full_actions = model.pressure_actions(298.15, 1e5, {1, 0, 0}, DensityBranch::Liquid,
            {PropertyObservable::ResidualEnthalpy}, {{temperature ? 1.0 : 0.0, temperature ? 0.0 : 1.0, {0, 0, 0}, {}}});
        const auto neutral_actions = HelmholtzModel(selected).pressure_actions(298.15, 1e5, {1}, DensityBranch::Liquid,
            {PropertyObservable::ResidualEnthalpy}, {{temperature ? 1.0 : 0.0, temperature ? 0.0 : 1.0, {0}, {}}});
        CHECK(full_actions.density_status == PropertyActionStatus::Available && neutral_actions.density_status == PropertyActionStatus::Available);
        CHECK(full_actions.results[0].action && neutral_actions.results[0].action);
        const auto suffix = temperature ? "/T" : "/P";
        agree(std::string("neutral-electrostatics/density-action") + suffix, temperature ? "d_rho_dT_at_P" : "d_rho_dP_at_T",
            temperature ? "mol/(m^3 K)" : "mol/(m^3 Pa)",
            full_actions.density_action[0], neutral_actions.density_action[0], 0, FIRST_ORDER, water_inputs,
            "supported pressure-closed density actions agree with neutral selection", "issue-231");
        agree(std::string("neutral-electrostatics/enthalpy-action") + suffix, temperature ? "d_h_res_dT_at_P" : "d_h_res_dP_at_T",
            temperature ? "J/(mol K)" : "J/(mol Pa)",
            *full_actions.results[0].action, *neutral_actions.results[0].action, 0, FIRST_ORDER, water_inputs,
            "supported caloric actions agree with neutral selection", "issue-231");
    }
    const auto no_solvent = model.pressure_actions(298.15, 1e5, {0, .5, .5}, DensityBranch::Liquid, {}, {{0, 0, {0, 0, 0}, {}}});
    truth("neutral-electrostatics/ionic-no-solvent", "positive-ion solvent-free states preserve the existing density failure",
        no_solvent.density_status != PropertyActionStatus::Available && !no_solvent.density_message.empty(),
        fields({field("T_K", 298.15), field("P_Pa", 1e5), field("composition", std::vector<double>{0, .5, .5}),
            field("cause", no_solvent.density_message)}), "issue-231");
    truth("neutral-electrostatics/charged-slots", "absent charged potentials, fugacities and parameter columns are unavailable",
        !neutral.residual_chemical_potential_over_rt[1] && !(*neutral.log_fugacity_coefficient)[1]
        && !neutral.log_fugacity_parameter_derivative[1] && neutral.component_property_status[1] == PropertyActionStatus::IonicInsertion
        && neutral.residual_chemical_potential_over_rt[0] && (*neutral.log_fugacity_coefficient)[0], water_inputs, "issue-231");
    std::array<double, 2> previous{1, 1};
    for (double salt : {0.0, 1e-4, 1e-6, 1e-8}) {
        const std::vector<double> x{1 - 2 * salt, salt, salt};
        const std::vector<double> point{298.15, neutral.molar_density, p.born_diameter_angstrom[1]};
        std::vector<detail::AD> variables(point.begin(), point.end());
        CppAD::Independent(variables);
        const auto record = detail::make_tape_parameters(p, active, {variables[2]});
        const auto terms = detail::supported_electrostatics(variables[0], variables[1],
            std::vector<detail::AD>(x.begin(), x.end()), x, record);
        CppAD::ADFun<double> tape(variables, std::vector<detail::AD>{terms.second, terms.first.debye_huckel});
        const auto values = tape.Forward(0, point), gradient = tape.Jacobian(point);
        const auto born_hessian = tape.Hessian(point, 0), dh_hessian = tape.Hessian(point, 1);
        const auto inputs = fields({field("salt_fraction", salt), field("T_K", point[0]), field("rho_mol_m3", point[1]),
            field("composition", x), field("Na_Born_diameter_Angstrom", point[2]), field("neutral_permittivity", neutral.bulk_relative_permittivity)});
        if (salt == 0) {
            bool zero = std::all_of(values.begin(), values.end(), [](double v) { return v == 0; });
            for (const auto* rows : {&gradient, &born_hessian, &dh_hessian})
                zero = zero && std::all_of(rows->begin(), rows->end(), [](double v) { return v == 0; });
            tape.Forward(0, point); tape.Forward(1, std::vector<double>{1, 1, 1}); tape.Forward(2, std::vector<double>{0, 0, 0});
            const auto third = tape.Forward(3, std::vector<double>{0, 0, 0});
            zero = zero && third[0] == 0 && third[1] == 0;
            truth("neutral-electrostatics/zero", "Born and DH and admitted mixed restriction actions are finite exact zeros",
                zero && neutral.born == 0 && neutral.debye_huckel == 0, inputs, "issue-231");
            continue;
        }
        const double eps = neutral.bulk_relative_permittivity;
        const double constant = ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
            / (4 * PI * VACUUM_PERMITTIVITY_F_PER_M * BOLTZMANN_J_PER_K * point[0]);
        double radial = 0, radial_first = 0;
        for (std::size_t i : {1u, 2u}) {
            const double diameter = p.born_diameter_angstrom[i] * ANGSTROM_M;
            const double outer = diameter * (1 + p.born_c_shell * (p.solvation_factor[0] - 1));
            const double term = (1 - 1 / eps) / outer + p.born_c_dielectric * (1 - 1 / p.ionic_region_relative_permittivity)
                * (1 / diameter - 1 / outer);
            radial += term;
            if (i == 1) radial_first = term;
        }
        const double born = -constant * salt * radial, ion_born = -constant * salt * radial_first;
        const double screening = std::sqrt(point[1] * AVOGADRO_PER_MOL * ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
            / (VACUUM_PERMITTIVITY_F_PER_M * eps * BOLTZMANN_J_PER_K * point[0]));
        const double dh = -ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C * screening * std::pow(2 * salt, 1.5)
            / (12 * PI * VACUUM_PERMITTIVITY_F_PER_M * eps * BOLTZMANN_J_PER_K * point[0]);
        for (std::size_t i = 0; i < 2; ++i) {
            const double leading = i == 0 ? born : dh;
            const double error = std::abs(values[i] / leading - 1);
            truth("neutral-electrostatics/limit/" + number(salt) + (i == 0 ? "/Born" : "/DH"),
                "the independent leading coefficient converges with the declared dilute path",
                error < previous[i] && (salt != 1e-8 || error <= 1e-3),
                inputs + ";" + fields({field("value", values[i]), field("leading", leading), field("normalized_error", error)}), "issue-231");
            previous[i] = error;
        }
        if (salt == 1e-8) {
            const std::array<double, 6> references{-born / point[0], 2 * born / (point[0] * point[0]),
                -1.5 * dh / point[0], 3.75 * dh / (point[0] * point[0]), dh / (2 * point[1]), -dh / (4 * point[1] * point[1])};
            const std::array<double, 6> actions{gradient[0], born_hessian[0], gradient[3], dh_hessian[0], gradient[4], dh_hessian[4]};
            const std::array<const char*, 6> names{"Born/dT", "Born/dT2", "DH/dT", "DH/dT2", "DH/drho", "DH/drho2"};
            const std::array<const char*, 6> quantities{"d_a_Born_dT", "d2_a_Born_dT2", "d_a_DH_dT", "d2_a_DH_dT2", "d_a_DH_drho", "d2_a_DH_drho2"};
            const std::array<const char*, 6> units{"1/K", "1/K^2", "1/K", "1/K^2", "m^3/mol", "m^6/mol^2"};
            for (std::size_t i = 0; i < references.size(); ++i)
                agree(std::string("neutral-electrostatics/coefficient-action/") + names[i], quantities[i], units[i],
                    actions[i], references[i], 0, 1e-3, inputs, "exact first and second T and density actions approach independent coefficient derivatives", "issue-231");
            for (std::size_t order = 1; order <= 2; ++order) {
                const double reference = order == 1 ? -ion_born / point[2] : 2 * ion_born / (point[2] * point[2]);
                agree("neutral-electrostatics/born-diameter/" + std::to_string(order),
                    order == 1 ? "d_a_Born_d_dBorn_Na" : "d2_a_Born_d_dBorn_Na2", order == 1 ? "1/Angstrom" : "1/Angstrom^2",
                    order == 1 ? gradient[2] : born_hessian[8], reference, 0, 1e-3, inputs,
                    "Born diameter first and second actions approach independent coefficient derivatives", "issue-231");
            }
        }
        const auto changed = model.state(298.15, neutral.molar_density, x);
        const auto repeated = model.state_at_pressure(298.15, 1e5, {1, 0, 0}, DensityBranch::Liquid);
        const auto fresh = HelmholtzModel(p, active).state(298.15, neutral.molar_density, x);
        agree("neutral-electrostatics/support-tape/" + number(salt), "h_res", "J/mol", changed.residual_enthalpy,
            fresh.residual_enthalpy, 0, 1e-12, inputs, "alternating supports preserve fresh-model calorics", "issue-231");
        agree("neutral-electrostatics/neutral-tape/" + number(salt), "h_res", "J/mol", repeated.residual_enthalpy,
            neutral.residual_enthalpy, 0, 1e-12, inputs, "ionic evaluation cannot contaminate the neutral tape", "issue-231");
    }
}

}  // namespace regression_checks
