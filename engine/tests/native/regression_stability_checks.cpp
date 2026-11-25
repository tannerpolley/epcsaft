#include "regression_checks_shared.hpp"
#include "water_isobutanol.hpp"

namespace regression_checks {
void check_stability_requirement() {
    using namespace epcsaft::equilibrium;
    HelmholtzModel model(water_isobutanol());
    Observation row; row.kind = ObservationKind::Stability; row.scale = .1;
    auto& p = row.problem;
    p.temperature = FreeValue::fixed(293.15); p.pressure = FreeValue::fixed(1e5);
    p.feed.kind = FeedKind::Amounts; p.feed.values = {.6, .4};
    p.solver_options.values[policy::row] = 1e-10;
    PhaseSpec liquid; liquid.name = "L"; liquid.kind = PhaseKind::Liquid;
    liquid.support = {0, 1}; liquid.support_declared = true; liquid.composition_guess = {.6, .4};
    p.phases = {liquid};
    auto initial = predict(model, row);
    CHECK(initial.stability && initial.stability->available);
    CHECK(*initial.stability->value < -*initial.stability->uncertainty);
    CHECK(!initial.stability->report.no_descent_in_trial_set);
    CHECK(initial.solved->residual_norm() <= 1e-10);
    auto supplied = liquid;
    supplied.packing_guess = model.packing_fraction_coefficient(293.15, liquid.composition_guess) * initial.solved->molar_densities.at(0);
    const auto bank = prepare_supplied_trials(model, 293.15, {supplied});
    const auto both = predict_validated(model, row, bank);
    CHECK(both.stability && both.stability->raw_minimum == initial.stability->raw_minimum);
    CHECK(both.supplied_available && both.supplied_stability->trials.size() == 1);
    CHECK(both.supplied_stability->evidence_counts({3}).at("evaluated") == 3);
    StabilityOptions options; options.prepared_trials = bank;
    const auto separate = stability_test(*initial.compiled, *initial.solved, options);
    CHECK(both.supplied_stability->trials[0].quantities == separate.trials.at(0).quantities);
    {
        auto uncertain = row;
        uncertain.stability_options.standard_uncertainty_over_rt = {.001, .001};
        auto vapor = liquid; vapor.kind = PhaseKind::Vapor;
        options.prepared_trials = prepare_supplied_trials(model, 293.15, {liquid, vapor});
        options.standard_uncertainty_over_rt = uncertain.stability_options.standard_uncertainty_over_rt;
        const auto rooted = predict_validated(model, uncertain, options.prepared_trials);
        const auto independent = stability_test(*initial.compiled, *initial.solved, options);
        CHECK(rooted.supplied_available && rooted.supplied_stability->supplied_minima().size() == 2);
        for (std::size_t i = 0; i < 2; ++i) {
            const auto& q = rooted.supplied_stability->trials.at(i).quantities;
            CHECK(q == independent.trials.at(i).quantities && q.at("standard_error_estimate_over_rt") > .0019);
            CHECK(q.at("density_side") == static_cast<double>(i ? DensityBranch::Vapor : DensityBranch::Liquid));
            std::cout << std::setprecision(17) << "supplied root " << q.at("density_side") << " density mol/m3 "
                << q.at("molar_density_mol_m3") << " TPD/RT " << q.at("tpd_over_rt") << " standard error "
                << q.at("standard_error_estimate_over_rt") << " total error " << q.at("descent_error_estimate_over_rt") << '\n';
        }
        auto water = liquid; water.support = {0}; water.composition_guess = {1., 0.};
        auto alcohol = liquid; alcohol.support = {1}; alcohol.composition_guess = {0., 1.};
        auto pure = row; pure.problem.feed.values = {1., 0.}; pure.problem.phases = {water};
        const auto mixed = predict_validated(model, pure, prepare_supplied_trials(model, 293.15, {alcohol, water}));
        CHECK(!mixed.supplied_available && mixed.supplied_stability->completed);
        const auto& points = mixed.supplied_stability->trials;
        CHECK(points.size() == 2 && points[0].evidence == "not_applicable" && points[1].evidence == "evaluated");
        CHECK(mixed.supplied_stability->unresolved_reasons().at("stability support outside fed/reaction union") == 1);
        CHECK(points[1].quantities == predict_validated(model, pure, prepare_supplied_trials(model, 293.15, {water}))
            .supplied_stability->trials.at(0).quantities);
    }
    auto changed = initial;
    changed.stability = std::make_shared<StabilityMeasure>(*initial.stability);
    auto& changed_trial = changed.stability->report.trials.at(*changed.stability->minimizing_trial);
    changed_trial.quantities["density_side"] = static_cast<double>(DensityBranch::Vapor);
    const auto branch_event = stability_transition(0, {}, initial, {}, changed);
    CHECK(initial.stability->minimizing_trial == changed.stability->minimizing_trial);
    CHECK(branch_event.kind.find("minimizer") != std::string::npos);
    std::cout << "{\"case\":\"same-index-root-change\",\"event\":"
        << recording_support::json_string(branch_event.kind) << "}\n";
    auto capped = row; capped.stability_options.trial_budget = 0;
    const auto incomplete = predict_validated(model, capped, bank);
    CHECK(incomplete.supplied_available && incomplete.supplied_stability->trials[0].quantities == separate.trials[0].quantities);
    CHECK(incomplete.stability && !incomplete.stability->available && !incomplete.value);
    auto tied = row; tied.problem.phases[0].composition_guess = {.999, .001};
    liquid.name = "L2"; liquid.composition_guess = {.3, .7}; tied.problem.phases.push_back(liquid);
    const auto tie = predict(model, tied);
    CHECK(tie.stability && !tie.stability->available && !tie.stability->event.empty() && !tie.value);
    CHECK(tie.stability->report.no_descent_in_trial_set && tie.stability->minima.size() >= 2);
    CHECK(differentiate(model, tied, tie).status != ObservationStatus::Available);
    std::cout << std::setprecision(17) << "distinct tie " << *tie.stability->raw_minimum << ' '
        << *tie.stability->uncertainty << ' ' << tie.stability->minima.size() << '\n';
    {
        auto reference = row; reference.kind = ObservationKind::EquilibriumPressure;
        validate(reference, model);
        const auto pressure = predict(model, reference);
        CHECK(pressure.status == ObservationStatus::Available);
        ReferenceCache references; references.pressures.emplace_back(reference, pressure);
        auto shared = tied; shared.pressure_reference = reference.problem; shared.phase = 1;
        CHECK(validate_declaration(shared, model));
        CHECK(recording_support::refusal([&] { validate(shared, model); }) == "observation phase index");
        CHECK(recording_support::refusal([&] { (void)predict(model, shared, {}, {}, &references, nullptr, nullptr, false,
            std::chrono::steady_clock::time_point::max(), {&shared}); }) == "observation phase index");
    }
    auto uncertain = row; uncertain.stability_options.standard_uncertainty_over_rt = {1., 1.};
    const auto transition = predict(model, uncertain, {}, {}, nullptr, &initial, &model, true);
    CHECK(transition.stability && transition.stability->report.no_descent_in_trial_set);
    const auto event = stability_transition(0, {}, initial, {}, transition);
    const auto action = differentiate(model, uncertain, transition);
    CHECK(event.kind.find("resolution") != std::string::npos);
    CHECK(action.status != ObservationStatus::Available);
    std::cout << "{\"case\":\"resolution-transition\",\"event\":" << recording_support::json_string(event.kind)
        << ",\"previous_classification\":" << recording_support::json_string(event.previous_classification)
        << ",\"classification\":" << recording_support::json_string(event.classification)
        << ",\"previous_tau\":" << *initial.stability->raw_minimum
        << ",\"tau\":" << *transition.stability->raw_minimum
        << ",\"previous_m\":" << *initial.stability->value
        << ",\"available\":" << (transition.stability->available ? "true" : "false")
        << ",\"action_status\":" << recording_support::json_string(observation_status_name(action.status)) << "}\n";
}

void check_nested_support_tangent_plane_action() {
    using namespace epcsaft::equilibrium;
    FitRequest request;
    auto& p = request.parameters = properties_support::parameters(2);
    p.segment_count = {1.6069, 1.6069};
    p.sigma_constant = {3.5206, 3.5206};
    p.dispersion_energy_over_k = {191.42, 191.42};
    p.molar_mass_kg_per_mol = {.03007, .03007};
    request.coordinates.coordinates = {
        {{{ParameterFamily::DispersionEnergyOverK, {0}}, {ParameterFamily::DispersionEnergyOverK, {1}}}, 191.42, 1.0, 150.0, 250.0, {}}
    };
    Observation row;
    row.kind = ObservationKind::Stability;
    auto& problem = row.problem;
    problem.temperature = FreeValue::fixed(200.0);
    problem.pressure = FreeValue::fixed(2e5);
    problem.feed.values = {1.0, 0.0};
    problem.reactions = {{"A=B", {-1.0, 1.0}, 1e-15, {}, {}, {}}};
    problem.solver_options.values[policy::row] = 1e-10;
    PhaseSpec liquid;
    liquid.name = "L";
    liquid.support = {0, 1};
    liquid.support_declared = true;
    problem.phases = {liquid};
    auto vapor = liquid;
    vapor.name = "V";
    vapor.kind = PhaseKind::Vapor;
    auto restricted = vapor;
    restricted.support = {0};
    row.stability_options.trials = std::vector{vapor, restricted};
    request.observations = {row};
    HelmholtzModel model(p, request.coordinates.active_parameters());
    const auto trial_set = [&](const Observation& observation, const HelmholtzModel& eos, double step) {
        const auto prediction = predict(eos, observation);
        if (!prediction.value) std::cout << prediction.message << std::endl;
        CHECK(prediction.stability && prediction.stability->available && prediction.value);
        const auto& measure = *prediction.stability;
        for (auto index : measure.minima) {
            const auto& trial = measure.report.trials[index];
            const auto& q = trial.quantities;
            CHECK(trial.resolved_descent && q.at("split_delta_g_over_rt") + q.at("split_error_estimate_over_rt") < 0.0);
            CHECK(q.at("root_density_side") == static_cast<double>(DensityBranch::Vapor));
            CHECK(std::abs(q.at("normalization_residual")) <= policy::COMPOSITION_SUM_TOLERANCE.value);
            CHECK(trial.composition[1] <= policy::ADDITION_MATERIAL_BALANCE.value);
            CHECK(q.at("reconstructed_stationarity_residual") <= trial.optimizer.convergence.criterion->value);
            std::cout << std::setprecision(17) << "trial " << step << ' ' << index << ' ' << trial.phase.support.size()
                << ' ' << q.at("root_density_side") << ' ' << trial.resolved_descent << ' ' << q.at("tpd_over_rt")
                << ' ' << q.at("plane_value_error_estimate_over_rt") << ' ' << q.at("reconstructed_stationarity_residual")
                << ' ' << q.at("pressure_row") << ' ' << q.at("split_delta_g_over_rt") << ' '
                << q.at("split_error_estimate_over_rt") << ' ' << trial.composition[1] << '\n';
        }
        return prediction;
    };
    const auto center = trial_set(row, model, 0.0);
    const auto direct = [&](const Observation& observation, const Prediction& prediction, const CoordinateSet& coordinates) {
        return differentiate_group(model, {&observation}, prediction, coordinates).front();
    };
    const auto action = direct(row, center, request.coordinates);
    CHECK(action.status == ObservationStatus::Available);
    const double expected = std::accumulate(action.value.begin(), action.value.end(), 0.0);
    CHECK(std::abs(expected - 0.0396268846633831) <= FIRST_ORDER);
    for (const auto& support : {vapor, restricted}) {
        auto single = row;
        single.stability_options.trials = std::vector{support};
        const auto prediction = trial_set(single, model, 0.0);
        const double resolution = *center.stability->uncertainty + *prediction.stability->uncertainty;
        CHECK(std::abs(*center.value - *prediction.value) <= resolution);
        const auto gradient = direct(single, prediction, request.coordinates);
        CHECK(gradient.status == ObservationStatus::Available);
        const double exact = std::accumulate(gradient.value.begin(), gradient.value.end(), 0.0);
        CHECK(std::abs(exact - expected) <= resolution);
        std::cout << "action " << support.support.size() << ' ' << exact << ' ' << expected << ' ' << resolution << '\n';
    }
    double prior = expected;
    for (double h : {0.005, 0.0025, 0.00125}) {
        const auto endpoint = [&](double u) {
            HelmholtzModel shifted(request.coordinates.apply(p, {u}), request.coordinates.active_parameters());
            return trial_set(row, shifted, u);
        };
        const auto plus = endpoint(h), minus = endpoint(-h);
        const double centered = (*plus.value - *minus.value) / (2.0 * h);
        CHECK(std::abs(centered - expected) <= FIRST_ORDER);
        std::cout << "centered " << h << ' ' << centered << ' ' << centered - expected << ' ' << centered - prior
            << ' ' << (*plus.stability->uncertainty + *minus.stability->uncertainty) / (2.0 * h) << '\n';
        prior = centered;
    }
    std::reverse(row.stability_options.trials->begin(), row.stability_options.trials->end());
    const auto reversed = trial_set(row, model, 0.0);
    CHECK(std::abs(*reversed.value - *center.value) <= *reversed.stability->uncertainty + *center.stability->uncertainty);
    CHECK(direct(row, reversed, request.coordinates).status == ObservationStatus::Available);

    const auto same_point = [&](bool admitted) {
        FitResult failures(request);
        ObservationEvaluator evaluator(request, failures, {0});
        ObservationEvaluation value, jacobian;
        CHECK(evaluator.evaluate({0.0}, false, value));
        CHECK(evaluator.evaluate({0.0}, true, jacobian) == admitted);
        if (admitted) {
            CHECK(std::abs(jacobian.optimizer_gradient[0] - expected) <= *center.stability->uncertainty);
        } else {
            CHECK(failures.failures.back().status == ObservationStatus::ActionUnavailable);
            CHECK(failures.failures.back().message.find("reduced-value gap") != std::string::npos);
            std::cout << failures.failures.back().message << '\n';
        }
    };
    same_point(true);
    request.controls.maximum_iterations = 1;
    request.controls.maximum_elapsed_time_seconds = 15.0;
    {
        core::Capture trace(true);
        const auto complete = fit(request);
        if (!complete.initial_cost) std::cout << complete.message << std::endl;
        CHECK(complete.initial_cost && complete.jacobian_evaluations > 0);
        const auto& events = trace.trace().events;
        const auto reused = std::count_if(events.begin(), events.end(), [&](const auto& event) {
            return std::string_view(event.name) == "reuse" && event.values == std::vector{3.0}
                && std::string_view(trace.trace().spans[event.span].name) == "fit_evaluation";
        });
        std::cout << "fit same-point Jacobian requests " << reused << '\n';
    }
    request.coordinates.coordinates = {{{}, std::log(1e-15), 1000.0, -40.0, -30.0, ReactionCoordinate{"A=B"}}};
    same_point(false);
}
}  // namespace regression_checks
