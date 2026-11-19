#include "regression_checks_shared.hpp"
#include "water_isobutanol.hpp"

namespace regression_checks {
void check_stability_pressure_reference() {
    using namespace epcsaft::equilibrium;
    constexpr double criterion = 1e-16;  // caller override from #85
    EquilibriumProblem reference;
    reference.temperature = FreeValue::fixed(584.75);
    reference.pressure = FreeValue::variable(1e6);
    reference.feed.kind = FeedKind::Amounts;
    reference.feed.values = {0.5, 0.5};
    reference.phases = {equilibrium_support::phase("L", PhaseKind::Liquid),
        equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0))};
    reference.solver_options.values[policy::row] = criterion;

    Observation parent;
    parent.kind = ObservationKind::Stability;
    parent.pressure_reference = reference;
    parent.problem = reference;
    parent.problem.pressure = FreeValue::fixed(1e7);
    parent.problem.phases.resize(1);
    parent.problem.solver_options = {};
    parent.pressure_ratio = 1.05;
    auto higher = parent;
    higher.pressure_ratio = 3.0;

    FitRequest request;
    request.parameters = properties_support::associating_binary();
    request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 204.0438}};
    request.observations = {parent, higher};
    request.weights = {1.0, 1.0};
    request.controls.maximum_iterations = 1;
    request.controls.maximum_elapsed_time_seconds = 30.0;
    request.controls.workers = 1;
    const auto result = fit(request);
    CHECK(result.acceptance == "pass");
    const HelmholtzModel replay(result.parameters, {}, {}, request.controls.solver_options, true);
    const auto cold = solve_compiled(compile_problem(replay, reference));
    CHECK(cold.has_state() && cold.convergence.residual <= criterion);
    for (std::size_t i = 0; i < request.observations.size(); ++i) {
        const auto& prediction = result.stability_predictions.at(i).value();
        CHECK(prediction.value == 0.0 && prediction.stability->report.completed);
        CHECK(prediction.pressure_reference->value == cold.pressure);
        CHECK(prediction.pressure_reference->convergence->residual <= criterion);
        CHECK(prediction.solved->pressure == request.observations[i].pressure_ratio * *cold.pressure);
    }
}

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
    const auto incomplete = predict(model, capped);
    CHECK(incomplete.stability && !incomplete.stability->available && !incomplete.value);
    auto tied = row; tied.problem.phases[0].composition_guess = {.999, .001};
    liquid.name = "L2"; liquid.composition_guess = {.3, .7}; tied.problem.phases.push_back(liquid);
    const auto tie = predict(model, tied);
    CHECK(tie.stability && !tie.stability->available && !tie.stability->event.empty() && !tie.value);
    CHECK(tie.stability->report.no_descent_in_trial_set && tie.stability->minima.size() >= 2);
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
}  // namespace regression_checks
