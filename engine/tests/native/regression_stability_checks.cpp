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
    // Stop an actual descending iterate in the second scheduled liquid trial, before convergence.
    stability_stop_trial_for_test = 2;
    Prediction interrupted;
    try { interrupted = predict(model, row); }
    catch (...) { stability_stop_trial_for_test.reset(); throw; }
    stability_stop_trial_for_test.reset();
    CHECK(interrupted.stability && !interrupted.stability->available && !interrupted.value);
    const auto& failed = interrupted.stability->report.trials.at(1);
    CHECK(!failed.optimizer.success && failed.optimizer.convergence.stop == core::Stop::budget_exhausted);
    CHECK(!failed.raw_coordinates.empty() && failed.resolved_descent);
    CHECK(failed.quantities.at("split_delta_g_over_rt") < -failed.quantities.at("split_error_estimate_over_rt"));
    std::cout << std::setprecision(17) << "{\"case\":\"failed-final-optimizer\",\"trial\":1,\"iteration\":"
        << failed.optimizer.iterations << ",\"success\":" << (failed.optimizer.success ? "true" : "false")
        << ",\"stop\":" << recording_support::json_string(core::stop_name(failed.optimizer.convergence.stop))
        << ",\"cause\":" << recording_support::json_string(failed.optimizer.convergence.cause)
        << ",\"stopping_residual\":" << failed.optimizer.convergence.residual << ",\"raw_coordinates\":[";
    for (std::size_t i = 0; i < failed.raw_coordinates.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << failed.raw_coordinates[i];
    }
    std::cout << "],\"composition\":[";
    for (std::size_t i = 0; i < failed.composition.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << failed.composition[i];
    }
    std::cout << "],\"quantities\":{";
    bool first = true;
    for (const auto& [name, value] : failed.quantities) { if (!first) std::cout << ','; first = false;
        std::cout << '"' << name << "\":" << value; }
    std::cout << "}}\n";
    auto tied = row; tied.problem.phases[0].composition_guess = {.999, .001};
    liquid.name = "L2"; liquid.composition_guess = {.3, .7}; tied.problem.phases.push_back(liquid);
    const auto tie = predict(model, tied);
    CHECK(tie.stability && !tie.stability->available && !tie.stability->event.empty() && !tie.value);
    CHECK(tie.stability->report.no_descent_in_trial_set && tie.stability->minima.size() >= 2);
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
