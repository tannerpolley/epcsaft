#include "regression_checks_shared.hpp"

#include "fixtures/association_stress.hpp"
#include "fixtures/hong_pathway_weight.hpp"

#include <filesystem>

namespace regression_checks {
namespace {
using namespace epcsaft::equilibrium;

Observation site_row(const EquilibriumProblem& problem, std::size_t phase, std::size_t site) {
    auto row = solved(ObservationKind::AssociationSiteFraction, problem, phase);
    row.site = site;
    return row;
}

EquilibriumResult center(const HelmholtzModel& model, const EquilibriumProblem& problem) {
    auto result = solve_equilibrium(model, problem);
    const auto context = fields({field("message", result.message), field("T_K", result.temperature.value_or(0)),
        field("P_Pa", result.pressure.value_or(0)), field("rho_mol_m3", result.molar_densities),
        field("amounts", result.amounts), field("x", result.mole_fractions)});
    CHECK(truth("equilibrium/admitted", result.success && result.phase_count == problem.phases.size()
        && result.topology_event == TopologyEvent::None, context));
    CHECK(agree("equilibrium/classified", "row norm", "1", result.convergence.residual, 0., 1e-10, 0., context));
    return result;
}

std::vector<double> phase_x(const EquilibriumResult& result, std::size_t phase) {
    const auto count = result.mole_fractions.size() / result.phase_count;
    const auto first = result.mole_fractions.begin() + phase * count;
    return {first, first + count};
}

double direct(const HelmholtzModel& model, const EquilibriumResult& state, const Observation& row) {
    return model.state(*state.temperature, state.molar_densities[row.phase], phase_x(state, row.phase)).site_fractions.at(*row.site);
}

double projection(const char* id, const HelmholtzModel& model, const Observation& row, const Prediction& prediction) {
    CHECK(prediction.status == ObservationStatus::Available && prediction.value && prediction.solved && prediction.unit == "1");
    const double reference = direct(model, *prediction.solved, row);
    CHECK(reference_scale_agrees(id, "site fraction", "1", *prediction.value, reference, 0., 32 * EPS,
        fields({field("phase", row.problem.phases[row.phase].name), field("site", *row.site)})));
    return reference;
}

EquilibriumProblem one_phase(double T, double P, const std::vector<double>& x, PhaseKind kind) {
    auto problem = equilibrium_support::flash(x.size());
    problem.temperature = FreeValue::fixed(T);
    problem.pressure = FreeValue::fixed(P);
    problem.feed.values = x;
    problem.phases.front().kind = kind;
    problem.solver_options.values[policy::row] = 1e-10;
    return problem;
}

void bubble_temperature() {
    using hong_pathway_weight::shifted;
    const auto parameters = properties_support::associating_binary();
    const ActiveParameter energy{ParameterFamily::AssociationEnergyOverK, {1, 4}};
    const HelmholtzModel model(parameters, {energy});
    auto problem = equilibrium_support::bubble_at(330);
    problem.solver_options.values[policy::row] = 1e-10;
    problem.pressure = FreeValue::fixed(*center(model, problem).pressure);
    problem.temperature = FreeValue::variable(330);
    const auto liquid = site_row(problem, 0, 0), vapor = site_row(problem, 1, 0);
    const auto low = predict_validated(model, liquid), high = predict_validated(model, vapor);
    const double xl = projection("bubble/liquid", model, liquid, low), xv = projection("bubble/vapor", model, vapor, high);
    CHECK(std::abs(xl - xv) > 32 * EPS * (xl + xv));
    auto omitted = vapor;
    omitted.site.reset();
    CHECK(recording_support::refusal([&] { (void)predict_validated(model, omitted); }).find("site selector") != std::string::npos);
    const auto action = differentiate(model, vapor, high);
    CHECK(action.status == ObservationStatus::Available && action.value.size() == 1);
    std::map<double, EquilibriumResult> states;
    const auto reference = [&](double step, bool hold_T) {
        const HelmholtzModel moved(shifted(parameters, energy, step));
        if (states.find(step) == states.end()) states.emplace(step, center(moved, problem));
        auto state = states.at(step);
        if (hold_T) state.temperature = high.solved->temperature;
        return direct(moved, state, vapor);
    };
    const auto rungs = ladder([&](double h) { return reference(h, false); }, 0.1);
    const double held = (reference(0.05, true) - reference(-0.05, true)) / 0.1;
    CHECK(ladder_agrees("bubble/energy", "dX/denergy", action.value[0], rungs, fields({field("h_K", 0.1)}), true));
    CHECK(truth("bubble/temperature-omission", std::abs(held - rungs.fine) > FIRST_ORDER * std::abs(rungs.fine),
        fields({field("held_slope", held), field("fine_slope", rungs.fine),
            field("gap", std::abs(held - rungs.fine)), field("allowance", FIRST_ORDER * std::abs(rungs.fine))})));

    SolvedStateObservable item{SolvedStateObservableKind::AssociationSiteFraction, 1, 0, {}, 0};
    SolvedStateActionDirection direction;
    direction.feed_amounts_mol = {0., 0.};
    direction.active_parameters = {1.};
    CHECK(recording_support::refusal([&] {
        (void)solved_state_actions(*high.compiled, *high.solved, {{item}, {direction, direction}});
    }).find("first actions only") != std::string::npos);
}

void pathway_mapping() {
    using namespace hong_pathway_weight;
    FitRequest request;
    request.parameters = hong(Form::Blend);
    request.controls = budgeted_controls();
    const auto active = coordinates(Form::Blend);
    const HelmholtzModel model(request.parameters, active);
    const auto& liquid = STATES.front();
    const ActiveParameter weight{ParameterFamily::PathwayWeight, {0}};
    request.coordinates.coordinates = {{{{weight.family, weight.indices, 1.}}, 1.7, 2., 0., 3., {}}};
    std::vector<Ladder> columns;
    for (const double density : {liquid.density, 42450.0}) {
        const auto problem = one_phase(liquid.temperature, model.state(liquid.temperature, density, liquid.x).pressure,
            liquid.x, PhaseKind::Liquid);
        const auto state = center(model, problem);
        auto row4 = site_row(problem, 0, 4), row9 = site_row(problem, 0, 9);
        row4.form = ResidualForm::Difference;
        row4.scale = 0.4;
        row4.observed = direct(model, state, row4);
        row9.observed = direct(model, state, row9);
        const auto first = request.observations.size();
        if (first == 0) request.observations.insert(request.observations.end(), {row4, row9});
        else request.observations.insert(request.observations.end(), {row9, row4});
        std::map<double, std::array<double, 2>> values;
        const auto weight_values = [&](double h) -> const auto& {
            if (values.find(h) == values.end()) {
                const HelmholtzModel moved(shifted(request.parameters, weight, h));
                const auto moved_state = center(moved, problem);
                values.emplace(h, std::array{direct(moved, moved_state, row4), direct(moved, moved_state, row9)});
                if (first == 0)
                    CHECK(direct(moved, state, row4) == row4.observed && direct(moved, state, row9) == row9.observed);
            }
            return values.at(h);
        };
        for (std::size_t i = first; i < first + 2; ++i) {
            const auto& row = request.observations[i];
            const auto pred = predict_validated(model, row);
            projection(("pathway/site-" + std::to_string(*row.site)).c_str(), model, row, pred);
            const double factor = 2. / (*row.site == 4 ? row.scale : row.observed);
            columns.push_back(ladder([&](double h) { return factor * weight_values(h)[*row.site == 4 ? 0 : 1]; }, 0.001));
            if (first != 0) continue;
            const auto derivative = differentiate(model, row, pred);
            CHECK(derivative.status == ObservationStatus::Available && derivative.value.size() == active.size());
            for (std::size_t c = 0; c < active.size(); ++c) {
                const auto& parameter = active[c];
                const auto slot = parameter.indices.front();
                const auto& fixed = detail::family_values(request.parameters, parameter.family).first;
                double step = 0.001 * std::abs(fixed[slot]);
                if (c == 13) step = 0.001;
                const auto rungs = ladder([&](double h) {
                    const HelmholtzModel moved(shifted(request.parameters, parameter, h));
                    return direct(moved, center(moved, problem), row);
                }, step);
                CHECK(ladder_agrees("pathway/site-" + std::to_string(*row.site) + "/coordinate-" + std::to_string(c),
                    "dX/dparameter", derivative.value[c], rungs, fields({field("slot", slot), field("h", step)}), true));
            }
        }
        CHECK(std::abs(columns[first].fine - columns[first + 1].fine) >
            FIRST_ORDER * (std::abs(columns[first].fine) + std::abs(columns[first + 1].fine)));
        auto density_row = solved(ObservationKind::PhaseDensity, problem);
        density_row.observed = state.molar_densities.front();
        request.observations.push_back(density_row);
        columns.emplace_back();  // This observation shares the solve, but has no new fitted derivative claim.
    }
    const auto check = [&](const std::string& path, const auto& jacobian) {
        CHECK(jacobian.size() == request.observations.size());
        for (std::size_t i = 0; i < jacobian.size(); ++i) if (request.observations[i].site)
            CHECK(ladder_agrees("pathway/" + path + "/row-" + std::to_string(i), "dr/doptimizer",
                jacobian[i], columns[i], fields({field("site", *request.observations[i].site), field("h", 0.001)}), true));
    };
    const auto evaluated = evaluated_columns(request, {0.});
    check("point", evaluated.residual_optimizer_gradient);
    request.controls.workers = 1;
    const auto serial = fit(request);
    request.controls.workers = 2;
    const auto worker = fit(request);
    CHECK(serial.usable && worker.usable && serial.optimizer == std::vector<double>{0.} && worker.optimizer == serial.optimizer);
    check("serial", serial.optimizer_jacobian);
    check("worker", worker.optimizer_jacobian);
    CHECK(serial.predictions == worker.predictions && serial.optimizer_jacobian == worker.optimizer_jacobian);
    CHECK(serial.rank_records.back().rank == 1 && worker.rank_records.back().rank == 1);
}

void reacting_actions() {
#include "fixtures/mea-published-state-19.inc"
    // The #356 flash uses the same packet-3 coefficients with both Born switches off.
    auto packet = model.parameters();
    packet.born_c_shell = packet.born_c_dielectric = 0.;
    problem.temperature = FreeValue::fixed(318.15);
    problem.phases[0].support.clear();
    problem.phases[0].support_declared = false;
    problem.pressure = FreeValue::fixed(3699.0803800848512);
    problem.feed.values[0] = 0.25557798094241313 * problem.feed.values[1];
    problem.source_reference.reset();
    problem.solver_options.values[policy::row] = 1e-10;
    const std::array ln_k{-605.348, -587.8, -948.2312, -25.1014, -298.6415};
    CHECK(problem.reactions.size() == ln_k.size());
    for (std::size_t i = 0; i < problem.reactions.size(); ++i) {
        auto& reaction = problem.reactions[i];
        reaction.reference.reset();
        reaction.correlation = ReactionLogPolynomial{ln_k[i], 0., 0., 0., 318.15, 250., 450., EOS_STANDARD_STATE_ID};
    }
    const ActiveParameter energy{ParameterFamily::DispersionEnergyOverK, {2}};
    const HelmholtzModel active(packet, {energy});
    const auto water = site_row(problem, 0, 2);
    const auto prediction = predict_validated(active, water);
    projection("reacting/water", active, water, prediction);
    CHECK(reference_scale_agrees("reacting/356-liquid-density", "density", "mol/m³",
        prediction.solved->molar_densities[0], 42120.183730573895, 0., 1e-10, "{}"));
    const auto action = differentiate(active, water, prediction);
    CHECK(action.status == ObservationStatus::Available && action.value.size() == 1);
    std::map<double, EquilibriumResult> states;
    const auto displaced_value = [&](double h, bool hold_rho, bool hold_x) {
        const HelmholtzModel moved(hong_pathway_weight::shifted(packet, energy, h));
        if (states.find(h) == states.end()) states.emplace(h, center(moved, problem));
        auto state = states.at(h);
        if (hold_rho) state.molar_densities = prediction.solved->molar_densities;
        if (hold_x) state.mole_fractions = prediction.solved->mole_fractions;
        return direct(moved, state, water);
    };
    const auto rungs = ladder([&](double h) { return displaced_value(h, false, false); }, 0.1);
    CHECK(ladder_agrees("reacting/dispersion", "dX/denergy", action.value[0], rungs, fields({field("h_K", 0.1)}), true));
    for (const bool hold_rho : {false, true}) {
        const double held = (displaced_value(0.05, hold_rho, !hold_rho) - displaced_value(-0.05, hold_rho, !hold_rho)) / 0.1;
        const double gap = std::abs(held - rungs.fine), allowance = FIRST_ORDER * std::abs(rungs.fine);
        CHECK(truth(hold_rho ? "reacting/density-omission" : "reacting/composition-omission", gap > allowance,
            fields({field("held_slope", held), field("fine_slope", rungs.fine), field("gap", gap), field("allowance", allowance)})));
    }
}

void declaration_refusals() {
    FitRequest request;
    request.controls = budgeted_controls();
    const ActiveParameter energy{ParameterFamily::DispersionEnergyOverK, {0}};
    request.coordinates.coordinates = {{{{energy.family, energy.indices, 1.}}, 191.42, 1., 150., 250., {}}};
    auto problem = one_phase(180., 1e6, {1., 0.}, PhaseKind::Liquid);
    auto& phase = problem.phases.front();
    phase.support = {0};
    phase.support_declared = true;
    const auto no_sites = properties_support::identical_ethane_pair();
    auto inactive = no_sites;
    inactive.site_component = {0, 1};
    inactive.site_multiplicity = {1, 1};
    inactive.association_active = {0, 1, 1, 0};
    inactive.association_energy_over_k = {0., 1804.17, 1804.17, 0.};
    inactive.association_volume = {0., 0.2039, 0.2039, 0.};
    const HelmholtzModel model(inactive, {energy});
    auto row = site_row(problem, 0, 0);
    const auto prediction = predict_validated(model, row);
    CHECK(projection("inactive/site-zero", model, row, prediction) == 1.);
    const auto derivative = differentiate(model, row, prediction);
    CHECK(truth("inactive/in-support-action", derivative.status == ObservationStatus::Available
        && derivative.value == std::vector<double>{0.}, fields({field("dX/denergy", derivative.value)})));
    for (const auto ordinal : {std::optional<std::size_t>{}, std::optional<std::size_t>{0}, {1}, {2}}) {
        request.parameters = ordinal == 0 ? no_sites : inactive;
        row.site = ordinal;
        request.observations = {solved(ObservationKind::PhaseDensity, problem), row};
        const auto result = fit(request);
        const std::string expected = !ordinal ? "selector" : *ordinal == 1 ? "phase support" : "ordinal";
        CHECK(truth("declaration/refused", result.status == "invalid_setup" && result.residual_evaluations == 0
            && result.jacobian_evaluations == 0 && result.message.find("observation 1:") == 0
            && result.message.find(expected) != std::string::npos,
            fields({field("selector", ordinal ? std::to_string(*ordinal) : "omitted"), field("message", result.message)})));
    }
}

void boundary_refusals() {
    using Status = SolvedStateActionStatus;
    for (const auto& [sites, q, status] : {std::tuple{2u, 1e26, Status::AssociationIllConditioned},
                                        std::tuple{128u, 1e27, Status::AssociationSingular}}) {
        const auto stressed = properties_support::association_stress(sites, q);
        auto weak = stressed;
        weak.association_energy_over_k[0] = 0.5;
        const HelmholtzModel regular(weak);
        const auto problem = one_phase(300., regular.state(300., 1000., {1.}).pressure, {1.}, PhaseKind::Vapor);
        const auto point = center(regular, problem);
        const auto row = site_row(problem, 0, 0);
        const ActiveParameter energy{ParameterFamily::DispersionEnergyOverK, {0}};
        const HelmholtzModel model(stressed, {energy});
        Prediction prediction;
        prediction.solved = point;
        prediction.compiled = compile_problem(regular, problem);
        // The weak model supplies the coordinates and compiled declaration for this boundary check.
        prediction.compiled->model = &model;
        project_prediction(model, row, prediction);
        const double exact = 2. / (1. + std::sqrt(1. + 4. * q));
        CHECK(prediction.value && prediction.residual);
        CHECK(reference_scale_agrees("refusal/stress-value", "site fraction", "1", *prediction.value, exact,
            0., 2e-12, fields({field("q", q), field("sites", sites)})));
        const auto action = differentiate(model, row, prediction);
        CHECK(truth("refusal/stress-action", action.status == ObservationStatus::ActionUnavailable
            && action.central_value == prediction.value && action.value.empty()
            && action.message.find(solved_state_status_name(status)) == 0,
            fields({field("cause", action.message), field("q", q), field("sites", sites)})));
        if (sites != 2) continue;

        weak.association_energy_over_k[0] = 3e5;
        const HelmholtzModel overflow(weak);
        std::optional<core::Convergence> failure;
        try {
            overflow.evaluate(*point.temperature, 1. / point.molar_densities[0], {1.}, DerivativeOrder::Value);
        } catch (const core::ConvergenceFailure& error) {
            failure = error.record;
        }
        CHECK(failure && failure->stop == core::Stop::nonfinite && failure->criterion
            && std::string_view(failure->criterion->name) == policy::ASSOCIATION_RESIDUAL.name);
        Prediction refused;
        refused.solved = point;
        project_prediction(overflow, row, refused);
        CHECK(truth("refusal/projection", !refused.value && !refused.residual && refused.convergence
            && core::describe(*refused.convergence) == core::describe(*failure),
            fields({field("cause", refused.message), field("convergence", core::describe(*failure))})));
    }
}

void fit_failure_records() {
    FitRequest request;
    auto& parameters = request.parameters = properties_support::associating_binary();
    // The fitted energy slope leaves the 330 K bubble invariant. The 329 K site
    // target forces a trial with overflowing association strength at its bound.
    parameters.temperature_law_reference_k = {330.};
    parameters.temperature_law_term_law = {0, 0};
    parameters.temperature_law_term_power = {0, 1};
    parameters.temperature_law_coefficient = {parameters.association_energy_over_k[1], 0.};
    parameters.temperature_law_coefficient_identity = {"energy", "slope"};
    parameters.temperature_law_slot = {{int(TemperatureLaw::AssociationEnergyOverK), 1, 0}};
    request.coordinates.coordinates = {{{{ParameterFamily::TemperatureLawCoefficient, {1}, 1.}},
        0., 1000., -1e6, 0., {}}};
    request.controls = budgeted_controls();
    request.controls.maximum_iterations = 3;
    // This test-only target below the rounding floor retains the bubble's
    // recovered density-scan diagnostic. The normal association target is unchanged.
    request.controls.solver_options.values[policy::association] = 2e-16;
    const HelmholtzModel model(parameters, request.coordinates.active_parameters(), {}, request.controls.solver_options);
    auto first = solved(ObservationKind::EquilibriumPressure, equilibrium_support::bubble_at(330.));
    const auto reference = predict_validated(model, first);
    CHECK(reference.value && reference.first_numerical_failure
        && reference.first_numerical_failure->stop == core::Stop::converged_to_floor);
    first.observed = *reference.value;
    auto cell = equilibrium_support::flash(2);
    cell.temperature = FreeValue::fixed(329.);
    cell.pressure = FreeValue::fixed(1e4);
    auto row = site_row(cell, 0, 0);
    row.observed = 1e-250;
    const auto initial = predict_validated(model, row);
    CHECK(initial.value);
    request.observations = {first, row};
    for (const int workers : {1, 2}) {
        request.controls.workers = workers;
        const auto result = fit(request);
        CHECK(result.optimizer == std::vector<double>{0.} && result.predictions.at(1) == *initial.value);
        CHECK(result.first_numerical_failure && !result.failures.empty());
        const auto& failure = result.failures.front();
        const HelmholtzModel trial(request.coordinates.apply(parameters, failure.optimizer), request.coordinates.active_parameters(),
            {}, request.controls.solver_options);
        const auto refused = predict(trial, row);
        CHECK(!refused.value && !refused.residual && refused.convergence && refused.convergence->criterion);
        CHECK(std::string_view(refused.convergence->criterion->name) == policy::ASSOCIATION_RESIDUAL.name
            && refused.convergence->stop == core::Stop::nonfinite);
        const auto earlier = core::describe(*reference.first_numerical_failure);
        const auto terminal = core::describe(*refused.convergence);
        CHECK(truth("fit-failure/records", failure.observation == 1 && !failure.jacobian && !failure.central_value
            && failure.status == refused.status && failure.message == refused.message && failure.convergence
            && core::describe(*failure.convergence) == terminal
            && failure.first_numerical_failure && core::describe(*failure.first_numerical_failure) == earlier
            && core::describe(*result.first_numerical_failure) == earlier && earlier != terminal,
            fields({field("workers", workers), field("row", failure.observation), field("first", earlier),
                field("terminal", terminal), field("optimizer", failure.optimizer)})));
    }

    // A finite site value can survive even when the scaled residual requested with
    // the initial Jacobian overflows. The worker must reply with that row's value.
    row.observed = 1.;
    row.form = ResidualForm::Difference;
    row.scale = std::numeric_limits<double>::min() / 16.;
    request.observations[1] = row;
    for (const int workers : {1, 2}) {
        request.controls.workers = workers;
        const auto result = fit(request);
        CHECK(result.status == "initial_evaluation_failed" && result.failures.size() == 1);
        const auto& failure = result.failures.front();
        CHECK(truth("fit-failure/central-value", failure.observation == 1 && failure.jacobian
            && failure.status == ObservationStatus::Nonfinite && failure.central_value == initial.value
            && result.residuals.empty() && result.optimizer_jacobian.empty(),
            fields({field("workers", workers), field("row", failure.observation), field("site_fraction", *initial.value)})));
    }
}

}  // namespace
}  // namespace regression_checks

int main(int argc, char** argv) {
    using namespace regression_checks;
    if (argc >= 2 && std::string_view(argv[1]) == "--fit-worker") {
#ifdef _WIN32
        serve_fit(std::stoull(argv[2]));
#else
        serve_fit(0);
#endif
    }
    set_fit_host({{std::filesystem::absolute(argv[0]).string(), "--fit-worker"}}, reinterpret_cast<const void*>(&bubble_temperature));
    return recording_support::run(argc, argv, {
        {{"native.regression.association-site-temperature", "temperature"}, bubble_temperature},
        {{"native.regression.association-site-mapping", "mapping"}, pathway_mapping},
        {{"native.regression.association-site-reaction", "reaction"}, reacting_actions},
        {{"native.regression.association-site-declaration", "declaration"}, declaration_refusals},
        {{"native.regression.association-site-refusals", "refusals"}, boundary_refusals},
        {{"native.regression.association-site-fit-failure", "fit-failure"}, fit_failure_records},
    });
}
