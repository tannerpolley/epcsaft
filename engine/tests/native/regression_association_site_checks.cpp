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

void kim_pure_water() {
#include "fixtures/kim_pure_water.inc"
    using Family = ParameterFamily;
    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = build_parameters(document);
    request.coordinates.coordinates = {
        {{{Family::TemperatureLawCoefficient, {0}}}, -0.06711345443725586, 0.01, -0.5, -0.06711345443725586, {}},
        {{{Family::TemperatureLawCoefficient, {1}}}, 0.0012092080101259242, 0.001, -0.01, 0.01, {}},
        {{{Family::Kij, {7, 56}}}, 0.25261865658522437, 0.01, -0.5, 0.5, {}},
        {{{Family::Kij, {15, 57}}}, -0.00956762026107625, 0.01, -0.5, 0.5, {}},
        {{{Family::KijReciprocalTemperatureSlope, {15, 57}}}, -15.929380742416088, 10., -1000., 1000., {}},
        {{{Family::Kij, {23, 58}}}, 1., 0.01, -1., 1., {}},
        {{{Family::SigmaConstant, {7}}, {Family::DebyeHuckelDiameter, {7}}, {Family::BornDiameter, {7}}},
            3.6552722619123172, 0.1, 3., 4.5, {}},
        {{{Family::DispersionEnergyOverK, {7}}}, 635.4895968646832, 10., 50., 950., {}},
    };
    const std::array temperatures = {273.15 + 60.02, 273.15 + 59.98, 273.15 + 79.99};
    const std::array observed_pressures = {19880., 19580., 47190.};
    for (std::size_t i = 0; i < temperatures.size(); ++i) {
        auto problem = pure_saturation(temperatures[i], observed_pressures[i]);
        problem.feed.values = {0., 1., 0., 0., 0., 0., 0., 0.};
        for (auto& phase : problem.phases) {
            phase.support = {1};
            phase.support_declared = true;
        }
        auto row = solved(ObservationKind::EquilibriumPressure, problem);
        row.observed = observed_pressures[i];
        request.observations.push_back(row);
    }
    const auto width = request.coordinates.coordinates.size();
    const auto evaluated = evaluated_columns(request, std::vector<double>(width, 0.));
    CHECK(evaluated.residual_optimizer_gradient.size() == temperatures.size() * width);
    bool exact = true;
    for (std::size_t i = 0; i < temperatures.size(); ++i) {
        CHECK(evaluated.statuses.at(i) == ObservationStatus::Available);
        for (std::size_t k = 0; k < width; ++k) {
            const double slope = evaluated.residual_optimizer_gradient[i * width + k];
            exact &= slope == 0.;
            recording_support::observe("Kim/pure-water/" + std::to_string(i) + "/" + std::to_string(k), "log-pressure slope",
                "1", slope, 0., "parameter regression", fields({field("T_K", temperatures[i]),
                    field("P_Pa", evaluated.predictions.at(i)), field("coordinate", k), field("exact_zero", slope == 0.)}));
        }
    }

    // At the liquid endpoint, MDEA insertion differentiates across the same site blocks.
    const auto amounts = request.observations[1].problem.feed.values;
    const HelmholtzModel model(request.parameters, {{Family::SigmaConstant, {6}}});
    const double temperature = temperatures[1];
    const double volume = 1. / model.state_at_pressure(temperature, evaluated.predictions.at(1), amounts,
        DensityBranch::Liquid).molar_density;
    const auto boundary = model.evaluate(temperature, volume, amounts, DerivativeOrder::Hessian);
    const std::size_t mdea = 6, insertion = 2 + mdea, width_with_parameter = boundary.gradient.size();
    const double mu = boundary.gradient.at(insertion);
    const auto perturbed = [&](double n) {
        auto inserted = amounts;
        inserted[mdea] = n;
        return model.evaluate(temperature, volume, inserted, DerivativeOrder::Value).value;
    };
    constexpr double h = 1e-5;
    CHECK(agree("Kim/MDEA-insertion", "residual chemical potential divided by RT", "1",
        mu, (4 * (perturbed(h) - boundary.value) - (perturbed(2 * h) - boundary.value)) / (2 * h), 0., FIRST_ORDER,
        fields({field("T_K", temperature), field("V_m3", volume), field("h_mol", h)})));
    const auto fd = ladder([&](double step) {
        auto parameters = request.parameters;
        parameters.sigma_constant[mdea] += step;
        return HelmholtzModel(parameters).evaluate(temperature, volume, amounts, DerivativeOrder::Gradient).gradient.at(insertion);
    }, 1e-3);
    CHECK(ladder_agrees("Kim/MDEA-insertion/diameter", "d(mu_res/RT)/d(sigma_MDEA)",  // reciprocal angstrom
        boundary.hessian.at(insertion * width_with_parameter + width_with_parameter - 1), fd,
        fields({field("T_K", temperature), field("V_m3", volume), field("h_angstrom", 1e-3)}), true));
    CHECK(truth("Kim/pure-water/exact-zero", exact, fields({field("rows", temperatures.size()), field("coordinates", width)})));
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
        CHECK(truth("refusal/projection", refused.status == ObservationStatus::CentralSolveUnavailable
            && !refused.value && !refused.residual && refused.convergence && refused.refusal_convergence
            && core::describe(*refused.convergence) == core::describe(*failure)
            && core::describe(*refused.refusal_convergence) == core::describe(*failure)
            && refused.refusal_convergence->evaluations == failure->evaluations,
            fields({field("status", observation_status_name(refused.status)), field("cause", refused.message),
                field("refusal_record", refused.refusal_convergence ? core::describe(*refused.refusal_convergence) : "absent"),
                field("mass_action_record", core::describe(*failure))})));
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
    auto solver_options = request.controls.solver_options;
    solver_options.nested_solve = true;
    const auto active = request.coordinates.active_parameters();
    const HelmholtzModel model(parameters, active, {}, solver_options);
    auto first = solved(ObservationKind::EquilibriumPressure, equilibrium_support::bubble_at(330.));
    const auto reference = predict_validated(model, first);
    CHECK(reference.value && reference.first_numerical_failure
        && reference.first_numerical_failure->stop == core::Stop::converged_to_floor);
    CHECK(agree("bubble/recovered-pressure", "pressure", "Pa", *reference.value,
        6023.3839326137504, 0., 1e-10, fields({field("T_K", 330.), field("feed", first.problem.feed.values)})));
    first.observed = *reference.value;
    request.weights = {1., 1.};

    {
        const HelmholtzModel retry_model(properties_support::associating_binary(), {}, {}, solver_options);
        core::Capture capture(true);
        const auto state = retry_model.state(330., 2918., {0.5, 0.5});
        const auto& spans = capture.trace().spans;
        const auto site_solve = std::find_if(spans.begin(), spans.end(),
            [](const auto& span) { return std::string_view(span.name) == "association_newton"; });
        CHECK(site_solve != spans.end() && site_solve->convergence
            && site_solve->convergence->stop == core::Stop::converged && site_solve->convergence->iterations == 57);
        CHECK(state.first_numerical_failure && state.first_numerical_failure->stop == core::Stop::converged_to_floor
            && state.first_numerical_failure->iterations == 50);
    }

    auto partial_pivot_bubble = first;
    partial_pivot_bubble.problem.feed.values = {0.25, 0.75};
    CHECK(predict_validated(model, partial_pivot_bubble).value);

    auto refused_bubble = first;
    refused_bubble.problem.feed.values = {0.2, 0.8};
    // The default pressure seed now recovers. A 10 kPa seed tests the terminal path.
    refused_bubble.problem.pressure = FreeValue::variable(1e4);
    const auto refused_reference = predict_validated(model, refused_bubble);
    CHECK(!refused_reference.value && refused_reference.first_numerical_failure
        && refused_reference.first_numerical_failure->stop == core::Stop::converged_to_floor
        && refused_reference.first_numerical_failure->criterion
        && std::string_view(refused_reference.first_numerical_failure->criterion->name) == policy::ASSOCIATION_RESIDUAL.name);

    auto cell = equilibrium_support::flash(2);
    cell.temperature = FreeValue::fixed(329.);
    // The fit-order checks need a central prediction without a recovered site failure.
    cell.pressure = FreeValue::fixed(9900.);
    const auto coalesced = equilibrium_support::bubble_at(700.);
    const auto failed_solve = solve_equilibrium(model, coalesced);
    CHECK(!failed_solve.success && !failed_solve.terminal_numerical_failure
        && failed_solve.convergence.stop == core::Stop::infeasible
        && failed_solve.convergence.cause == "phases_coalesced");

    struct Refusal {
        const char* stage;
        Observation row;
        ObservationStatus status;
        std::string message;
        std::optional<double> value;
        std::string record;
        bool initial_trial = true;
        double weight = 1.;
        bool diagnostic_after = false;
    };
    const auto describe = [](const std::optional<core::Convergence>& record) {
        return record ? core::describe(*record) : "absent";
    };
    auto stability = solved(ObservationKind::Stability, cell);
    stability.form = ResidualForm::Difference;
    stability.observed = 0.;
    stability.pressure_reference = coalesced;
    const auto pressure = predict_validated(model, solved(ObservationKind::EquilibriumPressure, coalesced));
    CHECK(!pressure.value && pressure.first_numerical_failure && pressure.first_numerical_failure->criterion);
    auto stability_request = request;
    stability_request.controls.workers = 1;
    stability_request.observations = {stability};
    stability_request.weights = {1.};
    const auto stability_fit = fit(stability_request);
    CHECK(stability_fit.failures.size() == 1);
    const auto& stability_failure = stability_fit.failures.front();
    CHECK(truth("fit-failure/pressure-reference", stability_failure.observation == 0
        && stability_failure.status == ObservationStatus::CentralSolveUnavailable
        && stability_failure.message == "pressure reference: " + pressure.message
        && describe(stability_failure.first_numerical_failure) == describe(pressure.first_numerical_failure)
        && describe(stability_fit.first_numerical_failure) == describe(pressure.first_numerical_failure)
        && describe(stability_failure.convergence) == core::describe(failed_solve.convergence),
        fields({field("first", describe(stability_failure.first_numerical_failure)),
            field("nested_first", describe(pressure.first_numerical_failure)), field("record", describe(stability_failure.convergence))})));

    const auto subnormal_scale = std::numeric_limits<double>::min() / 16.;
    auto large_mass = parameters;
    large_mass.molar_mass_kg_per_mol[0] = 1e308;
    const HelmholtzModel mass_model(large_mass, active, {}, solver_options);
    const auto terminal_parameters = request.coordinates.apply(parameters, {-1000.});
    const HelmholtzModel terminal_model(terminal_parameters, active, {}, solver_options);
    for (const auto& row : {solved(ObservationKind::PhaseDensity, cell), site_row(cell, 0, 0)}) {
        const bool site = row.site.has_value();
        const auto prediction = predict_validated(model, row);
        CHECK(prediction.value && prediction.solved && prediction.solved->success && prediction.convergence);
        const auto action = differentiate(model, row, prediction);
        CHECK(action.status == ObservationStatus::Available && action.value.size() == 1
            && std::isfinite(action.value.front()) && action.value.front() != 0.);

        auto terminal = row;
        terminal.observed = site ? 1e-250 : 1e250;
        const auto terminal_prediction = predict(terminal_model, terminal);
        CHECK(!terminal_prediction.value && terminal_prediction.convergence && terminal_prediction.convergence->criterion
            && terminal_prediction.convergence->stop == core::Stop::nonfinite
            && std::string_view(terminal_prediction.convergence->criterion->name) == policy::ASSOCIATION_RESIDUAL.name);
        auto driver = row;
        driver.problem = coalesced;
        const auto driver_prediction = predict_validated(model, driver);
        CHECK(!driver_prediction.value && driver_prediction.status == ObservationStatus::TopologyEvent);

        auto cost = row;
        cost.form = ResidualForm::Difference;
        cost.scale = 1e-160;
        cost.observed = 1.;
        const auto cost_prediction = predict_validated(model, cost);
        CHECK(cost_prediction.residual && std::isfinite(*cost_prediction.residual)
            && !std::isfinite(*cost_prediction.residual * *cost_prediction.residual));
        auto residual = cost;
        residual.scale = subnormal_scale;
        const auto residual_prediction = predict_validated(model, residual);
        CHECK(residual_prediction.value == prediction.value && !residual_prediction.residual);
        auto jacobian = residual;
        jacobian.observed = *prediction.value;
        CHECK(predict_validated(model, jacobian).residual == 0.);

        const double large_weight = std::numeric_limits<double>::max();
        CHECK(!std::isfinite(std::sqrt(large_weight) * *cost_prediction.residual));
        std::vector<Refusal> cases{
            {"driver", driver, ObservationStatus::TopologyEvent, driver_prediction.message, {}, core::describe(failed_solve.convergence)},
            {"residual", residual, ObservationStatus::Nonfinite, "the residual is not finite", prediction.value, "absent"},
            {"weighted-residual", cost, ObservationStatus::Nonfinite, "weighted residual not finite", prediction.value, "absent", true, large_weight},
            {"weighted-cost", cost, ObservationStatus::Nonfinite, "accumulated weighted cost not finite", prediction.value, "absent"},
            {"weighted-cost-prefix", cost, ObservationStatus::Nonfinite, "accumulated weighted cost not finite", prediction.value, "absent", true, 1., true},
            {"optimizer-jacobian", jacobian, ObservationStatus::Nonfinite, "optimizer Jacobian not finite", prediction.value, "absent"},
            {"terminal", terminal, terminal_prediction.status, terminal_prediction.message, {}, describe(terminal_prediction.convergence), false},
        };
        if (!site) {
            auto mass = row;
            mass.mass_basis = true;
            const auto projection = predict_validated(mass_model, mass);
            CHECK(projection.solved && projection.solved->success && !projection.value && !projection.residual);
            cases.push_back({"projection", mass, ObservationStatus::Nonfinite, "the predicted value is not finite", {}, "absent"});
        }

        const auto earlier = describe(reference.first_numerical_failure);
        for (const auto& item : cases) for (const int workers : {1, 2}) {
            auto trial = request;
            if (item.row.mass_basis) trial.parameters = large_mass;
            trial.controls.workers = workers;
            trial.observations = {first, item.row};
            trial.weights[1] = item.weight;
            std::optional<core::Convergence> fallback;
            if (item.diagnostic_after) {
                CHECK(!cost_prediction.first_numerical_failure && !action.first_numerical_failure);
                std::swap(trial.observations[0], trial.observations[1]);
                std::swap(trial.weights[0], trial.weights[1]);
                core::retain_failure(fallback, item.message, "fit evaluation");
            }
            const auto expected_first = item.diagnostic_after ? describe(fallback) : earlier;
            const auto result = fit(trial);
            CHECK(!result.failures.empty());
            const auto& failure = result.failures.front();
            CHECK(truth(std::string("fit-failure/") + item.stage, failure.observation == (item.diagnostic_after ? 0 : 1)
                && failure.status == item.status && failure.message == item.message && failure.central_value == item.value
                && describe(failure.convergence) == item.record && describe(failure.first_numerical_failure) == expected_first
                && describe(result.first_numerical_failure) == expected_first
                && (item.initial_trial ? failure.jacobian && result.status == "initial_evaluation_failed"
                    && result.residuals.empty() && result.optimizer_jacobian.empty()
                    : !failure.jacobian && failure.optimizer == std::vector<double>{-1000.}
                    && result.optimizer == std::vector<double>{0.} && result.predictions.at(1) == *prediction.value),
                fields({field("workers", workers), field("site_row", site), field("row", failure.observation),
                    field("cause", failure.message), field("first", describe(result.first_numerical_failure)),
                    field("record", describe(failure.convergence)), field("expected_record", item.record),
                    (failure.central_value ? field("central", *failure.central_value) : field("central", "absent")),
                    field("weighted_cost_residual", *cost_prediction.residual), field("first_action", action.value.front()),
                    field("optimizer", failure.optimizer)})));
        }
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
        {{"native.regression.kim-pure-water", "kim-pure-water"}, kim_pure_water},
        {{"native.regression.association-site-temperature", "temperature"}, bubble_temperature},
        {{"native.regression.association-site-mapping", "mapping"}, pathway_mapping},
        {{"native.regression.association-site-reaction", "reaction"}, reacting_actions},
        {{"native.regression.association-site-declaration", "declaration"}, declaration_refusals},
        {{"native.regression.association-site-refusals", "refusals"}, boundary_refusals},
        {{"native.regression.association-site-fit-failure", "fit-failure"}, fit_failure_records},
    });
}
