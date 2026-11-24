// R1 F1-F10, #86 affine columns and density actions, and R2 J1 prediction gradients.
#include "regression_checks_shared.hpp"

namespace regression_checks {

void check_f1() {
    const HelmholtzModel model(properties_support::ethane());
    for (const double observed : {0.0, -1.0}) {
        auto row = density(200.0, 217800.6949469248, {1.0}, DensityBranch::Liquid);
        row.observed = observed;
        CHECK(recording_support::refusal([&] { validate(row, model); }) == "observation positive LogRatio observation");
        FitRequest request; request.parameters = model.parameters(); request.observations = {row};
        FitResult failures;
        CHECK(recording_support::refusal([&] { ObservationEvaluator evaluator(request, failures, {0}); }) == "observation positive LogRatio observation");
    }
    const double mass = 0.03007;
    for (const auto& row : NIST) {
        const std::string point = "T" + std::to_string(static_cast<int>(row.T));
        const auto problem = pure_saturation(row.T, row.p_source);
        auto pressure = solved(ObservationKind::EquilibriumPressure, problem);
        pressure.observed = row.p_source;
        const auto p = predict(model, pressure);
        CHECK(p.status == ObservationStatus::Available && p.unit == "Pa" && p.residual.has_value());
        const auto inputs = fields({field("T_K", row.T), field("source_p_sat_Pa", row.p_source),
                                    field("source_rho_liq_kg_m3", row.rho_source), field("branch", "L free / V pinned 0")});
        agree("r1-f1/" + point + "/p_sat", "pressure", "Pa", *p.value, row.p_model, 0.0, 1.0e-9, inputs);
        agree("r1-f1/" + point + "/p_deviation", "relative_deviation", "1", (*p.value - row.p_source) / row.p_source,
              (row.p_model - row.p_source) / row.p_source, 1.0e-9, 0.0, inputs);
        auto rho_mass = solved(ObservationKind::PhaseDensity, problem, 0, {}, true);
        rho_mass.observed = row.rho_source;
        auto rho_molar = solved(ObservationKind::PhaseDensity, problem, 0, {}, false);
        const auto kg = predict(model, rho_mass), mol = predict(model, rho_molar);
        CHECK(kg.status == ObservationStatus::Available && kg.unit == "kg/m^3" && mol.unit == "mol/m^3");
        agree("r1-f1/" + point + "/rho_liq", "saturated_liquid_density", "kg/m^3", *kg.value, row.rho_model,
              0.0, 1.0e-9, inputs);
        agree("r1-f1/" + point + "/mass_chain", "saturated_liquid_density", "kg/m^3", *kg.value, *mol.value * mass,
              0.0, 1.0e-12, inputs);
    }
    const auto value = predict(model, density(200.0, 217800.6949469248, {1.0}, DensityBranch::Liquid));
    CHECK(value.status == ObservationStatus::Available);
    agree("r1-f1/T200/o1a_density", "molar_density", "mol/m^3", *value.value, 17541.25129083691, 0.0, 1.0e-9,
          fields({field("T_K", 200.0), field("P_Pa", 217800.6949469248), field("branch", "liquid")}));
    FitRequest request; request.parameters = model.parameters();
    request.observations = {density(200.0, 217800.6949469248, {1.0}, DensityBranch::Liquid)};
    FitResult failures;
    ObservationEvaluator evaluator(request, failures, {0});
    for (const bool columns : {false, true}) {
        ObservationEvaluation evaluated;
        CHECK(evaluator.evaluate({}, columns, evaluated));
        CHECK(evaluated.predictions.size() == 1 && evaluated.statuses == std::vector{ObservationStatus::Available});
        agree(columns ? "r1-f1/empty_coordinates/derivative" : "r1-f1/empty_coordinates/prediction",
              "molar_density", "mol/m^3", evaluated.predictions[0], 17541.251290836917, 0.0, 1.0e-9,
              fields({field("T_K", 200.0), field("P_Pa", 217800.6949469248)}));
        CHECK(evaluated.gradient.empty() && evaluated.optimizer_gradient.empty()
              && evaluated.residual_optimizer_gradient.empty() && evaluated.chart.map.empty());
    }
    const auto manufactured = observations_at(request, {}, {0.0});
    CHECK(manufactured.size() == 1);
    agree("r1-f1/empty_coordinates/manufactured", "molar_density", "mol/m^3", manufactured[0].observed,
          17541.251290836917, 0.0, 1.0e-9, "{}");
    // A neutral species always carries a positive molar mass by Properties admission; the missing-mass
    // path is reachable only for an ion without a recorded mass.
    auto ionic = properties_support::aqueous_electrolyte();
    ionic.molar_mass_kg_per_mol[1] = 0.0;
    const HelmholtzModel ions(ionic);
    const auto missing = predict(ions, density(298.15, 1.0e5, {0.96, 0.02, 0.02}, DensityBranch::Liquid, true));
    truth("r1-f1/missing_molar_mass",
          missing.status == ObservationStatus::MissingMolarMass && !missing.value.has_value(),
          fields({field("status", observation_status_name(missing.status))}));
}

// ---- R1 F2: affine coordinates, bounds and admission on all seven families.
CoordinateSet seven_families() {
    CoordinateSet set;
    set.coordinates = {
        {{{ParameterFamily::DispersionEnergyOverK, {0}}}, 221.0, 10.0, 150.0, 260.0},
        {{{ParameterFamily::SegmentCount, {0}}}, 1.9489, 0.1, 1.5, 6.0},
        {{{ParameterFamily::SigmaConstant, {0}}}, 3.7209, 0.1, 2.0, 4.2},
        {{{ParameterFamily::Kij, {1, 2}}}, 0.0, 1.0, -1.0, 1.0},
        {{{ParameterFamily::Lij, {1, 2}}}, 0.0, 1.0, -1.0, 1.0},
        {{{ParameterFamily::AssociationEnergyOverK, {1, 4}}}, 0.0, 1.0, 0.0, 5000.0},
        {{{ParameterFamily::AssociationVolume, {1, 4}}}, 0.0, 1.0, 0.0, 1.0},
    };
    return set;
}

void check_affine_start_signed();

void check_f2() {
    check_affine_start_signed();
    const auto parameters = properties_support::associating_binary();
    const auto set = seven_families();
    validate(set, parameters);
    const auto start = set.start(parameters);
    const auto physical = set.physical(start);
    const auto active = set.active_parameters();
    CHECK(physical.size() == 7 && active.size() == 7);
    for (std::size_t r = 0; r < 7; ++r) {
        const double truth_value = detail::family_values(parameters, active[r].family).first[active[r].indices.front()];
        const double ulp = 2.0 * EPS * std::max(std::abs(truth_value), std::abs(set.coordinates[r].origin));
        agree("r1-f2/roundtrip/" + std::to_string(r), "theta", "native", physical[r], truth_value, ulp, 0.0,
              fields({field("origin", set.coordinates[r].origin), field("scale", set.coordinates[r].scale),
                      field("u", start[r])}));
    }
    const auto applied = set.apply(parameters, start);
    bool identical = true;
    for (const auto& item : active) {
        const auto& lhs = detail::family_values(applied, item.family);
        const auto& rhs = detail::family_values(parameters, item.family);
        identical = identical && lhs.second == rhs.second && std::equal(lhs.first, lhs.first + lhs.second, rhs.first);
    }
    truth("r1-f2/apply_identity", identical, "");
    // Bounds map through the same affine transform and are inclusive.
    auto lower = start;
    lower[0] = (150.0 - 221.0) / 10.0;
    truth("r1-f2/bound_inclusive", set.within_bounds(lower), field("u0", lower[0]));
    lower[0] -= 1.0e-9;
    truth("r1-f2/bound_exclusive", !set.within_bounds(lower), field("u0", lower[0]));
    auto supplied_chart = set;
    supplied_chart.coordinates[0].lower = 230.0;
    auto supplied_point = start; supplied_point[0] = 0.9;
    const auto supplied = evaluate_coordinates(parameters, supplied_chart, supplied_point, {});
    CHECK(supplied.physical[0] == 230.0);
    CHECK(recording_support::refusal([&] { supplied_chart.start(parameters); }) == "coordinate start outside bounds");
    struct Rejection final { const char* phrase; std::function<void(CoordinateSet&)> mutate; };
    const std::vector<Rejection> rejections = {
        {"start outside bounds", [](CoordinateSet& s) { s.coordinates[0].lower = 230.0; }},
        {"scale", [](CoordinateSet& s) { s.coordinates[1].scale = 0.0; }},
        {"bounds", [](CoordinateSet& s) { s.coordinates[2].lower = 5.0; }},
        {"slot index", [](CoordinateSet& s) { s.coordinates[1].targets[0].indices = {99}; }},
        {"targets", [](CoordinateSet& s) { s.coordinates[1].targets.clear(); }},
        {"set is empty", [](CoordinateSet& s) { s.coordinates.clear(); }},
        {"diagonal parameter coordinate", [](CoordinateSet& s) { s.coordinates[3].targets[0].indices = {0}; }},
        {"reciprocal parameter coordinate", [](CoordinateSet& s) { s.coordinates[3].targets[0].indices = {1}; }},
    };
    for (const auto& item : rejections) {
        auto broken = set;
        item.mutate(broken);
        const auto message = recording_support::refusal([&] { validate(broken, parameters); });
        truth(std::string("r1-f2/reject/") + item.phrase, message.find(item.phrase) != std::string::npos,
              field("message", message));
    }
}

// ---- R1 F3 and F4: same-family sharing, reciprocal ownership, fixing and the active-list order.
void check_f3_f4() {
    auto parameters = properties_support::associating_binary();
    parameters.sigma_constant = {2.5, 2.5};
    CoordinateSet tied;
    tied.coordinates = {{{{ParameterFamily::SigmaConstant, {0}}, {ParameterFamily::SigmaConstant, {1}}}, 2.5, 0.1, 2.0, 4.0},
                        {{{ParameterFamily::Kij, {1, 2}}}, 0.0, 1.0, -1.0, 1.0}};
    validate(tied, parameters);
    const std::vector<std::size_t> tie = tied.tie();
    truth("r1-f3/tie_map",
          tie == std::vector<std::size_t>{0, 0, 1}, "");
    const std::vector<double> u = {0.3, -0.02};
    const auto applied = tied.apply(parameters, u);
    auto hand = parameters;
    hand.sigma_constant = {2.5 + 0.1 * 0.3, 2.5 + 0.1 * 0.3};
    hand.k_ij[1] = hand.k_ij[2] = -0.02;
    truth("r1-f3/shared_slots",
          applied.sigma_constant == hand.sigma_constant && applied.k_ij == hand.k_ij, field("u", u));
    const auto lhs = HelmholtzModel(applied).state(300.0, 5000.0, {0.4, 0.6});
    const auto rhs = HelmholtzModel(hand).state(300.0, 5000.0, {0.4, 0.6});
    agree("r1-f3/state_identity", "residual_helmholtz", "1", lhs.residual_helmholtz, rhs.residual_helmholtz, 0.0, 0.0,
          fields({field("T", 300.0), field("rho", 5000.0)}));
    auto half = parameters;
    half.k_ij[1] = -0.02;
    const auto message = recording_support::refusal([&] { (void)HelmholtzModel(half); });
    truth("r1-f3/one_slot_rejected",
          message.find("k_ij symmetry") != std::string::npos, field("message", message));
    auto disagree = parameters;
    disagree.sigma_constant[1] = 2.6;
    const auto start_message = recording_support::refusal([&] { (void)tied.start(disagree); });
    truth("r1-f3/tied_start_disagrees",
          start_message.find("tied values disagree") != std::string::npos, field("message", start_message));
    // F4: the model's active list is the coordinate set's list in order; a fixed family projects to zero.
    const HelmholtzModel model(applied, tied.active_parameters());
    const auto& active = model.active_parameters();
    bool ordered = active.size() == 3;
    for (std::size_t r = 0; ordered && r < 3; ++r)
        ordered = active[r].family == tied.active_parameters()[r].family && active[r].indices == tied.active_parameters()[r].indices;
    truth("r1-f4/active_order", ordered, "");
    const auto fixed = model.pressure_actions(300.0, 1.0e5, {0.4, 0.6}, DensityBranch::Liquid, {},
        {{0.0, 0.0, {0.0, 0.0}, {{ParameterFamily::DispersionEnergyOverK, {0}, 1.0}}}});
    CHECK(fixed.density_status == PropertyActionStatus::Available && fixed.density_action.size() == 1);
    agree("r1-f4/fixed_projects_to_zero", "density_action", "mol/m^3 per K", fixed.density_action[0], 0.0, 1.0e-12, 0.0,
          field("family", "DispersionEnergyOverK fixed"));
    auto outside = tied;
    outside.coordinates[0].targets[1].indices = {5};
    const auto slot_message = recording_support::refusal([&] { validate(outside, parameters); });
    truth("r1-f4/out_of_range_rejected",
          slot_message.find("slot index") != std::string::npos, field("message", slot_message));
}

// ---- R1 F5 and F6: pressure/temperature specifications and the pure saturation identity.
void check_f5_f6() {
    const HelmholtzModel model(properties_support::neutral_binary());
    const auto bubble = equilibrium_support::bubble_at(180.0);
    const auto direct = equilibrium::solve_equilibrium(model, bubble);
    const auto pressure = predict(model, solved(ObservationKind::EquilibriumPressure, bubble));
    CHECK(pressure.status == ObservationStatus::Available && direct.success);
    agree("r1-f5/bubble_pressure", "pressure", "Pa", *pressure.value, *direct.pressure, 0.0, 0.0,
          field("T_K", 180.0));
    // The pinned-zero free-T solve is solved cold from a guess above the answer (a fit's first trial
    // is always cold); a warm start from the O2 result would return the seed bit for bit and prove nothing.
    auto free_t = bubble;
    free_t.temperature = FreeValue::variable(190.0);
    free_t.pressure = FreeValue::fixed(*pressure.value);
    const auto temperature = predict(model, solved(ObservationKind::EquilibriumTemperature, free_t));
    truth("r1-f5/bubble_temperature_status",
          temperature.status == ObservationStatus::Available, field("status", observation_status_name(temperature.status)));
    agree("r1-f5/bubble_temperature", "temperature", "K", *temperature.value, 180.0, 1.0e-6, 0.0,
          fields({field("P_Pa", *pressure.value), field("T_guess", 190.0)}));
    // The cold start from below reaches the same bubble state; no refusal is admitted.
    free_t.temperature = FreeValue::variable(170.0);
    const auto below = predict(model, solved(ObservationKind::EquilibriumTemperature, free_t));
    truth("r1-f5/bubble_temperature_from_below_status",
          below.status == ObservationStatus::Available,
          fields({field("status", observation_status_name(below.status)), field("message", below.message)}));
    agree("r1-f5/bubble_temperature_from_below", "temperature", "K", *below.value, 180.0, 1e-6, 0.0,
          field("T_guess", 170.0));
    auto dew = bubble;
    dew.phases[0].amount = AmountSpec::pinned(0.0);
    dew.phases[1].amount = AmountSpec::free();
    const auto dew_pressure = predict(model, solved(ObservationKind::EquilibriumPressure, dew));
    CHECK(dew_pressure.status == ObservationStatus::Available);
    truth("r1-f5/dew_distinct",
          std::abs(*dew_pressure.value / *pressure.value - 1.0) > 1.0e-3,
          fields({field("P_bubble", *pressure.value), field("P_dew", *dew_pressure.value)}));
    // F6
    const HelmholtzModel ethane(properties_support::ethane());
    const auto up = predict(ethane, solved(ObservationKind::EquilibriumPressure, pure_saturation(200.0, 2.0e5)));
    const auto down = predict(ethane, solved(ObservationKind::EquilibriumPressure, pure_saturation(200.0, 2.0e5, true)));
    CHECK(up.status == ObservationStatus::Available && down.status == ObservationStatus::Available);
    agree("r1-f6/pure_bubble_dew", "pressure", "Pa", *up.value, *down.value, 0.0, 1.0e-8, field("T_K", 200.0));
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto x = predict(ethane, solved(ObservationKind::PhaseMoleFraction, pure_saturation(200.0, 2.0e5), phase, {0}));
        agree("r1-f6/pure_fraction/" + std::to_string(phase), "mole_fraction", "1", *x.value, 1.0, 0.0, 0.0,
              field("phase", phase));
    }
}

// ---- R1 F7 and F8: phase fractions on the E4-B tie line; aggregate compositions on May row 001.
void check_f7_f8() {
    const HelmholtzModel model(properties_support::neutral_binary());
    const auto tie = equilibrium_support::interior(model);
    std::vector<double> beta(2), fraction(4);
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto b = predict(model, solved(ObservationKind::PhaseFraction, tie.problem, phase), tie.start);
        CHECK(b.status == ObservationStatus::Available);
        beta[phase] = *b.value;
        for (std::size_t i = 0; i < 2; ++i)
            fraction[phase * 2 + i] = *predict(model, solved(ObservationKind::PhaseMoleFraction, tie.problem, phase, {i}), tie.start).value;
    }
    const auto inputs = fields({field("T_K", 180.0), field("P_Pa", 0.5 * tie.bubble_pressure), field("z", tie.problem.feed.values)});
    agree("r1-f7/fraction_sum", "sum_beta", "1", beta[0] + beta[1], 1.0, 1.0e-14, 0.0, inputs);
    for (std::size_t i = 0; i < 2; ++i)
        agree("r1-f7/lever/" + std::to_string(i), "feed_fraction", "1", beta[0] * fraction[i] + beta[1] * fraction[2 + i],
              tie.problem.feed.values[i], 1.0e-12, 0.0, inputs);
    auto scaled = tie.problem;
    for (auto& value : scaled.feed.values) value *= 3.0;
    const auto b3 = predict(model, solved(ObservationKind::PhaseFraction, scaled, 0), tie.start);
    CHECK(b3.status == ObservationStatus::Available);
    agree("r1-f7/fraction_homogeneous", "beta_L", "1", *b3.value, beta[0], 0.0, 1.0e-12, inputs);
    agree("r1-f7/amount_scales", "n_L", "mol", b3.solved->amounts[0], 3.0 * tie.start.amounts[0], 0.0, 1.0e-12, inputs);
    // F8: May 2015 row 001 (ref-may/may2015-ch4-c2h6-001/solve, E6 retained compositions).
    const std::vector<double> retained = {0.3025223259394309, 0.697477674060569, 0.6703563353246006, 0.3296436646753994};
    auto may = equilibrium_support::flash(2);
    may.phases = {equilibrium_support::phase("L", PhaseKind::Liquid), equilibrium_support::phase("V", PhaseKind::Vapor)};
    may.phases[0].composition_guess = {0.3099, 0.6901};
    may.phases[1].composition_guess = {0.6664, 0.3336};
    may.temperature = FreeValue::fixed(243.58);
    may.pressure = FreeValue::fixed(3.949e6);
    may.feed.values = {0.5 * (retained[0] + retained[2]), 0.5 * (retained[1] + retained[3])};
    may.solver_options.values[epcsaft::policy::row] = 1.0e-10;
    const auto may_inputs = fields({field("T_K", 243.58), field("P_Pa", 3.949e6), field("z", may.feed.values)});
    const auto both = predict(model, solved(ObservationKind::PhaseMoleFraction, may, 0, {0, 1}));
    CHECK(both.status == ObservationStatus::Available);
    agree("r1-f8/aggregate_sum", "x_CH4+x_C2H6", "1", *both.value, 1.0, 1.0e-15, 0.0, may_inputs);
    const auto x0 = predict(model, solved(ObservationKind::PhaseMoleFraction, may, 0, {0}), both.solved);
    // The start re-solves through the driver; two accepted solves of one problem agree to rounding.
    agree("r1-f8/single_species", "x_CH4", "1", *x0.value, both.solved->mole_fractions[0], 0.0, 1.0e-13, may_inputs);
    agree("r1-f8/retained_x", "x_CH4", "1", *x0.value, retained[0], 1.0e-9, 0.0, may_inputs);
    const auto y0 = predict(model, solved(ObservationKind::PhaseMoleFraction, may, 1, {0}), both.solved);
    agree("r1-f8/retained_y", "y_CH4", "1", *y0.value, retained[2], 1.0e-9, 0.0, may_inputs);
    // A mass basis is defined only for the two density kinds: a two-phase flash whose mole-basis
    // phase fraction differs from the physical mass fraction is refused by name, never reinterpreted.
    auto two = equilibrium_support::bubble_at(180.0);
    two.phases[1].amount = AmountSpec::free();
    two.pressure = FreeValue::fixed(1.0e6);
    const auto liquid_fraction = predict(model, solved(ObservationKind::PhaseFraction, two, 0));
    CHECK(liquid_fraction.status == ObservationStatus::Available && liquid_fraction.solved.has_value());
    const auto& state = *liquid_fraction.solved;
    const auto& masses = model.parameters().molar_mass_kg_per_mol;
    const double mass_liquid = state.amounts[0] * (state.mole_fractions[0] * masses[0] + state.mole_fractions[1] * masses[1]);
    const double mass_vapor = state.amounts[1] * (state.mole_fractions[2] * masses[0] + state.mole_fractions[3] * masses[1]);
    const auto two_inputs = fields({field("T_K", 180.0), field("P_Pa", 1.0e6), field("amounts", state.amounts)});
    agree("r1-f8/phase_fraction_mole_basis", "beta_L", "1", *liquid_fraction.value, 0.6738453209601837, 0.0, 1.0e-9, two_inputs);
    agree("r1-f8/phase_fraction_mass_control", "w_L", "1", mass_liquid / (mass_liquid + mass_vapor), 0.7580747775119783, 0.0, 1.0e-9,
          two_inputs);
    CHECK(std::abs(mass_liquid / (mass_liquid + mass_vapor) - *liquid_fraction.value) > 0.05);
    struct Rejection final { const char* phrase; Observation observation; };
    auto supported = may;
    supported.phases[0].support = {0};
    supported.phases[0].support_declared = true;
    const std::vector<Rejection> rejections = {
        {"species set is empty", solved(ObservationKind::PhaseMoleFraction, may, 0, {})},
        {"species repeated", solved(ObservationKind::PhaseMoleFraction, may, 0, {0, 0})},
        {"species outside the phase support", solved(ObservationKind::PhaseMoleFraction, supported, 0, {1})},
        {"phase index", solved(ObservationKind::PhaseFraction, may, 2)},
        {"mass basis is defined only for a density", solved(ObservationKind::PhaseFraction, two, 0, {}, true)},
        {"mass basis is defined only for a density", solved(ObservationKind::PhaseMoleFraction, two, 0, {0}, true)},
    };
    std::size_t rejection_index = 0;
    for (const auto& item : rejections) {
        const auto message = recording_support::refusal([&] { validate(item.observation, model); });
        truth("r1-f8/reject/" + std::to_string(rejection_index++) + "/" + item.phrase,
              message.find(item.phrase) != std::string::npos, field("message", message));
    }
}

// ---- R1 F9 and F10: preserved failure causes and the residual forms.
void check_f9_f10() {
    const HelmholtzModel model(properties_support::neutral_binary());
    auto unattainable = equilibrium_support::bubble_at(180.0);
    unattainable.solver_options.values[epcsaft::policy::row] = 1.0e-18;
    const auto floor = predict(model, solved(ObservationKind::EquilibriumPressure, unattainable));
    truth("r1-f9/unmet_tolerance",
          floor.status == ObservationStatus::CentralSolveUnavailable && !floor.value.has_value()
              && !floor.residual.has_value() && floor.requested_tolerance_met.has_value()
              && !*floor.requested_tolerance_met && floor.convergence.has_value()
              && std::isfinite(floor.convergence->residual),
          fields({field("status", observation_status_name(floor.status)), field("message", floor.message),
                  field("residual", floor.convergence ? floor.convergence->residual : -1.0)}));
    auto coincident = equilibrium_support::boundary();
    coincident.feed.values = {0.5, 0.5};
    coincident.temperature.value = 280.0;
    const auto topology = predict(model, solved(ObservationKind::EquilibriumPressure, coincident));
    // #565: observations keep the ordinary solve, without recovery: the convergence record is the default solve's.
    const auto ordinary = epcsaft::equilibrium::solve_compiled(epcsaft::equilibrium::compile_problem(model, coincident)).convergence;
    truth("r1-f9/topology_event",
          topology.status == ObservationStatus::TopologyEvent && !topology.value.has_value()
              && topology.message.find("phases_coincide") != std::string::npos && topology.convergence.has_value()
              && topology.convergence->iterations == ordinary.iterations && topology.convergence->evaluations == ordinary.evaluations
              && topology.convergence->stop == ordinary.stop && topology.convergence->residual == ordinary.residual,
          fields({field("status", observation_status_name(topology.status)), field("message", topology.message)}));
    const HelmholtzModel ethane(properties_support::ethane());
    // A declared branch selects among two stable roots; one stable root is returned as Unspecified with
    // its count, so the caller can see a branch that was not honoured. No root at all is the failure.
    const auto single = predict(ethane, density(200.0, 1.0e8, {1.0}, DensityBranch::Vapor));
    truth("r1-f9/single_root_branch",
          single.status == ObservationStatus::Available && single.branch == DensityBranch::Unspecified
              && single.stable_root_count == std::size_t{1},
          fields({field("status", observation_status_name(single.status)), field("rho", single.value.value_or(-1.0))}));
    const auto vapor = predict(ethane, density(200.0, 1.0e10, {1.0}, DensityBranch::Liquid));
    truth("r1-f9/density_unavailable",
          vapor.status == ObservationStatus::DensityUnavailable && !vapor.value.has_value() && !vapor.message.empty(),
          fields({field("status", observation_status_name(vapor.status)), field("message", vapor.message),
                  field("exception_shape", vapor.convergence.has_value() ? "ConvergenceFailure with record" : "runtime_error without record")}));
    // F10
    Observation difference = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    difference.form = ResidualForm::Difference;
    difference.observed = 1.0;
    difference.scale = 0.5;
    ObservationStatus status = ObservationStatus::Available;
    agree("r1-f10/difference", "residual", "1", *residual(difference, 3.0, status), 4.0, 0.0, 0.0,
          fields({field("y", 3.0), field("y_obs", 1.0), field("scale", 0.5)}));
    Observation ratio = difference;
    ratio.form = ResidualForm::LogRatio;
    agree("r1-f10/log_ratio", "residual", "1", *residual(ratio, std::exp(1.0), status), 1.0, 4.0 * EPS, 0.0,
          fields({field("y", std::exp(1.0)), field("y_obs", 1.0)}));
    ratio.observed = std::exp(1.0);
    agree("r1-f10/log_ratio_swapped", "residual", "1", *residual(ratio, 1.0, status), -1.0, 4.0 * EPS, 0.0,
          fields({field("y", 1.0), field("y_obs", std::exp(1.0))}));
    ratio.observed = 0.0;
    const auto rejected = residual(ratio, 1.0, status);
    truth("r1-f10/nonpositive_log",
          !rejected.has_value() && status == ObservationStatus::NonpositiveLogArgument,
          field("status", observation_status_name(status)));
    difference.scale = 0.0;
    const auto message = recording_support::refusal([&] { validate(difference, ethane); });
    truth("r1-f10/scale_rejected",
          message.find("scale") != std::string::npos, field("message", message));
    // The log ratio is finite across the dynamic range and keeps its accuracy near unity.
    ratio.observed = 1.0e-308;
    const double wide = *residual(ratio, 19936.48319859253, status);
    agree("r1-f10/log_ratio_wide_range", "residual", "1", wide, std::log(19936.48319859253) - std::log(1.0e-308), 0.0, 1.0e-12,
          fields({field("y", 19936.48319859253), field("y_obs", 1.0e-308)}));
    // Off an exactly divisible base the naive log(y / y_obs) loses five digits; the reference is the
    // 60-digit decimal ln y - ln y_obs of the rounded doubles.
    ratio.observed = 17541.251290836932;
    const double near = ratio.observed * (1.0 + 1.0e-12), naive = std::log(near / ratio.observed);
    agree("r1-f10/log_ratio_near_unity", "residual", "1", *residual(ratio, near, status), 1.0000617126413953e-12, 0.0, 1.0e-9,
          fields({field("y", near), field("y_obs", ratio.observed), field("naive_log", naive)}));
    CHECK(std::abs(naive - 1.0000617126413953e-12) > 1.0e-6 * 1.0000617126413953e-12);
    Observation overflow = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    overflow.form = ResidualForm::Difference;
    overflow.observed = -1.0e308;
    overflow.scale = 1.0e-308;
    status = ObservationStatus::Available;
    const auto overflowed = residual(overflow, 1.0, status);
    truth("r1-f10/difference_overflow",
          !overflowed.has_value() && status == ObservationStatus::Nonfinite, field("status", observation_status_name(status)));
    const auto kept = predict(ethane, overflow);
    truth("r1-f10/value_survives_nonfinite_residual",
          kept.status == ObservationStatus::Nonfinite && kept.value.has_value() && !kept.residual.has_value() && !kept.message.empty(),
          fields({field("status", observation_status_name(kept.status)), field("value", kept.value.value_or(-1.0))}));
    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = properties_support::ethane();
    request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 191.42, 10.0, 150.0, 260.0}};
    request.observations = {overflow, overflow, overflow};
    request.weights = {0.0, 1.0, 1.0};
    const auto refused = fit(request);
    truth("r1-f10/fit_nonfinite_typed",
          refused.status == "initial_evaluation_failed" && refused.failures.size() == 1 && refused.failures[0].observation == 1
              && refused.failures[0].status == ObservationStatus::Nonfinite && !refused.usable && std::isfinite(refused.initial_cost)
              && refused.statuses == std::vector<ObservationStatus>{ObservationStatus::NotEvaluated, ObservationStatus::Nonfinite, ObservationStatus::NotEvaluated},
          summary(refused));
    // A finite residual whose weighted form overflows is recorded by the initial-evaluation recorder with its row.
    Observation steep = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    steep.form = ResidualForm::Difference;
    steep.observed = *kept.value - 1.0e4;
    steep.scale = 1.0e-300;
    request.observations = {steep};
    request.weights = {1.0e40};
    const auto weighted = fit(request);
    truth("r1-f10/weighted_overflow_typed",
          weighted.status == "initial_evaluation_failed" && weighted.failures.size() == 1 && weighted.failures[0].observation == 0
              && weighted.failures[0].status == ObservationStatus::Nonfinite
              && weighted.initial_cost == 0.0,
          summary(weighted));
    // The coordinate chain can overflow with finite physical gradients and a finite objective;
    // refusal precedes publication of the initial successful-point cache and warm starts.
    steep.observed = *kept.value - 1.0;
    steep.scale = 1.0e-150;
    request.observations = {steep};
    request.weights.clear();
    request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 191.42, 1.0e200, 150.0, 260.0}};
    const auto steep_fit = fit(request);
    truth("r1-f10/jacobian_overflow_typed",
          steep_fit.status == "initial_evaluation_failed" && !steep_fit.usable && steep_fit.failures.size() == 1
              && steep_fit.failures[0].status == ObservationStatus::Nonfinite && steep_fit.failures[0].observation == 0
              && steep_fit.failures[0].jacobian,
          summary(steep_fit));
    // Neither the rejected Jacobian request at the committed point u nor the rejected trial at v is reused.
    {
        FitResult counts; ObservationEvaluator evaluator(request, counts, {0}); ObservationEvaluation committed, rejected, again;
        const std::vector<double> u = request.coordinates.start(request.parameters), v{1.0e-200};
        std::size_t step = 0;
        const auto at = [&](bool ok, std::size_t jacobians) {
            ++step;
            truth("r1-f9/no_reuse/step" + std::to_string(step),
                  ok && counts.residual_evaluations == 1 && counts.jacobian_evaluations == jacobians && counts.trial_failures == jacobians
                      && (jacobians == 0 || counts.failures.back().jacobian),
                  fields({field("residual_evaluations", counts.residual_evaluations), field("jacobian_evaluations", counts.jacobian_evaluations),
                          field("trial_failures", counts.trial_failures)}));
        };
        at(evaluator.evaluate(u, false, committed), 0);
        at(!evaluator.evaluate(u, true, rejected), 1);
        at(!evaluator.evaluate(v, true, rejected), 2);
        at(evaluator.evaluate(u, false, again) && again.predictions == committed.predictions && again.residuals == committed.residuals
               && again.weighted == committed.weighted && again.statuses == committed.statuses && again.cost == committed.cost
               && again.optimizer_gradient == committed.optimizer_gradient, 2);
        at(!evaluator.evaluate(u, true, rejected), 3);
        at(!evaluator.evaluate(v, true, rejected), 4);
    }
    // Three finite rows whose squares are representable but whose accumulated objective is not: the
    // failure is typed at the row where the running sum stops being finite, never an infinite cost.
    Observation steep_row = density(200.0, 1.0e5, {1.0}, DensityBranch::Liquid);
    steep_row.form = ResidualForm::Difference;
    steep_row.observed = *kept.value - 1.0;
    steep_row.scale = 1.0 / 1.2e154;
    request.observations = {steep_row, steep_row, steep_row};
    request.weights.clear();
    request.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 191.42, 10.0, 150.0, 260.0}};
    const auto accumulated = fit(request);
    truth("r1-f10/accumulated_cost_typed",
          accumulated.status == "initial_evaluation_failed" && !accumulated.usable && accumulated.failures.size() == 1
              && accumulated.failures[0].observation == 1 && accumulated.failures[0].status == ObservationStatus::Nonfinite
              && accumulated.initial_cost == 0.0
              && accumulated.statuses.size() == 3,
          fields({summary(accumulated), field("row_residual", *predict(ethane, steep_row).residual)}));
    // A refused start publishes no rows, so its zero costs cannot be read as evaluated residuals.
    for (const auto& [name, failed] : {std::pair{"refused", &refused}, {"weighted", &weighted}, {"accumulated", &accumulated}})
        truth(std::string("r1-f9/zero_cost_is_not_a_residual/") + name,
              failed->initial_cost == 0.0 && failed->final_cost == 0.0 && !failed->usable && failed->status == "initial_evaluation_failed"
                  && failed->termination == "not_run" && failed->predictions.empty() && failed->residuals.empty()
                  && failed->weighted_residuals.empty() && failed->covariance.outcome == FitCovarianceStatus::NotEvaluated && failed->acceptance == "not_evaluated",
              fields({summary(*failed), field("predictions", failed->predictions.size()), field("initial_cost", failed->initial_cost)}));
}

// ---- #86: first-order pressure-chart density actions against rebuilt-model closures.
Parameters coupled_lengths(double q) {
    auto p = properties_support::aqueous_electrolyte();
    p.sigma_constant[1] = q;
    p.packing_diameter_angstrom[1] = p.debye_huckel_diameter_angstrom[1] = 0.88 * q;
    p.born_diameter_angstrom[1] = 3.5;
    return p;
}
CoordinateSet coupled_coordinate() {
    return {{{{{ParameterFamily::PackingDiameter, {1}, 0.88}, {ParameterFamily::SigmaConstant, {1}, 1},
               {ParameterFamily::DebyeHuckelDiameter, {1}, 0.88}}, 3.0, 0.2, 2.4, 4.4}}};
}
HydrationHold coupled_hold() {
    Observation ob;
    ob.kind = ObservationKind::SolvationGibbsEnergy; ob.temperature = 298.15; ob.pressure = 1e5;
    ob.composition = {1, 0, 0}; ob.ion_indices = {1, 2}; ob.stoichiometry = {1, 1};
    ob.species = {1}; ob.branch = DensityBranch::Liquid;
    ob.observed = *predict(HelmholtzModel(coupled_lengths(3.2)), ob).value;
    return {ob, 1, 8, 1e-6, 1e-10, 50};
}
// The six activity and osmotic rows copy the hold's reference declaration.
std::vector<Observation> coupled_observations(const Observation& reference) {
    std::vector<Observation> rows;
    for (const auto kind : {ObservationKind::MeanIonicActivityCoefficient, ObservationKind::OsmoticCoefficient})
        for (const double m : {0.1, 1.0, 4.5}) {
            auto ob = reference;
            ob.kind = kind; ob.species.clear(); ob.molality = m;
            ob.form = ResidualForm::LogRatio; ob.observed = 1;
            rows.push_back(ob);
        }
    return rows;
}
void affine_difference(const std::string& id, double exact, const std::function<double(double)>& value,
                       double step, double scale) {
    // The step/2 ladder's coarse rung is this ladder's fine rung; only its finest rung is new.
    const auto a = ladder(value, step);
    const double half = step / 2;
    const Ladder b{a.fine, (value(0.5 * half) - value(-0.5 * half)) / half};
    const double reference = scale * (4 * b.fine - b.coarse) / 3;
    const double uncertainty = scale * std::abs((4 * b.fine - b.coarse - 4 * a.fine + a.coarse) / 3);
    const double tolerance = 1e-7 + 2e-5 * std::abs(exact);
    const auto inputs = fields({field("step", step), field("coarse", scale*a.coarse), field("middle", scale*a.fine),
        field("fine", scale*b.fine), field("uncertainty", uncertainty), field("tolerance", tolerance)});
    agree(id + "/refinement", "Richardson_uncertainty", "per optimizer unit", uncertainty, 0, tolerance/4, 0,
  inputs);
    agree(id, "optimizer_column", "per optimizer unit", exact, reference, tolerance, 0,
  inputs);
}
void check_affine_start_signed() {
    const auto p = coupled_lengths(3.2); const auto cs = coupled_coordinate();
    const auto v = cs.start(p); const auto mapped = cs.apply(p, v);
    for (const auto& target : cs.coordinates[0].targets) {
        const double old = detail::family_values(p,target.family).first[1], value = detail::family_values(mapped,target.family).first[1];
        agree("r1-f2/coupled-roundtrip/" + std::to_string(static_cast<int>(target.family)), "diameter", "angstrom",
            value, old, 8*EPS*std::max(std::abs(old),std::abs(value)), 0, field("v",v));
    }
    auto broken = p; broken.debye_huckel_diameter_angstrom[1] += 1e-8;
    const auto message = recording_support::refusal([&] { (void)cs.start(broken); });
    truth("r1-f2/coupled-disagreement",
        message.find("tied values disagree") != std::string::npos, field("message",message));
    FitRequest request; request.parameters = properties_support::neutral_binary();
    request.parameters.k_ij[1] = request.parameters.k_ij[2] = -0.088;
    request.coordinates = {{{{{ParameterFamily::Kij, {1,2}, -0.88}}, 0.1, 0.2, -1, 1}}};
    auto ob = density(200, 1e5, {0.4,0.6}, DensityBranch::Liquid); ob.form = ResidualForm::Difference;
    request.observations = {ob};
    const auto exact = evaluated_columns(request,{0}).optimizer_gradient[0];
    for (const double scale : {1.0, 3.0}) {
        request.observations[0].scale = scale;
        const auto columns = evaluated_columns(request,{0});
        const double residual = columns.residual_optimizer_gradient[0], expected = exact/scale;
        agree("r1-f2/residual-basis/"+std::to_string(scale), "density_residual_action", "per optimizer unit",
            residual,expected,8*EPS*std::max(std::abs(residual),std::abs(expected)),0,field("scale",scale));
    }
    const auto value = [&](double v) {
        auto physical = request.parameters; physical.k_ij[1] = physical.k_ij[2] = -0.88*(0.1+0.2*v);
        return *predict(HelmholtzModel(physical),ob).value;
    };
    affine_difference("r1-f2/signed-pair",exact,value,1e-3,1);
    truth("r1-f2/signed-signal",
        std::abs(exact)>1e-7+2e-5*std::abs(exact), field("exact",exact));
}
void check_coupled_columns() {
    FitRequest request; request.parameters = coupled_lengths(3.2); request.coordinates = coupled_coordinate();
    const auto hold = coupled_hold();
    request.observations = coupled_observations(hold.reference);
    for (const bool held : {false,true}) {
        request.holds = held ? std::vector<HydrationHold>{hold} : std::vector<HydrationHold>{};
        const auto exact = evaluated_columns(request,{1});
        for (std::size_t r = 0; r < 3; ++r)
            agree("affine-columns/map/"+std::to_string(held)+"/"+std::to_string(r), "D", "angstrom per v",
                exact.chart.map[r], r==1 ? 0.2 : 0.176, 1e-16,0,"");
        std::size_t root_row = 0;
        const auto rebuilt = [&](double shift) {
            auto physical = coupled_lengths(3.2+shift);
            if (held) {
                const CoordinateSet single{{{{{ParameterFamily::SigmaConstant,{1}}},0,1}}};
                const auto root = evaluate_coordinates(physical,single,{3.2+shift},{hold});
                agree("affine-columns/root/"+std::to_string(root_row++), "hydration_error", "J/mol", root.holds[0].energy_error,
                    0,1e-6,0,field("q",3.2+shift));
                physical = root.parameters;
            }
            return physical;
        };
        for (std::size_t row=0; row<request.observations.size(); ++row) {
            const auto value = [&](double shift) { return std::log(*predict(HelmholtzModel(rebuilt(shift)),request.observations[row]).value); };
            affine_difference("affine-columns/"+std::to_string(held)+"/"+std::to_string(row),
                exact.residual_optimizer_gradient[row],value,1e-3,0.2);
            const double residual = exact.residual_optimizer_gradient[row], expected = exact.optimizer_gradient[row]/exact.predictions[row];
            agree("affine-columns/residual-basis/"+std::to_string(held)+"/"+std::to_string(row), "log_ratio_action", "per optimizer unit",
                residual,expected,8*EPS*std::max(std::abs(residual),std::abs(expected)),0,field("prediction",exact.predictions[row]));
        }
        if (held) {
            const auto value = [&](double shift) { return rebuilt(shift).born_diameter_angstrom[1]; };
            affine_difference("affine-columns/Born",exact.chart.map[3],value,1e-3,0.2);
            agree("affine-columns/held-energy", "hydration_error", "J/mol",exact.chart.holds[0].energy_error,
                0,1e-6,0,"");
        }
    }
}
void check_coupled_recovery() {
    FitRequest request; request.parameters = coupled_lengths(3.0); request.coordinates = coupled_coordinate();
    request.controls = budgeted_controls();
    request.holds = {coupled_hold()}; request.observations = coupled_observations(request.holds[0].reference);
    const HelmholtzModel truth_model(coupled_lengths(3.2));
    for (auto& ob : request.observations) ob.observed = *predict(truth_model,ob).value;
    request.controls.solver_options.values[epcsaft::policy::function] = 1e-15; request.controls.solver_options.values[epcsaft::policy::gradient] = 1e-14;
    request.controls.solver_options.values[epcsaft::policy::parameter] = 1e-14;
    const auto result = fit(request);
    CHECK(result.usable);
    agree("affine-recovery/sigma", "sigma", "angstrom",result.parameters.sigma_constant[1],3.2,1e-6,0,
        summary(result));
    for (const auto family : {ParameterFamily::PackingDiameter,ParameterFamily::DebyeHuckelDiameter})
        agree("affine-recovery/"+std::to_string(static_cast<int>(family)), "proportional_diameter", "angstrom",
            detail::family_values(result.parameters,family).first[1],0.88*3.2,0.88e-6,0,"");
    agree("affine-recovery/stored-hydration", "energy_error", "J/mol",
        *predict(HelmholtzModel(result.parameters),request.holds[0].reference).value-request.holds[0].reference.observed,
        0,1e-6,0,"");
    for (std::size_t row=0;row<request.observations.size();++row)
        agree("affine-recovery/decade/"+std::to_string(row),"log10_ratio","decade",
            std::log10(*predict(HelmholtzModel(result.parameters),request.observations[row]).value/request.observations[row].observed),
            0,1e-8,0,"");
}

void check_d86() {
    const std::vector<ActiveParameter> active = {{ParameterFamily::DispersionEnergyOverK, {0}},
                                                 {ParameterFamily::SegmentCount, {0}},
                                                 {ParameterFamily::SigmaConstant, {0}}};
    const double T = 200.0, P = 217800.6949469248;
    const HelmholtzModel model(properties_support::ethane(), active);
    std::vector<PressureActionDirection> directions;
    for (const auto& item : active) directions.push_back({0.0, 0.0, {0.0}, {{item.family, item.indices, 1.0}}});
    const auto batch = model.pressure_actions(T, P, {1.0}, DensityBranch::Liquid, {}, directions);
    CHECK(batch.density_status == PropertyActionStatus::Available && batch.density_action.size() == 3
          && batch.results.empty() && batch.molar_density.has_value());
    const auto central = model.state_at_pressure(T, P, {1.0}, DensityBranch::Liquid);
    agree("d86/density_value", "molar_density", "mol/m^3", *batch.molar_density, central.molar_density, 0.0, 0.0,
          fields({field("T_K", T), field("P_Pa", P)}));
    const double steps[3] = {0.05, 5.0e-4, 5.0e-4};
    const char* names[3] = {"epsilon", "m", "sigma"};
    for (std::size_t k = 0; k < 3; ++k) {
        const auto rebuilt = [&](double shift) {
            auto p = properties_support::ethane();
            (k == 0 ? p.dispersion_energy_over_k : k == 1 ? p.segment_count : p.sigma_constant)[0] += shift;
            return p;
        };
        const auto closure = [&](double shift) {
            return HelmholtzModel(rebuilt(shift)).state_at_pressure(T, P, {1.0}, DensityBranch::Liquid).molar_density;
        };
        const auto inputs = fields({field("T_K", T), field("P_Pa", P), field("step", steps[k])});
        ladder_agrees(std::string("d86/parameter/") + names[k], "drho_dtheta", batch.density_action[k], ladder(closure, steps[k]),
                      inputs);
        // The implicit identity: drho/dtheta|P = -(dP/dtheta|rho) / (dP/drho), the explicit pressure slope from
        // rebuilt models at the resolved density; an explicit-only stand-in (zero) is off by the whole action.
        const auto explicit_pressure = [&](double shift) {
            return HelmholtzModel(rebuilt(shift)).state(T, central.molar_density, {1.0}).pressure;
        };
        const auto slope = ladder(explicit_pressure, steps[k]);
        agree(std::string("d86/implicit/") + names[k], "drho_dtheta", "native", batch.density_action[k],
              -slope.fine / central.pressure_density_derivative, FIRST_ORDER * std::max(1.0, std::abs(batch.density_action[k])), 0.0,
              inputs);
        truth(std::string("d86/omitted_term/") + names[k],
              std::abs(batch.density_action[k]) > 1.0e3 * FIRST_ORDER * std::max(1.0, std::abs(batch.density_action[k])),
              field("action", batch.density_action[k]));
    }
    // State directions: temperature, pressure and a zero-sum composition direction on the neutral binary.
    const auto thermal = model.pressure_actions(T, P, {1.0}, DensityBranch::Liquid, {}, {{1.0, 0.0, {0.0}, {}}, {0.0, 1.0, {0.0}, {}}});
    CHECK(thermal.density_status == PropertyActionStatus::Available && thermal.density_action.size() == 2);
    const auto over_t = [&](double shift) { return model.state_at_pressure(T + shift, P, {1.0}, DensityBranch::Liquid).molar_density; };
    ladder_agrees("d86/temperature", "drho_dT", thermal.density_action[0], ladder(over_t, 0.05),
                  fields({field("T_K", T), field("P_Pa", P)}));
    agree("d86/pressure", "drho_dP", "mol/(m^3 Pa)", thermal.density_action[1], 1.0 / central.pressure_density_derivative, 0.0, 1.0e-12,
          field("dP_drho", central.pressure_density_derivative));
    const HelmholtzModel binary(properties_support::neutral_binary());
    const double Tb = 180.0, Pb = 2.0 * *equilibrium::solve_equilibrium(binary, equilibrium_support::bubble_at(180.0)).pressure;
    const auto mixed = binary.pressure_actions(Tb, Pb, {0.5, 0.5}, DensityBranch::Liquid, {}, {{0.0, 0.0, {0.01, -0.01}, {}}});
    CHECK(mixed.density_status == PropertyActionStatus::Available);
    const auto over_x = [&](double shift) {
        return binary.state_at_pressure(Tb, Pb, {0.5 + 0.01 * shift, 0.5 - 0.01 * shift}, DensityBranch::Liquid).molar_density;
    };
    ladder_agrees("d86/composition", "drho_dx", mixed.density_action[0], ladder(over_x, 0.1),
                  fields({field("T_K", Tb), field("P_Pa", Pb), field("direction", std::vector<double>{0.01, -0.01})}));
    // Failure meaning and the association route.
    const auto absent = model.pressure_actions(T, 1.0e10, {1.0}, DensityBranch::Liquid, {}, directions);
    truth("d86/density_unavailable",
          absent.density_status == PropertyActionStatus::DensityUnavailable && absent.density_action.empty()
              && !absent.molar_density.has_value(), "");
    // The root search's cause travels with the batch; a record-less throw is an absent record, distinct
    // from the measured (finite, possibly zero) pressure residual of a resolved root.
    // An anchor outside the density domain throws with a record: the batch carries it, with its stop.
    const auto anchored = model.pressure_actions(T, P, {1.0}, DensityBranch::Liquid, {}, directions, 1.0e-6);
    truth("d86/density_failure_record",
          anchored.density_status == PropertyActionStatus::DensityUnavailable && anchored.density_convergence.has_value()
              && anchored.density_convergence->stop == core::Stop::left_domain && !anchored.density_message.empty(),
          fields({field("message", anchored.density_message), field("anchor", 1.0e-6)}));
    truth("d86/density_failure_cause",
          absent.density_message.find("no stable density root") != std::string::npos && !absent.density_convergence.has_value()
              && batch.density_message.empty() && !batch.density_convergence.has_value()
              && std::isfinite(batch.diagnostics.pressure_residual),
          fields({field("message", absent.density_message), field("resolved_pressure_residual", batch.diagnostics.pressure_residual)}));
    // Reuse the admitted state with a density-root tolerance below double precision for the action closure.
    const auto row = density(T, P, {1.0}, DensityBranch::Liquid);
    const auto prediction = predict_validated(model, row);
    CHECK(prediction.status == ObservationStatus::Available && prediction.state && prediction.value);
    policy::SolverOptions strict;
    strict.values[policy::density] = 1e-30;
    const HelmholtzModel parameter_model(model.parameters(), active, {}, strict);
    const HelmholtzModel input_model(model.parameters(), {}, {}, strict);
    ObservationInputDirection input;
    input.state.pressure_pa = 1.0;
    input.state.feed_amounts_mol = {0.0};
    for (const bool measured : {false, true}) {
        const auto expected = measured
            ? input_model.pressure_actions(T, P, {1.0}, row.branch, {}, {{0.0, 1.0, {0.0}, {}}}, prediction.state->molar_density)
            : parameter_model.pressure_actions(T, P, {1.0}, row.branch, {}, directions);
        CHECK(expected.density_status == PropertyActionStatus::DensityUnavailable && expected.density_convergence
            && expected.density_convergence->criterion && !expected.density_structure_refused);
        const auto gradient = measured ? differentiate_group(input_model, {&row}, prediction, {}, {input}).front()
                                       : differentiate(parameter_model, row, prediction);
        const auto record = core::describe(*expected.density_convergence);
        CHECK(truth(measured ? "d86/measured-input-refusal-record" : "d86/parameter-refusal-record",
            gradient.status == ObservationStatus::ActionUnavailable && gradient.value.empty() && gradient.refusal_convergence
            && core::describe(*gradient.refusal_convergence) == record
            && gradient.refusal_convergence->evaluations == expected.density_convergence->evaluations,
            fields({field("record", gradient.refusal_convergence ? core::describe(*gradient.refusal_convergence) : "absent"),
                field("closure_record", record)})));
    }
    {
        // The unanchored prediction succeeds at this tolerance, but the measured-input anchor stalls at the floor.
        policy::SolverOptions options;
        options.set(policy::density, 5e-15);
        const HelmholtzModel electrolyte(properties_support::aqueous_electrolyte(), {}, {}, options);
        auto gamma = density(325.0, 2e5, {1.0, 0.0, 0.0}, DensityBranch::Liquid);
        gamma.kind = ObservationKind::MeanIonicActivityCoefficient;
        gamma.molality = 5.0; gamma.ion_indices = {1, 2}; gamma.stoichiometry = {1, 1};
        const auto center = predict(electrolyte, gamma);
        CHECK(center.status == ObservationStatus::Available && center.activity && center.value);
        const auto x = formula_composition(gamma.composition, electrolyte.parameters(), gamma.ion_indices,
            gamma.stoichiometry, gamma.molality);
        const auto expected = electrolyte.pressure_actions(gamma.temperature, gamma.pressure, x, gamma.branch,
            {PropertyObservable::LogFugacityCoefficient}, {{0.0, 1.0, {0.0, 0.0, 0.0}, {}}}, center.activity->state.molar_density);
        CHECK(expected.density_status == PropertyActionStatus::DensityUnavailable && expected.density_convergence
            && expected.diagnostics.first_numerical_failure);
        const auto gradient = differentiate_group(electrolyte, {&gamma}, center, {}, {input}).front();
        const auto record = core::describe(*expected.density_convergence);
        CHECK(truth("d86/electrolyte-measured-input-refusal-record",
            gradient.status == ObservationStatus::ActionUnavailable && gradient.value.empty() && gradient.refusal_convergence
            && core::describe(*gradient.refusal_convergence) == record
            && gradient.refusal_convergence->evaluations == expected.density_convergence->evaluations
            && gradient.first_numerical_failure
            && core::describe(*gradient.first_numerical_failure) == core::describe(*expected.diagnostics.first_numerical_failure),
            fields({field("record", gradient.refusal_convergence ? core::describe(*gradient.refusal_convergence) : "absent"),
                field("closure_record", record)})));
    }
    const auto message = recording_support::refusal([&] { (void)model.density_actions(T, 1000.0, {1.0}, {}, {{0.0, 0.0, {0.0}, {}}}); });
    truth("d86/density_chart_rejected",
          message.find("observables") != std::string::npos, field("message", message));
    // A neutral insertion is its exact infinite-dilution limit (#73): f(0) = 2 f(eps) - f(2 eps) + O(eps^2).
    const auto inserted = [&](double eps) {
        return binary.pressure_actions(Tb, Pb, {1.0 - eps, eps}, DensityBranch::Liquid, {}, {{0.0, 0.0, {-0.01, 0.01}, {}}});
    };
    const auto support = inserted(0.0);
    const double limit = 2.0 * inserted(1.0e-7).density_action.at(0) - inserted(2.0e-7).density_action.at(0);
    truth("d86/support_change",
          support.density_status == PropertyActionStatus::Available && support.density_action.size() == 1
              && std::abs(support.density_action[0] - limit) <= 1.0e-9 * std::abs(limit),
          fields({field("status", static_cast<int>(support.density_status)), field("limit", limit),
                  field("drho", support.density_action.empty() ? 0.0 : support.density_action[0])}));
    // Unfalsified arms of the density status: AssociationResidual and the ill-conditioned/singular site causes
    // share the caloric gates but the P4 stress fixture throws in state() before a pressure-chart request forms;
    // PressureSingular is unreachable here because the root search admits only positive pressure slopes.
    const std::vector<ActiveParameter> site = {{ParameterFamily::AssociationEnergyOverK, {1, 4}}};
    const HelmholtzModel associating(properties_support::associating_binary(), site);
    const auto pair = associating.pressure_actions(300.0, 1.0e5, {0.4, 0.6}, DensityBranch::Liquid, {}, {{0.0, 0.0, {0.0, 0.0}, {{site[0].family, site[0].indices, 1.0}}}});
    CHECK(pair.density_status == PropertyActionStatus::Available);
    const auto over_site = [&](double shift) {
        auto p = properties_support::associating_binary();
        p.association_energy_over_k[1] += shift;
        p.association_energy_over_k[4] += shift;
        return HelmholtzModel(p).state_at_pressure(300.0, 1.0e5, {0.4, 0.6}, DensityBranch::Liquid).molar_density;
    };
    ladder_agrees("d86/association", "drho_depsAB", pair.density_action[0], ladder(over_site, 1.0),
                  fields({field("T_K", 300.0), field("P_Pa", 1.0e5)}));
}

// ---- R2 J1: prediction gradients against rebuilt-model central differences, then the residual chain.
CoordinateSet binary_coordinates() {
    CoordinateSet set;
    set.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 150.0, 10.0, 130.0, 170.0},
                       {{{ParameterFamily::Kij, {1, 2}}}, 0.0, 0.1, -0.2, 0.2},
                       {{{ParameterFamily::SegmentCount, {1}}}, 1.6, 0.1, 1.0, 2.5},
                       {{{ParameterFamily::KijReciprocalTemperatureSlope, {1, 2}}}, 0.0, 10.0, -100.0, 100.0}};
    return set;
}

void check_j1() {
    // k_ij(T) = k_ij + b (1/T - 1/T_ref) with b = 0 at T_ref = 200 K: the model is unchanged, dk/db is not.
    auto parameters = properties_support::neutral_binary();
    parameters.k_ij_reference_temperature_k[1] = parameters.k_ij_reference_temperature_k[2] = 200.0;
    const auto set = binary_coordinates();
    validate(set, parameters);
    const HelmholtzModel model(parameters, set.active_parameters());
    const auto tie = equilibrium_support::interior(model);
    const auto bubble = equilibrium_support::bubble_at(180.0);
    const auto bubble_solve = equilibrium::solve_equilibrium(model, bubble);
    const double P_bubble = *bubble_solve.pressure;
    auto free_t = bubble;
    free_t.temperature = FreeValue::variable(180.0);
    free_t.pressure = FreeValue::fixed(P_bubble);
    struct Case final { const char* name; Observation observation; };
    std::vector<Case> cases = {
        {"bubble_pressure", solved(ObservationKind::EquilibriumPressure, bubble)},
        {"bubble_temperature", solved(ObservationKind::EquilibriumTemperature, free_t)},
        {"vapor_methane", solved(ObservationKind::PhaseMoleFraction, bubble, 1, {0})},
        {"methane_partial_pressure", solved(ObservationKind::PartialPressure, bubble, 1, {0})},
        {"phase_fraction", solved(ObservationKind::PhaseFraction, tie.problem, 0)},
        {"liquid_mass_density", solved(ObservationKind::PhaseDensity, bubble, 0, {}, true)},
        {"density", density(180.0, 2.0 * P_bubble, {0.5, 0.5}, DensityBranch::Liquid, true)},
    };
    const auto u0 = set.start(parameters);
    for (const auto& item : cases) {
        const std::optional<EquilibriumResult> start = item.observation.kind == ObservationKind::PhaseFraction
            ? std::optional<EquilibriumResult>(tie.start)
            : item.observation.kind == ObservationKind::EquilibriumTemperature ? std::optional<EquilibriumResult>(bubble_solve)
            : std::nullopt;
        const auto prediction = predict(model, item.observation, start);
        CHECK(prediction.status == ObservationStatus::Available);
        const auto gradient = differentiate(model, item.observation, prediction);
        CHECK(gradient.status == ObservationStatus::Available && gradient.value.size() == 4);
        if (std::string(item.name) == "vapor_methane") {
            const double factor = 1.0 / 180.0 - 1.0 / 200.0;
            const double reference = factor * gradient.value[1];
            const double allowance = 1.0e-12 + 1.0e-8 * std::abs(reference);
            agree("r2-j1/reciprocal-temperature-factor", "dy_db", "1/K",
                gradient.value[3], reference, allowance, 0.0,
                fields({field("factor_per_K", factor), field("dy_dkij", gradient.value[1])}));
            CHECK(std::abs(reference) > allowance);
            auto at_reference = item.observation;
            at_reference.problem.temperature = FreeValue::fixed(200.0);
            const auto prediction_at_reference = predict(model, at_reference);
            const auto gradient_at_reference = differentiate(model, at_reference, prediction_at_reference);
            CHECK(prediction_at_reference.status == ObservationStatus::Available
                && gradient_at_reference.status == ObservationStatus::Available);
            agree("r2-j1/reciprocal-temperature-zero", "dy_db", "1/K",
                gradient_at_reference.value[3], 0.0, 1.0e-12, 0.0, field("T_K", 200.0));
        }
        for (std::size_t c = 0; c < 4; ++c) {
            const auto value = [&](double shift) {
                auto u = u0;
                u[c] += shift;
                const HelmholtzModel rebuilt(set.apply(parameters, u), set.active_parameters());
                const auto p = predict(rebuilt, item.observation, start);
                if (p.status != ObservationStatus::Available) throw std::runtime_error("rebuilt prediction failed");
                return *p.value;
            };
            const double scale = set.coordinates[c].scale;
            auto rungs = ladder(value, 0.005);
            rungs.coarse /= scale;
            rungs.fine /= scale;
            ladder_agrees(std::string("r2-j1/") + item.name + "/" + std::to_string(c), "dy_dtheta", gradient.value[c], rungs,
                          fields({field("value", *prediction.value), field("unit", prediction.unit), field("u_step", 0.005)}), c == 3);
        }
    }
    // The residual chain: optimizer_jacobian = sqrt(w) (dr/dy) (dy/dtheta) S diag(scale) with a tie and w != 1.
    auto tied = parameters;
    tied.sigma_constant = {3.6, 3.6};
    CoordinateSet shared;
    shared.coordinates = {{{{ParameterFamily::SigmaConstant, {0}}, {ParameterFamily::SigmaConstant, {1}}}, 3.6, 0.1, 3.0, 4.0},
                          {{{ParameterFamily::DispersionEnergyOverK, {0}}}, 150.0, 10.0, 130.0, 170.0}};
    FitRequest request;
    request.controls = budgeted_controls();
    request.parameters = tied;
    request.coordinates = shared;
    request.observations = {solved(ObservationKind::EquilibriumPressure, bubble),
                            solved(ObservationKind::PhaseMoleFraction, bubble, 1, {0})};
    request.observations[0].observed = 1.0e6;
    request.observations[1].form = ResidualForm::Difference;
    request.observations[1].observed = 0.9;
    request.observations[1].scale = 0.01;
    request.weights = {2.0, 0.5};
    request.controls.maximum_iterations = 1;
    const auto result = fit(request);
    truth("r2-j1/chain_setup", result.usable,
          fields({summary(result), field("message", result.message)}));
    CHECK(result.optimizer_jacobian.size() == 4 && result.physical_jacobian.size() == 6);
    for (std::size_t i = 0; i < 2; ++i) {
        const double slope = std::sqrt(request.weights[i]) * residual_slope(request.observations[i], result.predictions[i]);
        const double* j = result.physical_jacobian.data() + i * 3;
        const double hand[2] = {slope * (j[0] + j[1]) * 0.1, slope * j[2] * 10.0};
        // A tied sum may cancel (row 1: y_CH4 against a shared sigma), so the floor is the sum of the terms.
        const double magnitude[2] = {std::abs(slope) * (std::abs(j[0]) + std::abs(j[1])) * 0.1, std::abs(slope * j[2]) * 10.0};
        for (std::size_t c = 0; c < 2; ++c)
            agree("r2-j1/chain/" + std::to_string(i) + "/" + std::to_string(c), "dr_du", "1", result.optimizer_jacobian[i * 2 + c],
                  hand[c], 1.0e-13 * magnitude[c], 1.0e-13, fields({field("weight", request.weights[i]), field("form", i == 0 ? "LogRatio" : "Difference")}));
    }
    // Interleaved physical states, pooled species, fraction and both density bases keep their own rows.
    auto declared = tie.problem;
    for (std::size_t q = 0; q < declared.phases.size(); ++q) {
        declared.phases[q].composition_guess.assign(tie.start.mole_fractions.begin() + q * 2,
                                                    tie.start.mole_fractions.begin() + (q + 1) * 2);
        declared.phases[q].packing_guess = tie.start.packing_fractions[q];
    }
    FitRequest grouped;
    grouped.parameters = parameters;
    grouped.coordinates.coordinates = {set.coordinates.front()};
    grouped.observations = {solved(ObservationKind::PhaseFraction, declared, 0),
        solved(ObservationKind::EquilibriumPressure, equilibrium_support::bubble_at(190.0)),
        solved(ObservationKind::PhaseMoleFraction, declared, 0, {0, 1}),
        solved(ObservationKind::PartialPressure, declared, 1, {0}),
        solved(ObservationKind::PhaseDensity, declared, 0),
        solved(ObservationKind::PhaseDensity, declared, 0, {}, true),
        solved(ObservationKind::EquilibriumPressure, declared),
        solved(ObservationKind::PhaseFraction, declared, 1)};
    for (std::size_t i = 0; i < grouped.observations.size(); ++i) {
        auto& row = grouped.observations[i];
        row.observed = *predict(model, row).value;
        if (i % 2 == 0) { row.form = ResidualForm::Difference; row.scale = row.observed; }
        grouped.weights.push_back(1.0 + i);
    }
    grouped.controls.maximum_iterations = 1;
    grouped.controls.maximum_elapsed_time_seconds = 1.0e-6;
    const auto batch = fit(grouped);
    CHECK(batch.usable && batch.predictions.size() == grouped.observations.size());
    const auto point = grouped.coordinates.start(parameters);
    for (std::size_t i = 0; i < grouped.observations.size(); ++i) {
        const auto& row = grouped.observations[i];
        agree("r2-j1/group_value/" + std::to_string(i), "prediction", "native", batch.predictions[i], row.observed,
              1.0e-10 * row.observed, 1.0e-10, summary(batch));
        const auto value = [&](double shift) {
            auto u = point;
            u[0] += shift;
            const HelmholtzModel rebuilt(grouped.coordinates.apply(parameters, u));
            return *predict(rebuilt, row).value;
        };
        auto rungs = ladder(value, 0.005);
        rungs.coarse /= grouped.coordinates.coordinates.front().scale;
        rungs.fine /= grouped.coordinates.coordinates.front().scale;
        ladder_agrees("r2-j1/group_action/" + std::to_string(i), "dy_dtheta", batch.physical_jacobian[i], rungs,
                      fields({field("value", batch.predictions[i]), field("u_step", 0.005)}));
    }

}

}  // namespace regression_checks
