#include "regression_checks_shared.hpp"

namespace regression_checks {
void check_saturation_capped_addition() {
#include "fixtures/mea_capped_addition.inc"
  using namespace epcsaft::equilibrium;
  const HelmholtzModel model(parameters, coordinates.active_parameters(),
                             thermal);
  auto invalid = row;
  invalid.loading_cap = 0;
  CHECK(recording_support::refusal([&] { (void)validate_declaration(invalid, model); }) ==
        "observation positive basis and cap");

  const auto initial = compile_problem(model, row.problem);
  const AffineLoading departure{initial.n_feed, row.composition, 0.001, 0.01, 0,
                                row.pressure};
  const auto low =
      search_loading(model, row.final_problem, departure, {0.001, 0.01});
  CHECK(low.nodes.size() == 2);
  for (const auto &node : low.nodes) {
    CHECK(node.state.success && node.state.residual_norm() <= 1e-10 &&
          maximum_abs(node.state.classified_residuals) <= 1e-10);
    CHECK(node.state.iterations < 200);
  }
  agree("capped-addition/ordinary-loading-node", "pressure", "Pa",
        *low.nodes.back().state.pressure, 6467.867870101388, 1e-6, 0,
        fields({field("loading_mol_per_mol", low.nodes.back().loading)}));
  CHECK(low.complete && !low.accepted);
  CHECK(low.causes == std::vector<std::string>{"unresolved_slope_reversal"});

  const auto center = predict(model, row);
  CHECK(center.status == ObservationStatus::Available && center.value &&
        center.saturation);
  CHECK(center.saturation->accepted && center.cap_branch == "saturation_cap");
  CHECK(center.solved->residual_norm() <= 1e-10 &&
        center.final_solved->residual_norm() <= 1e-10);
  const double alpha = *center.endpoint_loading;
  const auto inputs =
      fields({field("T_K", row.temperature), field("P_Pa", row.pressure),
              field("basis_mol_MEA", row.basis_amount),
              field("cap_mol_per_mol", row.loading_cap),
              field("search_iterations", center.saturation->iterations)});
  agree("capped-addition/loading", "saturation_loading", "mol_CO2/mol_MEA",
        *center.saturation_loading, 1.0175305291914531, 1e-8, 0, inputs);
  agree("capped-addition/endpoint", "endpoint_loading", "mol_CO2/mol_MEA",
        alpha, *center.saturation_loading, 0, 0, inputs);

  // An output sigma changes the residual scale. A reaction-only step keeps the EOS reference preparation.
  const auto& reaction = row.problem.reactions.front();
  Coordinate intercept;
  intercept.reaction = ReactionCoordinate{reaction.name};
  intercept.origin = reaction.correlation->a;
  FitRequest request{parameters, thermal, {{intercept}}, {row}, {1.0}};
  request.uncertainty = MeasurementUncertainty{{}, {row.scale}, {}, {}};
  FitResult failures;
  ObservationEvaluator sequence(request, failures, {0}), fresh(request, failures, {0});
  ObservationEvaluation first, next, independent;
  CHECK(sequence.evaluate({0.0}, false, first));
  CHECK(sequence.evaluate({0.001}, true, next));
  CHECK(fresh.evaluate({0.001}, true, independent));
  const auto& nodes = next.saturation_predictions[0]->saturation->nodes;
  CHECK(nodes.back().prior_start == "accepted");
  CHECK(nodes.front().state.model_identity == first.saturation_predictions[0]->saturation->nodes.front().state.model_identity);
  CHECK(std::abs(next.predictions[0] / independent.predictions[0] - 1) <= 1e-12);
  CHECK(std::abs(next.residual_optimizer_gradient[0] / independent.residual_optimizer_gradient[0] - 1) <= 1e-12);

  // The dilute node is a poor start for the last loading. The ordinary continuation must still recover it.
  auto displaced = center;
  displaced.saturation = std::make_shared<SaturationSearch>(*center.saturation);
  displaced.saturation->nodes.back().state = low.nodes.front().state;
  const auto recovered = predict(model, row, {}, {}, nullptr, &displaced);
  CHECK(recovered.value && recovered.saturation->events == center.saturation->events);
  CHECK(std::abs(*recovered.value / *center.value - 1) <= 1e-12);
  CHECK(recovered.saturation->rejected_starts.size() == 1);
  const auto& rejected = recovered.saturation->rejected_starts.front();
  const auto& retried = recovered.saturation->nodes.back();
  CHECK(rejected.loading == retried.loading && rejected.driver_attempts > 0 && retried.driver_attempts > 0);
  CHECK(retried.prior_start.find("fallback:") == 0 && rejected.prior_start == retried.prior_start);

  auto bubble = row.final_problem;
  bubble.feed.values = row.problem.feed.values;
  bubble.feed.values[0] += row.basis_amount * alpha;
  const auto cold = solve_equilibrium(model, bubble);
  CHECK(cold.success && cold.pressure && cold.residual_norm() <= 1e-10);
  reference_scale_agrees("capped-addition/cold-pressure", "pressure", "Pa",
                         *cold.pressure, row.pressure, 0, 1e-9, inputs);
  const auto cold_compiled = compile_problem(model, bubble);
  SolvedStateActionDirection zero;
  zero.feed_amounts_mol.assign(9, 0);
  zero.active_parameters.assign(coordinates.active_parameters().size(), 0);
  const auto enthalpy = [&](const CompiledSystem &compiled,
                            const EquilibriumResult &solved) {
    const auto entry =
        solved_state_actions(
            compiled, solved,
            {{{SolvedStateObservableKind::SystemEnthalpy}}, {zero}})
            .results.front();
    CHECK(entry.status == SolvedStateActionStatus::Available && entry.value);
    return *entry.value;
  };
  const double initial_h = enthalpy(*center.compiled, *center.solved);
  const double final_h = enthalpy(cold_compiled, cold);
  ThermochemistryRecord feed_record;
  feed_record.standard_pressure = thermal.standard_pressure;
  feed_record.component_intervals = {thermal.component_intervals[0]};
  const auto feed_state = model.state_at_pressure(
      row.temperature, row.pressure, row.composition, DensityBranch::Vapor);
  const auto feed_ideal = evaluate_ideal_properties(
      feed_record, row.temperature, feed_state.molar_density, {1});
  CHECK(feed_ideal.value);
  const double feed_h =
      feed_ideal.value->enthalpy + feed_state.residual_enthalpy;
  const double added = row.basis_amount * alpha;
  const double reference = (final_h - initial_h) / added - feed_h;
  const auto heat_inputs =
      fields({field("H_initial_J", initial_h), field("H_endpoint_J", final_h),
              field("h_feed_J_per_mol", feed_h), field("added_mol", added),
              field("cold_residual", cold.residual_norm())});
  reference_scale_agrees("capped-addition/heat", "signed_addition_enthalpy", "J/mol_CO2",
      *center.value, reference, 0, 1e-6, heat_inputs);

  const auto exact =
      differentiate_group(model, {&row}, center, coordinates).front();
  CHECK(exact.status == ObservationStatus::Available &&
        exact.value.size() == 2);
  const auto root_compiled =
      compile_problem(model, *center.final_solved->declaration);
  std::vector<SolvedStateActionDirection> directions(2, zero);
  for (std::size_t k = 0; k < directions.size(); ++k)
    directions[k].active_parameters[k] = 1;
  const auto loading_actions = solved_state_first_actions(
      root_compiled, *center.final_solved,
      {{SolvedStateObservableKind::Loading}}, directions);
  for (const auto &action : loading_actions)
    CHECK(action.results.front().status == SolvedStateActionStatus::Available &&
          action.results.front().action);
  const auto start = coordinates.start(parameters, {row}, thermal);
  // Four retained central-difference rungs resolve second order; these are the
  // fine-pair Richardson allowances.
  const double heat_allowance[] = {
      4 * (-70.556718017905951 + 70.556719653541222) / 3,
      4 * (32.060942376847379 - 32.060940249357373) / 3};
  const double loading_allowance[] = {
      4 * (-0.0042851466573012864 + 0.0042851467169313651) / 3,
      4 * (0.0011072359793695696 - 0.0011072359179742364) / 3};
  constexpr double step = 0.005;
  for (std::size_t column = 0; column < 2; ++column) {
    std::array<std::array<double, 3>, 2> values{};
    for (std::size_t sign = 0; sign < 2; ++sign) {
      auto u = start;
      u[column] += (sign ? 1 : -1) * step;
      const HelmholtzModel perturbed(coordinates.apply(parameters, u),
                                     coordinates.active_parameters(), thermal);
      const auto prediction = predict(perturbed, row);
      CHECK(prediction.status == ObservationStatus::Available &&
            prediction.value && prediction.saturation->accepted);
      CHECK(prediction.cap_branch == "saturation_cap" &&
            row.loading_cap - *prediction.saturation_loading >= 0.15);
      if (column == 0 && sign == 0) {
        const auto continued = predict(perturbed, row, {}, {}, nullptr, &center);
        CHECK(continued.value && continued.cap_branch == prediction.cap_branch);
        CHECK(continued.saturation->events == prediction.saturation->events);
        CHECK(continued.saturation->nodes.size() == row.loading_grid.size());
        CHECK(std::any_of(continued.saturation->nodes.begin(), continued.saturation->nodes.end(),
                          [](const auto& node) { return node.prior_start == "accepted"; }));
        CHECK(std::abs(*continued.value / *prediction.value - 1) <= 1e-12);
      }
      const auto compiled =
          compile_problem(perturbed, *prediction.final_solved->declaration);
      const auto estimate =
          equilibrium_numerical_estimates(compiled, *prediction.final_solved);
      CHECK(estimate.loading.resolved &&
            row.loading_cap - *prediction.saturation_loading >
                estimate.loading.absolute);
      CHECK(prediction.solved->residual_norm() <= 1e-10 &&
            prediction.final_solved->residual_norm() <= 1e-10);
      values[sign] = {*prediction.value, *prediction.saturation_loading,
                      *prediction.endpoint_loading};
    }
    for (std::size_t quantity = 0; quantity < 3; ++quantity) {
      const double difference =
          (values[1][quantity] - values[0][quantity]) / (2 * step);
      const double scale = coordinates.coordinates[column].scale;
      const double action =
          scale * (quantity == 0
                       ? exact.value[column]
                       : *loading_actions[column].results.front().action);
      const double allowance =
          quantity == 0 ? heat_allowance[column] : loading_allowance[column];
      const auto column_inputs =
          fields({field("coordinate", column), field("step_u", step),
                  field("minus", values[0][quantity]),
                  field("plus", values[1][quantity]),
                  field("physical_action", action / scale),
                  field("coordinate_scale", scale)});
      const auto id = "capped-addition/column-" + std::to_string(column) + "/" +
                      std::to_string(quantity);
      agree(id, quantity == 0 ? "dQ_du" : "dloading_du",
            quantity == 0 ? "J/mol/u" : "mol/mol/u", action, difference,
            allowance, 0, column_inputs);
      CHECK(std::abs(action - difference) /
                std::max({1.0, std::abs(action), std::abs(difference)}) <=
            2e-5);
      if (quantity != 0)
        CHECK(std::abs(action) > 10 * loading_allowance[column]);
    }
  }
}
} // namespace regression_checks

namespace regression_checks {
void check_phase_bulk_relative_permittivity() {
    using namespace epcsaft::equilibrium;
    using Family = ParameterFamily;
    using Kind = SolvedStateObservableKind;
#include "fixtures/mea_capped_addition.inc"

    // The reciprocal-temperature column needs T distinct from the fixture's 313.15 K pivot.
    Observation epsilon;
    epsilon.kind = ObservationKind::PhaseBulkRelativePermittivity;
    epsilon.problem = row.problem;
    epsilon.problem.temperature = FreeValue::fixed(323.15);
    epsilon.problem.pressure = FreeValue::fixed(101325);
    epsilon.problem.feed.values[0] = 0.4;
    epsilon.observed = 60;
    epsilon.scale = 3.2;
    auto density = epsilon;
    density.kind = ObservationKind::PhaseDensity;
    density.observed = 55000;
    density.scale = 1000;
    const std::vector<Observation> rows{density, epsilon};

    coordinates.coordinates = {
        {{{Family::Kij, {22, 38}}}, 0, 0.001},
        {{{Family::KijReciprocalTemperatureSlope, {21, 29}}}, 0, 10},
        {{{Family::SigmaConstant, {1}}}, 0, 0.001},
        {{{Family::PackingDiameter, {4}}}, 0, 0.001},
        {{{Family::DebyeHuckelDiameter, {4}}}, 0, 0.001},
        {{{Family::BornDiameter, {4}}}, 0, 0.001},
        {{{Family::DispersionEnergyOverK, {1}}}, 0, 0.1},
        {{{Family::IonSuppressionCoefficient, {3}}}, 0, 0.005},
        {{{Family::IonSuppressionCoefficient, {4}}}, 0, 0.005},
    };
    const char* names[] = {"k_ij", "k_ij_reciprocal_temperature_slope", "segment_diameter", "packing_diameter",
                           "debye_huckel_diameter", "born_diameter", "dispersion_energy_over_k",
                           "ion_suppression_coefficient/meah", "ion_suppression_coefficient/carbamate"};
    for (const char* form : {"solvent-density-temperature-ion-specific-suppression",
                             "solvent-density-temperature-ion-specific-susceptibility-suppression"}) {
        const std::string record = std::string(form) + "/phase-permittivity/";
        auto declaration = document;
        declaration.permittivity_choice = form;
        for (std::size_t i = 3; i < 9; ++i)
            declaration.entries.push_back({"ion_suppression_coefficient", "ion_suppression_coefficient", {i},
                i == 3 ? 1.0 : i == 4 ? 2.0 : 0.0, "dimensionless"});
        parameters = build_parameters(declaration);
        const auto active = coordinates.active_parameters();
        const HelmholtzModel model(parameters, active, thermal);
        const auto center = predict_validated(model, epsilon);
        truth(record + "solved", center.status == ObservationStatus::Available && center.value && center.solved,
              fields({field("cause", center.message)}));
        CHECK(center.value && center.solved && center.compiled && center.unit == "1");
        const auto& state = *center.solved;
        const auto inputs = fields({field("form", form), field("T_K", *state.temperature),
            field("P_Pa", *state.pressure), field("loading", epsilon.problem.feed.values[0]),
            field("rho_mol_m3", state.molar_densities[0]), field("residual_norm", state.residual_norm())});
        CHECK(state.residual_norm() <= 1e-10);
        const double reference = model.state(*state.temperature, state.molar_densities[0], state.mole_fractions).bulk_relative_permittivity;
        CHECK(agree(record + "value", "bulk_relative_permittivity", "1", *center.value, reference,
                    32 * EPS * std::abs(reference), 0, inputs));

        const auto gradient = differentiate(model, epsilon, center);
        truth(record + "action-available", gradient.status == ObservationStatus::Available,
              fields({field("cause", gradient.message)}));
        CHECK(gradient.value.size() == active.size());
        const auto start = coordinates.start(parameters, rows, thermal);
        for (std::size_t column = 0; column < active.size(); ++column) {
            const bool suppression = active[column].family == Family::IonSuppressionCoefficient;
            const double step = coordinates.coordinates[column].scale;
            std::map<double, double> composition_frozen;
            const auto value = [&](double displacement) {
                auto optimizer = start;
                optimizer[column] += displacement / step;
                const HelmholtzModel shifted(coordinates.apply(parameters, optimizer), active, thermal);
                const auto prediction = predict(shifted, epsilon);
                CHECK(prediction.status == ObservationStatus::Available && prediction.value && prediction.solved);
                const auto& moved = *prediction.solved;
                CHECK(moved.residual_norm() <= 1e-10);
                if (suppression)
                    composition_frozen.emplace(displacement,
                        shifted.state(*moved.temperature, moved.molar_densities[0], state.mole_fractions).bulk_relative_permittivity);
                return *prediction.value;
            };
            const auto full = ladder(value, step);
            const auto id = record + names[column];
            CHECK(ladder_agrees(id, "d_epsilon_d_parameter", gradient.value[column], full,
                fields({field("form", form), field("step", step), field("T_K", *state.temperature)}), true));
            if (suppression) {
                const auto frozen = ladder([&](double displacement) { return composition_frozen.at(displacement); }, step);
                const double allowance = FIRST_ORDER * std::abs(full.fine);
                const auto power_inputs = fields({field("frozen_fine", frozen.fine), field("full_fine", full.fine),
                    field("frozen_coarse", frozen.coarse), field("full_coarse", full.coarse), field("fine_allowance", allowance)});
                truth(id + "/composition-power", std::abs(frozen.fine - full.fine) > allowance, power_inputs);
            } else CHECK(full.fine != 0);
        }

        FitRequest request;
        request.parameters = parameters;
        request.thermochemistry = thermal;
        request.coordinates = coordinates;
        request.observations = rows;
        request.controls = budgeted_controls();
        const auto evaluated = evaluated_columns(request, start);
        CHECK(evaluated.statuses == std::vector<ObservationStatus>(rows.size(), ObservationStatus::Available));
        for (std::size_t column = 0; column < active.size(); ++column) {
            const double expected = gradient.value[column] * coordinates.coordinates[column].scale / epsilon.scale;
            CHECK(agree(record + "grouped-fit/" + names[column], "dr_du", "1",
                evaluated.residual_optimizer_gradient.at(active.size() + column), expected,
                1e-13 * std::abs(expected), 0, inputs));
        }

        auto apparent_cp = epsilon;
        apparent_cp.kind = ObservationKind::EquilibriumIsobaricHeatCapacity;
        const auto cp_prediction = predict_validated(model, apparent_cp);
        CHECK(cp_prediction.status == ObservationStatus::Available && cp_prediction.value);
        const auto cp_gradient = differentiate(model, apparent_cp, cp_prediction);
        CHECK(cp_gradient.status == ObservationStatus::Available && cp_gradient.value.size() == active.size());
        CHECK(cp_gradient.value.front() != 0 && gradient.value.front() != 0);
        CHECK(same_state(epsilon, apparent_cp));
        const auto grouped = [&] {
            core::Capture capture(true);
            auto result = differentiate_group(model, {&epsilon, &apparent_cp}, center, coordinates);
            std::vector<double> widths;
            const auto& trace = capture.trace();
            for (const auto& event : trace.events)
                if (std::string_view(event.name) == "results" &&
                    std::string_view(trace.spans.at(event.span).name) == "solved_state_first_actions")
                    widths.push_back(event.values.at(0));
            CHECK(widths.size() == 2 * active.size() + 1);
            CHECK(static_cast<std::size_t>(std::count(widths.begin(), widths.end(), 1.0)) == active.size() + 1);
            CHECK(static_cast<std::size_t>(std::count(widths.begin(), widths.end(), 2.0)) == active.size());
            return result;
        }();
        CHECK(grouped.size() == 2);
        const PredictionGradient* separate[] = {&gradient, &cp_gradient};
        const char* quantities[] = {"bulk_relative_permittivity", "apparent_cp"};
        for (std::size_t o = 0; o < grouped.size(); ++o) {
            CHECK(grouped[o].status == ObservationStatus::Available && grouped[o].value.size() == active.size());
            for (std::size_t k = 0; k < active.size(); ++k) {
                const double expected = separate[o]->value[k];
                CHECK(agree(record + "grouped-cp/" + quantities[o] + "/" + names[k], quantities[o], "native",
                    grouped[o].value[k], expected, 1e-13 * std::abs(expected), 0, inputs));
            }
        }

        SolvedStateActionDirection direction;
        direction.feed_amounts_mol.resize(state.mole_fractions.size());
        direction.active_parameters.resize(active.size());
        bool order_refused = false;
        try {
            (void)solved_state_actions(*center.compiled, state, {{{Kind::PhaseBulkRelativePermittivity}}, {direction, direction}});
        } catch (const std::invalid_argument& error) {
            CHECK(std::string(error.what()) == "bulk relative permittivity supports first actions only");
            order_refused = true;
        }
        CHECK(order_refused);
    }
}
}  // namespace regression_checks

namespace regression_checks {
void check_merkley_saturation_loading() {
    using namespace epcsaft::equilibrium;
    using Family = ParameterFamily;
#include "fixtures/mdea-411-bubble.inc"
#include "fixtures/mdea-thermochemistry.inc"
    // Amine-Thermodynamics 152ad1bb, normalized MDEA-1, Merkley fit:126, trial 2.
    FitRequest request;
    request.parameters = model.parameters();
    request.thermochemistry = thermal;
    request.coordinates.coordinates = {
        {{{Family::TemperatureLawCoefficient, {0}}}, -0.06711345443725586, 0.01, -0.5, -0.06711345443725586},
        {{{Family::TemperatureLawCoefficient, {1}}}, 0.0012092080101259242, 0.001, -0.01, 0.01},
        {{{Family::Kij, {7, 56}}}, 0.25261865658522437, 0.01, -0.5, 0.5},
        {{{Family::Kij, {15, 57}}}, -0.00956762026107625, 0.01, -0.5, 0.5},
        {{{Family::KijReciprocalTemperatureSlope, {15, 57}}}, -15.929380742416088, 10, -1000, 1000},
        {{{Family::Kij, {23, 58}}}, 1, 0.01, -1, 1},
        {{{Family::SigmaConstant, {7}}, {Family::DebyeHuckelDiameter, {7}}, {Family::BornDiameter, {7}}},
         3.6552722619123172, 0.1, 3, 4.5},
        {{{Family::DispersionEnergyOverK, {7}}}, 635.4895968646832, 10, 50, 950},
    };
    auto& row = request.observations.emplace_back();
    row.kind = ObservationKind::SaturationCappedAdditionEnthalpy;
    row.temperature = 333.15;
    row.pressure = 1121000;
    row.observed = -55400;
    row.scale = 3000;
    row.branch = DensityBranch::Vapor;
    row.basis_amount = row.feed_amount = 1;
    row.composition = {1, 0, 0, 0, 0, 0, 0, 0};
    row.loading_cap = 0.941;
    row.loading_bounds = {0.05, 1.5};
    for (int i = 1; i <= 30; ++i) row.loading_grid.push_back(i * 0.05);
    bubble.temperature = FreeValue::fixed(row.temperature);
    bubble.pressure = FreeValue::variable(row.pressure);
    bubble.feed.values = {0.5, 9.9217608607804, 0, 0, 0, 0, 1, 0};
    bubble.phases.back().composition_guess = {0.9, 0.1, 1e-5};
    row.final_problem = bubble;
    row.problem = bubble;
    row.problem.phases.resize(1);
    row.problem.feed.values[0] = 0;
    row.problem.pressure = FreeValue::fixed(row.pressure);
    request.weights = {1};
    FitResult failures;
    ObservationEvaluator evaluator(request, failures, {0});
    ObservationEvaluation start, trial;
    CHECK(evaluator.evaluate(std::vector<double>(8, 0), false, start));
    const std::vector<double> trial_two = {0, -0.013512156522612428, -20.453703287524615, -3.9026605768041516,
        0.15355058202903896, 0, -0.5013581788020577, -9.360489808605179};
    const bool admitted = evaluator.evaluate(trial_two, false, trial);
    CHECK(truth("merkley-saturation-loading/admitted", admitted,
                fields({field("cause", admitted ? "" : failures.failures.at(0).message)})));
    const auto& result = *trial.saturation_predictions.at(0);
    CHECK(result.saturation && result.saturation->complete && result.saturation->accepted);
    const double tolerance = bubble.solver_options.value(policy::row);
    const auto inputs = fields({field("T_K", row.temperature), field("P_Pa", row.pressure), field("cap", row.loading_cap)});
    CHECK(result.endpoint_loading && reference_scale_agrees("merkley-saturation-loading/root", "endpoint_loading",
          "mol_CO2/mol_MDEA", *result.endpoint_loading, 0.908550536030039, 0, tolerance, inputs));
    CHECK(result.value && reference_scale_agrees("merkley-saturation-loading/heat", "signed_addition_enthalpy",
          "J/mol_CO2", *result.value, -58263.86274105281, 0, tolerance, inputs));
}
}  // namespace regression_checks
