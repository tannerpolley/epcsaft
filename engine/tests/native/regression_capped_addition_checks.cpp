#include "regression_checks_shared.hpp"

namespace regression_checks {
void check_saturation_capped_addition() {
#include "mea_capped_addition.inc"
  using namespace epcsaft::equilibrium;
  const HelmholtzModel model(parameters, coordinates.active_parameters(),
                             thermal);
  auto invalid = row;
  invalid.loading_cap = 0;
  CHECK(rejection([&] { (void)validate_declaration(invalid, model); }) ==
        "observation positive basis and cap");

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
        *center.saturation_loading, 1.0175305291914531, 1e-8, 0, inputs,
        "inverse loading reproduces the cold-bubble reference",
        "capped-addition");
  agree("capped-addition/endpoint", "endpoint_loading", "mol_CO2/mol_MEA",
        alpha, *center.saturation_loading, 0, 0, inputs,
        "the saturation-selected endpoint equals saturation loading",
        "capped-addition");

  auto bubble = row.final_problem;
  bubble.feed.values = row.problem.feed.values;
  bubble.feed.values[0] += row.basis_amount * alpha;
  const auto cold = solve_equilibrium(model, bubble);
  CHECK(cold.success && cold.pressure && cold.residual_norm() <= 1e-10);
  reference_scale_agrees("capped-addition/cold-pressure", "pressure", "Pa",
                         *cold.pressure, row.pressure, 0, 1e-9, inputs,
                         "ordinary cold bubble reproduces physical pressure",
                         "capped-addition");
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
  reference_scale_agrees(
      "capped-addition/heat", "signed_addition_enthalpy", "J/mol_CO2",
      *center.value, reference, 0, 1e-6, heat_inputs,
      "two referenced enthalpies and pure feed independently compose the "
      "signed heat",
      "capped-addition");

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
            allowance, 0, column_inputs,
            "native exact action agrees with the retained-ladder selected "
            "central difference",
            "capped-addition");
      CHECK(std::abs(action - difference) /
                std::max({1.0, std::abs(action), std::abs(difference)}) <=
            2e-5);
      if (quantity != 0)
        CHECK(std::abs(action) > 10 * loading_allowance[column]);
    }
  }
}
} // namespace regression_checks
