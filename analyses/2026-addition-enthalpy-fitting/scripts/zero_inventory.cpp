// The #231 zero-inventory readout: native calculations and comparisons on the frozen
// declarations scripts/zero_inventory.py exports (--native-input). It reuses the regression
// test program's helpers; scripts/zero_inventory.py --native-run compiles and runs it.
#include "regression_j2_inventory_checks.cpp"
#include "regression_reaction_checks.cpp"
using namespace regression_checks;

ParameterDocument read_parameter_document(std::istream& in) {
    epcsaft::ParameterDocument doc;
    std::size_t n, sites, entries;
    in >> n >> sites >> entries;
    doc.component_ids.resize(n); doc.site_component.resize(sites); doc.site_multiplicity.resize(sites);
    for (auto& id : doc.component_ids) in >> std::quoted(id);
    for (std::size_t i = 0; i < sites; ++i) in >> doc.site_component[i] >> doc.site_multiplicity[i];
    in >> std::quoted(doc.base_choice) >> std::quoted(doc.association_choice) >> std::quoted(doc.electrolyte_choice)
       >> std::quoted(doc.ionic_dispersion_choice) >> std::quoted(doc.permittivity_choice) >> std::quoted(doc.polar_choice);
    for (std::size_t i = 0; i < entries; ++i) {
        epcsaft::ParameterEntry entry;
        std::size_t width;
        in >> std::quoted(entry.identity) >> std::quoted(entry.field) >> width;
        entry.index.resize(width);
        for (auto& index : entry.index) in >> index;
        in >> entry.magnitude >> std::quoted(entry.unit);
        doc.entries.push_back(entry);
    }
    if (!in) throw std::runtime_error("incomplete parameter declarations");
    return doc;
}

// Frozen #231 inputs are declared by zero_inventory.py; calculations and comparisons stay native.
int zero_inventory_file(const char* filename) {
    std::ifstream in(filename);
    const auto parameters = build_parameters(read_parameter_document(in));
    EquilibriumProblem initial;
    initial.temperature = FreeValue::fixed(313.15);
    initial.pressure = FreeValue::fixed(1734000.0);
    initial.feed.kind = equilibrium::FeedKind::Amounts;
    initial.feed.values = {0.0, 1.0, 7.911062165006594, 0, 0, 0, 0, 0, 0};
    initial.phases = {equilibrium_support::phase("L", PhaseKind::Liquid)};
    initial.source_reference = NeutralReferenceDeclaration{2, {.125, .125, 0, .625, .25, .125, .125, .125, .125}, 1.0};
    initial.solver_options.values[epcsaft::policy::row] = 1e-10;
    for (std::size_t r = 0; r < 5; ++r) {
        equilibrium::ReactionSpec reaction;
        in >> std::quoted(reaction.name);
        reaction.stoichiometry.resize(9);
        for (auto& value : reaction.stoichiometry) in >> value;
        ReactionLogPolynomial law;
        std::string basis;
        double pressure;
        in >> law.a >> law.b >> law.c >> law.d >> law.reference_temperature >> law.temperature_min
           >> law.temperature_max >> std::quoted(law.standard_state_id) >> std::quoted(basis) >> pressure;
        reaction.correlation = law;
        reaction.reference = ReactionReference{pressure == 0 ? std::nullopt : std::optional<double>(pressure),
            basis == "RawMoleFractionInfiniteDilution" ? ReferenceSourceBasis::RawMoleFractionInfiniteDilution
                                                     : ReferenceSourceBasis::CommonMolalityInfiniteDilution};
        initial.reactions.push_back(std::move(reaction));
    }
    ThermochemistryRecord record;
    in >> record.standard_pressure;
    record.component_intervals.resize(9);
    for (std::size_t i = 0; i < 3; ++i) {
        std::size_t count;
        double lower, upper, formation, centre;
        in >> count;
        std::vector<double> coefficients(count);
        for (auto& value : coefficients) in >> value;
        in >> lower >> upper >> formation >> centre;
        IdealCorrelation correlation;
        correlation.source_pressure = record.standard_pressure;
        if (centre == 0) {
            IdealShomate shomate;
            std::copy(coefficients.begin(), coefficients.end(), shomate.coefficients.begin());
            shomate.formation_enthalpy = formation;
            correlation.form = shomate;
        } else {
            const auto& c = coefficients;
            correlation.form = IdealPolynomial{{c[0] + centre * (c[1] + centre * (c[2] + centre * c[3])),
                c[1] + centre * (2 * c[2] + centre * 3 * c[3]), c[2] + centre * 3 * c[3], c[3]}, centre, 0, 0};
        }
        record.component_intervals[i].push_back({lower, upper, true, true, correlation});
    }
    if (!in) throw std::runtime_error("incomplete frozen MEA input declarations");
    std::vector<ActiveParameter> active;
    for (const auto pair : {std::pair{4u, 2u}, {3u, 2u}, {5u, 2u}, {4u, 3u}})
        active.push_back({ParameterFamily::Kij, {9 * pair.second + pair.first, 9 * pair.first + pair.second}});
    active.push_back({ParameterFamily::KijReciprocalTemperatureSlope, {21, 29}});
    CoordinateSet coordinates;
    for (const auto& item : active) coordinates.coordinates.push_back({{{item.family, item.indices}}, 0, 1});
    coordinates.coordinates.push_back(reaction_coordinate("R2", ReactionCoefficient::B));
    coordinates.coordinates.push_back(reaction_coordinate("R2", ReactionCoefficient::C));
    const std::array<const char*, 7> coordinate_names{"k_MEACOO_water", "k_MEAH_water", "k_HCO3_water", "k_MEACOO_MEAH",
        "inverse_T_k_MEAH_water", "R2_b", "R2_c"};
    const std::array<const char*, 7> derivative_units{"J/mol", "J/mol", "J/mol", "J/mol", "J/(mol K)", "J/(mol K)", "J/mol"};
    HelmholtzModel model(parameters, active, record);
    auto boundary = initial;
    boundary.pressure = FreeValue::variable();
    auto vapor = equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0));
    vapor.support = {0, 1, 2}; vapor.support_declared = true; boundary.phases.push_back(vapor);
    const auto bubble = equilibrium::solve_equilibrium(model, boundary);
    check_inventory_boundary(model, boundary, bubble, "incipient/mea-zero", 1e-7);
    const auto face = equilibrium::compile_problem(model, boundary);
    const auto center = face.prepare(bubble.coordinates);
    truth("incipient/mea-zero/balance-rank", "three conserved balances and a fourteen-dimensional endpoint retain the two reaction motions",
  face.rank == 3 && face.rows.size() == 14 && face.removed_species == std::vector<std::size_t>({0, 4, 5, 6}),
  fields({field("removed_species", std::vector<double>{0, 4, 5, 6})}), "issue-269");
    for (const auto r : {0u, 4u}) {
        double affinity = 0, motion_error = 0;
        for (const auto i : face.problem.phases[0].support)
            affinity += boundary.reactions[r].stoichiometry[i] * (center.derivative[0].gradient[2 + i]
                + std::log(center.composition[0][i] / center.volume[0]) + center.standard.derivative[0][i]);
        for (const auto& row : face.B) {
            double contraction = 0;
            for (std::size_t i = 0; i < row.size(); ++i) contraction += row[i] * boundary.reactions[r].stoichiometry[i];
            motion_error = std::max(motion_error, std::abs(contraction));
        }
        agree("incipient/mea-zero/affinity/" + std::to_string(r), "reaction_affinity_over_RT", "1", affinity, 0, 1e-10, 0, "R1/R5", "surviving reactions meet chemical stationarity", "issue-269");
        agree("incipient/mea-zero/motion/" + std::to_string(r), "conserved_motion", "mol", motion_error, 0, 1e-12, 0, "R1/R5", "reaction motions preserve the face balances", "issue-269");
    }
    truth("incipient/mea-zero/separation", "the downstream pressure guard holds and surviving liquid ions remain positive",
  1734000 > *bubble.pressure * (1 + 1e-7) && bubble.mole_fractions[3] > 0 && bubble.mole_fractions[7] > 0 && bubble.mole_fractions[8] > 0,
  fields({field("bubble_pressure_Pa", *bubble.pressure), field("guard_pressure_Pa", *bubble.pressure * (1 + 1e-7)), field("observed_pressure_Pa", 1734000)}), "issue-269");
    auto control = boundary; control.phases[1].support = {1, 2};
    const auto explicit_face = equilibrium::solve_equilibrium(model, control);
    check_inventory_boundary(model, control, explicit_face, "incipient/mea-control", 1e-7);
    reference_scale_agrees("incipient/mea-zero/control-pressure", "boundary_pressure", "Pa", *bubble.pressure, *explicit_face.pressure, 1e-4, 1e-8, "same wheel", "automatic and explicit vapor faces agree", "issue-269");
    for (std::size_t i = 0; i < bubble.mole_fractions.size(); ++i)
        reference_scale_agrees("incipient/mea-zero/control-fraction/" + std::to_string(i), "mole_fraction", "1", bubble.mole_fractions[i], explicit_face.mole_fractions[i], 1e-12, 1e-8, "same wheel", "full chemistry fractions agree", "issue-269");
    for (std::size_t p = 0; p < 2; ++p)
        reference_scale_agrees("incipient/mea-zero/control-density/" + std::to_string(p), "molar_density", "mol/m^3", bubble.molar_densities[p], explicit_face.molar_densities[p], 1e-12, 1e-8, "same wheel", "mechanical branches agree", "issue-269");
    double previous = 1;
    for (const double alpha : {1e-3, 1e-4, 1e-6, 1e-8}) {
        auto loaded = boundary; loaded.feed.values[0] = alpha;
        const auto result = equilibrium::solve_equilibrium(model, loaded);
        const std::string id = "incipient/mea-continuity/" + number(alpha);
        check_inventory_boundary(model, loaded, result, id, 1e-7, false);
        const double discrepancy = std::abs(*result.pressure - *bubble.pressure) / *bubble.pressure;
        truth(id + "/positive-inventory", "positive carbon inventory retains species and approaches the exact face beyond numerical noise",
      equilibrium::compile_problem(model, loaded).removed_species.empty() && result.mole_fractions[9] > 0 && discrepancy <= previous + 1e-8,
      fields({field("alpha_mol_CO2_per_mol_MEA", alpha), field("relative_pressure_difference", discrepancy), field("vapor_CO2_fraction", result.mole_fractions[9])}), "issue-269");
        if (alpha == 1e-8) agree(id + "/limit", "relative_pressure_difference", "1", discrepancy, 0, 1e-7, 0, "unchanged positive inventory ladder", "the final diagnostic approaches the face pressure", "issue-269");
        previous = discrepancy;
    }
    for (const double alpha : {.11, .45}) {
        auto final = initial;
        final.feed.values[0] = alpha;
        const auto row = addition_observation(initial, final, alpha, 313.15, 1734000, {1, 0, 0, 0, 0, 0, 0, 0, 0});
        const auto prediction = predict(model, row);
        if (!prediction.value) throw std::runtime_error("zero-inventory addition: " + prediction.message);
        const std::string prefix = "zero-inventory/" + std::to_string(alpha);
        const auto inputs = fields({field("T_K", 313.15), field("P_Pa", 1734000.0), field("alpha_mol", alpha)});
        const auto compiled = equilibrium::compile_problem(model, initial);
        const auto solved_initial = equilibrium::solve_equilibrium(model, initial);
        const auto final_compiled = equilibrium::compile_problem(model, final);
        const auto solved_final = equilibrium::solve_equilibrium(model, final);
        CHECK(solved_initial.success && solved_final.success);
        for (const auto pair : {std::pair{"initial", &solved_initial}, std::pair{"final", &solved_final}}) {
            double classified = 0, chemical = 0;
            for (std::size_t k = 0; k < pair.second->rows.size(); ++k) {
                classified = std::max(classified, std::abs(pair.second->classified_residuals[k]));
                if (pair.second->rows[k].kind == equilibrium::RowKind::Chemical)
                    chemical = std::max(chemical, std::abs(pair.second->residuals[k]));
            }
            agree(prefix + "/" + pair.first + "/classified", "maximum_classified_row", "1", classified, 0,
          1e-10, 0, inputs, "both endpoint rows meet the unchanged request", "issue-231");
            agree(prefix + "/" + pair.first + "/stationarity", "maximum_raw_chemical_row", "1", chemical, 0,
          1e-10, 0, inputs, "both endpoints resolve raw chemical stationarity", "issue-231");
        }
        Eigen::MatrixXd state_balance(compiled.B.size(), 5);
        const std::array<std::size_t, 5> surviving{1, 2, 3, 7, 8};
        for (std::size_t k = 0; k < compiled.B.size(); ++k)
            for (std::size_t j = 0; j < surviving.size(); ++j) state_balance(k, j) = compiled.B[k][surviving[j]];
        const Eigen::FullPivLU<Eigen::MatrixXd> balance_rank(state_balance);
        truth(prefix + "/nullspace", "the surviving conserved balances have rank three and reaction dimension two",
      balance_rank.rank() == 3 && balance_rank.dimensionOfKernel() == 2, inputs, "issue-231");
        const auto state = model.state(313.15, solved_initial.molar_densities[0], solved_initial.mole_fractions);
        const auto prepared = compiled.prepare(solved_initial.coordinates);
        for (std::size_t r : {0u, 4u}) {
            double affinity = 0;
            Eigen::VectorXd motion(5);
            for (std::size_t j = 0; j < surviving.size(); ++j) {
                const auto i = surviving[j];
                motion(j) = initial.reactions[r].stoichiometry[i];
                affinity += motion(j) * (state.residual_chemical_potential_over_rt[i].value()
                    + std::log(solved_initial.mole_fractions[i] * solved_initial.molar_densities[0]) + prepared.standard.derivative[0][i]);
            }
            agree(prefix + "/reaction-span/" + std::to_string(r), "conserved_motion_norm", "1",
          (state_balance * motion).norm(), 0, 1e-12, 0, inputs, "R1 and R5 belong to the conserved face reaction space", "issue-231");
            agree(prefix + "/affinity/" + std::to_string(r), "reaction_affinity_over_RT", "1", affinity, 0, 1e-10, 0,
          inputs, "public chemical potentials satisfy R1 and R5 on the five-species face", "issue-231");
        }
        const double h_initial = endpoint_enthalpy(model, initial, compiled, solved_initial).enthalpy;
        const double h_final = endpoint_enthalpy(model, final, final_compiled, solved_final).enthalpy;
        auto feed_row = row;
        feed_row.kind = ObservationKind::FeedEnthalpy;
        const auto feed = predict(model, feed_row);
        CHECK(feed.value);
        reference_scale_agrees(prefix + "/feed-h", "h_feed", "J/mol", *feed.value, -393574.4494137983,
      1e-6, 1e-9, inputs, "the support-restricted feed calorics preserve the frozen CO2 result", "issue-231");
        const double direct = (h_final - h_initial) / alpha - *feed.value;
        const double allowance = 1e-9 * (std::abs(h_initial) + std::abs(h_final)) / alpha + 1e-6;
        agree(prefix + "/enthalpy", "delta_h_add", "J/mol", *prediction.value, direct, allowance, 0,
      inputs + ";" + fields({field("H_initial_J", h_initial), field("H_final_J", h_final), field("h_feed_J_per_mol", *feed.value)}),
      "separate public endpoint and feed enthalpies reproduce finite addition", "issue-231");
        const double conditioning = (std::abs(h_initial) + std::abs(h_final)) / (alpha * std::abs(direct));
        agree(prefix + "/conditioning", "enthalpy_subtraction_condition", "1", conditioning, 0, 1e4, 0, inputs,
      "endpoint subtraction conditioning is at most 1e4", "issue-231");
        if (alpha == .11) {
            agree(prefix + "/extrapolation", "delta_h_add", "J/mol", *prediction.value,
          -68169.42621165613, .01, 0, inputs, "exact zero endpoint agrees with the retained neighbouring-state trend", "issue-231");
            agree(prefix + "/initial-H", "H_initial", "J", h_initial, -2313691.892440736, .001, 0,
          inputs, "fresh-solvent enthalpy agrees with the retained extrapolation", "issue-231");
            const std::vector<std::size_t> order{8, 7, 6, 5, 4, 3, 2, 1, 0};
            auto permuted = initial;
            permuted.feed.values = equilibrium_support::permute_vector(initial.feed.values, order);
            permuted.source_reference->solvent_index = 6;
            permuted.source_reference->solute_molality_ratios = equilibrium_support::permute_vector(
                initial.source_reference->solute_molality_ratios, order);
            for (auto& reaction : permuted.reactions) {
                reaction.stoichiometry = equilibrium_support::permute_vector(reaction.stoichiometry, order);
                for (double& value : reaction.stoichiometry) value = -value;
                auto& law = *reaction.correlation;
                law.a = -law.a; law.b = -law.b; law.c = -law.c; law.d = -law.d;
            }
            std::reverse(permuted.reactions.begin(), permuted.reactions.end());
            auto thermal = record;
            thermal.component_intervals = equilibrium_support::permute_vector(record.component_intervals, order);
            HelmholtzModel permuted_model(equilibrium_support::permute_parameters(parameters, order), {}, thermal);
            const auto at = equilibrium::solve_equilibrium(permuted_model, permuted);
            CHECK(at.success);
            for (std::size_t i = 0; i < 9; ++i)
                reference_scale_agrees(prefix + "/permutation/" + std::to_string(i), "component_amount", "mol",
              at.amounts[0] * at.mole_fractions[i], solved_initial.amounts[0] * solved_initial.mole_fractions[order[i]],
              1e-12, 1e-8, inputs, "component/declaration permutation and reaction reversal preserve physical amounts", "issue-231");
        }
        const auto original_basis = equilibrium::balance_basis(9, initial.reactions);
        for (const auto endpoint : {std::tuple{"initial", &initial, &solved_initial}, std::tuple{"final", &final, &solved_final}}) {
            const auto& [name, problem, solved] = endpoint;
            double charge = 0.0;
            std::vector<double> delta(9);
            for (std::size_t i = 0; i < 9; ++i) {
                delta[i] = solved->amounts[0] * solved->mole_fractions[i] - problem->feed.values[i];
                charge += parameters.charge_number[i] * delta[i];
            }
            for (std::size_t k = 0; k < original_basis.size(); ++k)
                agree(prefix + "/" + name + "/balance/" + std::to_string(k), "original_inventory_error", "mol",
              std::inner_product(original_basis[k].begin(), original_basis[k].end(), delta.begin(), 0.0),
              0, 1e-7, 0, inputs, "original full-column conserved inventories close", "issue-231");
            agree(prefix + "/" + name + "/charge", "charge_error", "mol", charge, 0, 1e-7, 0, inputs,
          "endpoint electroneutrality cancels the charge gauge", "issue-231");
        }
        equilibrium::SolvedStateActionDirection direction;
        direction.feed_amounts_mol.assign(9, 0);
        direction.active_parameters.assign(active.size(), 0);
        direction.active_parameters[1] = 1;
        std::vector<equilibrium::SolvedStateObservable> observables;
        for (const std::size_t i : {0u, 4u, 5u, 6u}) {
            observables.push_back({equilibrium::SolvedStateObservableKind::PhaseComponentAmount, 0, i});
            observables.push_back({equilibrium::SolvedStateObservableKind::PhaseMoleFraction, 0, i});
        }
        auto carbon_direction = direction;
        carbon_direction.active_parameters.assign(active.size(), 0);
        carbon_direction.feed_amounts_mol[0] = 1;
        const auto face_actions = equilibrium::solved_state_first_actions(compiled, solved_initial, observables, {direction, carbon_direction});
        for (std::size_t i = 0; i < face_actions[0].results.size(); ++i) {
            const auto& value = face_actions[0].results[i];
            truth(prefix + "/absent/" + std::to_string(i), "forced absent amounts, fractions and fixed-face actions are exact zeros",
          value.value && value.action && *value.value == 0 && *value.action == 0, inputs, "issue-231");
            const auto& crossing = face_actions[1].results[i];
            truth(prefix + "/crossing/" + std::to_string(i), "carbon insertion refuses only the boundary-crossing action",
          crossing.value && *crossing.value == 0 && !crossing.action
          && crossing.status == equilibrium::SolvedStateActionStatus::SupportChange, inputs, "issue-231");
        }
        const auto gradients = differentiate_group(model, {&row}, prediction, coordinates);
        if (gradients[0].status != ObservationStatus::Available)
            throw std::runtime_error(prefix + "/gradient: " + observation_status_name(gradients[0].status) + ": " + gradients[0].message);
        CHECK(gradients[0].value.size() == 7);
        auto basis_coordinates = coordinates;
        basis_coordinates.coordinates.push_back(reaction_coordinate("R2", ReactionCoefficient::A));
        const auto basis_gradient = differentiate_group(model, {&row}, prediction, basis_coordinates).front();
        CHECK(basis_gradient.value.size() == 8);
        const double slope_term = 313.15 * basis_gradient.value[7] - 313.15 * 313.15 * basis_gradient.value[5];
        agree(prefix + "/c-at-Tref", "d_delta_h_dc", "J/mol", gradients[0].value[6], slope_term / 313.15,
      0, 1e-10, inputs, "the zero source-value c basis retains its nonzero caloric slope action", "issue-231");
        truth(prefix + "/c-omission", "omitting the c-at-Tref slope is resolved beyond the derivative criterion",
      std::abs(gradients[0].value[6]) > 1e3 * 5e-6 * std::max(1.0, std::abs(gradients[0].value[6])), inputs, "issue-231");
        auto initial_h_row = row;
        initial_h_row.kind = ObservationKind::SystemEnthalpy;
        auto final_h_row = initial_h_row;
        final_h_row.problem = final;
        const auto initial_h_prediction = predict(model, initial_h_row, solved_initial, compiled);
        const auto final_h_prediction = predict(model, final_h_row, solved_final, final_compiled);
        const auto initial_h_gradient = differentiate_group(model, {&initial_h_row}, initial_h_prediction, coordinates).front();
        const auto final_h_gradient = differentiate_group(model, {&final_h_row}, final_h_prediction, coordinates).front();
        const auto feed_gradient = differentiate(model, feed_row, feed);
        CHECK(initial_h_gradient.value.size() == 7 && final_h_gradient.value.size() == 7 && feed_gradient.value.size() >= 5);
        const auto start = coordinates.start(parameters, {row}, record);
        for (std::size_t column = 0; column < 7; ++column) {
            const double step = column < 4 || column == 6 ? 1e-4 : 1.0;
            const auto value = [&](double shift) {
                auto point = start;
                point[column] += shift;
                std::vector<Observation> changed{row};
                coordinates.apply(changed, point);
                HelmholtzModel rebuilt(coordinates.apply(parameters, point), {}, record);
                const auto at = predict(rebuilt, changed[0]);
                if (!at.value) throw std::runtime_error("zero-inventory derivative node: " + at.message);
                return *at.value;
            };
            const double plus = value(step), minus = value(-step);
            const double fine = (plus - minus) / (2 * step);
            const double coarse = (value(2 * step) - value(-2 * step)) / (4 * step);
            const double error = 5e-6 * std::max(1.0, std::abs(fine));
            const auto column_inputs = fields({field("alpha_mol", alpha), field("coordinate", coordinate_names[column]),
                field("derivative_order", 1), field("step", step), field("fine", fine), field("coarse", coarse),
                field("truncation_estimate", std::abs(fine - coarse) / 3), field("roundoff_estimate", EPS * (std::abs(plus) + std::abs(minus)) / (2 * step)),
                field("normalized_error", std::abs(gradients[0].value[column] - fine) / error)});
            agree(prefix + "/column/" + coordinate_names[column], "d_delta_h_d_coordinate", derivative_units[column],
          gradients[0].value[column], fine, error, 0, column_inputs, "exact fitted-coordinate action agrees with native centred re-solves", "issue-231");
            agree(prefix + "/coarse/" + coordinate_names[column], "d_delta_h_d_coordinate", derivative_units[column],
          gradients[0].value[column], coarse, 4 * error, 0, column_inputs, "coarse derivative rung retains the unchanged criterion", "issue-231");
        }
        auto shifted_record = record;
        std::get<IdealShomate>(shifted_record.component_intervals[0][0].correlation.form).formation_enthalpy += 1e4;
        HelmholtzModel shifted_model(parameters, active, shifted_record);
        const auto shifted = predict(shifted_model, row);
        CHECK(shifted.value);
        agree(prefix + "/anchor", "delta_h_add", "J/mol", *shifted.value, *prediction.value, allowance, 0,
      inputs, "CO2 feed anchor cancels across full caloric completion on unequal supports", "issue-231");
        const auto shifted_gradient = differentiate_group(shifted_model, {&row}, shifted, coordinates);
        for (std::size_t column = 0; column < 7; ++column)
            agree(prefix + "/anchor-column/" + coordinate_names[column], "d_delta_h_d_coordinate", derivative_units[column],
          shifted_gradient[0].value[column], gradients[0].value[column],
          1e-9 * ((std::abs(initial_h_gradient.value[column]) + std::abs(final_h_gradient.value[column])) / alpha
              + (column < feed_gradient.value.size() ? std::abs(feed_gradient.value[column]) : 0.0)), 0,
          inputs + ";" + fields({field("coordinate", coordinate_names[column]), field("derivative_order", 1)}),
          "CO2 anchor cancellation preserves exact first columns", "issue-231");
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    return zero_inventory_file(argv[1]);
}
