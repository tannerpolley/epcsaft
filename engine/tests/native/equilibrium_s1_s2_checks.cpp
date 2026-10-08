// S1 and S2: synthetic flash, bubble and reference-reaction solves.
#include "reactive_caloric_shared.hpp"

namespace equilibrium_checks {

void verify_s1(const epcsaft::HelmholtzModel& model,
               const epcsaft::equilibrium::EquilibriumProblem& problem,
               const epcsaft::equilibrium::EquilibriumResult& result,
               const char* case_id) {
    CHECK(result.success);
    const double limit = convergence_limit(result.convergence);
    CHECK(result.residuals.size() == result.rows.size());
    for (std::size_t row = 0; row < result.residuals.size(); ++row) {
        const auto inputs = recording_support::fields({
            recording_support::field("T", result.temperature.front()),
            recording_support::field("P", result.pressure.front()),
            recording_support::field("coordinates", result.coordinates),
            recording_support::field("row", row), recording_support::field("limit", limit),
        });
        recording_support::compare(
            (std::string("residual/") + case_id + "/" + std::to_string(row)).c_str(),
            "equilibrium_residual", "1", result.residuals[row], 0.0, inputs,
            "Relationship", "equilibrium-residual", "residual is within convergence limit",
            "reference_or_one", limit, 0.0, std::abs(result.residuals[row]) <= limit,
            "s1-residual", (std::string(case_id) + "/" + std::to_string(row)).c_str(),
            "r", "T,P,\\mathbf{z},row",
            "|r_i|\\leq\\ell", "row");
    }
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(result.convergence.residual == compiled.classified_residual_norm(result.coordinates));
    const std::size_t components = model.parameters().component_count();
    const double rt = epcsaft::GAS_CONSTANT_J_PER_MOL_K * result.temperature.front();
    const auto& charges = model.parameters().charge_number;
    for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
        const auto first = result.mole_fractions.begin()
            + static_cast<std::ptrdiff_t>(phase * components);
        const std::vector<double> x(first, first + static_cast<std::ptrdiff_t>(components));
        const double volume = 1.0 / result.molar_densities[phase];
        const auto derivative = model.evaluate(
            result.temperature.front(), volume, x, epcsaft::DerivativeOrder::Gradient);
        double gibbs = derivative.value + result.pressure.front() * volume / rt;
        double lambda_Bx = 0.0;
        double charge = 0.0;
        std::vector<double> chemical(components);
        for (std::size_t species = 0; species < components; ++species) {
            chemical[species] = derivative.gradient[2 + species]
                + std::log(std::max(x[species] / volume,
                                    std::numeric_limits<double>::min()))
                + compiled.standard_gibbs[species] / rt;
            gibbs += x[species] * (
                std::log(std::max(x[species], std::numeric_limits<double>::min()))
                - std::log(volume) - 1.0 + compiled.standard_gibbs[species] / rt);
            charge += charges[species] * x[species];
            for (std::size_t row = 0; row < compiled.rank; ++row)
                lambda_Bx += result.multipliers[row] * compiled.B[row][species] * x[species];
        }
        const auto phase_inputs = recording_support::fields({
            recording_support::field("T", result.temperature.front()),
            recording_support::field("P", result.pressure.front()),
            recording_support::field("x", x), recording_support::field("phase", phase),
            recording_support::field("limit", limit),
        });
        recording_support::compare(
            (std::string("gibbs-dual/") + case_id + "/" + std::to_string(phase)).c_str(),
            "gibbs_dual_residual", "1", gibbs, lambda_Bx, phase_inputs,
            "Relationship", "KKT", "phase Gibbs value equals balance dual",
            "reference_or_one", limit, 0.0, std::abs(gibbs - lambda_Bx) <= limit,
            "s1-gibbs-dual", (std::string(case_id) + "/" + std::to_string(phase)).c_str(),
            "g", "T,P,\\mathbf{x},phase",
            "g=\\lambda^TBx", "phase");
        bool ion_phase = false;
        for (const auto species : compiled.problem.phases[phase].support)
            ion_phase = ion_phase || charges[species] != 0;
        if (ion_phase) {
            recording_support::compare(
                (std::string("charge/") + case_id + "/" + std::to_string(phase)).c_str(),
                "phase_charge", "1", charge, 0.0, phase_inputs,
                "Relationship", "charge-balance", "ionic phase charge is zero",
                "reference_or_one", limit, 0.0, std::abs(charge) <= limit,
                "s1-charge", (std::string(case_id) + "/" + std::to_string(phase)).c_str(),
                "q", "T,P,\\mathbf{x},phase",
                "q=0", "phase");
        }
        for (const auto& reaction : compiled.problem.reactions) {
            double affinity = 0.0;
            for (std::size_t species = 0; species < components; ++species)
                affinity += reaction.stoichiometry[species] * chemical[species];
            recording_support::compare(
                (std::string("affinity/") + case_id + "/" + std::to_string(phase)
                 + "/" + reaction.name).c_str(),
                "reaction_affinity", "J/mol", affinity, 0.0, phase_inputs,
                "Relationship", "reaction-equilibrium", "reaction affinity is zero",
                "reference_or_one", limit, 0.0, std::abs(affinity) <= limit,
                "s1-affinity", (std::string(case_id) + "/" + std::to_string(phase)
                                 + "/" + reaction.name).c_str(),
                "A", "T,P,\\mathbf{x},phase",
                "A=0", "reaction");
        }
    }
    double feed_dual = 0.0;
    for (std::size_t row = 0; row < compiled.rank; ++row)
        feed_dual += result.multipliers[row] * compiled.balance[row];
    recording_support::compare(
        (std::string("feed-dual/") + case_id).c_str(), "objective", "1", result.objective, feed_dual,
        recording_support::fields({recording_support::field("T", result.temperature.front()),
                                   recording_support::field("P", result.pressure.front()),
                                   recording_support::field("coordinates", result.coordinates),
                                   recording_support::field("limit", limit)}),
        "Relationship", "KKT", "objective equals feed balance dual",
        "reference_or_one", limit, 0.0, std::abs(result.objective - feed_dual) <= limit,
        "s1-feed-dual", case_id, "g", "T,P,\\mathbf{z}", "g=\\lambda^Tb", "objective");
}

void check_s1() {
    epcsaft::HelmholtzModel binary_model(properties_support::neutral_binary());
    const auto liquid = equilibrium_support::flash(2);
    epcsaft::core::Capture capture;
    const auto liquid_result = epcsaft::equilibrium::solve_equilibrium(binary_model, liquid);
    verify_s1(binary_model, liquid, liquid_result, "liquid");
    const auto span = std::find_if(capture.trace().spans.begin(), capture.trace().spans.end(),
        [](const auto& item) { return std::string(item.name) == "solve_equilibrium"; });
    CHECK(span != capture.trace().spans.end() && span->convergence.has_value());
    if (span != capture.trace().spans.end() && span->convergence.has_value()) {
        const auto& record = *span->convergence;
        CHECK(record.residual == liquid_result.convergence.residual);
        CHECK(record.criterion == liquid_result.convergence.criterion);
        CHECK(record.floor == liquid_result.convergence.floor);
        CHECK(record.amplification == liquid_result.convergence.amplification);
        CHECK(record.contraction == liquid_result.convergence.contraction);
        CHECK(record.iterations == liquid_result.convergence.iterations);
        CHECK(record.evaluations == liquid_result.convergence.evaluations);
        CHECK(record.stop == liquid_result.convergence.stop);
    }
    auto vapor = equilibrium_support::flash(2);
    vapor.phases[0].kind = epcsaft::equilibrium::PhaseKind::Vapor;
    verify_s1(binary_model, vapor,
              epcsaft::equilibrium::solve_equilibrium(binary_model, vapor), "vapor");
    auto flash = equilibrium_support::flash(3);
    flash.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                    equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    flash.temperature.value = 180.0;
    epcsaft::HelmholtzModel ternary_model(properties_support::neutral_ternary());
    verify_s1(ternary_model, flash,
              epcsaft::equilibrium::solve_equilibrium(ternary_model, flash), "ternary");
    epcsaft::HelmholtzModel reactive_model(equilibrium_support::synthetic_reactive_ternary());
    const auto reactive = equilibrium_support::synthetic_reactive_problem();
    verify_s1(reactive_model, reactive,
              epcsaft::equilibrium::solve_equilibrium(reactive_model, reactive), "reactive");
    auto ethane = equilibrium_support::flash(1);
    ethane.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                     equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    ethane.temperature.value = 250.0;
    ethane.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    epcsaft::HelmholtzModel ethane_model(properties_support::ethane());
    verify_s1(ethane_model, ethane,
              epcsaft::equilibrium::solve_equilibrium(ethane_model, ethane), "ethane");
    auto salt = equilibrium_support::flash(3);
    salt.feed.values = {0.98, 0.01, 0.01};
    salt.temperature.value = 298.15;
    epcsaft::HelmholtzModel salt_model(properties_support::aqueous_electrolyte());
    verify_s1(salt_model, salt,
              epcsaft::equilibrium::solve_equilibrium(salt_model, salt), "salt");
}

const epcsaft::NeutralReferenceReactionDiagnostic::Derivative& temperature_parameter_column(
    const epcsaft::NeutralReferenceReactionDiagnostic& diagnostic, std::size_t coordinate) {
    for (const auto& column : diagnostic.lower_derivatives)
        if (column.temperature_order == 1 && column.pressure_order == 0 && column.parameter == coordinate) return column;
    throw std::logic_error("missing C_Ttheta column");
}

void check_s2() {
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto binary = equilibrium_support::flash(3);
    binary.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                     equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    binary.temperature.value = 180.0;
    binary.pressure.value = 1.0e5;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, binary);
    recording_support::compare(
        "binary-success", "success", "bool", result.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", binary.temperature.value),
                                   recording_support::field("P", binary.pressure.value),
                                   recording_support::field("feed", binary.feed.values)}),
        "Relationship", "stability", "binary solution succeeds", "reference_or_one", 0.0, 0.0,
        result.success, "s2-binary-success", "binary", "success", "T,P,feed",
        "success=1", "case");
    const bool binary_curvature_positive = std::isfinite(result.reduced_hessian_minimum)
        && result.reduced_hessian_minimum > 0.0;
    recording_support::compare_custom(
        "binary-curvature", "reduced_hessian_minimum", "1", result.reduced_hessian_minimum,
        0.0, recording_support::fields({recording_support::field("T", binary.temperature.value),
                                         recording_support::field("P", binary.pressure.value),
                                         recording_support::field("feed", binary.feed.values)}),
        "Relationship", "stability", "reduced Hessian minimum is positive",
        "isfinite(model) && model > reference", binary_curvature_positive,
        "s2-binary-curvature", "binary", "\\lambda_{min}(H)", "T,P,feed",
        "\\lambda_{min}(H)>0", "stability");
    CHECK(result.success && binary_curvature_positive);
    const auto reactive = epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::synthetic_reactive_ternary()),
        equilibrium_support::synthetic_reactive_problem());
    recording_support::compare(
        "reactive-success", "success", "bool", reactive.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", std::vector<double>{0.05, 0.475, 0.475})}),
        "Relationship", "stability", "reactive solution succeeds", "reference_or_one", 0.0, 0.0,
        reactive.success, "s2-reactive-success", "reactive", "success", "T,P,feed",
        "success=1", "case");
    const bool reactive_curvature_positive = std::isfinite(reactive.reduced_hessian_minimum)
        && reactive.reduced_hessian_minimum > 0.0;
    recording_support::compare_custom(
        "reactive-curvature", "reduced_hessian_minimum", "1", reactive.reduced_hessian_minimum,
        0.0, recording_support::fields({recording_support::field("T", 298.15),
                                         recording_support::field("P", 1.0e5),
                                         recording_support::field("feed", std::vector<double>{0.05, 0.475, 0.475})}),
        "Relationship", "stability", "reactive reduced Hessian minimum is positive",
        "isfinite(model) && model > reference", reactive_curvature_positive,
        "s2-reactive-curvature", "reactive", "\\lambda_{min}(H)", "T,P,feed",
        "\\lambda_{min}(H)>0", "stability");
    CHECK(reactive.success && reactive_curvature_positive);

    // #84: the native all-solute path must include the density factor, retain
    // the selected pressure, and reject a charged reaction before any solve.
    epcsaft::HelmholtzModel reference_model(properties_support::neutral_ternary());
    const auto ratios = [](std::vector<double> values) {
        epcsaft::NeutralReferenceDeclaration declaration;
        declaration.solvent_index = 0;
        declaration.standard_molality_mol_per_kg = 1.0;
        declaration.solute_molality_ratios = values;
        return declaration;
    };
    const std::vector<double> nu = {-1.0, 2.0, 0.0};
    const auto reference = epcsaft::ReactionReference{
        std::nullopt, epcsaft::ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    const auto first = epcsaft::reaction_reference(
        reference_model, 300.0, 1.0e5, ratios({0.0, 1.0, 1.0}), nu, reference);
    const auto second = epcsaft::reaction_reference(
        reference_model, 300.0, 1.0e5, ratios({0.0, 2.0, 1.0}), nu, reference);
    CHECK(first.status == epcsaft::NeutralReferenceStatus::Available
          && second.status == epcsaft::NeutralReferenceStatus::Available);
    CHECK(first.terminal_contraction.has_value() && second.terminal_contraction.has_value());
    const double path_gap = std::abs(*first.terminal_contraction - *second.terminal_contraction);
    recording_support::compare_custom(
        "neutral-reference/path-gap", "terminal_contraction_gap", "1", path_gap, 5.0e-5,
        recording_support::fields({recording_support::field("T", 300.0),
            recording_support::field("P", 1.0e5), recording_support::field("nu", nu)}),
        "Relationship", "independent-neutral-paths",
        "independent positive neutral solute paths agree at the finite terminal point",
        "isfinite(model) && model <= reference", path_gap <= 5.0e-5,
        "s2-neutral-reference", "path-gap", "|C_a-C_b|", "T,P,nu", "|C_a-C_b|<=5e-5", "path");
    const double molality = 1.0e-14;
    const double solvent_moles = 1.0 / reference_model.parameters().molar_mass_kg_per_mol[0];
    const double total = solvent_moles + 2.0 * molality;
    const std::vector<double> terminal_x = {
        solvent_moles / total, molality / total, molality / total};
    const auto terminal = reference_model.state_at_pressure(
        300.0, 1.0e5, terminal_x, epcsaft::DensityBranch::Liquid);
    CHECK(terminal.log_fugacity_coefficient.has_value());
    double identity = 0.0;
    for (std::size_t i = 0; i < nu.size(); ++i)
        identity += nu[i] * ((*terminal.log_fugacity_coefficient)[i].value()
            + std::log(1.0e5 / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0)));
    recording_support::compare(
        "neutral-reference/lambda-identity", "terminal_contraction", "1",
        *first.terminal_contraction, identity,
        recording_support::fields({recording_support::field("rho", terminal.molar_density),
            recording_support::field("observed_terminal_change", *first.observed_terminal_change)}),
        "Relationship", "fugacity-pressure-identity",
        "Lambda includes log(P/(RT rho0)) in the independent fugacity form",
        "reference_or_one", 0.0, 1.0e-9,
        close_with_tolerance(*first.terminal_contraction, identity, 1.0e-9),
        "s2-neutral-reference", "lambda-identity", "C", "rho,phi,T,P", "C=nu.Lambda", "identity");
    const auto fixed_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5, ratios({0.0, 1.0, 1.0}), nu,
        epcsaft::ReactionReference{1.0e5, reference.source_basis});
    const auto system_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5, ratios({0.0, 1.0, 1.0}), nu, reference);
    CHECK(fixed_pressure.status == epcsaft::NeutralReferenceStatus::Available
          && system_pressure.status == epcsaft::NeutralReferenceStatus::Available
          && fixed_pressure.reference_pressure_pa == 1.0e5
          && fixed_pressure.terminal_contraction.has_value()
          && system_pressure.terminal_contraction.has_value());
    CHECK(std::abs(*fixed_pressure.terminal_contraction - *first.terminal_contraction) <= 1.0e-12);
    const double pressure_gap = std::abs(
        *system_pressure.terminal_contraction - *fixed_pressure.terminal_contraction);
    recording_support::compare_custom(
        "neutral-reference/pressure-path", "pressure_contraction_gap", "1", pressure_gap, 1.0e-6,
        recording_support::fields({recording_support::field("system_pressure", 2.0e5),
            recording_support::field("reference_pressure", 1.0e5)}),
        "Relationship", "source-pressure-selection",
        "system and explicit source pressures produce distinct Delta-nu contraction values",
        "isfinite(model) && model > reference", pressure_gap > 1.0e-6,
        "s2-neutral-reference", "pressure-path", "|C_P-C_P0|", "T,P,P0", "|gap|>1e-6", "pressure");
    CHECK(system_pressure.terminal_contraction_pressure_derivative.has_value()
          && system_pressure.terminal_contraction_pressure_second_derivative.has_value()
          && fixed_pressure.terminal_contraction_pressure_derivative.has_value()
          && fixed_pressure.terminal_contraction_pressure_second_derivative.has_value()
          && *fixed_pressure.terminal_contraction_pressure_derivative == 0.0
          && *fixed_pressure.terminal_contraction_pressure_second_derivative == 0.0);
    const double pressure_step = 100.0;
    const auto lower_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5 - pressure_step, ratios({0.0, 1.0, 1.0}), nu, reference);
    const auto upper_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5 + pressure_step, ratios({0.0, 1.0, 1.0}), nu, reference);
    CHECK(lower_pressure.terminal_contraction.has_value()
          && upper_pressure.terminal_contraction.has_value());
    const double centered_pressure_derivative =
        (*upper_pressure.terminal_contraction - *lower_pressure.terminal_contraction)
        / (2.0 * pressure_step);
    const double centered_pressure_second_derivative =
        (*upper_pressure.terminal_contraction
            - 2.0 * *system_pressure.terminal_contraction
            + *lower_pressure.terminal_contraction)
        / (pressure_step * pressure_step);
    const auto jet_close = [](double left, double right) {
        return std::abs(left - right) <= 1.0e-3
            * std::max({std::abs(left), std::abs(right), 1.0e-12});
    };
    CHECK(jet_close(*system_pressure.terminal_contraction_pressure_derivative,
                    centered_pressure_derivative)
          && jet_close(*system_pressure.terminal_contraction_pressure_second_derivative,
                       centered_pressure_second_derivative));

    // Identical EOS species provide a finite nonideal counterproof: the density
    // term and residual chemical potential do not separate into an ideal -1/P law.
    auto identical_parameters = properties_support::parameters(2);
    identical_parameters.segment_count = {1.6069, 1.6069};
    identical_parameters.sigma_constant = {3.5206, 3.5206};
    identical_parameters.dispersion_energy_over_k = {191.42, 191.42};
    identical_parameters.molar_mass_kg_per_mol = {0.03007, 0.03007};
    epcsaft::HelmholtzModel identical_model(identical_parameters);
    const auto identical_declaration = ratios({0.0, 1.0});
    const std::vector<double> identical_nu = {-2.0, 1.0};
    const auto identical_reference = epcsaft::reaction_reference(
        identical_model, 300.0, 1.0e5, identical_declaration, identical_nu, reference);
    CHECK(identical_reference.status == epcsaft::NeutralReferenceStatus::Available
          && identical_reference.terminal_contraction_pressure_derivative.has_value()
          && identical_reference.terminal_contraction_pressure_second_derivative.has_value());
    const double identical_molality = 1.0e-14;
    const double identical_solvent_moles = 1.0 / identical_parameters.molar_mass_kg_per_mol[0];
    const double identical_total = identical_solvent_moles + identical_molality;
    const auto identical_terminal = identical_model.state_at_pressure(
        300.0, 1.0e5,
        {identical_solvent_moles / identical_total, identical_molality / identical_total},
        epcsaft::DensityBranch::Liquid);
    const double identical_cp = -1.0 / (
        identical_terminal.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0);
    const double identical_cpp = 1.0 / (
        identical_terminal.molar_density * identical_terminal.molar_density
        * epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0
        * identical_terminal.pressure_density_derivative);
    CHECK(jet_close(*identical_reference.terminal_contraction_pressure_derivative, identical_cp)
          && jet_close(*identical_reference.terminal_contraction_pressure_second_derivative,
                       identical_cpp));

    // #84 temperature slope, common native/Python case: aqueous NaCl with the formula
    // reaction H2O = Na+ + Cl- at the system pressure and at a fixed 2e5 Pa source.
    // dC/dT|_P against a centered same-branch terminal difference; floor: truncation
    // step^2 C^(3)/6 ~ (0.01/298)^2 |C'| and value roundoff/step are below 1e-7 relative.
    // The measured multiplier M = T|dD|/|dC| qualifies the observed slope criterion.
    {
        epcsaft::HelmholtzModel ionic(properties_support::aqueous_electrolyte());
        const auto declaration = ratios({0.0, 1.0, 1.0});
        const std::vector<std::vector<double>> rows = {{-1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}};
        const std::vector<epcsaft::ReactionReference> sources = {
            reference, epcsaft::ReactionReference{2.0e5, reference.source_basis}};
        const double temperature = 298.15, pressure = 1.0e5, step = 1.0e-2;
        const auto center = epcsaft::reaction_references(
            ionic, temperature, pressure, declaration, rows, sources, false, std::nullopt, 2);
        const auto lower = epcsaft::reaction_references(
            ionic, temperature - step, pressure, declaration, rows, sources, false, std::nullopt, 1);
        const auto upper = epcsaft::reaction_references(
            ionic, temperature + step, pressure, declaration, rows, sources, false, std::nullopt, 1);
        // #138 calorics: water carries a manufactured ideal-gas record (cp 33.59 J/(mol K),
        // h = -241826 J/mol at 298.15 K); Na+ and Cl- have none and are reaction-completed.
        epcsaft::ThermochemistryRecord water_record;
        water_record.standard_pressure = 1.0e5;
        water_record.component_intervals = {{{200.0, 500.0, true, true, {epcsaft::IdealPolynomial{
            {33.59}, 298.15, -241826.0, 188.83}, 1.0e5}}}, {}, {}};
        epcsaft::HelmholtzModel ionic_thermal(properties_support::aqueous_electrolyte(), {}, water_record);
        for (std::size_t row = 0; row < rows.size(); ++row) {
            CHECK(center[row].terminal_contraction_temperature_derivative.has_value()
                  && center[row].observed_temperature_derivative_change.has_value()
                  && center[row].observed_terminal_change.has_value()
                  && lower[row].terminal_contraction.has_value()
                  && upper[row].terminal_contraction.has_value());
            const double slope = *center[row].terminal_contraction_temperature_derivative;
            const double centered = (*upper[row].terminal_contraction
                - *lower[row].terminal_contraction) / (2.0 * step);
            const std::string label = row == 0 ? "system-pressure" : "fixed-pressure";
            const auto inputs = recording_support::fields({
                recording_support::field("T", temperature), recording_support::field("P", pressure),
                recording_support::field("P_ref", sources[row].reference_pressure_pa.value_or(pressure)),
                recording_support::field("step", step)});
            recording_support::compare(
                ("neutral-reference/temperature-slope/" + label).c_str(), "dC_dT", "1/K",
                slope, centered, inputs, "Relationship", "centered-terminal-difference",
                "dC/dT|_P carries the implicit density and association chains of the reference",
                "reference_scale", 0.0, 1.0e-6,
                std::abs(slope - centered) <= 1.0e-6 * std::abs(centered),
                "s2-neutral-reference", "temperature-slope", "dC/dT", "T,P,nu", "rel<=1e-6", "slope");
            const double change = *center[row].observed_temperature_derivative_change;
            recording_support::compare_custom(
                ("neutral-reference/temperature-limit/" + label).c_str(),
                "T_dC_dT_terminal_change", "1", change, 5.0e-5, inputs,
                "Relationship", "observed-slope-convergence",
                "the slope's terminal ladder change stays within the value criterion in RT units",
                "isfinite(model) && model <= reference", change <= 5.0e-5,
                "s2-neutral-reference", "temperature-limit", "T|dD|", "T,P,nu", "T|dD|<=5e-5", "limit");
            // The same manufactured formula-reaction solve as the Python test; its first-order
            // T and T/P actions of x_Na are the rows that the installed wheel must reproduce.
            auto problem = equilibrium_support::flash(3);
            problem.temperature = epcsaft::equilibrium::FreeValue::fixed(temperature);
            problem.pressure = epcsaft::equilibrium::FreeValue::fixed(pressure);
            problem.feed.values = {1.0, 1.0e-6, 1.0e-6};
            problem.source_reference = declaration;
            problem.reactions = {{"H2O=Na++Cl-", rows[row], std::nullopt, std::nullopt,
                epcsaft::ReactionLogPolynomial{std::log(1.0e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0,
                                               "analytic-source"}, sources[row]}};
            const auto solved = epcsaft::equilibrium::solve_equilibrium(ionic, problem);
            CHECK(solved.success);
            const auto compiled = epcsaft::equilibrium::compile_problem(ionic, problem);
            epcsaft::equilibrium::SolvedStateObservable sodium;
            sodium.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
            sodium.component = 1;
            for (const auto& [name, direction] : std::vector<std::pair<std::string,
                     epcsaft::equilibrium::SolvedStateActionDirection>>{
                     {"T", {1.0, 0.0, {0.0, 0.0, 0.0}, {}}}, {"T+P", {1.0, 1.0e5, {0.0, 0.0, 0.0}, {}}}}) {
                const auto action = epcsaft::equilibrium::solved_state_actions(
                    compiled, solved, {{sodium}, {direction}}).results[0];
                CHECK(action.status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                      && action.action.has_value());
                recording_support::observe(
                    ("reference-temperature/ionic/" + std::string(row == 0 ? "system" : "fixed")
                     + "/dx_Na/" + name).c_str(), "dx_Na", "1/K", action.action.value_or(0.0),
                    std::nullopt, inputs, "Relationship", "solved-state-first-action",
                    "native first-order referenced T/P action for the installed-wheel comparison");
            }
            recording_support::observe(
                ("neutral-reference/temperature-multiplier/" + label).c_str(), "M", "1",
                change / *center[row].observed_terminal_change, std::nullopt, inputs,
                "Relationship", "debye-huckel-sqrt-m-form",
                "measured slope-to-value terminal-change multiplier for this model and state");
            // Order two on the ionic ladder (Born, density and association third-derivative chains;
            // the Debye-Hueckel term vanishes at 1e-14 mol/kg and is not tested here) against the
            // centered difference of the order-one slope; floor: step^2 |d4C/dT4|/6 and slope
            // roundoff/step, below 1e-7 relative. Admission is sized to cp: R T^2 max|dD2|.
            CHECK(center[row].terminal_contraction_temperature_second_derivative.has_value()
                  && lower[row].terminal_contraction_temperature_derivative.has_value()
                  && upper[row].terminal_contraction_temperature_derivative.has_value());
            const double second = center[row].terminal_contraction_temperature_second_derivative.value_or(0.0);
            const double second_centered = (upper[row].terminal_contraction_temperature_derivative.value_or(0.0)
                - lower[row].terminal_contraction_temperature_derivative.value_or(0.0)) / (2.0 * step);
            recording_support::compare(
                ("neutral-reference/temperature-second-slope/" + label).c_str(), "d2C_dT2", "1/K^2",
                second, second_centered, inputs, "Relationship", "centered-slope-difference",
                "d2C/dT2|_P carries the Born, density and association chains through E3",
                "reference_scale", 0.0, 1.0e-6,
                std::abs(second - second_centered) <= 1.0e-6 * std::abs(second_centered),
                "s2-neutral-reference", "temperature-second-slope", "d2C/dT2", "T,P,nu", "rel<=1e-6", "slope");
            const double second_change = epcsaft::GAS_CONSTANT_J_PER_MOL_K
                * center[row].observed_temperature_second_derivative_change.value_or(1.0);
            recording_support::compare_custom(
                ("neutral-reference/temperature-second-limit/" + label).c_str(),
                "R_T2_d2C_dT2_terminal_change", "J/(mol K)", second_change, 1.0e-2, inputs,
                "Relationship", "observed-slope-convergence",
                "the second slope's terminal ladder change is negligible at the heat-capacity scale",
                "isfinite(model) && model <= reference", second_change <= 1.0e-2,
                "s2-neutral-reference", "temperature-second-limit", "R T^2|dD2|", "T,P,nu",
                "R T^2|dD2|<=0.01 J/(mol K)", "limit");
            const auto thermal_solved = epcsaft::equilibrium::solve_equilibrium(ionic_thermal, problem);
            CHECK(thermal_solved.success);
            epcsaft::equilibrium::SolvedStateObservable enthalpy, capacity;
            enthalpy.kind = capacity.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseProperty;
            enthalpy.property = epcsaft::PropertyObservable::TotalEnthalpy;
            capacity.property = epcsaft::PropertyObservable::TotalIsobaricHeatCapacity;
            const auto calorics = epcsaft::equilibrium::solved_state_actions(
                epcsaft::equilibrium::compile_problem(ionic_thermal, problem), thermal_solved,
                {{enthalpy, capacity}, {{0.0, 0.0, {0.0, 0.0, 0.0}, {}}}}).results;
            for (std::size_t i = 0; i < 2; ++i) {
                CHECK(calorics[i].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                      && calorics[i].value.has_value());
                recording_support::observe(
                    ("reference-calorics/ionic/" + std::string(row == 0 ? "system" : "fixed")
                     + (i == 0 ? "/h" : "/cp")).c_str(), i == 0 ? "h" : "cp", i == 0 ? "J/mol" : "J/(mol K)",
                    calorics[i].value.value_or(0.0), std::nullopt, inputs, "Relationship", "reference-species-calorics",
                    "native total caloric value with reaction-completed ions for the installed-wheel comparison");
            }
        }
    }
    // #333: at solvent permittivity 30 the 1e-10 to 1e-12 mol/kg NaCl contraction change
    // exceeds 5e-5; the change into the returned 1e-14 value admits it, as the ion witness does.
    {
        auto parameters = properties_support::aqueous_electrolyte();
        parameters.relative_permittivity[0] = 30.0;
        const auto low = epcsaft::reaction_reference(epcsaft::HelmholtzModel(parameters),
            298.15, 1.0e5, ratios({0.0, 1.0, 1.0}), {-1.0, 1.0, 1.0}, reference);
        CHECK(low.status == epcsaft::NeutralReferenceStatus::Available && low.terminal_contraction
              && *low.observed_terminal_change <= 5.0e-5);
    }
    // #138 order-two slope on the associating MEA-water fixture, solvent water, formula
    // reaction H2O = MEA: d2C/dT2|_P against the centered difference of dC/dT|_P; floor:
    // truncation step^2 |d4C/dT4|/6 and slope roundoff/step, both below 1e-7 relative.
    {
        epcsaft::HelmholtzModel associating(properties_support::associating_binary());
        epcsaft::NeutralReferenceDeclaration water;
        water.solvent_index = 1;
        water.solute_molality_ratios = {1.0, 0.0};
        water.standard_molality_mol_per_kg = 1.0;
        const std::vector<std::vector<double>> rows = {{1.0, -1.0}};
        const std::vector<epcsaft::ReactionReference> sources = {reference};
        const double temperature = 313.15, pressure = 1.0e5, step = 1.0e-2;
        const auto center = epcsaft::reaction_references(
            associating, temperature, pressure, water, rows, sources, false, std::nullopt, 2);
        const auto lower = epcsaft::reaction_references(
            associating, temperature - step, pressure, water, rows, sources, false, std::nullopt, 1);
        const auto upper = epcsaft::reaction_references(
            associating, temperature + step, pressure, water, rows, sources, false, std::nullopt, 1);
        CHECK(center[0].terminal_contraction_temperature_second_derivative.has_value()
              && lower[0].terminal_contraction_temperature_derivative.has_value()
              && upper[0].terminal_contraction_temperature_derivative.has_value());
        const double second = center[0].terminal_contraction_temperature_second_derivative.value_or(0.0);
        const double centered = (upper[0].terminal_contraction_temperature_derivative.value_or(0.0)
            - lower[0].terminal_contraction_temperature_derivative.value_or(0.0)) / (2.0 * step);
        const auto inputs = recording_support::fields({
            recording_support::field("T", temperature), recording_support::field("P", pressure),
            recording_support::field("step", step)});
        recording_support::compare(
            "neutral-reference/temperature-second-slope/associating", "d2C_dT2", "1/K^2",
            second, centered, inputs, "Relationship", "centered-slope-difference",
            "d2C/dT2|_P carries the implicit density and association chains through E3",
            "reference_scale", 0.0, 1.0e-6, std::abs(second - centered) <= 1.0e-6 * std::abs(centered),
            "s2-neutral-reference", "temperature-second-slope", "d2C/dT2", "T,P,nu", "rel<=1e-6", "slope");
        const double change = epcsaft::GAS_CONSTANT_J_PER_MOL_K
            * center[0].observed_temperature_second_derivative_change.value_or(1.0);
        recording_support::compare_custom(
            "neutral-reference/temperature-second-limit/associating", "R_T2_d2C_dT2_terminal_change", "J/(mol K)",
            change, 1.0e-2, inputs, "Relationship", "observed-slope-convergence",
            "the second slope's terminal ladder change is negligible at the heat-capacity scale",
            "isfinite(model) && model <= reference", change <= 1.0e-2,
            "s2-neutral-reference", "temperature-second-limit", "R T^2|dD2|", "T,P,nu",
            "R T^2|dD2|<=0.01 J/(mol K)", "limit");
    }

    // #210 F1(i): identical ethane keeps the reference species identical, so the
    // mixed reference slope is the pure-fluid residual enthalpy parameter action.
    {
        constexpr double temperature = 260.0, pressure = 3.0e6;
        const std::vector<epcsaft::ActiveParameter> active = {
            {epcsaft::ParameterFamily::DispersionEnergyOverK, {0}},
            {epcsaft::ParameterFamily::DispersionEnergyOverK, {1}}};
        const epcsaft::HelmholtzModel pair(identical_parameters, active);
        const epcsaft::HelmholtzModel pure(
            properties_support::ethane(),
            {{epcsaft::ParameterFamily::DispersionEnergyOverK, {0}}});
        const std::vector<epcsaft::ReactionReference> sources = {
            reference, {5.0e6, reference.source_basis}};
        for (std::size_t source = 0; source < sources.size(); ++source) {
            const auto diagnostic = epcsaft::reaction_references(
                pair, temperature, pressure, identical_declaration, {identical_nu},
                {sources[source]}, false, std::nullopt, 1, true).front();
            const double mixed = diagnostic.terminal_contraction_temperature_parameter_derivative.size() != 2
                ? std::numeric_limits<double>::quiet_NaN()
                : diagnostic.terminal_contraction_temperature_parameter_derivative[0].value_or(
                    std::numeric_limits<double>::quiet_NaN())
                    + diagnostic.terminal_contraction_temperature_parameter_derivative[1].value_or(
                        std::numeric_limits<double>::quiet_NaN());
            const double mixed_s0 = diagnostic.terminal_contraction_temperature_parameter_derivative.size() > 0
                ? diagnostic.terminal_contraction_temperature_parameter_derivative[0].value_or(
                    std::numeric_limits<double>::quiet_NaN())
                : std::numeric_limits<double>::quiet_NaN();
            const double mixed_s1 = diagnostic.terminal_contraction_temperature_parameter_derivative.size() > 1
                ? diagnostic.terminal_contraction_temperature_parameter_derivative[1].value_or(
                    std::numeric_limits<double>::quiet_NaN())
                : std::numeric_limits<double>::quiet_NaN();
            const double reference_pressure = sources[source].reference_pressure_pa.value_or(pressure);
            const auto pure_action = pure.pressure_actions(
                temperature, reference_pressure, {1.0}, epcsaft::DensityBranch::Liquid,
                {epcsaft::PropertyObservable::ResidualEnthalpy},
                {{0.0, 0.0, {0.0}, {{epcsaft::ParameterFamily::DispersionEnergyOverK, {0}, 1.0}}}});
            CHECK(diagnostic.status == epcsaft::NeutralReferenceStatus::Available
                  && std::isfinite(mixed) && pure_action.results.size() == 1
                  && pure_action.results[0].status == epcsaft::PropertyActionStatus::Available
                  && pure_action.results[0].action.has_value());
            const double scaled_mixed = epcsaft::GAS_CONSTANT_J_PER_MOL_K
                * temperature * temperature * mixed;
            const double pure_slope = pure_action.results.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : pure_action.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN());
            const bool agrees = std::isfinite(scaled_mixed) && std::isfinite(pure_slope)
                && std::abs(scaled_mixed - pure_slope) <= 1.0e-9 * std::abs(pure_slope);
            const auto inputs = recording_support::fields({
                recording_support::field("case", "N"),
                recording_support::field("source", source == 0 ? "system" : "fixed"),
                recording_support::field("T_K", temperature),
                recording_support::field("P_Pa", pressure),
                recording_support::field("P_ref_Pa", reference_pressure),
                recording_support::field("epsilon_over_k_K", 191.42),
                recording_support::field("C_Ttheta_s0", mixed_s0),
                recording_support::field("C_Ttheta_s1", mixed_s1),
                recording_support::field("C_Ttheta", mixed),
                recording_support::field("C_Ttheta_ladder_s0", std::vector<double>{
                    temperature_parameter_column(diagnostic, 0).ladder[0].value_or(std::nan("")),
                    temperature_parameter_column(diagnostic, 0).ladder[1].value_or(std::nan("")),
                    temperature_parameter_column(diagnostic, 0).ladder[2].value_or(std::nan(""))}),
                recording_support::field("C_Ttheta_ladder_s1", std::vector<double>{
                    temperature_parameter_column(diagnostic, 1).ladder[0].value_or(std::nan("")),
                    temperature_parameter_column(diagnostic, 1).ladder[1].value_or(std::nan("")),
                    temperature_parameter_column(diagnostic, 1).ladder[2].value_or(std::nan(""))}),
                recording_support::field("observed_C_Ttheta_admitted_change_s0", temperature_parameter_column(diagnostic, 0).scaled_changes[1]),
                recording_support::field("observed_C_Ttheta_admitted_change_s1", temperature_parameter_column(diagnostic, 1).scaled_changes[1]),
                recording_support::field("R_T2_C_Ttheta_J_per_mol_K", scaled_mixed),
                recording_support::field("pure_pressure_action_J_per_mol_K", pure_slope),
                recording_support::field("observed_C_T_admitted_change",
                    epcsaft::core::optional_number(diagnostic.observed_temperature_derivative_change)),
                recording_support::field("observed_terminal_change",
                    epcsaft::core::optional_number(diagnostic.observed_terminal_change))});
            const std::string id = std::string("neutral-reference/temperature-parameter/N/")
                + (source == 0 ? "system" : "fixed") + "/tied-epsilon/pure-identity";
            recording_support::compare_custom(
                id.c_str(), "R_T2_C_Ttheta", "J/(mol K)", scaled_mixed, pure_slope, inputs,
                "Relationship", "identical-species-gibbs-helmholtz",
                "the reference enthalpy parameter slope equals the pure-fluid pressure-chart action",
                "abs(model-reference)<=1e-9*abs(reference)", agrees,
                "s2-neutral-reference", "f1-N", "R T^2 C_{Ttheta}", "T,P,P_ref,epsilon",
                "R T^2 C_{Ttheta}=d h_res/d epsilon", "source");
        }
    }

    // #210 F1(ii): on aqueous NaCl, qualify both mixed-slope directions against
    // independent centered theta and T differences at h and 2h.
    {
        constexpr double temperature = 298.15, pressure = 1.0e5, temperature_step = 1.0e-2;
        const auto parameters = properties_support::aqueous_electrolyte();
        const epcsaft::NeutralReferenceDeclaration declaration = ratios({0.0, 1.0, 1.0});
        const std::vector<double> stoichiometry = {-1.0, 1.0, 1.0};
        const std::vector<epcsaft::ActiveParameter> active = {
            {epcsaft::ParameterFamily::DispersionEnergyOverK, {0}},
            {epcsaft::ParameterFamily::SigmaConstant, {1}}};
        const epcsaft::HelmholtzModel model(parameters, active);
        const std::vector<epcsaft::ReactionReference> sources = {
            reference, {2.0e5, reference.source_basis}};
        const double steps[] = {0.05, 5.0e-4};
        const char* names[] = {"water-epsilon", "sodium-sigma"};
        const char* units[] = {"K", "angstrom"};
        const auto shifted_parameters = [&](std::size_t coordinate, double shift) {
            auto shifted = properties_support::aqueous_electrolyte();
            if (coordinate == 0) shifted.dispersion_energy_over_k[0] += shift;
            else shifted.sigma_constant[1] += shift;
            return shifted;
        };
        const auto temperature_slope = [&](const epcsaft::HelmholtzModel& current,
                                           double at_temperature,
                                           const epcsaft::ReactionReference& source) {
            const auto diagnostic = epcsaft::reaction_references(
                current, at_temperature, pressure, declaration, {stoichiometry}, {source},
                false, std::nullopt, 1, false).front();
            CHECK(diagnostic.status == epcsaft::NeutralReferenceStatus::Available
                  && diagnostic.terminal_contraction_temperature_derivative.has_value());
            return diagnostic.terminal_contraction_temperature_derivative.value_or(
                std::numeric_limits<double>::quiet_NaN());
        };
        const auto parameter_slope = [&](double at_temperature,
                                         const epcsaft::ReactionReference& source,
                                         std::size_t coordinate) {
            const auto diagnostic = epcsaft::reaction_references(
                model, at_temperature, pressure, declaration, {stoichiometry}, {source},
                false, std::nullopt, 0, true).front();
            CHECK(diagnostic.status == epcsaft::NeutralReferenceStatus::Available
                  && diagnostic.terminal_contraction_parameter_derivative.size() == active.size()
                  && diagnostic.terminal_contraction_parameter_derivative[coordinate].has_value());
            return diagnostic.terminal_contraction_parameter_derivative.size() > coordinate
                ? diagnostic.terminal_contraction_parameter_derivative[coordinate].value_or(
                    std::numeric_limits<double>::quiet_NaN())
                : std::numeric_limits<double>::quiet_NaN();
        };
        for (std::size_t source = 0; source < sources.size(); ++source) {
            const auto exact_diagnostic = epcsaft::reaction_references(
                model, temperature, pressure, declaration, {stoichiometry}, {sources[source]},
                false, std::nullopt, 1, true).front();
            CHECK(exact_diagnostic.status == epcsaft::NeutralReferenceStatus::Available
                  && exact_diagnostic.terminal_contraction_temperature_parameter_derivative.size() == active.size());
            for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
                const double exact = exact_diagnostic.terminal_contraction_temperature_parameter_derivative.size() > coordinate
                    ? exact_diagnostic.terminal_contraction_temperature_parameter_derivative[coordinate].value_or(
                        std::numeric_limits<double>::quiet_NaN())
                    : std::numeric_limits<double>::quiet_NaN();
                const auto theta_difference = [&](double step) {
                    const epcsaft::HelmholtzModel plus(shifted_parameters(coordinate, step));
                    const epcsaft::HelmholtzModel minus(shifted_parameters(coordinate, -step));
                    return (temperature_slope(plus, temperature, sources[source])
                        - temperature_slope(minus, temperature, sources[source])) / (2.0 * step);
                };
                const auto temperature_difference = [&](double step) {
                    return (parameter_slope(temperature + step, sources[source], coordinate)
                        - parameter_slope(temperature - step, sources[source], coordinate)) / (2.0 * step);
                };
                const auto inputs = recording_support::fields({
                    recording_support::field("case", "I"),
                    recording_support::field("source", source == 0 ? "system" : "fixed"),
                    recording_support::field("parameter", names[coordinate]),
                    recording_support::field("parameter_unit", units[coordinate]),
                    recording_support::field("T_K", temperature),
                    recording_support::field("P_Pa", pressure),
                    recording_support::field("P_ref_Pa", sources[source].reference_pressure_pa.value_or(pressure)),
                    recording_support::field("C_Ttheta_ladder", std::vector<double>{
                        temperature_parameter_column(exact_diagnostic, coordinate).ladder[0].value_or(std::numeric_limits<double>::quiet_NaN()),
                        temperature_parameter_column(exact_diagnostic, coordinate).ladder[1].value_or(std::numeric_limits<double>::quiet_NaN()),
                        temperature_parameter_column(exact_diagnostic, coordinate).ladder[2].value_or(std::numeric_limits<double>::quiet_NaN())}),
                    recording_support::field("observed_C_Ttheta_admitted_change",
                        temperature_parameter_column(exact_diagnostic, coordinate).scaled_changes[1]),
                    recording_support::field("observed_C_T_admitted_change",
                        epcsaft::core::optional_number(exact_diagnostic.observed_temperature_derivative_change)),
                    recording_support::field("observed_terminal_change",
                        epcsaft::core::optional_number(exact_diagnostic.observed_terminal_change))});
                for (const auto& [route, coarse, fine] : std::array<std::tuple<const char*, double, double>, 2>{
                         std::tuple{"theta-difference", theta_difference(2.0 * steps[coordinate]),
                                    theta_difference(steps[coordinate])},
                         std::tuple{"temperature-difference", temperature_difference(2.0 * temperature_step),
                                    temperature_difference(temperature_step)}}) {
                    const double fine_error = std::abs(fine - exact);
                    const double coarse_error = std::abs(coarse - exact);
                    const bool agrees = std::isfinite(exact) && std::isfinite(coarse) && std::isfinite(fine)
                        && std::abs(exact) > 0.0 && fine_error <= 1.0e-6 * std::abs(exact)
                        && fine_error <= coarse_error / 2.0;
                    const std::string route_id = std::string("neutral-reference/temperature-parameter/I/")
                        + (source == 0 ? "system/" : "fixed/") + names[coordinate] + "/" + route;
                    const auto route_inputs = recording_support::fields({inputs,
                        recording_support::field("coarse_step", std::string(route) == "theta-difference"
                            ? 2.0 * steps[coordinate] : 2.0 * temperature_step),
                        recording_support::field("fine_step", std::string(route) == "theta-difference"
                            ? steps[coordinate] : temperature_step),
                        recording_support::field("coarse", coarse), recording_support::field("fine", fine),
                        recording_support::field("exact_C_Ttheta", exact),
                        recording_support::field("coarse_error", coarse_error),
                        recording_support::field("fine_error", fine_error)});
                    recording_support::compare_custom(
                        route_id.c_str(), "C_Ttheta", coordinate == 0 ? "K^-2" : "K^-1 angstrom^-1",
                        fine, exact, route_inputs, "Relationship", "centered-independent-difference",
                        "the exact mixed reference slope agrees with independent theta and T differences",
                        "abs(fine-exact)<=1e-6*abs(exact) && abs(fine-exact)<=abs(coarse-exact)/2",
                        agrees, "s2-neutral-reference", "f1-I", "C_{Ttheta}", "T,P,P_ref,theta,h",
                        "O(h^2) centered difference converges to C_{Ttheta}", route);
                }
            }
        }
    }

    // #210 F2: extensive nH actions against centered re-solves. The omission
    // diagnostic is the independently assembled direct reference-species term.
    {
        using ActionDirection = epcsaft::equilibrium::SolvedStateActionDirection;
        struct Ladder final { double action, coarse, fine, richardson, direct; };
        const auto enthalpy = [] {
            epcsaft::equilibrium::SolvedStateObservable item;
            item.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseProperty;
            item.property = epcsaft::PropertyObservable::TotalEnthalpy;
            return item;
        }();
        const auto amount = [] {
            epcsaft::equilibrium::SolvedStateObservable item;
            item.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseAmount;
            return item;
        }();
        const auto extensive = [&](const epcsaft::HelmholtzModel& model,
                                   const epcsaft::equilibrium::EquilibriumProblem& problem) {
            const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem);
            CHECK(solved.success);
            if (!solved.success) return std::numeric_limits<double>::quiet_NaN();
            const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
            ActionDirection zero;
            zero.feed_amounts_mol.assign(problem.feed.values.size(), 0.0);
            zero.active_parameters.assign(model.active_parameters().size(), 0.0);
            const auto values = epcsaft::equilibrium::solved_state_actions(
                compiled, solved, {{enthalpy, amount}, {zero}}).results;
            CHECK(values.size() == 2 && values[0].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                  && values[1].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                  && values[0].value.has_value() && values[1].value.has_value());
            if (values.size() != 2 || !values[0].value || !values[1].value)
                return std::numeric_limits<double>::quiet_NaN();
            return *values[0].value * *values[1].value;
        };
        const auto action_ladder = [&](const char* id, const epcsaft::HelmholtzModel& model,
                                       const epcsaft::equilibrium::EquilibriumProblem& problem,
                                       const ActionDirection& direction, double step,
                                       const auto& value_at_shift) {
            const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem);
            CHECK(solved.success);
            const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
            const auto batches = epcsaft::equilibrium::solved_state_first_actions(
                compiled, solved, {enthalpy, amount}, {direction});
            CHECK(solved.success && batches.size() == 1 && batches[0].results.size() == 2);
            const auto& h = batches[0].results[0];
            const auto& n = batches[0].results[1];
            CHECK(h.status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                  && n.status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                  && h.value.has_value() && h.action.has_value() && n.value.has_value() && n.action.has_value());
            const double exact = n.value.value_or(std::numeric_limits<double>::quiet_NaN())
                    * h.action.value_or(std::numeric_limits<double>::quiet_NaN())
                + h.value.value_or(std::numeric_limits<double>::quiet_NaN())
                    * n.action.value_or(std::numeric_limits<double>::quiet_NaN());
            const double coarse = (value_at_shift(2.0 * step) - value_at_shift(-2.0 * step)) / (4.0 * step);
            const double fine = (value_at_shift(step) - value_at_shift(-step)) / (2.0 * step);
            const double richardson = (4.0 * fine - coarse) / 3.0;
            const auto calorics = epcsaft::equilibrium::reference_calorics(
                compiled, problem.temperature.value, problem.pressure.value, false, false, true);
            CHECK(calorics.status == epcsaft::equilibrium::SolvedStateActionStatus::Available);
            double direct = 0.0;
            for (std::size_t component = 0; component < solved.mole_fractions.size(); ++component) {
                double slope = 0.0;
                for (std::size_t a = 0; a < direction.active_parameters.size(); ++a)
                    if (direction.active_parameters[a] != 0.0)
                        slope += direction.active_parameters[a] * (
                            reactive_caloric_checks::gibbs_contraction(calorics, component, 0, 0, a) -
                            problem.temperature.value * reactive_caloric_checks::gibbs_contraction(calorics, component, 1, 0, a));
                if (!direction.reaction_coefficients.empty())
                    slope -= epcsaft::GAS_CONSTANT_J_PER_MOL_K * direction.reaction_coefficients.front()[1]
                        * calorics.reaction_map[0][component];
                direct += solved.amounts[0] * solved.mole_fractions[component] * slope;
            }
            const double tolerance = 1.0e-6 * std::abs(richardson);
            const bool agrees = std::isfinite(exact) && std::isfinite(coarse) && std::isfinite(fine)
                && std::isfinite(richardson) && std::isfinite(direct) && tolerance > 0.0
                && std::abs(richardson - exact) <= tolerance
                && std::abs(fine - exact) <= std::abs(coarse - exact) / 2.0;
            const auto inputs = recording_support::fields({
                recording_support::field("coarse_step", 2.0 * step),
                recording_support::field("fine_step", step),
                recording_support::field("coarse", coarse), recording_support::field("fine", fine),
                recording_support::field("action", exact), recording_support::field("richardson", richardson),
                recording_support::field("direct_reference_species_term", direct),
                recording_support::field("criterion", tolerance)});
            recording_support::compare_custom(
                id, "d(nH)/dtheta", "J per native parameter unit", richardson, exact, inputs,
                "Relationship", "richardson-centered-resolve",
                "the extensive total enthalpy action equals rebuilt-model re-solves",
                "abs(Richardson-action)<=1e-6*abs(Richardson) && abs(fine-action)<=abs(coarse-action)/2",
                agrees, "s2-neutral-reference", id, "d(nH)/dtheta", "T,P,theta,h",
                "(4 fine-coarse)/3=d(nH)/dtheta", "coordinate");
            return Ladder{exact, coarse, fine, richardson, direct};
        };
        const auto discrimination = [&](const std::string& id, const Ladder& row) {
            const double tolerance = 1.0e-6 * std::abs(row.richardson);
            const double threshold = 10.0 * tolerance;
            const bool separated = threshold > 0.0 && std::abs(row.direct) >= threshold;
            const auto inputs = recording_support::fields({
                recording_support::field("action", row.action),
                recording_support::field("richardson", row.richardson),
                recording_support::field("criterion", tolerance),
                recording_support::field("omitted_term", row.direct),
                recording_support::field("omission_over_criterion",
                    tolerance > 0.0 ? std::abs(row.direct) / tolerance : std::numeric_limits<double>::quiet_NaN())});
            const std::string comparison = id + "/omission-discrimination";
            recording_support::compare_custom(
                comparison.c_str(), "omitted_direct_enthalpy_term", "J per native parameter unit",
                std::abs(row.direct), threshold,
                inputs, "Discrimination", "reference-species-parameter-chain",
                "omitting the direct referenced-species enthalpy slope exceeds ten action tolerances",
                "abs(omitted_term)>=10*1e-6*abs(Richardson)", separated,
                "s2-neutral-reference", id.c_str(), "n_sum_x_dh_dtheta", "n,x,theta",
                "n sum_i x_i dh_i/dtheta", "case");
        };

        // Case N reaction b at the system-pressure reference.
        epcsaft::ThermochemistryRecord neutral_record;
        neutral_record.standard_pressure = 1.0e5;
        neutral_record.component_intervals = {
            {{200.0, 400.0, true, true,
              {epcsaft::IdealPolynomial{{52.0, 0.12}, 300.0, 0.0, 0.0}, 1.0e5}}},
            {}};
        const epcsaft::HelmholtzModel neutral_model(identical_parameters, {}, neutral_record);
        const epcsaft::ReactionLogPolynomial neutral_law{
            std::log(0.3), 200.0, 3.0, -0.004, 300.0, 250.0, 350.0, "analytic-source"};
        const std::vector<epcsaft::ReactionReference> neutral_sources = {reference};
        for (std::size_t source = 0; source < neutral_sources.size(); ++source) {
            auto problem = equilibrium_support::flash(2);
            problem.temperature = epcsaft::equilibrium::FreeValue::fixed(260.0);
            problem.pressure = epcsaft::equilibrium::FreeValue::fixed(3.0e6);
            problem.feed.values = {2.0, 0.2};
            problem.source_reference = identical_declaration;
            problem.reactions = {{"2A=B", identical_nu, std::nullopt, std::nullopt,
                                  neutral_law, neutral_sources[source]}};
            ActionDirection b_direction;
            b_direction.feed_amounts_mol = {0.0, 0.0};
            b_direction.reaction_coefficients = {{{0.0, 1.0, 0.0, 0.0}}};
            const auto value_at_shift = [&](double shift) {
                auto shifted = problem;
                shifted.reactions[0].correlation->b += shift;
                return extensive(neutral_model, shifted);
            };
            const std::string id = std::string("reference-calorics/enthalpy-reaction/N/")
                + (source == 0 ? "system" : "fixed") + "/reaction-b/resolve";
            const auto row = action_ladder(id.c_str(), neutral_model, problem, b_direction, 1.0, value_at_shift);
            (void)row;
            // At Tref, phi_c=0 but phi_c'=1/T: a mixed chemical action cannot be an available zero.
            auto centered = problem;
            centered.reactions[0].correlation->reference_temperature = 260.0;
            const auto center = epcsaft::equilibrium::solve_equilibrium(neutral_model, centered);
            const auto compiled = epcsaft::equilibrium::compile_problem(neutral_model, centered);
            epcsaft::equilibrium::SolvedStateObservable fraction;
            fraction.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
            fraction.component = 1;
            ActionDirection c_direction;
            c_direction.feed_amounts_mol = {0.0, 0.0};
            c_direction.reaction_coefficients = {{{0.0, 0.0, 1.0, 0.0}}};
            const ActionDirection temperature{1.0, 0.0, {0.0, 0.0}, {}};
            const auto higher = epcsaft::equilibrium::solved_state_actions(
                compiled, center, {{fraction}, {temperature, c_direction}});
            CHECK(higher.results[0].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                  && higher.results[0].action.has_value());
            const auto x_at = [&](double t, double c) {
                auto changed = centered;
                changed.temperature.value += t;
                changed.reactions[0].correlation->c += c;
                const auto state = epcsaft::equilibrium::solve_equilibrium(neutral_model, changed, center);
                CHECK(state.has_state());
                return state.mole_fractions[1];
            };
            const auto difference = [&](double step) {
                return (x_at(step, step)-x_at(step,-step)-x_at(-step,step)+x_at(-step,-step))/(4*step*step);
            };
            const double exact=*higher.results[0].action;
            const reactive_caloric_checks::Quantity quantity{"d2x_dTdc", "1", 1e-11, fraction};
            CHECK(reactive_caloric_checks::scalar_ladder(
                "reference-calorics/parameter/reaction-c-at-Tref/mixed-action", quantity, exact, 0.2,
                recording_support::fields({recording_support::field("Tref_K",260.0),
                    recording_support::field("temperature_increment_K",1.0),recording_support::field("coefficient_c_increment",1.0)}),
                1e-6, difference));
        }

        // Case I water epsilon/k and reaction b, both system- and fixed-pressure references.
        constexpr double ionic_temperature = 298.15, ionic_pressure = 1.0e5;
        const auto ionic_parameters = properties_support::aqueous_electrolyte();
        const std::vector<epcsaft::ActiveParameter> ionic_active = {
            {epcsaft::ParameterFamily::DispersionEnergyOverK, {0}}};
        epcsaft::ThermochemistryRecord water_record;
        water_record.standard_pressure = 1.0e5;
        water_record.component_intervals = {
            {{200.0, 500.0, true, true,
              {epcsaft::IdealPolynomial{{33.59}, 298.15, -241826.0, 188.83}, 1.0e5}}},
            {}, {}};
        const epcsaft::HelmholtzModel ionic_model(ionic_parameters, ionic_active, water_record);
        const std::vector<epcsaft::ReactionReference> ionic_sources = {
            reference, {2.0e5, reference.source_basis}};
        for (std::size_t source = 0; source < ionic_sources.size(); ++source) {
            auto problem = equilibrium_support::flash(3);
            problem.temperature = epcsaft::equilibrium::FreeValue::fixed(ionic_temperature);
            problem.pressure = epcsaft::equilibrium::FreeValue::fixed(ionic_pressure);
            problem.feed.values = {1.0, 1.0e-6, 1.0e-6};
            problem.source_reference = ratios({0.0, 1.0, 1.0});
            problem.reactions = {{"H2O=Na++Cl-", {-1.0, 1.0, 1.0}, std::nullopt, std::nullopt,
                epcsaft::ReactionLogPolynomial{std::log(1.0e-8) + 2000.0 / ionic_temperature, -2000.0, 0.0, 0.0,
                 ionic_temperature, 250.0, 350.0, "analytic-source"}, ionic_sources[source]}};
            const ActionDirection epsilon_direction{
                0.0, 0.0, {0.0, 0.0, 0.0}, {1.0}, {}};
            const auto epsilon_value = [&](double shift) {
                auto shifted = ionic_parameters;
                shifted.dispersion_energy_over_k[0] += shift;
                const epcsaft::HelmholtzModel changed(shifted, ionic_active, water_record);
                return extensive(changed, problem);
            };
            const std::string epsilon_id = std::string("reference-calorics/enthalpy-parameter/I/")
                + (source == 0 ? "system" : "fixed") + "/water-epsilon/resolve";
            const auto epsilon_row = action_ladder(
                epsilon_id.c_str(), ionic_model, problem, epsilon_direction, 0.05, epsilon_value);
            discrimination(epsilon_id, epsilon_row);

            ActionDirection b_direction;
            b_direction.feed_amounts_mol = {0.0, 0.0, 0.0};
            b_direction.active_parameters = {0.0};
            b_direction.reaction_coefficients = {{{0.0, 1.0, 0.0, 0.0}}};
            const auto b_value = [&](double shift) {
                auto changed = problem;
                changed.reactions[0].correlation->b += shift;
                return extensive(ionic_model, changed);
            };
            const std::string b_id = std::string("reference-calorics/enthalpy-reaction/I/")
                + (source == 0 ? "system" : "fixed") + "/reaction-b/resolve";
            const auto b_row = action_ladder(b_id.c_str(), ionic_model, problem, b_direction, 1.0, b_value);
            discrimination(b_id, b_row);
        }
    }


    const auto charged = epcsaft::reaction_reference(
        epcsaft::HelmholtzModel(properties_support::aqueous_electrolyte()), 298.15, 1.0e5,
        ratios({0.0, 1.0, 1.0}), {-1.0, 1.0, 0.0}, reference);
    CHECK(charged.status == epcsaft::NeutralReferenceStatus::NonNeutralDeclaration);
    auto referenced = equilibrium_support::flash(3);
    referenced.feed.values = {2.0, 0.2, 0.2};
    referenced.source_reference = ratios({0.0, 1.0, 1.0});
    const auto law = epcsaft::ReactionLogPolynomial{
        std::log(0.3), 0.0, 0.0, 0.0, 300.0, 250.0, 350.0, "native-reference-source"};
    referenced.reactions = {{"A=2B", {-1.0, 2.0, 0.0}, std::nullopt, std::nullopt,
                             law, reference}};
    auto free_pressure_reference = referenced;
    free_pressure_reference.phases = {
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
            epcsaft::equilibrium::AmountSpec::pinned(0.0)),
    };
    free_pressure_reference.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto free_pressure_compiled = epcsaft::equilibrium::compile_problem(
        reference_model, free_pressure_reference);
    const auto free_pressure_point = free_pressure_compiled.generic_point();
    const auto free_pressure_prepared = free_pressure_compiled.prepare(free_pressure_point);
    CHECK(std::all_of(free_pressure_prepared.standard.pressure_derivative.begin(),
                      free_pressure_prepared.standard.pressure_derivative.end(),
                      [](double value) { return std::isfinite(value); }));
    const std::size_t control_width = 2 + reference_model.parameters().component_count();
    std::vector<double> pressure_control(control_width, 0.0);
    pressure_control[1] = 1.0;
    std::vector<std::vector<double>> pressure_state(3,
        std::vector<double>(free_pressure_prepared.point.size(), 0.0));
    pressure_state[0] = free_pressure_prepared.point;
    pressure_state[1][free_pressure_compiled.primal_size - 1] = 1.0;
    const auto pressure_series = free_pressure_compiled.row_series(
        free_pressure_prepared,
        pressure_state,
        pressure_control, 2);
    CHECK(pressure_series.energy == epcsaft::PropertyActionStatus::Available);
    const auto chemical = std::find_if(
        free_pressure_compiled.rows.begin(), free_pressure_compiled.rows.end(),
        [](const auto& row) { return row.kind == epcsaft::equilibrium::RowKind::Chemical; });
    CHECK(chemical != free_pressure_compiled.rows.end());
    if (chemical != free_pressure_compiled.rows.end()) {
        const std::size_t row = static_cast<std::size_t>(
            std::distance(free_pressure_compiled.rows.begin(), chemical));
        const double pressure_step = 100.0;
        auto plus = free_pressure_prepared.point;
        auto minus = free_pressure_prepared.point;
        plus[free_pressure_compiled.primal_size - 1] += pressure_step;
        minus[free_pressure_compiled.primal_size - 1] -= pressure_step;
        const auto plus_rows = free_pressure_compiled.evaluate_at(
            free_pressure_compiled.prepare(plus)).values;
        const auto minus_rows = free_pressure_compiled.evaluate_at(
            free_pressure_compiled.prepare(minus)).values;
        const auto central_rows = free_pressure_compiled.evaluate_at(free_pressure_prepared).values;
        const double centered_first = (plus_rows[row] - minus_rows[row])
            / (2.0 * pressure_step);
        const double centered_second_coefficient = (plus_rows[row]
            - 2.0 * central_rows[row] + minus_rows[row])
            / (2.0 * pressure_step * pressure_step);
        const auto close_pressure_jet = [](double actual, double expected) {
            return std::isfinite(actual) && std::isfinite(expected)
                && std::abs(actual - expected) <= 2.0e-3
                    * std::max({std::abs(actual), std::abs(expected), 1.0e-12});
        };
        CHECK(close_pressure_jet(pressure_series.rows[1][row], centered_first)
              && close_pressure_jet(pressure_series.rows[2][row], centered_second_coefficient));
    }
    bool pressure_third_rejected = false;
    try {
        (void)free_pressure_compiled.row_series(
            free_pressure_prepared,
            {free_pressure_prepared.point,
             std::vector<double>(free_pressure_prepared.point.size(), 0.0),
             std::vector<double>(free_pressure_prepared.point.size(), 0.0),
             std::vector<double>(free_pressure_prepared.point.size(), 0.0)},
            pressure_control, 3);
    } catch (const std::invalid_argument&) {
        pressure_third_rejected = true;
    }
    CHECK(pressure_third_rejected);
    const auto fixed_pressure_compiled = epcsaft::equilibrium::compile_problem(
        reference_model, referenced);
    const auto fixed_pressure_prepared = fixed_pressure_compiled.prepare(
        fixed_pressure_compiled.generic_point());
    bool fixed_pressure_third_rejected = false;
    try {
        (void)fixed_pressure_compiled.row_series(
            fixed_pressure_prepared,
            {fixed_pressure_prepared.point,
             std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
             std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
             std::vector<double>(fixed_pressure_prepared.point.size(), 0.0)},
            pressure_control, 3);
    } catch (const std::invalid_argument&) {
        fixed_pressure_third_rejected = true;
    }
    CHECK(fixed_pressure_third_rejected);
    std::fill(pressure_control.begin(), pressure_control.end(), 0.0);
    pressure_control[2] = 1.0;
    const auto fixed_feed_series = fixed_pressure_compiled.row_series(
        fixed_pressure_prepared,
        {fixed_pressure_prepared.point,
         std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
         std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
         std::vector<double>(fixed_pressure_prepared.point.size(), 0.0)},
        pressure_control, 3);
    CHECK(fixed_feed_series.energy == epcsaft::PropertyActionStatus::Available);
    epcsaft::ThermochemistryRecord thermal_record;
    thermal_record.standard_pressure = 1.0e5;
    for (std::size_t i = 0; i < 3; ++i)
        thermal_record.component_intervals.push_back(
            properties_support::p3_record(properties_support::p3_forms()[0].second)
                .component_intervals.front());
    epcsaft::HelmholtzModel thermal_reference_model(
        properties_support::neutral_ternary(), {}, thermal_record);
    const auto referenced_result = epcsaft::equilibrium::solve_equilibrium(
        thermal_reference_model, referenced);
    CHECK(referenced_result.success
          && referenced_result.reference_diagnostics.size() == 1
          && referenced_result.reference_diagnostics[0].terminal_contraction.has_value());
    const auto compiled_reference = epcsaft::equilibrium::compile_problem(
        thermal_reference_model, referenced);
    epcsaft::equilibrium::SolvedStateActionRequest action_request;
    epcsaft::equilibrium::SolvedStateObservable composition_observable;
    composition_observable.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
    composition_observable.component = 1;
    action_request.observables = {composition_observable};
    action_request.directions = {{1.0, 0.0, {0.0, 0.0, 0.0}, {}}};
    const auto temperature_action = epcsaft::equilibrium::solved_state_actions(
        compiled_reference, referenced_result, action_request);
    CHECK(temperature_action.results[0].status
              == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && temperature_action.results[0].action.has_value());
    composition_observable.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseProperty;
    composition_observable.property = epcsaft::PropertyObservable::IdealEnthalpy;
    const auto total_observable = [&] {
        auto result = composition_observable;
        result.property = epcsaft::PropertyObservable::TotalEnthalpy;
        return result;
    }();
    const auto residual_observable = [&] {
        auto result = composition_observable;
        result.property = epcsaft::PropertyObservable::ResidualEnthalpy;
        return result;
    }();
    action_request.observables = {composition_observable, total_observable, residual_observable};
    action_request.directions = {{0.0, 0.0, {0.0, 0.0, 0.0}, {}}};
    const auto caloric = epcsaft::equilibrium::solved_state_actions(
        compiled_reference, referenced_result, action_request);
    const std::vector<double> central_x(
        referenced_result.mole_fractions.begin(), referenced_result.mole_fractions.begin() + 3);
    const auto central_state = thermal_reference_model.state(
        referenced_result.temperature.front(), referenced_result.molar_densities[0], central_x);
    CHECK(central_state.ideal_enthalpy.has_value() && central_state.total_enthalpy.has_value()
          && std::isfinite(central_state.residual_enthalpy));
    CHECK(caloric.results.size() == 3
          && caloric.results[0].status
                 == epcsaft::equilibrium::SolvedStateActionStatus::ReferenceUnavailable
          && !caloric.results[0].value.has_value() && !caloric.results[0].action.has_value()
          && caloric.results[1].status
                 == epcsaft::equilibrium::SolvedStateActionStatus::ReferenceUnavailable
          && !caloric.results[1].value.has_value() && !caloric.results[1].action.has_value()
          && caloric.results[2].status
                 == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && caloric.results[2].value.has_value()
          && close_with_tolerance(*caloric.results[2].value, central_state.residual_enthalpy, 1.0e-9));
    const std::vector<epcsaft::equilibrium::SolvedStateObservable> source_observables = {
        temperature_action.results[0].observable, composition_observable};
    const std::vector<epcsaft::equilibrium::SolvedStateActionDirection> source_directions = {
        {1.0, 0.0, {0.0, 0.0, 0.0}, {}}, {0.0, 0.0, {0.01, 0.0, 0.0}, {}}};
    const auto source_batches = epcsaft::equilibrium::solved_state_first_actions(
        compiled_reference, referenced_result, source_observables, source_directions);
    bool source_match = source_batches.size() == source_directions.size();
    double source_worst = 0.0;
    for (std::size_t j = 0; j < source_batches.size() && source_match; ++j) {
        const auto separate = epcsaft::equilibrium::solved_state_actions(
            compiled_reference, referenced_result, {source_observables, {source_directions[j]}});
        source_match = source_batches[j].results.size() == separate.results.size();
        for (std::size_t i = 0; i < separate.results.size() && source_match; ++i) {
            const auto& actual = source_batches[j].results[i];
            const auto& expected = separate.results[i];
            source_match = actual.order == expected.order && actual.unit == expected.unit
                && actual.status == expected.status && actual.diagnostic_message == expected.diagnostic_message
                && actual.value.has_value() == expected.value.has_value()
                && actual.action.has_value() == expected.action.has_value();
            if (actual.value && expected.value)
                source_worst = std::max(source_worst, std::abs(*actual.value - *expected.value) / std::max(1.0, std::abs(*expected.value)));
            if (actual.action && expected.action)
                source_worst = std::max(source_worst, std::abs(*actual.action - *expected.action) / std::max(1.0, std::abs(*expected.action)));
        }
    }
    CHECK(source_match && source_worst <= 2.0e-12
          && source_batches[0].results[0].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && source_batches[0].results[0].value.has_value() && source_batches[0].results[0].action.has_value()
          && source_batches[1].results[0].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && source_batches[1].results[0].value.has_value() && source_batches[1].results[0].action.has_value()
          && source_batches[1].results[1].status == epcsaft::equilibrium::SolvedStateActionStatus::ReferenceUnavailable
          && source_batches[1].results[1].diagnostic_message.has_value());
}

}  // namespace equilibrium_checks
