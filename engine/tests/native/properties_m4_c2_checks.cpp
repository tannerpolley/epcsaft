// M4-M9 and C2: Hessians, pressure slopes, parameter derivatives, density, cache invariance and
// parameter coordinates.
#include "properties_checks_shared.hpp"

namespace properties_checks {

void hessian_reference() {
    const auto associating_parameters = properties_support::associating_binary();
    const std::vector<double> associating_amounts = {0.1385, 0.8615};
    epcsaft::HelmholtzModel associating(associating_parameters);
    const auto associating_result = associating.evaluate(
        376.88, 1.0e-4, associating_amounts, epcsaft::DerivativeOrder::Hessian
    );
    const auto associating_reference = properties_support::reference_tensor(
        associating_parameters, 376.88, 1.0e-4, associating_amounts,
        associating_result.site_fractions
    );
    const auto associating_inputs = recording_support::fields({
        model_inputs(associating_parameters, 376.88, 1.0e-4, associating_amounts),
        recording_support::field("site_fractions", associating_result.site_fractions),
        recording_support::field("reference_site_fractions", associating_result.site_fractions),
    });
    CHECK(associating_result.gradient.size() == associating_reference.gradient.size());
    for (std::size_t coordinate = 0;
         coordinate < associating_reference.gradient.size(); ++coordinate) {
        recording_support::compare(
            (std::string("associating-gradient/") + std::to_string(coordinate)).c_str(),
            "gradient", "1", associating_result.gradient[coordinate],
            associating_reference.gradient[coordinate],
            associating_inputs,
            "Relationship", "AD-reference", "implicit-site gradient matches reference tensor",
            "reference_or_one", 0.0, 1.0e-10,
            close(associating_result.gradient[coordinate], associating_reference.gradient[coordinate], 1e-10),
            "m4-associating-gradient", std::to_string(coordinate).c_str(), "\\nabla a", "T,V,\\mathbf{n},\\chi",
            "\\nabla a_{implicit}=\\nabla a_{reference}", "coordinate");
    }
    for (std::size_t entry = 0; entry < associating_reference.hessian.size(); ++entry) {
        const std::size_t width = associating_reference.gradient.size();
        const std::size_t row = entry / width;
        const std::size_t column = entry % width;
        const auto point = std::to_string(row) + "-" + std::to_string(column);
        recording_support::compare(
            (std::string("associating-hessian/") + point).c_str(), "hessian", "1",
            associating_result.hessian[entry], associating_reference.hessian[entry],
            associating_inputs,
            "Relationship", "AD-reference", "implicit-site Hessian matches reference tensor",
            "reference_or_one", 0.0, 1.0e-10,
            close(associating_result.hessian[entry], associating_reference.hessian[entry], 1e-10),
            "m4-associating-hessian", point.c_str(), "H^{res}", "T,V,\\mathbf{n},\\chi",
            "H_{ij}^{implicit}=H_{ij}^{reference}", "row,col");
    }

    const auto nine_parameters = properties_support::nine_species();
    const std::vector<double> nine_amounts = {
        0.795, 0.10, 0.05, 0.02, 0.02, 0.005, 0.001, 0.008, 0.001,
    };
    epcsaft::HelmholtzModel nine(
        nine_parameters,
        {{epcsaft::ParameterFamily::DispersionEnergyOverK, {0}}}
    );
    const auto nine_result = nine.evaluate(
        300.0, 1.0 / 1.2e4, nine_amounts, epcsaft::DerivativeOrder::Hessian
    );
    const auto nine_reference = properties_support::reference_tensor(
        nine_parameters, 300.0, 1.0 / 1.2e4, nine_amounts,
        nine_result.site_fractions,
        std::vector<epcsaft::ActiveParameter>{epcsaft::ActiveParameter{
            epcsaft::ParameterFamily::DispersionEnergyOverK, {0}
        }}
    );
    const auto nine_inputs = recording_support::fields({
        model_inputs(nine_parameters, 300.0, 1.0 / 1.2e4, nine_amounts),
        recording_support::field("site_fractions", nine_result.site_fractions),
        recording_support::field("reference_site_fractions", nine_result.site_fractions),
        recording_support::field("active_parameter", "dispersion_energy_over_k[0]"),
    });
    for (std::size_t coordinate = 0;
         coordinate < nine_reference.gradient.size(); ++coordinate) {
        recording_support::compare(
            (std::string("nine-gradient/") + std::to_string(coordinate)).c_str(),
            "gradient", "1", nine_result.gradient[coordinate], nine_reference.gradient[coordinate],
            nine_inputs,
            "Relationship", "AD-reference", "active-parameter gradient matches reference tensor",
            "reference_or_one", 0.0, 1.0e-10,
            close(nine_result.gradient[coordinate], nine_reference.gradient[coordinate], 1e-10),
            "m4-nine-gradient", std::to_string(coordinate).c_str(), "\\nabla a", "T,V,\\mathbf{n},k_{ij}",
            "\\nabla a_{active}=\\nabla a_{reference}", "coordinate");
    }
    for (std::size_t entry = 0; entry < nine_reference.hessian.size(); ++entry) {
        const std::size_t width = nine_reference.gradient.size();
        const std::size_t row = entry / width;
        const std::size_t column = entry % width;
        const auto point = std::to_string(row) + "-" + std::to_string(column);
        recording_support::compare(
            (std::string("nine-hessian/") + point).c_str(), "hessian", "1",
            nine_result.hessian[entry], nine_reference.hessian[entry],
            nine_inputs,
            "Relationship", "AD-reference", "active-parameter Hessian matches reference tensor",
            "reference_or_one", 0.0, 1.0e-10,
            close(nine_result.hessian[entry], nine_reference.hessian[entry], 1e-10),
            "m4-nine-hessian", point.c_str(), "H^{res}", "T,V,\\mathbf{n},k_{ij}",
            "H_{ij}^{active}=H_{ij}^{reference}", "row,col");
    }
}

void hessian_fd_and_pressure_slope() {
    // Compare the projected chart with the independent raw mass-action elimination,
    // including absent associating species, absent ions and an active parameter.
    auto absent_association = properties_support::associating_binary();
    // Manufactured topology: all four sites belong to the absent first species.
    absent_association.site_component = {0, 0, 0, 0};
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> projected_cases = {
        {properties_support::neutral_binary(), {0.4, 0.6}},
        {properties_support::associating_binary(), {0.0, 1.0}},
        {properties_support::aqueous_electrolyte(), {1.0, 0.0, 0.0}},
        {absent_association, {0.0, 1.0}},
        {properties_support::nine_species(), {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0}},
    };
    for (std::size_t index = 0; index < projected_cases.size(); ++index) {
        const auto& [parameters, amounts] = projected_cases[index];
        bool charged_support = false;
        for (std::size_t i = 0; i < amounts.size(); ++i)
            charged_support = charged_support || (amounts[i] > 0.0 && parameters.charge_number[i] != 0);
        const std::vector<epcsaft::ActiveParameter> active = {
            {epcsaft::ParameterFamily::DispersionEnergyOverK, {0}}};
        epcsaft::HelmholtzModel model(parameters, active);
        (void)model.evaluate(300.0, 1.0 / 4000.0, amounts, epcsaft::DerivativeOrder::Hessian);
        constexpr double temperature = 305.0, volume = 1.0 / 5000.0;
        const auto sites = model.evaluate(temperature, volume, amounts, epcsaft::DerivativeOrder::Value).site_fractions;
        const auto reference = properties_support::reference_tensor(parameters, temperature, volume, amounts, sites, active);
        const auto inputs = model_inputs(parameters, temperature, volume, amounts);
        const std::array orders = {epcsaft::DerivativeOrder::IsothermalPhaseHessian,
            epcsaft::DerivativeOrder::PhaseHessian, epcsaft::DerivativeOrder::Hessian};
        for (std::size_t mode = 0; mode < orders.size(); ++mode) {
            const auto actual = model.evaluate(temperature, volume, amounts, orders[mode]);
            const std::size_t width = reference.gradient.size();
            std::vector<std::size_t> selected;
            for (std::size_t i = 0; i < width; ++i) {
                // At no-ion support the charged insertion Hessian is outside the admitted restriction.
                if (!charged_support && i >= 2 && i < 2 + amounts.size()
                    && amounts[i - 2] == 0.0 && parameters.charge_number[i - 2] != 0) continue;
                if (mode == 2 || i == 1 || (i == 0 && mode == 1)
                    || (i >= 2 && i < 2 + amounts.size() && amounts[i - 2] > 0.0)) selected.push_back(i);
            }
            const auto compare = [&](double value, double expected, const std::string& entry) {
                const auto point = std::to_string(index) + "/" + std::to_string(mode) + "/" + entry;
                recording_support::compare(("projected/" + point).c_str(), "Helmholtz_derivative", "1",
                    value, expected, inputs, "Relationship", "raw-implicit-reference",
                    "selected derivatives and full controls preserve implicit mass-action chains",
                    "reference_or_one", 0.0, 1e-10, close(value, expected, 1e-10),
                    "m5-projected", point.c_str());
            };
            for (const auto row : selected) {
                compare(actual.gradient[row], reference.gradient[row], "g" + std::to_string(row));
                for (const auto column : selected)
                    compare(actual.hessian[row * width + column], reference.hessian[row * width + column],
                        std::to_string(row) + "-" + std::to_string(column));
            }
        }
    }
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> cases = {
        {properties_support::neutral_binary(), {0.4, 0.6}},
        {properties_support::associating_binary(), {0.1385, 0.8615}},
        {properties_support::aqueous_electrolyte(), {0.98, 0.01, 0.01}},
         {properties_support::nine_species(),
         {0.795, 0.10, 0.05, 0.02, 0.02, 0.005, 0.001, 0.008, 0.001}},
    };
    for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
        const auto& item = cases[case_index];
        epcsaft::HelmholtzModel model(item.first);
        const double temperature = 300.0;
        const double volume = 1.0 / 1.2e4;
        const auto actual = model.evaluate(
            temperature, volume, item.second, epcsaft::DerivativeOrder::Hessian
        );
        const std::size_t width = actual.gradient.size();
        const auto inputs = model_inputs(item.first, temperature, volume, item.second);
        for (std::size_t coordinate = 0; coordinate < item.second.size() + 2;
             ++coordinate) {
            const double step = coordinate == 0 ? temperature * 1e-5
                : coordinate == 1 ? volume * 1e-5 : 1e-5;
            auto plus = item.second;
            auto minus = item.second;
            double plus_temperature = temperature;
            double minus_temperature = temperature;
            double plus_volume = volume;
            double minus_volume = volume;
            if (coordinate == 0) {
                plus_temperature += step;
                minus_temperature -= step;
            } else if (coordinate == 1) {
                plus_volume += step;
                minus_volume -= step;
            } else {
                plus[coordinate - 2] += step;
                minus[coordinate - 2] -= step;
            }
            const auto plus_gradient = model.evaluate(
                plus_temperature, plus_volume, plus,
                epcsaft::DerivativeOrder::Gradient
            ).gradient;
            const auto minus_gradient = model.evaluate(
                minus_temperature, minus_volume, minus,
                epcsaft::DerivativeOrder::Gradient
            ).gradient;
            for (std::size_t row = 0; row < width; ++row) {
                const double finite_difference =
                    (plus_gradient[row] - minus_gradient[row]) / (2.0 * step);
                const auto point = std::to_string(case_index) + "/"
                    + std::to_string(row) + "-" + std::to_string(coordinate);
                recording_support::compare(
                    (std::string("hessian/") + point).c_str(), "hessian", "1",
                    actual.hessian[row * width + coordinate], finite_difference,
                    recording_support::fields({inputs, recording_support::field("row", row),
                                               recording_support::field("column", coordinate),
                                               recording_support::field("step", step),
                                               recording_support::field("T_plus", plus_temperature),
                                               recording_support::field("T_minus", minus_temperature),
                                               recording_support::field("V_plus", plus_volume),
                                               recording_support::field("V_minus", minus_volume),
                                               recording_support::field("n_plus", plus),
                                               recording_support::field("n_minus", minus)}),
                    "Relationship", "finite-difference", "Hessian=derivative of gradient",
                    "reference_or_one", 0.0, 1.0e-6,
                    close(actual.hessian[row * width + coordinate], finite_difference, 1e-6),
                    "m5-hessian", point.c_str(), "H^{res}", "T,V,\\mathbf{n},\\Delta z",
                    "H_{ij}=\\Delta a_i/\\Delta z_j", "row,col");
            }
        }
        const auto slope_state = model.state(temperature, 1.2e4, item.second);
        const auto plus_state = model.state(temperature, 1.2e4 * (1.0 + 1e-5), item.second);
        const auto minus_state = model.state(temperature, 1.2e4 * (1.0 - 1e-5), item.second);
        const double pressure_density_fd =
            (plus_state.pressure - minus_state.pressure) / (2.0 * 1.2e4 * 1e-5);
        recording_support::compare(
            (std::string("pressure-density-slope/") + std::to_string(case_index)).c_str(),
            "dP_drho", "Pa m^3/mol", slope_state.pressure_density_derivative,
            pressure_density_fd,
            recording_support::fields({inputs,
                                       recording_support::field("rho_plus", 1.2e4 * (1.0 + 1e-5)),
                                       recording_support::field("rho_minus", 1.2e4 * (1.0 - 1e-5))}),
            "Relationship", "finite-difference",
            "pressure density slope matches finite difference", "reference_or_one", 0.0, 1.0e-6,
            close(slope_state.pressure_density_derivative, pressure_density_fd, 1e-6),
            "m5-pressure-slope", std::to_string(case_index).c_str(), "P_\\rho", "T,V,\\mathbf{n}",
            "P_\\rho=\\Delta P/\\Delta\\rho", "rho");
        const auto pressure_slope = model.evaluate(
            temperature, volume, item.second, epcsaft::DerivativeOrder::PressureSlope
        );
        const double volume_step = volume * 1e-5;
        const auto plus_volume_gradient = model.evaluate(
            temperature, volume + volume_step, item.second,
            epcsaft::DerivativeOrder::Gradient
        );
        const auto minus_volume_gradient = model.evaluate(
            temperature, volume - volume_step, item.second,
            epcsaft::DerivativeOrder::Gradient
        );
        const double volume_second_fd =
            (plus_volume_gradient.gradient[1] - minus_volume_gradient.gradient[1])
                / (2.0 * volume_step);
        recording_support::compare(
            (std::string("volume-second-derivative/") + std::to_string(case_index)).c_str(),
            "a_VV", "1", pressure_slope.volume_second_derivative,
            volume_second_fd,
            recording_support::fields({inputs, recording_support::field("step_V", volume_step),
                                       recording_support::field("V_plus", volume + volume_step),
                                       recording_support::field("V_minus", volume - volume_step)}),
            "Relationship", "finite-difference", "a_VV=derivative of a_V",
            "reference_or_one", 0.0, 1.0e-6,
            close(pressure_slope.volume_second_derivative, volume_second_fd, 1e-6),
            "m5-volume-curvature", std::to_string(case_index).c_str(), "a_{VV}", "T,V,\\mathbf{n},\\Delta V",
            "a_{VV}=\\Delta a_V/\\Delta V", "volume");
    }
}

void same_sign_pair_is_inert() {
    const std::vector<double> amounts = {
        0.795, 0.10, 0.05, 0.02, 0.02, 0.005, 0.001, 0.008, 0.001,
    };
    auto excluded = properties_support::nine_species();
    auto changed = excluded;
    changed.k_ij[3 * 9 + 7] = changed.k_ij[7 * 9 + 3] = 0.5;
    const auto base_state = epcsaft::HelmholtzModel(excluded).state(
        300.0, 1.2e4, amounts
    );
    const auto changed_state = epcsaft::HelmholtzModel(changed).state(
        300.0, 1.2e4, amounts
    );
    const auto inputs = recording_support::fields({
        model_inputs(excluded, 300.0, 1.0 / 1.2e4, amounts),
        recording_support::field("pair", "3,7"),
        recording_support::field("kij_changed", 0.5),
    });
    recording_support::compare(
        "excluded-inert", "dispersion", "1", base_state.dispersion, changed_state.dispersion,
        inputs, "Relationship", "ionic-pair-selection", "same-sign excluded pair is inert",
        "reference_or_one", 0.0, 1.0e-14,
        close(base_state.dispersion, changed_state.dispersion, 1e-14),
        "m6-excluded", "same-sign", "a^{disp}", "T,\\rho,\\mathbf{n},k_{ij}",
        "\\Delta a^{disp}_{same-sign}=0", "pair");

    auto all_pairs = excluded;
    all_pairs.ionic_dispersion = epcsaft::IonicDispersion::AllIonPairs;
    auto all_changed = all_pairs;
    all_changed.k_ij[3 * 9 + 7] = all_changed.k_ij[7 * 9 + 3] = 0.5;
    const auto all_base_state = epcsaft::HelmholtzModel(all_pairs).state(
        300.0, 1.2e4, amounts
    );
    const auto all_changed_state = epcsaft::HelmholtzModel(all_changed).state(
        300.0, 1.2e4, amounts
    );
    const bool all_pairs_active = all_base_state.dispersion != all_changed_state.dispersion;
    recording_support::compare_custom(
        "all-pairs-active", "dispersion", "1", all_changed_state.dispersion,
        all_base_state.dispersion, inputs, "Relationship", "ionic-pair-selection",
        "all ion pairs respond to the changed interaction", "model != reference",
        all_pairs_active, "m6-all-pairs", "all-signs", "a^{disp}",
        "T,\\rho,\\mathbf{n},k_{ij}", "\\Delta a^{disp}_{all-pairs}\\ne0", "pair");
}

void parameter_derivative_and_helpers() {
    epcsaft::Parameters parameters = properties_support::neutral_binary();
    epcsaft::HelmholtzModel model(
        parameters, {{epcsaft::ParameterFamily::Kij, {1, 2}}}
    );
    const auto derivative = model.evaluate(300.0, 1.0 / 1.2e4, {0.4, 0.6},
                                           epcsaft::DerivativeOrder::Gradient);
    const double step = 1e-5;
    auto plus = parameters;
    auto minus = parameters;
    plus.k_ij[1] += step; plus.k_ij[2] += step;
    minus.k_ij[1] -= step; minus.k_ij[2] -= step;
    const double up = epcsaft::HelmholtzModel(plus).evaluate(
        300.0, 1.0 / 1.2e4, {0.4, 0.6}, epcsaft::DerivativeOrder::Value).value;
    const double down = epcsaft::HelmholtzModel(minus).evaluate(
        300.0, 1.0 / 1.2e4, {0.4, 0.6}, epcsaft::DerivativeOrder::Value).value;
    const auto inputs = model_inputs(parameters, 300.0, 1.0 / 1.2e4, {0.4, 0.6});
    const auto parameter_inputs = recording_support::fields({
        inputs,
        recording_support::field("parameter", "kij[1,2]"), recording_support::field("step", step),
        recording_support::field("kij", parameters.k_ij),
        recording_support::field("kij_plus", plus.k_ij),
        recording_support::field("kij_minus", minus.k_ij),
    });
    const double parameter_fd = (up - down) / (2.0 * step);
    recording_support::compare(
        "parameter-gradient", "gradient_kij", "1", derivative.gradient.back(), parameter_fd,
        parameter_inputs, "Relationship", "finite-difference", "parameter derivative=central difference",
        "reference_or_one", 0.0, 1.0e-6,
        close(derivative.gradient.back(), parameter_fd, 1e-6), "m7-parameter-gradient",
        "kij[1,2]", "\\partial a/\\partial k_{ij}", "T,\\rho,\\mathbf{x},\\Delta k",
        "a_k=\\Delta a/\\Delta k", "parameter");
    bool threw = false;
    try { auto bad = parameters; bad.k_ij.pop_back(); (void)epcsaft::HelmholtzModel(std::move(bad)); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    threw = false;
    try {
        auto bad = properties_support::aqueous_electrolyte();
        bad.ionic_dispersion = {};
        (void)epcsaft::HelmholtzModel(std::move(bad));
    } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    threw = false;
    try {
        auto bad = properties_support::neutral_binary();
        bad.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
        (void)epcsaft::HelmholtzModel(std::move(bad));
    } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    threw = false;
    try { auto bad = parameters; bad.k_ij[1] = 1.0; (void)epcsaft::HelmholtzModel(std::move(bad)); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    threw = false;
    try {
        auto bad = properties_support::associating_binary();
        bad.site_component[0] = 2;
        (void)epcsaft::HelmholtzModel(std::move(bad));
    } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

void solvent_activity_properties() {
    using namespace epcsaft;
    using F = ParameterFamily;
    for (const auto& p : {properties_support::associating_binary(), properties_support::nine_species()}) {
        const bool electrolyte = p.component_count() == 9;
        const std::size_t s = electrolyte ? 0 : 1;
        std::vector<double> solvent(p.component_count(), 0.0); solvent[s] = 1.0;
        const std::vector<double> x = electrolyte
            ? std::vector<double>{0.8,0.1,0.05,0.02,0.02,0.005,0,0.005,0}
            : std::vector<double>{0.05,0.95};
        const std::vector<ActiveParameter> active{{F::DispersionEnergyOverK,{s}}};
        HelmholtzModel model(p, active);
        const auto ref = solvent_reference(model, 298.15, 1e5, solvent);
        const auto result = solvent_activity(model, 298.15, 1e5, x, ref);
        CHECK(result.solvent_activity.has_value());
        const double identity = x[s] * std::exp(require_component_property(result.state.log_fugacity_coefficient->at(s))
            - require_component_property(ref.pure_solvent_state->log_fugacity_coefficient->at(s)));
        CHECK(close(*result.solvent_activity, identity, 1e-10));
        CHECK(close(std::get<double>(*result.osmotic_coefficient), -std::log(identity) / ((1.0-x[s])/x[s]), 1e-10));
        CHECK(result.state.pressure_density_derivative > 0.0);
        const double scale = p.dispersion_energy_over_k[s];
        for (double relative : {1e-5,0.5e-5}) {
            const double h = scale * relative;
            auto plus = p, minus = p;
            plus.dispersion_energy_over_k[s] += h; minus.dispersion_energy_over_k[s] -= h;
            HelmholtzModel positive(plus), negative(minus);
            const auto a = solvent_activity(positive, 298.15, 1e5, x, solvent_reference(positive, 298.15, 1e5, solvent));
            const auto b = solvent_activity(negative, 298.15, 1e5, x, solvent_reference(negative, 298.15, 1e5, solvent));
            CHECK(close(scale * *result.solvent_activity_parameter_derivative[0],
                scale * (*a.solvent_activity-*b.solvent_activity)/(2*h), 2e-6));
            CHECK(close(scale * *result.osmotic_parameter_derivative[0],
                scale * (std::get<double>(*a.osmotic_coefficient)-std::get<double>(*b.osmotic_coefficient))/(2*h), 2e-6));
        }
        const auto pure = solvent_activity(model, 298.15, 1e5, solvent, ref);
        CHECK(pure.solvent_activity == 1.0 && pure.solvent_activity_parameter_derivative[0] == 0.0);
        CHECK(pure.osmotic_cause == "OsmoticUndefinedPureSolvent" && !pure.osmotic_parameter_derivative[0]);
        CHECK(solvent_activity(model, 299.0, 1e5, x, ref).cause == "SolventReferenceMismatch");
        auto absent = x; absent[s] = 0.0; absent[1-s] += x[s];
        CHECK(solvent_activity(model, 298.15, 1e5, absent, ref).cause.find("PureSolventUnavailable") == 0);
        auto failed = ref; failed.pure_solvent_state.reset();
        const auto missing_reference = solvent_activity(model, 298.15, 1e5, x, failed);
        CHECK(missing_reference.cause.find("PureSolventUnavailable") == 0
            && missing_reference.status == PropertyActionStatus::CentralPropertyUnavailable
            && !missing_reference.solvent_activity && !missing_reference.osmotic_coefficient);
        CHECK(solvent_reference(model, 298.15, 1e5, x).solvent_cause.find("PureSolventUnavailable") == 0);
        auto missing = ref;
        missing.pure_solvent_state->log_fugacity_parameter_derivative[s].reset();
        missing.pure_solvent_state->component_property_status[s] = PropertyActionStatus::AssociationIllConditioned;
        const auto unavailable = solvent_activity(model, 298.15, 1e5, x, missing);
        CHECK(unavailable.solvent_activity.has_value() && !unavailable.solvent_activity_parameter_derivative[0]
            && !unavailable.osmotic_parameter_derivative[0]);
        CHECK(unavailable.parameter_status[0] == PropertyActionStatus::AssociationIllConditioned);
    }
    HelmholtzModel model(properties_support::aqueous_electrolyte(), {{F::DispersionEnergyOverK,{0}}});
    const std::vector<double> solvent{1,0,0};
    const auto ref = solvent_reference(model, 298.15, 1e5, solvent);
    for (double q : {0.5e-8,1e-8,2e-8}) {
        // Construct the declared particle ratio directly: formula_composition rounds the exact-floor case below it.
        const double D = q * model.parameters().molar_mass_kg_per_mol[0], water = 1/(1+D);
        const std::vector<double> x{water,D*water/2,D*water/2};
        const auto result = solvent_activity(model, 298.15, 1e5, x, ref);
        CHECK(result.solvent_activity.has_value());
        CHECK(std::holds_alternative<double>(*result.osmotic_coefficient) == (q >= 1e-8));
        CHECK(result.solvent_activity_parameter_derivative[0].has_value());
        CHECK(result.osmotic_parameter_derivative[0].has_value() == (q >= 1e-8));
    }
    CHECK(solvent_reference(model, 650.0, 100.0, solvent).solvent_cause.find("VaporSideState") == 0);
    auto bad = *ref.pure_solvent_state; bad.pressure_density_derivative = 0.0;
    const auto x = formula_composition(solvent, model.parameters(), {1,2}, {1,1}, 1.0);
    const auto density_refused = solvent_activity(model, 298.15, 1e5, x, ref, bad);
    CHECK(density_refused.cause.find("DensityUnavailable") == 0
        && !density_refused.osmotic_coefficient && density_refused.osmotic_cause.empty());
    auto vapor = bad; vapor.pressure_density_derivative = 1.0;
    vapor.density_branch = DensityBranch::Unspecified; vapor.negative_slope_density = 2*vapor.molar_density;
    const auto vapor_refused = solvent_activity(model, 298.15, 1e5, x, ref, vapor);
    CHECK(vapor_refused.cause.find("VaporSideState") == 0
        && !vapor_refused.osmotic_coefficient && vapor_refused.osmotic_cause.empty());
}

void density_and_reference() {
    solvent_activity_properties();
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    const auto reference_rows = std::array<std::tuple<double, double, double, double>, 2>{
        std::tuple{200.0, 217800.69494692178, 0.007241100995,
                   0.03007 / 523.97698},
        std::tuple{240.0, 969152.4807171321, 0.001745731329,
                   0.03007 / 465.30887}};
    for (std::size_t case_index = 0; case_index < reference_rows.size(); ++case_index) {
        const auto& row = reference_rows[case_index];
        const double temperature = std::get<0>(row);
        const double pressure = std::get<1>(row);
        const auto vapor = ethane.state_at_pressure(
            temperature, pressure, {1.0}, epcsaft::DensityBranch::Vapor
        );
        const auto liquid = ethane.state_at_pressure(
            temperature, pressure, {1.0}, epcsaft::DensityBranch::Liquid
        );
        const auto inputs = recording_support::fields({
            model_inputs(ethane.parameters(), temperature, 1.0 / vapor.molar_density, {1.0}),
            recording_support::field("P", pressure), recording_support::field("branch", "vapor/liquid"),
            recording_support::field("vapor_reference_m3_per_mol", std::get<2>(row)),
            recording_support::field("liquid_reference_m3_per_mol", std::get<3>(row)),
        });
        recording_support::compare(
            (std::string("vapor-pressure/") + std::to_string(case_index)).c_str(),
            "pressure", "Pa", vapor.pressure, pressure, inputs,
            "Relationship", "density-root", "vapor state reproduces requested pressure",
            "reference_or_one", 0.0, 1e-9, close(vapor.pressure, pressure, 1e-9),
            "m8-vapor-pressure", std::to_string(case_index).c_str(), "P", "T,P,\\mathbf{x}", "P_{vapor}=P", "branch");
        recording_support::compare(
            (std::string("liquid-pressure/") + std::to_string(case_index)).c_str(),
            "pressure", "Pa", liquid.pressure, pressure, inputs,
            "Relationship", "density-root", "liquid state reproduces requested pressure",
            "reference_or_one", 0.0, 1e-9, close(liquid.pressure, pressure, 1e-9),
            "m8-liquid-pressure", std::to_string(case_index).c_str(), "P", "T,P,\\mathbf{x}", "P_{liquid}=P", "branch");
        recording_support::compare(
            (std::string("density-order/") + std::to_string(case_index)).c_str(),
            "vapor_density_below_liquid", "bool",
            vapor.molar_density < liquid.molar_density ? 1.0 : 0.0, 1.0, inputs,
            "Relationship", "density-branches", "vapor density is less than liquid density",
            "absolute", 0.0, 0.0, vapor.molar_density < liquid.molar_density,
            "m8-vapor-order", std::to_string(case_index).c_str(), "\\rho_v<\\rho_l",
            "T,P,\\mathbf{x},\\rho_v,\\rho_l", "\\rho_v<\\rho_l", "case");
        CHECK(close(vapor.pressure, pressure, 1e-9));
        CHECK(close(liquid.pressure, pressure, 1e-9));
        recording_support::compare(
            (std::string("vapor-volume/") + std::to_string(case_index)).c_str(),
            "specific_volume", "m3/mol", 1.0 / vapor.molar_density, std::get<2>(row),
            inputs, "Relationship", "historical-reference", "historical vapor volume agrees within 5%",
            "reference_or_one", 0.0, 0.05,
            close(1.0 / vapor.molar_density, std::get<2>(row), 0.05),
            "m8-vapor-volume", std::to_string(case_index).c_str(), "v_v", "T,P,\\mathbf{x}",
            "v_v", "case");
        recording_support::compare(
            (std::string("liquid-volume/") + std::to_string(case_index)).c_str(),
            "specific_volume", "m3/mol", 1.0 / liquid.molar_density, std::get<3>(row),
            inputs, "Relationship", "historical-reference", "historical liquid volume agrees within 5%",
            "reference_or_one", 0.0, 0.05,
            close(1.0 / liquid.molar_density, std::get<3>(row), 0.05),
            "m8-liquid-volume", std::to_string(case_index).c_str(), "v_l", "T,P,\\mathbf{x}",
            "v_l", "case");
        bool ambiguous = false;
        try {
            ethane.state_at_pressure(
                temperature, pressure, {1.0}, epcsaft::DensityBranch::Unspecified
            );
        } catch (const std::runtime_error&) {
            ambiguous = true;
        }
        CHECK(ambiguous);
    }
    const auto vapor = ethane.state_at_pressure(320.0, 1e6, {1.0},
                                                epcsaft::DensityBranch::Vapor);
    const auto liquid = ethane.state_at_pressure(320.0, 1e6, {1.0},
                                                 epcsaft::DensityBranch::Liquid);
    const auto supercritical_inputs = recording_support::fields({
        model_inputs(ethane.parameters(), 320.0, 1.0 / vapor.molar_density, {1.0}),
        recording_support::field("P", 1e6), recording_support::field("branch", "vapor/liquid"),
        recording_support::field("vapor_density", vapor.molar_density),
        recording_support::field("liquid_density", liquid.molar_density),
    });
    recording_support::compare(
        "supercritical-pressure", "pressure", "Pa", vapor.pressure, 1e6,
        supercritical_inputs, "Relationship", "density-root",
        "supercritical vapor state reproduces requested pressure", "reference_or_one",
        0.0, 1e-9, close(vapor.pressure, 1e6, 1e-9), "m8-supercritical-pressure", "state",
        "P", "T,P,\\mathbf{x}", "P=P_{requested}", "state");
    recording_support::compare(
        "supercritical-density", "molar_density", "mol/m3", liquid.molar_density,
        vapor.molar_density, supercritical_inputs, "Relationship", "density-branches",
        "supercritical branches agree", "reference_or_one", 0.0, 1e-10,
        close(liquid.molar_density, vapor.molar_density, 1e-10), "m8-supercritical-density", "state",
        "\\rho", "T,P,\\mathbf{x}", "\\rho_l=\\rho_v", "state");
    const auto associating = epcsaft::HelmholtzModel(properties_support::associating_binary());
    const auto associating_liquid = associating.state_at_pressure(
        313.15, 1e5, {0.1385, 0.8615}, epcsaft::DensityBranch::Liquid
    );
    recording_support::compare(
        "associating-pressure", "pressure", "Pa", associating_liquid.pressure, 1e5,
        recording_support::fields({
            model_inputs(associating.parameters(), 313.15,
                         1.0 / associating_liquid.molar_density, {0.1385, 0.8615}),
            recording_support::field("P", 1e5), recording_support::field("branch", "liquid"),
        }), "Relationship", "density-root", "associating liquid reproduces requested pressure",
        "reference_or_one", 0.0, 1e-9, close(associating_liquid.pressure, 1e5, 1e-9),
        "m8-associating", "state", "P", "T,P,\\mathbf{x}", "P_{liquid}=P", "state");
    const auto nine = epcsaft::HelmholtzModel(properties_support::nine_species());
    const auto nine_state = nine.state_at_pressure(
        313.15, 1e5,
        {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
        epcsaft::DensityBranch::Liquid
    );
    const std::vector<double> nine_composition = {
        0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0
    };
    const double nine_packing = epcsaft::packing_fraction_from_geometry(
        nine_state.molar_density, nine_composition, nine.parameters(), epcsaft::segment_diameter(
            313.15, nine.parameters(), epcsaft::segment_sigma(313.15, nine.parameters())));
    const bool nine_packing_ok =
        nine_packing < epcsaft::policy::CLOSE_PACKING_FRACTION.value;
    recording_support::compare(
        "nine-packing", "packing_below_limit", "bool", nine_packing_ok ? 1.0 : 0.0,
        1.0, recording_support::fields({
            model_inputs(nine.parameters(), 313.15, 1.0 / nine_state.molar_density,
                         nine_composition),
            recording_support::field("packing_fraction", nine_packing),
            recording_support::field("packing_limit", epcsaft::policy::CLOSE_PACKING_FRACTION.value),
        }), "Relationship", "packing-bound", "phase packing fraction is below close packing",
        "absolute", 0.0, 0.0, nine_packing_ok, "m8-nine", "state", "\\eta<\\eta_{max}",
        "T,P,\\mathbf{x},\\rho", "\\eta<\\eta_{max}", "state");
    const auto aqueous = properties_support::aqueous_electrolyte();
    epcsaft::HelmholtzModel model(aqueous);
    const auto reference = epcsaft::infinite_dilution_reference(
        model, 298.15, 1e5, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}
    );
    const auto activity = epcsaft::electrolyte_activity(
        model, 298.15, 1e5, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}, 1e-8, reference
    );
    const auto pure = model.state_at_pressure(
        298.15, 1e5, {1.0, 0.0, 0.0}, epcsaft::DensityBranch::Liquid
    );
    const auto reference_inputs = recording_support::fields({
        model_inputs(model.parameters(), 298.15, 1.0 / pure.molar_density,
                     {1.0, 0.0, 0.0}),
        recording_support::field("P", 1e5),
        recording_support::field("ions", std::vector<double>{1.0, 2.0}),
        recording_support::field("stoichiometry", std::vector<double>{1.0, 1.0}),
        recording_support::field("activity_tolerance", 1e-8),
        recording_support::field("reference_molality", reference.molality),
    });
    recording_support::observe(
        "activity-finite", "mean_ionic_log_molality_scale", "1",
        activity.mean_ionic_log_molality_scale, std::nullopt, reference_inputs, "Relationship", "finite-value",
        "electrolyte activity scale is finite", "m8-activity-finite", "state", "\\log\\gamma",
        "T,P,\\mathbf{x}", "\\log\\gamma", "state");
    CHECK(!reference.convergence_error);
    CHECK(std::isfinite(activity.mean_ionic_log_molality_scale));
    for (double m : {0.001, 0.2, 1.0, 3.0, 7.0}) {
        const auto scalar = epcsaft::electrolyte_activity(model, 298.15, 1e5,
            {1.0, 0.0, 0.0}, {1, 2}, {1, 1}, m, reference);
        const auto full = epcsaft::electrolyte_activity(model, 298.15, 1e5,
            {1.0, 0.0, 0.0}, {1, 2}, {1, 1}, m, reference, std::nullopt, std::vector<double>{0.0, m, m});
        CHECK(std::abs(scalar.mean_ionic_log_molality_scale - full.mean_ionic_log_molality_scale) <= 1e-10);
    }
    std::vector<double> water(9, 0.0);
    water[0] = 1.0;
    const auto mixed_reference = epcsaft::infinite_dilution_reference(nine, 298.15, 1e5, water, {3, 4}, {1, 1});
    for (const auto& amounts : {std::array<double, 2>{0.2, 0.4}, std::array<double, 2>{0.7, 1.1},
                               std::array<double, 2>{1e-15, 0.4}}) {
        const double acid = amounts[0], background = amounts[1];
        std::vector<double> m(9, 0.0);
        m[3] = acid; m[7] = background; m[4] = acid + background;
        const auto full = epcsaft::electrolyte_activity(nine, 298.15, 1e5,
            water, {3, 4}, {1, 1}, acid, mixed_reference, std::nullopt, m);
        const auto& logs = *full.state.log_fugacity_coefficient;
        const double expected = 0.5 * (epcsaft::require_component_property(logs[3]) - mixed_reference.log_fugacity[0]
            + epcsaft::require_component_property(logs[4]) - mixed_reference.log_fugacity[1])
            - std::log(1.0 + nine.parameters().molar_mass_kg_per_mol[0] * 2.0 * (acid + background));
        CHECK(std::abs(full.mean_ionic_log_molality_scale - expected) <= 1e-10);
        CHECK(full.osmotic_coefficient && std::holds_alternative<double>(*full.osmotic_coefficient));
        m[1] = 0.15;  // Occupied molecular solute in addition to both electrolyte pairs.
        const auto multi = epcsaft::electrolyte_activity(nine, 298.15, 1e5, water, {3, 4}, {1, 1}, acid, mixed_reference, std::nullopt, m);
        epcsaft::NeutralReferenceDeclaration declaration{0, {0,1,1,4,1,1,1,1,1}, 1.0, {}};
        std::vector<double> nu(9, 0.0); nu[3] = nu[4] = 1.0;
        std::vector<std::vector<epcsaft::NeutralReferenceReactionDiagnostic::Derivative>> columns;
        const auto diagnostics = epcsaft::reaction_references(nine, 298.15, 1e5, declaration, {nu}, {{}}, false, std::nullopt, 0, false, 0, &columns);
        CHECK(diagnostics[0].status == epcsaft::NeutralReferenceStatus::Available);
        const double mass = nine.parameters().molar_mass_kg_per_mol[0], total_m = std::accumulate(m.begin(), m.end(), 0.0);
        std::vector<double> x(9); x[0] = 1.0 / (1.0 + mass * total_m);
        for (std::size_t i = 1; i < 9; ++i) x[i] = x[0] * mass * m[i];
        const auto read = epcsaft::species_activities(x, &multi.state, columns[0], &nine.parameters(), &declaration,
            epcsaft::ReferenceSourceBasis::CommonMolalityInfiniteDilution);
        const double mean = 0.5 * (*read.log_activities[3] - std::log(m[3]) + *read.log_activities[4] - std::log(m[4]));
        CHECK(std::abs(mean - multi.mean_ionic_log_molality_scale) <= 1e-10);
        const auto solvent = epcsaft::solvent_activity(nine, 298.15, 1e5, x, mixed_reference, multi.state);
        const auto& pure = *mixed_reference.pure_solvent_state;
        const double standard_gap = *columns[0][0].value - std::log(pure.molar_density) - *pure.residual_chemical_potential_over_rt[0];
        CHECK(solvent.solvent_activity && std::abs(*read.log_activities[0] - std::log(*solvent.solvent_activity) + standard_gap) <= 1e-10);
        std::printf("activity-multisolute,acid=%.17g,background=%.17g,mean=%.17g,legacy=%.17g,water=%.17g,standardGap=%.17g\n",
            acid, background, mean, multi.mean_ionic_log_molality_scale, *read.log_activities[0], standard_gap);
        m[1] = 0;
        const auto rejects = [&](const std::vector<double>& declaration, const auto& ref) {
            try {
                epcsaft::electrolyte_activity(nine, 298.15, 1e5, water, {3, 4}, {1, 1}, acid, ref, std::nullopt, declaration);
            } catch (const std::invalid_argument&) { return true; }
            return false;
        };
        m[4] = acid;
        CHECK(rejects(m, mixed_reference));
        m[4] = acid + background; m[3] = -acid;
        CHECK(rejects(m, mixed_reference));
        m[3] = acid; m[0] = 1.0;
        CHECK(rejects(m, mixed_reference));
        m[0] = 0.0;
        auto wrong = mixed_reference;
        wrong.pressure = 101325.0;
        CHECK(rejects(m, wrong));
        m.pop_back();
        CHECK(rejects(m, mixed_reference));
    }
    const double gamma_coefficient = epcsaft::ELEMENTARY_CHARGE_C
        * epcsaft::ELEMENTARY_CHARGE_C
        / (8.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M
           * pure.bulk_relative_permittivity * epcsaft::BOLTZMANN_J_PER_K * 298.15)
        * std::sqrt(
            2.0 * epcsaft::ELEMENTARY_CHARGE_C * epcsaft::ELEMENTARY_CHARGE_C
            * epcsaft::AVOGADRO_PER_MOL
            * pure.molar_density * model.parameters().molar_mass_kg_per_mol[0]
            / (epcsaft::VACUUM_PERMITTIVITY_F_PER_M
               * pure.bulk_relative_permittivity * epcsaft::BOLTZMANN_J_PER_K
               * 298.15)
        );
    const auto activity_scale_inputs = recording_support::fields({
        reference_inputs, recording_support::field("ionic_strength", 1e-8)
    });
    recording_support::compare(
        "activity-scale", "mean_ionic_log_molality_scale", "1",
        activity.mean_ionic_log_molality_scale, -gamma_coefficient * std::sqrt(1e-8),
        activity_scale_inputs,
        "Relationship", "infinite-dilution", "activity scale matches Debye limit",
        "reference_or_one", 0.0, 1e-3,
        close(activity.mean_ionic_log_molality_scale, -gamma_coefficient * std::sqrt(1e-8), 1e-3),
        "m8-activity-scale", "ionic-strength", "\\log\\gamma", "T,P,\\mathbf{x},I",
        "\\log\\gamma= -\\Gamma\\sqrt I", "I");
}

void cache_invariance() {
    epcsaft::HelmholtzModel model(properties_support::aqueous_electrolyte());
    (void)model.evaluate(298.15, 1.0e-3, {0.98, 0.01, 0.01},
                         epcsaft::DerivativeOrder::Hessian);
    (void)model.evaluate(287.15, 1.0e-3, {1.0, 0.0, 0.0},
                         epcsaft::DerivativeOrder::Value);
    const auto a_again = model.evaluate(
        323.15, 1.5e-3, {1.35, 0.12, 0.03}, epcsaft::DerivativeOrder::Hessian);
    (void)model.evaluate(305.15, 1.1e-3, {1.067, 0.033, 0.0},
                         epcsaft::DerivativeOrder::Value);
    const auto b_again = model.evaluate(
        315.15, 1.25e-3, {1.25, 0.0, 0.0}, epcsaft::DerivativeOrder::Hessian);
    const auto fresh_a = epcsaft::HelmholtzModel(
        properties_support::aqueous_electrolyte()).evaluate(
            323.15, 1.5e-3, {1.35, 0.12, 0.03}, epcsaft::DerivativeOrder::Hessian);
    const auto fresh_b = epcsaft::HelmholtzModel(
        properties_support::aqueous_electrolyte()).evaluate(
            315.15, 1.25e-3, {1.25, 0.0, 0.0}, epcsaft::DerivativeOrder::Hessian);
    const std::array<const epcsaft::HelmholtzDerivatives*, 2> retained = {
        &a_again, &b_again};
    const std::array<const epcsaft::HelmholtzDerivatives*, 2> fresh = {
        &fresh_a, &fresh_b};
    const std::array<const char*, 2> labels = {"a-reuse", "b-eviction-revisit"};
    const std::array<std::string, 2> inputs = {{
        model_inputs(model.parameters(), 323.15, 1.5e-3, {1.35, 0.12, 0.03}),
        model_inputs(model.parameters(), 315.15, 1.25e-3, {1.25, 0.0, 0.0})}};
    for (std::size_t case_index = 0; case_index < retained.size(); ++case_index) {
        const auto& actual = *retained[case_index];
        const auto& reference = *fresh[case_index];
        const bool value_ok = close(actual.value, reference.value, 1e-14);
        recording_support::compare(
            (std::string("value/") + labels[case_index]).c_str(), "residual_helmholtz", "1",
            actual.value, reference.value, inputs[case_index], "Relationship", "cache-invariance",
            "retained/evicted signature matches a fresh model", "reference_or_one", 0.0, 1e-14,
            value_ok, "m9-value", labels[case_index], "a^{res}", "T,V,\\mathbf{n}",
            "a_{retained}=a_{fresh}", "state");
        CHECK(actual.gradient.size() == reference.gradient.size());
        for (std::size_t i = 0; i < actual.gradient.size(); ++i) {
            const bool gradient_ok = close(actual.gradient[i], reference.gradient[i], 1e-14);
            recording_support::compare(
                (std::string("gradient/") + labels[case_index] + "/" + std::to_string(i)).c_str(),
                "gradient", "1", actual.gradient[i], reference.gradient[i], inputs[case_index],
                "Relationship", "cache-invariance", "retained/evicted signature matches a fresh model",
                "reference_or_one", 0.0, 1e-14, gradient_ok, "m9-gradient", labels[case_index],
                "\\nabla a", "T,V,\\mathbf{n}", "\\nabla a_{retained}=\\nabla a_{fresh}",
                "coordinate");
        }
        CHECK(actual.hessian.size() == reference.hessian.size());
        CHECK(!actual.gradient.empty());
        const std::size_t width = actual.gradient.size();
        for (std::size_t i = 0; i < actual.hessian.size(); ++i) {
            const auto point = std::string(labels[case_index]) + "/"
                + std::to_string(i / width) + "-" + std::to_string(i % width);
            const bool hessian_ok = close(actual.hessian[i], reference.hessian[i], 1e-14);
            recording_support::compare(
                (std::string("hessian/") + point).c_str(), "hessian", "1", actual.hessian[i],
                reference.hessian[i], inputs[case_index], "Relationship", "cache-invariance",
                "retained/evicted signature matches a fresh model", "reference_or_one", 0.0, 1e-14,
                hessian_ok, "m9-hessian", point.c_str(), "H^{res}", "T,V,\\mathbf{n}",
                "H_{retained}=H_{fresh}", "row,col");
        }
    }
    // One tape serves both sides of SCREENING_SERIES_SWITCH: recorded at kappa*d = 0.087,
    // reused at kappa*d = 2.0 against a fresh recording there. A series branch frozen at
    // recording diverges past kappa*d = 1 (value -9.6e5 against -34.6). No recording on
    // reuse shows the comparison exercises the reused tape.
    epcsaft::HelmholtzModel crossing(properties_support::aqueous_electrolyte());
    (void)crossing.evaluate(298.15, 1.0e-3, {0.98, 0.01, 0.01}, epcsaft::DerivativeOrder::Hessian);
    crossing.reset_counters();
    const auto reused = crossing.evaluate(298.15, 4.0e-5, {0.8, 0.1, 0.1}, epcsaft::DerivativeOrder::Hessian);
    const auto recorded = epcsaft::HelmholtzModel(properties_support::aqueous_electrolyte()).evaluate(
        298.15, 4.0e-5, {0.8, 0.1, 0.1}, epcsaft::DerivativeOrder::Hessian);
    CHECK(crossing.counters().recordings == 0);
    CHECK(close(reused.value, recorded.value, 1e-10));
    for (std::size_t i = 0; i < recorded.gradient.size(); ++i) CHECK(close(reused.gradient[i], recorded.gradient[i], 1e-10));
    for (std::size_t i = 0; i < recorded.hessian.size(); ++i) CHECK(close(reused.hessian[i], recorded.hessian[i], 1e-10));
}

void core_c2_parameter_coordinates() {
    auto rejects = [](epcsaft::Parameters parameters,
                      epcsaft::ActiveParameter active) {
        try {
            (void)epcsaft::HelmholtzModel(std::move(parameters), {active});
        } catch (const std::invalid_argument&) {
            return true;
        }
        return false;
    };
    auto binary = properties_support::neutral_binary();
    const auto kij = epcsaft::ParameterFamily::Kij;
    const auto lij = epcsaft::ParameterFamily::Lij;
    CHECK(rejects(binary, {epcsaft::ParameterFamily::SegmentCount, {0, 1}}));
    CHECK(rejects(binary, {epcsaft::ParameterFamily::SigmaConstant, {0, 1}}));
    CHECK(rejects(binary, {epcsaft::ParameterFamily::DispersionEnergyOverK, {0, 1}}));
    constexpr double scalar_temperature = 300.0;
    constexpr double scalar_volume = 1.0 / 1.2e4;
    const std::vector<double> scalar_amounts = {0.4, 0.6};
    const auto fixed_scalar = epcsaft::HelmholtzModel(binary).evaluate(
        scalar_temperature, scalar_volume, scalar_amounts, epcsaft::DerivativeOrder::Value);
    const std::vector<epcsaft::ParameterFamily> scalar_families = {
        epcsaft::ParameterFamily::SegmentCount,
        epcsaft::ParameterFamily::SigmaConstant,
        epcsaft::ParameterFamily::DispersionEnergyOverK
    };
    const std::array<const char*, 3> scalar_names = {
        "segment-count", "sigma", "dispersion-energy-over-k"
    };
    for (std::size_t i = 0; i < scalar_families.size(); ++i) {
        const auto family = scalar_families[i];
        const auto active_value = epcsaft::HelmholtzModel(
            binary, {{family, {0}}}).evaluate(
                scalar_temperature, scalar_volume, scalar_amounts,
                epcsaft::DerivativeOrder::Value).value;
        CHECK(close(active_value, fixed_scalar.value, 1e-12));
        auto overwritten = binary;
        switch (family) {
        case epcsaft::ParameterFamily::SegmentCount:
            overwritten.segment_count[1] = overwritten.segment_count[0];
            break;
        case epcsaft::ParameterFamily::SigmaConstant:
            overwritten.sigma_constant[1] = overwritten.sigma_constant[0];
            break;
        case epcsaft::ParameterFamily::DispersionEnergyOverK:
            overwritten.dispersion_energy_over_k[1] =
                overwritten.dispersion_energy_over_k[0];
            break;
        default:
            break;
        }
        const auto overwritten_value = epcsaft::HelmholtzModel(overwritten).evaluate(
            scalar_temperature, scalar_volume, scalar_amounts,
            epcsaft::DerivativeOrder::Value).value;
        CHECK(!close(overwritten_value, fixed_scalar.value, 1e-12));
        recording_support::compare(
            (std::string("c2-scalar-value/") + scalar_names[i]).c_str(),
            "residual_helmholtz", "1", active_value, fixed_scalar.value,
            recording_support::fields({
                model_inputs(binary, scalar_temperature, scalar_volume, scalar_amounts),
                recording_support::field("family", scalar_names[i]),
                recording_support::field("active_index", std::size_t{0})
            }), "Identity", "single-slot-coordinate", "fixed value is preserved",
            "reference_or_one", 0.0, 1e-12,
            close(active_value, fixed_scalar.value, 1e-12),
            "c2-scalar-values", scalar_names[i], "a^{res}", "T,V,\\mathbf{n}",
            "a_{active}=a_{fixed}", "scalar");
    }
    CHECK(rejects(binary, {kij, {1}}));
    CHECK(rejects(binary, {kij, {0}}));
    CHECK(rejects(binary, {kij, {1, 1}}));
    CHECK(rejects(binary, {kij, {1, 3}}));
    CHECK(rejects(binary, {lij, {1}}));
    CHECK(rejects(binary, {lij, {0}}));
    (void)epcsaft::HelmholtzModel(binary, {{kij, {1, 2}}, {lij, {1, 2}}});

    auto associating = properties_support::associating_binary();
    CHECK(!rejects(associating, {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}}));
    CHECK(rejects(associating, {epcsaft::ParameterFamily::AssociationEnergyOverK, {1}}));
    CHECK(rejects(associating, {epcsaft::ParameterFamily::AssociationEnergyOverK, {0}}));
    auto self = associating;
    self.association_active[0] = 1;
    self.association_energy_over_k[0] = 1200.0;
    self.association_volume[0] = 0.1;
    CHECK(!rejects(self, {epcsaft::ParameterFamily::AssociationEnergyOverK, {0}}));

    const auto association_active = std::vector<epcsaft::ActiveParameter>{
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}}
    };
    epcsaft::HelmholtzModel association_model(associating, association_active);
    const auto association_base = association_model.evaluate(
        376.88, 1.0e-4, {0.4, 0.6}, epcsaft::DerivativeOrder::Gradient);
    auto association_shifted = [&](double delta) {
        auto p = associating;
        p.association_energy_over_k[1] += delta;
        p.association_energy_over_k[4] += delta;
        return epcsaft::HelmholtzModel(std::move(p), association_active).evaluate(
            376.88, 1.0e-4, {0.4, 0.6}, epcsaft::DerivativeOrder::Value).value;
    };
    constexpr double association_h = 1e-5;
    const double association_fd =
        (association_shifted(association_h) - association_shifted(-association_h))
        / (2.0 * association_h);
    const double association_error = std::abs(
        association_base.gradient.back() - association_fd)
        / std::max(1.0, std::abs(association_fd));
    CHECK(association_error <= 5e-6);

    constexpr double temperature = 300.0;
    constexpr double volume = 1.0 / 1.2e4;
    const std::vector<double> amounts = {0.4, 0.6};
    binary.l_ij[1] = binary.l_ij[2] = 0.2;
    const auto active = std::vector<epcsaft::ActiveParameter>{{kij, {1, 2}}, {lij, {1, 2}}};
    epcsaft::HelmholtzModel model(binary, active);
    const auto base = model.evaluate(temperature, volume, amounts, epcsaft::DerivativeOrder::Hessian);
    const auto fixed = epcsaft::HelmholtzModel(binary).evaluate(
        temperature, volume, amounts, epcsaft::DerivativeOrder::Value);
    CHECK(close(base.value, fixed.value, 1e-12));
    const double h = 1e-5;
    auto shifted = [&](double dk, double dl, epcsaft::DerivativeOrder order) {
        auto p = binary;
        p.k_ij[1] += dk; p.k_ij[2] += dk;
        p.l_ij[1] += dl; p.l_ij[2] += dl;
        return epcsaft::HelmholtzModel(std::move(p), active).evaluate(
            temperature, volume, amounts, order);
    };
    const auto plus_k = shifted(h, 0.0, epcsaft::DerivativeOrder::Value);
    const auto minus_k = shifted(-h, 0.0, epcsaft::DerivativeOrder::Value);
    const auto plus_l = shifted(0.0, h, epcsaft::DerivativeOrder::Value);
    const auto minus_l = shifted(0.0, -h, epcsaft::DerivativeOrder::Value);
    const double gradient_k_fd = (plus_k.value - minus_k.value) / (2.0 * h);
    const double gradient_l_fd = (plus_l.value - minus_l.value) / (2.0 * h);
    CHECK(close(base.gradient[base.gradient.size() - 2], gradient_k_fd, 5e-6));
    CHECK(close(base.gradient.back(), gradient_l_fd, 5e-6));
    CHECK(std::abs(base.gradient[base.gradient.size() - 2]) > 1e-12);
    CHECK(std::abs(base.gradient.back()) > 1e-12);

    const auto plus_direction = shifted(0.7 * h, -0.4 * h, epcsaft::DerivativeOrder::Gradient);
    const auto minus_direction = shifted(-0.7 * h, 0.4 * h, epcsaft::DerivativeOrder::Gradient);
    const auto width = base.gradient.size();
    double error = 0.0;
    double scaling_error = 0.0;
    for (std::size_t row = 0; row < width; ++row) {
        const double fd = (plus_direction.gradient[row] - minus_direction.gradient[row]) / (2.0 * h);
        const double action = base.hessian[row * width + width - 2] * 0.7
            + base.hessian[row * width + width - 1] * -0.4;
        error = std::max(error, std::abs(fd - action) / std::max(1.0, std::abs(fd)));
    }
    CHECK(error <= 5e-6);

    const double scales[2] = {2.0, 3.0};
    for (std::size_t i = 0; i < 2; ++i) {
        const auto source_plus = shifted(i == 0 ? scales[i] * h : 0.0,
                                          i == 1 ? scales[i] * h : 0.0,
                                          epcsaft::DerivativeOrder::Value);
        const auto source_minus = shifted(i == 0 ? -scales[i] * h : 0.0,
                                           i == 1 ? -scales[i] * h : 0.0,
                                           epcsaft::DerivativeOrder::Value);
        const double source_gradient = (source_plus.value - source_minus.value) / (2.0 * h);
        const double gradient_error = std::abs(source_gradient
            - scales[i] * base.gradient[width - 2 + i])
            / std::max(1.0, std::abs(source_gradient));
        scaling_error = std::max(scaling_error, gradient_error);
        CHECK(gradient_error <= 5e-6);
        for (std::size_t j = 0; j < 2; ++j) {
            auto source_plus_gradient = shifted(
                i == 0 ? scales[i] * h : 0.0, i == 1 ? scales[i] * h : 0.0,
                epcsaft::DerivativeOrder::Gradient);
            auto source_minus_gradient = shifted(
                i == 0 ? -scales[i] * h : 0.0, i == 1 ? -scales[i] * h : 0.0,
                epcsaft::DerivativeOrder::Gradient);
            const double source_hessian =
                (source_plus_gradient.gradient[width - 2 + j]
                 - source_minus_gradient.gradient[width - 2 + j])
                * scales[j] / (2.0 * h);
            const double expected = base.hessian[(width - 2 + i) * width + width - 2 + j]
                * scales[i] * scales[j];
            const double hessian_error = std::abs(source_hessian - expected)
                / std::max(1.0, std::abs(source_hessian));
            scaling_error = std::max(scaling_error, hessian_error);
            CHECK(hessian_error <= 5e-6);
        }
    }
    recording_support::compare("c2-falsifier", "relative_error", "1", error, 0.0,
        recording_support::fields({model_inputs(binary, temperature, volume, amounts),
                                   recording_support::field("h", h),
                                   recording_support::field("v", std::vector<double>{0.7, -0.4})}), "Relationship",
        "centered-difference", "C2 Hessian action", "reference_or_one", 0.0, 5e-6,
        error <= 5e-6, "c2-falsifier", "H*v", "relative_error", "H*v", "T,V,n",
        "H*v=finite_difference");
    recording_support::compare("c2-association-gradient", "relative_error", "1",
        association_error, 0.0,
        recording_support::fields({model_inputs(associating, 376.88, 1.0e-4, {0.4, 0.6}),
                                   recording_support::field("h", association_h),
                                   recording_support::field("pair", "1,4")}),
        "Relationship", "centered-difference", "association gradient", "reference_or_one",
        0.0, 5e-6, association_error <= 5e-6, "c2-association-gradient", "distinct-site-pair",
        "relative_error", "dA/dtheta", "T,V,n", "dA/dtheta=finite_difference");
    recording_support::compare("c2-scaling", "relative_error", "1", scaling_error, 0.0,
        recording_support::fields({recording_support::field("h", h),
                                   recording_support::field("S", std::vector<double>{2.0, 3.0})}),
        "Relationship", "scaled-coordinate", "source gradient/Hessian scaling",
        "reference_or_one", 0.0, 5e-6, scaling_error <= 5e-6, "c2-scaling", "S=diag(2,3)",
        "relative_error", "Sg,S H S", "T,V,n", "source=S*native");
}

// #423: a sole cold root below the stored negative-slope sample is vapor-side and is no liquid reference.
void liquid_reference_branch_evidence() {
    const auto refused = [](const std::string& message) { return message.find("below the negative-slope sample") != std::string::npos; };
    const epcsaft::HelmholtzModel neutral(properties_support::neutral_binary());
    const auto ethane = [&](double T) {
        return epcsaft::reaction_reference(neutral, T, 100.0, {1, {1.0, 0.0}, 1.0, {}}, {1.0, -1.0}, {});
    };
    const auto gas = ethane(300.0);  // sole root 0.040091082 mol/m3, sample 7087.881992 mol/m3
    CHECK(gas.status == epcsaft::NeutralReferenceStatus::DensityUnavailable
          && gas.message.find("solvent species 1") != std::string::npos && refused(gas.message));
    CHECK(ethane(310.0).status == epcsaft::NeutralReferenceStatus::Available);  // no sample: sole root kept
    const auto pure = epcsaft::reaction_reference(neutral, 300.0, 100.0,
        {0, {0.0, 1.0}, 1.0, {epcsaft::SpeciesReference::PureLiquid, epcsaft::SpeciesReference::PureLiquid}}, {1.0, -1.0}, {});
    CHECK(pure.status == epcsaft::NeutralReferenceStatus::DensityUnavailable
          && pure.message.find("pure reference species 1: liquid branch not established by the finite search") != std::string::npos);
    // Delivered-review case: 2A <=> B, identical ethane EOS species, both referenced to pure gas.
    auto identical_ethane = properties_support::parameters(2);
    identical_ethane.segment_count = {1.6069, 1.6069};
    identical_ethane.sigma_constant = {3.5206, 3.5206};
    identical_ethane.dispersion_energy_over_k = {191.42, 191.42};
    identical_ethane.molar_mass_kg_per_mol = {0.03007, 0.03007};
    const auto pure_gas = epcsaft::reaction_reference(epcsaft::HelmholtzModel(identical_ethane), 300.0, 100.0,
        {0, {0.0, 1.0}, 1.0, {epcsaft::SpeciesReference::PureGas, epcsaft::SpeciesReference::PureGas}}, {-2.0, 1.0},
        {std::nullopt, epcsaft::ReferenceSourceBasis::RawMoleFractionInfiniteDilution});
    CHECK(pure_gas.status == epcsaft::NeutralReferenceStatus::Available && pure_gas.terminal_contraction
          && close(*pure_gas.terminal_contraction, 3.2166161799675423, 1e-10));
    epcsaft::State extreme; extreme.molar_density = 1; extreme.residual_chemical_potential_over_rt = {1000,-1000};
    const auto range = epcsaft::species_activities({0.5,0.5}, &extreme);
    CHECK(range.log_activities[0] && range.log_activities[1] && !range.activities[0] && !range.activities[1]);
    CHECK(!range.causes[0].empty() && !range.causes[1].empty());
    for (const auto choice : {epcsaft::SpeciesReference::PureGas, epcsaft::SpeciesReference::PureLiquid}) {
        const epcsaft::HelmholtzModel fluid(identical_ethane);
        const double pressure = choice == epcsaft::SpeciesReference::PureGas ? 100.0 : 1e7;
        const epcsaft::NeutralReferenceDeclaration declaration{0, {0,1}, 1.0, {choice,choice}};
        const epcsaft::ReactionReference source{pressure, epcsaft::ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
        std::vector<std::vector<epcsaft::NeutralReferenceReactionDiagnostic::Derivative>> columns;
        const auto contracted = epcsaft::reaction_references(fluid, 300, pressure, declaration, {{-2,1}}, {source}, false, std::nullopt, 0, false, 0, &columns);
        const auto finite = fluid.state_at_pressure(300, pressure, {0.4,0.6}, choice == epcsaft::SpeciesReference::PureGas ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Liquid);
        const auto read = epcsaft::species_activities({0.4,0.6}, &finite, columns[0]);
        std::printf("activity-pure,choice=%d,P=%.17g,cause=%s\n", int(choice), pressure, contracted[0].message.c_str());
        CHECK(contracted[0].terminal_contraction && read.log_activities[0] && read.log_activities[1]);
        const double expected = -2 * (std::log(0.4 * finite.molar_density) + *finite.residual_chemical_potential_over_rt[0])
            + std::log(0.6 * finite.molar_density) + *finite.residual_chemical_potential_over_rt[1] - *contracted[0].terminal_contraction;
        CHECK(std::abs(-2 * *read.log_activities[0] + *read.log_activities[1] - expected) <= 1e-10);
    }
    const epcsaft::HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    const auto message = [](const auto& call) {
        try { call(); } catch (const std::runtime_error& error) { return std::string(error.what()); }
        return std::string();
    };
    // 650 K, 100 Pa: sole root 0.018503460 mol/m3 below sample 17653.65 mol/m3 gave gamma = 1.000000196209869.
    CHECK(refused(message([&] { epcsaft::infinite_dilution_reference(aqueous, 650.0, 100.0, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}); })));
    CHECK(refused(message([&] { epcsaft::electrolyte_activity(aqueous, 650.0, 100.0, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}, 1e-6, [] { epcsaft::InfiniteDilutionReference r; r.ion_indices = {1,2}; return r; }()); })));
    // 630 K keeps its two-root liquid reference (main values).
    const auto reference = epcsaft::infinite_dilution_reference(aqueous, 630.0, 100.0, {1.0, 0.0, 0.0}, {1, 2}, {1, 1});
    const auto activity = epcsaft::electrolyte_activity(aqueous, 630.0, 100.0, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}, 1e-6, reference);
    CHECK(close(reference.pure_solvent_state->molar_density, 34616.636583442916, 1e-10));
    // gamma before #594 0.9996968293218973 times exp(-(10/9) mean 1e-12 -> 1e-14 rung change) = exp(-3.0408e-7).
    CHECK(close(activity.mean_ionic_molality_scale, 0.9996965253318902, 1e-10));
}

}  // namespace properties_checks
