#if defined(__linux__)
#include <quadmath.h>
#endif
#include <iostream>
// P1, F1, M2 and M3: residual properties, caloric ladders, derivatives and limits.
#include "properties_checks_shared.hpp"

namespace properties_checks {

struct P1Case final {
    const char* name;
    epcsaft::Parameters parameters;
    std::vector<double> composition;
    double temperature, pressure;
    epcsaft::DensityBranch branch;
    const char* branch_name;
};

std::array<P1Case, 5> p1_cases() {
    return {{
        {"ethane-vapor", properties_support::ethane(), {1.0}, 200.0, properties_support::ethane_two_root_pressure,
         epcsaft::DensityBranch::Vapor, "vapor"},
        {"ethane-liquid", properties_support::ethane(), {1.0}, 200.0, properties_support::ethane_two_root_pressure,
         epcsaft::DensityBranch::Liquid, "liquid"},
        {"associating-binary", properties_support::associating_binary(), {0.1385, 0.8615},
         313.15, 1.0e5, epcsaft::DensityBranch::Liquid, "liquid"},
        {"aqueous-electrolyte", properties_support::aqueous_electrolyte(), {0.98, 0.01, 0.01},
         298.15, 1.0e5, epcsaft::DensityBranch::Liquid, "liquid"},
        {"nine-species", properties_support::nine_species(),
         {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
         313.15, 1.0e5, epcsaft::DensityBranch::Liquid, "liquid"},
    }};
}

std::string p1_pressure_record(const char* prefix, const epcsaft::State& state,
                               double target_pressure, const char* requested_branch) {
    const auto* convergence = state.density_convergence.has_value()
        ? &*state.density_convergence : nullptr;
    const bool accepted = convergence != nullptr
        && convergence->criterion != nullptr
        && convergence->stop == epcsaft::core::Stop::converged
        && std::isfinite(convergence->residual)
        && convergence->residual <= convergence->criterion->value;
    CHECK(accepted);
    if (!accepted) return {};
    const auto name = [prefix](const char* suffix) {
        return std::string(prefix) + "_" + suffix;
    };
    return recording_support::fields({
        recording_support::field(name("requested_branch").c_str(), requested_branch),
        recording_support::field(name("target_pressure").c_str(), target_pressure),
        recording_support::field(name("achieved_pressure").c_str(), state.pressure),
        recording_support::field(name("density").c_str(), state.molar_density),
        recording_support::field(name("convergence_stop").c_str(),
                                 epcsaft::core::stop_name(convergence->stop)),
        recording_support::field(name("convergence_residual").c_str(), convergence->residual),
        recording_support::field(name("convergence_criterion").c_str(),
                                 convergence->criterion->value)
    });
}

void p1_compare(const std::string& id, const char* quantity, const char* unit,
                double model, double reference, const std::string& inputs, double tolerance,
                const char* group) {
    const bool cv = std::string(group) == "p1-cv";
    const bool cp = std::string(group) == "p1-cp";
    recording_support::compare(id, quantity, unit, model, reference, recording_support::Rule::ReferenceOrOne, 0.0,
        tolerance, cv ? "centered_temperature_derivative"
           : cp ? "independent_fixed_pressure_derivative" : "dilute_density_refinement",
        inputs);
}

void p1_observe(const std::string& id, const char* quantity, const char* unit, double value,
                const std::string& inputs, const char* method) {
    recording_support::observe(id, quantity, unit, value, std::nullopt, method, inputs);
}

void residual_internal_energy_identity_ladder() {
    // F1: the state's u, s and h against a centered temperature difference of a_res,
    // which is independent of the Engine's exact a_T and its caloric formulas.
    epcsaft::HelmholtzModel model(properties_support::ethane());
    constexpr double temperature = 300.0;
    constexpr double density = 1.0e4;
    constexpr double gas_constant = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const std::vector<double> composition = {1.0};
    const auto state = model.state(temperature, density, composition);
    const auto inputs = model_inputs(
        model.parameters(), temperature, 1.0 / density, composition
    );
    double a_t = 0.0;
    for (const double delta_temperature : {0.1, 0.03, 0.01, 0.003}) {
        const double plus = model.state(
            temperature + delta_temperature, density, composition
        ).residual_helmholtz;
        const double minus = model.state(
            temperature - delta_temperature, density, composition
        ).residual_helmholtz;
        a_t = (plus - minus) / (2.0 * delta_temperature);
        const double centered = -gas_constant * temperature * temperature * a_t;
        char point[64];
        std::snprintf(point, sizeof(point), "delta-%g", delta_temperature);
        char comparison_id[96];
        std::snprintf(comparison_id, sizeof(comparison_id), "temperature-ladder/%s", point);
        const auto point_inputs = recording_support::fields({
            inputs,
            recording_support::field("delta_T", delta_temperature),
            recording_support::field("T_plus", temperature + delta_temperature),
            recording_support::field("T_minus", temperature - delta_temperature),
        });
        recording_support::compare(comparison_id, "residual_internal_energy", "J/mol", state.residual_internal_energy,
            centered, recording_support::Rule::ReferenceOrOne, 0.0, 2.0e-7, "centered_temperature_derivative",
            point_inputs);
    }
    // a_T from the 0.003 K step: s = -R(a + T a_T) and h = RT(Z - 1) - RT^2 a_T.
    const double entropy = -gas_constant * (state.residual_helmholtz + temperature * a_t);
    const double enthalpy = gas_constant * temperature * (state.compressibility_factor - 1.0 - temperature * a_t);
    // The 0.003 K difference agrees to about 2e-10 (truncation-limited); 1e-8 keeps a 50x margin.
    recording_support::compare("residual_entropy", "residual_entropy", "J/mol/K", state.residual_entropy, entropy,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-8, "centered_temperature_derivative", inputs);
    recording_support::compare("residual_enthalpy", "residual_enthalpy", "J/mol", state.residual_enthalpy, enthalpy,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-8, "centered_temperature_derivative", inputs);
}

void residual_properties_and_caloric_ladders() {
    constexpr double gas_constant = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    constexpr std::array<double, 4> deltas = {0.1, 0.03, 0.01, 0.003};
    const auto cases = p1_cases();

    for (const auto& item : cases) {
        epcsaft::HelmholtzModel model(item.parameters);
        const auto pressure_state = model.state_at_pressure(item.temperature, item.pressure, item.composition, item.branch);
        const auto density_state = model.state(item.temperature, pressure_state.molar_density, item.composition);
        const double D = pressure_state.pressure_density_derivative / (gas_constant * item.temperature);
        CHECK(std::isfinite(D) && D > 0.0);
        CHECK(std::isfinite(density_state.residual_isochoric_heat_capacity));
        CHECK(pressure_state.residual_isobaric_heat_capacity.has_value());
        const double cp = *pressure_state.residual_isobaric_heat_capacity;
        CHECK(std::isfinite(cp));
        const auto inputs = recording_support::fields({model_inputs(item.parameters, item.temperature, 1.0 / pressure_state.molar_density, item.composition), p1_pressure_record("central", pressure_state, item.pressure, item.branch_name), recording_support::field("D", D)});
        double cv_coarse = 0.0;
        double cp_coarse = 0.0;
        for (std::size_t step = 0; step < deltas.size(); ++step) {
            const double delta = deltas[step];
            const auto cv_plus = model.state(item.temperature + delta, pressure_state.molar_density, item.composition);
            const auto cv_minus = model.state(item.temperature - delta, pressure_state.molar_density, item.composition);
            const double cv_reference = (cv_plus.residual_internal_energy - cv_minus.residual_internal_energy) / (2.0 * delta);
            const double tolerance = delta <= 0.01 ? 1.0e-8 : 1.0e-5;
            const std::string point = std::string(item.name) + "/" + recording_support::number(delta);
            p1_compare("cv-ladder/" + point, "residual_isochoric_heat_capacity", "J/mol/K", density_state.residual_isochoric_heat_capacity, cv_reference, recording_support::fields({inputs, recording_support::field("route", "fixed-density"), recording_support::field("delta_T", delta), recording_support::field("P_plus", cv_plus.pressure), recording_support::field("P_minus", cv_minus.pressure)}), tolerance, "p1-cv");
            const auto cp_plus = model.state_at_pressure(item.temperature + delta, item.pressure, item.composition, item.branch);
            const auto cp_minus = model.state_at_pressure(item.temperature - delta, item.pressure, item.composition, item.branch);
            CHECK(close(cp_plus.pressure, item.pressure, 1.0e-9));
            CHECK(close(cp_minus.pressure, item.pressure, 1.0e-9));
            const double cp_reference = (cp_plus.residual_enthalpy - cp_minus.residual_enthalpy) / (2.0 * delta);
            p1_compare("cp-ladder/" + point, "residual_isobaric_heat_capacity", "J/mol/K", cp, cp_reference, recording_support::fields({inputs, p1_pressure_record("plus", cp_plus, item.pressure, item.branch_name), p1_pressure_record("minus", cp_minus, item.pressure, item.branch_name), recording_support::field("route", "fixed-pressure"), recording_support::field("delta_T", delta)}), tolerance, "p1-cp");
            if (step == 2) {
                cv_coarse = cv_reference;
                cp_coarse = cp_reference;
            } else if (step == 3) {
                CHECK(close(cv_coarse, cv_reference, 2.0e-8));
                CHECK(close(cp_coarse, cp_reference, 2.0e-8));
            }
        }
    }

    const epcsaft::HelmholtzModel ethane(properties_support::ethane());
    for (double rho : {1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 1e-20, 1e-50, 1e-100}) {
        const auto state = ethane.state(300, rho, {1});
        const auto derivatives = ethane.evaluate(300, 1 / rho, {1}, epcsaft::DerivativeOrder::Hessian);
        const double reference = gas_constant * 300 * (-300 * derivatives.gradient[0] - derivatives.gradient[1] / rho);
        const auto inputs = model_inputs(ethane.parameters(), 300, 1 / rho, {1});
        recording_support::compare("dilute-enthalpy/" + recording_support::number(rho), "residual_enthalpy", "J/mol",
            state.residual_enthalpy, reference, recording_support::Rule::ReferenceScale, 0, 1e-10, "residual_enthalpy",
            inputs);
    }

    for (double rho : {1e-12, 1e-15}) {
        const auto state = ethane.state(300, rho, {1});
        const auto action = ethane.density_actions(300, rho, {1}, {epcsaft::PropertyObservable::ResidualEnthalpy}, {{1, 0, {0}, {}}});
        CHECK(action.results[0].value && std::abs(*action.results[0].value / state.residual_enthalpy - 1) <= 1e-10);
    }

    epcsaft::HelmholtzModel dilute_model(properties_support::neutral_binary());
    const std::vector<double> dilute_composition = {0.4, 0.6};
    const std::array<double, 4> dilute_densities = {10.0, 1.0, 0.1, 0.01};
    std::array<epcsaft::State, 4> dilute_states;
    for (std::size_t index = 0; index < dilute_densities.size(); ++index) {
        dilute_states[index] = dilute_model.state(300.0, dilute_densities[index], dilute_composition);
        const auto& state = dilute_states[index];
        CHECK(finite_state_without_cp(state));
        CHECK(state.residual_isobaric_heat_capacity.has_value());
        const double D = state.pressure_density_derivative / (gas_constant * 300.0);
        p1_observe("dilute-state/" + recording_support::number(dilute_densities[index]), "residual_helmholtz_over_rho", "m3/mol", state.residual_helmholtz / dilute_densities[index], recording_support::fields({model_inputs(dilute_model.parameters(), 300.0, 1.0 / dilute_densities[index], dilute_composition), recording_support::field("D", D), recording_support::field("fields_over_rho", std::vector<double>{state.residual_entropy / dilute_densities[index], state.residual_internal_energy / dilute_densities[index], state.residual_enthalpy / dilute_densities[index], state.residual_isochoric_heat_capacity / dilute_densities[index], *state.residual_isobaric_heat_capacity / dilute_densities[index]})}), "dilute_ladder");
    }
    const auto& coarse = dilute_states[2];
    const auto& fine = dilute_states[3];
    const auto dilute_inputs = recording_support::fields({model_inputs(dilute_model.parameters(), 300.0, 1.0 / dilute_densities[2], dilute_composition), recording_support::field("rho_coarse", dilute_densities[2]), recording_support::field("rho_fine", dilute_densities[3])});
    const std::array<std::pair<const char*, double epcsaft::State::*>, 4> dilute_fields = {{{"residual_entropy", &epcsaft::State::residual_entropy}, {"residual_internal_energy", &epcsaft::State::residual_internal_energy}, {"residual_enthalpy", &epcsaft::State::residual_enthalpy}, {"residual_isochoric_heat_capacity", &epcsaft::State::residual_isochoric_heat_capacity}}};
    for (const auto& [quantity, field] : dilute_fields) {
        p1_compare(std::string("dilute-limit/") + quantity, quantity,
                   std::string(quantity) == "residual_entropy"
                       || std::string(quantity) == "residual_isochoric_heat_capacity"
                       ? "J m3/(mol2 K)" : "J m3/mol2",
                   coarse.*field / dilute_densities[2], fine.*field / dilute_densities[3],
                   dilute_inputs, 5.0e-6, "p1-dilute-limit");
    }
    p1_compare("dilute-limit/residual_isobaric_heat_capacity", "residual_isobaric_heat_capacity", "J m3/(mol2 K)", *coarse.residual_isobaric_heat_capacity / dilute_densities[2], *fine.residual_isobaric_heat_capacity / dilute_densities[3], dilute_inputs, 5.0e-6, "p1-dilute-limit");

    epcsaft::HelmholtzModel conditioning_model(properties_support::ethane());
    const auto conditioned = conditioning_model.state(280.0, 3000.0, {1.0});
    const auto unavailable = conditioning_model.state(280.0, 4000.0, {1.0});
    const double conditioned_D = conditioned.pressure_density_derivative / (gas_constant * 280.0);
    const double unavailable_D = unavailable.pressure_density_derivative / (gas_constant * 280.0);
    CHECK(std::isfinite(conditioned_D) && conditioned_D > 0.0);
    CHECK(conditioned.residual_isobaric_heat_capacity.has_value());
    p1_observe("conditioning/near-zero-positive-D", "D", "1", conditioned_D, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 3000.0, {1.0}), recording_support::field("D", conditioned_D), recording_support::field("cp", *conditioned.residual_isobaric_heat_capacity)}), "conditioning_stress");
    p1_observe("conditioning/near-zero-positive-cp", "residual_isobaric_heat_capacity", "J/mol/K", *conditioned.residual_isobaric_heat_capacity, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 3000.0, {1.0}), recording_support::field("D", conditioned_D)}), "conditioning_stress");
    CHECK(std::isfinite(unavailable_D) && unavailable_D < 0.0);
    const bool unavailable_cp = !unavailable.residual_isobaric_heat_capacity.has_value();
    recording_support::compare_custom("conditioning/negative-D-cp-unavailable",
        "residual_isobaric_heat_capacity_available", "bool", unavailable_cp ? 0.0 : 1.0, 0.0, "absolute",
        unavailable_cp, "mechanical_conditioning",
        recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 4000.0, {1.0}), recording_support::field("D", unavailable_D), recording_support::field("unavailable_reason", "nonpositive mechanical denominator")}));
    const bool finite_other_outputs = finite_state_without_cp(unavailable);
    recording_support::compare_custom("conditioning/negative-D-other-fields-finite", "residual_properties_finite",
        "bool", finite_other_outputs ? 1.0 : 0.0, 1.0, "absolute", finite_other_outputs, "partial_availability",
        recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 4000.0, {1.0}), recording_support::field("D", unavailable_D)}));
}
void derivatives_and_limits() {
    {
#if defined(__linux__)
    // Gross2001 A.4/A.6/A.7 and the exact one-site mass-action root, in binary128.
    auto dilute = properties_support::ethane();
    dilute.segment_count[0] = 1.6069; dilute.sigma_constant[0] = 3.5206;
    dilute.dispersion_energy_over_k[0] = 191.42;
    dilute.site_component = {0}; dilute.site_multiplicity = {1}; dilute.association_active = {1};
    dilute.association_energy_over_k = {1500}; dilute.association_volume = {.03};
    const std::vector<double> pure_dilute{1.0};
    const auto dsigma = epcsaft::segment_sigma(300.0, dilute);
    const auto dd = epcsaft::segment_diameter(300.0, dilute, dsigma);
    using Wide = __float128;
    CHECK(FLT128_MANT_DIG == 113);
    for (double rho : {1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 1e-20, 1e-50, 1e-100, 1e-150, 1e-300}) {
        const auto delta = epcsaft::association_delta(300.0, rho, pure_dilute, dilute, dsigma, dd);
        const auto sites = epcsaft::solve_site_fractions(300, rho, pure_dilute, dilute, dsigma, dd);
        const auto actual = epcsaft::residual_helmholtz(300.0, rho, pure_dilute, dilute, dsigma, dd,
            sites, epcsaft::association_products(rho, pure_dilute, dilute, sites, delta));
        const Wide eta = Wide(epcsaft::PI) / 6 * Wide(rho) * Wide(epcsaft::NUMBER_DENSITY_PER_ANGSTROM3)
            * Wide(dilute.segment_count[0]) * Wide(dd[0]) * Wide(dd[0]) * Wide(dd[0]);
        const Wide u = 1 - eta;
        const Wide hc = Wide(dilute.segment_count[0]) * (4 * eta - 3 * eta * eta) / (u * u)
            - Wide(dilute.segment_count[0] - 1) * (log1pq(-eta / 2) - 3 * log1pq(-eta));
        const Wide K = Wide(rho) * Wide(epcsaft::NUMBER_DENSITY_PER_ANGSTROM3) * Wide(delta[0]);
        const Wide b = 2 * K / (1 + 2 * K + sqrtq(1 + 4 * K));
        const Wide assoc = log1pq(-b) + b / 2;
        const double hc_error = static_cast<double>(fabsq((Wide(actual.hard_chain) - hc) / hc));
        const double assoc_error = static_cast<double>(fabsq((Wide(actual.association) - assoc) / assoc));
        std::cout << std::setprecision(17) << "dilute," << rho << ',' << actual.hard_chain << ',' << static_cast<double>(hc)
            << ',' << hc_error << ',' << actual.association << ',' << static_cast<double>(assoc) << ',' << assoc_error << '\n';
        CHECK(hc_error <= 1e-10 && assoc_error <= 1e-10);
    }
#endif
    }
    for (const bool derived : {false, true}) {
        const epcsaft::HelmholtzModel dilute(properties_support::ethane());
        for (double rho : {1e-100, 1e-150, 1e-300}) {
            try {
                const auto state = dilute.state(300, rho, {1}, derived);
                CHECK(std::isfinite(state.pressure) && std::isfinite(state.residual_enthalpy));
                CHECK(std::abs(state.pressure / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300 * rho) - 1) <= 1e-10);
                CHECK(std::isfinite(state.pressure_density_derivative));
                std::cout << std::setprecision(17) << "state," << derived << ',' << rho << ',' << state.pressure << ','
                    << state.residual_enthalpy << ',' << state.pressure_density_derivative << '\n';
                CHECK(state.density_parameter_derivative.empty());
            } catch (const std::runtime_error& error) {
                std::cout << "refusal," << derived << ',' << rho << ',' << error.what() << '\n';
                CHECK(rho < 1e-100 && std::string(error.what()).find("nonfinite") != std::string::npos);
            }
        }
        CHECK(recording_support::refusal([&] { (void)dilute.state(std::numeric_limits<double>::denorm_min(), 100, {1}, derived); })
            .find("nonfinite pure thermal derivatives") != std::string::npos);
        const epcsaft::HelmholtzModel columns(properties_support::ethane(),
            {{epcsaft::ParameterFamily::SegmentCount, {0}}, {epcsaft::ParameterFamily::DispersionEnergyOverK, {0}}});
        bool unavailable_closure = false;
        for (double rho : {1000., 3000., 5000., 8000., 10000.}) {
            const auto state = columns.state(200, rho, {1}, derived);
            CHECK(state.permittivity_parameter_derivative.empty());
            CHECK(state.density_parameter_derivative.size() == 2 && state.log_fugacity_parameter_derivative.size() == 2);
            if (state.pressure_density_derivative > 0 && state.log_fugacity_coefficient) continue;
            unavailable_closure = true;
            for (std::size_t a = 0; a < 2; ++a) {
                CHECK(!state.density_parameter_derivative[a] && !state.log_fugacity_parameter_derivative[a]);
                const auto name = "density_parameter_derivative[" + std::to_string(a) + "]";
                CHECK(std::any_of(state.unavailable_properties.begin(), state.unavailable_properties.end(),
                    [&](const auto& entry) { return entry.first == name && entry.second != epcsaft::PropertyActionStatus::Available; }));
            }
        }
        CHECK(unavailable_closure);
        const auto action = dilute.density_actions(300, 1e-100, {1},
            {epcsaft::PropertyObservable::ResidualEnthalpy}, {{0, 1e-100, {0}, {}}}).results.front();
        CHECK((action.action && std::isfinite(*action.action)) || action.status != epcsaft::PropertyActionStatus::Available);
        std::cout << "action,1e-100," << static_cast<int>(action.status) << ',' << action.action.value_or(NAN) << '\n';
    }
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> cases = {
        {properties_support::neutral_binary(), {0.4, 0.6}},
        {properties_support::associating_binary(), {0.1385, 0.8615}},
        {properties_support::aqueous_electrolyte(), {0.98, 0.01, 0.01}},
        {properties_support::nine_species(), {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0, 0.005, 0}},
    };
    for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
        const auto& item = cases[case_index];
        epcsaft::HelmholtzModel model(item.first);
        const double temperature = 300.0;
        const double volume = 1.0 / 1.2e4;
        const auto value = model.evaluate(temperature, volume, item.second,
                                          epcsaft::DerivativeOrder::Value);
        const auto state = model.state(temperature, 1.2e4, item.second);
        const auto inputs = model_inputs(item.first, temperature, volume, item.second);
        recording_support::compare(std::string("value-state/") + std::to_string(case_index), "residual_helmholtz", "1",
            value.value, state.residual_helmholtz, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14,
            "state-evaluation", inputs);
        const double residual_sum = state.hard_chain + state.dispersion
            + state.association + state.debye_huckel + state.born;
        recording_support::compare(std::string("residual-sum/") + std::to_string(case_index), "residual_helmholtz",
            "1", state.residual_helmholtz, residual_sum, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12,
            "contribution-sum", inputs);
        const auto gradient = model.evaluate(temperature, volume, item.second,
                                             epcsaft::DerivativeOrder::Gradient);
        recording_support::observe(std::string("gradient-density/") + std::to_string(case_index), "a_rho", "1",
            gradient.gradient[1], std::nullopt, "finite-gradient", inputs);
        CHECK(std::isfinite(gradient.gradient[1]));
        const auto hessian = model.evaluate(temperature, volume, item.second,
                                            epcsaft::DerivativeOrder::Hessian);
        const auto phase = model.evaluate(temperature, volume, item.second, epcsaft::DerivativeOrder::IsothermalPhaseHessian);
        const auto selected = model.evaluate(temperature, volume, item.second, epcsaft::DerivativeOrder::IsothermalPhaseHessian, true);
        CHECK(selected.value == phase.value && selected.gradient == phase.gradient);
        CHECK(selected.gradient_floor == phase.gradient_floor && selected.site_fractions == phase.site_fractions);
        CHECK(selected.first_numerical_failure.has_value() == phase.first_numerical_failure.has_value());
        const auto width = selected.gradient.size();
        if (item.second.size() == 2) for (double scale : {0.3, 3.0}) {
            auto amounts = item.second; for (auto& n : amounts) n *= scale;
            const auto scaled = model.evaluate(temperature, volume * scale, amounts, epcsaft::DerivativeOrder::IsothermalPhaseHessian, true);
            for (std::size_t i = 1; i < width; ++i) for (std::size_t j = 1; j < width; ++j)
                CHECK(close(scale * scaled.hessian[i * width + j], selected.hessian[i * width + j], 1e-10));
        }

        for (std::size_t i = 1; i < 2 + item.second.size(); ++i) {
            double contraction = 0.0, magnitude = 0.0;
            for (std::size_t j = 1; j < 2 + item.second.size(); ++j) {
                CHECK(close(selected.hessian[i * width + j], phase.hessian[i * width + j], 1e-10));
                const double product = selected.hessian[i * width + j] * (j == 1 ? volume : item.second[j - 2]);
                contraction += product; magnitude += std::abs(product);
            }
            CHECK(std::abs(contraction) <= 1e-10 * std::max(1.0, magnitude));
        }
        for (std::size_t row = 0; row < gradient.gradient.size(); ++row) {
            for (std::size_t column = row + 1; column < gradient.gradient.size(); ++column) {
                const double left = hessian.hessian[row * gradient.gradient.size() + column];
                const double right = hessian.hessian[column * gradient.gradient.size() + row];
                const auto point = std::to_string(case_index) + "/"
                    + std::to_string(row) + "-" + std::to_string(column);
                recording_support::compare(std::string("hessian-symmetry/") + point, "residual_hessian", "1", left,
                    right, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "symmetry", inputs);
            }
        }
        double euler = volume * gradient.gradient[1];
        for (std::size_t component = 0; component < item.second.size(); ++component) {
            euler += item.second[component] * gradient.gradient[2 + component];
        }
        recording_support::compare(std::string("euler/") + std::to_string(case_index), "residual_helmholtz", "1",
            euler, value.value, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "Euler-homogeneity", inputs);
        const auto state_identity = model.state(temperature, 1.2e4, item.second);
        double chemical_identity = 0.0;
        for (std::size_t component = 0; component < item.second.size(); ++component) {
            chemical_identity += item.second[component]
                * state_identity.residual_chemical_potential_over_rt[component].value();
        }
        recording_support::compare(std::string("chemical-potential/") + std::to_string(case_index),
            "chemical_potential_identity", "1", chemical_identity,
            state_identity.residual_helmholtz + state_identity.compressibility_factor - 1.0,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "chemical-potential-identity", inputs);
    }
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const double close_density = model.molar_density_upper_bound(300.0, {0.4, 0.6})
        * epcsaft::policy::DENSITY_UPPER_FRACTION.value;
    const auto near_cap = model.evaluate(300.0, 1.0 / close_density, {0.4, 0.6},
                                         epcsaft::DerivativeOrder::Value);
    CHECK(std::isfinite(near_cap.value));
}

void association_volume_tangent() {
    using namespace epcsaft;
    auto p = properties_support::associating_binary();
    p.site_multiplicity.assign(4, 1);
    constexpr double t = 300.0, rho = 10000.0;
    const std::vector<double> x{0.4, 0.6};
    const HelmholtzModel symmetric(p, {{ParameterFamily::AssociationVolume, {1, 4}}});
    const HelmholtzModel asymmetric(p, {{ParameterFamily::AssociationVolume, {3, 12}}});
    const auto sigma = segment_sigma(t, p), diameter = segment_diameter(t, p, sigma);
    const auto strength = association_delta(t, rho, x, p, sigma, diameter);
    std::optional<AssociationAction> accepted;
    const auto sites = solve_site_fractions(t, rho, x, p, sigma, diameter, nullptr, nullptr, &strength,
        false, {}, nullptr, nullptr, nullptr, &accepted, &symmetric.association_plan());
    CHECK(accepted && accepted->strategy == AssociationPlan::BalancedBinary);
    const auto bonded = bonded_densities(x, p, sites, strength);
    const Eigen::FullPivLU<Eigen::MatrixXd> full(mass_action_jacobian(
        x, p, sites, rho * NUMBER_DENSITY_PER_ANGSTROM3, strength, bonded));
    for (const auto* model : {&symmetric, &asymmetric}) {
        const bool broken = model == &asymmetric;
        const std::vector<double> point{t, 1.0/rho, x[0], x[1], model->active_parameter_values()[0]};
        auto tape = detail::HelmholtzTape::record(*model, point, sites, {});
        const std::vector<std::vector<double>> path{point, {0, 0, 0, 0, 1}};
        std::vector<double> actual, reference, right;
        std::size_t solves = 0;
        tape->implicit_taylor(path, 1, &full, nullptr, nullptr, nullptr, nullptr, &reference);
        tape->implicit_taylor(path, 1, &full, &solves, nullptr, nullptr, nullptr, &actual, &right, &*accepted);
        CHECK((right[4] != right[5] || right[6] != right[7]) == broken);
        CHECK(solves == (broken ? 1 : 0));
        for (std::size_t i = 4; i < actual.size(); ++i)
            recording_support::compare(broken ? "cross-association-volume" : "self-association-volume",
                "site_volume_tangent", "1", actual[i], reference[i],
                recording_support::Rule::ReferenceOrOne, 0.0, 1e-10, "full implicit mass action");

        // The two signed insertions have the same center and tangent. Check every returned
        // coefficient against the negative insertion evaluated from order zero.
        auto inserted = path;
        inserted.resize(4, std::vector<double>(point.size(), 0.0));
        inserted[2][0] = -1.0;
        std::vector<std::vector<double>> expected_output, output;
        std::vector<double> expected_sites, expected_right;
        tape->implicit_taylor(inserted, 3, &full, nullptr, nullptr, nullptr, &expected_output, &expected_sites, &expected_right);
        inserted[2][0] = 1.0;
        tape->implicit_taylor(inserted, 3, &full, nullptr, nullptr, nullptr, &output, &actual, &right);
        inserted[2][0] = -1.0;
        tape->implicit_taylor(inserted, 3, &full, nullptr, nullptr, nullptr, &output, &actual, &right, nullptr, 2);
        CHECK(output == expected_output && actual == expected_sites && right == expected_right);
    }
}

void finite_difference_gradient() {
    epcsaft::HelmholtzModel model(
        properties_support::associating_binary(),
        {{epcsaft::ParameterFamily::Kij, {1, 2}}}
    );
    const double temperature = 376.88;
    const double volume = 1.0 / 1e4;
    const std::vector<double> amounts = {0.1385, 0.8615};
    const auto actual = model.evaluate(temperature, volume, amounts,
                                       epcsaft::DerivativeOrder::Gradient);
    const auto inputs = model_inputs(model.parameters(), temperature, volume, amounts);
    for (std::size_t coordinate = 0; coordinate < actual.gradient.size(); ++coordinate) {
        if (coordinate == amounts.size() + 2) {
            const double step = 1e-5;
            auto plus = properties_support::associating_binary();
            auto minus = plus;
            plus.k_ij[1] += step; plus.k_ij[2] += step;
            minus.k_ij[1] -= step; minus.k_ij[2] -= step;
            const double up = epcsaft::HelmholtzModel(plus).evaluate(
                temperature, volume, amounts, epcsaft::DerivativeOrder::Value
            ).value;
            const double down = epcsaft::HelmholtzModel(minus).evaluate(
                temperature, volume, amounts, epcsaft::DerivativeOrder::Value
            ).value;
            recording_support::compare("gradient/parameter", "gradient", "1", actual.gradient[coordinate],
                (up - down) / (2.0 * step), recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-6, "finite-difference",
                recording_support::fields({inputs, recording_support::field("coordinate", coordinate),
                                           recording_support::field("step", step),
                                           recording_support::field("parameter", "kij[1,2]"),
                                           recording_support::field("kij", model.parameters().k_ij),
                                           recording_support::field("kij_plus", plus.k_ij),
                                           recording_support::field("kij_minus", minus.k_ij)}));
            continue;
        }
        std::vector<double> plus = amounts;
        std::vector<double> minus = amounts;
        double t_plus = temperature;
        double t_minus = temperature;
        double v_plus = volume;
        double v_minus = volume;
        const double step = coordinate == 0 ? temperature * 1e-5
            : coordinate == 1 ? volume * 1e-5 : 1e-5;
        if (coordinate == 0) { t_plus += step; t_minus -= step; }
        else if (coordinate == 1) { v_plus += step; v_minus -= step; }
        else { plus[coordinate - 2] += step; minus[coordinate - 2] -= step; }
        const double up = model.evaluate(t_plus, v_plus, plus,
                                         epcsaft::DerivativeOrder::Value).value;
        const double down = model.evaluate(t_minus, v_minus, minus,
                                           epcsaft::DerivativeOrder::Value).value;
        const double finite_difference = (up - down) / (2.0 * step);
        recording_support::compare(std::string("gradient/") + std::to_string(coordinate), "gradient", "1",
            actual.gradient[coordinate], finite_difference, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-6,
            "finite-difference", recording_support::fields({inputs, recording_support::field("coordinate", coordinate),
                                       recording_support::field("step", step),
                                       recording_support::field("plus_T", t_plus),
                                       recording_support::field("minus_T", t_minus),
                                       recording_support::field("plus_V", v_plus),
                                       recording_support::field("minus_V", v_minus),
                                       recording_support::field("plus_n", plus),
                                       recording_support::field("minus_n", minus)}));
    }
    association_volume_tangent();
}

}  // namespace properties_checks
