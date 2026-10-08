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
        {"ethane-vapor", properties_support::ethane(), {1.0}, 200.0, 217800.69494692178,
         epcsaft::DensityBranch::Vapor, "vapor"},
        {"ethane-liquid", properties_support::ethane(), {1.0}, 200.0, 217800.69494692178,
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
                const char* group, const char* point) {
    const bool ok = close(model, reference, tolerance);
    const bool cv = std::string(group) == "p1-cv";
    const bool cp = std::string(group) == "p1-cp";
    recording_support::compare(
        id.c_str(), quantity, unit, model, reference, inputs, "Relationship",
        cv ? "centered_temperature_derivative"
           : cp ? "independent_fixed_pressure_derivative" : "dilute_density_refinement",
        cv ? "c_v^r=(du^r/dT)_(rho,x)"
           : cp ? "c_p^r=(dh^r/dT)_(P,x)" : "field/rho has a converged dilute limit",
        "reference_or_one", 0.0, tolerance, ok, group, point
    );
}

void p1_observe(const std::string& id, const char* quantity, const char* unit, double value,
                const std::string& inputs, const char* basis, const char* method,
                const char* claim, const char* group, const char* point) {
    recording_support::observe(
        id.c_str(), quantity, unit, value, std::nullopt, inputs, basis, method, claim,
        group, point
    );
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
        recording_support::compare(
            comparison_id, "residual_internal_energy", "J/mol", state.residual_internal_energy, centered,
            point_inputs, "Relationship", "centered_temperature_derivative",
            "u_res=-R*T^2*a_T_res", "reference_or_one", 0.0, 2.0e-7,
            close(state.residual_internal_energy, centered, 2.0e-7),
            "f1-ladder", point, "u^{res}", "T,\\rho,\\mathbf{x},\\Delta T",
            "u^{res}=-RT^2a_T^{res}", "T,rho,x,delta_T"
        );
    }
    // a_T from the 0.003 K step: s = -R(a + T a_T) and h = RT(Z - 1) - RT^2 a_T.
    const double entropy = -gas_constant * (state.residual_helmholtz + temperature * a_t);
    const double enthalpy = gas_constant * temperature * (state.compressibility_factor - 1.0 - temperature * a_t);
    // The 0.003 K difference agrees to about 2e-10 (truncation-limited); 1e-8 keeps a 50x margin.
    recording_support::compare(
        "residual_entropy", "residual_entropy", "J/mol/K", state.residual_entropy, entropy, inputs,
        "Relationship", "centered_temperature_derivative", "s_res=-R*(a_res+T*a_T_res)", "reference_or_one",
        0.0, 1.0e-8, close(state.residual_entropy, entropy, 1.0e-8),
        "f1-entropy", "s", "s^{res}", "T,\\rho,\\mathbf{x}", "s^{res}=-R(a^{res}+Ta_T^{res})");
    recording_support::compare(
        "residual_enthalpy", "residual_enthalpy", "J/mol", state.residual_enthalpy, enthalpy, inputs,
        "Relationship", "centered_temperature_derivative", "h_res=RT*(Z-1)-R*T^2*a_T_res", "reference_or_one",
        0.0, 1.0e-8, close(state.residual_enthalpy, enthalpy, 1.0e-8),
        "f1-enthalpy", "h", "h^{res}", "T,\\rho,\\mathbf{x}", "h^{res}=RT(Z-1)-RT^2a_T^{res}");
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
            p1_compare("cv-ladder/" + point, "residual_isochoric_heat_capacity", "J/mol/K", density_state.residual_isochoric_heat_capacity, cv_reference, recording_support::fields({inputs, recording_support::field("route", "fixed-density"), recording_support::field("delta_T", delta), recording_support::field("P_plus", cv_plus.pressure), recording_support::field("P_minus", cv_minus.pressure)}), tolerance, "p1-cv", point.c_str());
            const auto cp_plus = model.state_at_pressure(item.temperature + delta, item.pressure, item.composition, item.branch);
            const auto cp_minus = model.state_at_pressure(item.temperature - delta, item.pressure, item.composition, item.branch);
            CHECK(close(cp_plus.pressure, item.pressure, 1.0e-9));
            CHECK(close(cp_minus.pressure, item.pressure, 1.0e-9));
            const double cp_reference = (cp_plus.residual_enthalpy - cp_minus.residual_enthalpy) / (2.0 * delta);
            p1_compare("cp-ladder/" + point, "residual_isobaric_heat_capacity", "J/mol/K", cp, cp_reference, recording_support::fields({inputs, p1_pressure_record("plus", cp_plus, item.pressure, item.branch_name), p1_pressure_record("minus", cp_minus, item.pressure, item.branch_name), recording_support::field("route", "fixed-pressure"), recording_support::field("delta_T", delta)}), tolerance, "p1-cp", point.c_str());
            if (step == 2) {
                cv_coarse = cv_reference;
                cp_coarse = cp_reference;
            } else if (step == 3) {
                CHECK(close(cv_coarse, cv_reference, 2.0e-8));
                CHECK(close(cp_coarse, cp_reference, 2.0e-8));
            }
        }
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
        p1_observe("dilute-state/" + recording_support::number(dilute_densities[index]), "residual_helmholtz_over_rho", "m3/mol", state.residual_helmholtz / dilute_densities[index], recording_support::fields({model_inputs(dilute_model.parameters(), 300.0, 1.0 / dilute_densities[index], dilute_composition), recording_support::field("D", D), recording_support::field("fields_over_rho", std::vector<double>{state.residual_entropy / dilute_densities[index], state.residual_internal_energy / dilute_densities[index], state.residual_enthalpy / dilute_densities[index], state.residual_isochoric_heat_capacity / dilute_densities[index], *state.residual_isobaric_heat_capacity / dilute_densities[index]})}), "Observation", "dilute_ladder", "residual state over density remains finite", "p1-dilute", recording_support::number(dilute_densities[index]).c_str());
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
                   dilute_inputs, 5.0e-6, "p1-dilute-limit", quantity);
    }
    p1_compare("dilute-limit/residual_isobaric_heat_capacity", "residual_isobaric_heat_capacity", "J m3/(mol2 K)", *coarse.residual_isobaric_heat_capacity / dilute_densities[2], *fine.residual_isobaric_heat_capacity / dilute_densities[3], dilute_inputs, 5.0e-6, "p1-dilute-limit", "residual_isobaric_heat_capacity");

    epcsaft::HelmholtzModel conditioning_model(properties_support::ethane());
    const auto conditioned = conditioning_model.state(280.0, 3000.0, {1.0});
    const auto unavailable = conditioning_model.state(280.0, 4000.0, {1.0});
    const double conditioned_D = conditioned.pressure_density_derivative / (gas_constant * 280.0);
    const double unavailable_D = unavailable.pressure_density_derivative / (gas_constant * 280.0);
    CHECK(std::isfinite(conditioned_D) && conditioned_D > 0.0);
    CHECK(conditioned.residual_isobaric_heat_capacity.has_value());
    p1_observe("conditioning/near-zero-positive-D", "D", "1", conditioned_D, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 3000.0, {1.0}), recording_support::field("D", conditioned_D), recording_support::field("cp", *conditioned.residual_isobaric_heat_capacity)}), "Observation", "conditioning_stress", "near-zero positive D is conditioning stress", "p1-conditioning", "near-zero-positive-D");
    p1_observe("conditioning/near-zero-positive-cp", "residual_isobaric_heat_capacity", "J/mol/K", *conditioned.residual_isobaric_heat_capacity, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 3000.0, {1.0}), recording_support::field("D", conditioned_D)}), "Observation", "conditioning_stress", "finite c_p^r retained for positive D", "p1-conditioning", "near-zero-positive-cp");
    CHECK(std::isfinite(unavailable_D) && unavailable_D < 0.0);
    const bool unavailable_cp = !unavailable.residual_isobaric_heat_capacity.has_value();
    recording_support::compare_custom("conditioning/negative-D-cp-unavailable", "residual_isobaric_heat_capacity_available", "bool", unavailable_cp ? 0.0 : 1.0, 0.0, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 4000.0, {1.0}), recording_support::field("D", unavailable_D), recording_support::field("unavailable_reason", "nonpositive mechanical denominator")}), "Availability", "mechanical_conditioning", "c_p^r is unavailable when D<0", "absolute", unavailable_cp, "p1-conditioning", "negative-D-cp-unavailable");
    const bool finite_other_outputs = finite_state_without_cp(unavailable);
    recording_support::compare_custom("conditioning/negative-D-other-fields-finite", "residual_properties_finite", "bool", finite_other_outputs ? 1.0 : 0.0, 1.0, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 4000.0, {1.0}), recording_support::field("D", unavailable_D)}), "Availability", "partial_availability", "other residual outputs remain finite when c_p^r is unavailable", "absolute", finite_other_outputs, "p1-conditioning", "negative-D-other-fields-finite");
}
void derivatives_and_limits() {
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
        recording_support::compare(
            (std::string("value-state/") + std::to_string(case_index)).c_str(),
            "residual_helmholtz", "1", value.value,
            state.residual_helmholtz, inputs, "Relationship", "state-evaluation",
            "evaluate(value)=state(residual_helmholtz)", "reference_or_one", 0.0, 1.0e-14,
            close(value.value, state.residual_helmholtz, 1.0e-14),
            "m2-value-state", std::to_string(case_index).c_str(), "a^{res}",
            "T,V,\\mathbf{n}", "a_{eval}=a_{state}", "state");
        CHECK(close(value.value, state.residual_helmholtz, 1e-14));
        const double residual_sum = state.hard_chain + state.dispersion
            + state.association + state.debye_huckel + state.born;
        recording_support::compare(
            (std::string("residual-sum/") + std::to_string(case_index)).c_str(),
            "residual_helmholtz", "1", state.residual_helmholtz, residual_sum, inputs,
            "Relationship", "contribution-sum", "a_res=sum(contributions)",
            "reference_or_one", 0.0, 1.0e-12,
            close(state.residual_helmholtz, residual_sum, 1.0e-12),
            "m2-residual-sum", std::to_string(case_index).c_str(), "a^{res}",
            "T,V,\\mathbf{n}", "a^{res}=a^{hc}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}", "case");
        CHECK(close(state.residual_helmholtz, residual_sum, 1e-12));
        const auto gradient = model.evaluate(temperature, volume, item.second,
                                             epcsaft::DerivativeOrder::Gradient);
        recording_support::observe(
            (std::string("gradient-density/") + std::to_string(case_index)).c_str(),
            "a_rho", "1", gradient.gradient[1], std::nullopt, inputs, "Relationship",
            "finite-gradient", "density derivative is finite", "m2-gradient",
            std::to_string(case_index).c_str(), "a_\\rho", "T,V,\\mathbf{n}",
            "a_\\rho", "case");
        CHECK(std::isfinite(gradient.gradient[1]));
        const auto hessian = model.evaluate(temperature, volume, item.second,
                                            epcsaft::DerivativeOrder::Hessian);
        for (std::size_t row = 0; row < gradient.gradient.size(); ++row) {
            for (std::size_t column = row + 1; column < gradient.gradient.size(); ++column) {
                const double left = hessian.hessian[row * gradient.gradient.size() + column];
                const double right = hessian.hessian[column * gradient.gradient.size() + row];
                const auto point = std::to_string(case_index) + "/"
                    + std::to_string(row) + "-" + std::to_string(column);
                recording_support::compare(
                    (std::string("hessian-symmetry/") + point).c_str(),
                    "residual_hessian", "1", left, right, inputs, "Relationship",
                    "symmetry", "H_ij=H_ji", "reference_or_one", 0.0, 1.0e-12,
                    close(left, right, 1.0e-12), "m2-hessian-symmetry", point.c_str(),
                    "H^{res}", "T,V,\\mathbf{n}", "H_{ij}=H_{ji}", "row,col");
            }
        }
        double euler = volume * gradient.gradient[1];
        for (std::size_t component = 0; component < item.second.size(); ++component) {
            euler += item.second[component] * gradient.gradient[2 + component];
        }
        recording_support::compare(
            (std::string("euler/") + std::to_string(case_index)).c_str(),
            "residual_helmholtz", "1", euler, value.value, inputs,
            "Relationship", "Euler-homogeneity", "V*a_V+sum(n_i*a_n_i)=a",
            "reference_or_one", 0.0, 1.0e-12, close(euler, value.value, 1e-12),
            "m2-euler", std::to_string(case_index).c_str(), "a^{res}", "T,V,\\mathbf{n}",
            "Va_V+\\sum_i n_i a_{n_i}=a", "state");
        const auto state_identity = model.state(temperature, 1.2e4, item.second);
        double chemical_identity = 0.0;
        for (std::size_t component = 0; component < item.second.size(); ++component) {
            chemical_identity += item.second[component]
                * state_identity.residual_chemical_potential_over_rt[component].value();
        }
        recording_support::compare(
            (std::string("chemical-potential/") + std::to_string(case_index)).c_str(),
            "chemical_potential_identity", "1", chemical_identity,
            state_identity.residual_helmholtz + state_identity.compressibility_factor - 1.0,
            inputs, "Relationship", "chemical-potential-identity",
            "sum(x_i*mu_i/RT)=a+Z-1", "reference_or_one", 0.0, 1.0e-12,
            close(chemical_identity,
                  state_identity.residual_helmholtz + state_identity.compressibility_factor - 1.0,
            1.0e-12), "m2-chemical-potential", std::to_string(case_index).c_str(),
            "\\sum_i x_i\\mu_i", "T,V,\\mathbf{n}",
            "\\sum_i x_i\\mu_i/RT=a+Z-1", "component");
        CHECK(close(chemical_identity,
                    state_identity.residual_helmholtz
                    + state_identity.compressibility_factor - 1.0, 1e-12));
    }
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const double close_density = model.molar_density_upper_bound(300.0, {0.4, 0.6})
        * epcsaft::policy::DENSITY_UPPER_FRACTION.value;
    model.evaluate(300.0, 1.0 / close_density, {0.4, 0.6},
                   epcsaft::DerivativeOrder::Value);
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
            recording_support::compare(
                "gradient/parameter", "gradient", "1", actual.gradient[coordinate],
                (up - down) / (2.0 * step),
                recording_support::fields({inputs, recording_support::field("coordinate", coordinate),
                                           recording_support::field("step", step),
                                           recording_support::field("parameter", "kij[1,2]"),
                                           recording_support::field("kij", model.parameters().k_ij),
                                           recording_support::field("kij_plus", plus.k_ij),
                                           recording_support::field("kij_minus", minus.k_ij)}),
                "Relationship", "finite-difference", "gradient=central_difference",
                "reference_or_one", 0.0, 1.0e-6,
                close(actual.gradient[coordinate], (up - down) / (2.0 * step), 1e-6),
                "m3-gradient", "parameter", "\\partial a/\\partial k_{ij}",
                "T,V,\\mathbf{n},\\Delta k", "a_k=\\Delta a/\\Delta k", "parameter");
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
        recording_support::compare(
            (std::string("gradient/") + std::to_string(coordinate)).c_str(),
            "gradient", "1", actual.gradient[coordinate], finite_difference,
            recording_support::fields({inputs, recording_support::field("coordinate", coordinate),
                                       recording_support::field("step", step),
                                       recording_support::field("plus_T", t_plus),
                                       recording_support::field("minus_T", t_minus),
                                       recording_support::field("plus_V", v_plus),
                                       recording_support::field("minus_V", v_minus),
                                       recording_support::field("plus_n", plus),
                                       recording_support::field("minus_n", minus)}),
            "Relationship", "finite-difference", "gradient=central_difference",
            "reference_or_one", 0.0, 1.0e-6,
            close(actual.gradient[coordinate], finite_difference, 1e-6),
            "m3-gradient", std::to_string(coordinate).c_str(), "\\nabla a", "T,V,\\mathbf{n},\\Delta z",
            "a_i=\\Delta a/\\Delta z_i", "coordinate");
    }
}

}  // namespace properties_checks
