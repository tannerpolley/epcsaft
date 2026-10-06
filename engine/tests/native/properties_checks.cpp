#include "properties/helmholtz_model.hpp"
#include "properties_support.hpp"
#include "equilibrium_support.hpp"
#undef CHECK
#include "recording_support.hpp"
#include "contributions/association.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#define CHECK(condition) recording_support::check( \
    static_cast<bool>(condition), __FILE__, __LINE__, #condition)

namespace {

bool close(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

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

bool finite_state_without_cp(const epcsaft::State& state) {
    return std::isfinite(state.residual_helmholtz) && std::isfinite(state.pressure)
        && std::isfinite(state.residual_entropy) && std::isfinite(state.residual_internal_energy)
        && std::isfinite(state.residual_enthalpy)
        && std::isfinite(state.residual_isochoric_heat_capacity)
        && std::isfinite(state.molar_density)
        && std::isfinite(state.pressure_density_derivative);
}

void p1_group(const char* id, const char* prefix, std::size_t count,
              const char* definition, const char* coordinates) {
    recording_support::expect_group(
        recording_support::active_test_id(), id, prefix, count, definition, coordinates
    );
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
    CHECK(ok);
}

void p1_observe(const std::string& id, const char* quantity, const char* unit, double value,
                const std::string& inputs, const char* basis, const char* method,
                const char* claim, const char* group, const char* point) {
    recording_support::observe(
        id.c_str(), quantity, unit, value, std::nullopt, inputs, basis, method, claim,
        group, point
    );
}

std::string model_inputs(const epcsaft::Parameters& parameters, double temperature,
                         double volume, const std::vector<double>& amounts) {
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    const std::vector<double> composition = total == 0.0
        ? std::vector<double>{} : [&] {
            std::vector<double> result;
            result.reserve(amounts.size());
            for (const double amount : amounts) result.push_back(amount / total);
            return result;
        }();
    return recording_support::fields({
        recording_support::field("T", temperature),
        recording_support::field("V", volume),
        recording_support::field("rho", total / volume),
        recording_support::field("n", amounts),
        recording_support::field("x", composition),
        recording_support::field("segments", parameters.segment_count),
        recording_support::field("sigma_A", parameters.sigma_constant),
        recording_support::field("epsilon_over_k_K", parameters.dispersion_energy_over_k),
    });
}

void residual_internal_energy_identity_ladder() {
    // F1 previous-failure check: the temperature derivative is independent of
    // the state caloric fields, and the centered ladder checks its local error.
    epcsaft::HelmholtzModel model(properties_support::ethane());
    constexpr double temperature = 300.0;
    constexpr double density = 1.0e4;
    constexpr double gas_constant = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const std::vector<double> composition = {1.0};
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-a-res", "a_res", 1,
        "F1 actual state identities", "quantity");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-pressure", "pressure", 1,
        "F1 actual state identities", "quantity");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-compressibility", "compressibility_factor", 1,
        "F1 actual state identities", "quantity");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-entropy", "residual_entropy", 1,
        "F1 actual state identities", "quantity");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-enthalpy", "residual_enthalpy", 1,
        "F1 actual state identities", "quantity");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-ladder", "temperature-ladder/", 4,
        "F1 centered temperature derivative ladder", "delta_T");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-enthalpy-identity", "enthalpy-internal-energy-identity", 1,
        "F1 enthalpy identity", "quantity");
    recording_support::expect_group(
        recording_support::active_test_id(), "f1-internal-energy-identity", "internal-energy-identity", 1,
        "F1 internal energy identity", "quantity");
    const auto state = model.state(temperature, density, composition);
    const auto inputs = model_inputs(
        model.parameters(), temperature, 1.0 / density, composition
    );
    recording_support::observe(
        "a_res", "residual_helmholtz", "1", state.residual_helmholtz,
        state.hard_chain + state.dispersion + state.association
            + state.debye_huckel + state.born,
        inputs, "Relationship", "contribution_sum",
        "residual Helmholtz contribution identity", "f1-a-res", "a_res",
        "a^{res}", "T,\\rho,\\mathbf{x}", "a^{res}=a^{hc}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}"
    );
    recording_support::observe(
        "pressure", "pressure", "Pa", state.pressure,
        density * gas_constant * temperature * state.compressibility_factor,
        inputs, "Relationship", "pressure_identity",
        "P=rho*R*T*Z", "f1-pressure", "pressure", "P", "T,\\rho,\\mathbf{x}", "P=\\rho RTZ"
    );
    const auto derivative = model.evaluate(
        temperature, 1.0 / density, composition, epcsaft::DerivativeOrder::Gradient
    );
    recording_support::observe(
        "compressibility_factor", "compressibility_factor", "1",
        state.compressibility_factor, 1.0 - derivative.gradient[1] / density,
        inputs, "Relationship", "density_derivative_identity",
        "Z=1+rho*a_rho", "f1-compressibility", "Z", "Z", "T,\\rho,\\mathbf{x}", "Z=1+\\rho a_\\rho"
    );
    recording_support::observe(
        "residual_entropy", "residual_entropy", "J/mol/K", state.residual_entropy,
        -gas_constant * (derivative.value + temperature * derivative.gradient[0]),
        inputs, "Relationship", "temperature_derivative_identity",
        "s=-R*(a+T*a_T)", "f1-entropy", "s", "s^{res}", "T,\\rho,\\mathbf{x}", "s^{res}=-R(a^{res}+Ta_T^{res})"
    );
    recording_support::observe(
        "residual_enthalpy", "residual_enthalpy", "J/mol", state.residual_enthalpy,
        gas_constant * temperature * (-temperature * derivative.gradient[0]
        + state.compressibility_factor - 1.0),
        inputs, "Relationship", "enthalpy_identity",
        "h=u+RT*(Z-1)", "f1-enthalpy", "h", "h^{res}", "T,\\rho,\\mathbf{x}", "h^{res}=u^{res}+RT(Z-1)"
    );
    const double identity = -gas_constant * temperature * temperature
        * derivative.gradient.front();
    CHECK(std::isfinite(identity));
    for (const double delta_temperature : {0.1, 0.03, 0.01, 0.003}) {
        const double plus = model.state(
            temperature + delta_temperature, density, composition
        ).residual_helmholtz;
        const double minus = model.state(
            temperature - delta_temperature, density, composition
        ).residual_helmholtz;
        const double centered = -gas_constant * temperature * temperature
            * (plus - minus) / (2.0 * delta_temperature);
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
            comparison_id, "residual_internal_energy", "J/mol", centered, identity,
            point_inputs, "Relationship", "centered_temperature_derivative",
            "u_res=-R*T^2*a_T_res", "reference_or_one", 0.0, 2.0e-7,
            std::abs(centered - identity) <= 2.0e-7 * std::max(1.0, std::abs(identity)),
            "f1-ladder", point, "u^{res}", "T,\\rho,\\mathbf{x},\\Delta T",
            "u^{res}=-RT^2a_T^{res}", "T,rho,x,delta_T"
        );
        CHECK(close(centered, identity, 2.0e-7));
    }
    const double enthalpy_identity = state.residual_enthalpy
        - gas_constant * temperature * (state.compressibility_factor - 1.0);
    recording_support::compare(
        "enthalpy-internal-energy-identity", "residual_internal_energy", "J/mol",
        state.residual_internal_energy, enthalpy_identity, inputs, "Relationship",
        "enthalpy_identity", "u_res=h_res-R*T*(Z-1)", "reference_or_one", 0.0, 1.0e-10,
        close(state.residual_internal_energy, enthalpy_identity, 1.0e-10),
        "f1-enthalpy-identity", "enthalpy", "u^{res}", "T,\\rho,\\mathbf{x}", "u^{res}=h^{res}-RT(Z-1)"
    );
    recording_support::compare(
        "internal-energy-identity", "residual_internal_energy", "J/mol",
        state.residual_internal_energy, identity, inputs, "Relationship",
        "temperature_derivative_identity", "u_res=-R*T^2*a_T_res", "reference_or_one",
        0.0, 1.0e-10,
        close(state.residual_internal_energy, identity, 1.0e-10),
        "f1-internal-energy-identity", "derivative", "u^{res}", "T,\\rho,\\mathbf{x}", "u^{res}=-RT^2a_T^{res}"
    );
    CHECK(close(state.residual_internal_energy, identity, 1.0e-10));
    CHECK(close(state.residual_internal_energy, enthalpy_identity, 1.0e-10));
}

void residual_properties_and_caloric_ladders() {
    constexpr double gas_constant = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    constexpr std::array<double, 4> deltas = {0.1, 0.03, 0.01, 0.003};
    const auto cases = p1_cases();
    p1_group("p1-cv", "cv-ladder/", cases.size() * deltas.size(), "fixed-density residual isochoric heat-capacity ladder", "T,rho,x,delta_T");
    p1_group("p1-cp", "cp-ladder/", cases.size() * deltas.size(), "fixed-pressure residual isobaric heat-capacity ladder", "T,P,x,delta_T");
    p1_group("p1-dilute", "dilute-state/", 4, "neutral residual caloric dilute ladder", "T,rho,x");
    p1_group("p1-dilute-limit", "dilute-limit/", 5, "neutral residual caloric field over density limit", "T,rho,x");
    p1_group("p1-conditioning", "conditioning/", 4, "residual isobaric heat-capacity conditioning", "T,rho,x,D");

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
    CHECK(unavailable_cp);
    const bool finite_other_outputs = finite_state_without_cp(unavailable);
    recording_support::compare_custom("conditioning/negative-D-other-fields-finite", "residual_properties_finite", "bool", finite_other_outputs ? 1.0 : 0.0, 1.0, recording_support::fields({model_inputs(conditioning_model.parameters(), 280.0, 1.0 / 4000.0, {1.0}), recording_support::field("D", unavailable_D)}), "Availability", "partial_availability", "other residual outputs remain finite when c_p^r is unavailable", "absolute", finite_other_outputs, "p1-conditioning", "negative-D-other-fields-finite");
    CHECK(finite_other_outputs);
}
void derivatives_and_limits() {
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> cases = {
        {properties_support::neutral_binary(), {0.4, 0.6}},
        {properties_support::associating_binary(), {0.1385, 0.8615}},
        {properties_support::aqueous_electrolyte(), {0.98, 0.01, 0.01}},
        {properties_support::nine_species(), {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0, 0.005, 0}},
    };
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-value-state", "value-state/", cases.size(),
        "state/value evaluations for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-residual-sum", "residual-sum/", cases.size(),
        "residual contribution sum for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-gradient", "gradient-density/", cases.size(),
        "finite density gradients for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-hessian-symmetry", "hessian-symmetry/", 178,
        "symmetric Hessian entries for all parameter cases", "case,row,column");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-euler", "euler/", cases.size(),
        "Euler homogeneity for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-internal-energy", "internal-energy/", cases.size(),
        "enthalpy/internal-energy identity for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-internal-energy-entropy", "internal-energy-entropy/", cases.size(),
        "residual caloric identity for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-chemical-potential", "chemical-potential/", cases.size(),
        "chemical potential identity for all parameter cases", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-low-density", "low-density/", 4,
        "low-density residual Helmholtz limit", "rho");
    recording_support::expect_group(
        recording_support::active_test_id(), "m2-low-density-z", "low-density-z/", 4,
        "low-density compressibility limit", "rho");
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
            for (std::size_t column = 0; column < gradient.gradient.size(); ++column) {
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
                CHECK(close(left, right, 1e-12));
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
        CHECK(close(euler, value.value, 1e-12));
        const auto state_identity = model.state(temperature, 1.2e4, item.second);
        const double rt = epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature;
        const double internal_from_h = state_identity.residual_enthalpy
            - rt * (state_identity.compressibility_factor - 1.0);
        recording_support::compare(
            (std::string("internal-energy/") + std::to_string(case_index)).c_str(),
            "residual_internal_energy", "J/mol",
            state_identity.residual_internal_energy, internal_from_h, inputs,
            "Relationship", "caloric-identity", "u=h-RT(Z-1)", "reference_or_one",
            0.0, 1.0e-10,
            close(state_identity.residual_internal_energy, internal_from_h, 1.0e-10),
            "m2-internal-energy", std::to_string(case_index).c_str(), "u^{res}", "T,V,\\mathbf{n}",
            "u^{res}=h^{res}-RT(Z-1)", "state");
        CHECK(close(state_identity.residual_internal_energy,
                    state_identity.residual_enthalpy
                    - rt * (state_identity.compressibility_factor - 1.0), 1e-10));
        const double internal_from_entropy = rt * state_identity.residual_helmholtz
            + temperature * state_identity.residual_entropy;
        recording_support::compare(
            (std::string("internal-energy-entropy/") + std::to_string(case_index)).c_str(),
            "residual_internal_energy", "J/mol", state_identity.residual_internal_energy,
            internal_from_entropy, inputs, "Relationship", "caloric-identity",
            "u=RTa+Ts", "reference_or_one", 0.0, 1.0e-10,
            close(state_identity.residual_internal_energy, internal_from_entropy, 1.0e-10),
            "m2-internal-energy-entropy", std::to_string(case_index).c_str(), "u^{res}",
            "T,V,\\mathbf{n}", "u^{res}=RTa^{res}+Ts^{res}", "case");
        CHECK(close(state_identity.residual_internal_energy, internal_from_entropy, 1e-10));
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
    for (double density : {1e1, 1e0, 1e-1, 1e-2}) {
        const auto state = model.state(300.0, density, {0.4, 0.6});
        recording_support::compare(
            (std::string("low-density/") + recording_support::number(density)).c_str(),
            "residual_helmholtz", "1", state.residual_helmholtz, 0.0,
            model_inputs(model.parameters(), 300.0, 1.0 / density, {0.4, 0.6}),
            "Relationship", "limit", "residual Helmholtz vanishes at low density",
            "strict_absolute", 1e-2, 0.0,
            std::abs(state.residual_helmholtz) < 1e-2,
            "m2-low-density", recording_support::number(density).c_str(), "a^{res}",
            "T,\\rho,\\mathbf{x}", "\\lim_{\\rho\\to0}a^{res}=0", "rho");
        CHECK(std::abs(state.residual_helmholtz) < 1e-2);
        recording_support::compare(
            (std::string("low-density-z/") + recording_support::number(density)).c_str(),
            "compressibility_factor", "1", state.compressibility_factor, 1.0,
            model_inputs(model.parameters(), 300.0, 1.0 / density, {0.4, 0.6}),
            "Relationship", "limit", "compressibility factor approaches one",
            "reference_or_one", 0.0, 1e-2,
            close(state.compressibility_factor, 1.0, 1e-2),
            "m2-low-density-z", recording_support::number(density).c_str(), "Z",
            "T,\\rho,\\mathbf{x}", "\\lim_{\\rho\\to0}Z=1", "rho");
        CHECK(close(state.compressibility_factor, 1.0, 1e-2));
    }
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
    recording_support::expect_group(
        recording_support::active_test_id(), "m3-gradient", "gradient/", amounts.size() + 3,
        "temperature, volume, amounts and active parameter gradients", "coordinate");
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
            CHECK(close(actual.gradient[coordinate], (up - down) / (2.0 * step), 1e-6));
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
        CHECK(close(actual.gradient[coordinate], finite_difference, 1e-6));
    }
}

void hessian_reference() {
    const auto associating_parameters = properties_support::associating_binary();
    const std::vector<double> associating_amounts = {0.1385, 0.8615};
    recording_support::expect_group(
        recording_support::active_test_id(), "m4-associating-gradient", "associating-gradient/", 4,
        "associating implicit-site gradient", "coordinate");
    recording_support::expect_group(
        recording_support::active_test_id(), "m4-associating-hessian", "associating-hessian/", 16,
        "associating implicit-site Hessian", "row,column");
    recording_support::expect_group(
        recording_support::active_test_id(), "m4-nine-gradient", "nine-gradient/", 12,
        "active-parameter gradient", "coordinate");
    recording_support::expect_group(
        recording_support::active_test_id(), "m4-nine-hessian", "nine-hessian/", 144,
        "active-parameter Hessian", "row,column");
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
        CHECK(close(associating_result.gradient[coordinate],
                    associating_reference.gradient[coordinate], 1e-10));
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
        CHECK(close(associating_result.hessian[entry],
                    associating_reference.hessian[entry], 1e-10));
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
        CHECK(close(nine_result.gradient[coordinate],
                    nine_reference.gradient[coordinate], 1e-10));
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
        CHECK(close(nine_result.hessian[entry], nine_reference.hessian[entry], 1e-10));
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
    recording_support::expect_group(
        recording_support::active_test_id(), "m5-projected", "projected/", 514,
        "projected and subsequent full derivatives against raw implicit elimination", "case,order,row,column");
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
                CHECK(close(value, expected, 1e-10));
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
    recording_support::expect_group(
        recording_support::active_test_id(), "m5-hessian", "hessian/", 178,
        "Hessian finite differences for all parameter cases", "case,row,column");
    recording_support::expect_group(
        recording_support::active_test_id(), "m5-pressure-slope", "pressure-density-slope/", cases.size(),
        "pressure density slope finite differences", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m5-volume-curvature", "volume-second-derivative/", cases.size(),
        "volume curvature finite differences", "case");
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
                CHECK(close(actual.hessian[row * width + coordinate],
                            finite_difference,
                            1e-6));
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
        CHECK(close(slope_state.pressure_density_derivative, pressure_density_fd, 1e-6));
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
        CHECK(close(pressure_slope.volume_second_derivative, volume_second_fd, 1e-6));
    }
}

void same_sign_pair_is_inert() {
    const std::vector<double> amounts = {
        0.795, 0.10, 0.05, 0.02, 0.02, 0.005, 0.001, 0.008, 0.001,
    };
    auto excluded = properties_support::nine_species();
    auto changed = excluded;
    changed.k_ij[3 * 9 + 7] = changed.k_ij[7 * 9 + 3] = 0.5;
    recording_support::expect_group(
        recording_support::active_test_id(), "m6-excluded", "excluded-inert", 1,
        "same-sign excluded interaction", "pair");
    recording_support::expect_group(
        recording_support::active_test_id(), "m6-all-pairs", "all-pairs-active", 1,
        "same-sign interaction under all-pairs policy", "pair");
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
    CHECK(close(base_state.dispersion, changed_state.dispersion, 1e-14));

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
    CHECK(all_pairs_active);
}

void parameter_derivative_and_helpers() {
    epcsaft::Parameters parameters = properties_support::neutral_binary();
    epcsaft::HelmholtzModel model(
        parameters, {{epcsaft::ParameterFamily::Kij, {1, 2}}}
    );
    recording_support::expect_group(
        recording_support::active_test_id(), "m7-parameter-gradient", "parameter-gradient", 1,
        "active k_ij derivative", "parameter");
    recording_support::expect_group(
        recording_support::active_test_id(), "m7-packing", "packing-coefficient", 1,
        "packing coefficient geometry", "rho");
    recording_support::expect_group(
        recording_support::active_test_id(), "m7-phase-packing", "phase-packing", 1,
        "phase packing fraction geometry", "rho");
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
    CHECK(close(derivative.gradient.back(), parameter_fd, 1e-6));
    const double coefficient = model.packing_fraction_coefficient(300.0, {0.4, 0.6});
    const double packing = epcsaft::packing_fraction_from_geometry(
        1.2e4, {0.4, 0.6}, parameters,
        epcsaft::segment_diameter(
            300.0, parameters, epcsaft::segment_sigma(300.0, parameters)));
    recording_support::compare(
        "packing-coefficient", "packing_fraction", "1", coefficient * 1.2e4, packing,
        inputs, "Relationship", "geometry", "rho*packing_coefficient=packing_fraction",
        "reference_or_one", 0.0, 1.0e-14, close(coefficient * 1.2e4, packing, 1e-14),
        "m7-packing", "rho", "\\eta", "T,\\rho,\\mathbf{x}",
        "\\eta=\\rho c_{\\eta}", "rho");
    CHECK(close(coefficient * 1.2e4,
                packing,
                1e-14));
    const double phase_packing = model.phase_packing_fraction(
        300.0, {0.4, 0.6}, 1.0 / 1.2e4
    );
    recording_support::compare(
        "phase-packing", "packing_fraction", "1", phase_packing, coefficient * 1.2e4,
        inputs, "Relationship", "geometry", "phase packing fraction uses coefficient",
        "reference_or_one", 0.0, 1.0e-14,
        close(phase_packing, coefficient * 1.2e4, 1e-14),
        "m7-phase-packing", "rho", "\\eta", "T,\\rho,\\mathbf{x}",
        "\\eta_{phase}=\\rho c_{\\eta}", "rho");
    CHECK(close(phase_packing, coefficient * 1.2e4, 1e-14));
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

void density_and_reference() {
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    const auto reference_rows = std::array<std::tuple<double, double, double, double>, 2>{
        std::tuple{200.0, 217800.69494692178, 0.007241100995,
                   0.03007 / 523.97698},
        std::tuple{240.0, 969152.4807171321, 0.001745731329,
                   0.03007 / 465.30887}};
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-vapor-pressure", "vapor-pressure/",
        reference_rows.size(), "vapor pressure density roots", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-liquid-pressure", "liquid-pressure/",
        reference_rows.size(), "liquid pressure density roots", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-vapor-order", "density-order/",
        reference_rows.size(), "vapor density is below liquid density", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-vapor-volume", "vapor-volume/",
        reference_rows.size(), "historical vapor specific-volume observations", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-liquid-volume", "liquid-volume/",
        reference_rows.size(), "historical liquid specific-volume observations", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-supercritical-pressure", "supercritical-pressure", 1,
        "supercritical pressure root", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-supercritical-density", "supercritical-density", 1,
        "supercritical branches agree", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-anchored-liquid", "anchored-liquid", 1,
        "anchored liquid branch", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-anchored-vapor", "anchored-vapor", 1,
        "anchored vapor branch", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-associating", "associating-pressure", 1,
        "associating pressure root", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-nine", "nine-packing", 1,
        "nine-species packing bound", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-reference-finite", "reference-finite", 1,
        "infinite-dilution reference is finite", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-activity-finite", "activity-finite", 1,
        "electrolyte activity is finite", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-reference-error", "reference-error", 1,
        "infinite-dilution reference limit", "ionic-strength");
    recording_support::expect_group(
        recording_support::active_test_id(), "m8-activity-scale", "activity-scale", 1,
        "electrolyte activity limit", "ionic-strength");
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
        CHECK(vapor.molar_density < liquid.molar_density);
        CHECK(close(vapor.pressure, pressure, 1e-9));
        CHECK(close(liquid.pressure, pressure, 1e-9));
        recording_support::compare(
            (std::string("vapor-volume/") + std::to_string(case_index)).c_str(),
            "specific_volume", "m3/mol", 1.0 / vapor.molar_density, std::get<2>(row),
            inputs, "Relationship", "historical-reference", "historical vapor volume retained without validation claim",
            "reference_or_one", 0.0, 0.05,
            close(1.0 / vapor.molar_density, std::get<2>(row), 0.05),
            "m8-vapor-volume", std::to_string(case_index).c_str(), "v_v", "T,P,\\mathbf{x}",
            "v_v", "case");
        CHECK(close(1.0 / vapor.molar_density, std::get<2>(row), 0.05));
        recording_support::compare(
            (std::string("liquid-volume/") + std::to_string(case_index)).c_str(),
            "specific_volume", "m3/mol", 1.0 / liquid.molar_density, std::get<3>(row),
            inputs, "Relationship", "historical-reference", "historical liquid volume retained without validation claim",
            "reference_or_one", 0.0, 0.05,
            close(1.0 / liquid.molar_density, std::get<3>(row), 0.05),
            "m8-liquid-volume", std::to_string(case_index).c_str(), "v_l", "T,P,\\mathbf{x}",
            "v_l", "case");
        CHECK(close(1.0 / liquid.molar_density, std::get<3>(row), 0.05));
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
    CHECK(close(vapor.pressure, 1e6, 1e-9));
    recording_support::compare(
        "supercritical-density", "molar_density", "mol/m3", liquid.molar_density,
        vapor.molar_density, supercritical_inputs, "Relationship", "density-branches",
        "supercritical branches agree", "reference_or_one", 0.0, 1e-10,
        close(liquid.molar_density, vapor.molar_density, 1e-10), "m8-supercritical-density", "state",
        "\\rho", "T,P,\\mathbf{x}", "\\rho_l=\\rho_v", "state");
    CHECK(close(liquid.molar_density, vapor.molar_density, 1e-10));
    const auto cold_vapor = ethane.state_at_pressure(
        280.0, 1e6, {1.0}, epcsaft::DensityBranch::Vapor
    );
    const auto cold_liquid = ethane.state_at_pressure(
        280.0, 1e6, {1.0}, epcsaft::DensityBranch::Liquid
    );
    const auto anchored_liquid = ethane.state_at_pressure(
        280.0, 1e6, {1.0}, epcsaft::DensityBranch::Liquid,
        cold_liquid.molar_density
    );
    const auto anchored_vapor = ethane.state_at_pressure(
        280.0, 1e6, {1.0}, epcsaft::DensityBranch::Vapor,
        cold_vapor.molar_density
    );
    const auto anchored_inputs = recording_support::fields({
        model_inputs(ethane.parameters(), 280.0, 1.0 / cold_liquid.molar_density, {1.0}),
        recording_support::field("P", 1e6),
        recording_support::field("cold_liquid_density", cold_liquid.molar_density),
        recording_support::field("cold_vapor_density", cold_vapor.molar_density),
        recording_support::field("liquid_anchor", cold_liquid.molar_density),
        recording_support::field("vapor_anchor", cold_vapor.molar_density),
    });
    recording_support::compare(
        "anchored-liquid", "molar_density", "mol/m3", anchored_liquid.molar_density,
        cold_liquid.molar_density, anchored_inputs, "Relationship", "anchored-root",
        "anchored liquid branch reproduces cold liquid root", "reference_or_one", 0.0, 1e-10,
        close(anchored_liquid.molar_density, cold_liquid.molar_density, 1e-10),
        "m8-anchored-liquid", "state", "\\rho_l", "T,P,\\mathbf{x}", "\\rho_{anchored,l}=\\rho_l", "state");
    CHECK(close(anchored_liquid.molar_density, cold_liquid.molar_density, 1e-10));
    recording_support::compare(
        "anchored-vapor", "molar_density", "mol/m3", anchored_vapor.molar_density,
        cold_vapor.molar_density, anchored_inputs, "Relationship", "anchored-root",
        "anchored vapor branch reproduces cold vapor root", "reference_or_one", 0.0, 1e-10,
        close(anchored_vapor.molar_density, cold_vapor.molar_density, 1e-10),
        "m8-anchored-vapor", "state", "\\rho_v", "T,P,\\mathbf{x}", "\\rho_{anchored,v}=\\rho_v", "state");
    CHECK(close(anchored_vapor.molar_density, cold_vapor.molar_density, 1e-10));
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
    CHECK(close(associating_liquid.pressure, 1e5, 1e-9));
    const auto nine = epcsaft::HelmholtzModel(properties_support::nine_species());
    const auto nine_state = nine.state_at_pressure(
        313.15, 1e5,
        {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
        epcsaft::DensityBranch::Liquid
    );
    const std::vector<double> nine_composition = {
        0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0
    };
    const double nine_packing = nine.phase_packing_fraction(
        313.15, nine_composition, 1.0 / nine_state.molar_density
    );
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
    CHECK(nine_packing_ok);
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
        recording_support::field("reference_convergence_error", reference.convergence_error),
    });
    recording_support::observe(
        "reference-finite", "reference_convergence_error", "1",
        reference.convergence_error, std::nullopt, reference_inputs, "Relationship", "finite-value",
        "infinite-dilution reference convergence error is finite", "m8-reference-finite", "state",
        "\\epsilon", "T,P,\\mathbf{x}", "\\epsilon", "state");
    recording_support::observe(
        "activity-finite", "mean_ionic_log_molality_scale", "1",
        activity.mean_ionic_log_molality_scale, std::nullopt, reference_inputs, "Relationship", "finite-value",
        "electrolyte activity scale is finite", "m8-activity-finite", "state", "\\log\\gamma",
        "T,P,\\mathbf{x}", "\\log\\gamma", "state");
    CHECK(std::isfinite(reference.convergence_error));
    CHECK(std::isfinite(activity.mean_ionic_log_molality_scale));
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
    const auto reference_error_inputs = recording_support::fields({
        reference_inputs, recording_support::field("ionic_strength", 1e-10)
    });
    const auto activity_scale_inputs = recording_support::fields({
        reference_inputs, recording_support::field("ionic_strength", 1e-8)
    });
    recording_support::compare(
        "reference-error", "reference_convergence_error", "1", reference.convergence_error,
        gamma_coefficient * (std::sqrt(1e-10) - std::sqrt(1e-12)),
        reference_error_inputs,
        "Relationship", "infinite-dilution", "reference convergence error matches Debye limit",
        "reference_or_one", 0.0, 1e-3,
        close(reference.convergence_error,
              gamma_coefficient * (std::sqrt(1e-10) - std::sqrt(1e-12)), 1e-3),
        "m8-reference-error", "ionic-strength", "\\epsilon", "T,P,\\mathbf{x},I",
        "\\epsilon=\\Gamma(\\sqrt I-\\sqrt{I_0})", "I");
    CHECK(close(reference.convergence_error,
                gamma_coefficient * (std::sqrt(1e-10) - std::sqrt(1e-12)),
                1e-3));
    recording_support::compare(
        "activity-scale", "mean_ionic_log_molality_scale", "1",
        activity.mean_ionic_log_molality_scale,
        -gamma_coefficient * (std::sqrt(1e-8) - std::sqrt(1e-12)),
        activity_scale_inputs,
        "Relationship", "infinite-dilution", "activity scale matches Debye limit",
        "reference_or_one", 0.0, 1e-2,
        close(activity.mean_ionic_log_molality_scale,
              -gamma_coefficient * (std::sqrt(1e-8) - std::sqrt(1e-12)), 1e-2),
        "m8-activity-scale", "ionic-strength", "\\log\\gamma", "T,P,\\mathbf{x},I",
        "\\log\\gamma= -\\Gamma(\\sqrt I-\\sqrt{I_0})", "I");
    CHECK(close(activity.mean_ionic_log_molality_scale,
                -gamma_coefficient * (std::sqrt(1e-8) - std::sqrt(1e-12)),
                1e-2));
}

std::vector<double> p2_association_residuals(
    double temperature, double volume, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites
) {
    const double density = 1.0 / volume;
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const auto delta = epcsaft::association_delta(
        temperature, density, composition, parameters, sigma, diameter
    );
    return epcsaft::association_residuals(
        density, composition, parameters, sites, delta
    );
}

struct P2AssociationConditioning final {
    Eigen::MatrixXd jacobian;
    std::size_t rank = 0;
    double pivot_ratio = 1.0;
};

P2AssociationConditioning p2_association_conditioning(
    double temperature, double density, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites
) {
    const std::size_t count = sites.size();
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const auto delta = epcsaft::association_delta(
        temperature, density, composition, parameters, sigma, diameter
    );
    const double number_density = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(count), static_cast<Eigen::Index>(count)
    );
    for (std::size_t a = 0; a < count; ++a) {
        double bonded = 0.0;
        for (std::size_t b = 0; b < count; ++b) {
            const auto component = static_cast<std::size_t>(parameters.site_component[b]);
            const double mass = composition[component] * parameters.site_multiplicity[b];
            bonded += mass * sites[b] * delta[a * count + b];
        }
        for (std::size_t b = 0; b < count; ++b) {
            const auto component = static_cast<std::size_t>(parameters.site_component[b]);
            const double mass = composition[component] * parameters.site_multiplicity[b];
            jacobian(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
                (a == b ? 1.0 + number_density * bonded : 0.0)
                + sites[a] * number_density * mass * delta[a * count + b];
        }
    }
    P2AssociationConditioning result{std::move(jacobian)};
    if (count == 0) return result;
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(result.jacobian);
    factorization.setThreshold(epcsaft::policy::ASSOCIATION_JACOBIAN_RANK.value * count);
    result.rank = static_cast<std::size_t>(factorization.rank());
    const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
    result.pivot_ratio = pivots.minCoeff() / pivots.maxCoeff();
    return result;
}

enum class P2AssociationDirection { Temperature, Volume, AssociationEnergy };

struct P2AssociationAction final {
    double residual = 0.0;
    std::size_t rank = 0;
    double pivot_ratio = 1.0;
};

P2AssociationAction p2_association_action(
    double temperature, double volume, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites,
    P2AssociationDirection direction
) {
    const auto conditioning = p2_association_conditioning(
        temperature, 1.0 / volume, composition, parameters, sites
    );
    using AD = CppAD::AD<double>;
    const double coordinate = direction == P2AssociationDirection::Temperature
        ? temperature
        : direction == P2AssociationDirection::Volume
            ? volume : parameters.association_energy_over_k[1];
    std::vector<AD> independent = {AD(coordinate)};
    CppAD::Independent(independent);
    const AD temperature_ad = direction == P2AssociationDirection::Temperature
        ? independent.front() : AD(temperature);
    const AD volume_ad = direction == P2AssociationDirection::Volume
        ? independent.front() : AD(volume);
    std::vector<AD> composition_ad(composition.begin(), composition.end());
    std::vector<AD> sites_ad(sites.begin(), sites.end());
    std::vector<epcsaft::detail::ActiveParameter> active;
    if (direction == P2AssociationDirection::AssociationEnergy) {
        active.push_back({
            epcsaft::detail::ParameterFamily::AssociationEnergyOverK, {1, 4}
        });
    }
    const auto tape_parameters = epcsaft::detail::make_tape_parameters(
        parameters, active, independent
    );
    const auto sigma = epcsaft::segment_sigma(temperature_ad, tape_parameters);
    const auto diameter = epcsaft::segment_diameter(
        temperature_ad, tape_parameters, sigma
    );
    const auto density_ad = AD(1.0) / volume_ad;
    const auto delta = epcsaft::association_delta_values(
        temperature_ad, density_ad, composition_ad, tape_parameters, sigma, diameter,
        tape_parameters.association_energy_over_k, tape_parameters.association_volume
    );
    const auto residuals = epcsaft::association_residuals(
        density_ad, composition_ad, tape_parameters, sites_ad, delta
    );
    CppAD::ADFun<double> function(independent, residuals);
    const auto jacobian = function.Jacobian(std::vector<double>{coordinate});
    Eigen::VectorXd f_direction(static_cast<Eigen::Index>(sites.size()));
    for (std::size_t site = 0; site < sites.size(); ++site) {
        f_direction(static_cast<Eigen::Index>(site)) =
            jacobian[site];
    }
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(conditioning.jacobian);
    const Eigen::VectorXd x_direction = -factorization.solve(f_direction);
    const double residual = (conditioning.jacobian * x_direction + f_direction)
        .lpNorm<Eigen::Infinity>();
    return {residual, conditioning.rank, conditioning.pivot_ratio};
}

double p2_pressure_slope_floor_reference(
    double temperature, double volume, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites
) {
    const double density = 1.0 / volume;
    const auto conditioning = p2_association_conditioning(
        temperature, density, composition, parameters, sites
    );
    if (sites.empty()) return 0.0;
    const double step = volume * 1.0e-6;
    const auto plus = p2_association_residuals(
        temperature, volume + step, composition, parameters, sites
    );
    const auto minus = p2_association_residuals(
        temperature, volume - step, composition, parameters, sites
    );
    Eigen::VectorXd potential_cross(static_cast<Eigen::Index>(sites.size()));
    for (std::size_t site = 0; site < sites.size(); ++site) {
        const auto component = static_cast<std::size_t>(parameters.site_component[site]);
        const double mass = composition[component] * parameters.site_multiplicity[site];
        const double f_volume = (plus[site] - minus[site]) / (2.0 * step);
        potential_cross(static_cast<Eigen::Index>(site)) =
            -mass / sites[site] * f_volume;
    }
    const Eigen::VectorXd ones =
        Eigen::VectorXd::Ones(static_cast<Eigen::Index>(sites.size()));
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(conditioning.jacobian);
    // The floor propagates the achieved site residual, bounded below by eps per site.
    double achieved = std::numeric_limits<double>::epsilon() * static_cast<double>(sites.size());
    for (const double value : p2_association_residuals(temperature, volume, composition, parameters, sites))
        achieved = std::max(achieved, std::abs(value));
    return std::abs(potential_cross.dot(factorization.solve(ones))) * achieved;
}

struct P2PressureChain final {
    double measured = 0.0;
    double predicted = 0.0;
    double pressure_density_derivative = 0.0;
    double reduced_residual = 0.0;
    std::size_t stable_root_count = 0;
    epcsaft::DensityBranch branch = epcsaft::DensityBranch::Unspecified;
    std::size_t association_rank = 0;
    double association_pivot_ratio = 1.0;
};

P2PressureChain p2_fixed_pressure_chain(
    const epcsaft::HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, epcsaft::DensityBranch branch,
    double step, bool enthalpy
) {
    const auto central = model.state_at_pressure(
        temperature, pressure, composition, branch
    );
    const auto plus = model.state_at_pressure(
        temperature + step, pressure, composition, branch
    );
    const auto minus = model.state_at_pressure(
        temperature - step, pressure, composition, branch
    );
    const auto fixed_plus = model.state(
        temperature + step, central.molar_density, composition
    );
    const auto fixed_minus = model.state(
        temperature - step, central.molar_density, composition
    );
    const double pressure_temperature =
        (fixed_plus.pressure - fixed_minus.pressure) / (2.0 * step);
    const double density_temperature =
        -pressure_temperature / central.pressure_density_derivative;
    const double density_step = central.molar_density * 1.0e-5;
    const auto density_plus = model.state(
        temperature, central.molar_density + density_step, composition
    );
    const auto density_minus = model.state(
        temperature, central.molar_density - density_step, composition
    );
    const auto quantity = [enthalpy](const epcsaft::State& state) {
        return enthalpy ? state.residual_enthalpy : state.molar_density;
    };
    const double quantity_temperature = quantity(fixed_plus) - quantity(fixed_minus);
    const double quantity_density = quantity(density_plus) - quantity(density_minus);
    const double predicted = quantity_temperature / (2.0 * step)
        + quantity_density / (2.0 * density_step) * density_temperature;
    const double measured = (quantity(plus) - quantity(minus)) / (2.0 * step);
    const auto conditioning = p2_association_conditioning(
        temperature, central.molar_density, composition,
        model.parameters(), central.site_fractions
    );
    return {
        measured, predicted, central.pressure_density_derivative,
        central.density_convergence.has_value()
            ? central.density_convergence->residual : std::numeric_limits<double>::infinity(),
        central.stable_root_count.value_or(0), central.density_branch,
        conditioning.rank, conditioning.pivot_ratio
    };
}

void properties_p2() {
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-support-residual", "support-residual/", 4,
        "exact-support association residual", "support");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-support-hessian", "support-hessian/", 6,
        "association trace implicit Hessian", "trace");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-support-floor", "support-floor/", 2,
        "exact-support pressure-slope residual floor", "support");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-association-action", "association-action/", 3,
        "association implicit linear residual", "direction");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-root-residual", "root-residual/", 3,
        "density root pressure residual and dense-liquid closure", "branch");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-root-convergence", "root-convergence/", 2,
        "density root reduced convergence residual", "branch");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-root-coalescing", "root-coalescing/", 1,
        "density scan coalesces a root at a grid point", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-root-metadata", "root-metadata/", 2,
        "two-root branch and count metadata", "branch");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-single-metadata", "single-metadata/", 1,
        "single-root metadata is unbranched", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-anchor", "anchor/", 3,
        "anchored local continuation", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-anchor-invalid", "anchor-invalid/", 1,
        "invalid anchor is rejected at the trust boundary", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-density-selection", "density-selection/", 1,
        "unspecified multi-root selection is ambiguous", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-fixed-pressure", "fixed-pressure/", 9,
        "fixed-pressure density implicit derivative ladder", "delta_T");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-loop-vapor", "loop-vapor/", 13,
        "declared vapor branch is the smallest stable root of a sub-cell loop", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-loop-liquid", "loop-liquid/", 13,
        "declared liquid branch is the largest stable root of a sub-cell loop", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-loop-ambiguous", "loop-ambiguous/", 13,
        "unspecified request on a resolved loop fails closed", "case");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-loop-action", "loop-action/", 2,
        "pressure-chart density action on the selected loop branch", "branch");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-fold-resolution", "fold-resolution/", 1,
        "a pressure at a mechanical fold is explicitly unresolved", "case");

    const auto associating_parameters = properties_support::associating_binary();
    epcsaft::HelmholtzModel associating(associating_parameters);
    const auto record_zero_support = [&](const char* point,
                                         const epcsaft::HelmholtzModel& model,
                                         double temperature_value, double density_value,
                                         const std::vector<double>& composition) {
        const auto state = model.state(temperature_value, density_value, composition);
        const auto residuals = p2_association_residuals(
            temperature_value, 1.0 / density_value, composition,
            model.parameters(), state.site_fractions
        );
        double residual = 0.0;
        bool finite = finite_state_without_cp(state);
        for (const double value : state.site_fractions) finite = finite && std::isfinite(value);
        for (const double value : residuals) residual = std::max(residual, std::abs(value));
        recording_support::compare_custom(
            (std::string("support-residual/") + point).c_str(),
            "association_residual_infinity", "1", residual, 0.0,
            recording_support::fields({
                model_inputs(model.parameters(), temperature_value, 1.0 / density_value, composition),
                recording_support::field("support", point),
                recording_support::field("site_fractions", state.site_fractions),
            }), "Relationship", "mass-action",
            "exact neutral/zero-ion support remains finite and converged", "absolute",
            residual <= 2e-13 && finite, "p2-support-residual", point,
            "||F_a||_inf", "T,rho,x", "||F_a||_inf<=2e-13", "support");
        CHECK(finite && residual <= 2e-13);
    };
    const std::array<std::vector<double>, 2> supports = {{{1.0, 0.0}, {0.0, 1.0}}};
    for (std::size_t support = 0; support < supports.size(); ++support) {
        const auto& composition = supports[support];
        constexpr double temperature = 300.0;
        constexpr double density = 5000.0;
        const auto state = associating.state(temperature, density, composition);
        const auto residuals = p2_association_residuals(
            temperature, 1.0 / density, composition,
            associating.parameters(), state.site_fractions
        );
        double residual = 0.0;
        bool finite = finite_state_without_cp(state);
        for (const double value : state.site_fractions) finite = finite && std::isfinite(value);
        for (const double value : residuals) residual = std::max(residual, std::abs(value));
        const auto pressure_slope = associating.evaluate(
            temperature, 1.0 / density, composition, epcsaft::DerivativeOrder::PressureSlope
        );
        const double expected_floor = p2_pressure_slope_floor_reference(
            temperature, 1.0 / density, composition,
            associating.parameters(), state.site_fractions
        );
        // The PressureSlope adjoint gradient and the Gradient-order sensitivity gradient are
        // algebraically distinct formulas for the same reduced gradient.
        const auto gradient_order = associating.evaluate(
            temperature, 1.0 / density, composition, epcsaft::DerivativeOrder::Gradient);
        double gradient_error = 0.0;
        for (std::size_t i = 0; i < gradient_order.gradient.size(); ++i)
            gradient_error = std::max(gradient_error, std::abs(pressure_slope.gradient[i] - gradient_order.gradient[i])
                / std::max(1.0, std::abs(gradient_order.gradient[i])));
        CHECK(pressure_slope.gradient.size() == gradient_order.gradient.size() && gradient_error <= 1e-12);
        const auto inputs = recording_support::fields({
            model_inputs(associating.parameters(), temperature, 1.0 / density, composition),
            recording_support::field("support", support),
            recording_support::field("pressure_slope_gradient_identity_error", gradient_error),
            recording_support::field("site_fractions", state.site_fractions),
            recording_support::field("association_floor", pressure_slope.gradient_floor),
            recording_support::field("expected_floor", expected_floor),
        });
        const auto point = std::to_string(support);
        recording_support::compare_custom(
            (std::string("support-residual/") + point).c_str(),
            "association_residual_infinity", "1", residual, 0.0, inputs,
            "Relationship", "mass-action", "exact support remains finite and converged",
            "absolute", residual <= 2e-13 && finite,
            "p2-support-residual", point.c_str(), "||F_a||_inf", "T,rho,x",
            "||F_a||_inf<=2e-13", "support");
        CHECK(finite && residual <= 2e-13);
        const double floor_scale = std::max(expected_floor, 1.0e-300);
        const double floor_error = std::abs(pressure_slope.gradient_floor - expected_floor);
        const bool floor_ok = floor_error <= 1.0e-6 * floor_scale;
        recording_support::compare_custom(
            (std::string("support-floor/") + point).c_str(), "gradient_floor", "1",
            pressure_slope.gradient_floor, expected_floor, inputs, "Relationship",
            "mass-action-residual-bound", "PressureSlope floor uses U^{-1} one",
            "relative", floor_ok, "p2-support-floor", point.c_str(), "floor",
            "T,rho,x", "floor=bound(U^{-1}1)", "support");
        CHECK(floor_ok);
    }
    record_zero_support(
        "aqueous-zero-ion", epcsaft::HelmholtzModel(properties_support::aqueous_electrolyte()),
        298.15, 5000.0, {1.0, 0.0, 0.0}
    );
    record_zero_support(
        "nine-zero-ion", epcsaft::HelmholtzModel(properties_support::nine_species()),
        313.15, 5000.0, {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0}
    );

    const std::array<double, 5> association_trace = {
        1e-2, 1e-3, 1e-4, 1e-5, 0.0
    };
    std::array<double, 5> trace_to_zero_errors{};
    const auto zero_derivatives = associating.evaluate(
        300.0, 1.0 / 5000.0, {1.0, 0.0}, epcsaft::DerivativeOrder::Hessian
    );
    for (std::size_t trace_index = 0; trace_index < association_trace.size(); ++trace_index) {
        const double trace = association_trace[trace_index];
        const std::vector<double> composition = {1.0 - trace, trace};
        const auto state = associating.state(300.0, 5000.0, composition);
        const auto derivatives = associating.evaluate(
            300.0, 1.0 / 5000.0, composition, epcsaft::DerivativeOrder::Hessian
        );
        const auto reference = properties_support::reference_tensor(
            associating.parameters(), 300.0, 1.0 / 5000.0, composition,
            state.site_fractions
        );
        double hessian_error = 0.0;
        double trace_to_zero_error = 0.0;
        for (std::size_t index = 0; index < derivatives.hessian.size(); ++index) {
            hessian_error = std::max(
                hessian_error,
                std::abs(derivatives.hessian[index] - reference.hessian[index])
                    / std::max(1.0, std::abs(reference.hessian[index]))
            );
            trace_to_zero_error = std::max(
                trace_to_zero_error,
                std::abs(derivatives.hessian[index] - zero_derivatives.hessian[index])
                    / std::max(1.0, std::abs(zero_derivatives.hessian[index]))
            );
        }
        trace_to_zero_errors[trace_index] = trace_to_zero_error;
        const auto point = std::to_string(trace);
        recording_support::compare(
            (std::string("support-hessian/") + point).c_str(), "relative_hessian_error", "1",
            hessian_error, 0.0,
            recording_support::fields({
                model_inputs(associating.parameters(), 300.0, 1.0 / 5000.0, composition),
                recording_support::field("trace_amount", trace),
                recording_support::field("trace_to_zero_error", trace_to_zero_error),
                recording_support::field("scaled_trace_to_zero_error",
                    trace == 0.0 ? 0.0 : trace_to_zero_error / trace),
            }), "Relationship", "AD-reference",
            "association implicit Hessian remains support-reduced", "absolute", 1e-12, 0.0,
            hessian_error <= 1e-12, "p2-support-hessian", point.c_str(), "Hres",
            "T,rho,x", "Himplicit=Hreference", "trace");
        CHECK(hessian_error <= 1e-12);
    }
    const bool trace_refines = trace_to_zero_errors[0] >= trace_to_zero_errors[1]
        && trace_to_zero_errors[1] >= trace_to_zero_errors[2]
        && trace_to_zero_errors[2] >= trace_to_zero_errors[3]
        && trace_to_zero_errors[3] >= trace_to_zero_errors[4];
    recording_support::compare_custom(
        "support-hessian/trace-to-zero", "refinement", "bool",
        trace_refines ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("trace_errors", std::vector<double>(
                trace_to_zero_errors.begin(), trace_to_zero_errors.end())),
            recording_support::field("terminal_error", trace_to_zero_errors.back()),
            recording_support::field("scaled_terminal_error",
                trace_to_zero_errors[3] / association_trace[3]),
        }), "Relationship", "scaled-trace-ladder",
        "trace-to-zero Hessian error refines to the exact-zero row", "absolute", trace_refines,
        "p2-support-hessian", "trace-to-zero", "refinement", "trace",
        "error(trace)<=error(previous);error(0)=0", "trace");
    CHECK(trace_refines);

    const std::vector<double> action_composition = {0.4, 0.6};
    const std::array<P2AssociationDirection, 3> directions = {
        P2AssociationDirection::Temperature,
        P2AssociationDirection::Volume,
        P2AssociationDirection::AssociationEnergy
    };
    const std::array<const char*, 3> direction_names = {
        "temperature", "volume", "association-energy"
    };
    const std::array<double, 3> direction_steps = {
        1e-3, 1.0 / 5000.0 * 1e-3, 1e-3
    };
    const auto action_state = associating.state(300.0, 5000.0, action_composition);
    for (std::size_t index = 0; index < directions.size(); ++index) {
        std::vector<epcsaft::ActiveParameter> active;
        if (directions[index] == P2AssociationDirection::AssociationEnergy) {
            active.push_back({epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}});
        }
        const epcsaft::HelmholtzModel production_model(
            associating_parameters, active
        );
        const auto production = production_model.evaluate(
            300.0, 1.0 / 5000.0, action_composition, epcsaft::DerivativeOrder::Hessian
        );
        const auto reference = properties_support::reference_tensor(
            associating_parameters, 300.0, 1.0 / 5000.0, action_composition,
            action_state.site_fractions, active
        );
        const std::size_t coordinate = directions[index] == P2AssociationDirection::Temperature
            ? 0 : directions[index] == P2AssociationDirection::Volume
                ? 1 : 2 + action_composition.size();
        const double production_error = std::abs(
            production.hessian[coordinate * production.gradient.size() + coordinate]
            - reference.hessian[coordinate * reference.gradient.size() + coordinate]
        ) / std::max(1.0, std::abs(reference.hessian[coordinate * reference.gradient.size() + coordinate]));
        const auto action = p2_association_action(
            300.0, 1.0 / 5000.0, action_composition,
            associating.parameters(), action_state.site_fractions,
            directions[index]
        );
        const bool action_ok = production_error <= 1e-10
            && action.rank == action_state.site_fractions.size()
            && std::isfinite(action.pivot_ratio)
            && action.pivot_ratio > epcsaft::policy::ASSOCIATION_JACOBIAN_CONDITION.value
            && std::isfinite(action.residual) && action.residual <= 1e-10;
        const auto inputs = recording_support::fields({
            model_inputs(associating.parameters(), 300.0, 1.0 / 5000.0, action_composition),
            recording_support::field("direction", direction_names[index]),
            recording_support::field("production_hessian_coordinate", coordinate),
            recording_support::field("production_reference_error", production_error),
            recording_support::field("step", direction_steps[index]),
            recording_support::field("rank", action.rank),
            recording_support::field("pivot_ratio", action.pivot_ratio),
        });
        recording_support::compare_custom(
            (std::string("association-action/") + direction_names[index]).c_str(),
            "production_action_error", "1", production_error, 0.0, inputs, "Relationship",
            "production-vs-raw-residual", "production Hessian matches independent raw-residual elimination",
            "absolute", action_ok, "p2-association-action", direction_names[index],
            "H_q,||UX_q+F_q||inf", "T,rho,V,theta", "H_production=H_raw;||UX_q+F_q||inf<=1e-10", "direction");
        CHECK(action_ok);
    }

    const auto dense_state = associating.state_at_pressure(
        353.15, 1e5, {0.25, 0.75}, epcsaft::DensityBranch::Liquid);
    const double dense_pressure_closure = std::abs(dense_state.pressure - 1e5) / 1e5;
    recording_support::compare(
        "root-residual/associating-dense", "relative_pressure_residual", "1",
        dense_pressure_closure, 0.0, recording_support::fields({
            model_inputs(associating.parameters(), 353.15, 1.0 / dense_state.molar_density, {0.25, 0.75}),
            recording_support::field("P_star", 1e5),
            recording_support::field("P_published", dense_state.pressure),
        }), "Relationship", "density-publication",
        "dense-liquid publication closes the imposed pressure", "absolute", 1e-10, 0.0,
        dense_pressure_closure <= 1e-10, "p2-root-residual", "associating-dense",
        "|P(rho)-P_star|/P_star", "T,P,x,liquid", "|P(rho)-P_star|/P_star<=1e-10", "branch");
    CHECK(dense_pressure_closure <= 1e-10);

    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    constexpr double temperature = 200.0;
    constexpr double pressure = 217800.69494692178;
    const std::vector<double> composition = {1.0};
    const auto vapor_root = epcsaft::solve_density(
        ethane, temperature, pressure, composition, epcsaft::DensityBranch::Vapor, std::nullopt
    );
    const auto liquid_root = epcsaft::solve_density(
        ethane, temperature, pressure, composition, epcsaft::DensityBranch::Liquid, std::nullopt
    );
    const auto vapor_state = ethane.state_at_pressure(
        temperature, pressure, composition, epcsaft::DensityBranch::Vapor
    );
    const auto liquid_state = ethane.state_at_pressure(
        temperature, pressure, composition, epcsaft::DensityBranch::Liquid
    );
    const std::array<const epcsaft::DensityRoot*, 2> roots = {&vapor_root, &liquid_root};
    const std::array<const epcsaft::State*, 2> states = {&vapor_state, &liquid_state};
    const std::array<const char*, 2> names = {"vapor", "liquid"};
    const std::array<epcsaft::DensityBranch, 2> branches = {
        epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid
    };
    for (std::size_t index = 0; index < roots.size(); ++index) {
        const auto& root = *roots[index];
        const auto& state = *states[index];
        const auto root_state = ethane.state(temperature, root.molar_density, composition);
        const double root_pressure = root_state.pressure - pressure;
        const double pressure_error = std::abs(root.pressure_residual - root_pressure);
        const auto inputs = recording_support::fields({
            model_inputs(ethane.parameters(), temperature, 1.0 / root.molar_density, composition),
            recording_support::field("P_star", pressure),
            recording_support::field("requested_branch", names[index]),
            recording_support::field("stable_root_count", root.stable_root_count.value_or(0)),
            recording_support::field("density", root.molar_density),
            recording_support::field("root_pressure_residual_Pa", root.pressure_residual),
            recording_support::field("reduced_convergence_residual", root.convergence.residual),
            recording_support::field("pressure_density_derivative", root_state.pressure_density_derivative),
            recording_support::field("D", root_state.pressure_density_derivative
                / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature)),
            recording_support::field("stop", epcsaft::core::stop_name(root.convergence.stop)),
            recording_support::field("criterion", root.convergence.criterion->name),
            recording_support::field("association_rank", std::size_t{0}),
            recording_support::field("association_pivot_ratio", 1.0),
        });
        recording_support::compare(
            (std::string("root-residual/") + names[index]).c_str(), "pressure_residual", "Pa",
            root.pressure_residual, root_pressure, inputs, "Identity", "state-replay",
            "DensityRoot pressure residual is P(root)-P", "absolute", 1e-7, 0.0,
            pressure_error <= 1e-7, "p2-root-residual", names[index], "P_{root}-P",
            "T,P,x,branch", "rP=P(root)-P", "branch");
        CHECK(pressure_error <= 1e-7);
        recording_support::compare(
            (std::string("root-convergence/") + names[index]).c_str(), "reduced_residual", "1",
            root.convergence.residual, 0.0, inputs, "Availability", "density-policy",
            "density convergence uses the reduced pressure residual", "absolute", 1e-12, 0.0,
            root.convergence.residual <= 1e-12, "p2-root-convergence", names[index],
            "r_reduced", "T,P,x,branch", "r_reduced<=1e-12", "branch");
        CHECK(root.convergence.residual <= 1e-12);
        const double reuse_error = std::abs(state.molar_density - root.molar_density);
        const bool metadata = root.branch == branches[index]
            && root.stable_root_count == std::optional<std::size_t>{2}
            && state.density_branch == branches[index]
            && state.stable_root_count == std::optional<std::size_t>{2}
            && reuse_error <= 1e-12 * std::max(1.0, root.molar_density);
        const auto metadata_inputs = recording_support::fields({
            inputs,
            recording_support::field("state_branch", names[index]),
            recording_support::field("state_count", state.stable_root_count.value_or(0)),
            recording_support::field("state_density", state.molar_density),
            recording_support::field("density_reuse_error", reuse_error),
        });
        recording_support::compare(
            (std::string("root-metadata/") + names[index]).c_str(), "molar_density", "mol/m3",
            state.molar_density, root.molar_density, metadata_inputs, "Identity",
            "state-propagation", "selected root and State retain density, branch and count",
            "reference_or_one", 0.0, 1e-12, metadata, "p2-root-metadata", names[index],
            "rho", "T,P,x,branch", "rho_state=rho_root;branch=count", "branch");
        CHECK(metadata);
    }

    bool ambiguous = false;
    try {
        (void)epcsaft::solve_density(
            ethane, temperature, pressure, composition, epcsaft::DensityBranch::Unspecified,
            std::nullopt
        );
    } catch (const std::runtime_error& error) {
        ambiguous = std::string(error.what()).find("ambiguous") != std::string::npos;
    }
    recording_support::compare_custom(
        "density-selection/ambiguous", "selected", "bool", ambiguous ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("P_star", pressure),
            recording_support::field("requested_branch", "unspecified"),
            recording_support::field("stable_root_count", 2U),
            recording_support::field("outcome", "ambiguous"),
        }), "Availability", "bounded-scan", "unspecified multi-root request fails closed",
        "absolute", ambiguous, "p2-density-selection", "ambiguous", "selected",
        "T,P,x,branch,count", "no selected density", "case");
    CHECK(ambiguous);

    constexpr double coalescing_temperature = 320.0;
    const double coalescing_upper = ethane.molar_density_upper_bound(
        coalescing_temperature, composition
    ) * epcsaft::policy::DENSITY_UPPER_FRACTION.value;
    double coalescing_pressure = 1e6;
    double coalescing_density = 0.0;
    for (std::size_t iteration = 0; iteration < 40; ++iteration) {
        const double ideal_density = coalescing_pressure
            / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * coalescing_temperature);
        const double lower_density = ideal_density
            * epcsaft::policy::DENSITY_LOWER_FRACTION.value;
        coalescing_density = std::exp(
            std::log(lower_density) + 0.5
            * (std::log(coalescing_upper) - std::log(lower_density))
        );
        coalescing_pressure = ethane.state(
            coalescing_temperature, coalescing_density, composition
        ).pressure;
    }
    const auto coalescing_root = epcsaft::solve_density(
        ethane, coalescing_temperature, coalescing_pressure, composition,
        epcsaft::DensityBranch::Vapor,
        std::nullopt
    );
    const double coalescing_alignment = std::abs(
        coalescing_root.molar_density - coalescing_density
    ) / std::max(1.0, coalescing_density);
    const bool coalesced = close(
        coalescing_root.molar_density, coalescing_density, 1e-12
    ) && coalescing_root.stable_root_count == std::optional<std::size_t>{1};
    recording_support::compare_custom(
        "root-coalescing/grid", "grid_alignment", "1", coalescing_alignment, 0.0,
        recording_support::fields({
            recording_support::field("T", coalescing_temperature),
            recording_support::field("P_star", coalescing_pressure),
            recording_support::field("grid_index", std::size_t{20}),
            recording_support::field("grid_density", coalescing_density),
            recording_support::field("root_density", coalescing_root.molar_density),
            recording_support::field("stable_root_count",
                coalescing_root.stable_root_count.value_or(0)),
        }), "Relationship", "bounded-scan", "grid-point root is coalesced once",
        "absolute", coalesced, "p2-root-coalescing", "grid", "rho",
        "T,P,x,grid_index", "rho_root=rho_grid", "case");
    CHECK(coalesced);

    const auto supercritical_vapor = ethane.state_at_pressure(
        320.0, 1e6, composition, epcsaft::DensityBranch::Vapor
    );
    const auto supercritical_liquid = ethane.state_at_pressure(
        320.0, 1e6, composition, epcsaft::DensityBranch::Liquid
    );
    const auto supercritical_root = epcsaft::solve_density(
        ethane, 320.0, 1e6, composition, epcsaft::DensityBranch::Vapor, std::nullopt
    );
    const auto supercritical_inputs = recording_support::fields({
        model_inputs(ethane.parameters(), 320.0, 1.0 / supercritical_vapor.molar_density, composition),
        recording_support::field("P", 1e6),
        recording_support::field("vapor_density", supercritical_vapor.molar_density),
        recording_support::field("liquid_density", supercritical_liquid.molar_density),
        recording_support::field("root_density", supercritical_root.molar_density),
        recording_support::field("root_branch", "unspecified"),
        recording_support::field("root_count", supercritical_root.stable_root_count.value_or(0)),
        recording_support::field("state_branch", "unspecified"),
        recording_support::field("state_count", supercritical_vapor.stable_root_count.value_or(0)),
    });
    recording_support::compare(
        "single-metadata/root", "molar_density", "mol/m3", supercritical_root.molar_density,
        supercritical_vapor.molar_density, supercritical_inputs, "Relationship", "single-root",
        "one stable root is shared by DensityRoot and both State requests", "reference_or_one", 0.0, 1e-12,
        supercritical_root.branch == epcsaft::DensityBranch::Unspecified
            && supercritical_root.stable_root_count == std::optional<std::size_t>{1}
            && supercritical_vapor.density_branch == epcsaft::DensityBranch::Unspecified
            && supercritical_vapor.stable_root_count == std::optional<std::size_t>{1}
            && supercritical_liquid.density_branch == epcsaft::DensityBranch::Unspecified
            && supercritical_liquid.stable_root_count == std::optional<std::size_t>{1}
            && close(supercritical_root.molar_density, supercritical_vapor.molar_density, 1e-12),
        "p2-single-metadata", "root", "rho", "T,P,x", "rho_root=rho_state;branch=count", "state");

    const auto cold_vapor = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Vapor
    );
    const auto cold_liquid = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Liquid
    );
    const auto anchored_vapor = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Vapor,
        cold_vapor.molar_density
    );
    const auto anchored_liquid = ethane.state_at_pressure(
        280.0, 1e6, composition, epcsaft::DensityBranch::Liquid,
        cold_liquid.molar_density
    );
    const auto anchor_inputs = recording_support::fields({
        recording_support::field("T", 280.0), recording_support::field("P", 1e6),
        recording_support::field("requested_branch", "vapor/liquid"),
        recording_support::field("returned_branch", "vapor/liquid"),
        recording_support::field("stable_root_count", "absent"),
        recording_support::field("envelope_log_half_width",
            epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value),
        recording_support::field("scan_intervals",
            static_cast<std::size_t>(epcsaft::policy::DENSITY_ANCHOR_SCAN_INTERVALS.value)),
        recording_support::field("cold_vapor_density", cold_vapor.molar_density),
        recording_support::field("cold_liquid_density", cold_liquid.molar_density),
        recording_support::field("vapor_density", anchored_vapor.molar_density),
        recording_support::field("liquid_density", anchored_liquid.molar_density),
        recording_support::field("vapor_reduced_residual",
            anchored_vapor.density_convergence->residual),
        recording_support::field("liquid_reduced_residual",
            anchored_liquid.density_convergence->residual),
        recording_support::field("association_rank", std::size_t{0}),
        recording_support::field("association_pivot_ratio", 1.0),
    });
    recording_support::compare(
        "anchor/vapor", "molar_density", "mol/m3", anchored_vapor.molar_density,
        cold_vapor.molar_density, anchor_inputs, "Relationship", "local-continuation",
        "vapor anchor replays the local vapor root", "reference_or_one", 0.0, 1e-12,
        close(anchored_vapor.molar_density, cold_vapor.molar_density, 1e-12),
        "p2-anchor", "vapor", "rho_v", "T,P,x",
        "rho_anchor_v=rho_v", "case");
    CHECK(close(anchored_vapor.molar_density, cold_vapor.molar_density, 1e-12));
    recording_support::compare(
        "anchor/liquid", "molar_density", "mol/m3", anchored_liquid.molar_density,
        cold_liquid.molar_density, anchor_inputs, "Relationship", "local-continuation",
        "liquid anchor replays the local liquid root", "reference_or_one", 0.0, 1e-12,
        close(anchored_liquid.molar_density, cold_liquid.molar_density, 1e-12),
        "p2-anchor", "liquid", "rho_l", "T,P,x",
        "rho_anchor_l=rho_l", "case");
    CHECK(close(anchored_liquid.molar_density, cold_liquid.molar_density, 1e-12));

    const auto association_anchor = associating.state_at_pressure(
        376.88, 101330.0, {0.1385, 0.8615}, epcsaft::DensityBranch::Liquid
    );
    bool branch_unreachable = false;
    bool typed_failure = false;
    double failure_envelope = 0.0;
    std::size_t failure_intervals = 0;
    epcsaft::core::Stop failure_stop = epcsaft::core::Stop::stalled;
    double failure_residual = 0.0;
    try {
        (void)associating.state_at_pressure(
            376.88, 101330.0, {0.99, 0.01}, epcsaft::DensityBranch::Liquid,
            association_anchor.molar_density
        );
    } catch (const epcsaft::core::ConvergenceFailure& error) {
        typed_failure = true;
        branch_unreachable = error.record.stop == epcsaft::core::Stop::left_domain;
        failure_envelope = epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value;
        failure_intervals = static_cast<std::size_t>(
            epcsaft::policy::DENSITY_ANCHOR_SCAN_INTERVALS.value);
        failure_stop = error.record.stop;
        failure_residual = error.record.residual;
    } catch (const std::runtime_error&) {
        branch_unreachable = false;
    }
    const auto unreachable_inputs = recording_support::fields({
        recording_support::field("T", 376.88), recording_support::field("P", 101330.0),
        recording_support::field("anchor_x", std::vector<double>{0.1385, 0.8615}),
        recording_support::field("target_x", std::vector<double>{0.99, 0.01}),
        recording_support::field("anchor_density", association_anchor.molar_density),
        recording_support::field("typed_failure", typed_failure ? 1U : 0U),
        recording_support::field("stop", epcsaft::core::stop_name(failure_stop)),
        recording_support::field("residual", failure_residual),
        recording_support::field("envelope_log_half_width", failure_envelope),
        recording_support::field("scan_intervals", failure_intervals),
    });
    recording_support::compare_custom(
        "anchor/unreachable", "branch_reachable", "bool", branch_unreachable ? 1.0 : 0.0, 1.0,
        unreachable_inputs, "Availability", "local-continuation", "unreachable anchored branch fails closed",
        "absolute", branch_unreachable && typed_failure
            && failure_envelope == epcsaft::policy::DENSITY_ANCHOR_LOG_WINDOW.value
            && failure_intervals == static_cast<std::size_t>(
                epcsaft::policy::DENSITY_ANCHOR_SCAN_INTERVALS.value),
        "p2-anchor", "unreachable", "reachable",
        "T,P,x_anchor,x_target", "no cold fallback", "case");
    CHECK(branch_unreachable && typed_failure);

    bool invalid_anchor = false;
    try {
        (void)epcsaft::solve_density(
            ethane, temperature, pressure, composition, epcsaft::DensityBranch::Vapor, -1.0
        );
    } catch (const std::invalid_argument&) {
        invalid_anchor = true;
    }
    recording_support::compare_custom(
        "anchor-invalid/negative", "accepted", "bool", invalid_anchor ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("T", temperature),
            recording_support::field("P_star", pressure),
            recording_support::field("anchor_density", -1.0),
            recording_support::field("outcome", "invalid-anchor"),
        }), "Availability", "trust-boundary", "nonpositive anchor is rejected",
        "absolute", invalid_anchor, "p2-anchor-invalid", "negative", "accepted",
        "T,P,x,anchor", "anchor rejected", "case");
    CHECK(invalid_anchor);

    struct P2ChainCase final {
        const epcsaft::HelmholtzModel* model;
        double temperature;
        double pressure;
        std::vector<double> composition;
        epcsaft::DensityBranch branch;
        const char* branch_name;
        const char* point;
        const char* quantity;
        const char* unit;
        bool enthalpy;
    };
    const std::array<double, 3> delta_temperatures = {1e-2, 3e-3, 1e-3};
    const std::array<P2ChainCase, 3> chain_cases = {{
        {&ethane, 320.0, 1e6, composition, epcsaft::DensityBranch::Vapor,
            "vapor", "one-root-density", "density_temperature_derivative", "mol/m3/K", false},
        {&ethane, 200.0, pressure, composition, epcsaft::DensityBranch::Liquid,
            "liquid", "multi-root-enthalpy", "residual_enthalpy_temperature_derivative", "J/mol/K", true},
        {&associating, 313.15, 1e5, {0.1385, 0.8615}, epcsaft::DensityBranch::Liquid,
            "liquid", "associating-enthalpy", "residual_enthalpy_temperature_derivative", "J/mol/K", true},
    }};
    const auto record_chain = [&](const char* point, const char* quantity,
                                  const char* unit, double temperature_value,
                                  double pressure_value,
                                  const char* requested_branch,
                                  double delta_temperature,
                                  const P2PressureChain& chain) {
        const bool ok = close(chain.measured, chain.predicted, 5e-6)
            && chain.pressure_density_derivative > 0.0
            && std::isfinite(chain.reduced_residual);
        recording_support::compare(
            (std::string("fixed-pressure/") + point).c_str(), quantity, unit,
            chain.measured, chain.predicted,
            recording_support::fields({
                recording_support::field("T", temperature_value),
                recording_support::field("P_star", pressure_value),
                recording_support::field("requested_branch", requested_branch),
                recording_support::field("returned_branch",
                    chain.branch == epcsaft::DensityBranch::Unspecified
                        ? "unspecified" : requested_branch),
                recording_support::field("stable_root_count", chain.stable_root_count),
                recording_support::field("P_rho", chain.pressure_density_derivative),
                recording_support::field("D", chain.pressure_density_derivative
                    / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature_value)),
                recording_support::field("reduced_residual", chain.reduced_residual),
                recording_support::field("association_rank", chain.association_rank),
                recording_support::field("association_pivot_ratio",
                    chain.association_pivot_ratio),
                recording_support::field("delta_T", delta_temperature),
            }), "Relationship", "implicit-chain",
            "fixed-pressure property chain includes density response",
            "reference_or_one", 0.0, 5e-6, ok, "p2-fixed-pressure", point,
            quantity, "T,P,x,branch,delta_T", "Q_P=Q_T+Q_rho rho_P", "case");
        CHECK(ok);
    };
    for (const auto& item : chain_cases) {
        for (const double delta_temperature : delta_temperatures) {
            const auto chain = p2_fixed_pressure_chain(
                *item.model, item.temperature, item.pressure, item.composition,
                item.branch, delta_temperature, item.enthalpy
            );
            const std::string point = std::string(item.point) + "/"
                + std::to_string(delta_temperature);
            record_chain(point.c_str(), item.quantity, item.unit, item.temperature,
                item.pressure, item.branch_name, delta_temperature, chain);
        }
    }

    // #89: a van der Waals loop narrower than one 40-cell scan interval. The reference is
    // independent of solve_density: sign changes of P(T,rho)-P* on a 4000-point log grid at
    // fixed density, bisected to machine precision, classified by the analytic slope.
    // The 300 K state sits within 0.05 % of the liquid spinodal pressure: the liquid root lies
    // 1.1 % above the unstable one and the residual dip between them is far narrower than a cell.
    // The eighth state is 1e-4 K below the critical point: the negative-slope dip is about 0.002 wide
    // in log density, narrower than any fixed search width, and the pressure extrema differ by 0.021 Pa.
    // The ninth lies 7.4 Pa above it, near the vapor-side spinodal pressure: the vapor and unstable
    // roots are 3 mol/m3 apart (hence the denser reference grid) and the vapor root's residual
    // tolerance band is 0.05 mol/m3 wide, where a shared cell endpoint inside the band was once
    // reported as a third stable root. Near-tangent roots are compared at the solver's first-order
    // accuracy tol |d rho/d r| when that exceeds 1e-7 relative, and that accuracy must stay far
    // below the branch separation.
    // The final three inputs put the liquid root exactly on the former scan grid and one
    // pressure ULP on either side. A grid coincidence must not erase or duplicate a branch.
    struct LoopCase final {
        double temperature, pressure; std::size_t reference_points; const char* suffix = "";
        double reference_lower = 1000.0, reference_upper = 0.0;
    };
    const std::array<LoopCase, 13> loop_cases = {{
        {300.0, 3916073.568, 4000}, {305.0, 4758116.154, 4000}, {306.0, 4862786.389, 4000},
        {307.0, 4965901.440, 4000}, {308.0, 5067464.963, 4000}, {308.8, 5147500.0, 4000},
        {308.9, 5157549.110, 4000}, {308.9554032139963, 5163053.740696311, 4000},
        {308.95547769995335, 5163061.1390785752, 16000},
        {308.0, 5057915.461281232, 4000, "/grid"},
        {308.0, std::nextafter(5057915.461281232, 0.0), 4000, "/grid-below"},
        {308.0, std::nextafter(5057915.461281232, 1e7), 4000, "/grid-above"},
        // Residual-only stopping overshot its local sensitivity radius on the liquid branch.
        {308.9554732140194, 5163060.690898801, 1000, "/radius", 6380.0, 6400.0}}};
    for (const auto& item : loop_cases) {
        const double reference_upper = item.reference_upper > 0.0 ? item.reference_upper
            : ethane.molar_density_upper_bound(item.temperature, composition)
                * epcsaft::policy::DENSITY_UPPER_FRACTION.value;
        const auto residual_at = [&](double density) {
            return ethane.state(item.temperature, density, composition).pressure - item.pressure;
        };
        std::vector<std::pair<double, double>> references;
        double previous_density = item.reference_lower, previous_residual = residual_at(previous_density);
        for (std::size_t index = 1; index <= item.reference_points; ++index) {
            const double density = item.reference_lower * std::pow(reference_upper / item.reference_lower,
                                                     static_cast<double>(index) / static_cast<double>(item.reference_points));
            const double residual = residual_at(density);
            if (previous_residual * residual < 0.0) {
                double lower = previous_density, upper = density, lower_residual = previous_residual;
                for (std::size_t bisection = 0; bisection < 80; ++bisection) {
                    const double middle = 0.5 * (lower + upper), value = residual_at(middle);
                    if (lower_residual * value <= 0.0) upper = middle; else { lower = middle; lower_residual = value; }
                }
                const double root = 0.5 * (lower + upper);
                references.emplace_back(root, ethane.state(item.temperature, root, composition).pressure_density_derivative);
            }
            previous_density = density;
            previous_residual = residual;
        }
        CHECK(references.size() == 3 && references[0].second > 0.0 && references[1].second < 0.0
              && references[2].second > 0.0);
        const std::string label = std::to_string(item.temperature) + item.suffix;
        const auto loop_inputs = recording_support::fields({
            recording_support::field("T", item.temperature), recording_support::field("P_star", item.pressure),
            recording_support::field("reference_vapor", references[0].first),
            recording_support::field("reference_unstable", references[1].first),
            recording_support::field("reference_liquid", references[2].first),
            recording_support::field("loop_log_width", std::log(references[2].first / references[0].first)),
        });
        const std::array<std::pair<epcsaft::DensityBranch, const char*>, 2> branches = {{
            {epcsaft::DensityBranch::Vapor, "loop-vapor/"}, {epcsaft::DensityBranch::Liquid, "loop-liquid/"}}};
        for (std::size_t side = 0; side < 2; ++side) {
            const auto root = epcsaft::solve_density(
                ethane, item.temperature, item.pressure, composition, branches[side].first, std::nullopt);
            const double reference = references[side == 0 ? 0 : 2].first;
            const double accuracy = std::max(1e-7 * reference, epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                * root.convergence.amplification.value_or(0.0));
            CHECK(accuracy < 1e-2 * (references[2].first - references[0].first));
            const bool ok = std::abs(root.molar_density - reference) <= accuracy && root.branch == branches[side].first
                && root.stable_root_count == std::optional<std::size_t>{2} && root.pressure_density_derivative > 0.0
                && root.convergence.stop == epcsaft::core::Stop::converged;
            recording_support::compare((branches[side].second + label).c_str(), "molar_density", "mol/m3",
                root.molar_density, reference, recording_support::fields({loop_inputs,
                    recording_support::field("stable_root_count", root.stable_root_count.value_or(0)),
                    recording_support::field("pressure_density_derivative", root.pressure_density_derivative),
                    recording_support::field("first_order_accuracy", epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                        * root.convergence.amplification.value_or(0.0))}),
                "Relationship", "independent-bisection", "declared branch selects the extreme stable root",
                "absolute", accuracy, 0.0, ok, side == 0 ? "p2-loop-vapor" : "p2-loop-liquid",
                label.c_str(), "rho", "T,P,x,branch", "rho_branch=rho_ref", "case");
            CHECK(ok);
        }
        bool loop_ambiguous = false;
        try {
            (void)epcsaft::solve_density(ethane, item.temperature, item.pressure, composition,
                                         epcsaft::DensityBranch::Unspecified, std::nullopt);
        } catch (const std::runtime_error& error) {
            loop_ambiguous = std::string(error.what()).find("ambiguous") != std::string::npos;
        }
        recording_support::compare_custom(("loop-ambiguous/" + label).c_str(), "selected", "bool",
            loop_ambiguous ? 1.0 : 0.0, 1.0, loop_inputs, "Availability", "bounded-scan",
            "unspecified request on two stable roots fails closed", "absolute", loop_ambiguous,
            "p2-loop-ambiguous", label.c_str(), "selected", "T,P,x", "no selected density", "case");
        CHECK(loop_ambiguous);
        if (item.temperature != 308.8) continue;
        for (std::size_t side = 0; side < 2; ++side) {
            const auto batch = ethane.pressure_actions(item.temperature, item.pressure, composition,
                branches[side].first, {}, {epcsaft::PressureActionDirection{0.0, 1.0, {0.0}, {}}});
            const double reference = 1.0 / references[side == 0 ? 0 : 2].second;
            const bool ok = batch.density_status == epcsaft::PropertyActionStatus::Available
                && batch.density_action.size() == 1 && close(batch.density_action[0], reference, 1e-6);
            recording_support::compare((std::string("loop-action/") + (side == 0 ? "vapor" : "liquid")).c_str(),
                "density_action", "mol/(m3 Pa)", batch.density_action.empty() ? 0.0 : batch.density_action[0],
                reference, loop_inputs, "Identity", "reciprocal-slope",
                "d rho/dP on the selected branch is 1/P_rho at the independent root",
                "reference_or_one", 0.0, 1e-6, ok, "p2-loop-action", side == 0 ? "vapor" : "liquid",
                "d rho/dP", "T,P,x,branch", "rho_P=1/P_rho", "branch");
            CHECK(ok);
        }
    }
    // Independently locate the vapor-side fold at 308 K from its analytic fixed-density slope.
    // A pressure indistinguishable from this maximum cannot certify two nonsingular branches.
    double fold_lower = 5000.0, fold_upper = 6400.0;
    CHECK(ethane.state(308.0, fold_lower, composition).pressure_density_derivative > 0.0);
    CHECK(ethane.state(308.0, fold_upper, composition).pressure_density_derivative < 0.0);
    for (std::size_t i = 0; i < 80; ++i) {
        const double middle = 0.5 * (fold_lower + fold_upper);
        if (ethane.state(308.0, middle, composition).pressure_density_derivative > 0.0)
            fold_lower = middle;
        else fold_upper = middle;
    }
    const auto fold = ethane.state(308.0, 0.5 * (fold_lower + fold_upper), composition);
    CHECK(std::abs(fold.pressure_density_derivative) < 1e-8);
    bool unresolved = false;
    try {
        (void)epcsaft::solve_density(ethane, 308.0, fold.pressure, composition,
                                    epcsaft::DensityBranch::Vapor, std::nullopt);
    } catch (const std::runtime_error& error) {
        unresolved = std::string(error.what()).find("unresolved at spinodal") != std::string::npos;
    }
    recording_support::compare_custom("fold-resolution/vapor", "unresolved", "bool",
        unresolved ? 1.0 : 0.0, 1.0, recording_support::fields({
            recording_support::field("T", 308.0), recording_support::field("P_star", fold.pressure),
            recording_support::field("fold_density", fold.molar_density),
            recording_support::field("pressure_density_derivative", fold.pressure_density_derivative)}),
        "Availability", "independent-slope-bisection", "a mechanical fold is not a resolved stable branch",
        "absolute", unresolved, "p2-fold-resolution", "vapor", "unresolved", "T,P,x", "P_rho=0", "case");
    CHECK(unresolved);
}

void cache_invariance() {
    epcsaft::HelmholtzModel model(properties_support::aqueous_electrolyte());
    recording_support::expect_group(
        recording_support::active_test_id(), "m9-value", "value/", 2,
        "cached value invariance against a fresh model", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m9-gradient", "gradient/", 10,
        "cached gradient invariance against a fresh model", "coordinate");
    recording_support::expect_group(
        recording_support::active_test_id(), "m9-hessian", "hessian/", 50,
        "cached Hessian invariance against a fresh model", "row,column");
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
        CHECK(value_ok);
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
            CHECK(gradient_ok);
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
            CHECK(hessian_ok);
        }
    }
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
    recording_support::expect_group(recording_support::active_test_id(),
        "c2-scalar-values", "c2-scalar-value/", scalar_families.size(),
        "Core C2 scalar coordinate value identity", "family");
    recording_support::expect_group(recording_support::active_test_id(),
        "c2-falsifier", "c2-falsifier", 1, "Core C2 Hessian action", "row");
    recording_support::expect_group(recording_support::active_test_id(),
        "c2-association-gradient", "c2-association-gradient", 1,
        "Core C2 association gradient", "row");
    recording_support::expect_group(recording_support::active_test_id(),
        "c2-scaling", "c2-scaling", 1, "Core C2 source scaling", "row");
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

using properties_support::p3_forms;
using properties_support::p3_record;
using properties_support::p3_binary_record;

std::string p3_inputs(const char* form, double temperature = 320.0) {
    return recording_support::fields({
        recording_support::field("form", form),
        recording_support::field("T", temperature),
        recording_support::field("rho", 1000.0),
        recording_support::field("P0", 1.0e5),
    });
}

void p3_compare(const std::string& id, const char* quantity, const char* unit,
                double model, double reference, const std::string& inputs,
                const char* group, const char* point, double tolerance = 2.0e-12) {
    const bool ok = close(model, reference, tolerance);
    recording_support::compare(
        id.c_str(), quantity, unit, model, reference, inputs, "Relationship",
        "independent-reference", "P3 property/reference relation", "reference_or_one",
        0.0, tolerance, ok, group, point, quantity, "T,rho,P0", "model=reference", "state"
    );
    CHECK(ok);
}

void properties_thermochemistry_p3() {
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const double temperature = 320.0;
    const double density = 1000.0;
    const double standard_entropy = 50.0 + 30.0 * std::log(temperature / 300.0);
    const double expected_h = 1000.0 + 30.0 * (temperature - 300.0);
    const double expected_u = expected_h - R * temperature;
    const double expected_s = standard_entropy
        - R * std::log(density * R * temperature / 1.0e5);
    const auto forms = p3_forms();

    recording_support::expect_group(recording_support::active_test_id(),
        "p3-form-values", "p3-form-value/", forms.size() * 5,
        "all ideal forms match independent values", "form,quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-form-derivatives", "p3-form-derivative/", forms.size() * 8,
        "ideal h and s temperature identities", "form,quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-form-derivative-refinement", "p3-form-derivative-refinement/", forms.size() * 2,
        "final ideal derivative estimates refine", "form,quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-total", "p3-total/", forms.size() * 6,
        "production Ideal/Total state identities", "form,quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-boundaries", "p3-boundary/", 6,
        "endpoint and invalid interval behavior", "case");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-endpoints", "p3-endpoint/", forms.size() * 10,
        "all ideal forms match both owned endpoints", "form,endpoint,quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-mixing", "p3-mixing/", 5,
        "multicomponent ideal mixing and anchors", "quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-availability", "p3-availability/", 5,
        "typed thermochemistry availability causes", "case");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-reference", "p3-reference/", 13,
        "fixed standard-state conversion and pressure invariance", "quantity");
    recording_support::expect_group(recording_support::active_test_id(),
        "p3-reaction", "p3-reaction/", 34,
        "reaction basis and K(T) relations", "basis,quantity");

    for (const auto& [name, correlation] : forms) {
        const auto record = p3_record(correlation);
        const auto ideal = epcsaft::evaluate_ideal_properties(record, temperature, density, {1.0});
        CHECK(ideal.value.has_value() && ideal.status == epcsaft::ThermochemistryStatus::Available);
        const auto inputs = p3_inputs(name);
        p3_compare(std::string("p3-form-value/") + name + "/h", "ideal_enthalpy", "J/mol",
                   ideal.value->enthalpy, expected_h, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/u", "ideal_internal_energy", "J/mol",
                   ideal.value->internal_energy, expected_u, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/s", "ideal_entropy", "J/(mol K)",
                   ideal.value->entropy, expected_s, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/cp", "ideal_heat_capacity_cp", "J/(mol K)",
                   ideal.value->isobaric_heat_capacity, 30.0, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/cv", "ideal_heat_capacity_cv", "J/(mol K)",
                   ideal.value->isochoric_heat_capacity, 30.0 - R, inputs, "p3-form-values", name);

        constexpr std::array<double, 4> steps = {0.1, 0.03, 0.01, 0.003};
        double refined_dh[2] = {}, refined_ds[2] = {};
        for (std::size_t step = 0; step < steps.size(); ++step) {
            const double delta = steps[step];
            const auto plus = epcsaft::evaluate_ideal_correlation(correlation, temperature + delta, 1.0e5);
            const auto minus = epcsaft::evaluate_ideal_correlation(correlation, temperature - delta, 1.0e5);
            const double dh_dt = (plus.h - minus.h) / (2.0 * delta);
            const double t_ds_dt = temperature * (plus.s - minus.s) / (2.0 * delta);
            const auto point = std::string(name) + "/" + recording_support::number(delta);
            const double derivative_tolerance = delta <= 0.01 ? 1.0e-8 : 1.0e-5;
            p3_compare(std::string("p3-form-derivative/") + point + "/dh", "dh_dT", "J/(mol K)",
                       dh_dt, ideal.value->isobaric_heat_capacity, inputs, "p3-form-derivatives", point.c_str(), derivative_tolerance);
            p3_compare(std::string("p3-form-derivative/") + point + "/ds", "T_ds_dT", "J/(mol K)",
                       t_ds_dt, ideal.value->isobaric_heat_capacity, inputs, "p3-form-derivatives", point.c_str(), derivative_tolerance);
            if (step == 2) { refined_dh[0] = dh_dt; refined_ds[0] = t_ds_dt; }
            if (step == 3) { refined_dh[1] = dh_dt; refined_ds[1] = t_ds_dt; }
        }
        p3_compare(std::string("p3-form-derivative-refinement/") + name + "/dh", "dh_dT_refinement", "J/(mol K)",
                   refined_dh[1], refined_dh[0], inputs, "p3-form-derivative-refinement", name, 1.0e-8);
        p3_compare(std::string("p3-form-derivative-refinement/") + name + "/ds", "T_ds_dT_refinement", "J/(mol K)",
                   refined_ds[1], refined_ds[0], inputs, "p3-form-derivative-refinement", name, 1.0e-8);

        epcsaft::HelmholtzModel model(properties_support::ethane(), {}, record);
        const auto state = model.state(temperature, density, {1.0});
        CHECK(state.thermochemistry_status == epcsaft::ThermochemistryStatus::Available);
        CHECK(state.ideal_enthalpy.has_value() && state.ideal_internal_energy.has_value()
              && state.ideal_entropy.has_value() && state.ideal_isobaric_heat_capacity.has_value()
              && state.ideal_isochoric_heat_capacity.has_value()
              && state.total_enthalpy.has_value() && state.total_internal_energy.has_value()
              && state.total_entropy.has_value() && state.total_isochoric_heat_capacity.has_value());
        p3_compare(std::string("p3-total/") + name + "/h", "total_h_decomposition", "J/mol",
                   *state.total_enthalpy, *state.ideal_enthalpy + state.residual_enthalpy,
                   inputs, "p3-total", name);
        p3_compare(std::string("p3-total/") + name + "/u", "total_u_decomposition", "J/mol",
                   *state.total_internal_energy, *state.ideal_internal_energy + state.residual_internal_energy,
                   inputs, "p3-total", name);
        p3_compare(std::string("p3-total/") + name + "/s", "total_s_decomposition", "J/(mol K)",
                   *state.total_entropy, *state.ideal_entropy + state.residual_entropy,
                   inputs, "p3-total", name);
        p3_compare(std::string("p3-total/") + name + "/cp", "total_cp_decomposition", "J/(mol K)",
                   *state.total_isobaric_heat_capacity,
                   *state.ideal_isobaric_heat_capacity + *state.residual_isobaric_heat_capacity,
                   inputs, "p3-total", name);
        p3_compare(std::string("p3-total/") + name + "/cv", "total_cv_decomposition", "J/(mol K)",
                   *state.total_isochoric_heat_capacity,
                   *state.ideal_isochoric_heat_capacity + state.residual_isochoric_heat_capacity,
                   inputs, "p3-total", name);
        p3_compare(std::string("p3-total/") + name + "/pv", "total_h_u_pv", "J/mol",
                   *state.total_enthalpy,
                   *state.total_internal_energy + state.pressure / state.molar_density,
                   inputs, "p3-total", name);
    }

    for (const auto& [name, correlation] : forms) {
        for (double endpoint_temperature : {250.0, 350.0}) {
            const auto endpoint_state = epcsaft::evaluate_ideal_properties(
                p3_record(correlation), endpoint_temperature, density, {1.0}
            );
            CHECK(endpoint_state.value.has_value());
            const double h = 1000.0 + 30.0 * (endpoint_temperature - 300.0);
            const double u = h - R * endpoint_temperature;
            const double s = 50.0 + 30.0 * std::log(endpoint_temperature / 300.0)
                - R * std::log(density * R * endpoint_temperature / 1.0e5);
            const auto point = std::string(name) + "/" + recording_support::number(endpoint_temperature);
            const auto endpoint_inputs = p3_inputs(name, endpoint_temperature);
            p3_compare(std::string("p3-endpoint/") + point + "/h", "ideal_enthalpy", "J/mol",
                       endpoint_state.value->enthalpy, h, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/u", "ideal_internal_energy", "J/mol",
                       endpoint_state.value->internal_energy, u, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/s", "ideal_entropy", "J/(mol K)",
                       endpoint_state.value->entropy, s, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/cp", "ideal_heat_capacity_cp", "J/(mol K)",
                       endpoint_state.value->isobaric_heat_capacity, 30.0, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/cv", "ideal_heat_capacity_cv", "J/(mol K)",
                       endpoint_state.value->isochoric_heat_capacity, 30.0 - R, endpoint_inputs, "p3-endpoints", point.c_str());
        }
    }

    const auto binary = epcsaft::evaluate_ideal_properties(
        p3_binary_record(), temperature, density, {0.25, 0.75}
    );
    CHECK(binary.value.has_value());
    const double binary_h = 0.25 * expected_h + 0.75 * (2000.0 + 40.0 * (temperature - 300.0));
    const double binary_u = binary_h - R * temperature;
    const double binary_cp = 0.25 * 30.0 + 0.75 * 40.0;
    const double binary_s = 0.25 * (standard_entropy - R * std::log(
        0.25 * density * R * temperature / 1.0e5))
        + 0.75 * (70.0 + 40.0 * std::log(temperature / 300.0) - R * std::log(
            0.75 * density * R * temperature / 1.0e5));
    p3_compare("p3-mixing/h", "ideal_enthalpy", "J/mol", binary.value->enthalpy, binary_h,
               "x=0.25,0.75,T=320", "p3-mixing", "h");
    p3_compare("p3-mixing/u", "ideal_internal_energy", "J/mol", binary.value->internal_energy, binary_u,
               "x=0.25,0.75,T=320", "p3-mixing", "u");
    p3_compare("p3-mixing/s", "ideal_entropy", "J/(mol K)", binary.value->entropy, binary_s,
               "x=0.25,0.75,T=320", "p3-mixing", "s");
    p3_compare("p3-mixing/cp", "ideal_heat_capacity_cp", "J/(mol K)", binary.value->isobaric_heat_capacity,
               binary_cp, "x=0.25,0.75,T=320", "p3-mixing", "cp");
    p3_compare("p3-mixing/cv", "ideal_heat_capacity_cv", "J/(mol K)", binary.value->isochoric_heat_capacity,
               binary_cp - R, "x=0.25,0.75,T=320", "p3-mixing", "cv");

    const auto endpoint_record = p3_record(forms.front().second);
    const auto endpoint = epcsaft::evaluate_ideal_thermal(
        endpoint_record.component_intervals.front(), 250.0, endpoint_record.standard_pressure
    );
    const bool endpoint_value = std::isfinite(endpoint.value.h) && !endpoint.derivatives_available;
    recording_support::compare_custom("p3-boundary/endpoint", "endpoint_value_available", "bool",
        endpoint_value ? 1.0 : 0.0, 1.0, "T=250", "Availability", "endpoint-owner",
        "owned endpoint has value and no cross-interval derivative", "absolute", endpoint_value,
        "p3-boundaries", "endpoint");
    const auto outside = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(400.0, density, {1.0});
    const bool outside_partial = !outside.ideal_enthalpy.has_value()
        && std::isfinite(outside.residual_enthalpy);
    recording_support::compare_custom("p3-boundary/outside", "residual_only", "bool",
        outside_partial ? 1.0 : 0.0, 1.0, "T=400", "Availability", "interval-owner",
        "outside interval preserves residual-only state", "absolute", outside_partial,
        "p3-boundaries", "outside");
    bool gap_failed = false, overlap_failed = false, pressure_failed = false;
    try {
        auto gap = endpoint_record;
        gap.component_intervals.front().push_back(
            epcsaft::IdealInterval{351.0, 360.0, true, true, forms.front().second}
        );
        epcsaft::validate_thermochemistry(gap, 1);
    } catch (const std::invalid_argument&) { gap_failed = true; }
    try {
        auto overlap = endpoint_record;
        overlap.component_intervals.front().push_back(
            epcsaft::IdealInterval{250.0, 350.0, true, true, forms.front().second}
        );
        epcsaft::validate_thermochemistry(overlap, 1);
    } catch (const std::invalid_argument&) { overlap_failed = true; }
    try {
        auto invalid = endpoint_record;
        invalid.component_intervals.front().front().correlation.source_pressure = 0.0;
        epcsaft::validate_thermochemistry(invalid, 1);
    } catch (const std::invalid_argument&) { pressure_failed = true; }
    recording_support::compare_custom("p3-boundary/gap", "gap_rejected", "bool",
        gap_failed ? 1.0 : 0.0, 1.0, "351-360 after 250-350", "Availability", "interval-owner",
        "gapped intervals are rejected", "absolute", gap_failed, "p3-boundaries", "gap");
    recording_support::compare_custom("p3-boundary/overlap", "overlap_rejected", "bool",
        overlap_failed ? 1.0 : 0.0, 1.0, "250-350 plus 250-350", "Availability", "interval-owner",
        "overlapping intervals are rejected", "absolute", overlap_failed, "p3-boundaries", "overlap");
    recording_support::compare_custom("p3-boundary/pressure", "invalid_pressure_rejected", "bool",
        pressure_failed ? 1.0 : 0.0, 1.0, "Psrc=0", "Availability", "input-validation",
        "invalid source pressure is rejected", "absolute", pressure_failed, "p3-boundaries", "pressure");
    bool missing_component_failed = false;
    try {
        auto missing = endpoint_record;
        missing.component_intervals.clear();
        epcsaft::validate_thermochemistry(missing, 1);
    } catch (const std::invalid_argument&) { missing_component_failed = true; }
    recording_support::compare_custom("p3-boundary/missing-component", "missing_component_rejected", "bool",
        missing_component_failed ? 1.0 : 0.0, 1.0, "component thermal record absent", "Availability",
        "input-validation", "missing component thermochemistry is unavailable", "absolute",
        missing_component_failed, "p3-boundaries", "missing-component");

    const auto missing_record_state = epcsaft::HelmholtzModel(properties_support::ethane()).state(
        temperature, density, {1.0}
    );
    recording_support::compare_custom("p3-availability/missing-record", "thermochemistry_status", "enum",
        static_cast<double>(missing_record_state.thermochemistry_status), 0.0, "no record",
        "Availability", "status-owner", "missing record has a typed cause", "absolute",
        missing_record_state.thermochemistry_status == epcsaft::ThermochemistryStatus::MissingRecord,
        "p3-availability", "missing-record");
    recording_support::compare_custom("p3-availability/outside", "thermochemistry_status", "enum",
        static_cast<double>(outside.thermochemistry_status), 2.0, "T=400 outside [250,350]",
        "Availability", "status-owner", "outside interval has a typed cause", "absolute",
        outside.thermochemistry_status == epcsaft::ThermochemistryStatus::OutsideInterval,
        "p3-availability", "outside");
    auto nonfinite_record = endpoint_record;
    std::get<epcsaft::IdealPolynomial>(
        nonfinite_record.component_intervals.front().front().correlation.form
    ).coefficients[0] = 1.0e308;
    const auto nonfinite = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, nonfinite_record
    ).state(temperature, density, {1.0});
    recording_support::compare_custom("p3-availability/nonfinite", "thermochemistry_status", "enum",
        static_cast<double>(nonfinite.thermochemistry_status), 3.0, "finite record with nonfinite evaluation",
        "Availability", "status-owner", "nonfinite evaluation has a typed cause", "absolute",
        nonfinite.thermochemistry_status == epcsaft::ThermochemistryStatus::Nonfinite,
        "p3-availability", "nonfinite");
    const auto negative_d = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(280.0, 4000.0, {1.0});
    const bool negative_d_partial = negative_d.thermochemistry_status == epcsaft::ThermochemistryStatus::Available
        && !negative_d.total_isobaric_heat_capacity.has_value()
        && std::isfinite(negative_d.residual_enthalpy);
    recording_support::compare_custom("p3-availability/negative-D", "partial_state", "bool",
        negative_d_partial ? 1.0 : 0.0, 1.0, "T=280,rho=4000,D<0",
        "Availability", "partial-state", "negative-D preserves other fields and total cp unavailability",
        "absolute", negative_d_partial, "p3-availability", "negative-D");
    const auto available_state = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(temperature, density, {1.0});
    recording_support::compare_custom("p3-availability/available", "thermochemistry_status", "enum",
        static_cast<double>(available_state.thermochemistry_status), 1.0, "T=320 inside [250,350]",
        "Availability", "status-owner", "available record has a typed cause", "absolute",
        available_state.thermochemistry_status == epcsaft::ThermochemistryStatus::Available,
        "p3-availability", "available");

    auto high_pressure = endpoint_record;
    high_pressure.standard_pressure = 2.0e5;
    const auto low_state = epcsaft::evaluate_ideal_properties(endpoint_record, temperature, density, {1.0});
    const auto high_ideal = epcsaft::evaluate_ideal_properties(high_pressure, temperature, density, {1.0});
    CHECK(low_state.value.has_value() && high_ideal.value.has_value());
    p3_compare("p3-reference/h", "h_pressure_invariance", "J/mol",
               high_ideal.value->enthalpy, low_state.value->enthalpy, "P0=1e5 to 2e5",
               "p3-reference", "h");
    p3_compare("p3-reference/u", "u_pressure_invariance", "J/mol",
               high_ideal.value->internal_energy, low_state.value->internal_energy, "P0=1e5 to 2e5",
               "p3-reference", "u");
    p3_compare("p3-reference/s", "s_pressure_invariance", "J/(mol K)",
               high_ideal.value->entropy, low_state.value->entropy, "P0=1e5 to 2e5",
               "p3-reference", "s");
    p3_compare("p3-reference/cp", "cp_pressure_invariance", "J/(mol K)",
               high_ideal.value->isobaric_heat_capacity, low_state.value->isobaric_heat_capacity,
               "P0=1e5 to 2e5", "p3-reference", "cp");
    p3_compare("p3-reference/cv", "cv_pressure_invariance", "J/(mol K)",
               high_ideal.value->isochoric_heat_capacity, low_state.value->isochoric_heat_capacity,
               "P0=1e5 to 2e5", "p3-reference", "cv");
    const auto low_total = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(temperature, density, {1.0});
    const auto high_total = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, high_pressure
    ).state(temperature, density, {1.0});
    for (const auto& item : std::array<std::tuple<const char*, double, double, const char*>, 5>{
        std::tuple{"h", *low_total.total_enthalpy, *high_total.total_enthalpy, "J/mol"},
        {"u", *low_total.total_internal_energy, *high_total.total_internal_energy, "J/mol"},
        {"s", *low_total.total_entropy, *high_total.total_entropy, "J/(mol K)"},
        {"cp", *low_total.total_isobaric_heat_capacity, *high_total.total_isobaric_heat_capacity, "J/(mol K)"},
        {"cv", *low_total.total_isochoric_heat_capacity, *high_total.total_isochoric_heat_capacity, "J/(mol K)"}
    }) {
        const auto [quantity, low, high, unit] = item;
        const auto point = std::string("total-") + quantity;
        p3_compare(std::string("p3-reference/") + point, "total_pressure_invariance", unit,
                   high, low, "State P0=1e5 to 2e5", "p3-reference", point.c_str());
    }
    const auto converted = epcsaft::pressure_to_density_standard(
        1234.0, 5678.0, 21.0, 300.0, 1.0e5, 1000.0
    );
    const double shift = R * 300.0 * std::log(R * 300.0 * 1000.0 / 1.0e5);
    p3_compare("p3-reference/g", "g_density_conversion", "J/mol",
               converted.gibbs, 1234.0 + shift, "P0=1e5,rho0=1000",
               "p3-reference", "g");
    p3_compare("p3-reference/h0", "h_density_conversion", "J/mol",
               converted.enthalpy, 5678.0 - R * 300.0, "P0=1e5,rho0=1000",
               "p3-reference", "h0");
    p3_compare("p3-reference/s0", "s_density_conversion", "J/(mol K)",
               converted.entropy, 21.0 - R - shift / 300.0, "P0=1e5,rho0=1000",
               "p3-reference", "s0");

    const epcsaft::ReactionLogPolynomial correlation{
        2.0, -1200.0, 0.4, -0.003, 300.0, 293.15, 393.15, "P0"
    };
    const std::vector<double> stoichiometry = {1.0, -1.0};
    const auto direct = epcsaft::direct_activities({0.8, 1.2}, "P0");
    const auto x_activity = epcsaft::mole_fraction_activities({2.0, 3.0}, {0.4, 0.6}, "P0");
    const auto c_activity = epcsaft::molarity_activities({2.0, 3.0}, 500.0, {0.4, 0.6}, 100.0, "P0");
    const auto m_activity = epcsaft::molality_activities({2.0, 3.0}, {1.0, 2.0}, 1.0, 1.0, "P0");
    p3_compare("p3-reaction/direct", "ln_q", "1",
               epcsaft::reaction_ln_quotient(stoichiometry, direct, {}, correlation),
               std::log(0.8 / 1.2), "basis=activity", "p3-reaction", "direct");
    p3_compare("p3-reaction/x", "ln_q", "1",
               epcsaft::reaction_ln_quotient(stoichiometry, x_activity, {}, correlation),
               std::log(0.8 / 1.8), "basis=mole-fraction", "p3-reaction", "x");
    p3_compare("p3-reaction/c", "ln_q", "1",
               epcsaft::reaction_ln_quotient(stoichiometry, c_activity, {}, correlation),
               std::log(4.0 / 9.0), "basis=molarity,rho=500,c0=100 mol/m3", "p3-reaction", "c");
    p3_compare("p3-reaction/c-delta-nu", "ln_q", "1",
               epcsaft::reaction_ln_quotient({2.0, -1.0}, c_activity, {}, correlation),
               std::log(16.0 / 9.0), "basis=molarity,rho=500,Delta_nu=1", "p3-reaction", "c-delta-nu");
    p3_compare("p3-reaction/m", "ln_q", "1",
               epcsaft::reaction_ln_quotient(stoichiometry, m_activity, {}, correlation),
               std::log(2.0 / 6.0), "basis=molality,m0=1 mol/kg", "p3-reaction", "m");
    // Standard-scale identity (#84): for a solute in 1 kg of water, ln a^(m) - ln a^(x)
    // = ln((n_s + n_w)/(n_w M_w m0)) -> -ln(M_w m0) = 4.0165 at infinite dilution, so a reaction
    // with one net solute product (the R4/R5 pattern: solute + water -> solute + solute, water on
    // its own mole-fraction scale) shifts its quotient by exactly that amount, a factor 55.5 in K.
    {
        // Species order (MEACOO-, H2O, HCO3-, MEA) with distinct solute amounts; the solvent mass comes
        // from the water amount. The separate infinite-dilution value check below pins M_w m0.
        const double water_molar_mass = 0.01801528, n_water = 40.0, standard_molality = 1.0;
        const std::vector<double> amounts = {2.0e-6, n_water, 3.0e-6, 5.0e-6};
        const double total = amounts[0] + amounts[1] + amounts[2] + amounts[3];
        const auto x_scale = epcsaft::mole_fraction_activities({1.0, 1.0, 1.0, 1.0},
            {amounts[0] / total, amounts[1] / total, amounts[2] / total, amounts[3] / total}, "P0");
        const auto solutes = epcsaft::molality_activities({1.0, 1.0, 1.0, 1.0}, amounts,
            n_water * water_molar_mass, standard_molality, "P0");
        // Water stays on its own mole-fraction scale in the molality-basis set, as the source rows state.
        const auto m_scale = epcsaft::direct_activities(
            {solutes.values[0], x_scale.values[1], solutes.values[2], solutes.values[3]}, "P0");
        const std::vector<double> r4_pattern = {-1.0, -1.0, 1.0, 1.0};
        const double expected = std::log(total / (n_water * water_molar_mass * standard_molality));
        p3_compare("p3-reaction/molality-scale-solute", "ln_a_m_minus_ln_a_x", "1",
                   std::log(m_scale.values[2]) - std::log(x_scale.values[2]), expected,
                   "n_w=40 mol water, n_s=3e-6 mol", "p3-reaction", "molality-scale-solute");
        p3_compare("p3-reaction/molality-scale-r4", "ln_q_m_minus_ln_q_x", "1",
                   epcsaft::reaction_ln_quotient(r4_pattern, m_scale, {}, correlation)
                       - epcsaft::reaction_ln_quotient(r4_pattern, x_scale, {}, correlation),
                   expected, "R4 pattern, sum nu_solute = +1, quotient owner", "p3-reaction", "molality-scale-r4");
        CHECK(std::abs(expected - 4.016535) < 1.0e-5);
    }
    const auto neutral = epcsaft::direct_activities({1.0, 1.0}, "P0");
    p3_compare("p3-reaction/neutral-ionic", "ln_q", "1",
               epcsaft::reaction_ln_quotient(stoichiometry, neutral, {1, 1}, correlation),
               0.0, "charge-neutral ionic contraction", "p3-reaction", "neutral-ionic");
    bool gauge_failed = false;
    try {
        epcsaft::reaction_ln_quotient({1.0}, epcsaft::direct_activities({1.0}, "P0"), {1}, correlation);
    } catch (const std::invalid_argument& error) {
        gauge_failed = std::string(error.what()) == "ionic-gauge-missing";
    }
    recording_support::compare_custom("p3-reaction/gauge", "ionic_gauge_required", "bool",
        gauge_failed ? 1.0 : 0.0, 1.0, "charged nonneutral reaction", "Availability",
        "ionic-gauge", "nonneutral ionic contraction remains unavailable without a real gauge",
        "absolute", gauge_failed, "p3-reaction", "gauge");
    bool mismatch_failed = false;
    try {
        epcsaft::reaction_ln_quotient(stoichiometry, epcsaft::direct_activities({0.8, 1.2}, "rho0"), {}, correlation);
    } catch (const std::invalid_argument& error) {
        mismatch_failed = std::string(error.what()) == "reference-basis-mismatch";
    }
    recording_support::compare_custom("p3-reaction/mismatch", "basis_identity_required", "bool",
        mismatch_failed ? 1.0 : 0.0, 1.0, "activity id=rho0,K id=P0", "Availability",
        "standard-state-identity", "quotient and K standard-state identities must match",
        "absolute", mismatch_failed, "p3-reaction", "mismatch");

    const auto reaction = epcsaft::evaluate_reaction_properties(correlation, 350.0);
    const double expected_ln_k = 2.0 - 1200.0 / 350.0
        + 0.4 * std::log(350.0 / 300.0) - 0.003 * 350.0;
    const double expected_derivative = 1200.0 / (350.0 * 350.0) + 0.4 / 350.0 - 0.003;
    p3_compare("p3-reaction/ln-k", "ln_k", "1", reaction.ln_k, expected_ln_k,
               "T=350", "p3-reaction", "ln-k");
    CHECK(reaction.k.has_value());
    p3_compare("p3-reaction/k", "K", "1", *reaction.k, std::exp(expected_ln_k),
               "T=350", "p3-reaction", "k");
    p3_compare("p3-reaction/dk", "dln_k_dT", "1/K", reaction.derivative, expected_derivative,
               "T=350", "p3-reaction", "dk");
    p3_compare("p3-reaction/h", "reaction_enthalpy", "J/mol", reaction.enthalpy,
               R * 350.0 * 350.0 * expected_derivative, "T=350", "p3-reaction", "h");

    const auto log_correlation = epcsaft::normalize_log10_reciprocal_correlation(
        2677.91, 0.3869, 0.0004277, 293.15, 393.15, "P0", 300.0
    );
    const auto log_reaction = epcsaft::evaluate_reaction_properties(log_correlation, 300.0);
    p3_compare("p3-reaction/log10-ln", "ln_k", "1", log_reaction.ln_k,
               -std::log(10.0) * (2677.91 / 300.0 + 0.3869 + 0.0004277 * 300.0),
               "T=300", "p3-reaction", "log10-ln");
    CHECK(log_reaction.k.has_value());
    p3_compare("p3-reaction/log10-k", "K", "1", *log_reaction.k,
               std::exp(log_reaction.ln_k), "T=300", "p3-reaction", "log10-k");
    p3_compare("p3-reaction/log10-dk", "dln_k_dT", "1/K", log_reaction.derivative,
               std::log(10.0) * (2677.91 / (300.0 * 300.0) - 0.0004277),
               "T=300", "p3-reaction", "log10-dk");
    const auto legacy = epcsaft::normalize_legacy_reaction_correlation(
        2.0, -1200.0, 0.4, -0.003, 300.0, 293.15, 393.15, "P0"
    );
    p3_compare("p3-reaction/legacy", "ln_k", "1",
               epcsaft::evaluate_reaction_properties(legacy, 350.0).ln_k,
               2.0 - 1200.0 / 350.0 + 0.4 * std::log(350.0) - 0.003 * 350.0,
               "legacy ln(T/1 K)", "p3-reaction", "legacy");
    bool endpoints = false, finite_k = true, outside_k = false;
    try {
        const auto lower = epcsaft::evaluate_reaction_properties(correlation, correlation.temperature_min);
        const auto upper = epcsaft::evaluate_reaction_properties(correlation, correlation.temperature_max);
        endpoints = true;
        finite_k = lower.k.has_value() && upper.k.has_value();
    } catch (const std::exception&) {}
    try {
        epcsaft::evaluate_reaction_properties(correlation, 280.0);
    } catch (const std::invalid_argument&) { outside_k = true; }
    recording_support::compare_custom("p3-reaction/endpoints", "K_endpoint_available", "bool",
        endpoints ? 1.0 : 0.0, 1.0, "Tmin/Tmax inclusive", "Availability", "correlation-interval",
        "correlation endpoints are owned and available", "absolute", endpoints, "p3-reaction", "endpoints");
    recording_support::compare_custom("p3-reaction/finite-k", "K_finite", "bool",
        finite_k ? 1.0 : 0.0, 1.0, "Tmin/Tmax", "Availability", "finite-K",
        "finite ln K produces finite K", "absolute", finite_k, "p3-reaction", "finite-k");
    const auto overflow = epcsaft::evaluate_reaction_properties(
        {1000.0, 0.0, 0.0, 0.0, 300.0, 293.15, 393.15, "P0"}, 300.0
    );
    const auto underflow = epcsaft::evaluate_reaction_properties(
        {-1000.0, 0.0, 0.0, 0.0, 300.0, 293.15, 393.15, "P0"}, 300.0
    );
    const bool k_separate = !overflow.k.has_value() && !underflow.k.has_value()
        && std::isfinite(overflow.ln_k) && std::isfinite(underflow.ln_k)
        && std::isfinite(overflow.derivative) && std::isfinite(underflow.enthalpy);
    recording_support::compare_custom("p3-reaction/k-unavailable", "K_availability_separate", "bool",
        k_separate ? 1.0 : 0.0, 1.0, "ln K=+/-1000", "Availability", "finite-K",
        "finite ln K, derivative and enthalpy survive K overflow/underflow",
        "absolute", k_separate, "p3-reaction", "k-unavailable");
    recording_support::compare_custom("p3-reaction/outside", "K_outside_rejected", "bool",
        outside_k ? 1.0 : 0.0, 1.0, "T=280 outside [293.15,393.15]", "Availability",
        "correlation-interval", "out-of-domain K is unavailable", "absolute", outside_k,
        "p3-reaction", "outside");

    const std::vector<double> density_stoichiometry = {2.0, -1.0};
    const auto density_reaction = epcsaft::reaction_pressure_to_density(
        reaction, density_stoichiometry, "P0=1e5;rho0=1000", 350.0, 1.0e5, 1000.0
    );
    const double density_shift = std::log(R * 350.0 * 1000.0 / 1.0e5);
    p3_compare("p3-reaction/density-ln", "ln_k_density", "1", density_reaction.ln_k,
               reaction.ln_k - density_shift, "P0=1e5,rho0=1000,Delta_nu=1",
               "p3-reaction", "density-ln");
    CHECK(density_reaction.k.has_value());
    p3_compare("p3-reaction/density-k", "K_density", "1", *density_reaction.k,
               std::exp(reaction.ln_k - density_shift), "P0=1e5,rho0=1000,Delta_nu=1",
               "p3-reaction", "density-k");
    p3_compare("p3-reaction/density-dk", "dln_k_density_dT", "1/K", density_reaction.derivative,
               reaction.derivative - 1.0 / 350.0, "P0=1e5,rho0=1000,Delta_nu=1",
               "p3-reaction", "density-dk");
    p3_compare("p3-reaction/density-h", "reaction_enthalpy_density", "J/mol", density_reaction.enthalpy,
               reaction.enthalpy - R * 350.0, "P0=1e5,rho0=1000,Delta_nu=1",
               "p3-reaction", "density-h");
    const bool density_identity = density_reaction.standard_state_id == "P0=1e5;rho0=1000";
    recording_support::compare_custom("p3-reaction/density-identity", "standard_state_identity", "bool",
        density_identity ? 1.0 : 0.0, 1.0, "explicit P0/rho0 identity", "Availability",
        "standard-state-identity", "density conversion retains the caller's exact target identity",
        "absolute", density_identity, "p3-reaction", "density-identity");

    // Old-minus-new species standards at T0 = 350 K: Delta g = {500, 2000} and Delta h = {500, 2400}
    // J/mol give Delta s = (Delta h - Delta g)/T0; no heat-capacity difference is the constant shift.
    const auto converted_reaction = epcsaft::evaluate_reaction_properties(epcsaft::reaction_standard_state_shift(
        correlation, stoichiometry, {500.0, 2400.0}, {0.0, 400.0 / 350.0}, {0.0, 0.0}, {0.0, 0.0}, 350.0, "Q"), 350.0);
    const double gibbs_shift = -1500.0;
    p3_compare("p3-reaction/std-ln", "ln_k_standard", "1", converted_reaction.ln_k,
               reaction.ln_k + gibbs_shift / (R * 350.0), "P0 to Q",
               "p3-reaction", "std-ln");
    CHECK(converted_reaction.k.has_value());
    p3_compare("p3-reaction/std-k", "K_standard", "1", *converted_reaction.k,
               std::exp(reaction.ln_k + gibbs_shift / (R * 350.0)), "P0 to Q",
               "p3-reaction", "std-k");
    const double enthalpy_shift = -1900.0;
    p3_compare("p3-reaction/std-dk", "dln_k_standard_dT", "1/K", converted_reaction.derivative,
               reaction.derivative - enthalpy_shift / (R * 350.0 * 350.0), "P0 to Q",
               "p3-reaction", "std-dk");
    p3_compare("p3-reaction/std-h", "reaction_enthalpy_standard", "J/mol", converted_reaction.enthalpy,
               reaction.enthalpy - enthalpy_shift, "P0 to Q",
               "p3-reaction", "std-h");

    const double tref = correlation.reference_temperature;
    const double species_h_ref = R * (-correlation.b + correlation.c * tref + correlation.d * tref * tref);
    const auto species_h = [=](double t) {
        return species_h_ref + R * (correlation.c * (t - tref) + correlation.d * (t * t - tref * tref));
    };
    p3_compare("p3-reaction/species-vanthoff", "species_delta_h", "J/mol",
               species_h(350.0), reaction.enthalpy, "supplied h_i(Tref),cp_i(T)",
               "p3-reaction", "species-vanthoff");
    constexpr int integration_steps = 1000;
    const double integration_step = (350.0 - tref) / integration_steps;
    double species_integral = species_h(tref) / (R * tref * tref) + species_h(350.0) / (R * 350.0 * 350.0);
    for (int i = 1; i < integration_steps; ++i) {
        const double t = tref + i * integration_step;
        species_integral += (i % 2 == 0 ? 2.0 : 4.0) * species_h(t) / (R * t * t);
    }
    species_integral *= integration_step / 3.0;
    p3_compare("p3-reaction/species-integral", "ln_k_species_integral", "1",
               species_integral, reaction.ln_k
                   - epcsaft::evaluate_reaction_properties(correlation, tref).ln_k,
               "supplied h_i(Tref),cp_i(T),Simpson-1000", "p3-reaction", "species-integral", 1.0e-8);
    bool missing_standard_state = false;
    try {
        auto missing = correlation;
        missing.standard_state_id.clear();
        epcsaft::evaluate_reaction_properties(missing, 350.0);
    } catch (const std::invalid_argument&) { missing_standard_state = true; }
    recording_support::compare_custom("p3-reaction/missing-standard-state", "standard_state_required", "bool",
        missing_standard_state ? 1.0 : 0.0, 1.0, "correlation identity absent and no gauge substitute",
        "Availability", "input-validation", "missing source standard state remains unavailable",
        "absolute", missing_standard_state, "p3-reaction", "missing-standard-state");
}

double p4_state_value(const epcsaft::State& state, epcsaft::PropertyObservable observable) {
    const auto index = static_cast<std::size_t>(observable) % 5;
    const bool ideal = static_cast<std::size_t>(observable) < 5;
    const bool total = static_cast<std::size_t>(observable) >= 10;
    const auto residual = std::array<double, 5>{
        state.residual_enthalpy, state.residual_internal_energy, state.residual_entropy,
        state.residual_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.residual_isochoric_heat_capacity};
    const auto ideal_values = std::array<double, 5>{
        state.ideal_enthalpy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_internal_energy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_entropy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_isochoric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN())};
    const auto total_values = std::array<double, 5>{
        state.total_enthalpy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_internal_energy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_entropy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_isochoric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN())};
    return (ideal ? ideal_values : total ? total_values : residual)[index];
}

epcsaft::Parameters p4_shifted_parameters(
    const epcsaft::Parameters& base,
    const std::vector<epcsaft::ActiveParameter>& active,
    const std::vector<double>& shift
) {
    auto result = base;
    for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
        const auto& item = active[coordinate];
        auto apply = [&](std::vector<double>& values) {
            for (const auto index : item.indices) values[index] += shift[coordinate];
        };
        switch (item.family) {
        case epcsaft::ParameterFamily::AssociationEnergyOverK:
            apply(result.association_energy_over_k); break;
        case epcsaft::ParameterFamily::AssociationVolume:
            apply(result.association_volume); break;
        case epcsaft::ParameterFamily::Kij:
            apply(result.k_ij); break;
        case epcsaft::ParameterFamily::BornDiameter:
            apply(result.born_diameter_angstrom); break;
        case epcsaft::ParameterFamily::DebyeHuckelDiameter:
            apply(result.debye_huckel_diameter_angstrom); break;
        case epcsaft::ParameterFamily::PackingDiameter:
            apply(result.packing_diameter_angstrom); break;
        case epcsaft::ParameterFamily::IonSuppressionCoefficient:
            apply(result.ion_suppression_coefficient); break;
        default: throw std::invalid_argument("unsupported P4 reference parameter");
        }
    }
    return result;
}

std::vector<double> p4_take_signature(
    epcsaft::core::Capture& capture,
    const std::optional<epcsaft::ThermochemistryRecord>& thermochemistry,
    double temperature
) {
    std::vector<double> result;
    for (const auto& event : capture.trace().events)
        if (std::string(event.name) == "predicate_signature") result = event.values;
    if (thermochemistry.has_value())
        for (const auto& intervals : thermochemistry->component_intervals)
            for (const auto& interval : intervals)
                result.push_back(epcsaft::ideal_interval_contains(interval, temperature));
    capture.take();
    return result;
}

double p4_reference(
    const epcsaft::Parameters& parameters,
    const std::vector<epcsaft::ActiveParameter>& active,
    const std::optional<epcsaft::ThermochemistryRecord>& thermochemistry,
    double temperature, double molar_density, const std::vector<double>& composition,
    const std::vector<epcsaft::DensityActionDirection>& directions,
    epcsaft::PropertyObservable observable, double step,
    std::vector<double>* admission = nullptr
) {
    epcsaft::core::Capture capture;
    const epcsaft::HelmholtzModel central_model(parameters, {}, thermochemistry);
    (void)central_model.state(temperature, molar_density, composition);
    const auto central_signature = p4_take_signature(capture, thermochemistry, temperature);
    const std::size_t combinations = std::size_t{1} << directions.size();
    double result = 0.0;
    bool support_unchanged = true, predicate_unchanged = true;
    for (std::size_t mask = 0; mask < combinations; ++mask) {
        double signed_temperature = temperature;
        double signed_density = molar_density;
        std::vector<double> signed_composition = composition;
        std::vector<double> parameter_shift(active.size(), 0.0);
        int sign_product = 1;
        for (std::size_t index = 0; index < directions.size(); ++index) {
            const double sign = (mask & (std::size_t{1} << index)) != 0 ? 1.0 : -1.0;
            sign_product *= sign > 0.0 ? 1 : -1;
            signed_temperature += step * sign * directions[index].temperature;
            signed_density += step * sign * directions[index].molar_density;
            for (std::size_t component = 0; component < composition.size(); ++component)
                signed_composition[component] += step * sign * directions[index].composition[component];
            for (const auto& parameter : directions[index].parameter) {
                for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
                    if (active[coordinate].family == parameter.family
                        && active[coordinate].indices == parameter.indices) {
                        parameter_shift[coordinate] += step * sign * parameter.value;
                    }
                }
            }
        }
        const auto shifted = p4_shifted_parameters(parameters, active, parameter_shift);
        const epcsaft::HelmholtzModel model(shifted, {}, thermochemistry);
        const auto shifted_state = model.state(
            signed_temperature, signed_density, signed_composition);
        const bool predicate_match =
            p4_take_signature(capture, thermochemistry, signed_temperature)
                == central_signature;
        predicate_unchanged = predicate_unchanged && predicate_match;
        CHECK(predicate_match);
        for (std::size_t i = 0; i < composition.size(); ++i) {
            const bool support_match =
                (composition[i] > 0.0) == (signed_composition[i] > 0.0);
            support_unchanged = support_unchanged && support_match;
            CHECK(support_match);
        }
        result += sign_product * p4_state_value(shifted_state, observable);
    }
    if (admission != nullptr) admission->insert(admission->end(), {
        step, static_cast<double>(support_unchanged),
        static_cast<double>(predicate_unchanged)});
    return result / std::pow(2.0 * step, static_cast<int>(directions.size()));
}

double p4_pressure_reference(
    const epcsaft::Parameters& parameters,
    const std::vector<epcsaft::ActiveParameter>& active,
    const std::optional<epcsaft::ThermochemistryRecord>& thermochemistry,
    double temperature, double pressure, const std::vector<double>& composition,
    const std::vector<epcsaft::PressureActionDirection>& directions,
    epcsaft::DensityBranch branch, epcsaft::PropertyObservable observable, double step,
    std::vector<double>* admission = nullptr
) {
    epcsaft::core::Capture capture;
    const epcsaft::HelmholtzModel central_model(parameters, {}, thermochemistry);
    const auto central = central_model.state_at_pressure(temperature, pressure, composition, branch);
    const auto central_signature = p4_take_signature(capture, thermochemistry, temperature);
    double minimum_slope = std::numeric_limits<double>::infinity();
    double maximum_closure = 0.0, maximum_site_residual = 0.0, minimum_pivot_ratio = 1.0;
    double result = 0.0;
    const std::size_t combinations = std::size_t{1} << directions.size();
    for (std::size_t mask = 0; mask < combinations; ++mask) {
        double signed_temperature = temperature;
        double signed_pressure = pressure;
        std::vector<double> signed_composition = composition;
        std::vector<double> parameter_shift(active.size(), 0.0);
        const int sign_product = [&] {
            int product = 1;
            for (std::size_t index = 0; index < directions.size(); ++index) {
                const double sign = (mask & (std::size_t{1} << index)) != 0 ? 1.0 : -1.0;
                product *= sign > 0.0 ? 1 : -1;
                signed_temperature += step * sign * directions[index].temperature;
                signed_pressure += step * sign * directions[index].pressure;
                for (std::size_t component = 0; component < composition.size(); ++component)
                    signed_composition[component] += step * sign * directions[index].composition[component];
                for (const auto& parameter : directions[index].parameter) {
                    for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate)
                        if (active[coordinate].family == parameter.family
                            && active[coordinate].indices == parameter.indices)
                            parameter_shift[coordinate] += step * sign * parameter.value;
                }
            }
            return product;
        }();
        const auto shifted_parameters = p4_shifted_parameters(parameters, active, parameter_shift);
        const epcsaft::HelmholtzModel shifted(shifted_parameters, {}, thermochemistry);
        const auto shifted_state = shifted.state_at_pressure(
            signed_temperature, signed_pressure, signed_composition, branch);
        const auto shifted_signature = p4_take_signature(
            capture, thermochemistry, signed_temperature);
        const double closure = std::abs(shifted_state.pressure - signed_pressure) / signed_pressure;
        minimum_slope = std::min(minimum_slope, shifted_state.pressure_density_derivative);
        maximum_closure = std::max(maximum_closure, closure);
        // Density closure is owned in reduced units; relative pressure error is amplified by 1/Z.
        const double reduced_closure = std::abs(shifted_state.pressure - signed_pressure)
            / (shifted_state.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * signed_temperature);
        CHECK(std::isfinite(shifted_state.pressure_density_derivative)
            && shifted_state.pressure_density_derivative > 0.0
            && reduced_closure <= epcsaft::policy::DENSITY_ROOT_RESIDUAL.value);
        CHECK(shifted_state.density_branch == central.density_branch
            && shifted_state.stable_root_count == central.stable_root_count
            && shifted_signature == central_signature);
        for (std::size_t i = 0; i < composition.size(); ++i)
            CHECK((composition[i] > 0.0) == (signed_composition[i] > 0.0));
        const auto residuals = p2_association_residuals(signed_temperature,
            1.0 / shifted_state.molar_density, signed_composition,
            shifted_parameters, shifted_state.site_fractions);
        for (const auto residual : residuals) {
            CHECK(std::isfinite(residual) && std::abs(residual) <= 2.0e-13);
            maximum_site_residual = std::max(maximum_site_residual, std::abs(residual));
        }
        const auto conditioning = p2_association_conditioning(signed_temperature,
            shifted_state.molar_density, signed_composition,
            shifted_parameters, shifted_state.site_fractions);
        CHECK(conditioning.rank == shifted_state.site_fractions.size()
            && conditioning.pivot_ratio > epcsaft::policy::ASSOCIATION_JACOBIAN_CONDITION.value);
        minimum_pivot_ratio = std::min(minimum_pivot_ratio, conditioning.pivot_ratio);
        result += sign_product * p4_state_value(shifted_state, observable);
    }
    if (admission != nullptr) {
        admission->insert(admission->end(), {step, minimum_slope, maximum_closure,
            maximum_site_residual, minimum_pivot_ratio,
            static_cast<double>(central.density_branch),
            central.stable_root_count.has_value()
                ? static_cast<double>(*central.stable_root_count) : -1.0});
    }
    return result / std::pow(2.0 * step, static_cast<int>(directions.size()));
}

// Differentiate the preceding action once, instead of subtracting O(h^k) caloric values.
// ALG-TERMINATION: cancellations scale with the fixed chart-coordinate terms.
struct PropertyActionReference final {
    std::vector<double> coarse, middle, fine, value, uncertainty, scale;
};

template <class Direction>
PropertyActionReference property_action_reference(
    const epcsaft::Parameters& parameters, const std::vector<epcsaft::ActiveParameter>& active,
    const std::optional<epcsaft::ThermochemistryRecord>& record, double T, double mechanical,
    const std::vector<double>& x, const std::vector<epcsaft::PropertyObservable>& observables,
    const std::vector<Direction>& directions, epcsaft::DensityBranch branch = epcsaft::DensityBranch::Unspecified
) {
    using namespace epcsaft;
    const HelmholtzModel model(parameters, active, record);
    const auto batch = [&](const HelmholtzModel& m, double t, double q, const std::vector<double>& y,
                           const std::vector<Direction>& request) {
        if constexpr (std::is_same_v<Direction, DensityActionDirection>)
            return m.density_actions(t, q, y, observables, request);
        else return m.pressure_actions(t, q, y, branch, observables, request);
    };
    PropertyActionReference result;
    result.scale.assign(observables.size(), 0.0);
    // Omit only zero sub-directions; retain the entire multilinear coordinate split.
    std::vector<Direction> terms;
    const auto expand = [&](const auto& self, std::size_t j) -> void {
        if (j == directions.size()) {
            const auto actions = batch(model, T, mechanical, x, terms);
            for (std::size_t i = 0; i < observables.size(); ++i) {
                CHECK(actions.results[i].status == PropertyActionStatus::Available && actions.results[i].action);
                result.scale[i] += std::abs(*actions.results[i].action);
            }
            return;
        }
        for (int part = 0; part < 4; ++part) {
            Direction d{};
            d.composition.assign(x.size(), 0.0);
            if (part == 0) d.temperature = directions[j].temperature;
            if (part == 1) {
                if constexpr (std::is_same_v<Direction, DensityActionDirection>) d.molar_density = directions[j].molar_density;
                else d.pressure = directions[j].pressure;
            }
            if (part == 2) d.composition = directions[j].composition;
            if (part == 3) d.parameter = directions[j].parameter;
            bool nonzero = d.temperature != 0.0 || !d.parameter.empty();
            if constexpr (std::is_same_v<Direction, DensityActionDirection>) nonzero = nonzero || d.molar_density != 0.0;
            else nonzero = nonzero || d.pressure != 0.0;
            for (const auto v : d.composition) nonzero = nonzero || v != 0.0;
            if (nonzero) { terms.push_back(d); self(self, j + 1); terms.pop_back(); }
        }
    };
    expand(expand, 0);
    epcsaft::core::Capture capture;
    const auto central_state = [&] {
        if constexpr (std::is_same_v<Direction, DensityActionDirection>) return model.state(T, mechanical, x);
        else return model.state_at_pressure(T, mechanical, x, branch);
    }();
    const auto central_signature = p4_take_signature(capture, record, T);
    const auto& last = directions.back();
    const std::vector<Direction> lower(directions.begin(), directions.end() - 1);
    const auto values = [&](double h) {
        std::vector<double> shift(active.size(), 0.0), y = x;
        for (std::size_t i = 0; i < x.size(); ++i) y[i] += h * last.composition[i];
        for (const auto& parameter : last.parameter)
            for (std::size_t a = 0; a < active.size(); ++a)
                if (parameter.family == active[a].family && parameter.indices == active[a].indices) shift[a] += h * parameter.value;
        const HelmholtzModel shifted(p4_shifted_parameters(parameters, active, shift), active, record);
        const double t = T + h * last.temperature;
        double q = mechanical;
        if constexpr (std::is_same_v<Direction, DensityActionDirection>) q += h * last.molar_density;
        else q += h * last.pressure;
        const auto state = [&] {
            if constexpr (std::is_same_v<Direction, DensityActionDirection>) return shifted.state(t, q, y);
            else return shifted.state_at_pressure(t, q, y, branch);
        }();
        CHECK(p4_take_signature(capture, record, t) == central_signature);
        for (std::size_t i = 0; i < x.size(); ++i) CHECK((x[i] > 0.0) == (y[i] > 0.0));
        if constexpr (!std::is_same_v<Direction, DensityActionDirection>) {
            CHECK(state.density_branch == central_state.density_branch && state.stable_root_count == central_state.stable_root_count);
            CHECK(std::isfinite(state.pressure_density_derivative) && state.pressure_density_derivative > 0.0);
            CHECK(std::abs(state.pressure - q) / (state.molar_density * GAS_CONSTANT_J_PER_MOL_K * t)
                <= policy::DENSITY_ROOT_RESIDUAL.value);
            const auto residuals = p2_association_residuals(t, 1.0 / state.molar_density, y,
                shifted.parameters(), state.site_fractions);
            for (const auto residual : residuals) CHECK(std::isfinite(residual) && std::abs(residual) <= 2.0e-13);
            const auto conditioning = p2_association_conditioning(t, state.molar_density, y,
                shifted.parameters(), state.site_fractions);
            CHECK(conditioning.rank == state.site_fractions.size()
                && conditioning.pivot_ratio > policy::ASSOCIATION_JACOBIAN_CONDITION.value);
        }
        std::vector<double> values;
        if (lower.empty()) for (const auto obs : observables) values.push_back(p4_state_value(state, obs));
        else {
            const auto actions = batch(shifted, t, q, y, lower);
            for (const auto& action : actions.results) {
                CHECK(action.status == PropertyActionStatus::Available && action.action);
                values.push_back(*action.action);
            }
        }
        return values;
    };
    // Higher-action subtraction at 0.0025 is density-root rounding dominated;
    // these nodes resolve cp/cv without crossing support or thermal pieces.
    constexpr double steps[] = {0.1, 0.05, 0.025};
    std::vector<std::vector<double>> rungs;
    std::vector<double> rounding;
    for (const double h : steps) {
        const auto plus = values(h), minus = values(-h);
        std::vector<double> stencil;
        rounding.clear();
        for (std::size_t i = 0; i < observables.size(); ++i) {
            stencil.push_back((plus[i] - minus[i]) / (2.0 * h));
            rounding.push_back(32.0 * std::numeric_limits<double>::epsilon() * (std::abs(plus[i]) + std::abs(minus[i])) / (2.0 * h));
        }
        rungs.push_back(stencil);
    }
    result.coarse = rungs[0]; result.middle = rungs[1]; result.fine = rungs[2];
    for (std::size_t i = 0; i < observables.size(); ++i) {
        const double refined = (4.0 * rungs[2][i] - rungs[1][i]) / 3.0;
        const double coarser = (4.0 * rungs[1][i] - rungs[0][i]) / 3.0;
        result.value.push_back(refined);
        // Polynomial ideal cp/cv has identically zero higher chart terms.
        // An exactly constant preceding action subtracts exactly: its nonzero
        // value must not manufacture a rounding floor for a zero stencil.
        const bool exact_zero = result.scale[i] == 0.0 && rungs[0][i] == 0.0
            && rungs[1][i] == 0.0 && rungs[2][i] == 0.0;
        result.uncertainty.push_back(exact_zero ? 0.0 : std::abs(refined - coarser) + rounding[i]);
    }
    return result;
}

// EOS sections 6--7: charged-coordinate higher actions on a regular positive-ionic-strength state.
void electrolyte_parameter_actions() {
    using namespace epcsaft;
    auto parameters = properties_support::aqueous_electrolyte();
    parameters.born_c_shell = parameters.born_c_dielectric = 1.0;
    const std::vector<ActiveParameter> active = {
        {ParameterFamily::BornDiameter, {1}}, {ParameterFamily::DebyeHuckelDiameter, {1}},
        {ParameterFamily::PackingDiameter, {1}}, {ParameterFamily::IonSuppressionCoefficient, {1}}};
    const std::vector<double> x = {0.98, 0.01, 0.01};
    constexpr double T = 298.15, rho = 55000.0;
    const HelmholtzModel model(parameters, active);
    const auto observables = std::vector<PropertyObservable>{PropertyObservable::ResidualEnthalpy};
    for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
        const DensityActionDirection parameter{0.0, 0.0, {0.0, 0.0, 0.0},
            {{active[coordinate].family, active[coordinate].indices, 1.0}}};
        for (std::size_t order = 2; order <= 3; ++order) {
            std::vector<DensityActionDirection> directions(order, parameter);
            if (coordinate != 0) directions[1] = {1.0, 0.0, {0.0, 0.0, 0.0}, {}};
            const auto result = model.density_actions(T, rho, x, observables, directions).results.front();
            CHECK(result.status == PropertyActionStatus::Available && result.action);
            const std::string id = "electrolyte-coordinate/" + std::to_string(coordinate) + "/A" + std::to_string(order);
            double reference, uncertainty = 0.0, allowance;
            std::string inputs;
            if (coordinate == 0) {
                // Independently evaluated corrected radial charging formula, h_B = -K/d_B.
                const double epsilon_bulk = parameters.relative_permittivity[0]
                    / (1.0 + x[1] * parameters.ion_suppression_coefficient[1]
                        + x[2] * parameters.ion_suppression_coefficient[2]);
                double f_mix = 0.0;
                for (std::size_t i = 0; i < x.size(); ++i) f_mix += x[i] * parameters.solvation_factor[i];
                const double shell_ratio = 1.0 + (f_mix - 1.0) / std::abs(parameters.charge_number[1]);
                const double outer = (1.0 - 1.0 / epsilon_bulk) / shell_ratio;
                const double inner = (1.0 - 1.0 / parameters.ionic_region_relative_permittivity) * (1.0 - 1.0 / shell_ratio);
                const double prefactor = GAS_CONSTANT_J_PER_MOL_K * ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
                    * x[1] * parameters.charge_number[1] * parameters.charge_number[1]
                    / (4.0 * PI * VACUUM_PERMITTIVITY_F_PER_M * BOLTZMANN_J_PER_K * ANGSTROM_M);
                const double d = parameters.born_diameter_angstrom[1];
                const double derivative_factor = order == 2 ? -2.0 / std::pow(d, 3) : 6.0 / std::pow(d, 4);
                reference = derivative_factor * prefactor * (inner + outer);
                const double source_scale = std::abs(derivative_factor * prefactor) * (std::abs(inner) + std::abs(outer));
                allowance = 1.0e-10 * std::abs(reference) + 32.0 * std::numeric_limits<double>::epsilon() * source_scale;
                inputs = recording_support::fields({recording_support::field("K_J_A_per_mol", prefactor * (inner + outer)),
                    recording_support::field("d_B_A", d), recording_support::field("source_scale", source_scale)});
            } else {
                auto lower = directions;
                lower.pop_back();
                const auto stencil = [&](double h) {
                    const auto value = [&](double shift) {
                        std::vector<double> delta(active.size(), 0.0);
                        delta[coordinate] = order == 3 ? shift : 0.0;
                        const HelmholtzModel rebuilt(p4_shifted_parameters(parameters, active, delta), active);
                        const auto action = rebuilt.density_actions(order == 2 ? T + shift : T, rho, x, observables, lower).results.front();
                        CHECK(action.status == PropertyActionStatus::Available && action.action);
                        return *action.action;
                    };
                    return (value(h) - value(-h)) / (2.0 * h);
                };
                // A 0.01 suppression step loses the small T/parameter action
                // to cancellation in the preceding action; resolve it at larger nodes.
                const double step = coordinate == 3 && order == 3 ? 4.0 : 0.4;
                const double coarse = stencil(step), middle = stencil(step / 2.0), fine = stencil(step / 4.0);
                const double rich_coarse = (4.0 * middle - coarse) / 3.0;
                reference = (4.0 * fine - middle) / 3.0;
                uncertainty = std::abs(reference - rich_coarse);
                allowance = (order == 2 ? 5.0e-5 : 5.0e-4) * std::abs(reference);
                inputs = recording_support::fields({recording_support::field("coarse", coarse),
                    recording_support::field("middle", middle), recording_support::field("fine", fine),
                    recording_support::field("reference_uncertainty", uncertainty)});
            }
            const bool ok = std::abs(*result.action - reference) <= allowance
                && uncertainty <= allowance && std::abs(reference) > allowance;
            recording_support::compare_custom(id.c_str(), "enthalpy_parameter_action", "native", *result.action, reference,
                recording_support::fields({inputs, recording_support::field("T_K", T), recording_support::field("rho_mol_per_m3", rho),
                    recording_support::field("x", x), recording_support::field("allowance", allowance)}),
                "Relationship", coordinate == 0 ? "Born-charging-identity" : "preceding-action-Richardson-stencil",
                "nonzero second and third actions in charged EOS coordinates", "magnitude_scaled_with_resolved_reference", ok);
            CHECK(ok);
        }
    }
}

epcsaft::Parameters p4_association_stress(std::size_t sites, double q) {
    constexpr double temperature = 300.0;
    constexpr double density = 1000.0;
    auto parameters = properties_support::ethane();
    parameters.site_component.assign(sites, 0);
    parameters.site_multiplicity.assign(sites, 1);
    parameters.association_active.assign(sites * sites, 0);
    parameters.association_energy_over_k.assign(sites * sites, 0.0);
    parameters.association_volume.assign(sites * sites, 0.0);
    parameters.association_active[0] = 1;
    parameters.association_volume[0] = 0.1;
    parameters.association_energy_over_k[0] = temperature * std::log(2.0);
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const double unit_delta = epcsaft::association_delta(
        temperature, density, {1.0}, parameters, sigma, diameter)[0];
    const double number_density = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    parameters.association_energy_over_k[0] = temperature
        * std::log1p(q / (number_density * unit_delta));
    return parameters;
}

double p4_residual_enthalpy_reference(const epcsaft::Parameters& parameters) {
    constexpr double temperature = 300.0;
    constexpr double volume = 1.0e-3;
    constexpr double temperature_step = 1.0e-3;
    constexpr double volume_step = 1.0e-8;
    epcsaft::HelmholtzModel model(parameters);
    const auto potential = [&](double t, double v) {
        return model.evaluate(t, v, {1.0}, epcsaft::DerivativeOrder::Value).value;
    };
    const double temperature_derivative = (
        potential(temperature + temperature_step, volume)
        - potential(temperature - temperature_step, volume)) / (2.0 * temperature_step);
    const double volume_derivative = (
        potential(temperature, volume + volume_step)
        - potential(temperature, volume - volume_step)) / (2.0 * volume_step);
    return -epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature
        * (temperature * temperature_derivative + volume * volume_derivative);
}

void properties_p4() {
    const std::vector<epcsaft::ActiveParameter> active = {
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}},
        {epcsaft::ParameterFamily::AssociationVolume, {1, 4}},
    };
    epcsaft::HelmholtzModel model(
        properties_support::associating_binary(), active, p3_binary_record());
    const std::vector<double> composition = {0.25, 0.75};
    const auto observables = std::vector<epcsaft::PropertyObservable>{
        epcsaft::PropertyObservable::IdealEnthalpy,
        epcsaft::PropertyObservable::IdealInternalEnergy,
        epcsaft::PropertyObservable::IdealEntropy,
        epcsaft::PropertyObservable::IdealIsobaricHeatCapacity,
        epcsaft::PropertyObservable::IdealIsochoricHeatCapacity,
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualInternalEnergy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity,
        epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity,
        epcsaft::PropertyObservable::TotalEnthalpy,
        epcsaft::PropertyObservable::TotalInternalEnergy,
        epcsaft::PropertyObservable::TotalEntropy,
        epcsaft::PropertyObservable::TotalIsobaricHeatCapacity,
        epcsaft::PropertyObservable::TotalIsochoricHeatCapacity,
    };
    const std::array<epcsaft::DensityActionDirection, 3> directions = {{
        {0.2, 5.0, {-0.01, 0.01}, {
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}, 0.03},
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, 0.001},
        }},
        {-0.1, -3.0, {0.005, -0.005}, {
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}, -0.02},
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, -0.0005},
        }},
        {0.05, 2.0, {0.002, -0.002}, {
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}, 0.01},
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, 0.0002},
        }},
    }};
    recording_support::expect_group(recording_support::active_test_id(),
        "p4-direct", "direct/", observables.size() * directions.size() * 2,
        "P4 direct density directional property actions", "observable,order");
    recording_support::expect_group(recording_support::active_test_id(),
        "p4-pressure", "pressure/", observables.size() * directions.size() * 2,
        "P4 fixed-pressure directional property actions", "observable,order");
    recording_support::expect_group(recording_support::active_test_id(),
        "p4-supplemental", "supplemental/", 149,
        "P4 availability, selection and structural counter checks", "case");
    const auto state = model.state(313.15, 42441.09427799721, composition);
    const auto pressure_root = model.state_at_pressure(
        313.15, state.pressure, composition, epcsaft::DensityBranch::Unspecified);
    const std::array<epcsaft::PressureActionDirection, 3> pressure_directions_all = {{
        {directions[0].temperature, 15000.0, directions[0].composition, directions[0].parameter},
        {directions[1].temperature, -10000.0, directions[1].composition, directions[1].parameter},
        {directions[2].temperature, 6000.0, directions[2].composition, directions[2].parameter},
    }};
    for (std::size_t order = 1; order <= directions.size(); ++order) {
        const std::vector<epcsaft::DensityActionDirection> request(
            directions.begin(), directions.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = model.density_actions(
            313.15, 42441.09427799721, composition, observables, request);
        CHECK(batch.results.size() == observables.size());
        const auto reference = property_action_reference(model.parameters(), active, p3_binary_record(),
            313.15, 42441.09427799721, composition, observables, request);
        for (std::size_t index = 0; index < observables.size(); ++index) {
            const auto& result = batch.results[index];
            const bool available = result.status == epcsaft::PropertyActionStatus::Available
                && result.value.has_value() && result.action.has_value()
                && std::isfinite(*result.value) && std::isfinite(*result.action);
            const double fine = reference.value[index], middle = reference.middle[index], coarse = reference.coarse[index];
            const double tolerance = std::array<double, 3>{5.0e-6, 5.0e-5, 5.0e-4}[order - 1];
            const bool action_ok = available && std::abs(*result.action - fine)
                <= tolerance * reference.scale[index];
            const bool refinement_ok = reference.uncertainty[index] <= tolerance * reference.scale[index];
            recording_support::compare_custom(
                (std::string("direct/action/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "action", "native",
                available ? *result.action : std::numeric_limits<double>::quiet_NaN(), fine,
                recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order),
                    recording_support::field("central", p4_state_value(state, observables[index])),
                    recording_support::field("energy_order", result.energy_order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-refinement", "direct action matches independent refinement",
                "tolerance_times_chart_term_scale", action_ok, "p4-direct",
                (std::string("action/") + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                "action", "T,rho,x,theta", "available", "observable,order");
            recording_support::compare_custom(
                (std::string("direct/refinement/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "reference_refinement", "native",
                fine, middle,
                recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-ladder", "direct reference refines",
                "tolerance_times_chart_term_scale", refinement_ok, "p4-direct",
                (std::string("refinement/") + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                "reference", "T,rho,x,theta,h", "refined", "observable,order");
            CHECK(action_ok && refinement_ok);
        }
    }
    for (std::size_t order = 1; order <= directions.size(); ++order) {
        const std::vector<epcsaft::PressureActionDirection> request(
            pressure_directions_all.begin(),
            pressure_directions_all.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = model.pressure_actions(
            313.15, state.pressure, composition, pressure_root.density_branch,
            observables, request);
        const auto reference = property_action_reference(model.parameters(), active, p3_binary_record(),
            313.15, state.pressure, composition, observables, request, pressure_root.density_branch);
        for (std::size_t index = 0; index < observables.size(); ++index) {
            const auto& result = batch.results[index];
            const double fine = reference.value[index], middle = reference.middle[index], coarse = reference.coarse[index];
            const double tolerance = std::array<double, 3>{5.0e-6, 5.0e-5, 5.0e-4}[order - 1];
            const bool action_ok = result.status == epcsaft::PropertyActionStatus::Available
                && result.action.has_value() && std::abs(*result.action - fine)
                    <= tolerance * reference.scale[index];
            const bool refinement_ok = reference.uncertainty[index] <= tolerance * reference.scale[index];
            recording_support::compare_custom(
                (std::string("pressure/action/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "action", "native",
                result.action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order), recording_support::field("central",
                        p4_state_value(pressure_root, observables[index])),
                    recording_support::field("energy_order", result.energy_order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-refinement",
                "fixed-pressure action matches independent refinement",
                "tolerance_times_chart_term_scale", action_ok,
                "p4-pressure", (std::string("action/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "action", "T,P,x,theta", "available",
                "observable,order");
            recording_support::compare_custom(
                (std::string("pressure/refinement/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "reference_refinement", "native", fine,
                middle, recording_support::fields({recording_support::field("observable", index),
                    recording_support::field("order", order),
                    recording_support::field("reference_fine", fine),
                    recording_support::field("reference_middle", middle),
                    recording_support::field("reference_coarse", coarse),
                    recording_support::field("tolerance", tolerance),
                    recording_support::field("term_scale", reference.scale[index]),
                    recording_support::field("reference_uncertainty", reference.uncertainty[index]),
                    recording_support::field("allowance", tolerance * reference.scale[index])}),
                "Relationship", "central-difference-ladder", "pressure reference refines",
                "tolerance_times_chart_term_scale", refinement_ok, "p4-pressure",
                (std::string("refinement/") + std::to_string(index) + "/"
                    + std::to_string(order)).c_str(), "reference", "T,P,x,theta,h", "refined",
                "observable,order");
            CHECK(action_ok && refinement_ok);
        }
    }
    const std::array<epcsaft::PropertyObservable, 5> thermal_observables = {
        epcsaft::PropertyObservable::IdealEnthalpy,
        epcsaft::PropertyObservable::IdealInternalEnergy,
        epcsaft::PropertyObservable::IdealEntropy,
        epcsaft::PropertyObservable::IdealIsobaricHeatCapacity,
        epcsaft::PropertyObservable::IdealIsochoricHeatCapacity};
    const std::vector<epcsaft::PropertyObservable> thermal_request_observables(
        thermal_observables.begin(), thermal_observables.end());
    const double gas_constant = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const std::array<std::pair<const char*, epcsaft::IdealCorrelation>, 5> thermal_forms = {{
        {"polynomial", epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{30.0, 0.02, 1.0e-4, 1.0e-4}, 300.0, 1000.0, 50.0},
            1.0e5}},
        {"shomate", epcsaft::IdealCorrelation{
            epcsaft::IdealShomate{{30.0, 20.0, 10.0, 5.0, 2.0, 0.0, 0.0, 0.0}, 1000.0},
            1.0e5}},
        {"nasa7", epcsaft::IdealCorrelation{
            epcsaft::IdealNASA7{{30.0 / gas_constant, 0.02 / gas_constant,
                2.0e-5 / gas_constant, 1.0e-7 / gas_constant, 1.0e-9 / gas_constant,
                0.0, 0.0}}, 1.0e5}},
        {"nasa9", epcsaft::IdealCorrelation{
            epcsaft::IdealNASA9{{0.0, 0.0, 30.0 / gas_constant, 0.02 / gas_constant,
                2.0e-5 / gas_constant, 1.0e-7 / gas_constant, 1.0e-9 / gas_constant,
                0.0, 0.0}}, 1.0e5}},
        {"dippr107", epcsaft::IdealCorrelation{
            epcsaft::IdealDIPPR107{{0.0, 30.0, 250.0, 10.0, 600.0},
                300.0, 1000.0, 50.0}, 1.0e5}}
    }};
    const std::array<epcsaft::DensityActionDirection, 3> thermal_directions = {{
        {20.0, 0.0, {0.0}, {}}, {-10.0, 0.0, {0.0}, {}}, {5.0, 0.0, {0.0}, {}}}};
    for (const auto& [name, correlation] : thermal_forms) {
        const auto thermal_record = p3_record(correlation);
        epcsaft::HelmholtzModel thermal_model(
            properties_support::ethane(), {}, thermal_record);
        for (std::size_t order = 1; order <= thermal_directions.size(); ++order) {
            const std::vector<epcsaft::DensityActionDirection> request(
                thermal_directions.begin(), thermal_directions.begin() + static_cast<std::ptrdiff_t>(order));
            const auto batch = thermal_model.density_actions(
                320.0, 1000.0, {1.0}, thermal_request_observables, request);
            for (std::size_t index = 0; index < thermal_observables.size(); ++index) {
                const double fine = p4_reference(
                    thermal_model.parameters(), {}, thermal_record, 320.0, 1000.0,
                    {1.0}, request, thermal_observables[index],
                    std::array<double, 3>{1.0e-3, 3.0e-3, 1.0e-2}[order - 1]);
                const bool nonconstant_cp = order != 3 || index != 3 || std::abs(fine) > 5.0e-4;
                const bool ok = batch.results[index].status == epcsaft::PropertyActionStatus::Available
                    && batch.results[index].action.has_value()
                    && close(*batch.results[index].action, fine, 5.0e-4)
                    && nonconstant_cp;
                recording_support::compare_custom(
                    (std::string("supplemental/thermal-form/") + name + "/"
                        + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                    "action", "native",
                    batch.results[index].action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                    recording_support::fields({recording_support::field("form", name),
                        recording_support::field("observable", index), recording_support::field("order", order),
                        recording_support::field("reference", fine),
                        recording_support::field("nonconstant_cp", nonconstant_cp ? 1 : 0)}),
                    "Relationship", "central-difference-refinement", "thermal form action",
                    "relative", ok, "p4-supplemental", (std::string("thermal-form/") + name
                        + "/" + std::to_string(index) + "/" + std::to_string(order)).c_str(),
                    "action", "T,rho,x", "available", "form,observable,order");
                CHECK(ok);
            }
        }
    }
    const auto multiroot_directions = std::array<epcsaft::PressureActionDirection, 2>{
        epcsaft::PressureActionDirection{0.1, 2000.0, {0.0}, {}},
        epcsaft::PressureActionDirection{-0.05, -1000.0, {0.0}, {}}};
    epcsaft::HelmholtzModel ethane_model(properties_support::ethane());
    for (const auto& [site_count, q, expected_status, name]
         : std::array<std::tuple<std::size_t, double, epcsaft::PropertyActionStatus,
             const char*>, 2>{{
             {2, 1.0e26, epcsaft::PropertyActionStatus::AssociationIllConditioned,
              "association-ill-conditioned"},
             {128, 1.0e27, epcsaft::PropertyActionStatus::AssociationSingular,
              "association-singular"}}}) {
        const auto parameters = p4_association_stress(site_count, q);
        const auto sigma = epcsaft::segment_sigma(300.0, parameters);
        const auto diameter = epcsaft::segment_diameter(300.0, parameters, sigma);
        const auto sites = epcsaft::solve_site_fractions(
            300.0, 1000.0, {1.0}, parameters, sigma, diameter);
        const double exact_site = 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * q));
        const auto residuals = p2_association_residuals(
            300.0, 1.0e-3, {1.0}, parameters, sites);
        const auto conditioning = p2_association_conditioning(
            300.0, 1000.0, {1.0}, parameters, sites);
        epcsaft::HelmholtzModel probe(
            parameters, {}, p3_record(p3_forms()[0].second));
        const auto batch = probe.density_actions(300.0, 1000.0, {1.0},
            {epcsaft::PropertyObservable::ResidualEnthalpy,
             epcsaft::PropertyObservable::TotalEnthalpy,
             epcsaft::PropertyObservable::IdealEnthalpy},
            {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
        const double expected_residual_enthalpy = p4_residual_enthalpy_reference(parameters);
        const double maximum_residual = std::accumulate(residuals.begin(), residuals.end(), 0.0,
            [](double value, double residual) { return std::max(value, std::abs(residual)); });
        const bool inactive_exact = std::all_of(sites.begin() + 1, sites.end(),
            [](double value) { return value == 1.0; });
        const bool ok = batch.results.size() == 3 && sites.size() == site_count
            && close(sites[0], exact_site, 2.0e-12) && inactive_exact
            && maximum_residual <= epcsaft::policy::ASSOCIATION_RESIDUAL.value
            && batch.results[0].status == expected_status
            && batch.results[0].value.has_value() && !batch.results[0].action.has_value()
            && close(*batch.results[0].value, expected_residual_enthalpy, 2.0e-6)
            && batch.results[1].status == expected_status
            && batch.results[1].value.has_value() && !batch.results[1].action.has_value()
            && batch.results[2].status == epcsaft::PropertyActionStatus::Available
            && batch.results[2].value.has_value() && batch.results[2].action.has_value()
            && close(*batch.results[1].value,
                *batch.results[0].value + *batch.results[2].value, 2.0e-12)
            && batch.diagnostics.site_rank == conditioning.rank;
        recording_support::compare_custom(
            (std::string("supplemental/") + name).c_str(), "status_and_preservation", "bool",
            ok ? 1.0 : 0.0, 1.0,
            recording_support::fields({recording_support::field("site_count", site_count),
                recording_support::field("q", q),
                recording_support::field("site_fraction", sites[0]),
                recording_support::field("exact_site_fraction", exact_site),
                recording_support::field("maximum_primal_residual", maximum_residual),
                recording_support::field("rank", conditioning.rank),
                recording_support::field("pivot_ratio", conditioning.pivot_ratio),
                recording_support::field("residual_enthalpy", batch.results[0].value.value_or(
                    std::numeric_limits<double>::quiet_NaN())),
                recording_support::field("reference_residual_enthalpy", expected_residual_enthalpy),
                recording_support::field("central_ideal_enthalpy",
                    batch.results[2].value.value_or(std::numeric_limits<double>::quiet_NaN()))}),
            "Availability", "analytic-association-stress",
            "typed association failure preserves central residual, total and ideal values", "absolute", ok,
            "p4-supplemental", name, "status", "T,rho,x,q", "typed failure");
        CHECK(ok);
    }
    const auto missing_thermochemistry = ethane_model.density_actions(
        300.0, 1000.0, {1.0},
        {epcsaft::PropertyObservable::IdealEnthalpy},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool missing_ok = missing_thermochemistry.results[0].status
        == epcsaft::PropertyActionStatus::MissingThermochemistry;
    recording_support::compare_custom("supplemental/missing-thermochemistry", "status", "bool",
        missing_ok ? 1.0 : 0.0, 1.0, "ethane without source", "Availability", "thermo",
        "ideal source absence is typed", "absolute", missing_ok,
        "p4-supplemental", "missing-thermochemistry", "status", "thermo", "missing_thermochemistry");
    const auto partial = ethane_model.density_actions(
        300.0, 1000.0, {1.0},
        {epcsaft::PropertyObservable::IdealEnthalpy, epcsaft::PropertyObservable::ResidualEnthalpy},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool partial_ok = partial.results[0].status == epcsaft::PropertyActionStatus::MissingThermochemistry
        && partial.results[1].status == epcsaft::PropertyActionStatus::Available
        && partial.results[1].value.has_value();
    recording_support::compare_custom("supplemental/partial-availability", "sibling_preserved", "bool",
        partial_ok ? 1.0 : 0.0, 1.0, "ideal+residual without source", "Availability", "partial",
        "residual sibling remains available", "absolute", partial_ok,
        "p4-supplemental", "partial-availability", "sibling_preserved", "sibling", "partial");
    const auto negative_slope = ethane_model.density_actions(
        280.0, 4000.0, {1.0},
        {epcsaft::PropertyObservable::ResidualEnthalpy,
         epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity,
         epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool negative_ok = negative_slope.results[0].status == epcsaft::PropertyActionStatus::Available
        && negative_slope.results[1].status == epcsaft::PropertyActionStatus::PressureSingular
        && negative_slope.results[0].value.has_value()
        && negative_slope.results[2].status == epcsaft::PropertyActionStatus::Available
        && negative_slope.results[2].action.has_value();
    recording_support::compare_custom("supplemental/negative-slope", "status", "bool",
        negative_ok ? 1.0 : 0.0, 1.0, "ethane T=280,rho=4000", "Availability", "pressure",
        "negative pressure slope suppresses residual cp only", "absolute", negative_ok,
        "p4-supplemental", "negative-slope", "status", "pressure", "pressure_singular");
    const auto density_unavailable = ethane_model.pressure_actions(
        200.0, 217800.69494692178, {1.0}, epcsaft::DensityBranch::Unspecified,
        {epcsaft::PropertyObservable::ResidualEnthalpy},
        {epcsaft::PressureActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool density_unavailable_ok = density_unavailable.results.size() == 1
        && density_unavailable.results[0].status
            == epcsaft::PropertyActionStatus::DensityUnavailable;
    recording_support::compare_custom(
        "supplemental/density-unavailable", "status", "bool",
        density_unavailable_ok ? 1.0 : 0.0, 1.0,
        "ethane two-root state with unspecified branch", "Availability", "density-selection",
        "ambiguous fixed-pressure density remains unavailable", "absolute", density_unavailable_ok,
        "p4-supplemental", "density-unavailable", "status", "T,P,x,branch",
        "density_unavailable");
    CHECK(density_unavailable_ok);
    const std::vector<epcsaft::ActiveParameter> entropy_active = {
        {epcsaft::ParameterFamily::AssociationVolume, {1, 4}}
    };
    const epcsaft::PressureActionDirection entropy_direction{
        0.0, 0.0, {0.0, 0.0}, {
            {epcsaft::ParameterFamily::AssociationVolume, {1, 4}, 1.0}
        }};
    const auto entropy_batch = model.pressure_actions(
        313.15, state.pressure, composition, pressure_root.density_branch,
        {epcsaft::PropertyObservable::IdealEntropy}, {entropy_direction});
    const epcsaft::HelmholtzModel entropy_model(
        model.parameters(), entropy_active, p3_binary_record());
    const auto entropy_derivatives = entropy_model.evaluate(
        313.15, 1.0 / state.molar_density, composition, epcsaft::DerivativeOrder::Hessian);
    const std::size_t entropy_parameter_column = 2 + composition.size();
    const std::size_t entropy_point_width = entropy_parameter_column + entropy_active.size();
    const double pressure_parameter = -epcsaft::GAS_CONSTANT_J_PER_MOL_K * 313.15
        * entropy_derivatives.hessian[entropy_point_width + entropy_parameter_column];
    const double entropy_chain = epcsaft::GAS_CONSTANT_J_PER_MOL_K * pressure_parameter
        / (state.molar_density * state.pressure_density_derivative);
    const double entropy_fine = p4_pressure_reference(
        model.parameters(), entropy_active, p3_binary_record(), 313.15, state.pressure,
        composition, {entropy_direction}, pressure_root.density_branch,
        epcsaft::PropertyObservable::IdealEntropy, 1.0e-3);
    const double entropy_middle = p4_pressure_reference(
        model.parameters(), entropy_active, p3_binary_record(), 313.15, state.pressure,
        composition, {entropy_direction}, pressure_root.density_branch,
        epcsaft::PropertyObservable::IdealEntropy, 3.0e-3);
    const double entropy_coarse = p4_pressure_reference(
        model.parameters(), entropy_active, p3_binary_record(), 313.15, state.pressure,
        composition, {entropy_direction}, pressure_root.density_branch,
        epcsaft::PropertyObservable::IdealEntropy, 1.0e-2);
    const bool entropy_ok = entropy_batch.results.size() == 1
        && entropy_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && entropy_batch.results[0].action.has_value()
        && close(*entropy_batch.results[0].action, entropy_chain, 5.0e-6)
        && close(*entropy_batch.results[0].action, entropy_fine, 5.0e-6)
        && close(entropy_fine, entropy_middle, 1.0e-5);
    recording_support::compare_custom(
        "supplemental/parameter-ideal-entropy", "action", "J/(mol K)",
        entropy_batch.results.empty() ? std::numeric_limits<double>::quiet_NaN()
            : entropy_batch.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()),
        entropy_fine,
        recording_support::fields({recording_support::field("T", 313.15),
            recording_support::field("P", state.pressure),
            recording_support::field("rho", state.molar_density),
            recording_support::field("P_theta", pressure_parameter),
            recording_support::field("chain", entropy_chain),
            recording_support::field("reference_fine", entropy_fine),
            recording_support::field("reference_middle", entropy_middle),
            recording_support::field("reference_coarse", entropy_coarse)}),
        "Relationship", "implicit-pressure-chain", "parameter Ideal entropy retains P_theta/(rho P_rho)",
        "relative", entropy_ok, "p4-supplemental", "parameter-ideal-entropy", "action",
        "T,P,x,theta", "available", "reference_fine,reference_middle,reference_coarse");
    for (const auto branch : {epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid}) {
        const auto root = ethane_model.state_at_pressure(
            200.0, 217800.69494692178, {1.0}, branch);
        const std::vector<epcsaft::PropertyObservable> rows = {
            epcsaft::PropertyObservable::ResidualEnthalpy,
            epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity};
        for (std::size_t order = 1; order <= 2; ++order) {
            const std::vector<epcsaft::PressureActionDirection> request(
                multiroot_directions.begin(), multiroot_directions.begin() + static_cast<std::ptrdiff_t>(order));
            const auto branch_batch = ethane_model.pressure_actions(
                200.0, 217800.69494692178, {1.0}, branch, rows, request);
            for (std::size_t index = 0; index < rows.size(); ++index) {
                const double fine = p4_pressure_reference(
                    properties_support::ethane(), {}, std::nullopt, 200.0,
                    217800.69494692178, {1.0}, request, branch, rows[index],
                    std::array<double, 2>{3.0e-3, 1.0e-2}[order - 1]);
                const double coarse = p4_pressure_reference(
                    properties_support::ethane(), {}, std::nullopt, 200.0,
                    217800.69494692178, {1.0}, request, branch, rows[index],
                    std::array<double, 2>{1.0e-2, 3.0e-2}[order - 1]);
                const double tolerance = std::array<double, 2>{5.0e-6, 5.0e-5}[order - 1];
                const bool action_ok = branch_batch.results[index].status == epcsaft::PropertyActionStatus::Available
                    && branch_batch.results[index].action.has_value()
                    && close(*branch_batch.results[index].action, fine, tolerance);
                const bool refinement_ok = close(fine, coarse, 2.0 * tolerance);
                const std::string point = std::string(branch == epcsaft::DensityBranch::Vapor ? "vapor" : "liquid")
                    + "/" + std::to_string(index) + "/" + std::to_string(order);
                recording_support::compare_custom(("supplemental/multiroot/action/" + point).c_str(),
                    "action", "native",
                    branch_batch.results[index].action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                    recording_support::fields({recording_support::field("branch", point), recording_support::field("order", order)}),
                    "Relationship", "central-difference-refinement", "multiroot pressure action",
                    "relative", action_ok, "p4-supplemental", ("multiroot/action/" + point).c_str(),
                    "action", "T,P,x", "available", "branch,observable,order");
                recording_support::compare_custom(("supplemental/multiroot/refinement/" + point).c_str(),
                    "reference_refinement", "native", fine, coarse,
                    recording_support::fields({recording_support::field("branch", point), recording_support::field("order", order)}),
                    "Relationship", "central-difference-ladder", "multiroot reference refines",
                    "relative", refinement_ok, "p4-supplemental", ("multiroot/refinement/" + point).c_str(),
                    "reference", "T,P,x,h", "refined", "branch,observable,order");
                CHECK(action_ok && refinement_ok);
            }
        }
        CHECK(root.density_branch == branch);
    }
    const std::array<epcsaft::DensityActionDirection, 3> electrolyte_directions = {{
        {0.1, 1.0, {-0.02, 0.01, 0.01}, {}},
        {-0.05, -0.5, {0.01, -0.005, -0.005}, {}},
        {0.03, 0.2, {0.004, -0.002, -0.002}, {}}}};
    const std::vector<epcsaft::PropertyObservable> electrolyte_observables = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity};
    epcsaft::HelmholtzModel electrolyte_model(properties_support::aqueous_electrolyte());
    for (std::size_t order = 1; order <= electrolyte_directions.size(); ++order) {
        const std::vector<epcsaft::DensityActionDirection> request(
            electrolyte_directions.begin(), electrolyte_directions.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = electrolyte_model.density_actions(
            300.0, 1000.0, {0.98, 0.01, 0.01}, electrolyte_observables, request);
        for (std::size_t index = 0; index < electrolyte_observables.size(); ++index) {
            const double fine = p4_reference(
                electrolyte_model.parameters(), {}, std::nullopt, 300.0, 1000.0,
                {0.98, 0.01, 0.01}, request, electrolyte_observables[index],
                std::array<double, 3>{1.0e-3, 3.0e-3, 1.0e-2}[order - 1]);
            const bool ok = batch.results[index].status == epcsaft::PropertyActionStatus::Available
                && batch.results[index].action.has_value()
                && close(*batch.results[index].action, fine,
                    std::array<double, 3>{5.0e-6, 5.0e-5, 5.0e-4}[order - 1]);
            const std::string point = std::to_string(index) + "/" + std::to_string(order);
            recording_support::compare_custom(("supplemental/electrolyte/" + point).c_str(),
                "action", "native",
                batch.results[index].action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                recording_support::fields({recording_support::field("observable", index), recording_support::field("order", order)}),
                "Relationship", "central-difference-refinement", "charge-balanced electrolyte action",
                "relative", ok, "p4-supplemental", ("electrolyte/" + point).c_str(),
                "action", "T,rho,x", "available", "observable,order");
            CHECK(ok);
        }
    }
    epcsaft::ThermochemistryRecord electrolyte_record;
    electrolyte_record.standard_pressure = 1.0e5;
    for (std::size_t component = 0; component < 3; ++component) {
        electrolyte_record.component_intervals.push_back({
            epcsaft::IdealInterval{250.0, 350.0, true, true, p3_forms()[0].second}
        });
    }
    const std::vector<double> electrolyte_composition = {0.98, 0.01, 0.01};
    const epcsaft::DensityActionDirection electrolyte_total_direction{
        0.1, 1.0, {-0.02, 0.01, 0.01}, {}};
    epcsaft::HelmholtzModel electrolyte_with_ideal(
        properties_support::aqueous_electrolyte(), {}, electrolyte_record);
    const auto electrolyte_total = electrolyte_with_ideal.density_actions(
        300.0, 1000.0, electrolyte_composition,
        {epcsaft::PropertyObservable::TotalEnthalpy}, {electrolyte_total_direction});
    const std::optional<epcsaft::ThermochemistryRecord> electrolyte_record_optional = electrolyte_record;
    const double electrolyte_total_fine = p4_reference(
        electrolyte_with_ideal.parameters(), {}, electrolyte_record_optional, 300.0, 1000.0,
        electrolyte_composition, {electrolyte_total_direction},
        epcsaft::PropertyObservable::TotalEnthalpy, 1.0e-3);
    const double electrolyte_total_middle = p4_reference(
        electrolyte_with_ideal.parameters(), {}, electrolyte_record_optional, 300.0, 1000.0,
        electrolyte_composition, {electrolyte_total_direction},
        epcsaft::PropertyObservable::TotalEnthalpy, 3.0e-3);
    const double electrolyte_total_coarse = p4_reference(
        electrolyte_with_ideal.parameters(), {}, electrolyte_record_optional, 300.0, 1000.0,
        electrolyte_composition, {electrolyte_total_direction},
        epcsaft::PropertyObservable::TotalEnthalpy, 1.0e-2);
    const bool electrolyte_total_ok = electrolyte_total.results.size() == 1
        && electrolyte_total.results[0].status == epcsaft::PropertyActionStatus::Available
        && electrolyte_total.results[0].value.has_value()
        && electrolyte_total.results[0].action.has_value()
        && close(*electrolyte_total.results[0].action, electrolyte_total_fine, 5.0e-6)
        && close(electrolyte_total_fine, electrolyte_total_middle, 1.0e-5);
    recording_support::compare_custom(
        "supplemental/electrolyte-total", "action", "J/mol", electrolyte_total.results.empty()
            ? std::numeric_limits<double>::quiet_NaN()
            : electrolyte_total.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()), electrolyte_total_fine,
        recording_support::fields({recording_support::field("T", 300.0),
            recording_support::field("rho", 1000.0), recording_support::field("x", electrolyte_composition),
            recording_support::field("reference_fine", electrolyte_total_fine),
            recording_support::field("reference_middle", electrolyte_total_middle),
            recording_support::field("reference_coarse", electrolyte_total_coarse)}),
        "Relationship", "central-difference-refinement", "sourced Ideal plus residual total action",
        "relative", electrolyte_total_ok, "p4-supplemental", "electrolyte-total", "action",
        "T,rho,x", "available", "reference_fine,reference_middle,reference_coarse");
    const auto full_one = model.density_actions(
        313.15, 42441.09427799721, composition, observables, {directions[0]});
    const auto selected = model.density_actions(
        313.15, 42441.09427799721, composition,
        {observables[5], observables[13]}, {directions[0]});
    const auto result_matches = [](const epcsaft::PropertyActionResult& left,
                                   const epcsaft::PropertyActionResult& right) {
        if (left.observable != right.observable || left.status != right.status
            || left.energy_order != right.energy_order
            || left.value.has_value() != right.value.has_value()
            || left.action.has_value() != right.action.has_value()) return false;
        return (!left.value.has_value() || close(*left.value, *right.value, 2.0e-12))
            && (!left.action.has_value() || close(*left.action, *right.action, 2.0e-12));
    };
    const bool selected_h = selected.results.size() == 2 && full_one.results.size() == observables.size()
        && result_matches(selected.results[0], full_one.results[5]);
    const bool selected_cp = selected.results.size() == 2 && full_one.results.size() == observables.size()
        && result_matches(selected.results[1], full_one.results[13]);
    CHECK(selected_h && selected_cp);
    recording_support::compare_custom("supplemental/selected-h", "selected_full", "bool", selected_h ? 1.0 : 0.0,
        1.0, "residual_h", "Identity", "selected/full", "selected equals full", "absolute", selected_h,
        "p4-supplemental", "selected-h", "selected_full", "subset", "selected/full");
    recording_support::compare_custom("supplemental/selected-total-cp", "selected_full", "bool", selected_cp ? 1.0 : 0.0,
        1.0, "total_cp", "Identity", "selected/full", "selected equals full", "absolute", selected_cp,
        "p4-supplemental", "selected-total-cp", "selected_full", "subset", "selected/full");
    const auto full_missing = ethane_model.density_actions(
        300.0, 1000.0, {1.0}, observables,
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool failure_match = partial.results.size() == 2 && full_missing.results.size() == observables.size()
        && result_matches(partial.results[0], full_missing.results[0])
        && result_matches(partial.results[1], full_missing.results[5]);
    recording_support::compare_custom(
        "supplemental/selected-failure", "selected_full", "bool",
        failure_match ? 1.0 : 0.0, 1.0, "ideal+residual without source",
        "Identity", "selected/full", "selected preserves typed failure and valid sibling",
        "absolute", failure_match, "p4-supplemental", "selected-failure", "selected_full",
        "subset", "selected/full");
    const auto single_h = model.density_actions(
        313.15, 42441.09427799721, composition, {observables[5]}, {directions[0]});
    const auto single_cp = model.density_actions(
        313.15, 42441.09427799721, composition, {observables[13]}, {directions[0]});
    const bool output_selection_return = single_h.results.size() == 1
        && single_cp.results.size() == 1
        && result_matches(single_h.results[0], full_one.results[5])
        && result_matches(single_cp.results[0], full_one.results[13]);
    recording_support::compare_custom(
        "supplemental/output-selection-return", "selected_full", "bool",
        output_selection_return ? 1.0 : 0.0, 1.0,
        "residual-only then ideal+residual selection at the same state",
        "Identity", "selected/full", "returning to ideal output selection preserves values and actions",
        "absolute", output_selection_return, "p4-supplemental", "output-selection-return",
        "selected_full", "subset", "selected/full");
    epcsaft::HelmholtzModel cold_model(
        model.parameters(), active, p3_binary_record());
    const auto cold_batch = cold_model.density_actions(
        313.15, 42441.09427799721, composition,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity}, {directions[0]});
    const auto warm_batch = cold_model.density_actions(
        313.15, 42441.09427799721, composition,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity}, {directions[0]});
    const std::size_t full_paths = full_one.diagnostics.directional_paths;
    const std::size_t single_paths = single_h.diagnostics.directional_paths
        + single_cp.diagnostics.directional_paths;
    const bool work_reuse = full_one.results.size() == observables.size()
        && single_h.results.size() == 1 && single_cp.results.size() == 1
        && cold_batch.results.size() == 1 && warm_batch.results.size() == 1
        && cold_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && warm_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && full_paths > 0 && single_paths > full_paths
        && full_one.diagnostics.dense_hessian_calls == 0
        && cold_batch.diagnostics.recordings > 0
        && warm_batch.diagnostics.recordings == 0;
    recording_support::compare_custom(
        "supplemental/work-reuse", "directional_paths", "count", static_cast<double>(full_paths),
        static_cast<double>(single_paths),
        recording_support::fields({recording_support::field("full_paths", full_paths),
            recording_support::field("single_sum_paths", single_paths),
            recording_support::field("cold_recordings", cold_batch.diagnostics.recordings),
            recording_support::field("warm_recordings", warm_batch.diagnostics.recordings),
            recording_support::field("full_dense_hessian_calls",
                full_one.diagnostics.dense_hessian_calls)}),
        "Structure", "shared-batch-work", "full request shares overlapping primitive work",
        "absolute", work_reuse, "p4-supplemental", "work-reuse", "directional_paths",
        "full,selected", "paths,recordings,dense_hessian_calls");
    const auto state_values_match = [](
        const epcsaft::State& left, const epcsaft::State& right,
        const std::vector<epcsaft::PropertyObservable>& rows, double tolerance) {
        if (!close(left.molar_density, right.molar_density, tolerance)
            || !close(left.pressure, right.pressure, tolerance)
            || !close(left.pressure_density_derivative, right.pressure_density_derivative, tolerance))
            return false;
        for (const auto observable : rows) {
            const double left_value = p4_state_value(left, observable);
            const double right_value = p4_state_value(right, observable);
            if (std::isnan(left_value) && std::isnan(right_value)) continue;
            if (!std::isfinite(left_value) || !std::isfinite(right_value)
                || !close(left_value, right_value, tolerance)) return false;
        }
        return true;
    };
    const std::vector<double> pure_first = {1.0, 0.0};
    const std::vector<epcsaft::PropertyObservable> freshness_rows = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity,
        epcsaft::PropertyObservable::IdealEntropy};
    for (std::size_t order = 1; order <= 3; ++order) {
        const std::vector<epcsaft::DensityActionDirection> request(order, {0.2, 5.0, {0.0, 0.0}, {}});
        for (std::size_t point = 0; point < 3; ++point) {
            const bool mixed = point == 1;
            const auto& x = mixed ? composition : pure_first;
            const double t = point == 1 ? 320.0 : point == 2 ? 307.15 : 313.15;
            const double rho = point == 1 ? 6000.0 : point == 2 ? 4500.0 : 5000.0;
            const auto retained = model.density_actions(t, rho, x, freshness_rows, request);
            const auto fresh = epcsaft::HelmholtzModel(model.parameters(), active, p3_binary_record())
                .density_actions(t, rho, x, freshness_rows, request);
            for (std::size_t row = 0; row < freshness_rows.size(); ++row) {
                const auto& actual = retained.results[row];
                const auto& expected = fresh.results[row];
                const bool agrees = actual.status == epcsaft::PropertyActionStatus::Available
                    && expected.status == actual.status && actual.value && expected.value
                    && actual.action && expected.action
                    && close(*actual.value, *expected.value, 2.0e-12)
                    && close(*actual.action, *expected.action, 2.0e-12);
                const auto id = "supplemental/action-freshness/" + std::to_string(order)
                    + "/" + std::to_string(point) + "/" + std::to_string(row);
                recording_support::compare_custom(id.c_str(), "action", "native",
                    actual.action.value_or(std::numeric_limits<double>::quiet_NaN()),
                    expected.action.value_or(std::numeric_limits<double>::quiet_NaN()),
                    recording_support::fields({recording_support::field("T", t),
                        recording_support::field("rho", rho), recording_support::field("x", x),
                        recording_support::field("order", order), recording_support::field("row", row),
                        recording_support::field("value", actual.value.value_or(std::numeric_limits<double>::quiet_NaN())),
                        recording_support::field("fresh_value", expected.value.value_or(std::numeric_limits<double>::quiet_NaN()))}),
                    "Identity", "fresh-model-action", "retained action refreshes all state dependencies",
                    "relative", agrees, "p4-supplemental", id.c_str(), "action", "T,rho,x", "fresh-model");
                CHECK(agrees);
            }
        }
    }
    const auto retained_pure = model.state(313.15, 5000.0, pure_first);
    const auto retained_mixed = model.state(320.0, 6000.0, composition);
    const auto retained_pure_again = model.state(313.15, 5000.0, pure_first);
    const auto fresh_pure = epcsaft::HelmholtzModel(
        model.parameters(), active, p3_binary_record()).state(
            313.15, 5000.0, pure_first);
    const auto fresh_mixed = epcsaft::HelmholtzModel(
        model.parameters(), active, p3_binary_record()).state(320.0, 6000.0, composition);
    const auto fresh_pure_again = epcsaft::HelmholtzModel(
        model.parameters(), active, p3_binary_record()).state(
            313.15, 5000.0, pure_first);
    const bool retained_freshness = state_values_match(retained_pure, fresh_pure, observables, 2.0e-12)
        && state_values_match(retained_mixed, fresh_mixed, observables, 2.0e-12)
        && state_values_match(retained_pure_again, fresh_pure_again, observables, 2.0e-12)
        && close(p4_state_value(retained_pure, epcsaft::PropertyObservable::IdealEntropy),
                 p4_state_value(retained_pure_again, epcsaft::PropertyObservable::IdealEntropy), 2.0e-12);
    recording_support::compare_custom(
        "supplemental/retained-freshness", "state", "bool", retained_freshness ? 1.0 : 0.0,
        1.0, recording_support::fields({recording_support::field("pure_x", pure_first),
            recording_support::field("mixed_x", composition), recording_support::field("T_pure", 313.15),
            recording_support::field("T_mixed", 320.0),
            recording_support::field("rho_pure", 5000.0),
            recording_support::field("rho_mixed", 6000.0),
            recording_support::field("pure_entropy", p4_state_value(retained_pure,
                epcsaft::PropertyObservable::IdealEntropy)),
            recording_support::field("mixed_entropy", p4_state_value(retained_mixed,
                epcsaft::PropertyObservable::IdealEntropy)),
            recording_support::field("repeat_entropy", p4_state_value(retained_pure_again,
                epcsaft::PropertyObservable::IdealEntropy))}),
        "Identity", "retained-state-replay", "changing T and x refreshes ideal support and state values",
        "absolute", retained_freshness, "p4-supplemental", "retained-freshness", "state",
        "T,rho,x", "fresh-state", "ideal_entropy,all_properties");
    const std::vector<epcsaft::PropertyObservable> residual_rows = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualInternalEnergy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity,
        epcsaft::PropertyObservable::ResidualIsochoricHeatCapacity};
    bool branch_route_ok = true;
    std::vector<double> branch_densities;
    std::vector<double> branch_pressures;
    for (const auto branch : {epcsaft::DensityBranch::Vapor, epcsaft::DensityBranch::Liquid}) {
        const auto pressure_state = ethane_model.state_at_pressure(
            200.0, 217800.69494692178, {1.0}, branch);
        const auto direct_state = ethane_model.state(
            200.0, pressure_state.molar_density, {1.0});
        const auto fresh_state = epcsaft::HelmholtzModel(properties_support::ethane()).state_at_pressure(
            200.0, 217800.69494692178, {1.0}, branch);
        branch_densities.push_back(pressure_state.molar_density);
        branch_pressures.push_back(pressure_state.pressure);
        branch_route_ok = branch_route_ok
            && pressure_state.density_branch == branch
            && fresh_state.density_branch == branch
            && pressure_state.stable_root_count == fresh_state.stable_root_count
            && state_values_match(pressure_state, fresh_state, residual_rows, 2.0e-12)
            && state_values_match(pressure_state, direct_state, residual_rows, 1.0e-9);
    }
    recording_support::compare_custom(
        "supplemental/branch-route", "state", "bool", branch_route_ok ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 200.0),
            recording_support::field("P", 217800.69494692178),
            recording_support::field("rho", branch_densities),
            recording_support::field("achieved_P", branch_pressures)}),
        "Identity", "resolved-state-route", "same branch pressure and density routes agree",
        "absolute", branch_route_ok, "p4-supplemental", "branch-route", "state",
        "T,P,x,branch", "residual_properties,branch,root_count");
    const auto fixed = model.density_actions(
        313.15, 42441.09427799721, composition, {observables[5]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {0.0, 0.0},
            {{epcsaft::ParameterFamily::AssociationVolume, {0}, 1.0}}}});
    const bool fixed_zero = fixed.results[0].status == epcsaft::PropertyActionStatus::Available
        && fixed.results[0].action.has_value() && std::abs(*fixed.results[0].action) <= 1.0e-12;
    recording_support::compare_custom("supplemental/fixed-slot", "fixed_slot_zero", "bool",
        fixed_zero ? 1.0 : 0.0, 1.0, "association_volume[0]", "Availability", "fixed-slot",
        "known fixed slot is available structural zero", "absolute", fixed_zero,
        "p4-supplemental", "fixed-slot", "fixed_slot_zero", "slot", "fixed/zero");
    bool unknown_rejected = false;
    try {
        (void)model.density_actions(313.15, 42441.09427799721, composition, {observables[5]},
            {epcsaft::DensityActionDirection{0.0, 0.0, {0.0, 0.0},
                {{epcsaft::ParameterFamily::AssociationVolume, {999}, 1.0}}}});
    } catch (const std::invalid_argument&) { unknown_rejected = true; }
    recording_support::compare_custom("supplemental/unknown-slot", "unknown_rejected", "bool",
        unknown_rejected ? 1.0 : 0.0, 1.0, "association_volume[999]", "Availability", "unknown-slot",
        "unknown native slot rejects", "absolute", unknown_rejected,
        "p4-supplemental", "unknown-slot", "unknown_rejected", "slot", "invalid_argument");
    const auto rejects = [](auto&& call) {
        try { call(); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    const auto& valid_direction = directions[0];
    std::size_t malformed_rejections = 0;
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {}, {valid_direction}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0], observables[0]},
        {valid_direction}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]}, {}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        std::vector<epcsaft::DensityActionDirection>(4, valid_direction)); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{std::numeric_limits<double>::infinity(),
            0.0, {0.0, 0.0}, {}}}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {0.1, 0.0}, {}}}); });
    malformed_rejections += rejects([&] { model.density_actions(
        313.15, 42441.09427799721, composition, {observables[0]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {0.0, 0.0},
            {{epcsaft::ParameterFamily::Kij, {2, 1}, 1.0}}}}); });
    epcsaft::HelmholtzModel ionic_validation(properties_support::aqueous_electrolyte());
    malformed_rejections += rejects([&] { ionic_validation.density_actions(
        300.0, 1000.0, {0.98, 0.01, 0.01},
        {epcsaft::PropertyObservable::ResidualEnthalpy},
        {epcsaft::DensityActionDirection{0.0, 0.0, {-0.01, 0.01, 0.0}, {}}}); });
    const bool malformed_ok = malformed_rejections == 9;
    recording_support::compare_custom(
        "supplemental/malformed-requests", "invalid_argument", "count",
        static_cast<double>(malformed_rejections), 9.0,
        "empty,duplicate,order,dimension,nonfinite,balance,slot,charge",
        "Availability", "request-validation", "malformed requests fail at the trust boundary",
        "absolute", malformed_ok, "p4-supplemental", "malformed-requests",
        "invalid_argument", "request", "nine rejected cases");
    CHECK(malformed_ok);
    epcsaft::HelmholtzModel neutral(properties_support::neutral_binary());
    const auto support = neutral.density_actions(300.0, 1000.0, {1.0, 0.0}, {observables[5]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {-0.01, 0.01}, {}}});
    const bool support_change = support.results[0].status == epcsaft::PropertyActionStatus::SupportChange;
    recording_support::compare_custom("supplemental/support-change", "status", "bool",
        support_change ? 1.0 : 0.0, 1.0, "neutral x=(1,0), dx=(-.01,.01)", "Availability",
        "support", "support-changing action unavailable", "absolute", support_change,
        "p4-supplemental", "support-change", "status", "support", "support_change");
    epcsaft::HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    const auto ionic = aqueous.density_actions(300.0, 1000.0, {1.0, 0.0, 0.0}, {observables[5]},
        {epcsaft::DensityActionDirection{0.0, 0.0, {-0.02, 0.01, 0.01}, {}}});
    const bool ionic_insertion = ionic.results[0].status == epcsaft::PropertyActionStatus::IonicInsertion;
    recording_support::compare_custom("supplemental/ionic-insertion", "status", "bool",
        ionic_insertion ? 1.0 : 0.0, 1.0, "aqueous x=(1,0,0), dx=(-.02,.01,.01)", "Availability", "ionic",
        "ionic insertion unavailable", "absolute", ionic_insertion,
        "p4-supplemental", "ionic-insertion", "status", "ionic", "ionic_insertion");
    epcsaft::HelmholtzModel endpoint_model(
        properties_support::ethane(), {}, p3_record(p3_forms()[0].second));
    const auto endpoint_action = endpoint_model.density_actions(250.0, 1000.0, {1.0},
        {epcsaft::PropertyObservable::IdealEnthalpy},
        {epcsaft::DensityActionDirection{0.1, 0.0, {0.0}, {}}});
    const bool thermal_boundary = endpoint_action.results[0].status == epcsaft::PropertyActionStatus::ThermalBoundary;
    recording_support::compare_custom("supplemental/thermal-boundary", "status", "bool",
        thermal_boundary ? 1.0 : 0.0, 1.0, "T=250 endpoint, dT=.1", "Availability", "thermal",
        "thermal boundary action unavailable", "absolute", thermal_boundary,
        "p4-supplemental", "thermal-boundary", "status", "thermal", "thermal_boundary");
    model.reset_counters();
    const auto counter_batch = model.density_actions(313.15, 42441.09427799721, composition,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity},
        {directions[0], directions[1], directions[2]});
    const auto counter_values = model.counters();
    const auto pressure_counter_batch = model.pressure_actions(
        313.15, state.pressure, composition, pressure_root.density_branch,
        {epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity},
        {pressure_directions_all[0], pressure_directions_all[1], pressure_directions_all[2]});
    const bool counter_ok = counter_batch.results.size() == 1
        && pressure_counter_batch.results.size() == 1
        && counter_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && pressure_counter_batch.results[0].status == epcsaft::PropertyActionStatus::Available
        && counter_values.hessian_calls == 0 && counter_values.maximum_energy_order == 5
        && counter_values.directional_paths > 0
        && counter_batch.maximum_energy_order == 5
        && counter_batch.diagnostics.site_rank > 0
        && counter_batch.diagnostics.site_factorizations == 1
        && counter_batch.diagnostics.thermal_evaluations == 0
        && counter_batch.diagnostics.directional_paths > 0
        && counter_batch.diagnostics.directional_paths <= 93
        && counter_batch.diagnostics.forward_sweeps > 0
        && counter_batch.diagnostics.dense_hessian_calls == 0
        && pressure_counter_batch.diagnostics.pressure_factorizations == 1
        && pressure_counter_batch.diagnostics.pressure_solves > 0
        && pressure_counter_batch.diagnostics.directional_paths > 0
        && pressure_counter_batch.diagnostics.directional_paths <= 93
        && pressure_counter_batch.diagnostics.dense_hessian_calls == 0
        && pressure_counter_batch.diagnostics.pressure_density_derivative > 0.0
        && pressure_counter_batch.diagnostics.stable_root_count == pressure_root.stable_root_count
        && pressure_counter_batch.diagnostics.branch == pressure_root.density_branch;
    recording_support::compare_custom("supplemental/counters", "counter_contract", "bool",
        counter_ok ? 1.0 : 0.0, 1.0, "q3 residual cp", "Structure", "counters",
        "directional work reuses one U/P factorization without dense Hessian", "absolute", counter_ok,
        "p4-supplemental", "counters", "counter_contract", "counters", "E5/no-Hessian");
    CHECK(missing_ok && partial_ok && negative_ok && selected_h && selected_cp
        && failure_match && work_reuse && retained_freshness && branch_route_ok
        && electrolyte_total_ok && entropy_ok
        && output_selection_return
        && fixed_zero && unknown_rejected && support_change && ionic_insertion
        && thermal_boundary && counter_ok);
}

void neutral_property_pressure_hessian() {
    auto parameters = properties_support::ethane();
    parameters.segment_count[0] = 2.002;
    parameters.sigma_constant[0] = 3.6184;
    parameters.dispersion_energy_over_k[0] = 208.11;
    parameters.molar_mass_kg_per_mol[0] = 0.044096;
    epcsaft::HelmholtzModel model(parameters);
    constexpr double temperature = 400.0;
    constexpr double density = 60.0;
    const auto inputs = model_inputs(parameters, temperature, 1.0 / density, {1.0});
    const auto actual = model.neutral_property(
        temperature, density, epcsaft::NeutralProperty::Pressure,
        epcsaft::DerivativeOrder::Hessian
    );
    const auto value_after_hessian = model.neutral_property(
        temperature, density, epcsaft::NeutralProperty::Pressure,
        epcsaft::DerivativeOrder::Value
    );
    const auto state = model.state(temperature, density, {1.0});
    recording_support::expect_group(
        recording_support::active_test_id(), "issue10-pressure", "issue10-pressure/", 2,
        "neutral pressure agrees with state", "T,rho"
    );
    recording_support::expect_group(
        recording_support::active_test_id(), "issue10-hessian", "issue10-hessian/", 2,
        "neutral pressure Hessian action", "T,rho"
    );
    recording_support::compare(
        "issue10-pressure/value", "pressure", "Pa", actual.value, state.pressure, inputs,
        "Relationship", "state", "neutral pressure equals state pressure",
        "reference_or_one", 0.0, 1.0e-12,
        close(actual.value, state.pressure, 1.0e-12), "issue10-pressure", "value",
        "P", "T,rho", "P_neutral=P_state", "T,rho"
    );
    recording_support::compare(
        "issue10-pressure/value-after-hessian", "pressure", "Pa",
        value_after_hessian.value, actual.value, inputs,
        "Relationship", "order-transition", "value remains valid after Hessian",
        "reference_or_one", 0.0, 1.0e-12,
        close(value_after_hessian.value, actual.value, 1.0e-12),
        "issue10-pressure", "value-after-hessian", "P", "T,rho",
        "P_value=P_hessian", "T,rho"
    );
    const double delta = 1.0e-4;
    const auto plus = model.neutral_property(
        temperature + 400.0 * delta, density + 100.0 * delta,
        epcsaft::NeutralProperty::Pressure, epcsaft::DerivativeOrder::Gradient
    );
    const auto minus = model.neutral_property(
        temperature - 400.0 * delta, density - 100.0 * delta,
        epcsaft::NeutralProperty::Pressure, epcsaft::DerivativeOrder::Gradient
    );
    constexpr const char* units[] = {"Pa/K", "Pa m3/mol"};
    for (std::size_t row = 0; row < 2; ++row) {
        const double centered = (plus.gradient[row] - minus.gradient[row]) / (2.0 * delta);
        const double expected = actual.hessian[row] * 400.0
            + actual.hessian[row * 2 + 1] * 100.0;
        recording_support::compare(
            (std::string("issue10-hessian/") + std::to_string(row)).c_str(),
            "pressure_hessian_action", units[row], centered, expected,
            recording_support::fields({inputs, recording_support::field("delta", delta),
                                       recording_support::field("direction", row)}),
            "Relationship", "centered_gradient", "pressure Hessian action",
            "reference_or_one", 0.0, 1.0e-6,
            close(centered, expected, 1.0e-6), "issue10-hessian",
            std::to_string(row).c_str(), "H*d", "T,rho", "H*d=finite_difference", "T,rho"
        );
        CHECK(close(centered, expected, 1.0e-6));
    }
}

}  // namespace

// Issue 26: every native permittivity form on the aqueous 2B-water Na+/Cl- liquid.
// A density-dependent eps_r enters pressure, the closure and every pressure-chart chain.
void modular_permittivity_actions() {
    auto ion_suppression = properties_support::aqueous_electrolyte();
    auto component_mixing = ion_suppression;
    component_mixing.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    component_mixing.relative_permittivity = {78.09, 8.0, 8.0};
    auto density_temperature = ion_suppression;  // Yang2024 water d_V, d_T (full EOS 6.2)
    density_temperature.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    density_temperature.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    density_temperature.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    auto dipolar = ion_suppression;  // Rueben2024 PT2 records (full EOS 6.4)
    dipolar.permittivity_form = epcsaft::PermittivityForm::DipolarPerturbation;
    dipolar.dipole_scaling_debye2 = {5.199, 5.407, 2.563};
    dipolar.correlation_integral_parameter = {0.1276, 0.215, 0.138};
    const std::array<std::pair<const char*, epcsaft::Parameters>, 4> forms = {{
        {"ion-suppression", ion_suppression}, {"component-mixing", component_mixing},
        {"density-temperature", density_temperature}, {"dipolar", dipolar}}};
    const std::vector<epcsaft::ActiveParameter> active = {
        {epcsaft::ParameterFamily::Kij, {1, 3}},
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 2}}};
    const std::vector<epcsaft::PropertyObservable> observables = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity};
    const auto parameter = [](double kij, double energy) {
        return std::vector<epcsaft::ParameterActionDirection>{
            {epcsaft::ParameterFamily::Kij, {1, 3}, kij},
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 2}, energy}};
    };
    const std::array<std::vector<epcsaft::ParameterActionDirection>, 3> theta = {
        parameter(0.01, 20.0), parameter(-0.005, -10.0), parameter(0.002, 5.0)};
    const std::array<double, 3> temperature_direction = {1.0, -0.5, 0.2};
    const std::array<std::vector<double>, 3> composition_direction = {{
        {-0.02, 0.01, 0.01}, {0.01, -0.005, -0.005}, {0.004, -0.002, -0.002}}};
    const std::array<double, 3> density_direction = {200.0, -100.0, 50.0};
    const std::array<double, 3> pressure_direction = {3.0e5, -2.0e5, 1.0e5};  // every stencil P > 0
    const std::array<double, 3> tolerances = {5.0e-6, 5.0e-5, 5.0e-4};
    const std::vector<double> x = {0.96, 0.02, 0.02};
    // 2 MPa: at least 1 MPa from the ion-bearing vapor-root fold of every form (scan below).
    constexpr double temperature = 298.15, pressure = 2.0e6;
    const std::array<double, 12> scan = {
        2e3, 1e4, 2e4, 5e4, 9e4, 9.5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6};
    const char* test_id = recording_support::active_test_id();
    recording_support::expect_group(test_id, "issue26-pressure-identity", "issue26/pressure/",
        forms.size(), "residual pressure equals the density derivative of a_res", "form");
    recording_support::expect_group(test_id, "issue26-branch", "issue26/branch/",
        forms.size() * scan.size(), "declared liquid root persists across the vapor fold", "form,P");
    recording_support::expect_group(test_id, "issue26-actions", "issue26/action/",
        forms.size() * 2 * 3 * observables.size(),
        "liquid electrolyte actions in T, rho or P, x, k_ij and association energy", "form,chart,order,observable");
    for (const auto& [name, parameters] : forms) {
        const epcsaft::HelmholtzModel model(parameters, active);
        double previous_liquid = 0.0;
        for (const double scan_pressure : scan) {
            const auto liquid = model.state_at_pressure(
                temperature, scan_pressure, x, epcsaft::DensityBranch::Liquid);
            const auto vapor = model.state_at_pressure(
                temperature, scan_pressure, x, epcsaft::DensityBranch::Vapor);
            const int roots = liquid.stable_root_count.value_or(-1);
            // One stable root: both requests return it unlabelled. Two: liquid is the denser.
            // The liquid root rises monotonically in P, so no request switched branch.
            const bool ok = liquid.molar_density > previous_liquid
                && vapor.stable_root_count == liquid.stable_root_count
                && (roots == 1 ? vapor.molar_density == liquid.molar_density
                        && liquid.density_branch == epcsaft::DensityBranch::Unspecified
                    : roots == 2 && vapor.molar_density < 1.0e-2 * liquid.molar_density
                        && liquid.density_branch == epcsaft::DensityBranch::Liquid)
                && std::abs(liquid.pressure - scan_pressure) <= epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                    * liquid.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature;
            previous_liquid = liquid.molar_density;
            const std::string point = std::string(name) + "/" + recording_support::number(scan_pressure);
            recording_support::compare_custom(("issue26/branch/" + point).c_str(),
                "liquid_molar_density", "mol/m3", liquid.molar_density, vapor.molar_density,
                recording_support::fields({recording_support::field("form", name),
                    recording_support::field("P", scan_pressure),
                    recording_support::field("stable_roots", roots),
                    recording_support::field("epsilon_r_liquid", liquid.bulk_relative_permittivity),
                    recording_support::field("epsilon_r_vapor", vapor.bulk_relative_permittivity),
                    recording_support::field("liquid_pressure_residual", liquid.pressure - scan_pressure)}),
                "Relationship", "declared-branch-scan", "liquid root persists; sole roots are unlabelled",
                "custom", ok, "issue26-branch", point.c_str(), "\\rho_L", "T,P,\\mathbf{x}",
                "branch", "form,P");
            CHECK(ok);
        }
        const auto root = model.state_at_pressure(
            temperature, pressure, x, epcsaft::DensityBranch::Liquid);
        const double rho = root.molar_density;
        // Z - 1 = rho d(a_res)/d(rho) at fixed T, x: a dropped d(eps_r)/d(rho) fails it.
        const double h = 1.0e-5 * rho;
        const double slope = (model.state(temperature, rho + h, x).residual_helmholtz
            - model.state(temperature, rho - h, x).residual_helmholtz) / (2.0 * h);
        const double identity_reference = 1.0 + rho * slope;
        const bool identity_ok = std::abs(root.compressibility_factor - identity_reference)
            <= 1.0e-8 * std::abs(identity_reference - 1.0);
        recording_support::compare_custom(
            (std::string("issue26/pressure/") + name).c_str(), "compressibility_factor", "1",
            root.compressibility_factor, identity_reference,
            recording_support::fields({recording_support::field("form", name),
                recording_support::field("rho", rho),
                recording_support::field("epsilon_r", root.bulk_relative_permittivity),
                recording_support::field("relative_step", 1.0e-5)}),
            "Relationship", "central-difference", "pressure carries every density chain",
            "relative", identity_ok, "issue26-pressure-identity", name, "Z",
            "T,\\rho,\\mathbf{x}", "Z=1+\\rho\\,\\partial a^{res}/\\partial\\rho", "form");
        CHECK(identity_ok);
        for (const bool pressure_chart : {false, true}) {
            for (std::size_t order = 1; order <= 3; ++order) {
                std::vector<epcsaft::DensityActionDirection> density_request;
                std::vector<epcsaft::PressureActionDirection> pressure_request;
                for (std::size_t k = 0; k < order; ++k) {
                    density_request.push_back({temperature_direction[k], density_direction[k],
                        composition_direction[k], theta[k]});
                    pressure_request.push_back({temperature_direction[k], pressure_direction[k],
                        composition_direction[k], theta[k]});
                }
                const auto batch = pressure_chart
                    ? model.pressure_actions(temperature, pressure, x,
                        epcsaft::DensityBranch::Liquid, observables, pressure_request)
                    : model.density_actions(temperature, rho, x, observables, density_request);
                const auto action_reference = pressure_chart
                    ? property_action_reference(parameters, active, std::nullopt, temperature, pressure, x,
                        observables, pressure_request, epcsaft::DensityBranch::Liquid)
                    : property_action_reference(parameters, active, std::nullopt, temperature, rho, x,
                        observables, density_request);
                for (std::size_t index = 0; index < observables.size(); ++index) {
                    const double fine = action_reference.value[index];
                    const double middle = action_reference.middle[index];
                    const double allowance = tolerances[order - 1] * action_reference.scale[index];
                    const auto& result = batch.results[index];
                    const bool ok = result.status == epcsaft::PropertyActionStatus::Available
                        && result.action.has_value()
                        && std::abs(*result.action - fine) <= allowance
                        && action_reference.uncertainty[index] <= allowance;
                    const std::string point = std::string(name) + (pressure_chart ? "/P/" : "/rho/")
                        + std::to_string(order) + "/" + std::to_string(index);
                    recording_support::compare_custom(("issue26/action/" + point).c_str(),
                        "action", "native",
                        result.action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                        recording_support::fields({recording_support::field("form", name),
                            recording_support::field("chart", pressure_chart ? "T,P,x,theta" : "T,rho,x,theta"),
                            recording_support::field("order", order),
                            recording_support::field("observable", index),
                            recording_support::field("status", static_cast<int>(result.status)),
                            recording_support::field("reference_fine", fine),
                            recording_support::field("reference_middle", middle),
                            recording_support::field("tolerance", tolerances[order - 1]),
                            recording_support::field("term_scale", action_reference.scale[index]),
                            recording_support::field("reference_uncertainty", action_reference.uncertainty[index]),
                            recording_support::field("allowance", allowance)}),
                        "Relationship", "central-difference-refinement",
                        "liquid electrolyte action retains every permittivity chain",
                        "tolerance_times_chart_term_scale", ok, "issue26-actions", point.c_str(), "action",
                        "T,\\rho|P,\\mathbf{x},\\theta", "available", "form,chart,order,observable");
                    CHECK(ok);
                }
            }
        }
    }
}

// Issue 26: implicit-association E3/E4 energy actions against a centered difference of the
// next-lower exact order (u^T H v for E3, E3 for E4) along the last direction, Richardson h, h/2, h/4.
struct EnergyCase final {
    const char* name;
    epcsaft::Parameters parameters;
    std::vector<epcsaft::ActiveParameter> active;
    double temperature, molar_density;
    std::vector<double> composition;
    std::vector<std::vector<double>> directions;  // width 2 + C + A: T, V, n, theta
};

double shifted_energy(const EnergyCase& item, const std::vector<std::vector<double>>& lower,
                      const std::vector<double>& along, double t) {
    const std::size_t components = item.composition.size();
    std::vector<double> shift(item.active.size()), amounts(components);
    for (std::size_t a = 0; a < shift.size(); ++a) shift[a] = t * along[2 + components + a];
    for (std::size_t i = 0; i < components; ++i) amounts[i] = item.composition[i] + t * along[2 + i];
    const epcsaft::HelmholtzModel model(
        p4_shifted_parameters(item.parameters, item.active, shift), item.active);
    const double temperature = item.temperature + t * along[0];
    const double volume = 1.0 / item.molar_density + t * along[1];
    if (lower.size() == 3) return model.residual_energy_action(temperature, volume, amounts, lower).action.value();
    const auto derivatives = model.evaluate(temperature, volume, amounts, epcsaft::DerivativeOrder::Hessian);
    const std::size_t width = along.size();
    CHECK(derivatives.gradient.size() == width);
    double sum = 0.0;
    for (std::size_t i = 0; i < width; ++i)
        for (std::size_t j = 0; j < width; ++j)
            sum += lower[0][i] * derivatives.hessian[i * width + j] * lower[1][j];
    return sum;
}

void energy_action_ladders() {
    auto density_temperature = properties_support::aqueous_electrolyte();
    density_temperature.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    density_temperature.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    density_temperature.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    auto dipolar = properties_support::aqueous_electrolyte();
    dipolar.permittivity_form = epcsaft::PermittivityForm::DipolarPerturbation;
    dipolar.dipole_scaling_debye2 = {5.199, 5.407, 2.563};
    dipolar.correlation_integral_parameter = {0.1276, 0.215, 0.138};
    auto component_mixing = properties_support::aqueous_electrolyte();
    component_mixing.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    component_mixing.relative_permittivity = {78.09, 8.0, 8.0};
    const std::vector<epcsaft::ActiveParameter> aqueous_active = {
        {epcsaft::ParameterFamily::Kij, {1, 3}},
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 2}}};
    const double aqueous_volume = 1.0 / 56000.0;
    // Charge-neutral amount directions; V, T, k_ij and association-energy components nonzero.
    const std::vector<std::vector<double>> aqueous_directions = {
        {0.5, 1.0e-2 * aqueous_volume, -0.02, 0.01, 0.01, 0.01, 10.0},
        {-0.3, -0.5e-2 * aqueous_volume, 0.01, -0.005, -0.005, -0.004, -6.0},
        {0.2, 0.3e-2 * aqueous_volume, 0.004, -0.002, -0.002, 0.002, 4.0},
        {-0.1, 0.2e-2 * aqueous_volume, -0.006, 0.003, 0.003, -0.003, -3.0}};
    const std::vector<double> aqueous_x = {0.96, 0.02, 0.02};
    const double binary_volume = 1.0 / 42441.09427799721;
    const std::vector<EnergyCase> cases = {
        {"ion-suppression", properties_support::aqueous_electrolyte(), aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"component-mixing", component_mixing, aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"density-temperature", density_temperature, aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"dipolar", dipolar, aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"cross-associating-binary", properties_support::associating_binary(),
         {{epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}},
          {epcsaft::ParameterFamily::AssociationVolume, {1, 4}}},
         313.15, 42441.09427799721, {0.25, 0.75}, {
            {0.5, 1.0e-2 * binary_volume, -0.01, 0.01, 10.0, 1.0e-3},
            {-0.3, -0.5e-2 * binary_volume, 0.005, -0.004, -6.0, -5.0e-4},
            {0.2, 0.3e-2 * binary_volume, 0.003, -0.002, 4.0, 3.0e-4},
            {-0.1, 0.2e-2 * binary_volume, -0.004, 0.003, -3.0, 2.0e-4}}},
    };
    const char* test_id = recording_support::active_test_id();
    recording_support::expect_group(test_id, "issue26-energy", "issue26/energy/", cases.size() * 3,
        "E3 and E4 residual energy actions match the difference of the next lower exact order",
        "case,action");
    for (const auto& item : cases) {
        const auto& d = item.directions;
        const std::array<std::pair<const char*, std::vector<std::vector<double>>>, 3> requests = {{
            {"E3", {d[0], d[1], d[2]}}, {"E3-repeated", {d[0], d[0], d[2]}}, {"E4", {d[0], d[1], d[2], d[3]}}}};
        const epcsaft::HelmholtzModel model(item.parameters, item.active);
        for (const auto& [label, request] : requests) {
            const auto action = model.residual_energy_action(
                item.temperature, 1.0 / item.molar_density, item.composition, request);
            const std::vector<std::vector<double>> lower(request.begin(), request.end() - 1);
            const auto difference = [&](double h) {
                return (shifted_energy(item, lower, request.back(), h)
                    - shifted_energy(item, lower, request.back(), -h)) / (2.0 * h);
            };
            // Two Richardson extrapolations (h, h/2) and (h/2, h/4): their difference is an
            // asymptotic, conservative estimate of the finer one's error (about 15x it for an h^4
            // remainder), not a bound. It must resolve 1e-7 relative, a decade below the 1e-6
            // agreement criterion that this falsifier uses to expose a missing chain.
            const double d1 = difference(1.0e-2), d2 = difference(5.0e-3), d3 = difference(2.5e-3);
            const double r1 = (4.0 * d2 - d1) / 3.0, reference = (4.0 * d3 - d2) / 3.0;
            const double scale = std::abs(reference);
            const double uncertainty = std::abs(reference - r1);
            const bool resolved = uncertainty <= 1.0e-7 * scale;
            const bool ok = action.status == epcsaft::PropertyActionStatus::Available
                && action.action.has_value() && resolved
                && std::abs(*action.action - reference) <= 1.0e-6 * scale;
            const std::string point = std::string(item.name) + "/" + label;
            recording_support::compare_custom(("issue26/energy/" + point).c_str(), "energy_action",
                "extensive residual scalar", action.action.value_or(std::numeric_limits<double>::quiet_NaN()),
                reference, recording_support::fields({recording_support::field("case", item.name),
                    recording_support::field("action", label),
                    recording_support::field("status", static_cast<int>(action.status)),
                    recording_support::field("difference_h", d1),
                    recording_support::field("difference_h_over_2", d2),
                    recording_support::field("difference_h_over_4", d3),
                    recording_support::field("reference_uncertainty", uncertainty),
                    recording_support::field("association_residual", action.diagnostics.association_residual)}),
                "Relationship", "centered-difference-richardson",
                "implicit association and permittivity chains complete through fourth order",
                "relative", ok, "issue26-energy", point.c_str(), label, "T,V,\\mathbf{n},\\theta",
                "D^kA=\\partial_t D^{k-1}A", "case,action");
            CHECK(ok);
        }
    }
}

// Copied anchor compositions from 2026-born-unavailable-states/inputs, reordered
// CO2/MEA/water -> properties_support's water/MEA/CO2. EqIDs own the radial work.
void born_parameter_columns() {
    const std::array<std::vector<double>, 2> anchors = {{{0.8863734508600574, 0.07827110880066358, 3.8480681348098447e-08, 0.01818111244399814, 0.015766643324345936, 0.0003959814177106597, 0.0010068230563202352, 1.3460805650294977e-11, 4.841602761876763e-06},
        {0.8448866085505489, 0.004900981383765356, 0.0002933071241251476, 0.07818889597821514, 0.029096076931579815, 0.03617530531138919, 0.006458693289371563, 2.1372510655078036e-09, 1.2929375400277957e-07}}};
    for (std::size_t case_index = 0; case_index < anchors.size(); ++case_index)
    for (const double c : {0.0, 1.0}) for (const double epsilon : {2.0, 8.0}) {
        auto p = properties_support::nine_species();
        p.born_c_shell = 1.0; p.born_c_dielectric = c;
        p.ionic_region_relative_permittivity = epsilon;
        const auto& x = anchors[case_index];
        const epcsaft::HelmholtzModel model(p, {{epcsaft::ParameterFamily::BornCDielectric, {0}},
            {epcsaft::ParameterFamily::IonicRegionRelativePermittivity, {0}}});
        const double t = 313.15, rho = 45000.0;
        const auto state = model.state(t, rho, x);
        const auto derivatives = model.evaluate(t, 1.0 / rho, x, epcsaft::DerivativeOrder::Gradient);
        double solvation = 0.0, total = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) { solvation += x[i] * p.solvation_factor[i]; total += x[i]; }
        solvation /= total;
        double shell_sum = 0.0, correction = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) if (p.charge_number[i]) {
            const double z = p.charge_number[i], d = p.born_diameter_angstrom[i] * epcsaft::ANGSTROM_M;
            const double outer = d * (1.0 + (solvation - 1.0) / std::abs(z));
            shell_sum += x[i] * z * z / outer;
            correction += x[i] * z * z * (1.0 / d - 1.0 / outer);
        }
        const double prefactor = -epcsaft::ELEMENTARY_CHARGE_C * epcsaft::ELEMENTARY_CHARGE_C
            / (4.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M * epcsaft::BOLTZMANN_J_PER_K * t);
        const std::array<double, 3> reference = {
            prefactor * ((1.0 - 1.0 / state.bulk_relative_permittivity) * shell_sum + c * (1.0 - 1.0 / epsilon) * correction),
            prefactor * (1.0 - 1.0 / epsilon) * correction,
            prefactor * c * correction / (epsilon * epsilon)};
        const std::array<double, 3> actual = {state.born, derivatives.gradient[2 + x.size()], derivatives.gradient[3 + x.size()]};
        for (std::size_t column = 0; column < actual.size(); ++column) {
            const auto id = std::to_string(case_index + 1) + "/" + std::to_string(c) + "/" + std::to_string(epsilon) + "/" + std::to_string(column);
            const double atol = reference[column] == 0.0 ? 1e-12 : 0.0;
            const bool ok = std::abs(actual[column] - reference[column]) <= atol + 1e-12 * std::abs(reference[column]);
            CHECK(recording_support::compare(("born-parameter/" + id).c_str(), "Born radial work and first columns", "dimensionless",
                actual[column], reference[column], recording_support::fields({recording_support::field("epsilon", epsilon), recording_support::field("c", c)}),
                "Relationship", "born_radial_charging_work", "frozen-state Born value, c coefficient and dq/depsilon chain", "reference_scale",
                atol, 1e-12, ok, "born-parameter", id.c_str(), "a_Born", "c,epsilon", "radial work", "case,c,epsilon,column"));
        }
    }
}

namespace {
using PO = epcsaft::PropertyObservable;
using PD = epcsaft::PressureActionDirection;
using PF = epcsaft::ParameterFamily;
constexpr double reference_nan = std::numeric_limits<double>::quiet_NaN();
double factorial(unsigned n) { double v = 1; for (unsigned k=2; k<=n; ++k) v *= k; return v; }
double binomial(unsigned n, unsigned k) { return factorial(n)/(factorial(k)*factorial(n-k)); }

// Independent P4 generator: C_T=(h_res+RT)/(RT^2).
double pure_reference_thermal(const epcsaft::HelmholtzModel& model, double T, double P,
    unsigned a, unsigned b, bool theta, double rho) {
    const double ps=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).pressure_density_derivative
        *model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).molar_density;
    const auto h = [&](unsigned thermal) {
        std::vector<PD> d;
        for (unsigned j=0; j<thermal; ++j) { PD q; q.composition={0.0}; q.temperature=T; d.push_back(q); }
        for (unsigned j=0; j<b; ++j) { PD q; q.composition={0.0}; q.pressure=ps; d.push_back(q); }
        if (theta) { PD q; q.composition={0.0}; q.parameter={{PF::DispersionEnergyOverK,{0},191.42}}; d.push_back(q); }
        if (d.empty()) return model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).residual_enthalpy;
        // Heat capacity first/second actions independently supply h thermal derivatives.
        const auto observable = thermal ? PO::ResidualIsobaricHeatCapacity : PO::ResidualEnthalpy;
        if (thermal) d.erase(d.begin());
        if (d.empty()) return *model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).residual_isobaric_heat_capacity;
        return model.pressure_actions(T,P,{1},epcsaft::DensityBranch::Liquid,{observable},d,rho)
            .results.front().action.value_or(reference_nan) / (std::pow(T,int(thermal?thermal-1:0))*std::pow(ps,int(b))*(theta?191.42:1));
    };
    double value=0;
    for (unsigned j=0; j<a; ++j)
        value += binomial(a-1,j) * (j%2 ? -1:1) * factorial(j+1)
            * h(a-1-j)/std::pow(T,int(j+2))/epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    if (!b && !theta) value += ((a-1)%2 ? -1:1)*factorial(a-1)/std::pow(T,int(a));
    return value;
}

// Independent pressure-volume inverse through E4; no reference Taylor path.
double pure_reference_pressure(const epcsaft::HelmholtzModel& model, double T, double P,
    unsigned b, bool theta, double rho) {
    rho=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).molar_density;
    const double V=1/rho, rt=epcsaft::GAS_CONSTANT_J_PER_MOL_K*T;
    const auto e=model.evaluate(T,V,{1},epcsaft::DerivativeOrder::Hessian);
    const auto w=e.gradient.size();
    std::vector<double> v(w,0), k(w,0); v[1]=1; k[w-1]=1;
    const auto energy=[&](std::vector<std::vector<double>> d) {
        return model.residual_energy_action(T,V,{1},d).action.value_or(reference_nan);
    };
    const double pv=rt*(-1/(V*V)-e.hessian[w+1]);
    const double pvv=rt*(2/(V*V*V)-energy({v,v,v}));
    const double pvvv=rt*(-6/std::pow(V,4)-energy({v,v,v,v}));
    double value=V;
    if (theta) {
        const double vk=rt*e.hessian[w+w-1]/pv;
        const double pvk=-rt*energy({v,v,k})+pvv*vk;
        const double pvvk=-rt*energy({v,v,v,k})+pvvv*vk;
        if (b==1) value=vk;
        if (b==2) value=-pvk/(pv*pv);
        if (b==3) value=-pvvk/std::pow(pv,3)+3*pvv*pvk/std::pow(pv,4);
    } else {
        if (b==2) value=1/pv;
        if (b==3) value=-pvv/std::pow(pv,3);
        if (b==4) value=3*pvv*pvv/std::pow(pv,5)-pvvv/std::pow(pv,4);
    }
    return -value/rt;
}

void reference_f8() {
    epcsaft::HelmholtzModel model(properties_support::ethane(),{{PF::DispersionEnergyOverK,{0}}});
    const double T=260, R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    for (double P : {3e6,5e6}) {
        const auto state=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid);
        for (unsigned c=0; c<=1; ++c) for (unsigned degree=1; degree<=4; ++degree)
            for (unsigned b=0; b<=degree-c; ++b) {
                const unsigned a=degree-b-c;
                const double scale=R*std::pow(T,int(a))*std::pow(P,int(b))*(c?191.42:1);
                const auto exact=epcsaft::reference_partial(model,T,P,{1},{-1},a,b,c?std::optional<std::size_t>(0):std::nullopt,state.molar_density);
                const double oracle=(c && !a && !b) ? -model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,state.molar_density).log_fugacity_parameter_derivative[0].value() : a ? pure_reference_thermal(model,T,P,a,b,c,state.molar_density)
                    : pure_reference_pressure(model,T,P,b,c,state.molar_density);
                const double x=scale*exact.action.value_or(reference_nan), y=scale*oracle;
                const double target=1e-7+1e-9*std::abs(y);
                const bool ok=exact.status==epcsaft::PropertyActionStatus::Available && std::abs(x-y)<=target;
                const std::string id="F8/"+std::to_string(int(P))+"/"+std::to_string(a)+std::to_string(b)+std::to_string(c);
                recording_support::compare_custom(id.c_str(),"scaled reference partial","J/mol/K",x,y,
                    recording_support::fields({recording_support::field("T",T),recording_support::field("P_source",P),
                    recording_support::field("a",double(a)),recording_support::field("b",double(b)),recording_support::field("theta",double(c)),
                    recording_support::field("target",target),recording_support::field("error",std::abs(x-y)),recording_support::field("error_over_target",std::abs(x-y)/target),
                    recording_support::field("site_residual",exact.diagnostics.association_residual),
                    recording_support::field("pressure_series_residual",exact.diagnostics.pressure_residual)}),
                    "Relationship","pure-P4-identities","finite source-pressure reference derivatives","absolute+relative",ok,"F8");
                CHECK(ok);
            }
    }
}
}

namespace {
// E2/E3 Hessian-chain oracle at a literal fixed composition. No limiting selection.
double reference_lower(const epcsaft::HelmholtzModel& m, double T, double P,
    const std::vector<double>& x, unsigned a, unsigned b, double anchor) {
    const auto s=m.state_at_pressure(T,P,x,epcsaft::DensityBranch::Liquid,anchor);
    const double V=1/s.molar_density, R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const auto e=m.evaluate(T,V,x,epcsaft::DerivativeOrder::Hessian);
    const auto w=e.gradient.size();
    std::vector<double> nu(w,0), ev(w,0); nu[2]=-1;nu[3]=nu[4]=1;ev[1]=1;
    const double pv=R*T*(-1/(V*V)-e.hessian[w+1]);
    const double pt=R*(1/V-e.gradient[1])-R*T*e.hessian[w];
    double cv=-1/V;
    for (std::size_t i=2; i<5; ++i) cv+=nu[i]*e.hessian[i*w+1];
    std::vector<std::vector<double>> path;
    for (unsigned j=0; j<a+b; ++j) {
        std::vector<double> d(w,0);
        if (j<a) { d[0]=1;d[1]=-pt/pv; } else d[1]=1/pv;
        path.push_back(d);
    }
    if (a+b==1) {
        double val=cv*path[0][1];
        for (std::size_t i=2; i<5; ++i) val+=nu[i]*e.hessian[i*w]*path[0][0];
        return val;
    }
    const auto& u=path[0];const auto& v=path[1];
    const auto mixed=[&](const std::vector<double>& q) {
        // Size both noncomposition motions before half-sign polarization.
        const auto norm=[&](const auto& d) { return std::max(std::abs(d[0])/T,std::abs(d[1])/V); };
        const double alpha=norm(v)/norm(u);
        auto plus=v,minus=v;
        for(std::size_t i=0;i<w;++i){plus[i]+=alpha*u[i];minus[i]-=alpha*u[i];}
        const auto p=m.residual_energy_action(T,V,x,{q,plus,plus});
        const auto n=m.residual_energy_action(T,V,x,{q,minus,minus});
        return (p.action.value_or(reference_nan)-n.action.value_or(reference_nan))/(4*alpha);
    };
    const auto hv=[&](const auto& d){double v=0;for(std::size_t i=0;i<w;++i)v+=e.hessian[w+i]*d[i];return v;};
    const double p2=R*(u[0]*(-v[1]/(V*V)-hv(v))+v[0]*(-u[1]/(V*V)-hv(u)))
        +R*T*(2*u[1]*v[1]/std::pow(V,3)-mixed(ev));
    return mixed(nu)+u[1]*v[1]/(V*V)-cv*p2/pv;
}

void reference_f9() {
    const double T=298.15,R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    auto params=properties_support::aqueous_electrolyte();
    epcsaft::HelmholtzModel model(params,{{PF::DispersionEnergyOverK,{0}}});
    for(double P:{1e5,2e5}) for(double ratio:{1.0,2.0}) {
        std::optional<double> anchor;
        for(double molality:{1e-10,1e-12,1e-14}) {
            const double salt=molality*params.molar_mass_kg_per_mol[0]*ratio;
            const std::vector<double> x{1/(1+2*salt),salt/(1+2*salt),salt/(1+2*salt)};
            const auto state=model.state_at_pressure(T,P,x,epcsaft::DensityBranch::Liquid,anchor);
            anchor=state.molar_density;
            for(unsigned c=0;c<=1;++c) for(unsigned d=2;d<=4;++d) for(unsigned b=0;b<=d-c;++b) {
                const unsigned a=d-b-c;
                if(d==2 && !(c && b))continue;
                const auto exact=epcsaft::reference_partial(model,T,P,x,{-1,1,1},a,b,
                    c?std::optional<std::size_t>(0):std::nullopt,*anchor);
                const double scale=R*std::pow(T,int(a))*std::pow(P,int(b))*(c?353.95:1);
                const double value=scale*exact.action.value_or(reference_nan);
                const auto lower=[&](double step) {
                    auto p=params;
                    double t=T,pressure=P;
                    unsigned la=a,lb=b;
                    if(c)p.dispersion_energy_over_k[0]+=step;
                    else if(a){t+=step;--la;}else{pressure+=step;--lb;}
                    epcsaft::HelmholtzModel perturbed(p,{{PF::DispersionEnergyOverK,{0}}});
                    if(la+lb<=2)return reference_lower(perturbed,t,pressure,x,la,lb,*anchor);
                    return epcsaft::reference_partial(perturbed,t,pressure,x,{-1,1,1},la,lb,std::nullopt,*anchor).action.value_or(reference_nan);
                };
                double h=c?0.02:a?0.2:2e4;
                bool accepted=false;
                double coarse=scale*(lower(h)-lower(-h))/(2*h);
                for(unsigned halving=0;halving<=6;++halving) {
                    const double plus=lower(h/2),minus=lower(-h/2);
                    const double fine=scale*(plus-minus)/h;
                    const double rich=(4*fine-coarse)/3,refine=std::abs(fine-coarse)/3;
                    const double sigma=std::abs(value); // One reference-chart basis monomial.
                    const double target=1e-7+1e-6*std::max(std::abs(rich),sigma);
                    const bool resolved=refine<=target/4;
                    accepted=resolved&&std::abs(value-rich)<=target;
                    const std::string id="F9/partial/"+std::to_string(int(P))+"/"+std::to_string(int(ratio))+"/"+std::to_string(int(-std::log10(molality)))
                        +"/"+std::to_string(a)+std::to_string(b)+std::to_string(c);
                    const auto inputs=recording_support::fields({recording_support::field("T",T),recording_support::field("P_source",P),
                        recording_support::field("molality",molality),recording_support::field("ratio",ratio),
                        recording_support::field("a",double(a)),recording_support::field("b",double(b)),recording_support::field("theta",double(c)),
                        recording_support::field("h",h),recording_support::field("coarse",coarse),recording_support::field("fine",fine),
                        recording_support::field("refinement",refine),recording_support::field("target",target),recording_support::field("sigma",sigma),
                        recording_support::field("error",std::abs(value-rich)),recording_support::field("error_over_target",std::abs(value-rich)/target),recording_support::field("resolved",resolved?1:0),
                        recording_support::field("site_residual",exact.diagnostics.association_residual),
                        recording_support::field("density",state.molar_density),recording_support::field("pressure_root_residual",state.pressure-P),
                        recording_support::field("pressure_density_slope",state.pressure_density_derivative),
                        recording_support::field("site_pivot_ratio",exact.diagnostics.site_pivot_ratio),
                        recording_support::field("rounding_over_step",64*std::numeric_limits<double>::epsilon()*scale*(std::abs(plus)+std::abs(minus))/h),
                        recording_support::field("pressure_series_residual",exact.diagnostics.pressure_residual)});
                    if (!resolved && halving < 6)
                        recording_support::observe(("F9/refinement/"+id.substr(11)+"/"+std::to_string(halving)).c_str(),"scaled reference partial","J/mol/K",
                            value,rich,inputs,"Numerical","fixed-composition-Richardson","coarse step requires refinement");
                    else recording_support::compare_custom(id.c_str(),"scaled reference partial","J/mol/K",value,rich,inputs,
                        "Relationship","fixed-composition-Richardson","finite terminal ionic derivatives","absolute+relative",accepted,"F9");
                    if(resolved)break;
                    h/=2;coarse=fine;
                }
                CHECK(accepted);
            }
        }
    }
}
}
namespace {
void reference_limits() {
    const double T=298.15,P=1e5,R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const auto params=properties_support::aqueous_electrolyte();
    epcsaft::HelmholtzModel model(params,{{PF::DispersionEnergyOverK,{0}}});
    for(bool fixed:{false,true}) {
        const double src=fixed?2e5:P;
        epcsaft::NeutralReferenceDeclaration declaration{0,{0,1,1},1};
        epcsaft::ReactionReference reference{fixed?std::optional<double>(src):std::nullopt,
            epcsaft::ReferenceSourceBasis::CommonMolalityInfiniteDilution};
        const auto jet=epcsaft::reaction_references(model,T,P,declaration,{{-1,1,1}},{reference},true,std::nullopt,2,true,4)[0];
        CHECK(jet.status==epcsaft::NeutralReferenceStatus::Available);
        auto doubled=declaration;doubled.solute_molality_ratios={0,2,2};
        const auto jet2=epcsaft::reaction_references(model,T,P,doubled,{{-1,1,1}},{reference},true,std::nullopt,2,true,4)[0];
        CHECK(jet2.status==epcsaft::NeutralReferenceStatus::Available);
        for(std::size_t i=0;i<jet.derivatives.size();++i) {
            const auto& column=jet.derivatives[i];const auto& second=jet2.derivatives[i];
            CHECK(column.value && second.value);
            const double scale=R*std::pow(T,int(column.temperature_order))*std::pow(src,int(column.pressure_order))
                *(column.parameter?353.95:1);
            const double change=scale*std::abs(*column.value-*second.value);
            const std::string id="F9/doubled/"+std::to_string(fixed)+"/"+std::to_string(i);
            recording_support::compare_custom(id.c_str(),"doubled terminal concentration change","J/mol/K",change,0,
                recording_support::fields({recording_support::field("target",0.01)}),"Relationship","same-neutral-curve",
                "doubled terminal concentration, no different ionic-ratio claim","absolute",change<=0.01,"F9-limit");
            CHECK(change<=0.01);
            if(fixed && column.pressure_order) {
                const std::string zero="F8/system-P-zero/"+std::to_string(i);
                recording_support::compare_custom(zero.c_str(),"fixed-source system-pressure derivative","raw SI",*column.value,0,
                    "{}","Relationship","analytic-zero","fixed-source system-P terms vanish exactly","absolute",*column.value==0,"F8");
                CHECK(*column.value==0);
            }
        }

        for(const auto& column:jet.derivatives) {
            CHECK(column.value);
            for(unsigned k=0;k<2;++k) {
                const std::string id="F9/limit/"+std::to_string(fixed)+"/"+std::to_string(column.temperature_order)
                    +std::to_string(column.pressure_order)+std::to_string(column.parameter.has_value())+"/"+std::to_string(k);
                const bool ok=column.scaled_changes[k]<=0.01;
                recording_support::compare_custom(id.c_str(),"scaled consecutive reference change","J/mol/K",
                    column.scaled_changes[k],0,recording_support::fields({recording_support::field("target",0.01),
                    recording_support::field("a",double(column.temperature_order)),recording_support::field("b",double(column.pressure_order)),
                    recording_support::field("theta",column.parameter?1:0),recording_support::field("source_pressure",src)}),
                    "Relationship","terminal-ladder","independent consumed-column admission","absolute",ok,"F9-limit");CHECK(ok);
            }
        }
        // Existing equilibrium driver supplies only the case-center composition.
        // This independent verification completion never supplies runtime calorics.
        auto problem=equilibrium_support::flash(3);
        problem.temperature=epcsaft::equilibrium::FreeValue::fixed(T);
        problem.pressure=epcsaft::equilibrium::FreeValue::fixed(P);
        problem.feed.values={1,1e-6,1e-6};problem.source_reference=declaration;
        problem.reactions={{"H2O=Na++Cl-",{-1,1,1},1e-8,std::nullopt,std::nullopt,reference}};
        const auto solved=epcsaft::equilibrium::solve_equilibrium(model,problem);
        CHECK(solved.success);
        const auto& x=solved.mole_fractions;
        const auto center=model.state_at_pressure(T,P,x,epcsaft::DensityBranch::Liquid);
        PD dt,dp;dt.temperature=1;dp.pressure=1;dt.composition=dp.composition={0,0,0};
        const auto thermal=model.pressure_actions(T,P,x,epcsaft::DensityBranch::Liquid,{}, {dt},center.molar_density);
        const auto pressure=model.pressure_actions(T,P,x,epcsaft::DensityBranch::Liquid,{}, {dp},center.molar_density);
        const double v=1/center.molar_density, vt=-thermal.density_action[0]*v*v,vp=-pressure.density_action[0]*v*v;
        const double y=(x[1]+x[2])/2,anchor_weight=x[0]+y;
        std::array<double,3> cp{},cv{},u{};
        for(unsigned k=0;k<3;++k) {
            const double molality=k==0?1e-10:k==1?1e-12:1e-14;
            const double salt=molality*params.molar_mass_kg_per_mol[0];
            const std::vector<double> base{1/(1+2*salt),salt/(1+2*salt),salt/(1+2*salt)};
            const auto partial=[&](unsigned a,unsigned b) {
                return *epcsaft::reference_partial(model,T,src,base,{-1,1,1},a,b,std::nullopt).action;
            };
            const double ct=partial(1,0),ctt=partial(2,0);
            const double cpref=fixed?0:partial(0,1),ctpref=fixed?0:partial(1,1),cppref=fixed?0:partial(0,2);
            const double vg=v-y*R*T*cpref, vgt=vt-y*R*(cpref+T*ctpref),vgp=vp-y*R*T*cppref;
            const double h=anchor_weight*(-241826-R*T)+y*R*T*T*ct+R*T+center.residual_enthalpy;
            cp[k]=anchor_weight*(33.59-R)+y*R*(2*T*ct+T*T*ctt)+R+*center.residual_isobaric_heat_capacity;
            cv[k]=cp[k]+T*vgt*vgt/vgp;u[k]=h-P*vg;
            const std::string id="F9/completion/"+std::to_string(fixed)+"/"+std::to_string(k);
            recording_support::observe(id.c_str(),"completed frozen cv","J/mol/K",cv[k],std::nullopt,
                recording_support::fields({recording_support::field("cp",cp[k]),recording_support::field("u_G",u[k]),
                recording_support::field("v_G",vg),recording_support::field("v_G_T",vgt),recording_support::field("v_G_P",vgp),
                recording_support::field("quotient_magnitude",std::abs(T*vgt*vgt/vgp)),
                recording_support::field("x_water",x[0]),recording_support::field("x_Na",x[1]),recording_support::field("x_Cl",x[2]),
                recording_support::field("equilibrium_residual",solved.residual_norm()),
                recording_support::field("molality",molality),recording_support::field("source_pressure",src)}),
                "Numerical","independent-anchored-completion","case-center quotient conditioning");
            CHECK(vg>0 && vgp<0 && std::isfinite(cv[k]));
        }
        for(unsigned k=0;k<2;++k)for(bool is_cv:{false,true}) {
            const double change=std::abs((is_cv?cv:cp)[k+1]-(is_cv?cv:cp)[k]);
            const std::string id="F9/completed-change/"+std::to_string(fixed)+"/"+std::to_string(is_cv)+"/"+std::to_string(k);
            const bool ok=change<=0.01;
            recording_support::compare_custom(id.c_str(),is_cv?"completed cv change":"completed cp change","J/mol/K",change,0,
                recording_support::fields({recording_support::field("u_G_change",std::abs(u[k+1]-u[k])),recording_support::field("target",0.01)}),
                "Relationship","completed-ladder","completed heat-capacity convergence","absolute",ok,"F9-limit");CHECK(ok);
        }
    }
}
}
int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.properties.born-parameter", "Born_parameter_columns", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "born-parameter/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "BornParameter", "born-parameter/*:24"}, born_parameter_columns},
        {{"native.properties.reference-f8", "Reference_derivative_pure_identities", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "F8/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "F8", "F8/*:48"}, reference_f8, "F8"},
        {{"native.properties.reference-f9", "Reference_derivative_ionic_differences", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "F9/partial/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "F9", "F9/partial/*:204"}, reference_f9, "F9"},
        {{"native.properties.reference-limits", "Reference_derivative_terminal_changes", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "F9/limit/*,F9/completed-change/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "ReferenceLimits"}, reference_limits, "ReferenceLimits"},
        {{"native.properties.p1", "Residual_properties_and_caloric_ladders", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P1", "cv-ladder/*,cp-ladder/*,dilute-state/*,dilute-limit/*,conditioning/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "P1", "cv-ladder/*:20,cp-ladder/*:20,dilute-state/*:4,dilute-limit/*:5,conditioning/*:4"}, residual_properties_and_caloric_ladders, "P1"},
        {{"native.properties.f1", "Residual_internal_energy_identity_ladder", "Numerical", "Properties", "Automatic_fast", "Numerical", "Residual", "P1", "a_res,pressure,compressibility_factor,residual_entropy,residual_enthalpy,temperature-ladder/*,internal-energy-identity,enthalpy-internal-energy-identity", "engine/tests/native/properties_checks.cpp", "properties_checks", "F1", "a_res:1,pressure:1,compressibility_factor:1,residual_entropy:1,residual_enthalpy:1,temperature-ladder/*:4,internal-energy-identity:1,enthalpy-internal-energy-identity:1"}, residual_internal_energy_identity_ladder, "F1"},
        {{"native.properties.m2", "Derivatives_and_limits", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P1", "value-state/*,residual-sum/*,gradient-density/*,hessian-symmetry/*,euler/*,internal-energy/*,internal-energy-entropy/*,chemical-potential/*,low-density/*,low-density-z/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "M2", "value-state/*:4,residual-sum/*:4,gradient-density/*:4,hessian-symmetry/*:178,euler/*:4,internal-energy/*:4,internal-energy-entropy/*:4,chemical-potential/*:4,low-density/*:4,low-density-z/*:4"}, derivatives_and_limits},
        {{"native.properties.m3", "Finite_difference_gradient", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "gradient/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "M3", "gradient/*:5"}, finite_difference_gradient},
        {{"native.properties.m4", "Hessian_reference", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "associating-gradient/*,associating-hessian/*,nine-gradient/*,nine-hessian/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "M4", "associating-gradient/*:4,associating-hessian/*:16,nine-gradient/*:12,nine-hessian/*:144"}, hessian_reference},
        {{"native.properties.m5", "Hessian_and_pressure_slope", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "hessian/*,pressure-density-slope/*,volume-second-derivative/*,projected/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "M5", "hessian/*:178,pressure-density-slope/*:4,volume-second-derivative/*:4,projected/*:514"}, hessian_fd_and_pressure_slope},
        {{"native.properties.m6", "Same_sign_pair_inert", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "C3", "excluded-inert,all-pairs-active", "engine/tests/native/properties_checks.cpp", "properties_checks", "M6", "excluded-inert:1,all-pairs-active:1"}, same_sign_pair_is_inert},
        {{"native.properties.m7", "Parameter_derivative_and_helpers", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "parameter-gradient,packing-coefficient,phase-packing", "engine/tests/native/properties_checks.cpp", "properties_checks", "M7", "parameter-gradient:1,packing-coefficient:1,phase-packing:1"}, parameter_derivative_and_helpers},
        {{"native.properties.c2", "Core_C2_parameter_coordinates", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "CORE2", "c2-scalar-value/*,c2-falsifier,c2-association-gradient,c2-scaling", "engine/tests/native/properties_checks.cpp", "properties_checks", "C2", "c2-scalar-value/*:3,c2-falsifier:1,c2-association-gradient:1,c2-scaling:1"}, core_c2_parameter_coordinates},
        {{"native.properties.m8", "Density_and_reference", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P2", "vapor-pressure/*,liquid-pressure/*,density-order/*,vapor-volume/*,liquid-volume/*,supercritical-pressure,supercritical-density,anchored-liquid,anchored-vapor,associating-pressure,nine-packing,reference-finite,activity-finite,reference-error,activity-scale", "engine/tests/native/properties_checks.cpp", "properties_checks", "M8", "vapor-pressure/*:2,liquid-pressure/*:2,density-order/*:2,vapor-volume/*:2,liquid-volume/*:2,supercritical-pressure:1,supercritical-density:1,anchored-liquid:1,anchored-vapor:1,associating-pressure:1,nine-packing:1,reference-finite:1,activity-finite:1,reference-error:1,activity-scale:1"}, density_and_reference},
        {{"native.properties.p2", "Properties_P2_Bounded_Roots_and_Continuation", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P2", "support-residual/*,support-hessian/*,support-floor/*,association-action/*,root-residual/*,root-convergence/*,root-metadata/*,root-coalescing/*,single-metadata/*,density-selection/*,anchor/*,anchor-invalid/*,fixed-pressure/*,loop-vapor/*,loop-liquid/*,loop-ambiguous/*,loop-action/*,fold-resolution/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "P2", "support-residual/*:4,support-hessian/*:6,support-floor/*:2,association-action/*:3,root-residual/*:3,root-convergence/*:2,root-metadata/*:2,root-coalescing/*:1,single-metadata/*:1,density-selection/*:1,anchor/*:3,anchor-invalid/*:1,fixed-pressure/*:9,loop-vapor/*:13,loop-liquid/*:13,loop-ambiguous/*:13,loop-action/*:2,fold-resolution/*:1"}, properties_p2, "P2"},
        {{"native.properties.p3", "Properties_P3_Thermochemistry_and_References", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P3", "p3-form-value/*,p3-form-derivative/*,p3-form-derivative-refinement/*,p3-total/*,p3-boundary/*,p3-endpoint/*,p3-mixing/*,p3-availability/*,p3-reference/*,p3-reaction/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "P3", "p3-form-value/*:25,p3-form-derivative/*:40,p3-form-derivative-refinement/*:10,p3-total/*:30,p3-boundary/*:6,p3-endpoint/*:50,p3-mixing/*:5,p3-availability/*:5,p3-reference/*:13,p3-reaction/*:35"}, properties_thermochemistry_p3, "P3"},
        {{"native.properties.p4", "Properties_P4_Directional_Actions", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "direct/*,pressure/*,supplemental/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "P4", "direct/*:90,pressure/*:90,supplemental/*:149"}, properties_p4, "P4"},
        {{"native.properties.m9", "Cache_invariance", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "value/*,gradient/*,hessian/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "M9", "value/*:2,gradient/*:10,hessian/*:50"}, cache_invariance},
        {{"native.properties.issue26", "Modular_permittivity_actions", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "issue26/pressure/*,issue26/branch/*,issue26/action/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "Issue26", "issue26/pressure/*:4,issue26/branch/*:48,issue26/action/*:72"}, modular_permittivity_actions},
        {{"native.properties.electrolyte-coordinates", "Electrolyte_parameter_actions", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "electrolyte-coordinate/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "ElectrolyteCoordinates"}, electrolyte_parameter_actions},
        {{"native.properties.issue26e", "Energy_action_ladders", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "issue26/energy/*", "engine/tests/native/properties_checks.cpp", "properties_checks", "Issue26E", "issue26/energy/*:15"}, energy_action_ladders},
        {{"native.properties.issue10", "Neutral_property_pressure_hessian", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P4", "pressure,value,gradient,hessian", "engine/tests/native/properties_checks.cpp", "properties_checks", "Issue10", "issue10-pressure/:2,issue10-hessian/:2"}, neutral_property_pressure_hessian},
    });
}
