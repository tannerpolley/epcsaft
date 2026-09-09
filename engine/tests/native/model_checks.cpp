#include "model/helmholtz_model.hpp"
#include "model_support.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
    epcsaft::HelmholtzModel model(model_support::ethane());
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

void derivatives_and_limits() {
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> cases = {
        {model_support::neutral_binary(), {0.4, 0.6}},
        {model_support::associating_binary(), {0.1385, 0.8615}},
        {model_support::aqueous_electrolyte(), {0.98, 0.01, 0.01}},
        {model_support::nine_species(), {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0, 0.005, 0}},
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
                * state_identity.residual_chemical_potential_over_rt[component];
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
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
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
        model_support::associating_binary(),
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
            auto plus = model_support::associating_binary();
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
    const auto associating_parameters = model_support::associating_binary();
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
    const auto associating_reference = model_support::reference_tensor(
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

    const auto nine_parameters = model_support::nine_species();
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
    const auto nine_reference = model_support::reference_tensor(
        nine_parameters, 300.0, 1.0 / 1.2e4, nine_amounts,
        nine_result.site_fractions,
        epcsaft::ActiveParameter{
            epcsaft::ParameterFamily::DispersionEnergyOverK, {0}
        }
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
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> cases = {
        {model_support::neutral_binary(), {0.4, 0.6}},
        {model_support::associating_binary(), {0.1385, 0.8615}},
        {model_support::aqueous_electrolyte(), {0.98, 0.01, 0.01}},
         {model_support::nine_species(),
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
    auto excluded = model_support::nine_species();
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
    epcsaft::Parameters parameters = model_support::neutral_binary();
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
        auto bad = model_support::aqueous_electrolyte();
        bad.ionic_dispersion = {};
        (void)epcsaft::HelmholtzModel(std::move(bad));
    } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    threw = false;
    try {
        auto bad = model_support::neutral_binary();
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
        auto bad = model_support::associating_binary();
        bad.site_component[0] = 2;
        (void)epcsaft::HelmholtzModel(std::move(bad));
    } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
}

void density_and_reference() {
    epcsaft::HelmholtzModel ethane(model_support::ethane());
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
        cold_vapor.molar_density
    );
    const auto anchored_vapor = ethane.state_at_pressure(
        280.0, 1e6, {1.0}, epcsaft::DensityBranch::Vapor,
        cold_liquid.molar_density
    );
    const auto anchored_inputs = recording_support::fields({
        model_inputs(ethane.parameters(), 280.0, 1.0 / cold_liquid.molar_density, {1.0}),
        recording_support::field("P", 1e6),
        recording_support::field("cold_liquid_density", cold_liquid.molar_density),
        recording_support::field("cold_vapor_density", cold_vapor.molar_density),
        recording_support::field("liquid_anchor", cold_vapor.molar_density),
        recording_support::field("vapor_anchor", cold_liquid.molar_density),
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
    const auto associating = epcsaft::HelmholtzModel(model_support::associating_binary());
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
    const auto nine = epcsaft::HelmholtzModel(model_support::nine_species());
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
    const auto aqueous = model_support::aqueous_electrolyte();
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

void cache_invariance() {
    epcsaft::HelmholtzModel model(model_support::aqueous_electrolyte());
    recording_support::expect_group(
        recording_support::active_test_id(), "m9-value", "value", 1,
        "cached value invariance", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "m9-gradient", "gradient/", 5,
        "cached gradient invariance", "coordinate");
    recording_support::expect_group(
        recording_support::active_test_id(), "m9-hessian", "hessian/", 25,
        "cached Hessian invariance", "row,column");
    const auto a = model.evaluate(298.15, 1.0 / 1e3, {0.98, 0.01, 0.01},
                                  epcsaft::DerivativeOrder::Hessian);
    model.evaluate(298.15, 1.0 / 1e3, {1.0, 0.0, 0.0}, epcsaft::DerivativeOrder::Value);
    model.evaluate(298.15, 1.0 / 1e3, {1.0 - 2e-12, 1e-12, 1e-12}, epcsaft::DerivativeOrder::Value);
    const auto again = model.evaluate(298.15, 1.0 / 1e3, {0.98, 0.01, 0.01},
                                      epcsaft::DerivativeOrder::Hessian);
    const auto inputs = model_inputs(model.parameters(), 298.15,
                                     1.0 / 1e3, {0.98, 0.01, 0.01});
    recording_support::compare(
        "value", "residual_helmholtz", "1", a.value, again.value, inputs,
        "Relationship", "cache-invariance", "repeated evaluation preserves value",
        "reference_or_one", 0.0, 1e-14, close(a.value, again.value, 1e-14),
        "m9-value", "state", "a^{res}", "T,V,\\mathbf{n}", "a_{cached}=a", "state");
    CHECK(close(a.value, again.value, 1e-14));
    CHECK(a.gradient.size() == again.gradient.size());
    for (std::size_t i = 0; i < a.gradient.size(); ++i) {
        recording_support::compare(
            (std::string("gradient/") + std::to_string(i)).c_str(), "gradient", "1",
            a.gradient[i], again.gradient[i], inputs, "Relationship", "cache-invariance",
            "repeated evaluation preserves gradient", "reference_or_one", 0.0, 1e-14,
            close(a.gradient[i], again.gradient[i], 1e-14), "m9-gradient",
            std::to_string(i).c_str(), "\\nabla a", "T,V,\\mathbf{n}",
            "\\nabla a_{cached}=\\nabla a", "coordinate");
        CHECK(close(a.gradient[i], again.gradient[i], 1e-14));
    }
    CHECK(a.hessian.size() == again.hessian.size());
    for (std::size_t i = 0; i < a.hessian.size(); ++i) {
        const std::size_t width = a.gradient.size();
        const auto point = std::to_string(i / width) + "-" + std::to_string(i % width);
        recording_support::compare(
            (std::string("hessian/") + point).c_str(), "hessian", "1", a.hessian[i],
            again.hessian[i], inputs, "Relationship", "cache-invariance",
            "repeated evaluation preserves Hessian", "reference_or_one", 0.0, 1e-14,
            close(a.hessian[i], again.hessian[i], 1e-14), "m9-hessian", point.c_str(),
            "H^{res}", "T,V,\\mathbf{n}", "H_{cached}=H", "row,col");
        CHECK(close(a.hessian[i], again.hessian[i], 1e-14));
    }
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.model.f1", "Residual_internal_energy_identity_ladder", "Numerical", "Properties", "Automatic_fast", "Numerical", "Residual", "M1", "a_res,pressure,compressibility_factor,residual_entropy,residual_enthalpy,temperature-ladder/*,internal-energy-identity,enthalpy-internal-energy-identity", "engine/tests/native/model_checks.cpp", "model_checks", "F1", "a_res:1,pressure:1,compressibility_factor:1,residual_entropy:1,residual_enthalpy:1,temperature-ladder/*:4,internal-energy-identity:1,enthalpy-internal-energy-identity:1"}, residual_internal_energy_identity_ladder, "F1"},
        {{"native.model.m2", "Derivatives_and_limits", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M1", "value-state/*,residual-sum/*,gradient-density/*,hessian-symmetry/*,euler/*,internal-energy/*,internal-energy-entropy/*,chemical-potential/*,low-density/*,low-density-z/*", "engine/tests/native/model_checks.cpp", "model_checks", "M2", "value-state/*:4,residual-sum/*:4,gradient-density/*:4,hessian-symmetry/*:178,euler/*:4,internal-energy/*:4,internal-energy-entropy/*:4,chemical-potential/*:4,low-density/*:4,low-density-z/*:4"}, derivatives_and_limits},
        {{"native.model.m3", "Finite_difference_gradient", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M4", "gradient/*", "engine/tests/native/model_checks.cpp", "model_checks", "M3", "gradient/*:5"}, finite_difference_gradient},
        {{"native.model.m4", "Hessian_reference", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M4", "associating-gradient/*,associating-hessian/*,nine-gradient/*,nine-hessian/*", "engine/tests/native/model_checks.cpp", "model_checks", "M4", "associating-gradient/*:4,associating-hessian/*:16,nine-gradient/*:12,nine-hessian/*:144"}, hessian_reference},
        {{"native.model.m5", "Hessian_and_pressure_slope", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M4", "hessian/*,pressure-density-slope/*,volume-second-derivative/*", "engine/tests/native/model_checks.cpp", "model_checks", "M5", "hessian/*:178,pressure-density-slope/*:4,volume-second-derivative/*:4"}, hessian_fd_and_pressure_slope},
        {{"native.model.m6", "Same_sign_pair_inert", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P3", "excluded-inert,all-pairs-active", "engine/tests/native/model_checks.cpp", "model_checks", "M6", "excluded-inert:1,all-pairs-active:1"}, same_sign_pair_is_inert},
        {{"native.model.m7", "Parameter_derivative_and_helpers", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M4", "parameter-gradient,packing-coefficient,phase-packing", "engine/tests/native/model_checks.cpp", "model_checks", "M7", "parameter-gradient:1,packing-coefficient:1,phase-packing:1"}, parameter_derivative_and_helpers},
        {{"native.model.m8", "Density_and_reference", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M2", "vapor-pressure/*,liquid-pressure/*,density-order/*,vapor-volume/*,liquid-volume/*,supercritical-pressure,supercritical-density,anchored-liquid,anchored-vapor,associating-pressure,nine-packing,reference-finite,activity-finite,reference-error,activity-scale", "engine/tests/native/model_checks.cpp", "model_checks", "M8", "vapor-pressure/*:2,liquid-pressure/*:2,density-order/*:2,vapor-volume/*:2,liquid-volume/*:2,supercritical-pressure:1,supercritical-density:1,anchored-liquid:1,anchored-vapor:1,associating-pressure:1,nine-packing:1,reference-finite:1,activity-finite:1,reference-error:1,activity-scale:1"}, density_and_reference},
        {{"native.model.m9", "Cache_invariance", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "M4", "value,gradient/*,hessian/*", "engine/tests/native/model_checks.cpp", "model_checks", "M9", "gradient/*:5,hessian/*:25"}, cache_invariance},
    });
}
