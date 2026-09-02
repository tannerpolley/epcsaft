#include "model/helmholtz_model.hpp"
#include "model_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#define CHECK(condition) do { if (!(condition)) { std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition); std::abort(); } } while (false)

namespace {

bool close(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

void derivatives_and_limits() {
    const std::vector<std::pair<epcsaft::Parameters, std::vector<double>>> cases = {
        {model_support::neutral_binary(), {0.4, 0.6}},
        {model_support::associating_binary(), {0.1385, 0.8615}},
        {model_support::aqueous_electrolyte(), {0.98, 0.01, 0.01}},
        {model_support::nine_species(), {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0, 0.005, 0}},
    };
    for (const auto& item : cases) {
        epcsaft::HelmholtzModel model(item.first);
        const double temperature = 300.0;
        const double volume = 1.0 / 1.2e4;
        const auto value = model.evaluate(temperature, volume, item.second,
                                          epcsaft::DerivativeOrder::Value);
        const auto state = model.state(temperature, 1.2e4, item.second);
        CHECK(close(value.value, state.residual_helmholtz, 1e-14));
        CHECK(close(state.residual_helmholtz,
                    state.hard_chain + state.dispersion + state.association
                    + state.debye_huckel + state.born, 1e-12));
        const auto gradient = model.evaluate(temperature, volume, item.second,
                                             epcsaft::DerivativeOrder::Gradient);
        CHECK(std::isfinite(gradient.gradient[1]));
        const auto hessian = model.evaluate(temperature, volume, item.second,
                                            epcsaft::DerivativeOrder::Hessian);
        for (std::size_t row = 0; row < gradient.gradient.size(); ++row) {
            for (std::size_t column = 0; column < gradient.gradient.size(); ++column) {
                const double left = hessian.hessian[row * gradient.gradient.size() + column];
                const double right = hessian.hessian[column * gradient.gradient.size() + row];
                CHECK(close(left, right, 1e-12));
            }
        }
        double euler = volume * gradient.gradient[1];
        for (std::size_t component = 0; component < item.second.size(); ++component) {
            euler += item.second[component] * gradient.gradient[2 + component];
        }
        CHECK(close(euler, value.value, 1e-12));
        const auto state_identity = model.state(temperature, 1.2e4, item.second);
        double chemical_identity = 0.0;
        for (std::size_t component = 0; component < item.second.size(); ++component) {
            chemical_identity += item.second[component]
                * state_identity.residual_chemical_potential_over_rt[component];
        }
        CHECK(close(chemical_identity,
                    state_identity.residual_helmholtz
                    + state_identity.compressibility_factor - 1.0, 1e-12));
    }
    epcsaft::HelmholtzModel model(model_support::neutral_binary());
    for (double density : {1e1, 1e0, 1e-1, 1e-2}) {
        const auto state = model.state(300.0, density, {0.4, 0.6});
        CHECK(std::abs(state.residual_helmholtz) < 1e-2);
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
    const auto actual = model.evaluate(temperature, volume, amounts,
                                       epcsaft::DerivativeOrder::Gradient);
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
        CHECK(close(actual.gradient[coordinate], finite_difference, 1e-6));
    }
}

void hessian_reference() {
    const auto associating_parameters = model_support::associating_binary();
    const std::vector<double> associating_amounts = {0.1385, 0.8615};
    epcsaft::HelmholtzModel associating(associating_parameters);
    const auto associating_result = associating.evaluate(
        376.88, 1.0e-4, associating_amounts, epcsaft::DerivativeOrder::Hessian
    );
    const auto associating_reference = model_support::reference_tensor(
        associating_parameters, 376.88, 1.0e-4, associating_amounts,
        associating_result.site_fractions
    );
    CHECK(associating_result.gradient.size() == associating_reference.gradient.size());
    for (std::size_t coordinate = 0;
         coordinate < associating_reference.gradient.size(); ++coordinate) {
        CHECK(close(associating_result.gradient[coordinate],
                    associating_reference.gradient[coordinate], 1e-10));
    }
    for (std::size_t entry = 0; entry < associating_reference.hessian.size(); ++entry) {
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
    for (std::size_t coordinate = 0;
         coordinate < nine_reference.gradient.size(); ++coordinate) {
        CHECK(close(nine_result.gradient[coordinate],
                    nine_reference.gradient[coordinate], 1e-10));
    }
    for (std::size_t entry = 0; entry < nine_reference.hessian.size(); ++entry) {
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
    for (const auto& item : cases) {
        epcsaft::HelmholtzModel model(item.first);
        const double temperature = 300.0;
        const double volume = 1.0 / 1.2e4;
        const auto actual = model.evaluate(
            temperature, volume, item.second, epcsaft::DerivativeOrder::Hessian
        );
        const std::size_t width = actual.gradient.size();
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
                CHECK(close(actual.hessian[row * width + coordinate],
                            (plus_gradient[row] - minus_gradient[row]) / (2.0 * step),
                            1e-6));
            }
        }
        const auto slope_state = model.state(temperature, 1.2e4, item.second);
        const auto plus_state = model.state(temperature, 1.2e4 * (1.0 + 1e-5), item.second);
        const auto minus_state = model.state(temperature, 1.2e4 * (1.0 - 1e-5), item.second);
        CHECK(close(slope_state.pressure_density_derivative,
                    (plus_state.pressure - minus_state.pressure) / (2.0 * 1.2e4 * 1e-5),
                    1e-6));
    }
}

void parameter_derivative_and_helpers() {
    epcsaft::Parameters parameters = model_support::neutral_binary();
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
    CHECK(close(derivative.gradient.back(), (up - down) / (2.0 * step), 1e-6));
    const double coefficient = model.packing_fraction_coefficient(300.0, {0.4, 0.6});
    CHECK(close(coefficient * 1.2e4,
                epcsaft::packing_fraction_from_geometry(
                    1.2e4, {0.4, 0.6}, parameters,
                    epcsaft::segment_diameter(
                        300.0, parameters, epcsaft::segment_sigma(300.0, parameters))),
                1e-14));
    CHECK(close(model.phase_packing_fraction(300.0, {0.4, 0.6}, 1.0 / 1.2e4),
                coefficient * 1.2e4, 1e-14));
    bool threw = false;
    try { auto bad = parameters; bad.k_ij.pop_back(); (void)epcsaft::HelmholtzModel(std::move(bad)); }
    catch (const std::invalid_argument&) { threw = true; }
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
    for (const auto& row : std::array<std::tuple<double, double, double, double>, 2>{
             std::tuple{200.0, 217800.69494692178, 0.007241100995,
                        0.03007 / 523.97698},
             std::tuple{240.0, 969152.4807171321, 0.001745731329,
                        0.03007 / 465.30887}}) {
        const double temperature = std::get<0>(row);
        const double pressure = std::get<1>(row);
        const auto vapor = ethane.state_at_pressure(
            temperature, pressure, {1.0}, epcsaft::DensityBranch::Vapor
        );
        const auto liquid = ethane.state_at_pressure(
            temperature, pressure, {1.0}, epcsaft::DensityBranch::Liquid
        );
        CHECK(vapor.molar_density < liquid.molar_density);
        CHECK(close(vapor.pressure, pressure, 1e-9));
        CHECK(close(liquid.pressure, pressure, 1e-9));
        CHECK(close(1.0 / vapor.molar_density, std::get<2>(row), 0.05));
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
    CHECK(close(vapor.pressure, 1e6, 1e-9));
    CHECK(close(liquid.molar_density, vapor.molar_density, 1e-10));
    const auto nine = epcsaft::HelmholtzModel(model_support::nine_species());
    const auto nine_state = nine.state_at_pressure(
        313.15, 1e5,
        {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
        epcsaft::DensityBranch::Liquid
    );
    CHECK(nine.phase_packing_fraction(
              313.15,
              {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
              1.0 / nine_state.molar_density
          ) < epcsaft::policy::CLOSE_PACKING_FRACTION.value);
    const auto aqueous = model_support::aqueous_electrolyte();
    epcsaft::HelmholtzModel model(aqueous);
    const auto reference = epcsaft::infinite_dilution_reference(
        model, 298.15, 1e5, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}
    );
    const auto activity = epcsaft::electrolyte_activity(
        model, 298.15, 1e5, {1.0, 0.0, 0.0}, {1, 2}, {1, 1}, 1e-8, reference
    );
    CHECK(std::isfinite(reference.convergence_error));
    CHECK(std::isfinite(activity.mean_ionic_log_molality_scale));
    const auto pure = model.state_at_pressure(
        298.15, 1e5, {1.0, 0.0, 0.0}, epcsaft::DensityBranch::Liquid
    );
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
    CHECK(close(reference.convergence_error,
                gamma_coefficient * (std::sqrt(1e-10) - std::sqrt(1e-12)),
                1e-3));
    CHECK(close(activity.mean_ionic_log_molality_scale,
                -gamma_coefficient * (std::sqrt(1e-8) - std::sqrt(1e-12)),
                1e-2));
}

void cache_invariance() {
    epcsaft::HelmholtzModel model(model_support::aqueous_electrolyte());
    const auto a = model.evaluate(298.15, 1.0 / 1e3, {0.98, 0.01, 0.01},
                                  epcsaft::DerivativeOrder::Hessian);
    model.evaluate(298.15, 1.0 / 1e3, {1.0, 0.0, 0.0}, epcsaft::DerivativeOrder::Value);
    model.evaluate(298.15, 1.0 / 1e3, {1.0 - 2e-12, 1e-12, 1e-12}, epcsaft::DerivativeOrder::Value);
    const auto again = model.evaluate(298.15, 1.0 / 1e3, {0.98, 0.01, 0.01},
                                      epcsaft::DerivativeOrder::Hessian);
    CHECK(close(a.value, again.value, 1e-14));
    CHECK(a.gradient.size() == again.gradient.size());
    for (std::size_t i = 0; i < a.gradient.size(); ++i) CHECK(close(a.gradient[i], again.gradient[i], 1e-14));
    CHECK(a.hessian.size() == again.hessian.size());
    for (std::size_t i = 0; i < a.hessian.size(); ++i) CHECK(close(a.hessian[i], again.hessian[i], 1e-14));
}

}  // namespace

int main() {
    derivatives_and_limits();
    finite_difference_gradient();
    hessian_reference();
    hessian_fd_and_pressure_slope();
    parameter_derivative_and_helpers();
    density_and_reference();
    cache_invariance();
    std::printf("model checks passed\n");
}
