#pragma once

#include "properties/helmholtz_model.hpp"

#include <Eigen/Dense>
#include <cppad/cppad.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace properties_support {

inline epcsaft::Parameters parameters(std::size_t components) {
    epcsaft::Parameters p;
    p.segment_count.resize(components);
    p.sigma_constant.resize(components);
    p.dispersion_energy_over_k.resize(components);
    p.molar_mass_kg_per_mol.resize(components);
    p.charge_number.assign(components, 0);
    p.packing_diameter_angstrom.assign(components, 0.0);
    p.debye_huckel_diameter_angstrom.assign(components, 0.0);
    p.born_diameter_angstrom.assign(components, 0.0);
    p.solvation_factor.assign(components, 1.0);
    p.relative_permittivity.assign(components, 0.0);
    p.permittivity_reference_temperature_k.assign(components, 298.15);
    p.ion_suppression_coefficient.assign(components, 0.0);
    p.dielectric_volume_m3_per_mol.assign(components, 0.0);
    p.dielectric_temperature_k.assign(components, 0.0);
    p.dipole_scaling_debye2.assign(components, 0.0);
    p.polarizability_scaling_angstrom3.assign(components, 0.0);
    p.correlation_integral_parameter.assign(components, 0.0);
    p.k_ij.assign(components * components, 0.0);
    p.k_ij_reciprocal_temperature_slope.assign(components * components, 0.0);
    p.k_ij_reference_temperature_k.assign(components * components, 0.0);
    p.l_ij.assign(components * components, 0.0);
    return p;
}

inline epcsaft::Parameters ethane() {
    auto p = parameters(1);
    p.segment_count[0] = 1.6069;
    p.sigma_constant[0] = 3.5206;
    p.dispersion_energy_over_k[0] = 191.42;
    p.molar_mass_kg_per_mol[0] = 0.03007;
    return p;
}

inline epcsaft::Parameters neutral_binary() {
    auto p = parameters(2);
    p.segment_count = {1.0, 1.6069};
    p.sigma_constant = {3.7039, 3.5206};
    p.dispersion_energy_over_k = {150.03, 191.42};
    p.molar_mass_kg_per_mol = {0.016043, 0.03007};
    return p;
}

inline epcsaft::Parameters neutral_ternary() {
    auto p = parameters(3);
    p.segment_count = {1.0, 1.6069, 2.0020};
    p.sigma_constant = {3.7039, 3.5206, 3.6184};
    p.dispersion_energy_over_k = {150.03, 191.42, 208.11};
    p.molar_mass_kg_per_mol = {0.016043, 0.03007, 0.044097};
    return p;
}

inline epcsaft::Parameters associating_binary() {
    auto p = parameters(2);
    p.segment_count = {4.5354, 2.1945};
    p.sigma_constant = {2.6019, 2.229};
    p.dispersion_energy_over_k = {204.0438, 141.66};
    p.molar_mass_kg_per_mol = {0.06108, 0.01801528};
    p.k_ij[1] = p.k_ij[2] = -0.052;
    p.site_component = {0, 0, 1, 1};
    p.site_multiplicity = {2, 1, 2, 2};
    p.association_active.assign(16, 0);
    p.association_energy_over_k.assign(16, 0.0);
    p.association_volume.assign(16, 0.0);
    const auto set_edge = [&](std::size_t a, std::size_t b, double energy, double volume) {
        p.association_active[a * 4 + b] = 1;
        p.association_active[b * 4 + a] = 1;
        p.association_energy_over_k[a * 4 + b] = energy;
        p.association_energy_over_k[b * 4 + a] = energy;
        p.association_volume[a * 4 + b] = volume;
        p.association_volume[b * 4 + a] = volume;
    };
    set_edge(0, 1, 2383.4744, 0.118488);
    set_edge(0, 3, 2093.8222, 0.1540467874409226);
    set_edge(1, 2, 2093.8222, 0.1540467874409226);
    set_edge(2, 3, 1804.17, 0.2039);
    return p;
}

inline epcsaft::Parameters aqueous_electrolyte() {
    auto p = parameters(3);
    p.segment_count = {1.2047, 1.0, 1.0};
    p.sigma_constant = {2.7927, 2.8232, 2.756};
    p.sigma_term_component = {0, 0};
    p.sigma_term_amplitude = {10.11, -1.417};
    p.sigma_term_exponent = {-0.01775, -0.01146};
    p.dispersion_energy_over_k = {353.95, 230.0, 170.0};
    p.molar_mass_kg_per_mol = {0.0180153, 0.02298977, 0.0354527};
    p.charge_number = {0, 1, -1};
    p.packing_diameter_angstrom = {0.0, 2.484416, 2.42528};
    p.debye_huckel_diameter_angstrom = {0.0, 2.484416, 2.42528};
    p.born_diameter_angstrom = {0.0, 3.445, 4.1};
    p.solvation_factor = {1.5, 1.0, 1.0};
    p.relative_permittivity = {78.09, 0.0, 0.0};
    p.ion_suppression_coefficient = {0.0, 7.01, 7.01};
    p.electrolyte_form = epcsaft::ElectrolyteForm::DebyeHuckelBorn;
    p.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
    p.permittivity_form = epcsaft::PermittivityForm::IonSuppression;
    p.ionic_region_relative_permittivity = 8.0;
    p.site_component = {0, 0};
    p.site_multiplicity = {1, 1};
    p.association_active = {0, 1, 1, 0};
    p.association_energy_over_k = {0.0, 2425.7, 2425.7, 0.0};
    p.association_volume = {0.0, 0.04509, 0.04509, 0.0};
    p.k_ij[1] = p.k_ij[3] = -0.3;
    p.k_ij[5] = p.k_ij[7] = 0.8;
    return p;
}

inline epcsaft::Parameters nine_species() {
    auto p = parameters(9);
    p.segment_count = {1.2047, 3.0353, 2.079, 1, 1, 1, 1, 1, 1};
    p.sigma_constant = {2.7927, 3.0435, 2.7852, 3.48508556586, 3.53543525721,
                        2.9296, 2.4422, 3.4654, 2.0177};
    p.sigma_term_component = {0, 0};
    p.sigma_term_amplitude = {10.11, -1.417};
    p.sigma_term_exponent = {-0.01775, -0.01146};
    p.dispersion_energy_over_k = {353.95, 277.174, 173.44025, 232.687201645,
                                  453.265244384, 70.0, 249.26, 500.0, 650.0};
    p.molar_mass_kg_per_mol = {0.01801528, 0.06108, 0.04401, 0.06209, 0.10408,
                               0.0610168, 0.06001, 0.01902, 0.01701};
    p.charge_number = {0, 0, 0, 1, -1, -1, -2, 1, -1};
    p.packing_diameter_angstrom = {0, 0, 0, 3.0668752979568, 3.1111830263448,
                                   2.578048, 2.149136, 3.049552, 1.775576};
    p.debye_huckel_diameter_angstrom = p.packing_diameter_angstrom;
    p.born_diameter_angstrom = {0, 0, 0, 3.53322927146, 3.54107030822, 3, 3,
                                1.218, 3.08107689400404};
    p.solvation_factor = {1.5, 1, 1, 1, 1, 1, 1, 1, 1};
    p.relative_permittivity = {78.41043478578715, 32, 0, 0, 0, 0, 0, 0, 0};
    p.permittivity_term_component = {0, 0};
    p.permittivity_term_coefficient = {-0.36133766233691506, 0.00076555618295};
    p.permittivity_term_power = {1, 2};
    p.ion_suppression_coefficient = {0, 0, 0, 0, 0, 0, 0, 0, 0};  // revision 2 declares solvent-only; the physics layer spells that IonSuppression with alpha = 0
    p.k_ij[2 * 9 + 0] = p.k_ij[0 * 9 + 2] = 0.013262879176919628;
    p.electrolyte_form = epcsaft::ElectrolyteForm::DebyeHuckelBorn;
    p.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
    p.permittivity_form = epcsaft::PermittivityForm::IonSuppression;
    p.ionic_region_relative_permittivity = 8.0;
    p.site_component = {1, 1, 0, 0, 2, 2};
    p.site_multiplicity = {1, 1, 1, 1, 1, 1};
    p.association_active.assign(36, 0);
    p.association_energy_over_k.assign(36, 0.0);
    p.association_volume.assign(36, 0.0);
    const auto set_edge = [&](std::size_t a, std::size_t b, double energy, double volume) {
        p.association_active[a * 6 + b] = p.association_active[b * 6 + a] = 1;
        p.association_energy_over_k[a * 6 + b] = p.association_energy_over_k[b * 6 + a] = energy;
        p.association_volume[a * 6 + b] = p.association_volume[b * 6 + a] = volume;
    };
    set_edge(0, 1, 2586.3, 0.03747);
    set_edge(2, 3, 2425.7, 0.04509);
    set_edge(4, 3, 1212.85, 0.04509);
    set_edge(5, 2, 1212.85, 0.04509);
    set_edge(0, 3, 2506.0, 0.0411037990944876);
    set_edge(1, 2, 2506.0, 0.0411037990944876);
    p.k_ij[1] = p.k_ij[9] = -0.07352749874985018;
    p.k_ij[6 * 9] = p.k_ij[6] = -0.25;
    p.k_ij[7 * 9] = p.k_ij[7] = 0.25;
    p.k_ij[8 * 9] = p.k_ij[8] = -0.25;
    p.k_ij[3 * 9 + 4] = p.k_ij[4 * 9 + 3] = -0.00201813457644;
    return p;
}

struct ReferenceTensor final {
    std::vector<double> gradient;
    std::vector<double> hessian;
};

inline ReferenceTensor reference_tensor(
    const epcsaft::Parameters& parameters,
    double temperature,
    double volume,
    const std::vector<double>& amounts,
    const std::vector<double>& site_fractions,
    std::optional<epcsaft::ActiveParameter> active = std::nullopt,
    bool association_only = false
) {
    using AD = CppAD::AD<double>;
    const std::size_t components = parameters.component_count();
    const std::size_t sites_count = site_fractions.size();
    const std::size_t width = 2 + components + (active.has_value() ? 1 : 0);
    const std::size_t full_width = width + sites_count;
    std::vector<double> point;
    point.reserve(full_width);
    point.push_back(temperature);
    point.push_back(volume);
    point.insert(point.end(), amounts.begin(), amounts.end());
    if (active.has_value()) {
        const auto& values = active->family == epcsaft::ParameterFamily::Kij
            ? parameters.k_ij
            : active->family == epcsaft::ParameterFamily::DispersionEnergyOverK
            ? parameters.dispersion_energy_over_k
            : active->family == epcsaft::ParameterFamily::AssociationEnergyOverK
            ? parameters.association_energy_over_k
            : parameters.association_volume;
        point.push_back(values[active->indices.front()]);
    }
    point.insert(point.end(), site_fractions.begin(), site_fractions.end());
    std::vector<AD> independent(point.begin(), point.end());
    CppAD::Independent(independent);
    std::vector<AD> amount_values(
        independent.begin() + 2,
        independent.begin() + 2 + static_cast<std::ptrdiff_t>(components)
    );
    AD total = AD(0.0);
    for (const AD& amount : amount_values) total += amount;
    const AD density = total / independent[1];
    std::vector<AD> composition;
    for (const AD& amount : amount_values) composition.push_back(amount / total);
    const std::vector<AD> sites(independent.begin() + width, independent.end());
    std::vector<epcsaft::detail::ActiveParameter> active_set;
    if (active.has_value()) active_set.push_back(*active);
    const auto tape_parameters = epcsaft::detail::make_tape_parameters(
        parameters, active_set,
        std::vector<AD>(independent.begin() + 2 + components,
                        independent.begin() + width)
    );
    const auto sigma = epcsaft::segment_sigma(independent[0], tape_parameters);
    const auto diameter = epcsaft::segment_diameter(
        independent[0], tape_parameters, sigma
    );
    std::vector<double> reference_composition;
    double total_double = 0.0;
    for (double amount : amounts) total_double += amount;
    for (double amount : amounts) reference_composition.push_back(amount / total_double);
    const auto electrolyte = epcsaft::electrolyte_contribution(
        independent[0], density, temperature, total_double / volume,
        composition, reference_composition, parameters
    );
    const auto born = epcsaft::born_contribution(
        independent[0], composition, parameters,
        electrolyte.bulk_relative_permittivity
    );
    const auto terms = epcsaft::residual_helmholtz(
        independent[0], density, composition, tape_parameters, sigma, diameter,
        sites, electrolyte.debye_huckel, born
    );
    const auto delta = epcsaft::association_delta_values(
        independent[0], density, composition, tape_parameters, sigma, diameter,
        tape_parameters.association_energy_over_k,
        tape_parameters.association_volume
    );
    const auto residuals = epcsaft::association_residuals(
        density, composition, tape_parameters, sites, delta
    );
    AD value = total * (association_only ? terms.association : (
        terms.hard_chain + terms.dispersion + terms.debye_huckel
        + terms.born + terms.association
    ));
    std::vector<AD> dependents = {value};
    for (const AD& residual : residuals) dependents.push_back(residual);
    CppAD::ADFun<double> function(independent, dependents);
    const auto jacobian = function.Jacobian(point);
    Eigen::MatrixXd site_jacobian(
        static_cast<Eigen::Index>(sites_count),
        static_cast<Eigen::Index>(full_width)
    );
    for (std::size_t site = 0; site < sites_count; ++site) {
        for (std::size_t coordinate = 0; coordinate < full_width; ++coordinate) {
            site_jacobian(static_cast<Eigen::Index>(site),
                          static_cast<Eigen::Index>(coordinate)) =
                jacobian[(site + 1) * full_width + coordinate];
        }
    }
    const auto site_factorization = site_jacobian.rightCols(
        static_cast<Eigen::Index>(sites_count)
    );
    const Eigen::MatrixXd x_y = -Eigen::FullPivLU<Eigen::MatrixXd>(site_factorization)
        .solve(site_jacobian.leftCols(static_cast<Eigen::Index>(width)));
    const auto value_jacobian = function.Jacobian(point);
    Eigen::VectorXd value_x(static_cast<Eigen::Index>(sites_count));
    Eigen::VectorXd value_y(static_cast<Eigen::Index>(width));
    for (std::size_t site = 0; site < sites_count; ++site) {
        value_x(static_cast<Eigen::Index>(site)) =
            value_jacobian[width + site];
    }
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        value_y(static_cast<Eigen::Index>(coordinate)) = value_jacobian[coordinate];
    }
    ReferenceTensor result;
    result.gradient.resize(width);
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        result.gradient[coordinate] = value_y(static_cast<Eigen::Index>(coordinate))
            + value_x.dot(x_y.col(static_cast<Eigen::Index>(coordinate)));
    }
    result.hessian.assign(width * width, 0.0);
    std::vector<std::vector<double>> hessians;
    for (std::size_t dependent = 0; dependent <= sites_count; ++dependent) {
        hessians.push_back(function.Hessian(point, dependent));
    }
    for (std::size_t row = 0; row < width; ++row) {
        Eigen::VectorXd direction_row = Eigen::VectorXd::Zero(
            static_cast<Eigen::Index>(full_width)
        );
        for (std::size_t site = 0; site < sites_count; ++site) {
            direction_row(static_cast<Eigen::Index>(width + site)) =
                x_y(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(row));
        }
        direction_row(static_cast<Eigen::Index>(row)) = 1.0;
        for (std::size_t column = 0; column < width; ++column) {
            Eigen::VectorXd direction_column = Eigen::VectorXd::Zero(
                static_cast<Eigen::Index>(full_width)
            );
            for (std::size_t site = 0; site < sites_count; ++site) {
                direction_column(static_cast<Eigen::Index>(width + site)) =
                    x_y(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(column));
            }
            direction_column(static_cast<Eigen::Index>(column)) = 1.0;
            Eigen::VectorXd second_site(static_cast<Eigen::Index>(sites_count));
            for (std::size_t site = 0; site < sites_count; ++site) {
                double contraction = 0.0;
                for (std::size_t left = 0; left < full_width; ++left) {
                    for (std::size_t right = 0; right < full_width; ++right) {
                        contraction += direction_row(static_cast<Eigen::Index>(left))
                            * hessians[site + 1][left * full_width + right]
                            * direction_column(static_cast<Eigen::Index>(right));
                    }
                }
                second_site(static_cast<Eigen::Index>(site)) = contraction;
            }
            const Eigen::VectorXd x_second = -Eigen::FullPivLU<Eigen::MatrixXd>(
                site_factorization
            ).solve(second_site);
            double value_second = 0.0;
            for (std::size_t left = 0; left < full_width; ++left) {
                for (std::size_t right = 0; right < full_width; ++right) {
                    value_second += direction_row(static_cast<Eigen::Index>(left))
                        * hessians[0][left * full_width + right]
                        * direction_column(static_cast<Eigen::Index>(right));
                }
            }
            result.hessian[row * width + column] = value_second + value_x.dot(x_second);
        }
    }
    return result;
}

}  // namespace properties_support
