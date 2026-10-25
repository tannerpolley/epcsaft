#pragma once

#include "properties/helmholtz_model.hpp"
#include "core/physical_constants.hpp"

#include <Eigen/Dense>
#include <cppad/cppad.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
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

// P4 thermal test records: one nonconstant ideal form per correlation family (cp = 30 J/(mol K)
// with h(300 K) = 1000 J/mol and s(300 K) = 50 J/(mol K)), one component interval each.
inline std::vector<std::pair<const char*, epcsaft::IdealCorrelation>> p3_forms() {
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    return {
        {"polynomial", epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{30.0, 0.0}, 300.0, 1000.0, 50.0}, 1.0e5}},
        {"shomate", epcsaft::IdealCorrelation{
            epcsaft::IdealShomate{{30.0, 0.0, 0.0, 0.0, 0.0, -9.0,
                                   50.0 - 30.0 * std::log(0.3), 0.0}, 1000.0}, 1.0e5}},
        {"nasa7", epcsaft::IdealCorrelation{
            epcsaft::IdealNASA7{{30.0 / R, 0.0, 0.0, 0.0, 0.0,
                                 1000.0 / R - 30.0 * 300.0 / R,
                                 50.0 / R - 30.0 * std::log(300.0) / R}}, 1.0e5}},
        {"nasa9", epcsaft::IdealCorrelation{
            epcsaft::IdealNASA9{{0.0, 0.0, 30.0 / R, 0.0, 0.0, 0.0, 0.0,
                                 1000.0 / R - 30.0 * 300.0 / R,
                                 50.0 / R - 30.0 * std::log(300.0) / R}}, 1.0e5}},
        {"dippr107", epcsaft::IdealCorrelation{
            epcsaft::IdealDIPPR107{{0.0, 30.0, 0.0, 0.0, 0.0},
                                   300.0, 1000.0, 50.0}, 1.0e5}},
    };
}

inline epcsaft::ThermochemistryRecord p3_record(const epcsaft::IdealCorrelation& correlation) {
    epcsaft::ThermochemistryRecord result;
    result.standard_pressure = 1.0e5;
    result.component_intervals.push_back({
        epcsaft::IdealInterval{250.0, 350.0, true, true, correlation}
    });
    return result;
}

inline epcsaft::ThermochemistryRecord p3_binary_record() {
    epcsaft::ThermochemistryRecord result;
    result.standard_pressure = 1.0e5;
    result.component_intervals = {
        {epcsaft::IdealInterval{250.0, 350.0, true, true, p3_forms()[0].second}},
        {epcsaft::IdealInterval{250.0, 350.0, true, true,
            epcsaft::IdealCorrelation{epcsaft::IdealPolynomial{{40.0, 0.0}, 300.0, 2000.0, 70.0}, 1.0e5}}}
    };
    return result;
}

struct ReferenceTensor final {
    std::vector<double> gradient;
    std::vector<double> hessian;
    double scaled_condition = std::numeric_limits<double>::quiet_NaN();
    double first_sensitivity_residual = std::numeric_limits<double>::quiet_NaN();
    double second_sensitivity_residual = std::numeric_limits<double>::quiet_NaN();
    double sensitivity_error_proxy = std::numeric_limits<double>::quiet_NaN();
    std::size_t first_sensitivity_coordinate = 0;
    std::size_t second_sensitivity_row = 0;
    std::size_t second_sensitivity_column = 0;
};

inline ReferenceTensor reference_tensor(
    const epcsaft::Parameters& parameters,
    double temperature,
    double volume,
    const std::vector<double>& amounts,
    const std::vector<double>& site_fractions,
    std::vector<epcsaft::ActiveParameter> active = {},
    bool association_only = false
) {
    using AD = CppAD::AD<double>;
    const std::size_t components = parameters.component_count();
    const std::size_t sites_count = site_fractions.size();
    const std::size_t width = 2 + components + active.size();
    const std::size_t full_width = width + sites_count;
    std::vector<double> point;
    point.reserve(full_width);
    point.push_back(temperature);
    point.push_back(volume);
    point.insert(point.end(), amounts.begin(), amounts.end());
    for (const auto& item : active) {
        const auto& values = item.family == epcsaft::ParameterFamily::Kij
            ? parameters.k_ij
            : item.family == epcsaft::ParameterFamily::Lij ? parameters.l_ij
            : item.family == epcsaft::ParameterFamily::DispersionEnergyOverK
            ? parameters.dispersion_energy_over_k
            : item.family == epcsaft::ParameterFamily::AssociationEnergyOverK
            ? parameters.association_energy_over_k
            : parameters.association_volume;
        point.push_back(values[item.indices.front()]);
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
    active_set = active;
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
        independent[0], density,
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
    const Eigen::MatrixXd site_factorization = site_jacobian.rightCols(
        static_cast<Eigen::Index>(sites_count)
    );
    const Eigen::FullPivLU<Eigen::MatrixXd> factorization(site_factorization);
    const Eigen::MatrixXd x_y = -factorization.solve(
        site_jacobian.leftCols(static_cast<Eigen::Index>(width))
    );
    const Eigen::Map<const Eigen::VectorXd> value_x(
        jacobian.data() + width, static_cast<Eigen::Index>(sites_count));
    const Eigen::Map<const Eigen::VectorXd> value_y(
        jacobian.data(), static_cast<Eigen::Index>(width));
    ReferenceTensor result;
    const auto infinity_norm = [](const Eigen::MatrixXd& matrix) {
        if (matrix.size() == 0) return 0.0;
        if (!matrix.allFinite()) return std::numeric_limits<double>::quiet_NaN();
        const double result = matrix.cwiseAbs().rowwise().sum().maxCoeff();
        return std::isfinite(result) ? result : std::numeric_limits<double>::quiet_NaN();
    };
    const auto vector_infinity_norm = [](const Eigen::VectorXd& vector) {
        if (vector.size() == 0) return 0.0;
        if (!vector.allFinite()) return std::numeric_limits<double>::quiet_NaN();
        const double result = vector.cwiseAbs().maxCoeff();
        return std::isfinite(result) ? result : std::numeric_limits<double>::quiet_NaN();
    };
    const auto residual_ratio = [&](const Eigen::VectorXd& actual_residual,
                                    const Eigen::VectorXd& unknown,
                                    const Eigen::VectorXd& rhs,
                                    double matrix_norm) {
        const double numerator = vector_infinity_norm(actual_residual);
        const double unknown_norm = vector_infinity_norm(unknown);
        const double rhs_norm = vector_infinity_norm(rhs);
        const double denominator = matrix_norm * unknown_norm + rhs_norm;
        if (!std::isfinite(numerator) || !std::isfinite(unknown_norm)
            || !std::isfinite(rhs_norm) || !std::isfinite(matrix_norm)
            || !std::isfinite(denominator))
            return std::numeric_limits<double>::quiet_NaN();
        if (denominator == 0.0) return numerator == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
        return numerator / denominator;
    };
    const double total_amount = total_double;
    const auto scale = [&](std::size_t coordinate) {
        if (coordinate == 0) return temperature;
        if (coordinate == 1) return volume;
        if (coordinate < 2 + components) return total_amount;
        switch (active[coordinate - 2 - components].family) {
        case epcsaft::ParameterFamily::SegmentCount: return 1.0;
        case epcsaft::ParameterFamily::SigmaConstant: return 3.0;
        case epcsaft::ParameterFamily::DispersionEnergyOverK: return 1000.0;
        case epcsaft::ParameterFamily::Kij:
        case epcsaft::ParameterFamily::Lij: return 1.0;
        case epcsaft::ParameterFamily::AssociationEnergyOverK: return 1000.0;
        case epcsaft::ParameterFamily::AssociationVolume: return 0.1;
        }
        return 1.0;
    };
    const Eigen::MatrixXd site_diagonal = Eigen::Map<const Eigen::VectorXd>(
        site_fractions.data(), static_cast<Eigen::Index>(sites_count)).asDiagonal();
    const Eigen::MatrixXd scaled_jacobian = site_factorization * site_diagonal;
    const Eigen::FullPivLU<Eigen::MatrixXd> scaled_factorization(scaled_jacobian);
    if (sites_count != 0 && scaled_factorization.rank() == static_cast<Eigen::Index>(sites_count)) {
        const Eigen::MatrixXd inverse = scaled_factorization.solve(
            Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(sites_count))
        );
        const double matrix_norm = infinity_norm(scaled_jacobian);
        const double inverse_norm = infinity_norm(inverse);
        result.scaled_condition = std::isfinite(matrix_norm) && std::isfinite(inverse_norm)
            ? matrix_norm * inverse_norm : std::numeric_limits<double>::quiet_NaN();
    }
    const double matrix_norm = infinity_norm(scaled_jacobian);
    result.first_sensitivity_residual = 0.0;
    result.second_sensitivity_residual = 0.0;
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        Eigen::VectorXd unknown = x_y.col(static_cast<Eigen::Index>(coordinate)) * scale(coordinate);
        for (std::size_t site = 0; site < sites_count; ++site)
            unknown(static_cast<Eigen::Index>(site)) /= site_fractions[site];
        const Eigen::VectorXd rhs = site_jacobian.col(static_cast<Eigen::Index>(coordinate)) * scale(coordinate);
        const double residual = residual_ratio(
            (site_factorization * x_y.col(static_cast<Eigen::Index>(coordinate))
                + site_jacobian.col(static_cast<Eigen::Index>(coordinate))) * scale(coordinate),
            unknown, rhs, matrix_norm);
        if (!std::isfinite(residual)) {
            result.first_sensitivity_residual = std::numeric_limits<double>::quiet_NaN();
            result.first_sensitivity_coordinate = coordinate;
        }
        else if (std::isfinite(result.first_sensitivity_residual) && residual > result.first_sensitivity_residual) {
            result.first_sensitivity_residual = residual;
            result.first_sensitivity_coordinate = coordinate;
        }
    }
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
    std::vector<Eigen::VectorXd> directions(width, Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(full_width)));
    for (std::size_t coordinate = 0; coordinate < width; ++coordinate) {
        for (std::size_t site = 0; site < sites_count; ++site)
            directions[coordinate](static_cast<Eigen::Index>(width + site)) =
                x_y(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(coordinate));
        directions[coordinate](static_cast<Eigen::Index>(coordinate)) = 1.0;
    }
    for (std::size_t row = 0; row < width; ++row) {
        const auto& direction_row = directions[row];
        for (std::size_t column = 0; column < width; ++column) {
            const auto& direction_column = directions[column];
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
            const Eigen::VectorXd x_second = -factorization.solve(second_site);
            Eigen::VectorXd transformed_second = x_second * scale(row) * scale(column);
            for (std::size_t site = 0; site < sites_count; ++site)
                transformed_second(static_cast<Eigen::Index>(site)) /= site_fractions[site];
            const Eigen::VectorXd second_rhs = second_site * scale(row) * scale(column);
            const double sensitivity_residual = residual_ratio(
                (site_factorization * x_second + second_site) * scale(row) * scale(column),
                transformed_second, second_rhs, matrix_norm);
            if (!std::isfinite(sensitivity_residual)) {
                result.second_sensitivity_residual = std::numeric_limits<double>::quiet_NaN();
                result.second_sensitivity_row = row;
                result.second_sensitivity_column = column;
            }
            else if (std::isfinite(result.second_sensitivity_residual) && sensitivity_residual > result.second_sensitivity_residual) {
                result.second_sensitivity_residual = sensitivity_residual;
                result.second_sensitivity_row = row;
                result.second_sensitivity_column = column;
            }
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
    if (sites_count == 0) {
        result.scaled_condition = 1.0;
        result.first_sensitivity_residual = 0.0;
        result.second_sensitivity_residual = 0.0;
        result.sensitivity_error_proxy = 0.0;
    } else if (std::isfinite(result.scaled_condition)
        && std::isfinite(result.first_sensitivity_residual)
        && std::isfinite(result.second_sensitivity_residual)) {
        result.sensitivity_error_proxy = result.scaled_condition
            * (std::max(result.first_sensitivity_residual, result.second_sensitivity_residual)
               + static_cast<double>(sites_count) * std::numeric_limits<double>::epsilon());
    }
    return result;
}

}  // namespace properties_support
