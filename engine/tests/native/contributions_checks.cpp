// Properties of the Helmholtz terms that hold for every parameter set, checked
// against closed forms and exact symmetries rather than against stored numbers;
// they fail if the kernels are wrong at all.

#include "contributions/association.hpp"
#include "contributions/born.hpp"
#include "contributions/electrolyte.hpp"
#include "contributions/helmholtz.hpp"
#include "contributions/ideal.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Not assert: the release build defines NDEBUG, which would compile every check
// in this file away and leave it passing vacuously.
#define CHECK(condition) recording_support::check(                         \
    static_cast<bool>(condition), __FILE__, __LINE__, #condition)

namespace {

using epcsaft::Parameters;

bool close(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

template <typename Function>
double simpson(Function function, double lower, double upper, int intervals = 2000) {
    const double step = (upper - lower) / intervals;
    double sum = function(lower) + function(upper);
    for (int i = 1; i < intervals; ++i)
        sum += (i % 2 == 0 ? 2.0 : 4.0) * function(lower + i * step);
    return sum * step / 3.0;
}

// A two-component neutral mixture with no sites: the smallest parameter set the
// hard-chain and dispersion terms accept.
Parameters neutral_binary() {
    Parameters p;
    p.segment_count = {1.0, 2.457};
    p.sigma_constant = {3.7039, 3.6184};
    p.dispersion_energy_over_k = {150.03, 267.11};
    p.molar_mass_kg_per_mol = {0.016043, 0.058124};
    p.k_ij.assign(4, 0.0);
    p.k_ij_reciprocal_temperature_slope.assign(4, 0.0);
    p.k_ij_reference_temperature_k.assign(4, 0.0);
    p.l_ij.assign(4, 0.0);
    p.charge_number = {0, 0};
    p.packing_diameter_angstrom = {0.0, 0.0};
    p.relative_permittivity = {0.0, 0.0};
    p.ion_suppression_coefficient = {0.0, 0.0};
    p.dielectric_volume_m3_per_mol = {0.0, 0.0};
    p.dielectric_temperature_k = {0.0, 0.0};
    p.dipole_scaling_debye2 = {0.0, 0.0};
    p.polarizability_scaling_angstrom3 = {0.0, 0.0};
    p.correlation_integral_parameter = {0.0, 0.0};
    p.debye_huckel_diameter_angstrom = {0.0, 0.0};
    p.born_diameter_angstrom = {0.0, 0.0};
    p.solvation_factor = {0.0, 0.0};
    return p;
}

Parameters electrolyte_fixture() {
    Parameters p;
    p.segment_count = {1.0, 1.3, 1.0, 1.0};
    p.sigma_constant = {3.1, 3.7, 2.6, 3.4};
    p.dispersion_energy_over_k = {200.0, 240.0, 180.0, 220.0};
    p.molar_mass_kg_per_mol = {0.018, 0.046, 0.023, 0.071};
    p.k_ij.assign(16, 0.0);
    p.k_ij_reciprocal_temperature_slope.assign(16, 0.0);
    p.k_ij_reference_temperature_k.assign(16, 0.0);
    p.l_ij.assign(16, 0.0);
    p.charge_number = {0, 0, 1, -2};
    p.packing_diameter_angstrom = {0.0, 0.0, 1.7, 2.8};
    p.debye_huckel_diameter_angstrom = {0.0, 0.0, 2.1, 3.9};
    p.born_diameter_angstrom = {0.0, 0.0, 3.2, 5.1};
    p.solvation_factor = {1.2, 0.8, 1.1, 1.3};
    p.relative_permittivity = {70.0, 55.0, 12.0, 18.0};
    p.permittivity_term_component = {0, 1};
    p.permittivity_term_coefficient = {0.02, -0.01};
    p.permittivity_term_power = {1, 2};
    p.permittivity_reference_temperature_k = {300.0, 300.0, 300.0, 300.0};
    p.ion_suppression_coefficient = {0.0, 0.0, 5.0, 8.0};
    p.dielectric_volume_m3_per_mol = {1.5e-4, 2.5e-4, 0.0, 0.0};
    p.dielectric_temperature_k = {450.0, 600.0, 0.0, 0.0};
    p.dipole_scaling_debye2 = {1.0, 2.0, 0.0, 0.0};
    p.polarizability_scaling_angstrom3 = {1.1, 1.5, 0.0, 0.0};
    p.correlation_integral_parameter = {0.2, 0.5, 0.0, 0.0};
    p.electrolyte_form = epcsaft::ElectrolyteForm::DebyeHuckelBorn;
    p.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
    p.permittivity_form = epcsaft::PermittivityForm::IonSuppression;
    p.ionic_region_relative_permittivity = 15.0;
    p.born_c_shell = 1.0;
    p.born_c_dielectric = 1.0;
    return p;
}

std::string model_inputs(const Parameters& parameters, double temperature,
                         double molar_density, const std::vector<double>& composition) {
    return recording_support::fields({
        recording_support::field("T", temperature),
        recording_support::field("rho", molar_density),
        recording_support::field("x", composition),
        recording_support::field("segments", parameters.segment_count),
        recording_support::field("sigma_A", parameters.sigma_constant),
        recording_support::field("epsilon_over_k_K", parameters.dispersion_energy_over_k),
        recording_support::field("molar_mass_kg_per_mol", parameters.molar_mass_kg_per_mol),
        recording_support::field("charge", parameters.charge_number),
        recording_support::field("kij", parameters.k_ij),
        recording_support::field("site_component", parameters.site_component),
        recording_support::field("site_multiplicity", parameters.site_multiplicity),
        recording_support::field("association_active", parameters.association_active),
        recording_support::field("association_energy_over_k_K", parameters.association_energy_over_k),
        recording_support::field("association_volume", parameters.association_volume),
    });
}

epcsaft::HelmholtzResult<double> evaluate(
    const Parameters& p,
    double temperature,
    double molar_density,
    const std::vector<double>& x
) {
    const auto sigma = epcsaft::segment_sigma(temperature, p);
    const auto diameter = epcsaft::segment_diameter(temperature, p, sigma);
    const auto sites = epcsaft::solve_site_fractions(
        temperature, molar_density, x, p, sigma, diameter
    );
    return epcsaft::residual_helmholtz(
        temperature, molar_density, x, p, sigma, diameter, sites,
        epcsaft::association_products(molar_density, x, p, sites, epcsaft::association_delta(temperature, molar_density, x, p, sigma, diameter)));
}

// Carnahan-Starling. For one component the log term in the hard-sphere
// expression cancels exactly (zeta_2^3 / zeta_3^2 == zeta_0), leaving
// (4 eta - 3 eta^2) / (1 - eta)^2. A single segment makes the chain term vanish,
// so hard_chain is the hard-sphere free energy alone.
void hard_sphere_matches_carnahan_starling() {
    Parameters p;
    p.segment_count = {1.0};
    p.sigma_constant = {3.7039};
    p.dispersion_energy_over_k = {150.03};
    p.molar_mass_kg_per_mol = {0.016043};
    p.k_ij = {0.0};
    p.k_ij_reciprocal_temperature_slope = {0.0};
    p.k_ij_reference_temperature_k = {0.0};
    p.l_ij = {0.0};
    p.charge_number = {0};
    p.packing_diameter_angstrom = {0.0};
    p.relative_permittivity = {0.0};
    p.ion_suppression_coefficient = {0.0};
    p.dielectric_volume_m3_per_mol = {0.0};
    p.dielectric_temperature_k = {0.0};
    p.dipole_scaling_debye2 = {0.0};
    p.polarizability_scaling_angstrom3 = {0.0};
    p.correlation_integral_parameter = {0.0};
    p.debye_huckel_diameter_angstrom = {0.0};
    p.born_diameter_angstrom = {0.0};
    p.solvation_factor = {0.0};

    const double temperature = 250.0;
    const std::vector<double> x = {1.0};
    for (const double molar_density : {1.0e3, 8.0e3, 2.0e4}) {
        const auto sigma = epcsaft::segment_sigma(temperature, p);
        const auto diameter = epcsaft::segment_diameter(temperature, p, sigma);
        const double eta = epcsaft::packing_fraction_from_geometry(
            molar_density, x, p, diameter
        );
        const double expected =
            (4.0 * eta - 3.0 * eta * eta) / ((1.0 - eta) * (1.0 - eta));
        const auto result = evaluate(p, temperature, molar_density, x);
        const auto inputs = recording_support::fields({
            model_inputs(p, temperature, molar_density, x),
            recording_support::field("eta", eta),
        });
        const auto point = recording_support::number(molar_density);
        recording_support::compare(std::string("eta-") + point, "hard_chain", "1", result.hard_chain, expected,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "closed_form", inputs);
    }
}

void mixture_contributions_match_independent_equations() {
    const auto p = neutral_binary();
    const double temperature = 300.0;
    const double density = 1.2e4;
    const std::vector<double> x = {0.4, 0.6};
    const auto result = evaluate(p, temperature, density, x);
    const double expected_hard_chain = 3.2311322600846055;
    const double expected_dispersion = -6.390063657181498;
    const auto inputs = recording_support::fields({
        model_inputs(p, temperature, density, x),
        recording_support::field("independent_reference", "Gross2001 A.8-A.17"),
    });
    recording_support::compare("hard-chain", "hard_chain", "1", result.hard_chain, expected_hard_chain,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "independent_closed_form", inputs);
    recording_support::compare("dispersion", "dispersion", "1", result.dispersion, expected_dispersion,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "independent_closed_form", inputs);
}

// Every residual term is a departure from the ideal gas, so all of them vanish
// together as the density goes to zero.
void residual_vanishes_at_zero_density() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    double previous = 1.0;
    double previous_density = 0.0;
    std::size_t point_index = 0;
    for (const double molar_density : {1.0e1, 1.0e0, 1.0e-1, 1.0e-2}) {
        const auto result = evaluate(p, 300.0, molar_density, x);
        const double total = result.hard_chain + result.dispersion;
        const bool decreasing = std::abs(total) < std::abs(previous);
        recording_support::compare(std::string("density-monotonic/") + std::to_string(point_index), "abs_residual",
            "1", decreasing ? 1.0 : 0.0, 1.0, recording_support::Rule::Absolute, 0.0, 0.0, "density-limit-ladder",
            recording_support::fields({
                model_inputs(p, 300.0, molar_density, x),
                recording_support::field("previous_density", previous_density),
                recording_support::field("previous_abs_residual", std::abs(previous)),
                recording_support::field("current_abs_residual", std::abs(total)),
            }));
        previous = total;
        previous_density = molar_density;
        ++point_index;
    }
    recording_support::compare("zero-density-limit", "residual_helmholtz", "1", previous, 0.0,
        recording_support::Rule::StrictAbsolute, 1.0e-3, 0.0, "limit", model_inputs(p, 300.0, 1.0e-2, x));
    auto spheres = p;
    spheres.dispersion_energy_over_k.assign(2, 0.0);
    const std::vector<double> pure = {1.0, 0.0};
    const auto sigma = epcsaft::segment_sigma(300.0, spheres);
    const auto diameter = epcsaft::segment_diameter(300.0, spheres, sigma);
    epcsaft::Taylor<2> density;
    density.c[1] = 1.0;
    const auto terms = epcsaft::residual_helmholtz(300.0, density, pure, spheres, sigma, diameter,
        std::vector<epcsaft::Taylor<2>>{}, {});
    const double v = epcsaft::PI / 6.0 * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3
        * std::pow(0.88 * spheres.sigma_constant[0], 3);
    CHECK(std::abs(terms.hard_chain.c[1] / (4.0 * v) - 1.0) < 1e-14);
    CHECK(std::abs(2.0 * terms.hard_chain.c[2] / (10.0 * v * v) - 1.0) < 1e-14);
    spheres.ideal_gas = true;
    const auto ideal = epcsaft::residual_helmholtz(300.0, density, pure, spheres, sigma, diameter,
        std::vector<epcsaft::Taylor<2>>{}, {});
    CHECK(ideal.hard_chain.c[1] == 0.0 && ideal.hard_chain.c[2] == 0.0);
    auto associated = p;
    associated.site_component = {0, 0, 1};
    associated.site_multiplicity = {1, 2, 1};
    associated.association_active.assign(9, 1);
    associated.association_energy_over_k.assign(9, 1500.0);
    associated.association_volume.assign(9, 0.03);
    const auto ds = epcsaft::association_delta(300.0, density, x, associated, sigma, diameter);
    const auto sites = epcsaft::association_density_coefficients(density, x, associated, ds);
    const auto residuals = epcsaft::association_residuals(density, x, associated, sites, ds);
    for (std::size_t site = 0; site < sites.size(); ++site) {
        CHECK(residuals[site].c[0] == 0.0);
        for (std::size_t k = 1; k <= 2; ++k)
            CHECK(std::abs(residuals[site].c[k]) <= 1e-14 * std::abs(sites[site].c[k]));
    }

}

// Relabelling the components is not a physical change. This is the check that
// catches an index used where a component was meant, which is the failure a
// pair loop invites.
void relabelling_components_changes_nothing() {
    const auto p = neutral_binary();
    Parameters swapped = p;
    swapped.segment_count = {p.segment_count[1], p.segment_count[0]};
    swapped.sigma_constant = {p.sigma_constant[1], p.sigma_constant[0]};
    swapped.dispersion_energy_over_k = {
        p.dispersion_energy_over_k[1], p.dispersion_energy_over_k[0]
    };
    swapped.molar_mass_kg_per_mol = {
        p.molar_mass_kg_per_mol[1], p.molar_mass_kg_per_mol[0]
    };

    const double temperature = 300.0;
    const double molar_density = 1.2e4;
    const auto direct = evaluate(p, temperature, molar_density, {0.4, 0.6});
    const auto relabelled =
        evaluate(swapped, temperature, molar_density, {0.6, 0.4});
    const auto inputs = recording_support::fields({
        model_inputs(p, temperature, molar_density, {0.4, 0.6}),
        recording_support::field("x_relabelled", std::vector<double>{0.6, 0.4}),
        recording_support::field("segments_relabelled", swapped.segment_count),
        recording_support::field("sigma_relabelled_A", swapped.sigma_constant),
        recording_support::field("epsilon_relabelled_over_k_K", swapped.dispersion_energy_over_k),
        recording_support::field("molar_mass_relabelled_kg_per_mol", swapped.molar_mass_kg_per_mol),
    });
    recording_support::compare("hard-chain", "hard_chain", "1", direct.hard_chain, relabelled.hard_chain,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "permutation", inputs);
    recording_support::compare("dispersion", "dispersion", "1", direct.dispersion, relabelled.dispersion,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "permutation", inputs);
}

// One associating component, two sites, one bonding pair.
Parameters two_site_binary() {
    Parameters p = neutral_binary();
    p.site_component = {0, 0};
    p.site_multiplicity = {1, 1};
    p.association_active = {0, 1, 1, 0};
    p.association_energy_over_k = {0.0, 2425.7, 2425.7, 0.0};
    p.association_volume = {0.0, 0.0451, 0.0451, 0.0};
    return p;
}

// The solve has to land inside the domain and satisfy the mass-action residual to the policy tolerance.
void association_solve_satisfies_mass_action() {
    const Parameters p = two_site_binary();

    const double temperature = 300.0;
    const double molar_density = 1.5e4;
    const std::vector<double> x = {0.4, 0.6};
    const auto sigma = epcsaft::segment_sigma(temperature, p);
    const auto diameter = epcsaft::segment_diameter(temperature, p, sigma);
    const auto fractions = epcsaft::solve_site_fractions(
        temperature, molar_density, x, p, sigma, diameter
    );
    CHECK(fractions.size() == 2);
    for (std::size_t site = 0; site < fractions.size(); ++site) {
        const double value = fractions[site];
        const auto inputs = recording_support::fields({
            model_inputs(p, temperature, molar_density, x),
            recording_support::field("site_fractions", fractions),
            recording_support::field("site", site),
            recording_support::field("floor", epcsaft::policy::SITE_FRACTION_FLOOR.value),
            recording_support::field("ceiling", epcsaft::policy::SITE_FRACTION_CEILING.value),
        });
        const bool above_floor = value > epcsaft::policy::SITE_FRACTION_FLOOR.value;
        recording_support::compare(std::string("site-floor-") + std::to_string(site), "site_fraction_above_floor",
            "bool", above_floor ? 1.0 : 0.0, 1.0, recording_support::Rule::Absolute, 0.0, 0.0, "association-domain",
            inputs);
        const bool below_ceiling = value <= epcsaft::policy::SITE_FRACTION_CEILING.value;
        recording_support::compare(std::string("site-ceiling-") + std::to_string(site), "site_fraction_below_ceiling",
            "bool", below_ceiling ? 1.0 : 0.0, 1.0, recording_support::Rule::Absolute, 0.0, 0.0, "association-domain",
            inputs);
    }
    const auto delta = epcsaft::association_delta(
        temperature, molar_density, x, p, sigma, diameter
    );
    const auto residuals = epcsaft::association_residuals(
        molar_density, x, p, fractions, delta
    );
    for (std::size_t site = 0; site < residuals.size(); ++site) {
        const double residual = residuals[site];
        recording_support::compare(std::string("mass-action-") + std::to_string(site), "association_residual", "1",
            residual, 0.0, recording_support::Rule::ReferenceOrOne, epcsaft::policy::ASSOCIATION_RESIDUAL.value, 0.0,
            "mass_action", recording_support::fields({model_inputs(p, temperature, molar_density, x),
                                       recording_support::field("site_fractions", fractions),
                                       recording_support::field("site", site)}));
    }
    // Identical sites on one component are indistinguishable.
    const bool identical_sites = close(fractions[0], fractions[1], 1.0e-13);
    recording_support::compare("site-symmetry", "site_fraction_symmetry", "bool", identical_sites ? 1.0 : 0.0, 1.0,
        recording_support::Rule::Absolute, 0.0, 0.0, "site-symmetry",
        recording_support::fields({model_inputs(p, temperature, molar_density, x),
                                   recording_support::field("site_fractions", fractions)}));
}

// A captured rejected trial at 298.15 K; density in mol/m^3, strengths in Angstrom^3.
// The previous step bound exhausted 50 iterations at residual 5.480282894154698e-12.
void association_unit_ceiling() {
    constexpr double temperature = 298.15;
    constexpr double density = 52992.32016447701;
    const std::vector<double> mole_fractions = {
        0.9967664629120051, 9.581789719163911e-29, 0.003232311646483076,
        6.310494279431576e-07, 5.712187821800883e-10, 5.938208650913641e-07};
    Parameters parameters;
    parameters.site_component = {0, 0, 1, 1, 2, 2, 3, 3};
    parameters.site_multiplicity = std::vector<int>(8, 1);
    const std::vector<double> strength = {
        0.0, 15964.491628191587, 0.0, 1986249.8013975348, 0.0, 4550137.1076187035, 0.0, 0.0,
        15964.491628191587, 0.0, 1986249.8013975348, 0.0, 4550137.1076187035, 0.0, 0.0, 0.0,
        0.0, 1986249.8013975348, 0.0, 286578809.8455787, 0.0, 633297185.5581268, 0.0, 61842.874171620366,
        1986249.8013975348, 0.0, 286578809.8455787, 0.0, 633297185.5581268, 0.0, 61842.874171620366, 0.0,
        0.0, 4550137.1076187035, 0.0, 633297185.5581268, 0.0, 1412151858.732317, 0.0, 0.0,
        4550137.1076187035, 0.0, 633297185.5581268, 0.0, 1412151858.732317, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 61842.874171620366, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 61842.874171620366, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    // The capture kept no segment geometry. The 99.7 % species is water: its segment number and sigma
    // (data/packets/held-2014-strategy2) for all six species give packing 0.44, inside the fluid domain.
    parameters.segment_count = std::vector<double>(6, 1.2047);
    const std::vector<double> diameter(6, 2.7927);
    epcsaft::core::Convergence convergence;
    const auto fractions = epcsaft::solve_site_fractions(temperature, density, mole_fractions,
        parameters, {}, diameter, &convergence, nullptr, &strength);
    const auto residuals = epcsaft::association_residuals(density, mole_fractions,
        parameters, fractions, strength);
    double residual = 0.0;
    for (double value : residuals) {
        CHECK(std::isfinite(value));
        residual = std::max(residual, std::abs(value));
    }
    const double tolerance = epcsaft::policy::SolverOptions{}.value(epcsaft::policy::association);
    CHECK(residual <= tolerance);
    CHECK(convergence.iterations <= epcsaft::policy::ASSOCIATION_MAX_ITERATIONS.value);
    CHECK(std::all_of(fractions.begin(), fractions.end(),
        [](double value) { return std::isfinite(value) && value > 0.0 && value <= 1.0; }));
    std::printf("association unit ceiling: residual=%.17g tolerance=%.17g iterations=%zu\n",
        residual, tolerance, convergence.iterations);
}

// Closed forms here are independent reductions of the accepted mass-action
// equation. They cover multiplicity and support without adding another solver.
void association_equations_match_closed_forms_and_limits() {
    const double temperature = 300.0, density = 1.5e4;
    const std::vector<double> pure = {1.0, 0.0};
    const Parameters symmetric = two_site_binary();
    const auto sigma = epcsaft::segment_sigma(temperature, symmetric);
    const auto diameter = epcsaft::segment_diameter(temperature, symmetric, sigma);
    const auto fractions = epcsaft::solve_site_fractions(
        temperature, density, pure, symmetric, sigma, diameter
    );
    const auto delta = epcsaft::association_delta(
        temperature, density, pure, symmetric, sigma, diameter
    );
    const double k = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3 * delta[1];
    const double expected_x = 2.0 / (1.0 + std::sqrt(1.0 + 4.0 * k));
    const auto symmetric_result = evaluate(symmetric, temperature, density, pure);
    const double expected_association = 2.0 * (
        std::log(expected_x) - expected_x / 2.0 + 0.5
    );
    const auto inputs = [&](const Parameters& parameters, const std::vector<double>& x,
                            const std::string& extra) {
        return recording_support::fields({model_inputs(parameters, temperature, density, x),
            recording_support::field("l_ij", parameters.l_ij), extra});
    };
    const auto record = [&](const char* id, const char* quantity, double model,
                            double reference, const std::string& values, double tolerance = 1.0e-12,
                            const char* unit = "1") {
        recording_support::compare(id, quantity, unit, model, reference, recording_support::Rule::ReferenceOrOne, 0.0,
            tolerance, "independent_closed_form", values);
    };
    record("symmetric-x-0", "site_fraction", fractions[0], expected_x,
           inputs(symmetric, pure, recording_support::field("K", k)));
    record("symmetric-x-1", "site_fraction", fractions[1], expected_x,
           inputs(symmetric, pure, recording_support::field("K", k)));
    record("symmetric-helmholtz", "association_helmholtz", symmetric_result.association,
           expected_association, inputs(symmetric, pure, recording_support::field("K", k)));

    Parameters unequal = symmetric;
    unequal.site_multiplicity = {2, 3};
    const auto unequal_fractions = epcsaft::solve_site_fractions(
        temperature, density, pure, unequal, sigma, diameter
    );
    const auto unequal_result = evaluate(unequal, temperature, density, pure);
    const double discriminant = std::pow(1.0 + k * 5.0, 2) - 24.0 * k * k;
    const double bonds = 12.0 * k / (1.0 + 5.0 * k + std::sqrt(discriminant));
    const double expected_a = 1.0 - bonds / 2.0, expected_b = 1.0 - bonds / 3.0;
    const double expected_unequal_association = 2.0 * (
        std::log(expected_a) - expected_a / 2.0 + 0.5
    ) + 3.0 * (std::log(expected_b) - expected_b / 2.0 + 0.5);
    record("unequal-x-a", "site_fraction", unequal_fractions[0], expected_a,
           inputs(unequal, pure, recording_support::field("bonded_count", bonds)));
    record("unequal-x-b", "site_fraction", unequal_fractions[1], expected_b,
           inputs(unequal, pure, recording_support::field("bonded_count", bonds)));
    record("unequal-conservation", "bonded_count", 2.0 * (1.0 - unequal_fractions[0]),
           3.0 * (1.0 - unequal_fractions[1]),
           inputs(unequal, pure, recording_support::field("bonded_count", bonds)));
    record("unequal-helmholtz", "association_helmholtz", unequal_result.association,
           expected_unequal_association,
           inputs(unequal, pure, recording_support::field("bonded_count", bonds)));
    Parameters relabelled = unequal;
    relabelled.site_multiplicity = {3, 2};
    const auto relabelled_fractions = epcsaft::solve_site_fractions(
        temperature, density, pure, relabelled, sigma, diameter
    );
    const auto relabelled_result = evaluate(relabelled, temperature, density, pure);
    record("unequal-relabel-x", "site_fraction", relabelled_fractions[0], unequal_fractions[1],
           inputs(relabelled, pure, "site_order:BA"));
    record("unequal-relabel-helmholtz", "association_helmholtz",
           relabelled_result.association, unequal_result.association,
           inputs(relabelled, pure, "site_order:BA"));

    Parameters zero_strength = symmetric;
    zero_strength.association_energy_over_k = {0.0, 0.0, 0.0, 0.0};
    const auto zero_strength_fractions = epcsaft::solve_site_fractions(
        temperature, density, pure, zero_strength, sigma, diameter
    );
    const auto zero_strength_result = evaluate(zero_strength, temperature, density, pure);
    record("limit-zero-strength-x", "site_fraction", zero_strength_fractions[0], 1.0,
           inputs(zero_strength, pure, "limit:zero-strength"));
    record("limit-zero-strength-energy", "association_helmholtz",
           zero_strength_result.association, 0.0,
           inputs(zero_strength, pure, "limit:zero-strength"));
    // The smallest positive density: zero itself is outside the fluid domain.
    const double zero_density = std::numeric_limits<double>::min();
    const auto zero_density_fractions = epcsaft::solve_site_fractions(
        temperature, zero_density, pure, symmetric, sigma, diameter
    );
    record("limit-zero-density-x", "site_fraction", zero_density_fractions[0], 1.0,
           recording_support::fields({model_inputs(symmetric, temperature, zero_density, pure),
               recording_support::field("l_ij", symmetric.l_ij), "limit:zero-density"}));
    const auto zero_density_terms = epcsaft::residual_helmholtz(
        temperature, zero_density, pure, symmetric, sigma, diameter, zero_density_fractions,
        epcsaft::association_products(zero_density, pure, symmetric, zero_density_fractions, epcsaft::association_delta(temperature, zero_density, pure, symmetric, sigma, diameter)));
    record("limit-zero-density-energy", "association_helmholtz",
           zero_density_terms.association, 0.0,
           recording_support::fields({model_inputs(symmetric, temperature, zero_density, pure),
               recording_support::field("l_ij", symmetric.l_ij), "limit:zero-density"}));
    Parameters absent = symmetric;
    absent.site_component = {1};
    absent.site_multiplicity = {1};
    absent.association_active = {1};
    absent.association_energy_over_k = {2425.7};
    absent.association_volume = {0.0451};
    const std::vector<double> absent_x = {1.0, 0.0};
    const auto absent_sigma = epcsaft::segment_sigma(temperature, absent);
    const auto absent_diameter = epcsaft::segment_diameter(temperature, absent, absent_sigma);
    const auto absent_fractions = epcsaft::solve_site_fractions(
        temperature, density, absent_x, absent, absent_sigma, absent_diameter
    );
    const auto absent_result = evaluate(absent, temperature, density, absent_x);
    record("limit-absent-self-x", "site_fraction", absent_fractions[0], 1.0,
           inputs(absent, absent_x, "limit:absent-self-associated-component"));
    record("limit-absent-energy", "association_helmholtz", absent_result.association, 0.0,
           inputs(absent, absent_x, "limit:absent-self-associated-component"));

    Parameters cross = symmetric;
    cross.site_component = {0, 1};
    cross.site_multiplicity = {1, 2};
    cross.association_energy_over_k = {0.0, 2400.0, 2400.0, 0.0};
    cross.association_volume = {0.0, 0.04, 0.04, 0.0};
    cross.l_ij = {0.0, 0.12, 0.12, 0.0};
    const double cross_temperature = 320.0, cross_density = 9.0e3;
    const std::vector<double> cross_x = {0.25, 0.75};
    const auto cross_sigma = epcsaft::segment_sigma(cross_temperature, cross);
    const auto cross_diameter = epcsaft::segment_diameter(cross_temperature, cross, cross_sigma);
    const auto cross_delta = epcsaft::association_delta(
        cross_temperature, cross_density, cross_x, cross, cross_sigma, cross_diameter
    );
    const double rho_n = cross_density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    const double zeta_2 = epcsaft::PI / 6.0 * rho_n * (
        cross_x[0] * cross.segment_count[0] * std::pow(cross_diameter[0], 2)
        + cross_x[1] * cross.segment_count[1] * std::pow(cross_diameter[1], 2)
    );
    const double eta = epcsaft::PI / 6.0 * rho_n * (
        cross_x[0] * cross.segment_count[0] * std::pow(cross_diameter[0], 3)
        + cross_x[1] * cross.segment_count[1] * std::pow(cross_diameter[1], 3)
    );
    const double reduced = cross_diameter[0] * cross_diameter[1]
        / (cross_diameter[0] + cross_diameter[1]);
    const double contact = 1.0 / (1.0 - eta)
        + 3.0 * reduced * zeta_2 / std::pow(1.0 - eta, 2)
        + 2.0 * reduced * reduced * zeta_2 * zeta_2 / std::pow(1.0 - eta, 3);
    const double sigma_01 = 0.5 * (cross_sigma[0] + cross_sigma[1]) * (1.0 - cross.l_ij[1]);
    const auto expected_delta = [&](std::size_t pair) {
        return contact * (std::exp(cross.association_energy_over_k[pair] / cross_temperature) - 1.0)
            * std::pow(sigma_01, 3) * cross.association_volume[pair];
    };
    const auto cross_fractions = epcsaft::solve_site_fractions(
        cross_temperature, cross_density, cross_x, cross, cross_sigma, cross_diameter
    );
    const double cross_residual_a = cross_fractions[0] * (
        1.0 + rho_n * cross_x[1] * cross.site_multiplicity[1]
        * cross_fractions[1] * expected_delta(1)
    ) - 1.0;
    const double cross_residual_b = cross_fractions[1] * (
        1.0 + rho_n * cross_x[0] * cross.site_multiplicity[0]
        * cross_fractions[0] * expected_delta(2)
    ) - 1.0;
    const auto cross_inputs = recording_support::fields({
        model_inputs(cross, cross_temperature, cross_density, cross_x),
        recording_support::field("l_ij", cross.l_ij), recording_support::field("contact", contact),
        recording_support::field("sigma_01", sigma_01),
        recording_support::field("site_fractions", cross_fractions),
    });
    record("cross-strength-ab", "Delta_AB", cross_delta[1], expected_delta(1), cross_inputs, 1.0e-12,
           "angstrom^3");
    record("cross-strength-ba", "Delta_BA", cross_delta[2], expected_delta(2), cross_inputs, 1.0e-12,
           "angstrom^3");
    record("cross-residual-a", "association_residual", cross_residual_a, 0.0, cross_inputs, epcsaft::policy::ASSOCIATION_RESIDUAL.value);
    record("cross-residual-b", "association_residual", cross_residual_b, 0.0, cross_inputs, epcsaft::policy::ASSOCIATION_RESIDUAL.value);

    const std::vector<double> one_sided = {0.0, 1.0};
    const auto one_sided_delta = epcsaft::association_delta(
        cross_temperature, cross_density, one_sided, cross, cross_sigma, cross_diameter);
    epcsaft::core::Convergence one_sided_record;
    const auto pair = epcsaft::classify_association(cross, [](std::size_t) { return false; });
    const auto one_sided_fractions = epcsaft::solve_site_fractions(
        cross_temperature, cross_density, one_sided, cross, cross_sigma, cross_diameter, &one_sided_record,
        nullptr, &one_sided_delta, false, {}, nullptr, nullptr, nullptr, nullptr, &pair);
    CHECK(one_sided_record.iterations == 0);
    const double q = rho_n * cross.site_multiplicity[1] * one_sided_delta[1];
    const auto one_sided_inputs = recording_support::fields({
        model_inputs(cross, cross_temperature, cross_density, one_sided),
        recording_support::field("K_AB", q)});
    record("limit-absent-cross-x-a", "site_fraction", one_sided_fractions[0], 1.0 / (1.0 + q), one_sided_inputs);
    record("limit-absent-cross-x-b", "site_fraction", one_sided_fractions[1], 1.0, one_sided_inputs);
}

void kraska_model_iv() {
    using namespace epcsaft;
    auto p = two_site_binary();
    // Use the shuffled order (u, v, v, u), and make a, b, c, d = 1/2, 3, 2, 8.
    p.site_component = {0, 1, 1, 0};
    p.site_multiplicity = {1, 1, 1, 1};
    p.association_active = {0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
    p.association_energy_over_k.assign(16, 0);
    p.association_volume.assign(16, 0);
    const auto plan = classify_association(p, [](std::size_t) { return false; });
    const auto broken = classify_association(p, [](std::size_t edge) { return edge == 1; });
    CHECK(broken.strategy == AssociationPlan::General);
    const auto mismatch = [](Parameters copy) {
        copy.temperature_law_slot.push_back({static_cast<int>(TemperatureLaw::AssociationEnergyOverK), 1, 0});
        return classify_association(copy, [](std::size_t) { return false; });
    };
    CHECK(mismatch(p).strategy == AssociationPlan::General);
    auto duplicate = p;
    const auto energy = static_cast<int>(TemperatureLaw::AssociationEnergyOverK);
    duplicate.temperature_law_slot = {{energy, 1, 0}, {energy, 1, 1}, {energy, 11, 1}};
    duplicate.temperature_law_reference_k = {300.0, 300.0};
    duplicate.temperature_law_term_law = {0, 1};
    duplicate.temperature_law_term_power = {0, 0};
    duplicate.temperature_law_coefficient = {1100.0, 2100.0};
    CHECK(classify_association(duplicate, [](std::size_t) { return false; }).strategy == AssociationPlan::General);
    constexpr double t = 300, rho = 10000;
    const double density = rho * NUMBER_DENSITY_PER_ANGSTROM3;
    const std::vector<double> seed{0.56, 0.31, 0.28, 0.52}, reference{0.5, 0.25, 0.25, 0.5};
    const auto sigma = segment_sigma(t, p), diameter = segment_diameter(t, p, sigma);
    const auto compare = [](const std::string& id, double value, double reference, double tolerance) {
        recording_support::compare(id, "mass_action", "1", value, reference,
            recording_support::Rule::Absolute, tolerance, 0.0, "Kraska Model IV Eqs. 15–18");
    };
    // In the middle case ad = bc, so the quartic coefficient cancels.
    for (const auto& [name, a, b, c, d] : {std::tuple{"quartic", 0.5, 3.0, 2.0, 8.0},
             std::tuple{"degree-drop", 0.5, 3.0, 1.5, 9.0}, std::tuple{"zero-cross", 2.0, 0.0, 0.0, 12.0}}) {
        const double x0 = b != 0 ? c / (b + c) : 0.4;
        const std::vector<double> x{x0, 1 - x0};
        std::vector<double> strength(16);
        const auto bond = [&](int i, int j, double delta) { strength[4*i+j] = strength[4*j+i] = delta; };
        bond(0, 3, a / (density * x[0])); bond(1, 2, d / (density * x[1]));
        bond(0, 1, b / (density * x[1])); bond(2, 3, b / (density * x[1]));
        std::optional<AssociationAction> accepted;
        const auto solve = [&](const AssociationPlan* strategy, core::Convergence* work = nullptr) {
            return solve_site_fractions(t, rho, x, p, sigma, diameter, work, &seed, &strength,
                false, {}, nullptr, nullptr, nullptr, &accepted, strategy);
        };
        const auto general = solve(nullptr), reduced = solve(&plan);
        CHECK(accepted && accepted->strategy == AssociationPlan::BalancedBinary);
        const auto rows = association_residuals(rho, x, p, reduced, strength);
        const std::string id = std::string(name) + "/";
        for (std::size_t site = 0; site < seed.size(); ++site) {
            compare(id + "X/" + std::to_string(site), reduced[site], reference[site], 1e-12);
            compare(id + "Newton/" + std::to_string(site), reduced[site], general[site], 1e-12);
            compare(id + "residual/" + std::to_string(site), rows[site], 0.0, policy::ASSOCIATION_RESIDUAL.value);
        }
        if (b != 0) {
            const double v = reduced[1], r = 1 - v - d*v*v, z = c*v;
            compare(id + "polynomial", a*r*r + (1+b*v)*z*r - z*z, 0.0, 1e-12);
        }
        if (d == 8.0) {
            // Violate the supplied strength's declared symmetry to exercise the original-row check.
            strength[1] *= 1.1;
            core::Convergence general_work, fallback_work;
            const auto expected = solve(nullptr, &general_work), rejected = solve(&plan, &fallback_work);
            CHECK(!accepted);
            std::cout << "rejected-root evaluations: general=" << general_work.evaluations
                      << " fallback=" << fallback_work.evaluations << '\n';
            CHECK(fallback_work.evaluations == general_work.evaluations + 1);
            for (std::size_t site = 0; site < seed.size(); ++site)
                compare("rejected-root/" + std::to_string(site), rejected[site], expected[site], 0.0);
        }
    }
}

void huang_radosz_table_vii() {
    using namespace epcsaft;
    constexpr double t = 300, rho = 10000;
    const double delta = 1.0 / (rho * NUMBER_DENSITY_PER_ANGSTROM3);
    const auto p0 = two_site_binary();
    const auto sigma = segment_sigma(t, p0), diameter = segment_diameter(t, p0, sigma);
    const double golden = (std::sqrt(5.0) - 1) / 2.0;
    const double minority_4b = 2.0 / (3.0 + std::sqrt(13.0));
    const auto check = [&](const std::string& id, std::vector<int> m, const std::vector<double>& fractions) {
        auto p = p0;
        p.site_component.assign(m.size(), 0);
        p.site_multiplicity = std::move(m);
        p.association_active = fractions.size() == 1 ? std::vector<int>{1} : std::vector<int>{0, 1, 1, 0};
        std::vector<double> strength;
        for (int edge : p.association_active) strength.push_back(edge * delta);
        const auto plan = classify_association(p, [](std::size_t) { return false; });
        const auto got = solve_site_fractions(t, rho, {1.0, 0.0}, p, sigma, diameter, nullptr, nullptr,
            &strength, false, {}, nullptr, nullptr, nullptr, nullptr, &plan);
        for (std::size_t i = 0; i < got.size(); ++i)
            recording_support::compare(id + "/" + std::to_string(i), "site_fraction", "1", got[i], fractions[i],
                recording_support::Rule::Absolute, 1e-12, 0.0, "Huang–Radosz Table VII");
    };
    check("1A", {1}, {golden});
    check("2B", {1, 1}, {golden, golden});
    check("2A", {2}, {0.5});
    check("4C", {2, 2}, {0.5, 0.5});
    check("3B", {2, 1}, {1.0 / std::sqrt(2.0), std::sqrt(2.0) - 1.0});
    check("4B", {3, 1}, {1.0 / (1.0 + minority_4b), minority_4b});
}

// A term that is switched off contributes nothing, whatever else is in the
// parameter record.
void inactive_terms_contribute_nothing() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    const auto result = evaluate(p, 300.0, 1.2e4, x);
    const auto inputs = model_inputs(p, 300.0, 1.2e4, x);
    recording_support::compare("association", "association", "1", result.association, 0.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 0.0, "inactive-term", inputs);

    const double born = epcsaft::born_contribution(300.0, x, p, 78.4);
    recording_support::compare("born", "born", "1", born, 0.0, recording_support::Rule::ReferenceOrOne, 0.0, 0.0,
        "inactive-term", inputs);

    const auto electrolyte = epcsaft::electrolyte_contribution(
        300.0, 1.2e4, x, {0.4, 0.6}, p
    );
    recording_support::compare("debye-huckel", "debye_huckel", "1", electrolyte.screening, 0.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 0.0, "inactive-term", inputs);
}

// Bulow 2019 Eq. (7): a declared water component remains active at zero amount.
void component_permittivity_zero_solvent() {
    using AD = CppAD::AD<double>;
    auto p = electrolyte_fixture();
    p.charge_number = {0, 1, -1};
    p.relative_permittivity = {80.0, 11.0, 11.0};
    p.permittivity_term_component.clear();
    p.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    const auto record = [](const char* id, double value, double reference,
                           const std::string& inputs, double tolerance) {
        recording_support::compare(id, id, "1", value, reference, recording_support::Rule::Absolute, tolerance, 0.0,
            "exact_derivative_or_difference", inputs);
    };
    const auto value = [&](double t) {
        return epcsaft::bulk_relative_permittivity(298.15, 1000.0,
            std::vector<double>{t, (1.0 - t) / 2.0, (1.0 - t) / 2.0}, p);
    };
    for (double t : {0.0, 1.0e-3, 1.0e-2}) {
        std::vector<AD> u{t};
        CppAD::Independent(u);
        const std::vector<AD> x{u[0], (1.0 - u[0]) / 2.0, (1.0 - u[0]) / 2.0};
        CppAD::ADFun<double> tape(u, {epcsaft::bulk_relative_permittivity(AD(298.15), AD(1000.0), x, p)});
        const std::vector<double> path_point{t};
        const double first = tape.Jacobian(path_point)[0], second = tape.Hessian(path_point, 0)[0];
        const auto point = recording_support::field("t", t);
        record("path-first", first, 69.0, point, t == 0.0 ? 1.0e-12 : 1.0e-8);
        record("path-second", second, 0.0, point, t == 0.0 ? 1.0e-12 : 1.0e-6);
        for (double h : {t == 0.0 ? 1.0e-3 : t / 4.0, t == 0.0 ? 5.0e-4 : t / 8.0}) {
            const double f = value(t), plus = value(t + h), twice = value(t + 2.0 * h);
            const double other = value(t == 0.0 ? 3.0 * h : t - h);
            const double d1 = t == 0.0 ? (-3.0 * f + 4.0 * plus - twice) / (2.0 * h)
                                         : (plus - other) / (2.0 * h);
            const double d2 = t == 0.0 ? (2.0 * f - 5.0 * plus + 4.0 * twice - other) / (h * h)
                                         : (plus - 2.0 * f + other) / (h * h);
            const auto inputs = recording_support::fields({point, recording_support::field("h", h),
                recording_support::field("f", f), recording_support::field("f_plus", plus),
                recording_support::field("f_twice", twice), recording_support::field("f_other", other)});
            record("difference-first", d1, first, inputs, 1.0e-8);
            record("difference-second", d2, second, inputs, 1.0e-4);
        }
        const std::vector<double> point_x{t, (1.0 - t) / 2.0, (1.0 - t) / 2.0};
        for (bool amounts : {false, true}) {
            std::vector<AD> n(point_x.begin(), point_x.end());
            CppAD::Independent(n);
            const AD total = amounts ? n[0] + n[1] + n[2] : AD(1.0);
            const std::vector<AD> composition{n[0] / total, n[1] / total, n[2] / total};
            CppAD::ADFun<double> ambient(n, {epcsaft::bulk_relative_permittivity(
                AD(298.15), AD(1000.0), composition, p)});
            const auto gradient = ambient.Jacobian(point_x), hessian = ambient.Hessian(point_x, 0);
            const std::vector<double> a{80.0, 11.0, 11.0};
            for (std::size_t i = 0; i < 3; ++i) {
                const auto inputs = recording_support::fields({point,
                    recording_support::field("amounts", amounts), recording_support::field("i", i)});
                record("ambient-first", gradient[i], amounts ? a[i] - value(t) : a[i], inputs, 1.0e-8);
                for (std::size_t j = 0; j < 3; ++j)
                    record("ambient-second", hessian[3 * i + j], amounts ? 2.0 * value(t) - a[i] - a[j] : 0.0,
                           recording_support::fields({inputs, recording_support::field("j", j)}), 1.0e-6);
            }
        }
    }
}

// Bülow 2021 Eq. 19 and its explicitly solvent-normalized alternative.
void solvent_mole_permittivity() {
    using AD = CppAD::AD<double>;
    auto p = electrolyte_fixture();
    p.charge_number = {0, 0, 1, -1};
    p.permittivity_form = epcsaft::PermittivityForm::SolventMoleWeighting;
    p.permittivity_solvent = {1, 1, 0, 0};
    p.relative_permittivity = {80, 20, 0, 0};
    p.permittivity_reference_temperature_k = {300, 300, 0, 0};
    p.permittivity_term_component = {0, 0, 0, 1, 1, 1};
    p.permittivity_term_power = {1, 2, 3, 1, 2, 3};
    p.permittivity_term_coefficient = {-0.2, 0.001, 0.00001, -0.05, 0.0005, -0.000002};
    p.born_c_shell = p.born_c_dielectric = 0;
    p.born_diameter_angstrom = p.sigma_constant;
    auto check = [](const char* name, double value, double expected, double tolerance=1e-12) {
        recording_support::compare(name, name, "1", value, expected, recording_support::Rule::Absolute, tolerance, 0,
            "analytic", "T:300;x:0.7,0.1,0.1,0.1");
    };
    for (bool normalized : {false, true}) {
        p.permittivity_solvent_normalized = normalized;
        const double expected = normalized ? 72.5 : 58;
        const std::vector<double> x{0.7,0.1,0.1,0.1};
        for (double mass : {0.018,0.061}) for (double rho : {1000.,20000.}) {
            p.molar_mass_kg_per_mol = {mass,0.08-mass,0.023,0.035};
            check("value-mass-density",epcsaft::bulk_relative_permittivity(300.,rho,x,p),expected);
        }
        std::vector<AD> u{300,0.7,0.1,0.1,0.1,1000}; CppAD::Independent(u);
        std::vector<AD> ax{u[1],u[2],u[3],u[4]};
        const AD eps = epcsaft::bulk_relative_permittivity(u[0],u[5],ax,p);
        const AD born = epcsaft::born_contribution(u[0],ax,p,eps);
        AD radial=0; for (std::size_t i=2;i<4;++i) radial+=ax[i]/(p.sigma_constant[i]*epcsaft::ANGSTROM_M);
        const AD classic=-epcsaft::ELEMENTARY_CHARGE_C*epcsaft::ELEMENTARY_CHARGE_C
            /(4*epcsaft::PI*epcsaft::VACUUM_PERMITTIVITY_F_PER_M*epcsaft::BOLTZMANN_J_PER_K*u[0])*(1-1/eps)*radial;
        CppAD::ADFun<double> f(u,{eps,born-classic});
        const std::vector<double> point{300,0.7,0.1,0.1,0.1,1000};
        const auto j=f.Jacobian(point),h=f.Hessian(point,0),bj=f.Jacobian(point),bh=f.Hessian(point,1);
        check("thermal-first",j[0],normalized?-0.18125:-0.145);
        check("thermal-second",h[0],normalized?0.001875:0.0015);
        check("mixed-T-water",h[1],normalized?-0.0234375:-0.2);
        check("water-first",j[1],normalized?9.375:80);
        check("amine-first",j[2],normalized?-65.625:20);
        check("water-second",h[7],normalized?-23.4375:0);
        check("water-amine",h[8],normalized?70.3125:0);
        for (std::size_t i=0;i<6;++i) {check("density-mixed",h[30+i],0);check("classic-Born-gradient",bj[6+i],0,1e-10);}
        for (double v:bh) check("classic-Born-hessian",v,0,1e-10);
        check("classic-Born-value",f.Forward(0,point)[1],0,1e-10);
        auto derivatives=[&](double lambda) {
            auto z=point;z[1]+=lambda;z[3]-=lambda/2;z[4]-=lambda/2;
            std::vector<double> d(6,0);d[1]=1;d[3]=d[4]=-0.5;
            std::vector<double> result{f.Forward(0,z)[0]};double factorial=1;
            for (std::size_t k=1;k<=4;++k) {factorial*=k;result.push_back(f.Forward(k,k==1?d:std::vector<double>(6,0))[0]*factorial);}
            return result;
        };
        const auto ds=derivatives(0);
        check("composition-third",ds[3],normalized?87.890625:0,1e-10);
        check("composition-fourth",ds[4],normalized?-439.453125:0,1e-9);
        for (std::size_t k:{3U,4U}) {
            const double step=0.00025, half=step/2;
            const double coarse=(derivatives(step)[k-1]-derivatives(-step)[k-1])/(2*step);
            const double fine=(derivatives(half)[k-1]-derivatives(-half)[k-1])/(2*half);
            const double extrapolated=(4*fine-coarse)/3;
            const double refinement=std::abs(fine-coarse)/3;
            const double target=1e-6*std::max(1.,std::abs(ds[k]));
            std::cout << "contraction-refinement normalized=" << normalized << " order=" << k
                      << " h=" << step << " half=" << half << " coarse=" << std::setprecision(17) << coarse
                      << " fine=" << fine << " extrapolated=" << extrapolated
                      << " error=" << std::abs(extrapolated-ds[k]) << " refinement=" << refinement
                      << " target=" << target << std::endl;
            CHECK(refinement<=target/4);
            check("composition-Richardson",extrapolated,ds[k],target);
        }
        check("pure-water",epcsaft::bulk_relative_permittivity(300.,1000.,std::vector<double>{1,0,0,0},p),80);
        check("pure-amine",epcsaft::bulk_relative_permittivity(300.,1000.,std::vector<double>{0,1,0,0},p),20);
        check("salt-dilution",epcsaft::bulk_relative_permittivity(300.,1000.,std::vector<double>{0.35,0.05,0.3,0.3},p),normalized?72.5:29);
        CHECK((recording_support::refusal([&] {check("zero-solvent",epcsaft::bulk_relative_permittivity(300.,1000.,std::vector<double>{0,0,0.5,0.5},p),0);})
            .find("solvent denominator")!=std::string::npos)==normalized);
    }
}

// These references use only the accepted scalar equations, including numerical
// radial and finite-size integrals; they do not route through the contribution owners.
void electrolyte_born_and_dielectric_match_independent_equations() {
    solvent_mole_permittivity();
    component_permittivity_zero_solvent();
    const double temperature = 320.0, density = 9.0e3;
    const std::vector<double> x = {0.4, 0.3, 0.2, 0.1};
    Parameters p = electrolyte_fixture();
    const auto inputs = [&](const Parameters& parameters, double rho,
                            const std::vector<double>& composition,
                            const std::string& extra) {
        return recording_support::fields({
            model_inputs(parameters, temperature, rho, composition),
            recording_support::field("packing_diameter_A", parameters.packing_diameter_angstrom),
            recording_support::field("dh_diameter_A", parameters.debye_huckel_diameter_angstrom),
            recording_support::field("born_diameter_A", parameters.born_diameter_angstrom),
            recording_support::field("solvation", parameters.solvation_factor),
            recording_support::field("epsilon_r", parameters.relative_permittivity),
            recording_support::field("permittivity_term_component", parameters.permittivity_term_component),
            recording_support::field("permittivity_term_coefficient", parameters.permittivity_term_coefficient),
            recording_support::field("permittivity_term_power", parameters.permittivity_term_power),
            recording_support::field("permittivity_reference_temperature_K", parameters.permittivity_reference_temperature_k),
            recording_support::field("ion_suppression", parameters.ion_suppression_coefficient),
            recording_support::field("dielectric_volume_m3_per_mol", parameters.dielectric_volume_m3_per_mol),
            recording_support::field("dielectric_temperature_K", parameters.dielectric_temperature_k),
            recording_support::field("dipole_scaling_debye2", parameters.dipole_scaling_debye2),
            recording_support::field("polarizability_scaling_A3", parameters.polarizability_scaling_angstrom3),
            recording_support::field("correlation_integral_parameter", parameters.correlation_integral_parameter),
            recording_support::field("epsilon_inner", parameters.ionic_region_relative_permittivity),
            recording_support::field("born_c_shell", parameters.born_c_shell),
            recording_support::field("born_c_dielectric", parameters.born_c_dielectric),
            recording_support::field("permittivity_form", static_cast<int>(parameters.permittivity_form.require("permittivity_form"))),
            recording_support::field("electrolyte_form", static_cast<int>(parameters.electrolyte_form.require("electrolyte_form"))),
            extra,
        });
    };
    const auto record = [&](const char* id, const char* quantity, const char* unit,
                            double model, double reference, const std::string& values, double tolerance) {
        recording_support::compare(id, quantity, unit, model, reference, recording_support::Rule::ReferenceOrOne, 0.0,
            tolerance, "independent_equation_evaluation", values);
    };

    const auto chi_integral = [](double reduced, int intervals) {
        return 3.0 * simpson(
            [reduced](double t) { return t * t / (1.0 + reduced * t); },
            0.0, 1.0, intervals
        );
    };
    for (const double reduced : {0.0, 5.0e-4, 2.0e-3, 0.5}) {
        const bool use_series = reduced < epcsaft::policy::SCREENING_SERIES_SWITCH.value;
        const double coarse = chi_integral(reduced, 1000);
        const double reference = chi_integral(reduced, 4000);
        const double cancellation_bound = use_series ? 1.0e-12
            : std::max(1.0e-12, 64.0 * std::numeric_limits<double>::epsilon()
                / (reduced * reduced * reduced));
        const std::string point = recording_support::number(reduced);
        record((std::string("chi-value-") + point).c_str(), "chi", "1",
               epcsaft::finite_size_chi(reduced, use_series), reference,
               inputs(p, density, x, recording_support::fields({
                   recording_support::field("reduced_diameter", reduced),
                   recording_support::field("series_branch", static_cast<int>(use_series)),
                   recording_support::field("reference_panels", 4000),
                   recording_support::field("cancellation_bound", cancellation_bound),
               })), cancellation_bound);
        record((std::string("chi-refine-") + point).c_str(), "chi_integral", "1",
               coarse, reference, inputs(p, density, x, recording_support::fields({
                   recording_support::field("reduced_diameter", reduced),
                   recording_support::field("coarse_panels", 1000),
                   recording_support::field("reference_panels", 4000),
               })), 1.0e-11);
    }

    // Taped chi' = -3 int t^3/(1+yt)^2 and chi'' = 6 int t^4/(1+yt)^3 on both sides of the
    // series switch; the former closed form missed chi'' by 3.3x at 1e-3 and 2e-3 at 2e-3.
    for (const double reduced : {2.0e-3, 0.19, 0.21}) {
        const bool use_series = reduced < epcsaft::policy::SCREENING_SERIES_SWITCH.value;
        std::vector<CppAD::AD<double>> y = {reduced};
        CppAD::Independent(y);
        std::vector<CppAD::AD<double>> chi = {epcsaft::finite_size_chi(y[0], use_series)};
        CppAD::ADFun<double> tape(y, chi);
        tape.Forward(0, std::vector<double>{reduced});
        const double first = tape.Forward(1, std::vector<double>{1.0})[0];
        const double second = 2.0 * tape.Forward(2, std::vector<double>{0.0})[0];
        const auto integral = [reduced](int power, int order) {
            return simpson([=](double t) { return std::pow(t, power) / std::pow(1.0 + reduced * t, order); }, 0.0, 1.0, 4000);
        };
        const std::string point = recording_support::number(reduced);
        const auto fields = inputs(p, density, x, recording_support::fields({
            recording_support::field("reduced_diameter", reduced),
            recording_support::field("series_branch", static_cast<int>(use_series)),
            recording_support::field("reference_panels", 4000)}));
        record((std::string("chi-d1-") + point).c_str(), "dchi_dy", "1", first, -3.0 * integral(3, 2), fields, 1.0e-11);
        record((std::string("chi-d2-") + point).c_str(), "d2chi_dy2", "1", second, 6.0 * integral(4, 3), fields, 1.0e-9);
    }

    const auto component_permittivity = epcsaft::component_relative_permittivities(temperature, p);
    const double epsilon_0 = 70.0 + 0.02 * (temperature - 300.0);
    const double epsilon_1 = 55.0 - 0.01 * std::pow(temperature - 300.0, 2);
    record("dielectric-temperature-0", "epsilon_r", "1", component_permittivity[0], epsilon_0,
           inputs(p, density, x, "component:0"), 1.0e-12);
    record("dielectric-temperature-1", "epsilon_r", "1", component_permittivity[1], epsilon_1,
           inputs(p, density, x, "component:1"), 1.0e-12);
    const double solvent_mass = x[0] * p.molar_mass_kg_per_mol[0]
        + x[1] * p.molar_mass_kg_per_mol[1];
    const double salt_free = (x[0] * p.molar_mass_kg_per_mol[0] * epsilon_0
        + x[1] * p.molar_mass_kg_per_mol[1] * epsilon_1) / solvent_mass;
    const double expected_ion_suppression = salt_free
        / (1.0 + x[2] * p.ion_suppression_coefficient[2]
           + x[3] * p.ion_suppression_coefficient[3]);
    record("dielectric-ion-specific", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, x, p),
           expected_ion_suppression, inputs(p, density, x, "form:ion-specific"), 1.0e-12);
    const std::vector<double> single_solvent_x = {0.7, 0.0, 0.2, 0.1};
    const double expected_single = epsilon_0 / (1.0 + 0.2 * 5.0 + 0.1 * 8.0);
    record("dielectric-single-ion-specific", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, single_solvent_x, p),
           expected_single, inputs(p, density, single_solvent_x, "form:Zuber-single-solvent"), 1.0e-12);
    Parameters equal_suppression = p;
    equal_suppression.ion_suppression_coefficient[2] = 6.0;
    equal_suppression.ion_suppression_coefficient[3] = 6.0;
    record("dielectric-equal-suppression", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, x, equal_suppression),
           salt_free / (1.0 + 6.0 * (x[2] + x[3])),
           inputs(equal_suppression, density, x, "form:equal-ion-fraction"), 1.0e-12);
    Parameters component_mixing = p;
    component_mixing.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    record("dielectric-component-mixing", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, x, component_mixing),
           salt_free * (x[0] + x[1]) + x[2] * 12.0 + x[3] * 18.0,
           inputs(component_mixing, density, x, "form:component-mixing"), 1.0e-12);
    const std::vector<double> salt_free_x = {4.0 / 7.0, 3.0 / 7.0, 0.0, 0.0};
    record("dielectric-ion-salt-free", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, salt_free_x, p), salt_free,
           inputs(p, density, salt_free_x, "limit:salt-free-ion-suppression"), 1.0e-12);
    record("dielectric-component-salt-free", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(
               temperature, density, salt_free_x, component_mixing
           ), salt_free, inputs(component_mixing, density, salt_free_x,
                                "limit:salt-free-component-mixing"), 1.0e-12);
    Parameters density_temperature = p;
    density_temperature.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    const double mixed_volume = (x[0] * 1.5e-4 + x[1] * 2.5e-4) / (x[0] + x[1]);
    const double mixed_temperature = (x[0] * 450.0 + x[1] * 600.0) / (x[0] + x[1]);
    const double expected_density = 1.0 + density * (x[0] + x[1]) * mixed_volume
        * (mixed_temperature / temperature - 1.0);
    record("dielectric-density-temperature", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, x, density_temperature),
           expected_density, inputs(density_temperature, density, x, "form:density-temperature"), 1.0e-12);
    record("dielectric-density-dilute", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, 1.0e-6, x, density_temperature), 1.0,
           inputs(density_temperature, 1.0e-6, x, "limit:dilute"), 1.0e-10);
    Parameters dipolar = p;
    dipolar.permittivity_form = epcsaft::PermittivityForm::DipolarPerturbation;
    const auto dipolar_reference = [&](double rho) {
        const double number_density = rho * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
        double sum = 0.0, correlation = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            // Debye^2 in J A^3 from the SI Debye, epsilon_0 and the angstrom.
            sum += x[i] * (dipolar.dipole_scaling_debye2[i] * std::pow(epcsaft::DEBYE_C_M, 2)
                / (4.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M * 1.0e-30)
                / (epcsaft::BOLTZMANN_J_PER_K * temperature)
                + 3.0 * dipolar.polarizability_scaling_angstrom3[i]);
            correlation += x[i] * dipolar.correlation_integral_parameter[i];
        }
        const double y = 4.0 * epcsaft::PI / 9.0 * number_density * sum;
        const double integral = 1.0 + correlation * (std::exp(-y) - 1.0);
        return 1.0 + 3.0 * y * (1.0 + y + ((17.0 / 16.0) * integral - 1.0) * y * y);
    };
    record("dielectric-dipolar", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, density, x, dipolar),
           dipolar_reference(density), inputs(dipolar, density, x, "form:dipolar"), 1.0e-12);
    record("dielectric-dipolar-dilute", "epsilon_bulk", "1",
           epcsaft::bulk_relative_permittivity(temperature, 1.0e-6, x, dipolar), 1.0,
           inputs(dipolar, 1.0e-6, x, "limit:dilute"), 1.0e-10);

    const auto electrolyte = epcsaft::electrolyte_contribution(
        temperature, density, x, x, p
    );
    record("dh-bulk", "epsilon_bulk", "1", electrolyte.bulk_relative_permittivity,
           expected_ion_suppression, inputs(p, density, x, "case:multivalent"), 1.0e-12);
    const double charge_square_sum = x[2] + 4.0 * x[3];
    const double kappa = std::sqrt(density * epcsaft::AVOGADRO_PER_MOL
        * epcsaft::ELEMENTARY_CHARGE_C * epcsaft::ELEMENTARY_CHARGE_C * charge_square_sum
        / (epcsaft::VACUUM_PERMITTIVITY_F_PER_M * expected_ion_suppression
           * epcsaft::BOLTZMANN_J_PER_K * temperature));
    const double chi_cation_argument =
        kappa * p.debye_huckel_diameter_angstrom[2] * epcsaft::ANGSTROM_M;
    const double chi_anion_argument =
        kappa * p.debye_huckel_diameter_angstrom[3] * epcsaft::ANGSTROM_M;
    const double chi_quadrature_bound = std::max(
        0.4 * chi_cation_argument * chi_cation_argument / std::pow(4000.0, 4),
        0.4 * chi_anion_argument * chi_anion_argument / std::pow(4000.0, 4)
    );
    const double weighted_chi = x[2] * chi_integral(chi_cation_argument, 4000)
        + 4.0 * x[3] * chi_integral(chi_anion_argument, 4000);
    const double expected_dh = -kappa * epcsaft::ELEMENTARY_CHARGE_C
        * epcsaft::ELEMENTARY_CHARGE_C * weighted_chi
        / (12.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M
           * expected_ion_suppression * epcsaft::BOLTZMANN_J_PER_K * temperature);
    record("dh-finite-multivalent", "debye_huckel", "1", electrolyte.screening,
           expected_dh, inputs(p, density, x, recording_support::fields({
               "case:multivalent-unequal-dh-diameter",
               recording_support::field("chi_cation_argument", chi_cation_argument),
               recording_support::field("chi_anion_argument", chi_anion_argument),
               recording_support::field("chi_reference_panels", 4000),
               recording_support::field("chi_simpson_bound", chi_quadrature_bound),
           })), 1.0e-10);
    const double dilute_coefficient = -std::sqrt(
        epcsaft::AVOGADRO_PER_MOL * epcsaft::ELEMENTARY_CHARGE_C
        * epcsaft::ELEMENTARY_CHARGE_C * charge_square_sum
        / (epcsaft::VACUUM_PERMITTIVITY_F_PER_M * expected_ion_suppression
           * epcsaft::BOLTZMANN_J_PER_K * temperature)
    ) * epcsaft::ELEMENTARY_CHARGE_C * epcsaft::ELEMENTARY_CHARGE_C * charge_square_sum
        / (12.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M
           * expected_ion_suppression * epcsaft::BOLTZMANN_J_PER_K * temperature);
    for (const double rho : {1.0e-6, 1.0e-8}) {
        const auto dilute = epcsaft::electrolyte_contribution(
            temperature, rho, x, x, p
        );
        const double coefficient = dilute.screening / std::sqrt(rho);
        record((std::string("dh-dilute-") + recording_support::number(rho)).c_str(),
               "debye_huckel_coefficient_ratio", "1", coefficient / dilute_coefficient, 1.0,
               inputs(p, rho, x, recording_support::fields({
                   "limit:Debye-square-root",
                   recording_support::field("actual_coefficient", coefficient),
                   recording_support::field("analytic_coefficient", dilute_coefficient),
                   recording_support::field("relative_bound", 5.0e-5),
               })), 5.0e-5);
    }

    const auto born_reference = [&](const Parameters& parameters, double bulk) {
        double solvation = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) solvation += x[i] * parameters.solvation_factor[i];
        double radial = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            const int charge = parameters.charge_number[i];
            if (charge == 0) continue;
            const double inner = parameters.born_diameter_angstrom[i] * epcsaft::ANGSTROM_M;
            const double outer = inner * (1.0 + parameters.born_c_shell
                * (solvation - 1.0) / std::abs(charge));
            const double outer_region = simpson(
                [bulk](double) { return 1.0 - 1.0 / bulk; },
                0.0, 1.0 / outer, 2
            );
            const double inner_region = simpson(
                [&parameters](double) {
                    return parameters.born_c_dielectric
                        * (1.0 - 1.0 / parameters.ionic_region_relative_permittivity)
                        ;
                }, 1.0 / outer, 1.0 / inner, 2
            );
            radial += x[i] * charge * charge * (outer_region + inner_region);
        }
        return -epcsaft::ELEMENTARY_CHARGE_C * epcsaft::ELEMENTARY_CHARGE_C * radial
            / (4.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M
               * epcsaft::BOLTZMANN_J_PER_K * temperature);
    };
    const double born = epcsaft::born_contribution(
        temperature, x, p, electrolyte.bulk_relative_permittivity
    );
    record("born-two-region", "born", "1", born,
           born_reference(p, expected_ion_suppression),
           inputs(p, density, x, recording_support::fields({
               "case:two-region-c-shell=c-dielectric=1",
               recording_support::field("born_bulk_relative_permittivity", electrolyte.bulk_relative_permittivity),
               recording_support::field("reference_bulk_relative_permittivity", expected_ion_suppression),
               recording_support::field("radial_panels_per_region", 2),
           })), 1.0e-12);
    Parameters direct = p;
    direct.solvation_factor = {1.0, 1.0, 1.0, 1.0};
    const double direct_born = epcsaft::born_contribution(
        temperature, x, direct, electrolyte.bulk_relative_permittivity
    );
    record("born-direct-limit", "born", "1", direct_born,
           born_reference(direct, expected_ion_suppression),
           inputs(direct, density, x, recording_support::fields({
               "limit:outer=inner",
               recording_support::field("born_bulk_relative_permittivity", electrolyte.bulk_relative_permittivity),
           })), 1.0e-12);
    Parameters equal_permittivity = p;
    equal_permittivity.ionic_region_relative_permittivity = electrolyte.bulk_relative_permittivity;
    record("born-equal-permittivity", "born", "1",
           epcsaft::born_contribution(temperature, x, equal_permittivity,
                                      electrolyte.bulk_relative_permittivity),
           born_reference(equal_permittivity, expected_ion_suppression),
           inputs(equal_permittivity, density, x, recording_support::fields({
               "limit:epsilon-inner=epsilon-bulk",
               recording_support::field("born_bulk_relative_permittivity", electrolyte.bulk_relative_permittivity),
           })), 1.0e-12);
    Parameters zero_contrast = p;
    zero_contrast.ionic_region_relative_permittivity = 1.0;
    record("born-zero-contrast", "born", "1",
           epcsaft::born_contribution(temperature, x, zero_contrast, 1.0), 0.0,
           inputs(zero_contrast, density, x, recording_support::fields({
               "limit:epsilon-inner=epsilon-bulk=1",
               recording_support::field("born_bulk_relative_permittivity", 1.0),
           })), 1.0e-14);
    Parameters charge_scaled = direct;
    charge_scaled.charge_number = {0, 0, 2, -4};
    record("born-charge-square", "born", "1",
           epcsaft::born_contribution(temperature, x, charge_scaled,
                                      electrolyte.bulk_relative_permittivity),
           4.0 * direct_born, inputs(charge_scaled, density, x, recording_support::fields({
               "limit:charges-doubled",
               recording_support::field("born_bulk_relative_permittivity", electrolyte.bulk_relative_permittivity),
           })), 1.0e-12);
}

void same_sign_exclusion_is_a_limit() {
    auto excluded = neutral_binary();
    excluded.charge_number[1] = 1;
    excluded.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
    auto all_pairs = excluded;
    all_pairs.ionic_dispersion = epcsaft::IonicDispersion::AllIonPairs;

    const auto one_excluded = evaluate(excluded, 300.0, 1.2e4, {0.0, 1.0});
    const auto one_all = evaluate(all_pairs, 300.0, 1.2e4, {0.0, 1.0});
    const auto inputs = recording_support::fields({
        model_inputs(excluded, 300.0, 1.2e4, {0.0, 1.0}),
        recording_support::field("policy_excluded", "ExcludeSameSign"),
        recording_support::field("policy_all", "AllIonPairs"),
    });
    recording_support::observe("excluded-pure-ion", "dispersion", "1", one_excluded.dispersion, 0.0,
        "ionic-pair-selection", inputs);
    recording_support::observe("all-pairs-pure-ion", "dispersion", "1", one_all.dispersion, 0.0,
        "ionic-pair-selection", inputs);
    const bool excluded_zero = one_excluded.dispersion == 0.0;
    recording_support::compare("excluded-zero", "dispersion_zero", "bool", excluded_zero ? 1.0 : 0.0, 1.0,
        recording_support::Rule::Absolute, 0.0, 0.0, "ionic-pair-selection",
        recording_support::fields({inputs, recording_support::field("dispersion", one_excluded.dispersion)}));
    const bool all_negative = one_all.dispersion < 0.0;
    recording_support::compare("all-pairs-negative", "dispersion_negative", "bool", all_negative ? 1.0 : 0.0, 1.0,
        recording_support::Rule::Absolute, 0.0, 0.0, "ionic-pair-selection",
        recording_support::fields({inputs, recording_support::field("dispersion", one_all.dispersion)}));
    const auto half_excluded = evaluate(excluded, 300.0, 1.2e4, {0.5, 0.5});
    const auto half_all = evaluate(all_pairs, 300.0, 1.2e4, {0.5, 0.5});
    const bool half_different = half_excluded.dispersion != half_all.dispersion;
    recording_support::compare("half-mixture-different", "dispersion_policy_difference", "bool",
        half_different ? 1.0 : 0.0, 1.0, recording_support::Rule::Absolute, 0.0, 0.0, "ionic-pair-selection",
        recording_support::fields({model_inputs(excluded, 300.0, 1.2e4, {0.5, 0.5}),
                                   recording_support::field("excluded_dispersion", half_excluded.dispersion),
                                   recording_support::field("all_pairs_dispersion", half_all.dispersion)}));
}

// The packing fraction is linear in density and cubic in diameter; close packing
// is the ceiling the hard-sphere logarithm cannot cross.
void packing_fraction_is_bounded() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    const auto sigma = epcsaft::segment_sigma(300.0, p);
    const auto diameter = epcsaft::segment_diameter(300.0, p, sigma);
    const double low =
        epcsaft::packing_fraction_from_geometry(1.0e3, x, p, diameter);
    const double high =
        epcsaft::packing_fraction_from_geometry(2.0e3, x, p, diameter);
    recording_support::compare("density-scaling", "packing_fraction", "1", high, 2.0 * low,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "density-scaling",
        recording_support::fields({model_inputs(p, 300.0, 1.0e3, x),
                                   recording_support::field("rho_high", 2.0e3),
                                   recording_support::field("eta_low", low),
                                   recording_support::field("eta_high", high)}));
    const bool low_positive = low > 0.0;
    recording_support::compare("low-positive", "packing_fraction_positive", "bool", low_positive ? 1.0 : 0.0, 1.0,
        recording_support::Rule::Absolute, 0.0, 0.0, "geometry",
        recording_support::fields({model_inputs(p, 300.0, 1.0e3, x),
                                         recording_support::field("eta", low)}));
    const bool policy_lower = epcsaft::policy::CLOSE_PACKING_FRACTION.value > 0.74;
    recording_support::compare("policy-lower", "close_packing_fraction", "1", policy_lower ? 1.0 : 0.0, 1.0,
        recording_support::Rule::Absolute, 0.0, 0.0, "policy-bound",
        recording_support::field("policy", epcsaft::policy::CLOSE_PACKING_FRACTION.value));
    const bool policy_upper = epcsaft::policy::CLOSE_PACKING_FRACTION.value < 0.75;
    recording_support::compare("policy-upper", "close_packing_fraction", "1", policy_upper ? 1.0 : 0.0, 1.0,
        recording_support::Rule::Absolute, 0.0, 0.0, "policy-bound",
        recording_support::field("policy", epcsaft::policy::CLOSE_PACKING_FRACTION.value));
}

// The kernels are templated so the same source both evaluates a double and
// records a tape; the AD instantiation is checked here against a centred
// difference of the double path.
void kernels_record_a_tape() {
    using AD = CppAD::AD<double>;
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    const double temperature = 300.0;
    const double molar_density = 1.2e4;

    std::vector<AD> independent = {AD(temperature), AD(molar_density)};
    CppAD::Independent(independent);
    const std::vector<AD> composition = {AD(x[0]), AD(x[1])};
    const auto sigma = epcsaft::segment_sigma(independent[0], p);
    const auto diameter = epcsaft::segment_diameter(independent[0], p, sigma);
    const std::vector<AD> no_sites;
    const auto residual = epcsaft::residual_helmholtz(
        independent[0], independent[1], composition, p, sigma, diameter,
        no_sites, {});
    std::vector<AD> dependent = {residual.hard_chain + residual.dispersion};
    CppAD::ADFun<double> tape(independent, dependent);

    const std::vector<double> point = {temperature, molar_density};
    const auto gradient = tape.Jacobian(point);

    const double step = molar_density * 1.0e-6;
    const auto up = evaluate(p, temperature, molar_density + step, x);
    const auto down = evaluate(p, temperature, molar_density - step, x);
    const double difference =
        ((up.hard_chain + up.dispersion) - (down.hard_chain + down.dispersion))
        / (2.0 * step);
    recording_support::compare("density-gradient", "d_a_res_d_rho", "1", gradient[1], difference,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-7, "derivative",
        recording_support::fields({model_inputs(p, temperature, molar_density, x),
                                   recording_support::field("step_rho", step),
                                   recording_support::field("rho_plus", molar_density + step),
                                   recording_support::field("rho_minus", molar_density - step)}));

    // The tape reproduces the double path exactly at the recording point.
    const auto value = tape.Forward(0, point);
    const auto direct = evaluate(p, temperature, molar_density, x);
    recording_support::compare("tape-value", "a_res", "1", value[0], direct.hard_chain + direct.dispersion,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "AD-evaluation",
        model_inputs(p, temperature, molar_density, x));
}

struct IdealCase {
    const char* name;
    epcsaft::IdealCorrelation correlation;
    double cp, h, s;
    const char* inputs;
};
// Six ideal correlations, each with independent cp, h and s at 550 K and 2e5 Pa.
const std::vector<IdealCase> ideal_cases = {
    {"constant", {epcsaft::IdealPolynomial{{25.0}, 400.0, 1000.0, 2.0}, 1.0e5},
     25.0, 4750.0, 4.198196956319387, "cp:25;Tref:400;href:1000;sref:2;Psrc:1e5"},
    {"polynomial", {epcsaft::IdealPolynomial{{20.0, 0.03, -1.0e-4}, 400.0, 1000.0, 2.0}, 1.0e5},
     22.25, 4225.0, 3.0642238294077453, "cp:20,.03,-1e-4;Tref:400;href:1000;sref:2;Psrc:1e5"},
    {"shomate", {epcsaft::IdealShomate{{20.0, 3.0, -1.0, 0.2, 4.0, 2.0, 3.0, 4.0}, 1000.0}, 1.0e5},
     34.60391549586777, 3130.139706439393, -19.821614918023606,
     "A..H:20,3,-1,.2,4,2,3,4;formation_h:1000;Psrc:1e5"},
    {"nasa7", {epcsaft::IdealNASA7{{3.5, 1.0e-3, -2.0e-6, 4.0e-9, -1.0e-12, 2.0, 3.0}}, 1.0e5},
     33.415773296966506, 17034.45396995443, 206.5148184322377,
     "a0..a6:3.5,1e-3,-2e-6,4e-9,-1e-12,2,3;Psrc:1e5"},
    {"nasa9", {epcsaft::IdealNASA9{{1.0e5, -20.0, 3.5, 1.0e-3, -2.0e-6, 4.0e-9, -1.0e-12, 2.0, 3.0}}, 1.0e5},
     35.862011885447956, 14473.461900968241, 205.44287118559976,
     "a0..a8:1e5,-20,3.5,1e-3,-2e-6,4e-9,-1e-12,2,3;Psrc:1e5"},
    {"dippr107", {epcsaft::IdealDIPPR107{{15.0, 22.0, 700.0, 8.0, 500.0}, 400.0, 1000.0, 2.0}, 1.0e5},
     31.34268666475855, 5422.806964662188, 5.5912660099491704,
     "A..E:15,22,700,8,500;Tref:400;href:1000;sref:2;Psrc:1e5"},
};

void ideal_thermal_forms_match_integrals() {
    const double pressure = 2.0e5;
    const double lower = 450.0;
    const double upper = 650.0;
    for (const auto& item : ideal_cases) {
        const char* name = item.name;
        const auto& correlation = item.correlation;
        const auto start = epcsaft::evaluate_ideal_correlation(correlation, lower, pressure);
        const auto finish = epcsaft::evaluate_ideal_correlation(correlation, upper, pressure);
        const double coarse_h = simpson(
            [&](double temperature) {
                return epcsaft::evaluate_ideal_correlation(
                    correlation, temperature, pressure).cp;
            }, lower, upper, 1000);
        const double expected_h = simpson(
            [&](double temperature) {
                return epcsaft::evaluate_ideal_correlation(
                    correlation, temperature, pressure).cp;
            }, lower, upper);
        const double fine_h = simpson(
            [&](double temperature) {
                return epcsaft::evaluate_ideal_correlation(
                    correlation, temperature, pressure).cp;
            }, lower, upper, 4000);
        const double coarse_s = simpson(
            [&](double temperature) {
                const auto value = epcsaft::evaluate_ideal_correlation(
                    correlation, temperature, pressure);
                return value.cp / temperature;
            }, lower, upper, 1000);
        const double expected_s = simpson(
            [&](double temperature) {
                const auto value = epcsaft::evaluate_ideal_correlation(
                    correlation, temperature, pressure);
                return value.cp / temperature;
            }, lower, upper);
        const double fine_s = simpson(
            [&](double temperature) {
                // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage): clang-tidy 18 misreads a lambda's reference capture of a structured binding as uninitialized.
                const auto value = epcsaft::evaluate_ideal_correlation(
                    correlation, temperature, pressure);
                return value.cp / temperature;
            }, lower, upper, 4000);
        const auto inputs = recording_support::fields({
            recording_support::field("form", name),
            recording_support::field("T_lower", lower),
            recording_support::field("T_upper", upper),
            recording_support::field("pressure", pressure),
        });
        const std::string h_id = std::string("dh-") + name;
        const std::string s_id = std::string("ds-") + name;
        recording_support::compare(h_id, "h_increment", "J/mol", finish.h - start.h, expected_h,
            recording_support::Rule::ReferenceOrOne, 1.0e-8, 1.0e-8, "independent_Simpson_integral", inputs);
        CHECK(close(finish.h - start.h, expected_h, 1.0e-8));
        recording_support::compare(s_id, "s_increment", "J/mol/K", finish.s - start.s, expected_s,
            recording_support::Rule::ReferenceOrOne, 1.0e-8, 1.0e-8, "independent_Simpson_integral", inputs);
        CHECK(close(finish.s - start.s, expected_s, 1.0e-8));
        recording_support::compare(std::string("refine-dh-") + name, "simpson_refinement", "J/mol", expected_h, fine_h,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-8, "independent_Simpson_refinement", inputs);
        recording_support::compare(std::string("refine-ds-") + name, "simpson_refinement", "J/mol/K", expected_s,
            fine_s, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-8, "independent_Simpson_refinement", inputs);
        CHECK(close(coarse_h, expected_h, 1.0e-8));
        CHECK(close(coarse_s, expected_s, 1.0e-8));
    }
}

void ideal_thermal_absolute_values_and_ad_derivatives() {
    const double temperature = 550.0;
    const double pressure = 2.0e5;
    for (const auto& item : ideal_cases) {
        const auto value = epcsaft::evaluate_ideal_correlation(item.correlation, temperature, pressure);
        const std::string inputs = std::string(item.inputs) + ";T:550;P:2e5";
        const auto check = [&](const char* prefix, const char* quantity, const char* unit,
                               double model, double reference, double rtol) {
            const std::string id = std::string(prefix) + item.name;
            recording_support::compare(id, quantity, unit, model, reference, recording_support::Rule::ReferenceOrOne,
                0.0, rtol, "independent_equation_evaluation", inputs);
        };
        check("cp-", "cp", "J/mol/K", value.cp, item.cp, 1.0e-12);
        check("h-", "h", "J/mol", value.h, item.h, 1.0e-12);
        check("s-", "s", "J/mol/K", value.s, item.s, 1.0e-12);

        using AD = CppAD::AD<double>;
        std::vector<AD> independent = {AD(temperature)};
        CppAD::Independent(independent);
        const auto taped = epcsaft::evaluate_ideal_correlation(
            item.correlation, independent[0], pressure);
        CppAD::ADFun<double> tape(independent, {taped.h, taped.s});
        const auto derivative = tape.Jacobian(std::vector<double>{temperature});
        check("dh-", "dh_dT", "J/mol/K", derivative[0], item.cp, 1.0e-12);
        check("ds-", "T_ds_dT", "J/mol/K", temperature * derivative[1], item.cp, 1.0e-12);
    }
}

void dippr_limits_and_pressure_shift_are_exact() {
    const double T = 300.0;
    const auto zero = epcsaft::IdealCorrelation{
        epcsaft::IdealDIPPR107{{2.0, 10.0, 0.0, 7.0, 0.0}, 250.0, 11.0, 4.0}, 1.0e5};
    const auto value = epcsaft::evaluate_ideal_correlation(zero, T, 1.0e5);
    const double dt = T - 250.0;
    const double expected_s = 12.0 * std::log(T / 250.0);
    recording_support::compare("zero-cp", "cp", "J/mol/K", value.cp, 12.0, recording_support::Rule::ReferenceOrOne,
        0.0, 1.0e-14, "analytic_limit", "T:300;C:0;E:0");
    recording_support::compare("zero-h", "h", "J/mol", value.h, 11.0 + 12.0 * dt,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "analytic_limit", "T:300;C:0;E:0");
    recording_support::compare("zero-s", "s", "J/mol/K", value.s, 4.0 + expected_s,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "analytic_limit", "T:300;C:0;E:0");

    const auto signed_c = epcsaft::IdealCorrelation{
        epcsaft::IdealDIPPR107{{2.0, 10.0, -700.0, 7.0, 500.0}, 250.0, 11.0, 4.0}, 1.0e5};
    const auto signed_e = epcsaft::IdealCorrelation{
        epcsaft::IdealDIPPR107{{2.0, 10.0, 700.0, 7.0, -500.0}, 250.0, 11.0, 4.0}, 1.0e5};
    const auto left = epcsaft::evaluate_ideal_correlation(signed_c, T, 1.0e5);
    const auto right = epcsaft::evaluate_ideal_correlation(signed_e, T, 1.0e5);
    CHECK(close(left.h, right.h, 1.0e-14));
    CHECK(close(left.s, right.s, 1.0e-14));

    const auto huge = epcsaft::IdealCorrelation{
        epcsaft::IdealDIPPR107{{2.0, 10.0, 1.0e300, 7.0, 1.0e300}, 250.0, 11.0, 4.0}, 1.0e5};
    const auto huge_value = epcsaft::evaluate_ideal_correlation(huge, T, 1.0e5);
    const double huge_s = 4.0 + 2.0 * std::log(T / 250.0);
    recording_support::compare("large-cp", "cp", "J/mol/K", huge_value.cp, 2.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "stable_limit", "T:300;C_abs:1e300;E_abs:1e300");
    recording_support::compare("large-h", "h", "J/mol", huge_value.h, 111.0, recording_support::Rule::ReferenceOrOne,
        0.0, 1.0e-12, "stable_limit", "T:300;C_abs:1e300;E_abs:1e300");
    recording_support::compare("large-s", "s", "J/mol/K", huge_value.s, huge_s,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "stable_limit", "T:300;C_abs:1e300;E_abs:1e300");

    const auto shifted = epcsaft::IdealCorrelation{
        epcsaft::IdealPolynomial{{12.0}, 250.0, 11.0,
                                 4.0 - epcsaft::GAS_CONSTANT_J_PER_MOL_K * std::log(2.0)}, 2.0e5};
    const auto original = epcsaft::evaluate_ideal_correlation(zero, T, 2.0e5);
    const auto shifted_value = epcsaft::evaluate_ideal_correlation(shifted, T, 2.0e5);
    recording_support::compare("signed-cp", "cp", "J/mol/K", left.cp, right.cp,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "even_argument", "T:300;C:+/-700;E:+/-500");
    recording_support::compare("source-pressure", "s", "J/mol/K", original.s, shifted_value.s,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "standard_state_shift",
        "T:300;P:200000;Psrc:100000->200000;sentropy-anchor-shifted");
}

void ideal_intervals_require_owned_endpoints() {
    const auto left = epcsaft::IdealInterval{
        200.0, 300.0, true, false,
        {epcsaft::IdealPolynomial{{10.0}, 250.0, 0.0, 0.0}, 1.0e5}};
    const auto right = epcsaft::IdealInterval{
        300.0, 400.0, true, true,
        {epcsaft::IdealPolynomial{{20.0}, 350.0, 0.0, 0.0}, 1.0e5}};
    const std::vector<epcsaft::IdealInterval> intervals = {left, right};
    const auto owned = [](const std::vector<epcsaft::IdealInterval>& set, double T) {
        return epcsaft::evaluate_ideal_correlation(epcsaft::ideal_interval(set, T)->correlation, T, 1.0e5);
    };
    const auto boundary = owned(intervals, 300.0);
    const auto left_interior = owned(intervals, 250.0);
    const auto interior = owned(intervals, 350.0);
    const std::vector<epcsaft::IdealInterval> reverse_intervals = {
        {200.0, 300.0, true, true,
         {epcsaft::IdealPolynomial{{10.0}, 250.0, 0.0, 0.0}, 1.0e5}},
        {300.0, 400.0, false, true,
         {epcsaft::IdealPolynomial{{20.0}, 350.0, 0.0, 0.0}, 1.0e5}},
    };
    const auto reverse_boundary = owned(reverse_intervals, 300.0);
    recording_support::compare("interval-cp", "cp", "J/mol/K", boundary.cp, 20.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "owned_branch",
        "T:300;left-upper-open;right-lower-closed");
    recording_support::compare("interval-left-cp", "cp", "J/mol/K", left_interior.cp, 10.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "owned_branch", "T:250;left-interval");
    recording_support::compare("interval-reverse-cp", "cp", "J/mol/K", reverse_boundary.cp, 10.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "owned_branch",
        "T:300;left-upper-closed;right-lower-open");
    recording_support::compare("interval-interior-cp", "cp", "J/mol/K", interior.cp, 20.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "owned_branch", "T:350;right-interval");

    const auto gap = std::vector<epcsaft::IdealInterval>{
        {200.0, 300.0, true, true, left.correlation},
        {400.0, 500.0, true, true,
         {epcsaft::IdealPolynomial{{20.0}, 450.0, 0.0, 0.0}, 1.0e5}}};
    const auto overlap = std::vector<epcsaft::IdealInterval>{
        {200.0, 300.0, true, true, left.correlation},
        {300.0, 400.0, true, true, right.correlation}};
    const auto open_gap = std::vector<epcsaft::IdealInterval>{
        {200.0, 300.0, true, false, left.correlation},
        {300.0, 400.0, false, true, right.correlation}};
    CHECK(recording_support::refusal([&] { epcsaft::validate_ideal_intervals(gap); }).find("overlapping or gapped ideal intervals") != std::string::npos);
    CHECK(recording_support::refusal([&] { epcsaft::validate_ideal_intervals(overlap); }).find("overlapping or gapped ideal intervals") != std::string::npos);
    CHECK(recording_support::refusal([&] { epcsaft::validate_ideal_intervals(open_gap); }).find("overlapping or gapped ideal intervals") != std::string::npos);
}

void ideal_polynomial_tape_reuses_smooth_interval() {
    using AD = CppAD::AD<double>;
    const epcsaft::IdealCorrelation correlation{
        epcsaft::IdealPolynomial{{20.0, 0.03, -1.0e-4}, 400.0, 1000.0, 2.0}, 1.0e5};
    std::vector<AD> independent = {AD(450.0)};
    CppAD::Independent(independent);
    const auto value = epcsaft::evaluate_ideal_correlation(correlation, independent[0], 1.0e5);
    CppAD::ADFun<double> tape(independent, {value.cp, value.h, value.s});
    const auto first = tape.Forward(0, std::vector<double>{450.0});
    const auto second = tape.Forward(0, std::vector<double>{475.0});
    const auto direct = epcsaft::evaluate_ideal_correlation(correlation, 475.0, 1.0e5);
    recording_support::compare("tape-cp", "cp", "J/mol/K", second[0], direct.cp,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "CppAD_replay", "T_record:450;T_replay:475");
    recording_support::compare("tape-h", "h", "J/mol", second[1], direct.h, recording_support::Rule::ReferenceOrOne,
        0.0, 1.0e-14, "CppAD_replay", "T_record:450;T_replay:475");
    recording_support::compare("tape-s", "s", "J/mol/K", second[2], direct.s, recording_support::Rule::ReferenceOrOne,
        0.0, 1.0e-14, "CppAD_replay", "T_record:450;T_replay:475");
    CHECK(close(first[0], epcsaft::evaluate_ideal_correlation(correlation, 450.0, 1.0e5).cp, 1.0e-14));
}

void ideal_mechanical_owner_preserves_chart() {
    const std::vector<double> amounts = {0.2, 0.9, 0.0};
    const double volume = 2.3;
    const double rho0 = 1.7;
    const double expected = amounts[0] * (std::log(amounts[0] / (volume * rho0)) - 1.0)
        + amounts[1] * (std::log(amounts[1] / (volume * rho0)) - 1.0);
    const double actual = epcsaft::ideal_mechanical_helmholtz(amounts, volume, rho0);
    recording_support::compare("off-normalized", "a_ideal_mechanical", "1", actual, expected,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "closed_form",
        recording_support::fields({recording_support::field("amounts", amounts),
                                   recording_support::field("volume", volume),
                                   recording_support::field("rho0", rho0)}));
    CHECK(recording_support::refusal([&] {
        epcsaft::ideal_mechanical_helmholtz(std::vector<double>{0.2, -0.1}, volume, rho0);
    }).find("invalid ideal amount") != std::string::npos);

    using AD = CppAD::AD<double>;
    const double tiny = std::numeric_limits<double>::min() / 1024.0;
    std::vector<AD> independent = {AD(volume)};
    CppAD::Independent(independent);
    const auto taped = epcsaft::ideal_mechanical_helmholtz(
        std::vector<AD>{AD(tiny)}, independent[0], AD(rho0), true);
    CppAD::ADFun<double> tape(independent, {taped});
    const double direct = epcsaft::ideal_mechanical_helmholtz(
        std::vector<double>{tiny}, volume, rho0, true);
    const double replay = tape.Forward(0, std::vector<double>{volume})[0];
    recording_support::compare("clamp-tiny", "a_ideal_mechanical", "1", replay, direct,
        recording_support::Rule::ReferenceScale, 0.0, 1.0e-14, "scalar_policy",
        "n:DBL_MIN/1024;V:2.3;rho0:1.7;clamp:true");

    std::vector<AD> zero_independent = {AD(volume)};
    CppAD::Independent(zero_independent);
    const auto zero_taped = epcsaft::ideal_mechanical_helmholtz(
        std::vector<AD>{AD(0.0), AD(1.0)}, zero_independent[0], AD(rho0));
    CppAD::ADFun<double> zero_tape(zero_independent, {zero_taped});
    const auto zero_value = zero_tape.Forward(0, std::vector<double>{volume})[0];
    const auto zero_gradient = zero_tape.Jacobian(std::vector<double>{volume});
    CHECK(std::isfinite(zero_value) && std::isfinite(zero_gradient[0]));
}

void ideal_polynomial_entropy_is_stable_across_switches() {
    constexpr int degree = 60;
    const double reference = 400.0;
    std::vector<double> coefficients(degree + 1, 0.0);
    coefficients[degree] = std::pow(200.0, -degree);
    const epcsaft::IdealCorrelation correlation{
        epcsaft::IdealPolynomial{coefficients, reference, 0.0, 0.0}, 1.0e5};
    auto independent_integral = [&](double u) {
        const long double ul = u;
        long double power = std::pow(ul, degree + 1);
        long double sum = 0.0L;
        for (int j = 0; j != 400; ++j) {
            sum += power / static_cast<long double>(degree + j + 1);
            power *= -ul;
        }
        return static_cast<double>(static_cast<long double>(coefficients[degree])
                                   * std::pow(static_cast<long double>(reference), degree)
                                   * sum);
    };
    const std::vector<double> offsets = {-0.51, -0.50, -0.49, 0.49, 0.50, 0.51};
    for (const double u : offsets) {
        const double temperature = reference * (1.0 + u);
        const auto value = epcsaft::evaluate_ideal_correlation(correlation, temperature, 1.0e5);
        const double expected = independent_integral(u);
        const std::string point = recording_support::number(u);
        recording_support::compare(std::string("poly-") + point, "entropy_increment", "J/mol/K", value.s, expected,
            recording_support::Rule::ReferenceOrOne, 1.0e-14, 1.0e-12, "independent_long_double_series",
            recording_support::fields({recording_support::field("degree", degree),
                                       recording_support::field("Tref", reference),
                                       recording_support::field("u", u)}));
        CHECK(close(value.s, expected, 1.0e-12));
    }
    for (const auto& [stress_degree, stress_u] : std::vector<std::pair<int, double>>{
             {150, 0.999}, {200, 0.999}, {300, 0.999}, {600, 0.51},
             {600, -0.9}, {600, -0.95}, {600, -0.98}, {600, -0.999},
             {600, -std::exp(-1.0 / 600) * (1.0 - 1.0e-6)},
             {600, -std::exp(-1.0 / 600) * (1.0 + 1.0e-6)}, {60, 1.0}, {60, 1.01}, {2, 10.0}}) {
        std::vector<double> stress_coefficients(stress_degree + 1, 0.0);
        stress_coefficients[stress_degree] = std::pow(std::abs(stress_u), -stress_degree);
        const epcsaft::IdealCorrelation stress{
            epcsaft::IdealPolynomial{stress_coefficients, 1.0, 0.0, 0.0}, 1.0e5};
        const auto value = epcsaft::evaluate_ideal_correlation(
            stress, 1.0 + stress_u, 1.0e5);
        long double power = std::pow(static_cast<long double>(stress_u), stress_degree + 1);
        long double sum = 0.0L;
        if (stress_u > 0.0) {
            const long double ratio = static_cast<long double>(stress_u) / (1.0L + stress_u);
            power /= 1.0L + stress_u;
            long double term = 1.0L / (stress_degree + 1.0L);
            sum = term;
            for (int j = 0; j != 400; ++j) {
                term *= ratio * (j + 1.0L) / (stress_degree + j + 2.0L);
                sum += term;
            }
        } else {
            for (int j = 0; j != 100000; ++j) {
                sum += power / static_cast<long double>(stress_degree + j + 1);
                power *= -stress_u;
            }
        }
        const double expected = static_cast<double>(
            static_cast<long double>(stress_coefficients[stress_degree])
            * (stress_u > 0.0 ? power * sum : sum));
        const std::string point = std::to_string(stress_degree) + "/"
            + recording_support::number(stress_u);
        recording_support::compare(std::string("poly-stress-") + point, "entropy_increment", "J/mol/K", value.s,
            expected, recording_support::Rule::ReferenceOrOne, 1.0e-12, 1.0e-12, "independent_long_double_series",
            recording_support::fields({recording_support::field("degree", stress_degree),
                                       recording_support::field("Tref", 1.0),
                                       recording_support::field("u", stress_u)}));
        CHECK(close(value.s, expected, 1.0e-12));
    }
    const epcsaft::IdealCorrelation near_anchor{
        epcsaft::IdealPolynomial{{0.0, 0.0, 1.0}, reference, 0.0, 0.0}, 1.0e5};
    for (const double u : {-0.01, 0.01}) {
        const auto value = epcsaft::evaluate_ideal_correlation(
            near_anchor, reference * (1.0 + u), 1.0e5);
        long double power = static_cast<long double>(u) * u * u;
        long double expected_sum = 0.0L;
        for (int j = 0; j != 400; ++j) {
            expected_sum += power / static_cast<long double>(j + 3);
            power *= -u;
        }
        const double expected = static_cast<double>(reference * reference * expected_sum);
        CHECK(close(value.s, expected, 1.0e-12));
    }

    using AD = CppAD::AD<double>;
    std::vector<AD> independent = {AD(reference * 1.49)};
    CppAD::Independent(independent);
    const auto taped = epcsaft::evaluate_ideal_correlation(correlation, independent[0], 1.0e5);
    CppAD::ADFun<double> tape(independent, {taped.s});
    const double replay = tape.Forward(0, std::vector<double>{reference * 1.51})[0];
    const double direct = epcsaft::evaluate_ideal_correlation(
        correlation, reference * 1.51, 1.0e5).s;
    recording_support::compare("poly-tape-switch", "entropy_increment", "J/mol/K", replay, direct,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "CppAD_replay", "degree:60;Trecord:596;Treplay:604");

    std::vector<AD> boundary_independent = {AD(reference * 1.99)};
    CppAD::Independent(boundary_independent);
    const auto boundary_taped = epcsaft::evaluate_ideal_correlation(
        correlation, boundary_independent[0], 1.0e5);
    CppAD::ADFun<double> boundary_tape(boundary_independent, {boundary_taped.s});
    const double boundary_replay = boundary_tape.Forward(
        0, std::vector<double>{reference * 2.01})[0];
    const double boundary_direct = epcsaft::evaluate_ideal_correlation(
        correlation, reference * 2.01, 1.0e5).s;
    recording_support::compare("poly-tape-large-switch", "entropy_increment", "J/mol/K", boundary_replay,
        boundary_direct, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "CppAD_replay",
        "degree:60;Trecord:796;Treplay:804");

    std::vector<double> overflow_coefficients(151, 0.0);
    overflow_coefficients[150] = 1.0;
    const epcsaft::IdealCorrelation overflow_correlation{
        epcsaft::IdealPolynomial{overflow_coefficients, 400.0, 0.0, 0.0}, 1.0e5};
    const auto overflow_anchor = epcsaft::evaluate_ideal_correlation(
        overflow_correlation, 400.0, 1.0e5);
    const auto overflow_one = epcsaft::evaluate_ideal_correlation(
        overflow_correlation, 401.0, 1.0e5);
    long double overflow_power = std::pow(1.0L / 400.0L, 151);
    long double overflow_sum = 0.0L;
    for (int j = 0; j != 400; ++j) {
        overflow_sum += overflow_power / (151.0L + j);
        overflow_power *= -1.0L / 400.0L;
    }
    const double overflow_entropy = static_cast<double>(
        std::pow(400.0L, 150) * overflow_sum);
    for (const auto& [point, value, expected_cp, expected_h, expected_s] : std::vector<
             std::tuple<const char*, epcsaft::IdealThermalValue<double>, double, double, double>>{
             {"anchor", overflow_anchor, 0.0, 0.0, 0.0},
             {"T401", overflow_one, 1.0, 1.0 / 151.0, overflow_entropy}}) {
        const std::string inputs = std::string("a150:1;Tref:400;T:")
            + (point == std::string("anchor") ? "400" : "401");
        recording_support::compare(std::string("poly-overflow/") + point + "/cp", "cp", "J/mol/K", value.cp,
            expected_cp, recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "independent_equation_evaluation",
            inputs);
        recording_support::compare(std::string("poly-overflow/") + point + "/h", "h", "J/mol", value.h, expected_h,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "independent_equation_evaluation", inputs);
        recording_support::compare(std::string("poly-overflow/") + point + "/s", "s", "J/mol/K", value.s, expected_s,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "independent_equation_evaluation", inputs);
    }

    // A finite coefficient-weighted power must not form the overflowing raw power.
    overflow_coefficients[150] = 1.0e-300;
    const epcsaft::IdealCorrelation scaled_correlation{
        epcsaft::IdealPolynomial{overflow_coefficients, 400.0, 0.0, 0.0}, 1.0e5};
    const auto scaled = epcsaft::evaluate_ideal_correlation(scaled_correlation, 800.0, 1.0e5);
    const long double scaled_cp = 1.0e-300L * std::pow(400.0L, 150);
    long double scaled_term = 1.0L / 151, scaled_sum = scaled_term;
    for (int j = 0; j < 400; ++j) {
        scaled_term *= 0.5L * (j + 1) / (152 + j);
        scaled_sum += scaled_term;
    }
    for (const auto& [name, actual, expected] : std::vector<std::tuple<const char*, double, double>>{
             {"cp", scaled.cp, static_cast<double>(scaled_cp)},
             {"h", scaled.h, static_cast<double>(scaled_cp * 400 / 151)},
             {"s", scaled.s, static_cast<double>(scaled_cp * scaled_sum / 2)}}) {
        const std::string point = std::string("poly-scaled/") + name;
        recording_support::compare(point, name, std::string(name) == "h" ? "J/mol" : "J/mol/K", actual, expected,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "independent_long_double_series",
            "a150:1e-300;Tref:400;T:800");
    }

    // A nonconstant tape traverses both signs, both recurrences and the anchor.
    const epcsaft::IdealCorrelation quadratic{
        epcsaft::IdealPolynomial{{0.0, 0.0, 1.0}, 1.0, 0.0, 0.0}, 1.0e5};
    std::vector<AD> quadratic_record = {AD(1.5)};
    CppAD::Independent(quadratic_record);
    const auto quadratic_value = epcsaft::evaluate_ideal_correlation(
        quadratic, quadratic_record[0], 1.0e5);
    CppAD::ADFun<double> quadratic_tape(quadratic_record, {quadratic_value.h, quadratic_value.s});
    for (const double temperature : {0.1, 0.8, 0.4, 1.0, 1.2, 2.0, 10.0, 1.0e-20}) {
        const long double dt = static_cast<long double>(temperature) - 1;
        const double expected_s = static_cast<double>(dt*dt/2 - dt + std::log(static_cast<long double>(temperature)));
        const double cp = static_cast<double>(dt*dt), cp_prime = static_cast<double>(2*dt);
        const auto value = quadratic_tape.Forward(0, std::vector<double>{temperature});
        const auto first = quadratic_tape.Forward(1, std::vector<double>{1.0});
        const auto second = quadratic_tape.Forward(2, std::vector<double>{0.0});
        const auto third = quadratic_tape.Forward(3, std::vector<double>{0.0});
        for (const auto& [name, actual, expected] : std::vector<std::tuple<const char*, double, double>>{
                 {"s", value[1], expected_s}, {"dh", first[0], cp},
                 {"Tds", temperature*first[1], cp},
                 {"T2d2s", temperature*temperature*2*second[1], temperature*cp_prime-cp},
                 {"T3d3s", temperature*temperature*temperature*6*third[1],
                  2*temperature*temperature-2*temperature*cp_prime+2*cp}}) {
            const std::string point = "poly-replay/" + recording_support::number(temperature) + "/" + name;
            recording_support::compare(point, name, "J/mol/K", actual, expected,
                recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-12, "analytic_integral_and_derivatives",
                "a2:1;Tref:1;Trecord:1.5;T:" + recording_support::number(temperature));
        }
    }

    const epcsaft::IdealCorrelation constant_cp{
        epcsaft::IdealPolynomial{{30.0}, 400.0, 0.0, 0.0}, 1.0e5};
    std::vector<AD> positive_record = {AD(500.0)};
    CppAD::Independent(positive_record);
    const auto positive_taped = epcsaft::evaluate_ideal_correlation(
        constant_cp, positive_record[0], 1.0e5);
    CppAD::ADFun<double> positive_tape(positive_record, {positive_taped.s});
    for (const double temperature : {300.0, 200.0, 40.0}) {
        const double replay = positive_tape.Forward(0, std::vector<double>{temperature})[0];
        const double expected = 30.0 * std::log(temperature / 400.0);
        const std::string point = "negative-replay/" + std::to_string(temperature);
        recording_support::compare(point, "entropy", "J/mol/K", replay, expected,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "CppAD_replay",
            "cp:30;Tref:400;Trecord:500;Treplay:" + std::to_string(temperature));
    }
    std::vector<AD> negative_record = {AD(300.0)};
    CppAD::Independent(negative_record);
    const auto negative_taped = epcsaft::evaluate_ideal_correlation(
        constant_cp, negative_record[0], 1.0e5);
    CppAD::ADFun<double> negative_tape(negative_record, {negative_taped.s});
    for (const double temperature : {200.0, 40.0}) {
        const double replay = negative_tape.Forward(0, std::vector<double>{temperature})[0];
        const double expected = 30.0 * std::log(temperature / 400.0);
        const std::string point = "negative-negative/" + std::to_string(temperature);
        recording_support::compare(point, "entropy", "J/mol/K", replay, expected,
            recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "CppAD_replay",
            "cp:30;Tref:400;Trecord:300;Treplay:" + std::to_string(temperature));
    }
    std::vector<AD> derivative_record = {AD(300.0)};
    CppAD::Independent(derivative_record);
    const auto derivative_taped = epcsaft::evaluate_ideal_correlation(
        constant_cp, derivative_record[0], 1.0e5);
    CppAD::ADFun<double> derivative_tape(
        derivative_record, {derivative_taped.h, derivative_taped.s});
    const auto derivative = derivative_tape.Jacobian(std::vector<double>{200.0});
    recording_support::compare("negative-derivative/dh", "dh_dT", "J/mol/K", derivative[0], 30.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "CppAD_derivative", "cp:30;Tref:400;T:200");
    recording_support::compare("negative-derivative/ds", "T_ds_dT", "J/mol/K", 200.0 * derivative[1], 30.0,
        recording_support::Rule::ReferenceOrOne, 0.0, 1.0e-14, "CppAD_derivative", "cp:30;Tref:400;T:200");
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.contributions.p1", "P1"}, hard_sphere_matches_carnahan_starling},
        {{"native.contributions.c1_mixture", "C1"}, mixture_contributions_match_independent_equations},
        {{"native.contributions.p2", "P2"}, residual_vanishes_at_zero_density},
        {{"native.contributions.p3", "P3"}, relabelling_components_changes_nothing},
        {{"native.contributions.p4", "P4"}, association_solve_satisfies_mass_action},
        {{"native.contributions.association-unit-ceiling", "association-unit-ceiling"}, association_unit_ceiling},
        {{"native.contributions.c2", "C2"}, association_equations_match_closed_forms_and_limits},
        {{"native.contributions.kraska-model-iv", "Kraska"}, kraska_model_iv},
        {{"native.contributions.huang-radosz-table-vii", "HuangRadosz"}, huang_radosz_table_vii},
        {{"native.contributions.c3", "C3"}, electrolyte_born_and_dielectric_match_independent_equations},
        {{"native.contributions.p5", "P5"}, inactive_terms_contribute_nothing},
        {{"native.contributions.p6", "P6"}, same_sign_exclusion_is_a_limit},
        {{"native.contributions.p7", "P7"}, packing_fraction_is_bounded},
        {{"native.contributions.p8", "P8"}, kernels_record_a_tape},
        {{"native.contributions.c1_ideal", "C1"}, ideal_mechanical_owner_preserves_chart},
        {{"native.contributions.c1_thermal", "C1"}, ideal_thermal_forms_match_integrals},
        {{"native.contributions.c1_absolute", "C1"}, ideal_thermal_absolute_values_and_ad_derivatives},
        {{"native.contributions.c1_dippr", "C1"}, dippr_limits_and_pressure_shift_are_exact},
        {{"native.contributions.c1_poly_stability", "C1"}, ideal_polynomial_entropy_is_stable_across_switches},
        {{"native.contributions.c1_intervals", "C1"}, ideal_intervals_require_owned_endpoints},
        {{"native.contributions.c1_tape", "C1"}, ideal_polynomial_tape_reuses_smooth_interval},
    });
}
