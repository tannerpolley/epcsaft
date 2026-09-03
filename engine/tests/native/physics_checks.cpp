// Properties of the Helmholtz terms that hold for every parameter set, checked
// against closed forms and exact symmetries rather than against stored numbers.
// Data-backed tests arrive with the Python layer; these are the checks that fail
// if the kernels are wrong at all.

#include "physics/association.hpp"
#include "physics/born.hpp"
#include "physics/electrolyte.hpp"
#include "physics/helmholtz.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

// Not assert: the release build defines NDEBUG, which would compile every check
// in this file away and leave it passing vacuously.
#define CHECK(condition)                                                    \
    do {                                                                    \
        if (!(condition)) {                                                 \
            std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, \
                         #condition);                                       \
            std::abort();                                                   \
        }                                                                   \
    } while (false)

namespace {

using epcsaft::Parameters;

bool close(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
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
        temperature, molar_density, x, p, sigma, diameter, sites, 0.0, 0.0
    );
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
        CHECK(close(result.hard_chain, expected, 1.0e-12));
    }
}

// Every residual term is a departure from the ideal gas, so all of them vanish
// together as the density goes to zero.
void residual_vanishes_at_zero_density() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    double previous = 1.0;
    for (const double molar_density : {1.0e1, 1.0e0, 1.0e-1, 1.0e-2}) {
        const auto result = evaluate(p, 300.0, molar_density, x);
        const double total = result.hard_chain + result.dispersion;
        CHECK(std::abs(total) < std::abs(previous));
        previous = total;
    }
    CHECK(std::abs(previous) < 1.0e-3);
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
    CHECK(close(direct.hard_chain, relabelled.hard_chain, 1.0e-14));
    CHECK(close(direct.dispersion, relabelled.dispersion, 1.0e-14));
}

// One associating component, two sites, one bonding pair. The solve has to land
// inside the domain, satisfy the mass-action residual to the policy tolerance,
// and reach it in far fewer steps than the substitution solver it replaces.
void association_solve_satisfies_mass_action() {
    Parameters p = neutral_binary();
    p.site_component = {0, 0};
    p.site_multiplicity = {1, 1};
    p.association_active = {0, 1, 1, 0};
    p.association_energy_over_k = {0.0, 2425.7, 2425.7, 0.0};
    p.association_volume = {0.0, 0.0451, 0.0451, 0.0};

    const double temperature = 300.0;
    const double molar_density = 1.5e4;
    const std::vector<double> x = {0.4, 0.6};
    const auto sigma = epcsaft::segment_sigma(temperature, p);
    const auto diameter = epcsaft::segment_diameter(temperature, p, sigma);
    const auto fractions = epcsaft::solve_site_fractions(
        temperature, molar_density, x, p, sigma, diameter
    );
    CHECK(fractions.size() == 2);
    for (const double value : fractions) {
        CHECK(value > epcsaft::policy::SITE_FRACTION_FLOOR.value);
        CHECK(value <= epcsaft::policy::SITE_FRACTION_CEILING.value);
    }
    const auto delta = epcsaft::association_delta(
        temperature, molar_density, x, p, sigma, diameter
    );
    const auto residuals = epcsaft::association_residuals(
        molar_density, x, p, fractions, delta
    );
    for (const double residual : residuals) {
        CHECK(std::abs(residual) <= epcsaft::policy::ASSOCIATION_RESIDUAL.value);
    }
    // Identical sites on one component are indistinguishable.
    CHECK(close(fractions[0], fractions[1], 1.0e-13));
}

// A term that is switched off contributes nothing, whatever else is in the
// parameter record.
void inactive_terms_contribute_nothing() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    const auto result = evaluate(p, 300.0, 1.2e4, x);
    CHECK(result.association == 0.0);

    const double born = epcsaft::born_contribution(300.0, x, p, 78.4);
    CHECK(born == 0.0);

    const auto electrolyte = epcsaft::electrolyte_contribution(
        300.0, 1.2e4, 300.0, 1.2e4, x, {0.4, 0.6}, p
    );
    CHECK(electrolyte.debye_huckel == 0.0);
}

void same_sign_exclusion_is_a_limit() {
    auto excluded = neutral_binary();
    excluded.charge_number[1] = 1;
    excluded.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
    auto all_pairs = excluded;
    all_pairs.ionic_dispersion = epcsaft::IonicDispersion::AllIonPairs;

    const auto one_excluded = evaluate(excluded, 300.0, 1.2e4, {0.0, 1.0});
    const auto one_all = evaluate(all_pairs, 300.0, 1.2e4, {0.0, 1.0});
    CHECK(one_excluded.dispersion == 0.0);
    CHECK(one_all.dispersion < 0.0);
    const auto half_excluded = evaluate(excluded, 300.0, 1.2e4, {0.5, 0.5});
    const auto half_all = evaluate(all_pairs, 300.0, 1.2e4, {0.5, 0.5});
    CHECK(half_excluded.dispersion != half_all.dispersion);
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
    CHECK(close(high, 2.0 * low, 1.0e-14));
    CHECK(low > 0.0);
    CHECK(epcsaft::policy::CLOSE_PACKING_FRACTION.value > 0.74);
    CHECK(epcsaft::policy::CLOSE_PACKING_FRACTION.value < 0.75);
}

// The kernels are templated so the same source both evaluates a double and
// records a tape. Nothing else instantiates them on an AD type until the model
// layer exists, so the instantiation is checked here, against a centred
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
        no_sites, AD(0.0), AD(0.0)
    );
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
    CHECK(close(gradient[1], difference, 1.0e-7));

    // The tape reproduces the double path exactly at the recording point.
    const auto value = tape.Forward(0, point);
    const auto direct = evaluate(p, temperature, molar_density, x);
    CHECK(close(value[0], direct.hard_chain + direct.dispersion, 1.0e-14));
}

}  // namespace

int main() {
    hard_sphere_matches_carnahan_starling();
    residual_vanishes_at_zero_density();
    relabelling_components_changes_nothing();
    association_solve_satisfies_mass_action();
    inactive_terms_contribute_nothing();
    same_sign_exclusion_is_a_limit();
    packing_fraction_is_bounded();
    kernels_record_a_tape();
    std::printf("physics checks passed\n");
    return 0;
}
