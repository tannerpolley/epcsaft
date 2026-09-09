// Properties of the Helmholtz terms that hold for every parameter set, checked
// against closed forms and exact symmetries rather than against stored numbers.
// Data-backed tests arrive with the Python layer; these are the checks that fail
// if the kernels are wrong at all.

#include "contributions/association.hpp"
#include "contributions/born.hpp"
#include "contributions/electrolyte.hpp"
#include "contributions/helmholtz.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <string>
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

std::vector<double> as_double(const std::vector<int>& values) {
    return std::vector<double>(values.begin(), values.end());
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
        recording_support::field("charge", as_double(parameters.charge_number)),
        recording_support::field("kij", parameters.k_ij),
        recording_support::field("site_component", as_double(parameters.site_component)),
        recording_support::field("site_multiplicity", as_double(parameters.site_multiplicity)),
        recording_support::field("association_active", as_double(parameters.association_active)),
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
    recording_support::expect_group(
        recording_support::active_test_id(), "p1-hard-chain", "eta-", 3,
        "Carnahan-Starling hard-sphere values", "rho");
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
        recording_support::compare(
            (std::string("eta-") + point).c_str(), "hard_chain", "1",
            result.hard_chain, expected, inputs, "Relationship", "closed_form",
            "hard_chain=Carnahan_Starling(eta)", "reference_or_one", 0.0, 1.0e-12,
            close(result.hard_chain, expected, 1.0e-12), "p1-hard-chain", point.c_str(),
            "a^{hc}", "T,\\rho,\\mathbf{x}",
            "a^{hc}=(4\\eta-3\\eta^2)/(1-\\eta)^2", "rho");
        CHECK(close(result.hard_chain, expected, 1.0e-12));
    }
}

// Every residual term is a departure from the ideal gas, so all of them vanish
// together as the density goes to zero.
void residual_vanishes_at_zero_density() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-monotonic", "density-monotonic/", 4,
        "residual magnitude decreases along density ladder", "rho");
    recording_support::expect_group(
        recording_support::active_test_id(), "p2-limit", "zero-density-limit", 1,
        "zero-density residual limit", "rho");
    double previous = 1.0;
    double previous_density = 0.0;
    std::size_t point_index = 0;
    for (const double molar_density : {1.0e1, 1.0e0, 1.0e-1, 1.0e-2}) {
        const auto result = evaluate(p, 300.0, molar_density, x);
        const double total = result.hard_chain + result.dispersion;
        const bool decreasing = std::abs(total) < std::abs(previous);
        recording_support::compare(
            (std::string("density-monotonic/") + std::to_string(point_index)).c_str(),
            "abs_residual", "1", decreasing ? 1.0 : 0.0, 1.0,
            recording_support::fields({
                model_inputs(p, 300.0, molar_density, x),
                recording_support::field("previous_density", previous_density),
                recording_support::field("previous_abs_residual", std::abs(previous)),
                recording_support::field("current_abs_residual", std::abs(total)),
            }), "Relationship", "density-limit-ladder",
            "abs(a_res(rho)) decreases along the density ladder", "absolute", 0.0, 0.0,
            decreasing, "p2-monotonic", std::to_string(point_index).c_str(), "|a^{res}|",
            "T,\\rho,\\mathbf{x},a_{previous},a_{current}",
            "|a^{res}(\\rho_i)|<|a^{res}(\\rho_{i-1})|", "rho");
        CHECK(decreasing);
        previous = total;
        previous_density = molar_density;
        ++point_index;
    }
    recording_support::compare(
        "zero-density-limit", "residual_helmholtz", "1", previous, 0.0,
        model_inputs(p, 300.0, 1.0e-2, x),
        "Relationship", "limit", "a_res->0 as rho->0", "strict_absolute", 1.0e-3, 0.0,
        std::abs(previous) < 1.0e-3, "p2-limit", "rho-1e-2", "a^{res}",
        "T,\\rho,\\mathbf{x}", "\\lim_{\\rho\\to0}a^{res}=0", "rho");
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
    recording_support::expect_group(
        recording_support::active_test_id(), "p3-hard-chain", "hard-chain", 1,
        "component relabelling hard-chain invariance", "permutation");
    recording_support::expect_group(
        recording_support::active_test_id(), "p3-dispersion", "dispersion", 1,
        "component relabelling dispersion invariance", "permutation");
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
    recording_support::compare(
        "hard-chain", "hard_chain", "1", direct.hard_chain, relabelled.hard_chain,
        inputs, "Relationship", "permutation", "component relabelling invariance",
        "reference_or_one", 0.0, 1.0e-14,
        close(direct.hard_chain, relabelled.hard_chain, 1.0e-14),
        "p3-hard-chain", "hard-chain", "a^{hc}", "T,\\rho,\\mathbf{x}",
        "a^{hc}(x)=a^{hc}(\\pi x)", "x");
    recording_support::compare(
        "dispersion", "dispersion", "1", direct.dispersion, relabelled.dispersion,
        inputs, "Relationship", "permutation", "component relabelling invariance",
        "reference_or_one", 0.0, 1.0e-14,
        close(direct.dispersion, relabelled.dispersion, 1.0e-14),
        "p3-dispersion", "dispersion", "a^{disp}", "T,\\rho,\\mathbf{x}",
        "a^{disp}(x)=a^{disp}(\\pi x)", "x");
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
    recording_support::expect_group(
        recording_support::active_test_id(), "p4-mass-action", "mass-action-", 2,
        "association mass-action residuals", "site");
    recording_support::expect_group(
        recording_support::active_test_id(), "p4-site-floor", "site-floor-", 2,
        "association site fractions stay above floor", "site");
    recording_support::expect_group(
        recording_support::active_test_id(), "p4-site-ceiling", "site-ceiling-", 2,
        "association site fractions stay below ceiling", "site");
    recording_support::expect_group(
        recording_support::active_test_id(), "p4-site-symmetry", "site-symmetry", 1,
        "identical sites have identical fractions", "site");
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
        recording_support::compare(
            (std::string("site-floor-") + std::to_string(site)).c_str(),
            "site_fraction_above_floor", "bool", above_floor ? 1.0 : 0.0, 1.0,
            inputs, "Relationship", "association-domain", "site fraction is above floor",
            "absolute", 0.0, 0.0, above_floor, "p4-site-floor", std::to_string(site).c_str(),
            "\\chi", "T,\\rho,\\mathbf{x},\\chi_{floor}", "\\chi>\\chi_{floor}", "site");
        CHECK(above_floor);
        const bool below_ceiling = value <= epcsaft::policy::SITE_FRACTION_CEILING.value;
        recording_support::compare(
            (std::string("site-ceiling-") + std::to_string(site)).c_str(),
            "site_fraction_below_ceiling", "bool", below_ceiling ? 1.0 : 0.0, 1.0,
            inputs, "Relationship", "association-domain", "site fraction is below ceiling",
            "absolute", 0.0, 0.0, below_ceiling, "p4-site-ceiling", std::to_string(site).c_str(),
            "\\chi", "T,\\rho,\\mathbf{x},\\chi_{ceiling}", "\\chi\\le\\chi_{ceiling}", "site");
        CHECK(below_ceiling);
    }
    const auto delta = epcsaft::association_delta(
        temperature, molar_density, x, p, sigma, diameter
    );
    const auto residuals = epcsaft::association_residuals(
        molar_density, x, p, fractions, delta
    );
    for (std::size_t site = 0; site < residuals.size(); ++site) {
        const double residual = residuals[site];
        recording_support::compare(
            (std::string("mass-action-") + std::to_string(site)).c_str(),
            "association_residual", "1", residual, 0.0,
            recording_support::fields({model_inputs(p, temperature, molar_density, x),
                                       recording_support::field("site_fractions", fractions),
                                       recording_support::field("site", site)}),
            "Relationship", "mass_action", "association residual is zero",
            "reference_or_one", epcsaft::policy::ASSOCIATION_RESIDUAL.value, 0.0,
            std::abs(residual) <= epcsaft::policy::ASSOCIATION_RESIDUAL.value,
            "p4-mass-action", std::to_string(site).c_str(), "F", "T,\\rho,\\mathbf{x},\\chi",
            "f(\\chi)=0", "site");
        CHECK(std::abs(residual) <= epcsaft::policy::ASSOCIATION_RESIDUAL.value);
    }
    // Identical sites on one component are indistinguishable.
    const bool identical_sites = close(fractions[0], fractions[1], 1.0e-13);
    recording_support::compare(
        "site-symmetry", "site_fraction_symmetry", "bool",
        identical_sites ? 1.0 : 0.0, 1.0,
        recording_support::fields({model_inputs(p, temperature, molar_density, x),
                                   recording_support::field("site_fractions", fractions)}),
        "Relationship", "site-symmetry", "identical sites have identical fractions",
        "absolute", 0.0, 0.0, identical_sites, "p4-site-symmetry", "site",
        "\\chi", "T,\\rho,\\mathbf{x},\\chi", "\\chi_0=\\chi_1", "site");
    CHECK(identical_sites);
}

// A term that is switched off contributes nothing, whatever else is in the
// parameter record.
void inactive_terms_contribute_nothing() {
    const auto p = neutral_binary();
    const std::vector<double> x = {0.4, 0.6};
    recording_support::expect_group(
        recording_support::active_test_id(), "p5-association", "association", 1,
        "inactive association term", "term");
    recording_support::expect_group(
        recording_support::active_test_id(), "p5-born", "born", 1,
        "inactive Born term", "term");
    recording_support::expect_group(
        recording_support::active_test_id(), "p5-debye", "debye-huckel", 1,
        "inactive electrolyte term", "term");
    const auto result = evaluate(p, 300.0, 1.2e4, x);
    const auto inputs = model_inputs(p, 300.0, 1.2e4, x);
    recording_support::compare(
        "association", "association", "1", result.association, 0.0, inputs,
        "Relationship", "inactive-term", "inactive association term vanishes",
        "reference_or_one", 0.0, 0.0, result.association == 0.0,
        "p5-association", "inactive", "a^{assoc}", "T,\\rho,\\mathbf{x}",
        "a^{assoc}=0", "term");
    CHECK(result.association == 0.0);

    const double born = epcsaft::born_contribution(300.0, x, p, 78.4);
    recording_support::compare(
        "born", "born", "1", born, 0.0, inputs, "Relationship", "inactive-term",
        "inactive Born term vanishes", "reference_or_one", 0.0, 0.0, born == 0.0,
        "p5-born", "inactive", "a^{Born}", "T,\\mathbf{x}", "a^{Born}=0", "term");
    CHECK(born == 0.0);

    const auto electrolyte = epcsaft::electrolyte_contribution(
        300.0, 1.2e4, 300.0, 1.2e4, x, {0.4, 0.6}, p
    );
    recording_support::compare(
        "debye-huckel", "debye_huckel", "1", electrolyte.debye_huckel, 0.0, inputs,
        "Relationship", "inactive-term", "inactive electrolyte term vanishes",
        "reference_or_one", 0.0, 0.0, electrolyte.debye_huckel == 0.0,
        "p5-debye", "inactive", "a^{DH}", "T,\\rho,\\mathbf{x}", "a^{DH}=0", "term");
    CHECK(electrolyte.debye_huckel == 0.0);
}

void same_sign_exclusion_is_a_limit() {
    auto excluded = neutral_binary();
    excluded.charge_number[1] = 1;
    excluded.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;
    auto all_pairs = excluded;
    all_pairs.ionic_dispersion = epcsaft::IonicDispersion::AllIonPairs;
    recording_support::expect_group(
        recording_support::active_test_id(), "p6-excluded-observation", "excluded-pure-ion", 1,
        "same-sign pure-ion pair exclusion observation", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "p6-all-observation", "all-pairs-pure-ion", 1,
        "all-ion-pairs pure-ion observation", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "p6-excluded-verdict", "excluded-zero", 1,
        "excluded same-sign pair is zero", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "p6-negative-verdict", "all-pairs-negative", 1,
        "all-ion-pairs contribution is negative", "state");
    recording_support::expect_group(
        recording_support::active_test_id(), "p6-mixture-verdict", "half-mixture-different", 1,
        "ionic pair selection changes the half mixture", "state");

    const auto one_excluded = evaluate(excluded, 300.0, 1.2e4, {0.0, 1.0});
    const auto one_all = evaluate(all_pairs, 300.0, 1.2e4, {0.0, 1.0});
    const auto inputs = recording_support::fields({
        model_inputs(excluded, 300.0, 1.2e4, {0.0, 1.0}),
        recording_support::field("policy_excluded", "ExcludeSameSign"),
        recording_support::field("policy_all", "AllIonPairs"),
    });
    recording_support::observe(
        "excluded-pure-ion", "dispersion", "1", one_excluded.dispersion, 0.0, inputs,
        "Relationship", "ionic-pair-selection", "same-sign pair is excluded", "p6-excluded-observation");
    recording_support::observe(
        "all-pairs-pure-ion", "dispersion", "1", one_all.dispersion, 0.0, inputs,
        "Relationship", "ionic-pair-selection", "all ion pairs include the pair", "p6-all-observation");
    const bool excluded_zero = one_excluded.dispersion == 0.0;
    recording_support::compare(
        "excluded-zero", "dispersion_zero", "bool", excluded_zero ? 1.0 : 0.0, 1.0,
        recording_support::fields({inputs, recording_support::field("dispersion", one_excluded.dispersion)}),
        "Relationship", "ionic-pair-selection", "same-sign pair contributes zero",
        "absolute", 0.0, 0.0, excluded_zero, "p6-excluded-verdict", "state", "a^{disp}=0",
        "T,\\rho,\\mathbf{x},a^{disp}", "a^{disp}=0", "state");
    CHECK(excluded_zero);
    const bool all_negative = one_all.dispersion < 0.0;
    recording_support::compare(
        "all-pairs-negative", "dispersion_negative", "bool", all_negative ? 1.0 : 0.0, 1.0,
        recording_support::fields({inputs, recording_support::field("dispersion", one_all.dispersion)}),
        "Relationship", "ionic-pair-selection", "all-ion-pairs contribution is negative",
        "absolute", 0.0, 0.0, all_negative, "p6-negative-verdict", "state", "a^{disp}<0",
        "T,\\rho,\\mathbf{x},a^{disp}", "a^{disp}<0", "state");
    CHECK(all_negative);
    const auto half_excluded = evaluate(excluded, 300.0, 1.2e4, {0.5, 0.5});
    const auto half_all = evaluate(all_pairs, 300.0, 1.2e4, {0.5, 0.5});
    const bool half_different = half_excluded.dispersion != half_all.dispersion;
    recording_support::compare(
        "half-mixture-different", "dispersion_policy_difference", "bool",
        half_different ? 1.0 : 0.0, 1.0,
        recording_support::fields({model_inputs(excluded, 300.0, 1.2e4, {0.5, 0.5}),
                                   recording_support::field("excluded_dispersion", half_excluded.dispersion),
                                   recording_support::field("all_pairs_dispersion", half_all.dispersion)}),
        "Relationship", "ionic-pair-selection", "ionic pair policy changes the half mixture",
        "absolute", 0.0, 0.0, half_different, "p6-mixture-verdict", "state", "\\Delta a^{disp}",
        "T,\\rho,\\mathbf{x},a^{disp}_{excluded},a^{disp}_{all}",
        "a^{disp}_{excluded}\\ne a^{disp}_{all}", "state");
    CHECK(half_different);
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
    recording_support::expect_group(
        recording_support::active_test_id(), "p7-density", "density-scaling", 1,
        "packing fraction density scaling", "rho");
    recording_support::expect_group(
        recording_support::active_test_id(), "p7-low-positive", "low-positive", 1,
        "packing fraction positivity", "rho");
    recording_support::expect_group(
        recording_support::active_test_id(), "p7-policy-lower", "policy-lower", 1,
        "close packing policy lower bound", "policy");
    recording_support::expect_group(
        recording_support::active_test_id(), "p7-policy-upper", "policy-upper", 1,
        "close packing policy upper bound", "policy");
    recording_support::compare(
        "density-scaling", "packing_fraction", "1", high, 2.0 * low,
        recording_support::fields({model_inputs(p, 300.0, 1.0e3, x),
                                   recording_support::field("rho_high", 2.0e3),
                                   recording_support::field("eta_low", low),
                                   recording_support::field("eta_high", high)}),
        "Relationship", "density-scaling", "packing fraction is linear in density",
        "reference_or_one", 0.0, 1.0e-14, close(high, 2.0 * low, 1.0e-14),
        "p7-density", "rho", "\\eta", "T,\\rho,\\mathbf{x}",
        "\\eta(2\\rho)=2\\eta(\\rho)", "rho");
    CHECK(close(high, 2.0 * low, 1.0e-14));
    const bool low_positive = low > 0.0;
    recording_support::compare(
        "low-positive", "packing_fraction_positive", "bool", low_positive ? 1.0 : 0.0,
        1.0, recording_support::fields({model_inputs(p, 300.0, 1.0e3, x),
                                         recording_support::field("eta", low)}),
        "Relationship", "geometry", "packing fraction is positive", "absolute", 0.0, 0.0,
        low_positive, "p7-low-positive", "rho", "\\eta>0", "T,\\rho,\\mathbf{x},\\eta",
        "\\eta>0", "rho");
    CHECK(low_positive);
    const bool policy_lower = epcsaft::policy::CLOSE_PACKING_FRACTION.value > 0.74;
    recording_support::compare(
        "policy-lower", "close_packing_fraction", "1", policy_lower ? 1.0 : 0.0, 1.0,
        recording_support::field("policy", epcsaft::policy::CLOSE_PACKING_FRACTION.value),
        "Relationship", "policy-bound", "close packing policy is above 0.74", "absolute",
        0.0, 0.0, policy_lower, "p7-policy-lower", "policy", "\\eta_{max}>0.74", "policy",
        "\\eta_{max}>0.74", "policy");
    CHECK(policy_lower);
    const bool policy_upper = epcsaft::policy::CLOSE_PACKING_FRACTION.value < 0.75;
    recording_support::compare(
        "policy-upper", "close_packing_fraction", "1", policy_upper ? 1.0 : 0.0, 1.0,
        recording_support::field("policy", epcsaft::policy::CLOSE_PACKING_FRACTION.value),
        "Relationship", "policy-bound", "close packing policy is below 0.75", "absolute",
        0.0, 0.0, policy_upper, "p7-policy-upper", "policy", "\\eta_{max}<0.75", "policy",
        "\\eta_{max}<0.75", "policy");
    CHECK(policy_upper);
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
    recording_support::expect_group(
        recording_support::active_test_id(), "p8-gradient", "density-gradient", 1,
        "CppAD density gradient", "rho");
    recording_support::expect_group(
        recording_support::active_test_id(), "p8-value", "tape-value", 1,
        "CppAD value at recording point", "state");

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
    recording_support::compare(
        "density-gradient", "d_a_res_d_rho", "1", gradient[1], difference,
        recording_support::fields({model_inputs(p, temperature, molar_density, x),
                                   recording_support::field("step_rho", step),
                                   recording_support::field("rho_plus", molar_density + step),
                                   recording_support::field("rho_minus", molar_density - step)}),
        "Relationship", "derivative", "CppAD gradient matches centered difference",
        "reference_or_one", 0.0, 1.0e-7,
        close(gradient[1], difference, 1.0e-7), "p8-gradient", "rho", "a_\\rho",
        "T,\\rho,\\mathbf{x},\\Delta\\rho", "a_\\rho=d a/d\\rho", "rho");
    CHECK(close(gradient[1], difference, 1.0e-7));

    // The tape reproduces the double path exactly at the recording point.
    const auto value = tape.Forward(0, point);
    const auto direct = evaluate(p, temperature, molar_density, x);
    recording_support::compare(
        "tape-value", "a_res", "1", value[0], direct.hard_chain + direct.dispersion,
        model_inputs(p, temperature, molar_density, x),
        "Relationship", "AD-evaluation", "tape reproduces double path",
        "reference_or_one", 0.0, 1.0e-14,
        close(value[0], direct.hard_chain + direct.dispersion, 1.0e-14),
        "p8-value", "state", "a^{res}", "T,\\rho,\\mathbf{x}", "a_{AD}=a_{double}", "state");
    CHECK(close(value[0], direct.hard_chain + direct.dispersion, 1.0e-14));
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.contributions.p1", "Hard_sphere_Carnahan_Starling", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P1", "eta-*", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P1", "eta-*:3"}, hard_sphere_matches_carnahan_starling},
        {{"native.contributions.p2", "Residual_vanishes_at_zero_density", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P1", "density-monotonic/*,zero-density-limit", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P2", "density-monotonic/*:4,zero-density-limit:1"}, residual_vanishes_at_zero_density},
        {{"native.contributions.p3", "Component_relabelling_invariance", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P1", "hard-chain,dispersion", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P3", "hard-chain:1,dispersion:1"}, relabelling_components_changes_nothing},
        {{"native.contributions.p4", "Association_mass_action", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P2", "mass-action-*,site-floor-*,site-ceiling-*,site-symmetry", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P4", "mass-action-*:2,site-floor-*:2,site-ceiling-*:2,site-symmetry:1"}, association_solve_satisfies_mass_action},
        {{"native.contributions.p5", "Inactive_terms_vanish", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P3", "association,born,debye-huckel", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P5", "association:1,born:1,debye-huckel:1"}, inactive_terms_contribute_nothing},
        {{"native.contributions.p6", "Same_sign_ionic_exclusion", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P3", "excluded-pure-ion,all-pairs-pure-ion,excluded-zero,all-pairs-negative,half-mixture-different", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P6", "excluded-pure-ion:1,all-pairs-pure-ion:1,excluded-zero:1,all-pairs-negative:1,half-mixture-different:1"}, same_sign_exclusion_is_a_limit},
        {{"native.contributions.p7", "Packing_fraction_bound", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "P1", "density-scaling,low-positive,policy-lower,policy-upper", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P7", "density-scaling:1,low-positive:1,policy-lower:1,policy-upper:1"}, packing_fraction_is_bounded},
        {{"native.contributions.p8", "CppAD_tape_matches_double_path", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "M4", "density-gradient,tape-value", "engine/tests/native/contributions_checks.cpp", "contributions_checks", "P8", "density-gradient:1,tape-value:1"}, kernels_record_a_tape},
    });
}
