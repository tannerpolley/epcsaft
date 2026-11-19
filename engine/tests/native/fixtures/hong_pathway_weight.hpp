#pragma once
// Hong et al. 2027, Tables 1–3, and the declared two-block verification states.
#include "surface/parameter_builder.hpp"

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace hong_pathway_weight {
using epcsaft::ActiveParameter;
using epcsaft::DensityBranch;
using epcsaft::ParameterEntry;
using epcsaft::ParameterFamily;

constexpr double T_R = 298.15;

// Single, Symmetric (nu = 2, ordinary association) and Corrected (nu = 2 with omega_2) are one block each;
// Blend weighs a corrected symmetric block by w and a single-site block by 1 - w; Uncorrected drops the correction.
enum class Form { Single, Symmetric, Corrected, Blend, Uncorrected };

struct Pathway final {
    double k_slope, k_intercept, kappa, c, d;  // k = slope T + intercept (Eq. 12); eps/k = c T + d (Eq. 13)
    std::size_t nu;
};
constexpr Pathway SYMMETRIC{3.675e-3, -2.3981, 1.555e-3, 4.1160, 3546.55, 2};
constexpr Pathway SINGLE{5.939e-3, -1.9778, 6.5e-5, 12.5398, 1197.42, 1};

struct Weight final {
    double intercept, loading, temperature;  // b_0, b_alpha, b_T
};
constexpr Weight LOGISTIC{1.7, -10.0, -1.5};  // illustrative, shaped on Table S1, not fitted

inline ParameterEntry entry(std::string identity, std::string field, std::vector<std::size_t> index, double magnitude,
                     std::string unit) {
    return {std::move(identity), std::move(field), std::move(index), magnitude, std::move(unit)};
}

// q = a_0 + a_1 (T - T_r) for q = slope T + intercept.
inline void linear_law(std::vector<ParameterEntry>& e, const std::string& field, std::size_t law, std::size_t a, std::size_t b,
                double slope, double intercept) {
    const bool energy = field == "association_energy_over_k";
    const std::string name = "law/" + std::to_string(law);
    e.push_back(entry(name, field + "_temperature_law", {law, a, b}, T_R, "kelvin"));
    e.push_back(entry(name + "/a0", "temperature_law_coefficient", {law, 0}, intercept + slope * T_R, energy ? "kelvin" : "dimensionless"));
    e.push_back(entry(name + "/a1", "temperature_law_coefficient", {law, 1}, slope, energy ? "dimensionless" : "1 / kelvin"));
}

// Components MEA (2B), water (4C), CO2; per block sites MEA a, b, water a, b, CO2 r. order[p] is the component
// declared at position p (the same model in another component order).
using Order = std::array<std::size_t, 3>;
inline epcsaft::Parameters hong(Form form, Weight weight = LOGISTIC, Order order = {0, 1, 2}) {
    Order at{};
    for (std::size_t p = 0; p < 3; ++p) at[order[p]] = p;
    epcsaft::ParameterDocument document;
    for (const auto component : order) document.component_ids.push_back(std::array{"mea", "water", "co2"}[component]);
    document.base_choice = "pc-saft";
    document.association_choice = "general-site";
    document.electrolyte_choice = "none";
    document.permittivity_choice = "none";
    auto& e = document.entries;
    const std::array<double, 3> m = {3.1200, 2.1945, 2.0729};
    const std::array<double, 3> sigma = {3.0229, 2.2290, 2.7852};
    const std::array<double, 3> epsilon = {287.94, 141.66, 169.21};
    const std::array<double, 3> mass = {0.06108352, 0.01801528, 0.0440095};
    for (std::size_t i = 0; i < 3; ++i) {
        e.push_back(entry("", "molar_mass", {at[i]}, mass[i], "kilogram / mole"));
        e.push_back(entry("", "charge_number", {at[i]}, 0.0, "elementary-charge"));
        e.push_back(entry("", "segment_count", {at[i]}, m[i], "dimensionless"));
        e.push_back(entry("", "sigma_constant", {at[i]}, sigma[i], "angstrom"));
        e.push_back(entry("", "dispersion_energy_over_k", {at[i]}, epsilon[i], "kelvin"));
    }
    const bool blend = form == Form::Blend || form == Form::Uncorrected;
    const std::vector<Pathway> blocks = blend ? std::vector<Pathway>{SYMMETRIC, SINGLE}
        : std::vector<Pathway>{form == Form::Single ? SINGLE : SYMMETRIC};
    linear_law(e, "k_ij", 0, at[0], at[1], 9.0e-5, -0.0603);  // Table 2, x_MEA <= 0.20 set
    linear_law(e, "k_ij", 1, at[1], at[2], 2.46e-4, -0.0633);
    linear_law(e, "k_ij", 2, at[0], at[2], blocks[0].k_slope, blocks[0].k_intercept);
    std::size_t law = 3;
    if (blend) linear_law(e, "k_ij_complement", law++, at[0], at[2], SINGLE.k_slope, SINGLE.k_intercept);
    // Eq. 6: eps_MW = alpha_HB sqrt(eps_MEA eps_W), one power-0 law on every MEA-water edge; Eq. 5 kappa is fixed.
    const std::size_t cross = law++;
    e.push_back(entry("law/cross/a0", "temperature_law_coefficient", {cross, 0}, 1.0374 * std::sqrt(2355.03 * 1804.17), "kelvin"));
    const double kappa_cross = std::pow((std::cbrt(0.0313) + std::cbrt(0.2039)) / 2.0, 3.0);
    for (std::size_t block = 0; block < blocks.size(); ++block) {
        const std::size_t s = 5 * block;
        for (const std::size_t component : {0, 0, 1, 1, 2}) document.site_component.push_back(static_cast<int>(at[component]));
        for (const std::size_t nu : {std::size_t{1}, std::size_t{1}, std::size_t{2}, std::size_t{2}, blocks[block].nu})
            document.site_multiplicity.push_back(static_cast<int>(nu));
        e.push_back(entry("", "association_energy_over_k", {s, s + 1}, 2355.03, "kelvin"));
        e.push_back(entry("", "association_volume", {s, s + 1}, 0.0313, "dimensionless"));
        e.push_back(entry("", "association_energy_over_k", {s + 2, s + 3}, 1804.17, "kelvin"));
        e.push_back(entry("", "association_volume", {s + 2, s + 3}, 0.2039, "dimensionless"));
        for (const auto& [a, b] : {std::pair{s, s + 3}, std::pair{s + 1, s + 2}}) {
            e.push_back(entry("law/cross", "association_energy_over_k_temperature_law", {cross, a, b}, T_R, "kelvin"));
            e.push_back(entry("", "association_volume", {a, b}, kappa_cross, "dimensionless"));
        }
        linear_law(e, "association_energy_over_k", law++, s + 4, s + 1, blocks[block].c, blocks[block].d);
        e.push_back(entry("", "association_volume", {s + 4, s + 1}, blocks[block].kappa, "dimensionless"));
        for (std::size_t site = s; blend && site < s + 5; ++site)
            e.push_back(entry("", "site_pathway", {site}, static_cast<double>(block + 1), "dimensionless"));
    }
    if (form == Form::Corrected || form == Form::Blend) e.push_back(entry("", "site_correction", {4}, 1.0, "dimensionless"));
    e.push_back(entry("", "no_cross_association", {at[1], at[2]}, 0.0, "dimensionless"));
    if (blend) {
        e.push_back(entry("", "pathway_loading", {at[2], at[0]}, 0.0, "dimensionless"));
        e.push_back(entry("b0", "pathway_weight", {0}, weight.intercept, "dimensionless"));
        e.push_back(entry("b_alpha", "pathway_weight", {1}, weight.loading, "dimensionless"));
        e.push_back(entry("b_T", "pathway_weight", {2}, weight.temperature, "dimensionless"));
    }
    return epcsaft::build_parameters(document);
}

// The 16 blend coordinates: k_MW a0, a1; k_sym a0, a1; k_single a0, a1; eps_MW; eps_sym a0, a1; eps_single a0, a1;
// kappa_sym; kappa_single; b_0, b_alpha, b_T. A single-block model has k_MW, k_MEA-CO2, eps_MW, eps_r and kappa_r.
inline std::vector<ActiveParameter> coordinates(Form form, bool weights = true) {
    std::vector<ActiveParameter> result;
    const bool blend = form == Form::Blend || form == Form::Uncorrected;
    for (const std::size_t slot : blend ? std::vector<std::size_t>{0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12}
                                        : std::vector<std::size_t>{0, 1, 4, 5, 6, 7, 8})
        result.push_back({ParameterFamily::TemperatureLawCoefficient, {slot}});
    if (!blend) {
        result.push_back({ParameterFamily::AssociationVolume, {9, 21}});
        return result;
    }
    result.push_back({ParameterFamily::AssociationVolume, {14, 41}});
    result.push_back({ParameterFamily::AssociationVolume, {69, 96}});
    for (std::size_t k = 0; weights && k < 3; ++k) result.push_back({ParameterFamily::PathwayWeight, {k}});
    return result;
}
struct StatePoint final {
    const char* name;
    double temperature, density;
    std::vector<double> x;
    DensityBranch branch;
};
// The six component orders of the same model, for the double rounding estimate e_A.
constexpr std::array<Order, 6> ORDERS{{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
inline const std::vector<StatePoint> STATES = {
    {"L", 333.15, 42424.79, {0.10404183733215595, 0.8231288765353348, 0.07282928613250916}, DensityBranch::Liquid},
    {"V", 353.15, 34.49684, {1.0e-5, 0.6, 0.39999}, DensityBranch::Vapor},
    {"Z", 333.15, 43629.45, {0.11221432663480933, 0.8877856733651908, 0.0}, DensityBranch::Liquid}};

inline epcsaft::Parameters shifted(epcsaft::Parameters p, const ActiveParameter& coordinate, double step) {
    auto values = epcsaft::detail::family_values(p, coordinate.family);
    for (const auto index : coordinate.indices) values.first[index] += step;
    return p;
}

}  // namespace hong_pathway_weight
