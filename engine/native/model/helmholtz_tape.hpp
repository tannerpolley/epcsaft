#pragma once

#include "physics/association.hpp"
#include "physics/born.hpp"
#include "physics/electrolyte.hpp"
#include "physics/helmholtz.hpp"

#include <cppad/cppad.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace epcsaft::detail {

// Source: TRANSFER.md, Stage 2; the potential form keeps site fractions as
// dynamic parameters and exposes Phi plus association residuals on one tape.
using AD = CppAD::AD<double>;

struct TapeParameters final {
    std::vector<AD> segment_count;
    std::vector<AD> sigma_constant;
    std::vector<int> sigma_term_component;
    std::vector<double> sigma_term_amplitude;
    std::vector<double> sigma_term_exponent;
    std::vector<AD> dispersion_energy_over_k;
    std::vector<double> molar_mass_kg_per_mol;
    std::vector<AD> k_ij;
    std::vector<double> k_ij_reciprocal_temperature_slope;
    std::vector<double> k_ij_reference_temperature_k;
    std::vector<double> l_ij;
    std::vector<int> site_component;
    std::vector<int> site_multiplicity;
    std::vector<int> association_active;
    std::vector<AD> association_energy_over_k;
    std::vector<AD> association_volume;
    std::vector<int> charge_number;
    std::vector<double> packing_diameter_angstrom;
    std::vector<double> debye_huckel_diameter_angstrom;
    std::vector<double> born_diameter_angstrom;
    std::vector<double> solvation_factor;
    double ionic_region_relative_permittivity;
    double born_c_shell;
    double born_c_dielectric;
    ElectrolyteForm electrolyte_form;
    PermittivityForm permittivity_form;
    std::vector<double> relative_permittivity;
    std::vector<int> permittivity_term_component;
    std::vector<double> permittivity_term_coefficient;
    std::vector<int> permittivity_term_power;
    std::vector<double> permittivity_reference_temperature_k;
    std::vector<double> ion_suppression_coefficient;
    std::vector<double> dielectric_volume_m3_per_mol;
    std::vector<double> dielectric_temperature_k;
    std::vector<double> dipole_scaling_debye2;
    std::vector<double> polarizability_scaling_angstrom3;
    std::vector<double> correlation_integral_parameter;
};

enum class ParameterFamily {
    SegmentCount,
    SigmaConstant,
    DispersionEnergyOverK,
    Kij,
    AssociationEnergyOverK,
    AssociationVolume,
};

struct ActiveParameter final {
    ParameterFamily family;
    std::vector<std::size_t> indices;
};

using PredicateSignature = std::vector<unsigned char>;

inline std::vector<AD> active_vector(
    const std::vector<double>& values,
    const std::vector<ActiveParameter>& active,
    ParameterFamily family,
    const std::vector<AD>& coordinates
) {
    std::vector<AD> result(values.begin(), values.end());
    for (std::size_t parameter = 0; parameter < active.size(); ++parameter) {
        const auto& item = active[parameter];
        if (item.family != family) {
            continue;
        }
        for (const std::size_t index : item.indices) {
            result[index] = coordinates[parameter];
        }
    }
    return result;
}

inline TapeParameters make_tape_parameters(
    const Parameters& parameters,
    const std::vector<ActiveParameter>& active,
    const std::vector<AD>& coordinates
) {
    TapeParameters result;
    result.segment_count = active_vector(
        parameters.segment_count, active, ParameterFamily::SegmentCount,
        coordinates
    );
    result.sigma_constant = active_vector(
        parameters.sigma_constant, active, ParameterFamily::SigmaConstant,
        coordinates
    );
    result.dispersion_energy_over_k = active_vector(
        parameters.dispersion_energy_over_k,
        active,
        ParameterFamily::DispersionEnergyOverK,
        coordinates
    );
    result.k_ij = active_vector(
        parameters.k_ij, active, ParameterFamily::Kij, coordinates
    );
    result.association_energy_over_k = active_vector(
        parameters.association_energy_over_k,
        active,
        ParameterFamily::AssociationEnergyOverK,
        coordinates
    );
    result.association_volume = active_vector(
        parameters.association_volume,
        active,
        ParameterFamily::AssociationVolume,
        coordinates
    );
    result.sigma_term_component = parameters.sigma_term_component;
    result.sigma_term_amplitude = parameters.sigma_term_amplitude;
    result.sigma_term_exponent = parameters.sigma_term_exponent;
    result.molar_mass_kg_per_mol = parameters.molar_mass_kg_per_mol;
    result.k_ij_reciprocal_temperature_slope =
        parameters.k_ij_reciprocal_temperature_slope;
    result.k_ij_reference_temperature_k = parameters.k_ij_reference_temperature_k;
    result.l_ij = parameters.l_ij;
    result.site_component = parameters.site_component;
    result.site_multiplicity = parameters.site_multiplicity;
    result.association_active = parameters.association_active;
    result.charge_number = parameters.charge_number;
    result.packing_diameter_angstrom = parameters.packing_diameter_angstrom;
    result.debye_huckel_diameter_angstrom = parameters.debye_huckel_diameter_angstrom;
    result.born_diameter_angstrom = parameters.born_diameter_angstrom;
    result.solvation_factor = parameters.solvation_factor;
    result.ionic_region_relative_permittivity =
        parameters.ionic_region_relative_permittivity;
    result.born_c_shell = parameters.born_c_shell;
    result.born_c_dielectric = parameters.born_c_dielectric;
    result.electrolyte_form = parameters.electrolyte_form;
    result.permittivity_form = parameters.permittivity_form;
    result.relative_permittivity = parameters.relative_permittivity;
    result.permittivity_term_component = parameters.permittivity_term_component;
    result.permittivity_term_coefficient = parameters.permittivity_term_coefficient;
    result.permittivity_term_power = parameters.permittivity_term_power;
    result.permittivity_reference_temperature_k =
        parameters.permittivity_reference_temperature_k;
    result.ion_suppression_coefficient = parameters.ion_suppression_coefficient;
    result.dielectric_volume_m3_per_mol = parameters.dielectric_volume_m3_per_mol;
    result.dielectric_temperature_k = parameters.dielectric_temperature_k;
    result.dipole_scaling_debye2 = parameters.dipole_scaling_debye2;
    result.polarizability_scaling_angstrom3 =
        parameters.polarizability_scaling_angstrom3;
    result.correlation_integral_parameter = parameters.correlation_integral_parameter;
    return result;
}

struct HelmholtzTape final {
    CppAD::ADFun<double> function;
    PredicateSignature predicates;
    std::size_t site_count = 0;

    static std::unique_ptr<HelmholtzTape> record(
        const Parameters& parameters,
        const std::vector<ActiveParameter>& active,
        const std::vector<double>& point,
        const std::vector<double>& site_fractions,
        PredicateSignature predicates
    ) {
        const std::size_t components = parameters.component_count();
        double total_value = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            total_value += point[2 + component];
        }
        const double reference_density = total_value / point[1];
        std::vector<double> reference_composition;
        reference_composition.reserve(components);
        for (std::size_t component = 0; component < components; ++component) {
            reference_composition.push_back(point[2 + component] / total_value);
        }
        std::vector<AD> independent(point.begin(), point.end());
        std::vector<AD> dynamic(site_fractions.begin(), site_fractions.end());
        CppAD::Independent(independent, dynamic);

        std::vector<AD> amounts(
            independent.begin() + 2,
            independent.begin() + 2 + static_cast<std::ptrdiff_t>(components)
        );
        AD total = AD(0.0);
        for (const AD& amount : amounts) {
            total += amount;
        }
        AD density = total / independent[1];
        std::vector<AD> composition;
        composition.reserve(components);
        for (const AD& amount : amounts) {
            composition.push_back(amount / total);
        }
        const std::vector<AD> sites = dynamic;
        const auto tape_parameters = make_tape_parameters(
            parameters,
            active,
            std::vector<AD>(
                independent.begin() + 2 + static_cast<std::ptrdiff_t>(components),
                independent.end()
            )
        );
        const auto sigma = segment_sigma(independent[0], tape_parameters);
        const auto diameter = segment_diameter(
            independent[0], tape_parameters, sigma
        );
        const auto electrolyte = electrolyte_contribution(
            independent[0], density, point[0], reference_density,
            composition, reference_composition,
            parameters
        );
        const auto born = born_contribution(
            independent[0], composition, parameters,
            electrolyte.bulk_relative_permittivity
        );
        const auto contributions = residual_helmholtz(
            independent[0], density, composition, tape_parameters, sigma,
            diameter, sites, electrolyte.debye_huckel, born
        );
        const auto delta = association_delta_values(
            independent[0], density, composition, tape_parameters, sigma,
            diameter, tape_parameters.association_energy_over_k,
            tape_parameters.association_volume
        );
        const auto residuals = association_residuals(
            density, composition, tape_parameters, sites, delta
        );
        AD potential = contributions.association;
        for (std::size_t site = 0; site < sites.size(); ++site) {
            const std::size_t component = static_cast<std::size_t>(
                parameters.site_component[site]
            );
            potential -= AD(0.5) * composition[component]
                * parameters.site_multiplicity[site] * residuals[site];
        }
        AD phi = total * (
            contributions.hard_chain + contributions.dispersion
            + contributions.debye_huckel + contributions.born + potential
        );
        std::vector<AD> dependent = {phi};
        for (std::size_t site = 0; site < sites.size(); ++site) {
            const std::size_t component = static_cast<std::size_t>(
                parameters.site_component[site]
            );
            const AD multiplicity = AD(
                composition[component] * parameters.site_multiplicity[site]
            );
            dependent.push_back(
                -total * multiplicity / sites[site] * residuals[site]
            );
        }
        auto result = std::make_unique<HelmholtzTape>();
        result->function = CppAD::ADFun<double>(independent, dependent);
        result->predicates = std::move(predicates);
        result->site_count = sites.size();
        return result;
    }
};

}  // namespace epcsaft::detail
