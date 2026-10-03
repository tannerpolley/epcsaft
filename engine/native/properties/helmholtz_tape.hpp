#pragma once

#include "contributions/association.hpp"
#include "contributions/born.hpp"
#include "contributions/electrolyte.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/thermochemistry.hpp"

#include <cppad/cppad.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft::detail {

using AD = CppAD::AD<double>;

using TapeParameters = ParameterRecord<AD>;

template <class ParameterSet>
bool charged_support(const std::vector<double>& composition, const ParameterSet& parameters) {
    for (std::size_t i = 0; i < composition.size(); ++i)
        if (composition[i] > 0.0 && parameters.charge_number[i] != 0) return true;
    return false;
}

template <typename Scalar, class ParameterSet>
std::pair<ElectrolyteResult<Scalar>, Scalar> supported_electrostatics(const Scalar& temperature,
    const Scalar& density, double reference_temperature, double reference_density,
    const std::vector<Scalar>& composition, const std::vector<double>& reference,
    const ParameterSet& parameters) {
    const auto electrolyte = electrolyte_contribution(temperature, density, reference_temperature,
        reference_density, composition, reference, parameters);
    return {electrolyte, charged_support(reference, parameters)
        ? born_contribution(temperature, composition, parameters, electrolyte.bulk_relative_permittivity) : Scalar(0)};
}

enum class ParameterFamily {
    SegmentCount,
    SigmaConstant,
    DispersionEnergyOverK,
    Kij,
    Lij,
    AssociationEnergyOverK,
    AssociationVolume,
    BornDiameter, DebyeHuckelDiameter, PackingDiameter, IonSuppressionCoefficient,
    KijReciprocalTemperatureSlope,  // b_ij in k_ij(T) = k_ij + b_ij (1/T - 1/T_ref,ij), K
};
constexpr std::size_t PARAMETER_FAMILY_COUNT = static_cast<std::size_t>(ParameterFamily::KijReciprocalTemperatureSlope) + 1;

// The one family -> storage map, for Parameters (double) and TapeParameters (AD) alike.
template <class ParameterSet>
auto& family_values(ParameterSet& parameters, ParameterFamily family) {
    switch (family) {
    case ParameterFamily::PackingDiameter: return parameters.packing_diameter_angstrom;
    case ParameterFamily::DebyeHuckelDiameter: return parameters.debye_huckel_diameter_angstrom;
    case ParameterFamily::BornDiameter: return parameters.born_diameter_angstrom;
    case ParameterFamily::IonSuppressionCoefficient: return parameters.ion_suppression_coefficient;
    case ParameterFamily::SegmentCount: return parameters.segment_count;
    case ParameterFamily::SigmaConstant: return parameters.sigma_constant;
    case ParameterFamily::DispersionEnergyOverK: return parameters.dispersion_energy_over_k;
    case ParameterFamily::Kij: return parameters.k_ij;
    case ParameterFamily::Lij: return parameters.l_ij;
    case ParameterFamily::AssociationEnergyOverK: return parameters.association_energy_over_k;
    case ParameterFamily::AssociationVolume: return parameters.association_volume;
    case ParameterFamily::KijReciprocalTemperatureSlope: return parameters.k_ij_reciprocal_temperature_slope;
    }
    throw std::invalid_argument("invalid parameter family");
}

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
    for (std::size_t family = 0; family < PARAMETER_FAMILY_COUNT; ++family) {
        const auto item = static_cast<ParameterFamily>(family);
        family_values(result, item) = active_vector(family_values(parameters, item), active, item, coordinates);
    }
    static_cast<FixedParameters&>(result) = parameters;
    return result;
}

struct HelmholtzTape final {
    CppAD::ADFun<double> function;
    PredicateSignature predicates;
    std::size_t regular_count = 0, site_count = 0, ideal_offset = 0;
    bool residual_available = true, ideal_available = false;
    std::vector<double> site_values;

    static std::unique_ptr<HelmholtzTape> record(
        const Parameters& parameters,
        const std::vector<ActiveParameter>& active,
        const std::vector<double>& point,
        const std::vector<double>& site_fractions,
        PredicateSignature predicates,
        const std::vector<IdealCorrelation>& ideal_correlations = {},
        double ideal_pressure = 0.0,
        bool include_residual = true
    ) {
        const std::size_t components = parameters.component_count();
        const std::size_t point_size = point.size();
        const std::size_t sites_count = include_residual ? site_fractions.size() : 0;
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
        if (include_residual) {
            independent.insert(independent.end(), site_fractions.begin(), site_fractions.end());
        }
        CppAD::Independent(independent);

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
        std::vector<AD> dependent;
        if (include_residual) {
            const std::vector<AD> sites(
                independent.end() - static_cast<std::ptrdiff_t>(sites_count), independent.end());
            const auto tape_parameters = make_tape_parameters(
                parameters, active,
                std::vector<AD>(
                    independent.begin() + 2 + static_cast<std::ptrdiff_t>(components),
                    independent.begin() + static_cast<std::ptrdiff_t>(point_size)));
            const auto sigma = segment_sigma(independent[0], tape_parameters);
            const auto diameter = segment_diameter(independent[0], tape_parameters, sigma);
            const auto [electrolyte, born] = supported_electrostatics(
                independent[0], density, point[0], reference_density,
                composition, reference_composition, tape_parameters);
            const auto contributions = residual_helmholtz(
                independent[0], density, composition, tape_parameters, sigma,
                diameter, sites, electrolyte.debye_huckel, born);
            const auto delta = association_delta_values(
                independent[0], density, composition, tape_parameters, sigma,
                diameter, tape_parameters.association_energy_over_k,
                tape_parameters.association_volume);
            const auto residuals = association_residuals(
                density, composition, tape_parameters, sites, delta);
            AD potential = contributions.association;
            for (std::size_t site = 0; site < sites.size(); ++site) {
                const std::size_t component = static_cast<std::size_t>(
                    parameters.site_component[site]);
                potential -= AD(0.5) * composition[component]
                    * parameters.site_multiplicity[site] * residuals[site];
            }
            dependent.push_back(total * (
                contributions.hard_chain + contributions.dispersion
                + contributions.debye_huckel + contributions.born + potential
                + contributions.dipolar));
            dependent.insert(dependent.end(), residuals.begin(), residuals.end());
        }
        if (!ideal_correlations.empty()) {
            std::vector<bool> supported;
            supported.reserve(components);
            for (std::size_t i = 0; i < components; ++i) {
                supported.push_back(point[2 + i] > 0.0);
            }
            const auto ideal = evaluate_ideal_mixture(
                ideal_correlations, ideal_pressure, independent[0], density,
                composition, supported);
            dependent.insert(dependent.end(), {
                ideal.enthalpy, ideal.internal_energy, ideal.entropy,
                ideal.isobaric_heat_capacity, ideal.isochoric_heat_capacity});
        }
        auto result = std::make_unique<HelmholtzTape>();
        result->function = CppAD::ADFun<double>(independent, dependent);
        result->predicates = std::move(predicates);
        result->regular_count = point_size;
        result->site_count = sites_count;
        result->ideal_offset = include_residual ? 1 + sites_count : 0;
        result->residual_available = include_residual;
        result->ideal_available = !ideal_correlations.empty();
        result->site_values = include_residual ? site_fractions : std::vector<double>{};
        return result;
    }

    std::vector<double> implicit_taylor(
        const std::vector<std::vector<double>>& regular_coefficients,
        std::size_t order,
        const Eigen::FullPivLU<Eigen::MatrixXd>* factorization = nullptr,
        std::size_t* solve_count = nullptr,
        std::size_t* forward_count = nullptr,
        double* coefficient_residual = nullptr,
        std::vector<std::vector<double>>* output_prefix = nullptr
    ) {
        const std::size_t input_count = regular_count + site_count;
        std::vector<double> input(input_count), output;
        if (output_prefix != nullptr) output_prefix->clear();
        // Without U, callers may consume only the site-independent Ideal outputs.
        for (std::size_t current = 0; current <= order; ++current) {
            input.assign(input_count, 0.0);
            std::copy(regular_coefficients[current].begin(), regular_coefficients[current].end(), input.begin());
            if (!current) {
                std::copy(site_values.begin(), site_values.end(), input.begin() + regular_count);
                output = function.Forward(0, input);
                if (forward_count != nullptr) ++*forward_count;
                if (output_prefix != nullptr) output_prefix->push_back(output);
                continue;
            }
            const auto known = function.Forward(current, input);
            if (forward_count != nullptr) ++*forward_count;
            if (site_count && factorization != nullptr) {
                std::vector<double> rhs(site_count);
                for (std::size_t site = 0; site < site_count; ++site)
                    rhs[site] = -known[1 + site];
                const Eigen::VectorXd coefficient = factorization->solve(
                    Eigen::Map<const Eigen::VectorXd>(rhs.data(),
                        static_cast<Eigen::Index>(site_count)));
                for (std::size_t site = 0; site < site_count; ++site)
                    input[regular_count + site] = coefficient[static_cast<Eigen::Index>(site)];
                if (solve_count != nullptr) ++*solve_count;
            }
            output = known;
            if (site_count && factorization != nullptr) {
                output = function.Forward(current, input);
                if (forward_count != nullptr) ++*forward_count;
            }
            if (output_prefix != nullptr) output_prefix->push_back(output);
            if (coefficient_residual != nullptr && site_count != 0 && factorization != nullptr) {
                double scale = 1.0;
                for (std::size_t site = 0; site < site_count; ++site)
                    scale = std::max(scale, std::abs(known[1 + site]));
                for (std::size_t site = 0; site < site_count; ++site) {
                    const double residual = std::isfinite(known[1 + site])
                        && std::isfinite(output[1 + site])
                        ? std::abs(output[1 + site]) / scale
                        : std::numeric_limits<double>::infinity();
                    *coefficient_residual = std::max(*coefficient_residual, residual);
                }
            }
        }
        return output;
    }
};

}  // namespace epcsaft::detail
