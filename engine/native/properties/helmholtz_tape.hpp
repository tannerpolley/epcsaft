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
#include <cstring>
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

template <typename Scalar, class ParameterSet, typename Fixed>
std::pair<ElectrolyteResult<Scalar>, Scalar> supported_electrostatics(const Fixed& temperature,
    const Scalar& density,
    const std::vector<Fixed>& composition, const std::vector<double>& reference,
    const ParameterSet& parameters) {
    const auto electrolyte = electrolyte_contribution(temperature, density, composition, reference, parameters);
    return {electrolyte, charged_support(reference, parameters)
        ? born_contribution(temperature, composition, parameters, electrolyte.bulk_relative_permittivity) : Scalar(0)};
}

// ALG-ASSOCIATION: the stationary potential eliminates the highest implicit site coefficient.
template <typename Scalar, typename Fixed, typename ParameterSet>
Scalar stationary_association(Scalar potential, const std::vector<Fixed>& composition,
    const ParameterSet& parameters, const std::vector<Scalar>& residuals) {
    for (std::size_t site = 0; site < residuals.size(); ++site)
        potential -= Scalar(0.5) * composition[static_cast<std::size_t>(parameters.site_component[site])]
            * parameters.site_multiplicity[site] * residuals[site];
    return potential;
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
    BornCDielectric, IonicRegionRelativePermittivity,
    KijReciprocalTemperatureSlope,  // b_ij in k_ij(T) = k_ij + b_ij (1/T - 1/T_ref,ij), K
    DielectricVolume, DielectricTemperature,  // Schreckenberg d_V (m3/mol) and d_T (K) of a neutral solvent
    CationPairMixing, CationTripleMixing,
    TemperatureLawCoefficient,  // a_p or s of a declared temperature law, by coefficient slot
};
constexpr std::size_t PARAMETER_FAMILY_COUNT = static_cast<std::size_t>(ParameterFamily::TemperatureLawCoefficient) + 1;

// The one family -> (pointer, size) storage view, for double and AD alike.
template <class ParameterSet>
auto family_values(ParameterSet& parameters, ParameterFamily family, std::size_t size = 0) {
    const auto view = [size](auto& values) {
        if constexpr (!std::is_const_v<ParameterSet>) if (size) values.resize(size);
        return std::pair{values.data(), values.size()};
    };
    switch (family) {
    case ParameterFamily::BornCDielectric: return std::pair{&parameters.born_c_dielectric, std::size_t{1}};
    case ParameterFamily::IonicRegionRelativePermittivity: return std::pair{&parameters.ionic_region_relative_permittivity, std::size_t{1}};
    case ParameterFamily::PackingDiameter: return view(parameters.packing_diameter_angstrom);
    case ParameterFamily::DebyeHuckelDiameter: return view(parameters.debye_huckel_diameter_angstrom);
    case ParameterFamily::BornDiameter: return view(parameters.born_diameter_angstrom);
    case ParameterFamily::IonSuppressionCoefficient: return view(parameters.ion_suppression_coefficient);
    case ParameterFamily::SegmentCount: return view(parameters.segment_count);
    case ParameterFamily::SigmaConstant: return view(parameters.sigma_constant);
    case ParameterFamily::DispersionEnergyOverK: return view(parameters.dispersion_energy_over_k);
    case ParameterFamily::Kij: return view(parameters.k_ij);
    case ParameterFamily::Lij: return view(parameters.l_ij);
    case ParameterFamily::AssociationEnergyOverK: return view(parameters.association_energy_over_k);
    case ParameterFamily::AssociationVolume: return view(parameters.association_volume);
    case ParameterFamily::KijReciprocalTemperatureSlope: return view(parameters.k_ij_reciprocal_temperature_slope);
    case ParameterFamily::DielectricVolume: return view(parameters.dielectric_volume_m3_per_mol);
    case ParameterFamily::DielectricTemperature: return view(parameters.dielectric_temperature_k);
    case ParameterFamily::TemperatureLawCoefficient: return view(parameters.temperature_law_coefficient);
    case ParameterFamily::CationPairMixing: return view(parameters.cation_pair_mixing);
    case ParameterFamily::CationTripleMixing: return view(parameters.cation_triple_mixing);
    }
    throw std::invalid_argument("invalid parameter family");
}

struct ActiveParameter final {
    ParameterFamily family;
    std::vector<std::size_t> indices;
};

using PredicateSignature = std::vector<unsigned char>;

template<typename Scalar = AD>
inline ParameterRecord<Scalar> make_tape_parameters(
    const Parameters& parameters,
    const std::vector<ActiveParameter>& active,
    const std::vector<Scalar>& coordinates
) {
    ParameterRecord<Scalar> result;
    for (std::size_t family = 0; family < PARAMETER_FAMILY_COUNT; ++family) {
        const auto item = static_cast<ParameterFamily>(family);
        const auto source = family_values(parameters, item);
        const auto target = family_values(result, item, source.second);
        for (std::size_t i = 0; i < source.second; ++i) target.first[i] = Scalar(source.first[i]);
    }
    for (std::size_t i = 0; i < active.size(); ++i)
        for (const auto index : active[i].indices) family_values(result, active[i].family).first[index] = coordinates[i];
    static_cast<FixedParameters&>(result) = parameters;
    return result;
}

struct HelmholtzTape final {
    CppAD::ADFun<double> function;
    PredicateSignature predicates;
    std::size_t regular_count = 0, site_count = 0, ideal_offset = 0;
    bool residual_available = true, ideal_available = false;
    std::vector<double> site_values;
    std::vector<double> zero_point, zero_value;

    static std::unique_ptr<HelmholtzTape> record(
        const Parameters& parameters,
        const std::vector<ActiveParameter>& active,
        const std::vector<double>& point,
        const std::vector<double>& site_fractions,
        PredicateSignature predicates,
        const std::vector<IdealCorrelation>& ideal_correlations = {},
        double ideal_pressure = 0.0,
        bool include_residual = true,
        bool reference_insertion = false
    ) {
        const std::size_t components = parameters.component_count();
        const std::size_t point_size = point.size();
        const std::size_t sites_count = include_residual ? site_fractions.size() : 0;
        double total_value = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            total_value += point[2 + component];
        }
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
            auto [electrolyte, born] = supported_electrostatics(
                independent[0], density,
                composition, reference_composition, tape_parameters);
            // EqIDs born_residual_helmholtz, reference_single_insertion_taylor: DH is zero;
            // retain the Born first insertion only on the separate ion-free reference tape.
            if (reference_insertion) born = born_contribution(independent[0], composition,
                tape_parameters, electrolyte.bulk_relative_permittivity);
            const auto contributions = residual_helmholtz(
                independent[0], density, composition, tape_parameters, sigma,
                diameter, sites, electrolyte.debye_huckel, born);
            const auto delta = parameters.ideal_gas ? std::vector<AD>{} : association_delta_values(
                independent[0], density, composition, tape_parameters, sigma,
                diameter, tape_parameters.association_energy_over_k,
                tape_parameters.association_volume);
            const auto residuals = association_residuals(
                density, composition, tape_parameters, sites, delta);
            const auto potential = stationary_association(contributions.association, composition, parameters, residuals);
            dependent.push_back(total * (
                contributions.hard_chain + contributions.dispersion
                + contributions.debye_huckel + contributions.born + potential
                + contributions.dipolar));
            if (!parameters.cation_mixing.empty()) dependent.back() += total * contributions.cation_mixing;
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
        result->function.optimize("no_cumulative_sum_op");
        result->predicates = std::move(predicates);
        result->regular_count = point_size;
        result->site_count = sites_count;
        result->ideal_offset = include_residual ? 1 + sites_count : 0;
        result->residual_available = include_residual;
        result->ideal_available = !ideal_correlations.empty();
        result->site_values = include_residual ? site_fractions : std::vector<double>{};
        return result;
    }

    std::vector<double> forward_zero(const std::vector<double>& point, std::size_t* forward_count = nullptr) {
        if (zero_point.size() != point.size()
            || std::memcmp(zero_point.data(), point.data(), point.size() * sizeof(double)) != 0) {
            zero_point.clear();
            zero_value = function.Forward(0, point);
            zero_point = point;
            if (forward_count != nullptr) ++*forward_count;
        }
        return zero_value;
    }

    std::vector<double> implicit_taylor(
        const std::vector<std::vector<double>>& regular_coefficients,
        std::size_t order,
        const Eigen::FullPivLU<Eigen::MatrixXd>* factorization = nullptr,
        std::size_t* solve_count = nullptr,
        std::size_t* forward_count = nullptr,
        double* coefficient_residual = nullptr,
        std::vector<std::vector<double>>* output_prefix = nullptr,
        std::vector<double>* site_coefficients = nullptr, std::vector<double>* right_coefficients = nullptr
    ) {
        core::Scope scope("implicit_taylor"), evaluation("evaluate");
        scope.probe("derivative_order", static_cast<double>(order));
        const std::size_t input_count = regular_count + site_count;
        std::vector<double> input(input_count), output;
        if (output_prefix != nullptr) output_prefix->clear();
        if (site_coefficients) site_coefficients->clear();
        if (right_coefficients) right_coefficients->clear();
        // Without U, callers may consume only the site-independent Ideal outputs.
        for (std::size_t current = 0; current <= order; ++current) {
            input.assign(input_count, 0.0);
            std::copy(regular_coefficients[current].begin(), regular_coefficients[current].end(), input.begin());
            if (!current) {
                std::copy(site_values.begin(), site_values.end(), input.begin() + regular_count);
                output = forward_zero(input, forward_count);
                if (output_prefix != nullptr) output_prefix->push_back(output);
                if (site_coefficients) site_coefficients->insert(site_coefficients->end(), input.begin() + regular_count, input.end());
                if (right_coefficients) for (std::size_t a = 0; a < site_count; ++a) right_coefficients->push_back(-output[1+a]);
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
            if (site_coefficients) site_coefficients->insert(site_coefficients->end(), input.begin() + regular_count, input.end());
            if (right_coefficients) for (std::size_t a = 0; a < site_count; ++a) right_coefficients->push_back(-known[1+a]);
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
