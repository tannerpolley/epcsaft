#pragma once

#include "contributions/association.hpp"
#include "contributions/born.hpp"
#include "contributions/electrolyte.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_model.hpp"

#include <cppad/cppad.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft::detail {

using AD = CppAD::AD<double>;

struct RecordingOwner final {
    bool active = true;
    ~RecordingOwner() { if (active) AD::abort_recording(); }
};

using TapeParameters = ParameterRecord<AD>;

template <class ParameterSet>
bool charged_support(const std::vector<double>& composition, const ParameterSet& parameters) {
    for (std::size_t i = 0; i < composition.size(); ++i)
        if (composition[i] > 0.0 && parameters.charge_number[i] != 0) return true;
    return false;
}

template <typename Fixed>
struct SegmentGeometry final { std::vector<Fixed> sigma, diameter; };

// Segment sigma and the temperature-dependent diameter, admitted finite and positive (series by their lead value).
template <typename Fixed, class ParameterSet>
SegmentGeometry<Fixed> admitted_segment_geometry(const Fixed& temperature, const ParameterSet& parameters) {
    const auto require_positive = [](const std::vector<Fixed>& values, const char* message) {
        if (!std::all_of(values.begin(), values.end(), [](const Fixed& value) {
                const auto lead = scalar_lead(value);
                return wide_math::isfinite(lead) && lead > 0.0; }))
            throw std::invalid_argument(message);
    };
    SegmentGeometry<Fixed> result;
    result.sigma = segment_sigma(temperature, parameters);
    require_positive(result.sigma, "invalid sigma correlation");
    result.diameter = segment_diameter(temperature, parameters, result.sigma);
    require_positive(result.diameter, "invalid segment diameter");
    return result;
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
    const HelmholtzModel* quotient_model = nullptr;

    static std::unique_ptr<HelmholtzTape> record(
        const HelmholtzModel& model,
        const std::vector<double>& point,
        const std::vector<double>& site_fractions,
        PredicateSignature predicates,
        const std::vector<IdealCorrelation>& ideal_correlations = {},
        double ideal_pressure = 0.0,
        bool include_residual = true,
        bool reference_insertion = false
    ) {
        const auto& parameters = model.parameters();
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
        RecordingOwner recording;

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
                parameters, model.active_parameters(),
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
            const auto delta = parameters.ideal_gas ? std::vector<AD>{} : association_delta(
                independent[0], density, composition, tape_parameters, sigma,
                diameter);
            const auto products = association_products(density, composition, tape_parameters, sites, delta);
            const auto& residuals = products.second;
            dependent.push_back(extensive_residual_helmholtz(total, residual_helmholtz(independent[0], density, composition,
                tape_parameters, sigma, diameter, sites, products), electrolyte.debye_huckel, born, parameters));
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
        recording.active = false;
        result->function.optimize("no_cumulative_sum_op");
        result->predicates = std::move(predicates);
        if (include_residual && density_temperature_ion_suppression(parameters.permittivity_form))
            result->quotient_model = &model;
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
            if (quotient_model) {
                const auto& parameters = quotient_model->parameters();
                const auto& active = quotient_model->active_parameters();
                const auto end = point.begin() + 2 + parameters.component_count();
                std::vector<double> composition(point.begin() + 2, end);
                const double total = std::accumulate(composition.begin(), composition.end(), 0.0);
                for (auto& x : composition) x /= total;
                const auto evaluate = [&](const Parameters& mapped) {
                    return bulk_relative_permittivity(point[0], total / point[1], composition, mapped);
                };
                if (active.empty()) evaluate(parameters);
                else evaluate(make_tape_parameters<double>(parameters, active, {end, point.begin() + regular_count}));
            }
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
        core::CalculationWork::energy_evaluation();
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

// EqIDs packing_fraction and neutral_segment_diameter: the packing coefficient of composition x at T with its
// per-component terms and T derivative. No admission: callers own the packing domain.
template <typename Base>
PackingTermsValues<Base> packing_terms(const Parameters& model_parameters, Base temperature,
                                       const std::vector<Base>& composition) {
    const auto& parameters = [&]() -> decltype(auto) {
        if constexpr (std::is_same_v<Base, double>) return model_parameters;
        else return make_tape_parameters<Base>(model_parameters, {}, {});
    }();
    if (parameters.ideal_gas) {
        PackingTermsValues<Base> result;
        result.coefficient = 1.0; // q is log inverse unit-amount volume; composition and T partials vanish.
        result.component.assign(composition.size(), 0.0);
        return result;
    }
    const auto sigma = segment_sigma(temperature, parameters);
    const auto diameter = segment_diameter(temperature, parameters, sigma);
    const auto diameter_t = segment_diameter_temperature_derivative(temperature, parameters, sigma);
    PackingTermsValues<Base> result;
    result.component.resize(composition.size());
    const Base scale = (Base(PI) / Base(6.0)) * Base(NUMBER_DENSITY_PER_ANGSTROM3);
    for (std::size_t i = 0; i < composition.size(); ++i) {
        result.component[i] = scale * composition[i] * parameters.segment_count[i]
            * diameter[i] * diameter[i] * diameter[i];
        result.coefficient += result.component[i];
        result.temperature_derivative += scale * composition[i]
            * parameters.segment_count[i] * 3.0 * diameter[i] * diameter[i] * diameter_t[i];
    }
    return result;
}

}  // namespace epcsaft::detail
