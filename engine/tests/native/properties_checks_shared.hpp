#pragma once
// Includes and helpers shared by the properties_checks translation units.
#include "properties/helmholtz_model.hpp"
#include "properties_support.hpp"
#undef CHECK
#include "recording_support.hpp"
#include "contributions/association.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#define CHECK(condition) recording_support::check( \
    static_cast<bool>(condition), __FILE__, __LINE__, #condition)

namespace properties_checks {

using properties_support::p3_forms;
using properties_support::p3_record;
using properties_support::p3_binary_record;

inline bool close(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

inline bool finite_state_without_cp(const epcsaft::State& state) {
    return std::isfinite(state.residual_helmholtz) && std::isfinite(state.pressure)
        && std::isfinite(state.residual_entropy) && std::isfinite(state.residual_internal_energy)
        && std::isfinite(state.residual_enthalpy)
        && std::isfinite(state.residual_isochoric_heat_capacity)
        && std::isfinite(state.molar_density)
        && std::isfinite(state.pressure_density_derivative);
}

inline std::string model_inputs(const epcsaft::Parameters& parameters, double temperature,
                         double volume, const std::vector<double>& amounts) {
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    const std::vector<double> composition = total == 0.0
        ? std::vector<double>{} : [&] {
            std::vector<double> result;
            result.reserve(amounts.size());
            for (const double amount : amounts) result.push_back(amount / total);
            return result;
        }();
    return recording_support::fields({
        recording_support::field("T", temperature),
        recording_support::field("V", volume),
        recording_support::field("rho", total / volume),
        recording_support::field("n", amounts),
        recording_support::field("x", composition),
        recording_support::field("segments", parameters.segment_count),
        recording_support::field("sigma_A", parameters.sigma_constant),
        recording_support::field("epsilon_over_k_K", parameters.dispersion_energy_over_k),
    });
}

inline std::vector<double> p2_association_residuals(
    double temperature, double volume, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites
) {
    const double density = 1.0 / volume;
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const auto delta = epcsaft::association_delta(
        temperature, density, composition, parameters, sigma, diameter
    );
    return epcsaft::association_residuals(
        density, composition, parameters, sites, delta
    );
}

struct P2AssociationConditioning final {
    Eigen::MatrixXd jacobian;
    std::size_t rank = 0;
    double pivot_ratio = 1.0;
};

inline P2AssociationConditioning p2_association_conditioning(
    double temperature, double density, const std::vector<double>& composition,
    const epcsaft::Parameters& parameters, const std::vector<double>& sites
) {
    const std::size_t count = sites.size();
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const auto delta = epcsaft::association_delta(
        temperature, density, composition, parameters, sigma, diameter
    );
    const double number_density = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(count), static_cast<Eigen::Index>(count)
    );
    for (std::size_t a = 0; a < count; ++a) {
        double bonded = 0.0;
        for (std::size_t b = 0; b < count; ++b) {
            const auto component = static_cast<std::size_t>(parameters.site_component[b]);
            const double mass = composition[component] * parameters.site_multiplicity[b];
            bonded += mass * sites[b] * delta[a * count + b];
        }
        for (std::size_t b = 0; b < count; ++b) {
            const auto component = static_cast<std::size_t>(parameters.site_component[b]);
            const double mass = composition[component] * parameters.site_multiplicity[b];
            jacobian(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
                (a == b ? 1.0 + number_density * bonded : 0.0)
                + sites[a] * number_density * mass * delta[a * count + b];
        }
    }
    P2AssociationConditioning result{std::move(jacobian)};
    if (count == 0) return result;
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(result.jacobian);
    factorization.setThreshold(epcsaft::policy::ASSOCIATION_JACOBIAN_RANK.value * count);
    result.rank = static_cast<std::size_t>(factorization.rank());
    const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
    result.pivot_ratio = pivots.minCoeff() / pivots.maxCoeff();
    return result;
}

inline double p4_state_value(const epcsaft::State& state, epcsaft::PropertyObservable observable) {
    const auto index = static_cast<std::size_t>(observable) % 5;
    const bool ideal = static_cast<std::size_t>(observable) < 5;
    const bool total = static_cast<std::size_t>(observable) >= 10;
    const auto residual = std::array<double, 5>{
        state.residual_enthalpy, state.residual_internal_energy, state.residual_entropy,
        state.residual_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.residual_isochoric_heat_capacity};
    const auto ideal_values = std::array<double, 5>{
        state.ideal_enthalpy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_internal_energy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_entropy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.ideal_isochoric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN())};
    const auto total_values = std::array<double, 5>{
        state.total_enthalpy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_internal_energy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_entropy.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()),
        state.total_isochoric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN())};
    return (ideal ? ideal_values : total ? total_values : residual)[index];
}

inline epcsaft::Parameters p4_shifted_parameters(
    const epcsaft::Parameters& base,
    const std::vector<epcsaft::ActiveParameter>& active,
    const std::vector<double>& shift
) {
    auto result = base;
    for (std::size_t coordinate = 0; coordinate < active.size(); ++coordinate) {
        const auto& item = active[coordinate];
        auto apply = [&](std::vector<double>& values) {
            for (const auto index : item.indices) values[index] += shift[coordinate];
        };
        switch (item.family) {
        case epcsaft::ParameterFamily::AssociationEnergyOverK:
            apply(result.association_energy_over_k); break;
        case epcsaft::ParameterFamily::AssociationVolume:
            apply(result.association_volume); break;
        case epcsaft::ParameterFamily::Kij:
            apply(result.k_ij); break;
        case epcsaft::ParameterFamily::BornDiameter:
            apply(result.born_diameter_angstrom); break;
        case epcsaft::ParameterFamily::DebyeHuckelDiameter:
            apply(result.debye_huckel_diameter_angstrom); break;
        case epcsaft::ParameterFamily::PackingDiameter:
            apply(result.packing_diameter_angstrom); break;
        case epcsaft::ParameterFamily::IonSuppressionCoefficient:
            apply(result.ion_suppression_coefficient); break;
        default: throw std::invalid_argument("unsupported P4 reference parameter");
        }
    }
    return result;
}


inline std::vector<double> p4_take_signature(
    epcsaft::core::Capture& capture,
    const std::optional<epcsaft::ThermochemistryRecord>& thermochemistry,
    double temperature
) {
    std::vector<double> result;
    for (const auto& event : capture.trace().events)
        if (std::string(event.name) == "predicate_signature") result = event.values;
    if (thermochemistry.has_value())
        for (const auto& intervals : thermochemistry->component_intervals)
            for (const auto& interval : intervals)
                result.push_back(epcsaft::ideal_interval_contains(interval, temperature));
    capture.take();
    return result;
}

// Differentiate the preceding action once, instead of subtracting O(h^k) caloric values.
// ALG-TERMINATION: cancellations scale with the fixed chart-coordinate terms.
struct PropertyActionReference final {
    std::vector<double> coarse, middle, fine, value, uncertainty, scale;
};

template <class Direction>
PropertyActionReference property_action_reference(
    const epcsaft::Parameters& parameters, const std::vector<epcsaft::ActiveParameter>& active,
    const std::optional<epcsaft::ThermochemistryRecord>& record, double T, double mechanical,
    const std::vector<double>& x, const std::vector<epcsaft::PropertyObservable>& observables,
    const std::vector<Direction>& directions, epcsaft::DensityBranch branch = epcsaft::DensityBranch::Unspecified
) {
    using namespace epcsaft;
    const HelmholtzModel model(parameters, active, record);
    const auto batch = [&](const HelmholtzModel& m, double t, double q, const std::vector<double>& y,
                           const std::vector<Direction>& request) {
        if constexpr (std::is_same_v<Direction, DensityActionDirection>)
            return m.density_actions(t, q, y, observables, request);
        else return m.pressure_actions(t, q, y, branch, observables, request);
    };
    constexpr bool pressure_chart = !std::is_same_v<Direction, DensityActionDirection>;
    // One row per result (ln(phi_i) has one per component), then the pressure chart's density action.
    const auto rows = [&](const PropertyActionBatch& actions, std::size_t order) {
        std::vector<double> out;
        for (const auto& item : actions.results) {
            CHECK(item.status == PropertyActionStatus::Available && item.action);
            out.push_back(*item.action);
        }
        if constexpr (pressure_chart) {
            CHECK(actions.density_status == PropertyActionStatus::Available
                && (order == 1 ? actions.density_action.size() == 1 : actions.density_mixed_action.has_value()));
            out.push_back(order == 1 ? actions.density_action[0] : *actions.density_mixed_action);
        }
        return out;
    };
    PropertyActionReference result;
    // Omit only zero sub-directions; retain the entire multilinear coordinate split.
    std::vector<Direction> terms;
    const auto expand = [&](const auto& self, std::size_t j) -> void {
        if (j == directions.size()) {
            const auto actions = rows(batch(model, T, mechanical, x, terms), terms.size());
            result.scale.resize(actions.size(), 0.0);
            for (std::size_t i = 0; i < actions.size(); ++i) result.scale[i] += std::abs(actions[i]);
            return;
        }
        for (int part = 0; part < 4; ++part) {
            Direction d{};
            d.composition.assign(x.size(), 0.0);
            if (part == 0) d.temperature = directions[j].temperature;
            if (part == 1) {
                if constexpr (std::is_same_v<Direction, DensityActionDirection>) d.molar_density = directions[j].molar_density;
                else d.pressure = directions[j].pressure;
            }
            if (part == 2) d.composition = directions[j].composition;
            if (part == 3) d.parameter = directions[j].parameter;
            bool nonzero = d.temperature != 0.0 || !d.parameter.empty();
            if constexpr (std::is_same_v<Direction, DensityActionDirection>) nonzero = nonzero || d.molar_density != 0.0;
            else nonzero = nonzero || d.pressure != 0.0;
            for (const auto v : d.composition) nonzero = nonzero || v != 0.0;
            if (nonzero) { terms.push_back(d); self(self, j + 1); terms.pop_back(); }
        }
    };
    expand(expand, 0);
    epcsaft::core::Capture capture;
    const auto central_state = [&] {
        if constexpr (std::is_same_v<Direction, DensityActionDirection>) return model.state(T, mechanical, x);
        else return model.state_at_pressure(T, mechanical, x, branch);
    }();
    const auto central_signature = p4_take_signature(capture, record, T);
    const auto& last = directions.back();
    const std::vector<Direction> lower(directions.begin(), directions.end() - 1);
    const auto values = [&](double h) {
        std::vector<double> shift(active.size(), 0.0), y = x;
        for (std::size_t i = 0; i < x.size(); ++i) y[i] += h * last.composition[i];
        for (const auto& parameter : last.parameter)
            for (std::size_t a = 0; a < active.size(); ++a)
                if (parameter.family == active[a].family && parameter.indices == active[a].indices) shift[a] += h * parameter.value;
        const HelmholtzModel shifted(p4_shifted_parameters(parameters, active, shift), active, record);
        const double t = T + h * last.temperature;
        double q = mechanical;
        if constexpr (std::is_same_v<Direction, DensityActionDirection>) q += h * last.molar_density;
        else q += h * last.pressure;
        const auto state = [&] {
            if constexpr (std::is_same_v<Direction, DensityActionDirection>) return shifted.state(t, q, y);
            else return shifted.state_at_pressure(t, q, y, branch);
        }();
        CHECK(p4_take_signature(capture, record, t) == central_signature);
        for (std::size_t i = 0; i < x.size(); ++i) CHECK((x[i] > 0.0) == (y[i] > 0.0));
        if constexpr (!std::is_same_v<Direction, DensityActionDirection>) {
            CHECK(state.density_branch == central_state.density_branch && state.stable_root_count == central_state.stable_root_count);
            CHECK(std::isfinite(state.pressure_density_derivative) && state.pressure_density_derivative > 0.0);
            CHECK(std::abs(state.pressure - q) / (state.molar_density * GAS_CONSTANT_J_PER_MOL_K * t)
                <= policy::DENSITY_ROOT_RESIDUAL.value);
            const auto residuals = p2_association_residuals(t, 1.0 / state.molar_density, y,
                shifted.parameters(), state.site_fractions);
            for (const auto residual : residuals) CHECK(std::isfinite(residual) && std::abs(residual) <= 2.0e-13);
            const auto conditioning = p2_association_conditioning(t, state.molar_density, y,
                shifted.parameters(), state.site_fractions);
            CHECK(conditioning.rank == state.site_fractions.size()
                && conditioning.pivot_ratio > policy::ASSOCIATION_JACOBIAN_CONDITION.value);
        }
        if (!lower.empty()) return rows(batch(shifted, t, q, y, lower), lower.size());
        std::vector<double> values;
        for (const auto obs : observables) {
            if (obs != PropertyObservable::LogFugacityCoefficient) values.push_back(p4_state_value(state, obs));
            else for (const auto& value : *state.log_fugacity_coefficient) values.push_back(value.value());
        }
        if constexpr (pressure_chart) values.push_back(state.molar_density);
        return values;
    };
    // Higher-action subtraction at 0.0025 is density-root rounding dominated;
    // these nodes resolve cp/cv without crossing support or thermal pieces.
    constexpr double steps[] = {0.1, 0.05, 0.025};
    std::vector<std::vector<double>> rungs;
    std::vector<double> rounding;
    for (const double h : steps) {
        const auto plus = values(h), minus = values(-h);
        std::vector<double> stencil;
        rounding.clear();
        for (std::size_t i = 0; i < plus.size(); ++i) {
            stencil.push_back((plus[i] - minus[i]) / (2.0 * h));
            rounding.push_back(32.0 * std::numeric_limits<double>::epsilon() * (std::abs(plus[i]) + std::abs(minus[i])) / (2.0 * h));
        }
        rungs.push_back(stencil);
    }
    result.coarse = rungs[0]; result.middle = rungs[1]; result.fine = rungs[2];
    for (std::size_t i = 0; i < rungs[2].size(); ++i) {
        const double refined = (4.0 * rungs[2][i] - rungs[1][i]) / 3.0;
        const double coarser = (4.0 * rungs[1][i] - rungs[0][i]) / 3.0;
        result.value.push_back(refined);
        // Polynomial ideal cp/cv has identically zero higher chart terms.
        // An exactly constant preceding action subtracts exactly: its nonzero
        // value must not manufacture a rounding floor for a zero stencil.
        const bool exact_zero = result.scale[i] == 0.0 && rungs[0][i] == 0.0
            && rungs[1][i] == 0.0 && rungs[2][i] == 0.0;
        result.uncertainty.push_back(exact_zero ? 0.0 : std::abs(refined - coarser) + rounding[i]);
    }
    return result;
}

}  // namespace properties_checks
