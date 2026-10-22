#include "properties/thermochemistry.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace epcsaft {
namespace {

bool finite_values(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}
bool nonnegative_values(const std::vector<double>& values) {
    return finite_values(values) && std::all_of(values.begin(), values.end(), [](double value) { return value >= 0.0; });
}
void validate_identity(const std::string& identity) {
    if (identity.empty()) throw std::invalid_argument("missing standard-state identity");
}
void validate_composition(const std::vector<double>& composition) {
    if (!nonnegative_values(composition)) throw std::invalid_argument("invalid thermochemistry composition");
    if (std::abs(std::accumulate(composition.begin(), composition.end(), 0.0) - 1.0)
        > policy::COMPOSITION_SUM_TOLERANCE.value)
        throw std::invalid_argument("invalid thermochemistry composition sum");
}
void validate_correlation(const ReactionLogPolynomial& correlation) {
    const bool finite = std::isfinite(correlation.a) && std::isfinite(correlation.b)
        && std::isfinite(correlation.c) && std::isfinite(correlation.d)
        && std::isfinite(correlation.e) && std::isfinite(correlation.f)
        && std::isfinite(correlation.reference_temperature)
        && std::isfinite(correlation.temperature_min) && std::isfinite(correlation.temperature_max);
    if (!finite || correlation.reference_temperature <= 0.0 || correlation.temperature_min <= 0.0
        || correlation.temperature_max < correlation.temperature_min)
        throw std::invalid_argument("invalid reaction correlation");
    validate_identity(correlation.standard_state_id);
}
void validate_correlation_temperature(const ReactionLogPolynomial& correlation, double temperature) {
    validate_correlation(correlation);
    if (!std::isfinite(temperature) || temperature < correlation.temperature_min || temperature > correlation.temperature_max)
        throw std::invalid_argument("temperature outside reaction correlation");
}
ReactionProperties make_reaction_properties(
    double ln_k, double derivative, double enthalpy, const std::string& identity,
    double second_derivative, double third_derivative
) {
    validate_identity(identity);
    const double k = std::exp(ln_k);
    if (!std::isfinite(ln_k) || !std::isfinite(derivative) || !std::isfinite(enthalpy)
        || !std::isfinite(second_derivative) || !std::isfinite(third_derivative))
        throw std::domain_error("nonfinite reaction properties");
    return {ln_k, derivative, enthalpy, second_derivative, third_derivative,
            std::isfinite(k) && k > 0.0 ? std::optional<double>(k) : std::nullopt, identity};
}

}  // namespace

void validate_thermochemistry(const ThermochemistryRecord& record, std::size_t component_count) {
    if (!std::isfinite(record.standard_pressure) || record.standard_pressure <= 0.0)
        throw std::invalid_argument("invalid thermochemistry standard pressure");
    if (record.component_intervals.size() != component_count)
        throw std::invalid_argument("thermochemistry component order");
    // An empty list declares a species without an ideal-gas record (reaction-completed calorics).
    for (const auto& intervals : record.component_intervals)
        if (!intervals.empty()) validate_ideal_intervals(intervals);
}

IdealPropertiesEvaluation evaluate_ideal_properties(
    const ThermochemistryRecord& record, double temperature, double molar_density,
    const std::vector<double>& mole_fractions
) {
    if (record.component_intervals.size() != mole_fractions.size())
        throw std::invalid_argument("thermochemistry component order");
    if (!std::isfinite(temperature) || temperature <= 0.0 || !std::isfinite(molar_density) || molar_density <= 0.0)
        throw std::invalid_argument("invalid ideal state");
    validate_composition(mole_fractions);
    if (std::any_of(record.component_intervals.begin(), record.component_intervals.end(),
                    [](const auto& intervals) { return intervals.empty(); }))
        return {std::nullopt, ThermochemistryStatus::MissingRecord};
    if (std::any_of(record.component_intervals.begin(), record.component_intervals.end(),
                    [temperature](const auto& intervals) { return ideal_interval(intervals, temperature) == nullptr; }))
        return {std::nullopt, ThermochemistryStatus::OutsideInterval};

    try {
        std::vector<IdealCorrelation> correlations;
        std::vector<bool> supported;
        correlations.reserve(mole_fractions.size());
        supported.reserve(mole_fractions.size());
        for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
            correlations.push_back(ideal_interval(record.component_intervals[component], temperature)->correlation);
            supported.push_back(mole_fractions[component] > 0.0);
        }
        const auto value = evaluate_ideal_mixture(
            correlations, record.standard_pressure, temperature, molar_density,
            mole_fractions, supported);
        IdealProperties result{
            value.enthalpy, value.internal_energy, value.entropy,
            value.isobaric_heat_capacity, value.isochoric_heat_capacity};
        if (!std::isfinite(result.enthalpy) || !std::isfinite(result.internal_energy)
            || !std::isfinite(result.entropy) || !std::isfinite(result.isobaric_heat_capacity)
            || !std::isfinite(result.isochoric_heat_capacity))
            return {std::nullopt, ThermochemistryStatus::Nonfinite};
        return {result, ThermochemistryStatus::Available};
    } catch (const std::domain_error&) { return {std::nullopt, ThermochemistryStatus::Nonfinite}; }
}

ReactionProperties evaluate_reaction_properties(const ReactionLogPolynomial& correlation, double temperature) {
    validate_correlation_temperature(correlation, temperature);
    const double derivative = reaction_log_value(correlation, temperature, 1);
    const double enthalpy = GAS_CONSTANT_J_PER_MOL_K * temperature * temperature * derivative;
    return make_reaction_properties(reaction_log_value(correlation, temperature), derivative, enthalpy,
        correlation.standard_state_id, reaction_log_value(correlation, temperature, 2),
        reaction_log_value(correlation, temperature, 3));
}

ReactionThermochemistry reaction_thermochemistry_density(
    const ThermochemistryRecord& record, const std::vector<double>& stoichiometry, double temperature
) {
    core::Scope scope("reaction_thermochemistry_density");
    validate_thermochemistry(record, stoichiometry.size());
    if (!std::isfinite(temperature) || temperature <= 0.0 || !finite_values(stoichiometry))
        throw std::invalid_argument("invalid reaction thermochemistry request");
    ReactionThermochemistry result;
    Taylor<3> T(temperature), gp(0.0);
    T.c[1] = 1.0;
    double delta_nu = 0.0, hp = 0.0, cp = 0.0;
    for (std::size_t i = 0; i < stoichiometry.size(); ++i) {
        if (stoichiometry[i] == 0.0) continue;
        const auto* interval = ideal_interval(record.component_intervals[i], temperature);
        if (interval == nullptr) {
            result.status = record.component_intervals[i].empty() ? ThermochemistryStatus::MissingRecord : ThermochemistryStatus::OutsideInterval;
            result.cause = std::string(record.component_intervals[i].empty() ? "MissingRecord" : "OutsideInterval")
                + ": component " + std::to_string(i) + " at T=" + std::to_string(temperature) + " K";
            return result;
        }
        const auto thermal = evaluate_ideal_correlation(interval->correlation, T, record.standard_pressure);
        if (!std::isfinite(thermal.h.c[0]) || !std::isfinite(thermal.s.c[0]) || !std::isfinite(thermal.cp.c[0])) {
            result.status = ThermochemistryStatus::Nonfinite;
            result.cause = "Nonfinite: component " + std::to_string(i) + " at T=" + std::to_string(temperature) + " K";
            return result;
        }
        gp += stoichiometry[i] * (thermal.h - T * thermal.s);
        hp += stoichiometry[i] * thermal.h.c[0];
        cp += stoichiometry[i] * thermal.cp.c[0];
        delta_nu += stoichiometry[i];
    }
    const auto rt = GAS_CONSTANT_J_PER_MOL_K * T;
    // EqID: standard_pressure_density_conversion, rho0 = 1 mol/m3.
    const auto g = gp + delta_nu * rt * scalar_log(rt / record.standard_pressure);
    result.gibbs = g.c[0]; result.enthalpy = hp - delta_nu * rt.c[0];
    result.heat_capacity = cp - delta_nu * GAS_CONSTANT_J_PER_MOL_K;
    result.pressure_gibbs = gp.c[0]; result.ln_k_pressure = -gp.c[0] / rt.c[0]; result.ln_k_density = -g.c[0] / rt.c[0];
    result.gibbs_derivatives = {g.c[0], g.c[1], 2 * g.c[2], 6 * g.c[3]};
    if (!std::isfinite(result.enthalpy) || !std::isfinite(result.heat_capacity) || !std::isfinite(result.ln_k_pressure)
        || !std::isfinite(result.ln_k_density) || !std::all_of(result.gibbs_derivatives.begin(), result.gibbs_derivatives.end(),
            [](double v) { return std::isfinite(v); })) {
        result.status = ThermochemistryStatus::Nonfinite; result.cause = "Nonfinite reaction thermochemistry at T=" + std::to_string(temperature) + " K";
        return result;
    }
    for (auto item : {std::pair{result.ln_k_pressure, &result.k_pressure}, std::pair{result.ln_k_density, &result.k_density}}) {
        const double k = std::exp(item.first); if (std::isfinite(k) && k > 0.0) *item.second = k;
    }
    result.status = ThermochemistryStatus::Available;
    return result;
}

}  // namespace epcsaft
