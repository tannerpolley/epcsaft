#include "properties/thermochemistry.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

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
void require_sizes(const std::vector<double>& first, const std::vector<double>& second, const char* field) {
    if (first.size() != second.size()) throw std::invalid_argument(std::string("thermochemistry ") + field + " size");
}
void validate_positive_values(const std::vector<double>& values, const char* field) {
    if (!std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value) && value > 0.0; }))
        throw std::invalid_argument(std::string("invalid ") + field);
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
        && std::isfinite(correlation.reference_temperature)
        && std::isfinite(correlation.temperature_min) && std::isfinite(correlation.temperature_max);
    if (!finite || correlation.reference_temperature <= 0.0 || correlation.temperature_min <= 0.0
        || correlation.temperature_max < correlation.temperature_min
        || correlation.reference_temperature < correlation.temperature_min
        || correlation.reference_temperature > correlation.temperature_max)
        throw std::invalid_argument("invalid reaction correlation");
    validate_identity(correlation.standard_state_id);
}
void validate_correlation_temperature(const ReactionLogPolynomial& correlation, double temperature) {
    validate_correlation(correlation);
    if (!std::isfinite(temperature) || temperature < correlation.temperature_min || temperature > correlation.temperature_max)
        throw std::invalid_argument("temperature outside reaction correlation");
}
void validate_reaction_properties(const ReactionProperties& value) {
    if (!std::isfinite(value.ln_k) || !std::isfinite(value.derivative)
        || !std::isfinite(value.enthalpy)) throw std::invalid_argument("invalid reaction properties");
    if (value.k.has_value() && (!std::isfinite(*value.k) || *value.k <= 0.0))
        throw std::invalid_argument("invalid reaction properties");
    validate_identity(value.standard_state_id);
}
bool contains_temperature(const std::vector<IdealInterval>& intervals, double temperature) {
    return std::any_of(intervals.begin(), intervals.end(), [temperature](const IdealInterval& interval) {
        return ideal_interval_contains(interval, temperature);
    });
}
ReactionProperties make_reaction_properties(
    double ln_k, double derivative, double enthalpy, const std::string& identity
) {
    validate_identity(identity);
    const double k = std::exp(ln_k);
    if (!std::isfinite(ln_k) || !std::isfinite(derivative) || !std::isfinite(enthalpy))
        throw std::domain_error("nonfinite reaction properties");
    return {ln_k, derivative, enthalpy,
            std::isfinite(k) && k > 0.0 ? std::optional<double>(k) : std::nullopt, identity};
}

}  // namespace

void validate_thermochemistry(const ThermochemistryRecord& record, std::size_t component_count) {
    if (!std::isfinite(record.standard_pressure) || record.standard_pressure <= 0.0)
        throw std::invalid_argument("invalid thermochemistry standard pressure");
    if (record.component_intervals.size() != component_count)
        throw std::invalid_argument("thermochemistry component order");
    for (const auto& intervals : record.component_intervals) validate_ideal_intervals(intervals);
}

IdealPropertiesEvaluation evaluate_ideal_properties(
    const ThermochemistryRecord& record, double temperature, double molar_density,
    const std::vector<double>& mole_fractions
) {
    validate_thermochemistry(record, mole_fractions.size());
    if (!std::isfinite(temperature) || temperature <= 0.0 || !std::isfinite(molar_density) || molar_density <= 0.0)
        throw std::invalid_argument("invalid ideal state");
    validate_composition(mole_fractions);
    if (std::any_of(record.component_intervals.begin(), record.component_intervals.end(),
                    [temperature](const auto& intervals) { return !contains_temperature(intervals, temperature); }))
        return {std::nullopt, ThermochemistryStatus::OutsideInterval};

    try {
        std::vector<IdealCorrelation> correlations;
        std::vector<bool> supported;
        correlations.reserve(mole_fractions.size());
        supported.reserve(mole_fractions.size());
        for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
            correlations.push_back(ideal_interval(
                record.component_intervals[component], temperature).correlation);
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

DensityStandardState pressure_to_density_standard(
    double gibbs, double enthalpy, double entropy, double temperature, double pressure, double density
) {
    if (!std::isfinite(gibbs) || !std::isfinite(enthalpy) || !std::isfinite(entropy)
        || !std::isfinite(temperature) || temperature <= 0.0 || !std::isfinite(pressure) || pressure <= 0.0
        || !std::isfinite(density) || density <= 0.0) throw std::invalid_argument("invalid standard-state conversion");
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    const double shift = rt * std::log(rt * density / pressure);
    return {gibbs + shift, enthalpy - rt, entropy - GAS_CONSTANT_J_PER_MOL_K - shift / temperature};
}

ActivitySet direct_activities(const std::vector<double>& values, const std::string& identity) {
    validate_positive_values(values, "activities");
    validate_identity(identity);
    return {values, identity};
}
ActivitySet mole_fraction_activities(
    const std::vector<double>& coefficients, const std::vector<double>& mole_fractions,
    const std::string& identity
) {
    require_sizes(coefficients, mole_fractions, "mole fraction");
    validate_positive_values(coefficients, "activity coefficients");
    validate_composition(mole_fractions);
    validate_identity(identity);
    ActivitySet result{{}, identity};
    result.values.resize(coefficients.size());
    for (std::size_t i = 0; i < result.values.size(); ++i) result.values[i] = coefficients[i] * mole_fractions[i];
    return result;
}
ActivitySet molarity_activities(
    const std::vector<double>& coefficients, double molar_density, const std::vector<double>& mole_fractions,
    double standard, const std::string& identity
) {
    require_sizes(coefficients, mole_fractions, "molarity");
    validate_positive_values(coefficients, "activity coefficients");
    validate_identity(identity);
    validate_composition(mole_fractions);
    if (!std::isfinite(molar_density) || molar_density <= 0.0 || !std::isfinite(standard) || standard <= 0.0)
        throw std::invalid_argument("invalid molarity inputs");
    ActivitySet result{{}, identity};
    result.values.resize(coefficients.size());
    for (std::size_t i = 0; i < result.values.size(); ++i)
        result.values[i] = coefficients[i] * molar_density * mole_fractions[i] / standard;
    return result;
}
ActivitySet molality_activities(
    const std::vector<double>& coefficients, const std::vector<double>& amounts_mol, double solvent_mass_kg,
    double standard, const std::string& identity
) {
    require_sizes(coefficients, amounts_mol, "molality");
    validate_positive_values(coefficients, "activity coefficients");
    validate_identity(identity);
    if (!std::isfinite(solvent_mass_kg) || solvent_mass_kg <= 0.0
        || !std::isfinite(standard) || standard <= 0.0 || !nonnegative_values(amounts_mol))
        throw std::invalid_argument("invalid molality inputs");
    ActivitySet result{{}, identity};
    result.values.resize(coefficients.size());
    for (std::size_t i = 0; i < result.values.size(); ++i)
        result.values[i] = coefficients[i] * amounts_mol[i] / (solvent_mass_kg * standard);
    return result;
}

double reaction_ln_quotient(
    const std::vector<double>& stoichiometry, const ActivitySet& activities,
    const std::vector<int>& charges, const ReactionLogPolynomial& correlation
) {
    require_sizes(stoichiometry, activities.values, "reaction");
    validate_identity(activities.standard_state_id);
    validate_correlation(correlation);
    if (activities.standard_state_id != correlation.standard_state_id)
        throw std::invalid_argument("reference-basis-mismatch");
    if (!charges.empty() && charges.size() != stoichiometry.size()) throw std::invalid_argument("reaction charge size");
    if (!finite_values(stoichiometry)) throw std::invalid_argument("invalid stoichiometry");
    double charge = 0.0;
    bool ionic = false;
    for (std::size_t i = 0; i < stoichiometry.size(); ++i) {
        if (charges.empty()) continue;
        ionic = ionic || (stoichiometry[i] != 0.0 && charges[i] != 0);
        charge += stoichiometry[i] * static_cast<double>(charges[i]);
    }
    if (ionic && charge != 0.0) throw std::invalid_argument("ionic-gauge-missing");
    double result = 0.0;
    for (std::size_t i = 0; i < stoichiometry.size(); ++i) {
        if (stoichiometry[i] == 0.0) continue;
        if (!std::isfinite(activities.values[i]) || activities.values[i] <= 0.0)
            throw std::invalid_argument("reaction activity unavailable");
        result += stoichiometry[i] * std::log(activities.values[i]);
    }
    if (!std::isfinite(result)) throw std::domain_error("nonfinite reaction quotient");
    return result;
}

ReactionProperties evaluate_reaction_properties(const ReactionLogPolynomial& correlation, double temperature) {
    validate_correlation_temperature(correlation, temperature);
    const double derivative = -correlation.b / (temperature * temperature)
        + correlation.c / temperature + correlation.d;
    const double enthalpy = GAS_CONSTANT_J_PER_MOL_K * temperature * temperature * derivative;
    const double ln_k = correlation.a + correlation.b / temperature
        + correlation.c * std::log(temperature / correlation.reference_temperature) + correlation.d * temperature;
    return make_reaction_properties(ln_k, derivative, enthalpy, correlation.standard_state_id);
}

ReactionLogPolynomial normalize_legacy_reaction_correlation(
    double a_old, double b, double c, double d, double reference_temperature,
    double temperature_min, double temperature_max, const std::string& identity
) {
    return {a_old + c * std::log(reference_temperature), b, c, d, reference_temperature,
            temperature_min, temperature_max, identity};
}
ReactionLogPolynomial normalize_log10_reciprocal_correlation(
    double a, double b, double c, double temperature_min, double temperature_max,
    const std::string& identity, double reference_temperature
) {
    return {-std::log(10.0) * b, -std::log(10.0) * a, 0.0, -std::log(10.0) * c,
            reference_temperature, temperature_min, temperature_max, identity};
}
ReactionProperties reaction_pressure_to_density(
    const ReactionProperties& pressure_properties, const std::vector<double>& stoichiometry,
    const std::string& density_identity, double temperature, double pressure, double density
) {
    validate_reaction_properties(pressure_properties);
    if (!std::isfinite(temperature) || temperature <= 0.0 || !std::isfinite(pressure) || pressure <= 0.0
        || !std::isfinite(density) || density <= 0.0 || !finite_values(stoichiometry) || density_identity.empty())
        throw std::invalid_argument("invalid reaction standard-state conversion");
    const double delta_nu = std::accumulate(stoichiometry.begin(), stoichiometry.end(), 0.0);
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    const double shift = std::log(rt * density / pressure);
    return make_reaction_properties(
        pressure_properties.ln_k - delta_nu * shift,
        pressure_properties.derivative - delta_nu / temperature,
        pressure_properties.enthalpy - delta_nu * rt,
        density_identity
    );
}
ReactionProperties reaction_standard_state_conversion(
    const ReactionProperties& properties, const std::vector<double>& stoichiometry,
    const std::vector<double>& old_standard_gibbs_j_per_mol,
    const std::vector<double>& new_standard_gibbs_j_per_mol,
    const std::vector<double>& old_standard_enthalpy_j_per_mol,
    const std::vector<double>& new_standard_enthalpy_j_per_mol,
    const std::string& new_identity, double temperature
) {
    validate_reaction_properties(properties);
    require_sizes(stoichiometry, old_standard_gibbs_j_per_mol, "standard-state");
    require_sizes(stoichiometry, new_standard_gibbs_j_per_mol, "standard-state");
    require_sizes(stoichiometry, old_standard_enthalpy_j_per_mol, "standard-state");
    require_sizes(stoichiometry, new_standard_enthalpy_j_per_mol, "standard-state");
    validate_identity(new_identity);
    if (!std::isfinite(temperature) || temperature <= 0.0 || !finite_values(stoichiometry)
        || !finite_values(old_standard_gibbs_j_per_mol) || !finite_values(new_standard_gibbs_j_per_mol)
        || !finite_values(old_standard_enthalpy_j_per_mol) || !finite_values(new_standard_enthalpy_j_per_mol))
        throw std::invalid_argument("invalid standard-state conversion");
    double gibbs_shift = 0.0, enthalpy_shift = 0.0;
    for (std::size_t i = 0; i < stoichiometry.size(); ++i) {
        gibbs_shift += stoichiometry[i] * (old_standard_gibbs_j_per_mol[i] - new_standard_gibbs_j_per_mol[i]);
        enthalpy_shift += stoichiometry[i] * (old_standard_enthalpy_j_per_mol[i]
                                              - new_standard_enthalpy_j_per_mol[i]);
    }
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    return make_reaction_properties(
        properties.ln_k + gibbs_shift / rt, properties.derivative - enthalpy_shift / (rt * temperature),
        properties.enthalpy - enthalpy_shift, new_identity
    );
}

}  // namespace epcsaft
