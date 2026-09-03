#include "model/helmholtz_model.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace epcsaft {

std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const std::vector<double>& molar_masses,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality
) {
    // EqID: formula_unit_composition
    if (solvent_fractions.size() != molar_masses.size()
        || ion_indices[0] >= solvent_fractions.size()
        || ion_indices[1] >= solvent_fractions.size()
        || ion_indices[0] == ion_indices[1]
        || stoichiometry[0] <= 0 || stoichiometry[1] <= 0
        || !std::isfinite(molality) || molality <= 0.0) {
        throw std::invalid_argument("invalid reference composition inputs");
    }
    double fraction_sum = 0.0;
    double solvent_mass = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        if (!std::isfinite(solvent_fractions[component])
            || solvent_fractions[component] < 0.0
            || !std::isfinite(molar_masses[component])
            || (solvent_fractions[component] > 0.0 && molar_masses[component] <= 0.0)) {
            throw std::invalid_argument("invalid reference solvent fractions");
        }
        fraction_sum += solvent_fractions[component];
        solvent_mass += solvent_fractions[component] * molar_masses[component];
    }
    if (std::abs(fraction_sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value
        || solvent_fractions[ion_indices[0]] != 0.0
        || solvent_fractions[ion_indices[1]] != 0.0
        || !(solvent_mass > 0.0) || !std::isfinite(solvent_mass)) {
        throw std::invalid_argument("invalid reference solvent molar mass");
    }
    const double solvent_moles = 1.0 / solvent_mass;
    const double total_moles = solvent_moles
        + molality * static_cast<double>(stoichiometry[0] + stoichiometry[1]);
    std::vector<double> composition(solvent_fractions.size(), 0.0);
    for (std::size_t component = 0; component < composition.size(); ++component) {
        composition[component] = solvent_fractions[component] * solvent_moles / total_moles;
    }
    composition[ion_indices[0]] = stoichiometry[0] * molality / total_moles;
    composition[ion_indices[1]] = stoichiometry[1] * molality / total_moles;
    return composition;
}

InfiniteDilutionReference infinite_dilution_reference(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry
) {
    core::Scope scope("infinite_dilution_reference");
    std::array<double, 4> molalities = {
        policy::REFERENCE_MOLALITY_ONE.value,
        policy::REFERENCE_MOLALITY_TWO.value,
        policy::REFERENCE_MOLALITY_THREE.value,
        policy::REFERENCE_MOLALITY_FOUR.value,
    };
    std::optional<double> anchor;
    std::array<double, 2> previous{};
    std::array<double, 2> current{};
    double convergence_error = 0.0;
    double terminal_error = 0.0;
    bool observed = false;
    for (std::size_t index = 0; index < molalities.size(); ++index) {
        const double molality = molalities[index];
        const auto composition = formula_composition(
            solvent_fractions,
            model.parameters().molar_mass_kg_per_mol,
            ion_indices,
            stoichiometry,
            molality
        );
        const auto state = model.state_at_pressure(
            temperature, pressure, composition, DensityBranch::Liquid, anchor
        );
        if (!state.log_fugacity_coefficient.has_value()) {
            throw std::runtime_error("reference fugacity is unavailable");
        }
        current = {
            (*state.log_fugacity_coefficient)[ion_indices[0]],
            (*state.log_fugacity_coefficient)[ion_indices[1]],
        };
        scope.probe(
            "reference_ladder_point", {molality, current[0], current[1]},
            static_cast<std::uint32_t>(index)
        );
        if (observed) {
            // EqID: reference_limit_terminal_pair
            convergence_error = std::max(
                std::abs(current[0] - previous[0]),
                std::abs(current[1] - previous[1])
            );
            terminal_error = convergence_error;
        }
        previous = current;
        observed = true;
        anchor = state.molar_density;
        scope.probe("reference_convergence_error", convergence_error,
                    static_cast<std::uint32_t>(index));
    }
    if (convergence_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value) {
        const auto composition = formula_composition(
            solvent_fractions,
            model.parameters().molar_mass_kg_per_mol,
            ion_indices,
            stoichiometry,
            policy::REFERENCE_CONVERGENCE_WITNESS.value
        );
        const auto state = model.state_at_pressure(
            temperature, pressure, composition, DensityBranch::Liquid, anchor
        );
        if (!state.log_fugacity_coefficient.has_value()) {
            throw std::runtime_error("reference fugacity is unavailable");
        }
        const double witness_error = std::max(
            std::abs((*state.log_fugacity_coefficient)[ion_indices[0]] - current[0]),
            std::abs((*state.log_fugacity_coefficient)[ion_indices[1]] - current[1])
        );
        scope.probe("reference_witness_error", witness_error);
        if (witness_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value) {
            throw std::runtime_error("reference limit did not converge");
        }
    }
    std::optional<double> pure_solvent;
    std::size_t solvent_index = solvent_fractions.size();
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        if (solvent_fractions[component] == 0.0) continue;
        if (solvent_index != solvent_fractions.size()) {
            solvent_index = solvent_fractions.size();
            break;
        }
        solvent_index = component;
    }
    if (solvent_index < solvent_fractions.size()) {
        std::vector<double> pure(solvent_fractions.size(), 0.0);
        pure[solvent_index] = 1.0;
        const auto state = model.state_at_pressure(
            temperature, pressure, pure, DensityBranch::Liquid, std::nullopt
        );
        if (state.log_fugacity_coefficient.has_value()) {
            pure_solvent = (*state.log_fugacity_coefficient)[solvent_index];
        }
    }
    return {current, molalities.back(), terminal_error, pure_solvent};
}

ElectrolyteActivity electrolyte_activity(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality,
    const InfiniteDilutionReference& reference
) {
    const auto composition = formula_composition(
        solvent_fractions,
        model.parameters().molar_mass_kg_per_mol,
        ion_indices,
        stoichiometry,
        molality
    );
    const auto state = model.state_at_pressure(
        temperature, pressure, composition, DensityBranch::Liquid, std::nullopt
    );
    if (!state.log_fugacity_coefficient.has_value()) {
        throw std::runtime_error("reference fugacity is unavailable");
    }
    // EqID: single_ion_activity_coefficient
    const double first = (*state.log_fugacity_coefficient)[ion_indices[0]]
        - reference.log_fugacity[0];
    const double second = (*state.log_fugacity_coefficient)[ion_indices[1]]
        - reference.log_fugacity[1];
    const double coefficient_sum = stoichiometry[0] + stoichiometry[1];
    // EqID: mean_ionic_activity_coefficient
    const double mean_log_x = (
        stoichiometry[0] * first + stoichiometry[1] * second
    ) / coefficient_sum;
    double solvent_mass = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        solvent_mass += solvent_fractions[component]
            * model.parameters().molar_mass_kg_per_mol[component];
    }
    // EqID: molality_activity_conversion
    const double mean_log_m = mean_log_x
        - std::log(1.0 + solvent_mass * molality * coefficient_sum);
    ElectrolyteActivity result;
    result.single_ion = {std::exp(first), std::exp(second)};
    result.mean_ionic_mole_fraction_scale = std::exp(mean_log_x);
    result.mean_ionic_log_molality_scale = mean_log_m;
    result.mean_ionic_molality_scale = std::exp(mean_log_m);
    if (reference.pure_solvent_log_fugacity.has_value()) {
        std::size_t solvent_index = solvent_fractions.size();
        for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
            if (solvent_fractions[component] == 0.0) continue;
            if (solvent_index != solvent_fractions.size()) {
                solvent_index = solvent_fractions.size();
                break;
            }
            solvent_index = component;
        }
        if (solvent_index < solvent_fractions.size()) {
            const double log_activity = std::log(composition[solvent_index])
                + (*state.log_fugacity_coefficient)[solvent_index]
                - *reference.pure_solvent_log_fugacity;
            // EqID: practical_osmotic_coefficient
            result.osmotic_coefficient = -log_activity
                / (coefficient_sum * molality * solvent_mass);
        }
    }
    return result;
}

}  // namespace epcsaft
