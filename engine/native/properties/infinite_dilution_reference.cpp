#include "properties/helmholtz_model.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft {

namespace {

constexpr std::array<double, 5> NEUTRAL_REFERENCE_MOLALITIES = {
    policy::REFERENCE_MOLALITY_ONE.value,
    policy::REFERENCE_MOLALITY_TWO.value,
    policy::REFERENCE_MOLALITY_THREE.value,
    policy::REFERENCE_MOLALITY_FOUR.value,
    policy::REFERENCE_CONVERGENCE_WITNESS.value,
};

double validate_neutral_reference_declaration(
    const Parameters& parameters,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry
) {
    const std::size_t components = parameters.component_count();
    if (!declaration.solvent_index.has_value()
        || *declaration.solvent_index >= components
        || declaration.solute_molality_ratios.size() != components
        || !std::isfinite(declaration.standard_molality_mol_per_kg)
        || declaration.standard_molality_mol_per_kg <= 0.0
        || stoichiometry.size() != components) {
        throw std::invalid_argument("invalid neutral reference declaration");
    }
    const std::size_t solvent = *declaration.solvent_index;
    if (parameters.charge_number[solvent] != 0) {
        throw std::invalid_argument("neutral reference solvent is charged");
    }
    const auto& molar_masses = parameters.molar_mass_kg_per_mol;
    const double solvent_mass = molar_masses[solvent];
    if (!std::isfinite(solvent_mass) || solvent_mass <= 0.0) {
        throw std::invalid_argument("neutral reference solvent has invalid molar mass");
    }
    double charge = 0.0;
    double reaction_charge = 0.0;
    for (std::size_t component = 0; component < components; ++component) {
        const double ratio = declaration.solute_molality_ratios[component];
        if (!std::isfinite(ratio)
            || (component == solvent ? ratio != 0.0 : ratio <= 0.0)) {
            throw std::invalid_argument("invalid neutral reference solute ratio");
        }
        const double coefficient = stoichiometry[component];
        if (!std::isfinite(coefficient)) {
            throw std::invalid_argument("invalid neutral reference stoichiometry");
        }
        charge += ratio * static_cast<double>(parameters.charge_number[component]);
        reaction_charge += coefficient * static_cast<double>(parameters.charge_number[component]);
    }
    if (!std::isfinite(charge) || charge != 0.0) {
        throw std::invalid_argument("neutral reference solute ratios are not charge neutral");
    }
    if (!std::isfinite(reaction_charge) || reaction_charge != 0.0) {
        throw std::invalid_argument("neutral reference reaction is not charge neutral");
    }
    return solvent_mass;
}

std::vector<double> neutral_reference_composition(
    const NeutralReferenceDeclaration& declaration,
    double solvent_mass,
    double molality
) {
    if (!std::isfinite(molality) || molality <= 0.0) {
        throw std::invalid_argument("invalid neutral reference molality");
    }
    const std::size_t solvent = *declaration.solvent_index;
    const double solvent_moles = 1.0 / solvent_mass;
    double total_moles = solvent_moles;
    for (const double ratio : declaration.solute_molality_ratios)
        total_moles += molality * ratio;
    if (!std::isfinite(total_moles) || total_moles <= 0.0) {
        throw std::invalid_argument("invalid neutral reference composition total");
    }
    std::vector<double> composition(declaration.solute_molality_ratios.size(), 0.0);
    composition[solvent] = solvent_moles / total_moles;
    for (std::size_t component = 0; component < composition.size(); ++component) {
        if (component == solvent) continue;
        composition[component] = molality * declaration.solute_molality_ratios[component]
            / total_moles;
    }
    return composition;
}

}  // namespace

std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const Parameters& parameters,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality
) {
    // EqID: formula_unit_composition
    const auto& molar_masses = parameters.molar_mass_kg_per_mol;
    const auto& charges = parameters.charge_number;
    if (solvent_fractions.size() != molar_masses.size()
        || solvent_fractions.size() != charges.size()
        || ion_indices[0] >= solvent_fractions.size()
        || ion_indices[1] >= solvent_fractions.size()
        || ion_indices[0] == ion_indices[1]
        || stoichiometry[0] <= 0 || stoichiometry[1] <= 0
        || !std::isfinite(molality) || molality <= 0.0) {
        throw std::invalid_argument("invalid reference composition inputs");
    }
    if (charges[ion_indices[0]] == 0 || charges[ion_indices[1]] == 0
        || static_cast<std::int64_t>(stoichiometry[0]) * charges[ion_indices[0]]
            + static_cast<std::int64_t>(stoichiometry[1]) * charges[ion_indices[1]] != 0) {
        throw std::invalid_argument("reference formula must be charge neutral with two ions");
    }
    double fraction_sum = 0.0;
    double solvent_mass = 0.0;
    double solvent_charge = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        if (!std::isfinite(solvent_fractions[component])
            || solvent_fractions[component] < 0.0
            || !std::isfinite(molar_masses[component])
            || (solvent_fractions[component] > 0.0 && molar_masses[component] <= 0.0)) {
            throw std::invalid_argument("invalid reference solvent fractions");
        }
        fraction_sum += solvent_fractions[component];
        solvent_mass += solvent_fractions[component] * molar_masses[component];
        solvent_charge += solvent_fractions[component] * charges[component];
    }
    if (std::abs(fraction_sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value
        || solvent_fractions[ion_indices[0]] != 0.0
        || solvent_fractions[ion_indices[1]] != 0.0
        || !(solvent_mass > 0.0) || !std::isfinite(solvent_mass)
        || std::abs(solvent_charge) > policy::COMPOSITION_SUM_TOLERANCE.value) {
        throw std::invalid_argument("invalid reference solvent composition or molar mass");
    }
    const double solvent_moles = 1.0 / solvent_mass;
    const double total_moles = solvent_moles
        + molality * (static_cast<double>(stoichiometry[0]) + stoichiometry[1]);
    std::vector<double> composition(solvent_fractions.size(), 0.0);
    for (std::size_t component = 0; component < composition.size(); ++component) {
        composition[component] = solvent_fractions[component] * solvent_moles / total_moles;
    }
    composition[ion_indices[0]] = stoichiometry[0] * molality / total_moles;
    composition[ion_indices[1]] = stoichiometry[1] * molality / total_moles;
    return composition;
}

double reaction_reference_log_scale(
    const HelmholtzModel& model,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry
) {
    const double solvent_mass = validate_neutral_reference_declaration(
        model.parameters(), declaration, stoichiometry);
    double solute_stoichiometry = 0.0;
    const std::size_t solvent = *declaration.solvent_index;
    for (std::size_t component = 0; component < stoichiometry.size(); ++component)
        if (component != solvent) solute_stoichiometry += stoichiometry[component];
    const double scale = solute_stoichiometry
        * std::log(solvent_mass * declaration.standard_molality_mol_per_kg);
    if (!std::isfinite(scale)) throw std::invalid_argument("nonfinite neutral reference scale");
    return scale;
}

NeutralReferenceReactionDiagnostic reaction_reference(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry,
    const ReactionReference& reference
) {
    NeutralReferenceReactionDiagnostic result;
    try {
        const double solvent_mass = validate_neutral_reference_declaration(
            model.parameters(), declaration, stoichiometry);
        if (!std::isfinite(temperature) || temperature <= 0.0
            || !std::isfinite(system_pressure) || system_pressure <= 0.0) {
            result.status = NeutralReferenceStatus::InvalidDeclaration;
            result.message = "neutral reference temperature or system pressure is invalid";
            return result;
        }
        const double pressure = reference.reference_pressure_pa.value_or(system_pressure);
        if (!std::isfinite(pressure) || pressure <= 0.0) {
            result.status = NeutralReferenceStatus::InvalidDeclaration;
            result.message = "neutral reference pressure is invalid";
            return result;
        }
        switch (reference.source_basis) {
        case ReferenceSourceBasis::CommonMolalityInfiniteDilution:
        case ReferenceSourceBasis::RawMoleFractionInfiniteDilution:
            break;
        default:
            result.status = NeutralReferenceStatus::UnsupportedSourceBasis;
            result.message = "unsupported neutral reference source basis";
            return result;
        }
        result.reference_pressure_pa = pressure;
        std::optional<double> anchor;
        std::array<double, NEUTRAL_REFERENCE_MOLALITIES.size()> contractions{};
        for (std::size_t point = 0; point < NEUTRAL_REFERENCE_MOLALITIES.size(); ++point) {
            const double molality = NEUTRAL_REFERENCE_MOLALITIES[point];
            const auto composition = neutral_reference_composition(
                declaration, solvent_mass, molality);
            State state;
            try {
                state = model.state_at_pressure(
                    temperature, pressure, composition, DensityBranch::Liquid, anchor);
            } catch (const std::exception& error) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = std::string("neutral reference density unavailable: ")
                    + error.what();
                return result;
            }
            const bool cold_unique_root = point == 0
                && state.density_branch == DensityBranch::Unspecified
                && state.stable_root_count.value_or(0) == 1;
            if (state.density_branch != DensityBranch::Liquid && !cold_unique_root) {
                result.status = NeutralReferenceStatus::BranchChange;
                result.message = "neutral reference left the anchored liquid branch";
                return result;
            }
            if (!std::isfinite(state.molar_density) || state.molar_density <= 0.0
                || !std::isfinite(state.pressure_density_derivative)
                || state.pressure_density_derivative <= 0.0
                || state.residual_chemical_potential_over_rt.size()
                    != model.parameters().component_count()) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference state has invalid density response";
                return result;
            }
            const double log_density = std::log(state.molar_density);
            double contraction = 0.0;
            for (std::size_t component = 0;
                 component < state.residual_chemical_potential_over_rt.size(); ++component) {
                const double chemical = state.residual_chemical_potential_over_rt[component];
                if (!std::isfinite(chemical)) {
                    result.status = NeutralReferenceStatus::Nonfinite;
                    result.message = "neutral reference residual chemical potential is nonfinite";
                    return result;
                }
                contraction += stoichiometry[component] * (log_density + chemical);
            }
            if (!std::isfinite(contraction)) {
                result.status = NeutralReferenceStatus::Nonfinite;
                result.message = "neutral reference contraction is nonfinite";
                return result;
            }
            contractions[point] = contraction;
            anchor = state.molar_density;
        }
        const double observed_terminal_change = std::max(
            std::abs(contractions[4] - contractions[3]),
            std::abs(contractions[3] - contractions[2]));
        result.observed_terminal_change = observed_terminal_change;
        if (!std::isfinite(observed_terminal_change)) {
            result.status = NeutralReferenceStatus::Nonfinite;
            result.message = "neutral reference terminal change is nonfinite";
            return result;
        }
        if (observed_terminal_change > policy::REFERENCE_CONVERGENCE_TOLERANCE.value) {
            result.status = NeutralReferenceStatus::TerminalChangeExceeded;
            result.message = "neutral reference terminal change exceeds tolerance";
            return result;
        }
        result.terminal_contraction = contractions.back();
        result.status = NeutralReferenceStatus::Available;
        return result;
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        result.status = message.find("charged") != std::string::npos
                || message.find("charge neutral") != std::string::npos
            ? NeutralReferenceStatus::NonNeutralDeclaration
            : NeutralReferenceStatus::InvalidDeclaration;
        result.message = message;
        return result;
    } catch (const std::exception& error) {
        result.status = NeutralReferenceStatus::Nonfinite;
        result.message = error.what();
        return result;
    }
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
            model.parameters(),
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
            model.parameters(),
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
    core::Scope scope("electrolyte_activity");
    const auto composition = formula_composition(
        solvent_fractions,
        model.parameters(),
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
    const double coefficient_sum = static_cast<double>(stoichiometry[0]) + stoichiometry[1];
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
