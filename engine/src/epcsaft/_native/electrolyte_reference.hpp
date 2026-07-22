#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft::detail {

inline constexpr std::array<double, 4> REFERENCE_MOLALITIES = {
    1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12,
};
inline constexpr double REFERENCE_CONVERGENCE_TOLERANCE = 5.0e-5;

struct ReferenceObservation final {
    std::array<double, 2> log_fugacity;
    std::optional<std::array<double, 2>> derivative;
};

struct ReferenceSequenceResult final {
    std::array<double, 2> log_fugacity;
    std::optional<double> pure_solvent_log_fugacity_coefficient;
    std::optional<std::array<double, 2>> derivative;
    double molality;
    double convergence_error;
    std::optional<double> derivative_convergence_error;
};

inline std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const std::vector<double>& molar_masses_kg_per_mol,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality
) {
    if (solvent_fractions.size() != molar_masses_kg_per_mol.size()
        || ion_indices[0] >= solvent_fractions.size()
        || ion_indices[1] >= solvent_fractions.size()
        || ion_indices[0] == ion_indices[1]
        || stoichiometry[0] <= 0 || stoichiometry[1] <= 0
        || !std::isfinite(molality) || molality <= 0.0) {
        throw std::invalid_argument("invalid reference composition inputs");
    }
    double solvent_molar_mass = 0.0;
    double fraction_sum = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        if (!std::isfinite(solvent_fractions[component])
            || solvent_fractions[component] < 0.0) {
            throw std::invalid_argument("invalid reference solvent fractions");
        }
        fraction_sum += solvent_fractions[component];
        solvent_molar_mass += solvent_fractions[component]
            * molar_masses_kg_per_mol[component];
    }
    if (std::abs(fraction_sum - 1.0) > 1.0e-12
        || solvent_fractions[ion_indices[0]] != 0.0
        || solvent_fractions[ion_indices[1]] != 0.0
        || !std::isfinite(solvent_molar_mass) || solvent_molar_mass <= 0.0) {
        throw std::invalid_argument("invalid reference solvent molar mass");
    }
    const double solvent_moles = 1.0 / solvent_molar_mass;
    const double total_moles = solvent_moles
        + molality * static_cast<double>(stoichiometry[0] + stoichiometry[1]);
    std::vector<double> composition(solvent_fractions.size());
    for (std::size_t component = 0; component < composition.size(); ++component) {
        composition[component] = solvent_fractions[component] * solvent_moles
            / total_moles;
    }
    for (std::size_t ion = 0; ion < ion_indices.size(); ++ion) {
        composition[ion_indices[ion]] = stoichiometry[ion] * molality / total_moles;
    }
    return composition;
}

template <typename Evaluator, typename PureSolventEvaluator>
ReferenceSequenceResult electrolyte_reference_sequence(
    const std::vector<double>& solvent_fractions,
    const std::vector<double>& molar_masses_kg_per_mol,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    Evaluator evaluator,
    PureSolventEvaluator pure_solvent_evaluator
) {
    const auto pure_solvent_log_fugacity = pure_solvent_evaluator(solvent_fractions);
    if (pure_solvent_log_fugacity.has_value()
        && !std::isfinite(*pure_solvent_log_fugacity)) {
        throw std::runtime_error("nonfinite pure-solvent reference fugacity");
    }
    std::optional<ReferenceObservation> previous;
    std::optional<ReferenceObservation> current;
    double convergence_error = INFINITY;
    std::optional<double> derivative_convergence_error;

    for (double molality : REFERENCE_MOLALITIES) {
        current = evaluator(
            formula_composition(
                solvent_fractions,
                molar_masses_kg_per_mol,
                ion_indices,
                stoichiometry,
                molality
            )
        );
        for (double value : current->log_fugacity) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("nonfinite reference fugacity");
            }
        }
        if (previous.has_value()) {
            convergence_error = std::max(
                std::abs(current->log_fugacity[0] - previous->log_fugacity[0]),
                std::abs(current->log_fugacity[1] - previous->log_fugacity[1])
            );
            if (current->derivative.has_value() && previous->derivative.has_value()) {
                derivative_convergence_error = std::max(
                    std::abs(
                        (*current->derivative)[0] - (*previous->derivative)[0]
                    ),
                    std::abs(
                        (*current->derivative)[1] - (*previous->derivative)[1]
                    )
                );
            }
        }
        previous = current;
    }
    if (!current.has_value() || convergence_error > REFERENCE_CONVERGENCE_TOLERANCE) {
        throw std::runtime_error("reference limit did not converge");
    }
    if (current->derivative.has_value()) {
        if (!derivative_convergence_error.has_value()) {
            throw std::runtime_error("reference derivative limit unavailable");
        }
        for (double value : *current->derivative) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("nonfinite reference derivative");
            }
        }
        const double derivative_scale = std::max(
            1.0,
            std::max(
                std::abs((*current->derivative)[0]),
                std::abs((*current->derivative)[1])
            )
        );
        if (*derivative_convergence_error
            > REFERENCE_CONVERGENCE_TOLERANCE * derivative_scale) {
            throw std::runtime_error("reference derivative limit did not converge");
        }
    }
    return {
        current->log_fugacity,
        pure_solvent_log_fugacity,
        current->derivative,
        REFERENCE_MOLALITIES.back(),
        convergence_error,
        derivative_convergence_error,
    };
}

}  // namespace epcsaft::detail
