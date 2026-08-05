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

struct NoopReferenceCheckpoint final {
    // EqExclude: transport and buffer mechanics
    void operator()() const {}
};

inline constexpr std::array<double, 4> REFERENCE_MOLALITIES = {
    1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12,
};
// Keep the historical reported reference molality at 1e-12 while using one
// lower point as a convergence witness for associating organic solvents.
inline constexpr double REFERENCE_CONVERGENCE_WITNESS_MOLALITY = 1.0e-14;
inline constexpr double REFERENCE_CONVERGENCE_TOLERANCE = 5.0e-5;

struct ReferenceObservation final {
    std::array<double, 2> log_fugacity;
    std::optional<std::vector<double>> derivative;
};

struct ReferenceSequenceResult final {
    std::array<double, 2> log_fugacity;
    std::optional<std::vector<double>> derivative;
    double molality;
    double convergence_error;
    std::optional<double> derivative_convergence_error;
};

struct VectorReferenceSequenceResult final {
    std::vector<double> values;
    std::vector<double> convergence_errors;
    double molality;
    double convergence_error;
};

template <typename Evaluator, typename NeedsWitness>
// EqID: reference_molality_sequence
VectorReferenceSequenceResult reference_limit_sequence(
    Evaluator evaluator,
    NeedsWitness needs_witness
) {
    std::optional<std::vector<double>> previous;
    std::optional<std::vector<double>> current;
    std::vector<double> convergence_errors;
    double convergence_error = INFINITY;
    const auto observe = [&](double molality) {
        current = evaluator(molality);
        if (current->empty()
            || !std::all_of(current->begin(), current->end(), [](double value) {
                return std::isfinite(value);
            })) {
            throw std::runtime_error("nonfinite source-reference limit value");
        }
        if (previous.has_value()) {
            if (previous->size() != current->size()) {
                throw std::runtime_error(
                    "source-reference limit dimension changed"
                );
            }
            convergence_errors.assign(current->size(), 0.0);
            convergence_error = 0.0;
            for (std::size_t index = 0; index < current->size(); ++index) {
                convergence_errors[index] =
                    std::abs((*current)[index] - (*previous)[index]);
                convergence_error =
                    std::max(convergence_error, convergence_errors[index]);
            }
        }
        previous = current;
    };
    for (double molality : REFERENCE_MOLALITIES) {
        observe(molality);
    }
    if (!current.has_value()) {
        throw std::runtime_error("reference limit did not converge");
    }
    const auto terminal = *current;
    if (needs_witness(terminal, convergence_errors)) {
        observe(REFERENCE_CONVERGENCE_WITNESS_MOLALITY);
    }
    if (needs_witness(*current, convergence_errors)) {
        throw std::runtime_error("reference limit did not converge");
    }
    return {
        terminal,
        convergence_errors,
        REFERENCE_MOLALITIES.back(),
        convergence_error,
    };
}

template <typename Evaluator>
VectorReferenceSequenceResult reference_limit_sequence(Evaluator evaluator) {
    return reference_limit_sequence(
        std::move(evaluator),
        [](const std::vector<double>&, const std::vector<double>& errors) {
            return errors.empty()
                || *std::max_element(errors.begin(), errors.end())
                    > REFERENCE_CONVERGENCE_TOLERANCE;
        }
    );
}

// EqID: formula_unit_composition
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

template <typename Evaluator, typename Checkpoint = NoopReferenceCheckpoint>
// EqID: reference_molality_sequence
ReferenceSequenceResult electrolyte_reference_sequence(
    const std::vector<double>& solvent_fractions,
    const std::vector<double>& molar_masses_kg_per_mol,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    Evaluator evaluator,
    Checkpoint checkpoint = {}
) {
    std::optional<std::size_t> derivative_size;
    const auto sequence = reference_limit_sequence(
        [&](double molality) {
            checkpoint();
            const auto observation = evaluator(
                formula_composition(
                    solvent_fractions,
                    molar_masses_kg_per_mol,
                    ion_indices,
                    stoichiometry,
                    molality
                )
            );
            const std::size_t observed_derivative_size =
                observation.derivative.has_value()
                ? observation.derivative->size()
                : 0u;
            if (derivative_size.has_value()
                && *derivative_size != observed_derivative_size) {
                throw std::runtime_error("reference derivative size changed");
            }
            derivative_size = observed_derivative_size;
            std::vector<double> values{
                observation.log_fugacity[0], observation.log_fugacity[1]
            };
            if (observation.derivative.has_value()) {
                values.insert(
                    values.end(),
                    observation.derivative->begin(),
                    observation.derivative->end()
                );
            }
            return values;
        },
        [](const std::vector<double>&, const std::vector<double>& errors) {
            return errors.size() < 2
                || std::max(errors[0], errors[1])
                    > REFERENCE_CONVERGENCE_TOLERANCE;
        }
    );
    const double convergence_error = std::max(
        sequence.convergence_errors[0], sequence.convergence_errors[1]
    );
    std::optional<std::vector<double>> derivative;
    std::optional<double> derivative_convergence_error;
    if (sequence.values.size() > 2) {
        derivative = std::vector<double>(
            sequence.values.begin() + 2, sequence.values.end()
        );
        derivative_convergence_error = *std::max_element(
            sequence.convergence_errors.begin() + 2,
            sequence.convergence_errors.end()
        );
        double derivative_scale = 1.0;
        for (double value : *derivative) {
            derivative_scale = std::max(derivative_scale, std::abs(value));
        }
        if (*derivative_convergence_error
            > REFERENCE_CONVERGENCE_TOLERANCE * derivative_scale) {
            throw std::runtime_error("reference derivative limit did not converge");
        }
    }
    return {
        {sequence.values[0], sequence.values[1]},
        derivative,
        REFERENCE_MOLALITIES.back(),
        convergence_error,
        derivative_convergence_error,
    };
}

}  // namespace epcsaft::detail
