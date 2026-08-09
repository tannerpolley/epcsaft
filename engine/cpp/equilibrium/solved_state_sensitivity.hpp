#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "provider.hpp"

namespace epcsaft_engine::equilibrium {

inline constexpr double kSolvedStateSensitivityConditionNumberMax = 1.0e8;

struct SolvedPhaseState {
    std::vector<double> amounts_mol;
    double volume_m3 = 0.0;
};

struct SolvedStateSensitivity {
    std::string status = "unavailable";
    std::string failure_reason;
    std::vector<std::string> residual_identities;
    std::vector<std::string> state_identities;
    std::vector<std::string> output_identities;
    std::vector<std::string> parameter_identities;
    std::vector<double> residual_state_jacobian;
    std::vector<double> residual_parameter_jacobian;
    std::vector<double> state_parameter_jacobian;
    std::vector<double> output_parameter_jacobian;
    std::size_t rank = 0;
    std::size_t dimension = 0;
    double scaled_condition_number_2 = std::numeric_limits<double>::infinity();
    double linear_solve_backward_error = std::numeric_limits<double>::infinity();
    std::size_t reference_phase_index = 0;
};

[[nodiscard]] SolvedStateSensitivity evaluate_solved_state_sensitivity(
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases
);

}  // namespace epcsaft_engine::equilibrium
