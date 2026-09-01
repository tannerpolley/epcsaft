#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "eos.hpp"
#include "sensitivity_linear_solve.hpp"

namespace epcsaft_engine::equilibrium {

inline constexpr double kSolvedStateSensitivityConditionNumberMax =
    interfaces::kSolvedStateSensitivityConditionNumberMax;

[[nodiscard]] inline double solved_state_sensitivity_backward_error_max(
    std::size_t dimension
) {
    return interfaces::solved_state_sensitivity_backward_error_max(dimension);
}

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
    std::size_t eos_evaluation_count = 0;
    double scaled_condition_number_2 = std::numeric_limits<double>::infinity();
    double linear_solve_backward_error = std::numeric_limits<double>::infinity();
    std::size_t reference_phase_index = 0;
};

[[nodiscard]] bool solve_sensitivity_linear_system(
    SolvedStateSensitivity& result,
    std::size_t parameter_count
);

[[nodiscard]] SolvedStateSensitivity evaluate_solved_state_sensitivity(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases
);

[[nodiscard]] SolvedStateSensitivity evaluate_declared_state_sensitivity(
    const epcsaft::NativeModel& eos,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<SolvedPhaseState>& phases
);

}  // namespace epcsaft_engine::equilibrium
