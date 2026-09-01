#pragma once

#include "phase_stability_stage3.hpp"

#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {

struct PhaseStabilityLpCertificate {
    bool primal_feasible = false;
    bool dual_feasible = false;
    double primal_residual_inf = std::numeric_limits<double>::infinity();
    double dual_residual_inf = std::numeric_limits<double>::infinity();
    double complementarity_inf = std::numeric_limits<double>::infinity();
};

struct PhaseStabilityStage4CutEvidence {
    int id = -1;
    double intercept = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> slopes;
};

struct PhaseStabilityStage4Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    int major_iteration = -1;
    int upper_solve_count = -1;
    std::optional<double> upper_bound;
    std::optional<std::vector<double>> multipliers;
    std::vector<int> active_cut_ids;
    std::vector<PhaseStabilityStage4CutEvidence> cut_snapshot;
    std::optional<PhaseStabilityCoordinates> coordinate_snapshot;
    std::optional<PhaseStabilityLpCertificate> certificate;
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] PhaseStabilityStage4Result run_phase_stability_stage4(
    PhaseStabilityPersistentState& state
);

}  // namespace epcsaft_engine::equilibrium
