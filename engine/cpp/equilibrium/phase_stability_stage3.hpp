#pragma once

#include "phase_stability_stage2.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <optional>

namespace epcsaft_engine::equilibrium {

struct PhaseStabilityMPoint {
    std::uint64_t insertion_id = 0;
    std::vector<double> independent_modified_fractions;
    double volume = 0.0;
    double packing_fraction = std::numeric_limits<double>::quiet_NaN();
    double reduced_gibbs = 0.0;
    std::vector<double> reduced_gibbs_gradient;
    std::string origin;
};

struct PhaseStabilityPersistentState {
    int major_iteration = 0;
    int upper_solve_count = 0;
    double feed_reduced_gibbs = 0.0;
    double upper_bound = 0.0;
    double lower_value = -std::numeric_limits<double>::infinity();
    std::uint64_t next_start_ordinal = 0;
    bool step5_requires_new_member = false;
    std::optional<std::array<double, 2>> step5_volume_bounds;
    std::optional<PhaseStabilityCoordinates> coordinates;
    std::vector<double> feed;
    std::vector<double> multipliers;
    std::vector<PhaseStabilityMPoint> M;
};

using PhaseStabilityPressureRootEvaluator = std::function<PhaseStabilityPressureEnvelopeResult(
    const std::vector<double>&
)>;

struct PhaseStabilityStage3Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<PhaseStabilityPersistentState> state;
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] PhaseStabilityStage3Result run_phase_stability_stage3(
    const PhaseStabilityStage1Result& step1,
    const PhaseStabilityStage2Result& step2,
    const PhaseStabilityPressureRootEvaluator& pressure_roots
);

}  // namespace epcsaft_engine::equilibrium
