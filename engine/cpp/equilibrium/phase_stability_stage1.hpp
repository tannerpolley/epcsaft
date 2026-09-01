#pragma once

#include "phase_stability_coordinates.hpp"

#include <optional>

namespace epcsaft_engine::equilibrium {

struct PhaseStabilityStage1Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
    std::string provider_fingerprint;
    std::optional<PhaseStabilityCoordinates> coordinates;
    std::optional<std::vector<double>> independent_feed;
    std::optional<std::array<double, 2>> feed_volume_bounds;
    std::optional<PhaseStabilityVolumeBoundsEvaluator> volume_bounds;
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] PhaseStabilityStage1Result run_phase_stability_stage1(
    const std::vector<std::string>& component_ids,
    const std::vector<double>& charges,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& physical_feed,
    const PhaseStabilityPhysicalVolumeBoundsEvaluator& physical_volume_bounds,
    std::string provider_fingerprint,
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN()
);

}  // namespace epcsaft_engine::equilibrium
