#pragma once

#include "phase_stability_stage7.hpp"
#include "provider.hpp"

#include <cstddef>

namespace epcsaft_engine::equilibrium {

struct FlashInput {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> overall_mole_fractions;
};

struct PhaseStabilityThermodynamicAccess {
    const ProviderContext* provider = nullptr;
    std::vector<std::string> component_ids;
    std::vector<double> charges;
    std::string provider_fingerprint;
    PhaseStabilityStateEvaluator evaluate;
    PhaseStabilityPhysicalVolumeBoundsEvaluator volume_bounds_physical;
    PhaseStabilityPackingFractionEvaluator packing_fraction;
    double total_ion_mole_fraction_max =
        std::numeric_limits<double>::quiet_NaN();
};

struct PhaseStabilityPhase {
    std::uint64_t stable_id = 0;
    double phase_fraction = 0.0;
    std::vector<double> independent_modified_fractions;
    std::vector<double> physical_fractions_provider_order;
    double volume = 0.0;
    double packing_fraction = 0.0;
    double helmholtz_over_rt_reference_amount = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> chemical_potentials_over_rt;
    double reduced_gibbs = 0.0;
    std::vector<double> reduced_gibbs_gradient;
};

struct PhaseStabilityDiscoveryResult {
    std::string outcome = "indeterminate";
    std::string failure_stage;
    std::string failure_reason;
    std::string globality_certificate = "not_guaranteed";
    PhaseStabilityStage1Result step1;
    std::optional<PhaseStabilityStage2Result> step2;
    std::optional<PhaseStabilityStage3Result> step3;
    std::vector<PhaseStabilityStage4Result> step4_history;
    std::vector<PhaseStabilityStage5Result> step5_history;
    std::vector<PhaseStabilityStage6Result> step6_history;
    std::vector<PhaseStabilityStage7Result> step7_history;
    std::vector<PhaseStabilityPhase> phases;
    std::optional<double> total_free_energy_over_rt;
    std::optional<double> material_balance_inf;
    std::optional<double> electroneutrality_inf;
    std::optional<double> pressure_residual_inf;
    std::optional<double> chemical_potential_residual_inf;
    std::optional<double> kkt_stationarity_inf;
    std::optional<double> phase_distinctness;
    std::vector<PhaseStabilityStageTiming> step_timings;
};

[[nodiscard]] PhaseStabilityDiscoveryResult run_phase_stability(
    const PhaseStabilityThermodynamicAccess& thermodynamics,
    const FlashInput& input,
    const PhaseStabilityResourceProfile& resources
);

}  // namespace epcsaft_engine::equilibrium
