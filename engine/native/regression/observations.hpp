#pragma once

#include "equilibrium/problem.hpp"
#include "regression/coordinates.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::regression {

// The built-in observation inventory (#46). Values come from the shared owners only:
// Density and total c_p from the Properties pressure closure, other kinds from one declared
// specification of the single equilibrium driver. PartialPressure is P sum_{i in S} y_i.
enum class ObservationKind {
    Density, IsobaricHeatCapacity, PhaseDensity, EquilibriumPressure, EquilibriumTemperature, PhaseMoleFraction,
    PhaseFraction, PartialPressure, PhaseAmountRatio, LogEquilibriumConstant, AdditionEnthalpy, SystemEnthalpy, FeedEnthalpy
};
inline bool pressure_closed(ObservationKind kind) {
    return kind == ObservationKind::Density || kind == ObservationKind::IsobaricHeatCapacity || kind == ObservationKind::FeedEnthalpy;
}
enum class ResidualForm { Difference, LogRatio };
enum class ObservationStatus {
    Available, DensityUnavailable, CentralSolveUnavailable, TopologyEvent, MissingMolarMass,
    ThermochemistryUnavailable, NonpositiveLogArgument, ActionUnavailable, Nonfinite, InadmissibleTrial
};
const char* observation_status_name(ObservationStatus status) noexcept;

struct Observation final {
    ObservationKind kind = ObservationKind::Density;
    // Density and total c_p: T (K), P (Pa), normalized composition and the declared branch.
    double temperature = 0.0;
    double pressure = 0.0;
    std::vector<double> composition;
    DensityBranch branch = DensityBranch::Unspecified;
    // Equilibrium kinds: the declared problem; `phase` indexes its phases, `species` is a
    // nonempty set inside that phase's support (PhaseMoleFraction and PartialPressure sum the set).
    equilibrium::EquilibriumProblem problem;
    equilibrium::EquilibriumProblem final_problem;
    double feed_amount = 0.0;
    std::string reaction;  // LogEquilibriumConstant: named source record at problem.temperature.value.
    std::size_t phase = 0;
    std::vector<std::size_t> species;
    std::size_t denominator_phase = 0;
    std::vector<std::size_t> denominator_species;
    // Density and PhaseDensity: kg/m^3 through the fixed molar masses instead of mol/m^3.
    bool mass_basis = false;
    // The reported value in the prediction's unit, the residual form and, for Difference,
    // the positive nondimensionalization scale in that unit (not an uncertainty).
    double observed = 0.0;
    ResidualForm form = ResidualForm::Difference;
    double scale = 1.0;
};

struct Prediction final {
    std::optional<double> value;
    std::string unit;
    ObservationStatus status = ObservationStatus::Nonfinite;
    std::string message;
    std::optional<double> residual;  // the declared form, unweighted
    std::optional<core::Convergence> convergence;
    // Pressure closure: the branch the closure resolved (Unspecified when one stable root exists) and the count.
    std::optional<DensityBranch> branch;
    std::optional<std::size_t> stable_root_count;
    std::optional<bool> requested_tolerance_met;
    std::optional<equilibrium::EquilibriumResult> solved;  // the accepted central solve
    std::optional<equilibrium::CompiledSystem> compiled;    // the problem it solved, reused by differentiate
    std::optional<equilibrium::EquilibriumResult> final_solved;
    std::optional<equilibrium::CompiledSystem> final_compiled;
};

// d(value)/d(theta) in EOS-row order, then the requested reaction-coordinate order;
// value units per native coefficient unit, or empty with the blocking action status.
struct PredictionGradient final {
    std::vector<double> value;
    ObservationStatus status = ObservationStatus::ActionUnavailable;
    std::string message;
};

// Admission by name: finite inputs, phase/species indices inside the declared problem.
using equilibrium::same_reaction;
bool same_state(const Observation& left, const Observation& right);
void validate(const Observation& observation, const HelmholtzModel& model,
              const equilibrium::CompiledSystem* prepared = nullptr);
void project_prediction(const HelmholtzModel& model, const Observation& observation, Prediction& prediction);
Prediction predict(
    const HelmholtzModel& model, const Observation& observation,
    const std::optional<equilibrium::EquilibriumResult>& start = std::nullopt,
    std::optional<equilibrium::CompiledSystem> prepared = std::nullopt
);
PredictionGradient differentiate(
    const HelmholtzModel& model, const Observation& observation, const Prediction& prediction
);
std::vector<PredictionGradient> differentiate_group(
    const HelmholtzModel& model, const std::vector<const Observation*>& observations, const Prediction& prediction,
    const CoordinateSet& coordinates = {}
);
// r = (y - y_obs) / scale, or ln(y / y_obs) for positive y and y_obs; dr/dy = 1/scale or 1/y.
std::optional<double> residual(const Observation& observation, double value, ObservationStatus& status);
double residual_slope(const Observation& observation, double value);

}  // namespace epcsaft::regression
