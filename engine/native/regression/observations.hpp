#pragma once

#include "equilibrium/problem.hpp"
#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::regression {

// The built-in observation inventory (#46). Values come from the shared owners only:
// Density from the Properties pressure closure, every other kind from one declared
// specification of the single equilibrium driver.
enum class ObservationKind {
    Density, PhaseDensity, EquilibriumPressure, EquilibriumTemperature, PhaseMoleFraction,
    PhaseFraction
};
enum class ResidualForm { Difference, LogRatio };
enum class ObservationStatus {
    Available, DensityUnavailable, CentralSolveUnavailable, TopologyEvent, MissingMolarMass,
    NonpositiveLogArgument, ActionUnavailable, Nonfinite
};
const char* observation_status_name(ObservationStatus status) noexcept;

struct Observation final {
    ObservationKind kind = ObservationKind::Density;
    // Density: T (K), P (Pa), normalized composition and the declared branch.
    double temperature = 0.0;
    double pressure = 0.0;
    std::vector<double> composition;
    DensityBranch branch = DensityBranch::Unspecified;
    // Every other kind: the declared problem; `phase` indexes its phases, `species` is a
    // nonempty set inside that phase's support (PhaseMoleFraction sums the set).
    equilibrium::EquilibriumProblem problem;
    std::size_t phase = 0;
    std::vector<std::size_t> species;
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
    // Density: the branch the closure resolved (Unspecified when one stable root exists) and the count.
    std::optional<DensityBranch> branch;
    std::optional<std::size_t> stable_root_count;
    std::optional<bool> requested_tolerance_met;
    std::optional<equilibrium::EquilibriumResult> solved;  // the accepted central solve
};

// d(value)/d(theta) per active parameter of the model, in the value's unit per native
// parameter unit; empty with the blocking status when any required action is unavailable.
struct PredictionGradient final {
    std::vector<double> value;
    ObservationStatus status = ObservationStatus::ActionUnavailable;
    std::string message;
};

// Admission by name: finite inputs, phase/species indices inside the declared problem.
void validate(const Observation& observation, const HelmholtzModel& model);
Prediction predict(
    const HelmholtzModel& model, const Observation& observation,
    const std::optional<equilibrium::EquilibriumResult>& start = std::nullopt
);
PredictionGradient differentiate(
    const HelmholtzModel& model, const Observation& observation, const Prediction& prediction
);
// r = (y - y_obs) / scale, or ln(y / y_obs) for positive y and y_obs; dr/dy = 1/scale or 1/y.
std::optional<double> residual(const Observation& observation, double value, ObservationStatus& status);
double residual_slope(const Observation& observation, double value);

}  // namespace epcsaft::regression
