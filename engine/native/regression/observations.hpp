#pragma once

#include "equilibrium/problem.hpp"
#include "regression/coordinates.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::regression {

// The built-in observation inventory. Values come from the shared owners only:
// Density and total c_p from the Properties pressure closure, other kinds from one declared
// specification of the single equilibrium driver. PartialPressure is P sum_{i in S} y_i.
// MixtureOsmoticCoefficient is -ln a_w / (M_w sum_j m_j) of a full electroneutral liquid composition in one solvent.
enum class ObservationKind {
    MeanIonicActivityCoefficient, OsmoticCoefficient, MixtureOsmoticCoefficient, SolvationGibbsEnergy, TransferGibbsEnergy, BulkRelativePermittivity,
    Density, IsobaricHeatCapacity, EquilibriumIsobaricHeatCapacity, PhaseDensity, EquilibriumPressure, EquilibriumTemperature, PhaseMoleFraction,
    PhaseFraction, PartialPressure, PhaseAmountRatio, LogEquilibriumConstant, AdditionEnthalpy, SaturationCappedAdditionEnthalpy, SystemEnthalpy, FeedEnthalpy
};
inline bool addition(ObservationKind kind) {
    return kind == ObservationKind::AdditionEnthalpy || kind == ObservationKind::SaturationCappedAdditionEnthalpy;
}
inline bool pressure_closed(ObservationKind kind) {
    return kind == ObservationKind::Density || kind == ObservationKind::IsobaricHeatCapacity || kind == ObservationKind::FeedEnthalpy;
}
inline bool electrolyte_observation(ObservationKind kind) {
    return kind <= ObservationKind::BulkRelativePermittivity;
}
enum class ResidualForm { Difference, LogRatio };
enum class ObservationStatus {
    NotEvaluated, Available, ReferenceUnavailable, DensityUnavailable, CentralSolveUnavailable, TopologyEvent, MissingMolarMass,
    ThermochemistryUnavailable, NonpositiveLogArgument, ActionUnavailable, Nonfinite, InadmissibleTrial
};
const char* observation_status_name(ObservationStatus status) noexcept;

struct Observation final {
    ObservationKind kind = ObservationKind::Density;
    // Density and total c_p: T (K), P (Pa), normalized composition and the declared branch.
    double temperature = 0.0;
    double pressure = 0.0;
    std::vector<double> composition;
    std::array<std::size_t, 2> ion_indices{};
    std::array<int, 2> stoichiometry{};
    double molality = 0.0;
    std::optional<std::vector<double>> finite_solute_molalities;
    std::vector<double> reference_composition;
    DensityBranch branch = DensityBranch::Unspecified;
    // Equilibrium kinds: the declared problem; `phase` indexes its phases, `species` is a
    // nonempty set inside that phase's support (PhaseMoleFraction and PartialPressure sum the set).
    equilibrium::EquilibriumProblem problem;
    equilibrium::EquilibriumProblem final_problem;
    double feed_amount = 0.0;
    double basis_amount = 0, loading_cap = 0;
    std::array<double, 2> loading_bounds{};
    std::vector<double> loading_grid;
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

struct HydrationHold final {
    Observation reference;  // SolvationGibbsEnergy, target in observed (J/mol)
    double lower = 0.0, upper = 0.0;
    double energy_tolerance = 0.0, diameter_tolerance = 0.0;
    int maximum_iterations = 0;
};
// correction is the last Newton step |H/H_d| (Angstrom); uncertainty bounds the achieved diameter
// error, min(correction, bracket width), the quantity the diameter criterion accepts.
struct HoldResult final {
    double energy_error = 0.0, correction = 0.0, uncertainty = 0.0, derivative = 0.0;
    std::vector<double> parameter_derivative;
    int iterations = 0;
    core::NumericalRank rank;
};
struct EvaluatedCoordinates final {
    Parameters parameters;
    std::vector<ActiveParameter> active;
    std::vector<double> physical, map;  // independent EOS, held Born, then reaction rows; row-major D
    std::vector<HoldResult> holds;
    std::map<std::string, ReactionLogPolynomial> reaction_laws;
};
EvaluatedCoordinates evaluate_coordinates(const Parameters& parameters, const CoordinateSet& set,
    const std::vector<double>& optimizer, const std::vector<HydrationHold>& holds,
    const std::optional<ThermochemistryRecord>& thermochemistry = {}, const policy::SolverOptions& solver_options = {});

struct Prediction final {
    policy::SolverOptions solver_options;
    std::optional<InfiniteDilutionReference> reference, second_reference;
    std::optional<ElectrolyteActivity> activity;
    std::optional<State> state, reference_state;
    std::optional<double> value;
    std::string unit;
    ObservationStatus status = ObservationStatus::Nonfinite;
    std::string message;
    std::optional<double> residual;  // the declared form, unweighted
    std::optional<core::Convergence> convergence;
    std::optional<core::Convergence> first_numerical_failure;
    // Pressure closure: the branch the closure resolved (Unspecified when one stable root exists) and the count.
    std::optional<DensityBranch> branch;
    std::optional<std::size_t> stable_root_count;
    std::optional<bool> requested_tolerance_met;
    std::optional<equilibrium::EquilibriumResult> solved;  // the accepted central solve
    std::optional<equilibrium::CompiledSystem> compiled;    // the problem it solved, reused by differentiate
    std::optional<equilibrium::EquilibriumResult> final_solved;
    std::optional<equilibrium::CompiledSystem> final_compiled;
    std::shared_ptr<equilibrium::SaturationSearch> saturation;
    std::optional<double> saturation_loading, endpoint_loading, added_amount, enthalpy_difference;
    std::string cap_branch;
    std::vector<double> saturation_loading_gradient, endpoint_loading_gradient;
};

// d(value)/d(theta) in EOS-row order, then the requested reaction-coordinate order;
// value units per native coefficient unit, or empty with the blocking action status.
struct PredictionGradient final {
    std::vector<double> value;
    ObservationStatus status = ObservationStatus::ActionUnavailable;
    std::string message;
    std::vector<ElectrolyteInputAction> input_actions;
    std::vector<ResidualEnergyAction> parameter_actions{};
    std::vector<equilibrium::SolvedStateActionDiagnostics> state_diagnostics;
    std::optional<core::Convergence> first_numerical_failure;
    std::vector<double> saturation_loading_gradient, endpoint_loading_gradient;
};
struct ObservationInputDirection final {
    equilibrium::SolvedStateActionDirection state;
    double molality = 0.0;
};

// Admission by name: finite inputs, phase/species indices inside the declared problem.
using equilibrium::same_reaction;
bool same_state(const Observation& left, const Observation& right);
void validate(const Observation& observation, const HelmholtzModel& model,
              const equilibrium::CompiledSystem* prepared = nullptr);
void project_prediction(const HelmholtzModel& model, const Observation& observation, Prediction& prediction);
using ReferenceCache = std::vector<std::pair<Observation, InfiniteDilutionReference>>;
Prediction predict(
    const HelmholtzModel& model, const Observation& observation,
    const std::optional<equilibrium::EquilibriumResult>& start = std::nullopt,
    std::optional<equilibrium::CompiledSystem> prepared = std::nullopt, ReferenceCache* references = nullptr,
    const Prediction* previous = nullptr, const HelmholtzModel* previous_model = nullptr,
    bool same_parameters = false,
    std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max()
);
PredictionGradient differentiate(
    const HelmholtzModel& model, const Observation& observation, const Prediction& prediction
);
std::vector<PredictionGradient> differentiate_group(
    const HelmholtzModel& model, const std::vector<const Observation*>& observations, const Prediction& prediction,
    const CoordinateSet& coordinates = {},
    const std::vector<ObservationInputDirection>& inputs = {}
);
// r = (y - y_obs) / scale, or ln(y / y_obs) for positive y and y_obs; dr/dy = 1/scale or 1/y.
std::optional<double> residual(const Observation& observation, double value, ObservationStatus& status);
double residual_slope(const Observation& observation, double value);

}  // namespace epcsaft::regression
