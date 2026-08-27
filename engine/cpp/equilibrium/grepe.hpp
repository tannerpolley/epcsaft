#pragma once

#include "chemical_equilibrium.hpp"
#include "phase_equilibrium.hpp"
#include "provider.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_engine::equilibrium {

struct GrepeIntensive {
    bool solved = false;
    double value = 0.0;
    std::array<double, 2> bounds{};
    std::vector<double> starts;
};

struct GrepePhase {
    std::string identity;
    PhaseRole role = PhaseRole::Liquid;
    bool finite = true;
    bool ideal_gas = false;
    std::string model_reference_id;
    std::string support_name;
    std::vector<std::size_t> support_indices;
    std::array<double, 2> packing_fraction_interval{};
    std::optional<ReactiveFinitePhaseStart> start;
};

struct GrepeDirectionSystem {
    std::vector<std::string> phase_ids;
    std::vector<std::string> species_ids;
    std::vector<double> charges;
    DenseMatrix balance_matrix;
    DenseMatrix direction_matrix;
    std::vector<double> driving_constants;
};

struct GrepeSourceStandardState {
    std::string id;
    std::vector<double> log_activity_scale_factors;
    double reference_pressure_pa = 0.0;
    SourceReferenceRequest reference;
};

struct GrepeFeedCoordinate {
    std::string identity;
    std::array<double, 2> bounds{};
};

struct GrepeFeed {
    bool affine = false;
    std::vector<double> amounts_or_base_mol;
    DenseMatrix component_coordinate_matrix;
    std::vector<std::string> coordinate_ids;
    std::vector<GrepeFeedCoordinate> coordinates;
};

struct GrepeIntensiveBoundary {
    std::string identity;
    std::string coordinate_identity;
    std::string phase_identity;
    std::string selector;
    std::vector<double> coefficients;
    double target = 0.0;
    double initial = 0.0;
};

struct GrepeContinuationPhase {
    std::string identity;
    std::vector<double> mole_fractions;
    double molar_volume_m3_per_mol = 0.0;
};

struct GrepeRequest {
    std::string identity;
    GrepeIntensive temperature;
    GrepeIntensive pressure;
    std::vector<GrepePhase> phases;
    ReactionSystemInput reaction;
    double trace_floor = 0.0;
    std::vector<std::string> reaction_phase_ids;
    std::vector<ReactiveFinitePhaseOutput> outputs;
    std::optional<GrepeDirectionSystem> directions;
    std::optional<GrepeSourceStandardState> source;
    std::optional<GrepeFeed> feed;
    std::optional<GrepeIntensiveBoundary> boundary;
    std::vector<GrepeContinuationPhase> continuation_phases;
    double continuation_pressure_pa = 0.0;
    bool perform_derivative_audits = false;
    bool declared_starts_only = false;
};

struct GrepeResult {
    std::string formulation;
    std::string certificate_scope;
    ReactiveFinitePhaseResult equilibrium;
    std::string failure_code;
    double temperature_k = std::numeric_limits<double>::quiet_NaN();
    double pressure_pa = std::numeric_limits<double>::quiet_NaN();
    double source_reference_representation_residual_inf_norm = 0.0;
    std::vector<std::pair<std::string, double>> evidence;
    std::string branch_kind;
    std::vector<double> branch_coordinates;
    std::vector<std::string> phase_mechanical_classes;
};

[[nodiscard]] GrepeResult solve_grepe(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    const GrepeRequest& request
);

[[nodiscard]] GrepeResult evaluate_grepe_local(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
);

[[nodiscard]] GrepeResult finalize_grepe_local(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
);

[[nodiscard]] GrepeResult initialize_grepe_local(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
);

}  // namespace epcsaft_engine::equilibrium
