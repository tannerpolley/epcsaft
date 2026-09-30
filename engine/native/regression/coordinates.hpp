#pragma once

#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft::regression {

struct Observation;
enum class ReactionCoefficient { A, B, C, D };
struct ReactionCoordinate final {
    std::string name;
    ReactionCoefficient coefficient = ReactionCoefficient::A;
};

// An affine optimizer coordinate theta = origin + scale * u, with inclusive native-unit
// bounds: either tied EOS slots of one family, or one named reaction coefficient.
struct Coordinate final {
    ParameterFamily family = ParameterFamily::SegmentCount;
    std::vector<std::vector<std::size_t>> slots;
    double origin = 0.0;
    double scale = 1.0;
    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();
    std::optional<ReactionCoordinate> reaction;
};

struct CoordinateSet final {
    std::vector<Coordinate> coordinates;

    // Physical rows: EOS slots in coordinate order, then reaction coordinates in order.
    // active_parameters() contains only the EOS rows; tie() maps every row to optimizer order.
    std::vector<ActiveParameter> active_parameters() const;
    std::vector<std::size_t> tie() const;
    // Optimizer values of the parameters' current physical values; tied slots must agree.
    std::vector<double> start(const Parameters& parameters, const std::vector<Observation>& observations = {},
                              const std::optional<ThermochemistryRecord>& thermochemistry = {}) const;
    // theta per physical row from u.
    std::vector<double> physical(const std::vector<double>& optimizer) const;
    bool within_bounds(const std::vector<double>& optimizer) const;
    // Optimizer-space bounds (lower, upper): (bound - origin) / scale.
    std::pair<std::vector<double>, std::vector<double>> optimizer_bounds() const;
    // A copy of `parameters` with every slot (both mirrored entries of a pair) rewritten.
    Parameters apply(Parameters parameters, const std::vector<double>& optimizer) const;
    void apply(std::vector<Observation>& observations, const std::vector<double>& optimizer) const;
};

// Reject invalid charts, conflicting named reaction declarations and constants fixed by other inputs.
void validate(const CoordinateSet& set, const Parameters& parameters, const std::vector<Observation>& observations = {},
              const std::optional<ThermochemistryRecord>& thermochemistry = {});

}  // namespace epcsaft::regression
