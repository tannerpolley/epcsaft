#pragma once
#include "core/export.hpp"

#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <array>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft::regression {

struct Observation;
enum class ReactionCoefficient { A, B, C, D, EnthalpyShift };
struct ReactionCoordinate final {
    std::string name;
    ReactionCoefficient coefficient = ReactionCoefficient::A;
    double pivot_temperature = 0.0;  // K, required for EnthalpyShift
};
EPC_API std::array<double, 4> reaction_direction(const ReactionCoordinate& coordinate);
// -log10 K = A/T + B + C*T; refuses a nonzero logarithmic coefficient.
EPC_API std::array<double, 3> reaction_log10_reciprocal(const ReactionLogPolynomial& law);

struct CoordinateTarget final {
    ParameterFamily family;
    std::vector<std::size_t> indices;
    double coefficient = 1.0;
};

// Base value q = origin + scale * v, with inclusive bounds on q.
// EOS targets are coefficient*q; reaction coordinates have no EOS targets.
struct Coordinate final {
    std::vector<CoordinateTarget> targets;
    double origin = 0.0;
    double scale = 1.0;
    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();
    std::optional<ReactionCoordinate> reaction;
};

struct CoordinateSet final {
    std::vector<Coordinate> coordinates;

    // Physical rows: EOS slots in coordinate order, then reaction coordinates in order.
    EPC_API // active_parameters() contains only the EOS rows; tie() maps every row to optimizer order.
    EPC_API std::vector<ActiveParameter> active_parameters() const;
    EPC_API std::vector<std::size_t> tie() const;
    // Optimizer values of the parameters' current physical values; tied slots must agree.
    EPC_API std::vector<double> start(const Parameters& parameters, const std::vector<Observation>& observations = {},
                              const std::optional<ThermochemistryRecord>& thermochemistry = {}) const;
    // Mapped target value per physical row from v.
    EPC_API std::vector<double> physical(const std::vector<double>& optimizer) const;
    EPC_API bool within_bounds(const std::vector<double>& optimizer) const;
    // Optimizer-space bounds (lower, upper): (bound - origin) / scale.
    EPC_API std::pair<std::vector<double>, std::vector<double>> optimizer_bounds() const;
    // A copy of `parameters` with every slot (both mirrored entries of a pair) rewritten.
    EPC_API Parameters apply(Parameters parameters, const std::vector<double>& optimizer) const;
    EPC_API std::map<std::string, ReactionLogPolynomial> apply(std::vector<Observation>& observations, const std::vector<double>& optimizer) const;
};

// Reject invalid charts, conflicting named reaction declarations and constants fixed by other inputs.
EPC_API void validate(const CoordinateSet& set, const Parameters& parameters, const std::vector<Observation>& observations = {},
              const std::optional<ThermochemistryRecord>& thermochemistry = {});

}  // namespace epcsaft::regression
