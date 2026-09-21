#pragma once

#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <limits>
#include <vector>

namespace epcsaft::regression {

// One optimizer coordinate u. Every slot it drives is one Properties active-parameter
// coordinate of the same family (a component index, or a reciprocal pair in the owner's
// index convention); tied slots share the physical value theta = origin + scale * u.
// Bounds are inclusive and in the family's native unit.
struct Coordinate final {
    ParameterFamily family = ParameterFamily::SegmentCount;
    std::vector<std::vector<std::size_t>> slots;
    double origin = 0.0;
    double scale = 1.0;
    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();
};

struct CoordinateSet final {
    std::vector<Coordinate> coordinates;

    // The ordered physical coordinate list handed to HelmholtzModel: physical row r is driven
    // by optimizer column tie()[r]. Every derivative with respect to the model's active
    // parameters is indexed by these rows, in this order.
    std::vector<ActiveParameter> active_parameters() const;
    std::vector<std::size_t> tie() const;
    // Optimizer values of the parameters' current physical values; tied slots must agree.
    std::vector<double> start(const Parameters& parameters) const;
    // theta per physical row from u.
    std::vector<double> physical(const std::vector<double>& optimizer) const;
    bool within_bounds(const std::vector<double>& optimizer) const;
    // A copy of `parameters` with every slot (both mirrored entries of a pair) rewritten.
    Parameters apply(Parameters parameters, const std::vector<double>& optimizer) const;
};

// Rejects by name: no coordinates, empty slots, nonfinite or nonpositive scale, lower above
// upper, a slot outside its family array, a start outside its bounds or disagreeing tied
// values; then constructs the model so Properties applies its own pair/association admission.
void validate(const CoordinateSet& set, const Parameters& parameters);

}  // namespace epcsaft::regression
