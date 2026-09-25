#include "regression/coordinates.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace epcsaft::regression {
namespace {

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("coordinate ") + what);
}

using detail::family_values;

}  // namespace

std::vector<ActiveParameter> CoordinateSet::active_parameters() const {
    std::vector<ActiveParameter> result;
    for (const auto& coordinate : coordinates)
        for (const auto& slot : coordinate.slots) result.push_back({coordinate.family, slot});
    return result;
}

std::vector<std::size_t> CoordinateSet::tie() const {
    std::vector<std::size_t> result;
    for (std::size_t column = 0; column < coordinates.size(); ++column)
        result.insert(result.end(), coordinates[column].slots.size(), column);
    return result;
}

std::vector<double> CoordinateSet::start(const Parameters& parameters) const {
    std::vector<double> result;
    for (const auto& coordinate : coordinates) {
        const auto& values = family_values(parameters, coordinate.family);
        require(!coordinate.slots.empty() && !coordinate.slots.front().empty(), "slots");
        const double first = values.at(coordinate.slots.front().front());
        for (const auto& slot : coordinate.slots)
            for (const auto index : slot) require(values.at(index) == first, "tied values disagree");
        result.push_back((first - coordinate.origin) / coordinate.scale);
    }
    return result;
}

std::vector<double> CoordinateSet::physical(const std::vector<double>& optimizer) const {
    require(optimizer.size() == coordinates.size(), "optimizer width");
    std::vector<double> result;
    for (std::size_t column = 0; column < coordinates.size(); ++column) {
        const auto& coordinate = coordinates[column];
        result.insert(result.end(), coordinate.slots.size(),
                      std::clamp(coordinate.origin + coordinate.scale * optimizer[column], coordinate.lower, coordinate.upper));
    }
    return result;
}

// Checked in optimizer space against optimizer_bounds, the values the optimizer is handed: rounded
// (theta - origin) / scale is monotone in theta, so a physical value inside its bounds is inside these.
bool CoordinateSet::within_bounds(const std::vector<double>& optimizer) const {
    if (optimizer.size() != coordinates.size()) return false;
    const auto [lower, upper] = optimizer_bounds();
    for (std::size_t column = 0; column < coordinates.size(); ++column)
        if (!std::isfinite(optimizer[column]) || optimizer[column] < lower[column] || optimizer[column] > upper[column]) return false;
    return true;
}

std::pair<std::vector<double>, std::vector<double>> CoordinateSet::optimizer_bounds() const {
    std::vector<double> lower, upper;
    for (const auto& c : coordinates) {
        lower.push_back((c.lower - c.origin) / c.scale);
        upper.push_back((c.upper - c.origin) / c.scale);
    }
    return {lower, upper};
}

Parameters CoordinateSet::apply(Parameters parameters, const std::vector<double>& optimizer) const {
    require(optimizer.size() == coordinates.size(), "optimizer width");
    for (std::size_t column = 0; column < coordinates.size(); ++column) {
        const auto& coordinate = coordinates[column];
        auto& values = family_values(parameters, coordinate.family);
        // Clamped: at a bound the affine map can round one ulp past it.
        const double value = std::clamp(coordinate.origin + coordinate.scale * optimizer[column], coordinate.lower, coordinate.upper);
        for (const auto& slot : coordinate.slots)
            for (const auto index : slot) values.at(index) = value;
    }
    return parameters;
}

void validate(const CoordinateSet& set, const Parameters& parameters) {
    require(!set.coordinates.empty(), "set is empty");
    for (const auto& coordinate : set.coordinates) {
        require(!coordinate.slots.empty(), "slots");
        require(std::isfinite(coordinate.scale) && coordinate.scale > 0.0, "scale");
        require(std::isfinite(coordinate.origin), "origin");
        require(!std::isnan(coordinate.lower) && !std::isnan(coordinate.upper)
                && coordinate.lower <= coordinate.upper, "bounds");
        const auto& values = family_values(parameters, coordinate.family);
        for (const auto& slot : coordinate.slots) {
            require(!slot.empty(), "slot is empty");
            for (const auto index : slot) require(index < values.size(), "slot index");
        }
    }
    require(set.within_bounds(set.start(parameters)), "start outside bounds");
    (void)HelmholtzModel(parameters, set.active_parameters());
}

}  // namespace epcsaft::regression
