#include "regression/coordinates.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace epcsaft::regression {
namespace {

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("coordinate ") + what);
}

template <class ParameterSet>
auto& family_values(ParameterSet& parameters, ParameterFamily family) {
    switch (family) {
    case ParameterFamily::SegmentCount: return parameters.segment_count;
    case ParameterFamily::SigmaConstant: return parameters.sigma_constant;
    case ParameterFamily::DispersionEnergyOverK: return parameters.dispersion_energy_over_k;
    case ParameterFamily::Kij: return parameters.k_ij;
    case ParameterFamily::Lij: return parameters.l_ij;
    case ParameterFamily::AssociationEnergyOverK: return parameters.association_energy_over_k;
    case ParameterFamily::AssociationVolume: return parameters.association_volume;
    }
    throw std::invalid_argument("coordinate family");
}

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
                      coordinate.origin + coordinate.scale * optimizer[column]);
    }
    return result;
}

bool CoordinateSet::within_bounds(const std::vector<double>& optimizer) const {
    if (optimizer.size() != coordinates.size()) return false;
    for (std::size_t column = 0; column < coordinates.size(); ++column) {
        const auto& coordinate = coordinates[column];
        const double value = coordinate.origin + coordinate.scale * optimizer[column];
        if (!std::isfinite(value) || value < coordinate.lower || value > coordinate.upper) return false;
    }
    return true;
}

Parameters CoordinateSet::apply(Parameters parameters, const std::vector<double>& optimizer) const {
    require(optimizer.size() == coordinates.size(), "optimizer width");
    for (std::size_t column = 0; column < coordinates.size(); ++column) {
        const auto& coordinate = coordinates[column];
        auto& values = family_values(parameters, coordinate.family);
        const double value = coordinate.origin + coordinate.scale * optimizer[column];
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
