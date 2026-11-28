#pragma once

#include "core/bitwise_equal.hpp"
#include "equilibrium/problem.hpp"

#include <algorithm>
#include <set>
#include <utility>

namespace epcsaft::equilibrium {

struct PhaseAddition final {
    PhaseKind kind;
    std::vector<std::size_t> support;
    const std::vector<double>* seed;
    std::set<std::size_t> omitted_species;

    bool operator==(const PhaseAddition& other) const {
        return kind == other.kind && support == other.support && omitted_species == other.omitted_species
            && core::bitwise_equal(*seed, *other.seed);
    }
};

inline bool remember_phase_addition(std::vector<PhaseAddition>& attempted, PhaseAddition addition) {
    if (std::find(attempted.begin(), attempted.end(), addition) != attempted.end()) return false;
    attempted.push_back(std::move(addition));
    return true;
}

}  // namespace epcsaft::equilibrium
