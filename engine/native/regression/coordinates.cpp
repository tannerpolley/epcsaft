#include "regression/coordinates.hpp"
#include "regression/observations.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <algorithm>
#include <array>
#include <set>
#include <cmath>
#include <stdexcept>
#include <string>

namespace epcsaft::regression {
namespace {

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("coordinate ") + what);
}

double& coefficient(ReactionLogPolynomial& law, ReactionCoefficient k) {
    constexpr std::array members = {&ReactionLogPolynomial::a, &ReactionLogPolynomial::b,
                                   &ReactionLogPolynomial::c, &ReactionLogPolynomial::d};
    return law.*members.at(static_cast<std::size_t>(k));
}

using detail::family_values;

}  // namespace

std::vector<ActiveParameter> CoordinateSet::active_parameters() const {
    std::vector<ActiveParameter> result;
    for (const auto& coordinate : coordinates)
        if (!coordinate.reaction) for (const auto& slot : coordinate.slots) result.push_back({coordinate.family, slot});
    return result;
}

std::vector<std::size_t> CoordinateSet::tie() const {
    std::vector<std::size_t> result;
    for (std::size_t column = 0; column < coordinates.size(); ++column)
        result.insert(result.end(), coordinates[column].slots.size(), column);
    for (std::size_t column = 0; column < coordinates.size(); ++column)
        if (coordinates[column].reaction) result.push_back(column);
    return result;
}

std::vector<double> CoordinateSet::start(const Parameters& parameters, const std::vector<Observation>& observations,
                                         const std::optional<ThermochemistryRecord>& thermochemistry) const {
    std::vector<double> result;
    std::set<std::pair<std::string, ReactionCoefficient>> seen;
    for (const auto& coordinate : coordinates) {
        double first = 0.0;
        if (coordinate.reaction) {
            const auto& target = *coordinate.reaction;
            const auto check = [&](bool ok, const char* cause) {
                if (!ok) throw std::invalid_argument("reaction coordinate " + target.name + ": " + cause);
            };
            check(static_cast<std::size_t>(target.coefficient) < 4, "unknown coefficient");
            check(!target.name.empty() && seen.emplace(target.name, target.coefficient).second, "empty or duplicate target");
            const equilibrium::ReactionSpec* prior = nullptr;
            const equilibrium::EquilibriumProblem* prior_problem = nullptr;
            for (const auto& observation : observations) {
                if (pressure_closed(observation.kind)) continue;
                const auto& problem = observation.problem;
                std::size_t matches = 0;
                for (const auto& reaction : problem.reactions) {
                    if (reaction.name != target.name) continue;
                    check(++matches == 1, "repeated name in problem");
                    check(!reaction.affinity_j_per_mol, "affinity target");
                    check(reaction.correlation || reaction.equilibrium_constant, "missing ln K target");
                    check(problem.standard_gibbs_j_per_mol.empty(), "supplied standard Gibbs values");
                    check(equilibrium::reaction_rank(parameters.component_count(), problem.reactions) == problem.reactions.size(), "dependent reaction set");
                    check(!reaction.reference || !problem.temperature.free, "referenced reaction at free temperature");
                    check(reaction.correlation || target.coefficient == ReactionCoefficient::A, "scalar K admits only A");
                    check(reaction.reference || !thermochemistry || !reaction_thermochemistry_density(*thermochemistry,
                          reaction.stoichiometry, problem.temperature.value), "K fixed by thermochemistry");
                    if (prior) {
                        check(same_reaction(*prior, reaction), "tied declarations or starts disagree");
                        if (reaction.reference) {
                            check(problem.source_reference && prior_problem->source_reference, "missing source declaration");
                            const auto& x = *problem.source_reference;
                            const auto& y = *prior_problem->source_reference;
                            check(x.solvent_index == y.solvent_index && (reaction.reference->source_basis != ReferenceSourceBasis::CommonMolalityInfiniteDilution
                                  || x.standard_molality_mol_per_kg == y.standard_molality_mol_per_kg), "tied source bases disagree");
                        }
                    } else {
                        prior = &reaction;
                        prior_problem = &problem;
                        auto law = reaction.correlation.value_or(ReactionLogPolynomial{});
                        first = reaction.correlation ? coefficient(law, target.coefficient) : std::log(reaction.equilibrium_constant.value());
                    }
                }
            }
            check(prior, "unmatched name");
        } else {
            const auto& values = family_values(parameters, coordinate.family);
            require(!coordinate.slots.empty() && !coordinate.slots.front().empty(), "slots");
            first = values.at(coordinate.slots.front().front());
            for (const auto& slot : coordinate.slots)
                for (const auto index : slot) require(values.at(index) == first, "tied values disagree");
        }
        result.push_back((first - coordinate.origin) / coordinate.scale);
    }
    return result;
}

std::vector<double> CoordinateSet::physical(const std::vector<double>& optimizer) const {
    require(optimizer.size() == coordinates.size(), "optimizer width");
    std::vector<double> result;
    for (const auto column : tie()) {
        const auto& coordinate = coordinates[column];
        result.push_back(std::clamp(coordinate.origin + coordinate.scale * optimizer[column], coordinate.lower, coordinate.upper));
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
    const auto physical_values = physical(optimizer);
    const auto active = active_parameters();
    for (std::size_t row = 0; row < active.size(); ++row) {
        auto& values = family_values(parameters, active[row].family);
        for (const auto index : active[row].indices) values.at(index) = physical_values[row];
    }
    return parameters;
}

void CoordinateSet::apply(std::vector<Observation>& observations, const std::vector<double>& optimizer) const {
    const auto values = physical(optimizer);
    const auto order = tie();
    for (std::size_t row = 0; row < order.size(); ++row) {
        const auto& target = coordinates[order[row]].reaction;
        if (!target) continue;
        for (auto& observation : observations)
            for (auto& reaction : observation.problem.reactions) {
                if (reaction.name != target->name) continue;
                if (reaction.correlation) coefficient(*reaction.correlation, target->coefficient) = values[row];
                else reaction.equilibrium_constant = std::exp(values[row]);
            }
    }
}

void validate(const CoordinateSet& set, const Parameters& parameters, const std::vector<Observation>& observations,
              const std::optional<ThermochemistryRecord>& thermochemistry) {
    require(!set.coordinates.empty(), "set is empty");
    for (const auto& coordinate : set.coordinates) {
        require(coordinate.reaction ? coordinate.slots.empty() : !coordinate.slots.empty(), "slots");
        require(std::isfinite(coordinate.scale) && coordinate.scale > 0.0, "scale");
        require(std::isfinite(coordinate.origin), "origin");
        require(!std::isnan(coordinate.lower) && !std::isnan(coordinate.upper)
                && coordinate.lower <= coordinate.upper, "bounds");
        if (coordinate.reaction) continue;
        const auto& values = family_values(parameters, coordinate.family);
        for (const auto& slot : coordinate.slots) {
            require(!slot.empty(), "slot is empty");
            for (const auto index : slot) require(index < values.size(), "slot index");
        }
    }
    require(set.within_bounds(set.start(parameters, observations, thermochemistry)), "start outside bounds");
    (void)HelmholtzModel(parameters, set.active_parameters());
}

}  // namespace epcsaft::regression
