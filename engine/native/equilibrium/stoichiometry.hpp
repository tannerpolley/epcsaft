#pragma once

#include "equilibrium/problem.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::equilibrium {

bool exact_zero_dot(const std::vector<double>& left, const std::vector<double>& right);
std::vector<double> reconstruct_conserved_functional(std::vector<double> candidate,
    const std::vector<ReactionSpec>& reactions);
void restrict_inventory_balances(CompiledSystem& compiled, const std::vector<bool>& active);

std::vector<double> resolve_feed(
    const HelmholtzModel& model, const FeedSpec& feed
);

std::vector<double> resolve_feed(const IdealMoleFraction& model, const FeedSpec& feed);

std::vector<std::vector<double>> balance_basis(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

std::vector<double> standard_state(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied, const ThermochemistryRecord* thermochemistry = nullptr
);

// series[k][i] is d^k g_i / dT^k of the minimum-norm standard state (one gauge for every
// order) or of a supplied constant vector; order <= 3. A supplied vector must satisfy the
// reaction targets at this temperature by the inferred-vector criterion.
std::vector<std::vector<double>> standard_state_series(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied, std::size_t order, const ThermochemistryRecord* thermochemistry = nullptr
);

// series[k][r] is d^k t_r / dT^k of the reaction target t_r = nu_r . g in J/mol:
// affinity, -R T ln K, or -R T ln K(T) from the P3 correlation on the EOS density basis.
std::vector<std::vector<double>> reaction_target_series(
    double temperature, const std::vector<ReactionSpec>& reactions, std::size_t order,
    const ThermochemistryRecord* thermochemistry = nullptr, const std::string& standard_state_id = EOS_STANDARD_STATE_ID
);

// A constant supplied Gibbs vector admits no temperature-dependent reaction target: the first
// reaction whose |T^k t^(k)|, k = 1..3, exceeds sqrt(eps) max(1, |t|), if any.
std::optional<std::size_t> thermal_reference_missing(
    double temperature, const std::vector<ReactionSpec>& reactions, const std::vector<double>& supplied, const ThermochemistryRecord* thermochemistry = nullptr
);

// True when the declared reaction space conserves total moles: rank([B; e]) == rank(B).
bool conserves_total_moles(const std::vector<std::vector<double>>& basis);

// participant[i]: species i has a nonzero coefficient in a declared reaction. With the fed species
// it forms ALG-STABILITY's U.
std::vector<bool> reaction_participants(std::size_t components, const std::vector<ReactionSpec>& reactions);

// ALG-EQUILIBRIUM phase support: an omitted support is the candidates, zero-charge only in a vapor;
// a declared support is sorted and refused (core::DeclarationError) if empty, repeated, out of range,
// or, in a vapor, holding a charged component outside DensityTemperature permittivity.
std::vector<std::size_t> resolve_support(const PhaseSpec& phase, const std::vector<std::size_t>& candidates,
    const std::vector<int>& z, const Parameters* parameters = nullptr);

std::vector<double> matvec(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector
);

double sum(const std::vector<double>& values);

std::size_t reaction_rank(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

void validate_reaction_set(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

std::string reaction_name(const ReactionSpec& reaction, std::size_t index);

bool all_numeric(const std::vector<double>& values);

}  // namespace epcsaft::equilibrium
