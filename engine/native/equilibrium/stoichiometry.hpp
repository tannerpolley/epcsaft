#pragma once

#include "equilibrium/problem.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft::equilibrium {

std::vector<double> resolve_feed(
    const HelmholtzModel& model, const FeedSpec& feed
);

std::vector<std::vector<double>> balance_basis(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

std::vector<double> standard_state(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied
);

std::vector<double> matvec(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector
);

double vector_norm(const std::vector<double>& values);

double feed_total(const std::vector<double>& feed);

std::vector<double> normalized_feed(const std::vector<double>& feed);

std::vector<double> reaction_affinities(
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& standard_gibbs
);

std::vector<std::vector<double>> reaction_projector(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

std::size_t reaction_rank(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

void validate_reaction_set(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

std::string reaction_name(const ReactionSpec& reaction, std::size_t index);

struct StoichiometryReport final {
    std::size_t components = 0;
    std::size_t reactions = 0;
    std::size_t rank = 0;
    bool independent = true;
    bool finite = true;
};

StoichiometryReport inspect_stoichiometry(
    std::size_t components, const std::vector<ReactionSpec>& reactions
);

std::vector<double> reaction_residuals(
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& standard_gibbs,
    const std::vector<double>& requested_affinities
);

std::vector<double> balance_residual(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& amounts,
    const std::vector<double>& target
);

std::vector<std::vector<double>> gram_matrix(
    const std::vector<std::vector<double>>& basis
);

std::vector<double> transpose_matvec(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector
);

bool is_orthonormal(
    const std::vector<std::vector<double>>& basis, double tolerance
);

std::vector<double> reaction_affinity_targets(
    double temperature, const std::vector<ReactionSpec>& reactions
);

std::vector<double> project_to_balance(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& composition
);

std::vector<double> reconstruct_from_balance(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& coordinates
);

double l1_norm(const std::vector<double>& values);

bool all_finite(const std::vector<double>& values);

double feed_mass_kg(
    const HelmholtzModel& model, const std::vector<double>& amounts
);

std::vector<double> scale_amounts(
    const std::vector<double>& amounts, double total
);

bool is_composition(const std::vector<double>& values, double tolerance);

double projection_error(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& composition
);

double maximum_abs(const std::vector<double>& values);

double basis_trace(const std::vector<std::vector<double>>& basis);

bool same_subspace(
    const std::vector<std::vector<double>>& left,
    const std::vector<std::vector<double>>& right,
    double tolerance
);

double balance_total(const std::vector<double>& amounts);

bool basis_matches_components(
    const std::vector<std::vector<double>>& basis, std::size_t components
);

bool contains_reaction_name(
    const std::vector<ReactionSpec>& reactions, const std::string& name
);

}  // namespace epcsaft::equilibrium

