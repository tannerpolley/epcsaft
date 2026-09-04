#pragma once

#include "core/convergence.hpp"
#include "model/helmholtz_model.hpp"

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <array>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::equilibrium {

enum class PhaseKind { Liquid, Vapor };
enum class AmountKind { Free, Pinned };
enum class FeedKind { MoleFractions, Amounts, MassFractions, Molalities };

inline const char* phase_kind_name(PhaseKind kind) noexcept {
    return kind == PhaseKind::Liquid ? "liquid" : "vapor";
}

inline const char* amount_kind_name(AmountKind kind) noexcept {
    return kind == AmountKind::Free ? "free" : "pinned";
}

inline const char* feed_kind_name(FeedKind kind) noexcept {
    switch (kind) {
    case FeedKind::MoleFractions: return "mole_fractions";
    case FeedKind::Amounts: return "amounts";
    case FeedKind::MassFractions: return "mass_fractions";
    case FeedKind::Molalities: return "molalities";
    }
    return "unknown";
}

struct AmountSpec final {
    AmountKind kind = AmountKind::Free;
    double fraction = 0.0;
    static AmountSpec free() { return {}; }
    static AmountSpec pinned(double value) { return {AmountKind::Pinned, value}; }
};

struct FreeValue final {
    bool free = false;
    double value = 0.0;
    bool absent = false;
    static FreeValue fixed(double value) { return {false, value, false}; }
    static FreeValue variable() { return {true, 0.0, true}; }
    static FreeValue variable(double guess) { return {true, guess, false}; }

    double seed_or(double fallback) const noexcept {
        return free && absent ? fallback : value;
    }

    bool is_positive() const noexcept {
        return !absent && std::isfinite(value) && value > 0.0;
    }
};

struct PhaseSpec final {
    std::string name;
    PhaseKind kind = PhaseKind::Liquid;
    std::vector<std::size_t> support;
    AmountSpec amount;
    std::vector<double> composition_guess;
    std::optional<double> packing_guess;

    bool is_pinned() const noexcept { return amount.kind == AmountKind::Pinned; }
    bool is_free() const noexcept { return amount.kind == AmountKind::Free; }
};

struct FeedSpec final {
    FeedKind kind = FeedKind::Amounts;
    std::vector<double> values;
    std::vector<std::size_t> solvent;
    double solvent_mass_kg = 1.0;

    bool has_solvent() const noexcept { return !solvent.empty(); }
};

struct ReactionSpec final {
    std::string name;
    std::vector<double> stoichiometry;
    double equilibrium_constant = 1.0;
    std::optional<double> affinity_j_per_mol;
};

struct EquilibriumProblem final {
    std::vector<PhaseSpec> phases;
    FreeValue temperature;
    FreeValue pressure;
    FeedSpec feed;
    std::vector<ReactionSpec> reactions;
    std::vector<double> standard_gibbs_j_per_mol;
    bool topology_declared = true;

    std::size_t component_count() const noexcept {
        return feed.values.size();
    }
};

inline std::vector<int> species_z(const Parameters& parameters) {
    return parameters.charge_number;
}

enum class RowKind { Chemical, Pressure, Moiety, Electro, Normalization };

inline const char* row_kind_name(RowKind kind) noexcept {
    switch (kind) {
    case RowKind::Chemical: return "chemical_potential_stationarity";
    case RowKind::Pressure: return "pressure_equality";
    case RowKind::Moiety: return "material_balance";
    case RowKind::Electro: return "phase_electroneutrality";
    case RowKind::Normalization: return "phase_normalization";
    }
    return "unknown";
}

struct Row final {
    RowKind kind = RowKind::Chemical;
    std::size_t phase = 0;
    std::size_t index = 0;
    const char* eqid = "chemical_potential_stationarity";

    bool is_intensive() const noexcept {
        return kind == RowKind::Chemical || kind == RowKind::Pressure;
    }
};

struct PointEvaluation final {
    std::vector<double> values;
    std::vector<double> jacobian;
    double objective = 0.0;
    double floor = 0.0;

    double residual_norm() const noexcept {
        double result = 0.0;
        for (const double value : values) result = std::max(result, std::abs(value));
        return result;
    }
};

struct CompiledSystem final {
    const HelmholtzModel* model = nullptr;
    EquilibriumProblem problem;
    std::vector<double> n_feed;
    std::vector<double> n_feed_normalised;
    std::vector<double> feed_composition;
    std::vector<double> balance;
    std::vector<double> standard_gibbs;
    std::vector<std::vector<double>> B;
    std::vector<Row> rows;
    std::size_t pinned_count = 0;
    std::size_t free_intensive_count = 0;
    std::size_t rank = 0;
    std::size_t charge_row_count = 0;
    std::size_t primal_size = 0;
    long rows_minus_unknowns = 0;
    std::size_t generic_nullity = 0;
    bool feasible = false;
    std::vector<double> F_y;
    std::vector<std::size_t> removed_species;
    std::vector<std::string> row_names;
    double feed_total = 0.0;
    bool has_ions = false;
    bool has_reaction = false;
    std::size_t jacobian_rank = 0;

    PointEvaluation evaluate(const std::vector<double>& point) const;
    std::vector<double> generic_point() const;
    std::vector<double> jacobian(const std::vector<double>& point) const;
    std::vector<double> lever_jacobian(const std::vector<double>& point) const;
    std::vector<double> objective_gradient(const std::vector<double>& point) const;
    std::vector<double> hessian(
        const std::vector<double>& point,
        const std::vector<double>& constraint_multipliers = {}
    ) const;

    std::size_t unknown_count() const noexcept {
        return primal_size + rank + charge_row_count;
    }

    std::size_t constraint_count() const noexcept { return rows.size(); }

    bool is_square() const noexcept { return rows_minus_unknowns == 0; }
};

CompiledSystem compile_problem(
    const HelmholtzModel& model, const EquilibriumProblem& problem
);

struct EquilibriumResult final {
    bool success = false;
    std::string message;
    std::vector<double> temperature;
    std::vector<double> pressure;
    std::vector<double> amounts;
    std::vector<double> mole_fractions;
    std::vector<double> molar_densities;
    std::vector<double> packing_fractions;
    std::vector<double> multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<double> residuals;
    std::vector<double> sensitivity;
    std::vector<double> coordinates;
    std::size_t phase_count = 0;
    std::size_t iterations = 0;
    std::size_t evaluations = 0;
    std::size_t hessian_calls = 0;
    double objective = 0.0;
    double reduced_hessian_minimum = 0.0;
    core::Convergence convergence;
    std::vector<Row> rows;

    bool has_state() const noexcept {
        return success && !temperature.empty() && !pressure.empty()
            && !amounts.empty() && !mole_fractions.empty();
    }

    double residual_norm() const noexcept {
        double result = 0.0;
        for (const double value : residuals) result = std::max(result, std::abs(value));
        return result;
    }
};

EquilibriumResult solve_equilibrium(
    const HelmholtzModel& model,
    const EquilibriumProblem& problem,
    const std::optional<EquilibriumResult>& start = std::nullopt
);

}  // namespace epcsaft::equilibrium
