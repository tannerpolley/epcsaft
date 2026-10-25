#pragma once
// Includes and helpers shared by the equilibrium_checks translation units.
#include "quad_base.hpp"
#include "properties/helmholtz_tape.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "properties/thermochemistry.hpp"
#include "equilibrium_support.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "recording_support.hpp"
#include "binary128_support.hpp"
#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <memory>
#include <limits>
#include <string>
#include <stdexcept>
#include <tuple>
#include <map>
#include <numeric>

namespace equilibrium_checks {

// EqID reaction_activity_quotient: compare with fresh raw rows, then remove every stored diagnostic.
template<class Model>
void check_activity_defect(const Model& model, const epcsaft::equilibrium::CompiledSystem& compiled,
    epcsaft::equilibrium::EquilibriumResult candidate, double log_k) {
    using namespace epcsaft::equilibrium;
    candidate.success = false; candidate.requested_tolerance_met = false;
    const auto raw = compiled.evaluate(candidate.coordinates).values;
    const auto& nu = candidate.declaration->reactions[0].stoichiometry;
    double defect = 0;
    for (std::size_t j = 0; j < raw.size(); ++j) {
        const auto& row = compiled.rows[j];
        if (row.kind == RowKind::Chemical && row.phase == 0) defect += nu[row.index] * raw[j];
        if (row.kind == RowKind::ReactionMassAction && row.index == 0) defect = raw[j];
    }
    const auto read = reaction_activities(model, candidate)[0];
    CHECK(read.log_quotient && std::abs(*read.log_quotient - log_k) > 1e-4);
    CHECK(std::abs(*read.log_quotient - log_k - defect) <= 1e-10);
    candidate.residuals.assign(candidate.residuals.size(), NAN); candidate.classified_residuals = {INFINITY};
    candidate.multipliers = {NAN}; candidate.constraint_multipliers = {INFINITY};
    CHECK(reaction_activities(model, candidate)[0].species.log_activities == read.species.log_activities);
    auto declaration = std::make_shared<EquilibriumProblem>(*candidate.declaration);
    if (declaration->reactions[0].correlation) declaration->reactions[0].correlation->a += 0.7;
    else declaration->reactions[0].equilibrium_constant = 2.0;
    candidate.declaration = declaration;
    const auto changed = reaction_activities(model, candidate)[0];
    CHECK(changed.basis == read.basis && changed.log_quotient == read.log_quotient && changed.species.activities == read.species.activities);
    std::printf("activity-defect,lnQ=%.17g,lnK=%.17g,raw=%.17g,error=%.17g,poison=0,target=0\n",
        *read.log_quotient, log_k, defect, *read.log_quotient - log_k - defect);
}

inline double convergence_limit(const epcsaft::core::Convergence& convergence) {
    CHECK(convergence.criterion != nullptr);
    return std::max(convergence.criterion->value, convergence.floor.value_or(0.0));
}

inline bool close(double left, double right) {
    return std::abs(left - right) <= 1.0e-6 * std::max(1.0, std::abs(right));
}

inline bool close_with_tolerance(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

inline double centered_step(double value) {
    return std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(std::abs(value), 1.0);
}

// Old normalized correlation ln K = a + 1200/T + 0.4 ln T + 0.002 T; P3 carries the
// implicit 1 K reference through its legacy normalization at 300 K on [250, 450] K.
// The retained legacy offset a = 0.3 gives K ~ 870-1010 (dilute A, x_A ~ 1e-3); the
// owner-selected qualification offset makes K(350 K) = 2 with unchanged b, c, d.
constexpr double DILUTE_ISOMERISATION_A = 0.3;
inline double moderate_isomerisation_a() {
    return std::log(2.0) - (1200.0 / 350.0 + 0.4 * std::log(350.0) + 0.002 * 350.0);
}
inline epcsaft::ReactionLogPolynomial legacy_isomerisation(double sign, double a) {
    return {sign * a + sign * 0.4 * std::log(300.0), sign * 1200.0, sign * 0.4, sign * 0.002, 300.0, 250.0, 450.0,
        epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
}

inline void check_finite_difference(const epcsaft::equilibrium::CompiledSystem& compiled,
                             const std::vector<double>& point) {
    const auto analytic = compiled.jacobian(point);
    const auto base = compiled.evaluate(point).values;
    const auto context = equilibrium_support::weighted_hessian_context(compiled);
    CHECK(analytic.size() == base.size() * point.size());
    for (std::size_t column = 0; column < point.size(); ++column) {
        auto plus = point;
        auto minus = point;
        // Centred differences balance truncation and rounding at eps^(1/3), not eps^(1/2):
        // the rows cancel O(10) potentials, so the sqrt(eps) step carries ~5e-7 rounding noise.
        const double delta = std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(std::abs(point[column]), 1.0);
        plus[column] += delta;
        minus[column] -= delta;
        const auto up = compiled.evaluate(plus).values;
        const auto down = compiled.evaluate(minus).values;
        for (std::size_t row = 0; row < base.size(); ++row) {
            const double centered = (up[row] - down[row]) / (2.0 * delta);
            const double model_value = analytic[row * point.size() + column];
            const bool verdict = close(model_value, centered);
            const auto point_id = std::to_string(row) + "-" + std::to_string(column);
            recording_support::compare(
                (std::string("jacobian/") + point_id).c_str(),
                "equilibrium_jacobian", "1", model_value, centered,
                recording_support::fields({
                    context,
                    recording_support::field("point", point),
                    recording_support::field("row", row),
                    recording_support::field("column", column),
                    recording_support::field("delta", delta),
                }),
                "Relationship", "finite-difference", "equilibrium Jacobian matches centered residual derivative",
                "reference_or_one", 0.0, 1.0e-6, verdict,
                "f5-jacobian", point_id.c_str(), "J", "T,P,\\mathbf{x},\\Delta z",
                "J_{ij}=\\Delta r_i/\\Delta z_j", "row,column");
        }
    }
}

inline std::vector<std::size_t> callback_rows(
    const epcsaft::equilibrium::CompiledSystem& compiled
) {
    std::vector<std::size_t> result;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.feasibility_cell()
            || (compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Chemical
                && compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Pressure)) {
            result.push_back(row);
        }
    }
    return result;
}

inline std::vector<double> weighted_gradient(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers,
    double objective_factor
) {
    auto result = compiled.objective_gradient(point);
    const auto jacobian = compiled.jacobian(point);
    const auto rows = callback_rows(compiled);
    if (result.size() != point.size()
        || jacobian.size() != compiled.rows.size() * point.size()
        || multipliers.size() != rows.size()) {
        throw std::runtime_error("weighted-gradient shape mismatch");
    }
    for (double& value : result) value *= objective_factor;
    for (std::size_t multiplier = 0; multiplier < rows.size(); ++multiplier) {
        for (std::size_t column = 0; column < point.size(); ++column) {
            result[column] += multipliers[multiplier]
                * jacobian[rows[multiplier] * point.size() + column];
        }
    }
    return result;
}

// ALG-EQUILIBRIUM: the minimizer's callbacks measure each free amount in its start amount, floored at
// PHASE_DISAPPEARANCE_AMOUNT N_feed and fixed for the solve: physical = unit * callback coordinate (#243).
inline std::vector<double> callback_unit(
    const epcsaft::equilibrium::CompiledSystem& compiled, const std::vector<double>& start
) {
    std::vector<double> unit(start.size(), 1.0);
    if (compiled.feasibility_cell()) return unit;
    const auto offsets = compiled.phase_offsets();
    for (std::size_t p = 0; p < offsets.size(); ++p) {
        const auto& phase = compiled.problem.phases[p];
        if (!phase.is_free()) continue;
        const std::size_t i = offsets[p] + phase.support.size();
        unit[i] = std::max(start[i], epcsaft::policy::PHASE_DISAPPEARANCE_AMOUNT.value * compiled.feed_total);
    }
    return unit;
}

inline double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) result = std::max(result, std::abs(value));
    return result;
}

inline void check_weighted_hessian_gradient_ladder(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers,
    const char* label,
    double tolerance
) {
    const auto rows = callback_rows(compiled);
    epcsaft::equilibrium::GibbsMinimizer nlp(compiled, point);
    const auto unit = callback_unit(compiled, point);
    auto scaled = point;  // the same physical point in the callback coordinates
    for (std::size_t i = 0; i < scaled.size(); ++i) scaled[i] /= unit[i];
    Ipopt::Index variables = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
    Ipopt::TNLP::IndexStyleEnum index_style;
    CHECK(nlp.get_nlp_info(variables, constraints, jacobian_entries, hessian_entries, index_style));
    CHECK(variables >= 0 && constraints >= 0);
    CHECK(point.size() == static_cast<std::size_t>(variables));
    CHECK(multipliers.size() == rows.size());
    CHECK(multipliers.size() == static_cast<std::size_t>(constraints));

    const auto objective_factors = {0.0, 0.5, 1.0};
    const auto step_factors = {10.0, 3.0, 1.0};
    const auto context = equilibrium_support::weighted_hessian_context(compiled);

    std::vector<Ipopt::Index> hessian_rows(static_cast<std::size_t>(hessian_entries));
    std::vector<Ipopt::Index> hessian_columns(static_cast<std::size_t>(hessian_entries));
    CHECK(nlp.eval_h(variables, scaled.data(), true, 1.0, constraints,
                     multipliers.data(), true, hessian_entries,
                     hessian_rows.data(), hessian_columns.data(), nullptr));
    for (const double objective_factor : objective_factors) {
        std::vector<Ipopt::Number> hessian(
            static_cast<std::size_t>(hessian_entries), 0.0
        );
        CHECK(nlp.eval_h(variables, scaled.data(), false, objective_factor,
                         constraints, multipliers.data(), false,
                         hessian_entries, hessian_rows.data(),
                         hessian_columns.data(), hessian.data()));
        for (const double step_factor : step_factors) {
            for (std::size_t column = 0; column < point.size(); ++column) {
                auto plus = point;
                auto minus = point;
                const double step = step_factor * centered_step(point[column]);
                plus[column] += step;
                minus[column] -= step;
                const auto grad_plus = weighted_gradient(
                    compiled, plus, multipliers, objective_factor
                );
                const auto grad_minus = weighted_gradient(
                    compiled, minus, multipliers, objective_factor
                );
                for (std::size_t row = column; row < point.size(); ++row) {
                    const double centered = (grad_plus[row] - grad_minus[row]) / (2.0 * step);
                    const std::size_t index = row * (row + 1) / 2 + column;
                    const double restored = hessian[index] / (unit[row] * unit[column]);  // S H S back to H
                    const auto comparison_id = std::string("hessian/") + label
                        + "/objective-" + recording_support::number(objective_factor)
                        + "/step-" + recording_support::number(step_factor)
                        + "/" + std::to_string(row) + "-" + std::to_string(column);
                    const auto point_id = std::string("case-") + label
                        + "/objective-" + recording_support::number(objective_factor)
                        + "/step-" + recording_support::number(step_factor)
                        + "/row-" + std::to_string(row)
                        + "/column-" + std::to_string(column);
                    const auto coordinates = recording_support::fields({
                        recording_support::field("case", label),
                        recording_support::field("objective_factor", objective_factor),
                        recording_support::field("step_factor", step_factor),
                        recording_support::field("row", row),
                        recording_support::field("column", column),
                    });
                    const auto inputs = recording_support::fields({
                        context,
                        recording_support::field("case", label),
                        recording_support::field("objective_factor", objective_factor),
                        recording_support::field("step_factor", step_factor),
                        recording_support::field("point", point),
                        recording_support::field("multipliers", multipliers),
                        recording_support::field("row", row),
                        recording_support::field("column", column),
                        recording_support::field("step", step),
                        recording_support::field("variables", point.size()),
                        recording_support::field("constraints", multipliers.size()),
                    });
                    const bool pressure_entry = compiled.problem.pressure.free
                        && (row == compiled.primal_size - 1 || column == compiled.primal_size - 1);
                    const double coordinate_scale = pressure_entry ? compiled.problem.pressure.value : 1.0;
                    const double absolute = pressure_entry ? 1.0e-12 / coordinate_scale : 0.0;
                    const double allowance = pressure_entry ? absolute + tolerance * std::abs(centered)
                        : tolerance * std::max(1.0, std::abs(centered));
                    const bool verdict = std::abs(restored - centered) <= allowance;
                    recording_support::compare(
                        comparison_id.c_str(), "weighted_lagrangian_hessian", "coordinate",
                        restored, centered, inputs, "Relationship",
                        "weighted_hessian_gradient_finite_difference",
                        "nabla2_L=d(weighted_gradient)/d(coordinate)", pressure_entry ? "reference_scale" : "reference_or_one",
                        absolute, tolerance, verdict, "f5-hessian-ladder", point_id.c_str(),
                        "\\nabla^2\\mathcal L", "\\mathbf{z},\\lambda,\\Delta z",
                        "\\nabla^2\\mathcal L_{ij}=\\Delta(\\nabla\\mathcal L)_i/\\Delta z_j",
                        coordinates.c_str());
                }
            }
        }
    }
}
// E2 and E5 share the identical-ethane pair and its ideal-gas records.
namespace e2 {

inline epcsaft::Parameters identical_ethane_pair() {
    auto result = properties_support::parameters(2);
    const auto ethane = properties_support::ethane();
    for (std::size_t i = 0; i < 2; ++i) {
        result.segment_count[i] = ethane.segment_count[0];
        result.sigma_constant[i] = ethane.sigma_constant[0];
        result.dispersion_energy_over_k[i] = ethane.dispersion_energy_over_k[0];
        result.molar_mass_kg_per_mol[i] = ethane.molar_mass_kg_per_mol[0];
    }
    return result;
}

inline epcsaft::ThermochemistryRecord species_record(double cp_a, double cp_b) {
    epcsaft::ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    const double anchor_entropy = std::log(350.0 / 300.0);
    record.component_intervals = {
        {epcsaft::IdealInterval{250.0, 450.0, true, true, epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{cp_a}, 350.0, 1500.0, 30.0 * anchor_entropy}, 1.0e5}}},
        {epcsaft::IdealInterval{250.0, 450.0, true, true, epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{cp_b}, 350.0, 1000.0, 20.0 * anchor_entropy}, 1.0e5}}}};
    return record;
}

}  // namespace e2

}  // namespace equilibrium_checks
