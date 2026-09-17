#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "properties/thermochemistry.hpp"
#include "equilibrium_support.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <stdexcept>

namespace {

double convergence_limit(const epcsaft::core::Convergence& convergence) {
    CHECK(convergence.criterion != nullptr);
    return std::max(convergence.criterion->value, convergence.floor.value_or(0.0));
}

bool close(double left, double right) {
    return std::abs(left - right) <= 1.0e-6 * std::max(1.0, std::abs(right));
}

bool close_with_tolerance(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

double centered_step(double value) {
    return std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(std::abs(value), 1.0);
}

std::string weighted_hessian_context(
    const epcsaft::equilibrium::CompiledSystem& compiled
) {
    const auto& problem = compiled.problem;
    std::string phase_definition;
    for (const auto& phase : problem.phases) {
        if (!phase_definition.empty()) phase_definition += ',';
        phase_definition += phase.name + ":"
            + epcsaft::equilibrium::phase_kind_name(phase.kind) + ":"
            + epcsaft::equilibrium::amount_kind_name(phase.amount.kind) + ":"
            + recording_support::number(phase.amount.fraction)
            + ":support=" + std::to_string(phase.support.size()) + ":guess=";
        for (std::size_t index = 0; index < phase.composition_guess.size(); ++index) {
            if (index != 0) phase_definition += ',';
            phase_definition += recording_support::number(phase.composition_guess[index]);
        }
        if (phase.packing_guess.has_value())
            phase_definition += ":packing=" + recording_support::number(*phase.packing_guess);
    }
    std::string reaction_definition;
    for (const auto& reaction : problem.reactions) {
        if (!reaction_definition.empty()) reaction_definition += ',';
        reaction_definition += reaction.name + ":K="
            + (reaction.equilibrium_constant.has_value()
               ? recording_support::number(*reaction.equilibrium_constant)
               : reaction.affinity_j_per_mol.has_value() ? "A=" + recording_support::number(*reaction.affinity_j_per_mol)
               : std::string("correlation")) + ":nu=";
        for (std::size_t index = 0; index < reaction.stoichiometry.size(); ++index) {
            if (index != 0) reaction_definition += ',';
            reaction_definition += recording_support::number(reaction.stoichiometry[index]);
        }
    }
    const auto* model = compiled.model;
    if (model == nullptr) {
        return recording_support::fields({
            recording_support::field("compiled_model", "missing"),
            recording_support::field("temperature", problem.temperature.value),
            recording_support::field("pressure", problem.pressure.value),
            recording_support::field("feed_kind",
                                     epcsaft::equilibrium::feed_kind_name(problem.feed.kind)),
            recording_support::field("feed", problem.feed.values),
            recording_support::field("phase_definition", phase_definition),
            recording_support::field("reaction_definition", reaction_definition),
        });
    }
    const auto& parameters = model->parameters();
    std::string row_names;
    for (const auto& name : compiled.row_names) {
        if (!row_names.empty()) row_names += ',';
        row_names += name;
    }
    return recording_support::fields({
        recording_support::field("temperature", problem.temperature.value),
        recording_support::field("temperature_free", problem.temperature.free ? 1U : 0U),
        recording_support::field("temperature_absent", problem.temperature.absent ? 1U : 0U),
        recording_support::field("pressure", problem.pressure.value),
        recording_support::field("pressure_free", problem.pressure.free ? 1U : 0U),
        recording_support::field("pressure_absent", problem.pressure.absent ? 1U : 0U),
        recording_support::field("feed_kind",
                                 epcsaft::equilibrium::feed_kind_name(problem.feed.kind)),
        recording_support::field("feed", problem.feed.values),
        recording_support::field("feed_solvent_count", problem.feed.solvent.size()),
        recording_support::field("feed_solvent_mass_kg", problem.feed.solvent_mass_kg),
        recording_support::field("phase_definition", phase_definition),
        recording_support::field("reaction_definition", reaction_definition),
        recording_support::field("standard_gibbs", problem.standard_gibbs_j_per_mol),
        recording_support::field("topology_declared", problem.topology_declared ? 1U : 0U),
        recording_support::field("compiled_primal_size", compiled.primal_size),
        recording_support::field("compiled_rank", compiled.rank),
        recording_support::field("compiled_pinned_count", compiled.pinned_count),
        recording_support::field("compiled_free_intensive_count", compiled.free_intensive_count),
        recording_support::field("compiled_charge_row_count", compiled.charge_row_count),
        recording_support::field("compiled_rows", compiled.rows.size()),
        recording_support::field("compiled_rows_minus_unknowns", compiled.rows_minus_unknowns),
        recording_support::field("compiled_generic_nullity", compiled.generic_nullity),
        recording_support::field("compiled_jacobian_rank", compiled.jacobian_rank),
        recording_support::field("compiled_n_feed", compiled.n_feed),
        recording_support::field("compiled_feed_composition", compiled.feed_composition),
        recording_support::field("compiled_balance", compiled.balance),
        recording_support::field("compiled_feed_total", compiled.feed_total),
        recording_support::field("compiled_feasible", compiled.feasible ? 1U : 0U),
        recording_support::field("compiled_has_ions", compiled.has_ions ? 1U : 0U),
        recording_support::field("compiled_has_reaction", compiled.has_reaction ? 1U : 0U),
        recording_support::field("compiled_row_names", row_names),
        recording_support::field("model_component_count", parameters.component_count()),
        recording_support::field("model_segment_count", parameters.segment_count),
        recording_support::field("model_sigma_constant", parameters.sigma_constant),
        recording_support::field("model_sigma_term_component", parameters.sigma_term_component),
        recording_support::field("model_sigma_term_amplitude", parameters.sigma_term_amplitude),
        recording_support::field("model_sigma_term_exponent", parameters.sigma_term_exponent),
        recording_support::field("model_dispersion_energy_over_k", parameters.dispersion_energy_over_k),
        recording_support::field("model_molar_mass_kg_per_mol", parameters.molar_mass_kg_per_mol),
        recording_support::field("model_k_ij", parameters.k_ij),
        recording_support::field("model_k_ij_reciprocal_temperature_slope",
                                 parameters.k_ij_reciprocal_temperature_slope),
        recording_support::field("model_k_ij_reference_temperature_k",
                                 parameters.k_ij_reference_temperature_k),
        recording_support::field("model_l_ij", parameters.l_ij),
        recording_support::field("model_site_component", parameters.site_component),
        recording_support::field("model_site_multiplicity", parameters.site_multiplicity),
        recording_support::field("model_association_active", parameters.association_active),
        recording_support::field("model_association_energy_over_k", parameters.association_energy_over_k),
        recording_support::field("model_association_volume", parameters.association_volume),
        recording_support::field("model_charge_number", parameters.charge_number),
        recording_support::field("model_packing_diameter_angstrom", parameters.packing_diameter_angstrom),
        recording_support::field("model_debye_huckel_diameter_angstrom", parameters.debye_huckel_diameter_angstrom),
        recording_support::field("model_born_diameter_angstrom", parameters.born_diameter_angstrom),
        recording_support::field("model_solvation_factor", parameters.solvation_factor),
        recording_support::field("model_ionic_region_relative_permittivity",
                                 parameters.ionic_region_relative_permittivity),
        recording_support::field("model_born_c_shell", parameters.born_c_shell),
        recording_support::field("model_born_c_dielectric", parameters.born_c_dielectric),
        recording_support::field("model_electrolyte_form_declared",
                                 parameters.electrolyte_form.declared() ? 1U : 0U),
        recording_support::field("model_electrolyte_form",
                                 parameters.electrolyte_form.declared()
                                     ? static_cast<int>(parameters.electrolyte_form.require("electrolyte_form"))
                                     : -1),
        recording_support::field("model_ionic_dispersion_declared",
                                 parameters.ionic_dispersion.declared() ? 1U : 0U),
        recording_support::field("model_ionic_dispersion",
                                 parameters.ionic_dispersion.declared()
                                     ? static_cast<int>(parameters.ionic_dispersion.require("ionic_dispersion"))
                                     : -1),
        recording_support::field("model_relative_permittivity", parameters.relative_permittivity),
        recording_support::field("model_permittivity_term_component",
                                 parameters.permittivity_term_component),
        recording_support::field("model_permittivity_term_coefficient",
                                 parameters.permittivity_term_coefficient),
        recording_support::field("model_permittivity_term_power",
                                 parameters.permittivity_term_power),
        recording_support::field("model_permittivity_reference_temperature_k",
                                 parameters.permittivity_reference_temperature_k),
        recording_support::field("model_ion_suppression_coefficient", parameters.ion_suppression_coefficient),
        recording_support::field("model_permittivity_form_declared",
                                 parameters.permittivity_form.declared() ? 1U : 0U),
        recording_support::field("model_permittivity_form",
                                 parameters.permittivity_form.declared()
                                     ? static_cast<int>(parameters.permittivity_form.require("permittivity_form"))
                                     : -1),
        recording_support::field("model_dielectric_volume_m3_per_mol",
                                 parameters.dielectric_volume_m3_per_mol),
        recording_support::field("model_dielectric_temperature_k",
                                 parameters.dielectric_temperature_k),
        recording_support::field("model_dipole_scaling_debye2", parameters.dipole_scaling_debye2),
        recording_support::field("model_polarizability_scaling_angstrom3",
                                 parameters.polarizability_scaling_angstrom3),
        recording_support::field("model_correlation_integral_parameter",
                                 parameters.correlation_integral_parameter),
    });
}

// Old normalized correlation ln K = a + 1200/T + 0.4 ln T + 0.002 T; P3 carries the
// implicit 1 K reference through its legacy normalization at 300 K on [250, 450] K.
// The retained legacy offset a = 0.3 gives K ~ 870-1010 (dilute A, x_A ~ 1e-3); the
// owner-selected qualification offset makes K(350 K) = 2 with unchanged b, c, d.
constexpr double DILUTE_ISOMERISATION_A = 0.3;
double moderate_isomerisation_a() {
    return std::log(2.0) - (1200.0 / 350.0 + 0.4 * std::log(350.0) + 0.002 * 350.0);
}
epcsaft::ReactionLogPolynomial legacy_isomerisation(double sign, double a) {
    return epcsaft::normalize_legacy_reaction_correlation(
        sign * a, sign * 1200.0, sign * 0.4, sign * 0.002, 300.0, 250.0, 450.0,
        epcsaft::equilibrium::EOS_STANDARD_STATE_ID);
}

void check_finite_difference(const epcsaft::equilibrium::CompiledSystem& compiled,
                             const std::vector<double>& point) {
    const auto analytic = compiled.jacobian(point);
    const auto base = compiled.evaluate(point).values;
    const auto context = weighted_hessian_context(compiled);
    CHECK(analytic.size() == base.size() * point.size());
    for (std::size_t column = 0; column < point.size(); ++column) {
        auto plus = point;
        auto minus = point;
        const double delta = centered_step(point[column]);
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
            CHECK(verdict);
        }
    }
}

std::vector<std::size_t> callback_rows(
    const epcsaft::equilibrium::CompiledSystem& compiled
) {
    std::vector<std::size_t> result;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.pinned_count != 0
            || (compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Chemical
                && compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Pressure)) {
            result.push_back(row);
        }
    }
    return result;
}

std::vector<double> weighted_gradient(
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

std::vector<double> kkt_gradient(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers
) {
    const auto objective = compiled.objective_gradient(point);
    const auto jacobian = compiled.jacobian(point);
    const auto rows = callback_rows(compiled);
    if (objective.size() != point.size()
        || jacobian.size() != compiled.rows.size() * point.size()
        || multipliers.size() != rows.size()) {
        throw std::runtime_error("KKT gradient shape mismatch");
    }
    auto result = objective;
    for (std::size_t multiplier = 0; multiplier < rows.size(); ++multiplier) {
        for (std::size_t column = 0; column < point.size(); ++column) {
            result[column] += multipliers[multiplier]
                * jacobian[rows[multiplier] * point.size() + column];
        }
    }
    return result;
}

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) result = std::max(result, std::abs(value));
    return result;
}

void expect_weighted_hessian_gradient_ladder(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers,
    const char* label
) {
    epcsaft::equilibrium::GibbsMinimizer nlp(compiled, point);
    const auto shape = nlp.shape();
    const auto objective_factors = {0.0, 0.5, 1.0};
    const auto step_factors = {10.0, 3.0, 1.0};
    const std::size_t lower_triangle_count = point.size() * (point.size() + 1) / 2;
    const auto context = weighted_hessian_context(compiled);
    const auto definition = recording_support::fields({
        context,
        recording_support::field("case", label),
        recording_support::field("point", point),
        recording_support::field("multipliers", multipliers),
        recording_support::field("variables", point.size()),
        recording_support::field("constraints", multipliers.size()),
        recording_support::field("matrix", "lower_triangle"),
        recording_support::field("objective_factors", std::vector<double>{0.0, 0.5, 1.0}),
        recording_support::field("step_factors", std::vector<double>{10.0, 3.0, 1.0}),
        recording_support::field("hessian_entries", static_cast<std::size_t>(shape.hessian_entries)),
    });
    for (const double objective_factor : objective_factors) {
        for (const double step_factor : step_factors) {
            const auto prefix = std::string("hessian/") + label
                + "/objective-" + recording_support::number(objective_factor)
                + "/step-" + recording_support::number(step_factor) + "/";
            const auto group_definition = recording_support::fields({
                definition,
                recording_support::field("objective_factor", objective_factor),
                recording_support::field("step_factor", step_factor),
                recording_support::field("expected_rows", lower_triangle_count),
            });
            recording_support::expect_group(
                "native.equilibrium.f5", "f5-hessian-ladder", prefix.c_str(),
                lower_triangle_count, group_definition,
                recording_support::field(
                    "coordinate_schema", "case,objective_factor,step_factor,row,column"
                ).c_str()
            );
        }
    }
}

bool check_weighted_hessian_gradient_ladder(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers,
    const char* label,
    double tolerance,
    std::array<std::size_t, 3>* false_verdicts = nullptr
) {
    const auto rows = callback_rows(compiled);
    epcsaft::equilibrium::GibbsMinimizer nlp(compiled, point);
    const auto shape = nlp.shape();
    CHECK(shape.variables >= 0 && shape.constraints >= 0);
    CHECK(point.size() == static_cast<std::size_t>(shape.variables));
    CHECK(multipliers.size() == rows.size());
    CHECK(multipliers.size() == static_cast<std::size_t>(shape.constraints));
    if (point.size() != static_cast<std::size_t>(shape.variables)
        || multipliers.size() != static_cast<std::size_t>(shape.constraints)) return false;

    const auto objective_factors = {0.0, 0.5, 1.0};
    const auto step_factors = {10.0, 3.0, 1.0};
    const auto context = weighted_hessian_context(compiled);

    std::vector<Ipopt::Index> hessian_rows(static_cast<std::size_t>(shape.hessian_entries));
    std::vector<Ipopt::Index> hessian_columns(static_cast<std::size_t>(shape.hessian_entries));
    CHECK(nlp.eval_h(shape.variables, point.data(), true, 1.0, shape.constraints,
                     multipliers.data(), true, shape.hessian_entries,
                     hessian_rows.data(), hessian_columns.data(), nullptr));
    bool all_verdicts = true;
    std::array<std::size_t, 3> false_counts{};
    std::size_t objective_index = 0;
    for (const double objective_factor : objective_factors) {
        std::vector<Ipopt::Number> hessian(
            static_cast<std::size_t>(shape.hessian_entries), 0.0
        );
        CHECK(nlp.eval_h(shape.variables, point.data(), false, objective_factor,
                         shape.constraints, multipliers.data(), false,
                         shape.hessian_entries, hessian_rows.data(),
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
                    const bool verdict = close_with_tolerance(hessian[index], centered, tolerance);
                    recording_support::compare(
                        comparison_id.c_str(), "weighted_lagrangian_hessian", "coordinate",
                        hessian[index], centered, inputs, "Relationship",
                        "weighted_hessian_gradient_finite_difference",
                        "nabla2_L=d(weighted_gradient)/d(coordinate)", "reference_or_one",
                        0.0, tolerance, verdict, "f5-hessian-ladder", point_id.c_str(),
                        "\\nabla^2\\mathcal L", "\\mathbf{z},\\lambda,\\Delta z",
                        "\\nabla^2\\mathcal L_{ij}=\\Delta(\\nabla\\mathcal L)_i/\\Delta z_j",
                        coordinates.c_str());
                    all_verdicts = all_verdicts && verdict;
                    if (!verdict) ++false_counts[objective_index];
                }
            }
        }
        ++objective_index;
    }
    if (false_verdicts != nullptr) *false_verdicts = false_counts;
    return all_verdicts;
}

void verify_s1(const epcsaft::HelmholtzModel& model,
               const epcsaft::equilibrium::EquilibriumProblem& problem,
               const epcsaft::equilibrium::EquilibriumResult& result,
               const char* case_id) {
    CHECK(result.success);
    if (!result.success) return;
    const double limit = convergence_limit(result.convergence);
    CHECK(result.residuals.size() == result.rows.size());
    for (std::size_t row = 0; row < result.residuals.size(); ++row) {
        const auto inputs = recording_support::fields({
            recording_support::field("T", result.temperature.front()),
            recording_support::field("P", result.pressure.front()),
            recording_support::field("coordinates", result.coordinates),
            recording_support::field("row", row), recording_support::field("limit", limit),
        });
        recording_support::compare(
            (std::string("residual/") + case_id + "/" + std::to_string(row)).c_str(),
            "equilibrium_residual", "1", result.residuals[row], 0.0, inputs,
            "Relationship", "equilibrium-residual", "residual is within convergence limit",
            "reference_or_one", limit, 0.0, std::abs(result.residuals[row]) <= limit,
            "s1-residual", (std::string(case_id) + "/" + std::to_string(row)).c_str(),
            "r", "T,P,\\mathbf{z},row",
            "|r_i|\\leq\\ell", "row");
        CHECK(std::abs(result.residuals[row]) <= limit);
    }
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    CHECK(result.convergence.residual == compiled.evaluate(result.coordinates).residual_norm());
    const std::size_t components = model.parameters().component_count();
    const double rt = epcsaft::GAS_CONSTANT_J_PER_MOL_K * result.temperature.front();
    const auto charges = epcsaft::equilibrium::species_z(model.parameters());
    for (std::size_t phase = 0; phase < result.phase_count; ++phase) {
        const auto first = result.mole_fractions.begin()
            + static_cast<std::ptrdiff_t>(phase * components);
        const std::vector<double> x(first, first + static_cast<std::ptrdiff_t>(components));
        const double volume = 1.0 / result.molar_densities[phase];
        const auto derivative = model.evaluate(
            result.temperature.front(), volume, x, epcsaft::DerivativeOrder::Gradient);
        double gibbs = derivative.value + result.pressure.front() * volume / rt;
        double lambda_Bx = 0.0;
        double charge = 0.0;
        std::vector<double> chemical(components);
        for (std::size_t species = 0; species < components; ++species) {
            chemical[species] = derivative.gradient[2 + species]
                + std::log(std::max(x[species] / volume,
                                    std::numeric_limits<double>::min()))
                + compiled.standard_gibbs[species] / rt;
            gibbs += x[species] * (
                std::log(std::max(x[species], std::numeric_limits<double>::min()))
                - std::log(volume) - 1.0 + compiled.standard_gibbs[species] / rt);
            charge += charges[species] * x[species];
            for (std::size_t row = 0; row < compiled.rank; ++row)
                lambda_Bx += result.multipliers[row] * compiled.B[row][species] * x[species];
        }
        const auto phase_inputs = recording_support::fields({
            recording_support::field("T", result.temperature.front()),
            recording_support::field("P", result.pressure.front()),
            recording_support::field("x", x), recording_support::field("phase", phase),
            recording_support::field("limit", limit),
        });
        recording_support::compare(
            (std::string("gibbs-dual/") + case_id + "/" + std::to_string(phase)).c_str(),
            "gibbs_dual_residual", "1", gibbs, lambda_Bx, phase_inputs,
            "Relationship", "KKT", "phase Gibbs value equals balance dual",
            "reference_or_one", limit, 0.0, std::abs(gibbs - lambda_Bx) <= limit,
            "s1-gibbs-dual", (std::string(case_id) + "/" + std::to_string(phase)).c_str(),
            "g", "T,P,\\mathbf{x},phase",
            "g=\\lambda^TBx", "phase");
        CHECK(std::abs(gibbs - lambda_Bx) <= limit);
        bool ion_phase = false;
        for (const auto species : compiled.problem.phases[phase].support)
            ion_phase = ion_phase || charges[species] != 0;
        if (ion_phase) {
            recording_support::compare(
                (std::string("charge/") + case_id + "/" + std::to_string(phase)).c_str(),
                "phase_charge", "1", charge, 0.0, phase_inputs,
                "Relationship", "charge-balance", "ionic phase charge is zero",
                "reference_or_one", limit, 0.0, std::abs(charge) <= limit,
                "s1-charge", (std::string(case_id) + "/" + std::to_string(phase)).c_str(),
                "q", "T,P,\\mathbf{x},phase",
                "q=0", "phase");
            CHECK(std::abs(charge) <= limit);
        }
        for (const auto& reaction : compiled.problem.reactions) {
            double affinity = 0.0;
            for (std::size_t species = 0; species < components; ++species)
                affinity += reaction.stoichiometry[species] * chemical[species];
            recording_support::compare(
                (std::string("affinity/") + case_id + "/" + std::to_string(phase)
                 + "/" + reaction.name).c_str(),
                "reaction_affinity", "J/mol", affinity, 0.0, phase_inputs,
                "Relationship", "reaction-equilibrium", "reaction affinity is zero",
                "reference_or_one", limit, 0.0, std::abs(affinity) <= limit,
                "s1-affinity", (std::string(case_id) + "/" + std::to_string(phase)
                                 + "/" + reaction.name).c_str(),
                "A", "T,P,\\mathbf{x},phase",
                "A=0", "reaction");
            CHECK(std::abs(affinity) <= limit);
        }
    }
    double feed_dual = 0.0;
    for (std::size_t row = 0; row < compiled.rank; ++row)
        feed_dual += result.multipliers[row] * compiled.balance[row];
    recording_support::compare(
        (std::string("feed-dual/") + case_id).c_str(), "objective", "1", result.objective, feed_dual,
        recording_support::fields({recording_support::field("T", result.temperature.front()),
                                   recording_support::field("P", result.pressure.front()),
                                   recording_support::field("coordinates", result.coordinates),
                                   recording_support::field("limit", limit)}),
        "Relationship", "KKT", "objective equals feed balance dual",
        "reference_or_one", limit, 0.0, std::abs(result.objective - feed_dual) <= limit,
        "s1-feed-dual", case_id, "g", "T,P,\\mathbf{z}", "g=\\lambda^Tb", "objective");
    CHECK(std::abs(result.objective - feed_dual) <= limit);
}

void check_s1() {
    epcsaft::HelmholtzModel binary_model(properties_support::neutral_binary());
    const auto liquid = equilibrium_support::flash(2);
    epcsaft::core::Capture capture;
    const auto liquid_result = epcsaft::equilibrium::solve_equilibrium(binary_model, liquid);
    verify_s1(binary_model, liquid, liquid_result, "liquid");
    const auto span = std::find_if(capture.trace().spans.begin(), capture.trace().spans.end(),
        [](const auto& item) { return std::string(item.name) == "solve_equilibrium"; });
    CHECK(span != capture.trace().spans.end() && span->convergence.has_value());
    if (span != capture.trace().spans.end() && span->convergence.has_value()) {
        const auto& record = *span->convergence;
        CHECK(record.residual == liquid_result.convergence.residual);
        CHECK(record.criterion == liquid_result.convergence.criterion);
        CHECK(record.floor == liquid_result.convergence.floor);
        CHECK(record.amplification == liquid_result.convergence.amplification);
        CHECK(record.contraction == liquid_result.convergence.contraction);
        CHECK(record.iterations == liquid_result.convergence.iterations);
        CHECK(record.evaluations == liquid_result.convergence.evaluations);
        CHECK(record.stop == liquid_result.convergence.stop);
    }
    auto vapor = equilibrium_support::flash(2);
    vapor.phases[0].kind = epcsaft::equilibrium::PhaseKind::Vapor;
    verify_s1(binary_model, vapor,
              epcsaft::equilibrium::solve_equilibrium(binary_model, vapor), "vapor");
    auto flash = equilibrium_support::flash(3);
    flash.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                    equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    flash.temperature.value = 180.0;
    epcsaft::HelmholtzModel ternary_model(properties_support::neutral_ternary());
    verify_s1(ternary_model, flash,
              epcsaft::equilibrium::solve_equilibrium(ternary_model, flash), "ternary");
    epcsaft::HelmholtzModel reactive_model(equilibrium_support::synthetic_reactive_ternary());
    const auto reactive = equilibrium_support::synthetic_reactive_problem();
    verify_s1(reactive_model, reactive,
              epcsaft::equilibrium::solve_equilibrium(reactive_model, reactive), "reactive");
    auto ethane = equilibrium_support::flash(1);
    ethane.phases = {equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
                     equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
                                                epcsaft::equilibrium::AmountSpec::pinned(0.0))};
    ethane.temperature.value = 250.0;
    ethane.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    epcsaft::HelmholtzModel ethane_model(properties_support::ethane());
    verify_s1(ethane_model, ethane,
              epcsaft::equilibrium::solve_equilibrium(ethane_model, ethane), "ethane");
    auto salt = equilibrium_support::flash(3);
    salt.feed.values = {0.98, 0.01, 0.01};
    salt.temperature.value = 298.15;
    epcsaft::HelmholtzModel salt_model(properties_support::aqueous_electrolyte());
    verify_s1(salt_model, salt,
              epcsaft::equilibrium::solve_equilibrium(salt_model, salt), "salt");
}

void check_s2() {
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto binary = equilibrium_support::flash(3);
    binary.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                     equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    binary.temperature.value = 180.0;
    binary.pressure.value = 1.0e5;
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, binary);
    recording_support::compare(
        "binary-success", "success", "bool", result.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", binary.temperature.value),
                                   recording_support::field("P", binary.pressure.value),
                                   recording_support::field("feed", binary.feed.values)}),
        "Relationship", "stability", "binary solution succeeds", "reference_or_one", 0.0, 0.0,
        result.success, "s2-binary-success", "binary", "success", "T,P,feed",
        "success=1", "case");
    const bool binary_curvature_positive = std::isfinite(result.reduced_hessian_minimum)
        && result.reduced_hessian_minimum > 0.0;
    recording_support::compare_custom(
        "binary-curvature", "reduced_hessian_minimum", "1", result.reduced_hessian_minimum,
        0.0, recording_support::fields({recording_support::field("T", binary.temperature.value),
                                         recording_support::field("P", binary.pressure.value),
                                         recording_support::field("feed", binary.feed.values)}),
        "Relationship", "stability", "reduced Hessian minimum is positive",
        "isfinite(model) && model > reference", binary_curvature_positive,
        "s2-binary-curvature", "binary", "\\lambda_{min}(H)", "T,P,feed",
        "\\lambda_{min}(H)>0", "stability");
    CHECK(result.success && binary_curvature_positive);
    const auto reactive = epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::synthetic_reactive_ternary()),
        equilibrium_support::synthetic_reactive_problem());
    recording_support::compare(
        "reactive-success", "success", "bool", reactive.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", std::vector<double>{0.05, 0.475, 0.475})}),
        "Relationship", "stability", "reactive solution succeeds", "reference_or_one", 0.0, 0.0,
        reactive.success, "s2-reactive-success", "reactive", "success", "T,P,feed",
        "success=1", "case");
    const bool reactive_curvature_positive = std::isfinite(reactive.reduced_hessian_minimum)
        && reactive.reduced_hessian_minimum > 0.0;
    recording_support::compare_custom(
        "reactive-curvature", "reduced_hessian_minimum", "1", reactive.reduced_hessian_minimum,
        0.0, recording_support::fields({recording_support::field("T", 298.15),
                                         recording_support::field("P", 1.0e5),
                                         recording_support::field("feed", std::vector<double>{0.05, 0.475, 0.475})}),
        "Relationship", "stability", "reactive reduced Hessian minimum is positive",
        "isfinite(model) && model > reference", reactive_curvature_positive,
        "s2-reactive-curvature", "reactive", "\\lambda_{min}(H)", "T,P,feed",
        "\\lambda_{min}(H)>0", "stability");
    CHECK(reactive.success && reactive_curvature_positive);
}

void check_s6() {
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto two = equilibrium_support::synthetic_reactive_problem();
    two.reactions.clear();
    auto one = two;
    one.phases.resize(1);
    const auto one_result = epcsaft::equilibrium::solve_equilibrium(model, one);
    const auto two_result = epcsaft::equilibrium::solve_equilibrium(model, two);
    recording_support::compare(
        "one-success", "success", "bool", one_result.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", two.feed.values),
                                   recording_support::field("phase_count", 1U)}),
        "Relationship", "solve", "one-phase solve succeeds", "reference_or_one", 0.0, 0.0,
        one_result.success, "s6-one-success", "one", "success", "T,P,feed,phases", "success=1", "case");
    recording_support::compare(
        "two-success", "success", "bool", two_result.success ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", two.feed.values),
                                   recording_support::field("phase_count", 2U)}),
        "Relationship", "solve", "two-phase solve succeeds", "reference_or_one", 0.0, 0.0,
        two_result.success, "s6-two-success", "two", "success", "T,P,feed,phases", "success=1", "case");
    const double objective_gap = one_result.objective - two_result.objective;
    const double objective_limit = two_result.convergence.criterion == nullptr
        ? 0.0 : convergence_limit(two_result.convergence);
    const bool objective_gap_ok = std::isfinite(objective_gap) && std::isfinite(objective_limit)
        && objective_gap > objective_limit;
    recording_support::compare_custom(
        "objective-gap", "objective_gap", "1", objective_gap, objective_limit,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", two.feed.values),
                                   recording_support::field("one_objective", one_result.objective),
                                   recording_support::field("two_objective", two_result.objective),
                                   recording_support::field("convergence_limit", objective_limit)}),
        "Relationship", "phase-split", "two-phase objective gap exceeds convergence limit",
        "isfinite(model) && isfinite(reference) && model > reference", objective_gap_ok,
        "s6-objective-gap", "two-minus-one", "\\Delta g", "T,P,feed,objective,limit",
        "\\Delta g>\\ell", "objective");
    CHECK(one_result.success && two_result.success);
    if (!one_result.success || !two_result.success) return;
    CHECK(two_result.objective < one_result.objective
          - convergence_limit(two_result.convergence));
}

void check_s7() {
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto problem = equilibrium_support::synthetic_reactive_problem();
    problem.reactions.clear();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(first.success);
    if (!first.success) return;
    for (const double fraction : {0.3, 0.7}) {
        for (std::size_t i = 0; i < 3; ++i)
            problem.feed.values[i] = (1.0 - fraction) * first.mole_fractions[i]
                + fraction * first.mole_fractions[3 + i];
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem, first);
        CHECK(result.success);
        if (!result.success) continue;
        for (std::size_t i = 0; i < 6; ++i) {
            const double difference = std::abs(
                result.mole_fractions[i] - first.mole_fractions[i]);
            const auto point = recording_support::number(fraction) + "/" + std::to_string(i);
            recording_support::compare(
                (std::string("continuation/") + point).c_str(),
                "mole_fraction_difference", "1", difference, 0.0,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("feed", problem.feed.values),
                                           recording_support::field("fraction", fraction),
                                           recording_support::field("component", i),
                                           recording_support::field("reference_solution", first.mole_fractions),
                                           recording_support::field("solution", result.mole_fractions)}),
                "Relationship", "continuation", "continuation preserves the solved composition",
                "reference_or_one", 1.0e-8, 0.0, difference <= 1.0e-8,
                "s7-continuation", point.c_str(), "\\Delta x_i",
                "T,P,feed,\\alpha,component", "\\Delta x_i\\leq10^{-8}", "component");
            CHECK(difference <= 1.0e-8);
        }
    }
}

void check_c1_amount_scale() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto normal = equilibrium_support::flash(2);
    auto tiny = normal;
    normal.feed.values = {0.5, 0.5};
    tiny.feed.values = {5.0e-10, 5.0e-10};
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, normal);
    const auto second = epcsaft::equilibrium::solve_equilibrium(model, tiny);
    CHECK(first.success && second.success);
    if (!first.success || !second.success) return;
    CHECK(first.amounts.size() == second.amounts.size());
    const double first_total = normal.feed.values[0] + normal.feed.values[1];
    const double second_total = tiny.feed.values[0] + tiny.feed.values[1];
    for (std::size_t i = 0; i < first.amounts.size(); ++i) {
        const double first_normalized = first.amounts[i] / first_total;
        const double second_normalized = second.amounts[i] / second_total;
        CHECK(std::abs(first_normalized - 1.0) <= 1.0e-7);
        CHECK(std::abs(second_normalized - 1.0) <= 1.0e-7);
        const auto normal_inputs = recording_support::fields({
            recording_support::field("T", 300.0),
            recording_support::field("P", 1.0e5),
            recording_support::field("feed", normal.feed.values),
            recording_support::field("solver_status", first.solver_status),
            recording_support::field("convergence_stop",
                epcsaft::core::stop_name(first.convergence.stop)),
            recording_support::field("residual", first.convergence.residual),
            recording_support::field("phase", i),
        });
        recording_support::compare(
            (std::string("c1-amount-scale/normal/") + std::to_string(i)).c_str(),
            "phase_amount/feed_total", "1", first_normalized, 1.0, normal_inputs,
            "Relationship", "extensive-scaling", "normalized phase amount remains one",
            "strict_absolute", 1.0e-7, 0.0, std::abs(first_normalized - 1.0) <= 1.0e-7,
            "c1-amount-scale", "normal", "n/feed_total", "T,P,feed,phase",
            "n_phase/feed_total=1", "phase");
        const auto tiny_inputs = recording_support::fields({
            recording_support::field("T", 300.0),
            recording_support::field("P", 1.0e5),
            recording_support::field("feed", tiny.feed.values),
            recording_support::field("solver_status", second.solver_status),
            recording_support::field("convergence_stop",
                epcsaft::core::stop_name(second.convergence.stop)),
            recording_support::field("residual", second.convergence.residual),
            recording_support::field("phase", i),
        });
        recording_support::compare(
            (std::string("c1-amount-scale/tiny/") + std::to_string(i)).c_str(),
            "phase_amount/feed_total", "1", second_normalized, 1.0, tiny_inputs,
            "Relationship", "extensive-scaling", "normalized phase amount remains one",
            "strict_absolute", 1.0e-7, 0.0, std::abs(second_normalized - 1.0) <= 1.0e-7,
            "c1-amount-scale", "tiny", "n/feed_total", "T,P,feed,phase",
            "n_phase/feed_total=1", "phase");
    }
}

void check_c1_multiplier_scale() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto small = equilibrium_support::flash(2);
    auto large = small;
    large.feed.values = {5.0, 5.0};
    const auto small_result = epcsaft::equilibrium::solve_equilibrium(model, small);
    const auto large_result = epcsaft::equilibrium::solve_equilibrium(model, large);
    CHECK(small_result.success && large_result.success);
    if (!small_result.success || !large_result.success) return;

    const auto small_compiled = epcsaft::equilibrium::compile_problem(model, small);
    const auto large_compiled = epcsaft::equilibrium::compile_problem(model, large);
    const auto small_rows = callback_rows(small_compiled);
    const auto large_rows = callback_rows(large_compiled);
    CHECK(small_rows == large_rows);
    CHECK(small_result.constraint_multipliers.size() == small_rows.size());
    CHECK(large_result.constraint_multipliers.size() == large_rows.size());
    if (small_rows != large_rows
        || small_result.constraint_multipliers.size() != small_rows.size()
        || large_result.constraint_multipliers.size() != large_rows.size()) {
        return;
    }

    const std::vector<double> small_point(
        small_result.coordinates.begin(),
        small_result.coordinates.begin()
            + static_cast<std::ptrdiff_t>(small_compiled.primal_size));
    const std::vector<double> large_point(
        large_result.coordinates.begin(),
        large_result.coordinates.begin()
            + static_cast<std::ptrdiff_t>(large_compiled.primal_size));
    const auto small_objective = small_compiled.objective_gradient(small_point);
    const auto large_objective = large_compiled.objective_gradient(large_point);
    const auto small_kkt = kkt_gradient(
        small_compiled, small_point, small_result.constraint_multipliers);
    const auto large_kkt = kkt_gradient(
        large_compiled, large_point, large_result.constraint_multipliers);
    const double small_kkt_scale = std::max(1.0, maximum_abs(small_objective));
    const double large_kkt_scale = std::max(1.0, maximum_abs(large_objective));
    const double small_kkt_limit =
        100.0 * epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value * small_kkt_scale;
    const double large_kkt_limit =
        100.0 * epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value * large_kkt_scale;
    const double small_kkt_residual = maximum_abs(small_kkt);
    const double large_kkt_residual = maximum_abs(large_kkt);
    const auto kkt_inputs = recording_support::fields({
        recording_support::field("T", 300.0),
        recording_support::field("P", 1.0e5),
        recording_support::field("feed_small", small.feed.values),
        recording_support::field("feed_large", large.feed.values),
        recording_support::field("feed_total_small", small_compiled.feed_total),
        recording_support::field("feed_total_large", large_compiled.feed_total),
        recording_support::field("small_solver_status", small_result.solver_status),
        recording_support::field("large_solver_status", large_result.solver_status),
        recording_support::field("small_convergence_stop",
            epcsaft::core::stop_name(small_result.convergence.stop)),
        recording_support::field("large_convergence_stop",
            epcsaft::core::stop_name(large_result.convergence.stop)),
        recording_support::field("active_constraint_count", small_rows.size()),
        recording_support::field("small_objective_gradient", small_objective),
        recording_support::field("large_objective_gradient", large_objective),
        recording_support::field("small_constraint_multipliers",
            small_result.constraint_multipliers),
        recording_support::field("large_constraint_multipliers",
            large_result.constraint_multipliers),
        recording_support::field("small_kkt_scale", small_kkt_scale),
        recording_support::field("large_kkt_scale", large_kkt_scale),
    });
    recording_support::compare(
        "c1-multiplier-scale/kkt-small", "kkt_residual", "1", small_kkt_residual, 0.0,
        kkt_inputs, "Relationship", "KKT", "small-feed returned multipliers satisfy KKT",
        "reference_or_one", small_kkt_limit, 0.0,
        small_kkt_residual <= small_kkt_limit, "c1-multiplier-scale", "small",
        "||r_KKT||", "T,P,feed,lambda,scale", "||r_KKT||<=100*row_residual*scale", "feed");
    recording_support::compare(
        "c1-multiplier-scale/kkt-large", "kkt_residual", "1", large_kkt_residual, 0.0,
        kkt_inputs, "Relationship", "KKT", "large-feed returned multipliers satisfy KKT",
        "reference_or_one", large_kkt_limit, 0.0,
        large_kkt_residual <= large_kkt_limit, "c1-multiplier-scale", "large",
        "||r_KKT||", "T,P,feed,lambda,scale", "||r_KKT||<=100*row_residual*scale", "feed");
    CHECK(small_kkt_residual <= small_kkt_limit);
    CHECK(large_kkt_residual <= large_kkt_limit);

    const double feed_scale = large_compiled.feed_total / small_compiled.feed_total;
    double multiplier_difference = 0.0;
    for (std::size_t i = 0; i < small_result.constraint_multipliers.size(); ++i) {
        const double expected = feed_scale * small_result.constraint_multipliers[i];
        multiplier_difference = std::max(
            multiplier_difference, std::abs(large_result.constraint_multipliers[i] - expected)
                / std::max(1.0, std::abs(expected)));
    }
    recording_support::compare(
        "c1-multiplier-scale/multipliers", "relative_multiplier_difference", "1",
        multiplier_difference, 0.0, kkt_inputs, "Relationship", "extensive-scaling",
        "constraint multipliers scale with feed amount", "reference_or_one", 1.0e-7, 0.0,
        multiplier_difference <= 1.0e-7, "c1-multiplier-scale", "multipliers",
        "delta_lambda", "T,P,feed,lambda", "lambda_large=10*lambda_small", "feed");
    CHECK(multiplier_difference <= 1.0e-7);

    double intensive_difference = 0.0;
    intensive_difference = std::max(
        intensive_difference,
        std::abs(small_result.temperature.front() - large_result.temperature.front())
            / std::max(1.0, std::abs(small_result.temperature.front())));
    intensive_difference = std::max(
        intensive_difference,
        std::abs(small_result.pressure.front() - large_result.pressure.front())
            / std::max(1.0, std::abs(small_result.pressure.front())));
    for (std::size_t i = 0; i < small_result.mole_fractions.size(); ++i) {
        intensive_difference = std::max(
            intensive_difference,
            std::abs(small_result.mole_fractions[i] - large_result.mole_fractions[i]));
    }
    for (std::size_t i = 0; i < small_result.molar_densities.size(); ++i) {
        intensive_difference = std::max(
            intensive_difference,
            std::abs(small_result.molar_densities[i] - large_result.molar_densities[i])
                / std::max(1.0, std::abs(small_result.molar_densities[i])));
        intensive_difference = std::max(
            intensive_difference,
            std::abs(small_result.packing_fractions[i] - large_result.packing_fractions[i]));
    }
    recording_support::compare(
        "c1-multiplier-scale/intensives", "relative_intensive_difference", "1",
        intensive_difference, 0.0, kkt_inputs, "Relationship", "extensive-scaling",
        "intensive state is invariant under feed scaling", "reference_or_one", 1.0e-8, 0.0,
        intensive_difference <= 1.0e-8, "c1-multiplier-scale", "intensives",
        "delta_state", "T,P,feed", "state_large=state_small", "feed");
    CHECK(intensive_difference <= 1.0e-8);
    CHECK(small_result.reduced_hessian_minimum > 0.0);
    CHECK(large_result.reduced_hessian_minimum > 0.0);

    auto double_scaled = large_result.constraint_multipliers;
    for (double& value : double_scaled) value *= feed_scale;
    const double double_scaled_kkt = maximum_abs(
        kkt_gradient(large_compiled, large_point, double_scaled));
    recording_support::compare_custom(
        "c1-multiplier-scale/double-scaling", "kkt_residual", "1", double_scaled_kkt,
        large_kkt_limit, kkt_inputs, "Relationship", "KKT-regression-control",
        "extra feed scaling violates the KKT residual", "model > reference",
        double_scaled_kkt > large_kkt_limit, "c1-multiplier-scale", "double-scaling",
        "||r_KKT||", "T,P,feed,lambda", "double-scaled lambda rejected", "feed");
    CHECK(double_scaled_kkt > large_kkt_limit);

    auto tiny = small;
    tiny.feed.values = {5.0e-10, 5.0e-10};
    const auto tiny_result = epcsaft::equilibrium::solve_equilibrium(model, tiny);
    CHECK(tiny_result.success);
    if (!tiny_result.success) return;
    const auto tiny_compiled = epcsaft::equilibrium::compile_problem(model, tiny);
    const auto tiny_rows = callback_rows(tiny_compiled);
    CHECK(tiny_rows == small_rows);
    CHECK(tiny_result.constraint_multipliers.size() == tiny_rows.size());
    if (tiny_rows != small_rows
        || tiny_result.constraint_multipliers.size() != tiny_rows.size()) return;

    const auto tiny_point = std::vector<double>(
        tiny_result.coordinates.begin(),
        tiny_result.coordinates.begin()
            + static_cast<std::ptrdiff_t>(tiny_compiled.primal_size));
    const auto tiny_objective = tiny_compiled.objective_gradient(tiny_point);
    const auto tiny_kkt = kkt_gradient(
        tiny_compiled, tiny_point, tiny_result.constraint_multipliers);
    const double tiny_kkt_scale = maximum_abs(tiny_objective);
    CHECK(tiny_kkt_scale > 0.0);
    const double tiny_kkt_limit =
        100.0 * epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value * tiny_kkt_scale;
    const double tiny_kkt_residual = maximum_abs(tiny_kkt);
    const auto tiny_inputs = recording_support::fields({
        recording_support::field("T", 300.0),
        recording_support::field("P", 1.0e5),
        recording_support::field("feed", tiny.feed.values),
        recording_support::field("feed_total", tiny_compiled.feed_total),
        recording_support::field("solver_status", tiny_result.solver_status),
        recording_support::field("convergence_stop",
            epcsaft::core::stop_name(tiny_result.convergence.stop)),
        recording_support::field("residual", tiny_result.convergence.residual),
        recording_support::field("convergence_limit", convergence_limit(tiny_result.convergence)),
        recording_support::field("objective_gradient", tiny_objective),
        recording_support::field("constraint_multipliers",
            tiny_result.constraint_multipliers),
        recording_support::field("kkt_scale", tiny_kkt_scale),
        recording_support::field("kkt_limit", tiny_kkt_limit),
        recording_support::field("reduced_hessian_minimum",
            tiny_result.reduced_hessian_minimum),
    });
    const bool tiny_status_ok = tiny_result.solver_status == "STOP_AT_TINY_STEP";
    recording_support::compare(
        "c1-multiplier-scale/tiny-status", "tiny_step_status", "bool",
        tiny_status_ok ? 1.0 : 0.0, 1.0, tiny_inputs, "Relationship", "Ipopt",
        "tiny-feed case uses the recorded tiny-step termination", "absolute", 0.0, 0.0,
        tiny_status_ok, "c1-multiplier-scale", "tiny-status", "status", "feed",
        "status=STOP_AT_TINY_STEP", "feed");
    const bool tiny_stop_ok = tiny_result.convergence.stop == epcsaft::core::Stop::converged;
    recording_support::compare(
        "c1-multiplier-scale/tiny-convergence", "convergence_stop", "bool",
        tiny_stop_ok ? 1.0 : 0.0, 1.0, tiny_inputs, "Relationship", "equilibrium-residual",
        "tiny-step is accepted only after residual convergence", "absolute", 0.0, 0.0,
        tiny_stop_ok, "c1-multiplier-scale", "tiny-convergence", "stop", "feed,residual",
        "stop=converged", "feed");
    recording_support::compare(
        "c1-multiplier-scale/tiny-kkt", "kkt_residual", "1", tiny_kkt_residual, 0.0,
        tiny_inputs, "Relationship", "KKT", "tiny-feed returned multipliers satisfy KKT",
        "strict_absolute", tiny_kkt_limit, 0.0, tiny_kkt_residual <= tiny_kkt_limit,
        "c1-multiplier-scale", "tiny", "||r_KKT||", "T,P,feed,lambda,scale",
        "||r_KKT||<=100*row_residual*||grad f||", "feed");
    const bool tiny_reduced_hessian_ok =
        std::isfinite(tiny_result.reduced_hessian_minimum)
        && tiny_result.reduced_hessian_minimum > 0.0;
    recording_support::compare_custom(
        "c1-multiplier-scale/tiny-reduced-hessian", "reduced_hessian_minimum", "1",
        tiny_result.reduced_hessian_minimum, 0.0, tiny_inputs, "Relationship", "KKT",
        "tiny-feed reduced Hessian remains positive", "model > 0", tiny_reduced_hessian_ok,
        "c1-multiplier-scale", "tiny-reduced-hessian", "lambda_min", "feed,lambda",
        "lambda_min>0", "feed");
    CHECK(tiny_status_ok);
    CHECK(tiny_stop_ok);
    CHECK(tiny_kkt_residual <= tiny_kkt_limit);
    CHECK(tiny_reduced_hessian_ok);

    const double tiny_feed_scale = tiny_compiled.feed_total / small_compiled.feed_total;
    double tiny_multiplier_difference = 0.0;
    for (std::size_t i = 0; i < tiny_result.constraint_multipliers.size(); ++i) {
        const double expected = tiny_feed_scale * small_result.constraint_multipliers[i];
        const double denominator = std::max(std::abs(expected), tiny_kkt_limit);
        tiny_multiplier_difference = std::max(
            tiny_multiplier_difference,
            std::abs(tiny_result.constraint_multipliers[i] - expected) / denominator);
    }
    recording_support::compare(
        "c1-multiplier-scale/tiny-multipliers", "relative_multiplier_difference", "1",
        tiny_multiplier_difference, 0.0, tiny_inputs, "Relationship", "extensive-scaling",
        "tiny-feed multipliers scale with feed amount", "reference_or_one", 1.0e-6, 0.0,
        tiny_multiplier_difference <= 1.0e-6, "c1-multiplier-scale", "tiny-multipliers",
        "delta_lambda", "feed,lambda,scale", "lambda_tiny=1e-9*lambda_small", "feed");
    CHECK(tiny_multiplier_difference <= 1.0e-6);

    double tiny_intensive_difference = 0.0;
    tiny_intensive_difference = std::max(
        tiny_intensive_difference,
        std::abs(tiny_result.temperature.front() - small_result.temperature.front())
            / std::max(1.0, std::abs(small_result.temperature.front())));
    tiny_intensive_difference = std::max(
        tiny_intensive_difference,
        std::abs(tiny_result.pressure.front() - small_result.pressure.front())
            / std::max(1.0, std::abs(small_result.pressure.front())));
    for (std::size_t i = 0; i < small_result.mole_fractions.size(); ++i) {
        tiny_intensive_difference = std::max(
            tiny_intensive_difference,
            std::abs(tiny_result.mole_fractions[i] - small_result.mole_fractions[i]));
    }
    for (std::size_t i = 0; i < small_result.molar_densities.size(); ++i) {
        tiny_intensive_difference = std::max(
            tiny_intensive_difference,
            std::abs(tiny_result.molar_densities[i] - small_result.molar_densities[i])
                / std::max(1.0, std::abs(small_result.molar_densities[i])));
        tiny_intensive_difference = std::max(
            tiny_intensive_difference,
            std::abs(tiny_result.packing_fractions[i] - small_result.packing_fractions[i]));
    }
    recording_support::compare(
        "c1-multiplier-scale/tiny-intensives", "relative_intensive_difference", "1",
        tiny_intensive_difference, 0.0, tiny_inputs, "Relationship", "extensive-scaling",
        "tiny-feed intensive state is invariant", "reference_or_one", 1.0e-8, 0.0,
        tiny_intensive_difference <= 1.0e-8, "c1-multiplier-scale", "tiny-intensives",
        "delta_state", "feed,T,P", "state_tiny=state_small", "feed");
    CHECK(tiny_intensive_difference <= 1.0e-8);
}

void check_c1_ideal_caller_chart() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto problem = equilibrium_support::flash(2);
    problem.feed.values = {0.2, 0.9};
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    auto point = compiled.generic_point();
    CHECK(point.size() >= 2);
    if (point.size() < 2) return;
    point[0] = std::log(0.2);
    point[1] = std::log(0.9);
    const auto base = compiled.evaluate(point);
    const auto analytic = compiled.objective_gradient(point);
    const auto inputs = [&](std::size_t column, double step) {
        return recording_support::fields({
            recording_support::field("T", problem.temperature.value),
            recording_support::field("P", problem.pressure.value),
            recording_support::field("feed", problem.feed.values),
            recording_support::field("x", std::vector<double>{0.2, 0.9}),
            recording_support::field("sum_x", 1.1),
            recording_support::field("coordinate", column),
            recording_support::field("step", step),
            recording_support::field("objective", base.objective),
        });
    };
    recording_support::expect_group(
        recording_support::active_test_id(), "c1-equilibrium-caller", "objective-gradient/", 6,
        "off-normalized CompiledSystem objective gradient matches objective value ladder", "coordinate,step");
    for (const std::size_t column : {std::size_t(0), std::size_t(1)}) {
        for (const double step : {1.0e-2, 3.0e-3, 1.0e-3}) {
            auto plus = point;
            auto minus = point;
            plus[column] += step;
            minus[column] -= step;
            const double finite_difference =
                (compiled.evaluate(plus).objective - compiled.evaluate(minus).objective)
                / (2.0 * step);
            const double tolerance = step == 1.0e-2 ? 1.0e-4
                : step == 3.0e-3 ? 1.0e-5 : 1.0e-6;
            const bool verdict = close_with_tolerance(analytic[column], finite_difference, tolerance);
            const std::string point_id = std::to_string(column) + "/" + recording_support::number(step);
            recording_support::compare(
                (std::string("objective-gradient/") + point_id).c_str(), "objective_gradient", "1",
                analytic[column], finite_difference, inputs(column, step), "Relationship",
                "centered_objective_value_ladder", "off-normalized caller gradient matches objective",
                "reference_or_one", 0.0, tolerance, verdict, "c1-equilibrium-caller", point_id.c_str(),
                "dG/dz", "T,P,feed,x,sum_x,coordinate,step", "dG/dz=\u0394G/(2\u0394z)", "coordinate,step");
            CHECK(verdict);
        }
    }
}

void check_s8() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto small = equilibrium_support::flash(2);
    auto large = small;
    large.feed.values = {5.0, 5.0};
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, small);
    const auto second = epcsaft::equilibrium::solve_equilibrium(model, large);
    CHECK(first.success && second.success);
    if (!first.success || !second.success) return;
    const bool same_intensives = first.temperature == second.temperature
        && first.pressure == second.pressure;
    const bool same_amount_shape = first.amounts.size() == second.amounts.size();
    if (!same_amount_shape) {
        CHECK(same_intensives);
        CHECK(same_amount_shape);
        return;
    }
    for (std::size_t i = 0; i < first.amounts.size(); ++i) {
        const double ratio = second.amounts[i] / first.amounts[i];
        recording_support::compare(
            (std::string("amount-scale/") + std::to_string(i)).c_str(), "phase_amount_ratio", "1",
            ratio, 10.0,
            recording_support::fields({recording_support::field("T", first.temperature.front()),
                                       recording_support::field("P", first.pressure.front()),
                                       recording_support::field("feed_small", small.feed.values),
                                       recording_support::field("feed_large", large.feed.values),
                                       recording_support::field("phase", i)}),
            "Relationship", "extensive-scaling", "amounts scale with feed by ten",
            "strict_absolute", 1.0e-7, 0.0, std::abs(ratio - 10.0) < 1.0e-7,
        "s8-amount-scale", std::to_string(i).c_str(), "n", "T,P,feed,phase",
        "n_{large}/n_{small}=10", "phase");
        CHECK(std::abs(ratio - 10.0) < 1.0e-7);
    }
    CHECK(same_intensives);
    CHECK(same_amount_shape);
}

void check_s10() {
    epcsaft::HelmholtzModel model(properties_support::neutral_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
                      equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    problem.temperature.value = 180.0;
    const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(solved.success);
    if (!solved.success) return;
    const auto flash = epcsaft::equilibrium::compile_problem(model, problem);
    const auto flash_point = solved.coordinates;
    check_finite_difference(flash, flash_point);
    std::vector<double> primal(flash_point.begin(),
                               flash_point.begin() + static_cast<std::ptrdiff_t>(flash.primal_size));
    const auto active_rows = callback_rows(flash);
    std::vector<double> lambda = solved.constraint_multipliers;
    CHECK(lambda.size() == active_rows.size());
    if (lambda.size() != active_rows.size()) return;
    for (std::size_t row = 0; row < lambda.size(); ++row) {
        if (std::abs(lambda[row]) < 1.0e-12)
            lambda[row] = 0.25 + 0.1 * static_cast<double>(row);
    }
    // F5 previous-failure check: vary objective factors and nonzero callback
    // multipliers independently, then differentiate the complete weighted gradient.
    auto pinned_problem = equilibrium_support::boundary();
    pinned_problem.temperature.value = 180.0;
    epcsaft::HelmholtzModel pinned_model(properties_support::neutral_binary());
    const auto pinned = epcsaft::equilibrium::compile_problem(
        pinned_model, pinned_problem
    );
    const auto pinned_point = pinned.generic_point();
    const auto pinned_rows = callback_rows(pinned);
    std::vector<double> pinned_lambda(pinned_rows.size());
    for (std::size_t row = 0; row < pinned_lambda.size(); ++row)
        pinned_lambda[row] = 0.35 + 0.1 * static_cast<double>(row);
    // Publish both complete sweeps before either ladder starts evaluating
    // weighted gradients, so a regular-case failure cannot hide pinned rows.
    expect_weighted_hessian_gradient_ladder(
        flash, primal, lambda, "regular"
    );
    expect_weighted_hessian_gradient_ladder(
        pinned, pinned_point, pinned_lambda, "pinned"
    );
    std::array<std::size_t, 3> regular_false_verdicts{};
    std::array<std::size_t, 3> pinned_false_verdicts{};
    const bool regular_ok = check_weighted_hessian_gradient_ladder(
        flash, primal, lambda, "regular", 1.0e-6, &regular_false_verdicts
    );
    const bool pinned_ok = check_weighted_hessian_gradient_ladder(
        pinned, pinned_point, pinned_lambda, "pinned", 1.0e-5, &pinned_false_verdicts
    );
    // Free intensive values must be differentiated at the callback point,
    // including reaction offsets and simultaneous temperature/pressure freedom.
    bool intensive_curvature_ok = true;
    bool correlated_jacobian_ok = true;
    for (int variant = 0; variant < 6; ++variant) {
        auto request = pinned_problem;
        request.temperature = epcsaft::equilibrium::FreeValue::variable(180.0);
        request.pressure = epcsaft::equilibrium::FreeValue::fixed(100000.0);
        const char* label = "pinned-free-T";
        if (variant == 1) {
            request.pressure = epcsaft::equilibrium::FreeValue::variable(100000.0);
            request.phases[1].amount = epcsaft::equilibrium::AmountSpec::pinned(0.2);
            auto third_phase = request.phases[1];
            third_phase.name = "L2";
            third_phase.kind = epcsaft::equilibrium::PhaseKind::Liquid;
            third_phase.amount = epcsaft::equilibrium::AmountSpec::pinned(0.3);
            request.phases.push_back(third_phase);
            label = "pinned-free-TP";
        } else if (variant == 2) {
            request.reactions = {{"A-to-B", {-1.0, 1.0}, std::nullopt, 850.0}};
            label = "pinned-reactive-T";
        } else if (variant == 5) {
            // Non-affine correlated standard state: its exact s', s'' enter the T column and
            // curvature; the callback point below evaluates at 313 K inside [250, 450] K.
            request.temperature.value = 300.0;
            request.reactions = {{"A-to-B", {-1.0, 1.0}, std::nullopt, std::nullopt, legacy_isomerisation(1.0, moderate_isomerisation_a())}};
            label = "pinned-correlated-T";
        } else if (variant >= 3) {
            request.temperature.value = 300.0;
            label = variant == 3 ? "pinned-association-T" : "pinned-ionic-T";
            if (variant == 4) request.feed.values = {0.98, 0.01, 0.01};
        }
        epcsaft::HelmholtzModel action_model(variant == 3
            ? properties_support::associating_binary() : variant == 4
            ? properties_support::aqueous_electrolyte() : properties_support::neutral_binary());
        const auto compiled = epcsaft::equilibrium::compile_problem(action_model, request);
        auto callback_point = compiled.generic_point();
        callback_point[compiled.primal_size - compiled.free_intensive_count] = variant >= 3 ? 313.0 : 193.0;
        if (variant == 1) callback_point[compiled.primal_size - 1] = 120000.0;
        std::vector<double> weights(callback_rows(compiled).size(), 0.0);
        for (std::size_t row = 0; row < weights.size(); ++row)
            weights[row] = 0.35 + 0.1 * static_cast<double>(row);
        expect_weighted_hessian_gradient_ladder(compiled, callback_point, weights, label);
        intensive_curvature_ok = check_weighted_hessian_gradient_ladder(
            compiled, callback_point, weights, label, 1.0e-5) && intensive_curvature_ok;
        if (variant == 5) {
            // The pinned ladder criterion (1e-5) is looser than the affine-shortcut defect at this
            // point (~2e-6), so the T,T curvature entry is compared directly at 1e-9 scaled: the
            // falsifier flips if the retired s'' = -2 s'/T shortcut returns.
            epcsaft::equilibrium::GibbsMinimizer nlp(compiled, callback_point);
            const auto shape = nlp.shape();
            const std::size_t t_index = compiled.primal_size - compiled.free_intensive_count;
            std::vector<Ipopt::Number> hessian(static_cast<std::size_t>(shape.hessian_entries), 0.0);
            CHECK(nlp.eval_h(shape.variables, callback_point.data(), true, 0.0, shape.constraints,
                             weights.data(), true, shape.hessian_entries, nullptr, nullptr, hessian.data()));
            const double h_tt = hessian[t_index * (t_index + 1) / 2 + t_index];
            const double tt_step = centered_step(callback_point[t_index]);
            auto plus = callback_point, minus = callback_point;
            plus[t_index] += tt_step; minus[t_index] -= tt_step;
            const double h_tt_centered = (weighted_gradient(compiled, plus, weights, 0.0)[t_index]
                - weighted_gradient(compiled, minus, weights, 0.0)[t_index]) / (2.0 * tt_step);
            const bool tt_ok = close_with_tolerance(h_tt, h_tt_centered, 1.0e-9);
            recording_support::compare("correlated-jacobian/hessian-TT", "weighted_lagrangian_hessian", "coordinate",
                h_tt, h_tt_centered, recording_support::fields({weighted_hessian_context(compiled),
                    recording_support::field("point", callback_point), recording_support::field("step", tt_step)}),
                "Relationship", "weighted_hessian_gradient_finite_difference",
                "correlated standard-state T,T curvature is exact beyond the affine shortcut",
                "reference_or_one", 0.0, 1.0e-9, tt_ok, "f5-correlated-jacobian", "hessian-TT",
                "\\nabla^2\\mathcal L_{TT}", "\\mathbf{z},\\lambda,\\Delta T", "exact s''", "entry");
            correlated_jacobian_ok = correlated_jacobian_ok && tt_ok;
            // Rows, Jacobian temperature column and lever column against centered differences.
            std::vector<double> full = callback_point;
            full.resize(compiled.unknown_count(), 0.0);
            const auto jacobian = compiled.jacobian(full);
            const auto levers = compiled.lever_jacobian(full);
            const std::size_t lever_width = levers.size() / compiled.rows.size();
            const std::size_t column = compiled.primal_size - compiled.free_intensive_count;
            const double step = 1.0e-3;
            auto warmer = full, colder = full;
            warmer[column] += step; colder[column] -= step;
            const auto warm = compiled.evaluate(warmer).values, cold = compiled.evaluate(colder).values;
            for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
                const double centered = (warm[row] - cold[row]) / (2.0 * step);
                const double analytic = jacobian[row * full.size() + column];
                const bool verdict = close(analytic, centered)
                    && (compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Chemical
                        || close(levers[row * lever_width], centered));
                const auto id = std::string("correlated-jacobian/") + std::to_string(row);
                recording_support::compare(id.c_str(), "equilibrium_jacobian", "1", analytic, centered,
                    recording_support::fields({weighted_hessian_context(compiled),
                        recording_support::field("point", full), recording_support::field("row", row),
                        recording_support::field("lever_column", levers[row * lever_width]),
                        recording_support::field("delta_T", step)}),
                    "Relationship", "finite-difference",
                    "correlated standard-state temperature column matches centered rows",
                    "reference_or_one", 0.0, 1.0e-6, verdict, "f5-correlated-jacobian", id.c_str(),
                    "J_{r,T}", "T,P,\\mathbf{x},\\Delta T", "J_{r,T}=\\Delta r/\\Delta T", "row");
                correlated_jacobian_ok = correlated_jacobian_ok && verdict;
            }
        }
    }
    const std::array<double, 3> objective_factors = {0.0, 0.5, 1.0};
    for (std::size_t i = 0; i < objective_factors.size(); ++i) {
        const auto inputs = recording_support::fields({
            recording_support::field("objective_factor", objective_factors[i]),
            recording_support::field("regular_false_verdict_count", regular_false_verdicts[i]),
            recording_support::field("pinned_false_verdict_count", pinned_false_verdicts[i]),
        });
        recording_support::compare(
            (std::string("f5-summary/regular/")
                + recording_support::number(objective_factors[i])).c_str(),
            "false_verdict_count", "count", static_cast<double>(regular_false_verdicts[i]),
            0.0, inputs, "Relationship", "hessian-ladder-summary",
            "regular Hessian ladder has no false verdicts", "absolute", 0.0, 0.0,
            regular_false_verdicts[i] == 0, "f5-hessian-summary", "regular", "count",
            "objective_factor", "false_verdict_count=0", "objective_factor");
        recording_support::compare(
            (std::string("f5-summary/pinned/")
                + recording_support::number(objective_factors[i])).c_str(),
            "false_verdict_count", "count", static_cast<double>(pinned_false_verdicts[i]),
            0.0, inputs, "Relationship", "hessian-ladder-summary",
            "pinned Hessian ladder has no false verdicts", "absolute", 0.0, 0.0,
            pinned_false_verdicts[i] == 0, "f5-hessian-summary", "pinned", "count",
            "objective_factor", "false_verdict_count=0", "objective_factor");
    }
    CHECK(regular_ok);
    CHECK(pinned_ok);
    CHECK(intensive_curvature_ok);
    CHECK(correlated_jacobian_ok);
}

void check_s11() {
    const std::vector<epcsaft::equilibrium::ReactionSpec> reactions = {
        {"isomerisation", {-1.0, 1.0, 0.0}, 2.0, std::nullopt}
    };
    const auto standard = epcsaft::equilibrium::standard_state(300.0, 3, reactions, {});
    const auto affinity = epcsaft::equilibrium::reaction_affinities(reactions, standard);
    const double affinity_reference = -epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0 * std::log(2.0);
    recording_support::compare(
        "standard-affinity", "reaction_affinity", "J/mol", affinity.front(), affinity_reference,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("stoichiometry", std::vector<double>{-1.0, 1.0, 0.0}),
                                   recording_support::field("K", 2.0)}),
        "Relationship", "reaction-standard-state", "standard affinity follows K=2",
        "strict_absolute", 1e-10, 0.0,
        std::abs(affinity.front() - affinity_reference) < 1e-10,
        "s11-standard-affinity", "reaction", "A", "T,\\nu,K", "A=-RT\\ln K", "reaction");
    CHECK(affinity.size() == 1);
    CHECK(std::isfinite(affinity.front()));
    CHECK(std::abs(affinity.front() + epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0
                   * std::log(2.0)) < 1.0e-10);
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto problem = equilibrium_support::flash(3);
    problem.phases.front().kind = epcsaft::equilibrium::PhaseKind::Vapor;
    problem.phases.front().composition_guess = {2.0 / 9.0, 4.0 / 9.0, 1.0 / 3.0};
    problem.reactions = reactions;
    double previous = std::numeric_limits<double>::infinity();
    for (const double pressure : {1.0e4, 1.0e3, 1.0e2}) {
        problem.pressure.value = pressure;
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
        CHECK(result.success && result.mole_fractions.size() == 3);
        if (!result.success || result.mole_fractions.size() != 3) {
            std::fprintf(stderr, "S11 solve: %s\n", result.message.c_str());
            return;
        }
        const double gap = std::abs(std::log(result.mole_fractions[1] / result.mole_fractions[0])
            + affinity.front() / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * result.temperature.front()));
        const auto pressure_id = recording_support::number(pressure);
        const auto previous_gap = std::isfinite(previous)
            ? recording_support::number(previous) : std::string("none");
        const auto gap_inputs = recording_support::fields({
            recording_support::field("T", result.temperature.front()),
            recording_support::field("P", pressure),
            recording_support::field("x", result.mole_fractions),
            recording_support::field("affinity", affinity.front()),
            recording_support::field("previous_gap", previous_gap),
        });
        const bool gap_contracts = std::isfinite(gap) && std::isfinite(previous)
            && gap <= previous;
        if (std::isfinite(previous)) {
            recording_support::compare_custom(
                (std::string("pressure-gap/") + pressure_id).c_str(),
                "reaction_equilibrium_gap", "1", gap, previous, gap_inputs,
                "Relationship", "reaction-equilibrium",
                "reaction gap contracts as pressure decreases",
                "isfinite(model) && isfinite(reference) && model <= reference", gap_contracts,
                "s11-gap", pressure_id.c_str(), "\\Delta", "T,P,\\mathbf{x},A,\\Delta_{previous}",
                "\\Delta(P)\\leq\\Delta(P_{previous})", "pressure");
        } else {
            recording_support::observe(
                (std::string("pressure-gap/") + pressure_id).c_str(),
                "reaction_equilibrium_gap", "1", gap, std::nullopt, gap_inputs,
                "Relationship", "reaction-equilibrium",
                "first reaction gap establishes the contraction baseline",
                "s11-gap", pressure_id.c_str(), "\\Delta", "T,P,\\mathbf{x},A",
                "\\Delta(P_{first})", "pressure");
        }
        CHECK(std::isfinite(gap) && (std::isinf(previous) || gap_contracts));
        previous = gap;
        if (pressure == 1.0e2) {
            recording_support::compare(
                "final-gap", "reaction_equilibrium_gap", "1", gap, 0.0,
                recording_support::fields({recording_support::field("T", result.temperature.front()),
                                           recording_support::field("P", pressure),
                                           recording_support::field("x", result.mole_fractions),
                                           recording_support::field("affinity", affinity.front())}),
                "Relationship", "reaction-equilibrium", "final reaction gap is below one millith",
                "reference_or_one", 1.0e-3, 0.0, gap <= 1.0e-3,
                "s11-final-gap", "final", "\\Delta", "T,P,\\mathbf{x},A", "\\Delta\\leq10^{-3}", "pressure");
            CHECK(gap <= 1.0e-3);
        }
    }
}

void check_s13() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto problem = equilibrium_support::flash(2);
    const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
    const double limit = result.convergence.criterion == nullptr
        ? 0.0 : convergence_limit(result.convergence);
    recording_support::compare(
        "phase-count", "phase_count", "count", static_cast<double>(result.phase_count), 1.0,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values),
                                   recording_support::field("limit", limit)}),
        "Relationship", "single-phase", "single-phase solve returns one phase",
        "reference_or_one", 0.0, 0.0, result.phase_count == 1,
        "s13-phase-count", "phase", "N_{phase}", "T,P,feed", "N_{phase}=1", "phase");
    CHECK(result.success && result.phase_count == 1);
    CHECK(result.mole_fractions.size() == 2);
    recording_support::compare(
        "composition-0", "mole_fraction", "1", result.mole_fractions[0], problem.feed.values[0],
        recording_support::fields({recording_support::field("T", 300.0), recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values), recording_support::field("limit", limit),
                                   recording_support::field("component", 0U)}),
        "Relationship", "feed-balance", "phase composition follows feed", "reference_or_one", limit, 0.0,
        std::abs(result.mole_fractions[0] - problem.feed.values[0]) <= limit,
        "s13-composition", "component-0", "x_0", "T,P,feed,component", "x_0=x_{feed,0}", "component");
    recording_support::compare(
        "composition-1", "mole_fraction", "1", result.mole_fractions[1], problem.feed.values[1],
        recording_support::fields({recording_support::field("T", 300.0), recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values), recording_support::field("limit", limit),
                                   recording_support::field("component", 1U)}),
        "Relationship", "feed-balance", "phase composition follows feed", "reference_or_one", limit, 0.0,
        std::abs(result.mole_fractions[1] - problem.feed.values[1]) <= limit,
        "s13-composition", "component-1", "x_1", "T,P,feed,component", "x_1=x_{feed,1}", "component");
    CHECK(std::abs(result.mole_fractions[0] - problem.feed.values[0]) <= limit);
    CHECK(std::abs(result.mole_fractions[1] - problem.feed.values[1]) <= limit);
    CHECK(result.molar_densities.size() == 1 && result.molar_densities.front() > 0.0);
    const auto root = epcsaft::solve_density(
        model, 300.0, 1.0e5, problem.feed.values, epcsaft::DensityBranch::Liquid, std::nullopt);
    const double density_limit = std::max(limit, convergence_limit(root.convergence));
    const double specific_volume_difference =
        std::abs(1.0 / result.molar_densities.front() - 1.0 / root.molar_density);
    recording_support::compare(
        "density-root", "specific_volume_difference", "m3/mol", specific_volume_difference, 0.0,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values),
                                   recording_support::field("density_limit", density_limit),
                                   recording_support::field("rho_result", result.molar_densities.front()),
                                   recording_support::field("rho_root", root.molar_density)}),
        "Relationship", "density-root", "equilibrium density matches direct density root",
        "reference_or_one", density_limit, 0.0, specific_volume_difference <= density_limit,
        "s13-density-root", "liquid", "\\Delta v", "T,P,feed,branch", "\\Delta v\\leq\\ell", "density");
    CHECK(specific_volume_difference <= density_limit);
}

void check_s15() {
    epcsaft::HelmholtzModel model(properties_support::aqueous_electrolyte());
    auto problem = equilibrium_support::flash(3);
    problem.phases = {
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
    };
    problem.phases[1].support = {0};
    problem.temperature.value = 298.15;
    problem.pressure.value = 1.0e5;
    problem.feed.kind = epcsaft::equilibrium::FeedKind::Molalities;
    problem.feed.values = {0.0, 1.0, 1.0};
    problem.feed.solvent = {0};
    problem.feed.solvent_mass_kg = 1.0;
    const auto compiled = epcsaft::equilibrium::compile_problem(model, problem);
    epcsaft::core::Capture deficit_capture;
    const auto failure = epcsaft::equilibrium::solve_equilibrium(model, problem);
    recording_support::compare(
        "vanished-phase", "success", "bool", failure.success ? 1.0 : 0.0, 0.0,
        recording_support::fields({recording_support::field("T", 298.15),
                                   recording_support::field("P", 1.0e5),
                                   recording_support::field("feed", problem.feed.values),
                                   recording_support::field("vanished_phase", "V")}),
        "Relationship", "boundary-failure", "vanished vapor phase is rejected",
        "reference_or_one", 0.0, 0.0, !failure.success, "s15-vanished-phase", "V",
        "success", "T,P,feed,phase", "success=0", "phase");
    const bool failure_status_ok = !failure.success && !failure.has_state();
    const bool failure_message_ok = failure.message.find("phase V vanished") != std::string::npos;
    double expected_deficit = 1.0;
    bool driving_row = false;
    const double limit = failure.convergence.criterion == nullptr
        ? 0.0 : convergence_limit(failure.convergence);
    for (std::size_t row = 0; row < failure.rows.size(); ++row)
        if (failure.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical
            && failure.rows[row].phase == 1) {
            expected_deficit -= std::exp(-failure.residuals[row]);
            driving_row = driving_row || std::abs(failure.residuals[row]) > limit;
        }
    const bool exact_failure_message = failure.message
        == "phase V vanished; Michelsen deficit " + std::to_string(expected_deficit);
    const auto& trace = deficit_capture.trace();
    const bool trace_ok = !trace.spans.empty()
        && std::string(trace.spans.front().name) == "solve_equilibrium";
    const auto deficit = std::find_if(trace.events.rbegin(), trace.events.rend(), [](const auto& event) {
            return event.span == 0 && std::string(event.name) == "michelsen_deficit";
        });
    const auto distance = std::find_if(trace.events.rbegin(), trace.events.rend(), [](const auto& event) {
            return event.span == 0 && std::string(event.name) == "phase_distance";
        });
    const bool deficit_present = deficit != trace.events.rend() && deficit->values.size() == 1;
    const bool distance_present = distance != trace.events.rend() && distance->values.size() == 1;
    bool deficit_ok = false;
    bool distance_ok = false;
    if (deficit_present) {
        deficit_ok = std::abs(deficit->values.front() - expected_deficit)
            <= 1.0e-12 * std::max(1.0, std::abs(expected_deficit));
        recording_support::compare(
            "deficit", "michelsen_deficit", "1", deficit->values.front(), expected_deficit,
            recording_support::fields({recording_support::field("T", 298.15),
                                       recording_support::field("P", 1.0e5),
                                       recording_support::field("feed", problem.feed.values),
                                       recording_support::field("expected_deficit", expected_deficit)}),
            "Relationship", "Michelsen-deficit", "trace deficit equals recomputed phase deficit",
            "reference_or_one", 0.0, 1e-12,
            deficit_ok,
            "s15-deficit", "trace", "D", "T,P,feed", "D=1-\\sum_i\\exp(-r_i)", "deficit");
    }
    if (distance_present) {
        const double distance_reference = compiled.feed_composition[1]
            + compiled.feed_composition[2];
        distance_ok = std::isfinite(distance->values.front())
            && std::isfinite(distance_reference)
            && distance->values.front() >= distance_reference;
        recording_support::compare_custom(
            "distance", "phase_distance", "1", distance->values.front(), distance_reference,
            recording_support::fields({recording_support::field("T", 298.15),
                                       recording_support::field("P", 1.0e5),
                                       recording_support::field("feed", problem.feed.values),
                                       recording_support::field("reference_distance", distance_reference)}),
            "Relationship", "boundary-failure", "phase distance reaches the vanished-phase feed amount",
            "isfinite(model) && isfinite(reference) && model >= reference", distance_ok,
            "s15-distance", "trace", "d", "T,P,feed", "d\\geq n_{vanished}", "trace");
    }
    CHECK(failure_status_ok);
    CHECK(failure_message_ok);
    CHECK(exact_failure_message);
    CHECK(trace_ok);
    CHECK(deficit_present);
    CHECK(distance_present);
    CHECK(deficit_ok);
    CHECK(distance_ok);
    CHECK(driving_row);
    auto symmetric = equilibrium_support::synthetic_reactive_problem();
    symmetric.phases[0].composition_guess.clear();
    symmetric.phases[1].composition_guess.clear();
    bool refused = false;
    try { (void)epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::synthetic_reactive_ternary()), symmetric); }
    catch (const std::runtime_error& error) {
        refused = std::string(error.what()).find("start_symmetric") != std::string::npos;
    }
    CHECK(refused);
}

void check_c3_availability() {
    epcsaft::HelmholtzModel failure_model(properties_support::aqueous_electrolyte());
    auto failure_problem = equilibrium_support::flash(3);
    failure_problem.phases = {
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor),
    };
    failure_problem.phases[1].support = {0};
    failure_problem.temperature.value = 298.15;
    failure_problem.pressure.value = 1.0e5;
    failure_problem.feed.kind = epcsaft::equilibrium::FeedKind::Molalities;
    failure_problem.feed.values = {0.0, 1.0, 1.0};
    failure_problem.feed.solvent = {0};
    failure_problem.feed.solvent_mass_kg = 1.0;
    const auto failure = epcsaft::equilibrium::solve_equilibrium(failure_model, failure_problem);
    const bool retained = !failure.solver_status.empty()
        && failure.convergence.criterion != nullptr
        && std::isfinite(failure.convergence.residual)
        && failure.rows.size() == failure.residuals.size()
        && std::all_of(failure.residuals.begin(), failure.residuals.end(),
                       [](double value) { return std::isfinite(value); })
        && !failure.message.empty();
    const bool cleared = !failure.success && !failure.has_state()
        && failure.temperature.empty() && failure.pressure.empty()
        && failure.amounts.empty() && failure.mole_fractions.empty()
        && failure.molar_densities.empty() && failure.coordinates.empty();
    const bool failure_reason = failure.message.find("phase V vanished") != std::string::npos;
    const bool tiny_step_stalled = failure.solver_status == "STOP_AT_TINY_STEP"
        && failure.convergence.stop == epcsaft::core::Stop::stalled;

    epcsaft::HelmholtzModel success_model(properties_support::neutral_binary());
    auto success_problem = equilibrium_support::flash(2);
    success_problem.reactions.clear();
    success_problem.phases.resize(1);
    const auto success = epcsaft::equilibrium::solve_equilibrium(success_model, success_problem);
    const bool success_state = success.success && success.has_state()
        && success.convergence.criterion != nullptr
        && std::isfinite(success.convergence.residual)
        && !success.solver_status.empty();

    recording_support::expect_group(
        recording_support::active_test_id(), "c3-availability", "c3-availability", 2,
        "C3 equilibrium availability", "check");
    recording_support::compare(
        "c3-availability", "failed-result", "1",
        (retained && cleared && failure_reason && tiny_step_stalled) ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("solver_status", failure.solver_status),
            recording_support::field("stop", epcsaft::core::stop_name(failure.convergence.stop)),
            recording_support::field("residual", failure.convergence.residual),
            recording_support::field("row_count", failure.rows.size()),
            recording_support::field("residual_count", failure.residuals.size()),
            recording_support::field("message", failure.message),
            recording_support::field("T", failure_problem.temperature.value),
            recording_support::field("P", failure_problem.pressure.value),
            recording_support::field("feed", failure_problem.feed.values),
        }),
        "Relationship", "unavailable-result", "failed solve retains diagnostics and clears state",
        "absolute", 0.0, 0.0,
        retained && cleared && failure_reason && tiny_step_stalled, "native.equilibrium.c3", "C3",
        "failed-result", "availability", "T,P,feed", "has_state=0 with diagnostics retained");
    recording_support::compare(
        "c3-availability", "success-state", "1", success_state ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("solver_status", success.solver_status),
            recording_support::field("stop", epcsaft::core::stop_name(success.convergence.stop)),
            recording_support::field("residual", success.convergence.residual),
            recording_support::field("T", success_problem.temperature.value),
            recording_support::field("P", success_problem.pressure.value),
            recording_support::field("feed", success_problem.feed.values),
        }),
        "Relationship", "available-result", "successful solve exposes state diagnostics",
        "absolute", 0.0, 0.0, success_state, "native.equilibrium.c3", "C3",
        "success-state", "availability", "T,P,feed", "has_state=1 with convergence retained");
    CHECK(retained && cleared && failure_reason && tiny_step_stalled);
    CHECK(success_state);
}

void check_s16() {
    epcsaft::HelmholtzModel model(equilibrium_support::synthetic_reactive_ternary());
    auto problem = equilibrium_support::synthetic_reactive_problem();
    const auto first = epcsaft::equilibrium::solve_equilibrium(model, problem);
    CHECK(first.success);
    if (!first.success) return;
    for (const double fraction : {0.25, 0.5, 0.75}) {
        for (std::size_t i = 0; i < 3; ++i)
            problem.feed.values[i] = (1.0 - fraction) * first.mole_fractions[i]
                + fraction * first.mole_fractions[3 + i];
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem, first);
        CHECK(result.success);
        if (!result.success) continue;
        const double limit = convergence_limit(result.convergence);
        for (std::size_t i = 0; i < 6; ++i) {
            const double difference = std::abs(
                result.mole_fractions[i] - first.mole_fractions[i]);
            const auto point = recording_support::number(fraction) + "/" + std::to_string(i);
            recording_support::compare(
                (std::string("continuation/") + point).c_str(),
                "mole_fraction_difference", "1", difference, 0.0,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("feed", problem.feed.values),
                                           recording_support::field("fraction", fraction),
                                           recording_support::field("component", i),
                                           recording_support::field("limit", limit),
                                           recording_support::field("solution", result.mole_fractions)}),
                "Relationship", "continuation", "continuation composition stays within convergence limit",
                "reference_or_one", limit, 0.0, difference <= limit,
                "s16-composition", point.c_str(), "\\Delta x_i",
                "T,P,feed,\\alpha,component", "\\Delta x_i\\leq\\ell", "component");
            CHECK(difference <= limit);
        }
        for (std::size_t phase = 0; phase < 2; ++phase) {
            double affinity = 0.0;
            for (std::size_t row = 0; row < result.rows.size(); ++row)
                if (result.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical
                    && result.rows[row].phase == phase)
                    affinity += problem.reactions[0].stoichiometry[result.rows[row].index]
                        * result.residuals[row];
            recording_support::compare(
                (std::string("affinity/") + recording_support::number(fraction) + "/"
                 + std::to_string(phase)).c_str(),
                "reaction_affinity", "1", affinity, 0.0,
                recording_support::fields({recording_support::field("T", problem.temperature.value),
                                           recording_support::field("P", problem.pressure.value),
                                           recording_support::field("feed", problem.feed.values),
                                           recording_support::field("fraction", fraction),
                                           recording_support::field("phase", phase),
                                           recording_support::field("limit", limit)}),
                "Relationship", "reaction-equilibrium", "continuation reaction affinity is zero",
                "reference_or_one", limit, 0.0, std::abs(affinity) <= limit,
                "s16-affinity", (std::to_string(fraction) + "-" + std::to_string(phase)).c_str(),
                "A", "T,P,feed,\\alpha,phase", "A=0", "fraction,phase");
            CHECK(std::abs(affinity) <= limit);
        }
    }
}

// ------------------------------------------------------------------ E2 solves

namespace e2 {

using epcsaft::equilibrium::AmountSpec;
using epcsaft::equilibrium::FreeValue;
using epcsaft::equilibrium::PhaseKind;

bool record_bool(const char* id, const char* group, const char* claim, bool verdict,
                 const std::string& inputs) {
    recording_support::compare(id, "verdict", "bool", verdict ? 1.0 : 0.0, 1.0, inputs,
        "Relationship", "equilibrium", claim, "absolute", 0.0, 0.0, verdict, group, id,
        "verdict", "inputs", claim, "case");
    return verdict;
}

double largest_row(const epcsaft::equilibrium::EquilibriumResult& result) {
    return result.residual_norm();
}

std::string solve_inputs(const epcsaft::equilibrium::EquilibriumProblem& problem,
                         const epcsaft::equilibrium::EquilibriumResult& result) {
    return recording_support::fields({
        recording_support::field("feed", problem.feed.values),
        recording_support::field("T", problem.temperature.value),
        recording_support::field("P", problem.pressure.value),
        recording_support::field("solver_status", result.solver_status),
        recording_support::field("stop", epcsaft::core::stop_name(result.convergence.stop)),
        recording_support::field("residual", result.convergence.residual),
        recording_support::field("x", result.mole_fractions),
        recording_support::field("amounts", result.amounts),
        recording_support::field("rho", result.molar_densities)});
}

bool not_run_record(const epcsaft::equilibrium::EquilibriumResult& result,
                    const epcsaft::equilibrium::CompiledSystem& compiled) {
    return !result.success && result.solver_status == "NOT_RUN"
        && result.message == "rank-deficient declared equations; nonlinear solve not run"
        && result.rows.size() == compiled.rows.size()
        && result.unknown_count == compiled.unknown_count()
        && result.jacobian_rank == compiled.jacobian_rank
        && result.generic_nullity == compiled.generic_nullity && result.generic_nullity > 0
        && result.residuals.empty() && result.iterations == 0 && result.evaluations == 0
        && result.hessian_calls == 0 && result.convergence.stop == epcsaft::core::Stop::infeasible
        && std::isnan(result.convergence.residual) && result.convergence.criterion == nullptr
        && !result.requested_tolerance_met.has_value() && !result.convergence.floor.has_value()
        && !result.convergence.amplification.has_value() && !result.has_state()
        && result.coordinates.empty() && result.multipliers.empty();
}

epcsaft::Parameters identical_ethane_pair() {
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

epcsaft::ThermochemistryRecord species_record(double cp_a, double cp_b) {
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

void check_e2_thermal() {
    using namespace e2;
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    epcsaft::HelmholtzModel pair(identical_ethane_pair());
    bool isomerisation_ok = true;
    // Qualification rows use the moderate offset; the dilute legacy offset is retained below as
    // a documented limitation of the unpinned Gibbs route (owner decision, 2026-09-16).
    struct Offset final { double a; bool qualified; };
    for (const Offset offset : {Offset{moderate_isomerisation_a(), true}, Offset{DILUTE_ISOMERISATION_A, false}}) {
    for (const double sign : offset.qualified ? std::vector<double>{1.0, -1.0} : std::vector<double>{1.0}) {
        for (const double temperature : {330.0, 350.0, 370.0}) {
            auto problem = equilibrium_support::flash(2);
            problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
            problem.feed.values = {0.5, 0.5};
            problem.temperature = FreeValue::fixed(temperature);
            problem.pressure = FreeValue::fixed(2.0e5);
            problem.reactions = {{"A=B", {-sign, sign}, std::nullopt, std::nullopt, legacy_isomerisation(sign, offset.a)}};
            const auto result = epcsaft::equilibrium::solve_equilibrium(pair, problem);
            const double ln_k = offset.a + 1200.0 / temperature + 0.4 * std::log(temperature) + 0.002 * temperature;
            const double expected = std::exp(ln_k) / (1.0 + std::exp(ln_k));
            const double observed = result.success ? result.mole_fractions[1] : std::numeric_limits<double>::quiet_NaN();
            const auto ethane_state = pair.state_at_pressure(temperature, 2.0e5, {1.0, 0.0}, epcsaft::DensityBranch::Vapor);
            const double density_gap = result.success
                ? std::abs(result.molar_densities[0] / ethane_state.molar_density - 1.0)
                : std::numeric_limits<double>::quiet_NaN();
            const auto id = std::string("e2-2/isomerisation/") + (offset.qualified ? "" : "dilute-limit/")
                + (sign > 0.0 ? "forward/" : "reversed/") + recording_support::number(temperature);
            // Unpinned chemical rows are recovered from the KKT gradient N x_i mu_i, so the
            // dilute participant's row is bounded only by row_tolerance / (N x_A).
            const double dilute_bound = result.convergence.criterion->value / (1.0 * (1.0 - expected));
            const auto inputs = recording_support::fields({solve_inputs(problem, result),
                recording_support::field("a", offset.a), recording_support::field("ln_K", ln_k),
                recording_support::field("x_B_expected", expected),
                recording_support::field("rho_ethane", ethane_state.molar_density),
                recording_support::field("density_gap", density_gap),
                recording_support::field("dilute_row_bound", dilute_bound)});
            const bool verdict = result.success && std::abs(observed - expected) <= 1.0e-10
                && largest_row(result) <= result.convergence.criterion->value && density_gap <= 1.0e-9;
            recording_support::compare(id.c_str(), "mole_fraction", "1", observed, expected, inputs,
                "Relationship", "reaction-equilibrium", offset.qualified
                    ? "identical-species isomerisation follows x_B=K/(1+K) with the EOS density"
                    : "documented limitation: dilute unpinned participant row exceeds the criterion",
                "strict_absolute", 1.0e-10, 0.0, verdict, "e2-2-thermal", id.c_str(), "x_B", "T,P,K(T)", "x_B=K/(1+K)", "case");
            if (offset.qualified) {
                isomerisation_ok = isomerisation_ok && verdict;
            } else {
                const bool bounded = std::isfinite(result.convergence.residual)
                    && result.convergence.residual <= dilute_bound;
                recording_support::compare_custom((id + "/row-bound").c_str(), "row_residual", "1",
                    result.convergence.residual, dilute_bound, inputs, "Relationship", "conditioning",
                    "dilute recovered chemical row stays within row_tolerance/(N x_A)",
                    "isfinite(model) && model <= reference", bounded, "e2-2-thermal", id.c_str(),
                    "|r|", "T,P,K(T)", "|r|\\leq tol/(N x_A)", "case");
                // The documented limit is asserted; the unavailable x_B row above is not.
                isomerisation_ok = isomerisation_ok && bounded;
            }
        }
    }
    }
    // A -> 2B with matching species thermochemistry and an explicitly preconverted correlation.
    const double tref = 300.0;
    const epcsaft::ReactionLogPolynomial pressure_law{-10.0 / R, 3000.0 / R, 10.0 / R, 0.0, tref, 250.0, 450.0, "pressure-basis-P0"};
    const epcsaft::ReactionLogPolynomial density_law{
        pressure_law.a - std::log(R * tref * 1.0 / 1.0e5), pressure_law.b, pressure_law.c - 1.0, 0.0,
        tref, 250.0, 450.0, epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
    const auto p3_pressure = epcsaft::evaluate_reaction_properties(pressure_law, 350.0);
    const auto p3_density = epcsaft::reaction_pressure_to_density(
        p3_pressure, {-1.0, 2.0}, epcsaft::equilibrium::EOS_STANDARD_STATE_ID, 350.0, 1.0e5, 1.0);
    const auto p3_preconverted = epcsaft::evaluate_reaction_properties(density_law, 350.0);
    const double cp_density = R * (2.0 * 350.0 * p3_density.derivative + 350.0 * 350.0 * p3_density.second_derivative);
    const auto p3_inputs = recording_support::fields({recording_support::field("T", 350.0),
        recording_support::field("a_rho", density_law.a), recording_support::field("b_rho", density_law.b),
        recording_support::field("c_rho", density_law.c), recording_support::field("ln_K_preconverted", p3_preconverted.ln_k)});
    recording_support::compare("e2-2/a-to-2b/ln_K_rho", "ln_K", "1", p3_density.ln_k, 3.5505788308151547, p3_inputs,
        "Relationship", "standard-state-conversion", "P3 pressure-to-density conversion reproduces the preconverted correlation",
        "strict_absolute", 1.0e-12, 0.0, std::abs(p3_density.ln_k - 3.5505788308151547) <= 1.0e-12
            && std::abs(p3_preconverted.ln_k - p3_density.ln_k) <= 1.0e-12,
        "e2-2-thermal", "ln_K_rho", "\\ln K_\\rho", "T", "\\ln K_\\rho=3.5505788308151547", "case");
    recording_support::compare("e2-2/a-to-2b/dh_rho", "reaction_enthalpy", "J/mol", p3_density.enthalpy, -2410.061916353634, p3_inputs,
        "Relationship", "standard-state-conversion", "density-basis reaction enthalpy", "strict_absolute", 1.0e-9, 0.0,
        std::abs(p3_density.enthalpy + 2410.061916353634) <= 1.0e-9, "e2-2-thermal", "dh_rho", "\\Delta h_\\rho", "T", "\\Delta h_\\rho", "case");
    recording_support::compare("e2-2/a-to-2b/dcp_rho", "reaction_heat_capacity", "J/(mol K)", cp_density, 1.6855373818467605, p3_inputs,
        "Relationship", "standard-state-conversion", "density-basis reaction heat capacity from P3 derivatives", "strict_absolute", 1.0e-12, 0.0,
        std::abs(cp_density - 1.6855373818467605) <= 1.0e-12, "e2-2-thermal", "dcp_rho", "\\Delta c_{p,\\rho}", "T", "\\Delta c_{p,\\rho}", "case");
    CHECK(std::abs(p3_density.ln_k - 3.5505788308151547) <= 1.0e-12 && std::abs(p3_preconverted.ln_k - p3_density.ln_k) <= 1.0e-12);
    CHECK(std::abs(p3_density.enthalpy + 2410.061916353634) <= 1.0e-9 && std::abs(cp_density - 1.6855373818467605) <= 1.0e-12);
    epcsaft::HelmholtzModel thermal(identical_ethane_pair(), {}, species_record(30.0, 20.0));
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    const auto pure = ethane.state_at_pressure(350.0, 2.0e5, {1.0}, epcsaft::DensityBranch::Vapor);
    const double compressibility = 2.0e5 / (pure.molar_density * R * 350.0);
    const double phi = std::exp(pure.residual_chemical_potential_over_rt[0]) / compressibility;
    const double a = std::exp(p3_pressure.ln_k) * 1.0e5 / (phi * 2.0e5);
    const double y = 2.0 * a / (a + std::sqrt(a * a + 4.0 * a));
    const double total = 1.5 / (2.0 - y);
    auto doubling = equilibrium_support::flash(2);
    doubling.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    doubling.phases[0].composition_guess = {1.0 - y, y};
    doubling.feed.values = {0.5, 0.5};
    doubling.temperature = FreeValue::fixed(350.0);
    doubling.pressure = FreeValue::fixed(2.0e5);
    doubling.reactions = {{"A=2B", {-1.0, 2.0}, std::nullopt, std::nullopt, density_law}};
    const auto doubling_result = epcsaft::equilibrium::solve_equilibrium(thermal, doubling);
    const double composition_gap = doubling_result.success ? std::abs(doubling_result.mole_fractions[1] - y) : std::numeric_limits<double>::quiet_NaN();
    const double amount_gap = doubling_result.success
        ? std::max(std::abs(doubling_result.amounts[0] * doubling_result.mole_fractions[0] - total * (1.0 - y)),
                   std::abs(doubling_result.amounts[0] * doubling_result.mole_fractions[1] - total * y))
        : std::numeric_limits<double>::quiet_NaN();
    const bool doubling_ok = doubling_result.success && composition_gap <= 1.0e-8 && amount_gap <= 1.0e-8
        && largest_row(doubling_result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value;
    recording_support::compare("e2-2/a-to-2b/x_B", "mole_fraction", "1",
        doubling_result.success ? doubling_result.mole_fractions[1] : std::numeric_limits<double>::quiet_NaN(), y,
        recording_support::fields({solve_inputs(doubling, doubling_result), recording_support::field("Z", compressibility),
            recording_support::field("phi", phi), recording_support::field("A", a), recording_support::field("N_expected", total),
            recording_support::field("amount_gap", amount_gap)}),
        "Relationship", "reaction-equilibrium", "A=2B composition follows the independent pure-state relation y^2/(1-y)=K_P P0/(phi P)",
        "strict_absolute", 1.0e-8, 0.0, doubling_ok, "e2-2-thermal", "a-to-2b", "x_B", "T,P,K_P,\\phi", "y^2/(1-y)=A", "case");
    CHECK(doubling_ok);
    auto unconverted = doubling;
    unconverted.reactions[0].correlation = pressure_law;
    std::string mismatch;
    try { (void)epcsaft::equilibrium::compile_problem(thermal, unconverted); }
    catch (const std::invalid_argument& error) { mismatch = error.what(); }
    CHECK(record_bool("e2-2/a-to-2b/reference-basis-mismatch", "e2-2-thermal", "reference-basis-mismatch",
        mismatch.find("reference-basis-mismatch") != std::string::npos,
        recording_support::fields({recording_support::field("identity", pressure_law.standard_state_id), recording_support::field("message", mismatch)})));
    epcsaft::HelmholtzModel mismatched(identical_ethane_pair(), {}, species_record(30.0, 21.0));
    std::string inconsistency;
    try { (void)epcsaft::equilibrium::compile_problem(mismatched, doubling); }
    catch (const std::invalid_argument& error) { inconsistency = error.what(); }
    CHECK(record_bool("e2-2/a-to-2b/heat-capacity-inconsistency", "e2-2-thermal", "reference inconsistency",
        inconsistency.find("reference inconsistency") != std::string::npos,
        recording_support::fields({recording_support::field("cp_B", 21.0), recording_support::field("message", inconsistency)})));
    // Solved charged reaction with a manufactured consistent standard state at an admitted EOS state.
    const std::vector<double> x_star = {0.98, 0.01, 0.01};
    const std::vector<double> nu = {2.0, -1.0, -1.0};
    const std::vector<double> n_star = {1.96, 0.02, 0.02};
    for (const std::vector<std::size_t>& order : {std::vector<std::size_t>{0, 1, 2}, std::vector<std::size_t>{2, 0, 1}}) {
        epcsaft::HelmholtzModel salt(equilibrium_support::permute_parameters(properties_support::aqueous_electrolyte(), order));
        const auto x = equilibrium_support::permute_vector(x_star, order);
        const auto state = salt.state_at_pressure(298.15, 1.0e5, x, epcsaft::DensityBranch::Liquid);
        const double rt = R * 298.15;
        std::vector<double> gibbs(3);
        for (std::size_t i = 0; i < 3; ++i)
            gibbs[i] = -rt * (std::log(state.molar_density * x[i]) + state.residual_chemical_potential_over_rt[i]);
        const auto stoichiometry = equilibrium_support::permute_vector(nu, order);
        double contraction = 0.0;
        for (std::size_t i = 0; i < 3; ++i) contraction += stoichiometry[i] * gibbs[i];
        auto problem = equilibrium_support::flash(3);
        problem.feed.values = equilibrium_support::permute_vector(std::vector<double>{0.0, 1.0, 1.0}, order);
        problem.temperature = FreeValue::fixed(298.15);
        problem.phases[0].support = {0, 1, 2};
        problem.phases[0].composition_guess = x;
        problem.reactions = {{"2W=Na+Cl", stoichiometry, std::exp(-contraction / rt), std::nullopt}};
        problem.standard_gibbs_j_per_mol = gibbs;
        const auto result = epcsaft::equilibrium::solve_equilibrium(salt, problem);
        const auto expected = equilibrium_support::permute_vector(n_star, order);
        double amount_gap = std::numeric_limits<double>::quiet_NaN(), charge = std::numeric_limits<double>::quiet_NaN();
        if (result.success) {
            amount_gap = 0.0; charge = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                amount_gap = std::max(amount_gap, std::abs(result.amounts[0] * result.mole_fractions[i] - expected[i]));
                charge += salt.parameters().charge_number[i] * result.mole_fractions[i];
            }
        }
        // Rows use the suite's solved-row comparator max(criterion, floor): the association
        // gradient floor limits this ionic state (converged_to_floor), as in S1.
        const double limit = result.convergence.criterion == nullptr ? 0.0 : convergence_limit(result.convergence);
        const bool ok = result.success && amount_gap <= 1.0e-8 && std::abs(charge) <= 1.0e-10 && largest_row(result) <= limit;
        const auto id = std::string("e2-2/charged-solve/") + (order[0] == 0 ? "identity" : "permuted");
        recording_support::compare(id.c_str(), "amount_gap", "mol", ok ? amount_gap : std::numeric_limits<double>::quiet_NaN(), 0.0,
            recording_support::fields({solve_inputs(problem, result), recording_support::field("n_star", expected),
                recording_support::field("g_manufactured", gibbs), recording_support::field("K", std::exp(-contraction / rt)),
                recording_support::field("charge", charge), recording_support::field("row_limit", limit),
                recording_support::field("floor", epcsaft::core::optional_number(result.convergence.floor)),
                recording_support::field("order", std::vector<double>(order.begin(), order.end()))}),
            "Relationship", "reaction-equilibrium", "shared driver reproduces the manufactured charged solution n*",
            "strict_absolute", 1.0e-8, 0.0, ok, "e2-2-thermal", id.c_str(), "\\Delta n", "T,P,feed,g,K", "\\Delta n\\leq10^{-8}", "case");
        CHECK(ok);
    }
    CHECK(isomerisation_ok);
}

void check_e2_solve() {
    using namespace e2;
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    auto one = equilibrium_support::flash(2);
    one.feed.values = {0.5, 0.5};
    one.phases[0].support = {0, 1};
    one.phases[0].amount = AmountSpec::pinned(1.0);
    one.temperature = FreeValue::variable(300.0);
    auto two = equilibrium_support::boundary();
    two.feed.values = {0.5, 0.5};
    two.phases[0].amount = AmountSpec::pinned(0.5);
    two.phases[1].amount = AmountSpec::pinned(0.5);
    two.temperature = FreeValue::variable(300.0);
    for (const auto* item : {&one, &two}) {
        const auto compiled = epcsaft::equilibrium::compile_problem(model, *item);
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, *item);
        const auto id = std::string("e2-6/all-pinned-not-run/") + std::to_string(item->phases.size());
        CHECK(record_bool(id.c_str(), "e2-6-solve", "nonreactive all-pinned declaration returns the exact pre-solve NOT_RUN record",
            compiled.is_square() && compiled.generic_nullity > 0 && not_run_record(result, compiled),
            recording_support::fields({recording_support::field("phases", item->phases.size()),
                recording_support::field("rows", compiled.rows.size()), recording_support::field("unknowns", compiled.unknown_count()),
                recording_support::field("jacobian_rank", compiled.jacobian_rank), recording_support::field("generic_nullity", compiled.generic_nullity),
                recording_support::field("message", result.message)})));
    }
    // Positive pin in both orientations, then the count-deficient fixed-T/fixed-P form.
    epcsaft::HelmholtzModel ternary(properties_support::neutral_ternary());
    auto bubble = equilibrium_support::flash(3);
    bubble.phases = {equilibrium_support::phase("L", PhaseKind::Liquid),
                     equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.0))};
    bubble.temperature.value = 180.0;
    bubble.pressure = FreeValue::variable(1.0e5);
    const auto bubble_result = epcsaft::equilibrium::solve_equilibrium(ternary, bubble);
    CHECK(bubble_result.success);
    for (const bool free_pressure : {true, false}) {
        auto pinned = bubble;
        pinned.phases[1].amount = AmountSpec::pinned(0.01);
        pinned.temperature = free_pressure ? FreeValue::fixed(180.0) : FreeValue::variable(180.0);
        pinned.pressure = free_pressure ? FreeValue::variable(bubble_result.pressure.front())
                                        : FreeValue::fixed(bubble_result.pressure.front());
        const auto compiled = epcsaft::equilibrium::compile_problem(ternary, pinned);
        const auto result = epcsaft::equilibrium::solve_equilibrium(ternary, pinned, bubble_result);
        const auto id = std::string("e2-6/positive-pin/") + (free_pressure ? "free-P" : "free-T");
        const bool ok = result.success && compiled.pinned_count == 1 && compiled.free_intensive_count == 1
            && largest_row(result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value
            && result.amounts[1] == 0.01 * compiled.feed_total;
        recording_support::compare(id.c_str(), "row_residual", "1", largest_row(result), 0.0,
            recording_support::fields({solve_inputs(pinned, result), recording_support::field("pin", 0.01)}),
            "Relationship", "equilibrium-residual", "one positive pin with one free intensive solves through the shared driver",
            "strict_absolute", epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0, ok, "e2-6-solve", id.c_str(),
            "|r|", "T,P,feed,pin", "|r|\\leq10^{-10}", "case");
        CHECK(ok);
    }
    auto overdetermined = bubble;
    overdetermined.phases[1].amount = AmountSpec::pinned(0.01);
    overdetermined.pressure = FreeValue::fixed(bubble_result.pressure.front());
    std::string message;
    try { (void)epcsaft::equilibrium::compile_problem(ternary, overdetermined); }
    catch (const std::invalid_argument& error) { message = error.what(); }
    CHECK(record_bool("e2-6/positive-pin/fixed-T-fixed-P", "e2-6-solve", "positive pin with fixed T and P is a named structural rejection",
        message.find("pinned") != std::string::npos && message.find("free T or P") != std::string::npos,
        recording_support::fields({recording_support::field("message", message)})));
    // E2-7: explicit kind under a renamed phase, the single pure vapor phase, and the caller criterion.
    auto named = equilibrium_support::flash(2);
    auto renamed = named;
    renamed.phases[0].name = "liquid-phase";
    const auto named_result = epcsaft::equilibrium::solve_equilibrium(model, named);
    const auto renamed_result = epcsaft::equilibrium::solve_equilibrium(model, renamed);
    CHECK(record_bool("e2-7/rename", "e2-7-criterion", "renaming an explicitly typed phase changes nothing physical",
        named_result.success && renamed_result.success && named_result.mole_fractions == renamed_result.mole_fractions
            && named_result.molar_densities == renamed_result.molar_densities && named_result.rows.size() == renamed_result.rows.size(),
        solve_inputs(renamed, renamed_result)));
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    auto vapor = equilibrium_support::flash(1);
    vapor.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    vapor.temperature.value = 320.0;
    const auto vapor_result = epcsaft::equilibrium::solve_equilibrium(ethane, vapor);
    const auto vapor_root = epcsaft::solve_density(ethane, 320.0, 1.0e5, {1.0}, epcsaft::DensityBranch::Vapor, std::nullopt);
    recording_support::compare("e2-7/pure-vapor", "row_residual", "1", largest_row(vapor_result), 0.0,
        recording_support::fields({solve_inputs(vapor, vapor_result), recording_support::field("root_density", vapor_root.molar_density),
            recording_support::field("stable_root_count", vapor_root.stable_root_count.value_or(0))}),
        "Relationship", "equilibrium-residual", "a single declared pure phase solves from its unique density root",
        "strict_absolute", epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 0.0,
        vapor_result.success && largest_row(vapor_result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value,
        "e2-7-criterion", "pure-vapor", "|r|", "T,P", "|r|\\leq10^{-10}", "case");
    CHECK(vapor_result.success && largest_row(vapor_result) <= epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    for (const double tolerance : {epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value, 1.0e-8, 1.0e-12}) {
        auto request = equilibrium_support::flash(2);
        request.feed.values = {0.5, 0.5};
        request.row_tolerance = tolerance;
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, request);
        const bool default_request = tolerance == epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value;
        const bool criterion_ok = result.convergence.criterion != nullptr
            && result.convergence.criterion->value == tolerance
            && (default_request ? result.convergence.criterion == &epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL
                                : std::string(result.convergence.criterion->name) == "EQUILIBRIUM_REQUESTED_ROW_RESIDUAL");
        const bool flag_ok = std::isfinite(result.convergence.residual)
            ? (result.requested_tolerance_met.has_value()
               && *result.requested_tolerance_met == (result.convergence.residual <= tolerance))
            : !result.requested_tolerance_met.has_value();
        const bool warning_ok = !result.success || !result.requested_tolerance_met.has_value()
            || *result.requested_tolerance_met || result.message.find("requested row tolerance") != std::string::npos;
        const auto id = std::string("e2-7/tolerance/") + recording_support::number(tolerance);
        recording_support::compare(id.c_str(), "requested_tolerance_met", "bool",
            result.requested_tolerance_met.has_value() ? (*result.requested_tolerance_met ? 1.0 : 0.0) : std::numeric_limits<double>::quiet_NaN(),
            std::isfinite(result.convergence.residual) ? (result.convergence.residual <= tolerance ? 1.0 : 0.0) : std::numeric_limits<double>::quiet_NaN(),
            recording_support::fields({solve_inputs(request, result), recording_support::field("row_tolerance", tolerance),
                recording_support::field("criterion", result.convergence.criterion == nullptr ? "absent" : result.convergence.criterion->name),
                recording_support::field("floor", epcsaft::core::optional_number(result.convergence.floor)),
                recording_support::field("success", result.success ? 1U : 0U), recording_support::field("message", result.message)}),
            "Relationship", "caller-criterion", "requested_tolerance_met is exactly recomputed_residual <= row_tolerance",
            "absolute", 0.0, 0.0, criterion_ok && flag_ok && warning_ok && (!default_request || result.success),
            "e2-7-criterion", id.c_str(), "met", "tolerance,residual", "met=(|r|\\leq tol)", "tolerance");
        CHECK(criterion_ok && flag_ok && warning_ok && (!default_request || result.success));
    }
    // Owned caller criterion survives copies through the trace after the source record is gone,
    // on the existing finite-floor fixture (aqueous electrolyte, association gradient floor).
    std::optional<epcsaft::core::Convergence> copied;
    {
        epcsaft::core::Capture capture;
        epcsaft::HelmholtzModel salt(properties_support::aqueous_electrolyte());
        auto request = equilibrium_support::flash(3);
        request.feed.values = {0.98, 0.01, 0.01};
        request.temperature.value = 298.15;
        request.row_tolerance = 1.0e-8;
        {
            const auto result = epcsaft::equilibrium::solve_equilibrium(salt, request);
            CHECK(result.success && result.convergence.floor.has_value());
        }
        for (const auto& span : capture.trace().spans)
            if (std::string(span.name) == "solve_equilibrium" && span.convergence.has_value()) copied = span.convergence;
    }
    const bool lifetime_ok = copied.has_value() && copied->criterion != nullptr && copied->criterion_owner
        && copied->floor.has_value() && std::isfinite(*copied->floor)
        && copied->criterion == copied->criterion_owner.get() && copied->criterion->value == 1.0e-8
        && std::string(copied->criterion->name) == "EQUILIBRIUM_REQUESTED_ROW_RESIDUAL"
        && epcsaft::core::classify(copied->residual, *copied->criterion, copied->floor, copied->contraction, false,
                                   copied->amplification) == copied->stop;
    CHECK(record_bool("e2-7/owned-criterion-lifetime", "e2-7-criterion",
        "a copied convergence record classifies identically after the producing result is destroyed", lifetime_ok,
        recording_support::fields({recording_support::field("residual", copied ? copied->residual : std::numeric_limits<double>::quiet_NaN()),
            recording_support::field("floor", copied ? epcsaft::core::optional_number(copied->floor) : "absent"),
            recording_support::field("stop", copied ? epcsaft::core::stop_name(copied->stop) : "absent")})));
}

void check_e2_reuse() {
    using namespace e2;
    struct Case final { const char* label; epcsaft::Parameters parameters; epcsaft::equilibrium::EquilibriumProblem problem; bool solved; };
    auto pinned_problem = equilibrium_support::boundary();
    pinned_problem.temperature.value = 180.0;
    auto free_t = pinned_problem;
    free_t.temperature = FreeValue::variable(180.0);
    free_t.pressure = FreeValue::fixed(1.0e5);
    auto reactive = free_t;
    reactive.reactions = {{"A-to-B", {-1.0, 1.0}, std::nullopt, 850.0}};
    auto associating = free_t;
    associating.temperature.value = 300.0;
    auto regular = equilibrium_support::flash(3);
    regular.phases = {equilibrium_support::phase("V", PhaseKind::Vapor), equilibrium_support::phase("L", PhaseKind::Liquid)};
    regular.temperature.value = 180.0;
    const std::vector<Case> cases = {
        {"regular", properties_support::neutral_ternary(), regular, true},
        {"pinned", properties_support::neutral_binary(), pinned_problem, false},
        {"pinned-free-T", properties_support::neutral_binary(), free_t, false},
        {"pinned-reactive-T", properties_support::neutral_binary(), reactive, false},
        {"pinned-association-T", properties_support::associating_binary(), associating, false},
    };
    const double criterion = 1.0e-9;
    for (const auto& item : cases) {
        const epcsaft::HelmholtzModel model(item.parameters);
        const auto compiled = epcsaft::equilibrium::compile_problem(model, item.problem);
        std::vector<double> a;
        if (item.solved) {
            const auto solved = epcsaft::equilibrium::solve_equilibrium(model, item.problem);
            CHECK(solved.success);
            a.assign(solved.coordinates.begin(), solved.coordinates.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size));
        } else {
            a = compiled.generic_point();
            a[compiled.primal_size - compiled.free_intensive_count] = item.problem.temperature.free ? 193.0 : a[compiled.primal_size - compiled.free_intensive_count];
        }
        epcsaft::equilibrium::GibbsMinimizer nlp(compiled, a);
        const auto shape = nlp.shape();
        const auto n = shape.variables;
        const auto m = shape.constraints;
        CHECK(static_cast<std::size_t>(n) == a.size());
        std::vector<double> lambda_one(static_cast<std::size_t>(m)), lambda_two(static_cast<std::size_t>(m));
        for (std::size_t row = 0; row < lambda_one.size(); ++row) {
            lambda_one[row] = 0.35 + 0.1 * static_cast<double>(row);
            lambda_two[row] = -0.2 + 0.07 * static_cast<double>(row);
        }
        auto b = a;
        for (double& value : b) value += 1.0e-3;
        std::vector<Ipopt::Index> rows_index(static_cast<std::size_t>(shape.hessian_entries)), cols_index(static_cast<std::size_t>(shape.hessian_entries));
        CHECK(nlp.eval_h(n, a.data(), true, 1.0, m, lambda_one.data(), true, shape.hessian_entries, rows_index.data(), cols_index.data(), nullptr));
        const auto callback = [&](const std::vector<double>& point, double objective_factor, const std::vector<double>& lambda) {
            Ipopt::Number f = 0.0;
            std::vector<Ipopt::Number> g(static_cast<std::size_t>(m)), grad(static_cast<std::size_t>(n));
            std::vector<Ipopt::Number> jac(static_cast<std::size_t>(shape.jacobian_entries)), hess(static_cast<std::size_t>(shape.hessian_entries));
            CHECK(nlp.eval_f(n, point.data(), true, f));
            CHECK(nlp.eval_g(n, point.data(), false, m, g.data()));
            CHECK(nlp.eval_grad_f(n, point.data(), false, grad.data()));
            CHECK(nlp.eval_jac_g(n, point.data(), false, m, shape.jacobian_entries, nullptr, nullptr, jac.data()));
            CHECK(nlp.eval_h(n, point.data(), false, objective_factor, m, lambda.data(), true, shape.hessian_entries, nullptr, nullptr, hess.data()));
            std::vector<double> flat = {f};
            flat.insert(flat.end(), g.begin(), g.end());
            flat.insert(flat.end(), grad.begin(), grad.end());
            flat.insert(flat.end(), jac.begin(), jac.end());
            flat.insert(flat.end(), hess.begin(), hess.end());
            return flat;
        };
        const auto fresh = [&](const std::vector<double>& point, double objective_factor, const std::vector<double>& lambda) {
            const auto evaluation = compiled.evaluate(point);
            const auto gradient = compiled.objective_gradient(point);
            const auto jacobian = compiled.jacobian(point);
            const auto objective = compiled.hessian(point, {});
            const auto weighted = compiled.hessian(point, lambda);
            const auto active = callback_rows(compiled);
            std::vector<double> flat = {compiled.pinned_count == 0 ? evaluation.objective : 0.0};
            for (const std::size_t row : active) flat.push_back(evaluation.values[row]);
            flat.insert(flat.end(), gradient.begin(), gradient.end());
            for (const std::size_t row : active)
                for (Ipopt::Index column = 0; column < n; ++column) flat.push_back(jacobian[row * static_cast<std::size_t>(n) + static_cast<std::size_t>(column)]);
            for (Ipopt::Index row = 0; row < n; ++row)
                for (Ipopt::Index column = 0; column <= row; ++column) {
                    const std::size_t index = static_cast<std::size_t>(row) * static_cast<std::size_t>(n) + static_cast<std::size_t>(column);
                    flat.push_back(objective_factor * objective[index] + (weighted[index] - objective[index]));
                }
            return flat;
        };
        const auto compare_flat = [&](const char* step, const std::vector<double>& left, const std::vector<double>& right) {
            CHECK(left.size() == right.size());
            double worst = 0.0;
            for (std::size_t i = 0; i < left.size(); ++i)
                worst = std::max(worst, std::abs(left[i] - right[i]) / std::max(1.0, std::abs(right[i])));
            const auto id = std::string("e2-reuse/") + item.label + "/" + step;
            recording_support::compare(id.c_str(), "callback_difference", "1", worst, 0.0,
                recording_support::fields({recording_support::field("case", item.label), recording_support::field("step", step),
                    recording_support::field("entries", left.size()), recording_support::field("point", a)}),
                "Relationship", "current-point-reuse", "cached callback values equal fresh evaluation at the same point",
                "reference_or_one", criterion, 0.0, worst <= criterion, "e2-reuse", id.c_str(), "\\Delta", "z,\\lambda,\\sigma",
                "\\Delta\\leq10^{-9}", "case,step");
            CHECK(worst <= criterion);
        };
        const auto first = callback(a, 0.5, lambda_one);
        compare_flat("a-first", first, fresh(a, 0.5, lambda_one));
        compare_flat("a-changed-multipliers", callback(a, 1.0, lambda_two), fresh(a, 1.0, lambda_two));
        compare_flat("a-zero-objective", callback(a, 0.0, lambda_one), fresh(a, 0.0, lambda_one));
        compare_flat("b", callback(b, 0.5, lambda_one), fresh(b, 0.5, lambda_one));
        compare_flat("a-return", callback(a, 0.5, lambda_one), first);
        const auto& work = nlp.work();
        const bool work_ok = work.prepared_points == 3 && work.row_evaluations == 3 && work.jacobians == 3
            && work.gradients == 3 && work.objective_hessians == 3 && work.constraint_hessians == 5;
        const auto id = std::string("e2-reuse/") + item.label + "/work";
        recording_support::compare(id.c_str(), "prepared_points", "count", static_cast<double>(work.prepared_points), 3.0,
            recording_support::fields({recording_support::field("case", item.label),
                recording_support::field("row_evaluations", work.row_evaluations), recording_support::field("jacobians", work.jacobians),
                recording_support::field("gradients", work.gradients), recording_support::field("objective_hessians", work.objective_hessians),
                recording_support::field("constraint_hessians", work.constraint_hessians),
                recording_support::field("callback_calls", 5U * 5U)}),
            "Relationship", "current-point-reuse", "five callback rounds at three distinct points prepare each point once",
            "absolute", 0.0, 0.0, work_ok, "e2-reuse", id.c_str(), "N_{prepare}", "case", "N_{prepare}=3", "case");
        CHECK(work_ok);
    }
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.equilibrium.f5", "Weighted_equilibrium_curvature", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "hessian/regular/*,hessian/pinned/*,jacobian/*,correlated-jacobian/*,f5-summary/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "F5"}, check_s10, "S10"},
        {{"native.equilibrium.s1", "Equilibrium_S1", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "residual/*,gibbs-dual/*,charge/*,affinity/*,feed-dual/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S1"}, check_s1},
        {{"native.equilibrium.s2", "Equilibrium_S2", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "binary-success,binary-curvature,reactive-success,reactive-curvature", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S2"}, check_s2},
        {{"native.equilibrium.s6", "Equilibrium_S6", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "one-success,two-success,objective-gap", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S6"}, check_s6},
        {{"native.equilibrium.s7", "Equilibrium_S7", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "continuation/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S7", "continuation/*:12"}, check_s7},
        {{"native.equilibrium.s8", "Equilibrium_S8", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "amount-scale/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S8"}, check_s8},
        {{"native.equilibrium.c1", "Core_C1_amount_scale", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "CORE1", "c1-amount-scale/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "C1"}, check_c1_amount_scale},
        {{"native.equilibrium.c1-multipliers", "Core_C1_multiplier_scale", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "CORE1", "c1-multiplier-scale/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "C1"}, check_c1_multiplier_scale},
        {{"native.equilibrium.c1-ideal-caller", "Contributions_C1_ideal_caller_chart", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "C1", "objective-gradient/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "C1", "objective-gradient/*:6"}, check_c1_ideal_caller_chart},
        {{"native.equilibrium.s11", "Equilibrium_S11", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "standard-affinity,pressure-gap/*,final-gap", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S11", "pressure-gap/*:3"}, check_s11},
        {{"native.equilibrium.s13", "Equilibrium_S13", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "phase-count,composition-*,density-root", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S13", "composition-*:2"}, check_s13},
        {{"native.equilibrium.s15", "Equilibrium_S15", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "vanished-phase,deficit,distance", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S15"}, check_s15},
        {{"native.equilibrium.c3", "Core_C3_availability", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "CORE3", "c3-availability", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "C3", "c3-availability:2"}, check_c3_availability},
        {{"native.equilibrium.s16", "Equilibrium_S16", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "continuation/*,affinity/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S16", "continuation/*:18,affinity/*:6"}, check_s16},
        {{"native.equilibrium.e2-thermal", "Equilibrium_E2_thermal", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-2/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E2-2"}, check_e2_thermal},
        {{"native.equilibrium.e2-solve", "Equilibrium_E2_solve", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-6/*,e2-7/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E2-6"}, check_e2_solve},
        {{"native.equilibrium.e2-reuse", "Equilibrium_E2_reuse", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E2", "e2-reuse/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E2-R"}, check_e2_reuse},
    });
}
