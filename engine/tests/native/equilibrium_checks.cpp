#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
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
            + recording_support::number(reaction.equilibrium_constant) + ":nu=";
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
    epcsaft::HelmholtzModel reactive_model(equilibrium_support::ascani_reactive_ternary());
    const auto reactive = equilibrium_support::ascani_reactive_problem();
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
        epcsaft::HelmholtzModel(equilibrium_support::ascani_reactive_ternary()),
        equilibrium_support::ascani_reactive_problem());
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
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto two = equilibrium_support::ascani_reactive_problem();
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
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto problem = equilibrium_support::ascani_reactive_problem();
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
        recording_support::compare_custom(
            (std::string("f5-summary/pinned/")
                + recording_support::number(objective_factors[i])).c_str(),
            "false_verdict_count", "count", static_cast<double>(pinned_false_verdicts[i]),
            0.0, inputs, "Relationship", "hessian-ladder-summary",
            "pinned Hessian ladder false-verdict count is retained", "model > 0",
            pinned_false_verdicts[i] > 0, "f5-hessian-summary", "pinned", "count",
            "objective_factor", "false_verdict_count>0", "objective_factor");
    }
    CHECK(regular_ok);
    CHECK(pinned_ok);
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
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
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
    auto symmetric = equilibrium_support::ascani_reactive_problem();
    symmetric.phases[0].composition_guess.clear();
    symmetric.phases[1].composition_guess.clear();
    bool refused = false;
    try { (void)epcsaft::equilibrium::solve_equilibrium(
        epcsaft::HelmholtzModel(equilibrium_support::ascani_reactive_ternary()), symmetric); }
    catch (const std::runtime_error& error) {
        refused = std::string(error.what()).find("start_symmetric") != std::string::npos;
    }
    CHECK(refused);
}

void check_s16() {
    epcsaft::HelmholtzModel model(equilibrium_support::ascani_reactive_ternary());
    auto problem = equilibrium_support::ascani_reactive_problem();
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

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.equilibrium.f5", "Weighted_equilibrium_curvature", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "hessian/regular/*,hessian/pinned/*,jacobian/*,f5-summary/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "F5"}, check_s10, "S10"},
        {{"native.equilibrium.s1", "Equilibrium_S1", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "residual/*,gibbs-dual/*,charge/*,affinity/*,feed-dual/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S1"}, check_s1},
        {{"native.equilibrium.s2", "Equilibrium_S2", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "binary-success,binary-curvature,reactive-success,reactive-curvature", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S2"}, check_s2},
        {{"native.equilibrium.s6", "Equilibrium_S6", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "one-success,two-success,objective-gap", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S6"}, check_s6},
        {{"native.equilibrium.s7", "Equilibrium_S7", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "continuation/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S7", "continuation/*:12"}, check_s7},
        {{"native.equilibrium.s8", "Equilibrium_S8", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "amount-scale/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S8"}, check_s8},
        {{"native.equilibrium.c1", "Core_C1_amount_scale", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "C1", "c1-amount-scale/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "C1"}, check_c1_amount_scale},
        {{"native.equilibrium.c1-multipliers", "Core_C1_multiplier_scale", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "C1", "c1-multiplier-scale/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "C1"}, check_c1_multiplier_scale},
        {{"native.equilibrium.s11", "Equilibrium_S11", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "standard-affinity,pressure-gap/*,final-gap", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S11", "pressure-gap/*:3"}, check_s11},
        {{"native.equilibrium.s13", "Equilibrium_S13", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "phase-count,composition-*,density-root", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S13", "composition-*:2"}, check_s13},
        {{"native.equilibrium.s15", "Equilibrium_S15", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "vanished-phase,deficit,distance", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S15"}, check_s15},
        {{"native.equilibrium.s16", "Equilibrium_S16", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "continuation/*,affinity/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S16", "continuation/*:18,affinity/*:6"}, check_s16},
    });
}
