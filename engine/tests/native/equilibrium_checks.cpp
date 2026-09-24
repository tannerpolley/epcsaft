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
#include <functional>
#include <iterator>
#include <memory>
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
    CHECK(result.convergence.residual == compiled.classified_residual_norm(result.coordinates));
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

    // #84: the native all-solute path must include the density factor, retain
    // the selected pressure, and reject a charged reaction before any solve.
    epcsaft::HelmholtzModel reference_model(properties_support::neutral_ternary());
    const auto ratios = [](std::vector<double> values) {
        epcsaft::NeutralReferenceDeclaration declaration;
        declaration.solvent_index = 0;
        declaration.standard_molality_mol_per_kg = 1.0;
        declaration.solute_molality_ratios = values;
        return declaration;
    };
    const std::vector<double> nu = {-1.0, 2.0, 0.0};
    const auto reference = epcsaft::ReactionReference{
        std::nullopt, epcsaft::ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    const auto first = epcsaft::reaction_reference(
        reference_model, 300.0, 1.0e5, ratios({0.0, 1.0, 1.0}), nu, reference);
    const auto second = epcsaft::reaction_reference(
        reference_model, 300.0, 1.0e5, ratios({0.0, 2.0, 1.0}), nu, reference);
    CHECK(first.status == epcsaft::NeutralReferenceStatus::Available
          && second.status == epcsaft::NeutralReferenceStatus::Available);
    CHECK(first.terminal_contraction.has_value() && second.terminal_contraction.has_value());
    const double path_gap = std::abs(*first.terminal_contraction - *second.terminal_contraction);
    recording_support::compare_custom(
        "neutral-reference/path-gap", "terminal_contraction_gap", "1", path_gap, 5.0e-5,
        recording_support::fields({recording_support::field("T", 300.0),
            recording_support::field("P", 1.0e5), recording_support::field("nu", nu)}),
        "Relationship", "independent-neutral-paths",
        "independent positive neutral solute paths agree at the finite terminal point",
        "isfinite(model) && model <= reference", path_gap <= 5.0e-5,
        "s2-neutral-reference", "path-gap", "|C_a-C_b|", "T,P,nu", "|C_a-C_b|<=5e-5", "path");
    const double molality = 1.0e-14;
    const double solvent_moles = 1.0 / reference_model.parameters().molar_mass_kg_per_mol[0];
    const double total = solvent_moles + 2.0 * molality;
    const std::vector<double> terminal_x = {
        solvent_moles / total, molality / total, molality / total};
    const auto terminal = reference_model.state_at_pressure(
        300.0, 1.0e5, terminal_x, epcsaft::DensityBranch::Liquid);
    CHECK(terminal.log_fugacity_coefficient.has_value());
    double identity = 0.0;
    for (std::size_t i = 0; i < nu.size(); ++i)
        identity += nu[i] * ((*terminal.log_fugacity_coefficient)[i]
            + std::log(1.0e5 / (epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0)));
    recording_support::compare(
        "neutral-reference/lambda-identity", "terminal_contraction", "1",
        *first.terminal_contraction, identity,
        recording_support::fields({recording_support::field("rho", terminal.molar_density),
            recording_support::field("observed_terminal_change", *first.observed_terminal_change)}),
        "Relationship", "fugacity-pressure-identity",
        "Lambda includes log(P/(RT rho0)) in the independent fugacity form",
        "reference_or_one", 0.0, 1.0e-9,
        close_with_tolerance(*first.terminal_contraction, identity, 1.0e-9),
        "s2-neutral-reference", "lambda-identity", "C", "rho,phi,T,P", "C=nu.Lambda", "identity");
    const auto fixed_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5, ratios({0.0, 1.0, 1.0}), nu,
        epcsaft::ReactionReference{1.0e5, reference.source_basis});
    const auto system_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5, ratios({0.0, 1.0, 1.0}), nu, reference);
    CHECK(fixed_pressure.status == epcsaft::NeutralReferenceStatus::Available
          && system_pressure.status == epcsaft::NeutralReferenceStatus::Available
          && fixed_pressure.reference_pressure_pa == 1.0e5
          && fixed_pressure.terminal_contraction.has_value()
          && system_pressure.terminal_contraction.has_value());
    CHECK(std::abs(*fixed_pressure.terminal_contraction - *first.terminal_contraction) <= 1.0e-12);
    const double pressure_gap = std::abs(
        *system_pressure.terminal_contraction - *fixed_pressure.terminal_contraction);
    recording_support::compare_custom(
        "neutral-reference/pressure-path", "pressure_contraction_gap", "1", pressure_gap, 1.0e-6,
        recording_support::fields({recording_support::field("system_pressure", 2.0e5),
            recording_support::field("reference_pressure", 1.0e5)}),
        "Relationship", "source-pressure-selection",
        "system and explicit source pressures produce distinct Delta-nu contraction values",
        "isfinite(model) && model > reference", pressure_gap > 1.0e-6,
        "s2-neutral-reference", "pressure-path", "|C_P-C_P0|", "T,P,P0", "|gap|>1e-6", "pressure");
    CHECK(system_pressure.terminal_contraction_pressure_derivative.has_value()
          && system_pressure.terminal_contraction_pressure_second_derivative.has_value()
          && fixed_pressure.terminal_contraction_pressure_derivative.has_value()
          && fixed_pressure.terminal_contraction_pressure_second_derivative.has_value()
          && *fixed_pressure.terminal_contraction_pressure_derivative == 0.0
          && *fixed_pressure.terminal_contraction_pressure_second_derivative == 0.0);
    const double pressure_step = 100.0;
    const auto lower_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5 - pressure_step, ratios({0.0, 1.0, 1.0}), nu, reference);
    const auto upper_pressure = epcsaft::reaction_reference(
        reference_model, 300.0, 2.0e5 + pressure_step, ratios({0.0, 1.0, 1.0}), nu, reference);
    CHECK(lower_pressure.terminal_contraction.has_value()
          && upper_pressure.terminal_contraction.has_value());
    const double centered_pressure_derivative =
        (*upper_pressure.terminal_contraction - *lower_pressure.terminal_contraction)
        / (2.0 * pressure_step);
    const double centered_pressure_second_derivative =
        (*upper_pressure.terminal_contraction
            - 2.0 * *system_pressure.terminal_contraction
            + *lower_pressure.terminal_contraction)
        / (pressure_step * pressure_step);
    const auto jet_close = [](double left, double right) {
        return std::abs(left - right) <= 1.0e-3
            * std::max({std::abs(left), std::abs(right), 1.0e-12});
    };
    CHECK(jet_close(*system_pressure.terminal_contraction_pressure_derivative,
                    centered_pressure_derivative)
          && jet_close(*system_pressure.terminal_contraction_pressure_second_derivative,
                       centered_pressure_second_derivative));

    // Identical EOS species provide a finite nonideal counterproof: the density
    // term and residual chemical potential do not separate into an ideal -1/P law.
    auto identical_parameters = properties_support::parameters(2);
    identical_parameters.segment_count = {1.6069, 1.6069};
    identical_parameters.sigma_constant = {3.5206, 3.5206};
    identical_parameters.dispersion_energy_over_k = {191.42, 191.42};
    identical_parameters.molar_mass_kg_per_mol = {0.03007, 0.03007};
    epcsaft::HelmholtzModel identical_model(identical_parameters);
    const auto identical_declaration = ratios({0.0, 1.0});
    const std::vector<double> identical_nu = {-2.0, 1.0};
    const auto identical_reference = epcsaft::reaction_reference(
        identical_model, 300.0, 1.0e5, identical_declaration, identical_nu, reference);
    CHECK(identical_reference.status == epcsaft::NeutralReferenceStatus::Available
          && identical_reference.terminal_contraction_pressure_derivative.has_value()
          && identical_reference.terminal_contraction_pressure_second_derivative.has_value());
    const double identical_molality = 1.0e-14;
    const double identical_solvent_moles = 1.0 / identical_parameters.molar_mass_kg_per_mol[0];
    const double identical_total = identical_solvent_moles + identical_molality;
    const auto identical_terminal = identical_model.state_at_pressure(
        300.0, 1.0e5,
        {identical_solvent_moles / identical_total, identical_molality / identical_total},
        epcsaft::DensityBranch::Liquid);
    const double identical_cp = -1.0 / (
        identical_terminal.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0);
    const double identical_cpp = 1.0 / (
        identical_terminal.molar_density * identical_terminal.molar_density
        * epcsaft::GAS_CONSTANT_J_PER_MOL_K * 300.0
        * identical_terminal.pressure_density_derivative);
    CHECK(jet_close(*identical_reference.terminal_contraction_pressure_derivative, identical_cp)
          && jet_close(*identical_reference.terminal_contraction_pressure_second_derivative,
                       identical_cpp));

    // #84 temperature slope, common native/Python case: aqueous NaCl with the formula
    // reaction H2O = Na+ + Cl- at the system pressure and at a fixed 2e5 Pa source.
    // dC/dT|_P against a centered same-branch terminal difference; floor: truncation
    // step^2 C^(3)/6 ~ (0.01/298)^2 |C'| and value roundoff/step are below 1e-7 relative.
    // The measured multiplier M = T|dD|/|dC| qualifies the observed slope criterion.
    {
        epcsaft::HelmholtzModel ionic(properties_support::aqueous_electrolyte());
        const auto declaration = ratios({0.0, 1.0, 1.0});
        const std::vector<std::vector<double>> rows = {{-1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}};
        const std::vector<epcsaft::ReactionReference> sources = {
            reference, epcsaft::ReactionReference{2.0e5, reference.source_basis}};
        const double temperature = 298.15, pressure = 1.0e5, step = 1.0e-2;
        const auto center = epcsaft::reaction_references(
            ionic, temperature, pressure, declaration, rows, sources, false, std::nullopt, 2);
        const auto lower = epcsaft::reaction_references(
            ionic, temperature - step, pressure, declaration, rows, sources, false, std::nullopt, 1);
        const auto upper = epcsaft::reaction_references(
            ionic, temperature + step, pressure, declaration, rows, sources, false, std::nullopt, 1);
        // #138 calorics: water carries a manufactured ideal-gas record (cp 33.59 J/(mol K),
        // h = -241826 J/mol at 298.15 K); Na+ and Cl- have none and are reaction-completed.
        epcsaft::ThermochemistryRecord water_record;
        water_record.standard_pressure = 1.0e5;
        water_record.component_intervals = {{{200.0, 500.0, true, true, {epcsaft::IdealPolynomial{
            {33.59}, 298.15, -241826.0, 188.83}, 1.0e5}}}, {}, {}};
        epcsaft::HelmholtzModel ionic_thermal(properties_support::aqueous_electrolyte(), {}, water_record);
        for (std::size_t row = 0; row < rows.size(); ++row) {
            CHECK(center[row].terminal_contraction_temperature_derivative.has_value()
                  && center[row].observed_temperature_derivative_change.has_value()
                  && center[row].observed_terminal_change.has_value()
                  && lower[row].terminal_contraction.has_value()
                  && upper[row].terminal_contraction.has_value());
            const double slope = *center[row].terminal_contraction_temperature_derivative;
            const double centered = (*upper[row].terminal_contraction
                - *lower[row].terminal_contraction) / (2.0 * step);
            const std::string label = row == 0 ? "system-pressure" : "fixed-pressure";
            const auto inputs = recording_support::fields({
                recording_support::field("T", temperature), recording_support::field("P", pressure),
                recording_support::field("P_ref", sources[row].reference_pressure_pa.value_or(pressure)),
                recording_support::field("step", step)});
            recording_support::compare(
                ("neutral-reference/temperature-slope/" + label).c_str(), "dC_dT", "1/K",
                slope, centered, inputs, "Relationship", "centered-terminal-difference",
                "dC/dT|_P carries the implicit density and association chains of the reference",
                "reference_scale", 0.0, 1.0e-6,
                std::abs(slope - centered) <= 1.0e-6 * std::abs(centered),
                "s2-neutral-reference", "temperature-slope", "dC/dT", "T,P,nu", "rel<=1e-6", "slope");
            const double change = *center[row].observed_temperature_derivative_change;
            recording_support::compare_custom(
                ("neutral-reference/temperature-limit/" + label).c_str(),
                "T_dC_dT_terminal_change", "1", change, 5.0e-5, inputs,
                "Relationship", "observed-slope-convergence",
                "the slope's terminal ladder change stays within the value criterion in RT units",
                "isfinite(model) && model <= reference", change <= 5.0e-5,
                "s2-neutral-reference", "temperature-limit", "T|dD|", "T,P,nu", "T|dD|<=5e-5", "limit");
            // The same manufactured formula-reaction solve as the Python test; its first-order
            // T and T/P actions of x_Na are the rows that the installed wheel must reproduce.
            auto problem = equilibrium_support::flash(3);
            problem.temperature = epcsaft::equilibrium::FreeValue::fixed(temperature);
            problem.pressure = epcsaft::equilibrium::FreeValue::fixed(pressure);
            problem.feed.values = {1.0, 1.0e-6, 1.0e-6};
            problem.source_reference = declaration;
            problem.reactions = {{"H2O=Na++Cl-", rows[row], std::nullopt, std::nullopt,
                epcsaft::ReactionLogPolynomial{std::log(1.0e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0,
                                               "analytic-source"}, sources[row]}};
            const auto solved = epcsaft::equilibrium::solve_equilibrium(ionic, problem);
            CHECK(solved.success);
            const auto compiled = epcsaft::equilibrium::compile_problem(ionic, problem);
            epcsaft::equilibrium::SolvedStateObservable sodium;
            sodium.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
            sodium.component = 1;
            for (const auto& [name, direction] : std::vector<std::pair<std::string,
                     epcsaft::equilibrium::SolvedStateActionDirection>>{
                     {"T", {1.0, 0.0, {0.0, 0.0, 0.0}, {}}}, {"T+P", {1.0, 1.0e5, {0.0, 0.0, 0.0}, {}}}}) {
                const auto action = epcsaft::equilibrium::solved_state_actions(
                    compiled, solved, {{sodium}, {direction}}).results[0];
                CHECK(action.status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                      && action.action.has_value());
                recording_support::observe(
                    ("reference-temperature/ionic/" + std::string(row == 0 ? "system" : "fixed")
                     + "/dx_Na/" + name).c_str(), "dx_Na", "1/K", action.action.value_or(0.0),
                    std::nullopt, inputs, "Relationship", "solved-state-first-action",
                    "native first-order referenced T/P action for the installed-wheel comparison");
            }
            recording_support::observe(
                ("neutral-reference/temperature-multiplier/" + label).c_str(), "M", "1",
                change / *center[row].observed_terminal_change, std::nullopt, inputs,
                "Relationship", "debye-huckel-sqrt-m-form",
                "measured slope-to-value terminal-change multiplier for this model and state");
            // Order two on the ionic ladder (Debye-Hueckel and Born third derivatives) against the
            // centered difference of the order-one slope; floor: step^2 |d4C/dT4|/6 and slope
            // roundoff/step, below 1e-7 relative. Admission is sized to cp: R T^2 max|dD2|.
            CHECK(center[row].terminal_contraction_temperature_second_derivative.has_value()
                  && lower[row].terminal_contraction_temperature_derivative.has_value()
                  && upper[row].terminal_contraction_temperature_derivative.has_value());
            const double second = center[row].terminal_contraction_temperature_second_derivative.value_or(0.0);
            const double second_centered = (upper[row].terminal_contraction_temperature_derivative.value_or(0.0)
                - lower[row].terminal_contraction_temperature_derivative.value_or(0.0)) / (2.0 * step);
            recording_support::compare(
                ("neutral-reference/temperature-second-slope/" + label).c_str(), "d2C_dT2", "1/K^2",
                second, second_centered, inputs, "Relationship", "centered-slope-difference",
                "d2C/dT2|_P carries the ionic, density and association chains through E3",
                "reference_scale", 0.0, 1.0e-6,
                std::abs(second - second_centered) <= 1.0e-6 * std::abs(second_centered),
                "s2-neutral-reference", "temperature-second-slope", "d2C/dT2", "T,P,nu", "rel<=1e-6", "slope");
            const double second_change = epcsaft::GAS_CONSTANT_J_PER_MOL_K
                * center[row].observed_temperature_second_derivative_change.value_or(1.0);
            recording_support::compare_custom(
                ("neutral-reference/temperature-second-limit/" + label).c_str(),
                "R_T2_d2C_dT2_terminal_change", "J/(mol K)", second_change, 1.0e-2, inputs,
                "Relationship", "observed-slope-convergence",
                "the second slope's terminal ladder change is negligible at the heat-capacity scale",
                "isfinite(model) && model <= reference", second_change <= 1.0e-2,
                "s2-neutral-reference", "temperature-second-limit", "R T^2|dD2|", "T,P,nu",
                "R T^2|dD2|<=0.01 J/(mol K)", "limit");
            const auto thermal_solved = epcsaft::equilibrium::solve_equilibrium(ionic_thermal, problem);
            CHECK(thermal_solved.success);
            epcsaft::equilibrium::SolvedStateObservable enthalpy, capacity;
            enthalpy.kind = capacity.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseProperty;
            enthalpy.property = epcsaft::PropertyObservable::TotalEnthalpy;
            capacity.property = epcsaft::PropertyObservable::TotalIsobaricHeatCapacity;
            const auto calorics = epcsaft::equilibrium::solved_state_actions(
                epcsaft::equilibrium::compile_problem(ionic_thermal, problem), thermal_solved,
                {{enthalpy, capacity}, {{0.0, 0.0, {0.0, 0.0, 0.0}, {}}}}).results;
            for (std::size_t i = 0; i < 2; ++i) {
                CHECK(calorics[i].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
                      && calorics[i].value.has_value());
                recording_support::observe(
                    ("reference-calorics/ionic/" + std::string(row == 0 ? "system" : "fixed")
                     + (i == 0 ? "/h" : "/cp")).c_str(), i == 0 ? "h" : "cp", i == 0 ? "J/mol" : "J/(mol K)",
                    calorics[i].value.value_or(0.0), std::nullopt, inputs, "Relationship", "reference-species-calorics",
                    "native total caloric value with reaction-completed ions for the installed-wheel comparison");
            }
        }
    }
    // #138 order-two slope on the associating MEA-water fixture, solvent water, formula
    // reaction H2O = MEA: d2C/dT2|_P against the centered difference of dC/dT|_P; floor:
    // truncation step^2 |d4C/dT4|/6 and slope roundoff/step, both below 1e-7 relative.
    {
        epcsaft::HelmholtzModel associating(properties_support::associating_binary());
        epcsaft::NeutralReferenceDeclaration water;
        water.solvent_index = 1;
        water.solute_molality_ratios = {1.0, 0.0};
        water.standard_molality_mol_per_kg = 1.0;
        const std::vector<std::vector<double>> rows = {{1.0, -1.0}};
        const std::vector<epcsaft::ReactionReference> sources = {reference};
        const double temperature = 313.15, pressure = 1.0e5, step = 1.0e-2;
        const auto center = epcsaft::reaction_references(
            associating, temperature, pressure, water, rows, sources, false, std::nullopt, 2);
        const auto lower = epcsaft::reaction_references(
            associating, temperature - step, pressure, water, rows, sources, false, std::nullopt, 1);
        const auto upper = epcsaft::reaction_references(
            associating, temperature + step, pressure, water, rows, sources, false, std::nullopt, 1);
        CHECK(center[0].terminal_contraction_temperature_second_derivative.has_value()
              && lower[0].terminal_contraction_temperature_derivative.has_value()
              && upper[0].terminal_contraction_temperature_derivative.has_value());
        const double second = center[0].terminal_contraction_temperature_second_derivative.value_or(0.0);
        const double centered = (upper[0].terminal_contraction_temperature_derivative.value_or(0.0)
            - lower[0].terminal_contraction_temperature_derivative.value_or(0.0)) / (2.0 * step);
        const auto inputs = recording_support::fields({
            recording_support::field("T", temperature), recording_support::field("P", pressure),
            recording_support::field("step", step)});
        recording_support::compare(
            "neutral-reference/temperature-second-slope/associating", "d2C_dT2", "1/K^2",
            second, centered, inputs, "Relationship", "centered-slope-difference",
            "d2C/dT2|_P carries the implicit density and association chains through E3",
            "reference_scale", 0.0, 1.0e-6, std::abs(second - centered) <= 1.0e-6 * std::abs(centered),
            "s2-neutral-reference", "temperature-second-slope", "d2C/dT2", "T,P,nu", "rel<=1e-6", "slope");
        const double change = epcsaft::GAS_CONSTANT_J_PER_MOL_K
            * center[0].observed_temperature_second_derivative_change.value_or(1.0);
        recording_support::compare_custom(
            "neutral-reference/temperature-second-limit/associating", "R_T2_d2C_dT2_terminal_change", "J/(mol K)",
            change, 1.0e-2, inputs, "Relationship", "observed-slope-convergence",
            "the second slope's terminal ladder change is negligible at the heat-capacity scale",
            "isfinite(model) && model <= reference", change <= 1.0e-2,
            "s2-neutral-reference", "temperature-second-limit", "R T^2|dD2|", "T,P,nu",
            "R T^2|dD2|<=0.01 J/(mol K)", "limit");
    }


    const auto charged = epcsaft::reaction_reference(
        epcsaft::HelmholtzModel(properties_support::aqueous_electrolyte()), 298.15, 1.0e5,
        ratios({0.0, 1.0, 1.0}), {-1.0, 1.0, 0.0}, reference);
    CHECK(charged.status == epcsaft::NeutralReferenceStatus::NonNeutralDeclaration);
    auto referenced = equilibrium_support::flash(3);
    referenced.feed.values = {2.0, 0.2, 0.2};
    referenced.source_reference = ratios({0.0, 1.0, 1.0});
    const auto law = epcsaft::ReactionLogPolynomial{
        std::log(0.3), 0.0, 0.0, 0.0, 300.0, 250.0, 350.0, "native-reference-source"};
    referenced.reactions = {{"A=2B", {-1.0, 2.0, 0.0}, std::nullopt, std::nullopt,
                             law, reference}};
    auto free_pressure_reference = referenced;
    free_pressure_reference.phases = {
        equilibrium_support::phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        equilibrium_support::phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
            epcsaft::equilibrium::AmountSpec::pinned(0.0)),
    };
    free_pressure_reference.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    const auto free_pressure_compiled = epcsaft::equilibrium::compile_problem(
        reference_model, free_pressure_reference);
    const auto free_pressure_point = free_pressure_compiled.generic_point();
    const auto free_pressure_prepared = free_pressure_compiled.prepare(free_pressure_point);
    CHECK(std::all_of(free_pressure_prepared.standard.pressure_derivative.begin(),
                      free_pressure_prepared.standard.pressure_derivative.end(),
                      [](double value) { return std::isfinite(value); }));
    const std::size_t control_width = 2 + reference_model.parameters().component_count();
    std::vector<double> pressure_control(control_width, 0.0);
    pressure_control[1] = 1.0;
    std::vector<std::vector<double>> pressure_state(3,
        std::vector<double>(free_pressure_prepared.point.size(), 0.0));
    pressure_state[0] = free_pressure_prepared.point;
    pressure_state[1][free_pressure_compiled.primal_size - 1] = 1.0;
    const auto pressure_series = free_pressure_compiled.row_series(
        free_pressure_prepared,
        pressure_state,
        pressure_control, 2);
    CHECK(pressure_series.energy == epcsaft::PropertyActionStatus::Available);
    const auto chemical = std::find_if(
        free_pressure_compiled.rows.begin(), free_pressure_compiled.rows.end(),
        [](const auto& row) { return row.kind == epcsaft::equilibrium::RowKind::Chemical; });
    CHECK(chemical != free_pressure_compiled.rows.end());
    if (chemical != free_pressure_compiled.rows.end()) {
        const std::size_t row = static_cast<std::size_t>(
            std::distance(free_pressure_compiled.rows.begin(), chemical));
        const double pressure_step = 100.0;
        auto plus = free_pressure_prepared.point;
        auto minus = free_pressure_prepared.point;
        plus[free_pressure_compiled.primal_size - 1] += pressure_step;
        minus[free_pressure_compiled.primal_size - 1] -= pressure_step;
        const auto plus_rows = free_pressure_compiled.evaluate_at(
            free_pressure_compiled.prepare(plus)).values;
        const auto minus_rows = free_pressure_compiled.evaluate_at(
            free_pressure_compiled.prepare(minus)).values;
        const auto central_rows = free_pressure_compiled.evaluate_at(free_pressure_prepared).values;
        const double centered_first = (plus_rows[row] - minus_rows[row])
            / (2.0 * pressure_step);
        const double centered_second_coefficient = (plus_rows[row]
            - 2.0 * central_rows[row] + minus_rows[row])
            / (2.0 * pressure_step * pressure_step);
        const auto close_pressure_jet = [](double actual, double expected) {
            return std::isfinite(actual) && std::isfinite(expected)
                && std::abs(actual - expected) <= 2.0e-3
                    * std::max({std::abs(actual), std::abs(expected), 1.0e-12});
        };
        CHECK(close_pressure_jet(pressure_series.rows[1][row], centered_first)
              && close_pressure_jet(pressure_series.rows[2][row], centered_second_coefficient));
    }
    bool pressure_third_rejected = false;
    try {
        (void)free_pressure_compiled.row_series(
            free_pressure_prepared,
            {free_pressure_prepared.point,
             std::vector<double>(free_pressure_prepared.point.size(), 0.0),
             std::vector<double>(free_pressure_prepared.point.size(), 0.0),
             std::vector<double>(free_pressure_prepared.point.size(), 0.0)},
            pressure_control, 3);
    } catch (const std::invalid_argument&) {
        pressure_third_rejected = true;
    }
    CHECK(pressure_third_rejected);
    const auto fixed_pressure_compiled = epcsaft::equilibrium::compile_problem(
        reference_model, referenced);
    const auto fixed_pressure_prepared = fixed_pressure_compiled.prepare(
        fixed_pressure_compiled.generic_point());
    bool fixed_pressure_third_rejected = false;
    try {
        (void)fixed_pressure_compiled.row_series(
            fixed_pressure_prepared,
            {fixed_pressure_prepared.point,
             std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
             std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
             std::vector<double>(fixed_pressure_prepared.point.size(), 0.0)},
            pressure_control, 3);
    } catch (const std::invalid_argument&) {
        fixed_pressure_third_rejected = true;
    }
    CHECK(fixed_pressure_third_rejected);
    std::fill(pressure_control.begin(), pressure_control.end(), 0.0);
    pressure_control[2] = 1.0;
    const auto fixed_feed_series = fixed_pressure_compiled.row_series(
        fixed_pressure_prepared,
        {fixed_pressure_prepared.point,
         std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
         std::vector<double>(fixed_pressure_prepared.point.size(), 0.0),
         std::vector<double>(fixed_pressure_prepared.point.size(), 0.0)},
        pressure_control, 3);
    CHECK(fixed_feed_series.energy == epcsaft::PropertyActionStatus::Available);
    epcsaft::ThermochemistryRecord thermal_record;
    thermal_record.standard_pressure = 1.0e5;
    for (std::size_t i = 0; i < 3; ++i)
        thermal_record.component_intervals.push_back(
            properties_support::p3_record(properties_support::p3_forms()[0].second)
                .component_intervals.front());
    epcsaft::HelmholtzModel thermal_reference_model(
        properties_support::neutral_ternary(), {}, thermal_record);
    const auto referenced_result = epcsaft::equilibrium::solve_equilibrium(
        thermal_reference_model, referenced);
    CHECK(referenced_result.success
          && referenced_result.reference_diagnostics.size() == 1
          && referenced_result.reference_diagnostics[0].terminal_contraction.has_value());
    const auto compiled_reference = epcsaft::equilibrium::compile_problem(
        thermal_reference_model, referenced);
    epcsaft::equilibrium::SolvedStateActionRequest action_request;
    epcsaft::equilibrium::SolvedStateObservable composition_observable;
    composition_observable.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseMoleFraction;
    composition_observable.component = 1;
    action_request.observables = {composition_observable};
    action_request.directions = {{1.0, 0.0, {0.0, 0.0, 0.0}, {}}};
    const auto temperature_action = epcsaft::equilibrium::solved_state_actions(
        compiled_reference, referenced_result, action_request);
    CHECK(temperature_action.results[0].status
              == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && temperature_action.results[0].action.has_value());
    composition_observable.kind = epcsaft::equilibrium::SolvedStateObservableKind::PhaseProperty;
    composition_observable.property = epcsaft::PropertyObservable::IdealEnthalpy;
    const auto total_observable = [&] {
        auto result = composition_observable;
        result.property = epcsaft::PropertyObservable::TotalEnthalpy;
        return result;
    }();
    const auto residual_observable = [&] {
        auto result = composition_observable;
        result.property = epcsaft::PropertyObservable::ResidualEnthalpy;
        return result;
    }();
    action_request.observables = {composition_observable, total_observable, residual_observable};
    action_request.directions = {{0.0, 0.0, {0.0, 0.0, 0.0}, {}}};
    const auto caloric = epcsaft::equilibrium::solved_state_actions(
        compiled_reference, referenced_result, action_request);
    const std::vector<double> central_x(
        referenced_result.mole_fractions.begin(), referenced_result.mole_fractions.begin() + 3);
    const auto central_state = thermal_reference_model.state(
        referenced_result.temperature.front(), referenced_result.molar_densities[0], central_x);
    CHECK(central_state.ideal_enthalpy.has_value() && central_state.total_enthalpy.has_value()
          && std::isfinite(central_state.residual_enthalpy));
    CHECK(caloric.results.size() == 3
          && caloric.results[0].status
                 == epcsaft::equilibrium::SolvedStateActionStatus::ReferenceUnavailable
          && !caloric.results[0].value.has_value() && !caloric.results[0].action.has_value()
          && caloric.results[1].status
                 == epcsaft::equilibrium::SolvedStateActionStatus::ReferenceUnavailable
          && !caloric.results[1].value.has_value() && !caloric.results[1].action.has_value()
          && caloric.results[2].status
                 == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && caloric.results[2].value.has_value()
          && close_with_tolerance(*caloric.results[2].value, central_state.residual_enthalpy, 1.0e-9));
    const std::vector<epcsaft::equilibrium::SolvedStateObservable> source_observables = {
        temperature_action.results[0].observable, composition_observable};
    const std::vector<epcsaft::equilibrium::SolvedStateActionDirection> source_directions = {
        {1.0, 0.0, {0.0, 0.0, 0.0}, {}}, {0.0, 0.0, {0.01, 0.0, 0.0}, {}}};
    const auto source_batches = epcsaft::equilibrium::solved_state_first_actions(
        compiled_reference, referenced_result, source_observables, source_directions);
    bool source_match = source_batches.size() == source_directions.size();
    double source_worst = 0.0;
    std::size_t source_factorizations = 0;
    for (std::size_t j = 0; j < source_batches.size() && source_match; ++j) {
        const auto separate = epcsaft::equilibrium::solved_state_actions(
            compiled_reference, referenced_result, {source_observables, {source_directions[j]}});
        source_factorizations += source_batches[j].diagnostics.state_factorizations;
        source_match = source_batches[j].results.size() == separate.results.size();
        for (std::size_t i = 0; i < separate.results.size() && source_match; ++i) {
            const auto& actual = source_batches[j].results[i];
            const auto& expected = separate.results[i];
            source_match = actual.order == expected.order && actual.unit == expected.unit
                && actual.status == expected.status && actual.diagnostic_message == expected.diagnostic_message
                && actual.value.has_value() == expected.value.has_value()
                && actual.action.has_value() == expected.action.has_value();
            if (actual.value && expected.value)
                source_worst = std::max(source_worst, std::abs(*actual.value - *expected.value) / std::max(1.0, std::abs(*expected.value)));
            if (actual.action && expected.action)
                source_worst = std::max(source_worst, std::abs(*actual.action - *expected.action) / std::max(1.0, std::abs(*expected.action)));
        }
    }
    CHECK(source_match && source_worst <= 2.0e-12 && source_factorizations == 1
          && source_batches[0].results[0].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && source_batches[0].results[0].value.has_value() && source_batches[0].results[0].action.has_value()
          && source_batches[1].results[0].status == epcsaft::equilibrium::SolvedStateActionStatus::Available
          && source_batches[1].results[0].value.has_value() && source_batches[1].results[0].action.has_value()
          && source_batches[1].results[1].status == epcsaft::equilibrium::SolvedStateActionStatus::ReferenceUnavailable
          && source_batches[1].results[1].diagnostic_message.has_value());
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
    // Since the #87 correction the returned constraint multipliers are formed at the returned point
    // from the row multipliers, so these rows verify that reconstruction's KKT residual and feed
    // scaling, no longer the optimizer's own duals.
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
    // The vanished phase's KKT-scaled stationarity rows vanish with its amount, so the
    // classified residual meets the association floor while the topology predicate fails:
    // unavailable by topology, with the tiny-step solver status retained.
    const bool tiny_step_stalled = failure.solver_status == "STOP_AT_TINY_STEP"
        && (failure.convergence.stop == epcsaft::core::Stop::stalled
            || failure.convergence.stop == epcsaft::core::Stop::converged_to_floor);

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

// The classifying residual: KKT-scaled stationarity rows for free-amount cells, unscaled
// rows for pinned cells (owner decision 2026-09-16); the unscaled rows stay diagnostics.
double largest_row(const epcsaft::equilibrium::EquilibriumResult& result) {
    return result.convergence.residual;
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

// ------------------------------------------------------------------ E5 solved-state actions

namespace e5 {

using epcsaft::PropertyObservable;
using epcsaft::equilibrium::AmountSpec;
using epcsaft::equilibrium::CompiledSystem;
using epcsaft::equilibrium::EquilibriumProblem;
using epcsaft::equilibrium::EquilibriumResult;
using epcsaft::equilibrium::FreeValue;
using epcsaft::equilibrium::PhaseKind;
using epcsaft::equilibrium::SolvedStateActionBatch;
using epcsaft::equilibrium::SolvedStateActionDirection;
using epcsaft::equilibrium::SolvedStateActionRequest;
using epcsaft::equilibrium::SolvedStateActionStatus;
using epcsaft::equilibrium::SolvedStateObservable;
using Kind = epcsaft::equilibrium::SolvedStateObservableKind;

constexpr double LADDER[3][3] = {{1.0e-2, 3.0e-3, 1.0e-3}, {3.0e-2, 1.0e-2, 3.0e-3}, {1.0e-1, 3.0e-2, 1.0e-2}};
constexpr double CRITERION[3] = {5.0e-6, 5.0e-5, 5.0e-4};
constexpr double STRUCTURAL_ZERO = 1.0e-12;

SolvedStateObservable obs(Kind kind, std::size_t phase = 0, std::size_t component = 0,
                          PropertyObservable property = PropertyObservable::IdealEnthalpy) {
    SolvedStateObservable result;
    result.kind = kind; result.phase = phase; result.component = component; result.property = property;
    return result;
}

SolvedStateActionDirection direction(double T, double P, std::vector<double> feed, std::vector<double> active = {}) {
    SolvedStateActionDirection result;
    result.temperature_k = T; result.pressure_pa = P;
    result.feed_amounts_mol = std::move(feed); result.active_parameters = std::move(active);
    return result;
}

// One qualification fixture: a model recipe, a problem, three directions and observables.
struct Fixture final {
    std::string label;
    epcsaft::Parameters parameters;
    std::vector<epcsaft::ActiveParameter> active;
    std::optional<epcsaft::ThermochemistryRecord> record;
    EquilibriumProblem problem;
    std::vector<SolvedStateActionDirection> directions;
    std::vector<SolvedStateObservable> observables;
    bool scaled = true;
    // Corner branch admission for a pure declared phase (case 6) and an analytic corner
    // reference for selected observables (the A->2B pure-state relation).
    std::optional<epcsaft::DensityBranch> corner_branch;
    std::function<std::optional<double>(const SolvedStateObservable&, double)> analytic;

    // A fresh model with every active parameter shifted (reciprocal pairs move together).
    epcsaft::HelmholtzModel model(const std::vector<double>& shift = {}) const {
        auto p = parameters;
        for (std::size_t a = 0; a < active.size() && a < shift.size(); ++a) {
            std::vector<double>& family = active[a].family == epcsaft::ParameterFamily::SegmentCount ? p.segment_count
                : active[a].family == epcsaft::ParameterFamily::SigmaConstant ? p.sigma_constant
                : active[a].family == epcsaft::ParameterFamily::DispersionEnergyOverK ? p.dispersion_energy_over_k
                : active[a].family == epcsaft::ParameterFamily::Kij ? p.k_ij
                : active[a].family == epcsaft::ParameterFamily::Lij ? p.l_ij
                : active[a].family == epcsaft::ParameterFamily::AssociationEnergyOverK ? p.association_energy_over_k
                : p.association_volume;
            for (const std::size_t index : active[a].indices) family[index] += shift[a];
        }
        return epcsaft::HelmholtzModel(p, active, record);
    }
};

// Every geometry, phase-pressure and residual-potential identity on the supported components, plus T and P.
std::vector<SolvedStateObservable> geometry_observables(const EquilibriumProblem& problem, std::size_t components) {
    std::vector<SolvedStateObservable> result = {obs(Kind::Temperature), obs(Kind::Pressure)};
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        for (const Kind kind : {Kind::PhaseAmount, Kind::PhaseVolume, Kind::PhaseMolarDensity, Kind::PhasePackingFraction, Kind::PhasePressure})
            result.push_back(obs(kind, phase));
        for (std::size_t i = 0; i < components; ++i) {
            result.push_back(obs(Kind::PhaseComponentAmount, phase, i));
            result.push_back(obs(Kind::PhaseMoleFraction, phase, i));
            const auto& support = problem.phases[phase].support;
            if (support.empty() || std::find(support.begin(), support.end(), i) != support.end()) {
                result.push_back(obs(Kind::PhaseResidualChemicalPotentialOverRt, phase, i));
                result.push_back(obs(Kind::PhaseLogActivity, phase, i));
            }
        }
    }
    return result;
}

// The physical observable map of an accepted result through the public owners; never E5.
std::optional<double> mapped(const SolvedStateObservable& item, const epcsaft::HelmholtzModel& model,
                             const EquilibriumResult& result) {
    const std::size_t components = model.parameters().component_count();
    const auto x = [&](std::size_t phase) {
        return std::vector<double>(result.mole_fractions.begin() + static_cast<std::ptrdiff_t>(phase * components),
                                   result.mole_fractions.begin() + static_cast<std::ptrdiff_t>((phase + 1) * components));
    };
    const auto state = [&](std::size_t phase) {
        return model.state(result.temperature.front(), result.molar_densities[phase], x(phase));
    };
    switch (item.kind) {
    case Kind::Temperature: return result.temperature.front();
    case Kind::Pressure: return result.pressure.front();
    case Kind::PhaseAmount: return result.amounts[item.phase];
    case Kind::PhaseComponentAmount: return result.amounts[item.phase] * x(item.phase)[item.component];
    case Kind::PhaseMoleFraction: return x(item.phase)[item.component];
    case Kind::PhaseVolume: return result.amounts[item.phase] / result.molar_densities[item.phase];
    case Kind::PhaseMolarDensity: return result.molar_densities[item.phase];
    case Kind::PhasePackingFraction: return result.packing_fractions[item.phase];
    case Kind::PhasePressure: return state(item.phase).pressure;
    case Kind::PhaseResidualChemicalPotentialOverRt: return state(item.phase).residual_chemical_potential_over_rt[item.component];
    case Kind::PhaseLogActivity:  // ln(x rho / (1 mol/m^3)) + mu_res/RT through the public State
        return std::log(x(item.phase)[item.component] * result.molar_densities[item.phase])
            + state(item.phase).residual_chemical_potential_over_rt[item.component];
    case Kind::PhaseProperty: {
        const auto s = state(item.phase);
        switch (item.property) {
        case PropertyObservable::IdealEnthalpy: return s.ideal_enthalpy;
        case PropertyObservable::IdealInternalEnergy: return s.ideal_internal_energy;
        case PropertyObservable::IdealEntropy: return s.ideal_entropy;
        case PropertyObservable::IdealIsobaricHeatCapacity: return s.ideal_isobaric_heat_capacity;
        case PropertyObservable::IdealIsochoricHeatCapacity: return s.ideal_isochoric_heat_capacity;
        case PropertyObservable::ResidualEnthalpy: return s.residual_enthalpy;
        case PropertyObservable::ResidualInternalEnergy: return s.residual_internal_energy;
        case PropertyObservable::ResidualEntropy: return s.residual_entropy;
        case PropertyObservable::ResidualIsobaricHeatCapacity: return s.residual_isobaric_heat_capacity;
        case PropertyObservable::ResidualIsochoricHeatCapacity: return s.residual_isochoric_heat_capacity;
        case PropertyObservable::TotalEnthalpy: return s.total_enthalpy;
        case PropertyObservable::TotalInternalEnergy: return s.total_internal_energy;
        case PropertyObservable::TotalEntropy: return s.total_entropy;
        case PropertyObservable::TotalIsobaricHeatCapacity: return s.total_isobaric_heat_capacity;
        case PropertyObservable::TotalIsochoricHeatCapacity: return s.total_isochoric_heat_capacity;
        }
        return std::nullopt;
    }
    }
    return std::nullopt;
}

// The sole test-only stencil scaling: max(|q_i|, one physical unit) per external coordinate.
double magnitude(const SolvedStateActionDirection& d, const CompiledSystem& compiled, const std::vector<double>& active_values) {
    double result = std::max(std::abs(d.temperature_k) / std::max(compiled.problem.temperature.value, 1.0),
                             std::abs(d.pressure_pa) / std::max(compiled.problem.pressure.value, 1.0));
    for (std::size_t i = 0; i < d.feed_amounts_mol.size(); ++i)
        result = std::max(result, std::abs(d.feed_amounts_mol[i]) / std::max(std::abs(compiled.n_feed[i]), 1.0));
    for (std::size_t a = 0; a < d.active_parameters.size(); ++a)
        result = std::max(result, std::abs(d.active_parameters[a]) / std::max(std::abs(active_values[a]), 1.0));
    return result;
}

std::vector<double> active_values(const Fixture& fixture) {
    std::vector<double> result;
    for (const auto& item : fixture.active) {
        const auto& p = fixture.parameters;
        const std::vector<double>& family = item.family == epcsaft::ParameterFamily::SegmentCount ? p.segment_count
            : item.family == epcsaft::ParameterFamily::SigmaConstant ? p.sigma_constant
            : item.family == epcsaft::ParameterFamily::DispersionEnergyOverK ? p.dispersion_energy_over_k
            : item.family == epcsaft::ParameterFamily::Kij ? p.k_ij : item.family == epcsaft::ParameterFamily::Lij ? p.l_ij
            : item.family == epcsaft::ParameterFamily::AssociationEnergyOverK ? p.association_energy_over_k
            : p.association_volume;
        result.push_back(family[item.indices.front()]);
    }
    return result;
}

// Independent reference: a fresh nonlinear solve at every stencil corner, mapped physically.
struct Reference final {
    std::vector<std::optional<double>> value;   // per observable, absent when a corner was rejected
    std::vector<std::string> rejected;          // corner failure names
};

Reference reference(const Fixture& fixture, const CompiledSystem& compiled, const EquilibriumResult& central,
                    std::size_t order, double multiplier) {
    const std::size_t components = compiled.model->parameters().component_count();
    const auto values = active_values(fixture);
    std::vector<double> scale(order, 1.0);
    for (std::size_t j = 0; j < order; ++j)
        if (fixture.scaled) scale[j] = magnitude(fixture.directions[j], compiled, values);
    Reference result;
    result.value.assign(fixture.observables.size(), 0.0);
    double product = 1.0;
    for (const double m : scale) product *= m;
    if (product == 0.0) return result;
    const double limit = std::min(fixture.problem.row_tolerance, epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    for (std::size_t corner = 0; corner < (std::size_t(1) << order); ++corner) {
        auto problem = fixture.problem;
        std::vector<double> shift(fixture.active.size(), 0.0);
        double weight = 1.0;
        for (std::size_t j = 0; j < order; ++j) {
            const double sign = corner & (std::size_t(1) << j) ? -1.0 : 1.0;
            weight *= sign;
            const double h = sign * multiplier / scale[j];
            const auto& d = fixture.directions[j];
            if (!problem.temperature.free) problem.temperature.value += h * d.temperature_k;
            if (!problem.pressure.free) problem.pressure.value += h * d.pressure_pa;
            for (std::size_t i = 0; i < components; ++i) problem.feed.values[i] += h * d.feed_amounts_mol[i];
            for (std::size_t a = 0; a < shift.size(); ++a) shift[a] += h * d.active_parameters[a];
        }
        const auto model = fixture.model(shift);
        const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem, central);
        std::string failure;
        if (!solved.success) failure = "corner solve: " + solved.message;
        else if (solved.phase_count != central.phase_count) failure = "corner phase count";
        else if (!solved.requested_tolerance_met.value_or(false)) failure = "corner row tolerance";
        else if (!(solved.convergence.residual <= limit)) failure = "corner classified rows";
        else if (fixture.corner_branch.has_value()) {
            const auto root = epcsaft::solve_density(model, problem.temperature.value, problem.pressure.value, {1.0}, *fixture.corner_branch, std::nullopt);
            if (!(std::abs(std::log(solved.molar_densities[0] / root.molar_density)) <= 1.0e-8)) failure = "corner branch";
        }
        if (!failure.empty()) {
            result.rejected.push_back(failure);
            for (auto& value : result.value) value.reset();
            return result;
        }
        for (std::size_t i = 0; i < fixture.observables.size(); ++i) {
            const auto analytic = fixture.analytic ? fixture.analytic(fixture.observables[i], problem.temperature.value) : std::nullopt;
            const auto value = analytic.has_value() ? analytic : mapped(fixture.observables[i], model, solved);
            if (value.has_value() && result.value[i].has_value()) *result.value[i] += weight * *value;
            else result.value[i].reset();
        }
    }
    for (auto& value : result.value)
        if (value.has_value()) *value *= product / std::pow(2.0 * multiplier, static_cast<double>(order));
    return result;
}

std::string status_field(const SolvedStateActionBatch& batch, std::size_t index) {
    return epcsaft::equilibrium::solved_state_status_name(batch.results[index].status);
}

// Qualify one fixture through the requested orders; every observable/order/step is retained.
bool qualify(const Fixture& fixture, const char* group, std::size_t orders = 3, std::size_t first_order = 1) {
    const auto model = fixture.model();
    const auto compiled = epcsaft::equilibrium::compile_problem(model, fixture.problem);
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    bool ok = central.success;
    CHECK(central.success);
    if (!central.success) return false;
    for (std::size_t order = first_order; order <= orders; ++order) {
        SolvedStateActionRequest request;
        request.observables = fixture.observables;
        request.directions.assign(fixture.directions.begin(), fixture.directions.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = epcsaft::equilibrium::solved_state_actions(compiled, central, request);
        std::vector<Reference> ladder;
        for (const double multiplier : LADDER[order - 1]) ladder.push_back(reference(fixture, compiled, central, order, multiplier));
        for (std::size_t i = 0; i < fixture.observables.size(); ++i) {
            const auto& item = batch.results[i];
            const auto id = std::string(group) + "/" + fixture.label + "/" + epcsaft::equilibrium::solved_state_observable_id(item.observable)
                + "/A" + std::to_string(order);
            const auto inputs = recording_support::fields({
                recording_support::field("case", fixture.label), recording_support::field("order", order),
                recording_support::field("status", status_field(batch, i)), recording_support::field("unit", item.unit),
                recording_support::field("value", epcsaft::core::optional_number(item.value)),
                recording_support::field("scaled", fixture.scaled ? 1U : 0U),
                recording_support::field("ladder", std::vector<double>(LADDER[order - 1], LADDER[order - 1] + 3)),
                recording_support::field("rank", batch.diagnostics.rank.has_value() ? std::to_string(*batch.diagnostics.rank) : "absent"),
                recording_support::field("condition", epcsaft::core::optional_number(batch.diagnostics.condition_estimate)),
                recording_support::field("linear_residual", epcsaft::core::optional_number(batch.diagnostics.maximum_linear_residual)),
                recording_support::field("central_residual", epcsaft::core::optional_number(batch.diagnostics.central_residual)),
                recording_support::field("unscaled_residual", central.residual_norm()),
                recording_support::field("energy_order", batch.diagnostics.maximum_energy_order),
                recording_support::field("rejected_corners", ladder.back().rejected.empty() ? "none" : ladder.back().rejected.front())});
            const double action = item.action.value_or(std::numeric_limits<double>::quiet_NaN());
            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (std::size_t step = 0; step + 1 < ladder.size(); ++step) {
                recording_support::observe((id + "/step" + std::to_string(step)).c_str(), "action", item.unit.c_str(), action,
                    ladder[step].value[i], inputs, "Relationship", "centered-solve-ladder",
                    "coarser reference step retained for the rate record", group, id.c_str());
            }
            // Richardson reference with rate certification (owner decision 2026-09-16): the
            // centered references scale as h^2, so the last two steps extrapolate to the
            // reference; the three steps certify the rate, or already agree within 2c.
            const auto& coarse = ladder[0].value[i];
            const auto& mid = ladder[1].value[i];
            const auto& fine = ladder[2].value[i];
            const double h1 = LADDER[order - 1][0], h2 = LADDER[order - 1][1], h3 = LADDER[order - 1][2];
            const bool structural = fine.has_value() && mid.has_value() && *fine == 0.0 && *mid == 0.0;
            const bool complete = fine.has_value() && mid.has_value() && coarse.has_value();
            // A finest step whose difference grows again is rounding-bounded (the corner
            // solves agree to ~1e-15 relative); its ladder extrapolates from the two coarser
            // steps and certifies nothing beyond their agreement (Execution review, 2026-09-16).
            const bool noise_bounded = complete && std::abs(*mid - *fine) > std::abs(*coarse - *mid);
            std::optional<double> extrapolated;
            if (complete && noise_bounded) extrapolated = *mid + (*mid - *coarse) / (h1 * h1 / (h2 * h2) - 1.0);
            else if (fine.has_value() && mid.has_value()) extrapolated = *fine + (*fine - *mid) / (h2 * h2 / (h3 * h3) - 1.0);
            const double criterion = structural ? STRUCTURAL_ZERO : CRITERION[order - 1] * std::max(1.0, std::abs(extrapolated.value_or(0.0)));
            const double expected_rate = (h1 * h1 - h2 * h2) / (h2 * h2 - h3 * h3);
            bool refined = complete && (noise_bounded ? std::abs(*coarse - *mid) : std::abs(*fine - *mid)) <= 2.0 * criterion;
            double rate = nan;
            if (complete && !refined && !noise_bounded) {
                rate = (*coarse - *mid) / (*mid - *fine);
                refined = rate >= 0.5 * expected_rate && rate <= 2.0 * expected_rate;
            }
            const bool verdict = item.status == SolvedStateActionStatus::Available && item.action.has_value() && extrapolated.has_value()
                && std::abs(action - *extrapolated) <= criterion;
            // The best-step relative agreement is the falsifying strength independent of the
            // absolute floor: a row with |R| <= c would pass a zero action on the criterion alone.
            double best_step = nan;
            for (const auto* step : {&coarse, &mid, &fine})
                if (step->has_value() && std::isfinite(action))
                    best_step = std::min(std::isfinite(best_step) ? best_step : std::numeric_limits<double>::infinity(),
                                         std::abs(action - **step) / std::max(std::abs(**step), std::numeric_limits<double>::min()));
            const auto ladder_inputs = recording_support::fields({inputs, recording_support::field("reference_steps",
                std::vector<double>{coarse.value_or(nan), mid.value_or(nan), fine.value_or(nan)}),
                recording_support::field("rate", rate), recording_support::field("expected_rate", expected_rate),
                recording_support::field("noise_bounded", noise_bounded ? 1U : 0U),
                recording_support::field("absolute_floor", !structural && std::abs(extrapolated.value_or(0.0)) <= CRITERION[order - 1] ? 1U : 0U),
                recording_support::field("best_step_relative_agreement", best_step)});
            recording_support::compare(id.c_str(), "action", item.unit.c_str(), action, extrapolated.value_or(nan),
                ladder_inputs, "Relationship", "centered-solve-ladder", "solved-state action agrees with the h^2-extrapolated fresh-solve reference",
                structural ? "strict_absolute" : "reference_or_one", structural ? STRUCTURAL_ZERO : 0.0, structural ? 0.0 : CRITERION[order - 1],
                verdict, group, id.c_str(), "D^k g", "case,observable,order", "|D^k g-R|\\leq c", "observable,order");
            // A ladder certified by step agreement records that agreement, not an uncomputed rate.
            const bool by_rate = std::isfinite(rate);
            const double spread = complete ? std::abs(noise_bounded ? *coarse - *mid : *fine - *mid) : nan;
            recording_support::compare_custom((id + "/rate").c_str(), by_rate ? "difference_ratio" : "step_agreement",
                by_rate ? "1" : item.unit.c_str(),
                by_rate ? rate : spread, by_rate ? expected_rate : 2.0 * criterion, ladder_inputs,
                "Relationship", "centered-solve-ladder",
                "the reference ladder certifies second-order convergence, or its last two informative steps agree within twice the criterion",
                "|fine-mid| <= 2c || expected/2 <= rate <= 2 expected || (noise_bounded && |coarse-mid| <= 2c)", refined,
                group, id.c_str(), by_rate ? "(\\Delta_1-\\Delta_2)/(\\Delta_2-\\Delta_3)" : "|R_i-R_{i+1}|", "case,observable,order",
                by_rate ? "\\approx(h_1^2-h_2^2)/(h_2^2-h_3^2)" : "\\leq 2c", "observable,order");
            ok = ok && verdict && refined;
        }
    }
    return ok;
}

bool record_bool(const char* id, const char* group, const char* claim, bool verdict, const std::string& inputs) {
    recording_support::compare(id, "verdict", "bool", verdict ? 1.0 : 0.0, 1.0, inputs, "Relationship", "solved-state-action",
        claim, "absolute", 0.0, 0.0, verdict, group, id, "verdict", "inputs", claim, "case");
    return verdict;
}

SolvedStateActionBatch actions(const Fixture& fixture, const epcsaft::HelmholtzModel& model, const EquilibriumResult& central,
                               std::vector<SolvedStateObservable> observables, std::size_t order) {
    SolvedStateActionRequest request;
    request.observables = std::move(observables);
    request.directions.assign(fixture.directions.begin(), fixture.directions.begin() + static_cast<std::ptrdiff_t>(order));
    return epcsaft::equilibrium::solved_state_actions(epcsaft::equilibrium::compile_problem(model, fixture.problem), central, request);
}

Fixture neutral_case() {
    Fixture f;
    f.label = "neutral-ternary";
    f.parameters = properties_support::neutral_ternary();
    f.active = {{epcsaft::ParameterFamily::Kij, {1, 3}}, {epcsaft::ParameterFamily::SegmentCount, {0}}};
    f.problem = equilibrium_support::flash(3);
    f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor), equilibrium_support::phase("L", PhaseKind::Liquid)};
    f.problem.temperature.value = 180.0;
    f.directions = {direction(0.2, 2000.0, {-0.01, 0.006, 0.004}, {0.001, 0.01}),
                    direction(-0.1, -1000.0, {0.003, -0.004, 0.001}, {-0.0005, 0.004}),
                    direction(0.05, 500.0, {0.002, 0.001, -0.003}, {0.0002, -0.002})};
    f.observables = geometry_observables(f.problem, 3);
    return f;
}

Fixture pinned_case(double fraction) {
    Fixture f;
    f.label = fraction > 0.0 ? "positive-pin" : "zero-pin";
    f.parameters = properties_support::neutral_binary();
    f.problem = equilibrium_support::boundary();
    f.problem.phases[1].amount = AmountSpec::pinned(fraction);
    f.problem.temperature = FreeValue::fixed(180.0);
    f.problem.feed.values = {0.5, 0.5};
    f.directions = {direction(0.2, 0.0, {0.01, 0.02}), direction(-0.1, 0.0, {0.005, -0.003}), direction(0.05, 0.0, {-0.002, 0.006})};
    f.observables = geometry_observables(f.problem, 2);
    f.observables.push_back(obs(Kind::PhaseProperty, 1, 0, PropertyObservable::ResidualEnthalpy));
    f.observables.push_back(obs(Kind::PhaseProperty, 1, 0, PropertyObservable::ResidualIsobaricHeatCapacity));
    return f;
}

}  // namespace e5

void check_e5_neutral() {
    bool all_ok = true;
    using namespace e5;
    const auto fixture = neutral_case();
    all_ok = (qualify(fixture, "e5-1")) && all_ok;
    // Permutations: values and actions follow physical identity, not vector position.
    const auto model = fixture.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, fixture.problem);
    auto independent_directions = fixture.directions;
    independent_directions.push_back(direction(0.0, 0.0, {0.0, 0.0, 0.0}, {0.0, 0.0}));
    const auto independent = epcsaft::equilibrium::solved_state_first_actions(
        compiled, central, fixture.observables, independent_directions);
    bool same_a1 = independent.size() == independent_directions.size();
    double worst_a1 = 0.0;
    std::size_t factorizations = 0;
    for (std::size_t j = 0; j < independent.size() && same_a1; ++j) {
        const auto separate = epcsaft::equilibrium::solved_state_actions(
            compiled, central, {fixture.observables, {independent_directions[j]}});
        factorizations += independent[j].diagnostics.state_factorizations;
        same_a1 = independent[j].results.size() == separate.results.size();
        for (std::size_t i = 0; i < separate.results.size() && same_a1; ++i) {
            const auto& actual = independent[j].results[i];
            const auto& expected = separate.results[i];
            same_a1 = actual.order == 1 && actual.unit == expected.unit && actual.status == expected.status
                && actual.diagnostic_message == expected.diagnostic_message
                && actual.value.has_value() == expected.value.has_value()
                && actual.action.has_value() == expected.action.has_value();
            if (actual.value && expected.value)
                worst_a1 = std::max(worst_a1, std::abs(*actual.value - *expected.value) / std::max(1.0, std::abs(*expected.value)));
            if (actual.action && expected.action)
                worst_a1 = std::max(worst_a1, std::abs(*actual.action - *expected.action) / std::max(1.0, std::abs(*expected.action)));
        }
    }
    all_ok = (record_bool("e5-1/independent-A1", "e5-1", "four independent A1 responses, including zero, equal separate actions with one center factorization",
        same_a1 && worst_a1 <= 2.0e-12 && factorizations == 1,
        recording_support::fields({recording_support::field("maximum_scaled_change", worst_a1),
            recording_support::field("shared_factorizations", factorizations),
            recording_support::field("directions", independent_directions.size())}))) && all_ok;
    const auto reference_batch = actions(fixture, model, central, fixture.observables, 2);
    Fixture swapped = fixture;
    swapped.label = "phase-swap";
    std::swap(swapped.problem.phases[0], swapped.problem.phases[1]);
    std::swap(swapped.active[0], swapped.active[1]);
    for (auto& d : swapped.directions) std::swap(d.active_parameters[0], d.active_parameters[1]);
    for (auto& item : swapped.observables) if (item.kind != Kind::Temperature && item.kind != Kind::Pressure) item.phase = 1 - item.phase;
    const auto swapped_model = swapped.model();
    const auto swapped_batch = actions(swapped, swapped_model, epcsaft::equilibrium::solve_equilibrium(swapped_model, swapped.problem), swapped.observables, 2);
    Fixture permuted = fixture;
    permuted.label = "component-permutation";
    const std::vector<std::size_t> order = {2, 0, 1};
    permuted.parameters = equilibrium_support::permute_parameters(fixture.parameters, order);
    permuted.problem.feed.values = equilibrium_support::permute_vector(fixture.problem.feed.values, order);
    permuted.active = {{epcsaft::ParameterFamily::Kij, {5, 7}}, {epcsaft::ParameterFamily::SegmentCount, {1}}};
    for (auto& d : permuted.directions) d.feed_amounts_mol = equilibrium_support::permute_vector(d.feed_amounts_mol, order);
    for (auto& item : permuted.observables)
        if (item.kind == Kind::PhaseComponentAmount || item.kind == Kind::PhaseMoleFraction
            || item.kind == Kind::PhaseResidualChemicalPotentialOverRt || item.kind == Kind::PhaseLogActivity)
            item.component = static_cast<std::size_t>(std::find(order.begin(), order.end(), item.component) - order.begin());
    const auto permuted_model = permuted.model();
    const auto permuted_batch = actions(permuted, permuted_model, epcsaft::equilibrium::solve_equilibrium(permuted_model, permuted.problem), permuted.observables, 2);
    for (const auto* pair : {&swapped_batch, &permuted_batch}) {
        double worst = 0.0;
        bool statuses = true;
        for (std::size_t i = 0; i < reference_batch.results.size(); ++i) {
            const auto& a = reference_batch.results[i];
            const auto& b = pair->results[i];
            statuses = statuses && a.status == b.status && a.action.has_value() && b.action.has_value() && a.value.has_value() && b.value.has_value();
            if (statuses) worst = std::max(worst, std::max(std::abs(*a.action - *b.action) / std::max(1.0, std::abs(*a.action)),
                                                              std::abs(*a.value - *b.value) / std::max(1.0, std::abs(*a.value))));
        }
        const std::string label = pair == &swapped_batch ? "phase-swap" : "component-permutation";
        recording_support::compare(("e5-9/" + label).c_str(), "identity_difference", "1", worst, 0.0,
            recording_support::fields({recording_support::field("case", label), recording_support::field("results", reference_batch.results.size())}),
            "Relationship", "permutation", "A2 values and actions follow phase/component/active identity under permutation",
            "reference_or_one", 1.0e-8, 0.0, statuses && worst <= 1.0e-8, "e5-9", label.c_str(), "\\Delta", "identity", "\\Delta\\leq10^{-8}", "case");
        all_ok = (statuses && worst <= 1.0e-8) && all_ok;
    }
    // Warm reuse: a repeated identical request records nothing new on the shared tape.
    model.reset_counters();
    const auto again = actions(fixture, model, central, fixture.observables, 2);
    const auto counters = model.counters();
    all_ok = (record_bool("e5-9/warm-reuse", "e5-9", "a repeated request reuses the recorded tape and reproduces the batch",
        counters.recordings == 0 && again.results.size() == reference_batch.results.size(),
        recording_support::fields({recording_support::field("recordings", counters.recordings),
            recording_support::field("forward_sweeps", counters.forward), recording_support::field("directional_paths", counters.directional_paths)}))) && all_ok;
    CHECK(all_ok);
}

void check_e5_pins() {
    bool all_ok = true;
    using namespace e5;
    const auto positive = pinned_case(0.01);
    all_ok = (qualify(positive, "e5-2")) && all_ok;
    const auto model = positive.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, positive.problem);
    const std::vector<SolvedStateObservable> requested = {obs(Kind::PhaseAmount, 1), obs(Kind::PhaseVolume, 1), obs(Kind::PhaseMolarDensity, 1)};
    std::vector<SolvedStateActionBatch> batches;
    for (std::size_t order = 1; order <= 3; ++order) batches.push_back(actions(positive, model, central, requested, order));
    // Pinned amount: N_V = 0.01 N_feed, so A1 = 0.01 sum(dfeed) and A2/A3 vanish exactly.
    const double expected = 0.01 * (positive.directions[0].feed_amounts_mol[0] + positive.directions[0].feed_amounts_mol[1]);
    const bool amount_ok = std::abs(batches[0].results[0].action.value_or(1.0) - expected) <= STRUCTURAL_ZERO
        && std::abs(batches[1].results[0].action.value_or(1.0)) <= STRUCTURAL_ZERO
        && std::abs(batches[2].results[0].action.value_or(1.0)) <= STRUCTURAL_ZERO;
    recording_support::compare("e5-2/positive-pin/pinned-amount-A1", "action", "mol", batches[0].results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()),
        expected, recording_support::fields({recording_support::field("A2", epcsaft::core::optional_number(batches[1].results[0].action)),
            recording_support::field("A3", epcsaft::core::optional_number(batches[2].results[0].action))}),
        "Relationship", "pinned-amount", "pinned amount A1 is the fraction times the feed direction sum; A2 and A3 vanish",
        "strict_absolute", STRUCTURAL_ZERO, 0.0, amount_ok, "e5-2", "pinned-amount", "D N_V", "dfeed", "D N_V=0.01\\sum dn", "case");
    all_ok = (amount_ok) && all_ok;
    // Extensive volume follows the product rule V = N / rho from the returned N and rho actions.
    const auto& a1 = batches[0].results;
    const auto& a2 = batches[1].results;
    const double N = *a1[0].value, rho = *a1[2].value, Na = *a1[0].action, ra = *a1[2].action;
    const auto b1 = actions(positive, model, central, requested, 1);
    (void)b1;
    const double v1 = Na / rho - N * ra / (rho * rho);
    SolvedStateActionRequest second;
    second.observables = requested;
    second.directions = {positive.directions[1]};
    const auto other = epcsaft::equilibrium::solved_state_actions(epcsaft::equilibrium::compile_problem(model, positive.problem), central, second);
    const double Nb = *other.results[0].action, rb = *other.results[2].action, Nab = *a2[0].action, rab = *a2[2].action;
    const double v2 = Nab / rho - (Na * rb + Nb * ra) / (rho * rho) - N * rab / (rho * rho) + 2.0 * N * ra * rb / (rho * rho * rho);
    const bool volume_ok = std::abs(*a1[1].action - v1) <= 1.0e-10 * std::max(1.0, std::abs(v1))
        && std::abs(*a2[1].action - v2) <= 1.0e-10 * std::max(1.0, std::abs(v2));
    recording_support::compare("e5-2/positive-pin/volume-product-rule", "action", "m^3", *a2[1].action, v2,
        recording_support::fields({recording_support::field("A1_returned", *a1[1].action), recording_support::field("A1_product_rule", v1)}),
        "Relationship", "product-rule", "extensive volume actions equal the product rule of the returned amount and density actions",
        "reference_or_one", 0.0, 1.0e-10, volume_ok, "e5-2", "volume-product-rule", "D^2 V", "N,\\rho", "V=N/\\rho", "case");
    all_ok = (volume_ok) && all_ok;
    // Exactly zero pin: extensive rows are structural zero, intensive rows follow the ladders.
    auto zero = pinned_case(0.0);
    std::vector<SolvedStateObservable> intensive;
    for (const auto& item : zero.observables) {
        const bool extensive = item.phase == 1 && (item.kind == Kind::PhaseAmount || item.kind == Kind::PhaseComponentAmount || item.kind == Kind::PhaseVolume);
        if (!extensive) intensive.push_back(item);
    }
    const auto zero_model = zero.model();
    const auto zero_central = epcsaft::equilibrium::solve_equilibrium(zero_model, zero.problem);
    bool structural = zero_central.success;
    for (std::size_t order = 1; order <= 3 && structural; ++order) {
        const auto batch = actions(zero, zero_model, zero_central, {obs(Kind::PhaseAmount, 1), obs(Kind::PhaseComponentAmount, 1, 0), obs(Kind::PhaseVolume, 1)}, order);
        for (const auto& item : batch.results)
            structural = structural && item.status == SolvedStateActionStatus::Available && std::abs(item.value.value_or(1.0)) <= STRUCTURAL_ZERO
                && std::abs(item.action.value_or(1.0)) <= STRUCTURAL_ZERO;
    }
    all_ok = (record_bool("e5-2/zero-pin/structural-zero", "e5-2", "a zero pinned amount, its component amounts and its extensive volume are structural zeros through A3",
        structural, recording_support::fields({recording_support::field("success", zero_central.success ? 1U : 0U)}))) && all_ok;
    zero.observables = intensive;
    all_ok = (qualify(zero, "e5-2")) && all_ok;
    CHECK(all_ok);
}

namespace e5 {

Fixture reactive_case() {
    Fixture f;
    f.label = "reactive-two-liquid";
    f.parameters = equilibrium_support::synthetic_reactive_ternary();
    f.problem = equilibrium_support::synthetic_reactive_problem();
    f.directions = {direction(0.2, 1000.0, {0.001, -0.001, 0.0}), direction(-0.1, -500.0, {0.0, 0.001, -0.001}),
                    direction(0.05, 250.0, {-0.001, 0.0, 0.001})};
    f.observables = geometry_observables(f.problem, 3);
    return f;
}

Fixture isomerisation_case(double total, double a) {
    Fixture f;
    f.label = total == 1.0 ? "isomerisation" : "isomerisation-two-mol";
    f.parameters = e2::identical_ethane_pair();
    f.problem = equilibrium_support::flash(2);
    f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    f.problem.feed.values = {0.5 * total, 0.5 * total};
    f.problem.temperature = FreeValue::fixed(350.0);
    f.problem.pressure = FreeValue::fixed(2.0e5);
    f.problem.reactions = {{"A=B", {-1.0, 1.0}, std::nullopt, std::nullopt, legacy_isomerisation(1.0, a)}};
    f.directions = {direction(1.0, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0})};
    f.observables = {obs(Kind::PhaseMoleFraction, 0, 1), obs(Kind::PhaseComponentAmount, 0, 1)};
    f.scaled = false;
    return f;
}

}  // namespace e5

void check_e5_reactive() {
    bool all_ok = true;
    using namespace e5;
    const auto fixture = reactive_case();
    all_ok = (qualify(fixture, "e5-3")) && all_ok;
    {
        // At least one nonzero mixed action and the retained central affinities of both phases.
        const auto model = fixture.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
        const auto batch = actions(fixture, model, central, fixture.observables, 2);
        double largest = 0.0;
        for (const auto& item : batch.results) if (item.action.has_value()) largest = std::max(largest, std::abs(*item.action));
        std::vector<double> affinity(2, 0.0);
        for (std::size_t row = 0; row < central.rows.size(); ++row)
            if (central.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical)
                affinity[central.rows[row].phase] += fixture.problem.reactions[0].stoichiometry[central.rows[row].index] * central.residuals[row];
        all_ok = (record_bool("e5-3/reactive-two-liquid/mixed-nonzero", "e5-3", "the reactive two-liquid state carries a nonzero A2 mixed action with zero reaction affinity in both phases",
            largest > 1.0e-6 && std::abs(affinity[0]) <= 1.0e-9 && std::abs(affinity[1]) <= 1.0e-9,
            recording_support::fields({recording_support::field("largest_mixed_action", largest), recording_support::field("affinity", affinity)}))) && all_ok;
    }
    // Exact nonconstant-correlation isomerisation: p = K/(1+K) with analytic thermal derivatives.
    const double a = moderate_isomerisation_a(), T = 350.0;
    const double L1 = -1200.0 / (T * T) + 0.4 / T + 0.002, L2 = 2.0 * 1200.0 / (T * T * T) - 0.4 / (T * T), L3 = -6.0 * 1200.0 / (T * T * T * T) + 2.0 * 0.4 / (T * T * T);
    const double ln_k = a + 1200.0 / T + 0.4 * std::log(T) + 0.002 * T, p = std::exp(ln_k) / (1.0 + std::exp(ln_k));
    const double analytic[3] = {p * (1.0 - p) * L1, p * (1.0 - p) * (L2 + (1.0 - 2.0 * p) * L1 * L1),
        p * (1.0 - p) * (L3 + 3.0 * (1.0 - 2.0 * p) * L1 * L2 + (1.0 - 6.0 * p + 6.0 * p * p) * L1 * L1 * L1)};
    const double expected[3] = {-1.4784580498866213e-3, 8.4348805281749936e-6, -2.7304974653104138e-9};
    for (const double total : {1.0, 2.0}) {
        auto fixture_t = isomerisation_case(total, a);
        const auto model = fixture_t.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture_t.problem);
        all_ok = (central.success) && all_ok;
        for (std::size_t order = 1; order <= 3; ++order) {
            for (const bool rescaled : {false, true}) {
                auto directions = fixture_t;
                if (rescaled) {
                    directions.directions = order == 2 ? std::vector<SolvedStateActionDirection>{direction(1.0e-6, 0.0, {0.0, 0.0}), direction(1.0e6, 0.0, {0.0, 0.0})}
                        : std::vector<SolvedStateActionDirection>{direction(1.0e-6, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0}), direction(1.0e6, 0.0, {0.0, 0.0})};
                    if (order == 1) continue;
                }
                const auto batch = actions(directions, model, central, fixture_t.observables, order);
                const double bar = 2.0e-10 * std::abs(analytic[order - 1]) + 1.0e-14;
                const bool fraction_ok = batch.results[0].action.has_value() && std::abs(*batch.results[0].action - analytic[order - 1]) <= bar
                    && std::abs(analytic[order - 1] - expected[order - 1]) <= bar;
                const bool amount_ok = batch.results[1].action.has_value() && std::abs(*batch.results[1].action - total * analytic[order - 1]) <= total * bar;
                const auto id = std::string("e5-3/") + fixture_t.label + (rescaled ? "/rescaled" : "") + "/A" + std::to_string(order);
                recording_support::compare(id.c_str(), "action", "1/K^k", batch.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()), analytic[order - 1],
                    recording_support::fields({recording_support::field("total", total), recording_support::field("order", order), recording_support::field("rescaled", rescaled ? 1U : 0U),
                        recording_support::field("amount_action", epcsaft::core::optional_number(batch.results[1].action)), recording_support::field("bar", bar),
                        recording_support::field("status", status_field(batch, 0))}),
                    "Relationship", "analytic-correlation", "x_B thermal actions equal the analytic logistic derivatives of ln K(T); the B amount scales with the total",
                    "strict_absolute", bar, 0.0, fraction_ok && amount_ok, "e5-3", id.c_str(), "d^k x_B/dT^k", "T,K(T)", "p^{(k)}=p(1-p)[\\ldots]", "case");
                all_ok = (fraction_ok && amount_ok) && all_ok;
            }
        }
        if (total == 2.0) {
            // Mixed amount/temperature chain: n_B = N p(T) is degree one in the feed, x_B degree zero.
            auto mixed = fixture_t;
            mixed.directions = {direction(0.0, 0.0, {0.05, 0.05}), direction(1.0, 0.0, {0.0, 0.0})};
            const auto batch = actions(mixed, model, central, fixture_t.observables, 2);
            const bool ok = batch.results[0].action.has_value() && std::abs(*batch.results[0].action) <= 1.0e-12
                && batch.results[1].action.has_value() && std::abs(*batch.results[1].action - 0.1 * analytic[0]) <= 2.0e-10 * std::abs(analytic[0]) + 1.0e-14;
            recording_support::compare("e5-3/isomerisation-two-mol/amount-temperature", "action", "mol/K", batch.results[1].action.value_or(std::numeric_limits<double>::quiet_NaN()), 0.1 * analytic[0],
                recording_support::fields({recording_support::field("fraction_mixed", epcsaft::core::optional_number(batch.results[0].action))}),
                "Relationship", "analytic-correlation", "the amount/temperature mixed action is sum(dfeed) p'(T) while the composition mixed action vanishes",
                "strict_absolute", 2.0e-10 * std::abs(analytic[0]) + 1.0e-14, 0.0, ok, "e5-3", "amount-temperature", "D^2 n_B", "dn,dT", "D^2 n_B=\\sum dn\\,p'", "case");
            all_ok = (ok) && all_ok;
        }
    }
    {
        // Zero directions give exact zero actions after admission.
        auto zero = isomerisation_case(1.0, a);
        zero.directions = {direction(0.0, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0}), direction(0.0, 0.0, {0.0, 0.0})};
        const auto model = zero.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, zero.problem);
        const auto batch = actions(zero, model, central, zero.observables, 3);
        all_ok = (record_bool("e5-3/isomerisation/zero-direction", "e5-3", "a zero direction yields an available exact zero action with the central value retained",
            batch.results[0].status == SolvedStateActionStatus::Available && batch.results[0].action == 0.0 && batch.results[0].value.has_value(),
            recording_support::fields({recording_support::field("status", status_field(batch, 0)), recording_support::field("value", epcsaft::core::optional_number(batch.results[0].value))}))) && all_ok;
    }
    {
        // Fixed-T scalar K = 2: a conserved constant Gibbs shift leaves actions unchanged, but a
        // point-only supplied vector cannot supply the correlated thermal reference.
        auto scalar = isomerisation_case(1.0, a);
        scalar.label = "scalar-k";
        scalar.problem.reactions = {{"A=B", {-1.0, 1.0}, 2.0, std::nullopt}};
        scalar.directions = {direction(0.0, 1.0e4, {0.01, -0.01}), direction(0.0, -5.0e3, {0.0, 0.005}), direction(0.0, 2.0e3, {0.003, 0.0})};
        scalar.observables = geometry_observables(scalar.problem, 2);
        auto shifted = scalar;
        shifted.label = "scalar-k-shifted";
        // nu.g = -RT ln K with a conserved constant offset of 1500 J/mol on both species.
        shifted.problem.standard_gibbs_j_per_mol = {1500.0 + epcsaft::GAS_CONSTANT_J_PER_MOL_K * 350.0 * std::log(2.0) / 2.0,
                                                    1500.0 - epcsaft::GAS_CONSTANT_J_PER_MOL_K * 350.0 * std::log(2.0) / 2.0};
        const auto model = scalar.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, scalar.problem);
        const auto shifted_central = epcsaft::equilibrium::solve_equilibrium(model, shifted.problem);
        const auto batch = actions(scalar, model, central, scalar.observables, 2);
        const auto shifted_batch = actions(shifted, model, shifted_central, shifted.observables, 2);
        double worst = 0.0;
        bool available = true;
        for (std::size_t i = 0; i < batch.results.size(); ++i) {
            const auto& x = batch.results[i]; const auto& y = shifted_batch.results[i];
            available = available && x.action.has_value() && y.action.has_value();
            if (available) worst = std::max(worst, std::abs(*x.action - *y.action) / std::max(1.0, std::abs(*x.action)));
        }
        recording_support::compare("e5-3/scalar-k/gauge", "action_difference", "1", worst, 0.0,
            recording_support::fields({recording_support::field("results", batch.results.size())}),
            "Relationship", "gauge", "a conserved constant standard-state shift leaves fixed-T A2 actions unchanged",
            "reference_or_one", 1.0e-8, 0.0, available && worst <= 1.0e-8, "e5-3", "gauge", "\\Delta", "g", "\\Delta\\leq10^{-8}", "case");
        all_ok = (available && worst <= 1.0e-8) && all_ok;
        auto thermal = shifted;
        thermal.directions = {direction(1.0, 0.0, {0.0, 0.0})};
        const auto thermal_batch = actions(thermal, model, shifted_central, thermal.observables, 1);
        bool typed = true;
        for (const auto& item : thermal_batch.results) {
            const bool fixed_input = item.observable.kind == Kind::Temperature || item.observable.kind == Kind::Pressure;
            typed = typed && item.value.has_value()
                && (fixed_input ? item.status == SolvedStateActionStatus::Available && item.action.has_value()
                                : item.status == SolvedStateActionStatus::MissingThermochemistry && !item.action.has_value());
        }
        all_ok = (record_bool("e5-3/scalar-k/point-only-thermal", "e5-3", "a point-only supplied Gibbs vector returns missing_thermochemistry for thermal state actions and keeps fixed-input identities",
            typed, recording_support::fields({recording_support::field("status", status_field(thermal_batch, 2))}))) && all_ok;
        const auto compiled = epcsaft::equilibrium::compile_problem(model, shifted.problem);
        const auto nonthermal = epcsaft::equilibrium::solved_state_actions(
            compiled, shifted_central, {shifted.observables, {scalar.directions[0]}});
        const auto independent = epcsaft::equilibrium::solved_state_first_actions(
            compiled, shifted_central, shifted.observables, {thermal.directions[0], scalar.directions[0]});
        bool isolated = independent.size() == 2;
        double worst_a1 = 0.0;
        for (std::size_t j = 0; j < independent.size() && isolated; ++j) {
            const auto& expected = j == 0 ? thermal_batch : nonthermal;
            isolated = independent[j].results.size() == expected.results.size();
            for (std::size_t i = 0; i < expected.results.size() && isolated; ++i) {
                const auto& actual = independent[j].results[i];
                const auto& single = expected.results[i];
                isolated = actual.status == single.status && actual.value.has_value() == single.value.has_value()
                    && actual.action.has_value() == single.action.has_value();
                if (actual.action && single.action)
                    worst_a1 = std::max(worst_a1, std::abs(*actual.action - *single.action) / std::max(1.0, std::abs(*single.action)));
            }
        }
        all_ok = (record_bool("e5-3/scalar-k/independent-thermal", "e5-3", "a refused thermal A1 does not suppress a following available feed A1 at the same center",
            isolated && worst_a1 <= 2.0e-12 && independent[0].results[2].status == SolvedStateActionStatus::MissingThermochemistry
                && independent[1].results[2].status == SolvedStateActionStatus::Available,
            recording_support::fields({recording_support::field("maximum_scaled_change", worst_a1)}))) && all_ok;
    }
    {
        // A -> 2B with matching species thermochemistry: thermal A1-A3 of composition, amounts and calorics.
        const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K, tref = 300.0;
        const epcsaft::ReactionLogPolynomial density_law{-10.0 / R - std::log(R * tref / 1.0e5), 3000.0 / R, 10.0 / R - 1.0, 0.0,
            tref, 250.0, 450.0, epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
        Fixture f;
        f.label = "a-to-2b";
        f.parameters = e2::identical_ethane_pair();
        f.record = e2::species_record(30.0, 20.0);
        f.problem = equilibrium_support::flash(2);
        f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
        f.problem.feed.values = {0.5, 0.5};
        f.problem.temperature = FreeValue::fixed(350.0);
        f.problem.pressure = FreeValue::fixed(2.0e5);
        f.problem.reactions = {{"A=2B", {-1.0, 2.0}, std::nullopt, std::nullopt, density_law}};
        f.directions = {direction(20.0, 0.0, {0.0, 0.0}), direction(-10.0, 0.0, {0.0, 0.0}), direction(5.0, 0.0, {0.0, 0.0})};
        f.observables = {obs(Kind::PhaseMoleFraction, 0, 1), obs(Kind::PhaseComponentAmount, 0, 0), obs(Kind::PhaseComponentAmount, 0, 1), obs(Kind::PhaseAmount, 0),
            obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity)};
        f.scaled = false;
        // E2's independent pure-state relation y^2/(1-y) = K_P P0 / (phi P) at each corner
        // temperature supplies the composition and amount references; calorics use fresh solves.
        const epcsaft::ReactionLogPolynomial pressure_law{-10.0 / R, 3000.0 / R, 10.0 / R, 0.0, tref, 250.0, 450.0, "pressure-basis-P0"};
        const auto ethane = std::make_shared<epcsaft::HelmholtzModel>(properties_support::ethane());
        f.analytic = [ethane, pressure_law, R](const SolvedStateObservable& item, double T) -> std::optional<double> {
            if (item.kind != Kind::PhaseMoleFraction && item.kind != Kind::PhaseComponentAmount && item.kind != Kind::PhaseAmount) return std::nullopt;
            const auto pure = ethane->state_at_pressure(T, 2.0e5, {1.0}, epcsaft::DensityBranch::Vapor);
            const double compressibility = 2.0e5 / (pure.molar_density * R * T);
            const double phi = std::exp(pure.residual_chemical_potential_over_rt[0]) / compressibility;
            const double a = std::exp(epcsaft::evaluate_reaction_properties(pressure_law, T).ln_k) * 1.0e5 / (phi * 2.0e5);
            const double y = 2.0 * a / (a + std::sqrt(a * a + 4.0 * a));
            const double total = 1.5 / (2.0 - y);
            if (item.kind == Kind::PhaseMoleFraction) return item.component == 1 ? y : 1.0 - y;
            if (item.kind == Kind::PhaseAmount) return total;
            return item.component == 1 ? total * y : total * (1.0 - y);
        };
        all_ok = (qualify(f, "e5-3")) && all_ok;
    }
    CHECK(all_ok);
}

namespace e5 {

Fixture associating_case() {
    Fixture f;
    f.label = "associating";
    f.parameters = properties_support::associating_binary();
    f.active = {{epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}}, {epcsaft::ParameterFamily::AssociationVolume, {1, 4}}};
    f.record = properties_support::p3_binary_record();
    f.problem = equilibrium_support::flash(2);
    f.problem.feed.values = {0.25, 0.75};
    f.problem.temperature = FreeValue::fixed(313.15);
    f.problem.pressure = FreeValue::fixed(epcsaft::HelmholtzModel(f.parameters).state(313.15, 42441.09427799721, {0.25, 0.75}).pressure);
    f.directions = {direction(0.2, 15000.0, {-0.01, 0.01}, {0.03, 0.001}), direction(-0.1, -10000.0, {0.005, -0.005}, {-0.02, -0.0005}),
                    direction(0.05, 6000.0, {0.002, -0.002}, {0.01, 0.0002})};
    f.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy),
                     obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity), obs(Kind::PhaseAmount, 0), obs(Kind::PhaseMolarDensity, 0),
                     obs(Kind::PhaseLogActivity, 0, 0), obs(Kind::PhaseLogActivity, 0, 1)};
    f.scaled = false;  // owner decision 2026-09-16: the scaled A3 corners leave the feed domain
    return f;
}

}  // namespace e5

void check_e5_associating() {
    bool all_ok = true;
    using namespace e5;
    const auto fixture = associating_case();
    all_ok = (qualify(fixture, "e5-4")) && all_ok;
    const auto model = fixture.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    // Selected results equal the corresponding full results.
    const auto full = actions(fixture, model, central, fixture.observables, 3);
    double worst = 0.0;
    bool ok = true;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto selected = actions(fixture, model, central, {fixture.observables[i]}, 3);
        ok = ok && selected.results[0].action.has_value() && full.results[i].action.has_value();
        if (ok) worst = std::max(worst, std::abs(*selected.results[0].action - *full.results[i].action) / std::max(1.0, std::abs(*full.results[i].action)));
    }
    recording_support::compare("e5-4/associating/selected-vs-full", "action_difference", "1", worst, 0.0,
        recording_support::fields({recording_support::field("site_rank", full.results[0].property_diagnostics.has_value() ? full.results[0].property_diagnostics->site_rank : 0),
            recording_support::field("energy_order", full.diagnostics.maximum_energy_order)}),
        "Relationship", "selection", "selected A3 property results equal their full-batch values", "reference_or_one", 2.0e-12, 0.0,
        ok && worst <= 2.0e-12, "e5-4", "selected-vs-full", "\\Delta", "selection", "\\Delta\\leq2\\times10^{-12}", "case");
    all_ok = (ok && worst <= 2.0e-12) && all_ok;
    // Mixed-action symmetry: the direction order changes neither A2 nor A3.
    auto rotated = fixture;
    rotated.directions = {fixture.directions[1], fixture.directions[2], fixture.directions[0]};
    auto swapped_pair = fixture;
    swapped_pair.directions = {fixture.directions[1], fixture.directions[0]};
    const auto natural_two = actions(fixture, model, central, fixture.observables, 2);
    const auto swapped_two = actions(swapped_pair, model, central, fixture.observables, 2);
    const auto rotated_three = actions(rotated, model, central, fixture.observables, 3);
    double asymmetry = 0.0;
    bool symmetric = true;
    for (std::size_t i = 0; i < fixture.observables.size(); ++i) {
        symmetric = symmetric && swapped_two.results[i].action.has_value() && natural_two.results[i].action.has_value()
            && rotated_three.results[i].action.has_value() && full.results[i].action.has_value();
        if (!symmetric) break;
        asymmetry = std::max({asymmetry, std::abs(*swapped_two.results[i].action - *natural_two.results[i].action) / std::max(1.0, std::abs(*natural_two.results[i].action)),
            std::abs(*rotated_three.results[i].action - *full.results[i].action) / std::max(1.0, std::abs(*full.results[i].action))});
    }
    recording_support::compare("e5-4/associating/mixed-symmetry", "action_difference", "1", asymmetry, 0.0,
        recording_support::fields({recording_support::field("results", fixture.observables.size())}),
        "Relationship", "symmetry", "A2 and A3 mixed actions are invariant under direction-order permutation", "reference_or_one", 1.0e-10, 0.0,
        symmetric && asymmetry <= 1.0e-10, "e5-4", "mixed-symmetry", "\\Delta", "order", "\\Delta\\leq10^{-10}", "case");
    all_ok = (symmetric && asymmetry <= 1.0e-10) && all_ok;
    // Radial identities: molar h has Dh[r] = 0 and D2h[v, r] = -Dh[v]; extensive H = N h has DH[r] = H, D2H[v, r] = 0.
    auto radial = fixture;
    radial.directions = {direction(0.0, 0.0, {0.02, -0.01}, {0.0, 0.0}), direction(0.0, 0.0, {0.25, 0.75}, {0.0, 0.0})};
    const std::vector<SolvedStateObservable> requested = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy), obs(Kind::PhaseAmount, 0)};
    const auto dv = actions(radial, model, central, requested, 1);
    auto only_r = radial;
    only_r.directions = {radial.directions[1]};
    const auto dr = actions(only_r, model, central, requested, 1);
    const auto dvr = actions(radial, model, central, requested, 2);
    const double h = *dv.results[0].value, N = *dv.results[1].value;
    const double hv = *dv.results[0].action, hr = *dr.results[0].action, hvr = *dvr.results[0].action;
    const double Nv = *dv.results[1].action, Nr = *dr.results[1].action, Nvr = *dvr.results[1].action;
    const double tolerance = 1.0e-8 * std::max(1.0, std::abs(hv));
    const bool molar_ok = std::abs(hr) <= tolerance && std::abs(hvr + hv) <= tolerance;
    const double H_r = Nr * h + N * hr, H_vr = Nvr * h + Nv * hr + Nr * hv + N * hvr;
    const bool extensive_ok = std::abs(H_r - N * h) <= 1.0e-8 * std::max(1.0, std::abs(N * h)) && std::abs(H_vr) <= 1.0e-8 * std::max(1.0, std::abs(N * hv));
    recording_support::compare("e5-4/associating/radial-molar", "action", "J/mol", hvr, -hv,
        recording_support::fields({recording_support::field("Dh_r", hr), recording_support::field("Dh_v", hv), recording_support::field("H_r", H_r),
            recording_support::field("H", N * h), recording_support::field("H_vr", H_vr)}),
        "Relationship", "homogeneity", "molar enthalpy obeys Dh[r]=0 and D2h[v,r]=-Dh[v]; the extensive product H=N h obeys DH[r]=H and D2H[v,r]=0",
        "reference_or_one", 0.0, 1.0e-8, molar_ok && extensive_ok, "e5-4", "radial", "D^2 h[v,r]", "v,r", "D^2h[v,r]=-Dh[v]", "case");
    all_ok = (molar_ok && extensive_ok) && all_ok;
    CHECK(all_ok);
}

namespace e5 {

Fixture electrolyte_case(bool with_record) {
    Fixture f;
    f.label = with_record ? "electrolyte" : "electrolyte-residual-only";
    f.parameters = properties_support::aqueous_electrolyte();
    if (with_record) {
        epcsaft::ThermochemistryRecord record;
        record.standard_pressure = 1.0e5;
        for (std::size_t i = 0; i < 3; ++i) record.component_intervals.push_back(properties_support::p3_record(properties_support::p3_forms()[0].second).component_intervals.front());
        f.record = record;
    }
    f.problem = equilibrium_support::flash(3);
    f.problem.feed.values = {0.98, 0.01, 0.01};
    f.problem.temperature = FreeValue::fixed(298.15);
    f.directions = {direction(0.1, 1000.0, {-0.002, 0.001, 0.001}), direction(-0.05, -500.0, {0.001, -0.0005, -0.0005}),
                    direction(0.02, 250.0, {-0.0004, 0.0002, 0.0002})};
    f.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualIsobaricHeatCapacity)};
    if (with_record) {
        f.observables.push_back(obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy));
        f.observables.push_back(obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity));
    }
    f.scaled = false;
    return f;
}

Fixture two_root_case(PhaseKind kind) {
    Fixture f;
    f.label = kind == PhaseKind::Vapor ? "two-root-vapor" : "two-root-liquid";
    f.parameters = properties_support::ethane();
    f.problem = equilibrium_support::flash(1);
    f.problem.phases = {equilibrium_support::phase("P", kind)};
    f.problem.temperature = FreeValue::fixed(200.0);
    f.problem.pressure = FreeValue::fixed(217800.69494692178);
    f.directions = {direction(0.1, 2000.0, {0.0}), direction(-0.05, -1000.0, {0.0}), direction(0.03, 500.0, {0.0})};
    f.observables = {obs(Kind::PhasePressure, 0), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy),
                     obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualIsobaricHeatCapacity), obs(Kind::PhaseMolarDensity, 0)};
    f.corner_branch = kind == PhaseKind::Vapor ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Liquid;
    return f;
}

Fixture thermal_case(double temperature, bool with_record) {
    Fixture f;
    f.label = "thermal-" + recording_support::number(temperature) + (with_record ? "" : "-no-record");
    f.parameters = properties_support::ethane();
    if (with_record) f.record = properties_support::p3_record(properties_support::p3_forms()[1].second);
    f.problem = equilibrium_support::flash(1);
    f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    f.problem.temperature = FreeValue::fixed(temperature);
    f.directions = {direction(20.0, 0.0, {0.0}), direction(-10.0, 0.0, {0.0}), direction(5.0, 0.0, {0.0})};
    f.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::IdealEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy),
                     obs(Kind::PhaseProperty, 0, 0, PropertyObservable::IdealIsobaricHeatCapacity), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity)};
    f.scaled = false;
    return f;
}

}  // namespace e5

void check_e5_electrolyte() {
    bool all_ok = true;
    using namespace e5;
    auto with_record = electrolyte_case(true);
    for (std::size_t i = 0; i < 3; ++i) with_record.observables.push_back(obs(Kind::PhaseLogActivity, 0, i));  // solvent and both ions
    all_ok = (qualify(with_record, "e5-5")) && all_ok;
    // Residual-only without a record: residual siblings stay available, total rows are typed.
    auto residual_only = electrolyte_case(false);
    residual_only.observables.push_back(obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy));
    const auto model = residual_only.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, residual_only.problem);
    const auto batch = actions(residual_only, model, central, residual_only.observables, 2);
    const bool ok = batch.results[0].status == SolvedStateActionStatus::Available && batch.results[0].action.has_value()
        && batch.results[1].status == SolvedStateActionStatus::Available && batch.results[1].action.has_value()
        && batch.results[2].status == SolvedStateActionStatus::MissingThermochemistry && !batch.results[2].action.has_value()
        && !batch.results[2].value.has_value() && batch.results[2].diagnostic_message.has_value();
    all_ok = (record_bool("e5-5/electrolyte-residual-only/sibling", "e5-5", "a missing thermal record types the total enthalpy row while residual siblings stay available",
        ok, recording_support::fields({recording_support::field("status_total", status_field(batch, 2)), recording_support::field("status_residual", status_field(batch, 0)),
            recording_support::field("message", batch.results[2].diagnostic_message.value_or("absent"))}))) && all_ok;
    residual_only.observables.pop_back();
    all_ok = (qualify(residual_only, "e5-5")) && all_ok;
    CHECK(all_ok);
}

void check_e5_roots() {
    bool all_ok = true;
    using namespace e5;
    for (const PhaseKind kind : {PhaseKind::Vapor, PhaseKind::Liquid}) {
        const auto fixture = two_root_case(kind);
        const auto model = fixture.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
        const auto root = epcsaft::solve_density(model, 200.0, 217800.69494692178, {1.0},
            kind == PhaseKind::Vapor ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Liquid, std::nullopt);
        const bool branch_ok = central.success && std::abs(std::log(central.molar_densities[0] / root.molar_density)) <= 1.0e-8
            && root.stable_root_count.value_or(0) == 2;
        all_ok = (record_bool((std::string("e5-6/") + fixture.label + "/central-root").c_str(), "e5-6",
            "the declared branch selects its own root of the two stable ethane roots", branch_ok,
            recording_support::fields({recording_support::field("rho", central.success ? central.molar_densities[0] : 0.0),
                recording_support::field("rho_branch", root.molar_density), recording_support::field("stable_roots", root.stable_root_count.value_or(0))}))) && all_ok;
        all_ok = (qualify(fixture, "e5-6")) && all_ok;
    }
    CHECK(all_ok);
}

void check_e5_thermal() {
    bool all_ok = true;
    using namespace e5;
    all_ok = (qualify(thermal_case(320.0, true), "e5-7")) && all_ok;
    const auto endpoint = thermal_case(250.0, true);
    const auto model = endpoint.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, endpoint.problem);
    const auto batch = actions(endpoint, model, central, endpoint.observables, 1);
    bool boundary = central.success;
    for (const auto& item : batch.results)
        boundary = boundary && item.value.has_value() && item.status == SolvedStateActionStatus::ThermalBoundary && !item.action.has_value();
    all_ok = (record_bool("e5-7/thermal-250/boundary", "e5-7", "at the closed interval endpoint central values are retained and thermal actions are typed thermal_boundary",
        boundary, recording_support::fields({recording_support::field("status", status_field(batch, 0)), recording_support::field("value", epcsaft::core::optional_number(batch.results[0].value))}))) && all_ok;
    const auto outside = thermal_case(249.0, true);
    const auto outside_model = outside.model();
    const auto outside_central = epcsaft::equilibrium::solve_equilibrium(outside_model, outside.problem);
    const auto outside_batch = actions(outside, outside_model, outside_central, outside.observables, 1);
    bool interval = outside_central.success;
    for (const auto& item : outside_batch.results) interval = interval && item.status == SolvedStateActionStatus::OutsideInterval && !item.action.has_value();
    all_ok = (record_bool("e5-7/thermal-249/outside-interval", "e5-7", "below the closed interval dependent ideal and total rows are typed outside_interval",
        interval, recording_support::fields({recording_support::field("status", status_field(outside_batch, 0))}))) && all_ok;
    const auto missing = thermal_case(320.0, false);
    const auto missing_model = missing.model();
    const auto missing_central = epcsaft::equilibrium::solve_equilibrium(missing_model, missing.problem);
    const auto missing_batch = actions(missing, missing_model, missing_central, missing.observables, 1);
    bool typed = missing_central.success;
    for (const auto& item : missing_batch.results) typed = typed && item.status == SolvedStateActionStatus::MissingThermochemistry;
    all_ok = (record_bool("e5-7/thermal-no-record/missing", "e5-7", "without a thermal record dependent rows are typed missing_thermochemistry",
        typed, recording_support::fields({recording_support::field("status", status_field(missing_batch, 0))}))) && all_ok;
    CHECK(all_ok);
}

void check_e5_refusals() {
    bool all_ok = true;
    using namespace e5;
    const auto fixture = neutral_case();
    const auto model = fixture.model();
    const auto compiled = epcsaft::equilibrium::compile_problem(model, fixture.problem);
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    all_ok = (central.success) && all_ok;
    const auto invalid = [&](const char* id, const char* phrase, const SolvedStateActionRequest& request, const CompiledSystem& system, const EquilibriumResult& result) {
        std::string message;
        try { (void)epcsaft::equilibrium::solved_state_actions(system, result, request); }
        catch (const std::invalid_argument& error) { message = error.what(); }
        all_ok = (record_bool((std::string("e5-8/invalid/") + id).c_str(), "e5-8", phrase, message.find(phrase) != std::string::npos,
            recording_support::fields({recording_support::field("message", message)}))) && all_ok;
    };
    SolvedStateActionRequest base;
    base.observables = {obs(Kind::PhaseAmount, 0)};
    base.directions = {fixture.directions[0]};
    auto empty = base; empty.observables.clear();
    invalid("no-observables", "no solved-state observables", empty, compiled, central);
    auto four = base; four.directions.assign(4, fixture.directions[0]);
    invalid("order-four", "solved-state action order", four, compiled, central);
    auto short_feed = base; short_feed.directions[0].feed_amounts_mol = {0.0, 0.0};
    invalid("feed-size", "solved-state action direction", short_feed, compiled, central);
    auto duplicate = base; duplicate.observables.push_back(obs(Kind::PhaseAmount, 0));
    invalid("duplicate", "duplicate solved-state observable", duplicate, compiled, central);
    auto ordinal = base; ordinal.observables = {obs(Kind::PhaseAmount, 5)};
    invalid("phase-ordinal", "solved-state observable phase ordinal", ordinal, compiled, central);
    auto unsupported = base; unsupported.observables = {obs(Kind::PhaseResidualChemicalPotentialOverRt, 0, 2)};
    auto restricted = fixture.problem; restricted.phases[0].support = {0, 1};
    const auto restricted_compiled = epcsaft::equilibrium::compile_problem(model, restricted);
    invalid("unsupported-potential", "residual chemical potential outside the phase support", unsupported, restricted_compiled, epcsaft::equilibrium::solve_equilibrium(model, restricted));
    unsupported.observables = {obs(Kind::PhaseLogActivity, 0, 2)};
    invalid("unsupported-log-activity", "log activity outside the phase support", unsupported, restricted_compiled, epcsaft::equilibrium::solve_equilibrium(model, restricted));
    auto free_t = base;
    const auto pinned = pinned_case(0.01);
    const auto pinned_model = pinned.model();
    const auto pinned_compiled = epcsaft::equilibrium::compile_problem(pinned_model, pinned.problem);
    const auto pinned_central = epcsaft::equilibrium::solve_equilibrium(pinned_model, pinned.problem);
    free_t.directions = {direction(0.0, 1.0, {0.0, 0.0})};
    invalid("free-pressure-control", "free pressure is a solved output", free_t, pinned_compiled, pinned_central);
    auto foreign = free_t; foreign.directions = {direction(0.0, 0.0, {0.0, 0.0})};
    invalid("mismatched-central", "central result does not match the compiled system", foreign, pinned_compiled, central);
    // Typed refusals with retained central values: a same-layout result of another problem
    // does not satisfy this problem's defining rows.
    auto other = fixture.problem; other.temperature.value = 181.0;
    const auto other_compiled = epcsaft::equilibrium::compile_problem(model, other);
    const auto other_batch = epcsaft::equilibrium::solved_state_actions(other_compiled, central, base);
    const double other_rows = other_compiled.evaluate(central.coordinates).residual_norm();
    all_ok = (record_bool("e5-8/typed/defining-equations", "e5-8", "an accepted result of another same-layout problem is typed defining_equations_not_satisfied with its mapped value retained",
        other_rows > 1.0e-10 && other_batch.results[0].status == SolvedStateActionStatus::DefiningEquationsNotSatisfied
            && other_batch.results[0].value.has_value() && !other_batch.results[0].action.has_value(),
        recording_support::fields({recording_support::field("central_rows", other_rows), recording_support::field("status", status_field(other_batch, 0))}))) && all_ok;
    auto exact = equilibrium_support::flash(2);
    exact.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, AmountSpec::pinned(0.5)), equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.5))};
    exact.temperature = FreeValue::variable(300.0); exact.pressure = FreeValue::variable(1.0e5);
    const auto& binary_model = pinned_model;
    const auto exact_compiled = epcsaft::equilibrium::compile_problem(binary_model, exact);
    const auto exact_central = epcsaft::equilibrium::solve_equilibrium(binary_model, exact);
    SolvedStateActionRequest exact_request;
    exact_request.observables = {obs(Kind::PhaseAmount, 0)};
    exact_request.directions = {direction(0.0, 0.0, {0.01, 0.0})};
    const auto not_run = epcsaft::equilibrium::solved_state_actions(exact_compiled, exact_central, exact_request);
    all_ok = (record_bool("e5-8/typed/central-unavailable", "e5-8", "the pre-solve NOT_RUN record is typed central_solve_unavailable",
        !exact_central.success && not_run.results[0].status == SolvedStateActionStatus::CentralSolveUnavailable && !not_run.results[0].value.has_value(),
        recording_support::fields({recording_support::field("solver_status", exact_central.solver_status)}))) && all_ok;
    // A zero-feed non-participant species leaves a rank-deficient balance row (E2), so no
    // admitted central state can carry a removed species: support_change and ionic_insertion
    // are unforced defensive statuses, retained as claim limits rather than physical cases.
    auto sparse = equilibrium_support::flash(3); sparse.feed.values = {0.5, 0.5, 0.0}; sparse.temperature.value = 180.0; sparse.pressure.value = 2.0e6;
    const auto sparse_compiled = epcsaft::equilibrium::compile_problem(model, sparse);
    const auto sparse_central = epcsaft::equilibrium::solve_equilibrium(model, sparse);
    SolvedStateActionRequest insertion;
    insertion.observables = {obs(Kind::PhaseAmount, 0), obs(Kind::Temperature)};
    insertion.directions = {direction(0.0, 0.0, {0.0, 0.0, 0.01}, {0.0, 0.0})};
    const auto inserted = epcsaft::equilibrium::solved_state_actions(sparse_compiled, sparse_central, insertion);
    all_ok = (record_bool("e5-8/typed/removed-species", "e5-8", "a feed with a removed neutral species has no admitted central state (NOT_RUN rank record); E5 types it central_solve_unavailable, so support_change is unforced",
        !sparse_central.success && sparse_central.solver_status == "NOT_RUN" && !sparse_compiled.removed_species.empty()
            && inserted.results[0].status == SolvedStateActionStatus::CentralSolveUnavailable && inserted.results[1].status == SolvedStateActionStatus::CentralSolveUnavailable,
        recording_support::fields({recording_support::field("status", status_field(inserted, 0)), recording_support::field("message", sparse_central.message),
            recording_support::field("removed_species", sparse_compiled.removed_species.size())}))) && all_ok;
    epcsaft::HelmholtzModel salt(properties_support::aqueous_electrolyte());
    auto water = equilibrium_support::flash(3); water.feed.values = {1.0, 0.0, 0.0}; water.temperature = FreeValue::fixed(298.15);
    const auto water_compiled = epcsaft::equilibrium::compile_problem(salt, water);
    const auto water_central = epcsaft::equilibrium::solve_equilibrium(salt, water);
    SolvedStateActionRequest ionic_request;
    ionic_request.observables = {obs(Kind::PhaseMolarDensity, 0)};
    ionic_request.directions = {direction(0.0, 0.0, {-0.02, 0.01, 0.01})};
    const auto ionic_batch = epcsaft::equilibrium::solved_state_actions(water_compiled, water_central, ionic_request);
    all_ok = (record_bool("e5-8/typed/removed-ions", "e5-8", "a neutral-limit feed with removed ions has no admitted central state (NOT_RUN rank record); E5 types it central_solve_unavailable, so ionic_insertion is unforced",
        !water_central.success && water_central.solver_status == "NOT_RUN" && ionic_batch.results[0].status == SolvedStateActionStatus::CentralSolveUnavailable,
        recording_support::fields({recording_support::field("status", status_field(ionic_batch, 0)), recording_support::field("message", water_central.message)}))) && all_ok;
    // The matrix owner on a rank-deficient Jacobian: the exact all-pinned system loses one
    // rank at its generic point; nonfinite and excessive-residual paths remain unexecuted.
    const auto inspected = epcsaft::equilibrium::inspect_sensitivity(exact_compiled, exact_compiled.generic_point(), {});
    all_ok = (record_bool("e5-8/matrix/rank-deficient", "e5-8", "the shared matrix owner reports the rank-deficient all-pinned Jacobian as singular with its rank below the dimension",
        !inspected.nonsingular && inspected.rank + 1 == exact_compiled.generic_point().size() && inspected.rank == exact_compiled.jacobian_rank,
        recording_support::fields({recording_support::field("rank", inspected.rank), recording_support::field("dimension", exact_compiled.generic_point().size()),
            recording_support::field("condition", inspected.condition_estimate)}))) && all_ok;
    // Central preservation: an unavailable property action keeps the central value and its available siblings.
    auto partial = electrolyte_case(false);
    partial.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy), obs(Kind::PhaseMolarDensity, 0), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy)};
    const auto partial_model = partial.model();
    const auto partial_central = epcsaft::equilibrium::solve_equilibrium(partial_model, partial.problem);
    const auto partial_batch = actions(partial, partial_model, partial_central, partial.observables, 1);
    all_ok = (record_bool("e5-8/typed/sibling-preservation", "e5-8", "an unavailable ideal dependency leaves the density and residual siblings available with finite values",
        partial_batch.results[0].status == SolvedStateActionStatus::MissingThermochemistry && partial_batch.results[1].status == SolvedStateActionStatus::Available
            && partial_batch.results[1].action.has_value() && partial_batch.results[2].status == SolvedStateActionStatus::Available && partial_batch.results[2].action.has_value(),
        recording_support::fields({recording_support::field("statuses", status_field(partial_batch, 0) + "," + status_field(partial_batch, 1) + "," + status_field(partial_batch, 2))}))) && all_ok;
    CHECK(all_ok);
}

// S3 (#53): exact identities of solved-state actions on manufactured systems. Nonzero references
// use the retained CRITERION, zeros the retained STRUCTURAL_ZERO; no stencil and no new threshold.
void check_e5_identities() {
    bool all_ok = true;
    using namespace e5;
    const auto within = [&](const std::string& id, const char* claim, double value, double reference, double tolerance,
                            const SolvedStateActionBatch& batch, const std::string& stated) {
        const bool verdict = std::isfinite(value) && std::abs(value - reference) <= tolerance;
        recording_support::compare(id.c_str(), "action_identity", "1", value, reference,
            recording_support::fields({stated, recording_support::field("tolerance", tolerance),
                recording_support::field("condition_estimate", batch.diagnostics.condition_estimate.value_or(0.0)),
                recording_support::field("linear_residual", batch.diagnostics.maximum_linear_residual.value_or(0.0))}),
            "Relationship", "exact-identity", claim, "absolute", tolerance, 0.0, verdict, "e5-10");
        all_ok = verdict && all_ok;
    };
    const auto action = [](const SolvedStateActionBatch& batch, std::size_t i) { return batch.results[i].action.value_or(std::nan("")); };

    // H3: Euler homogeneity in the feed amounts, radial and mixed with a non-radial direction.
    auto neutral = neutral_case();
    const auto neutral_model = neutral.model();
    const auto neutral_central = epcsaft::equilibrium::solve_equilibrium(neutral_model, neutral.problem);
    neutral.directions = {direction(0.0, 0.0, neutral.problem.feed.values, {0.0, 0.0}), direction(0.0, 0.0, {-0.01, 0.006, 0.004}, {0.0, 0.0})};
    std::vector<SolvedStateObservable> list;
    std::vector<bool> extensive;
    for (std::size_t phase = 0; phase < 2; ++phase) {
        list.push_back(obs(Kind::PhaseAmount, phase)); extensive.push_back(true);
        list.push_back(obs(Kind::PhaseMolarDensity, phase)); extensive.push_back(false);
        for (std::size_t i = 0; i < 3; ++i) {
            list.push_back(obs(Kind::PhaseComponentAmount, phase, i)); extensive.push_back(true);
            list.push_back(obs(Kind::PhaseMoleFraction, phase, i)); extensive.push_back(false);
            list.push_back(obs(Kind::PhaseLogActivity, phase, i)); extensive.push_back(false);
        }
    }
    const auto euler = recording_support::fields({recording_support::field("n_feed", neutral.directions[0].feed_amounts_mol),
        recording_support::field("v", neutral.directions[1].feed_amounts_mol)});
    const auto radial = actions(neutral, neutral_model, neutral_central, list, 1);
    const auto mixed = actions(neutral, neutral_model, neutral_central, list, 2);
    std::swap(neutral.directions[0], neutral.directions[1]);
    const auto skew = actions(neutral, neutral_model, neutral_central, list, 1);
    double wrong_sign = 0.0;
    for (std::size_t i = 0; i < list.size(); ++i) {
        const std::string id = "e5-10/euler/" + epcsaft::equilibrium::solved_state_observable_id(list[i]);
        const double value = radial.results[i].value.value_or(std::nan("")), first = action(skew, i);
        if (extensive[i]) {
            within(id + "/radial", "degree one: D y[n] = y", action(radial, i), value, CRITERION[0] * std::abs(value), radial, euler);
            within(id + "/mixed", "degree one: D2 y[n, v] = 0", action(mixed, i), 0.0, STRUCTURAL_ZERO * std::max(1.0, std::abs(first)), mixed, euler);
        } else {
            within(id + "/radial", "degree zero: D y[n] = 0", action(radial, i), 0.0, STRUCTURAL_ZERO * std::max(1.0, std::abs(value)), radial, euler);
            within(id + "/mixed", "degree zero: D2 y[n, v] = -D y[v]", action(mixed, i), -first, CRITERION[1] * std::abs(first), mixed, euler);
            wrong_sign = std::max(wrong_sign, std::abs(action(mixed, i) - first) / std::abs(first));
        }
    }
    all_ok = (record_bool("e5-10/euler/wrong-sign", "e5-10", "the mixed identity with the wrong sign misses by a relative two",
        wrong_sign > 1.0, recording_support::fields({recording_support::field("relative_defect", wrong_sign)}))) && all_ok;

    // H4: in each liquid of A + B = C, sum_i nu_i ln a_i = ln K(T) along every admitted control path.
    const std::vector<double> nu = {-1.0, -1.0, 1.0};
    const double b = 1500.0, t0 = 298.15;  // labelled manufactured EOS-basis ln K = b (1/T - 1/T0), K(T0) = 1
    auto reactive = reactive_case();
    reactive.problem.standard_gibbs_j_per_mol.clear();
    reactive.problem.reactions[0].equilibrium_constant.reset();
    reactive.problem.reactions[0].correlation = epcsaft::ReactionLogPolynomial{-b / t0, b, 0.0, 0.0, t0, 250.0, 450.0,
        epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
    const auto reactive_model = reactive.model();
    const auto reactive_central = epcsaft::equilibrium::solve_equilibrium(reactive_model, reactive.problem);
    all_ok = reactive_central.success && all_ok;
    list.clear();
    for (std::size_t phase = 0; phase < 2; ++phase) for (std::size_t i = 0; i < 3; ++i) list.push_back(obs(Kind::PhaseLogActivity, phase, i));
    const auto manifold = [&](const char* label, std::vector<SolvedStateActionDirection> directions, double reference, bool values) {
        reactive.directions = std::move(directions);
        const std::size_t order = reactive.directions.size();
        const auto batch = actions(reactive, reactive_model, reactive_central, list, order);
        std::vector<double> d_t, d_p;
        for (const auto& item : reactive.directions) { d_t.push_back(item.temperature_k); d_p.push_back(item.pressure_pa); }
        const auto stated = recording_support::fields({recording_support::field("b", b), recording_support::field("T0", t0),
            recording_support::field("dT", d_t), recording_support::field("dP", d_p), recording_support::field("feed_direction", reactive.directions[0].feed_amounts_mol)});
        for (std::size_t phase = 0; phase < 2; ++phase) {
            double sum = 0.0, scale = 0.0, perturbed = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                const double term = values ? batch.results[3 * phase + i].value.value_or(std::nan("")) : action(batch, 3 * phase + i);
                sum += nu[i] * term; scale += std::abs(nu[i] * term); perturbed += (i == 2 ? 2.0 : nu[i]) * term;
            }
            const std::string id = std::string("e5-10/manifold/") + label + "/phase" + std::to_string(phase);
            // The central identity holds to the retained affinity bound of e5-3; actions to the retained criterion.
            const double tolerance = values ? 1.0e-9 : CRITERION[order - 1] * scale;
            within(id, "sum_i nu_i ln a_i follows ln K(T) in each phase", sum, reference, tolerance, batch, stated);
            all_ok = (record_bool((id + "/perturbed-nu").c_str(), "e5-10", "a wrong stoichiometric coefficient fails the identity's own criterion",
                std::abs(perturbed - reference) > tolerance, recording_support::fields({stated, recording_support::field("defect", perturbed - reference),
                    recording_support::field("tolerance", tolerance), recording_support::field("defect_over_tolerance", std::abs(perturbed - reference) / tolerance),
                    recording_support::field("term_scale", scale)}))) && all_ok;
        }
        return batch;
    };
    // Two phases, three species and one reaction leave no intensive freedom at fixed T and P, so a feed
    // direction moves only the phase amounts: the manifold identity would read 0 = 0 there. Pressure and
    // temperature directions do move the intensive state; K does not depend on pressure.
    const auto v = direction(0.0, 0.0, {0.001, -0.001, 0.0}), dP = direction(0.0, 1.0e5, {0.0, 0.0, 0.0}), dT = direction(1.0, 0.0, {0.0, 0.0, 0.0});
    const auto invariant = manifold("value", {v}, 0.0, true);
    manifold("pressure", {dP}, 0.0, false);
    manifold("pressure-temperature", {dP, dT}, 0.0, false);
    manifold("temperature", {dT}, -b / (t0 * t0), false);
    manifold("temperature-temperature", {dT, dT}, 2.0 * b / (t0 * t0 * t0), false);
    // Third order is where the weight j of the series logarithm first matters.
    manifold("temperature-temperature-temperature", {dT, dT, dT}, -6.0 * b / (t0 * t0 * t0 * t0), false);
    manifold("pressure-temperature-temperature", {dP, dT, dT}, 0.0, false);
    for (std::size_t i = 0; i < list.size(); ++i) {
        within("e5-10/phase-rule/" + epcsaft::equilibrium::solved_state_observable_id(list[i]), "a feed direction leaves the invariant intensive state unchanged",
               action(invariant, i), 0.0, STRUCTURAL_ZERO * std::max(1.0, std::abs(invariant.results[i].value.value_or(std::nan("")))), invariant,
               recording_support::fields({recording_support::field("feed_direction", v.feed_amounts_mol)}));
    }
    CHECK(all_ok);
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.equilibrium.f5", "Weighted_equilibrium_curvature", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "hessian/regular/*,hessian/pinned/*,jacobian/*,correlated-jacobian/*,f5-summary/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "F5"}, check_s10, "S10"},
        {{"native.equilibrium.s1", "Equilibrium_S1", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "residual/*,gibbs-dual/*,charge/*,affinity/*,feed-dual/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S1"}, check_s1},
        {{"native.equilibrium.s2", "Equilibrium_S2", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E1", "binary-success,binary-curvature,reactive-success,reactive-curvature,neutral-reference/*,reference-temperature/*,reference-calorics/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "S2"}, check_s2},
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
        {{"native.equilibrium.e5-neutral", "Equilibrium_E5_neutral", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-1/*,e5-9/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-1"}, check_e5_neutral},
        {{"native.equilibrium.e5-pins", "Equilibrium_E5_pins", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-2/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-2"}, check_e5_pins},
        {{"native.equilibrium.e5-reactive", "Equilibrium_E5_reactive", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-3/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-3"}, check_e5_reactive},
        {{"native.equilibrium.e5-associating", "Equilibrium_E5_associating", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-4/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-4"}, check_e5_associating},
        {{"native.equilibrium.e5-electrolyte", "Equilibrium_E5_electrolyte", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-5/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-5"}, check_e5_electrolyte},
        {{"native.equilibrium.e5-roots", "Equilibrium_E5_roots", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-6/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-6"}, check_e5_roots},
        {{"native.equilibrium.e5-thermal", "Equilibrium_E5_thermal", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-7/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-7"}, check_e5_thermal},
        {{"native.equilibrium.e5-refusals", "Equilibrium_E5_refusals", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "E5", "e5-8/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-8"}, check_e5_refusals},
        {{"native.equilibrium.e5-identities", "Equilibrium_E5_identities", "Numerical", "Equilibrium", "Explicit", "Numerical", nullptr, "S3", "e5-10/*", "engine/tests/native/equilibrium_checks.cpp", "equilibrium_checks", "E5-10"}, check_e5_identities},
    });
}
