#include "properties/helmholtz_model.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft {

namespace {

std::array<double, 2> ion_log_fugacity(const State& state, const std::array<std::size_t, 2>& ions) {
    if (!state.log_fugacity_coefficient.has_value())
        throw std::runtime_error("reference fugacity is unavailable");
    return {require_component_property((*state.log_fugacity_coefficient)[ions[0]]),
            require_component_property((*state.log_fugacity_coefficient)[ions[1]])};
}

constexpr std::array<double, 5> NEUTRAL_REFERENCE_MOLALITIES = {
    policy::REFERENCE_MOLALITY_ONE.value,
    policy::REFERENCE_MOLALITY_TWO.value,
    policy::REFERENCE_MOLALITY_THREE.value,
    policy::REFERENCE_MOLALITY_FOUR.value,
    policy::REFERENCE_CONVERGENCE_WITNESS.value,
};

double validate_neutral_reference_declaration(
    const Parameters& parameters,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry
) {
    const std::size_t components = parameters.component_count();
    if (!declaration.solvent_index.has_value()
        || *declaration.solvent_index >= components
        || declaration.solute_molality_ratios.size() != components
        || !std::isfinite(declaration.standard_molality_mol_per_kg)
        || declaration.standard_molality_mol_per_kg <= 0.0
        || stoichiometry.size() != components) {
        throw std::invalid_argument("invalid neutral reference declaration");
    }
    const std::size_t solvent = *declaration.solvent_index;
    if (parameters.charge_number[solvent] != 0) {
        throw std::invalid_argument("neutral reference solvent is charged");
    }
    const auto& molar_masses = parameters.molar_mass_kg_per_mol;
    const double solvent_mass = molar_masses[solvent];
    if (!std::isfinite(solvent_mass) || solvent_mass <= 0.0) {
        throw std::invalid_argument("neutral reference solvent has invalid molar mass");
    }
    double charge = 0.0;
    double reaction_charge = 0.0;
    for (std::size_t component = 0; component < components; ++component) {
        const double ratio = declaration.solute_molality_ratios[component];
        if (!std::isfinite(ratio)
            || (component == solvent ? ratio != 0.0 : ratio <= 0.0)) {
            throw std::invalid_argument("invalid neutral reference solute ratio");
        }
        const double coefficient = stoichiometry[component];
        if (!std::isfinite(coefficient)) {
            throw std::invalid_argument("invalid neutral reference stoichiometry");
        }
        charge += ratio * static_cast<double>(parameters.charge_number[component]);
        reaction_charge += coefficient * static_cast<double>(parameters.charge_number[component]);
    }
    if (!std::isfinite(charge) || charge != 0.0) {
        throw std::invalid_argument("neutral reference solute ratios are not charge neutral");
    }
    if (!std::isfinite(reaction_charge) || reaction_charge != 0.0) {
        throw std::invalid_argument("neutral reference reaction is not charge neutral");
    }
    return solvent_mass;
}

std::vector<double> neutral_reference_composition(
    const NeutralReferenceDeclaration& declaration,
    double solvent_mass,
    double molality
) {
    if (!std::isfinite(molality) || molality <= 0.0) {
        throw std::invalid_argument("invalid neutral reference molality");
    }
    const std::size_t solvent = *declaration.solvent_index;
    const double solvent_moles = 1.0 / solvent_mass;
    double total_moles = solvent_moles;
    for (const double ratio : declaration.solute_molality_ratios)
        total_moles += molality * ratio;
    if (!std::isfinite(total_moles) || total_moles <= 0.0) {
        throw std::invalid_argument("invalid neutral reference composition total");
    }
    std::vector<double> composition(declaration.solute_molality_ratios.size(), 0.0);
    composition[solvent] = solvent_moles / total_moles;
    for (std::size_t component = 0; component < composition.size(); ++component) {
        if (component == solvent) continue;
        composition[component] = molality * declaration.solute_molality_ratios[component]
            / total_moles;
    }
    return composition;
}

ResidualEnergyAction reference_partial_from_state(const HelmholtzModel& model, double temperature,
    const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    std::size_t temperature_order, std::size_t pressure_order,
    std::optional<std::size_t> parameter, const State& state) {
    const std::size_t order = temperature_order + pressure_order + std::size_t(parameter.has_value());
    ResidualEnergyAction result;
    const double pressure_scale = state.pressure_density_derivative * state.molar_density;
    const double theta_scale = parameter ? std::max(1.0, std::abs(model.active_parameter_values()[*parameter])) : 1.0;
    const double scale = std::pow(temperature, int(temperature_order)) * std::pow(pressure_scale, int(pressure_order)) * theta_scale;
    double contraction = 0.0;
    for (std::size_t mask = 0; mask < (std::size_t(1) << (order - 1)); ++mask) {
        PressureActionDirection direction;
        double weight = 1.0;
        for (std::size_t j = 0; j < order; ++j) {
            const double sign = j && (mask & (std::size_t(1) << (j-1))) ? -1.0 : 1.0;
            weight *= sign;
            if (j < temperature_order) direction.temperature += sign * temperature;
            else if (j < temperature_order + pressure_order) direction.pressure += sign * pressure_scale;
            else {
                const auto& active = model.active_parameters()[*parameter];
                direction.parameter.push_back({active.family, active.indices, sign * theta_scale});
            }
        }
        const auto series = model.reference_taylor(temperature, 1.0 / state.molar_density,
            composition, stoichiometry, direction, order);
        result.status = series.status;
        result.diagnostics.site_factorizations += series.diagnostics.site_factorizations;
        result.diagnostics.site_solves += series.diagnostics.site_solves;
        result.diagnostics.pressure_factorizations += series.diagnostics.pressure_factorizations;
        result.diagnostics.pressure_solves += series.diagnostics.pressure_solves;
        result.diagnostics.pressure_density_derivative = series.diagnostics.pressure_density_derivative;
        result.diagnostics.forward_sweeps += series.diagnostics.forward_sweeps;
        result.diagnostics.directional_paths += series.diagnostics.directional_paths;
        result.diagnostics.association_residual = std::max(result.diagnostics.association_residual, series.diagnostics.association_residual);
        result.diagnostics.pressure_residual = std::max(result.diagnostics.pressure_residual, series.diagnostics.pressure_residual);
        result.diagnostics.site_rank = series.diagnostics.site_rank;
        result.diagnostics.site_pivot_ratio = std::max(result.diagnostics.site_pivot_ratio, series.diagnostics.site_pivot_ratio);
        if (series.status != PropertyActionStatus::Available) return result;
        contraction += weight * series.contraction[order];
    }
    result.action = contraction / (double(std::size_t(1) << (order-1)) * scale);
    if (!std::isfinite(*result.action)) { result.action.reset(); result.status = PropertyActionStatus::Nonfinite; }
    return result;
}

struct NeutralReferenceTrace final {
    double pressure = 0.0;
    double solvent_mass = 0.0;
    std::array<std::vector<double>, NEUTRAL_REFERENCE_MOLALITIES.size()> composition;
    std::array<State, NEUTRAL_REFERENCE_MOLALITIES.size()> state;
    std::optional<HelmholtzDerivatives> terminal_derivatives;
    std::size_t temperature_order = 0;
    bool parameter_jets = false;
    double pressure_v = 0.0;
    double pressure_vv = 0.0;
};

struct NeutralReferenceTraceResult final {
    std::optional<NeutralReferenceTrace> trace;
    NeutralReferenceStatus status = NeutralReferenceStatus::Nonfinite;
    std::string message;
};

NeutralReferenceTraceResult build_neutral_reference_trace(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const NeutralReferenceDeclaration& declaration,
    bool pressure_jets,
    std::optional<double> solvent_anchor,
    std::size_t temperature_order,
    bool parameter_jets
) {
    core::Scope scope("neutral_reference_trace");
    NeutralReferenceTraceResult result;
    try {
        const std::size_t components = model.parameters().component_count();
        const double solvent_mass = validate_neutral_reference_declaration(
            model.parameters(), declaration, std::vector<double>(components, 0.0));
        if (!std::isfinite(temperature) || temperature <= 0.0
            || !std::isfinite(pressure) || pressure <= 0.0) {
            result.status = NeutralReferenceStatus::InvalidDeclaration;
            result.message = "neutral reference temperature or pressure is invalid";
            return result;
        }
        NeutralReferenceTrace trace;
        trace.temperature_order = temperature_order;
        trace.parameter_jets = parameter_jets;
        trace.pressure = pressure;
        trace.solvent_mass = solvent_mass;
        // A known liquid solvent density is tried first as a locally unique anchor; any failure
        // falls back to the cold scan and its branch admission.
        std::optional<double> anchor = solvent_anchor;
        for (std::size_t point = 0; point < NEUTRAL_REFERENCE_MOLALITIES.size(); ++point) {
            trace.composition[point] = neutral_reference_composition(
                declaration, solvent_mass, NEUTRAL_REFERENCE_MOLALITIES[point]);
            try {
                try {
                    trace.state[point] = model.state_at_pressure(
                        temperature, pressure, trace.composition[point], DensityBranch::Liquid, anchor);
                } catch (const std::exception& error) {
                    if (point != 0 || !anchor) throw;
                    if (scope.diagnostics_on()) scope.probe("anchor_cold_retry",
                        {*anchor, static_cast<double>(point)}, 0, error.what());
                    trace.state[point] = model.state_at_pressure(
                        temperature, pressure, trace.composition[point], DensityBranch::Liquid, std::nullopt);
                }
            } catch (const std::exception& error) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = std::string("neutral reference density unavailable: ")
                    + error.what();
                return result;
            }
            const bool cold_unique_root = point == 0
                && trace.state[point].density_branch == DensityBranch::Unspecified
                && trace.state[point].stable_root_count.value_or(0) == 1;
            if (trace.state[point].density_branch != DensityBranch::Liquid && !cold_unique_root) {
                result.status = NeutralReferenceStatus::BranchChange;
                result.message = "neutral reference left the anchored liquid branch";
                return result;
            }
            if (!std::isfinite(trace.state[point].molar_density)
                || trace.state[point].molar_density <= 0.0
                || !std::isfinite(trace.state[point].pressure_density_derivative)
                || trace.state[point].pressure_density_derivative <= 0.0
                || trace.state[point].residual_chemical_potential_over_rt.size() != components) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference state has invalid density response";
                return result;
            }
            anchor = trace.state[point].molar_density;
        }
        if (pressure_jets) {
            const double volume = 1.0 / trace.state.back().molar_density;
            const auto derivatives = model.evaluate(
                temperature, volume, trace.composition.back(), DerivativeOrder::Hessian);
            const std::size_t width = derivatives.gradient.size();
            if (width != 2 + components + model.active_parameters().size()
                || derivatives.hessian.size() != width * width) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference pressure jet has invalid Hessian shape";
                return result;
            }
            std::vector<double> volume_direction(width, 0.0);
            volume_direction[1] = 1.0;
            const auto pressure_third = model.residual_energy_action(
                temperature, volume, trace.composition.back(),
                {volume_direction, volume_direction, volume_direction});
            if (!pressure_third.action.has_value()) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference pressure jet unavailable from residual E3";
                return result;
            }
            const double phi_vv = derivatives.hessian[width + 1];
            const double pressure_v = GAS_CONSTANT_J_PER_MOL_K * temperature
                * (-1.0 / (volume * volume) - phi_vv);
            const double pressure_vv = GAS_CONSTANT_J_PER_MOL_K * temperature
                * (2.0 / (volume * volume * volume) - *pressure_third.action);
            if (!std::isfinite(pressure_v) || pressure_v >= 0.0
                || !std::isfinite(pressure_vv)) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference pressure response is nonfinite";
                return result;
            }
            trace.terminal_derivatives = derivatives;
            trace.pressure_v = pressure_v;
            trace.pressure_vv = pressure_vv;
        }
        result.trace = std::move(trace);
        result.status = NeutralReferenceStatus::Available;
        return result;
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        result.status = message.find("charged") != std::string::npos
                || message.find("charge neutral") != std::string::npos
            ? NeutralReferenceStatus::NonNeutralDeclaration
            : NeutralReferenceStatus::InvalidDeclaration;
        result.message = message;
        return result;
    } catch (const std::exception& error) {
        result.status = NeutralReferenceStatus::Nonfinite;
        result.message = error.what();
        return result;
    }
}

bool supported_reference_basis(ReferenceSourceBasis basis) {
    return basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution
        || basis == ReferenceSourceBasis::RawMoleFractionInfiniteDilution;
}

NeutralReferenceReactionDiagnostic neutral_reference_from_trace(
    const HelmholtzModel& model,
    double temperature,
    const std::vector<double>& stoichiometry,
    const ReactionReference& reference,
    bool pressure_jets,
    const NeutralReferenceTrace& trace
) {
    // EqID: reference_single_insertion_taylor
    // Raw lower T/P/EOS partials reuse each resolved terminal state; admission stays below.
    NeutralReferenceReactionDiagnostic result;
    result.reference_pressure_pa = trace.pressure;
    try {
        std::array<double, NEUTRAL_REFERENCE_MOLALITIES.size()> contractions{};
        for (std::size_t point = 0; point < contractions.size(); ++point) {
            const auto& state = trace.state[point];
            const double log_density = std::log(state.molar_density / 1.0);
            double contraction = 0.0;
            for (std::size_t component = 0; component < stoichiometry.size(); ++component) {
                if (stoichiometry[component] == 0.0) continue;
                const double chemical = require_component_property(state.residual_chemical_potential_over_rt[component]);
                if (!std::isfinite(chemical)) {
                    result.status = NeutralReferenceStatus::Nonfinite;
                    result.message = "neutral reference residual chemical potential is nonfinite";
                    return result;
                }
                contraction += stoichiometry[component] * (log_density + chemical);
            }
            if (!std::isfinite(contraction)) {
                result.status = NeutralReferenceStatus::Nonfinite;
                result.message = "neutral reference contraction is nonfinite";
                return result;
            }
            contractions[point] = contraction;
        }
        result.observed_terminal_change = std::max(
            std::abs(contractions[4] - contractions[3]),
            std::abs(contractions[3] - contractions[2]));
        if (!std::isfinite(*result.observed_terminal_change)) {
            result.status = NeutralReferenceStatus::Nonfinite;
            result.message = "neutral reference terminal change is nonfinite";
            return result;
        }
        if (*result.observed_terminal_change > policy::REFERENCE_CONVERGENCE_TOLERANCE.value) {
            result.status = NeutralReferenceStatus::TerminalChangeExceeded;
            result.message = "neutral reference terminal change exceeds tolerance";
            return result;
        }
        result.terminal_contraction = contractions.back();
        const auto admit = [&](std::size_t a, std::size_t b, std::optional<std::size_t> parameter,
                               double scale, double tolerance, std::optional<double>* observed = nullptr) -> std::optional<double> {
            auto& column = result.lower_derivatives.emplace_back();
            column.temperature_order = a; column.pressure_order = b; column.parameter = parameter; column.criterion = tolerance;
            for (std::size_t point = 0; point < 3; ++point) {
                const auto rung = NEUTRAL_REFERENCE_MOLALITIES.size() - 3 + point;
                try {
                    const auto partial = reference_partial_from_state(model, temperature, trace.composition[rung],
                        stoichiometry, a, b, parameter, trace.state[rung]);
                    column.ladder[point] = partial.action;
                    if (!partial.action) {
                        column.action_status = partial.status;
                        column.message = "reference action status " + std::to_string(int(partial.status));
                    }
                } catch (const std::exception& error) { column.message = error.what(); }
            }
            const auto& ladder = column.ladder;
            if (!ladder[0] || !ladder[1] || !ladder[2]) return std::nullopt;
            if (!std::isfinite(*ladder[0] + *ladder[1] + *ladder[2])) {
                column.message = "nonfinite lower reference ladder"; return std::nullopt;
            }
            for (std::size_t k = 0; k < 2; ++k) column.scaled_changes[k] = scale * std::abs(*ladder[k+1] - *ladder[k]);
            const double change = std::max(column.scaled_changes[0], column.scaled_changes[1]);
            if (observed) *observed = change;
            if (change <= tolerance) column.value = ladder[2];
            else column.message = "scaled terminal change " + std::to_string(change) + " exceeds " + std::to_string(tolerance);
            return column.value;
        };
        const double tolerance = policy::REFERENCE_CONVERGENCE_TOLERANCE.value;
        if (trace.temperature_order > 0 || trace.parameter_jets) {
            result.terminal_contraction_temperature_derivative = admit(1, 0, std::nullopt, temperature, tolerance,
                &result.observed_temperature_derivative_change);
            if (trace.temperature_order > 1)
                result.terminal_contraction_temperature_second_derivative = admit(2, 0, std::nullopt, temperature * temperature,
                    policy::REFERENCE_HEAT_CAPACITY_TOLERANCE.value / GAS_CONSTANT_J_PER_MOL_K,
                    &result.observed_temperature_second_derivative_change);
            if (pressure_jets)
                result.terminal_contraction_temperature_pressure_derivative = admit(1, 1, std::nullopt, temperature * trace.pressure, tolerance);
        }
        const auto values = model.active_parameter_values();
        for (std::size_t a = 0; trace.parameter_jets && a < values.size(); ++a) {
            const double scale = std::max(std::abs(values[a]), 1.0);
            result.terminal_contraction_parameter_derivative.push_back(admit(0, 0, a, scale, tolerance));
            if (trace.temperature_order > 0) {
                result.observed_temperature_parameter_change.emplace_back();
                result.terminal_contraction_temperature_parameter_derivative.push_back(admit(1, 0, a, temperature * scale, tolerance,
                    &result.observed_temperature_parameter_change.back()));
                result.temperature_parameter_ladder.push_back(result.lower_derivatives.back().ladder);
            }
        }
        if (reference.reference_pressure_pa.has_value()) {
            result.terminal_contraction_pressure_derivative = 0.0;
            result.terminal_contraction_pressure_second_derivative = 0.0;
            result.terminal_contraction_temperature_pressure_derivative = 0.0;
        } else if (!pressure_jets) {
            // The reference value and terminal-change gate are available from the
            // ladder alone.  Defer E2/E3 pressure contractions until a pressure
            // direction is actually requested by the equilibrium caller.
        } else {
            if (!trace.terminal_derivatives.has_value()) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference pressure jet is missing";
                return result;
            }
            const auto& derivatives = *trace.terminal_derivatives;
            const std::size_t width = derivatives.gradient.size();
            std::vector<double> volume_direction(width, 0.0);
            std::vector<double> stoichiometry_direction(width, 0.0);
            volume_direction[1] = 1.0;
            for (std::size_t component = 0; component < stoichiometry.size(); ++component)
                stoichiometry_direction[2 + component] = stoichiometry[component];
            const double volume = 1.0 / trace.state.back().molar_density;
            const auto contraction_third = model.residual_energy_action(
                temperature, volume, trace.composition.back(),
                {stoichiometry_direction, volume_direction, volume_direction});
            if (!contraction_third.action.has_value()) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference pressure jet unavailable from residual E3";
                return result;
            }
            double delta_nu = 0.0;
            double chemical_volume = 0.0;
            for (std::size_t component = 0; component < stoichiometry.size(); ++component) {
                delta_nu += stoichiometry[component];
                chemical_volume += stoichiometry[component]
                    * derivatives.hessian[(2 + component) * width + 1];
            }
            const double contraction_v = -delta_nu / volume + chemical_volume;
            const double contraction_vv = delta_nu / (volume * volume)
                + *contraction_third.action;
            const double volume_pressure = 1.0 / trace.pressure_v;
            const double volume_pressure_pressure = -trace.pressure_vv
                / (trace.pressure_v * trace.pressure_v * trace.pressure_v);
            const double contraction_pressure = contraction_v * volume_pressure;
            const double contraction_pressure_pressure = contraction_vv
                    * volume_pressure * volume_pressure
                + contraction_v * volume_pressure_pressure;
            if (!std::isfinite(contraction_pressure)
                || !std::isfinite(contraction_pressure_pressure)) {
                result.status = NeutralReferenceStatus::DensityUnavailable;
                result.message = "neutral reference pressure jet is nonfinite";
                return result;
            }
            result.terminal_contraction_pressure_derivative = contraction_pressure;
            result.terminal_contraction_pressure_second_derivative =
                contraction_pressure_pressure;
        }
        result.status = NeutralReferenceStatus::Available;
        return result;
    } catch (const std::invalid_argument& error) {
        const std::string message = error.what();
        result.status = message.find("charged") != std::string::npos
                || message.find("charge neutral") != std::string::npos
            ? NeutralReferenceStatus::NonNeutralDeclaration
            : NeutralReferenceStatus::InvalidDeclaration;
        result.message = message;
        return result;
    } catch (const std::exception& error) {
        result.status = NeutralReferenceStatus::Nonfinite;
        result.message = error.what();
        return result;
    }
}

}  // namespace

std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const Parameters& parameters,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality
) {
    // EqID: formula_unit_composition
    const auto& molar_masses = parameters.molar_mass_kg_per_mol;
    const auto& charges = parameters.charge_number;
    if (solvent_fractions.size() != molar_masses.size()
        || solvent_fractions.size() != charges.size()
        || ion_indices[0] >= solvent_fractions.size()
        || ion_indices[1] >= solvent_fractions.size()
        || ion_indices[0] == ion_indices[1]
        || stoichiometry[0] <= 0 || stoichiometry[1] <= 0
        || !std::isfinite(molality) || molality <= 0.0) {
        throw std::invalid_argument("invalid reference composition inputs");
    }
    if (charges[ion_indices[0]] == 0 || charges[ion_indices[1]] == 0
        || static_cast<std::int64_t>(stoichiometry[0]) * charges[ion_indices[0]]
            + static_cast<std::int64_t>(stoichiometry[1]) * charges[ion_indices[1]] != 0) {
        throw std::invalid_argument("reference formula must be charge neutral with two ions");
    }
    double fraction_sum = 0.0;
    double solvent_mass = 0.0;
    double solvent_charge = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        if (!std::isfinite(solvent_fractions[component])
            || solvent_fractions[component] < 0.0
            || !std::isfinite(molar_masses[component])
            || (solvent_fractions[component] > 0.0 && molar_masses[component] <= 0.0)) {
            throw std::invalid_argument("invalid reference solvent fractions");
        }
        fraction_sum += solvent_fractions[component];
        solvent_mass += solvent_fractions[component] * molar_masses[component];
        solvent_charge += solvent_fractions[component] * charges[component];
    }
    if (std::abs(fraction_sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value
        || solvent_fractions[ion_indices[0]] != 0.0
        || solvent_fractions[ion_indices[1]] != 0.0
        || !(solvent_mass > 0.0) || !std::isfinite(solvent_mass)
        || std::abs(solvent_charge) > policy::COMPOSITION_SUM_TOLERANCE.value) {
        throw std::invalid_argument("invalid reference solvent composition or molar mass");
    }
    const double solvent_moles = 1.0 / solvent_mass;
    const double total_moles = solvent_moles
        + molality * (static_cast<double>(stoichiometry[0]) + stoichiometry[1]);
    std::vector<double> composition(solvent_fractions.size(), 0.0);
    for (std::size_t component = 0; component < composition.size(); ++component) {
        composition[component] = solvent_fractions[component] * solvent_moles / total_moles;
    }
    composition[ion_indices[0]] = stoichiometry[0] * molality / total_moles;
    composition[ion_indices[1]] = stoichiometry[1] * molality / total_moles;
    return composition;
}

ResidualEnergyAction reference_partial(const HelmholtzModel& model, double temperature,
    double pressure, const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    std::size_t temperature_order, std::size_t pressure_order,
    std::optional<std::size_t> parameter, std::optional<double> anchor) {
    const std::size_t order = temperature_order + pressure_order + std::size_t(parameter.has_value());
    if (order == 0 || order > 4 || (parameter && *parameter >= model.active_parameters().size()))
        throw std::invalid_argument("reference partial order or coordinate");
    const auto state = model.state_at_pressure(temperature, pressure, composition, DensityBranch::Liquid, anchor);
    return reference_partial_from_state(model, temperature, composition, stoichiometry,
        temperature_order, pressure_order, parameter, state);
}

double reaction_reference_log_scale(
    const HelmholtzModel& model,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry
) {
    const double solvent_mass = validate_neutral_reference_declaration(
        model.parameters(), declaration, stoichiometry);
    double solute_stoichiometry = 0.0;
    const std::size_t solvent = *declaration.solvent_index;
    for (std::size_t component = 0; component < stoichiometry.size(); ++component)
        if (component != solvent) solute_stoichiometry += stoichiometry[component];
    const double scale = solute_stoichiometry
        * std::log(solvent_mass * declaration.standard_molality_mol_per_kg);
    if (!std::isfinite(scale)) throw std::invalid_argument("nonfinite neutral reference scale");
    return scale;
}

std::vector<NeutralReferenceReactionDiagnostic> reaction_references(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<std::vector<double>>& stoichiometries,
    const std::vector<ReactionReference>& references,
    bool pressure_jets,
    std::optional<double> solvent_anchor,
    std::size_t temperature_order,
    bool parameter_jets,
    std::size_t reference_order
) {
    core::Scope scope("reaction_references");
    if (reference_order > 4) throw std::invalid_argument("reference derivative order");
    if (scope.on()) {
        scope.probe("reference_inputs", {temperature, system_pressure, double(pressure_jets),
            double(temperature_order), double(parameter_jets), double(solvent_anchor.has_value()), solvent_anchor.value_or(0.0)});
        scope.probe("reference_declaration", {double(declaration.solvent_index.has_value()),
            double(declaration.solvent_index.value_or(0)), declaration.standard_molality_mol_per_kg});
        scope.probe("reference_ratios", declaration.solute_molality_ratios);
        scope.observe([&] { scope.probe("active_parameter_values", model.active_parameter_values()); });
    }
    if (stoichiometries.size() != references.size())
        throw std::invalid_argument("neutral reference batch size");
    struct CachedTrace final {
        double pressure = 0.0;
        bool pressure_jets = false;
        NeutralReferenceTraceResult result;
    };
    std::vector<CachedTrace> cache;
    std::vector<NeutralReferenceReactionDiagnostic> result(stoichiometries.size());
    if (!std::isfinite(temperature) || temperature <= 0.0
        || !std::isfinite(system_pressure) || system_pressure <= 0.0) {
        for (auto& diagnostic : result) {
            diagnostic.status = NeutralReferenceStatus::InvalidDeclaration;
            diagnostic.message = "neutral reference temperature or system pressure is invalid";
        }
        return result;
    }
    for (std::size_t row = 0; row < stoichiometries.size(); ++row) {
        auto& diagnostic = result[row];
        const auto& reference = references[row];
        if (!supported_reference_basis(reference.source_basis)) {
            diagnostic.status = NeutralReferenceStatus::UnsupportedSourceBasis;
            diagnostic.message = "unsupported neutral reference source basis";
            continue;
        }
        try {
            validate_neutral_reference_declaration(
                model.parameters(), declaration, stoichiometries[row]);
        } catch (const std::invalid_argument& error) {
            const std::string message = error.what();
            diagnostic.status = message.find("charged") != std::string::npos
                    || message.find("charge neutral") != std::string::npos
                ? NeutralReferenceStatus::NonNeutralDeclaration
                : NeutralReferenceStatus::InvalidDeclaration;
            diagnostic.message = message;
            continue;
        }
        const double pressure = reference.reference_pressure_pa.value_or(system_pressure);
        const bool need_jets = pressure_jets && !reference.reference_pressure_pa.has_value();
        auto found = std::find_if(cache.begin(), cache.end(), [&](const auto& item) {
            return item.pressure == pressure && (!need_jets || item.pressure_jets);
        });
        if (scope.on()) scope.probe("reference_row", {double(row), pressure,
            double(reference.reference_pressure_pa.has_value()), double(reference.source_basis), double(need_jets), double(found != cache.end())});
        if (found == cache.end()) {
            CachedTrace item;
            item.pressure = pressure;
            item.pressure_jets = need_jets;
            item.result = build_neutral_reference_trace(model, temperature, pressure,
                declaration, need_jets, solvent_anchor, temperature_order, parameter_jets);
            cache.push_back(std::move(item));
            found = std::prev(cache.end());
        }
        if (!found->result.trace.has_value()) {
            diagnostic.reference_pressure_pa = pressure;
            diagnostic.status = found->result.status;
            diagnostic.message = found->result.message;
            continue;
        }
        diagnostic = neutral_reference_from_trace(
            model, temperature, stoichiometries[row], reference,
            need_jets, *found->result.trace);
        diagnostic.solvent_density = found->result.trace->state.front().molar_density;
        // Additional consumed columns, with independent admission on both terminal changes.
        const auto& trace = *found->result.trace;
        const auto values = model.active_parameter_values();
        for (std::size_t c = 0; c <= std::size_t(parameter_jets); ++c)
            for (std::size_t d = 2; d <= reference_order; ++d)
                for (std::size_t b = 0; b <= d-c; ++b) {
                    const std::size_t a = d-b-c;
                    if (d == 2 && !(c && b)) continue;  // existing lower slopes keep their admitted owner
                    for (std::size_t coordinate = 0; coordinate < (c ? values.size() : 1); ++coordinate) {
                        NeutralReferenceReactionDiagnostic::Derivative column;
                        column.temperature_order = a; column.pressure_order = b;
                        if (c) column.parameter = coordinate;
                        if (b && reference.reference_pressure_pa) {
                            column.ladder = {0.0, 0.0, 0.0}; column.value = 0.0;
                        } else {
                            for (std::size_t k = 0; k < 3; ++k) {
                                const std::size_t rung = NEUTRAL_REFERENCE_MOLALITIES.size()-3+k;
                                try {
                                    const auto partial = reference_partial(model, temperature, pressure,
                                        trace.composition[rung], stoichiometries[row], a, b, column.parameter,
                                        trace.state[rung].molar_density);
                                    column.ladder[k] = partial.action;
                                    if (!partial.action) column.message = "reference action status " + std::to_string(int(partial.status));
                                } catch (const std::exception& error) { column.message = error.what(); }
                            }
                            const double scale = GAS_CONSTANT_J_PER_MOL_K * std::pow(temperature, int(a))
                                * std::pow(pressure, int(b)) * (c ? std::max(std::abs(values[coordinate]), 1.0) : 1.0);
                            if (column.ladder[0] && column.ladder[1] && column.ladder[2]) {
                                bool admitted = true;
                                for (std::size_t k = 0; k < 2; ++k) {
                                    column.scaled_changes[k] = scale * std::abs(*column.ladder[k+1]-*column.ladder[k]);
                                    admitted &= std::isfinite(column.scaled_changes[k])
                                        && column.scaled_changes[k] <= policy::REFERENCE_HEAT_CAPACITY_TOLERANCE.value;
                                }
                                if (admitted) column.value = column.ladder[2];
                                else column.message = "reference scaled terminal change exceeds heat-capacity tolerance";
                            }
                        }
                        if (!column.value) column.message = "reaction " + std::to_string(row) + " C_"
                            + std::to_string(a) + std::to_string(b) + (c ? " coordinate " + std::to_string(coordinate) : "")
                            + ": " + column.message;
                        diagnostic.derivatives.push_back(std::move(column));
                    }
                }

    }
    return result;
}

NeutralReferenceReactionDiagnostic reaction_reference(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry,
    const ReactionReference& reference
) {
    return reaction_references(
        model, temperature, system_pressure, declaration,
        {stoichiometry}, {reference}).front();
}

InfiniteDilutionReference infinite_dilution_reference(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry
) {
    core::Scope scope("infinite_dilution_reference");
    std::array<double, 4> molalities = {
        policy::REFERENCE_MOLALITY_ONE.value,
        policy::REFERENCE_MOLALITY_TWO.value,
        policy::REFERENCE_MOLALITY_THREE.value,
        policy::REFERENCE_MOLALITY_FOUR.value,
    };
    std::optional<double> anchor;
    std::array<double, 2> previous{};
    std::array<double, 2> current{};
    double convergence_error = 0.0;
    double terminal_error = 0.0;
    const std::size_t columns = model.active_parameters().size();
    std::vector<std::array<double, 2>> previous_columns(columns), current_columns(columns);
    std::vector<double> column_error(columns, 0.0);
    const auto values = model.active_parameter_values();
    const auto read_columns = [&](const State& state, auto& out) {
        for (std::size_t a = 0; a < columns; ++a)
            for (std::size_t i = 0; i < 2; ++i)
                out[a][i] = state.log_fugacity_parameter_derivative.empty() ? NAN
                    : require_component_property(state.log_fugacity_parameter_derivative[ion_indices[i] * columns + a]);
    };
    bool observed = false;
    for (std::size_t index = 0; index < molalities.size(); ++index) {
        const double molality = molalities[index];
        const auto composition = formula_composition(
            solvent_fractions,
            model.parameters(),
            ion_indices,
            stoichiometry,
            molality
        );
        const auto state = model.state_at_pressure(
            temperature, pressure, composition, DensityBranch::Liquid, anchor
        );
        current = ion_log_fugacity(state, ion_indices);
        scope.probe(
            "reference_ladder_point", {molality, current[0], current[1]},
            static_cast<std::uint32_t>(index)
        );
        if (observed) {
            // EqID: reference_limit_terminal_pair
            convergence_error = std::max(
                std::abs(current[0] - previous[0]),
                std::abs(current[1] - previous[1])
            );
            terminal_error = convergence_error;
        }
        read_columns(state, current_columns);
        if (observed) for (std::size_t a = 0; a < columns; ++a)
            column_error[a] = std::max(1.0, std::abs(values[a])) * std::max(
                std::abs(current_columns[a][0] - previous_columns[a][0]),
                std::abs(current_columns[a][1] - previous_columns[a][1]));
        previous_columns = current_columns;
        previous = current;
        observed = true;
        anchor = state.molar_density;
        scope.probe("reference_convergence_error", convergence_error,
                    static_cast<std::uint32_t>(index));
    }
    if (convergence_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value
        || std::any_of(column_error.begin(), column_error.end(), [](double e) {
            return !std::isfinite(e) || e > policy::REFERENCE_CONVERGENCE_TOLERANCE.value; })) {
        const auto composition = formula_composition(
            solvent_fractions,
            model.parameters(),
            ion_indices,
            stoichiometry,
            policy::REFERENCE_CONVERGENCE_WITNESS.value
        );
        const auto state = model.state_at_pressure(
            temperature, pressure, composition, DensityBranch::Liquid, anchor
        );
        const auto witness_log_fugacity = ion_log_fugacity(state, ion_indices);
        const double witness_error = std::max(
            std::abs(witness_log_fugacity[0] - current[0]),
            std::abs(witness_log_fugacity[1] - current[1])
        );
        scope.probe("reference_witness_error", witness_error);
        if (convergence_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value
            && witness_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value) {
            throw std::runtime_error("reference limit did not converge");
        }
        std::vector<std::array<double, 2>> witness(columns);
        read_columns(state, witness);
        for (std::size_t a = 0; a < columns; ++a)
            if (!std::isfinite(column_error[a]) || column_error[a] > policy::REFERENCE_CONVERGENCE_TOLERANCE.value)
                column_error[a] = std::max(1.0, std::abs(values[a])) * std::max(
                    std::abs(current_columns[a][0] - witness[a][0]), std::abs(current_columns[a][1] - witness[a][1]));
    }
    std::optional<double> pure_solvent;
    std::optional<State> pure_state;
    if (std::count_if(solvent_fractions.begin(), solvent_fractions.end(),
                     [](double x) { return x > 0.0; }) == 1) {
        const auto solvent = std::distance(solvent_fractions.begin(),
            std::max_element(solvent_fractions.begin(), solvent_fractions.end()));
        std::vector<double> pure(solvent_fractions.size(), 0.0);
        pure[solvent] = 1.0;
        pure_state = model.state_at_pressure(
            temperature, pressure, pure, DensityBranch::Liquid, std::nullopt);
        if (pure_state->log_fugacity_coefficient.has_value()) {
            pure_solvent = (*pure_state->log_fugacity_coefficient)[solvent];
        }
    }
    InfiniteDilutionReference result{current, molalities.back(), terminal_error, pure_solvent, std::move(pure_state), {}, column_error};
    result.parameter_derivative.resize(columns);
    for (std::size_t a = 0; a < columns; ++a)
        if (std::isfinite(column_error[a]) && column_error[a] <= policy::REFERENCE_CONVERGENCE_TOLERANCE.value)
            result.parameter_derivative[a] = {current_columns[a][0], current_columns[a][1]};
    return result;
}

ElectrolyteActivity electrolyte_activity(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality,
    const InfiniteDilutionReference& reference
) {
    core::Scope scope("electrolyte_activity");
    const auto composition = formula_composition(
        solvent_fractions,
        model.parameters(),
        ion_indices,
        stoichiometry,
        molality
    );
    const auto state = model.state_at_pressure(
        temperature, pressure, composition, DensityBranch::Liquid, std::nullopt
    );
    const auto log_fugacity = ion_log_fugacity(state, ion_indices);
    // EqID: single_ion_activity_coefficient
    const double first = log_fugacity[0] - reference.log_fugacity[0];
    const double second = log_fugacity[1] - reference.log_fugacity[1];
    const double coefficient_sum = static_cast<double>(stoichiometry[0]) + stoichiometry[1];
    // EqID: mean_ionic_activity_coefficient
    const double mean_log_x = (
        stoichiometry[0] * first + stoichiometry[1] * second
    ) / coefficient_sum;
    double solvent_mass = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        solvent_mass += solvent_fractions[component]
            * model.parameters().molar_mass_kg_per_mol[component];
    }
    // EqID: molality_activity_conversion
    const double mean_log_m = mean_log_x
        - std::log(1.0 + solvent_mass * molality * coefficient_sum);
    ElectrolyteActivity result;
    result.single_ion = {std::exp(first), std::exp(second)};
    result.mean_ionic_mole_fraction_scale = std::exp(mean_log_x);
    result.mean_ionic_log_molality_scale = mean_log_m;
    result.mean_ionic_molality_scale = std::exp(mean_log_m);
    if (reference.pure_solvent_state.has_value()) {
        if (molality < policy::OSMOTIC_MOLALITY_FLOOR.value) {
            result.osmotic_coefficient = core::Stop::converged_to_floor;
        } else {
            const auto solvent = std::distance(solvent_fractions.begin(),
                std::max_element(solvent_fractions.begin(), solvent_fractions.end()));
            const auto& pure = *reference.pure_solvent_state;
            const double s = coefficient_sum * molality * solvent_mass;
            // Extensive Helmholtz Euler identity, with the prescribed pressure in Z.
            long double chemical_difference = static_cast<long double>(state.residual_helmholtz) - pure.residual_helmholtz;
            chemical_difference += pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature)
                * ((static_cast<long double>(pure.molar_density) - state.molar_density) / (state.molar_density * static_cast<long double>(pure.molar_density)));
            for (const auto ion : ion_indices) chemical_difference -= composition[ion]
                * (static_cast<long double>(require_component_property(state.residual_chemical_potential_over_rt[ion])) - require_component_property(pure.residual_chemical_potential_over_rt[solvent]));
            chemical_difference /= composition[solvent];
            const double log_activity = -std::log1p(s) + static_cast<double>(chemical_difference)
                + std::log1p((state.molar_density - pure.molar_density) / pure.molar_density);
            // EqID: practical_osmotic_coefficient
            result.osmotic_coefficient = -log_activity / s;
        }
    }
    const std::size_t columns = model.active_parameters().size();
    for (std::size_t a = 0; a < columns; ++a) {
        // EqID: electrolyte_activity_parameter_derivative
        const auto& ref = reference.parameter_derivative.at(a);
        const auto g = [&](std::size_t i) { return require_component_property(state.log_fugacity_parameter_derivative.at(i * columns + a)); };
        result.mean_ionic_parameter_derivative.push_back(ref[0] && ref[1]
            ? result.mean_ionic_molality_scale * (stoichiometry[0] * (g(ion_indices[0]) - *ref[0])
                + stoichiometry[1] * (g(ion_indices[1]) - *ref[1])) / coefficient_sum : NAN);
        if (result.osmotic_coefficient && std::holds_alternative<double>(*result.osmotic_coefficient)) {
            const auto solvent = std::distance(solvent_fractions.begin(), std::max_element(solvent_fractions.begin(), solvent_fractions.end()));
            result.osmotic_parameter_derivative.push_back(-(g(solvent)
                - require_component_property(reference.pure_solvent_state->log_fugacity_parameter_derivative.at(solvent * columns + a)))
                / (coefficient_sum * molality * solvent_mass));
        }
    }
    return result;
}

}  // namespace epcsaft
