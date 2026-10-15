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

State properties::close_density(const HelmholtzModel& model, double temperature, double pressure,
                                const std::vector<double>& composition, DensityBranch branch, std::optional<double> prior_density) {
    if (prior_density) {
        core::Scope scope("density_anchor");
        try {
            auto state = model.state_at_pressure(temperature, pressure, composition, branch, prior_density);
            if (properties::regular_density(state.molar_density, state.pressure_density_derivative, state.density_convergence)) return state;
            if (state.density_convergence) scope.finish(*state.density_convergence);
            scope.probe("anchor_failure", {*prior_density}, 0, "target density regularity unavailable or exceeds DENSITY_ANCHOR_REGULARITY");
        } catch (const std::runtime_error& error) {
            if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error)) scope.finish(failure->record);
            scope.probe("anchor_failure", {*prior_density}, 0, error.what());
        }
    }
    return model.state_at_pressure(temperature, pressure, composition, branch);
}

namespace {

std::array<double, 2> ion_log_fugacity(const State& state, const std::array<std::size_t, 2>& ions) {
    if (!state.log_fugacity_coefficient.has_value())
        throw std::runtime_error("reference fugacity is unavailable");
    return {require_component_property((*state.log_fugacity_coefficient)[ions[0]]),
            require_component_property((*state.log_fugacity_coefficient)[ions[1]])};
}

// EqID: reference_molality_sequence; the fifth rung is the convergence witness.
constexpr std::array<double, 5> REFERENCE_MOLALITIES = {
    policy::REFERENCE_MOLALITY_ONE.value,
    policy::REFERENCE_MOLALITY_TWO.value,
    policy::REFERENCE_MOLALITY_THREE.value,
    policy::REFERENCE_MOLALITY_FOUR.value,
    policy::REFERENCE_CONVERGENCE_WITNESS.value,
};

// A declaration whose solvent, pure reference species, solute ratios or reaction carry net charge.
struct NonNeutralDeclaration final : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

NeutralReferenceStatus declaration_status(const std::invalid_argument& error) {
    return dynamic_cast<const NonNeutralDeclaration*>(&error) != nullptr
        ? NeutralReferenceStatus::NonNeutralDeclaration : NeutralReferenceStatus::InvalidDeclaration;
}

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
    if (!declaration.species_reference.empty()) {
        if (declaration.species_reference.size() != components)
            throw std::invalid_argument("species reference size");
        for (std::size_t i = 0; i < components; ++i) {
            const auto choice = declaration.species_reference[i];
            if (choice != SpeciesReference::InfiniteDilution && choice != SpeciesReference::PureLiquid
                && choice != SpeciesReference::PureGas)
                throw std::invalid_argument("species reference convention");
            if (choice != SpeciesReference::InfiniteDilution && parameters.charge_number[i] != 0)
                throw NonNeutralDeclaration("pure reference species is charged");
        }
    }
    const std::size_t solvent = *declaration.solvent_index;
    if (parameters.charge_number[solvent] != 0) {
        throw NonNeutralDeclaration("neutral reference solvent is charged");
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
        throw NonNeutralDeclaration("neutral reference solute ratios are not charge neutral");
    }
    if (!std::isfinite(reaction_charge) || reaction_charge != 0.0) {
        throw NonNeutralDeclaration("neutral reference reaction is not charge neutral");
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

// Reference-column scale of active coordinate a: max(1, |theta|); d_V (m3/mol) uses a 1 L/mol floor.
double parameter_scale(const HelmholtzModel& model, const std::vector<double>& values, std::size_t a) {
    return std::max(model.active_parameters()[a].family == ParameterFamily::DielectricVolume ? 1.0e-3 : 1.0, std::abs(values[a]));
}

}  // namespace

ResidualEnergyAction reference_partial_from_state(const HelmholtzModel& model, double temperature,
    const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    std::size_t temperature_order, std::size_t pressure_order,
    std::optional<std::size_t> parameter, const State& state) {
    const std::size_t order = temperature_order + pressure_order + std::size_t(parameter.has_value());
    ResidualEnergyAction result;
    const double pressure_scale = state.pressure_density_derivative * state.molar_density;
    const double theta_scale = parameter ? parameter_scale(model, model.active_parameter_values(), *parameter) : 1.0;
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

namespace {

struct NeutralReferenceTrace final {
    double pressure = 0.0;
    double solvent_mass = 0.0;
    std::size_t solvent_index = 0;
    std::array<std::vector<double>, REFERENCE_MOLALITIES.size()> composition;
    std::array<State, REFERENCE_MOLALITIES.size()> state;
    std::vector<std::optional<State>> pure_state;
    std::size_t temperature_order = 0;
    bool parameter_jets = false;
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
        trace.solvent_index = *declaration.solvent_index;
        // A known liquid solvent density is tried first as a locally unique anchor; any failure
        // falls back to the cold scan and its branch admission.
        std::optional<double> anchor = solvent_anchor;
        for (std::size_t point = 0; point < REFERENCE_MOLALITIES.size(); ++point) {
            trace.composition[point] = neutral_reference_composition(
                declaration, solvent_mass, REFERENCE_MOLALITIES[point]);
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
        std::vector<std::size_t> pure_indices;
        for (std::size_t i = 0; i < declaration.species_reference.size(); ++i) {
            const auto choice = declaration.species_reference[i];
            if (choice != SpeciesReference::InfiniteDilution
                && !(i == *declaration.solvent_index && choice == SpeciesReference::PureLiquid)) pure_indices.push_back(i);
        }
        if (!pure_indices.empty()) {
            trace.pure_state.resize(components);
            for (const auto i : pure_indices) {
                const auto choice = declaration.species_reference[i];
                std::vector<double> composition(components, 0.0);
                composition[i] = 1.0;
                try {
                    auto state = model.state_at_pressure(temperature, pressure, composition,
                        choice == SpeciesReference::PureGas ? DensityBranch::Vapor : DensityBranch::Liquid);
                    if (!(state.pressure_density_derivative > 0.0))
                        throw std::runtime_error("nonpositive pure-reference pressure slope");
                    if (state.density_branch == DensityBranch::Unspecified) {
                        // One stable root serves the declared branch unless it lies past the mechanical
                        // loop. The ALG-DENSITY scan grid finds the folds; a stable root is already sided
                        // against an unrefined fold sign bracket, so no refinement or budget is requested.
                        const double low = std::log(policy::DENSITY_LOWER_FRACTION.value * pressure
                            / (GAS_CONSTANT_J_PER_MOL_K * temperature)), span = std::log(policy::DENSITY_UPPER_FRACTION.value
                            * model.molar_density_upper_bound(temperature, composition)) - low;
                        const auto n = static_cast<std::size_t>(policy::DENSITY_BRACKET_INTERVALS.value);
                        std::vector<double> nodes(n + 1);
                        for (std::size_t k = 0; k <= n; ++k) nodes[k] = std::exp(low + static_cast<double>(k) / n * span);
                        constexpr double none = std::numeric_limits<double>::max();
                        const auto loop = model.density_diagnostics(temperature, pressure, composition, nodes, none, none,
                            std::numeric_limits<std::size_t>::max(), none, true, std::nullopt, none, DerivativeOrder::PressureSlope);
                        const bool gas = choice == SpeciesReference::PureGas;
                        if (!loop.unresolved.empty() || (!loop.spinodals.empty() && (gas
                            ? state.molar_density > loop.spinodals.front().upper
                            : state.molar_density < loop.spinodals.back().lower)))
                            throw std::runtime_error(std::string(gas ? "vapor" : "liquid")
                                + " branch absent: the single stable root lies across the mechanical loop");
                    }
                    trace.pure_state[i] = std::move(state);
                } catch (const std::exception& error) {
                    result.status = NeutralReferenceStatus::DensityUnavailable;
                    result.message = "pure reference species " + std::to_string(i) + ": " + error.what();
                    return result;
                }
            }
        }
        result.trace = std::move(trace);
        result.status = NeutralReferenceStatus::Available;
        return result;
    } catch (const std::invalid_argument& error) {
        result.status = declaration_status(error);
        result.message = error.what();
        return result;
    } catch (const std::exception& error) {
        result.status = NeutralReferenceStatus::Nonfinite;
        result.message = error.what();
        return result;
    }
}

ResidualEnergyAction reference_partial_from_trace(const HelmholtzModel& model, double temperature,
    const NeutralReferenceTrace& trace, std::size_t rung, const std::vector<double>& stoichiometry,
    std::size_t a, std::size_t b, std::optional<std::size_t> parameter) {
    if (trace.pure_state.empty())
        return reference_partial_from_state(model, temperature, trace.composition[rung],
            stoichiometry, a, b, parameter, trace.state[rung]);
    auto dilute = stoichiometry;
    for (std::size_t i = 0; i < dilute.size(); ++i)
        if (trace.pure_state[i]) dilute[i] = 0.0;
    ResidualEnergyAction result;
    result.status = PropertyActionStatus::Available;
    result.action = 0.0;
    if (std::any_of(dilute.begin(), dilute.end(), [](double v) { return v != 0.0; })) {
        result = reference_partial_from_state(model, temperature, trace.composition[rung],
            dilute, a, b, parameter, trace.state[rung]);
        if (!result.action) return result;
    }
    for (std::size_t i = 0; i < dilute.size(); ++i) {
        if (!trace.pure_state[i] || stoichiometry[i] == 0.0) continue;
        std::vector<double> composition(dilute.size(), 0.0), insertion(dilute.size(), 0.0);
        composition[i] = 1.0; insertion[i] = stoichiometry[i];
        const auto partial = reference_partial_from_state(model, temperature, composition,
            insertion, a, b, parameter, *trace.pure_state[i]);
        if (!partial.action) return partial;
        *result.action += *partial.action;
    }
    if (!std::isfinite(*result.action)) {
        result.action.reset(); result.status = PropertyActionStatus::Nonfinite;
    }
    return result;
}

// EqID: reference_limit_terminal_pair: the 1e-10, 1e-12, 1e-14 mol/kg rungs carry the
// Debye-Hueckel tail f0 + c sqrt(m) + O(m); each rung pair, a factor 100 in m, extrapolates
// sqrt(m) away as (10 f(m/100) - f(m)) / 9. The column returns the terminal pair's limit and
// is admitted on the scaled change between the two extrapolations, O(m) ~ 1e-11 relative to the
// tail. scaled_changes keep the raw rung changes. Returns the admitted change when the ladder is
// complete and finite.
std::optional<double> admit_ladder(NeutralReferenceReactionDiagnostic::Derivative& column, double scale, double tolerance) {
    column.criterion = tolerance;
    const auto& ladder = column.ladder;
    if (!ladder[0] || !ladder[1] || !ladder[2]) return std::nullopt;
    if (!std::isfinite(*ladder[0] + *ladder[1] + *ladder[2])) {
        column.message = "nonfinite reference ladder"; return std::nullopt;
    }
    for (std::size_t k = 0; k < 2; ++k) column.scaled_changes[k] = scale * std::abs(*ladder[k+1] - *ladder[k]);
    const auto limit = [&](std::size_t k) { return (10.0 * *ladder[k+1] - *ladder[k]) / 9.0; };
    const double change = scale * std::abs(limit(1) - limit(0));
    if (change <= tolerance) column.value = limit(1);
    else column.message = "scaled terminal change " + std::to_string(change) + " exceeds " + std::to_string(tolerance);
    return change;
}

// Scale is T^a P^b times parameter_scale. Columns that enter ln K and its first slopes admit at
// REFERENCE_CONVERGENCE_TOLERANCE; d2C/dT2 and the temperature-bearing higher columns (reference
// order >= 2) enter the reaction heat capacity and its derivatives, so they admit at
// REFERENCE_HEAT_CAPACITY_TOLERANCE / R. Pure-pressure higher columns are volume and compressibility
// quantities that share that tolerance by convention; their scaled changes are 1e-21 to 2e-17 in the
// reference-limits case, against 1.2e-3.
NeutralReferenceReactionDiagnostic::Derivative reference_column(const HelmholtzModel& model, double temperature,
    const NeutralReferenceTrace& trace, const std::vector<double>& stoichiometry, std::size_t a, std::size_t b,
    std::optional<std::size_t> parameter, double scale, double tolerance, std::optional<double>* observed = nullptr) {
    core::Scope scope("reference_column");
    NeutralReferenceReactionDiagnostic::Derivative column;
    column.temperature_order = a; column.pressure_order = b; column.parameter = parameter;
    for (std::size_t point = 0; point < 3; ++point) {
        try {
            const auto partial = reference_partial_from_trace(model, temperature, trace,
                REFERENCE_MOLALITIES.size() - 3 + point, stoichiometry, a, b, parameter);
            column.ladder[point] = partial.action;
            if (!partial.action) {
                column.action_status = partial.status;
                column.message = "reference action status " + std::to_string(int(partial.status));
            }
        } catch (const std::exception& error) { column.message = error.what(); }
    }
    const auto change = admit_ladder(column, scale, tolerance);
    if (observed && change) *observed = change;
    return column;
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
    result.solver_options = model.solver_options();
    result.reference_pressure_pa = trace.pressure;
    try {
        std::array<double, REFERENCE_MOLALITIES.size()> contractions{};
        for (std::size_t point = 0; point < contractions.size(); ++point) {
            const auto& state = trace.state[point];
            const double log_density = std::log(state.molar_density / 1.0);
            double contraction = 0.0;
            for (std::size_t component = 0; component < stoichiometry.size(); ++component) {
                if (stoichiometry[component] == 0.0) continue;
                const auto& source = !trace.pure_state.empty() && trace.pure_state[component]
                    ? *trace.pure_state[component] : state;
                if (const auto cause = &source == &state ? properties::vapor_side(state) : std::nullopt) {
                    result.status = NeutralReferenceStatus::DensityUnavailable;
                    result.message = "neutral reference solvent species " + std::to_string(trace.solvent_index) + ": " + *cause;
                    return result;
                }
                const double component_log_density = &source == &state ? log_density : std::log(source.molar_density / 1.0);
                const double chemical = require_component_property(source.residual_chemical_potential_over_rt[component]);
                if (!std::isfinite(chemical)) {
                    result.status = NeutralReferenceStatus::Nonfinite;
                    result.message = "neutral reference residual chemical potential is nonfinite";
                    return result;
                }
                contraction += stoichiometry[component] * (component_log_density + chemical);
            }
            if (!std::isfinite(contraction)) {
                result.status = NeutralReferenceStatus::Nonfinite;
                result.message = "neutral reference contraction is nonfinite";
                return result;
            }
            contractions[point] = contraction;
        }
        NeutralReferenceReactionDiagnostic::Derivative terminal;
        terminal.ladder = {contractions[2], contractions[3], contractions[4]};
        result.observed_terminal_change = admit_ladder(terminal, 1.0, policy::REFERENCE_CONVERGENCE_TOLERANCE.value);
        if (!terminal.value) {
            result.status = NeutralReferenceStatus::TerminalChangeExceeded;
            result.message = "neutral reference " + terminal.message;
            return result;
        }
        result.terminal_contraction = terminal.value;
        const auto admit = [&](std::size_t a, std::size_t b, std::optional<std::size_t> parameter,
                               double scale, double tolerance, std::optional<double>* observed = nullptr) {
            return result.lower_derivatives.emplace_back(reference_column(model, temperature, trace, stoichiometry,
                a, b, parameter, scale, tolerance, observed)).value;
        };
        const double tolerance = policy::REFERENCE_CONVERGENCE_TOLERANCE.value;
        if (trace.temperature_order > 0) {
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
            const double scale = parameter_scale(model, values, a);
            result.terminal_contraction_parameter_derivative.push_back(admit(0, 0, a, scale, tolerance));
            if (trace.temperature_order > 0)
                result.terminal_contraction_temperature_parameter_derivative.push_back(admit(1, 0, a, temperature * scale, tolerance));
        }
        if (reference.reference_pressure_pa.has_value()) {
            result.terminal_contraction_pressure_derivative = 0.0;
            result.terminal_contraction_pressure_second_derivative = 0.0;
            result.terminal_contraction_temperature_pressure_derivative = 0.0;
        } else if (!pressure_jets) {
            // The reference value and terminal-change gate are available from the
            // ladder alone.  Defer the pressure columns until a pressure
            // direction is actually requested by the equilibrium caller.
        } else {
            // The pressure slopes come from the same sqrt(m) limit as the value and its T slopes.
            result.terminal_contraction_pressure_derivative = admit(0, 1, std::nullopt, trace.pressure, tolerance);
            result.terminal_contraction_pressure_second_derivative = admit(0, 2, std::nullopt,
                trace.pressure * trace.pressure, policy::REFERENCE_HEAT_CAPACITY_TOLERANCE.value / GAS_CONSTANT_J_PER_MOL_K);
        }
        result.status = NeutralReferenceStatus::Available;
        return result;
    } catch (const std::invalid_argument& error) {
        result.status = declaration_status(error);
        result.message = error.what();
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
        throw std::invalid_argument("invalid reference composition inputs: need one solvent fraction per component, two "
            "distinct ions, positive stoichiometry and finite molality > 0 (got molality " + core::optional_number(molality) + ")");
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

std::optional<std::vector<NeutralReferenceReactionDiagnostic>> HelmholtzModel::stored_references(const std::vector<double>& inputs) const {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& [stored, result] : references_) if (stored == inputs) return result;
    return std::nullopt;
}

// 256 entries keep a fit evaluation's repeated start-pressure requests (largest reuse gap 130 stores,
// #449); a bubble solve's moving pressure only cycles the FIFO.
void HelmholtzModel::store_references(std::vector<double> inputs, const std::vector<NeutralReferenceReactionDiagnostic>& result) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (references_.size() == 256) references_.pop_front();
    references_.emplace_back(std::move(inputs), result);
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
    // Every input of this pure function of the model, segment sizes first so the flat key is unambiguous.
    std::vector<double> inputs{temperature, system_pressure, double(pressure_jets), double(solvent_anchor.has_value()),
        solvent_anchor.value_or(0.0), double(temperature_order), double(parameter_jets), double(reference_order),
        double(declaration.solvent_index.has_value()), double(declaration.solvent_index.value_or(0)),
        declaration.standard_molality_mol_per_kg, double(declaration.solute_molality_ratios.size()),
        double(declaration.species_reference.size()), double(stoichiometries.size())};
    inputs.insert(inputs.end(), declaration.solute_molality_ratios.begin(), declaration.solute_molality_ratios.end());
    for (const auto species : declaration.species_reference) inputs.push_back(double(species));
    for (std::size_t row = 0; row < stoichiometries.size(); ++row) {
        inputs.push_back(double(stoichiometries[row].size()));
        inputs.insert(inputs.end(), stoichiometries[row].begin(), stoichiometries[row].end());
        inputs.insert(inputs.end(), {double(references[row].reference_pressure_pa.has_value()),
            references[row].reference_pressure_pa.value_or(0.0), double(references[row].source_basis)});
    }
    if (auto stored = model.stored_references(inputs)) return std::move(*stored);
    scope.probe("reference_cache_miss", 1.0);
    struct CachedTrace final {
        double pressure = 0.0;
        NeutralReferenceTraceResult result;
    };
    std::vector<CachedTrace> cache;
    std::vector<NeutralReferenceReactionDiagnostic> result(stoichiometries.size());
    for (auto& diagnostic : result) diagnostic.solver_options = model.solver_options();
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
        try {
            validate_neutral_reference_declaration(
                model.parameters(), declaration, stoichiometries[row]);
        } catch (const std::invalid_argument& error) {
            diagnostic.status = declaration_status(error);
            diagnostic.message = error.what();
            continue;
        }
        const double pressure = reference.reference_pressure_pa.value_or(system_pressure);
        const bool need_jets = pressure_jets && !reference.reference_pressure_pa.has_value();
        auto found = std::find_if(cache.begin(), cache.end(), [&](const auto& item) {
            return item.pressure == pressure;
        });
        if (scope.on()) scope.probe("reference_row", {double(row), pressure,
            double(reference.reference_pressure_pa.has_value()), double(reference.source_basis), double(need_jets), double(found != cache.end())});
        if (found == cache.end()) {
            CachedTrace item;
            item.pressure = pressure;
            item.result = build_neutral_reference_trace(model, temperature, pressure,
                declaration, solvent_anchor, temperature_order, parameter_jets);
            cache.push_back(std::move(item));
            found = std::prev(cache.end());
        }
        if (!found->result.trace.has_value()) {
            diagnostic.reference_pressure_pa = pressure;
            diagnostic.status = found->result.status;
            diagnostic.message = found->result.message;
            continue;
        }
        if (!found->result.trace->pure_state.empty() && reference.source_basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution) {
            diagnostic.status = NeutralReferenceStatus::UnsupportedSourceBasis;
            diagnostic.message = "mixed species references require raw mole fractions; common molality is not qualified";
            continue;
        }
        diagnostic = neutral_reference_from_trace(
            model, temperature, stoichiometries[row], reference,
            need_jets, *found->result.trace);
        diagnostic.solvent_density = found->result.trace->state.front().molar_density;
        // Additional consumed columns, each the admitted sqrt(m) limit of its ladder.
        const auto& trace = *found->result.trace;
        const auto values = model.active_parameter_values();
        for (std::size_t c = 0; c <= std::size_t(parameter_jets); ++c)
            for (std::size_t d = 2; d <= reference_order; ++d)
                for (std::size_t b = 0; b <= d-c; ++b) {
                    const std::size_t a = d-b-c;
                    if (d == 2 && !(c && b)) continue;  // existing lower slopes keep their admitted owner
                    for (std::size_t coordinate = 0; coordinate < (c ? values.size() : 1); ++coordinate) {
                        const auto parameter = c ? std::optional<std::size_t>(coordinate) : std::nullopt;
                        NeutralReferenceReactionDiagnostic::Derivative column;
                        if (b && reference.reference_pressure_pa) {
                            column.temperature_order = a; column.pressure_order = b; column.parameter = parameter;
                            column.ladder = {0.0, 0.0, 0.0}; column.value = 0.0;
                        } else {
                            column = reference_column(model, temperature, trace, stoichiometries[row], a, b, parameter,
                                std::pow(temperature, int(a)) * std::pow(pressure, int(b))
                                    * (c ? parameter_scale(model, values, coordinate) : 1.0),
                                policy::REFERENCE_HEAT_CAPACITY_TOLERANCE.value / GAS_CONSTANT_J_PER_MOL_K);
                        }
                        if (!column.value) column.message = "reaction " + std::to_string(row) + " C_"
                            + std::to_string(a) + std::to_string(b) + (c ? " coordinate " + std::to_string(coordinate) : "")
                            + ": " + column.message;
                        diagnostic.derivatives.push_back(std::move(column));
                    }
                }

    }
    model.store_references(std::move(inputs), result);
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
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& solvent_fractions, const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    std::optional<double> first_rung_density, std::optional<double> pure_solvent_density
) {
    core::Scope scope("infinite_dilution_reference");
    InfiniteDilutionReference result;
    result.solver_options = model.solver_options();
    std::optional<double> anchor;
    std::array<double, 2> previous{};
    const std::size_t columns = model.active_parameters().size();
    std::vector<std::array<double, 2>> previous_columns(columns);
    result.parameter_convergence_error.assign(columns, 0.0);
    result.parameter_derivative.resize(columns);
    result.parameter_difference.resize(columns); result.witness_parameter_difference.resize(columns);
    const auto values = model.active_parameter_values();
    for (std::size_t index = 0; index < REFERENCE_MOLALITIES.size(); ++index) {
        const bool witness = index + 1 == REFERENCE_MOLALITIES.size();
        const auto failed_column = [&](double e) {
            return !std::isfinite(e) || e > policy::REFERENCE_CONVERGENCE_TOLERANCE.value;
        };
        if (witness && !(result.convergence_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value
            || std::any_of(result.parameter_convergence_error.begin(), result.parameter_convergence_error.end(), failed_column))) break;
        const auto composition = formula_composition(solvent_fractions, model.parameters(),
            ion_indices, stoichiometry, REFERENCE_MOLALITIES[index]);
        auto state = index == 0 ? properties::close_density(model, temperature, pressure, composition, DensityBranch::Liquid, first_rung_density)
            : model.state_at_pressure(temperature, pressure, composition, DensityBranch::Liquid, anchor);
        if (const auto cause = properties::vapor_side(state)) throw std::runtime_error(*cause);
        const auto current = ion_log_fugacity(state, ion_indices);
        // EqID: reference_limit_terminal_pair; the witness never replaces terminal values or columns.
        const double error = std::max(std::abs(current[0] - previous[0]), std::abs(current[1] - previous[1]));
        if (index) for (std::size_t i = 0; i < 2; ++i)
            (witness ? result.witness_difference : result.value_difference)[i] = current[i] - previous[i];
        if (witness) {
            scope.probe("reference_witness_error", error);
            if (result.convergence_error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value
                && error > policy::REFERENCE_CONVERGENCE_TOLERANCE.value)
                throw std::runtime_error("reference limit did not converge");
        } else {
            if (index) result.convergence_error = error;
            result.log_fugacity = current;
            result.molality = REFERENCE_MOLALITIES[index];
            scope.probe("reference_ladder_point", {REFERENCE_MOLALITIES[index], current[0], current[1]}, static_cast<std::uint32_t>(index));
        }
        for (std::size_t a = 0; a < columns; ++a) {
            std::array<double, 2> column;
            for (std::size_t i = 0; i < 2; ++i) column[i] = state.log_fugacity_parameter_derivative.empty() ? NAN
                : require_component_property(state.log_fugacity_parameter_derivative[ion_indices[i] * columns + a]);
            if (index) for (std::size_t i = 0; i < 2; ++i)
                (witness ? result.witness_parameter_difference : result.parameter_difference)[a][i] = column[i] - previous_columns[a][i];
            if (index && (!witness || failed_column(result.parameter_convergence_error[a])))
                result.parameter_convergence_error[a] = parameter_scale(model, values, a)
                    * std::max(std::abs(column[0] - previous_columns[a][0]), std::abs(column[1] - previous_columns[a][1]));
            if (!witness) previous_columns[a] = column;
        }
        result.ladder_states.push_back(state);
        if (!witness) {
            previous = current;
            anchor = state.molar_density;
            scope.probe("reference_convergence_error", result.convergence_error, static_cast<std::uint32_t>(index));
            if (index == 0) result.first_rung_state = std::move(state);
        }
    }
    for (std::size_t a = 0; a < columns; ++a)
        if (std::isfinite(result.parameter_convergence_error[a])
            && result.parameter_convergence_error[a] <= policy::REFERENCE_CONVERGENCE_TOLERANCE.value)
            result.parameter_derivative[a] = {previous_columns[a][0], previous_columns[a][1]};
    if (std::count_if(solvent_fractions.begin(), solvent_fractions.end(), [](double x) { return x > 0.0; }) == 1) {
        const auto solvent = std::distance(solvent_fractions.begin(), std::max_element(solvent_fractions.begin(), solvent_fractions.end()));
        std::vector<double> pure(solvent_fractions.size(), 0.0);
        pure[solvent] = 1.0;
        result.pure_solvent_state = properties::close_density(model, temperature, pressure, pure, DensityBranch::Liquid, pure_solvent_density);
        if (const auto cause = properties::vapor_side(*result.pure_solvent_state)) throw std::runtime_error(*cause);
        if (result.pure_solvent_state->log_fugacity_coefficient)
            result.pure_solvent_log_fugacity = (*result.pure_solvent_state->log_fugacity_coefficient)[solvent];
    }
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
    const InfiniteDilutionReference& reference, std::optional<double> finite_density
) {
    core::Scope scope("electrolyte_activity");
    if (!(reference.solver_options == model.solver_options()))
        throw std::invalid_argument("electrolyte reference solver options do not match the model");
    const auto composition = formula_composition(
        solvent_fractions,
        model.parameters(),
        ion_indices,
        stoichiometry,
        molality
    );
    ElectrolyteActivity result;
    result.solver_options = model.solver_options();
    result.state = properties::close_density(model, temperature, pressure, composition, DensityBranch::Liquid, finite_density);
    if (const auto cause = properties::vapor_side(result.state)) throw std::runtime_error(*cause);
    const auto& state = result.state;
    const auto log_fugacity = ion_log_fugacity(state, ion_indices);
    // EqID: single_ion_activity_coefficient
    const double first = log_fugacity[0] - reference.log_fugacity[0];
    const double second = log_fugacity[1] - reference.log_fugacity[1];
    const double coefficient_sum = static_cast<double>(stoichiometry[0]) + stoichiometry[1];
    // EqID: mean_ionic_activity_coefficient
    const double mean_log_x = (
        stoichiometry[0] * first + stoichiometry[1] * second
    ) / coefficient_sum;
    const auto solvent = std::distance(solvent_fractions.begin(),
        std::max_element(solvent_fractions.begin(), solvent_fractions.end()));
    double solvent_mass = 0.0;
    for (std::size_t component = 0; component < solvent_fractions.size(); ++component) {
        solvent_mass += solvent_fractions[component]
            * model.parameters().molar_mass_kg_per_mol[component];
    }
    // EqID: molality_activity_conversion
    const double mean_log_m = mean_log_x
        - std::log(1.0 + solvent_mass * molality * coefficient_sum);
    result.single_ion = {std::exp(first), std::exp(second)};
    result.mean_ionic_mole_fraction_scale = std::exp(mean_log_x);
    result.mean_ionic_log_molality_scale = mean_log_m;
    result.mean_ionic_molality_scale = std::exp(mean_log_m);
    if (reference.pure_solvent_state.has_value()) {
        if (molality < policy::OSMOTIC_MOLALITY_FLOOR.value) {
            result.osmotic_coefficient = core::Stop::converged_to_floor;
        } else {
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
            result.osmotic_parameter_derivative.push_back(-(g(solvent)
                - require_component_property(reference.pure_solvent_state->log_fugacity_parameter_derivative.at(solvent * columns + a)))
                / (coefficient_sum * molality * solvent_mass));
        }
    }
    return result;
}

ElectrolyteInputAction electrolyte_input_action(const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry,
    double molality, const InfiniteDilutionReference& reference, const ElectrolyteActivity& activity,
    const ElectrolyteInputDirection& direction, double standard_deviation) {
    ElectrolyteInputAction result;
    if (!(reference.solver_options == model.solver_options()) || !(activity.solver_options == model.solver_options()))
        throw std::invalid_argument("electrolyte input solver options do not match the model");
    result.solver_options = model.solver_options();
    const auto& parameters = model.parameters();
    const std::size_t count = parameters.component_count();
    const auto sidx = static_cast<std::size_t>(std::max_element(solvent.begin(), solvent.end()) - solvent.begin());
    if (std::count(solvent.begin(), solvent.end(), 1.0) != 1) {
        result.message = "input actions require one pure solvent"; return result;
    }
    const double mass = parameters.molar_mass_kg_per_mol[sidx];
    const double nu = double(stoichiometry[0] + stoichiometry[1]), s = nu * mass * molality;
    const double ds = nu * mass * direction.molality;
    auto composition = formula_composition(solvent, parameters, ions, stoichiometry, molality);
    PressureActionDirection finite{direction.temperature, direction.pressure, std::vector<double>(count, 0.0), {}};
    finite.composition[sidx] = -nu * mass * direction.molality / ((1.0 + s) * (1.0 + s));
    for (std::size_t i = 0; i < 2; ++i)
        finite.composition[ions[i]] = stoichiometry[i] * mass * direction.molality / ((1.0 + s) * (1.0 + s));
    const auto actions = [&](const std::vector<double>& x, const State& state, const PressureActionDirection& d,
                             const std::vector<std::size_t>& consumed) {
        const auto batch = model.pressure_actions(temperature, pressure, x, DensityBranch::Liquid,
            {PropertyObservable::LogFugacityCoefficient}, {d}, state.molar_density);
        result.diagnostics.push_back(batch.diagnostics);
        std::vector<double> columns(count, 0.0);
        if (batch.density_status != PropertyActionStatus::Available) {
            result.status = batch.density_status; result.message = batch.density_message; return std::vector<double>{};
        }
        for (const auto i : consumed) {
            const auto found = std::find_if(batch.results.begin(), batch.results.end(), [&](const auto& r) { return r.component == i; });
            if (found == batch.results.end() || found->status != PropertyActionStatus::Available || !found->action) {
                result.status = found == batch.results.end() ? PropertyActionStatus::CentralPropertyUnavailable : found->status;
                result.message = "unavailable electrolyte log-fugacity input action"; return std::vector<double>{};
            }
            columns[i] = *found->action;
        }
        return columns;
    };
    const auto dg = actions(composition, activity.state, finite, {sidx, ions[0], ions[1]});
    if (dg.empty()) return result;
    double pure_column = 0.0;
    if (direction.temperature != 0.0 || direction.pressure != 0.0) {
        const PressureActionDirection d{direction.temperature, direction.pressure, std::vector<double>(count, 0.0), {}};
        std::array<double, 2> previous{};
        for (std::size_t k = 0; k < REFERENCE_MOLALITIES.size(); ++k) {
            const bool witness = k + 1 == REFERENCE_MOLALITIES.size();
            const double error = standard_deviation * std::max(std::abs(result.reference_difference[0]), std::abs(result.reference_difference[1]));
            if (witness && std::isfinite(error) && error <= policy::REFERENCE_CONVERGENCE_TOLERANCE.value) break;
            if (witness) result.witness_evaluated = true;
            const auto x = formula_composition(solvent, parameters, ions, stoichiometry, REFERENCE_MOLALITIES[k]);
            const auto extra = k < reference.ladder_states.size() ? std::optional<State>{}
                : std::optional<State>{model.state_at_pressure(temperature, pressure, x, DensityBranch::Liquid, reference.ladder_states.back().molar_density)};
            const auto& state = extra ? *extra : reference.ladder_states[k];
            if (const auto cause = properties::vapor_side(state)) { result.message = *cause; return result; }
            const auto column = actions(x, state, d, {ions[0], ions[1]});
            if (column.empty()) return result;
            for (std::size_t i = 0; i < 2; ++i) {
                if (witness) result.witness_difference[i] = column[ions[i]] - previous[i];
                else {
                    if (k) result.reference_difference[i] = column[ions[i]] - previous[i];
                    result.reference_column[i] = previous[i] = column[ions[i]];
                }
            }
            if (witness && standard_deviation * std::max(std::abs(result.witness_difference[0]), std::abs(result.witness_difference[1]))
                > policy::REFERENCE_CONVERGENCE_TOLERANCE.value) {
                result.message = "reference input action did not converge at the witness"; return result;
            }
        }
        if (!reference.pure_solvent_state) { result.message = "pure solvent input action unavailable"; return result; }
        const auto pure = actions(solvent, *reference.pure_solvent_state, d, {sidx});
        if (pure.empty()) return result;
        pure_column = pure[sidx];
    }
    double mean = -ds / (1.0 + s);
    for (std::size_t i = 0; i < 2; ++i) mean += stoichiometry[i] * (dg[ions[i]] - result.reference_column[i]) / nu;
    result.mean_ionic = activity.mean_ionic_molality_scale * mean;
    if (activity.osmotic_coefficient && std::holds_alternative<double>(*activity.osmotic_coefficient))
        result.osmotic = -(finite.composition[sidx] / composition[sidx] + dg[sidx] - pure_column
            + std::get<double>(*activity.osmotic_coefficient) * ds) / s;
    result.status = std::isfinite(result.mean_ionic) ? PropertyActionStatus::Available : PropertyActionStatus::Nonfinite;
    return result;
}

}  // namespace epcsaft
