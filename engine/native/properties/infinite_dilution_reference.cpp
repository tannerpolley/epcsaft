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
#include <tuple>
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
        } catch (const DensityStructureError&) { throw;
        } catch (const std::runtime_error& error) {
            if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error)) scope.finish(failure->record);
            scope.probe("anchor_failure", {*prior_density}, 0, error.what());
        }
    }
    return model.state_at_pressure(temperature, pressure, composition, branch);
}

namespace {

// EqIDs single_ion_activity_coefficient, reaction_activity_quotient: retain log differences.
double activity_log(double finite, double reference, double log_concentration = 0.0) {
    return log_concentration == 0.0 ? finite - reference : (finite - reference) + log_concentration;
}

std::array<double, 2> ion_log_fugacity(const State& state, const std::array<std::size_t, 2>& ions) {
    if (!state.log_fugacity_coefficient) throw std::runtime_error("reference fugacity is unavailable");
    return {require_component_property((*state.log_fugacity_coefficient)[ions[0]]),
        require_component_property((*state.log_fugacity_coefficient)[ions[1]])};
}

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

// Reference-column scale of active coordinate a: max(1, |theta|); d_V (m3/mol) uses a 1 L/mol floor.
double parameter_scale(const HelmholtzModel& model, const std::vector<double>& values, std::size_t a) {
    return std::max(model.active_parameters()[a].family == ParameterFamily::DielectricVolume ? 1.0e-3 : 1.0, std::abs(values[a]));
}

}  // namespace

ResidualEnergyAction reference_partial_from_state(const HelmholtzModel& model, double temperature,
    const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    std::size_t temperature_order, std::size_t pressure_order,
    std::optional<std::size_t> parameter, const State& state, bool infinite_dilution) {
    const std::size_t order = temperature_order + pressure_order + std::size_t(parameter.has_value());
    ResidualEnergyAction result;
    const double pressure_scale = state.pressure_density_derivative * state.molar_density;
    const double theta_scale = parameter ? parameter_scale(model, model.active_parameter_values(), *parameter) : 1.0;
    const double scale = std::pow(temperature, int(temperature_order)) * std::pow(pressure_scale, int(pressure_order)) * theta_scale;
    double contraction = 0.0;
    const std::size_t paths = order ? std::size_t(1) << (order - 1) : 1;
    for (std::size_t mask = 0; mask < paths; ++mask) {
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
            composition, stoichiometry, direction, order, infinite_dilution);
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
    result.action = contraction / (double(paths) * scale);
    if (!std::isfinite(*result.action)) { result.action.reset(); result.status = PropertyActionStatus::Nonfinite; }
    return result;
}

struct NeutralReferenceTrace final {
    double pressure = 0.0;
    std::size_t solvent_index = 0;
    std::vector<double> composition;
    State state;
    std::vector<std::optional<State>> pure_state;
    std::size_t temperature_order = 0;
    bool parameter_jets = false;
};

void reference_errors(const HelmholtzModel& model, double temperature, const std::vector<double>& nu,
    NeutralReferenceReactionDiagnostic& diagnostic) {
    try {
        if (!diagnostic.trace) throw std::runtime_error("reference state unavailable");
        const auto& trace = *diagnostic.trace;
        diagnostic.numerical_error = {0.0, 0.0};
        std::array<WiderScalar, 2> replayed_sum{};
        const auto accumulate = [&](const State& state, const std::vector<double>& composition,
            const std::vector<double>& insertion, bool dilute) {
            ReferenceTaylorError constant_error, error;
            PressureActionDirection direction; direction.temperature = temperature;
            const auto constant = model.reference_taylor(temperature, 1 / state.molar_density, composition, insertion, {}, 0, dilute, &constant_error);
            const auto value = model.reference_taylor(temperature, 1 / state.molar_density, composition, insertion, direction, 1, dilute, &error);
            if (constant.status != PropertyActionStatus::Available || value.status != PropertyActionStatus::Available)
                throw std::runtime_error("reference error column unavailable");
            const auto cp = reference_partial_from_state(model, temperature, composition, insertion, 0, 1, {}, state, dilute);
            const auto ctp = reference_partial_from_state(model, temperature, composition, insertion, 1, 1, {}, state, dilute);
            if (!cp.action || !ctp.action) throw std::runtime_error("reference pressure sensitivity unavailable");
            double b0 = wide_math::abs(error.pressure[0] - trace.pressure) + error.pressure_floor[0]
                + std::numeric_limits<double>::epsilon() * wide_math::abs(trace.pressure);
            if (state.density_convergence && state.density_convergence->floor)
                b0 += *state.density_convergence->floor * GAS_CONSTANT_J_PER_MOL_K * temperature * state.molar_density;
            const double b1 = wide_math::abs(error.pressure[1]) + error.pressure_floor[1];
            diagnostic.numerical_error[0] += constant_error.fixed[0] + wide_math::abs(*cp.action) * b0;
            diagnostic.numerical_error[1] += (error.fixed[1] + wide_math::abs(*cp.action) * b1) / wide_math::abs(error.temperature_coefficient)
                + wide_math::abs(*ctp.action) * b0;
            replayed_sum[0] += constant.contraction[0];
            replayed_sum[1] += static_cast<WiderScalar>(value.contraction[1]) / temperature;
        };
        auto dilute = nu;
        for (std::size_t i = 0; i < trace.pure_state.size(); ++i) if (trace.pure_state[i]) dilute[i] = 0;
        if (std::any_of(dilute.begin(), dilute.end(), [](double v) { return v != 0; }))
            accumulate(trace.state, trace.composition, dilute, true);
        for (std::size_t i = 0; i < trace.pure_state.size(); ++i) if (trace.pure_state[i] && nu[i] != 0) {
            std::vector<double> composition(nu.size(), 0), insertion(nu.size(), 0);
            composition[i] = 1; insertion[i] = nu[i];
            accumulate(*trace.pure_state[i], composition, insertion, false);
        }
        if (!diagnostic.terminal_contraction || !diagnostic.terminal_contraction_temperature_derivative)
            throw std::runtime_error("reference consumed center unavailable");
        diagnostic.numerical_error[0] = static_cast<double>(diagnostic.numerical_error[0] + wide_math::abs(static_cast<WiderScalar>(*diagnostic.terminal_contraction) - replayed_sum[0]));
        diagnostic.numerical_error[1] = static_cast<double>(diagnostic.numerical_error[1] + wide_math::abs(static_cast<WiderScalar>(*diagnostic.terminal_contraction_temperature_derivative) - replayed_sum[1]));
        if (diagnostic.observed_terminal_change) diagnostic.numerical_error[0] += wide_math::abs(*diagnostic.observed_terminal_change);
        if (diagnostic.observed_temperature_derivative_change) diagnostic.numerical_error[1] += wide_math::abs(*diagnostic.observed_temperature_derivative_change);
        for (double e : diagnostic.numerical_error) if (!wide_math::isfinite(e) || e < 0) throw std::runtime_error("nonfinite reference error");
        diagnostic.numerical_error_cause.clear();
    } catch (const std::exception& error) {
        diagnostic.numerical_error = {NAN, NAN}; diagnostic.numerical_error_cause = error.what();
    }
}

namespace {

struct NeutralReferenceTraceResult final {
    std::shared_ptr<const NeutralReferenceTrace> trace;
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
        validate_neutral_reference_declaration(
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
        trace.solvent_index = *declaration.solvent_index;
        trace.composition.assign(components, 0.0);
        trace.composition[trace.solvent_index] = 1.0;
        try {
            trace.state = properties::close_density(model, temperature, pressure,
                trace.composition, DensityBranch::Liquid, solvent_anchor);
        } catch (const std::exception& error) {
            result.status = NeutralReferenceStatus::DensityUnavailable;
            result.message = std::string("neutral reference density unavailable: ") + error.what();
            return result;
        }
        const bool unique_root = trace.state.density_branch == DensityBranch::Unspecified
            && trace.state.stable_root_count.value_or(0) == 1;
        if (trace.state.density_branch != DensityBranch::Liquid && !unique_root) {
            result.status = NeutralReferenceStatus::BranchChange;
            result.message = "neutral reference left the liquid branch";
            return result;
        }
        if (!properties::regular_density(trace.state.molar_density,
            trace.state.pressure_density_derivative, trace.state.density_convergence)) {
            result.status = NeutralReferenceStatus::DensityUnavailable;
            result.message = "neutral reference state has invalid density response";
            return result;
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
                                + " branch not established by the finite search: " + (loop.unresolved.empty()
                                    ? "single stable root across the mechanical loop" : loop.unresolved.front().cause));
                    }
                    trace.pure_state[i] = std::move(state);
                } catch (const std::exception& error) {
                    result.status = NeutralReferenceStatus::DensityUnavailable;
                    result.message = "pure reference species " + std::to_string(i) + ": " + error.what();
                    return result;
                }
            }
        }
        result.trace = std::make_shared<const NeutralReferenceTrace>(std::move(trace));
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
    const NeutralReferenceTrace& trace, const std::vector<double>& stoichiometry,
    std::size_t a, std::size_t b, std::optional<std::size_t> parameter) {
    if (trace.pure_state.empty())
        return reference_partial_from_state(model, temperature, trace.composition,
            stoichiometry, a, b, parameter, trace.state, true);
    auto dilute = stoichiometry;
    for (std::size_t i = 0; i < dilute.size(); ++i)
        if (trace.pure_state[i]) dilute[i] = 0.0;
    ResidualEnergyAction result;
    result.status = PropertyActionStatus::Available;
    result.action = 0.0;
    if (std::any_of(dilute.begin(), dilute.end(), [](double v) { return v != 0.0; })) {
        result = reference_partial_from_state(model, temperature, trace.composition,
            dilute, a, b, parameter, trace.state, true);
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

// EqID: reference_single_insertion_taylor; each requested column uses the same solvent closure.
NeutralReferenceReactionDiagnostic::Derivative reference_column(const HelmholtzModel& model, double temperature,
    const NeutralReferenceTrace& trace, const std::vector<double>& stoichiometry, std::size_t a, std::size_t b,
    std::optional<std::size_t> parameter) {
    core::Scope scope("reference_column");
    NeutralReferenceReactionDiagnostic::Derivative column;
    column.temperature_order = a; column.pressure_order = b; column.parameter = parameter;
    try {
        const auto partial = reference_partial_from_trace(model, temperature, trace, stoichiometry, a, b, parameter);
        column.value = partial.action;
        column.diagnostics = partial.diagnostics;
        column.action_status = partial.status;
        if (!column.value) column.message = "reference action status " + std::to_string(int(partial.status));
    } catch (const std::exception& error) { column.message = error.what(); }
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
    // EqID: reference_single_insertion_taylor; raw partials at the ion-free solvent state.
    NeutralReferenceReactionDiagnostic result;
    result.solver_options = model.solver_options();
    result.reference_pressure_pa = trace.pressure;
    try {
        bool solvent_consumed = false;
        for (std::size_t i = 0; i < stoichiometry.size(); ++i)
            solvent_consumed = solvent_consumed || (stoichiometry[i] != 0.0
                && (trace.pure_state.empty() || !trace.pure_state[i]));
        if (const auto cause = solvent_consumed ? properties::vapor_side(trace.state) : std::nullopt) {
            result.status = NeutralReferenceStatus::DensityUnavailable;
            result.message = "neutral reference solvent species " + std::to_string(trace.solvent_index) + ": " + *cause;
            return result;
        }
        const auto value = reference_partial_from_trace(model, temperature, trace, stoichiometry, 0, 0, std::nullopt);
        result.terminal_contraction = value.action;
        result.diagnostics = value.diagnostics;
        result.density_convergence = trace.state.density_convergence;
        if (!value.action) {
            result.status = NeutralReferenceStatus::Nonfinite;
            result.message = "neutral reference insertion unavailable: " + std::to_string(int(value.status));
            return result;
        }
        const auto admit = [&](std::size_t a, std::size_t b, std::optional<std::size_t> parameter) {
            return result.lower_derivatives.emplace_back(reference_column(model, temperature, trace, stoichiometry,
                a, b, parameter)).value;
        };
        if (trace.temperature_order > 0) {
            result.terminal_contraction_temperature_derivative = admit(1, 0, std::nullopt);
            if (trace.temperature_order > 1)
                result.terminal_contraction_temperature_second_derivative = admit(2, 0, std::nullopt);
            if (pressure_jets)
                result.terminal_contraction_temperature_pressure_derivative = admit(1, 1, std::nullopt);
        }
        for (std::size_t a = 0; trace.parameter_jets && a < model.active_parameters().size(); ++a) {
            result.terminal_contraction_parameter_derivative.push_back(admit(0, 0, a));
            if (trace.temperature_order > 0)
                result.terminal_contraction_temperature_parameter_derivative.push_back(admit(1, 0, a));
        }
        if (reference.reference_pressure_pa.has_value()) {
            result.terminal_contraction_pressure_derivative = 0.0;
            result.terminal_contraction_pressure_second_derivative = 0.0;
            result.terminal_contraction_temperature_pressure_derivative = 0.0;
        } else if (pressure_jets) {
            result.terminal_contraction_pressure_derivative = admit(0, 1, std::nullopt);
            result.terminal_contraction_pressure_second_derivative = admit(0, 2, std::nullopt);
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
    double molality, const std::optional<std::vector<double>>& finite_solute_molalities, double* declared_solute_ratio
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
    if (declared_solute_ratio) *declared_solute_ratio = solvent_mass * molality * (static_cast<double>(stoichiometry[0]) + stoichiometry[1]);
    const double solvent_moles = 1.0 / solvent_mass;
    const double total_moles = solvent_moles
        + molality * (static_cast<double>(stoichiometry[0]) + stoichiometry[1]);
    std::vector<double> composition(solvent_fractions.size(), 0.0);
    for (std::size_t component = 0; component < composition.size(); ++component) {
        composition[component] = solvent_fractions[component] * solvent_moles / total_moles;
    }
    composition[ion_indices[0]] = stoichiometry[0] * molality / total_moles;
    composition[ion_indices[1]] = stoichiometry[1] * molality / total_moles;
    if (finite_solute_molalities) {
        const auto& m = *finite_solute_molalities;
        if (m.size() != composition.size()
            || std::count(solvent_fractions.begin(), solvent_fractions.end(), 1.0) != 1)
            throw std::invalid_argument("full solute molalities require one pure solvent and one entry per component");
        const auto solvent = std::size_t(std::max_element(solvent_fractions.begin(), solvent_fractions.end()) - solvent_fractions.begin());
        if (parameters.charge_number[solvent] != 0 || m[solvent] != 0.0)
            throw std::invalid_argument("full solute molalities require an uncharged solvent with zero solute molality");
        double total = 0.0, charge = 0.0;
        for (std::size_t i = 0; i < m.size(); ++i) {
            if (!std::isfinite(m[i]) || m[i] < 0.0)
                throw std::invalid_argument("full solute molalities must be finite and nonnegative");
            total += m[i];
            charge += m[i] * parameters.charge_number[i];
        }
        const double mass = parameters.molar_mass_kg_per_mol[solvent];
        if (declared_solute_ratio) *declared_solute_ratio = mass * total;
        const double denominator = 1.0 + mass * total;
        if (!std::isfinite(denominator) || !std::isfinite(charge)
            || m[ion_indices[0]] <= 0.0 || m[ion_indices[1]] <= 0.0
            || std::abs(mass * charge / denominator) > policy::COMPOSITION_SUM_TOLERANCE.value
            || std::abs(mass * (m[ion_indices[0]] - stoichiometry[0] * molality) / denominator)
                > policy::COMPOSITION_SUM_TOLERANCE.value)
            throw std::invalid_argument("full solute molalities must be neutral and match the selected formula cation molality");
        composition[solvent] = 1.0 / denominator;
        for (std::size_t i = 0; i < m.size(); ++i)
            if (i != solvent) composition[i] = mass * m[i] / denominator;
    }
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

std::optional<std::vector<NeutralReferenceReactionDiagnostic>> HelmholtzModel::stored_references(
    const std::vector<double>& inputs, std::size_t first, std::size_t last) const {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& [stored, result] : references_) {
        bool covers = stored.size() == inputs.size();
        for (std::size_t i = 0; covers && i < inputs.size(); ++i)
            covers = first <= i && i < last ? inputs[i] <= stored[i] : inputs[i] == stored[i];
        if (covers) return result;
    }
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
    std::size_t reference_order,
    std::vector<std::vector<NeutralReferenceReactionDiagnostic::Derivative>>* species_values
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
    std::vector<double> inputs{temperature, system_pressure, double(solvent_anchor.has_value()),
        solvent_anchor.value_or(0.0), double(pressure_jets), double(temperature_order), double(parameter_jets),
        double(reference_order), double(declaration.solvent_index.has_value()), double(declaration.solvent_index.value_or(0)),
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
    if (!species_values) if (auto stored = model.stored_references(inputs, 4, 8)) return std::move(*stored);
    if (species_values) species_values->resize(stoichiometries.size());
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
        if (!found->result.trace) {
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
        diagnostic.trace = found->result.trace;
        diagnostic.solvent_density = found->result.trace->state.molar_density;
        // Additional consumed columns at the same solvent state.
        const auto& trace = *found->result.trace;
        if (species_values && diagnostic.status == NeutralReferenceStatus::Available) {
            auto& columns = (*species_values)[row];
            columns.resize(stoichiometries[row].size());
            for (std::size_t i = 0; i < columns.size(); ++i) {
                std::vector<double> insertion(columns.size(), 0.0);
                insertion[i] = 1.0;
                columns[i] = reference_column(model, temperature, trace, insertion, 0, 0, std::nullopt);
            }
        }
        const auto values = model.active_parameter_values();
        for (std::size_t c = 0; c <= std::size_t(parameter_jets); ++c)
            for (std::size_t d = 2; d <= reference_order; ++d)
                for (std::size_t b = 0; b <= (pressure_jets ? d-c : 0); ++b) {
                    const std::size_t a = d-b-c;
                    if (d == 2 && !(c && b)) continue;  // existing lower slopes keep their admitted owner
                    for (std::size_t coordinate = 0; coordinate < (c ? values.size() : 1); ++coordinate) {
                        const auto parameter = c ? std::optional<std::size_t>(coordinate) : std::nullopt;
                        NeutralReferenceReactionDiagnostic::Derivative column;
                        if (b && reference.reference_pressure_pa) {
                            column.temperature_order = a; column.pressure_order = b; column.parameter = parameter;
                            column.value = 0.0;
                        } else {
                            column = reference_column(model, temperature, trace, stoichiometries[row], a, b, parameter);
                        }
                        if (!column.value) column.message = "reaction " + std::to_string(row) + " C_"
                            + std::to_string(a) + std::to_string(b) + (c ? " coordinate " + std::to_string(coordinate) : "")
                            + ": " + column.message;
                        diagnostic.derivatives.push_back(std::move(column));
                    }
                }

    }
    if (species_values) for (std::size_t r = 0; r < result.size(); ++r) {
        auto& columns = (*species_values)[r];
        if (columns.empty()) {
            columns.resize(stoichiometries[r].size());
            for (auto& column : columns) column.message = result[r].message;
        }
    }
    if (!species_values) model.store_references(std::move(inputs), result);
    return result;
}

SpeciesActivities species_activities(const std::vector<double>& composition, const State* state,
    const std::vector<NeutralReferenceReactionDiagnostic::Derivative>& references,
    const Parameters* parameters, const NeutralReferenceDeclaration* declaration, ReferenceSourceBasis basis) {
    SpeciesActivities result;
    if (declaration) {
        result.solvent_index = declaration->solvent_index;
        result.standard_molality_mol_per_kg = declaration->standard_molality_mol_per_kg;
    }
    result.log_activities.resize(composition.size()); result.activities.resize(composition.size());
    result.causes.resize(composition.size());
    for (std::size_t i = 0; i < composition.size(); ++i) {
        try {
            if (!references.empty() && !references.at(i).value)
                throw std::runtime_error("species reference unavailable: " + references[i].message);
            const double standard = references.empty() ? 0.0 : *references[i].value;
            double log_concentration = std::log(composition[i]);
            if (declaration && basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution && i != *declaration->solvent_index)
                log_concentration -= std::log(parameters->molar_mass_kg_per_mol[*declaration->solvent_index]
                    * declaration->standard_molality_mol_per_kg);
            const double value = composition[i] == 0.0 ? -INFINITY : state ? activity_log(std::log(state->molar_density / 1.0)
                + require_component_property(state->residual_chemical_potential_over_rt.at(i)), standard, log_concentration)
                : log_concentration;
            if (std::isnan(value) || value == INFINITY) throw std::runtime_error("nonfinite log activity");
            result.log_activities[i] = value;
            const double activity = state ? std::exp(value) : composition[i];
            if (!std::isfinite(activity) || (activity == 0.0 && composition[i] > 0.0))
                result.causes[i] = "activity exponentiation overflow or underflow";
            else result.activities[i] = activity;
        } catch (const std::exception& error) { result.causes[i] = error.what(); }
    }
    return result;
}

std::vector<HenryConstant> henry_constants(const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& solvent_fractions, const std::vector<std::size_t>& solute_indices) {
    core::Scope scope("henry_constants");
    const auto& parameters = model.parameters();
    const std::size_t components = parameters.component_count();
    const auto solvent = std::size_t(std::find(solvent_fractions.begin(), solvent_fractions.end(), 1.0)
        - solvent_fractions.begin());
    auto solutes = solute_indices;
    std::sort(solutes.begin(), solutes.end());
    if (solvent_fractions.size() != components || solvent == components
        || std::count(solvent_fractions.begin(), solvent_fractions.end(), 0.0) != std::ptrdiff_t(components - 1)
        || solutes.empty() || solutes.back() >= components
        || std::adjacent_find(solutes.begin(), solutes.end()) != solutes.end()
        || std::binary_search(solutes.begin(), solutes.end(), solvent))
        throw std::invalid_argument("henry constants need one pure solvent and distinct known solutes");
    std::vector<HenryConstant> result(solutes.size());
    for (std::size_t row = 0; row < solutes.size(); ++row) result[row].component_index = solutes[row];
    // A charge-balanced pair of ion slots passes the ratio validation; Henry admits all-neutral models only.
    if (std::any_of(parameters.charge_number.begin(), parameters.charge_number.end(), [](int z) { return z != 0; })) {
        for (auto& henry : result) {
            henry.status = NeutralReferenceStatus::NonNeutralDeclaration;
            henry.message = "henry constants require an all-neutral model";
        }
        return result;
    }
    NeutralReferenceDeclaration declaration{solvent, std::vector<double>(components, 1.0), 1.0, {}};
    declaration.solute_molality_ratios[solvent] = 0.0;
    std::vector<std::vector<double>> rows(solutes.size(), std::vector<double>(components, 0.0));
    for (std::size_t row = 0; row < solutes.size(); ++row) rows[row][solutes[row]] = 1.0;
    const auto references = reaction_references(model, temperature, pressure, declaration, rows,
        std::vector<ReactionReference>(rows.size()), /*pressure_jets=*/false);
    for (std::size_t row = 0; row < solutes.size(); ++row) {
        auto& henry = result[row];
        const auto& reference = references[row];
        henry.raw_log_reference = reference.terminal_contraction;
        henry.observed_terminal_change = reference.observed_terminal_change;
        henry.status = reference.status;
        henry.message = reference.message;
        if (!reference.terminal_contraction) continue;
        // EqID: henry_constant_px
        const double value = GAS_CONSTANT_J_PER_MOL_K * temperature * std::exp(*reference.terminal_contraction);
        if (std::isfinite(value) && value > 0.0) henry.henry_constant_px = value;
        else {
            henry.status = NeutralReferenceStatus::Nonfinite;
            henry.message = "henry constant transform is nonfinite or zero";
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

InfiniteDilutionReference solvent_reference(const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& fractions, std::optional<double> density) {
    InfiniteDilutionReference result;
    result.solver_options = model.solver_options();
    result.model = &model; result.temperature = temperature; result.pressure = pressure;
    result.solvent_fractions = fractions;
    const auto& parameters = model.parameters();
    double sum = 0.0;
    std::size_t solvent = fractions.size(), occupied = 0;
    for (std::size_t i = 0; i < fractions.size(); ++i) {
        if (!std::isfinite(fractions[i]) || fractions[i] < 0.0) {
            result.solvent_cause = "PureSolventUnavailable: invalid solvent fractions"; return result;
        }
        sum += fractions[i];
        if (fractions[i] > 0.0) { solvent = i; ++occupied; }
    }
    if (fractions.size() != parameters.component_count() || occupied != 1
        || std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value
        || parameters.charge_number[solvent] != 0
        || !std::isfinite(parameters.molar_mass_kg_per_mol[solvent]) || parameters.molar_mass_kg_per_mol[solvent] <= 0.0) {
        result.solvent_cause = "PureSolventUnavailable: require one neutral solvent with positive molar mass"; return result;
    }
    std::vector<double> pure(fractions.size(), 0.0);
    pure[solvent] = 1.0;
    try {
        result.pure_solvent_state = properties::close_density(model, temperature, pressure, pure, DensityBranch::Liquid, density);
        const auto& state = *result.pure_solvent_state;
        if (const auto cause = properties::vapor_side(state)) result.solvent_cause = "VaporSideState: pure solvent: " + *cause;
        else if (!properties::regular_density(state.molar_density, state.pressure_density_derivative, state.density_convergence))
            result.solvent_cause = "DensityUnavailable: pure solvent regularity";
        else {
            result.solvent_status = PropertyActionStatus::Available;
            if (state.log_fugacity_coefficient) result.pure_solvent_log_fugacity = (*state.log_fugacity_coefficient)[solvent];
            return result;
        }
    } catch (const std::runtime_error& error) {
        if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error)) result.solvent_convergence = failure->record;
        result.solvent_cause = "DensityUnavailable: pure solvent: " + std::string(error.what());
    }
    result.solvent_status = PropertyActionStatus::DensityUnavailable;
    return result;
}

SolventActivity solvent_activity(const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, const InfiniteDilutionReference& reference, std::optional<State> finite,
    std::optional<double> formula_solute_ratio) {
    SolventActivity result;
    const std::size_t columns = model.active_parameters().size();
    result.solvent_activity_parameter_derivative.resize(columns);
    result.osmotic_parameter_derivative.resize(columns);
    result.parameter_status.assign(columns, PropertyActionStatus::CentralPropertyUnavailable);
    result.parameter_cause.resize(columns);
    if (reference.model && (reference.model != &model || reference.temperature != temperature || reference.pressure != pressure
        || !(reference.solver_options == model.solver_options()))) {
        result.cause = "SolventReferenceMismatch"; return result;
    }
    if (reference.solvent_status != PropertyActionStatus::Available || !reference.pure_solvent_state) {
        if (reference.solvent_status != PropertyActionStatus::Available) result.status = reference.solvent_status;
        result.cause = reference.solvent_cause.empty() ? "PureSolventUnavailable: missing pure reference" : reference.solvent_cause;
        result.state.density_convergence = reference.solvent_convergence; return result;
    }
    const auto solvent = std::size_t(std::max_element(reference.solvent_fractions.begin(), reference.solvent_fractions.end()) - reference.solvent_fractions.begin());
    if (composition.size() != model.parameters().component_count()) throw std::invalid_argument("invalid composition size");
    double sum = 0.0, charge = 0.0, solute = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        if (!std::isfinite(composition[i]) || composition[i] < 0.0) throw std::invalid_argument("invalid composition fraction");
        sum += composition[i]; charge += composition[i] * model.parameters().charge_number[i];
        if (i != solvent) solute += composition[i];
    }
    if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value
        || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value) throw std::invalid_argument("invalid normalized neutral composition");
    if (composition[solvent] <= 0.0) { result.cause = "PureSolventUnavailable: zero solvent fraction"; return result; }
    const auto& pure = *reference.pure_solvent_state;
    // Preserve the one-salt declaration's particle ratio without a round trip through mole fractions.
    const double D = formula_solute_ratio.value_or(solute / composition[solvent]);
    if (D == 0.0) {
        result.osmotic_cause = "OsmoticUndefinedPureSolvent";
        result.osmotic_coefficient = core::Stop::converged_to_floor;
        result.state = pure; result.solvent_activity = 1.0; result.status = PropertyActionStatus::Available;
        result.parameter_status.assign(columns, PropertyActionStatus::Available);
        result.solvent_activity_parameter_derivative.assign(columns, 0.0); return result;
    }
    try {
        result.state = finite ? std::move(*finite) : properties::close_density(model, temperature, pressure, composition, DensityBranch::Liquid, std::nullopt);
    } catch (const std::runtime_error& error) {
        result.status = PropertyActionStatus::DensityUnavailable;
        result.cause = "DensityUnavailable: mixture: " + std::string(error.what());
        if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error)) result.state.density_convergence = failure->record;
        return result;
    }
    const auto& state = result.state;
    if (const auto cause = properties::vapor_side(state)) { result.cause = "VaporSideState: mixture: " + *cause; return result; }
    if (!properties::regular_density(state.molar_density, state.pressure_density_derivative, state.density_convergence)) {
        result.status = PropertyActionStatus::DensityUnavailable; result.cause = "DensityUnavailable: mixture regularity"; return result;
    }
    // EqID: practical_osmotic_coefficient; extensive Helmholtz Euler identity at common pressure.
    WiderScalar chemical_difference = static_cast<WiderScalar>(state.residual_helmholtz) - pure.residual_helmholtz;
    chemical_difference += pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature)
        * ((static_cast<WiderScalar>(pure.molar_density) - state.molar_density) / (state.molar_density * static_cast<WiderScalar>(pure.molar_density)));
    for (std::size_t i = 0; i < composition.size(); ++i) if (i != solvent && composition[i] > 0.0)
        chemical_difference -= composition[i] * (static_cast<WiderScalar>(require_component_property(state.residual_chemical_potential_over_rt[i]))
            - require_component_property(pure.residual_chemical_potential_over_rt[solvent]));
    chemical_difference /= composition[solvent];
    const double L = -wide_math::log1p(D) + static_cast<double>(chemical_difference)
        + wide_math::log1p((state.molar_density - pure.molar_density) / pure.molar_density);
    const double activity = wide_math::exp(L);
    if (!wide_math::isfinite(L) || !wide_math::isfinite(activity)) { result.status = PropertyActionStatus::Nonfinite; result.cause = "Nonfinite: solvent activity"; return result; }
    result.solvent_activity = activity; result.status = PropertyActionStatus::Available;
    result.parameter_status.assign(columns, PropertyActionStatus::Available);
    const bool admitted = D / model.parameters().molar_mass_kg_per_mol[solvent] >= policy::OSMOTIC_MOLALITY_FLOOR.value;
    if (admitted) result.osmotic_coefficient = -L / D;
    else {
        result.osmotic_coefficient = core::Stop::converged_to_floor;
        result.osmotic_cause = "OsmoticMolalityFloor";
    }
    // EqID: electrolyte_activity_parameter_derivative; physical coordinates at fixed T/P/x.
    for (std::size_t a = 0; a < columns; ++a) {
        const auto index = solvent * columns + a;
        const auto f = state.log_fugacity_parameter_derivative.at(index), p = pure.log_fugacity_parameter_derivative.at(index);
        if (f && p && std::isfinite(*f - *p)) {
            result.solvent_activity_parameter_derivative[a] = *result.solvent_activity * (*f - *p);
            if (admitted) result.osmotic_parameter_derivative[a] = -(*f - *p) / D;
        } else {
            result.parameter_status[a] = !f ? state.component_property_status[solvent] : !p ? pure.component_property_status[solvent] : PropertyActionStatus::Nonfinite;
            result.parameter_cause[a] = !f ? "mixture solvent fugacity column unavailable" : "pure solvent fugacity column unavailable or nonfinite";
        }
    }
    return result;
}

InfiniteDilutionReference infinite_dilution_reference(
    const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& solvent_fractions, const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    std::optional<double> solvent_density
) {
    core::Scope scope("infinite_dilution_reference");
    // Reuse the finite-composition declaration checks without evaluating a finite state.
    (void)formula_composition(solvent_fractions, model.parameters(), ion_indices, stoichiometry, 1.0);
    for (std::size_t i = 0; i < solvent_fractions.size(); ++i)
        if (solvent_fractions[i] > 0.0 && model.parameters().charge_number[i] != 0)
            throw IonicReferenceBackground("IonicReferenceBackground: infinite dilution requires an ion-free solvent");
    InfiniteDilutionReference result;
    result.solver_options = model.solver_options();
    result.model = &model; result.temperature = temperature; result.pressure = pressure;
    result.solvent_fractions = solvent_fractions; result.ion_indices = ion_indices;
    result.solvent_state = properties::close_density(model, temperature, pressure,
        solvent_fractions, DensityBranch::Liquid, solvent_density);
    const auto& state = *result.solvent_state;
    if (const auto cause = properties::vapor_side(state)) throw std::runtime_error(*cause);
    if (!properties::regular_density(state.molar_density, state.pressure_density_derivative, state.density_convergence))
        throw std::runtime_error("DensityUnavailable: infinite-dilution solvent regularity");
    const std::size_t columns = model.active_parameters().size();
    result.parameter_derivative.resize(columns);
    for (std::size_t i = 0; i < 2; ++i) {
        std::vector<double> insertion(solvent_fractions.size(), 0.0);
        insertion[ion_indices[i]] = 1.0;
        const auto value = reference_partial_from_state(model, temperature, solvent_fractions,
            insertion, 0, 0, std::nullopt, state, true);
        if (!value.action) throw std::runtime_error("infinite-dilution insertion unavailable: " + std::to_string(int(value.status)));
        // L_i = ln(rho/rho0) + b_i; ln(phi_i) = L_i - ln(P/(RT rho0)).
        result.log_fugacity[i] = *value.action - std::log(pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature));
        result.diagnostics[i] = value.diagnostics;
        for (std::size_t a = 0; a < columns; ++a)
            result.parameter_derivative[a][i] = reference_partial_from_state(model, temperature, solvent_fractions,
                insertion, 0, 0, a, state, true);
    }
    if (std::count_if(solvent_fractions.begin(), solvent_fractions.end(), [](double x) { return x > 0.0; }) == 1) {
        const auto solvent = std::size_t(std::max_element(solvent_fractions.begin(), solvent_fractions.end()) - solvent_fractions.begin());
        result.pure_solvent_state = state;
        result.solvent_status = PropertyActionStatus::Available;
        if (state.log_fugacity_coefficient) result.pure_solvent_log_fugacity = (*state.log_fugacity_coefficient)[solvent];
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
    const InfiniteDilutionReference& reference, std::optional<double> finite_density,
    const std::optional<std::vector<double>>& finite_solute_molalities
) {
    core::Scope scope("electrolyte_activity");
    if (!reference.ion_indices) throw std::invalid_argument("IonicReferenceUnavailable: solvent-only reference");
    if (!(reference.solver_options == model.solver_options()))
        throw std::invalid_argument("electrolyte reference solver options do not match the model");
    double solute_ratio = 0.0;
    auto composition = formula_composition(
        solvent_fractions,
        model.parameters(),
        ion_indices,
        stoichiometry,
        molality, finite_solute_molalities, &solute_ratio
    );
    if (finite_solute_molalities) {
        if (reference.model != &model || reference.temperature != temperature || reference.pressure != pressure
            || reference.solvent_fractions != solvent_fractions)
            throw std::invalid_argument("electrolyte reference model, temperature, pressure or solvent mismatch");
    }
    ElectrolyteActivity result;
    result.solver_options = model.solver_options();
    result.state = properties::close_density(model, temperature, pressure, composition, DensityBranch::Liquid, finite_density);
    if (const auto cause = properties::vapor_side(result.state)) throw std::runtime_error(*cause);
    const auto& state = result.state;
    const auto log_fugacity = ion_log_fugacity(state, ion_indices);
    // EqID: single_ion_activity_coefficient
    const double first = activity_log(log_fugacity[0], reference.log_fugacity[0]);
    const double second = activity_log(log_fugacity[1], reference.log_fugacity[1]);
    const double coefficient_sum = static_cast<double>(stoichiometry[0]) + stoichiometry[1];
    // EqID: mean_ionic_activity_coefficient
    const double mean_log_x = (
        stoichiometry[0] * first + stoichiometry[1] * second
    ) / coefficient_sum;
    // EqID: molality_activity_conversion
    const double mean_log_m = activity_log(mean_log_x,
        std::log(1.0 + solute_ratio));
    result.single_ion = {std::exp(first), std::exp(second)};
    result.mean_ionic_mole_fraction_scale = std::exp(mean_log_x);
    result.mean_ionic_log_molality_scale = mean_log_m;
    result.mean_ionic_molality_scale = std::exp(mean_log_m);
    if (reference.pure_solvent_state && !finite_solute_molalities && molality < policy::OSMOTIC_MOLALITY_FLOOR.value)
        result.osmotic_coefficient = core::Stop::converged_to_floor;
    else if (reference.pure_solvent_state) {
        auto water = solvent_activity(model, temperature, pressure, composition, reference, std::move(result.state),
            solute_ratio);
        result.state = std::move(water.state);
        if (!water.solvent_activity) throw std::runtime_error(water.cause);
        result.osmotic_coefficient = water.osmotic_coefficient;
        for (const auto& column : water.osmotic_parameter_derivative)
            result.osmotic_parameter_derivative.push_back(require_component_property(column));
    }
    const std::size_t columns = model.active_parameters().size();
    for (std::size_t a = 0; a < columns; ++a) {
        // EqID: electrolyte_activity_parameter_derivative
        const auto& ref = reference.parameter_derivative.at(a);
        const auto g = [&](std::size_t i) { return require_component_property(state.log_fugacity_parameter_derivative.at(i * columns + a)); };
        auto action = !ref[0].action ? ref[0] : ref[1];
        if (ref[0].action && ref[1].action) {
            action.action = result.mean_ionic_molality_scale * (stoichiometry[0] * (g(ion_indices[0]) - *ref[0].action)
                + stoichiometry[1] * (g(ion_indices[1]) - *ref[1].action)) / coefficient_sum;
            if (!std::isfinite(*action.action)) { action.action.reset(); action.status = PropertyActionStatus::Nonfinite; }
        }
        result.mean_ionic_parameter_derivative.push_back(action);
    }
    return result;
}

ElectrolyteInputAction electrolyte_input_action(const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry,
    double molality, const InfiniteDilutionReference& reference, const ElectrolyteActivity& activity,
    const ElectrolyteInputDirection& direction) {
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
        const PressureActionDirection d{direction.temperature, direction.pressure, {}, {}};
        for (std::size_t i = 0; i < 2; ++i) {
            std::vector<double> insertion(count, 0.0); insertion[ions[i]] = 1.0;
            const auto series = model.reference_taylor(temperature, 1.0 / reference.solvent_state->molar_density,
                solvent, insertion, d, 1, true);
            result.diagnostics.push_back(series.diagnostics);
            if (series.status != PropertyActionStatus::Available) {
                result.status = series.status; result.message = "reference input action unavailable"; return result;
            }
            result.reference_column[i] = series.contraction[1] - direction.pressure / pressure + direction.temperature / temperature;
        }
        if (!reference.pure_solvent_state) { result.message = "pure solvent input action unavailable"; return result; }
        std::vector<double> insertion(count, 0.0); insertion[sidx] = 1.0;
        const auto pure = model.reference_taylor(temperature, 1.0 / reference.pure_solvent_state->molar_density,
            solvent, insertion, d, 1, true);
        result.diagnostics.push_back(pure.diagnostics);
        if (pure.status != PropertyActionStatus::Available) {
            result.status = pure.status; result.message = "pure solvent input action unavailable"; return result;
        }
        pure_column = pure.contraction[1] - direction.pressure / pressure + direction.temperature / temperature;
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
