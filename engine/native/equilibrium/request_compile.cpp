#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"
#include "core/trace.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_tape.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

std::vector<Layout> layouts(const CompiledSystem& compiled) {
    std::vector<Layout> result;
    std::size_t offset = 0;
    for (const auto& phase : compiled.problem.phases) {
        Layout item;
        item.x.resize(phase.support.size());
        for (std::size_t i = 0; i < phase.support.size(); ++i) item.x[i] = offset++;
        if (phase.amount.kind == AmountKind::Free) item.amount = offset++;
        item.volume = offset++;
        result.push_back(std::move(item));
    }
    return result;
}

std::vector<std::size_t> complete_support(
    std::size_t components, const std::vector<std::size_t>& support
) {
    if (support.empty()) {
        std::vector<std::size_t> all(components);
        std::iota(all.begin(), all.end(), 0);
        return all;
    }
    std::vector<std::size_t> result = support;
    std::sort(result.begin(), result.end());
    if (std::adjacent_find(result.begin(), result.end()) != result.end()) {
        throw std::invalid_argument("duplicate support index");
    }
    if (result.back() >= components) throw std::invalid_argument("support index out of range");
    return result;
}

double phase_amount(const PhaseSpec& phase, const std::vector<double>& feed) {
    if (phase.amount.kind == AmountKind::Free) return 1.0;
    return phase.amount.fraction * vector_norm(feed);
}

std::size_t count_rows(const CompiledSystem& compiled, RowKind kind) {
    std::size_t result = 0;
    for (const auto& row : compiled.rows) result += row.kind == kind;
    return result;
}

void validate_intensive_seed(const FreeValue& value, const char* name) {
    if (!value.free && !value.is_positive()) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
    if (value.free && !value.absent && !value.is_positive()) {
        throw std::invalid_argument(std::string("invalid ") + name + " guess");
    }
}

PackingTerms packing_terms(const HelmholtzModel& model, double temperature,
                           const std::vector<double>& composition) {
    const auto& parameters = model.parameters();
    const auto sigma = segment_sigma(temperature, parameters);
    const auto diameter = segment_diameter(temperature, parameters, sigma);
    std::vector<double> sigma_t(sigma.size(), 0.0);
    for (std::size_t term = 0; term < parameters.sigma_term_component.size(); ++term) {
        const auto component = static_cast<std::size_t>(parameters.sigma_term_component[term]);
        sigma_t[component] += parameters.sigma_term_amplitude[term]
            * parameters.sigma_term_exponent[term] * std::exp(
                parameters.sigma_term_exponent[term] * temperature);
    }
    PackingTerms result;
    result.component.resize(composition.size());
    const double scale = (PI / 6.0) * NUMBER_DENSITY_PER_ANGSTROM3;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        double diameter_t = 0.0;
        if (parameters.charge_number[i] == 0) {
            const double exponential = std::exp(-3.0
                * parameters.dispersion_energy_over_k[i] / temperature);
            diameter_t = sigma_t[i] * (1.0 - 0.12 * exponential)
                - sigma[i] * 0.36 * parameters.dispersion_energy_over_k[i]
                    * exponential / (temperature * temperature);
        }
        result.component[i] = scale * composition[i] * parameters.segment_count[i]
            * diameter[i] * diameter[i] * diameter[i];
        result.coefficient += result.component[i];
        result.temperature_derivative += scale * composition[i]
            * parameters.segment_count[i] * 3.0 * diameter[i] * diameter[i] * diameter_t;
    }
    return result;
}

// Compare supplied species thermochemistry with every reaction target on the density
// basis: [t, T t', T^2 t''] against [g, g - h, -T cp] at the same J/mol scale.
void check_reference_consistency(
    const ThermochemistryRecord& record, double temperature,
    const std::vector<ReactionSpec>& reactions
) {
    const auto targets = reaction_target_series(temperature, reactions, 2);
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        const auto species = reaction_thermochemistry_density(
            record, reactions[row].stoichiometry, temperature);
        if (!species.has_value()) continue;
        const double t0 = targets[0][row];
        const double t1 = temperature * targets[1][row];
        const double t2 = temperature * temperature * targets[2][row];
        const double d0 = t0 - species->gibbs;
        const double d1 = t1 - (species->gibbs - species->enthalpy);
        const double d2 = t2 + temperature * species->heat_capacity;
        const double scale = std::sqrt(std::numeric_limits<double>::epsilon())
            * std::max(1.0, std::sqrt(t0 * t0 + t1 * t1 + t2 * t2));
        if (std::sqrt(d0 * d0 + d1 * d1 + d2 * d2) > scale) {
            throw std::invalid_argument("reference inconsistency in reaction "
                                        + reaction_name(reactions[row], row));
        }
    }
}

bool has_reaction_reference(const EquilibriumProblem& problem) {
    return std::any_of(problem.reactions.begin(), problem.reactions.end(),
        [](const ReactionSpec& reaction) { return reaction.reference.has_value(); });
}

std::vector<ReactionSpec> central_reference_reactions(
    const HelmholtzModel &model, const EquilibriumProblem &request, double temperature, double pressure,
    std::vector<NeutralReferenceReactionDiagnostic> &diagnostics, bool pressure_jets,
    const std::vector<NeutralReferenceReactionDiagnostic> *fixed_pressure = nullptr,
    std::size_t temperature_order = 0, bool parameter_jets = false, std::size_t reference_order = 0) {
    if (!has_reaction_reference(request)) return request.reactions;
    if (!request.source_reference.has_value()) {
        throw std::invalid_argument("reference_unavailable: source reference declaration is missing");
    }
    if (!request.standard_gibbs_j_per_mol.empty()) {
        throw std::invalid_argument(
            "reference_unavailable: supplied standard Gibbs conflicts with source reference");
    }
    std::vector<ReactionSpec> result = request.reactions;
    std::vector<NeutralReferenceReactionDiagnostic> reused;
    if (fixed_pressure != nullptr && temperature_order == 0 && !parameter_jets) {
        // Fixed-pressure references depend only on the compiled model, temperature,
        // declaration and stoichiometry; rebuild only system-pressure traces.
        for (const auto& diagnostic : *fixed_pressure)
            if (request.reactions[diagnostic.reaction_index].reference->reference_pressure_pa)
                reused.push_back(diagnostic);
    }
    diagnostics.clear();
    diagnostics.reserve(request.reactions.size());
    std::vector<std::size_t> reference_rows;
    std::vector<std::vector<double>> stoichiometries;
    std::vector<ReactionReference> references;
    for (std::size_t row = 0; row < request.reactions.size(); ++row) {
        const auto& reaction = request.reactions[row];
        if (!reaction.reference.has_value()) continue;
        if (std::any_of(reused.begin(), reused.end(),
                [row](const auto& item) { return item.reaction_index == row; })) continue;
        if (reaction.affinity_j_per_mol.has_value()
            || (!reaction.correlation.has_value() && !reaction.equilibrium_constant.has_value())) {
            throw std::invalid_argument(
                "reference_unavailable: reaction " + reaction_name(reaction, row)
                + " needs a scalar K or correlation source");
        }
        reference_rows.push_back(row);
        stoichiometries.push_back(reaction.stoichiometry);
        references.push_back(*reaction.reference);
    }
    std::optional<double> solvent_anchor;
    if (fixed_pressure != nullptr && !fixed_pressure->empty() && fixed_pressure->front().solvent_density > 0.0)
        solvent_anchor = fixed_pressure->front().solvent_density;
    const auto reference_results = reaction_references(
        model, temperature, pressure, *request.source_reference, stoichiometries, references, pressure_jets,
        solvent_anchor, temperature_order, parameter_jets, reference_order);
    std::vector<NeutralReferenceReactionDiagnostic> resolved = reused;
    for (std::size_t reference_row = 0; reference_row < reference_rows.size(); ++reference_row) {
        resolved.push_back(reference_results[reference_row]);
        resolved.back().reaction_index = reference_rows[reference_row];
    }
    std::sort(resolved.begin(), resolved.end(), [](const auto& a, const auto& b) {
        return a.reaction_index < b.reaction_index;
    });
    for (const auto& reference : resolved) {
        const std::size_t row = reference.reaction_index;
        const auto& reaction = request.reactions[row];
        diagnostics.push_back(reference);
        if (reference.status != NeutralReferenceStatus::Available
            || !reference.terminal_contraction.has_value()) {
            throw std::invalid_argument(
                "reference_unavailable: reaction " + reaction_name(reaction, row)
                + (reference.message.empty() ? "" : " " + reference.message));
        }
        double ln_k = 0.0;
        if (reaction.correlation.has_value()) {
            try {
                ln_k = evaluate_reaction_properties(
                    *reaction.correlation, temperature).ln_k;
            } catch (const std::exception& error) {
                throw std::invalid_argument(
                    "reference_unavailable: reaction " + reaction_name(reaction, row)
                    + " " + error.what());
            }
        } else {
            ln_k = std::log(*reaction.equilibrium_constant);
        }
        double scale = 0.0;
        if (reaction.reference->source_basis
                == ReferenceSourceBasis::CommonMolalityInfiniteDilution) {
            try {
                scale = reaction_reference_log_scale(
                    model, *request.source_reference, reaction.stoichiometry);
            } catch (const std::exception& error) {
                throw std::invalid_argument(
                    "reference_unavailable: reaction " + reaction_name(reaction, row)
                    + " " + error.what());
            }
        }
        const double ln_k_density = ln_k + scale + *reference.terminal_contraction;
        const double target = -GAS_CONSTANT_J_PER_MOL_K * temperature * ln_k_density;
        if (!std::isfinite(target)) {
            throw std::invalid_argument(
                "reference_unavailable: reaction " + reaction_name(reaction, row)
                + " has a nonfinite density-basis target");
        }
        result[row].equilibrium_constant.reset();
        result[row].correlation.reset();
        result[row].affinity_j_per_mol = target;
        result[row].reference.reset();
    }
    return result;
}

// The same minimum-norm reaction-to-species map for dimensionless standard-term slopes.
std::vector<double> species_standard_slopes(double temperature, std::size_t components,
    std::vector<ReactionSpec> reactions, const std::vector<double>& values) {
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        reactions[row].equilibrium_constant.reset();
        reactions[row].correlation.reset();
        reactions[row].reference.reset();
        reactions[row].affinity_j_per_mol = values[row];
    }
    return standard_state(temperature, components, reactions, {});
}

StandardTerm standard_term(
    const CompiledSystem& compiled, double temperature, double pressure, std::size_t order,
    bool pressure_jets, bool temperature_jets = false, const std::vector<double>& active_control = {}
) {
    const auto& problem = compiled.problem;
    const std::size_t components = compiled.model->parameters().component_count();
    if (!compiled.reference_diagnostics.empty()) {
        if (compiled.model->thermochemistry().has_value()) {
            std::vector<ReactionSpec> ordinary_reactions;
            ordinary_reactions.reserve(problem.reactions.size());
            for (const auto& reaction : problem.reactions)
                if (!reaction.reference.has_value()) ordinary_reactions.push_back(reaction);
            if (!ordinary_reactions.empty()) {
                check_reference_consistency(*compiled.model->thermochemistry(), temperature,
                                            ordinary_reactions);
            }
        }
        auto diagnostics = compiled.reference_diagnostics;
        auto reactions = problem.reactions;
        auto gibbs = compiled.standard_gibbs;
        const bool parameter_jets = std::any_of(active_control.begin(), active_control.end(),
            [](double value) { return value != 0.0; });
        // Fixed-T/P references do not change with phase coordinates or feed directions.
        const std::size_t thermal = temperature_jets ? order : 0;
        if (problem.pressure.free || pressure_jets || temperature_jets || parameter_jets) {
            reactions = central_reference_reactions(
                *compiled.model, problem, temperature, pressure, diagnostics, pressure_jets,
                &compiled.reference_diagnostics,
                thermal, parameter_jets,
                parameter_jets && order > 1 && (temperature_jets || pressure_jets) ? 2 : 0);
            gibbs = standard_state(temperature, components, reactions, {});
        }
        StandardTerm result;
        result.derivative.assign(order + 1, std::vector<double>(components, 0.0));
        result.pressure_derivative.assign(components, 0.0);
        result.pressure_second_derivative.assign(components, 0.0);
        result.reference_diagnostics = diagnostics;
        const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
        for (std::size_t species = 0; species < components; ++species)
            result.derivative[0][species] = gibbs[species] / rt;
        // One slope of s_r = g_r/(RT) per row, mapped to species by the same linear minimum-norm
        // owner as s. A referenced row has s_r = -ln K_rho, so its slope is -(source slope + slope
        // of C); an unadmitted slope of C refuses the direction.
        const auto species_slope = [&](std::vector<double> values, const auto& slope_of, const std::string& name,
                                       std::size_t source_order, std::size_t a = 0, std::size_t b = 0,
                                       std::optional<std::size_t> parameter = std::nullopt) {
            if (source_order) a = source_order;
            for (const auto& diagnostic : diagnostics) {
                const std::size_t row = diagnostic.reaction_index;
                const std::optional<double> slope = slope_of(diagnostic);
                if (!slope) {
                    const auto cause = diagnostic.derivative_cause(a, b, parameter);
                    throw std::invalid_argument("reference_unavailable: reaction " + reaction_name(problem.reactions[row], row)
                        + " C_" + std::to_string(a) + std::to_string(b)
                        + (parameter ? " coordinate " + std::to_string(*parameter) : "") + " (" + name + "): " + cause);
                }
                const auto& source = problem.reactions[row].correlation;
                const auto law = source ? evaluate_reaction_properties(*source, temperature) : ReactionProperties{};
                values[row] = -*slope - (source_order == 1 ? law.derivative : source_order == 2 ? law.second_derivative : 0.0);
            }
            return species_standard_slopes(temperature, components, reactions, values);
        };
        using Diagnostic = NeutralReferenceReactionDiagnostic;
        const std::vector<double> zero(reactions.size(), 0.0);
        // EqIDs: standard_state_transfer_temperature_(second_)derivative. Unreferenced rows keep
        // s' = (t' - t/T)/(RT) and s'' = t''/(RT) - 2 s'/T.
        std::vector<double> first(reactions.size(), 0.0), second(reactions.size(), 0.0);
        const auto targets = thermal > 0 ? reaction_target_series(temperature, reactions, thermal)
                                         : std::vector<std::vector<double>>{};
        for (std::size_t row = 0; thermal > 0 && row < reactions.size(); ++row) {
            first[row] = (targets[1][row] - targets[0][row] / temperature) / rt;
            if (thermal > 1) second[row] = targets[2][row] / rt - 2.0 * first[row] / temperature;
        }
        if (thermal > 0) result.derivative[1] = species_slope(first, [](const Diagnostic& d) {
            return d.terminal_contraction_temperature_derivative; }, "temperature derivative", 1);
        if (thermal > 1) result.derivative[2] = species_slope(second, [](const Diagnostic& d) {
            return d.terminal_contraction_temperature_second_derivative; }, "second temperature derivative", 2);
        // ln K_src and the molality scale are parameter free: ds_r/dtheta_a = -dC_r/dtheta_a
        // (EqID: standard_state_transfer_parameter_derivative).
        if (parameter_jets) {
            result.parameter_derivative.assign(active_control.size(), std::vector<double>(components, 0.0));
            if (order > 1 && temperature_jets)
                result.temperature_parameter_derivative.assign(active_control.size(), std::vector<double>(components, 0.0));
            if (order > 1 && pressure_jets)
                result.pressure_parameter_derivative.assign(active_control.size(), std::vector<double>(components, 0.0));
        }
        for (std::size_t a = 0; parameter_jets && a < active_control.size(); ++a) {
            if (active_control[a] == 0.0) continue;
            result.parameter_derivative[a] = species_slope(zero, [a](const Diagnostic& d) {
                return a < d.terminal_contraction_parameter_derivative.size()
                    ? d.terminal_contraction_parameter_derivative[a] : std::nullopt; }, "parameter derivative", 0, 0, 0, a);
            if (order > 1) {
                if (temperature_jets) result.temperature_parameter_derivative[a] = species_slope(
                    zero,
                    [a](const Diagnostic &d) {
                        return a < d.terminal_contraction_temperature_parameter_derivative.size()
                                   ? d.terminal_contraction_temperature_parameter_derivative[a]
                                   : std::nullopt;
                    },
                    "T parameter derivative", 0, 1, 0, a);
                if (pressure_jets)
                    result.pressure_parameter_derivative[a] = species_slope(
                        zero,
                        [a](const Diagnostic &d) -> std::optional<double> {
                            for (const auto &column : d.derivatives)
                                if (column.temperature_order == 0 && column.pressure_order == 1 &&
                                    column.parameter == a)
                                    return column.value;
                            return std::nullopt;
                        },
                        "P parameter derivative", 0, 0, 1, a);
            }
        }
        // EqIDs: neutral_reaction_reference_pressure_jet, standard_state_transfer_temperature_pressure_derivative.
        if (pressure_jets) {
            result.pressure_derivative = species_slope(zero, [](const Diagnostic& d) {
                return d.terminal_contraction_pressure_derivative; }, "pressure derivative", 0, 0, 1);
            result.pressure_second_derivative = species_slope(zero, [](const Diagnostic& d) {
                return d.terminal_contraction_pressure_second_derivative; }, "second pressure derivative", 0, 0, 2);
            if (thermal > 1) result.temperature_pressure_derivative = species_slope(zero, [](const Diagnostic& d) {
                return d.terminal_contraction_temperature_pressure_derivative; }, "temperature-pressure derivative", 0, 1, 1);
        }
        return result;
    }
    const auto gibbs = standard_state_series(
        temperature, components, problem.reactions, problem.standard_gibbs_j_per_mol, order);
    if (compiled.model->thermochemistry().has_value() && !problem.reactions.empty()) {
        check_reference_consistency(*compiled.model->thermochemistry(), temperature,
                                    problem.reactions);
    }
    // s = g u / R with u = 1/T; Leibniz over u' = -1/T^2, u'' = 2/T^3, u''' = -6/T^4.
    const double T = temperature;
    const double u[4] = {1.0 / T, -1.0 / (T * T), 2.0 / (T * T * T), -6.0 / (T * T * T * T)};
    const double binomial[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
    StandardTerm result;
    result.derivative.assign(order + 1, std::vector<double>(components, 0.0));
    result.pressure_derivative.assign(components, 0.0);
    result.pressure_second_derivative.assign(components, 0.0);
    for (std::size_t k = 0; k <= order; ++k) {
        for (std::size_t species = 0; species < components; ++species) {
            double value = 0.0;
            for (std::size_t j = 0; j <= k; ++j) value += binomial[k][j] * gibbs[j][species] * u[k - j];
            result.derivative[k][species] = value / GAS_CONSTANT_J_PER_MOL_K;
        }
    }
    return result;
}

// PH/PS consumes complete records at every callback temperature.
void admit_calorics(const CompiledSystem& compiled, double temperature) {
    if (!compiled.has_caloric_condition()) return;
    const auto& p = compiled.problem;
    if (p.enthalpy_j && p.entropy_j_per_k)
        throw std::invalid_argument("one extensive caloric condition required");
    if (!p.temperature.free || p.pressure.free || compiled.pinned_count
        || !p.topology_declared || p.discovery)
        throw std::invalid_argument("PH/PS requires declared free amounts, free T and fixed P");
    if (compiled.has_ions || p.source_reference || has_reaction_reference(p))
        throw std::invalid_argument("charged/source-referenced PH/PS unavailable");
    if (!std::isfinite(p.enthalpy_j ? *p.enthalpy_j : *p.entropy_j_per_k))
        throw std::invalid_argument("nonfinite caloric target");
    const auto& record = compiled.model->thermochemistry();
    if (!record) throw std::invalid_argument("missing thermochemistry for caloric condition");
    for (const auto& phase : p.phases)
        for (auto i : phase.support) {
            const auto& intervals = record->component_intervals[i];
            if (intervals.empty()) throw std::invalid_argument("missing caloric record: " + std::to_string(i));
            if (std::none_of(intervals.begin(), intervals.end(), [&](const auto& v) {
                    return ideal_interval_contains(v, temperature); }))
                throw std::invalid_argument("caloric record outside interval: species " + std::to_string(i)
                    + " at " + std::to_string(temperature) + " K");
        }
}

// Critical and envelope cells: one condition, neutral, nonreactive, T and P free (issue #224 B).
void admit_conditions(const CompiledSystem& compiled) {
    const auto& p = compiled.problem;
    if (!p.critical && !p.envelope) return;
    if ((p.critical && p.envelope) || compiled.has_caloric_condition() || !p.temperature.free || !p.pressure.free
        || p.discovery || compiled.has_ions || compiled.has_reaction || p.source_reference)
        throw std::invalid_argument("critical and envelope cells require one neutral nonreactive condition with free T and P");
    const auto& support = p.phases[0].support;
    if (p.critical) {
        const auto& c = *p.critical;
        if (p.phases.size() != 1 || !p.phases[0].is_free() || support.size() != compiled.n_feed.size()
            || !(std::isfinite(c.density_guess_mol_per_m3) && c.density_guess_mol_per_m3 > 0.0)
            || (!c.mode_guess.empty() && c.mode_guess.size() != support.size())
            || std::any_of(compiled.feed_composition.begin(), compiled.feed_composition.end(), [](double z) { return !(z > 0.0); }))
            throw std::invalid_argument("critical cell requires one free phase on a positive full support and a density seed");
        return;
    }
    const auto& e = *p.envelope;
    if (p.phases.size() != 2 || compiled.pinned_count != 1 || !std::isfinite(e.target) || p.phases[1].support != support
        || std::none_of(p.phases.begin(), p.phases.end(), [](const PhaseSpec& x) { return x.is_pinned() && x.amount.fraction == 0.0; })
        || (e.kind == EnvelopeCoordinate::LogK && std::find(support.begin(), support.end(), e.component) == support.end()))
        throw std::invalid_argument("envelope requires a free feed phase and a zero-pinned incipient phase on one support");
}

// Caloric, critical and envelope rows: analytic in record_geometry, completed by row_series.
bool condition_row(RowKind kind) { return kind >= RowKind::Enthalpy; }

template<class Scalar>
IdealMixtureValue<Scalar> caloric_ideal(const CompiledSystem& compiled, std::size_t phase,
    const Scalar& T, const Scalar& volume, const std::vector<Scalar>& x, double center) {
    const auto& record = *compiled.model->thermochemistry();
    std::vector<IdealCorrelation> correlations;
    std::vector<Scalar> fractions;
    for (auto i : compiled.problem.phases[phase].support) {
        const auto& intervals = record.component_intervals[i];
        const auto selected = std::find_if(intervals.begin(), intervals.end(), [&](const auto& v) {
            return ideal_interval_contains(v, center); });
        correlations.push_back(selected->correlation);
        fractions.push_back(x[i]);
    }
    return evaluate_ideal_mixture(correlations, record.standard_pressure, T, Scalar(1.0) / volume,
        fractions, std::vector<bool>(fractions.size(), true));
}

std::vector<double> caloric_hessian(const CompiledSystem&, const PhasePoint&, double);
std::vector<double> critical_hessian(const CompiledSystem&, const PhasePoint&, const std::vector<double>&);

std::vector<double> pinned_hessian(
    const CompiledSystem& compiled, const PhasePoint& prepared,
    const std::vector<double>& multipliers
) {
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    // The pinned objective is constant; only optimizer constraint weights matter.
    if (multipliers.empty()) return result;
    if (multipliers.size() != compiled.rows.size())
        throw std::invalid_argument("constraint multiplier count");
    const auto& problem = compiled.problem;
    const auto& model = *compiled.model;
    const auto& parameters = model.parameters();
    const std::size_t components = parameters.component_count();
    const std::size_t physical = 2 + components;
    const auto& shape = prepared.shape;
    const std::size_t temperature_index = compiled.primal_size - compiled.free_intensive_count;
    const double temperature = prepared.temperature;
    const auto& standard = prepared.standard;
    const auto& charges = parameters.charge_number;
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        std::vector<double> chemical(components, 0.0);
        double pressure_weight = 0.0;
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            const auto& descriptor = compiled.rows[row];
            if (descriptor.phase != phase) continue;
            if (descriptor.kind == RowKind::Chemical)
                chemical[descriptor.index] += multipliers[row];
            if (descriptor.kind == RowKind::Pressure) pressure_weight += multipliers[row];
        }
        using AD = CppAD::AD<double>;
        std::vector<AD> q(point.begin(), point.end());
        CppAD::Independent(q);
        const AD T = problem.temperature.free ? q[temperature_index] : AD(temperature);
        const AD P = problem.pressure.free ? q[compiled.primal_size - 1]
                                          : AD(problem.pressure.value);
        const AD delta_P = P - AD(prepared.pressure);
        std::vector<AD> x(components, AD(0.0));
        for (std::size_t local = 0; local < shape[phase].x.size(); ++local)
            x[problem.phases[phase].support[local]] = CppAD::exp(q[shape[phase].x[local]]);
        const auto diameter = segment_diameter(T, parameters, segment_sigma(T, parameters));
        const AD volume = packing_fraction_from_geometry(AD(1.0), x, parameters, diameter)
            / CppAD::exp(q[shape[phase].volume]);
        const AD amount = shape[phase].amount.has_value() ? q[*shape[phase].amount]
            : AD(phase_amount(problem.phases[phase], compiled.n_feed));
        AD simple = pressure_weight * (AD(1.0) - P * volume / (GAS_CONSTANT_J_PER_MOL_K * T));
        for (std::size_t local = 0; local < shape[phase].x.size(); ++local) {
            const std::size_t species = problem.phases[phase].support[local];
            const AD delta_T = T - temperature;
            // Exact second-order Taylor term of the standard state in the free temperature.
            const AD standard_value = standard.derivative[0][species]
                + standard.derivative[1][species] * delta_T
                + 0.5 * standard.derivative[2][species] * delta_T * delta_T
                + standard.pressure_derivative[species] * delta_P
                + 0.5 * standard.pressure_second_derivative[species]
                    * delta_P * delta_P;
            simple += chemical[species] * (q[shape[phase].x[local]]
                - CppAD::log(volume) + standard_value);
        }
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            const auto& descriptor = compiled.rows[row];
            if (descriptor.kind == RowKind::Moiety) {
                for (std::size_t species = 0; species < components; ++species)
                    simple += multipliers[row] * amount * compiled.B[descriptor.index][species]
                        * x[species] / compiled.feed_total;
            } else if (descriptor.phase == phase && descriptor.kind == RowKind::Electro) {
                for (std::size_t species = 0; species < components; ++species)
                    simple += multipliers[row] * charges[species] * x[species];
            } else if (descriptor.phase == phase && descriptor.kind == RowKind::Normalization) {
                for (const auto& fraction : x) simple += multipliers[row] * fraction;
            }
        }
        std::vector<AD> output = {T, volume};
        output.insert(output.end(), x.begin(), x.end());
        output.push_back(simple);
        CppAD::ADFun<double> geometry(q, output);
        const auto values = geometry.Forward(0, point);
        const std::vector<double> amounts(values.begin() + 2, values.begin() + physical);
        const auto& derivative = prepared.derivative[phase];
        const std::size_t width = derivative.gradient.size();
        std::vector<double> contraction(width, 0.0);
        contraction[1] = -pressure_weight * values[1];
        std::copy(chemical.begin(), chemical.end(), contraction.begin() + 2);
        std::vector<double> map_weights(physical + 1, 0.0);
        map_weights.back() = 1.0;
        for (std::size_t i = 0; i < physical; ++i) {
            for (std::size_t k = 0; k < width; ++k)
                map_weights[i] += derivative.hessian[i * width + k] * contraction[k];
        }
        map_weights[1] -= pressure_weight * derivative.gradient[1];
        const auto mapped = geometry.Hessian(point, map_weights);
        const auto jacobian = geometry.Jacobian(point);
        std::vector<std::vector<double>> directions(n, std::vector<double>(width, 0.0));
        std::vector<double> volume_hessian_action(n, 0.0);
        for (std::size_t column = 0; column < n; ++column) {
            for (std::size_t i = 0; i < physical; ++i) {
                directions[column][i] = jacobian[i * n + column];
                volume_hessian_action[column] += derivative.hessian[width + i]
                    * directions[column][i];
            }
        }
        // Trilinear contraction: a coordinate that does not move this phase contributes exactly zero.
        // ponytail: quadratic calls over moving pairs; contract a physical basis if cost requires it.
        std::vector<bool> moves(n);
        for (std::size_t c = 0; c < n; ++c)
            moves[c] = std::any_of(directions[c].begin(), directions[c].end(), [](double v) { return v != 0.0; });
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t column = 0; column <= row; ++column) {
                double action = 0.0;
                if (moves[row] && moves[column]) {
                    const auto energy = model.residual_energy_action(temperature, values[1], amounts,
                        {contraction, directions[row], directions[column]});
                    if (!energy.action.has_value()) throw std::runtime_error(std::string("third directional energy "
                        "contraction unavailable: ") + solved_state_status_name(solved_state_status(energy.status)));
                    action = *energy.action;
                }
                const double value = mapped[row * n + column] + action
                    - pressure_weight * (directions[row][1] * volume_hessian_action[column]
                        + directions[column][1] * volume_hessian_action[row]);
                result[row * n + column] += value;
                if (row != column) result[column * n + row] += value;
            }
        }
    }
    if (compiled.has_caloric_condition()) {
        const auto caloric = caloric_hessian(compiled, prepared, multipliers.back());
        for (std::size_t i = 0; i < result.size(); ++i) result[i] += caloric[i];
    }
    if (problem.critical) {
        const auto critical = critical_hessian(compiled, prepared, multipliers);
        for (std::size_t i = 0; i < result.size(); ++i) result[i] += critical[i];
    }
    // EqID envelope_local_condition: d2 ln q/dq2 = -1/q^2 for q = T or P; a ln K difference is affine.
    if (problem.envelope && problem.envelope->kind != EnvelopeCoordinate::LogK) {
        const std::size_t i = problem.envelope->kind == EnvelopeCoordinate::LogTemperature ? temperature_index
            : compiled.primal_size - 1;
        result[i * n + i] -= multipliers.back() / (point[i] * point[i]);
    }
    return result;
}

// The geometry map and every analytic row term as one CppAD function of (y, controls) with
// controls [T, P, n_feed, active parameters]: outputs per phase [T, v, x (C), active (A),
// N, eta], then the analytic part of each row (the residual gradient terms are added by the
// model's contractions in row_series).
CppAD::ADFun<double> record_geometry(
    const CompiledSystem& compiled, const PhasePoint& prepared, const StandardTerm& standard
) {
    using AD = CppAD::AD<double>;
    const auto& problem = compiled.problem;
    const auto& model = *compiled.model;
    const auto& parameters = model.parameters();
    const std::size_t components = parameters.component_count();
    const std::size_t n = prepared.point.size();
    std::vector<AD> u(prepared.point.begin(), prepared.point.end());
    u.push_back(prepared.temperature);
    u.push_back(prepared.pressure);
    for (const double amount : compiled.n_feed) u.push_back(amount);
    const auto theta_center = model.active_parameter_values();
    for (const double value : theta_center) u.push_back(value);
    if (compiled.has_caloric_condition()) {
        u.push_back(problem.enthalpy_j.value_or(0.0));
        u.push_back(problem.entropy_j_per_k.value_or(0.0));
    }
    CppAD::Independent(u);
    const AD T = problem.temperature.free ? u[compiled.primal_size - compiled.free_intensive_count] : u[n];
    const AD P = problem.pressure.free ? u[compiled.primal_size - 1] : u[n + 1];
    AD feed_total = AD(0.0);
    for (std::size_t i = 0; i < components; ++i) feed_total += u[n + 2 + i];
    const std::vector<AD> theta(u.begin() + static_cast<std::ptrdiff_t>(n + 2 + components), u.begin() + static_cast<std::ptrdiff_t>(n + 2 + components + theta_center.size()));
    const auto tape_parameters = detail::make_tape_parameters(parameters, model.active_parameters(), theta);
    const auto sigma = segment_sigma(T, tape_parameters);
    const auto diameter = segment_diameter(T, tape_parameters, sigma);
    const auto& charges = parameters.charge_number;
    std::vector<std::optional<std::size_t>> phi(problem.phases.size());
    std::size_t phi_offset = 0;
    for (const auto& row : compiled.rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = compiled.primal_size + compiled.rank + phi_offset++;
    std::vector<AD> output;
    std::vector<AD> row_values(compiled.rows.size(), AD(0.0));
    const AD delta_T = T - prepared.temperature;
    const AD delta_P = P - prepared.pressure;
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        const auto& layout = prepared.shape[phase];
        std::vector<AD> x(components, AD(0.0));
        for (std::size_t local = 0; local < spec.support.size(); ++local)
            x[spec.support[local]] = CppAD::exp(u[layout.x[local]]);
        const AD eta = CppAD::exp(u[layout.volume]);
        const AD volume = packing_fraction_from_geometry(AD(1.0), x, tape_parameters, diameter) / eta;
        const AD amount = layout.amount.has_value() ? u[*layout.amount] : AD(spec.amount.fraction) * feed_total;
        output.push_back(T);
        output.push_back(volume);
        output.insert(output.end(), x.begin(), x.end());
        output.insert(output.end(), theta.begin(), theta.end());
        output.push_back(amount);
        output.push_back(eta);
        if (compiled.has_caloric_condition()) {
            const bool h = problem.enthalpy_j.has_value();
            const AD denominator = feed_total * GAS_CONSTANT_J_PER_MOL_K * (h ? T : AD(1.0));
            const auto ideal = caloric_ideal(compiled, phase, T, volume, x, prepared.temperature);
            row_values.back() += amount * (h ? ideal.enthalpy : ideal.entropy) / denominator;
            const AD weight = amount * GAS_CONSTANT_J_PER_MOL_K / denominator;
            output.push_back(h ? AD(0.0) : -weight);
            output.push_back(h ? -weight * T * T : -weight * T);
            output.push_back(h ? -weight * T * volume : AD(0.0));
        }
        if (problem.critical) {
            // EqID critical_equations at n = x: D = diag(sqrt(z)) of the moving feed, w = D u, ideal
            // Q_ii = 1/n_i and F_iii = -1/n_i^2; row_series adds the residual F_nn w and F_nnn[w,w,w].
            const std::size_t mode = compiled.rows.size() - components - 2;
            row_values[mode + components] = AD(-1.0);
            for (std::size_t i = 0; i < components; ++i) {
                const AD u_i = u[compiled.mode_offset() + i], d = CppAD::sqrt(u[n + 2 + i] / feed_total), w = d * u_i;
                output.push_back(d);
                output.push_back(w);
                row_values[mode + i] = d * w / x[i];
                row_values[mode + components] += u_i * u_i;
                row_values[mode + components + 1] -= w * w * w / (x[i] * x[i]);
            }
        }
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            const auto& item = compiled.rows[row];
            if (condition_row(item.kind)) continue;
            if (item.kind == RowKind::Moiety) {
                for (std::size_t species = 0; species < components; ++species)
                    row_values[row] += amount * compiled.B[item.index][species] * x[species];
                continue;
            }
            if (item.phase != phase) continue;
            if (item.kind == RowKind::Chemical) {
                const std::size_t species = item.index;
                AD value = CppAD::log(x[species]) - CppAD::log(volume);
                AD power = AD(1.0);
                // Exact Taylor term of the standard state in the free or driven temperature.
                for (std::size_t k = 0; k < standard.derivative.size(); ++k) {
                    value += standard.derivative[k][species] * power;
                    power *= delta_T / AD(static_cast<double>(k + 1));
                }
                value += standard.pressure_derivative[species] * delta_P
                    + 0.5 * standard.pressure_second_derivative[species]
                        * delta_P * delta_P;
                if (!standard.temperature_pressure_derivative.empty())
                    value += standard.temperature_pressure_derivative[species] * delta_T * delta_P;
                for (std::size_t a = 0; a < standard.parameter_derivative.size(); ++a) {
                    AD slope = standard.parameter_derivative[a][species];
                    if (!standard.temperature_parameter_derivative.empty())
                        slope += standard.temperature_parameter_derivative[a][species] * delta_T;
                    if (!standard.pressure_parameter_derivative.empty())
                        slope += standard.pressure_parameter_derivative[a][species] * delta_P;
                    value += slope * (theta[a] - theta_center[a]);
                }
                for (std::size_t balance = 0; balance < compiled.rank; ++balance)
                    value -= u[compiled.primal_size + balance] * compiled.B[balance][species];
                if (phi[phase].has_value()) value -= AD(charges[species]) * u[*phi[phase]];
                row_values[row] = value;
            } else if (item.kind == RowKind::Pressure) {
                row_values[row] = AD(1.0) - P * volume / (GAS_CONSTANT_J_PER_MOL_K * T);
            } else if (item.kind == RowKind::Electro) {
                for (std::size_t species = 0; species < components; ++species)
                    row_values[row] += AD(charges[species]) * x[species];
            } else {
                for (const std::size_t species : spec.support) row_values[row] += x[species];
                row_values[row] -= AD(1.0);
            }
        }
    }
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.rows[row].kind != RowKind::Moiety) continue;
        for (std::size_t species = 0; species < components; ++species)
            row_values[row] -= compiled.B[compiled.rows[row].index][species] * u[n + 2 + species];
        row_values[row] /= AD(compiled.feed_total);
    }
    if (compiled.has_caloric_condition()) {
        const bool h = problem.enthalpy_j.has_value();
        row_values.back() -= u[n + 2 + components + theta_center.size() + (h ? 0 : 1)]
            / (feed_total * GAS_CONSTANT_J_PER_MOL_K * (h ? T : AD(1.0)));
    }
    if (problem.envelope) {
        // EqID envelope_local_condition; both phases share one support and the incipient phase is the zero pin.
        const auto& e = *problem.envelope;
        const std::size_t pin = problem.phases[0].is_pinned() ? 0 : 1;
        const auto& support = problem.phases[0].support;
        const std::size_t local = static_cast<std::size_t>(std::find(support.begin(), support.end(), e.component) - support.begin());
        row_values.back() = (e.kind == EnvelopeCoordinate::LogTemperature ? CppAD::log(T)
            : e.kind == EnvelopeCoordinate::LogPressure ? CppAD::log(P)
            : u[prepared.shape[pin].x[local]] - u[prepared.shape[1 - pin].x[local]]) - e.target;
    }
    output.insert(output.end(), row_values.begin(), row_values.end());
    return CppAD::ADFun<double>(u, output);
}

std::vector<double> caloric_hessian(const CompiledSystem& compiled, const PhasePoint& prepared, double weight) {
    const std::size_t n = prepared.point.size();
    std::vector<double> result(n * n, 0.0);
    if (weight == 0.0) return result;
    const std::size_t width = 2 + compiled.n_feed.size() + compiled.model->active_parameters().size();
    const std::size_t stride = width + 5;
    auto geometry = record_geometry(compiled, prepared, prepared.standard);
    using AD = CppAD::AD<double>;
    std::vector<AD> q(prepared.point.begin(), prepared.point.end());
    CppAD::Independent(q);
    std::vector<AD> input = q;
    input.push_back(prepared.temperature);
    input.push_back(prepared.pressure);
    for (auto v : compiled.n_feed) input.push_back(v);
    for (auto v : compiled.model->active_parameter_values()) input.push_back(v);
    input.push_back(compiled.problem.enthalpy_j.value_or(0.0));
    input.push_back(compiled.problem.entropy_j_per_k.value_or(0.0));
    auto active_geometry = geometry.base2ad();
    const auto mapped = active_geometry.Forward(0, input);
    AD value = mapped.back();
    for (std::size_t p = 0; p < prepared.shape.size(); ++p) {
        const auto& d = prepared.derivative[p];
        std::vector<double> center = {prepared.temperature, prepared.volume[p]};
        center.insert(center.end(), prepared.composition[p].begin(), prepared.composition[p].end());
        const auto theta = compiled.model->active_parameter_values();
        center.insert(center.end(), theta.begin(), theta.end());
        std::vector<AD> delta(width);
        for (std::size_t i = 0; i < width; ++i) delta[i] = mapped[p * stride + i] - center[i];
        AD energy = d.value, gt = d.gradient[0], gv = d.gradient[1];
        for (std::size_t i = 0; i < width; ++i) {
            energy += d.gradient[i] * delta[i];
            gt += d.hessian[i] * delta[i];
            gv += d.hessian[width + i] * delta[i];
            for (std::size_t j = 0; j < width; ++j)
                energy += 0.5 * d.hessian[i * width + j] * delta[i] * delta[j];
        }
        value += mapped[p * stride + width + 2] * energy
            + mapped[p * stride + width + 3] * gt + mapped[p * stride + width + 4] * gv;
    }
    CppAD::ADFun<double> pullback(q, std::vector<AD>{value});
    result = pullback.Hessian(prepared.point, 0);
    std::vector<double> center(prepared.point);
    center.push_back(prepared.temperature);
    center.push_back(prepared.pressure);
    center.insert(center.end(), compiled.n_feed.begin(), compiled.n_feed.end());
    const auto theta = compiled.model->active_parameter_values();
    center.insert(center.end(), theta.begin(), theta.end());
    center.push_back(compiled.problem.enthalpy_j.value_or(0.0));
    center.push_back(compiled.problem.entropy_j_per_k.value_or(0.0));
    const auto values = geometry.Forward(0, center);
    const auto jac = geometry.Jacobian(center);
    for (std::size_t p = 0; p < prepared.shape.size(); ++p) {
        std::vector<double> b(width, 0.0);
        b[0] = values[p * stride + width + 3];
        b[1] = values[p * stride + width + 4];
        std::vector<std::vector<double>> direction(n, std::vector<double>(width));
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t i = 0; i < width; ++i) direction[a][i] = jac[(p * stride + i) * center.size() + a];
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t c = 0; c <= a; ++c) {
                const auto action = compiled.model->residual_energy_action(prepared.temperature,
                    prepared.volume[p], prepared.composition[p], {b, direction[a], direction[c]});
                if (!action.action) throw std::domain_error(std::string("caloric constraint curvature: ") + solved_state_status_name(solved_state_status(action.status)));
                result[a * n + c] += *action.action;
                if (a != c) result[c * n + a] += *action.action;
            }
    }
    for (auto& v : result) v *= weight;
    return result;
}

// Weighted critical-row curvature by exact polarization of second Taylor coefficients from the one row
// owner: H_ab = c2(e_a + e_b) - c2(e_a) - c2(e_b), c2(d) = sum_r m_r [t^2] r(y + t d) (E4 and E5).
// ponytail: primal(primal+1)/2 order-2 series, critical cells only; a pullback if C1-C3 cost matters.
std::vector<double> critical_hessian(const CompiledSystem& compiled, const PhasePoint& prepared,
                                     const std::vector<double>& multipliers) {
    const std::size_t n = prepared.point.size(), first = compiled.rows.size() - compiled.mode_count - 2;
    std::vector<double> result(n * n, 0.0), diagonal(compiled.primal_size, 0.0);
    if (std::all_of(multipliers.begin() + static_cast<std::ptrdiff_t>(first), multipliers.end(),
                    [](double m) { return m == 0.0; })) return result;
    std::vector<std::vector<double>> state(3, std::vector<double>(n, 0.0));
    state[0] = prepared.point;
    const std::vector<double> control(2 + compiled.n_feed.size() + compiled.model->active_parameters().size(), 0.0);
    const auto c2 = [&] {
        const auto rows = compiled.row_series(prepared, state, control, 2).rows[2];
        return std::inner_product(rows.begin() + static_cast<std::ptrdiff_t>(first), rows.end(),
                                  multipliers.begin() + static_cast<std::ptrdiff_t>(first), 0.0);
    };
    for (std::size_t a = 0; a < compiled.primal_size; ++a) {
        state[1][a] = 1.0;
        diagonal[a] = c2();
        result[a * n + a] = 2.0 * diagonal[a];
        for (std::size_t b = 0; b < a; ++b) {
            state[1][b] = 1.0;
            result[a * n + b] = result[b * n + a] = c2() - diagonal[a] - diagonal[b];
            state[1][b] = 0.0;
        }
        state[1][a] = 0.0;
    }
    return result;
}

// Reaction-aware feed-relative pin admission: individual [0, 1] domain always; the exact
// sum rules only when the declared reaction space conserves total moles.
void admit_pins(CompiledSystem& compiled) {
    double pin_sum = 0.0;
    std::size_t pinned = 0;
    for (const auto& phase : compiled.problem.phases) {
        if (phase.amount.kind != AmountKind::Pinned) continue;
        const double fraction = phase.amount.fraction;
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
            throw std::invalid_argument("invalid pinned phase fraction in phase " + phase.name);
        }
        pin_sum += fraction;
        ++pinned;
    }
    if (!compiled.conserves_total_moles || pinned == 0) return;
    if (pin_sum > 1.0) {
        throw std::invalid_argument("pinned phase fractions exceed one: S_pin="
                                    + core::optional_number(pin_sum));
    }
    if (pinned == compiled.problem.phases.size() && pin_sum != 1.0) {
        throw std::invalid_argument("all-pinned phase fractions must sum to one: S_pin="
                                    + core::optional_number(pin_sum));
    }
}

}  // namespace

// EqID: reference_species_calorics
// Supplied ideal-gas records fix h and cp of their species; record-less species solve
// N_D x_D = Delta_r x - N_S x_S on the density basis for x = h, cp (min-norm), where only
// the ionic charge gauge may remain free: it cancels in every electroneutral phase.
ReferenceCalorics reference_calorics(const CompiledSystem &compiled, double temperature, double pressure,
                                     bool heat_capacity, bool pressure_slope, bool parameter_slopes,
                                     std::size_t derivative_order) {
    derivative_order = std::max(derivative_order, heat_capacity ? std::size_t(2) : std::size_t(1));
    using AD = CppAD::AD<double>;
    ReferenceCalorics result;
    const auto &model = *compiled.model;
    const auto &problem = compiled.problem;
    const auto &record = model.thermochemistry();
    const std::size_t components = model.parameters().component_count(), rows = problem.reactions.size();
    const auto refuse = [&](SolvedStateActionStatus status, const std::string &message) {
        result.status = status;
        result.message = message;
        return result;
    };
    if (!record)
        return refuse(SolvedStateActionStatus::MissingThermochemistry,
                      "reference calorics need record anchors");
    std::vector<std::size_t> derived;
    std::vector<AD> T = {temperature}, output;
    std::vector<std::optional<IdealCorrelation>> selected(components);
    for (std::size_t i = 0; i < components; ++i) {
        const auto& intervals = record->component_intervals[i];
        if (intervals.empty()) {
            derived.push_back(i);
            continue;
        }
        const auto found = std::find_if(intervals.begin(), intervals.end(),
            [&](const IdealInterval& interval) { return ideal_interval_contains(interval, temperature); });
        if (found == intervals.end())
            return refuse(SolvedStateActionStatus::OutsideInterval, "record excludes temperature");
        result.thermal_boundary |= temperature == found->lower || temperature == found->upper;
        selected[i] = found->correlation;
    }
    // One exact correlation tape supplies record g_rho and reaction source g through the consumed order.
    CppAD::Independent(T);
    const AD rt = GAS_CONSTANT_J_PER_MOL_K * T[0];
    for (const auto &correlation : selected) {
        AD g = 0.0;
        if (correlation) {
            const auto value = evaluate_ideal_correlation(*correlation, T[0], record->standard_pressure);
            g = value.h - T[0] * value.s + rt * CppAD::log(rt / record->standard_pressure);
        }
        output.push_back(g);
    }
    for (const auto &reaction : problem.reactions) {
        AD g = 0.0;
        if (reaction.affinity_j_per_mol)
            g = *reaction.affinity_j_per_mol;
        else {
            AD L = reaction.equilibrium_constant ? std::log(*reaction.equilibrium_constant) : 0.0;
            if (reaction.correlation) {
                const auto &law = *reaction.correlation;
                (void)evaluate_reaction_properties(law, temperature);
                L = reaction_log_value(law, T[0]);
            }
            if (reaction.reference &&
                reaction.reference->source_basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution)
                L += reaction_reference_log_scale(model, *problem.source_reference, reaction.stoichiometry);
            g = -rt * L;
        }
        output.push_back(g);
    }
    CppAD::ADFun<double> sources(T, output);
    std::vector<std::vector<double>> thermal(derivative_order + 1);
    double factorial = 1.0;
    for (std::size_t a = 0; a <= derivative_order; ++a) {
        thermal[a] = sources.Forward(a, std::vector<double>{a == 0 ? temperature : a == 1 ? 1.0 : 0.0});
        if (a)
            factorial *= double(a);
        for (auto &value : thermal[a])
            value *= factorial;
    }
    Eigen::MatrixXd matrix(rows, derived.size());
    Eigen::VectorXd charge(derived.size());
    for (std::size_t j = 0; j < derived.size(); ++j) {
        charge(j) = model.parameters().charge_number[derived[j]];
        for (std::size_t r = 0; r < rows; ++r)
            matrix(r, j) = problem.reactions[r].stoichiometry[derived[j]];
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition;
    if (!derived.empty())
        decomposition.compute(matrix);
    const auto nullity =
        derived.empty() ? Eigen::Index(0) : Eigen::Index(derived.size()) - decomposition.rank();
    if (nullity > 1 ||
        (nullity == 1 &&
         (charge.norm() == 0.0 ||
          (matrix * charge).norm() > std::sqrt(std::numeric_limits<double>::epsilon()) * charge.norm())))
        return refuse(SolvedStateActionStatus::ReferenceUnavailable,
                      "record-less calorics undetermined beyond charge gauge");
    result.reaction_map.assign(rows, std::vector<double>(components, 0.0));
    for (std::size_t r = 0; r < rows; ++r) {
        Eigen::VectorXd unit = Eigen::VectorXd::Zero(rows);
        unit(r) = 1.0;
        if (!derived.empty()) {
            const Eigen::VectorXd image = decomposition.solve(unit);
            for (std::size_t j = 0; j < derived.size(); ++j)
                result.reaction_map[r][derived[j]] = image(j);
        }
    }
    result.reaction_defect.assign(rows, std::vector<double>(rows, 0.0));
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t s = 0; s < rows; ++s) {
            result.reaction_defect[r][s] = r == s ? 1.0 : 0.0;
            for (std::size_t i = 0; i < components; ++i)
                result.reaction_defect[r][s] -= problem.reactions[r].stoichiometry[i] * result.reaction_map[s][i];
        }
    std::vector<std::vector<double>> nu;
    std::vector<ReactionReference> references;
    std::vector<std::size_t> reference_rows;
    for (std::size_t r = 0; r < rows; ++r)
        if (problem.reactions[r].reference) {
            nu.push_back(problem.reactions[r].stoichiometry);
            references.push_back(*problem.reactions[r].reference);
            reference_rows.push_back(r);
        }
    const auto diagnostics =
        references.empty()
            ? std::vector<NeutralReferenceReactionDiagnostic>{}
            : reaction_references(model, temperature, pressure, *problem.source_reference, nu, references,
                                  pressure_slope, std::nullopt, std::min(derivative_order, std::size_t(2)),
                                  parameter_slopes, derivative_order);
    const auto partial = [](const NeutralReferenceReactionDiagnostic &d, std::size_t a, std::size_t b,
                            std::optional<std::size_t> p) -> std::optional<double> {
        if (p && a + b == 0)
            return d.terminal_contraction_parameter_derivative[*p];
        if (p && a == 1 && b == 0)
            return d.terminal_contraction_temperature_parameter_derivative[*p];
        if (!p && a + b == 0)
            return d.terminal_contraction;
        if (!p && a == 1 && b == 0)
            return d.terminal_contraction_temperature_derivative;
        if (!p && a == 2 && b == 0)
            return d.terminal_contraction_temperature_second_derivative;
        if (!p && a == 0 && b == 1)
            return d.terminal_contraction_pressure_derivative;
        if (!p && a == 0 && b == 2)
            return d.terminal_contraction_pressure_second_derivative;
        if (!p && a == 1 && b == 1)
            return d.terminal_contraction_temperature_pressure_derivative;
        for (const auto &column : d.derivatives)
            if (column.temperature_order == a && column.pressure_order == b && column.parameter == p)
                return column.value;
        return std::nullopt;
    };
    for (std::size_t c = 0; c <= std::size_t(parameter_slopes); ++c)
        for (std::size_t a = 0; a + c <= derivative_order; ++a)
            for (std::size_t b = 0; a + b + c <= derivative_order; ++b) {
                if (b && !pressure_slope)
                    continue;
                for (std::size_t p = 0; p < (c ? model.active_parameters().size() : 1); ++p) {
                    ReferenceCalorics::Derivative entry{
                        a, b, c ? std::optional<std::size_t>(p) : std::nullopt, {}};
                    auto &column = entry.column;
                    column.values.assign(components, 0.0);
                    Eigen::VectorXd right(rows), scale(rows);
                    for (std::size_t i = 0; i < components; ++i)
                        if (!b && !c)
                            column.values[i] = thermal[a][i];
                    for (std::size_t r = 0; r < rows; ++r)
                        right(r) = !b && !c ? thermal[a][components + r] : 0.0;
                    for (std::size_t r = 0; r < diagnostics.size(); ++r) {
                        const auto upper = partial(diagnostics[r], a, b, entry.parameter);
                        const auto lower = a ? partial(diagnostics[r], a - 1, b, entry.parameter)
                                             : std::optional<double>(0.0);
                        if (!upper || !lower) {
                            column.status = SolvedStateActionStatus::ReferenceUnavailable;
                            const auto missing_a = !upper ? a : a - 1;
                            column.message = "reference_unavailable: reaction " +
                                reaction_name(problem.reactions[reference_rows[r]], reference_rows[r]) +
                                " C_" + std::to_string(missing_a) + std::to_string(b) +
                                (c ? " parameter " + std::to_string(p) : "") + ": " +
                                diagnostics[r].derivative_cause(missing_a, b, entry.parameter);
                        } else
                            right(reference_rows[r]) -=
                                GAS_CONSTANT_J_PER_MOL_K * (temperature * *upper + double(a) * *lower);
                    }
                    scale = right.cwiseAbs();
                    for (std::size_t r = 0; r < rows; ++r)
                        for (std::size_t i = 0; i < components; ++i)
                            if (selected[i]) {
                                const double term = problem.reactions[r].stoichiometry[i] * column.values[i];
                                right(r) -= term;
                                scale(r) += std::abs(term);
                            }
                    if (column.status == SolvedStateActionStatus::Available) {
                        const Eigen::VectorXd solution =
                            derived.empty() ? Eigen::VectorXd() : Eigen::VectorXd(decomposition.solve(right));
                        const Eigen::VectorXd residual =
                            derived.empty() ? right : Eigen::VectorXd(matrix * solution - right);
                        for (Eigen::Index r = 0; r < residual.size(); ++r)
                            if (!(std::abs(residual(r)) <= std::sqrt(std::numeric_limits<double>::epsilon()) *
                                                               std::max(1.0, scale(r)))) {
                                column.status = SolvedStateActionStatus::ReferenceUnavailable;
                                column.message = "supplied records inconsistent with reaction " +
                                                 reaction_name(problem.reactions[r], r);
                            }
                        for (std::size_t j = 0; j < derived.size(); ++j)
                            column.values[derived[j]] = solution(j);
                    }
                    result.gibbs.push_back(std::move(entry));
                }
            }
    return result;
}

static CompiledSystem compile_with_basis(
    const HelmholtzModel& model, const EquilibriumProblem& request,
    const std::vector<std::vector<double>>* trial_basis, core::Scope& scope
) {
    if (!request.topology_declared) {
        throw std::invalid_argument("phase count unknown; declare phases");
    }
    if (request.phases.empty()) throw std::invalid_argument("no phases declared");
    validate_intensive_seed(request.temperature, "temperature");
    validate_intensive_seed(request.pressure, "pressure");
    if (!std::isfinite(request.row_tolerance) || request.row_tolerance <= 0.0) {
        throw std::invalid_argument("invalid row tolerance");
    }
    const std::size_t components = model.parameters().component_count();
    CompiledSystem compiled;
    compiled.model = &model;
    compiled.problem = request;
    compiled.n_feed = resolve_feed(model, request.feed);
    compiled.feed_total = feed_total(compiled.n_feed);
    compiled.has_reaction = !request.reactions.empty();
    const auto& z = model.parameters().charge_number;
    double net_z = 0.0;
    for (std::size_t i = 0; i < components; ++i) net_z += z[i] * compiled.n_feed[i];
    if (std::abs(net_z) > std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(1.0, vector_norm(compiled.n_feed))) {
        throw std::invalid_argument("feed carries net z");
    }
    for (std::size_t row = 0; row < request.reactions.size(); ++row) {
        const auto& nu = request.reactions[row].stoichiometry;
        if (nu.size() != components) throw std::invalid_argument("stoichiometry size");
        double total = 0.0;
        for (std::size_t i = 0; i < components; ++i) total += z[i] * nu[i];
        if (std::abs(total) > std::sqrt(std::numeric_limits<double>::epsilon())
            * std::max(1.0, vector_norm(nu))) {
            throw std::invalid_argument("reaction " + reaction_name(request.reactions[row], row)
                                        + " violates z");
        }
    }
    compiled.B = trial_basis ? *trial_basis : balance_basis(components, request.reactions);
    compiled.rank = compiled.B.size();
    compiled.conserves_total_moles = conserves_total_moles(compiled.B);
    const double temperature_guess = request.temperature.seed_or(300.0);
    const double pressure_guess = request.pressure.seed_or(1.0e5);
    const auto central_reactions = central_reference_reactions(
        model, request, temperature_guess, pressure_guess,
        compiled.reference_diagnostics, false);
    compiled.standard_gibbs = standard_state(
        temperature_guess, components, central_reactions,
        has_reaction_reference(request) ? std::vector<double>{}
                                        : request.standard_gibbs_j_per_mol
    );
    if (request.temperature.free && !request.standard_gibbs_j_per_mol.empty()) {
        // A constant supplied vector admits a free temperature only for constant targets.
        const auto targets = reaction_target_series(temperature_guess, request.reactions, 3);
        for (std::size_t row = 0; row < request.reactions.size(); ++row) {
            const double scale = std::sqrt(std::numeric_limits<double>::epsilon())
                * std::max(1.0, std::abs(targets[0][row]));
            double power = 1.0;
            for (std::size_t k = 1; k <= 3; ++k) {
                power *= temperature_guess;
                if (std::abs(targets[k][row]) * power > scale) {
                    throw std::invalid_argument(
                        "temperature-dependent standard state required for reaction "
                        + reaction_name(request.reactions[row], row));
                }
            }
        }
    }
    compiled.feed_composition = normalized_feed(compiled.n_feed);
    compiled.balance = matvec(compiled.B, compiled.n_feed);
    std::vector<bool> covered(components, false);
    std::vector<bool> participant(components, false);
    for (const auto& item : request.reactions) {
        for (std::size_t i = 0; i < components; ++i)
            participant[i] = participant[i] || item.stoichiometry[i] != 0.0;
    }
    for (auto& phase : compiled.problem.phases) {
        if (phase.support_declared && phase.support.empty())
            throw std::invalid_argument("empty support in phase " + phase.name);
        phase.support = complete_support(components, phase.support);
        compiled.declared_supports.push_back(phase.support);
    }
    auto inventory = compiled;
    if (!trial_basis) {
        auto& support = inventory.problem.phases.front().support;
        for (std::size_t i = 0; i < components; ++i)
            if (participant[i] && std::find(support.begin(), support.end(), i) == support.end()) support.push_back(i);
    }
    std::vector<std::size_t> indices;
    for (const auto& phase : inventory.problem.phases)
        indices.insert(indices.end(), phase.support.begin(), phase.support.end());
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    const auto bounds = inventory_bounds(inventory, indices, false, 0.0);
    compiled.inventory_upper_bounds.assign(components, 1.0);
    std::vector<bool> absent(components, false);
    const auto certify = [&](std::size_t i, std::vector<double> candidate) {
        candidate = reconstruct_conserved_functional(std::move(candidate), request.reactions);
        return !candidate.empty() && candidate[i] > 0.0
            && std::all_of(candidate.begin(), candidate.end(), [](double value) { return value >= 0.0; })
            && exact_zero_dot(candidate, compiled.n_feed);
    };
    for (const auto& bound : bounds) {
        if (bound.status == "infeasible_inventory") throw std::invalid_argument("infeasible_inventory: " + bound.termination);
        if (bound.status == "optimal") compiled.inventory_upper_bounds[bound.species] = *bound.feed_normalized_amount;
        if (!trial_basis && bound.status == "optimal" && *bound.maximum == 0.0 && !bound.conserved_functional.empty())
            absent[bound.species] = certify(bound.species, bound.conserved_functional);
    }
    if (!trial_basis) for (std::size_t i = 0; i < components; ++i) {
        for (const auto& row : compiled.B) absent[i] = absent[i] || certify(i, row);
        if (absent[i]) compiled.removed_species.push_back(i);
    }
    for (auto& phase : compiled.problem.phases) {
        phase.support_declared = true;
        const auto support_before = phase.support;
        phase.support.erase(std::remove_if(
            phase.support.begin(), phase.support.end(), [&](std::size_t index) {
                return absent[index] || (compiled.n_feed[index] == 0.0 && !participant[index]);
            }), phase.support.end());
        if (phase.is_pinned() && phase.amount.fraction == 0.0) {
            const auto removed = std::find_if(support_before.begin(), support_before.end(),
                [&](std::size_t index) { return absent[index]; });
            if (removed != support_before.end()) {
                std::string reason;
                if (phase.support.empty()) reason = "empty effective support";
                else if (std::any_of(phase.support.begin(), phase.support.end(),
                    [&](std::size_t index) { return z[index] != 0; }))
                    reason = "effective incipient support must contain only uncharged species";
                else if (request.phases.size() != 2 || request.phases[0].kind == request.phases[1].kind)
                    reason = "requires exactly two phases of opposite liquid/vapor kinds";
                else if (request.temperature.free || !request.temperature.is_positive() || !request.pressure.free)
                    reason = "requires fixed positive T and free P";
                else if (std::count_if(request.phases.begin(), request.phases.end(),
                    [](const PhaseSpec& item) { return item.is_free(); }) != 1)
                    reason = "requires one zero pin and one free amount";
                if (!reason.empty()) throw std::invalid_argument("zero_inventory_incipient_support: component "
                    + std::to_string(*removed) + " in phase " + phase.name + ": " + reason);
            }
        }
        restrict_composition_guess(phase, support_before);
        for (const std::size_t index : support_before) {
            if (std::find(phase.support.begin(), phase.support.end(), index) == phase.support.end()
                && std::find(compiled.removed_species.begin(), compiled.removed_species.end(), index)
                    == compiled.removed_species.end()) {
                compiled.removed_species.push_back(index);
            }
        }
        if (phase.support.empty()) throw std::invalid_argument("empty support in phase " + phase.name);
        for (const std::size_t index : phase.support) covered[index] = true;
        bool positive = false;
        bool negative = false;
        for (const std::size_t index : phase.support) {
            positive = positive || z[index] > 0;
            negative = negative || z[index] < 0;
        }
        compiled.has_ions = compiled.has_ions || positive || negative;
        if (positive != negative && (positive || negative)) {
            throw std::invalid_argument("support z sign in phase " + phase.name);
        }
    }
    for (std::size_t i = 0; i < components; ++i) {
        if (!covered[i] && !absent[i] && (compiled.n_feed[i] != 0.0 || participant[i])) {
            throw std::invalid_argument("species " + std::to_string(i) + " has no support");
        }
    }
    if (!trial_basis && !compiled.removed_species.empty()) restrict_inventory_balances(compiled, covered);
    admit_pins(compiled);
    compiled.pinned_count = 0;
    for (const auto& phase : compiled.problem.phases) {
        compiled.pinned_count += phase.amount.kind == AmountKind::Pinned;
    }
    compiled.free_intensive_count = request.temperature.free + request.pressure.free;
    const std::size_t caloric_count = compiled.has_caloric_condition();
    admit_conditions(compiled);
    if (compiled.free_intensive_count != compiled.pinned_count + compiled.condition_count()) {
        if (caloric_count) throw std::invalid_argument("PH/PS requires declared free amounts, free T and fixed P");
        if (compiled.free_intensive_count < compiled.pinned_count) {
            std::string phase_name;
            for (const auto& phase : compiled.problem.phases) {
                if (phase.amount.kind == AmountKind::Pinned) {
                    phase_name = phase.name;
                    break;
                }
            }
            throw std::invalid_argument(
                "pinned phase " + phase_name + " requires a free T or P"
            );
        }
        throw std::invalid_argument(
            "free temperature or pressure requires a pinned phase amount"
        );
    }
    if (compiled.problem.phases.size() - compiled.pinned_count > compiled.rank + caloric_count) {
        throw std::invalid_argument(
            "phase count exceeds balance rank for phase "
            + compiled.problem.phases.back().name
            + "; pin a phase or free T or P"
        );
    }
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        for (const std::size_t index : compiled.problem.phases[phase].support) {
        compiled.rows.push_back({RowKind::Chemical, phase, index,
                                 "chemical_potential_stationarity"});
            compiled.row_names.emplace_back(row_kind_name(RowKind::Chemical));
        }
        compiled.rows.push_back({RowKind::Pressure, phase, 0, "pressure_equality"});
        compiled.row_names.emplace_back(row_kind_name(RowKind::Pressure));
    }
    std::vector<std::size_t> ion_phases;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        for (const std::size_t index : compiled.problem.phases[phase].support) {
            if (z[index] != 0) {
                ion_phases.push_back(phase);
                break;
            }
        }
    }
    std::optional<std::size_t> rowless;
    for (const std::size_t phase : ion_phases) {
        if (compiled.problem.phases[phase].amount.kind == AmountKind::Free
            || compiled.problem.phases[phase].amount.fraction > 0.0) {
            rowless = phase;
            break;
        }
    }
    for (const std::size_t phase : ion_phases) {
        if (rowless.has_value() && phase == *rowless) continue;
        compiled.rows.push_back({RowKind::Electro, phase, 0, "phase_electroneutrality"});
        compiled.row_names.emplace_back(row_kind_name(RowKind::Electro));
        ++compiled.charge_row_count;
    }
    for (std::size_t row = 0; row < compiled.B.size(); ++row) {
        compiled.rows.push_back({RowKind::Moiety, 0, row, "material_balance"});
        compiled.row_names.emplace_back(row_kind_name(RowKind::Moiety));
    }
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        compiled.rows.push_back({RowKind::Normalization, phase, 0, "phase_normalization"});
        compiled.row_names.emplace_back(row_kind_name(RowKind::Normalization));
    }
    if (caloric_count) {
        const auto kind = request.enthalpy_j ? RowKind::Enthalpy : RowKind::Entropy;
        compiled.rows.push_back({kind, 0, 0, row_kind_name(kind)});
        compiled.row_names.emplace_back(row_kind_name(kind));
    }
    // Critical rows M u (C), u'u - 1 and c; or one envelope row. Their positions are the last rows.
    compiled.mode_count = request.critical ? components : 0;
    for (std::size_t i = 0; i < compiled.mode_count + 2 * request.critical.has_value() + request.envelope.has_value(); ++i) {
        const auto kind = request.envelope ? RowKind::EnvelopeCondition : i < components ? RowKind::CriticalMode
            : i == components ? RowKind::CriticalNormalization : RowKind::CriticalCubic;
        compiled.rows.push_back({kind, 0, i, row_kind_name(kind)});
        compiled.row_names.emplace_back(row_kind_name(kind));
    }
    const auto shape = layouts(compiled);
    compiled.primal_size = 0;
    for (const auto& item : shape) compiled.primal_size = std::max(compiled.primal_size, item.volume + 1);
    compiled.primal_size += compiled.mode_count + compiled.free_intensive_count;
    compiled.rows_minus_unknowns = static_cast<long>(compiled.pinned_count + compiled.condition_count())
        - static_cast<long>(compiled.free_intensive_count);
    compiled.generic_nullity = 1;
    compiled.feasible = false;
    scope.probe("chemical_rows", static_cast<double>(count_rows(compiled, RowKind::Chemical)));
    scope.probe("pressure_rows", static_cast<double>(count_rows(compiled, RowKind::Pressure)));
    scope.probe("balance_rows", static_cast<double>(count_rows(compiled, RowKind::Moiety)));
    scope.probe("normalization_rows",
                static_cast<double>(count_rows(compiled, RowKind::Normalization)));
    scope.probe("removed_species", static_cast<double>(compiled.removed_species.size()));
    scope.probe("conserves_total_moles", compiled.conserves_total_moles ? 1.0 : 0.0);
    if (compiled.is_square()) {
        const auto point = compiled.generic_point();
        compiled.F_y = compiled.jacobian(point);
        Eigen::MatrixXd matrix(
            static_cast<Eigen::Index>(compiled.rows.size()),
            static_cast<Eigen::Index>(point.size()));
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            for (std::size_t column = 0; column < point.size(); ++column) {
                matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                    compiled.F_y[row * point.size() + column];
            }
        }
        compiled.jacobian_rank = static_cast<std::size_t>(matrix.fullPivLu().rank());
        compiled.generic_nullity = point.size() > compiled.jacobian_rank
            ? point.size() - compiled.jacobian_rank : 0;
        compiled.feasible = compiled.generic_nullity == 0;
    }
    if (compiled.row_names.size() != compiled.rows.size()) {
        throw std::logic_error("row labels and rows disagree");
    }
    for (const auto& row : compiled.rows) {
        if (row.eqid == nullptr || std::string(row.eqid).empty()) {
            throw std::logic_error("equilibrium row has no EqID");
        }
    }
    scope.probe("jacobian_rank", static_cast<double>(compiled.jacobian_rank));
    scope.probe("unknowns", static_cast<double>(compiled.unknown_count()));
    scope.probe("constraints", static_cast<double>(compiled.constraint_count()));
    // The compiled record is the sole source for both the driver and callers.
    scope.probe("rows", static_cast<double>(compiled.rows.size()));
    scope.probe("rank", static_cast<double>(compiled.rank));
    return compiled;
}

void restrict_composition_guess(PhaseSpec& phase, const std::vector<std::size_t>& support) {
    const auto guess = std::move(phase.composition_guess);
    phase.composition_guess.clear();
    if (guess.size() != support.size()) return;
    for (std::size_t i = 0; i < support.size(); ++i)
        if (std::find(phase.support.begin(), phase.support.end(), support[i]) != phase.support.end())
            phase.composition_guess.push_back(guess[i]);
}

CompiledSystem compile_problem(const HelmholtzModel& model, const EquilibriumProblem& request) {
    core::Scope scope("compile_problem");
    return compile_with_basis(model, request, nullptr, scope);
}

// The reduced trial space is internal: the public no-reaction compiler retains species balances.
CompiledSystem compile_stability_trial(const HelmholtzModel& model, double T, double P,
    const PhaseSpec& phase, const std::vector<double>& w, const std::vector<double>& c,
    std::vector<std::vector<double>>& transformation) {
    core::Scope scope("compile_problem");
    const auto& support = phase.support;
    const auto& z = model.parameters().charge_number;
    const bool has_charge = std::any_of(support.begin(), support.end(), [&](auto i) { return z[i] != 0; });
    Eigen::MatrixXd physical(has_charge ? 2 : 1, support.size());
    physical.row(0).setOnes();
    if (has_charge) for (std::size_t i = 0; i < support.size(); ++i) physical(1, i) = z[support[i]];
    const Eigen::MatrixXd kernel = physical.fullPivLu().kernel();
    std::vector<ReactionSpec> directions;
    if (support.size() > static_cast<std::size_t>(physical.rows())) {
        for (Eigen::Index j = 0; j < kernel.cols(); ++j) {
            ReactionSpec direction;
            direction.stoichiometry.assign(kernel.col(j).data(), kernel.col(j).data() + kernel.rows());
            // Only the basis owner's input validation consumes this zero target.
            // These directions never enter the trial's reaction declaration.
            direction.affinity_j_per_mol = 0.0;
            directions.push_back(std::move(direction));
        }
    }
    const auto reduced = balance_basis(support.size(), directions);
    std::vector<std::vector<double>> B(reduced.size(), std::vector<double>(w.size(), 0.0));
    transformation.clear();
    for (std::size_t k = 0; k < B.size(); ++k) {
        for (std::size_t i = 0; i < support.size(); ++i) B[k][support[i]] = reduced[k][i];
        const Eigen::VectorXd coefficients = physical.transpose().colPivHouseholderQr().solve(
            Eigen::Map<const Eigen::VectorXd>(reduced[k].data(), reduced[k].size()));
        transformation.emplace_back(coefficients.data(), coefficients.data() + coefficients.size());
    }
    EquilibriumProblem request;
    request.temperature = FreeValue::fixed(T); request.pressure = FreeValue::fixed(P);
    request.feed.values = w;
    request.phases = {phase};
    request.phases[0].amount = AmountSpec::free();
    request.row_tolerance = policy::EQUILIBRIUM_ROW_RESIDUAL.value;
    for (double coefficient : c) request.standard_gibbs_j_per_mol.push_back(coefficient * GAS_CONSTANT_J_PER_MOL_K * T);
    return compile_with_basis(model, request, &B, scope);
}

std::vector<double> CompiledSystem::generic_point() const {
    core::Scope scope("equilibrium_generic_point");
    const auto shape = layouts(*this);
    std::vector<double> point(primal_size + rank + charge_row_count, 0.0);
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        double total = 0.0;
        for (std::size_t i = 0; i < shape[phase].x.size(); ++i) {
            const double value = 1.0 + 0.37 * static_cast<double>(i) + 0.11 * static_cast<double>(phase)
                + (phase > 1 ? 0.05 * std::pow(static_cast<double>(i), phase) : 0.0);  // 3+ phases not coplanar
            point[shape[phase].x[i]] = value;
            total += value;
        }
        for (const std::size_t coordinate : shape[phase].x)
            point[coordinate] = std::log(point[coordinate] / total);
        if (shape[phase].amount.has_value()) point[*shape[phase].amount] = 1.0;
        point[shape[phase].volume] = std::log(0.25 + 0.05 * static_cast<double>(phase));
    }
    for (std::size_t i = 0; i < mode_count; ++i)  // a generic mode direction, aligned with no axis
        point[mode_offset() + i] = (1.0 + 0.29 * static_cast<double>(i)) / static_cast<double>(mode_count + 1);
    std::size_t offset = primal_size - free_intensive_count;
    if (problem.temperature.free) point[offset++] = problem.temperature.seed_or(300.0);
    if (problem.pressure.free) point[offset++] = problem.pressure.seed_or(1.0e5);
    return point;
}

std::vector<std::size_t> CompiledSystem::phase_offsets() const {
    std::vector<std::size_t> result;
    for (const auto& item : layouts(*this)) result.push_back(item.x.front());
    return result;
}

PhasePoint CompiledSystem::prepare(const std::vector<double>& point, bool state_only) const {
    core::Scope scope("equilibrium_prepare");
    if (point.size() < primal_size) {
        throw std::invalid_argument("equilibrium point is too short");
    }
    const std::size_t components = model->parameters().component_count();
    PhasePoint result;
    result.point = point;
    result.shape = layouts(*this);
    std::size_t offset = primal_size - free_intensive_count;
    result.temperature = problem.temperature.free ? point[offset++] : problem.temperature.value;
    result.pressure = problem.pressure.free ? point[offset++] : problem.pressure.value;
    admit_calorics(*this, result.temperature);
    for (std::size_t phase = 0; phase < result.shape.size(); ++phase) {
        const auto& layout = result.shape[phase];
        std::vector<double> values(components, 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < layout.x.size(); ++i) {
            const double value = std::exp(point[layout.x[i]]);
            if (!std::isfinite(value) || value <= 0.0) throw std::invalid_argument("invalid composition");
            values[problem.phases[phase].support[i]] = value;
            total += value;
        }
        if (!(total > 0.0)) throw std::invalid_argument("composition has zero total");
        result.packing.push_back(packing_terms(*model, result.temperature, values));
        result.eta.push_back(std::exp(point[layout.volume]));
        if (!admitted_packing_fraction(result.eta.back()))
            throw std::invalid_argument("phase packing fraction outside the fluid domain");
        result.volume.push_back(result.packing.back().coefficient / result.eta.back());
        result.amount.push_back(layout.amount.has_value()
            ? point[*layout.amount] : phase_amount(problem.phases[phase], n_feed));
        const auto order = !state_only ? DerivativeOrder::Hessian : problem.temperature.free
            ? DerivativeOrder::PhaseHessian : DerivativeOrder::IsothermalPhaseHessian;
        result.derivative.push_back(model->evaluate(
            result.temperature, result.volume.back(), values, order));
        scope.probe("gradient_floor", result.derivative.back().gradient_floor,
                    static_cast<std::uint32_t>(phase));
        result.composition.push_back(std::move(values));
    }
    result.standard = standard_term(
        *this, result.temperature, result.pressure, 2, problem.pressure.free, problem.temperature.free);
    return result;
}

PointEvaluation CompiledSystem::evaluate(const std::vector<double>& point) const {
    return evaluate_at(prepare(point));
}

PointEvaluation CompiledSystem::evaluate_at(const PhasePoint& prepared) const {
    core::Scope scope("equilibrium_evaluate");
    const auto& point = prepared.point;
    const auto& shape = prepared.shape;
    const std::size_t components = model->parameters().component_count();
    const auto& z = model->parameters().charge_number;
    const auto& x = prepared.composition;
    const auto& amount = prepared.amount;
    const auto& volume = prepared.volume;
    std::vector<double> lambda(rank, 0.0);
    for (std::size_t row = 0; row < rank && primal_size + row < point.size(); ++row) {
        lambda[row] = point[primal_size + row];
    }
    std::vector<double> phi(charge_row_count, 0.0);
    for (std::size_t row = 0; row < charge_row_count
         && primal_size + rank + row < point.size(); ++row) {
        phi[row] = point[primal_size + rank + row];
    }
    const double rt = GAS_CONSTANT_J_PER_MOL_K * prepared.temperature;
    const auto& standard = prepared.standard.derivative[0];
    PointEvaluation result;
    for (const auto& derivative : prepared.derivative) result.floor = std::max(result.floor, derivative.gradient_floor);
    std::vector<std::optional<std::size_t>> phase_phi(shape.size());
    std::size_t phi_index = 0;
    for (const auto& item : rows) {
        if (item.kind == RowKind::Electro) phase_phi[item.phase] = phi_index++;
    }
    std::vector<double> series;  // every condition row from one zeroth coefficient
    for (const auto& row : rows) {
        double value = 0.0;
        if (row.kind == RowKind::Chemical) {
            const auto& derivative = prepared.derivative[row.phase];
            const std::size_t species = row.index;
            double constraint = 0.0;
            for (std::size_t balance_row = 0; balance_row < rank; ++balance_row) {
                constraint += lambda[balance_row] * B[balance_row][species];
            }
            const double electro = z[species] * (
                phase_phi[row.phase].has_value() ? phi[*phase_phi[row.phase]] : 0.0);
            value = std::log(std::max(x[row.phase][species] / volume[row.phase],
                                      std::numeric_limits<double>::min()))
                + derivative.gradient[2 + species] + standard[species]
                - constraint - electro;
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = prepared.derivative[row.phase];
            value = volume[row.phase] * (
                1.0 / volume[row.phase] - derivative.gradient[1] - prepared.pressure / rt);
        } else if (row.kind == RowKind::Moiety) {
            for (std::size_t phase = 0; phase < shape.size(); ++phase) {
                for (std::size_t species = 0; species < components; ++species) {
                    value += amount[phase] * B[row.index][species] * x[phase][species];
                }
            }
            value = (value - balance[row.index]) / feed_total;
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t species = 0; species < components; ++species) value += z[species] * x[row.phase][species];
        } else if (condition_row(row.kind)) {
            if (series.empty())
                series = row_series(prepared, {point}, std::vector<double>(2 + components + model->active_parameters().size()), 0).rows[0];
            value = series[result.values.size()];
        } else {
            for (const std::size_t species : problem.phases[row.phase].support) value += x[row.phase][species];
            value -= 1.0;
        }
        result.values.push_back(value);
    }
    scope.probe("row_values", result.values);
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        const auto& derivative = prepared.derivative[phase];
        double ideal = ideal_mechanical_helmholtz(x[phase], volume[phase], 1.0);
        result.objective += amount[phase] * (
            derivative.value + ideal + prepared.pressure * volume[phase] / rt
            + std::inner_product(standard.begin(), standard.end(), x[phase].begin(), 0.0)
        );
    }
    return result;
}

std::vector<double> CompiledSystem::classified_residuals(
    const PhasePoint& prepared, const std::vector<double>& values
) const {
    std::vector<double> result = values;
    if (feasibility_cell() || values.size() != rows.size()) return result;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        const auto& item = rows[row];
        if (item.kind == RowKind::Chemical) {
            result[row] *= prepared.amount[item.phase] * prepared.composition[item.phase][item.index] / feed_total;
        } else if (item.kind == RowKind::Pressure) {
            result[row] *= prepared.amount[item.phase] / feed_total;
        }
    }
    return result;
}

double CompiledSystem::classified_residual_norm(const std::vector<double>& point) const {
    const auto prepared = prepare(point);
    double result = 0.0;
    for (const double value : classified_residuals(prepared, evaluate_at(prepared).values)) {
        if (!std::isfinite(value)) return std::numeric_limits<double>::infinity();
        result = std::max(result, std::abs(value));
    }
    return result;
}

std::vector<double> CompiledSystem::jacobian(const std::vector<double>& point) const {
    return jacobian_at(prepare(point));
}

std::vector<double> CompiledSystem::jacobian_at(const PhasePoint& prepared) const {
    core::Scope scope("equilibrium_jacobian");
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    const auto& compositions = prepared.composition;
    const auto& amount = prepared.amount;
    const auto& volumes = prepared.volume;
    const auto& packings = prepared.packing;
    const auto& derivatives = prepared.derivative;
    const auto& shape = prepared.shape;
    std::size_t offset = primal_size - free_intensive_count;
    const std::optional<std::size_t> temperature = problem.temperature.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const std::optional<std::size_t> pressure = problem.pressure.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const double T = prepared.temperature;
    const double P = prepared.pressure;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto& standard_temperature_derivative = prepared.standard.derivative[1];
    const auto& z = model->parameters().charge_number;
    std::vector<std::optional<std::size_t>> phi(problem.phases.size());
    std::size_t phi_offset = 0;
    for (const auto& row : rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = rank + phi_offset++;
    std::vector<double> result(rows.size() * n, 0.0);
    auto at = [&](std::size_t row, std::size_t column) -> double& {
        return result[row * n + column];
    };
    std::vector<std::vector<double>> columns;  // one first coefficient per coordinate serves every condition row
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const auto& row = rows[row_index];
        if (row.kind == RowKind::Chemical) {
            const auto& derivative = derivatives[row.phase];
            const std::size_t width = derivative.gradient.size();
            const auto& phase = problem.phases[row.phase];
            const std::size_t species = row.index;
            const double physical_v = derivative.hessian[(2 + species) * width + 1]
                - 1.0 / volumes[row.phase];
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const double physical_x = derivative.hessian[
                    (2 + species) * width + (2 + component)]
                    + (component == species ? 1.0 / compositions[row.phase][species] : 0.0);
                const double dv = packings[row.phase].component[component] / prepared.eta[row.phase];
                at(row_index, shape[row.phase].x[local]) =
                    physical_x * compositions[row.phase][component] + physical_v * dv;
            }
            at(row_index, shape[row.phase].volume) = -volumes[row.phase] * physical_v;
            if (temperature.has_value()) at(row_index, *temperature) =
                derivative.hessian[(2 + species) * width]
                + physical_v * packings[row.phase].temperature_derivative / prepared.eta[row.phase]
                + standard_temperature_derivative[species];
            if (pressure.has_value()) at(row_index, *pressure) =
                prepared.standard.pressure_derivative[species];
            for (std::size_t balance_row = 0; balance_row < rank; ++balance_row)
                at(row_index, primal_size + balance_row) = -B[balance_row][species];
            if (phi[row.phase].has_value()) at(row_index, primal_size + *phi[row.phase]) = -z[species];
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = derivatives[row.phase];
            const std::size_t width = derivative.gradient.size();
            const auto& phase = problem.phases[row.phase];
            const double physical_v = -derivative.gradient[1]
                - volumes[row.phase] * derivative.hessian[width + 1] - P / rt;
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const double dv = packings[row.phase].component[component] / prepared.eta[row.phase];
                at(row_index, shape[row.phase].x[local]) =
                    -volumes[row.phase] * derivative.hessian[width + 2 + component]
                        * compositions[row.phase][component] + physical_v * dv;
            }
            at(row_index, shape[row.phase].volume) = -volumes[row.phase] * physical_v;
            if (temperature.has_value()) at(row_index, *temperature) = -derivative.hessian[width]
                * volumes[row.phase] + P * volumes[row.phase]
                    / (GAS_CONSTANT_J_PER_MOL_K * T * T)
                + physical_v * packings[row.phase].temperature_derivative / prepared.eta[row.phase];
            if (pressure.has_value()) at(row_index, *pressure) = -volumes[row.phase] / rt;
        } else if (row.kind == RowKind::Moiety) {
            for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
                for (std::size_t local = 0; local < problem.phases[phase].support.size(); ++local) {
                    const std::size_t species = problem.phases[phase].support[local];
                    at(row_index, shape[phase].x[local]) = amount[phase] * B[row.index][species]
                        * compositions[phase][species] / feed_total;
                }
                if (shape[phase].amount.has_value()) {
                    for (const std::size_t species : problem.phases[phase].support)
                        at(row_index, *shape[phase].amount) += B[row.index][species]
                            * compositions[phase][species] / feed_total;
                }
            }
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, shape[row.phase].x[local]) =
                    z[problem.phases[row.phase].support[local]]
                    * compositions[row.phase][problem.phases[row.phase].support[local]];
        } else if (condition_row(row.kind)) {
            std::vector<std::vector<double>> state = {point, std::vector<double>(n, 0.0)};
            const std::vector<double> control(2 + n_feed.size() + model->active_parameters().size(), 0.0);
            for (std::size_t c = columns.size(); c < n; ++c) {
                state[1][c] = 1.0;
                columns.push_back(row_series(prepared, state, control, 1).rows[1]);
                state[1][c] = 0.0;
            }
            for (std::size_t c = 0; c < n; ++c) at(row_index, c) = columns[c][row_index];
        } else {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, shape[row.phase].x[local]) =
                    compositions[row.phase][problem.phases[row.phase].support[local]];
        }
    }
    scope.probe("rows", static_cast<double>(rows.size()));
    scope.probe("columns", static_cast<double>(n));
    return result;
}

std::vector<double> CompiledSystem::lever_jacobian(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_lever_jacobian");
    if (!reference_diagnostics.empty()) {
        throw std::invalid_argument(
            "reference_unavailable: source reference T/P/parameter lever derivatives are not admitted");
    }
    auto prepared = prepare(point);
    prepared.point.resize(std::max(prepared.point.size(), unknown_count()), 0.0);
    const std::size_t width = 2 + model->parameters().component_count() + model->active_parameters().size();
    std::vector<double> result(rows.size() * width), control(width);
    std::vector<std::vector<double>> state = {prepared.point, std::vector<double>(prepared.point.size())};
    for (std::size_t column = 0; column < width; ++column) {
        control[column] = 1.0;
        // The callable retains free T/P columns: drive their state coordinate rather than the ignored control.
        if (column == 0 && problem.temperature.free)
            state[1][primal_size - free_intensive_count] = 1.0;
        if (column == 1 && problem.pressure.free)
            state[1][primal_size - 1] = 1.0;
        const auto series = row_series(prepared, state, control, 1);
        for (std::size_t row = 0; row < rows.size(); ++row)
            result[row * width + column] = series.rows[1][row];
        control[column] = 0.0;
        std::fill(state[1].begin(), state[1].end(), 0.0);
    }
    return result;
}

std::vector<double> CompiledSystem::objective_gradient(const std::vector<double>& point) const {
    return objective_gradient_at(prepare(point));
}

std::vector<double> CompiledSystem::objective_gradient_at(const PhasePoint& prepared) const {
    core::Scope scope("equilibrium_objective_gradient");
    const auto& point = prepared.point;
    std::vector<double> result(point.size(), 0.0);
    if (feasibility_cell()) return result;
    const std::size_t components = model->parameters().component_count();
    const double T = prepared.temperature;
    const double P = prepared.pressure;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto& standard_value = prepared.standard.derivative[0];
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        const auto& layout = prepared.shape[phase];
        const auto& x = prepared.composition[phase];
        const double amount = prepared.amount[phase];
        const auto& packing = prepared.packing[phase];
        const double eta = prepared.eta[phase];
        const double volume = prepared.volume[phase];
        const auto& derivative = prepared.derivative[phase];
        double standard = 0.0;
        double ideal = ideal_mechanical_helmholtz(x, volume, 1.0, false);
        for (std::size_t species = 0; species < components; ++species) {
            standard += standard_value[species] * x[species];
        }
        const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        const double bracket = derivative.value + ideal + P * volume / rt + standard;
        const double physical_v = derivative.gradient[1] - sum_x / volume + P / rt;
        for (std::size_t local = 0; local < spec.support.size(); ++local) {
            const std::size_t species = spec.support[local];
            const double chemical =
                derivative.gradient[2 + species] + std::log(std::max(
                    x[species] / volume, std::numeric_limits<double>::min()))
                + standard_value[species];
            result[layout.x[local]] = amount * (
                chemical * x[species] + physical_v * packing.component[species] / eta);
        }
        result[layout.volume] = -amount * volume * physical_v;
        result[*layout.amount] = bracket;
    }
    return result;
}

std::vector<double> CompiledSystem::hessian(
    const std::vector<double>& point,
    const std::vector<double>& constraint_multipliers
) const {
    return hessian_at(prepare(point), constraint_multipliers);
}

std::vector<double> CompiledSystem::hessian_at(
    const PhasePoint& prepared, const std::vector<double>& constraint_multipliers
) const {
    core::Scope scope("equilibrium_hessian");
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (feasibility_cell()) return constraint_hessian_at(prepared, constraint_multipliers);
    const double P = prepared.pressure;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * prepared.temperature;
    const auto& standard_value = prepared.standard.derivative[0];
    auto put = [&](std::size_t row, std::size_t column, double value) {
        result[row * n + column] = value;
        result[column * n + row] = value;
    };
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        const auto& layout = prepared.shape[phase];
        const auto& x = prepared.composition[phase];
        const std::size_t volume = layout.volume;
        const double phase_amount_value = prepared.amount[phase];
        const auto& packing = prepared.packing[phase];
        const double eta = prepared.eta[phase];
        const double physical_volume = prepared.volume[phase];
        const auto& derivative = prepared.derivative[phase];
        const std::size_t width = derivative.gradient.size();
        const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        const double physical_v_gradient = phase_amount_value * (
            derivative.gradient[1] - sum_x / physical_volume + P / rt);
        const double physical_vv = phase_amount_value * (
            derivative.hessian[width + 1] + sum_x / (physical_volume * physical_volume));
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const std::size_t si = spec.support[i];
            const double dv_i = packing.component[si] / eta;
            const double physical_x_gradient = phase_amount_value * (
                derivative.gradient[2 + si]
                + std::log(std::max(x[si] / physical_volume,
                                    std::numeric_limits<double>::min()))
                + standard_value[si]);
            for (std::size_t j = 0; j < spec.support.size(); ++j) {
                const std::size_t sj = spec.support[j];
                const double dv_j = packing.component[sj] / eta;
                const double h_xx = phase_amount_value * (
                    derivative.hessian[(2 + si) * width + (2 + sj)]
                    + (si == sj ? 1.0 / x[si] : 0.0));
                const double h_xv_i = phase_amount_value * (
                    derivative.hessian[(2 + si) * width + 1] - 1.0 / physical_volume);
                const double h_xv_j = phase_amount_value * (
                    derivative.hessian[(2 + sj) * width + 1] - 1.0 / physical_volume);
                put(layout.x[i], layout.x[j],
                    x[si] * x[sj] * h_xx + x[si] * h_xv_i * dv_j
                    + dv_i * h_xv_j * x[sj] + dv_i * physical_vv * dv_j
                    + (si == sj ? physical_x_gradient * x[si]
                                      + physical_v_gradient * dv_i : 0.0));
            }
            const double h_xv = phase_amount_value * (
                derivative.hessian[(2 + si) * width + 1] - 1.0 / physical_volume);
            put(layout.x[i], volume,
                -physical_volume * (x[si] * h_xv + dv_i * physical_vv)
                - physical_v_gradient * dv_i);
        }
        put(volume, volume, physical_volume * physical_volume * physical_vv
            + physical_v_gradient * physical_volume);
        if (layout.amount.has_value()) {
            for (std::size_t i = 0; i < spec.support.size(); ++i) {
                const std::size_t species = spec.support[i];
                const double mu = derivative.gradient[2 + species]
                    + std::log(std::max(x[species] / physical_volume,
                                        std::numeric_limits<double>::min()))
                    + standard_value[species];
                const double pressure_gradient = derivative.gradient[1]
                    - sum_x / physical_volume + P / rt;
                put(*layout.amount, layout.x[i], mu * x[species]
                    + pressure_gradient * packing.component[species] / eta);
            }
            put(*layout.amount, volume, -physical_volume * (
                derivative.gradient[1] - sum_x / physical_volume + P / rt));
        }
    }
    if (!constraint_multipliers.empty()) {
        const auto weighted = constraint_hessian_at(prepared, constraint_multipliers);
        for (std::size_t index = 0; index < result.size(); ++index) result[index] += weighted[index];
    }
    scope.probe("entries", static_cast<double>(result.size()));
    return result;
}

RowSeries CompiledSystem::row_series(
    const PhasePoint& prepared, const std::vector<std::vector<double>>& state,
    const std::vector<double>& control, std::size_t order
) const {
    core::Scope scope("equilibrium_row_series");
    const std::size_t n = prepared.point.size();
    const std::size_t components = model->parameters().component_count();
    const std::size_t width = 2 + components + model->active_parameters().size();
    if (state.size() != order + 1 ||
        (control.size() != width && control.size() != width + 4 * problem.reactions.size() &&
         control.size() != width + 4 * problem.reactions.size() + 2) ||
        order > (problem.critical ? 2u : 3u))  // a third critical coefficient needs E6
        throw std::invalid_argument("row series shape");
    for (const auto& item : state)
        if (item.size() != n) throw std::invalid_argument("row series state");
    const bool pressure_jets = problem.pressure.free || control[1] != 0.0;
    const bool temperature_jets = problem.temperature.free || control[0] != 0.0;
    const std::vector<double> active_control(control.begin() + static_cast<std::ptrdiff_t>(2 + components), control.begin() + static_cast<std::ptrdiff_t>(width));
    const bool parameter_jets = std::any_of(active_control.begin(), active_control.end(),
        [](double value) { return value != 0.0; });
    if (!reference_diagnostics.empty() &&
        (order > 2 && (parameter_jets || pressure_jets || temperature_jets))) {
        throw std::invalid_argument(
            "reference_unavailable: row series pressure/T/parameter order is unsupported");
    }
    // The prepared term carries the current free-T/P jets. Parameter jets enter the rows only via
    // sum_a ds/dtheta_a (theta_a - theta_a,center): jets along every active parameter serve every direction
    // exactly. If any slope is refused, each direction resolves its own, leaving sibling directions intact.
    const bool reuse = !reference_diagnostics.empty() && control[0] == 0.0 && control[1] == 0.0;
    const auto parameter_cache_ready = [&] {
        return prepared.parameter_standard.has_value()
            && !(order > 1 && pressure_jets && prepared.parameter_standard->pressure_parameter_derivative.empty())
            && !(order > 1 && temperature_jets && prepared.parameter_standard->temperature_parameter_derivative.empty());
    };
    if (reuse && parameter_jets && !parameter_cache_ready()) {
        try {
            prepared.parameter_standard = standard_term(*this, prepared.temperature, prepared.pressure, 2,
                pressure_jets, temperature_jets, std::vector<double>(active_control.size(), 1.0));
        } catch (const std::invalid_argument&) {}
    }
    const auto& standard = reference_diagnostics.empty()
        ? standard_term(*this, prepared.temperature, prepared.pressure, order, false)
        : !reuse || (parameter_jets && !parameter_cache_ready())
            ? standard_term(*this, prepared.temperature, prepared.pressure, order, pressure_jets, temperature_jets, active_control)
        : parameter_jets ? *prepared.parameter_standard : prepared.standard;
    auto geometry = record_geometry(*this, prepared, standard);
    const std::size_t stride = width + 2 + (has_caloric_condition() ? 3 : 0) + (problem.critical ? 2 * components : 0);
    std::vector<std::vector<double>> output(order + 1);
    for (std::size_t k = 0; k <= order; ++k) {
        std::vector<double> input(n + width + (has_caloric_condition() ? 2 : 0), 0.0);
        std::copy(state[k].begin(), state[k].end(), input.begin());
        if (k == 0) {
            input[n] = prepared.temperature;
            input[n + 1] = prepared.pressure;
            std::copy(n_feed.begin(), n_feed.end(), input.begin() + static_cast<std::ptrdiff_t>(n + 2));
            const auto values = model->active_parameter_values();
            std::copy(values.begin(), values.end(), input.begin() + static_cast<std::ptrdiff_t>(n + 2 + components));
        } else if (k == 1) {
            std::copy_n(control.begin(), width, input.begin() + static_cast<std::ptrdiff_t>(n));
        }
        if (has_caloric_condition()) {
            if (k == 0) {
                input[n + width] = problem.enthalpy_j.value_or(0.0);
                input[n + width + 1] = problem.entropy_j_per_k.value_or(0.0);
            } else if (k == 1 && control.size() == width + 4 * problem.reactions.size() + 2) {
                input[n + width] = control[control.size() - 2];
                input[n + width + 1] = control.back();
            }
        }
        output[k] = geometry.Forward(k, input);
    }
    RowSeries result;
    result.rows.assign(order + 1, std::vector<double>(rows.size(), 0.0));
    result.temperature.assign(order + 1, 0.0);
    result.pressure.assign(order + 1, 0.0);
    result.phase.resize(problem.phases.size());
    const std::size_t row_offset = problem.phases.size() * stride;
    const std::size_t temperature_index = primal_size - free_intensive_count;
    for (std::size_t k = 0; k <= order; ++k) {
        for (std::size_t row = 0; row < rows.size(); ++row) result.rows[k][row] = output[k][row_offset + row];
        result.temperature[k] = k == 0 ? prepared.temperature
            : problem.temperature.free ? state[k][temperature_index] : k == 1 ? control[0] : 0.0;
        result.pressure[k] = k == 0 ? prepared.pressure
            : problem.pressure.free ? state[k][primal_size - 1] : k == 1 ? control[1] : 0.0;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            auto& series = result.phase[phase];
            const auto begin = output[k].begin() + static_cast<std::ptrdiff_t>(phase * stride);
            series.point.emplace_back(begin, begin + static_cast<std::ptrdiff_t>(width));
            series.amount.push_back(begin[static_cast<std::ptrdiff_t>(width)]);
            series.packing.push_back(begin[static_cast<std::ptrdiff_t>(width + 1)]);
        }
    }
    if (order > 0 && control.size() > width && !problem.reactions.empty()) {
        using AD = CppAD::AD<double>;
        std::vector<AD> q = {prepared.temperature, 0.0}, functions;
        CppAD::Independent(q);
        for (const auto &reaction : problem.reactions) {
            const auto r = std::size_t(&reaction - problem.reactions.data());
            const auto start = width + 4 * r;
            const double aa = control[start], bb = control[start + 1], cc = control[start + 2],
                         dd = control[start + 3];
            if (!reaction.correlation && (bb != 0.0 || cc != 0.0 || dd != 0.0))
                throw std::invalid_argument("b/c/d action needs a reaction law");
            const ReactionLogPolynomial law{aa, bb, cc, dd,
                reaction.correlation ? reaction.correlation->reference_temperature : 1.0, 0, 0, {}};
            const AD L = reaction_log_value(law, q[0]);
            const AD LT = reaction_log_value(law, q[0], 1);
            const AD LTT = reaction_log_value(law, q[0], 2);
            functions.push_back(q[1] * L);
            functions.push_back(q[1] * LT);
            functions.push_back(q[1] * LTT);
        }
        CppAD::ADFun<double> source(q, functions);
        for (auto &family : result.reaction_shift)
            family.assign(problem.reactions.size(), std::vector<double>(order + 1, 0.0));
        for (std::size_t k = 0; k <= order; ++k) {
            const auto values =
                source.Forward(k, std::vector<double>{result.temperature[k], k == 1 ? 1.0 : 0.0});
            std::vector<double> targets(problem.reactions.size());
            for (std::size_t r = 0; r < targets.size(); ++r) {
                for (std::size_t j = 0; j < 3; ++j)
                    result.reaction_shift[j][r][k] = values[3 * r + j];
                targets[r] = -values[3 * r];
            }
            const auto species =
                species_standard_slopes(prepared.temperature, components, problem.reactions, targets);
            for (std::size_t row = 0; row < rows.size(); ++row)
                if (rows[row].kind == RowKind::Chemical)
                    result.rows[k][row] += species[rows[row].index];
        }
    }
    // Residual gradient coefficients: the Hessian action of each coefficient plus the model's
    // third and fourth contractions of the lower coefficients (Faa di Bruno through order 3).
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        auto& series = result.phase[phase];
        const auto& derivative = prepared.derivative[phase];
        series.gradient.assign(order + 1, std::vector<double>(width, 0.0));
        series.gradient[0] = derivative.gradient;
        std::vector<std::size_t> needed = {1};
        if (has_caloric_condition()) needed.push_back(0);
        for (const std::size_t species : problem.phases[phase].support) needed.push_back(2 + species);
        const auto contraction = [&](std::vector<std::vector<double>> directions, std::optional<std::size_t> component) {
            if (component) {
                std::vector<double> unit(width, 0.0);
                unit[*component] = 1.0;
                directions.push_back(unit);
            }
            const auto energy = model->residual_energy_action(
                prepared.temperature, prepared.volume[phase], prepared.composition[phase], directions);
            result.maximum_energy_order = std::max(result.maximum_energy_order, directions.size());
            auto& total = result.energy_diagnostics;
            total.forward_sweeps += energy.diagnostics.forward_sweeps;
            total.directional_paths += energy.diagnostics.directional_paths;
            total.site_solves += energy.diagnostics.site_solves;
            total.site_factorizations += energy.diagnostics.site_factorizations;
            total.recordings += energy.diagnostics.recordings;
            total.site_rank = energy.diagnostics.site_rank;
            total.site_pivot_ratio = std::max(total.site_pivot_ratio, energy.diagnostics.site_pivot_ratio);
            total.association_residual = std::max(total.association_residual, energy.diagnostics.association_residual);
            if (energy.status != PropertyActionStatus::Available && result.energy == PropertyActionStatus::Available)
                result.energy = energy.status;
            return energy.action.value_or(std::numeric_limits<double>::quiet_NaN());
        };
        for (std::size_t k = 1; k <= order; ++k) {
            for (std::size_t i = 0; i < width; ++i)
                for (std::size_t j = 0; j < width; ++j)
                    series.gradient[k][i] += derivative.hessian[i * width + j] * series.point[k][j];
            if (k == 2) {
                for (const std::size_t j : needed)
                    series.gradient[2][j] += 0.5 * contraction({series.point[1], series.point[1]}, j);
            } else if (k == 3) {
                for (const std::size_t j : needed) {
                    series.gradient[3][j] += contraction({series.point[1], series.point[2]}, j)
                        + contraction({series.point[1], series.point[1], series.point[1]}, j) / 6.0;
                }
            }
        }
        if (has_caloric_condition()) {
            std::vector<double> energy(order + 1, 0.0);
            energy[0] = derivative.value;
            for (std::size_t k = 1; k <= order; ++k)
                energy[k] = std::inner_product(derivative.gradient.begin(), derivative.gradient.end(), series.point[k].begin(), 0.0);
            for (std::size_t k = 2; k <= order; ++k)
                for (std::size_t i = 0; i < width; ++i)
                    for (std::size_t j = 0; j < width; ++j)
                        energy[k] += derivative.hessian[i * width + j] * series.point[1][i]
                            * series.point[k - 1][j] * (k == 2 ? 0.5 : 1.0);
            if (order == 3)
                energy[3] += contraction({series.point[1], series.point[1], series.point[1]}, std::nullopt) / 6.0;
            for (std::size_t k = 0; k <= order; ++k)
                for (std::size_t j = 0; j <= k; ++j) {
                    const auto offset = phase * stride + width + 2;
                    result.rows[k].back() += output[j][offset] * energy[k - j]
                        + output[j][offset + 1] * series.gradient[k - j][0]
                        + output[j][offset + 2] * series.gradient[k - j][1];
                }
        }
        if (problem.critical) {
            // Residual parts of EqID critical_equations through order 2: d_i [F_nn(p) w]_i (E2-E4) and
            // F_nnn(p)[w,w,w] (E3-E5) along p(t), w(t), d(t); a zero direction contributes exactly zero.
            const std::size_t offset = phase * stride + width + 2, mode = rows.size() - components - 2;
            std::vector<std::vector<double>> d(order + 1, std::vector<double>(components)),
                w(order + 1, std::vector<double>(width, 0.0));
            for (std::size_t k = 0; k <= order; ++k)
                for (std::size_t i = 0; i < components; ++i) {
                    d[k][i] = output[k][offset + 2 * i];
                    w[k][2 + i] = output[k][offset + 2 * i + 1];
                }
            const auto& p = series.point;
            const auto F = [&](const std::vector<std::vector<double>>& v, std::optional<std::size_t> j = std::nullopt) {
                const bool zero = std::any_of(v.begin(), v.end(), [](const std::vector<double>& a) {
                    return std::all_of(a.begin(), a.end(), [](double x) { return x == 0.0; }); });
                return zero ? 0.0 : contraction(v, j);
            };
            for (std::size_t i = 0; i < components; ++i) {
                std::vector<double> g(order + 1, 0.0);
                for (std::size_t k = 0; k <= order; ++k)
                    for (std::size_t j = 0; j < width; ++j) g[k] += derivative.hessian[(2 + i) * width + j] * w[k][j];
                if (order > 0) g[1] += F({p[1], w[0]}, 2 + i);
                if (order > 1) g[2] += F({p[1], w[1]}, 2 + i) + F({p[2], w[0]}, 2 + i) + 0.5 * F({p[1], p[1], w[0]}, 2 + i);
                for (std::size_t k = 0; k <= order; ++k)
                    for (std::size_t j = 0; j <= k; ++j) result.rows[k][mode + i] += d[j][i] * g[k - j];
            }
            std::array<double, 3> c{F({w[0], w[0], w[0]}), 0.0, 0.0};
            if (order > 0) c[1] = 3.0 * F({w[0], w[0], w[1]}) + F({w[0], w[0], w[0], p[1]});
            if (order > 1) c[2] = 3.0 * F({w[0], w[0], w[2]}) + 3.0 * F({w[0], w[1], w[1]}) + 3.0 * F({w[0], w[0], w[1], p[1]})
                + F({w[0], w[0], w[0], p[2]}) + 0.5 * F({w[0], w[0], w[0], p[1], p[1]});
            for (std::size_t k = 0; k <= order; ++k) result.rows[k][mode + components + 1] += c[k];
        }
        for (std::size_t k = 0; k <= order; ++k) {
            for (std::size_t row = 0; row < rows.size(); ++row) {
                if (rows[row].phase != phase) continue;
                if (rows[row].kind == RowKind::Chemical) {
                    result.rows[k][row] += series.gradient[k][2 + rows[row].index];
                } else if (rows[row].kind == RowKind::Pressure) {
                    for (std::size_t j = 0; j <= k; ++j)
                        result.rows[k][row] -= series.point[j][1] * series.gradient[k - j][1];
                }
            }
        }
    }
    scope.probe("order", static_cast<double>(order));
    scope.probe("energy_order", static_cast<double>(result.maximum_energy_order));
    return result;
}

std::vector<double> CompiledSystem::constraint_hessian_at(
    const PhasePoint& prepared, const std::vector<double>& constraint_multipliers
) const {
    if (feasibility_cell()) return pinned_hessian(*this, prepared, constraint_multipliers);
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (constraint_multipliers.empty()) return result;
    std::size_t multiplier = 0;
    for (const auto& row : rows) {
        if (row.kind == RowKind::Chemical || row.kind == RowKind::Pressure) continue;
        if (multiplier >= constraint_multipliers.size()) {
            throw std::invalid_argument("constraint multiplier count");
        }
        const double weight = constraint_multipliers[multiplier++];
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            if (row.kind != RowKind::Moiety && row.phase != phase) continue;
            const auto& spec = problem.phases[phase];
            const auto& layout = prepared.shape[phase];
            const std::size_t amount = *layout.amount;
            for (std::size_t local = 0; local < spec.support.size(); ++local) {
                const std::size_t species = spec.support[local];
                const std::size_t coordinate = layout.x[local];
                const double x = prepared.composition[phase][species];
                double coefficient = 1.0;
                if (row.kind == RowKind::Electro)
                    coefficient = model->parameters().charge_number[species];
                if (row.kind == RowKind::Moiety)
                    coefficient = B[row.index][species] / feed_total;
                result[coordinate * n + coordinate] += weight * coefficient * x
                    * (row.kind == RowKind::Moiety ? point[amount] : 1.0);
                if (row.kind == RowKind::Moiety) {
                    result[coordinate * n + amount] += weight * coefficient * x;
                    result[amount * n + coordinate] += weight * coefficient * x;
                }
            }
        }
    }
    if (multiplier != constraint_multipliers.size())
        throw std::invalid_argument("constraint multiplier count");
    return result;
}

}  // namespace epcsaft::equilibrium
