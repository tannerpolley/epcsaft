// CompiledSystem coordinate labels, row values, uncertainty, Jacobians and Hessians, with the standard and
// reference terms they read.
#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"
#include "core/trace.hpp"
#include "core/linear_error.hpp"
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

// EqIDs ideal_mole_fraction_reaction, ideal_molarity_reaction and ideal_reaction_amount_coordinates.
CppAD::ADFun<double> ideal_rows(const CompiledSystem& c, const PhasePoint& p) {
    using AD = CppAD::AD<double>;
    std::vector<AD> y(p.point.begin(), p.point.end()), rows;
    CppAD::Independent(y);
    const double shift = *std::max_element(p.point.begin(), p.point.end());
    AD total = 0.0;
    for (const auto& value : y) total += CppAD::exp(value - shift);
    const AD log_total = shift + CppAD::log(total);
    const auto& support = c.problem.phases[0].support;
    for (std::size_t r = 0; r < c.problem.reactions.size(); ++r) {
        AD value = -c.ideal_log_k[r];
        for (std::size_t j = 0; j < y.size(); ++j)
            value += c.problem.reactions[r].stoichiometry[support[j]] * (c.ideal->activity_basis == "molarity"
                ? c.ideal->molarity_log_activity(AD(y[j] + std::log(c.feed_total)), support[j]) : AD(y[j] - log_total));
        rows.push_back(value);
    }
    for (std::size_t k = 0; k < c.rank; ++k) {
        AD value = -c.balance[k] / c.feed_total;
        for (std::size_t j = 0; j < y.size(); ++j) value += c.B[k][support[j]] * CppAD::exp(y[j]);
        rows.push_back(value);
    }
    return CppAD::ADFun<double>(y, rows);
}

// Compare supplied species thermochemistry with every reaction target on the density
// basis: [t, T t', T^2 t''] against [g, g - h, -T cp] at the same J/mol scale.
void check_reference_consistency(
    const ThermochemistryRecord& record, double temperature,
    const std::vector<ReactionSpec>& reactions
) {
    const auto targets = reaction_target_series(temperature, reactions, 2, &record);
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        if (reactions[row].uses_species_thermochemistry()) continue;
        const auto species = reaction_thermochemistry_density(
            record, reactions[row].stoichiometry, temperature);
        if (species.status != ThermochemistryStatus::Available) continue;
        const double t0 = targets[0][row];
        const double t1 = temperature * targets[1][row];
        const double t2 = temperature * temperature * targets[2][row];
        const double d0 = t0 - species.gibbs;
        const double d1 = t1 - (species.gibbs - species.enthalpy);
        const double d2 = t2 + temperature * species.heat_capacity;
        const double scale = std::sqrt(std::numeric_limits<double>::epsilon())
            * std::max(1.0, std::sqrt(t0 * t0 + t1 * t1 + t2 * t2));
        if (std::sqrt(d0 * d0 + d1 * d1 + d2 * d2) > scale) {
            throw std::invalid_argument("reference inconsistency in reaction "
                                        + reaction_name(reactions[row], row));
        }
    }
}

// The compiled pure-solvent liquid density starts every later reference ladder of the problem.
std::optional<double> reference_solvent_anchor(const std::vector<NeutralReferenceReactionDiagnostic>& d) {
    return d.empty() || !(d.front().solvent_density > 0.0) ? std::nullopt : std::optional(d.front().solvent_density);
}

// The source part of ln K^rho = ln_k + scale + C of a reaction or solid law at T: ln K_src with its T slope and curvature
// from a correlation, or the log of a scalar K, and the common-molality scale when its source basis needs one.
SolidLogK source_log_k(const HelmholtzModel& model, const std::optional<NeutralReferenceDeclaration>& declaration,
                       const ReactionSpec& reaction, double temperature) {
    SolidLogK result;
    if (reaction.correlation) {
        const auto law = evaluate_reaction_properties(*reaction.correlation, temperature);
        result.ln_k = law.ln_k, result.slope = law.derivative, result.curvature = law.second_derivative;
    } else result.ln_k = std::log(*reaction.equilibrium_constant);
    if (reaction.reference && reaction.reference->source_basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution)
        result.scale = reaction_reference_log_scale(model, *declaration, reaction.stoichiometry);
    return result;
}

static Eigen::VectorXd map_error(const Eigen::MatrixXd& n, const Eigen::MatrixXd& w, const Eigen::VectorXd& b,
    const Eigen::VectorXd& s, const Eigen::VectorXd& eb) {
    const double eps = std::numeric_limits<double>::epsilon(), gamma = n.cols() * eps / (1 - n.cols() * eps);
    const Eigen::VectorXd defect = core::compensated_defect(n, s, b), application = core::compensated_defect(w, b, s);
    return w.cwiseAbs() * (eb + gamma * (b.cwiseAbs() + n.cwiseAbs() * s.cwiseAbs()))
        + (w * defect).cwiseAbs() + application.cwiseAbs() + eps * (s.cwiseAbs() + w.cwiseAbs() * b.cwiseAbs());
}

// The ln K^rho_s(T, P0) = ln K_src + scale + C_s columns an action moves, from the solid_log_k owner: dC/dtheta_a for every
// moved coordinate; with T jets of order two L'' + C_TT and d2C/dT dtheta_a; d2C/dtheta_a dtheta_b over the pairs. A
// requested column the reference owner did not admit refuses with its cause.
std::vector<SolidLogK> solid_columns(const CompiledSystem& c, double T, double P, std::size_t order, bool temperature_jets,
                                     const std::vector<std::size_t>& parameters, ReferenceJets& jets) {
    const auto& solids = c.problem.solids;
    const std::size_t thermal = order > 1 && temperature_jets ? 2 : 0;
    if (solids.empty() || (parameters.empty() && !thermal)) return {};
    jets.temperature_order = thermal;
    const auto& pairs = jets.pairs;
    auto laws = solid_log_k(*c.model, T, P, c.problem.source_reference, solids, jets);
    for (std::size_t s = 0; s < solids.size(); ++s) {
        const auto need = [&](std::size_t a, std::optional<std::size_t> p = {}, std::optional<std::size_t> q = {}) {
            if (const auto column = laws[s].reference.column(a, 0, p, q); column.state != ReferenceColumnState::Available)
                throw ReferenceUnavailable("solid " + solids[s].name + " " + *column.cause);
        };
        if (thermal) need(2);
        for (const auto a : parameters) {
            need(0, a);
            if (thermal) need(1, a);
        }
        for (std::size_t i = 0; i < pairs.size(); ++i)
            for (std::size_t j = i; j < pairs.size(); ++j) need(0, pairs[i], pairs[j]);
    }
    return laws;
}

StandardTerm standard_term(
    const CompiledSystem& compiled, double temperature, double pressure, std::size_t order,
    bool pressure_jets, bool temperature_jets = false, const std::vector<std::size_t>& parameters = {},
    StandardTerm::NumericalError* error = nullptr, const std::vector<std::size_t>& pairs = {},
    const std::vector<std::size_t>& parameter_union = {}, const StandardTerm* consumed_center = nullptr
) {
    const auto& problem = compiled.problem;
    const std::size_t components = compiled.model->parameters().component_count();
    const std::size_t thermal = temperature_jets ? order : 0;
    ReferenceJets jets{thermal, pressure_jets ? order : 0, parameter_union.empty() ? parameters : parameter_union, 0, pairs};
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
        std::vector<std::vector<double>> inverse;
        // Fixed-T/P references do not change with phase coordinates or feed directions.
        if (consumed_center || problem.pressure.free || pressure_jets || temperature_jets || !parameters.empty()) {
            jets.reference_order = !parameters.empty() && order > 1 && (temperature_jets || pressure_jets) ? 2 : 0;
            reactions = central_reference_reactions(
                *compiled.model, problem, temperature, pressure, diagnostics, jets,
                consumed_center ? &consumed_center->reference_diagnostics : &compiled.reference_diagnostics, consumed_center != nullptr);
            gibbs = standard_state_series(temperature, components, reactions, {}, 0, nullptr, error ? &inverse : nullptr).front();
        }
        StandardTerm result;
        result.derivative.assign(order + 1, std::vector<double>(components, 0.0));
        result.pressure_derivative.assign(components, 0.0);
        result.pressure_second_derivative.assign(components, 0.0);
        result.reference_diagnostics = diagnostics;
        const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
        for (std::size_t species = 0; species < components; ++species)
            result.derivative[0][species] = gibbs[species] / rt;
        if (consumed_center) result.derivative[0] = consumed_center->derivative[0];
        // One slope of s_r = g_r/(RT) per row, mapped to species by the same linear minimum-norm
        // owner as s, factored once for every column. A referenced row has s_r = -ln K_rho, so its slope
        // is -(source slope + slope of C); an unadmitted slope of C refuses the direction. Only pure T
        // slopes move the source.
        std::optional<SpeciesSlopeMap> slopes;
        const auto species_slope = [&](std::vector<double> values, const std::string& name, std::size_t a,
                                       std::size_t b = 0, std::optional<std::size_t> parameter = std::nullopt,
                                       std::optional<std::size_t> second = std::nullopt) {
            const std::size_t source_order = !parameter && b == 0 ? a : 0;
            for (const auto& diagnostic : diagnostics) {
                const std::size_t row = diagnostic.reaction_index;
                const ReferenceColumn slope = diagnostic.column(a, b, parameter, second);
                if (slope.state != ReferenceColumnState::Available)
                    throw ReferenceUnavailable("reaction " + reaction_name(problem.reactions[row], row)
                        + " C_" + std::to_string(a) + std::to_string(b)
                        + (parameter ? " coordinate " + std::to_string(*parameter) : "")
                        + (second ? "," + std::to_string(*second) : "") + " (" + name + "): " + *slope.cause);
                const auto& source = problem.reactions[row].correlation;
                const auto law = source ? evaluate_reaction_properties(*source, temperature) : ReactionProperties{};
                values[row] = -slope.value - (source_order == 1 ? law.derivative : source_order == 2 ? law.second_derivative : 0.0);
            }
            if (!slopes) slopes.emplace(temperature, components, reactions);
            return (*slopes)(values);
        };
        const std::vector<double> zero(reactions.size(), 0.0);
        // EqIDs: standard_state_transfer_temperature_(second_)derivative. Unreferenced rows keep
        // s' = (t' - t/T)/(RT) and s'' = t''/(RT) - 2 s'/T.
        std::vector<double> first(reactions.size(), 0.0), second(reactions.size(), 0.0);
        const auto targets = thermal || error ? reaction_target_series(temperature, reactions, thermal)
                                         : std::vector<std::vector<double>>{};
        for (std::size_t row = 0; thermal > 0 && row < reactions.size(); ++row) {
            first[row] = (targets[1][row] - targets[0][row] / temperature) / rt;
            if (thermal > 1) second[row] = targets[2][row] / rt - 2.0 * first[row] / temperature;
        }
        if (thermal > 0) result.derivative[1] = species_slope(first, "temperature derivative", 1);
        if (thermal > 1) result.derivative[2] = species_slope(second, "second temperature derivative", 2);
        // ln K_src and the molality scale are parameter free: ds_r/dtheta_a = -dC_r/dtheta_a
        // (EqID: standard_state_transfer_parameter_derivative).
        if (!parameters.empty()) {
            result.parameter_derivative.assign(compiled.model->active_parameters().size(), std::vector<double>(components, 0.0));
            if (order > 1 && temperature_jets)
                result.temperature_parameter_derivative.assign(compiled.model->active_parameters().size(), std::vector<double>(components, 0.0));
            if (order > 1 && pressure_jets)
                result.pressure_parameter_derivative.assign(compiled.model->active_parameters().size(), std::vector<double>(components, 0.0));
        }
        for (const auto a : parameters) {
            result.parameter_derivative[a] = species_slope(zero, "parameter derivative", 0, 0, a);
            if (order > 1 && temperature_jets)
                result.temperature_parameter_derivative[a] = species_slope(zero, "T parameter derivative", 1, 0, a);
            if (order > 1 && pressure_jets)
                result.pressure_parameter_derivative[a] = species_slope(zero, "P parameter derivative", 0, 1, a);
        }
        // d2s/dtheta_a dtheta_b = -W d2C/dtheta_a dtheta_b, stored symmetric over the requested pairs.
        result.pair_coordinates = pairs;
        result.parameter_pair_derivative.assign(pairs.size(), std::vector<std::vector<double>>(pairs.size()));
        for (std::size_t i = 0; i < pairs.size(); ++i)
            for (std::size_t j = i; j < pairs.size(); ++j)
                result.parameter_pair_derivative[i][j] = result.parameter_pair_derivative[j][i]
                    = species_slope(zero, "parameter pair derivative", 0, 0, pairs[i], pairs[j]);
        // EqIDs: neutral_reaction_reference_pressure_jet, standard_state_transfer_temperature_pressure_derivative.
        if (pressure_jets && order) result.pressure_derivative = species_slope(zero, "pressure derivative", 0, 1);
        if (pressure_jets && order > 1)
            result.pressure_second_derivative = species_slope(zero, "second pressure derivative", 0, 2);
        if (pressure_jets && thermal > 1)
            result.temperature_pressure_derivative = species_slope(zero, "temperature-pressure derivative", 1, 1);
        if (error) {
            Eigen::MatrixXd n(reactions.size(), components), w(components, reactions.size());
            for (std::size_t r = 0; r < reactions.size(); ++r) for (std::size_t i = 0; i < components; ++i) {
                n(r, i) = reactions[r].stoichiometry[i]; w(i, r) = inverse[i][r];
            }
            for (auto& diagnostic : diagnostics) {
                reference_errors(*compiled.model, temperature, problem.reactions[diagnostic.reaction_index].stoichiometry, diagnostic, order > 0);
                if (!diagnostic.numerical_error_cause.empty()) error->cause = diagnostic.numerical_error_cause;
            }
            for (std::size_t k = 0; k <= std::min(order, std::size_t(1)); ++k) {
                Eigen::VectorXd b(reactions.size()), eb(reactions.size()), species(components);
                for (std::size_t r = 0; r < reactions.size(); ++r) {
                    b(r) = k ? first[r] : targets[0][r] / rt;
                    eb(r) = std::numeric_limits<double>::epsilon() * std::abs(b(r));
                }
                for (const auto& diagnostic : diagnostics) {
                    const auto r = diagnostic.reaction_index;
                    const auto& reaction = problem.reactions[r];
                    const auto law = source_log_k(*compiled.model, problem.source_reference, reaction, temperature);
                    const double source = k ? law.slope.value_or(0.0) : law.ln_k;
                    const double scale = k ? 0.0 : law.scale;
                    const double c = k ? *diagnostic.terminal_contraction_temperature_derivative : *diagnostic.terminal_contraction;
                    b(r) = -source - scale - c;
                    eb(r) = diagnostic.numerical_error[k] + std::numeric_limits<double>::epsilon() * (std::abs(source) + std::abs(scale) + std::abs(c));
                }
                for (std::size_t i = 0; i < components; ++i) species(i) = result.derivative[k][i];
                const auto mapped = map_error(n, w, b, species, eb);
                error->absolute[k].assign(mapped.data(), mapped.data() + mapped.size());
            }
        }
        result.solid = solid_columns(compiled, temperature, pressure, order, temperature_jets, parameters, jets);
        return result;
    }
    const auto gibbs = standard_state_series(
        temperature, components, problem.reactions, problem.standard_gibbs_j_per_mol, order, compiled.model->thermochemistry() ? &*compiled.model->thermochemistry() : nullptr);
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
    result.solid = solid_columns(compiled, temperature, pressure, order, temperature_jets, parameters, jets);
    return result;
}

// PH/PS consumes complete records at every callback temperature.
void admit_calorics(const CompiledSystem& compiled, double temperature) {
    if (!compiled.has_caloric_condition()) return;
    const auto& p = compiled.problem;
    if (int(p.enthalpy_j.has_value()) + int(p.entropy_j_per_k.has_value()) + int(p.internal_energy_j.has_value()) != 1)
        throw std::invalid_argument("one extensive caloric condition required");
    if (!p.temperature.free || (p.pressure.free != p.internal_energy_j.has_value()) || compiled.pinned_count
        || !p.topology_declared || p.discovery)
        throw std::invalid_argument("PH/PS requires declared free amounts, free T and fixed P");
    if (!std::isfinite(p.internal_energy_j ? *p.internal_energy_j : p.enthalpy_j ? *p.enthalpy_j : *p.entropy_j_per_k))
        throw std::invalid_argument("nonfinite caloric target");
    const auto& record = compiled.model->thermochemistry();
    if (!record) throw std::invalid_argument("missing thermochemistry for caloric condition");
    if (compiled.anchored_calorics()) return;  // caloric_reference admits records, intervals and completion
    for (const auto& phase : p.phases)
        for (auto i : phase.support) {
            const auto& intervals = record->component_intervals[i];
            if (intervals.empty()) throw std::invalid_argument("missing caloric record: " + std::to_string(i));
            const auto selected = std::find_if(intervals.begin(), intervals.end(), [&](const auto& v) {
                return ideal_interval_contains(v, temperature); });
            if (selected == intervals.end())
                throw std::invalid_argument("caloric record outside interval: species " + std::to_string(i)
                    + " at " + std::to_string(temperature) + " K");
            if (p.internal_energy_j && (temperature == selected->lower || temperature == selected->upper))
                throw std::invalid_argument("caloric derivative at interval endpoint: species " + std::to_string(i));
        }
}

// The anchored g of a PH/PS cell at (T, P): every pure-T column through the row and curvature order (3 with a source
// reference, whose row series stop at order 2; else 4) must exist, else the callback is refused with its cause.
std::shared_ptr<const ReferenceCalorics> caloric_reference(const CompiledSystem& c, double T, double P,
                                                           bool pressure = false, bool parameters = false) {
    ReferenceJets jets;
    jets.temperature_order = c.reference_diagnostics.empty() ? 4 : 3;
    jets.pressure_order = pressure ? jets.temperature_order : 0;
    if (parameters) jets.request_every_parameter(*c.model);
    auto result = std::make_shared<ReferenceCalorics>(reference_calorics(c, T, P, jets));
    if (result->status != SolvedStateActionStatus::Available)
        throw std::invalid_argument("caloric reference at " + std::to_string(T) + " K: " + result->message);
    for (const auto& e : result->gibbs)
        if (!e.pressure_order && !e.parameter && e.column.status != SolvedStateActionStatus::Available)
            throw std::invalid_argument("caloric reference at " + std::to_string(T) + " K: " + e.column.message);
    return result;
}

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

std::vector<double> caloric_hessian(const CompiledSystem&, const PhasePoint&, double, bool = false, std::vector<double>* = nullptr);
std::vector<double> critical_hessian(const CompiledSystem&, const PhasePoint&, const std::vector<double>&);
CppAD::ADFun<double> record_geometry(const CompiledSystem&, const PhasePoint&, const StandardTerm&, bool = false, bool = false);
std::vector<double> controls(const CompiledSystem&, const PhasePoint&);

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
    if (std::all_of(multipliers.begin(), multipliers.end(), [](double value) { return value == 0.0; })) return result;
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
            : AD(compiled.phase_amount(problem.phases[phase]));
        const auto feed = compiled.feed_at(q);
        AD total = 0.0;
        for (const auto& value : feed) total += value;
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
                        * x[species] / (problem.loading ? total : AD(compiled.feed_total));
                if (problem.loading && phase == 0)
                    for (std::size_t species = 0; species < components; ++species)
                        simple -= multipliers[row] * compiled.B[descriptor.index][species] * feed[species] / total;
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
        const auto [status, triangle] = model.contracted_energy_hessian(temperature, values[1], amounts, contraction, directions);
        if (status != PropertyActionStatus::Available) throw std::runtime_error(std::string("third directional energy "
            "contraction unavailable: ") + solved_state_status_name(solved_state_status(status)));
        std::size_t entry = 0;
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t column = 0; column <= row; ++column) {
                const double action = triangle[entry++];
                const double value = mapped[row * n + column] + action
                    - pressure_weight * (directions[row][1] * volume_hessian_action[column]
                        + directions[column][1] * volume_hessian_action[row]);
                result[row * n + column] += value;
                if (row != column) result[column * n + row] += value;
            }
        }
    }
    if (problem.volume_m3) {
        auto geometry = record_geometry(compiled, prepared, prepared.standard);
        std::vector<double> weights(geometry.Range(), 0.0);
        for (std::size_t r = 0; r < compiled.rows.size(); ++r)
            if (compiled.rows[r].kind == RowKind::Volume)
                weights[geometry.Range() - compiled.rows.size() - 2 * compiled.solids.size() + r] = multipliers[r];
        const auto h = geometry.Hessian(controls(compiled, prepared), weights);
        const auto width = geometry.Domain();
        for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) result[i*n+j] += h[i*width+j];
    }
    if (compiled.has_caloric_condition()) {
        const auto caloric = caloric_hessian(compiled, prepared, multipliers.back());
        for (std::size_t i = 0; i < result.size(); ++i) result[i] += caloric[i];
    }
    if (problem.critical || problem.spinodal) {
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

// The record_geometry inputs at the prepared point: y, T, P, n_feed, active parameters, then H and S
// in a caloric cell.
std::vector<double> controls(const CompiledSystem& compiled, const PhasePoint& prepared) {
    std::vector<double> result = prepared.point;
    result.push_back(prepared.temperature);
    result.push_back(prepared.pressure);
    result.insert(result.end(), compiled.n_feed.begin(), compiled.n_feed.end());
    const auto theta = compiled.model->active_parameter_values();
    result.insert(result.end(), theta.begin(), theta.end());
    if (compiled.has_caloric_condition()) {
        result.push_back(compiled.problem.enthalpy_j.value_or(0.0));
        result.push_back(compiled.problem.entropy_j_per_k.value_or(0.0));
    }
    if (compiled.problem.volume_m3) result.push_back(*compiled.problem.volume_m3);
    if (compiled.problem.internal_energy_j) result.push_back(*compiled.problem.internal_energy_j);
    return result;
}

// The geometry map and every analytic row term as one CppAD function of (y, controls) with
// controls [T, P, n_feed, active parameters]: outputs per phase [T, v, x (C), active (A),
// N, eta], then the analytic part of each row (the residual gradient terms are added by the
// model's contractions in row_series).
CppAD::ADFun<double> record_geometry(
    const CompiledSystem& compiled, const PhasePoint& prepared, const StandardTerm& standard, bool entropy, bool plane
) {
    using AD = CppAD::AD<double>;
    const auto& problem = compiled.problem;
    const auto& model = *compiled.model;
    const auto& parameters = model.parameters();
    const std::size_t components = parameters.component_count();
    const std::size_t n = prepared.point.size();
    const auto center = controls(compiled, prepared);
    std::vector<AD> u(center.begin(), center.end());
    const std::vector<double> theta_center(center.begin() + static_cast<std::ptrdiff_t>(n + 2 + components),
        center.begin() + static_cast<std::ptrdiff_t>(n + compiled.control_width()));
    std::vector<double> chemical(components, 0.0);
    if (entropy) {
        const auto reference = reference_calorics(compiled, prepared.temperature, prepared.pressure, {});
        if (reference.status != SolvedStateActionStatus::Available) throw std::invalid_argument(reference.message);
        for (std::size_t p = 0; p < prepared.shape.size(); ++p)
            for (auto i : problem.phases[p].support)
                chemical[i] = std::log(prepared.composition[p][i] / prepared.volume[p])
                    + prepared.derivative[p].gradient[2+i] + reference.gibbs.front().column.values[i]
                    / (GAS_CONSTANT_J_PER_MOL_K * prepared.temperature);
    }
    CppAD::Independent(u);
    const AD T = problem.temperature.free ? u[compiled.primal_size - compiled.free_intensive_count] : u[n];
    const AD P = problem.pressure.free ? u[compiled.primal_size - 1] : u[n + 1];
    AD feed_total = AD(0.0);
    auto feed = compiled.feed_at(u);
    for (std::size_t i = 0; i < components; ++i) {
        feed[i] += u[n + 2 + i] - compiled.n_feed[i];
        feed_total += feed[i];
    }
    const std::vector<AD> theta(u.begin() + static_cast<std::ptrdiff_t>(n + 2 + components), u.begin() + static_cast<std::ptrdiff_t>(n + 2 + components + theta_center.size()));
    const auto tape_parameters = detail::make_tape_parameters(parameters, model.active_parameters(),
        prepared.parameter_actions ? theta : std::vector<AD>{});
    const auto sigma = segment_sigma(T, tape_parameters);
    const auto diameter = segment_diameter(T, tape_parameters, sigma);
    const auto& charges = parameters.charge_number;
    std::vector<AD> output;
    std::vector<AD> row_values(compiled.rows.size(), AD(0.0));
    AD material_potential = 0.0;
    const AD delta_T = T - prepared.temperature;
    const AD delta_P = P - prepared.pressure;
    // Adds the exact Taylor term of species g/RT in the free or driven temperature, pressure and parameters.
    const auto add_standard = [&](AD value, std::size_t species) {
        AD power = AD(1.0);
        for (std::size_t k = 0; k < standard.derivative.size(); ++k) {
            value += standard.derivative[k][species] * power;
            power *= delta_T / AD(static_cast<double>(k + 1));
        }
        value += standard.pressure_derivative[species] * delta_P
            + 0.5 * standard.pressure_second_derivative[species] * delta_P * delta_P;
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
        const auto& pairs = standard.pair_coordinates;
        for (std::size_t i = 0; i < pairs.size(); ++i)
            for (std::size_t j = 0; j < pairs.size(); ++j)
                value += 0.5 * standard.parameter_pair_derivative[i][j][species]
                    * (theta[pairs[i]] - theta_center[pairs[i]]) * (theta[pairs[j]] - theta_center[pairs[j]]);
        return value;
    };
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        const auto& layout = prepared.shape[phase];
        std::vector<AD> x(components, AD(0.0));
        for (std::size_t local = 0; local < spec.support.size(); ++local)
            x[spec.support[local]] = CppAD::exp(u[layout.x[local]]);
        const AD eta = CppAD::exp(u[layout.volume]);
        const AD volume = (tape_parameters.ideal_gas ? AD(1.0) : packing_fraction_from_geometry(AD(1.0), x, tape_parameters, diameter)) / eta;
        const AD amount = layout.amount.has_value() ? u[*layout.amount] : AD(spec.amount.fraction) * feed_total;
        output.push_back(T);
        output.push_back(volume);
        output.insert(output.end(), x.begin(), x.end());
        output.insert(output.end(), theta.begin(), theta.end());
        output.push_back(amount);
        output.push_back(tape_parameters.ideal_gas ? AD(0.0) : eta);
        if (problem.volume_m3) for (std::size_t r = 0; r < compiled.rows.size(); ++r)
            if (compiled.rows[r].kind == RowKind::Volume)
                row_values[r] += amount * volume / u[u.size() - 1 - problem.internal_energy_j.has_value()];
        if (compiled.has_caloric_condition()) {
            const bool h = problem.enthalpy_j.has_value(), internal = problem.internal_energy_j.has_value();
            const AD denominator = feed_total * GAS_CONSTANT_J_PER_MOL_K * ((h || internal) ? T : AD(1.0));
            std::vector<AD> d_theta(theta.size());
            for (std::size_t a = 0; a < theta.size(); ++a) d_theta[a] = theta[a] - theta_center[a];
            const auto ideal = prepared.calorics ? anchored_ideal(*prepared.calorics, spec.support, x, volume, T, delta_T, delta_P, d_theta)
                                                 : caloric_ideal(compiled, phase, T, volume, x, prepared.temperature);
            row_values.back() += amount * (internal ? ideal.internal_energy : h ? ideal.enthalpy : ideal.entropy) / denominator;
            const AD weight = amount * GAS_CONSTANT_J_PER_MOL_K / denominator;
            output.push_back((h || internal) ? AD(0.0) : -weight);
            output.push_back((h || internal) ? -weight * T * T : -weight * T);
            output.push_back(h ? -weight * T * volume : AD(0.0));
            if (entropy) {
                const double rt = GAS_CONSTANT_J_PER_MOL_K * prepared.temperature;
                AD mu = 0.0;
                for (auto i : spec.support) mu += chemical[i] * x[i];
                row_values.back() -= amount * ideal.internal_energy / denominator;
                row_values.back() += amount * (-ideal.entropy / GAS_CONSTANT_J_PER_MOL_K
                    + ideal.internal_energy / rt + prepared.pressure * volume / rt);
                material_potential -= amount * mu;
                output[output.size()-3] = amount;
                output[output.size()-2] = amount * T * (1.0 - T / prepared.temperature);
            }
        }
        if (problem.critical || problem.spinodal) {
            // EqID critical_equations at n = x: D = diag(sqrt(z)) of the moving feed, w = D u, ideal
            // Q_ii = 1/n_i and F_iii = -1/n_i^2; row_series adds the residual F_nn w and F_nnn[w,w,w].
            const std::size_t mode = compiled.mode_row();
            row_values[mode + components] = AD(-1.0);
            for (std::size_t i = 0; i < components; ++i) {
                const AD u_i = u[compiled.mode_offset() + i], d = CppAD::sqrt(u[n + 2 + i] / feed_total), w = d * u_i;
                output.push_back(d);
                output.push_back(w);
                row_values[mode + i] = d * w / x[i];
                row_values[mode + components] += u_i * u_i;
                if (problem.critical) row_values[mode + components + 1] -= w * w * w / (x[i] * x[i]);
            }
        }
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            const auto& item = compiled.rows[row];
            if (condition_row(item.kind) || item.kind == RowKind::SolidAffinity) continue;
            if (item.kind == RowKind::Moiety) {
                for (std::size_t species = 0; species < components; ++species)
                    row_values[row] += amount * compiled.B[item.index][species] * x[species];
                continue;
            }
            if (item.phase != phase) continue;
            if (item.kind == RowKind::Chemical) {
                const std::size_t species = item.index;
                AD value = add_standard(CppAD::log(x[species]) - CppAD::log(volume), species);
                for (std::size_t balance = 0; balance < compiled.rank; ++balance)
                    value -= u[compiled.primal_size + balance] * compiled.B[balance][species];
                if (const auto& charge = compiled.charge_row[phase])
                    value -= AD(charges[species]) * u[compiled.primal_size + compiled.rank + *charge];
                row_values[row] = value;
            } else if (item.kind == RowKind::Pressure) {
                row_values[row] = (tape_parameters.ideal_gas ? std::accumulate(x.begin(), x.end(), AD(0.0)) : AD(1.0))
                    - P * volume / (GAS_CONSTANT_J_PER_MOL_K * T);
            } else if (item.kind == RowKind::Electro) {
                for (std::size_t species = 0; species < components; ++species)
                    row_values[row] += AD(charges[species]) * x[species];
            } else {
                for (const std::size_t species : spec.support) row_values[row] += x[species];
                row_values[row] -= AD(1.0);
            }
        }
    }
    // EqIDs solid_standard_gibbs and solid_saturation_index: the affinity subtracts the second-order Taylor term of
    // ln K^rho in (T, theta), K0 + (L' + C_T) dT + C_a dtheta_a + (L'' + C_TT) dT^2/2 + C_Ta dT dtheta_a + C_ab dtheta_a dtheta_b/2,
    // whose columns are those standard.solid carries; V_s (P - P0)/RT is exact.
    std::vector<AD> solid_output;
    for (std::size_t s = 0; s < compiled.solids.size(); ++s) {
        const auto& solid = compiled.solids[s];
        const auto& spec = problem.solids[s];
        AD affinity = -(solid.log_k + solid.log_k_slope.value_or(0.0) * delta_T + spec.molar_volume_m3_per_mol
            * (P - *spec.reference->reference_pressure_pa) / (GAS_CONSTANT_J_PER_MOL_K * T));
        if (!standard.solid.empty()) {  // solid_columns refused every moved column it lacks; an unmoved one is zero
            const auto& C = standard.solid[s];
            const auto column = [&C](std::size_t a, auto... p) {
                const ReferenceColumn c = C.reference.column(a, 0, p...);
                return c.state == ReferenceColumnState::Available ? c.value : 0.0;
            };
            affinity -= 0.5 * (C.curvature.value_or(0.0) + column(2, std::nullopt)) * delta_T * delta_T;
            for (std::size_t a = 0; a < theta.size(); ++a) {
                const AD d = theta[a] - theta_center[a];
                affinity -= (column(0, a) + column(1, a) * delta_T) * d;
                for (std::size_t b = 0; b <= a; ++b) affinity -= (b < a ? 1.0 : 0.5) * column(0, b, a) * (theta[b] - theta_center[b]) * d;
            }
        }
        for (std::size_t i = 0; i < components; ++i)
            if (spec.stoichiometry[i] != 0.0) affinity -= spec.stoichiometry[i] * add_standard(AD(0.0), i);
        for (std::size_t j = 0; j < compiled.rank; ++j) affinity += u[compiled.primal_size + j] * solid.balance[j];
        solid_output.push_back(solid.column ? u[*solid.column] : AD(0.0));
        solid_output.push_back(affinity);
    }
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        const auto& item = compiled.rows[row];
        if (item.kind == RowKind::SolidAffinity)
            row_values[row] = compiled.solids[item.index].present ? -solid_output[2 * item.index + 1] : solid_output[2 * item.index];
        if (item.kind != RowKind::Moiety) continue;
        for (std::size_t species = 0; species < components; ++species)
            row_values[row] -= compiled.B[item.index][species] * feed[species];
        for (std::size_t s = 0; s < compiled.solids.size(); ++s)
            row_values[row] += compiled.solids[s].balance[item.index] * solid_output[2 * s];
        row_values[row] /= problem.loading ? feed_total : AD(compiled.feed_total);
    }
    if (problem.volume_m3) for (std::size_t r = 0; r < compiled.rows.size(); ++r)
        if (compiled.rows[r].kind == RowKind::Volume) row_values[r] -= 1.0;
    if (problem.loading) row_values.back() = (P - problem.loading->pressure) / problem.loading->pressure;
    if (compiled.has_caloric_condition()) {
        const bool h = problem.enthalpy_j.has_value(), internal = problem.internal_energy_j.has_value();
        if (!entropy) row_values.back() -= u[internal ? u.size()-1 : n + compiled.control_width() + (h ? 0 : 1)]
            / (feed_total * GAS_CONSTANT_J_PER_MOL_K * ((h || internal) ? T : AD(1.0)));
    }
    for (std::size_t row = 0; row < compiled.rows.size(); ++row)
        if (compiled.rows[row].kind == RowKind::AzeotropeCondition) {  // EqID azeotrope_condition: ln K below at 0
            const std::size_t pin = problem.phases[0].is_pinned() ? 0 : 1, local = compiled.rows[row].index;
            row_values[row] = u[prepared.shape[pin].x[local]] - u[prepared.shape[1 - pin].x[local]];
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
    output.insert(output.end(), solid_output.begin(), solid_output.end());
    if (plane) for (std::size_t i = 0; i < components; ++i) {
        AD value = add_standard(AD(0.0), i);
        for (std::size_t k = 0; k < compiled.rank; ++k) value -= compiled.B[k][i] * u[compiled.primal_size + k];
        output.push_back(value);
    }
    if (entropy) output.push_back(material_potential);
    return CppAD::ADFun<double>(u, output);
}

std::vector<double> caloric_hessian(const CompiledSystem& compiled, const PhasePoint& prepared, double weight, bool entropy, std::vector<double>* magnitude) {
    const std::size_t n = prepared.point.size();
    std::vector<double> result(n * n, 0.0);
    if (weight == 0.0) return result;
    const std::size_t width = compiled.control_width();
    const std::size_t stride = width + 5;
    auto geometry = record_geometry(compiled, prepared, prepared.standard, entropy);
    using AD = CppAD::AD<double>;
    std::vector<AD> q(prepared.point.begin(), prepared.point.end());
    CppAD::Independent(q);
    const auto center = controls(compiled, prepared);
    std::vector<AD> input = q;
    input.insert(input.end(), center.begin() + static_cast<std::ptrdiff_t>(n), center.end());
    const std::vector<double> theta(center.begin() + static_cast<std::ptrdiff_t>(n + 2 + compiled.n_feed.size()),
                                    center.begin() + static_cast<std::ptrdiff_t>(n + width));
    auto active_geometry = geometry.base2ad();
    const auto mapped = active_geometry.Forward(0, input);
    AD value = entropy ? mapped[mapped.size()-2] + mapped.back() : mapped.back();
    for (std::size_t p = 0; p < prepared.shape.size(); ++p) {
        const auto& d = prepared.derivative[p];
        std::vector<double> phase_center = {prepared.temperature, prepared.volume[p]};
        phase_center.insert(phase_center.end(), prepared.composition[p].begin(), prepared.composition[p].end());
        phase_center.insert(phase_center.end(), theta.begin(), theta.end());
        std::vector<AD> delta(width);
        for (std::size_t i = 0; i < width; ++i) delta[i] = mapped[p * stride + i] - phase_center[i];
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
    CppAD::ADFun<double> pullback(q, std::vector<AD>{value, entropy ? mapped.back() : AD(0.0)});
    result = pullback.Hessian(prepared.point, 0);
    if (entropy) {
        const auto material = pullback.Hessian(prepared.point, 1);
        if (magnitude) { magnitude->resize(result.size());
            for (std::size_t i=0;i<result.size();++i) (*magnitude)[i]=std::abs(result[i]-material[i])+std::abs(material[i]); }
        auto gradient = pullback.Jacobian(prepared.point);
        for (std::size_t p = 0; p < prepared.shape.size(); ++p) {
            const auto& layout = prepared.shape[p];
            double normal = 0.0, total = 0.0;
            for (std::size_t i = 0; i < layout.x.size(); ++i) {
                normal -= gradient[layout.x[i]];
                total += prepared.composition[p][compiled.problem.phases[p].support[i]];
            }
            for (std::size_t i = 0; i < layout.x.size(); ++i) {
                const auto j = layout.x[i];
                const double term = normal * prepared.composition[p][compiled.problem.phases[p].support[i]] / total;
                gradient[j] += term; result[j*n+j] += term;
                if (magnitude) (*magnitude)[j*n+j] += std::abs(term);
            }
        }
        double defect = 0.0;
        for (std::size_t j = 0; j + 1 < compiled.primal_size; ++j)
            defect = std::max(defect, std::abs(gradient[j]) * (j == compiled.primal_size-2 ? prepared.temperature : 1.0) / compiled.feed_total);
        if (defect > compiled.problem.solver_options.value(policy::row))
            throw std::domain_error("entropy Lagrangian stationarity defect: " + std::to_string(defect));
    }
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
    const std::size_t n = prepared.point.size(), first = compiled.mode_row();
    std::vector<double> result(n * n, 0.0), diagonal(compiled.primal_size, 0.0);
    if (std::all_of(multipliers.begin() + static_cast<std::ptrdiff_t>(first), multipliers.end(),
                    [](double m) { return m == 0.0; })) return result;
    std::vector<std::vector<double>> state(3, std::vector<double>(n, 0.0));
    state[0] = prepared.point;
    const std::vector<double> control(compiled.control_width(), 0.0);
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

}  // namespace

// Reads jets.temperature_order, pressure_order, parameters, reference_order and pairs.
std::vector<ReactionSpec> central_reference_reactions(
    const HelmholtzModel &model, const EquilibriumProblem &request, double temperature, double pressure,
    std::vector<NeutralReferenceReactionDiagnostic> &diagnostics, const ReferenceJets& jets,
    const std::vector<NeutralReferenceReactionDiagnostic> *fixed_pressure, bool consumed_center) {
    if (!request.has_reaction_reference()) return request.reactions;
    if (!request.source_reference.has_value()) {
        throw ReferenceUnavailable("source reference declaration is missing");
    }
    if (!request.standard_gibbs_j_per_mol.empty()) {
        throw ReferenceUnavailable("supplied standard Gibbs conflicts with source reference");
    }
    const auto unavailable = [&](std::size_t row, const std::string& tail) {
        return ReferenceUnavailable("reaction " + reaction_name(request.reactions[row], row) + tail);
    };
    std::vector<ReactionSpec> result = request.reactions;
    std::vector<NeutralReferenceReactionDiagnostic> reused;
    if (fixed_pressure != nullptr && jets.temperature_order == 0 && jets.parameters.empty()) {
        // Fixed-pressure references depend only on the compiled model, temperature,
        // declaration and stoichiometry; rebuild only system-pressure traces.
        for (const auto& diagnostic : *fixed_pressure)
            if (consumed_center || request.reactions[diagnostic.reaction_index].reference->reference_pressure_pa)
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
        if (consumed_center && std::none_of(reused.begin(), reused.end(),
                [row](const auto& item) { return item.reaction_index == row; }))
            throw unavailable(row, " consumed reference unavailable");
        if (std::any_of(reused.begin(), reused.end(),
                [row](const auto& item) { return item.reaction_index == row; })) continue;
        if (reaction.affinity_j_per_mol.has_value()
            || (!reaction.correlation.has_value() && !reaction.equilibrium_constant.has_value()))
            throw unavailable(row, " needs a scalar K or correlation source");
        reference_rows.push_back(row);
        stoichiometries.push_back(reaction.stoichiometry);
        references.push_back(*reaction.reference);
    }
    const auto solvent_anchor = fixed_pressure ? reference_solvent_anchor(*fixed_pressure) : std::nullopt;
    const auto reference_results = consumed_center ? std::vector<NeutralReferenceReactionDiagnostic>{} : reaction_references(
        model, temperature, pressure, *request.source_reference, stoichiometries, references, jets.pressure_order,
        solvent_anchor, jets.temperature_order, jets.parameters, jets.reference_order, nullptr, jets.pairs, jets.temperature_order > 1);
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
        if (reference.status != NeutralReferenceStatus::Available || !reference.terminal_contraction.has_value())
            throw unavailable(row, reference.message.empty() ? "" : " " + reference.message);
        SolidLogK source;
        try {
            source = source_log_k(model, request.source_reference, reaction, temperature);
        } catch (const std::exception& error) {
            throw unavailable(row, std::string(" ") + error.what());
        }
        const double ln_k_density = source.ln_k + source.scale + *reference.terminal_contraction;
        const double target = -GAS_CONSTANT_J_PER_MOL_K * temperature * ln_k_density;
        if (!std::isfinite(target)) throw unavailable(row, " has a nonfinite density-basis target");
        result[row].equilibrium_constant.reset();
        result[row].correlation.reset();
        result[row].affinity_j_per_mol = target;
        result[row].reference.reset();
    }
    return result;
}

void admit_solid_declarations(const HelmholtzModel& model, const std::optional<NeutralReferenceDeclaration>& declaration,
                              const std::vector<SolidSpec>& solids) {
    const auto& z = model.parameters().charge_number;
    for (const auto& solid : solids) {
        const auto refuse = [&](const std::string& cause) { return std::invalid_argument("solid " + solid.name + ": " + cause); };
        const auto& v = solid.stoichiometry;
        if (v.size() != z.size() || std::all_of(v.begin(), v.end(), [](double x) { return x == 0.0; })
            || std::abs(std::inner_product(v.begin(), v.end(), z.begin(), 0.0)) > std::sqrt(std::numeric_limits<double>::epsilon()))
            throw refuse("stoichiometry must be nonzero, charge neutral and of component width");
        if (solid.equilibrium_constant.has_value() == solid.correlation.has_value()
            || (solid.equilibrium_constant && !(std::isfinite(*solid.equilibrium_constant) && *solid.equilibrium_constant > 0.0)))
            throw refuse("needs exactly one of a positive equilibrium constant or a correlation");
        if (!solid.reference || !solid.reference->reference_pressure_pa || !declaration)
            throw refuse("needs a neutral reference and reference_pressure_pa, the source standard pressure");
        if (!(std::isfinite(solid.molar_volume_m3_per_mol) && solid.molar_volume_m3_per_mol > 0.0))
            throw refuse("molar volume must be positive and finite");
    }
}

std::vector<SolidLogK> solid_log_k(const HelmholtzModel& model, double temperature, double pressure,
    const std::optional<NeutralReferenceDeclaration>& declaration, const std::vector<SolidSpec>& solids,
    const ReferenceJets& jets) {
    admit_solid_declarations(model, declaration, solids);
    std::vector<std::vector<double>> nu;
    std::vector<ReactionReference> references;
    for (const auto& solid : solids) nu.push_back(solid.stoichiometry), references.push_back(*solid.reference);
    const auto diagnostics = reaction_references(model, temperature, pressure, *declaration, nu, references, 0, std::nullopt,
        jets.temperature_order, jets.parameters, 0, nullptr, jets.pairs, jets.temperature_order > 1);
    std::vector<SolidLogK> result;
    for (std::size_t s = 0; s < solids.size(); ++s) {
        const auto& solid = solids[s];
        if (diagnostics[s].status != NeutralReferenceStatus::Available || !diagnostics[s].terminal_contraction)
            throw ReferenceUnavailable("solid " + solid.name + " " + diagnostics[s].message);
        result.push_back(source_log_k(model, declaration, solid_law(solid), temperature));
        result.back().reference = diagnostics[s];
    }
    return result;
}

bool CompiledSystem::anchored_calorics() const {
    if (!has_caloric_condition() || ideal || !model->thermochemistry()) return false;
    if (problem.has_reaction_reference()) return true;
    for (const auto& phase : problem.phases)
        for (const auto i : phase.support)
            if (model->thermochemistry()->component_intervals[i].empty()) return true;
    return false;
}

// EqID: reference_species_calorics
// Supplied ideal-gas records fix h and cp of their species; record-less species solve
// N_D x_D = Delta_r x - N_S x_S on the density basis for x = h, cp (min-norm), where only
// the ionic charge gauge may remain free: it cancels in every electroneutral phase.
ReferenceCalorics reference_calorics(const CompiledSystem &compiled, double temperature, double pressure,
                                     const ReferenceJets &jets) {
    if (compiled.ideal) throw std::invalid_argument("ideal reference calorics unavailable");
    const std::size_t derivative_order = std::max(jets.temperature_order, std::size_t(1));
    const bool pressure_slope = jets.pressure_order, parameter_slopes = !jets.parameters.empty(), numerical_errors = jets.errors;
    const auto &pair_coordinates = jets.pairs;
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
    if (model.parameters().ideal_gas && derivative_order > 3)
        return refuse(SolvedStateActionStatus::EnergyOrderUnavailable, "ideal-gas fourth-order calorics are deferred");
    for (std::size_t row = 0; row < problem.reactions.size(); ++row)
        if (problem.reactions[row].uses_species_thermochemistry())
            return refuse(SolvedStateActionStatus::ReferenceUnavailable, "species-derived calorics deferred for reaction " + reaction_name(problem.reactions[row], row));
    if (!record)
        return refuse(SolvedStateActionStatus::MissingThermochemistry,
                      "reference calorics need record anchors");
    std::vector<std::size_t> derived;
    bool supplied_ion = false;
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
        supplied_ion |= model.parameters().charge_number[i] != 0.0;
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
            const auto source = source_log_k(model, problem.source_reference, reaction, temperature);
            AD L = reaction.correlation ? reaction_log_value(*reaction.correlation, T[0]) : AD(source.ln_k);
            if (reaction.reference) L += source.scale;
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
    // Row sums of |nu_D|: nonzero on rows with record-less species; their maximum is the operator inf-norm of nu_D.
    const Eigen::VectorXd derived_weight = matrix.cwiseAbs().rowwise().sum();
    const double matrix_norm = derived_weight.lpNorm<Eigen::Infinity>();
    const auto nullity =
        derived.empty() ? Eigen::Index(0) : Eigen::Index(derived.size()) - decomposition.rank();
    if (nullity > 1 ||
        (nullity == 1 &&
         (supplied_ion || charge.norm() == 0.0 ||
          (matrix * charge).norm() > std::sqrt(std::numeric_limits<double>::epsilon()) * charge.norm())))
        return refuse(SolvedStateActionStatus::ReferenceUnavailable,
                      "undetermined neutral calorics: record-less freedom exceeds full charge gauge; rank " + std::to_string(decomposition.rank())
                      + "/" + std::to_string(derived.size()) + "; threshold " + core::optional_number(decomposition.threshold()));
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
    auto diagnostics =
        references.empty()
            ? std::vector<NeutralReferenceReactionDiagnostic>{}
            : reaction_references(model, temperature, pressure, *problem.source_reference, nu, references,
                                  jets.pressure_order, reference_solvent_anchor(compiled.reference_diagnostics),
                                  std::min(derivative_order, std::size_t(2)),
                                  jets.parameters, derivative_order, nullptr, pair_coordinates);
    if (numerical_errors && !derived.empty()) for (std::size_t r = 0; r < diagnostics.size(); ++r) {
        reference_errors(model, temperature, problem.reactions[reference_rows[r]].stoichiometry, diagnostics[r]);
        if (!diagnostics[r].numerical_error_cause.empty()) result.numerical_error_cause = diagnostics[r].numerical_error_cause;
    }
    // Parameter slots by count c: none, one coordinate, a requested coordinate pair.
    std::array<std::vector<std::array<std::optional<std::size_t>, 2>>, 3> slots;
    slots[0].push_back({});
    for (std::size_t p = 0; p < model.active_parameters().size(); ++p)
        slots[1].push_back({p, std::nullopt});
    for (std::size_t i = 0; i < pair_coordinates.size(); ++i)
        for (std::size_t j = i; j < pair_coordinates.size(); ++j)
            slots[2].push_back({pair_coordinates[i], pair_coordinates[j]});
    for (std::size_t c = 0; c <= (pair_coordinates.empty() ? std::size_t(parameter_slopes) : 2); ++c)
        for (std::size_t a = 0; a + c <= derivative_order; ++a)
            for (std::size_t b = 0; a + b + c <= derivative_order; ++b) {
                if (b && !pressure_slope)
                    continue;
                for (const auto &slot : slots[c]) {
                    ReferenceCalorics::Derivative entry{a, b, slot[0], slot[1], {}};
                    auto &column = entry.column;
                    column.values.assign(components, 0.0);
                    Eigen::VectorXd right(rows), scale(rows);
                    for (std::size_t i = 0; i < components; ++i)
                        if (!b && !c)
                            column.values[i] = thermal[a][i];
                    for (std::size_t r = 0; r < rows; ++r)
                        right(r) = !b && !c ? thermal[a][components + r] : 0.0;
                    for (std::size_t r = 0; r < diagnostics.size(); ++r) {
                        const auto upper = diagnostics[r].column(a, b, entry.parameter, entry.second_parameter);
                        const auto lower = a ? diagnostics[r].column(a - 1, b, entry.parameter, entry.second_parameter)
                                             : ReferenceColumn{ReferenceColumnState::Available};
                        const bool has_upper = upper.state == ReferenceColumnState::Available;
                        if (!has_upper || lower.state != ReferenceColumnState::Available) {
                            column.status = SolvedStateActionStatus::ReferenceUnavailable;
                            column.message = ReferenceUnavailable::label(
                                "reaction " + reaction_name(problem.reactions[reference_rows[r]], reference_rows[r]) +
                                " C_" + std::to_string(has_upper ? a - 1 : a) + std::to_string(b) +
                                (c ? " parameter " + std::to_string(*entry.parameter) : "") + ": " +
                                *(has_upper ? lower : upper).cause);
                        } else
                            right(reference_rows[r]) -=
                                GAS_CONSTANT_J_PER_MOL_K * (temperature * upper.value + double(a) * lower.value);
                    }
                    Eigen::VectorXd source_error = Eigen::VectorXd::Zero(numerical_errors ? rows : 0);
                    if (numerical_errors && !b && !c && a <= 1) for (std::size_t r = 0; r < diagnostics.size(); ++r) {
                        const auto& e = diagnostics[r].numerical_error;
                        source_error(reference_rows[r]) = GAS_CONSTANT_J_PER_MOL_K * (temperature * e[a] + (a ? e[0] : 0.0));
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
                        // Completion roundoff follows ||nu_D|| ||x||, so only rows with record-less species carry it.
                        const double bound = matrix_norm * solution.lpNorm<Eigen::Infinity>();
                        for (Eigen::Index r = 0; r < residual.size(); ++r) {
                            const double criterion = std::sqrt(std::numeric_limits<double>::epsilon())
                                * std::max(1.0, scale(r) + (derived_weight(r) > 0.0 ? bound : 0.0));
                            if (!(std::abs(residual(r)) <= criterion)) {
                                column.status = SolvedStateActionStatus::ReferenceUnavailable;
                                column.message = "supplied records inconsistent with reaction " +
                                                 reaction_name(problem.reactions[r], r)
                                    + "; residual " + core::optional_number(std::abs(residual(r)))
                                    + "; criterion " + core::optional_number(criterion);
                                if (!derived.empty()) column.message += "; reference_completion rank " + std::to_string(decomposition.rank())
                                    + "/" + std::to_string(derived.size()) + "; threshold " + core::optional_number(decomposition.threshold());
                            }
                        }
                        for (std::size_t j = 0; j < derived.size(); ++j)
                            column.values[derived[j]] = solution(j);
                        if (numerical_errors && !b && !c && a <= 1) {
                            Eigen::MatrixXd inverse(derived.size(), rows);
                            for (std::size_t r = 0; r < rows; ++r) for (std::size_t j = 0; j < derived.size(); ++j)
                                inverse(j, r) = result.reaction_map[r][derived[j]];
                            source_error += std::numeric_limits<double>::epsilon() * (scale + right.cwiseAbs());
                            const auto mapped = map_error(matrix, inverse, right, solution, source_error);
                            column.numerical_error.assign(components, 0.0);
                            for (std::size_t j = 0; j < derived.size(); ++j) column.numerical_error[derived[j]] = mapped(j);
                        }
                    }
                    result.gibbs.push_back(std::move(entry));
                }
            }
    return result;
}

std::vector<std::size_t> CompiledSystem::phase_offsets() const {
    if (ideal) return {0};
    std::vector<std::size_t> result;
    for (const auto& item : shape) result.push_back(item.x.front());
    return result;
}

std::vector<double> CompiledSystem::entropy_hessian_at(const PhasePoint& prepared, std::vector<double>* magnitude) const {
    return caloric_hessian(*this, prepared, 1.0, true, magnitude);
}

// The layout owners of generic_point name each column: phase compositions, free amounts and volume slots (log packing,
// or log molar density for an ideal-gas model), free solids, critical mode, loading, free T/P, then multipliers.
std::vector<std::string> CompiledSystem::coordinate_labels() const {
    core::Scope scope("coordinate_labels");
    const auto id = [](std::size_t component) { return "component " + std::to_string(component); };  // model order
    std::vector<std::string> labels(unknown_count());
    if (ideal) {
        for (std::size_t j = 0; j < primal_size; ++j) labels[j] = "ln(n[" + id(problem.phases[0].support[j]) + "]/N_feed)";
        return labels;
    }
    for (std::size_t p = 0; p < shape.size(); ++p) {
        const auto& name = problem.phases[p].name;
        for (std::size_t i = 0; i < shape[p].x.size(); ++i)
            labels[shape[p].x[i]] = "ln x[" + name + "," + id(problem.phases[p].support[i]) + "]";
        if (shape[p].amount) labels[*shape[p].amount] = "n[" + name + "] / mol";
        labels[shape[p].volume] = model->parameters().ideal_gas ? "ln rho[" + name + "] / (mol/m3)" : "ln eta[" + name + "]";
    }
    for (std::size_t s = 0; s < solids.size(); ++s)
        if (solids[s].column) labels[*solids[s].column] = "n_solid[" + problem.solids[s].name + "] / mol";
    for (std::size_t k = 0; k < mode_count; ++k) labels[mode_offset() + k] = "u[" + id(k) + "]";
    if (problem.loading) labels[loading_column()] = "loading";
    std::size_t offset = primal_size - free_intensive_count;
    if (problem.temperature.free) labels[offset++] = "T / K";
    if (problem.pressure.free) labels[offset++] = "P / Pa";
    for (std::size_t k = 0; k < rank; ++k) labels[primal_size + k] = "lambda[moiety " + std::to_string(k) + "]";
    for (std::size_t r = 0; r < charge_row_count; ++r) labels[primal_size + rank + r] = "phi[charge row " + std::to_string(r) + "]";
    return labels;
}

std::vector<core::DifferenceCheck> CompiledSystem::jacobian_difference(const std::vector<double>& point, double relative_step,
    double tau, double kappa) const {
    core::Scope scope("jacobian_difference");
    core::admit_difference_parameters(relative_step, tau, kappa);
    if (point.size() != unknown_count()) throw std::invalid_argument("jacobian difference needs every solved-chart coordinate");
    const auto exact = jacobian(point);
    const auto labels = coordinate_labels();
    const std::size_t n = point.size();
    std::vector<core::DifferenceCheck> out;
    for (std::size_t j = 0; j < n; ++j) {
        auto plus = point, minus = point;
        plus[j] += relative_step * std::max(1.0, std::abs(point[j]));
        minus[j] -= relative_step * std::max(1.0, std::abs(point[j]));
        const auto up = evaluate(plus), down = evaluate(minus);
        for (std::size_t i = 0; i < rows.size(); ++i)
            out.push_back(core::difference_check(std::string(row_kind_name(rows[i].kind)) + "[" + std::to_string(rows[i].phase) + ","
                + std::to_string(rows[i].index) + "] / " + labels[j], exact[i * n + j], up.values[i], down.values[i], up.floor, down.floor,
                0.5 * (plus[j] - minus[j]), tau, kappa));
    }
    return out;
}

PhasePoint CompiledSystem::prepare(const std::vector<double>& point, bool state_only,
    std::optional<core::Convergence>* first, std::optional<std::size_t> chemical_phase, const std::vector<bool>& site_phases, bool parameter_actions) const {
    return prepare_values(point, state_only, first, chemical_phase, site_phases, parameter_actions);
}

template <typename Base>
PhasePointValues<Base> CompiledSystem::prepare_values(const std::vector<Base>& point, bool state_only,
    std::optional<core::Convergence>* first_numerical_failure, std::optional<std::size_t> chemical_phase,
    const std::vector<bool>& site_phases, bool parameter_actions, bool gradient_only, bool rank_only) const {
    core::Scope scope("equilibrium_prepare");
    if (point.size() < primal_size) {
        throw std::invalid_argument("equilibrium point is too short");
    }
    if constexpr (std::is_same_v<Base, double>) if (ideal) {
        if (point.size() != primal_size) throw std::invalid_argument("ideal amount coordinate count");
        PhasePointValues<Base> result; result.point = point;
        result.temperature = problem.temperature.value; result.pressure = problem.pressure.value;
        std::vector<Base> amounts(n_feed.size(), 0.0);
        for (std::size_t j = 0; j < point.size(); ++j) {
            const Base amount = feed_total * wide_math::exp(point[j]);
            if (!wide_math::isfinite(amount) || amount <= 0.0) throw std::invalid_argument("invalid ideal amount");
            amounts[problem.phases[0].support[j]] = amount;
        }
        const Base total = sum(amounts);
        if (!wide_math::isfinite(total) || total <= 0.0) throw std::invalid_argument("invalid ideal total amount");
        for (Base& amount : amounts) amount /= total;
        result.composition.push_back(std::move(amounts)); result.amount.push_back(total);
        return result;
    }
    const std::size_t components = model->parameters().component_count();
    PhasePointValues<Base> result;
    result.point = point;
    result.shape = shape;
    std::size_t offset = primal_size - free_intensive_count;
    result.temperature = problem.temperature.free ? point[offset++] : problem.temperature.value;
    result.pressure = problem.pressure.free ? point[offset++] : problem.pressure.value;
    if constexpr (std::is_same_v<Base, double>) admit_calorics(*this, result.temperature);
    if constexpr (std::is_same_v<Base, double>)
        if (!state_only && !chemical_phase && !parameter_actions && model->active_parameters().empty()) result.state_model = model->identity();
    result.composition.resize(result.shape.size()); result.packing.resize(result.shape.size());
    result.eta.resize(result.shape.size()); result.volume.resize(result.shape.size());
    result.amount.resize(result.shape.size()); result.derivative.resize(result.shape.size());
    for (std::size_t phase = 0; phase < result.shape.size(); ++phase) {
        if (chemical_phase && phase != *chemical_phase) continue;
        core::AssociationObservations phase_observations(phase);
        const auto& layout = result.shape[phase];
        std::vector<Base> values(components, 0.0);
        Base total = 0.0;
        for (std::size_t i = 0; i < layout.x.size(); ++i) {
            const Base value = wide_math::exp(point[layout.x[i]]);
            if (!wide_math::isfinite(value) || value <= 0.0) throw std::invalid_argument("invalid composition in phase " + problem.phases[phase].name);
            values[problem.phases[phase].support[i]] = value;
            total += value;
        }
        if (!(total > 0.0)) throw std::invalid_argument("composition in phase " + problem.phases[phase].name + " has zero total");
        result.packing[phase] = detail::packing_terms(model->parameters(), result.temperature, values);
        result.eta[phase] = wide_math::exp(point[layout.volume]);
        if (!model->parameters().ideal_gas && log_packing_unbounded() && !admitted_packing_fraction(double(result.eta[phase])))
            throw std::invalid_argument("phase " + problem.phases[phase].name + " packing fraction outside the fluid domain");
        result.volume[phase] = result.packing[phase].coefficient / result.eta[phase];
        result.amount[phase] = layout.amount.has_value()
            ? point[*layout.amount] : phase_amount(problem.phases[phase]);
        const auto order = !state_only ? (parameter_actions ? DerivativeOrder::Hessian : DerivativeOrder::FixedParameterHessian)
            : problem.temperature.free ? (gradient_only ? DerivativeOrder::PhaseGradient : DerivativeOrder::PhaseHessian)
                                       : (gradient_only ? DerivativeOrder::IsothermalPhaseGradient : DerivativeOrder::IsothermalPhaseHessian);
        if constexpr (std::is_same_v<Base, double>) {
        if (state_only && one_component_two_density_eligible(*this)) {
            result.pure.push_back(model->pure_isothermal(result.temperature, result.volume[phase], total));
            result.derivative[phase] = result.pure.back().derivative;
        } else result.derivative[phase] = model->evaluate(result.temperature, result.volume[phase], values, order,
            phase_calculation == PhaseCalculation::Declared, !site_phases.empty() && site_phases[phase], rank_only);
        } else result.derivative[phase] = model->isothermal(result.temperature, result.volume[phase], values);
        if (first_numerical_failure) core::retain_first(*first_numerical_failure, result.derivative[phase].first_numerical_failure);
        scope.probe("gradient_floor", double(result.derivative[phase].gradient_floor),
                    static_cast<std::uint32_t>(phase));
        result.composition[phase] = std::move(values);
    }
    result.parameter_actions = parameter_actions;
    result.standard = fixed_standard ? *fixed_standard : standard_term(
        *this, double(result.temperature), double(result.pressure), 2, problem.pressure.free, problem.temperature.free);
    if constexpr (std::is_same_v<Base, double>) if (anchored_calorics())
        result.calorics = caloric_reference(*this, result.temperature, result.pressure);
    return result;
}

template PhasePointValues<double> CompiledSystem::prepare_values(const std::vector<double>&, bool, std::optional<core::Convergence>*, std::optional<std::size_t>, const std::vector<bool>&, bool, bool, bool) const;
template PhasePointValues<WiderScalar> CompiledSystem::prepare_values(const std::vector<WiderScalar>&, bool, std::optional<core::Convergence>*, std::optional<std::size_t>, const std::vector<bool>&, bool, bool, bool) const;

PointEvaluation CompiledSystem::evaluate(const std::vector<double>& point) const {
    return evaluate_at(prepare(point));
}

PointEvaluation CompiledSystem::evaluate_at(const PhasePoint& at, std::optional<std::size_t> chemical_phase) const { return evaluate_at<double>(at, chemical_phase); }

template <typename Base>
Base CompiledSystem::chemical_row(const PhasePointValues<Base>& prepared, std::size_t phase, std::size_t species, const Base& fraction) const {
    const auto& point = prepared.point;
    Base constraint = 0.0;
    for (std::size_t k = 0; k < rank; ++k)
        constraint += (primal_size + k < point.size() ? point[primal_size + k] : Base(0.0)) * B[k][species];
    const auto charge = charge_row[phase];
    const Base electro = model->parameters().charge_number[species]
        * (charge && primal_size + rank + *charge < point.size() ? point[primal_size + rank + *charge] : Base(0.0));
    return wide_math::log(std::max(fraction / prepared.volume[phase], std::numeric_limits<Base>::min()))
        + prepared.derivative[phase].gradient[2 + species] + prepared.standard.derivative[0][species] - constraint - electro;
}
template double CompiledSystem::chemical_row(const PhasePointValues<double>&, std::size_t, std::size_t, const double&) const;

template <typename Base>
PointEvaluationValues<Base> CompiledSystem::evaluate_at(const PhasePointValues<Base>& prepared, std::optional<std::size_t> chemical_phase) const {
    core::CalculationWork::equilibrium_entry();
    core::Scope scope("equilibrium_evaluate");
    if constexpr (std::is_same_v<Base, double>) if (ideal) return {ideal_rows(*this, prepared).Forward(0, prepared.point), 0.0, 0.0};
    const auto& point = prepared.point;
    const auto balance = [&] { if constexpr (std::is_same_v<Base, double>) if (problem.loading) return matvec(B, feed_at(point)); return this->balance; }();
    const double feed_total = [&] { if constexpr (std::is_same_v<Base, double>) return feed_total_at(prepared); else return this->feed_total; }();
    const std::size_t components = model->parameters().component_count();
    const auto& z = model->parameters().charge_number;
    const auto& x = prepared.composition;
    const auto& amount = prepared.amount;
    const auto& volume = prepared.volume;
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * prepared.temperature;
    const auto& standard = prepared.standard.derivative[0];
    PointEvaluationValues<Base> result;
    if (!chemical_phase) for (const auto& derivative : prepared.derivative) result.floor = std::max(result.floor, derivative.gradient_floor);
    std::vector<Base> series;  // every condition row from one zeroth coefficient
    for (const auto& row : rows) {
        if (chemical_phase && (row.kind != RowKind::Chemical || row.phase != *chemical_phase)) {
            result.values.push_back(0.0); continue;
        }
        Base value = 0.0;
        if (row.kind == RowKind::Chemical) {
            value = chemical_row(prepared, row.phase, row.index, x[row.phase][row.index]);
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = prepared.derivative[row.phase];
            value = volume[row.phase] * (
                (model->parameters().ideal_gas ? sum(x[row.phase]) : 1.0) / volume[row.phase] - derivative.gradient[1] - prepared.pressure / rt);
        } else if (row.kind == RowKind::Moiety) {
            for (std::size_t phase = 0; phase < shape.size(); ++phase) {
                for (std::size_t species = 0; species < components; ++species) {
                    value += amount[phase] * B[row.index][species] * x[phase][species];
                }
            }
            for (const auto& solid : solids)
                if (solid.column) value += solid.balance[row.index] * point[*solid.column];
            value = (value - balance[row.index]) / feed_total;
        } else if (row.kind == RowKind::PrescribedComposition) {
            value = point[shape[0].x[row.index]] - Base(std::log((*problem.prescribed_liquid_composition)[problem.phases[0].support[row.index]]));
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t species = 0; species < components; ++species) value += z[species] * x[row.phase][species];
        } else if (row.kind == RowKind::SolidAffinity) {
            if constexpr (std::is_same_v<Base, double>) { const auto& solid = solids[row.index];
            value = solid.present ? -solid_affinity(prepared, row.index) : point[*solid.column]; }
        } else if (row.kind == RowKind::Volume) {
            value = -1.0;
            for (std::size_t p = 0; p < amount.size(); ++p) value += amount[p] * volume[p] / *problem.volume_m3;
        } else if (condition_row(row.kind)) {
            if constexpr (std::is_same_v<Base, double>) { if (series.empty())
                series = row_series(prepared, {point}, std::vector<Base>(control_width()), 0).rows[0];
            value = series[result.values.size()]; }
        } else {
            for (const std::size_t species : problem.phases[row.phase].support) value += x[row.phase][species];
            value -= 1.0;
        }
        result.values.push_back(value);
    }
    if constexpr (std::is_same_v<Base, double>) scope.probe("row_values", result.values);
    if (problem.prescribed_liquid_composition) {
        result.objective = NAN;
        return result;
    }
    if (chemical_phase) return result;
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        const auto& derivative = prepared.derivative[phase];
        Base ideal = ideal_mechanical_helmholtz(x[phase], volume[phase], Base(1));
        result.objective += amount[phase] * (
            derivative.value + ideal + (problem.volume_m3 ? Base(0) : prepared.pressure * volume[phase] / rt)
            #ifdef EPCSAFT_SOFTWARE_WIDER
            + std::inner_product(standard.begin(), standard.end(), x[phase].begin(), 0.0,
                [](double sum, Base product) { return static_cast<double>(sum + product); }, std::multiplies<Base>{})
#else
            + std::inner_product(standard.begin(), standard.end(), x[phase].begin(), 0.0)
#endif
        );
    }
    if constexpr (std::is_same_v<Base, double>) for (std::size_t s = 0; s < solids.size(); ++s)
        if (solids[s].column) result.objective += point[*solids[s].column] * solid_gibbs(prepared, s);
    if constexpr (std::is_same_v<Base, double>) if (problem.internal_energy_j) {
        auto entropy = *this;
        entropy.problem.internal_energy_j.reset(); entropy.problem.entropy_j_per_k = 0.0;
        result.objective = -entropy.row_series(prepared, {point}, std::vector<double>(control_width()), 0).rows[0].back() * feed_total_at(prepared);
    }
    return result;
}

template PointEvaluationValues<WiderScalar> CompiledSystem::evaluate_at(const PhasePointValues<WiderScalar>&, std::optional<std::size_t>) const;

RowUncertainty CompiledSystem::row_uncertainty(const PhasePoint& at) const {
    const auto balance = problem.loading ? matvec(B, feed_at(at.point)) : this->balance;
    const double feed_total = feed_total_at(at);
    RowUncertainty out;
    out.terms.assign(rows.size(), 0.0); out.inner = out.terms;
    if (ideal) {
        const auto& support = problem.phases[0].support;
        const double shift = *std::max_element(at.point.begin(), at.point.end());
        double total = 0.0;
        for (double y : at.point) total += std::exp(y - shift);
        const double log_total = shift + std::log(total);
        for (std::size_t r = 0; r < problem.reactions.size(); ++r) {
            out.terms[r] = std::abs(ideal_log_k[r]);
            for (std::size_t j = 0; j < support.size(); ++j) {
                const double activity = ideal->activity_basis == "molarity"
                    ? ideal->molarity_log_activity(at.point[j] + std::log(feed_total), support[j]) : at.point[j] - log_total;
                out.terms[r] += std::abs(problem.reactions[r].stoichiometry[support[j]] * activity);
            }
        }
        for (std::size_t k = 0; k < rank; ++k) {
            auto& terms = out.terms[problem.reactions.size() + k];
            terms = std::abs(balance[k] / feed_total);
            for (std::size_t j = 0; j < support.size(); ++j) terms += std::abs(B[k][support[j]] * std::exp(at.point[j]));
        }
        return out;
    }
    const auto components = model->parameters().component_count(), width = 2 + components + model->active_parameters().size();
    const double rt = GAS_CONSTANT_J_PER_MOL_K * at.temperature;
    std::vector<std::vector<double>> gradient_error(at.shape.size(), std::vector<double>(width, 0.0));
    std::vector<double> energy_error(at.shape.size(), 0.0);
    for (std::size_t p = 0; p < at.shape.size(); ++p) {
        std::vector<double> point{at.temperature, at.volume[p]};
        point.insert(point.end(), at.composition[p].begin(), at.composition[p].end());
        for (std::size_t j = 0; j < point.size(); ++j) {
            if (j >= 2 && at.composition[p][j - 2] == 0.0) { gradient_error[p][j] = NAN; continue; }
            std::vector<double> unit(point.size(), 0.0); unit[j] = 1.0;
            const auto error = model->residual_coefficient_errors({point, unit}, at.derivative[p].site_fractions);
            energy_error[p] = error[0]; gradient_error[p][j] = error[1];
        }
    }
    StandardTerm::NumericalError standard_error;
    if (!reference_diagnostics.empty())
        standard_term(*this, at.temperature, at.pressure, 0, false, false, {}, &standard_error, {}, {}, &at.standard);
    if (!standard_error.cause.empty()) { out.cause = standard_error.cause; return out; }
    else if (!standard_error.absolute[0].empty()) for (auto& error : gradient_error)
        for (std::size_t j = 0; j < components; ++j) error[2 + j] += standard_error.absolute[0][j];
    out = row_uncertainty_at(at, gradient_error);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i]; const auto p = row.phase, species = row.index;
        double& terms = out.terms[i]; double& inner = out.inner[i];
        if ((row.kind == RowKind::CriticalNormalization || row.kind == RowKind::SpinodalNormalization)) {
            terms = 1;
            for (std::size_t j = 0; j < mode_count; ++j) terms += std::pow(at.point[mode_offset() + j], 2);
        } else if (row.kind == RowKind::AzeotropeCondition) {
            const auto local = problem.phases[0].support[species];
            terms = wide_math::abs(wide_math::log(at.composition[0][local])) + wide_math::abs(wide_math::log(at.composition[1][local]));
        } else if (row.kind == RowKind::SaturationPressure) {
            terms = (wide_math::abs(at.pressure) + problem.loading->pressure) / problem.loading->pressure;
        } else if (row.kind == RowKind::EnvelopeCondition) {
            const auto& e = *problem.envelope;
            terms = wide_math::abs(e.target);
            if (e.kind != EnvelopeCoordinate::LogK) terms += wide_math::abs(wide_math::log(e.kind == EnvelopeCoordinate::LogPressure ? at.pressure : at.temperature));
            else terms += wide_math::abs(wide_math::log(at.composition[0][e.component])) + wide_math::abs(wide_math::log(at.composition[1][e.component]));
        } else if (row.kind == RowKind::Enthalpy || row.kind == RowKind::Entropy || row.kind == RowKind::InternalEnergy) {
            const bool enthalpy = row.kind == RowKind::Enthalpy, internal = row.kind == RowKind::InternalEnergy;
            const double denominator = feed_total * GAS_CONSTANT_J_PER_MOL_K * ((enthalpy || internal) ? at.temperature : 1.0);
            terms = wide_math::abs(internal ? *problem.internal_energy_j : enthalpy ? *problem.enthalpy_j : *problem.entropy_j_per_k) / denominator;
            // Anchored g and g_T carry the reference and completion errors of their owner.
            ReferenceJets jets;
            jets.temperature_order = 1;
            jets.errors = true;
            const auto errors = at.calorics ? reference_calorics(*this, at.temperature, at.pressure, jets) : ReferenceCalorics{};
            if (!errors.numerical_error_cause.empty()) { out.cause = errors.numerical_error_cause; continue; }
            for (std::size_t a = 0; a < at.shape.size(); ++a) {
                const auto& support = problem.phases[a].support;
                const auto ideal = at.calorics ? anchored_ideal(*at.calorics, support, at.composition[a], at.volume[a], at.temperature, 0.0)
                                               : caloric_ideal(*this, a, at.temperature, at.volume[a], at.composition[a], at.temperature);
                const double weight = at.amount[a] / feed_total;
                for (std::size_t i = 0; at.calorics && i < support.size(); ++i) {
                    const auto& e0 = errors.gibbs[0].column.numerical_error, &e1 = errors.gibbs[1].column.numerical_error;
                    inner += at.amount[a] * at.composition[a][support[i]] * (enthalpy ? e0[support[i]] + at.temperature * e1[support[i]] : e1[support[i]]) / denominator;
                }
                terms += wide_math::abs(at.amount[a] * (internal ? ideal.internal_energy : enthalpy ? ideal.enthalpy : ideal.entropy)) / denominator;
                terms += weight * (at.temperature * wide_math::abs(at.derivative[a].gradient[0])
                    + (enthalpy ? wide_math::abs(at.volume[a] * at.derivative[a].gradient[1]) : internal ? 0.0 : wide_math::abs(at.derivative[a].value)));
                inner += weight * (at.temperature * gradient_error[a][0]
                    + (enthalpy ? wide_math::abs(at.volume[a]) * gradient_error[a][1] : internal ? 0.0 : energy_error[a]));
            }
        } else if ((row.kind == RowKind::CriticalMode || row.kind == RowKind::SpinodalMode) || row.kind == RowKind::CriticalCubic) {
            const auto& derivative = at.derivative[p];
            std::vector<double> point{at.temperature, at.volume[p]}, w(width, 0.0);
            point.insert(point.end(), at.composition[p].begin(), at.composition[p].end());
            const auto active = model->active_parameter_values(); point.insert(point.end(), active.begin(), active.end());
            for (std::size_t j = 0; j < components; ++j) w[2 + j] = wide_math::sqrt(n_feed[j] / feed_total) * at.point[mode_offset() + j];
            const auto coefficient_error = [&](const std::vector<double>& direction, std::size_t order) {
                std::vector<std::vector<double>> path(order + 1, std::vector<double>(width, 0.0));
                path[0] = point; path[1] = direction;
                return model->residual_coefficient_errors(path, derivative.site_fractions)[order];
            };
            if ((row.kind == RowKind::CriticalMode || row.kind == RowKind::SpinodalMode)) {
                const double d = wide_math::sqrt(n_feed[species] / feed_total);
                terms = wide_math::abs(d * w[2 + species] / at.composition[p][species]);
                for (std::size_t j = 0; j < components; ++j)
                    terms += wide_math::abs(d * derivative.hessian[(2 + species) * width + 2 + j] * w[2 + j]);
                // Polarization of E2 supplies this requested mixed contraction without a tensor.
                std::vector<double> unit(width, 0.0), sum = w; unit[2 + species] = 1.0; sum[2 + species] += 1.0;
                inner = wide_math::abs(d) * (coefficient_error(sum, 2) + coefficient_error(unit, 2) + coefficient_error(w, 2));
            } else {
                const auto nominal = model->residual_energy_action(at.temperature, at.volume[p], at.composition[p], {w, w, w});
                if (!nominal.action) { out.cause = "critical_contraction_unavailable"; continue; }
                terms = wide_math::abs(*nominal.action);
                for (std::size_t j = 0; j < components; ++j)
                    terms += wide_math::abs(std::pow(w[2 + j], 3) / std::pow(at.composition[p][j], 2));
                inner = 6 * coefficient_error(w, 3);
            }
        } else if (row.kind == RowKind::SolidAffinity) {
            const auto& solid = solids[species]; const auto& spec = problem.solids[species];
            if (!solid.present) { terms = wide_math::abs(at.point[*solid.column]); continue; }
            terms = wide_math::abs(solid.log_k) + wide_math::abs(spec.molar_volume_m3_per_mol * (at.pressure - *spec.reference->reference_pressure_pa) / rt);
            for (std::size_t j = 0; j < components; ++j) terms += wide_math::abs(spec.stoichiometry[j] * at.standard.derivative[0][j]);
            for (std::size_t j = 0; j < rank; ++j) terms += wide_math::abs(at.point[primal_size + j] * solid.balance[j]);
            auto reference = *solid.reference;
            reference_errors(*model, at.temperature, spec.stoichiometry, reference);
            if (!reference.numerical_error_cause.empty() || !standard_error.cause.empty()) {
                inner = NAN; out.cause = !reference.numerical_error_cause.empty() ? reference.numerical_error_cause : standard_error.cause;
            } else {
                const auto source = source_log_k(*model, problem.source_reference, solid_law(spec), at.temperature);
                inner = reference.numerical_error[0] + std::numeric_limits<double>::epsilon()
                    * (wide_math::abs(source.ln_k) + wide_math::abs(source.scale) + wide_math::abs(*reference.terminal_contraction));
                for (std::size_t j = 0; j < components; ++j) if (!standard_error.absolute[0].empty())
                    inner += wide_math::abs(spec.stoichiometry[j]) * standard_error.absolute[0][j];
            }

        }
    }
    return out;
}

template <typename Base>
RowUncertaintyValues<Base> CompiledSystem::row_uncertainty_at(const PhasePointValues<Base>& at,
    const std::vector<std::vector<Base>>& gradient_error) const {
    const auto balance = [&] { if constexpr (std::is_same_v<Base, double>) if (problem.loading) return matvec(B, feed_at(at.point)); return this->balance; }();
    const double feed_total = [&] { if constexpr (std::is_same_v<Base, double>) return feed_total_at(at); else return this->feed_total; }();
    RowUncertaintyValues<Base> out;
    out.terms.assign(rows.size(), 0); out.inner = out.terms;
    const auto components = model->parameters().component_count();
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * at.temperature;
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i]; const auto p = row.phase, species = row.index;
        Base& terms = out.terms[i]; Base& inner = out.inner[i];
        if (row.kind == RowKind::Chemical) {
            terms = wide_math::abs(wide_math::log(at.composition[p][species] / at.volume[p]))
                + wide_math::abs(at.derivative[p].gradient[2 + species]) + wide_math::abs(at.standard.derivative[0][species]);
            inner = gradient_error[p][2 + species];
            for (std::size_t k = 0; k < rank; ++k) terms += wide_math::abs(at.point[primal_size + k] * B[k][species]);
            if (charge_row[p]) terms += wide_math::abs(model->parameters().charge_number[species] * at.point[primal_size + rank + *charge_row[p]]);
        } else if (row.kind == RowKind::Pressure) {
            terms = (model->parameters().ideal_gas ? wide_math::abs(sum(at.composition[p])) : 1.0)
                + wide_math::abs(at.volume[p] * at.derivative[p].gradient[1]) + wide_math::abs(at.pressure * at.volume[p] / rt);
            inner = wide_math::abs(at.volume[p]) * gradient_error[p][1];
        } else if (row.kind == RowKind::Moiety) {
            terms = wide_math::abs(balance[species]);
            for (std::size_t a = 0; a < at.shape.size(); ++a)
                for (std::size_t j = 0; j < components; ++j) terms += wide_math::abs(at.amount[a] * B[species][j] * at.composition[a][j]);
            for (const auto& solid : solids) if (solid.column) terms += wide_math::abs(solid.balance[species] * at.point[*solid.column]);
            terms /= feed_total;
        } else if (row.kind == RowKind::PrescribedComposition) {
            terms = wide_math::abs(at.point[shape[0].x[row.index]])
                + std::abs(std::log((*problem.prescribed_liquid_composition)[problem.phases[0].support[row.index]]));
        } else if (row.kind == RowKind::Volume) {
            terms = 1.0;
            for (std::size_t a = 0; a < at.amount.size(); ++a)
                terms += wide_math::abs(at.amount[a] * at.volume[a] / *problem.volume_m3);
        } else if (row.kind == RowKind::Electro || row.kind == RowKind::Normalization) {
            terms = row.kind == RowKind::Normalization ? 1.0 : 0.0;
            for (std::size_t j = 0; j < components; ++j)
                terms += wide_math::abs(at.composition[p][j] * (row.kind == RowKind::Electro ? model->parameters().charge_number[j] : 1));
        }
    }
    return out;
}
template RowUncertaintyValues<double> CompiledSystem::row_uncertainty_at(const PhasePoint&, const std::vector<std::vector<double>>&) const;
template RowUncertaintyValues<WiderScalar> CompiledSystem::row_uncertainty_at(const PhasePointValues<WiderScalar>&, const std::vector<std::vector<WiderScalar>>&) const;

double CompiledSystem::solid_gibbs(const PhasePoint& prepared, std::size_t solid) const {
    const auto& spec = problem.solids[solid];
    double value = solids[solid].log_k + spec.molar_volume_m3_per_mol * (prepared.pressure - *spec.reference->reference_pressure_pa)
        / (GAS_CONSTANT_J_PER_MOL_K * prepared.temperature);
    for (std::size_t i = 0; i < spec.stoichiometry.size(); ++i) value += spec.stoichiometry[i] * prepared.standard.derivative[0][i];
    return value;
}

double CompiledSystem::solid_affinity(const PhasePoint& prepared, std::size_t solid) const {
    double value = -solid_gibbs(prepared, solid);
    for (std::size_t j = 0; j < rank && primal_size + j < prepared.point.size(); ++j)
        value += prepared.point[primal_size + j] * solids[solid].balance[j];
    return value;
}

CompiledSystem classify_solids(CompiledSystem compiled, const std::vector<double>& point) {
    for (auto& solid : compiled.solids)
        solid.present = solid.column && *solid.column < point.size() && point[*solid.column] > 0.0;
    return compiled;
}

std::vector<double> CompiledSystem::classified_residuals(
    const PhasePoint& prepared, const std::vector<double>& values
) const {
    std::vector<double> result = trace_moiety_relative(prepared, values);
    const auto numerical = log_charge_rows.empty() ? result : numerical_rows(prepared, values);
    for (const auto& [row, phase] : log_charge_rows)
        result[row] = std::copysign(std::max(wide_math::abs(result[row]), wide_math::abs(numerical[row])), numerical[row]);
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

std::vector<double> CompiledSystem::trace_moiety_relative(const PhasePoint& at, std::vector<double> values, std::size_t width) const {
    return trace_moiety_relative<double>(at, std::move(values), width);
}
template <typename Base>
std::vector<Base> CompiledSystem::trace_moiety_relative(
    const PhasePointValues<Base>& at, std::vector<Base> values, std::size_t width, const PhasePointValues<Base>* sizes
) const {
    const auto& prepared = sizes ? *sizes : at;
    const auto n_feed = feed_at(prepared.point);
    const double feed_total = [&] { if constexpr (std::is_same_v<Base, double>) return feed_total_at(prepared); else return this->feed_total; }();
    const Base trace = policy::EQUILIBRIUM_TRACE_MOIETY_FRACTION.value;
    if (values.size() != rows.size() * width) return values;
    for (std::size_t row = 0; row < rows.size(); ++row) {
        if (rows[row].kind != RowKind::Moiety) continue;
        const auto& b = B[rows[row].index];
        Base inventory = 0.0, terms = 0.0;
        for (std::size_t i = 0; i < b.size(); ++i) inventory += wide_math::abs(b[i]) * n_feed[i];
        if (inventory / feed_total > trace) continue;
        for (std::size_t p = 0; p < prepared.amount.size(); ++p)
            for (std::size_t i = 0; i < b.size(); ++i) terms += wide_math::abs(b[i]) * prepared.amount[p] * prepared.composition[p][i];
        for (const auto& solid : solids)
            if (solid.column) terms += wide_math::abs(solid.balance[rows[row].index] * prepared.point[*solid.column]);
        const Base size = std::max(inventory, terms) / feed_total;
        if (!(size > 0.0) || size > trace) continue;
        for (std::size_t j = 0; j < width; ++j) values[row * width + j] /= size;
    }
    return values;
}

template std::vector<WiderScalar> CompiledSystem::trace_moiety_relative(const PhasePointValues<WiderScalar>&, std::vector<WiderScalar>, std::size_t,
    const PhasePointValues<WiderScalar>*) const;

// Numerical charge representation only; physical equations and Ipopt derivatives remain unchanged.
std::vector<double> CompiledSystem::numerical_rows(
    const PhasePoint& prepared, std::vector<double> values, std::size_t width, const PhasePoint* sizes
) const {
    values = trace_moiety_relative<double>(prepared, std::move(values), width, sizes);
    for (const auto& [row, phase] : log_charge_rows) {
        const auto& z = model->parameters().charge_number;
        const auto& shape = prepared.shape[phase];
        const auto& support = problem.phases[phase].support;
        double logarithm[2];
        for (int sign = 0; sign < 2; ++sign) {
            double largest = -INFINITY, sum = 0.0;
            for (std::size_t j = 0; j < support.size(); ++j)
                if (z[support[j]] && (z[support[j]] > 0) == (sign == 0))
                    largest = std::max(largest, prepared.point[shape.x[j]] + std::log(std::abs(z[support[j]])));
            for (std::size_t j = 0; j < support.size(); ++j)
                if (z[support[j]] && (z[support[j]] > 0) == (sign == 0))
                    sum += std::exp(prepared.point[shape.x[j]] + std::log(std::abs(z[support[j]])) - largest);
            logarithm[sign] = largest + std::log(sum);
        }
        if (width == 1) values[row] = logarithm[0] - logarithm[1];
        else {
            std::fill_n(values.begin() + row * width, width, 0.0);
            for (std::size_t j = 0; j < support.size(); ++j) {
                const int charge = z[support[j]];
                if (charge) values[row * width + shape.x[j]] = std::copysign(std::exp(
                    prepared.point[shape.x[j]] + std::log(std::abs(charge)) - logarithm[charge < 0]), charge);
            }
        }
    }
    return values;
}

double CompiledSystem::classified_residual_norm(const std::vector<double>& point) const {
    const auto prepared = prepare(point);
    return maximum_abs(classified_residuals(prepared, evaluate_at(prepared).values));
}

std::vector<double> CompiledSystem::jacobian(const std::vector<double>& point) const {
    return jacobian_at(prepare(point));
}

std::vector<double> CompiledSystem::jacobian_at(const PhasePoint& at, bool constraints_only) const { return jacobian_at<double>(at, constraints_only); }

template <typename Base>
std::vector<Base> CompiledSystem::jacobian_at(const PhasePointValues<Base>& prepared, bool constraints_only) const {
    core::CalculationWork::equilibrium_entry();
    core::Scope scope("equilibrium_jacobian");
    if constexpr (std::is_same_v<Base, double>) {
    if (problem.loading) {
        const auto n = prepared.point.size();
        std::vector<double> out(rows.size() * n), unit(n);
        for (std::size_t j = 0; j < n; ++j) {
            unit[j] = 1;
            const auto series = row_series(prepared, {prepared.point, unit}, std::vector<double>(control_width()), 1);
            if (series.energy != PropertyActionStatus::Available) throw std::runtime_error("affine row action unavailable");
            for (std::size_t i = 0; i < rows.size(); ++i) out[i * n + j] = series.rows[1][i];
            unit[j] = 0;
        }
        return out;
    }
    if (ideal) return ideal_rows(*this, prepared).Jacobian(prepared.point);
    }
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
    const Base T = prepared.temperature;
    const Base P = prepared.pressure;
    const Base rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto& standard_temperature_derivative = prepared.standard.derivative[1];
    const auto& z = model->parameters().charge_number;
    std::vector<Base> result(rows.size() * n, 0.0);
    auto at = [&](std::size_t row, std::size_t column) -> Base& {
        return result[row * n + column];
    };
    std::vector<std::vector<Base>> columns;  // one first coefficient per coordinate serves every condition row
    for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
        const auto& row = rows[row_index];
        if (constraints_only && stationarity_row(row.kind)) continue;
        if (row.kind == RowKind::Chemical) {
            const auto& derivative = derivatives[row.phase];
            const std::size_t width = derivative.gradient.size();
            const auto& phase = problem.phases[row.phase];
            const std::size_t species = row.index;
            const Base physical_v = derivative.hessian[(2 + species) * width + 1]
                - 1.0 / volumes[row.phase];
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const Base physical_x = derivative.hessian[
                    (2 + species) * width + (2 + component)]
                    + (component == species ? 1.0 / compositions[row.phase][species] : 0.0);
                const Base dv = packings[row.phase].component[component] / prepared.eta[row.phase];
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
            if (charge_row[row.phase].has_value()) at(row_index, primal_size + rank + *charge_row[row.phase]) = -z[species];
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = derivatives[row.phase];
            const std::size_t width = derivative.gradient.size();
            const auto& phase = problem.phases[row.phase];
            const Base physical_v = -derivative.gradient[1]
                - volumes[row.phase] * derivative.hessian[width + 1] - P / rt;
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const Base dv = packings[row.phase].component[component] / prepared.eta[row.phase];
                at(row_index, shape[row.phase].x[local]) =
                    -volumes[row.phase] * derivative.hessian[width + 2 + component]
                        * compositions[row.phase][component] + physical_v * dv
                        + (model->parameters().ideal_gas ? compositions[row.phase][component] : 0.0);
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
            for (const auto& solid : solids)
                if (solid.column) at(row_index, *solid.column) = solid.balance[row.index] / feed_total;
        } else if (row.kind == RowKind::SolidAffinity) {
            const auto& solid = solids[row.index];
            if (!solid.present) at(row_index, *solid.column) = 1.0;
            for (std::size_t j = 0; solid.present && primal_size + j < n; ++j) at(row_index, primal_size + j) = -solid.balance[j];
        } else if (row.kind == RowKind::PrescribedComposition) {
            at(row_index, shape[0].x[row.index]) = 1;
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, shape[row.phase].x[local]) =
                    z[problem.phases[row.phase].support[local]]
                    * compositions[row.phase][problem.phases[row.phase].support[local]];
        } else if (row.kind == RowKind::Volume) {
            for (std::size_t p = 0; p < shape.size(); ++p) {
                const auto scale = amount[p] / *problem.volume_m3;
                at(row_index, *shape[p].amount) = volumes[p] / *problem.volume_m3;
                at(row_index, shape[p].volume) = -scale * volumes[p];
                for (std::size_t i = 0; i < shape[p].x.size(); ++i)
                    at(row_index, shape[p].x[i]) = scale * packings[p].component[problem.phases[p].support[i]] / prepared.eta[p];
                if (temperature) at(row_index, *temperature) += scale * packings[p].temperature_derivative / prepared.eta[p];
            }
        } else if (condition_row(row.kind)) {
            if constexpr (std::is_same_v<Base, double>) { std::vector<std::vector<Base>> state = {point, std::vector<Base>(n, 0.0)};
            const std::vector<Base> control(control_width(), 0.0);
            for (std::size_t c = columns.size(); c < n; ++c) {
                state[1][c] = 1.0;
                columns.push_back(row_series(prepared, state, control, 1).rows[1]);
                state[1][c] = 0.0;
            }
            for (std::size_t c = 0; c < n; ++c) at(row_index, c) = columns[c][row_index]; }
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

template std::vector<WiderScalar> CompiledSystem::jacobian_at(const PhasePointValues<WiderScalar>&, bool) const;

std::vector<double> CompiledSystem::lever_jacobian(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_lever_jacobian");
    if (problem.prescribed_liquid_composition) throw std::invalid_argument("feed actions unavailable for prescribed liquid composition");
    if (!reference_diagnostics.empty()) {
        throw ReferenceUnavailable("source reference T/P/parameter lever derivatives are not admitted");
    }
    auto prepared = prepare(point, false, nullptr, {}, {}, true);
    prepared.point.resize(std::max(prepared.point.size(), unknown_count()), 0.0);
    const std::size_t width = control_width();
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
    for (std::size_t s = 0; s < solids.size(); ++s)
        if (solids[s].column) result[*solids[s].column] = solid_gibbs(prepared, s);
    return result;
}

std::vector<double> CompiledSystem::hessian(
    const std::vector<double>& point,
    const std::vector<double>& constraint_multipliers
) const {
    return hessian_at(prepare(point), constraint_multipliers);
}

std::vector<double> CompiledSystem::hessian_at(
    const PhasePoint& prepared, const std::vector<double>& constraint_multipliers, bool physical
) const {
    core::CalculationWork::equilibrium_entry();
    core::Scope scope("equilibrium_hessian");
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (feasibility_cell() && !physical) return constraint_hessian_at(prepared, constraint_multipliers);
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
        // Per-mole chemical potentials and pressure gradient; the extensive ones are N_p times them.
        const double pressure_gradient = derivative.gradient[1] - sum_x / physical_volume + P / rt;
        const double physical_v_gradient = phase_amount_value * pressure_gradient;
        const double physical_vv = phase_amount_value * (
            derivative.hessian[width + 1] + sum_x / (physical_volume * physical_volume));
        std::vector<double> mu(spec.support.size()), h_xv(spec.support.size());
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const std::size_t si = spec.support[i];
            mu[i] = derivative.gradient[2 + si]
                + std::log(std::max(x[si] / physical_volume, std::numeric_limits<double>::min()))
                + standard_value[si];
            h_xv[i] = phase_amount_value * (derivative.hessian[(2 + si) * width + 1] - 1.0 / physical_volume);
        }
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const std::size_t si = spec.support[i];
            const double dv_i = packing.component[si] / eta;
            for (std::size_t j = 0; j < spec.support.size(); ++j) {
                const std::size_t sj = spec.support[j];
                const double dv_j = packing.component[sj] / eta;
                // Keep the ideal diagonal analytic and apply residual curvature before trace fractions multiply.
                const double h_xx = phase_amount_value * (x[si]
                    * (x[sj] * derivative.hessian[(2 + si) * width + (2 + sj)]))
                    + (si == sj ? phase_amount_value * x[si] : 0.0);
                put(layout.x[i], layout.x[j],
                    h_xx + x[si] * h_xv[i] * dv_j
                    + dv_i * h_xv[j] * x[sj] + dv_i * physical_vv * dv_j
                    + (si == sj ? phase_amount_value * mu[i] * x[si]
                                      + physical_v_gradient * dv_i : 0.0));
            }
            put(layout.x[i], volume,
                -physical_volume * (x[si] * h_xv[i] + dv_i * physical_vv)
                - physical_v_gradient * dv_i);
        }
        put(volume, volume, physical_volume * physical_volume * physical_vv
            + physical_v_gradient * physical_volume);
        if (layout.amount.has_value()) {
            for (std::size_t i = 0; i < spec.support.size(); ++i)
                put(*layout.amount, layout.x[i], mu[i] * x[spec.support[i]]
                    + pressure_gradient * packing.component[spec.support[i]] / eta);
            put(*layout.amount, volume, -physical_volume * pressure_gradient);
        }
    }
    if (!constraint_multipliers.empty()) {
        const auto weighted = constraint_hessian_at(prepared, constraint_multipliers, physical);
        for (std::size_t index = 0; index < result.size(); ++index) result[index] += weighted[index];
    }
    scope.probe("entries", static_cast<double>(result.size()));
    return result;
}

RowSeries CompiledSystem::row_series(
    const PhasePoint& prepared, const std::vector<std::vector<double>>& state,
    const std::vector<double>& control, std::size_t order, bool plane, bool parameter_pairs, const std::vector<std::size_t>& parameter_union
) const {
    core::CalculationWork::equilibrium_entry();
    core::Scope scope("equilibrium_row_series");
    if (ideal) throw std::invalid_argument("ideal solved-state actions unavailable");
    const std::size_t n = prepared.point.size();
    const std::size_t components = model->parameters().component_count();
    const std::size_t width = control_width(), extended = extended_control_width();
    if (state.size() != order + 1 ||
        (control.size() != width && control.size() != extended - 2 - problem.volume_m3.has_value() - problem.internal_energy_j.has_value() && control.size() != extended) ||
        order > (mode_count ? 2u : 3u))  // a third critical coefficient needs E6
        throw std::invalid_argument("row series shape");
    for (const auto& item : state)
        if (item.size() != n) throw std::invalid_argument("row series state");
    const bool pressure_jets = problem.pressure.free || control[1] != 0.0;
    const bool temperature_jets = problem.temperature.free || control[0] != 0.0;
    const auto active_control = control.begin() + static_cast<std::ptrdiff_t>(2 + components);
    std::vector<std::size_t> parameters;
    for (std::size_t a = 0; a < width - 2 - components; ++a)
        if (active_control[a] != 0.0) parameters.push_back(a);
    // Polarised paths need parameter pairs only when two directions move parameters.
    const auto& pairs = parameter_pairs ? parameters : std::vector<std::size_t>{};
    for (std::size_t s = 0; control[0] != 0.0 && s < solids.size(); ++s)
        if (!solids[s].log_k_slope)
            throw ReferenceUnavailable("solid " + problem.solids[s].name + " has no admitted ln K temperature slope");
    if (!reference_diagnostics.empty() && order > 2 && (!parameters.empty() || pressure_jets || temperature_jets))
        throw ReferenceUnavailable("row series pressure/T/parameter order is unsupported");
    // A batch shares requested columns. Only fixed-intensive reaction directions also share mapped slopes.
    const bool reuse = !reference_diagnostics.empty() && control[0] == 0.0 && control[1] == 0.0;
    const bool shared = !parameters.empty() && order == 1 && !parameter_pairs
        && !parameter_union.empty() && parameter_union.back() < width - 2 - components
        && std::is_sorted(parameter_union.begin(), parameter_union.end())
        && std::includes(parameter_union.begin(), parameter_union.end(), parameters.begin(), parameters.end());
    auto& cached = prepared.parameter_standard;
    if (reuse && shared && (!cached || cached->positions != parameter_union)) {
        core::Scope scope("parameter_standard");
        try {
            cached = ParameterStandard{parameter_union, standard_term(*this, prepared.temperature, prepared.pressure, order,
                pressure_jets, temperature_jets, parameter_union)};
        } catch (const std::invalid_argument& error) {
            cached = ParameterStandard{parameter_union, ParameterStandardRefused{}};
            scope.probe("parameter_standard_refused", {}, 0, error.what());
        }
    }
    const auto* mapped = reuse && shared ? std::get_if<StandardTerm>(&cached->outcome) : nullptr;
    const auto& standard = mapped ? *mapped : reuse && parameters.empty() ? prepared.standard
        : standard_term(*this, prepared.temperature, prepared.pressure, order, pressure_jets, temperature_jets,
                        parameters, nullptr, pairs, shared ? parameter_union : parameters);
    if (prepared.calorics) {
        // Pressure and parameter columns join the anchored g at this centre on first need; a motion that reads an
        // unavailable column is refused with its cause.
        const auto& gibbs = prepared.calorics->gibbs;
        const bool has_pressure = std::any_of(gibbs.begin(), gibbs.end(), [](const auto& e) { return e.pressure_order > 0; });
        const bool has_parameters = std::any_of(gibbs.begin(), gibbs.end(), [](const auto& e) { return e.parameter.has_value(); });
        if ((control[1] != 0.0 && !has_pressure) || (!parameters.empty() && !has_parameters))
            prepared.calorics = caloric_reference(*this, prepared.temperature, prepared.pressure,
                                                  has_pressure || control[1] != 0.0, has_parameters || !parameters.empty());
        for (const auto& e : prepared.calorics->gibbs)
            if (e.column.status != SolvedStateActionStatus::Available
                && ((e.pressure_order && control[1] != 0.0) || (e.parameter && active_control[*e.parameter] != 0.0)))
                throw ReferenceUnavailable("caloric " + e.column.message);
    }
    auto geometry = record_geometry(*this, prepared, standard, false, plane);
    const std::size_t stride = width + 2 + (has_caloric_condition() ? 3 : 0) + (mode_count ? 2 * components : 0);
    std::vector<std::vector<double>> output(order + 1);
    for (std::size_t k = 0; k <= order; ++k) {
        auto input = k == 0 ? controls(*this, prepared) : std::vector<double>(geometry.Domain(), 0.0);
        std::copy(state[k].begin(), state[k].end(), input.begin());
        if (k == 1) {
            std::copy_n(control.begin(), width, input.begin() + static_cast<std::ptrdiff_t>(n));
            if (control.size() == extended) {
                const auto target = extended - 2 - problem.volume_m3.has_value() - problem.internal_energy_j.has_value();
                if (problem.volume_m3) input[input.size()-1-problem.internal_energy_j.has_value()] = control[target+2];
                if (problem.internal_energy_j) input.back() = control.back();
                if (has_caloric_condition()) { input[n+width] = control[target]; input[n+width+1] = control[target+1]; }
            }
        }
        output[k] = geometry.Forward(k, input);
    }
    RowSeries result;
    result.rows.assign(order + 1, std::vector<double>(rows.size(), 0.0));
    result.temperature.assign(order + 1, 0.0);
    result.pressure.assign(order + 1, 0.0);
    if (problem.loading) for (const auto& coefficient : state) result.loading.push_back(coefficient[loading_column()]);
    result.phase.resize(problem.phases.size());
    if (plane) result.standard_minus_plane.assign(components, std::vector<double>(order + 1));
    const std::size_t row_offset = problem.phases.size() * stride;
    const std::size_t temperature_index = primal_size - free_intensive_count;
    result.solid_amount.assign(solids.size(), std::vector<double>(order + 1));
    result.solid_affinity = result.solid_amount;
    for (std::size_t k = 0; k <= order; ++k) {
        if (plane) for (std::size_t i = 0; i < components; ++i)
            result.standard_minus_plane[i][k] = output[k][row_offset + rows.size() + 2 * solids.size() + i];
        for (std::size_t row = 0; row < rows.size(); ++row) result.rows[k][row] = output[k][row_offset + row];
        for (std::size_t s = 0; s < solids.size(); ++s) {
            result.solid_amount[s][k] = output[k][row_offset + rows.size() + 2 * s];
            result.solid_affinity[s][k] = output[k][row_offset + rows.size() + 2 * s + 1];
        }
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
        SpeciesSlopeMap slopes(prepared.temperature, components, problem.reactions);
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
            const auto species = slopes(targets);
            if (plane) for (std::size_t i = 0; i < components; ++i) result.standard_minus_plane[i][k] += species[i];
            for (std::size_t row = 0; row < rows.size(); ++row)
                if (rows[row].kind == RowKind::Chemical)
                    result.rows[k][row] += species[rows[row].index];
        }
    }
    // EqID solid_saturation_index: a source ln K_s increment lowers A_s/RT linearly in t (solids fix T).
    for (std::size_t s = 0; order > 0 && control.size() > width && s < solids.size(); ++s) {
        const double shift = control[width + 4 * problem.reactions.size() + s];
        result.solid_affinity[s][1] -= shift;
        for (std::size_t row = 0; solids[s].present && row < rows.size(); ++row)
            if (rows[row].kind == RowKind::SolidAffinity && rows[row].index == s) result.rows[1][row] += shift;
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
        // Central inputs are fixed by the prepared state; a key preserves raw direction bits and shape.
        const auto action_key = [&](const std::vector<std::vector<double>>& directions) {
            std::string key(reinterpret_cast<const char*>(&phase), sizeof(phase));
            for (const auto& direction : directions) {
                const auto size = direction.size();
                key.append(reinterpret_cast<const char*>(&size), sizeof(size));
                key.append(reinterpret_cast<const char*>(direction.data()), size * sizeof(double));
            }
            return key;
        };
        const auto with_unit = [&](std::vector<std::vector<double>> directions, std::optional<std::size_t> component) {
            if (component) directions.emplace_back(width, 0.0).at(*component) = 1.0;
            return directions;
        };
        // D3F[v, v, e_j] for every needed j from one model call, kept with the prepared state as the contraction of each j.
        const auto contract_together = [&](const std::vector<double>& v) {
            if (prepared.energy_actions.count(action_key(with_unit({v, v}, needed.front())))) return;
            const auto energy = model->repeated_direction_energy_gradient(prepared.temperature, prepared.volume[phase],
                prepared.composition[phase], v, needed);
            for (std::size_t i = 0; i < needed.size(); ++i) prepared.energy_actions.emplace(action_key(with_unit({v, v}, needed[i])),
                ResidualEnergyAction{energy.gradient.empty() ? std::nullopt : std::optional(energy.gradient[i]), energy.status, energy.diagnostics});
        };
        const auto contraction = [&](std::vector<std::vector<double>> directions, std::optional<std::size_t> component) {
            directions = with_unit(std::move(directions), component);
            auto key = action_key(directions);
            const auto found = prepared.energy_actions.find(key);
            const bool reused = found != prepared.energy_actions.end();
            const auto& energy = reused ? found->second : prepared.energy_actions.emplace(
                std::move(key), model->residual_energy_action(prepared.temperature,
                    prepared.volume[phase], prepared.composition[phase], directions)).first->second;
            result.maximum_energy_order = std::max(result.maximum_energy_order, directions.size());
            if (energy.status != PropertyActionStatus::Available && result.energy == PropertyActionStatus::Available)
                result.energy = energy.status;
            return energy.action.value_or(std::numeric_limits<double>::quiet_NaN());
        };
        for (std::size_t k = 1; k <= order; ++k) {
            for (std::size_t i = 0; i < width; ++i)
                for (std::size_t j = 0; j < width; ++j)
                    series.gradient[k][i] += derivative.hessian[i * width + j] * series.point[k][j];
            if (k == 2) {
                contract_together(series.point[1]);
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
        if (mode_count) {
            // Residual parts of EqID critical_equations through order 2: d_i [F_nn(p) w]_i (E2-E4) and
            // F_nnn(p)[w,w,w] (E3-E5) along p(t), w(t), d(t); a zero direction contributes exactly zero.
            const std::size_t offset = phase * stride + width + 2, mode = mode_row();
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
            if (problem.critical) {
            std::array<double, 3> c{F({w[0], w[0], w[0]}), 0.0, 0.0};
            if (order > 0) c[1] = 3.0 * F({w[0], w[0], w[1]}) + F({w[0], w[0], w[0], p[1]});
            if (order > 1) c[2] = 3.0 * F({w[0], w[0], w[2]}) + 3.0 * F({w[0], w[1], w[1]}) + 3.0 * F({w[0], w[0], w[1], p[1]})
                + F({w[0], w[0], w[0], p[2]}) + 0.5 * F({w[0], w[0], w[0], p[1], p[1]});
            for (std::size_t k = 0; k <= order; ++k) result.rows[k][mode + components + 1] += c[k];
            }
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
    const PhasePoint& prepared, const std::vector<double>& constraint_multipliers, bool physical
) const {
    core::CalculationWork::equilibrium_entry();
    if (ideal) {
        if (!constraint_multipliers.empty() && constraint_multipliers.size() != rows.size())
            throw std::invalid_argument("constraint multiplier count");
        return ideal_rows(*this, prepared).Hessian(prepared.point,
            constraint_multipliers.empty() ? std::vector<double>(rows.size(), 0.0) : constraint_multipliers);
    }
    if (feasibility_cell() && !physical) return pinned_hessian(*this, prepared, constraint_multipliers);
    const auto& point = prepared.point;
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (constraint_multipliers.empty()) return result;
    std::size_t multiplier = 0;
    for (const auto& row : rows) {
        if (stationarity_row(row.kind)) continue;
        if (multiplier >= constraint_multipliers.size()) {
            throw std::invalid_argument("constraint multiplier count");
        }
        const double weight = constraint_multipliers[multiplier++];
        if (row.kind == RowKind::Volume) continue;
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
