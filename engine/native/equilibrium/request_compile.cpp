#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"
#include "core/trace.hpp"
#include "contributions/helmholtz.hpp"
#include "properties/helmholtz_tape.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
    result.parameter_derivative.resize(model.active_parameters().size());
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
    for (std::size_t parameter = 0; parameter < model.active_parameters().size(); ++parameter) {
        const auto& item = model.active_parameters()[parameter];
        for (const std::size_t i : item.indices) {
            if (item.family == ParameterFamily::SegmentCount) {
                result.parameter_derivative[parameter] +=
                    scale * composition[i] * diameter[i] * diameter[i] * diameter[i];
            } else if (item.family == ParameterFamily::SigmaConstant
                       && parameters.charge_number[i] == 0) {
                result.parameter_derivative[parameter] += scale * composition[i]
                    * parameters.segment_count[i] * 3.0 * diameter[i] * diameter[i]
                    * diameter[i] / sigma[i];
            } else if (item.family == ParameterFamily::DispersionEnergyOverK
                       && parameters.charge_number[i] == 0) {
                result.parameter_derivative[parameter] += scale * composition[i]
                    * parameters.segment_count[i] * 3.0 * diameter[i] * diameter[i]
                    * sigma[i] * 0.36 / temperature * std::exp(
                        -3.0 * parameters.dispersion_energy_over_k[i] / temperature);
            }
        }
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

StandardTerm standard_term(const CompiledSystem& compiled, double temperature, std::size_t order) {
    const auto& problem = compiled.problem;
    const std::size_t components = compiled.model->parameters().component_count();
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
    for (std::size_t k = 0; k <= order; ++k) {
        for (std::size_t species = 0; species < components; ++species) {
            double value = 0.0;
            for (std::size_t j = 0; j <= k; ++j) value += binomial[k][j] * gibbs[j][species] * u[k - j];
            result.derivative[k][species] = value / GAS_CONSTANT_J_PER_MOL_K;
        }
    }
    return result;
}

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
    const auto charges = species_z(parameters);
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
                + 0.5 * standard.derivative[2][species] * delta_T * delta_T;
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
        // ponytail: quadratic contraction calls; batch them if measured solve cost requires it.
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t column = 0; column <= row; ++column) {
                const auto energy = model.residual_energy_action(temperature, values[1], amounts,
                    {contraction, directions[row], directions[column]});
                if (!energy.action.has_value()) {
                    throw std::runtime_error(std::string("third directional energy contraction unavailable: ")
                        + solved_state_status_name(solved_state_status(energy.status)));
                }
                const double value = mapped[row * n + column] + *energy.action
                    - pressure_weight * (directions[row][1] * volume_hessian_action[column]
                        + directions[column][1] * volume_hessian_action[row]);
                result[row * n + column] += value;
                if (row != column) result[column * n + row] += value;
            }
        }
    }
    return result;
}

std::vector<double> active_values(const HelmholtzModel& model) {
    const auto& p = model.parameters();
    std::vector<double> result;
    for (const auto& item : model.active_parameters()) {
        const std::vector<double>& family = item.family == ParameterFamily::SegmentCount ? p.segment_count
            : item.family == ParameterFamily::SigmaConstant ? p.sigma_constant
            : item.family == ParameterFamily::DispersionEnergyOverK ? p.dispersion_energy_over_k
            : item.family == ParameterFamily::Kij ? p.k_ij : item.family == ParameterFamily::Lij ? p.l_ij
            : item.family == ParameterFamily::AssociationEnergyOverK ? p.association_energy_over_k
            : p.association_volume;
        result.push_back(family[item.indices.front()]);
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
    for (const double value : active_values(model)) u.push_back(value);
    CppAD::Independent(u);
    const AD T = problem.temperature.free ? u[compiled.primal_size - compiled.free_intensive_count] : u[n];
    const AD P = problem.pressure.free ? u[compiled.primal_size - 1] : u[n + 1];
    AD feed_total = AD(0.0);
    for (std::size_t i = 0; i < components; ++i) feed_total += u[n + 2 + i];
    const std::vector<AD> theta(u.begin() + static_cast<std::ptrdiff_t>(n + 2 + components), u.end());
    const auto tape_parameters = detail::make_tape_parameters(parameters, model.active_parameters(), theta);
    const auto sigma = segment_sigma(T, tape_parameters);
    const auto diameter = segment_diameter(T, tape_parameters, sigma);
    const auto charges = species_z(parameters);
    std::vector<std::optional<std::size_t>> phi(problem.phases.size());
    std::size_t phi_offset = 0;
    for (const auto& row : compiled.rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = compiled.primal_size + compiled.rank + phi_offset++;
    std::vector<AD> output;
    std::vector<AD> row_values(compiled.rows.size(), AD(0.0));
    const AD delta_T = T - prepared.temperature;
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
        for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
            const auto& item = compiled.rows[row];
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
    output.insert(output.end(), row_values.begin(), row_values.end());
    return CppAD::ADFun<double>(u, output);
}

std::string represented(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return buffer;
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
                                    + represented(pin_sum));
    }
    if (pinned == compiled.problem.phases.size() && pin_sum != 1.0) {
        throw std::invalid_argument("all-pinned phase fractions must sum to one: S_pin="
                                    + represented(pin_sum));
    }
}

}  // namespace

CompiledSystem compile_problem(
    const HelmholtzModel& model, const EquilibriumProblem& request
) {
    core::Scope scope("compile_problem");
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
    const auto z = species_z(model.parameters());
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
    compiled.B = balance_basis(components, request.reactions);
    compiled.rank = compiled.B.size();
    compiled.conserves_total_moles = conserves_total_moles(compiled.B);
    const double temperature_guess = request.temperature.seed_or(300.0);
    compiled.standard_gibbs = standard_state(
        temperature_guess, components, request.reactions, request.standard_gibbs_j_per_mol
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
    admit_pins(compiled);
    std::vector<bool> covered(components, false);
    std::vector<bool> participant(components, false);
    for (const auto& item : request.reactions) {
        for (std::size_t i = 0; i < components; ++i)
            participant[i] = participant[i] || item.stoichiometry[i] != 0.0;
    }
    for (auto& phase : compiled.problem.phases) {
        if (phase.support_declared && phase.support.empty()) {
            throw std::invalid_argument("empty support in phase " + phase.name);
        }
        phase.support = complete_support(components, phase.support);
        phase.support_declared = true;
        const auto support_before = phase.support;
        phase.support.erase(std::remove_if(
            phase.support.begin(), phase.support.end(), [&](std::size_t index) {
                return compiled.n_feed[index] == 0.0 && !participant[index];
            }), phase.support.end());
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
        if (!covered[i] && (compiled.n_feed[i] != 0.0 || participant[i])) {
            throw std::invalid_argument("species " + std::to_string(i) + " has no support");
        }
    }
    compiled.pinned_count = 0;
    for (const auto& phase : compiled.problem.phases) {
        compiled.pinned_count += phase.amount.kind == AmountKind::Pinned;
    }
    compiled.free_intensive_count = request.temperature.free + request.pressure.free;
    if (compiled.free_intensive_count != compiled.pinned_count) {
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
    if (compiled.problem.phases.size() - compiled.pinned_count > compiled.rank) {
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
    const auto shape = layouts(compiled);
    compiled.primal_size = 0;
    for (const auto& item : shape) compiled.primal_size = std::max(compiled.primal_size, item.volume + 1);
    compiled.primal_size += compiled.free_intensive_count;
    compiled.rows_minus_unknowns = static_cast<long>(compiled.pinned_count)
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

std::vector<double> CompiledSystem::generic_point() const {
    core::Scope scope("equilibrium_generic_point");
    const auto shape = layouts(*this);
    std::vector<double> point(primal_size + rank + charge_row_count, 0.0);
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        double total = 0.0;
        for (std::size_t i = 0; i < shape[phase].x.size(); ++i) {
            const double value = 1.0 + 0.37 * static_cast<double>(i)
                + 0.11 * static_cast<double>(phase);
            point[shape[phase].x[i]] = value;
            total += value;
        }
        for (const std::size_t coordinate : shape[phase].x)
            point[coordinate] = std::log(point[coordinate] / total);
        if (shape[phase].amount.has_value()) point[*shape[phase].amount] = 1.0;
        point[shape[phase].volume] = std::log(0.25 + 0.05 * static_cast<double>(phase));
    }
    std::size_t offset = primal_size - free_intensive_count;
    if (problem.temperature.free) point[offset++] = problem.temperature.seed_or(300.0);
    if (problem.pressure.free) point[offset++] = problem.pressure.seed_or(1.0e5);
    return point;
}

PhasePoint CompiledSystem::prepare(const std::vector<double>& point) const {
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
        result.volume.push_back(result.packing.back().coefficient / result.eta.back());
        result.amount.push_back(layout.amount.has_value()
            ? point[*layout.amount] : phase_amount(problem.phases[phase], n_feed));
        result.derivative.push_back(model->evaluate(
            result.temperature, result.volume.back(), values, DerivativeOrder::Hessian));
        scope.probe("gradient_floor", result.derivative.back().gradient_floor,
                    static_cast<std::uint32_t>(phase));
        result.composition.push_back(std::move(values));
    }
    result.standard = standard_term(*this, result.temperature, 2);
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
    const auto z = species_z(model->parameters());
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
    if (pinned_count != 0 || values.size() != rows.size()) return result;
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
    const auto z = species_z(model->parameters());
    std::vector<std::optional<std::size_t>> phi(problem.phases.size());
    std::size_t phi_offset = 0;
    for (const auto& row : rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = rank + phi_offset++;
    std::vector<double> result(rows.size() * n, 0.0);
    auto at = [&](std::size_t row, std::size_t column) -> double& {
        return result[row * n + column];
    };
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
    const auto prepared = prepare(point);
    const std::size_t components = model->parameters().component_count();
    const double T = prepared.temperature;
    const double P = prepared.pressure;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto& standard_temperature_derivative = prepared.standard.derivative[1];
    std::size_t levers = 0;
    std::vector<double> result;
    std::vector<double> pinned_amount_feed_term(rank, 0.0);
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& composition = prepared.composition[phase];
        const auto& packing = prepared.packing[phase];
        const double eta = prepared.eta[phase];
        const double volume = prepared.volume[phase];
        const auto& derivative = prepared.derivative[phase];
        const std::size_t width = derivative.gradient.size();
        if (result.empty()) {
            levers = width;
            result.resize(rows.size() * levers, 0.0);
        }
        if (problem.phases[phase].amount.kind == AmountKind::Pinned) {
            for (std::size_t balance_row = 0; balance_row < rank; ++balance_row) {
                for (std::size_t species = 0; species < components; ++species) {
                    pinned_amount_feed_term[balance_row] +=
                        problem.phases[phase].amount.fraction * B[balance_row][species]
                        * composition[species] / feed_total;
                }
            }
        }
        for (std::size_t row = 0; row < rows.size(); ++row) {
            if (rows[row].phase != phase) continue;
            if (rows[row].kind == RowKind::Chemical) {
                const double physical_v = derivative.hessian[
                    (2 + rows[row].index) * width + 1] - 1.0 / volume;
                result[levers * row] = derivative.hessian[(2 + rows[row].index) * width]
                    + physical_v * packing.temperature_derivative / eta
                    + standard_temperature_derivative[rows[row].index];
                for (std::size_t active = 2 + components; active < width; ++active) {
                    const double active_parameter_column = derivative.hessian[
                        (2 + rows[row].index) * width + active];
                    // EqID: packing_coefficient_parameter_derivative
                    const double packing_parameter_term = packing.parameter_derivative[
                        active - 2 - components] / eta;
                    result[levers * row + active] = active_parameter_column
                        + physical_v * packing_parameter_term;
                }
            } else if (rows[row].kind == RowKind::Pressure) {
                const double physical_v = -derivative.gradient[1]
                    - volume * derivative.hessian[width + 1] - P / rt;
                result[levers * row] = -volume * derivative.hessian[width]
                    + P * volume / (GAS_CONSTANT_J_PER_MOL_K * T * T)
                    + physical_v * packing.temperature_derivative / eta;
                result[levers * row + 1] = -packing.coefficient / eta / rt;
                for (std::size_t active = 2 + components; active < width; ++active) {
                    const double active_parameter_column =
                        -volume * derivative.hessian[width + active];
                    const double packing_parameter_term = packing.parameter_derivative[
                        active - 2 - components] / eta;
                    result[levers * row + active] = active_parameter_column
                        + physical_v * packing_parameter_term;
                }
            }
        }
    }
    for (std::size_t row = 0; row < rows.size(); ++row) {
        if (rows[row].kind == RowKind::Moiety) {
            for (std::size_t species = 0; species < components; ++species)
                result[levers * row + 2 + species] = -B[rows[row].index][species]
                    / feed_total + pinned_amount_feed_term[rows[row].index];
        }
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
    if (pinned_count != 0) return result;
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
    if (pinned_count != 0) return constraint_hessian_at(prepared, constraint_multipliers);
    const std::size_t components = model->parameters().component_count();
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
    if (state.size() != order + 1 || control.size() != width || order > 3)
        throw std::invalid_argument("row series shape");
    for (const auto& item : state)
        if (item.size() != n) throw std::invalid_argument("row series state");
    auto geometry = record_geometry(*this, prepared, standard_term(*this, prepared.temperature, 3));
    const std::size_t stride = width + 2;
    std::vector<std::vector<double>> output(order + 1);
    for (std::size_t k = 0; k <= order; ++k) {
        std::vector<double> input(n + width, 0.0);
        std::copy(state[k].begin(), state[k].end(), input.begin());
        if (k == 0) {
            input[n] = prepared.temperature;
            input[n + 1] = prepared.pressure;
            std::copy(n_feed.begin(), n_feed.end(), input.begin() + static_cast<std::ptrdiff_t>(n + 2));
            const auto values = active_values(*model);
            std::copy(values.begin(), values.end(), input.begin() + static_cast<std::ptrdiff_t>(n + 2 + components));
        } else if (k == 1) {
            std::copy(control.begin(), control.end(), input.begin() + static_cast<std::ptrdiff_t>(n));
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
    // Residual gradient coefficients: the Hessian action of each coefficient plus the model's
    // third and fourth contractions of the lower coefficients (Faa di Bruno through order 3).
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        auto& series = result.phase[phase];
        const auto& derivative = prepared.derivative[phase];
        series.gradient.assign(order + 1, std::vector<double>(width, 0.0));
        series.gradient[0] = derivative.gradient;
        std::vector<std::size_t> needed = {1};
        for (const std::size_t species : problem.phases[phase].support) needed.push_back(2 + species);
        const auto contraction = [&](std::vector<std::vector<double>> directions, std::size_t component) {
            std::vector<double> unit(width, 0.0);
            unit[component] = 1.0;
            directions.push_back(unit);
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
    if (pinned_count != 0) return pinned_hessian(*this, prepared, constraint_multipliers);
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
                    coefficient = species_z(model->parameters())[species];
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
