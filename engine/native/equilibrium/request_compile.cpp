#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"
#include "core/trace.hpp"
#include "contributions/helmholtz.hpp"

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

struct Layout final {
    std::vector<std::size_t> x;
    std::optional<std::size_t> amount;
    std::size_t volume = 0;
};

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
    if (!std::isfinite(phase.amount.fraction) || phase.amount.fraction < 0.0
        || phase.amount.fraction > 1.0) {
        throw std::invalid_argument("invalid phase amount");
    }
    const double total = vector_norm(feed);
    return phase.amount.fraction * total;
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

struct PackingTerms final {
    double coefficient = 0.0;
    double temperature_derivative = 0.0;
    std::vector<double> component;
    std::vector<double> parameter_derivative;
};

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

struct StandardTerm final {
    std::vector<double> value;
    std::vector<double> temperature_derivative;
};

StandardTerm standard_term(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied
) {
    const auto gibbs = standard_state(temperature, components, reactions, supplied);
    std::vector<double> gibbs_temperature_derivative(components, 0.0);
    if (supplied.empty() && !reactions.empty()) {
        auto derivative_reactions = reactions;
        for (auto& reaction : derivative_reactions) {
            if (reaction.affinity_j_per_mol.has_value()) reaction.affinity_j_per_mol = 0.0;
        }
        gibbs_temperature_derivative = standard_state(
            1.0, components, derivative_reactions, {}
        );
    }
    StandardTerm result;
    result.value.resize(components);
    result.temperature_derivative.resize(components);
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    for (std::size_t species = 0; species < components; ++species) {
        result.value[species] = gibbs[species] / rt;
        result.temperature_derivative[species] =
            gibbs_temperature_derivative[species] / rt
            - gibbs[species] / (GAS_CONSTANT_J_PER_MOL_K * temperature * temperature);
    }
    return result;
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
            const std::string name = request.reactions[row].name.empty()
                ? std::to_string(row) : request.reactions[row].name;
            throw std::invalid_argument("reaction " + name + " violates z");
        }
    }
    compiled.B = balance_basis(components, request.reactions);
    compiled.rank = compiled.B.size();
    const double temperature_guess = request.temperature.seed_or(300.0);
    compiled.standard_gibbs = standard_state(
        temperature_guess, components, request.reactions, request.standard_gibbs_j_per_mol
    );
    compiled.feed_composition = normalized_feed(compiled.n_feed);
    compiled.balance = matvec(compiled.B, compiled.n_feed);
    std::vector<bool> covered(components, false);
    std::vector<bool> produced(components, false);
    for (const auto& item : request.reactions) {
        for (std::size_t i = 0; i < components; ++i) produced[i] = produced[i] || item.stoichiometry[i] > 0.0;
    }
    for (auto& phase : compiled.problem.phases) {
        phase.support = complete_support(components, phase.support);
        const auto support_before = phase.support;
        phase.support.erase(std::remove_if(
            phase.support.begin(), phase.support.end(), [&](std::size_t index) {
                return compiled.n_feed[index] == 0.0 && !produced[index];
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
        if (!covered[i] && (compiled.n_feed[i] != 0.0 || produced[i])) {
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
    bool all_pinned = compiled.pinned_count == compiled.problem.phases.size();
    if (all_pinned) {
        throw std::invalid_argument(
            "no free phase carries the feed; phase "
            + compiled.problem.phases.front().name + " is pinned"
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
    std::size_t offset = 0;
    for (const auto& phase : problem.phases) {
        offset += phase.support.size() + 1;
        if (phase.amount.kind == AmountKind::Free) ++offset;
    }
    if (problem.temperature.free) point[offset++] = problem.temperature.seed_or(300.0);
    if (problem.pressure.free) point[offset++] = problem.pressure.seed_or(1.0e5);
    return point;
}

PointEvaluation CompiledSystem::evaluate(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_evaluate");
    const auto shape = layouts(*this);
    if (point.size() < primal_size) {
        throw std::invalid_argument("equilibrium point is too short");
    }
    const std::size_t components = model->parameters().component_count();
    const auto z = species_z(model->parameters());
    std::size_t offset = primal_size - free_intensive_count;
    const double temperature = problem.temperature.free ? point[offset++] : problem.temperature.value;
    const double pressure = problem.pressure.free ? point[offset++] : problem.pressure.value;
    std::vector<std::vector<double>> x;
    std::vector<double> amount;
    std::vector<double> volume;
    x.reserve(shape.size());
    amount.reserve(shape.size());
    volume.reserve(shape.size());
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        std::vector<double> values(components, 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < shape[phase].x.size(); ++i) {
            const double value = std::exp(point[shape[phase].x[i]]);
            if (!std::isfinite(value) || value <= 0.0) throw std::invalid_argument("invalid composition");
            values[problem.phases[phase].support[i]] = value;
            total += value;
        }
        if (!(total > 0.0)) throw std::invalid_argument("composition has zero total");
        const auto packing = packing_terms(*model, temperature, values);
        x.push_back(std::move(values));
        amount.push_back(shape[phase].amount.has_value()
            ? point[*shape[phase].amount] : phase_amount(problem.phases[phase], n_feed));
        volume.push_back(packing.coefficient / std::exp(point[shape[phase].volume]));
    }
    std::vector<double> lambda(rank, 0.0);
    for (std::size_t row = 0; row < rank && primal_size + row < point.size(); ++row) {
        lambda[row] = point[primal_size + row];
    }
    std::vector<double> phi(charge_row_count, 0.0);
    for (std::size_t row = 0; row < charge_row_count
         && primal_size + rank + row < point.size(); ++row) {
        phi[row] = point[primal_size + rank + row];
    }
    const double rt = GAS_CONSTANT_J_PER_MOL_K * temperature;
    const auto standard = standard_term(
        temperature, components, problem.reactions, problem.standard_gibbs_j_per_mol
    );
    std::vector<std::optional<HelmholtzDerivatives>> derivatives(shape.size());
    for (std::size_t phase = 0; phase < shape.size(); ++phase) {
        derivatives[phase] = model->evaluate(temperature, volume[phase], x[phase], DerivativeOrder::Hessian);
        scope.probe("gradient_floor", derivatives[phase]->gradient_floor,
                    static_cast<std::uint32_t>(phase));
    }
    PointEvaluation result;
    for (const auto& derivative : derivatives) result.floor = std::max(result.floor, derivative->gradient_floor);
    std::vector<std::optional<std::size_t>> phase_phi(shape.size());
    std::size_t phi_index = 0;
    for (const auto& item : rows) {
        if (item.kind == RowKind::Electro) phase_phi[item.phase] = phi_index++;
    }
    for (const auto& row : rows) {
        double value = 0.0;
        if (row.kind == RowKind::Chemical) {
            const auto& derivative = *derivatives[row.phase];
            const std::size_t species = row.index;
            double constraint = 0.0;
            for (std::size_t balance_row = 0; balance_row < rank; ++balance_row) {
                constraint += lambda[balance_row] * B[balance_row][species];
            }
            const double electro = z[species] * (
                phase_phi[row.phase].has_value() ? phi[*phase_phi[row.phase]] : 0.0);
            value = std::log(std::max(x[row.phase][species] / volume[row.phase],
                                      std::numeric_limits<double>::min()))
                + derivative.gradient[2 + species] + standard.value[species]
                - constraint - electro;
        } else if (row.kind == RowKind::Pressure) {
            const auto& derivative = *derivatives[row.phase];
            value = volume[row.phase] * (
                1.0 / volume[row.phase] - derivative.gradient[1] - pressure / rt);
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
        const auto& derivative = *derivatives[phase];
        double ideal = ideal_mechanical_helmholtz(x[phase], volume[phase], 1.0);
        result.objective += amount[phase] * (
            derivative.value + ideal + pressure * volume[phase] / rt
            + std::inner_product(standard.value.begin(), standard.value.end(),
                                 x[phase].begin(), 0.0)
        );
    }
    return result;
}

std::vector<double> CompiledSystem::jacobian(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_jacobian");
    const auto base = evaluate(point);
    const std::size_t n = point.size();
    const std::size_t components = model->parameters().component_count();
    std::vector<std::vector<double>> compositions;
    std::vector<double> amount;
    std::vector<double> volumes;
    std::vector<PackingTerms> packings;
    std::vector<std::size_t> offsets;
    std::size_t offset = 0;
    for (const auto& phase : problem.phases) {
        offsets.push_back(offset);
        std::vector<double> composition(components, 0.0);
        for (std::size_t local = 0; local < phase.support.size(); ++local)
            composition[phase.support[local]] = std::exp(point[offset + local]);
        offset += phase.support.size();
        amount.push_back(phase.amount.kind == AmountKind::Free
            ? point[offset++] : phase_amount(phase, n_feed));
        const auto packing = packing_terms(*model,
            problem.temperature.free ? point[primal_size - free_intensive_count]
                                     : problem.temperature.value,
            composition);
        volumes.push_back(packing.coefficient / std::exp(point[offset++]));
        packings.push_back(packing);
        compositions.push_back(std::move(composition));
    }
    const std::optional<std::size_t> temperature = problem.temperature.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const std::optional<std::size_t> pressure = problem.pressure.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const double T = temperature.has_value() ? point[*temperature] : problem.temperature.value;
    const double P = pressure.has_value() ? point[*pressure] : problem.pressure.value;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto standard = standard_term(
        T, components, problem.reactions, problem.standard_gibbs_j_per_mol
    );
    std::vector<HelmholtzDerivatives> derivatives;
    derivatives.reserve(problem.phases.size());
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase)
        derivatives.push_back(model->evaluate(T, volumes[phase], compositions[phase], DerivativeOrder::Hessian));
    const auto z = species_z(model->parameters());
    std::vector<std::optional<std::size_t>> phi(problem.phases.size());
    std::size_t phi_offset = 0;
    for (const auto& row : rows)
        if (row.kind == RowKind::Electro) phi[row.phase] = rank + phi_offset++;
    std::vector<double> result(base.values.size() * n, 0.0);
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
            for (std::size_t local = 0; local < phase.support.size(); ++local) {
                const std::size_t component = phase.support[local];
                const double physical_x = derivative.hessian[
                    (2 + species) * width + (2 + component)]
                    + (component == species ? 1.0 / compositions[row.phase][species] : 0.0);
                const double physical_v = derivative.hessian[(2 + species) * width + 1]
                    - 1.0 / volumes[row.phase];
                const double dv = packings[row.phase].component[component]
                    / std::exp(point[offsets[row.phase] + phase.support.size()
                        + (phase.amount.kind == AmountKind::Free ? 1 : 0)]);
                at(row_index, offsets[row.phase] + local) =
                    physical_x * compositions[row.phase][component] + physical_v * dv;
            }
            const std::size_t volume = offsets[row.phase] + phase.support.size()
                + (phase.amount.kind == AmountKind::Free ? 1 : 0);
            const double physical_v = derivative.hessian[(2 + species) * width + 1]
                - 1.0 / volumes[row.phase];
            at(row_index, volume) = -volumes[row.phase] * physical_v;
            if (temperature.has_value()) at(row_index, *temperature) =
                derivative.hessian[(2 + species) * width]
                + physical_v * packings[row.phase].temperature_derivative
                    / std::exp(point[volume])
                + standard.temperature_derivative[species];
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
                const double dv = packings[row.phase].component[component]
                    / std::exp(point[offsets[row.phase] + phase.support.size()
                        + (phase.amount.kind == AmountKind::Free ? 1 : 0)]);
                at(row_index, offsets[row.phase] + local) =
                    -volumes[row.phase] * derivative.hessian[width + 2 + component]
                        * compositions[row.phase][component] + physical_v * dv;
            }
            const std::size_t volume = offsets[row.phase] + phase.support.size()
                + (phase.amount.kind == AmountKind::Free ? 1 : 0);
            at(row_index, volume) = -volumes[row.phase] * physical_v;
            if (temperature.has_value()) at(row_index, *temperature) = -derivative.hessian[width]
                * volumes[row.phase] + P * volumes[row.phase]
                    / (GAS_CONSTANT_J_PER_MOL_K * T * T)
                + physical_v * packings[row.phase].temperature_derivative
                    / std::exp(point[volume]);
            if (pressure.has_value()) at(row_index, *pressure) = -volumes[row.phase] / rt;
        } else if (row.kind == RowKind::Moiety) {
            for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
                for (std::size_t local = 0; local < problem.phases[phase].support.size(); ++local) {
                    const std::size_t species = problem.phases[phase].support[local];
                    at(row_index, offsets[phase] + local) = amount[phase] * B[row.index][species]
                        * compositions[phase][species] / feed_total;
                }
                if (problem.phases[phase].amount.kind == AmountKind::Free) {
                    const std::size_t column = offsets[phase] + problem.phases[phase].support.size();
                    for (const std::size_t species : problem.phases[phase].support)
                        at(row_index, column) += B[row.index][species]
                            * compositions[phase][species] / feed_total;
                }
            }
        } else if (row.kind == RowKind::Electro) {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, offsets[row.phase] + local) =
                    z[problem.phases[row.phase].support[local]]
                    * compositions[row.phase][problem.phases[row.phase].support[local]];
        } else {
            for (std::size_t local = 0; local < problem.phases[row.phase].support.size(); ++local)
                at(row_index, offsets[row.phase] + local) =
                    compositions[row.phase][problem.phases[row.phase].support[local]];
        }
    }
    scope.probe("rows", static_cast<double>(base.values.size()));
    scope.probe("columns", static_cast<double>(n));
    return result;
}

std::vector<double> CompiledSystem::lever_jacobian(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_lever_jacobian");
    std::size_t offset = 0;
    std::vector<std::size_t> offsets;
    for (const auto& phase : problem.phases) {
        offsets.push_back(offset);
        offset += phase.support.size();
        if (phase.amount.kind == AmountKind::Free) ++offset;
        ++offset;
    }
    const std::size_t components = model->parameters().component_count();
    const double T = problem.temperature.free ? point[offset++] : problem.temperature.value;
    const double P = problem.pressure.free ? point[offset++] : problem.pressure.value;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto standard = standard_term(
        T, components, problem.reactions, problem.standard_gibbs_j_per_mol
    );
    std::size_t levers = 0;
    std::vector<double> result;
    std::vector<double> pinned_amount_feed_term(rank, 0.0);
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        std::vector<double> composition(components, 0.0);
        for (std::size_t local = 0; local < problem.phases[phase].support.size(); ++local)
            composition[problem.phases[phase].support[local]] =
                std::exp(point[offsets[phase] + local]);
        std::size_t cursor = offsets[phase] + problem.phases[phase].support.size();
        if (problem.phases[phase].amount.kind == AmountKind::Free) ++cursor;
        const auto packing = packing_terms(*model, T, composition);
        const double volume = packing.coefficient / std::exp(point[cursor]);
        const auto derivative = model->evaluate(T, volume, composition, DerivativeOrder::Hessian);
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
                const double standard_term_temperature_derivative =
                    standard.temperature_derivative[rows[row].index];
                result[levers * row] = derivative.hessian[(2 + rows[row].index) * width]
                    + physical_v * packing.temperature_derivative / std::exp(point[cursor])
                    + standard_term_temperature_derivative;
                for (std::size_t active = 2 + components; active < width; ++active) {
                    const double active_parameter_column = derivative.hessian[
                        (2 + rows[row].index) * width + active];
                    // EqID: packing_coefficient_parameter_derivative
                    const double packing_parameter_term = packing.parameter_derivative[
                        active - 2 - components] / std::exp(point[cursor]);
                    result[levers * row + active] = active_parameter_column
                        + physical_v * packing_parameter_term;
                }
            } else if (rows[row].kind == RowKind::Pressure) {
                const double physical_v = -derivative.gradient[1]
                    - volume * derivative.hessian[width + 1] - P / rt;
                result[levers * row] = -volume * derivative.hessian[width]
                    + P * volume / (GAS_CONSTANT_J_PER_MOL_K * T * T)
                    + physical_v * packing.temperature_derivative / std::exp(point[cursor]);
                for (std::size_t active = 2 + components; active < width; ++active) {
                    const double active_parameter_column =
                        -volume * derivative.hessian[width + active];
                    const double packing_parameter_term = packing.parameter_derivative[
                        active - 2 - components] / std::exp(point[cursor]);
                    result[levers * row + active] = active_parameter_column
                        + physical_v * packing_parameter_term;
                }
            }
        }
    }
    for (std::size_t row = 0; row < rows.size(); ++row) {
        if (rows[row].kind == RowKind::Pressure) {
            const auto& phase = problem.phases[rows[row].phase];
            const std::size_t cursor = offsets[rows[row].phase] + phase.support.size()
                + (phase.amount.kind == AmountKind::Free ? 1 : 0);
            std::vector<double> composition(components, 0.0);
            for (std::size_t local = 0; local < phase.support.size(); ++local)
                composition[phase.support[local]] = std::exp(
                    point[offsets[rows[row].phase] + local]);
            const auto packing = packing_terms(*model, T, composition);
            result[levers * row + 1] = -packing.coefficient
                / std::exp(point[cursor]) / rt;
        } else if (rows[row].kind == RowKind::Moiety) {
            for (std::size_t species = 0; species < components; ++species)
                result[levers * row + 2 + species] = -B[rows[row].index][species]
                    / feed_total + pinned_amount_feed_term[rows[row].index];
        }
    }
    return result;
}

std::vector<double> CompiledSystem::objective_gradient(const std::vector<double>& point) const {
    core::Scope scope("equilibrium_objective_gradient");
    std::vector<double> result(point.size(), 0.0);
    if (pinned_count != 0) return result;
    const std::size_t components = model->parameters().component_count();
    std::size_t offset = 0;
    std::vector<std::size_t> offsets;
    for (const auto& phase : problem.phases) {
        offsets.push_back(offset);
        offset += phase.support.size() + (phase.amount.kind == AmountKind::Free ? 1 : 0) + 1;
    }
    const std::optional<std::size_t> temperature = problem.temperature.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const std::optional<std::size_t> pressure = problem.pressure.free
        ? std::optional<std::size_t>(offset++) : std::nullopt;
    const double T = temperature.has_value() ? point[*temperature] : problem.temperature.value;
    const double P = pressure.has_value() ? point[*pressure] : problem.pressure.value;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto standard_at_temperature = standard_term(
        T, components, problem.reactions, problem.standard_gibbs_j_per_mol
    );
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        std::vector<double> x(components, 0.0);
        for (std::size_t local = 0; local < spec.support.size(); ++local)
            x[spec.support[local]] = std::exp(point[offsets[phase] + local]);
        std::size_t cursor = offsets[phase] + spec.support.size();
        const double amount = point[cursor++];
        const auto packing = packing_terms(*model, T, x);
        const double volume = packing.coefficient / std::exp(point[cursor]);
        const auto derivative = model->evaluate(T, volume, x, DerivativeOrder::Hessian);
        double standard = 0.0;
        double standard_term_temperature_derivative = 0.0;
        double ideal = ideal_mechanical_helmholtz(x, volume, 1.0, false);
        for (std::size_t species = 0; species < components; ++species) {
            standard += standard_at_temperature.value[species] * x[species];
            standard_term_temperature_derivative +=
                standard_at_temperature.temperature_derivative[species] * x[species];
        }
        const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        const double bracket = derivative.value + ideal + P * volume / rt + standard;
        for (std::size_t local = 0; local < spec.support.size(); ++local) {
            const std::size_t species = spec.support[local];
            const double chemical =
                derivative.gradient[2 + species] + std::log(std::max(
                    x[species] / volume, std::numeric_limits<double>::min()))
                + standard_at_temperature.value[species];
            const double physical_v = derivative.gradient[1] - sum_x / volume + P / rt;
            result[offsets[phase] + local] = amount * (
                chemical * x[species] + physical_v * packing.component[species]
                    / std::exp(point[cursor]));
        }
        result[cursor] = -amount * volume
            * (derivative.gradient[1] - sum_x / volume + P / rt);
        if (temperature.has_value()) result[*temperature] += amount * (
            derivative.gradient[0] - P * volume
                / (GAS_CONSTANT_J_PER_MOL_K * T * T)
            + standard_term_temperature_derivative);
        if (pressure.has_value()) result[*pressure] += amount * volume / rt;
        result[offsets[phase] + spec.support.size()] = bracket;
    }
    return result;
}

std::vector<double> CompiledSystem::hessian(
    const std::vector<double>& point,
    const std::vector<double>& constraint_multipliers
) const {
    core::Scope scope("equilibrium_hessian");
    const std::size_t n = point.size();
    std::vector<double> result(n * n, 0.0);
    if (pinned_count != 0) return result;
    const std::size_t components = model->parameters().component_count();
    std::size_t offset = 0;
    std::vector<std::size_t> offsets;
    for (const auto& phase : problem.phases) {
        offsets.push_back(offset);
        offset += phase.support.size() + 1 + (phase.amount.kind == AmountKind::Free ? 1 : 0);
    }
    const double T = problem.temperature.value;
    const double P = problem.pressure.value;
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const auto standard_at_temperature = standard_term(
        T, components, problem.reactions, problem.standard_gibbs_j_per_mol
    );
    auto put = [&](std::size_t row, std::size_t column, double value) {
        result[row * n + column] = value;
        result[column * n + row] = value;
    };
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        const auto& spec = problem.phases[phase];
        std::vector<double> x(components, 0.0);
        for (std::size_t local = 0; local < spec.support.size(); ++local)
            x[spec.support[local]] = std::exp(point[offsets[phase] + local]);
        std::size_t cursor = offsets[phase] + spec.support.size();
        const std::optional<std::size_t> amount = spec.amount.kind == AmountKind::Free
            ? std::optional<std::size_t>(cursor++) : std::nullopt;
        const std::size_t volume = cursor;
        const double phase_amount_value = amount.has_value() ? point[*amount]
            : phase_amount(spec, n_feed);
        const auto packing = packing_terms(*model, T, x);
        const double eta = std::exp(point[volume]);
        const double physical_volume = packing.coefficient / eta;
        const auto derivative = model->evaluate(T, physical_volume, x, DerivativeOrder::Hessian);
        const std::size_t width = derivative.gradient.size();
        const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
        const double physical_v_gradient = phase_amount_value * (
            derivative.gradient[1] - sum_x / physical_volume + P / rt);
        const double physical_vv = phase_amount_value * (
            derivative.hessian[width + 1] + sum_x / (physical_volume * physical_volume));
        double standard = 0.0;
        for (std::size_t species = 0; species < components; ++species)
            standard += standard_at_temperature.value[species] * x[species];
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const std::size_t si = spec.support[i];
            const double dv_i = packing.component[si] / eta;
            const double physical_x_gradient = phase_amount_value * (
                derivative.gradient[2 + si]
                + std::log(std::max(x[si] / physical_volume,
                                    std::numeric_limits<double>::min()))
                + standard_at_temperature.value[si]);
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
                put(offsets[phase] + i, offsets[phase] + j,
                    x[si] * x[sj] * h_xx + x[si] * h_xv_i * dv_j
                    + dv_i * h_xv_j * x[sj] + dv_i * physical_vv * dv_j
                    + (si == sj ? physical_x_gradient * x[si]
                                      + physical_v_gradient * dv_i : 0.0));
            }
            const double h_xv = phase_amount_value * (
                derivative.hessian[(2 + si) * width + 1] - 1.0 / physical_volume);
            put(offsets[phase] + i, volume,
                -physical_volume * (x[si] * h_xv + dv_i * physical_vv)
                - physical_v_gradient * dv_i);
        }
        put(volume, volume, physical_volume * physical_volume * physical_vv
            + physical_v_gradient * physical_volume);
        if (amount.has_value()) {
            for (std::size_t i = 0; i < spec.support.size(); ++i) {
                const std::size_t species = spec.support[i];
                const double mu = derivative.gradient[2 + species]
                    + std::log(std::max(x[species] / physical_volume,
                                        std::numeric_limits<double>::min()))
                    + standard_at_temperature.value[species];
                const double pressure_gradient = derivative.gradient[1]
                    - sum_x / physical_volume + P / rt;
                put(*amount, offsets[phase] + i, mu * x[species]
                    + pressure_gradient * packing.component[species] / eta);
            }
            put(*amount, volume, -physical_volume * (
                derivative.gradient[1] - sum_x / physical_volume + P / rt));
        }
    }
    if (!constraint_multipliers.empty()) {
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
                const std::size_t amount = offsets[phase] + spec.support.size();
                for (std::size_t local = 0; local < spec.support.size(); ++local) {
                    const std::size_t species = spec.support[local];
                    const std::size_t coordinate = offsets[phase] + local;
                    const double x = std::exp(point[coordinate]);
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
    }
    scope.probe("entries", static_cast<double>(result.size()));
    return result;
}

}  // namespace epcsaft::equilibrium
