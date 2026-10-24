#include "properties/helmholtz_model.hpp"

#include "core/trace.hpp"
#include "contributions/association.hpp"
#include "contributions/born.hpp"
#include "contributions/electrolyte.hpp"
#include "contributions/helmholtz.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <limits>
#include <numeric>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace epcsaft {

HelmholtzModel::~HelmholtzModel() = default;
namespace {

bool finite_vector(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

bool positive_vector(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value) && value > 0.0;
    });
}

template <typename Value>
bool symmetric(const std::vector<Value>& values, std::size_t size) {
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = row + 1; column < size; ++column) {
            if (values[row * size + column] != values[column * size + row]) return false;
        }
    }
    return true;
}

}  // namespace

Eigen::FullPivLU<Eigen::MatrixXd> properties::site_factorization(
    const std::vector<double>& sites, const std::vector<double>& composition, const Parameters& parameters,
    double number_density, const std::vector<double>& delta, PropertyActionStatus* status
) {
    Eigen::FullPivLU<Eigen::MatrixXd> factorization;
    site_factorization(factorization, sites, composition, parameters, number_density, delta, status);
    return factorization;
}

void properties::site_factorization(
    Eigen::FullPivLU<Eigen::MatrixXd>& factorization,
    const std::vector<double>& sites,
    const std::vector<double>& composition,
    const Parameters& parameters,
    double number_density,
    const std::vector<double>& delta,
    PropertyActionStatus* status
) {
    const std::size_t sites_count = sites.size();
    if (sites_count == 0) { factorization = {}; return; }  // no association sites: no mass-action system to factor
    factorization.compute(mass_action_jacobian(composition, parameters, sites,
        number_density, delta, bonded_densities(composition, parameters, sites, delta)));
    factorization.setThreshold(policy::ASSOCIATION_JACOBIAN_RANK.value * sites_count);
    auto scope = core::Scope::detailed("association_derivative");
    const auto failure = [&](const char* reason) {
        if (!scope) return;
        scope->observe([&] {
            const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
            scope->probe("association_derivative_failure", {static_cast<double>(sites_count),
                static_cast<double>(factorization.rank()), pivots.maxCoeff(), pivots.minCoeff(),
                pivots.minCoeff() == 0.0 ? std::numeric_limits<double>::infinity()
                    : pivots.maxCoeff() / pivots.minCoeff()}, 0, reason);
        });
    };
    if (factorization.rank() < static_cast<Eigen::Index>(sites_count)) {
        failure("association Jacobian is singular");
        if (status != nullptr) {
            *status = PropertyActionStatus::AssociationSingular;
            return;
        }
        throw std::domain_error("association Jacobian is singular");
    }
    if (!admitted_pivot_ratio(factorization)) {
        failure("association Jacobian is ill-conditioned");
        if (status != nullptr) {
            *status = PropertyActionStatus::AssociationIllConditioned;
            return;
        }
        throw std::runtime_error("association Jacobian is ill-conditioned");
    }
}

namespace {

// One central association factorization with its rank and pivot diagnostics; a failure
// is reported through `status` and the factorization is still returned.
std::unique_ptr<Eigen::FullPivLU<Eigen::MatrixXd>> association_factorization(
    const Parameters& parameters, double temperature, double density, const std::vector<double>& sites,
    const std::vector<double>& composition, PropertyActionStatus& status, PropertyActionDiagnostics& diagnostics
) {
    if (sites.empty()) return nullptr;
    const auto sigma = segment_sigma(temperature, parameters), diameter = segment_diameter(temperature, parameters, sigma);
    auto result = std::make_unique<Eigen::FullPivLU<Eigen::MatrixXd>>(properties::site_factorization(
        sites, composition, parameters, density * NUMBER_DENSITY_PER_ANGSTROM3,
        association_delta(temperature, density, composition, parameters, sigma, diameter), &status));
    diagnostics.site_factorizations = 1;
    diagnostics.site_rank = static_cast<std::size_t>(result->rank());
    const auto pivots = result->matrixLU().diagonal().cwiseAbs();
    if (pivots.size() != 0 && pivots.minCoeff() > 0.0) diagnostics.site_pivot_ratio = pivots.maxCoeff() / pivots.minCoeff();
    return result;
}

void require(bool condition, const char* field) {
    if (!condition) {
        throw std::invalid_argument(std::string("invalid ") + field);
    }
}

// Each law slot in double: tape replay skips the domain refusal (thermal-rule volumes: the double site solve).
void require_law_domain(double T, const Parameters& parameters) {
    const auto laws = temperature_law_values(T, parameters);
    const std::size_t sites = parameters.site_component.size(), components = parameters.segment_count.size();
    for (const auto& [kind, slot, law] : parameters.temperature_law_slot)  // kinds 3, 4 own site-pair slots
        temperature_law_value(T, parameters, laws, TemperatureLaw(kind), slot, kind > 2 ? sites : components);
}

using Jet = std::vector<double>;

Jet jet_add(const Jet& left, const Jet& right, double sign = 1.0) {
    Jet result(left.size());
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = left[i] + sign * right[i];
    return result;
}
Jet jet_scale(const Jet& value, double scale) {
    Jet result = value;
    for (double& item : result) item *= scale;
    return result;
}
Jet jet_multiply(const Jet& left, const Jet& right) {
    Jet result(left.size(), 0.0);
    for (std::size_t mask = 0; mask < result.size(); ++mask) {
        for (std::size_t subset = mask;; subset = (subset - 1) & mask) {
            result[mask] += left[subset] * right[mask ^ subset];
            if (subset == 0) break;
        }
    }
    return result;
}
Jet jet_inverse(const Jet& value) {
    if (value.empty() || !std::isfinite(value.front()) || value.front() == 0.0)
        throw std::domain_error("nonfinite jet inverse");
    Jet result(value.size(), 0.0);
    result.front() = 1.0 / value.front();
    for (std::size_t mask = 1; mask < value.size(); ++mask) {
        double known = 0.0;
        for (std::size_t subset = mask; subset != 0; subset = (subset - 1) & mask)
            known += value[subset] * result[mask ^ subset];
        result[mask] = -known * result.front();
    }
    return result;
}
// log(v0 (1 + u)), v0 > 0, with nilpotent u: the series ends at the jet's direction count.
Jet jet_log(const Jet& value) {
    Jet u = jet_scale(value, 1.0 / value.front()), power, result(value.size(), 0.0);
    u.front() = 0.0;
    power = u;
    for (std::size_t k = 1; (std::size_t(1) << (k - 1)) < value.size(); ++k) {
        result = jet_add(result, power, (k % 2 ? 1.0 : -1.0) / static_cast<double>(k));
        power = jet_multiply(power, u);
    }
    result.front() = std::log(value.front());
    return result;
}

struct ActionDirection final {
    double temperature = 0.0, molar_density = 0.0;
    std::vector<double> composition, parameter;
};

std::vector<std::vector<double>> action_coefficients(
    const std::vector<double>& point, double molar_density,
    const std::vector<ActionDirection>& directions, std::size_t extra_temperature,
    std::size_t extra_volume, std::size_t order
) {
    std::vector<std::vector<double>> result(order + 1, std::vector<double>(point.size()));
    result.front() = point;
    const std::size_t components = directions.empty()
        ? point.size() - 2 : directions.front().composition.size();
    const std::size_t parameters = directions.empty()
        ? point.size() - 2 - components : directions.front().parameter.size();
    double temperature = 0.0, density = 0.0;
    std::vector<double> composition(components), parameter(parameters);
    for (const auto& direction : directions) {
        temperature += direction.temperature;
        density += direction.molar_density;
        for (std::size_t i = 0; i < components; ++i) composition[i] += direction.composition[i];
        for (std::size_t i = 0; i < parameters; ++i) parameter[i] += direction.parameter[i];
    }
    if (order != 0) {
        result[1][0] = temperature + static_cast<double>(extra_temperature) * point[0];
        result[1][1] = -density / (molar_density * molar_density)
            + static_cast<double>(extra_volume) / molar_density;
        for (std::size_t i = 0; i < components; ++i) result[1][2 + i] = composition[i];
        for (std::size_t i = 0; i < parameters; ++i)
            result[1][2 + components + i] = parameter[i];
    }
    for (std::size_t k = 2; k <= order; ++k)
        result[k][1] = (k & 1 ? -1.0 : 1.0)
            * std::pow(density, static_cast<int>(k))
            / std::pow(molar_density, static_cast<int>(k + 1));
    return result;
}

struct ActionContext final {
    detail::HelmholtzTape& tape;
    const std::vector<double>& point;
    double molar_density;
    const std::vector<ActionDirection>& directions;
    const Eigen::FullPivLU<Eigen::MatrixXd>* factorization;
    std::vector<std::size_t> unit_ids = {};  // e_i amount directions at fixed T and V, for mu_i
    struct CachedPath final {
        std::size_t order = 0;
        std::vector<std::vector<double>> outputs;
    };
    std::map<std::tuple<std::size_t, std::size_t, std::size_t>, CachedPath> paths;
    std::size_t forward_sweeps = 0, solves = 0;
    std::size_t maximum_order = 0;
    bool planning = false;
    double coefficient_residual = 0.0;
    std::vector<double> planned_output = std::vector<double>(tape.function.Range());

    void run() {
        for (auto& [key, entry] : paths) {
            if (!entry.outputs.empty()) continue;
            const auto [mask, temperature, volume] = key;
            std::vector<ActionDirection> selected;
            for (std::size_t i = 0; i < directions.size(); ++i)
                if (mask & (std::size_t(1) << i)) selected.push_back(directions[i]);
            tape.implicit_taylor(
                action_coefficients(point, molar_density,
                    selected, temperature, volume, entry.order),
                entry.order, factorization, &solves, &forward_sweeps,
                &coefficient_residual, &entry.outputs);
        }
        planning = false;
    }

    const std::vector<double>& path(
        std::size_t mask, std::size_t extra_temperature, std::size_t extra_volume,
        std::size_t required_order
    ) {
        const auto key = std::make_tuple(mask, extra_temperature, extra_volume);
        auto [entry, inserted] = paths.try_emplace(key);
        if (planning) {
            const auto order = mask == 0 && extra_temperature + extra_volume != 0
                ? std::max(required_order, maximum_order) : required_order;
            if (entry->second.order < order && !entry->second.outputs.empty())
                throw std::logic_error("action path order was not planned");
            entry->second.order = std::max(entry->second.order, order);
            return planned_output;
        }
        return entry->second.outputs.at(required_order);
    }

    Jet extract(
        const std::vector<std::size_t>& ids, std::size_t extra_temperature,
        std::size_t extra_volume, std::size_t output
    ) {
        const std::size_t combinations = std::size_t(1) << ids.size();
        const double temperature_scale = extra_temperature == 0 ? 1.0
            : extra_temperature == 1 ? point[0] : point[0] * point[0];
        Jet result(combinations, 0.0);
        for (std::size_t target = 0; target < combinations; ++target) {
            std::vector<std::size_t> selected_ids;
            for (std::size_t i = 0; i < ids.size(); ++i)
                if (target & (std::size_t(1) << i)) selected_ids.push_back(ids[i]);
            const std::size_t labels = selected_ids.size()
                + extra_temperature + extra_volume;
            double value = 0.0;
            const std::size_t subsets = std::size_t(1) << labels;
            for (std::size_t selected = 1; selected < subsets; ++selected) {
                std::size_t mask = 0, temperature = 0, volume = 0;
                for (std::size_t i = 0; i < selected_ids.size(); ++i)
                    if (selected & (std::size_t(1) << i))
                        mask |= std::size_t(1) << selected_ids[i];
                for (std::size_t i = selected_ids.size(); i < labels; ++i) {
                    if (!(selected & (std::size_t(1) << i))) continue;
                    if (i < selected_ids.size() + extra_temperature) ++temperature;
                    else ++volume;
                }
                const double sign = std::bitset<std::numeric_limits<std::size_t>::digits>(
                    selected).count() % 2
                    == labels % 2 ? 1.0 : -1.0;
                value += sign * path(mask, temperature, volume, labels)[output];
            }
            if (labels == 0) value = path(0, 0, 0, 0)[output];
            result[target] = (extra_volume
                ? value * std::pow(molar_density, static_cast<int>(extra_volume)) : value)
                / temperature_scale;
        }
        return result;
    }
};

Jet temperature_jet(const ActionContext& context, const std::vector<std::size_t>& ids) {
    Jet result(std::size_t(1) << ids.size(), 0.0);
    result.front() = context.point[0];
    for (std::size_t i = 0; i < ids.size(); ++i)
        result[std::size_t(1) << i] = context.directions[ids[i]].temperature;
    return result;
}

Jet volume_jet(const ActionContext& context, const std::vector<std::size_t>& ids) {
    Jet result(std::size_t(1) << ids.size(), 0.0);
    result.front() = context.point[1];
    for (std::size_t mask = 1; mask < result.size(); ++mask) {
        const std::size_t order =
            std::bitset<std::numeric_limits<std::size_t>::digits>(mask).count();
        double density_product = 1.0;
        for (std::size_t i = 0; i < ids.size(); ++i)
            if (mask & (std::size_t(1) << i))
                density_product *= context.directions[ids[i]].molar_density;
        result[mask] = (order & 1 ? -1.0 : 1.0)
            * std::tgamma(static_cast<double>(order + 1)) * density_product
            / std::pow(context.molar_density, static_cast<int>(order + 1));
    }
    return result;
}

struct PropertyJets final {
    std::array<Jet, 5> ideal, residual;
    std::vector<Jet> chemical_potential;  // residual mu_i / RT; ln Z is added at fixed P
    Jet pressure, volume_temperature, volume_pressure;
    double pressure_slope = 0.0;
};

PropertyJets property_jets(
    ActionContext& context, const std::vector<std::size_t>& ids,
    const std::array<bool, 5>& need_ideal,
    const std::array<bool, 5>& need_residual, bool need_pressure, bool need_fugacity = false
) {
    PropertyJets result;
    result.chemical_potential.resize(context.unit_ids.size());
    for (std::size_t i = 0; i < 5; ++i) {
        if (need_ideal[i] && context.tape.ideal_available)
            result.ideal[i] = context.extract(ids, 0, 0, context.tape.ideal_offset + i);
    }
    const bool h = need_residual[0], u = need_residual[1], s = need_residual[2];
    const bool cp = need_residual[3], cv = need_residual[4];
    if (!std::any_of(need_residual.begin(), need_residual.end(), [](bool value) { return value; })
        && !need_pressure && !need_fugacity) return result;
    const double gas_constant = GAS_CONSTANT_J_PER_MOL_K;
    const auto temperature = temperature_jet(context, ids);
    const auto volume = volume_jet(context, ids);
    const auto temperature_square = jet_multiply(temperature, temperature);
    Jet phi, phi_t, phi_v, phi_tt, phi_tv, phi_vv, inverse_volume;
    if (s) phi = context.extract(ids, 0, 0, 0);
    if (h || u || s || cp || cv)
        phi_t = context.extract(ids, 1, 0, 0);
    if (h || cp || need_pressure) phi_v = context.extract(ids, 0, 1, 0);
    if (cv || cp) phi_tt = context.extract(ids, 2, 0, 0);
    if (cp) {
        phi_tv = context.extract(ids, 1, 1, 0);
        phi_vv = context.extract(ids, 0, 2, 0);
    }
    if (cp || need_pressure) inverse_volume = jet_inverse(volume);
    Jet internal;
    if (h || u) {
        internal = jet_scale(jet_multiply(temperature_square, phi_t), -gas_constant);
    }
    if (u) result.residual[1] = internal;
    if (h) result.residual[0] = jet_add(internal, jet_scale(
        jet_multiply(jet_multiply(temperature, volume), phi_v), -gas_constant));
    if (s) result.residual[2] = jet_scale(jet_add(
        phi, jet_multiply(temperature, phi_t)), -gas_constant);
    if (cv || cp) result.residual[4] = jet_scale(jet_add(
        jet_scale(jet_multiply(temperature, phi_t), 2.0),
        jet_multiply(temperature_square, phi_tt)), -gas_constant);
    if (cp) {
        const Jet pressure_volume = jet_scale(jet_multiply(temperature,
            jet_add(jet_multiply(inverse_volume, inverse_volume), phi_vv)), -gas_constant);
        result.pressure_slope = -pressure_volume.front() * context.point[1] * context.point[1];
        const Jet pressure_temperature = jet_scale(jet_add(jet_add(
            inverse_volume, phi_v, -1.0),
            jet_multiply(temperature, phi_tv), -1.0), gas_constant);
        result.volume_pressure = jet_inverse(pressure_volume);
        result.volume_temperature =
            jet_scale(jet_multiply(pressure_temperature, result.volume_pressure), -1.0);
        result.residual[3] = jet_add(result.residual[4], jet_scale(
            jet_multiply(temperature, jet_multiply(jet_multiply(
                pressure_temperature, pressure_temperature), jet_inverse(pressure_volume))),
            -1.0));
        result.residual[3].front() -= gas_constant;
    }
    if (need_pressure) result.pressure = jet_scale(jet_multiply(temperature,
        jet_add(inverse_volume, phi_v, -1.0)), gas_constant);
    // EqID: log_fugacity_coefficient -- the d(A_res/RT)/dn_i term at T, V, n_j.
    auto unit = ids;
    unit.push_back(0);
    for (std::size_t k = 0; need_fugacity && k < context.unit_ids.size(); ++k) {
        unit.back() = context.unit_ids[k];
        const Jet mu = context.extract(unit, 0, 0, 0);
        result.chemical_potential[k].assign(mu.begin() + static_cast<std::ptrdiff_t>(temperature.size()), mu.end());
    }
    return result;
}

bool is_ideal(PropertyObservable observable) { return static_cast<std::size_t>(observable) < 5; }
bool is_fugacity(PropertyObservable observable) { return observable == PropertyObservable::LogFugacityCoefficient; }
bool is_volume_partial(PropertyObservable observable) {
    return static_cast<std::size_t>(observable) >= 15 && !is_fugacity(observable);
}
bool is_total(PropertyObservable observable) {
    return static_cast<std::size_t>(observable) >= 10 && !is_volume_partial(observable) && !is_fugacity(observable);
}
std::size_t property_index(PropertyObservable observable) {
    return is_volume_partial(observable) ? 3 : is_fugacity(observable) ? 0 : static_cast<std::size_t>(observable) % 5;
}
bool is_heat_capacity(PropertyObservable observable) {
    const auto index = property_index(observable);
    return !is_fugacity(observable) && (index == 3 || index == 4);
}
bool zero_direction(const ActionDirection& direction) {
    return direction.temperature == 0.0 && direction.molar_density == 0.0
        && std::all_of(direction.composition.begin(), direction.composition.end(),
            [](double value) { return value == 0.0; })
        && std::all_of(direction.parameter.begin(), direction.parameter.end(),
            [](double value) { return value == 0.0; });
}
// The tape's independent point (T, V, n, theta).
std::vector<double> action_point(
    const HelmholtzModel& model, double temperature, double volume,
    const std::vector<double>& amounts
) {
    std::vector<double> point = {temperature, volume};
    point.insert(point.end(), amounts.begin(), amounts.end());
    for (const double value : model.active_parameter_values()) point.push_back(value);
    return point;
}

void validate_parameter_coordinate(
    const HelmholtzModel& model, ParameterFamily family,
    const std::vector<std::size_t>& indices, bool active
) {
    const auto& parameters = model.parameters();
    const auto& values = detail::family_values(parameters, family);
    require(!indices.empty(), "parameter indices");
    for (const auto index : indices) require(index < values.second, "parameter index");
    const bool association = family == ParameterFamily::AssociationEnergyOverK
        || family == ParameterFamily::AssociationVolume;
    const bool slope = family == ParameterFamily::KijReciprocalTemperatureSlope;
    const bool pair = association || slope || family == ParameterFamily::Kij || family == ParameterFamily::Lij;
    if (!pair) {
        require(indices.size() == 1, "scalar parameter coordinate");
        if (family == ParameterFamily::BornDiameter || family == ParameterFamily::DebyeHuckelDiameter
            || family == ParameterFamily::PackingDiameter || family == ParameterFamily::IonSuppressionCoefficient)
            require(parameters.charge_number[indices.front()] != 0, "electrolyte coordinate needs charged slot");
        if (family == ParameterFamily::DielectricVolume || family == ParameterFamily::DielectricTemperature)
            require(parameters.charge_number[indices.front()] == 0 && parameters.permittivity_form.is(PermittivityForm::DensityTemperature),
                "dielectric coordinate needs a neutral density-temperature slot");
        return;
    }
    const auto size = association ? parameters.site_count() : parameters.component_count();
    const auto first = indices.front(), row = first / size, column = first % size;
    if (row == column) {
        require(indices.size() == 1 && (association || !active), "diagonal parameter coordinate");
    } else {
        require(indices.size() == 2 && first < indices[1]
            && indices[1] == column * size + row, "reciprocal parameter coordinate");
    }
    // A law owns its slot and a thermal combining rule its edge, so their stored constants are unconsumed.
    const auto owner = family == ParameterFamily::Lij ? TemperatureLaw::Lij
        : family == ParameterFamily::AssociationEnergyOverK ? TemperatureLaw::AssociationEnergyOverK
        : family == ParameterFamily::AssociationVolume ? TemperatureLaw::AssociationVolume : TemperatureLaw::Kij;
    for (const auto index : indices) {
        require(!active || !association || parameters.association_active[index] != 0, "inactive association coordinate");
        require(!active || !slope || parameters.k_ij_reference_temperature_k[index] > 0.0, "k_ij slope without reference temperature");
        for (const auto& [target, slot, law] : parameters.temperature_law_slot)
            require(!active || target != static_cast<int>(owner) || static_cast<std::size_t>(slot) != index, "unconsumed law-owned coordinate");
        for (const auto& rule : parameters.association_thermal_rule)
            require(!active || !association || static_cast<std::size_t>(rule[0]) != index, "unconsumed law-owned coordinate");
    }
}

void validate_action_request(
    const HelmholtzModel& model,
    const std::vector<PropertyObservable>& observables,
    const std::vector<DensityActionDirection>& directions
) {
    require(directions.size() >= 1 && directions.size() <= 3, "action order");
    std::vector<bool> seen(18, false);
    for (const auto observable : observables) {
        const auto index = static_cast<std::size_t>(observable);
        require(index < seen.size() && !seen[index], "observables");
        seen[index] = true;
    }
    const auto components = model.parameters().component_count();
    for (const auto& direction : directions) {
        require(direction.composition.size() == components,
            "direction composition");
        require(std::isfinite(direction.temperature)
            && std::isfinite(direction.molar_density), "direction");
        require(finite_vector(direction.composition), "direction composition");
        double balance = 0.0, charge = 0.0;
        for (std::size_t i = 0; i < components; ++i) {
            balance += direction.composition[i];
            charge += direction.composition[i] * model.parameters().charge_number[i];
        }
        require(std::abs(balance) <= policy::COMPOSITION_SUM_TOLERANCE.value,
            "direction composition balance");
        require(std::abs(charge) <= policy::COMPOSITION_SUM_TOLERANCE.value,
            "direction charge balance");
        for (const auto& parameter : direction.parameter) {
            require(std::isfinite(parameter.value), "direction parameter");
            validate_parameter_coordinate(model, parameter.family, parameter.indices, false);
        }
    }
}

std::vector<double> project_parameters(const HelmholtzModel& model, const std::vector<ParameterActionDirection>& requested) {
    std::vector<double> result(model.active_parameters().size(), 0.0);
    for (const auto& item : requested) {
        for (std::size_t coordinate = 0; coordinate < model.active_parameters().size(); ++coordinate) {
            const auto& active = model.active_parameters()[coordinate];
            if (active.family == item.family && active.indices == item.indices) {
                result[coordinate] += item.value;
                break;
            }
        }
    }
    return result;
}

struct ActionPlan final {
    std::array<bool, 5> ideal{}, residual{};
    bool pressure = false, fugacity = false;
};

std::vector<IdealCorrelation> selected_ideal_correlations(
    const std::optional<ThermochemistryRecord>& record, double temperature,
    bool& endpoint
) {
    endpoint = false;
    if (!record.has_value()) return {};
    std::vector<IdealCorrelation> result;
    result.reserve(record->component_intervals.size());
    for (const auto& intervals : record->component_intervals) {
        const auto* selected = ideal_interval(intervals, temperature);
        if (selected == nullptr) return {};
        endpoint = endpoint || temperature == selected->lower || temperature == selected->upper;
        result.push_back(selected->correlation);
    }
    return result;
}

// {ideal and total entropy, every other result}: a neutral insertion follows its infinite-dilution limit (#73)
// except through x ln x and in a model with a multipole term, whose A2-only branch with every moment absent
// omits A3 beyond second composition order (full EOS 4.1); an ionic insertion is singular throughout.
std::array<PropertyActionStatus, 2> support_status(
    const HelmholtzModel& model, const std::vector<double>& composition,
    const std::vector<ActionDirection>& directions
) {
    bool support = false, ionic = false;
    for (const auto& direction : directions) {
        for (std::size_t i = 0; i < composition.size(); ++i) {
            if (composition[i] != 0.0 || direction.composition[i] == 0.0) continue;
            (model.parameters().charge_number[i] == 0 ? support : ionic) = true;
        }
    }
    const auto status = [&](bool singular) { return ionic ? PropertyActionStatus::IonicInsertion
        : singular ? PropertyActionStatus::SupportChange : PropertyActionStatus::Available; };
    const auto& parameters = model.parameters();
    const bool multipole = !parameters.dipole_count.empty() || !parameters.quadrupole_count.empty();
    return {status(support), status(support && multipole)};
}

bool ionic_insertion_value(const Parameters& p, const std::vector<double>& x, std::size_t component, bool charged) {
    if (!charged && p.charge_number[component] != 0) return true;
    if (p.cation_mixing.empty()) return false;
    const auto [s, i, j, k] = p.cation_mixing.front();
    return x[s] == 0.0 && ((component == i && x[i] == 0.0 && x[j] > 0.0)
        || (component == j && x[j] == 0.0 && x[i] > 0.0));
}

PropertyActionBatch unavailable_batch(
    const std::vector<PropertyObservable>& observables, PropertyActionStatus status,
    std::size_t order, std::size_t components
) {
    PropertyActionBatch result;
    result.maximum_energy_order = order + 2;
    for (const auto observable : observables)
        for (std::size_t i = 0; i < (is_fugacity(observable) ? components : 1); ++i)
            result.results.push_back({observable,
                is_fugacity(observable) ? std::optional<std::size_t>(i) : std::nullopt,
                std::nullopt, std::nullopt, status, order + (is_heat_capacity(observable) ? 2 : 1)});
    return result;
}

PropertyActionBatch evaluate_actions(
    ActionContext& context, const std::vector<std::size_t>& ids,
    const std::vector<PropertyObservable>& observables, const ActionPlan& plan,
    std::array<PropertyActionStatus, 2> support, PropertyActionStatus ideal_status,
    PropertyActionStatus association_status, bool thermal_boundary,
    bool temperature_action, bool pressure_chart, double pressure_slope
) {
    const bool association_available = association_status == PropertyActionStatus::Available;
    auto residual = plan.residual;
    if (!context.tape.residual_available) {
        residual.fill(false);
    } else if (!association_available) {
        residual[3] = false;
        residual[4] = false;
    }
    const auto jets = property_jets(context, ids, plan.ideal, residual, false,
        plan.fugacity && context.tape.residual_available && association_available);
    const bool zero_action = std::any_of(ids.begin(), ids.end(), [&](std::size_t id) {
        return zero_direction(context.directions[id]);
    });
    PropertyActionBatch result;
    result.results.reserve(observables.size());
    for (const auto observable : observables)
    for (std::size_t component = 0; component < (is_fugacity(observable) ? context.unit_ids.size() : 1); ++component) {
        const auto index = property_index(observable);
        const bool ideal = is_ideal(observable), total = is_total(observable), fugacity = is_fugacity(observable);
        PropertyActionResult item{observable, std::nullopt, std::nullopt, std::nullopt, PropertyActionStatus::Available,
            ids.size() + (is_heat_capacity(observable) ? 2 : 1)};
        if (fugacity) item.component = component;
        const Jet& residual_jet = fugacity ? jets.chemical_potential[component] : jets.residual[index];
        result.maximum_energy_order = std::max(result.maximum_energy_order, item.energy_order);
        if ((ideal || total) && !context.tape.ideal_available)
            item.status = ideal_status;
        else if ((ideal || total) && thermal_boundary && temperature_action)
            item.status = PropertyActionStatus::ThermalBoundary;
        else if (const auto s = support[(ideal || total) && index == 2 ? 0 : 1]; s != PropertyActionStatus::Available)
            item.status = s;
        else if ((!ideal || (pressure_chart && index == 2)) && !association_available)
            item.status = association_status;
        else if (!ideal && !context.tape.residual_available)
            item.status = PropertyActionStatus::CentralPropertyUnavailable;
        const double slope = pressure_chart ? pressure_slope
            : plan.pressure ? jets.pressure_slope : pressure_slope;
        if (item.status == PropertyActionStatus::Available
            && (pressure_chart || (!ideal && index == 3))
            && (!std::isfinite(slope) || slope <= 0.0))
            item.status = PropertyActionStatus::PressureSingular;

        const bool ideal_value = context.tape.ideal_available && (!total || context.tape.residual_available)
            && jets.ideal[index].size() != 0;
        const bool residual_value = context.tape.residual_available && residual_jet.size() != 0;
        if (ideal && ideal_value) item.value = jets.ideal[index].front();
        else if (!ideal && !total && residual_value) item.value = residual_jet.front();
        else if (total && ideal_value && residual_value)
            item.value = jets.ideal[index].front() + jets.residual[index].front();
        if (item.status == PropertyActionStatus::PressureSingular && index == 3)
            item.value.reset();

        if (item.status == PropertyActionStatus::Available && !zero_action) {
            if (ideal && ideal_value) item.action = jets.ideal[index].back();
            else if (!ideal && !total && residual_value) item.action = residual_jet.back();
            else if (total && ideal_value && residual_value)
                item.action = jets.ideal[index].back() + jets.residual[index].back();
        } else if (item.status == PropertyActionStatus::Available) {
            item.action = 0.0;
        }
        if (is_volume_partial(observable)) {
            const auto &volume = observable == PropertyObservable::MolarVolumeTemperatureDerivative
                                     ? jets.volume_temperature
                                     : jets.volume_pressure;
            item.value.reset();
            item.action.reset();
            if (!volume.empty()) {
                item.value = volume.front();
                if (item.status == PropertyActionStatus::Available)
                    item.action = zero_action ? 0.0 : volume.back();
            }
        }
        if (item.value.has_value() && !std::isfinite(*item.value)) {
            item.value.reset();
            item.action.reset();
            item.status = PropertyActionStatus::Nonfinite;
        }
        if (item.action.has_value() && !std::isfinite(*item.action)) {
            item.action.reset();
            item.status = PropertyActionStatus::Nonfinite;
        }
        result.results.push_back(std::move(item));
    }
    return result;
}
}  // namespace

void validate_cation_mixing(const Parameters& p) {
    const auto size = p.cation_mixing.size();
    require(size <= 1 && p.cation_pair_mixing.size() == size && p.cation_triple_mixing.size() == size, "cation mixing record");
    if (!size) return;
    const auto [s, i, j, k] = p.cation_mixing.front();
    for (const int index : {s, i, j, k}) require(index >= 0 && static_cast<std::size_t>(index) < p.component_count(), "cation mixing index");
    require(i < j && p.charge_number[s] == 0 && p.charge_number[i] > 0 && p.charge_number[j] > 0
        && p.charge_number[k] < 0 && p.electrolyte_form.declared(), "cation mixing participants");
    require(std::isfinite(p.molar_mass_kg_per_mol[s]) && p.molar_mass_kg_per_mol[s] > 0.0
        && finite_vector(p.cation_pair_mixing) && finite_vector(p.cation_triple_mixing), "cation mixing values");
}

HelmholtzModel::HelmholtzModel(
    Parameters parameters,
    std::vector<detail::ActiveParameter> active,
    std::optional<ThermochemistryRecord> thermochemistry, policy::SolverOptions solver_options
)
    : solver_options_(solver_options.resolved({})), parameters_(std::move(parameters)), active_(std::move(active)),
      thermochemistry_(std::move(thermochemistry)) {
    const std::size_t components = parameters_.component_count();
    const std::size_t sites = parameters_.site_count();
    const std::size_t pairs = components * components;
    const std::size_t site_pairs = sites * sites;
    if (parameters_.ideal_gas) {
        require(components > 0 && positive_vector(parameters_.molar_mass_kg_per_mol), "ideal-gas molar masses");
        require(parameters_.charge_number.size() == components && std::all_of(parameters_.charge_number.begin(),
            parameters_.charge_number.end(), [](int z) { return z == 0; }), "ideal-gas neutral species");
        require(active_.empty() && parameters_.segment_count.empty() && parameters_.sigma_constant.empty()
            && parameters_.dispersion_energy_over_k.empty() && sites == 0 && !parameters_.electrolyte_form.declared()
            && !parameters_.permittivity_form.declared(), "ideal-gas unconsumed residual inputs");
        if (thermochemistry_) validate_thermochemistry(*thermochemistry_, components);
        return;
    }
    require(components > 0, "segment_count");
    require(parameters_.sigma_constant.size() == components, "sigma_constant");
    require(parameters_.dispersion_energy_over_k.size() == components,
            "dispersion_energy_over_k");
    require(parameters_.molar_mass_kg_per_mol.size() == components, "molar_mass");
    require(parameters_.charge_number.size() == components, "charge_number");
    validate_cation_mixing(parameters_);
    require(parameters_.molecular_ion.empty() || parameters_.molecular_ion.size() == components, "charged_molecule_ion");
    require(parameters_.relative_permittivity.size() == components,
            "relative_permittivity");
    require(parameters_.packing_diameter_angstrom.size() == components,
            "packing_diameter");
    require(parameters_.debye_huckel_diameter_angstrom.size() == components,
            "debye_huckel_diameter");
    require(parameters_.born_diameter_angstrom.size() == components, "born_diameter");
    require(parameters_.solvation_factor.size() == components, "solvation_factor");
    require(parameters_.permittivity_reference_temperature_k.size() == components,
            "permittivity_reference_temperature");
    require(parameters_.ion_suppression_coefficient.size() == components,
            "ion_suppression_coefficient");
    require(parameters_.dielectric_volume_m3_per_mol.size() == components,
            "dielectric_volume");
    require(parameters_.dielectric_temperature_k.size() == components,
            "dielectric_temperature");
    require(parameters_.dipole_scaling_debye2.size() == components, "dipole_scaling");
    require(parameters_.polarizability_scaling_angstrom3.size() == components,
            "polarizability_scaling");
    require(parameters_.correlation_integral_parameter.size() == components,
            "correlation_integral");
    require(parameters_.k_ij.size() == pairs, "k_ij");
    require(parameters_.k_ij_reciprocal_temperature_slope.size() == pairs,
            "k_ij_temperature_slope");
    require(parameters_.k_ij_reference_temperature_k.size() == pairs,
            "k_ij_reference_temperature");
    require(parameters_.l_ij.size() == pairs, "l_ij");
    require(parameters_.site_multiplicity.size() == sites, "site_multiplicity");
    require(parameters_.association_active.size() == site_pairs, "association_active");
    require(parameters_.association_energy_over_k.size() == site_pairs,
            "association_energy_over_k");
    require(parameters_.association_volume.size() == site_pairs, "association_volume");
    const std::size_t sigma_terms = parameters_.sigma_term_component.size();
    require(parameters_.sigma_term_amplitude.size() == sigma_terms, "sigma_amplitude");
    require(parameters_.sigma_term_exponent.size() == sigma_terms, "sigma_exponent");
    const std::size_t permittivity_terms =
        parameters_.permittivity_term_component.size();
    require(parameters_.permittivity_term_coefficient.size() == permittivity_terms,
            "permittivity_coefficient");
    require(parameters_.permittivity_term_power.size() == permittivity_terms,
            "permittivity_power");

    require(positive_vector(parameters_.segment_count), "segment_count");
    require(positive_vector(parameters_.sigma_constant), "sigma_constant");
    require(std::all_of(
        parameters_.dispersion_energy_over_k.begin(),
        parameters_.dispersion_energy_over_k.end(),
        [](double value) { return std::isfinite(value) && value >= 0.0; }
    ), "dispersion_energy_over_k");
    for (std::size_t component = 0; component < components; ++component) {
        const double mass = parameters_.molar_mass_kg_per_mol[component];
        require(std::isfinite(mass) && mass >= 0.0, "molar_mass");
        require(parameters_.charge_number[component] != 0 || mass > 0.0, "molar_mass");
    }
    require(finite_vector(parameters_.k_ij), "k_ij");
    require(finite_vector(parameters_.l_ij), "l_ij");
    require(finite_vector(parameters_.k_ij_reciprocal_temperature_slope),
            "k_ij_temperature_slope");
    require(finite_vector(parameters_.k_ij_reference_temperature_k),
            "k_ij_reference_temperature");
    require(finite_vector(parameters_.sigma_term_amplitude), "sigma_amplitude");
    require(finite_vector(parameters_.sigma_term_exponent), "sigma_exponent");
    require(finite_vector(parameters_.molar_mass_kg_per_mol), "molar_mass");
    require(finite_vector(parameters_.relative_permittivity), "relative_permittivity");
    require(finite_vector(parameters_.permittivity_term_coefficient),
            "permittivity_coefficient");
    require(finite_vector(parameters_.permittivity_reference_temperature_k),
            "permittivity_reference_temperature");
    require(finite_vector(parameters_.ion_suppression_coefficient),
            "ion_suppression_coefficient");
    require(finite_vector(parameters_.dielectric_volume_m3_per_mol), "dielectric_volume");
    require(finite_vector(parameters_.dielectric_temperature_k), "dielectric_temperature");
    require(finite_vector(parameters_.dipole_scaling_debye2), "dipole_scaling");
    require(finite_vector(parameters_.polarizability_scaling_angstrom3),
            "polarizability_scaling");
    require(finite_vector(parameters_.correlation_integral_parameter),
            "correlation_integral");
    require(finite_vector(parameters_.association_energy_over_k),
            "association_energy_over_k");
    require(finite_vector(parameters_.association_volume), "association_volume");
    require(symmetric(parameters_.k_ij, components), "k_ij symmetry");
    require(symmetric(parameters_.l_ij, components), "l_ij symmetry");
    require(symmetric(parameters_.association_active, sites), "association_active symmetry");
    require(symmetric(parameters_.association_energy_over_k, sites),
            "association_energy_over_k symmetry");
    require(symmetric(parameters_.association_volume, sites),
            "association_volume symmetry");
    // build_parameters admits law structure; coefficients are the active values a fit may change.
    const std::size_t terms = parameters_.temperature_law_term_law.size();
    require(parameters_.temperature_law_coefficient.size() == terms && parameters_.temperature_law_term_power.size() == terms
        && finite_vector(parameters_.temperature_law_coefficient), "temperature_law_coefficient");
    for (const int component : parameters_.site_component) {
        require(component >= 0 && static_cast<std::size_t>(component) < components,
                "site_component");
    }
    for (const int active : parameters_.association_active) {
        require(active >= 0 && active <= 2, "association_active");
    }
    for (std::size_t pair = 0; pair < site_pairs; ++pair) {
        if (parameters_.association_active[pair] != 0) {
            require(parameters_.association_energy_over_k[pair] >= 0.0,
                    "association_energy_over_k");
            require(parameters_.association_volume[pair] >= 0.0,
                    "association_volume");
        }
    }
    for (const int component : parameters_.sigma_term_component) {
        require(component >= 0 && static_cast<std::size_t>(component) < components,
                "sigma_term_component");
    }
    for (const int component : parameters_.permittivity_term_component) {
        require(component >= 0 && static_cast<std::size_t>(component) < components,
                "permittivity_term_component");
    }
    for (std::size_t component = 0; component < components; ++component) {
        if (parameters_.charge_number[component] == 0) {
            continue;
        }
        require(molecular_ion(parameters_, component) || parameters_.packing_diameter_angstrom[component] > 0.0,
                "packing_diameter");
        require(parameters_.debye_huckel_diameter_angstrom[component] > 0.0,
                "debye_huckel_diameter");
        if (parameters_.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn)) {
            require(parameters_.born_diameter_angstrom[component] > 0.0,
                    "born_diameter");
        }
    }
    if (parameters_.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn)) {
        require(parameters_.ionic_region_relative_permittivity > 0.0,
                "ionic_region_relative_permittivity");
    }
    require(std::isfinite(parameters_.ionic_region_relative_permittivity),
            "ionic_region_relative_permittivity");
    require(std::isfinite(parameters_.born_c_shell), "born_c_shell");
    require(std::isfinite(parameters_.born_c_dielectric), "born_c_dielectric");
    // ComponentMixing stays defined on ions alone (Bulow 2019 Eq. 7). IonSuppression and DensityTemperature
    // average neutral-solvent records; Rueben 2024 fit DipolarPerturbation ion records to water-salt mixtures.
    if (parameters_.permittivity_form.declared()
        && !parameters_.permittivity_form.is(PermittivityForm::ComponentMixing)) {
        require(std::any_of(
            parameters_.charge_number.begin(), parameters_.charge_number.end(),
            [&](int charge) { return charge == 0; }
        ), "dielectric solvent");
    }
    const bool has_ion = std::any_of(
        parameters_.charge_number.begin(), parameters_.charge_number.end(),
        [](int charge) { return charge != 0; }
    );
    require(has_ion == parameters_.ionic_dispersion.declared(), "ionic_dispersion");
    require(has_ion == parameters_.electrolyte_form.declared(), "electrolyte_form");
    require(parameters_.electrolyte_form.declared() == parameters_.permittivity_form.declared(), "permittivity_form");
    if (thermochemistry_.has_value()) {
        validate_thermochemistry(*thermochemistry_, components);
    }

    std::set<std::pair<ParameterFamily, std::size_t>> seen;
    for (const auto& item : active_) {
        validate_parameter_coordinate(*this, item.family, item.indices, true);
        for (const auto index : item.indices) {
            require(seen.emplace(item.family, index).second,
                    "active_parameter duplicate");
        }
    }
}

void HelmholtzModel::validate_composition(
    const std::vector<double>& composition
) const {
    if (composition.size() != parameters_.component_count())
        throw std::invalid_argument("invalid composition size " + std::to_string(composition.size())
                                    + "; expected " + std::to_string(parameters_.component_count()));
    require(finite_vector(composition), "composition");
    require(std::all_of(composition.begin(), composition.end(), [](double value) {
        return value >= 0.0;
    }), "composition");
    double total = 0.0;
    for (double value : composition) total += value;
    if (!(std::abs(total - 1.0) <= policy::COMPOSITION_SUM_TOLERANCE.value))
        throw std::invalid_argument("invalid composition sum " + core::optional_number(total) + "; expected 1 within "
                                    + core::optional_number(policy::COMPOSITION_SUM_TOLERANCE.value));
}

void HelmholtzModel::validate_electroneutrality(const std::vector<double>& mole_fractions) const {
    validate_composition(mole_fractions);
    double charge = 0.0;
    for (std::size_t i = 0; i < mole_fractions.size(); ++i) charge += parameters_.charge_number[i] * mole_fractions[i];
    if (!(std::abs(charge) <= policy::COMPOSITION_SUM_TOLERANCE.value))
        throw std::invalid_argument("non-electroneutral composition: charge sum " + core::optional_number(charge)
                                    + "; expected 0 within " + core::optional_number(policy::COMPOSITION_SUM_TOLERANCE.value));
}

void HelmholtzModel::validate_point(
    double temperature,
    double volume,
    const std::vector<double>& amounts
) const {
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    require(std::isfinite(volume) && volume > 0.0, "volume");
    require(amounts.size() == parameters_.component_count(), "amounts size");
    require(finite_vector(amounts), "amounts");
    require(std::all_of(amounts.begin(), amounts.end(), [](double value) {
        return value >= 0.0;
    }), "amounts");
    double total = 0.0;
    for (double amount : amounts) total += amount;
    require(std::isfinite(total) && total > 0.0, "total amount");
}

detail::PredicateSignature HelmholtzModel::predicate_signature(
    double temperature, const std::vector<double>& composition, bool ideal
) const {
    detail::PredicateSignature result;
    for (const auto fraction : composition) result.push_back(fraction > 0.0);
    if (ideal && thermochemistry_.has_value()) {
        for (const auto& intervals : thermochemistry_->component_intervals) {
            for (const auto& interval : intervals)
                result.push_back(ideal_interval_contains(interval, temperature));
        }
    }
    return result;
}

HelmholtzDerivatives HelmholtzModel::evaluate(
    double temperature,
    double volume,
    const std::vector<double>& amounts,
    DerivativeOrder order
) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return evaluate_locked(temperature, volume, amounts, order);
}

HelmholtzDerivatives HelmholtzModel::evaluate_locked(
    double temperature,
    double volume,
    const std::vector<double>& amounts,
    DerivativeOrder order
) const {
    core::Scope entry("helmholtz_evaluate_locked"), scope("evaluate");
    entry.probe("requested_order", static_cast<double>(order));
    validate_point(temperature, volume, amounts);
    const std::size_t components = parameters_.component_count();
    const std::size_t point_size = 2 + components + active_.size();
    double total = 0.0;
    for (double amount : amounts) total += amount;
    const double density = total / volume;
    std::vector<double> composition;
    composition.reserve(components);
    for (double amount : amounts) composition.push_back(amount / total);
    const auto sigma = segment_sigma(temperature, parameters_);
    require(positive_vector(sigma), "sigma correlation");
    const auto diameter = segment_diameter(temperature, parameters_, sigma);
    require(positive_vector(diameter), "segment diameter");
    require_law_domain(temperature, parameters_);
    core::Convergence site_record;
    const auto sites = solve_site_fractions(
        temperature, density, composition, parameters_, sigma, diameter, &site_record, nullptr, nullptr, false, solver_options_
    );
    // The gradient floor propagates the achieved site residual, not the site tolerance.
    const double site_residual = std::max(site_record.residual, site_record.floor.value_or(0.0));
    scope.observe([&] {
        std::vector<double> declared_choices;
        if (parameters_.electrolyte_form.declared()) {
            declared_choices.push_back(0.0);
            declared_choices.push_back(static_cast<double>(
                parameters_.electrolyte_form.require("electrolyte_form")
            ));
        }
        if (parameters_.permittivity_form.declared()) {
            declared_choices.push_back(1.0);
            declared_choices.push_back(static_cast<double>(
                parameters_.permittivity_form.require("permittivity_form")
            ));
        }
        scope.probe("declared_choices", declared_choices);
    });
    const auto point = action_point(*this, temperature, volume, amounts);
    const auto signature = predicate_signature(temperature, composition);
    scope.observe([&] {
        scope.probe("predicate_signature", std::vector<double>(signature.begin(), signature.end()));
    });
    if (ensure_tape(point, sites, signature, {}, 0.0, true)) scope.probe("recording", 1.0);
    std::vector<double> full_point = point;
    full_point.insert(full_point.end(), sites.begin(), sites.end());
    const auto values = tape_->forward_zero(full_point, &counters_.forward);
    tape_->site_values = sites;
    HelmholtzDerivatives result;
    result.site_fractions = sites;
    result.value = values.front();
    scope.probe("residual_helmholtz_extensive", result.value);
    if (order == DerivativeOrder::Value) return result;

    const std::size_t sites_count = sites.size();
    const std::size_t input_count = point_size + sites_count;
    const std::size_t output_count = 1 + sites_count + (tape_->ideal_available ? 5 : 0);
    if (sites_count == 0 && order == DerivativeOrder::Gradient) {
        std::vector<double> weights(output_count, 0.0);
        weights.front() = 1.0;
        result.gradient = tape_->function.Reverse(1, weights);
        ++counters_.reverse;
        scope.probe("gradient", result.gradient);
        return result;
    }
    const auto delta = parameters_.ideal_gas ? std::vector<double>{} : association_delta(
        temperature, density, composition, parameters_, sigma, diameter);
    auto factorization = properties::site_factorization(
        sites, composition, parameters_, density * NUMBER_DENSITY_PER_ANGSTROM3, delta);
    if (sites_count != 0 && order == DerivativeOrder::PressureSlope) {
        // Density closure consumes the reduced gradient, a_VV and the floor only: one adjoint for
        // the site multipliers and one volume tangent replace the full site Jacobian.
        std::vector<double> weights(output_count, 0.0), direction(input_count, 0.0);
        weights.front() = 1.0;
        direction[1] = 1.0;
        const auto raw_gradient = tape_->function.Reverse(1, weights);
        const auto first = tape_->function.Forward(1, direction);
        Eigen::VectorXd phi_sites(static_cast<Eigen::Index>(sites_count)), volume_rows = phi_sites, cross = phi_sites;
        for (std::size_t site = 0; site < sites_count; ++site) {
            const auto s = static_cast<Eigen::Index>(site);
            phi_sites(s) = raw_gradient[point_size + site];
            volume_rows(s) = -first[1 + site];
            cross(s) = -total * composition[static_cast<std::size_t>(parameters_.site_component[site])]
                * parameters_.site_multiplicity[site] / sites[site] * first[1 + site];
        }
        const Eigen::VectorXd multipliers = factorization.transpose().solve(phi_sites);
        const Eigen::VectorXd volume_sites = factorization.solve(volume_rows);
        result.gradient_floor = std::abs(cross.dot(factorization.solve(
            Eigen::VectorXd::Ones(static_cast<Eigen::Index>(sites_count))))) * site_residual;
        for (std::size_t site = 0; site < sites_count; ++site) {
            weights[1 + site] = -multipliers(static_cast<Eigen::Index>(site));
            direction[point_size + site] = volume_sites(static_cast<Eigen::Index>(site));
        }
        tape_->function.Forward(1, direction);
        const auto second = tape_->function.Reverse(2, weights);
        counters_.forward += 2; counters_.reverse += 2;
        result.gradient.assign(point_size, 0.0);
        for (std::size_t row = 0; row < point_size; ++row) result.gradient[row] = second[2 * row];
        result.volume_second_derivative = second[2 * 1 + 1];
        for (std::size_t site = 0; site < sites_count; ++site)
            result.volume_second_derivative += volume_sites(static_cast<Eigen::Index>(site)) * second[2 * (point_size + site) + 1];
        scope.probe("gradient", result.gradient);
        scope.probe("volume_second_derivative", result.volume_second_derivative);
        scope.probe("gradient_floor", result.gradient_floor);
        return result;
    }
    // Requested chart: phase Hessians drop fixed controls, absent species and active parameters.
    std::vector<std::size_t> coordinates;
    if (order == DerivativeOrder::PressureSlope) {
        coordinates = {1};
    } else if (order == DerivativeOrder::PhaseHessian || order == DerivativeOrder::IsothermalPhaseHessian) {
        if (order == DerivativeOrder::PhaseHessian) coordinates.push_back(0);
        coordinates.push_back(1);
        for (std::size_t i = 0; i < components; ++i)
            if (amounts[i] > 0.0) coordinates.push_back(2 + i);
    } else {
        coordinates.resize(point_size);
        std::iota(coordinates.begin(), coordinates.end(), std::size_t{0});
    }
    const std::size_t width = coordinates.size();
    // The sites-free PressureSlope keeps every gradient row it has always returned.
    std::vector<std::size_t> rows = coordinates;
    if (order == DerivativeOrder::PressureSlope) {
        rows.resize(point_size);
        std::iota(rows.begin(), rows.end(), std::size_t{0});
    }
    std::vector<double> weights(output_count, 0.0);
    weights.front() = 1.0;
    Eigen::MatrixXd site_sensitivity(static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(width));
    result.gradient.assign(point_size, 0.0);
    if (sites_count != 0) {
        const auto raw_gradient = tape_->function.Reverse(1, weights);
        ++counters_.reverse;
        // Site partials of the requested columns from the cheaper of forward and reverse sweeps.
        Eigen::MatrixXd partials(static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(width));
        if (width <= sites_count) {
            for (std::size_t local = 0; local < width; ++local) {
                std::vector<double> direction(input_count, 0.0);
                direction[coordinates[local]] = 1.0;
                const auto first = tape_->function.Forward(1, direction);
                ++counters_.forward;
                for (std::size_t site = 0; site < sites_count; ++site)
                    partials(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(local)) = first[1 + site];
            }
        } else {
            for (std::size_t site = 0; site < sites_count; ++site) {
                std::vector<double> site_weight(output_count, 0.0);
                site_weight[1 + site] = 1.0;
                const auto partial = tape_->function.Reverse(1, site_weight);
                ++counters_.reverse;
                for (std::size_t local = 0; local < width; ++local)
                    partials(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(local)) = partial[coordinates[local]];
            }
        }
        site_sensitivity = factorization.solve(-partials);
        for (std::size_t local = 0; local < width; ++local) {
            const auto column = coordinates[local];
            result.gradient[column] = raw_gradient[column];
            for (std::size_t site = 0; site < sites_count; ++site)
                result.gradient[column] += raw_gradient[point_size + site]
                    * site_sensitivity(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(local));
        }
        scope.probe("gradient", result.gradient);
        if (order == DerivativeOrder::Gradient) return result;
        const Eigen::VectorXd multipliers = factorization.transpose().solve(
            Eigen::Map<const Eigen::VectorXd>(raw_gradient.data() + point_size, static_cast<Eigen::Index>(sites_count)));
        for (std::size_t site = 0; site < sites_count; ++site)
            weights[1 + site] = -multipliers(static_cast<Eigen::Index>(site));
    }
    result.hessian.assign(point_size * point_size, 0.0);
    if (order != DerivativeOrder::PressureSlope) ++counters_.hessian_calls;
    for (std::size_t offset = 0; offset < width; ++offset) {
        const std::size_t column = coordinates[offset];
        if (order == DerivativeOrder::StateHessian && column >= 2 && column < 2 + components) continue;
        std::vector<double> direction(input_count, 0.0);
        direction[column] = 1.0;
        for (std::size_t site = 0; site < sites_count; ++site)
            direction[point_size + site] =
                site_sensitivity(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(offset));
        tape_->function.Forward(1, direction);
        ++counters_.forward;
        const auto second = tape_->function.Reverse(2, weights);
        ++counters_.reverse;
        // Contract each row with that row coordinate's site tangent: O(sites) per entry.
        for (std::size_t local = 0; local < rows.size(); ++local) {
            const std::size_t row = rows[local];
            if (sites_count == 0 && offset == 0) result.gradient[row] = second[2 * row];
            double value = second[2 * row + 1];
            for (std::size_t site = 0; site < sites_count; ++site)
                value += site_sensitivity(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(local))
                    * second[2 * (point_size + site) + 1];
            result.hessian[row * point_size + column] = value;
        }
        if (sites_count == 0 && offset == 0) scope.probe("gradient", result.gradient);
    }
    if (order == DerivativeOrder::PressureSlope) {
        result.volume_second_derivative = result.hessian[1 * point_size + 1];
        scope.probe("volume_second_derivative", result.volume_second_derivative);
        scope.probe("gradient_floor", result.gradient_floor);
        return result;
    }
    if (sites_count != 0) {
        // Floor over the requested chart only: it bounds exactly the gradient entries computed here.
        double maximum = 0.0;
        for (std::size_t local = 0; local < width; ++local) {
            double norm = 0.0;
            for (std::size_t site = 0; site < sites_count; ++site) {
                const auto component = static_cast<std::size_t>(parameters_.site_component[site]);
                const double diagonal = total * composition[component]
                    * parameters_.site_multiplicity[site] / sites[site];
                norm += std::abs(site_sensitivity(static_cast<Eigen::Index>(site),
                    static_cast<Eigen::Index>(local)) * diagonal);
            }
            maximum = std::max(maximum, norm);
        }
        result.gradient_floor = maximum * site_residual;
    }
    result.volume_second_derivative = result.hessian[1 * point_size + 1];
    scope.probe("hessian", result.hessian);
    scope.probe("volume_second_derivative", result.volume_second_derivative);
    scope.probe("site_fractions", result.site_fractions);
    scope.probe("gradient_floor", result.gradient_floor);
    return result;
}

ResidualEnergyAction HelmholtzModel::residual_energy_action(double temperature, double volume,
    const std::vector<double>& amounts, const std::vector<std::vector<double>>& directions) const {
    validate_point(temperature, volume, amounts);
    const std::size_t components = parameters_.component_count(), order = directions.size();
    const std::size_t width = 2 + components + active_.size();
    require(order >= 3 && order <= 5, "energy action order");
    for (const auto& direction : directions)
        require(direction.size() == width && finite_vector(direction), "direction");
    double total = 0.0; for (double amount : amounts) total += amount;
    ResidualEnergyAction result;
    std::vector<ActionDirection> inserted;
    for (const auto& direction : directions)
        inserted.push_back({0.0, 0.0, {direction.begin() + 2, direction.begin() + 2 + static_cast<std::ptrdiff_t>(components)}, {}});
    if ((result.status = support_status(*this, amounts, inserted)[0]) != PropertyActionStatus::Available) return result;
    const auto point = action_point(*this, temperature, volume, amounts);
    std::vector<double> norms(order);
    std::vector<std::vector<double>> normalized(directions);
    for (std::size_t direction = 0; direction < order; ++direction) {
        for (std::size_t i = 0; i < width; ++i) {
            const double scale = i == 0 ? temperature : i == 1 ? volume
                : i < 2 + components ? total : std::max(1.0, std::abs(point[i]));
            norms[direction] = std::max(norms[direction], std::abs(directions[direction][i]) / scale);
        }
        for (double& value : normalized[direction]) value /= std::max(norms[direction], std::numeric_limits<double>::min());
    }
    std::lock_guard<std::mutex> lock(mutex_);
    const auto refreshed = evaluate_locked(temperature, volume, amounts, DerivativeOrder::Value);
    const double density = total / volume;
    std::vector<double> composition(amounts); for (double& amount : composition) amount /= total;
    result.diagnostics.recordings = ensure_tape(point, refreshed.site_fractions,
        predicate_signature(temperature, composition), {}, 0.0, true) ? 1 : 0;
    result.status = PropertyActionStatus::Available;
    const auto factorization = association_factorization(parameters_, temperature, density,
        refreshed.site_fractions, composition, result.status, result.diagnostics);
    if (result.status != PropertyActionStatus::Available) return result;
    if (std::any_of(norms.begin(), norms.end(), [](double norm) { return norm == 0.0; })) { result.action = 0.0; return result; }
    std::optional<std::pair<std::size_t, std::size_t>> repeated;
    if (order == 3) {
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = i + 1; j < 3; ++j)
                if (normalized[i] == normalized[j] && normalized[3 - i - j] != normalized[i])
                    repeated = {{3 - i - j, i}};
    }
    double coefficient = 0.0;
    if (repeated.has_value()) {
        // EqID: repeated_direction_energy_third
        // [t^5](F(p+v*t+u*t^3)-F(p+v*t-u*t^3)) = D3F[u,v,v].
        // The sparse paths avoid subtracting the singular pure-u third derivative
        // from an otherwise regular mixed response at dilute ionic states.
        for (const double sign : {1.0, -1.0}) {
            std::vector<std::vector<double>> series(6, std::vector<double>(width, 0.0));
            series[0] = point;
            series[1] = normalized[repeated->second];
            for (std::size_t i = 0; i < width; ++i)
                series[3][i] = sign * normalized[repeated->first][i];
            coefficient += sign * tape_->implicit_taylor(series, 5, factorization.get(),
                &result.diagnostics.site_solves, &result.diagnostics.forward_sweeps,
                &result.diagnostics.association_residual).front();
            ++result.diagnostics.directional_paths;
        }
    } else for (std::size_t mask = 0; mask < (std::size_t(1) << (order - 1)); ++mask) {
        // Half of the 2^q linear sign paths, weighted by their sign products.
        std::vector<std::vector<double>> series(order + 1, std::vector<double>(width, 0.0));
        series.front() = point;
        double weight = 1.0;
        for (std::size_t direction = 0; direction < order; ++direction) {
            const double sign = direction != 0 && (mask & (std::size_t(1) << (direction - 1))) ? -1.0 : 1.0;
            weight *= sign;
            for (std::size_t i = 0; i < width; ++i) series[1][i] += sign * normalized[direction][i];
        }
        coefficient += weight * tape_->implicit_taylor(series, order, factorization.get(), &result.diagnostics.site_solves,
            &result.diagnostics.forward_sweeps, &result.diagnostics.association_residual).front();
        ++result.diagnostics.directional_paths;
    }
    counters_.directional_paths += result.diagnostics.directional_paths; counters_.forward += result.diagnostics.forward_sweeps;
    counters_.maximum_energy_order = std::max(counters_.maximum_energy_order, order);
    double action = repeated.has_value() ? coefficient
        : coefficient / static_cast<double>(std::size_t(1) << (order - 1));
    for (const double norm : norms) action *= norm;
    if (result.diagnostics.association_residual > policy::ASSOCIATION_ACTION_RESIDUAL.value) result.status = PropertyActionStatus::AssociationResidual;
    else if (!std::isfinite(action)) result.status = PropertyActionStatus::Nonfinite;
    else result.action = action;
    return result;
}

ReferenceTaylor HelmholtzModel::reference_taylor(double temperature, double volume,
    const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    const PressureActionDirection& direction, std::size_t order, bool infinite_dilution) const {
    validate_point(temperature, volume, composition);
    validate_composition(composition);
    require(order <= 4 && stoichiometry.size() == composition.size()
        && finite_vector(stoichiometry) && std::isfinite(direction.temperature)
        && std::isfinite(direction.pressure) && direction.composition.empty(), "reference path");
    ReferenceTaylor result;
    std::vector<ActionDirection> insertion{{0.0, 0.0, stoichiometry, {}}};
    require(!infinite_dilution || !detail::charged_support(composition, parameters_), "ion-free reference solvent required");
    result.status = infinite_dilution ? PropertyActionStatus::Available : support_status(*this, composition, insertion)[0];
    if (result.status != PropertyActionStatus::Available) return result;
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<double> key{temperature, volume};
    key.insert(key.end(), composition.begin(), composition.end());
    auto found = std::find_if(centrals_.begin(), centrals_.end(), [&](const auto& entry) { return entry.first == key; });
    if (found == centrals_.end()) {
        if (centrals_.size() == 16) centrals_.pop_front();
        centrals_.emplace_back(std::move(key), evaluate_locked(temperature, volume, composition, DerivativeOrder::Hessian));
        found = std::prev(centrals_.end());
    }
    const auto& central = found->second;  // valid: centrals_ is not touched again while the lock is held
    const std::size_t width = central.gradient.size(), components = composition.size();
    const auto point = action_point(*this, temperature, volume, composition);
    std::vector<double> motion(width, 0.0);
    motion[0] = direction.temperature;
    for (const auto& parameter : direction.parameter) {
        require(std::isfinite(parameter.value), "reference parameter direction");
        const auto found = std::find_if(active_.begin(), active_.end(), [&](const auto& a) {
            return a.family == parameter.family && a.indices == parameter.indices;
        });
        require(found != active_.end(), "reference active parameter");
        motion[2 + components + std::size_t(found - active_.begin())] += parameter.value;
    }
    const double gas = GAS_CONSTANT_J_PER_MOL_K;
    const double pressure = gas * temperature * (1.0 / volume - central.gradient[1]);
    const double pv = gas * temperature * (-1.0 / (volume * volume) - central.hessian[width + 1]);
    result.diagnostics.pressure_density_derivative = -pv * volume * volume;
    if (!std::isfinite(pv) || pv >= 0.0) { result.status = PropertyActionStatus::PressureSingular; return result; }
    double norm = std::abs(direction.pressure) / (-pv * volume);
    for (std::size_t i = 0; i < width; ++i) {
        const double scale = i == 0 ? temperature : i == 1 ? volume
            : i < 2 + components ? 1.0 : std::max(1.0, std::abs(point[i]));
        norm = std::max(norm, std::abs(motion[i]) / scale);
    }
    if (norm == 0.0) norm = 1.0;
    for (double& v : motion) v /= norm;
    const auto signature = predicate_signature(temperature, composition);
    if (infinite_dilution) {
        const auto electro = electrolyte_contribution(temperature, 1.0 / volume, composition, composition, parameters_);
        (void)born_contribution(temperature, composition, parameters_, electro.bulk_relative_permittivity);
        if (!reference_tape_ || reference_tape_->predicates != signature) {
            reference_tape_ = detail::HelmholtzTape::record(parameters_, active_, point,
                central.site_fractions, signature, {}, 0.0, true, true);
            ++counters_.recordings;
            result.diagnostics.recordings = 1;
        }
    } else result.diagnostics.recordings = ensure_tape(point, central.site_fractions, signature, {}, 0.0, true) ? 1 : 0;
    auto& tape = infinite_dilution ? *reference_tape_ : *tape_;
    tape.site_values = central.site_fractions;
    const auto factorization = association_factorization(parameters_, temperature, 1.0 / volume,
        central.site_fractions, composition, result.status, result.diagnostics);
    if (result.status != PropertyActionStatus::Available) return result;
    const std::size_t insertion_degree = order + 1, degree = 2 * order + 1;
    std::vector<std::vector<double>> path(degree + 1, std::vector<double>(width, 0.0));
    path[0] = point; path[1] = motion;
    // EqID: reference_single_insertion_taylor. Only one insertion reaches degree 2r+1.
    const auto gradient_series = [&](const std::vector<double>& u, std::size_t last) {
        std::vector<double> gradient(last + 1, 0.0);
        for (double sign : {1.0, -1.0}) {
            auto inserted = path;
            for (std::size_t i = 0; i < width; ++i) inserted[insertion_degree][i] += sign * u[i];
            std::vector<std::vector<double>> outputs;
            tape.implicit_taylor(inserted, insertion_degree + last, factorization.get(),
                &result.diagnostics.site_solves, &result.diagnostics.forward_sweeps,
                &result.diagnostics.association_residual, &outputs);
            for (std::size_t k = 0; k <= last; ++k)
                gradient[k] += 0.5 * sign * outputs[insertion_degree + k][0];
            ++result.diagnostics.directional_paths;
        }
        return gradient;
    };
    const auto inverse_log = [&]() {
        std::pair<std::vector<double>, std::vector<double>> out;
        auto& [inverse, logarithm] = out;
        inverse.assign(order + 1, 0.0); logarithm.assign(order + 1, 0.0);
        inverse[0] = 1.0 / volume; logarithm[0] = std::log(volume);
        for (std::size_t k = 1; k <= order; ++k) {
            logarithm[k] = double(k) * path[k][1];
            for (std::size_t j = 1; j <= k; ++j) inverse[k] -= path[j][1] * inverse[k-j] / volume;
            for (std::size_t j = 1; j < k; ++j) logarithm[k] -= double(j) * logarithm[j] * path[k-j][1];
            logarithm[k] /= double(k) * volume;
        }
        return out;
    };
    std::vector<double> u(width, 0.0); u[1] = volume;
    result.diagnostics.pressure_factorizations = 1;
    for (std::size_t k = 1; k <= order; ++k) {
        const auto g = gradient_series(u, k);
        const auto inverse = inverse_log().first;
        const double known = gas * (temperature * (inverse[k] - g[k] / volume)
            + motion[0] * (inverse[k-1] - g[k-1] / volume));
        path[k][1] = ((k == 1 ? direction.pressure / norm : 0.0) - known) / pv;
        ++result.diagnostics.pressure_solves;
    }
    const auto g = gradient_series(u, order);
    const auto [inverse, logarithm] = inverse_log();
    for (std::size_t k = 1; k <= order; ++k) {
        const double residual = gas * (temperature * (inverse[k] - g[k] / volume)
            + motion[0] * (inverse[k-1] - g[k-1] / volume)) - (k == 1 ? direction.pressure / norm : 0.0);
        result.diagnostics.pressure_residual = std::max(result.diagnostics.pressure_residual,
            std::abs(residual) / std::max(pressure, std::abs(pv * path[k][1])));
    }
    u.assign(width, 0.0);
    std::copy(stoichiometry.begin(), stoichiometry.end(), u.begin() + 2);
    result.contraction = gradient_series(u, order);
    const double delta = std::accumulate(stoichiometry.begin(), stoichiometry.end(), 0.0);
    result.volume.resize(order + 1);
    for (std::size_t k = 0; k <= order; ++k) {
        result.contraction[k] = (result.contraction[k] - delta * logarithm[k]) * std::pow(norm, int(k));
        result.volume[k] = path[k][1] * std::pow(norm, int(k));
    }
    if (result.diagnostics.association_residual > policy::ASSOCIATION_ACTION_RESIDUAL.value)
        result.status = PropertyActionStatus::AssociationResidual;
    else if (!finite_vector(result.contraction) || !finite_vector(result.volume)) result.status = PropertyActionStatus::Nonfinite;
    counters_.directional_paths += result.diagnostics.directional_paths;
    counters_.forward += result.diagnostics.forward_sweeps;
    counters_.maximum_energy_order = std::max(counters_.maximum_energy_order, order + 1);
    return result;
}

// ponytail: retain two signatures for alternating phases; widen only for measured multi-signature churn.
bool HelmholtzModel::ensure_tape(const std::vector<double>& point, const std::vector<double>& sites,
    const detail::PredicateSignature& signature, const std::vector<IdealCorrelation>& correlations,
    double ideal_pressure, bool residual) const {
    const auto matches = [&](const auto& tape) {
        return tape != nullptr && tape->predicates == signature && tape->residual_available == residual
            && tape->ideal_available == !correlations.empty();
    };
    if (matches(tape_)) return false;
    if (matches(previous_tape_)) {
        tape_.swap(previous_tape_);
        return false;
    }
    previous_tape_ = detail::HelmholtzTape::record(parameters_, active_, point, sites, signature, correlations, ideal_pressure, residual);
    tape_.swap(previous_tape_);
    ++counters_.recordings;
    return true;
}

PropertyActionBatch HelmholtzModel::directional_actions(
    double temperature, double molar_density, const std::vector<double>& mole_fractions,
    const std::vector<PropertyObservable>& observables,
    const std::vector<DensityActionDirection>& directions,
    const DensityRoot* root, const std::vector<double>& pressure_directions, double pressure
) const {
    ActionPlan plan;
    for (const auto observable : observables) {
        const auto index = property_index(observable);
        plan.fugacity = plan.fugacity || is_fugacity(observable);
        if (is_fugacity(observable)) continue;
        if (is_ideal(observable) || is_total(observable)) plan.ideal[index] = true;
        if (!is_ideal(observable)) plan.residual[index] = true;
        if (!is_ideal(observable) && index == 3) plan.pressure = true;
    }
    require(root != nullptr || !plan.fugacity, "log fugacity coefficient needs the pressure chart");
    std::vector<ActionDirection> all;
    all.reserve(directions.size());
    for (const auto& direction : directions)
        all.push_back({direction.temperature, direction.molar_density,
            direction.composition, project_parameters(*this, direction.parameter)});
    const bool pressure_chain = root != nullptr && (observables.empty() || plan.ideal[2] || plan.fugacity
        || std::any_of(plan.residual.begin(), plan.residual.end(), [](bool value) { return value; }));
    const auto point = action_point(*this, temperature, 1.0 / molar_density, mole_fractions);
    // Linear polarization mixes each direction with the relative-T, relative-V and 1 mol labels.
    // Actions are multilinear, so scale every direction to those label sizes, on the pressure chart before and
    // after its density lift (dV/V ~ 1e-9 per Pa in a liquid), and rescale the results.
    std::vector<double> scales(directions.size(), 1.0);
    const auto normalize = [&](ActionDirection& d) {
        double s = std::max(std::abs(d.temperature) / temperature, std::abs(d.molar_density) / molar_density);
        for (const double v : d.composition) s = std::max(s, std::abs(v));
        for (std::size_t a = 0; a < d.parameter.size(); ++a)
            s = std::max(s, std::abs(d.parameter[a]) / std::max(std::abs(point[2 + mole_fractions.size() + a]), 1.0));
        if (!(s > 0.0) || !std::isfinite(s)) return 1.0;
        s = std::exp2(std::round(std::log2(s)));  // a power of two rescales exactly
        d.temperature /= s;
        d.molar_density /= s;
        for (double& v : d.composition) v /= s;
        for (double& v : d.parameter) v /= s;
        return s;
    };
    const bool need_ideal = std::any_of(plan.ideal.begin(), plan.ideal.end(), [](bool value) { return value; });
    bool need_residual = pressure_chain || plan.fugacity || std::any_of(plan.residual.begin(), plan.residual.end(), [](bool value) { return value; });
    bool endpoint = false;
    const auto correlations = need_ideal
        ? selected_ideal_correlations(thermochemistry_, temperature, endpoint)
        : std::vector<IdealCorrelation>{};
    const bool ideal_available = need_ideal && correlations.size() == parameters_.component_count();
    const bool record_missing = !thermochemistry_.has_value() || std::any_of(
        thermochemistry_->component_intervals.begin(), thermochemistry_->component_intervals.end(),
        [](const auto& intervals) { return intervals.empty(); });
    const auto ideal_status = !need_ideal || ideal_available ? PropertyActionStatus::Available
        : !record_missing ? PropertyActionStatus::OutsideInterval
        : PropertyActionStatus::MissingThermochemistry;
    PropertyActionStatus association_status = PropertyActionStatus::Available;
    std::vector<double> sites, sigma, diameter;
    if (need_residual) {
        sigma = segment_sigma(temperature, parameters_);
        require(positive_vector(sigma), "sigma correlation");
        diameter = segment_diameter(temperature, parameters_, sigma);
        require(positive_vector(diameter), "segment diameter");
        require_law_domain(temperature, parameters_);
        try {
            sites = solve_site_fractions(temperature, molar_density, mole_fractions,
                parameters_, sigma, diameter, nullptr, nullptr, nullptr, false, solver_options_);
        } catch (const core::ConvergenceFailure&) {
            need_residual = false;
            association_status = PropertyActionStatus::CentralPropertyUnavailable;
        }
    }
    const auto signature = predicate_signature(temperature, mole_fractions, ideal_available);
    PropertyActionBatch result;
    result.diagnostics.recordings = ensure_tape(point, sites, signature,
        ideal_available ? correlations : std::vector<IdealCorrelation>{},
        ideal_available ? thermochemistry_->standard_pressure : 0.0, need_residual) ? 1 : 0;
    tape_->site_values = sites;
    const auto factorization = need_residual ? association_factorization(parameters_, temperature, molar_density,
        tape_->site_values, mole_fractions, association_status, result.diagnostics) : nullptr;
    ActionContext context{*tape_, point, molar_density, all,
        association_status == PropertyActionStatus::Available ? factorization.get() : nullptr};
    context.maximum_order = directions.size()
        + (plan.residual[3] || plan.residual[4] ? 2 : 1);
    const auto result_support = support_status(*this, mole_fractions, all);
    // Direction subsets are bitmasks: base, lifted, radial and unit directions must fit one word.
    require(2 * directions.size() + 1 + (plan.fugacity ? mole_fractions.size() : 0) < 64, "direction count");
    for (std::size_t i = 0; plan.fugacity && i < mole_fractions.size(); ++i) {
        context.unit_ids.push_back(all.size());
        all.push_back({0.0, 0.0, std::vector<double>(mole_fractions.size()), std::vector<double>(active_.size())});
        all.back().composition[i] = 1.0;
    }
    const bool temperature_action = std::any_of(all.begin(), all.end(),
        [](const auto& direction) { return direction.temperature != 0.0; });
    auto diagnostics = result.diagnostics;
    for (std::size_t i = 0; i < directions.size(); ++i) scales[i] = normalize(all[i]);
    if (!pressure_chain || association_status != PropertyActionStatus::Available) {
        std::vector<std::size_t> ids(directions.size());
        for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = i;
        context.planning = true;
        (void)evaluate_actions(context, ids, observables, plan, result_support,
            ideal_status, association_status, endpoint, temperature_action,
            root != nullptr, root != nullptr ? root->pressure_density_derivative : 0.0);
        context.run();
        result = evaluate_actions(context, ids, observables, plan, result_support,
            ideal_status, association_status, endpoint, temperature_action,
            root != nullptr, root != nullptr ? root->pressure_density_derivative : 0.0);
        if (pressure_chain) result.density_status = association_status;
    } else {
        const double inverse_pressure_slope = 1.0 / root->pressure_density_derivative;
        std::size_t pressure_solves = 0;
        const auto density_solve = [&](double rhs) {
            ++pressure_solves;
            return rhs * inverse_pressure_slope;
        };
        const std::array<bool, 5> none{};
        const auto pressure_of = [&](const std::vector<std::size_t>& ids) {
            return property_jets(context, ids, none, none, true).pressure.back();
        };
        std::vector<std::size_t> base_ids(directions.size());
        for (std::size_t i = 0; i < base_ids.size(); ++i) base_ids[i] = i;
        context.planning = true;
        for (const auto id : base_ids) (void)pressure_of({id});
        context.run();
        std::vector<std::size_t> lifted_ids;
        std::vector<ActionDirection> lifted;
        lifted.reserve(directions.size());
        for (std::size_t i = 0; i < directions.size(); ++i) {
            const double density = density_solve(pressure_directions[i] / scales[i] - pressure_of({base_ids[i]}));
            lifted.push_back({all[i].temperature, density, all[i].composition, all[i].parameter});
            scales[i] *= normalize(lifted.back());
            all.push_back(lifted.back());
            lifted_ids.push_back(all.size() - 1);
        }
        const std::size_t radial_id = all.size();
        all.push_back({0.0, molar_density,
            std::vector<double>(mole_fractions.size()), std::vector<double>(active_.size())});
        const auto status = support_status(*this, mole_fractions, lifted);
        const bool temperature_action = std::any_of(lifted.begin(), lifted.end(),
            [](const auto& direction) { return direction.temperature != 0.0; });
        context.planning = true;
        (void)property_jets(context, lifted_ids, plan.ideal, plan.residual, true, plan.fugacity);
        if (directions.size() >= 2 && !observables.empty())
            (void)property_jets(context, {radial_id}, plan.ideal, plan.residual, false, plan.fugacity);
        if (directions.size() == 3 && !observables.empty()) {
            for (const auto id : lifted_ids)
                (void)property_jets(context, {id, radial_id}, plan.ideal, plan.residual, true, plan.fugacity);
        }
        context.run();
        result = evaluate_actions(context, lifted_ids, observables, plan, status,
            ideal_status, association_status, endpoint, temperature_action, true,
            root->pressure_density_derivative);
        result.molar_density = molar_density;
        result.density_status = status[1];
        for (const auto& direction : lifted) result.density_action.push_back(direction.molar_density);
        if (!finite_vector(result.density_action)) result.density_status = PropertyActionStatus::Nonfinite;
        std::array<double, 3> density_pairs{};
        if (directions.size() >= 2 && !observables.empty()) {
            density_pairs[0] = density_solve(-pressure_of({lifted_ids[0], lifted_ids[1]}));
        }
        if (directions.size() == 3 && !observables.empty()) {
            density_pairs[1] = density_solve(-pressure_of({lifted_ids[0], lifted_ids[2]}));
            density_pairs[2] = density_solve(-pressure_of({lifted_ids[1], lifted_ids[2]}));
        }
        PropertyActionBatch radial;
        std::vector<PropertyActionBatch> mixed;
        double density_123 = 0.0;
        if (directions.size() >= 2 && !observables.empty()) {
            radial = evaluate_actions(context, {radial_id}, observables, plan,
                {}, ideal_status, association_status,
                endpoint, false, true, root->pressure_density_derivative);
            if (directions.size() == 3) {
                const double pressure_123 = pressure_of(lifted_ids);
                density_123 = -pressure_123;
                for (std::size_t i = 0; i < 3; ++i) {
                    const std::array<std::size_t, 2> ids = {lifted_ids[i], radial_id};
                    const double pressure_ir = pressure_of({ids[0], ids[1]});
                    density_123 -= pressure_ir * density_pairs[2 - i]
                        / molar_density;
                    mixed.push_back(evaluate_actions(context, {ids[0], ids[1]}, observables,
                        plan, status, ideal_status, association_status, endpoint,
                        temperature_action, true, root->pressure_density_derivative));
                }
                density_123 = density_solve(density_123);
            }
            result.density_mixed_action = directions.size() == 2 ? density_pairs[0] : density_123;
        }
        diagnostics.pressure_factorizations = 1;
        diagnostics.pressure_solves = pressure_solves;
        for (std::size_t row = 0; row < result.results.size() && directions.size() >= 2; ++row) {
            if (result.results[row].status != PropertyActionStatus::Available) continue;
            const auto usable = [&](const PropertyActionBatch& batch) {
                return batch.results[row].status == PropertyActionStatus::Available
                    && batch.results[row].action.has_value();
            };
            if (!usable(radial)) {
                result.results[row].action.reset();
                result.results[row].status = radial.results[row].status;
                continue;
            }
            double action = *result.results[row].action
                + *radial.results[row].action * density_pairs[0] / molar_density;
            if (directions.size() == 3) {
                bool available = true;
                for (const auto& batch : mixed) available = available && usable(batch);
                if (!available) {
                    result.results[row].action.reset();
                    for (const auto& batch : mixed) {
                        if (!usable(batch)) { result.results[row].status = batch.results[row].status; break; }
                    }
                    continue;
                }
                action = *result.results[row].action + (
                    mixed[0].results[row].action.value() * density_pairs[2]
                    + mixed[1].results[row].action.value() * density_pairs[1]
                    + mixed[2].results[row].action.value() * density_pairs[0]
                    + radial.results[row].action.value() * density_123)
                    / molar_density;
            }
            result.results[row].action = std::isfinite(action)
                ? std::optional<double>(action) : std::nullopt;
            if (!result.results[row].action.has_value())
                result.results[row].status = PropertyActionStatus::Nonfinite;
        }
        // ln(phi_i) = mu_i - ln Z with Z = P / (R T rho) on the pressure-closed path: 1 - V phi_V
        // cancels in a dense liquid, while the total density actions above are well conditioned.
        Jet p(std::size_t(1) << directions.size(), 0.0), t = p, rho = p;
        p.front() = pressure; t.front() = temperature; rho.front() = molar_density;
        for (std::size_t i = 0; i < directions.size(); ++i) {
            p[std::size_t(1) << i] = pressure_directions[i] / scales[i];
            t[std::size_t(1) << i] = lifted[i].temperature;
            rho[std::size_t(1) << i] = lifted[i].molar_density;
        }
        for (std::size_t k = 0; k < directions.size() * (directions.size() - 1) / 2; ++k) rho[k == 0 ? 3 : 4 + k] = density_pairs[k];
        if (directions.size() == 3) rho[7] = density_123;
        const Jet log_z = jet_log(jet_multiply(p, jet_inverse(jet_multiply(t, rho))));
        for (auto& item : result.results)
            if (is_fugacity(item.observable) && item.action) *item.action -= log_z.back();
    }
    const bool charged = detail::charged_support(mole_fractions, parameters_);
    for (auto& item : result.results)  // fixed-pressure value, as state_at_pressure, which has no absent-ion value
        if (is_fugacity(item.observable) && ionic_insertion_value(parameters_, mole_fractions, *item.component, charged))
            item = {item.observable, item.component, std::nullopt, std::nullopt, PropertyActionStatus::IonicInsertion, item.energy_order};
        else if (is_fugacity(item.observable) && item.value)
            *item.value -= std::log(pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature * molar_density));
    const double product = std::accumulate(scales.begin(), scales.end(), 1.0, std::multiplies<>());
    for (auto& item : result.results) if (item.action) item.action = *item.action * product;
    for (std::size_t i = 0; i < result.density_action.size(); ++i) result.density_action[i] *= scales[i];
    if (result.density_mixed_action) result.density_mixed_action = *result.density_mixed_action * product;
    result.diagnostics = diagnostics;
    result.diagnostics.thermal_evaluations = ideal_available ? correlations.size() : 0;
    result.diagnostics.directional_paths = context.paths.size();
    result.diagnostics.forward_sweeps = context.forward_sweeps;
    result.diagnostics.site_solves = context.solves;
    result.diagnostics.association_residual = context.coefficient_residual;
    if (association_status == PropertyActionStatus::Available
        && context.coefficient_residual > policy::ASSOCIATION_ACTION_RESIDUAL.value) {
        if (pressure_chain) result.density_status = PropertyActionStatus::AssociationResidual;
        for (auto& item : result.results) {
            if ((!is_ideal(item.observable) || (pressure_chain && property_index(item.observable) == 2))
                && item.status == PropertyActionStatus::Available) {
                item.action.reset();
                item.status = PropertyActionStatus::AssociationResidual;
            }
        }
    }
    if (result.density_mixed_action && !std::isfinite(*result.density_mixed_action))
        result.density_status = PropertyActionStatus::Nonfinite;
    if (result.density_status != PropertyActionStatus::Available) {
        result.density_action.clear();
        result.density_mixed_action.reset();
    }
    if (root != nullptr) {
        result.diagnostics.pressure_density_derivative = root->pressure_density_derivative;
        result.diagnostics.stability_determinant = root->pressure_density_derivative
            / (GAS_CONSTANT_J_PER_MOL_K * temperature);
        result.diagnostics.pressure_residual = root->pressure_residual;
        result.diagnostics.branch = root->branch;
        result.diagnostics.stable_root_count = root->stable_root_count;
    }
    counters_.directional_paths += context.paths.size();
    counters_.forward += context.forward_sweeps;
    counters_.maximum_energy_order = std::max(counters_.maximum_energy_order,
        result.maximum_energy_order);
    result.solver_options = solver_options_;
    return result;
}

PropertyActionBatch HelmholtzModel::density_actions(
    double temperature, double molar_density, const std::vector<double>& mole_fractions,
    const std::vector<PropertyObservable>& observables,
    const std::vector<DensityActionDirection>& directions
) const {
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    require(std::isfinite(molar_density) && molar_density > 0.0, "molar_density");
    validate_composition(mole_fractions);
    require(!observables.empty(), "observables");
    validate_action_request(*this, observables, directions);
    std::lock_guard<std::mutex> lock(mutex_);
    return directional_actions(temperature, molar_density, mole_fractions, observables, directions);
}

PropertyActionBatch HelmholtzModel::pressure_actions(
    double temperature, double pressure, const std::vector<double>& mole_fractions,
    DensityBranch branch, const std::vector<PropertyObservable>& observables,
    const std::vector<PressureActionDirection>& directions,
    std::optional<double> anchor_density
) const {
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    require(std::isfinite(pressure) && pressure > 0.0, "pressure");
    validate_composition(mole_fractions);
    require(!directions.empty() && directions.size() <= 3, "action order");
    std::vector<DensityActionDirection> validation;
    validation.reserve(directions.size());
    std::vector<double> pressure_directions;
    for (const auto& direction : directions) {
        require(std::isfinite(direction.pressure), "pressure direction");
        pressure_directions.push_back(direction.pressure);
        validation.push_back({direction.temperature, 0.0, direction.composition, direction.parameter});
    }
    validate_action_request(*this, observables, validation);
    DensityRoot root;
    try {
        root = solve_density(*this, temperature, pressure, mole_fractions, branch, anchor_density);
    } catch (const std::runtime_error& error) {
        auto batch = unavailable_batch(observables, PropertyActionStatus::DensityUnavailable, directions.size(), mole_fractions.size());
        batch.density_status = PropertyActionStatus::DensityUnavailable;
        batch.solver_options = solver_options_;
        batch.density_message = error.what();
        batch.density_structure_refused = dynamic_cast<const DensityStructureError*>(&error) != nullptr;
        if (const auto* failure = dynamic_cast<const core::ConvergenceFailure*>(&error))
            batch.density_convergence = failure->record;
        return batch;
    }
    std::lock_guard<std::mutex> lock(mutex_);
    return directional_actions(temperature, root.molar_density, mole_fractions,
        observables, validation, &root, pressure_directions, pressure);
}

NeutralPropertyDerivatives HelmholtzModel::neutral_property(
    double temperature,
    double molar_density,
    NeutralProperty property,
    DerivativeOrder order
) const {
    core::Scope scope("neutral_property");
    std::lock_guard<std::mutex> lock(mutex_);
    require(order == DerivativeOrder::Value || order == DerivativeOrder::Gradient
        || order == DerivativeOrder::Hessian, "neutral property derivative order");
    require(
        property == NeutralProperty::ResidualHelmholtz || property == NeutralProperty::Pressure,
        "neutral property"
    );
    require(parameters_.component_count() == 1, "neutral property components");
    require(parameters_.site_count() == 0, "neutral property association sites");
    require(active_.empty(), "neutral property active parameters");
    require(parameters_.charge_number.front() == 0, "neutral property charge");
    require(!parameters_.electrolyte_form.declared(), "neutral property electrolyte");
    require(!parameters_.permittivity_form.declared(), "neutral property permittivity");
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    require(std::isfinite(molar_density) && molar_density > 0.0, "molar_density");

    const double volume = 1.0 / molar_density;
    const std::vector<double> amounts = {1.0};
    const auto energy_order = property == NeutralProperty::ResidualHelmholtz
        ? order
        : (order == DerivativeOrder::Value ? DerivativeOrder::Gradient : DerivativeOrder::Hessian);
    const auto energy = evaluate_locked(
        temperature, volume, amounts, energy_order
    );
    const double v1 = -1.0 / (molar_density * molar_density);
    const double v2 = 2.0 / (molar_density * molar_density * molar_density);

    NeutralPropertyDerivatives result;
    if (property == NeutralProperty::ResidualHelmholtz) {
        result.value = energy.value;
        if (order == DerivativeOrder::Value) return result;
        const auto& phi_gradient = energy.gradient;
        const double phi_t = phi_gradient[0];
        const double phi_v = phi_gradient[1];
        result.gradient = {phi_t, phi_v * v1};
        if (order == DerivativeOrder::Gradient) return result;
        const auto& phi_hessian = energy.hessian;
        const double phi_tt = phi_hessian[0];
        const double phi_tv = phi_hessian[1];
        const double phi_vv = phi_hessian[4];
        result.hessian = {
            phi_tt, phi_tv * v1,
            phi_tv * v1, phi_vv * v1 * v1 + phi_v * v2,
        };
        return result;
    }

    const auto& phi_gradient = energy.gradient;
    const auto& phi_hessian = energy.hessian;
    const double phi_v = phi_gradient[1];

    result.value = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (molar_density - phi_v);
    if (order == DerivativeOrder::Value) return result;
    const double phi_tv = phi_hessian[1];
    const double phi_vv = phi_hessian[4];
    const double p_t = GAS_CONSTANT_J_PER_MOL_K * (molar_density - phi_v)
        - GAS_CONSTANT_J_PER_MOL_K * temperature * phi_tv;
    const double p_rho = GAS_CONSTANT_J_PER_MOL_K * temperature * (1.0 - phi_vv * v1);
    result.gradient = {p_t, p_rho};
    if (order == DerivativeOrder::Gradient) return result;

    std::vector<double> point = {temperature, volume, 1.0};
    const std::vector<double> weight = {1.0};
    auto third = [&](std::size_t direction) {
        core::Scope scope("pressure_third_derivative"), evaluation("evaluate");
        scope.probe("direction", static_cast<double>(direction));
        tape_->forward_zero(point, &counters_.forward);
        std::vector<double> first(3, 0.0);
        first[direction] = 1.0;
        tape_->function.Forward(1, first);
        ++counters_.forward;
        tape_->function.Forward(2, std::vector<double>(3, 0.0));
        ++counters_.forward;
        const auto reverse = tape_->function.Reverse(3, weight);
        ++counters_.reverse;
        return reverse;
    };
    const auto third_t = third(0);
    const auto third_v = third(1);
    const double phi_ttv = 2.0 * third_t[5];
    const double phi_tvv = 2.0 * third_v[2];
    const double phi_vvv = 2.0 * third_v[5];
    const double p_tt = -2.0 * GAS_CONSTANT_J_PER_MOL_K * phi_tv
        - GAS_CONSTANT_J_PER_MOL_K * temperature * phi_ttv;
    const double p_t_rho = GAS_CONSTANT_J_PER_MOL_K * (1.0 - phi_vv * v1)
        - GAS_CONSTANT_J_PER_MOL_K * temperature * phi_tvv * v1;
    const double p_rho_rho = -GAS_CONSTANT_J_PER_MOL_K * temperature
        * (phi_vvv * v1 * v1 + phi_vv * v2);
    result.hessian = {p_tt, p_t_rho, p_t_rho, p_rho_rho};
    return result;
}

State HelmholtzModel::state(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    bool derived_properties
) const {
    core::Scope scope("state");
    validate_composition(mole_fractions);
    require(std::isfinite(molar_density) && molar_density > 0.0, "molar_density");
    const double volume = 1.0 / molar_density;
    const auto derivatives = evaluate(temperature, volume, mole_fractions,
        core::current_trace() ? DerivativeOrder::Hessian : DerivativeOrder::StateHessian);
    const auto& sites = derivatives.site_fractions;
    const auto [electro, terms] = [&] {
        core::Scope scope("state_residual_recompute"), evaluation("evaluate");
        const auto sigma = segment_sigma(temperature, parameters_);
        const auto diameter = segment_diameter(temperature, parameters_, sigma);
        const auto [electro, born] = detail::supported_electrostatics(
            temperature, molar_density,
            mole_fractions, mole_fractions, parameters_
        );
        const auto terms = residual_helmholtz(
            temperature, molar_density, mole_fractions, parameters_, sigma, diameter,
            sites, electro.debye_huckel, born
        );
        return std::pair{electro, terms};
    }();
    // EqID: residual_chemical_potential
    const double phi_v = derivatives.gradient[1];
    const double phi_t = derivatives.gradient[0];
    const double phi_vv = derivatives.hessian[1 * derivatives.gradient.size() + 1];
    const double phi_tv = derivatives.hessian[1];
    // EqID: compressibility_factor
    const double z = 1.0 - volume * phi_v;
    // EqID: pressure
    const double pressure = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (1.0 / volume - phi_v);
    const double dpressure_dt = pressure / temperature
        - GAS_CONSTANT_J_PER_MOL_K * temperature * phi_tv;
    const double dpressure_dv = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (-1.0 / (volume * volume) - phi_vv);
    const double dpressure_drho = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (1.0 + volume * volume * phi_vv);
    const double cv_over_r = -2.0 * temperature * phi_t
        - temperature * temperature * derivatives.hessian[0];
    State result;
    result.solver_options = solver_options_;
    result.hard_chain = terms.hard_chain;
    result.dispersion = terms.dispersion;
    result.association = terms.association;
    result.debye_huckel = terms.debye_huckel;
    result.born = terms.born;
    result.dipolar = terms.dipolar;
    result.cation_mixing = terms.cation_mixing;
    result.bulk_relative_permittivity = electro.bulk_relative_permittivity;
    result.residual_helmholtz = derivatives.value;
    result.compressibility_factor = z;
    result.pressure = pressure;
    result.residual_chemical_potential_over_rt.resize(mole_fractions.size());
    result.component_property_status.assign(mole_fractions.size(), PropertyActionStatus::Available);
    const bool charged = detail::charged_support(mole_fractions, parameters_);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        if (ionic_insertion_value(parameters_, mole_fractions, component, charged)) {
            result.component_property_status[component] = PropertyActionStatus::IonicInsertion;
            continue;
        }
        result.residual_chemical_potential_over_rt[component] =
            derivatives.gradient[2 + component];
    }
    if (z > 0.0 && std::isfinite(z)) {
        result.log_fugacity_coefficient.emplace(mole_fractions.size());
        for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
            if (!result.residual_chemical_potential_over_rt[component]) continue;
            (*result.log_fugacity_coefficient)[component] =
                *result.residual_chemical_potential_over_rt[component] - std::log(z);
        }
    }
    if (result.log_fugacity_coefficient && std::isfinite(dpressure_drho) && dpressure_drho > 0.0) {
        const std::size_t columns = active_.size(), width = derivatives.gradient.size();
        result.log_fugacity_parameter_derivative.resize(mole_fractions.size() * columns);
        result.density_parameter_derivative.resize(columns);
        for (std::size_t a = 0; a < columns; ++a) {
            // EqID: log_fugacity_coefficient_parameter_derivative
            const std::size_t index = 2 + mole_fractions.size() + a;
            const double va = -derivatives.hessian[width + index] / (1.0 / (volume * volume) + phi_vv);
            result.density_parameter_derivative[a] = -molar_density * molar_density * va;
            for (std::size_t i = 0; i < mole_fractions.size(); ++i) if (result.residual_chemical_potential_over_rt[i])
                result.log_fugacity_parameter_derivative[i * columns + a] = derivatives.hessian[(2 + i) * width + index]
                    + (derivatives.hessian[(2 + i) * width + 1] - 1.0 / volume) * va;
        }
        if (!active_.empty() && parameters_.permittivity_form.declared()) {
            std::vector<detail::AD> input{molar_density};
            const auto values = active_parameter_values();
            input.insert(input.end(), values.begin(), values.end());
            CppAD::Independent(input);
            const auto record = detail::make_tape_parameters(parameters_, active_, {input.begin() + 1, input.end()});
            const auto eps = bulk_relative_permittivity(detail::AD(temperature), input[0],
                std::vector<detail::AD>(mole_fractions.begin(), mole_fractions.end()), record);
            CppAD::ADFun<double> tape(input, std::vector<detail::AD>{eps});
            std::vector<double> point{molar_density};
            point.insert(point.end(), values.begin(), values.end());
            const auto jacobian = tape.Jacobian(point);
            for (std::size_t a = 0; a < columns; ++a)
                result.permittivity_parameter_derivative.push_back(jacobian[1 + a] + jacobian[0] * result.density_parameter_derivative[a]);
        }
    }
    // EqID: residual_entropy
    result.residual_entropy = GAS_CONSTANT_J_PER_MOL_K * (-derivatives.value - temperature * phi_t);
    // EqID: residual_internal_energy
    result.residual_internal_energy = -GAS_CONSTANT_J_PER_MOL_K * temperature * temperature * phi_t;
    // EqID: residual_enthalpy
    result.residual_enthalpy = GAS_CONSTANT_J_PER_MOL_K * temperature
        * (-temperature * phi_t + z - 1.0);
    // EqID: residual_heat_capacity_cv
    result.residual_isochoric_heat_capacity = GAS_CONSTANT_J_PER_MOL_K * cv_over_r;
    if (std::isfinite(dpressure_drho) && dpressure_drho > 0.0) {
        // EqID: residual_heat_capacity_cp
        const double cp = result.residual_isochoric_heat_capacity
            - temperature * dpressure_dt * dpressure_dt
                / dpressure_dv - GAS_CONSTANT_J_PER_MOL_K;
        if (std::isfinite(cp)) result.residual_isobaric_heat_capacity = cp;
    }
    if (thermochemistry_.has_value()) {
        const auto ideal = evaluate_ideal_properties(
            *thermochemistry_, temperature, molar_density, mole_fractions
        );
        result.thermochemistry_status = ideal.status;
        if (ideal.value.has_value()) {
            result.ideal_enthalpy = ideal.value->enthalpy;
            result.ideal_internal_energy = ideal.value->internal_energy;
            result.ideal_entropy = ideal.value->entropy;
            result.ideal_isobaric_heat_capacity = ideal.value->isobaric_heat_capacity;
            result.ideal_isochoric_heat_capacity = ideal.value->isochoric_heat_capacity;
            result.total_enthalpy = ideal.value->enthalpy + result.residual_enthalpy;
            result.total_internal_energy = ideal.value->internal_energy
                + result.residual_internal_energy;
            result.total_entropy = ideal.value->entropy + result.residual_entropy;
            // EqID: total_heat_capacity_cv
            result.total_isochoric_heat_capacity = ideal.value->isochoric_heat_capacity
                + result.residual_isochoric_heat_capacity;
            if (result.residual_isobaric_heat_capacity.has_value()) {
                // EqID: total_heat_capacity_cp
                result.total_isobaric_heat_capacity = ideal.value->isobaric_heat_capacity
                    + *result.residual_isobaric_heat_capacity;
            }
        }
    }
    result.molar_density = molar_density;
    result.pressure_density_derivative = dpressure_drho;
    result.site_fractions = sites;
    if (!derived_properties) return result;
    // Fixed-composition response and partial molar properties from the T/V Hessian columns above.
    const auto finite = [](double value) { return std::isfinite(value) ? std::optional(value) : std::nullopt; };
    const std::size_t components = mole_fractions.size(), width = derivatives.gradient.size();
    const bool slope = std::isfinite(dpressure_drho) && dpressure_drho > 0.0;
    const double volume_t = -dpressure_dt / dpressure_dv;  // (dv/dT)_P,x
    if (slope) {
        // EqID: isothermal_compressibility
        result.isothermal_compressibility = finite(1.0 / (molar_density * dpressure_drho));
        // EqID: isobaric_expansivity
        result.isobaric_expansivity = finite(dpressure_dt / (molar_density * dpressure_drho));
    }
    const auto& cp = result.total_isobaric_heat_capacity;  // available only with a positive slope
    const auto& cv = result.total_isochoric_heat_capacity;
    if (cp && *cp > 0.0) {
        // EqID: joule_thomson_coefficient
        result.joule_thomson_coefficient = finite((temperature * volume_t - volume) / *cp);
        if (cv && *cv > 0.0 && result.isothermal_compressibility) {
            // EqID: heat_capacity_ratio
            result.heat_capacity_ratio = finite(*cp / *cv);
            // EqID: isentropic_compressibility
            result.isentropic_compressibility = finite(*result.isothermal_compressibility * *cv / *cp);
            double mass = 0.0;
            bool masses = true;
            for (std::size_t i = 0; i < components; ++i) if (mole_fractions[i] > 0.0) {
                masses = masses && parameters_.molar_mass_kg_per_mol[i] > 0.0;
                mass += mole_fractions[i] * parameters_.molar_mass_kg_per_mol[i];
            }
            // EqID: speed_of_sound
            if (masses && result.heat_capacity_ratio)
                result.speed_of_sound = finite(std::sqrt(*result.heat_capacity_ratio * dpressure_drho / mass));
        }
    }
    auto& volumes = result.partial_molar_volume.emplace(components);
    auto& enthalpies = result.partial_molar_enthalpy.emplace(components);
    auto& entropies = result.partial_molar_entropy.emplace(components);
    const bool ideal = result.thermochemistry_status == ThermochemistryStatus::Available && result.ideal_enthalpy;
    for (std::size_t i = 0; slope && i < components; ++i) {
        if (!result.residual_chemical_potential_over_rt[i]) continue;  // IonicInsertion
        const double phi_nt = derivatives.hessian[(2 + i) * width], phi_nv = derivatives.hessian[(2 + i) * width + 1];
        // EqID: partial_molar_volume
        volumes[i] = finite(-GAS_CONSTANT_J_PER_MOL_K * temperature * (1.0 / volume - phi_nv) / dpressure_dv);
        if (!ideal) continue;
        const auto thermal = evaluate_ideal_correlation(
            ideal_interval(thermochemistry_->component_intervals[i], temperature)->correlation,
            temperature, thermochemistry_->standard_pressure);
        // d_i = (d ln phi_i / dT)_P,x on the density chart
        const double d = phi_nt + phi_nv * volume_t + 1.0 / temperature - volume_t / volume;
        // EqID: partial_molar_enthalpy
        enthalpies[i] = finite(thermal.h - GAS_CONSTANT_J_PER_MOL_K * temperature * temperature * d);
        if (mole_fractions[i] == 0.0) continue;  // ln x_i diverges
        const double log_ideal = std::log(mole_fractions[i] * molar_density * GAS_CONSTANT_J_PER_MOL_K * temperature
            / thermochemistry_->standard_pressure) + *result.residual_chemical_potential_over_rt[i];
        // EqID: partial_molar_entropy
        entropies[i] = finite(thermal.s - GAS_CONSTANT_J_PER_MOL_K * (log_ideal + temperature * d));
    }
    return result;
}

DensityBranch State::returned_density_branch() const {
    return properties::density_side(density_branch, molar_density, negative_slope_density);
}
bool State::density_branch_mismatch() const {
    const auto returned = returned_density_branch();
    return requested_density_branch != DensityBranch::Unspecified
        && returned != DensityBranch::Unspecified && requested_density_branch != returned;
}

State HelmholtzModel::state_at_pressure(
    double temperature,
    double pressure,
    const std::vector<double>& mole_fractions,
    DensityBranch branch,
    std::optional<double> anchor_density,
    bool derived_properties
) const {
    core::Scope scope("state_at_pressure");
    const auto root = solve_density(
        *this, temperature, pressure, mole_fractions, branch, anchor_density
    );
    State result = state(temperature, root.molar_density, mole_fractions, derived_properties);
    if (result.log_fugacity_coefficient) {
        // Fixed-pressure fugacity uses the prescribed P; EOS pressure retains the root residual.
        const double log_z = std::log(pressure / (GAS_CONSTANT_J_PER_MOL_K * temperature * root.molar_density));
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) if (result.residual_chemical_potential_over_rt[i])
            (*result.log_fugacity_coefficient)[i] = *result.residual_chemical_potential_over_rt[i] - log_z;
    }
    result.density_convergence = root.convergence;
    result.density_branch = root.branch;
    result.requested_density_branch = branch;
    result.stable_root_count = root.stable_root_count;
    result.negative_slope_density = root.negative_slope_density;
    scope.finish(root.convergence);
    return result;
}

double HelmholtzModel::packing_fraction_coefficient(
    double temperature,
    const std::vector<double>& mole_fractions
) const {
    validate_composition(mole_fractions);
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    if (parameters_.ideal_gas) return 0.0;
    const auto sigma = segment_sigma(temperature, parameters_);
    const auto diameter = segment_diameter(temperature, parameters_, sigma);
    require(positive_vector(sigma), "sigma correlation");
    require(positive_vector(diameter), "segment diameter");
    return packing_fraction_from_geometry(1.0, mole_fractions, parameters_, diameter);
}

double HelmholtzModel::molar_density_upper_bound(
    double temperature,
    const std::vector<double>& mole_fractions
) const {
    core::Scope scope("molar_density_upper_bound");
    const double coefficient = packing_fraction_coefficient(temperature, mole_fractions);
    require(std::isfinite(coefficient) && coefficient > 0.0, "packing coefficient");
    // EqID: density_upper_bound
    return policy::CLOSE_PACKING_FRACTION.value / coefficient;
}

std::vector<double> HelmholtzModel::active_parameter_values() const {
    std::vector<double> result;
    for (const auto& item : active_) result.push_back(detail::family_values(parameters_, item.family).first[item.indices.front()]);
    return result;
}

SweepCounters HelmholtzModel::counters() const {
    core::Scope scope("counters");
    std::lock_guard<std::mutex> lock(mutex_);
    return counters_;
}

void HelmholtzModel::reset_counters() const {
    core::Scope scope("reset_counters");
    std::lock_guard<std::mutex> lock(mutex_);
    counters_ = {};
}

}  // namespace epcsaft
