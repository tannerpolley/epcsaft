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
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace epcsaft {
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

bool symmetric(const std::vector<double>& values, std::size_t size) {
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = row + 1; column < size; ++column) {
            if (values[row * size + column] != values[column * size + row]) {
                return false;
            }
        }
    }
    return true;
}

bool symmetric(const std::vector<int>& values, std::size_t size) {
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = row + 1; column < size; ++column) {
            if (values[row * size + column] != values[column * size + row]) {
                return false;
            }
        }
    }
    return true;
}

Eigen::FullPivLU<Eigen::MatrixXd> site_factorization(
    const std::vector<double>& sites,
    const std::vector<double>& composition,
    const Parameters& parameters,
    double number_density,
    const std::vector<double>& delta,
    PropertyActionStatus* status = nullptr
) {
    const std::size_t sites_count = sites.size();
    Eigen::MatrixXd mass_action_jacobian = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(sites_count));
    for (std::size_t a = 0; a < sites_count; ++a) {
        double bonded = 0.0;
        for (std::size_t b = 0; b < sites_count; ++b) {
            const double mass_b = composition[parameters.site_component[b]]
                * parameters.site_multiplicity[b];
            const double delta_ab = delta[a * sites_count + b];
            bonded += mass_b * sites[b] * delta_ab;
            mass_action_jacobian(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) =
                sites[a] * number_density * mass_b * delta_ab;
        }
        mass_action_jacobian(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(a)) +=
            1.0 + number_density * bonded;
    }
    Eigen::FullPivLU<Eigen::MatrixXd> factorization(mass_action_jacobian);
    factorization.setThreshold(policy::ASSOCIATION_JACOBIAN_RANK.value * sites_count);
    if (sites_count != 0) {
        if (factorization.rank() < static_cast<Eigen::Index>(sites_count)) {
            if (status != nullptr) {
                *status = PropertyActionStatus::AssociationSingular;
                return factorization;
            }
            throw std::domain_error("association Jacobian is singular");
        }
        const auto pivots = factorization.matrixLU().diagonal().cwiseAbs();
        if (!std::isfinite(pivots.maxCoeff())
            || pivots.minCoeff() <= policy::ASSOCIATION_JACOBIAN_CONDITION.value
                * pivots.maxCoeff()) {
            if (status != nullptr) {
                *status = PropertyActionStatus::AssociationIllConditioned;
                return factorization;
            }
            throw std::runtime_error("association Jacobian is ill-conditioned");
        }
    }
    return factorization;
}

void require(bool condition, const char* field) {
    if (!condition) {
        throw std::invalid_argument(std::string("invalid ") + field);
    }
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
        result[1][0] = temperature + static_cast<double>(extra_temperature);
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
            result[target] = extra_volume
                ? value * std::pow(molar_density, static_cast<int>(extra_volume)) : value;
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
    Jet pressure;
    double pressure_slope = 0.0;
};

PropertyJets property_jets(
    ActionContext& context, const std::vector<std::size_t>& ids,
    const std::array<bool, 5>& need_ideal,
    const std::array<bool, 5>& need_residual, bool need_pressure
) {
    PropertyJets result;
    for (std::size_t i = 0; i < 5; ++i) {
        if (need_ideal[i] && context.tape.ideal_available)
            result.ideal[i] = context.extract(ids, 0, 0, context.tape.ideal_offset + i);
    }
    const bool h = need_residual[0], u = need_residual[1], s = need_residual[2];
    const bool cp = need_residual[3], cv = need_residual[4];
    if (!std::any_of(need_residual.begin(), need_residual.end(), [](bool value) { return value; })
        && !need_pressure) return result;
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
        result.residual[3] = jet_add(result.residual[4], jet_scale(
            jet_multiply(temperature, jet_multiply(jet_multiply(
                pressure_temperature, pressure_temperature), jet_inverse(pressure_volume))),
            -1.0));
        result.residual[3].front() -= gas_constant;
    }
    if (need_pressure) result.pressure = jet_scale(jet_multiply(temperature,
        jet_add(inverse_volume, phi_v, -1.0)), gas_constant);
    return result;
}

bool is_ideal(PropertyObservable observable) { return static_cast<std::size_t>(observable) < 5; }
bool is_total(PropertyObservable observable) { return static_cast<std::size_t>(observable) >= 10; }
std::size_t property_index(PropertyObservable observable) { return static_cast<std::size_t>(observable) % 5; }
bool is_heat_capacity(PropertyObservable observable) {
    const auto index = property_index(observable);
    return index == 3 || index == 4;
}
bool zero_direction(const ActionDirection& direction) {
    return direction.temperature == 0.0 && direction.molar_density == 0.0
        && std::all_of(direction.composition.begin(), direction.composition.end(),
            [](double value) { return value == 0.0; })
        && std::all_of(direction.parameter.begin(), direction.parameter.end(),
            [](double value) { return value == 0.0; });
}
const std::vector<double>& parameter_values(const Parameters& parameters, ParameterFamily family) {
    switch (family) {
    case ParameterFamily::SegmentCount: return parameters.segment_count;
    case ParameterFamily::SigmaConstant: return parameters.sigma_constant;
    case ParameterFamily::DispersionEnergyOverK: return parameters.dispersion_energy_over_k;
    case ParameterFamily::Kij: return parameters.k_ij;
    case ParameterFamily::Lij: return parameters.l_ij;
    case ParameterFamily::AssociationEnergyOverK: return parameters.association_energy_over_k;
    case ParameterFamily::AssociationVolume: return parameters.association_volume;
    }
    throw std::invalid_argument("invalid parameter family");
}
std::vector<double> action_point(
    const HelmholtzModel& model, double temperature, double rho,
    const std::vector<double>& composition
) {
    std::vector<double> point = {temperature, 1.0 / rho};
    point.insert(point.end(), composition.begin(), composition.end());
    for (const auto& item : model.active_parameters())
        point.push_back(parameter_values(model.parameters(), item.family)[item.indices.front()]);
    return point;
}

void validate_parameter_coordinate(
    const HelmholtzModel& model, ParameterFamily family,
    const std::vector<std::size_t>& indices, bool active
) {
    const auto& parameters = model.parameters();
    const auto& values = parameter_values(parameters, family);
    require(!indices.empty(), "parameter indices");
    for (const auto index : indices) require(index < values.size(), "parameter index");
    const bool association = family == ParameterFamily::AssociationEnergyOverK
        || family == ParameterFamily::AssociationVolume;
    const bool pair = association || family == ParameterFamily::Kij || family == ParameterFamily::Lij;
    if (!pair) {
        require(indices.size() == 1, "scalar parameter coordinate");
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
    if (active && association) {
        for (const auto index : indices)
            require(parameters.association_active[index] != 0, "inactive association coordinate");
    }
}

void validate_action_request(
    const HelmholtzModel& model,
    const std::vector<PropertyObservable>& observables,
    const std::vector<DensityActionDirection>& directions
) {
    require(!observables.empty(), "observables"); require(directions.size() >= 1 && directions.size() <= 3, "action order");
    std::vector<bool> seen(15, false);
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
    bool pressure = false;
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
        const IdealInterval* selected = nullptr;
        for (const auto& interval : intervals) {
            if (ideal_interval_contains(interval, temperature)) selected = &interval;
        }
        if (selected == nullptr) return {};
        endpoint = endpoint || temperature == selected->lower || temperature == selected->upper;
        result.push_back(selected->correlation);
    }
    return result;
}

PropertyActionStatus support_status(
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
    return support ? PropertyActionStatus::SupportChange
        : ionic ? PropertyActionStatus::IonicInsertion : PropertyActionStatus::Available;
}

PropertyActionBatch unavailable_batch(
    const std::vector<PropertyObservable>& observables, PropertyActionStatus status,
    std::size_t order
) {
    PropertyActionBatch result;
    result.maximum_energy_order = order + 2;
    for (const auto observable : observables) result.results.push_back({
        observable, std::nullopt, std::nullopt, status,
        order + (is_heat_capacity(observable) ? 2 : 1)});
    return result;
}

PropertyActionBatch evaluate_actions(
    ActionContext& context, const std::vector<std::size_t>& ids,
    const std::vector<PropertyObservable>& observables, const ActionPlan& plan,
    PropertyActionStatus support, PropertyActionStatus ideal_status,
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
    const auto jets = property_jets(context, ids, plan.ideal, residual, false);
    const bool zero_action = std::any_of(ids.begin(), ids.end(), [&](std::size_t id) {
        return zero_direction(context.directions[id]);
    });
    PropertyActionBatch result;
    result.results.reserve(observables.size());
    for (const auto observable : observables) {
        const auto index = property_index(observable);
        const bool ideal = is_ideal(observable), total = is_total(observable);
        PropertyActionResult item{
            observable, std::nullopt, std::nullopt, PropertyActionStatus::Available,
            ids.size() + (is_heat_capacity(observable) ? 2 : 1)};
        result.maximum_energy_order = std::max(result.maximum_energy_order, item.energy_order);
        if ((ideal || total) && !context.tape.ideal_available)
            item.status = ideal_status;
        else if ((ideal || total) && thermal_boundary && temperature_action)
            item.status = PropertyActionStatus::ThermalBoundary;
        else if (support != PropertyActionStatus::Available) item.status = support;
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
        const bool residual_value = context.tape.residual_available
            && jets.residual[index].size() != 0;
        if (ideal && ideal_value) item.value = jets.ideal[index].front();
        else if (!ideal && !total && residual_value) item.value = jets.residual[index].front();
        else if (total && ideal_value && residual_value)
            item.value = jets.ideal[index].front() + jets.residual[index].front();
        if (item.status == PropertyActionStatus::PressureSingular && index == 3)
            item.value.reset();

        if (item.status == PropertyActionStatus::Available && !zero_action) {
            if (ideal && ideal_value) item.action = jets.ideal[index].back();
            else if (!ideal && !total && residual_value) item.action = jets.residual[index].back();
            else if (total && ideal_value && residual_value)
                item.action = jets.ideal[index].back() + jets.residual[index].back();
        } else if (item.status == PropertyActionStatus::Available) {
            item.action = 0.0;
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

HelmholtzModel::HelmholtzModel(
    Parameters parameters,
    std::vector<detail::ActiveParameter> active,
    std::optional<ThermochemistryRecord> thermochemistry
)
    : parameters_(std::move(parameters)), active_(std::move(active)),
      thermochemistry_(std::move(thermochemistry)) {
    const std::size_t components = parameters_.component_count();
    const std::size_t sites = parameters_.site_count();
    const std::size_t pairs = components * components;
    const std::size_t site_pairs = sites * sites;
    require(components > 0, "segment_count");
    require(parameters_.sigma_constant.size() == components, "sigma_constant");
    require(parameters_.dispersion_energy_over_k.size() == components,
            "dispersion_energy_over_k");
    require(parameters_.molar_mass_kg_per_mol.size() == components, "molar_mass");
    require(parameters_.charge_number.size() == components, "charge_number");
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
    for (const int component : parameters_.site_component) {
        require(component >= 0 && static_cast<std::size_t>(component) < components,
                "site_component");
    }
    for (const int active : parameters_.association_active) {
        require(active == 0 || active == 1, "association_active");
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
        require(parameters_.packing_diameter_angstrom[component] > 0.0,
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
    if (parameters_.permittivity_form.declared()) {
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
    require(composition.size() == parameters_.component_count(), "composition size");
    require(finite_vector(composition), "composition");
    require(std::all_of(composition.begin(), composition.end(), [](double value) {
        return value >= 0.0;
    }), "composition");
    double total = 0.0;
    for (double value : composition) total += value;
    require(std::abs(total - 1.0) <= policy::COMPOSITION_SUM_TOLERANCE.value,
            "composition sum");
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
    double temperature,
    double density,
    const std::vector<double>& composition, bool residual, bool ideal
) const {
    detail::PredicateSignature result;
    for (const auto fraction : composition) result.push_back(fraction > 0.0);
    if (residual) {
    double charge_square_sum = 0.0;
    for (std::size_t component = 0; component < composition.size(); ++component) {
        if (!dielectric_solvent_component(component, parameters_)) {
            const double charge = parameters_.charge_number[component];
            charge_square_sum += composition[component] * charge * charge;
        }
    }
    result.push_back(charge_square_sum == 0.0 ? 1 : 0);
    if (charge_square_sum != 0.0) {
        const double permittivity = bulk_relative_permittivity(
            temperature, density, composition, parameters_
        );
        const double kappa = std::sqrt(
            density * AVOGADRO_PER_MOL * ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
            * charge_square_sum /
            (VACUUM_PERMITTIVITY_F_PER_M * permittivity
             * BOLTZMANN_J_PER_K * temperature)
        );
        for (std::size_t component = 0; component < composition.size(); ++component) {
            if (parameters_.charge_number[component] == 0) continue;
            result.push_back(
                kappa * ANGSTROM_M * parameters_.debye_huckel_diameter_angstrom[component]
                    < policy::SCREENING_SERIES_SWITCH.value ? 1 : 0
            );
        }
    }
    }
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
    core::Scope scope("evaluate");
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
    const auto sites = solve_site_fractions(
        temperature, density, composition, parameters_, sigma, diameter
    );
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
    std::vector<double> point = {temperature, volume};
    point.insert(point.end(), amounts.begin(), amounts.end());
    for (const auto& item : active_) {
        point.push_back(parameter_values(parameters_, item.family)[item.indices.front()]);
    }
    const auto signature = predicate_signature(temperature, density, composition);
    scope.probe(
        "predicate_signature",
        std::vector<double>(signature.begin(), signature.end())
    );
    if (tape_ == nullptr || tape_->predicates != signature
        || !tape_->residual_available) {
        tape_ = detail::HelmholtzTape::record(
            parameters_, active_, point, sites, signature, {}, 0.0, true
        );
        ++counters_.recordings;
        scope.probe("recording", 1.0);
    }
    std::vector<double> full_point = point;
    full_point.insert(full_point.end(), sites.begin(), sites.end());
    const auto values = tape_->function.Forward(0, full_point);
    ++counters_.forward;
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
    const auto jacobian = tape_->function.Jacobian(full_point);
    const auto delta = association_delta(
        temperature, density, composition, parameters_, sigma, diameter);
    auto factorization = site_factorization(
        sites, composition, parameters_, density * NUMBER_DENSITY_PER_ANGSTROM3, delta);
    Eigen::MatrixXd site_sensitivity = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(sites_count), static_cast<Eigen::Index>(point_size)
    );
    if (sites_count != 0) {
        for (std::size_t column = 0; column < point_size; ++column) {
            Eigen::VectorXd right_hand_side(static_cast<Eigen::Index>(sites_count));
            for (std::size_t row = 0; row < sites_count; ++row) {
                right_hand_side(static_cast<Eigen::Index>(row)) =
                    -jacobian[(1 + row) * input_count + column];
            }
            site_sensitivity.col(static_cast<Eigen::Index>(column)) =
                factorization.solve(right_hand_side);
        }
    }
    std::vector<double> gradient_weights(output_count, 0.0);
    gradient_weights.front() = 1.0;
    tape_->function.Forward(0, full_point);
    ++counters_.forward;
    const auto raw_gradient = tape_->function.Reverse(1, gradient_weights);
    ++counters_.reverse;
    result.gradient.resize(point_size);
    for (std::size_t column = 0; column < point_size; ++column) {
        result.gradient[column] = raw_gradient[column];
        for (std::size_t site = 0; site < sites_count; ++site) {
            result.gradient[column] += raw_gradient[point_size + site]
                * site_sensitivity(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(column));
        }
    }
    if (sites_count != 0 && order == DerivativeOrder::PressureSlope) {
        std::vector<double> volume_direction(input_count, 0.0);
        volume_direction[1] = 1.0;
        const auto first = tape_->function.Forward(1, volume_direction);
        ++counters_.forward;
        Eigen::VectorXd potential_cross(static_cast<Eigen::Index>(sites_count));
        for (std::size_t site = 0; site < sites_count; ++site) {
            const std::size_t component = static_cast<std::size_t>(
                parameters_.site_component[site]
            );
            const double mass = composition[component]
                * parameters_.site_multiplicity[site];
            potential_cross(static_cast<Eigen::Index>(site)) =
                -total * mass / sites[site] * first[1 + site];
        }
        const Eigen::VectorXd ones =
            Eigen::VectorXd::Ones(static_cast<Eigen::Index>(sites_count));
        result.gradient_floor = std::abs(
            potential_cross.dot(factorization.solve(ones))
        ) * policy::ASSOCIATION_RESIDUAL.value;
    }
    scope.probe("gradient", result.gradient);
    if (order == DerivativeOrder::Gradient) return result;

    std::vector<double> weights(output_count, 0.0);
    weights.front() = 1.0;
    if (sites_count != 0) {
        Eigen::VectorXd phi_sites(static_cast<Eigen::Index>(sites_count));
        for (std::size_t site = 0; site < sites_count; ++site) {
            phi_sites(static_cast<Eigen::Index>(site)) =
                jacobian[point_size + site];
        }
        const Eigen::VectorXd multipliers = factorization.transpose().solve(phi_sites);
        for (std::size_t site = 0; site < sites_count; ++site) {
            weights[1 + site] = -multipliers(static_cast<Eigen::Index>(site));
        }
    }
    const std::size_t requested_columns =
        order == DerivativeOrder::PressureSlope ? 1 : point_size;
    result.hessian.assign(point_size * point_size, 0.0);
    if (order != DerivativeOrder::PressureSlope) ++counters_.hessian_calls;
    for (std::size_t offset = 0; offset < requested_columns; ++offset) {
        const std::size_t column = order == DerivativeOrder::PressureSlope ? 1 : offset;
        std::vector<double> direction(input_count, 0.0);
        direction[column] = 1.0;
        for (std::size_t site = 0; site < sites_count; ++site) {
            direction[point_size + site] =
                site_sensitivity(static_cast<Eigen::Index>(site), static_cast<Eigen::Index>(column));
        }
        tape_->function.Forward(1, direction);
        ++counters_.forward;
        const auto second = tape_->function.Reverse(2, weights);
        ++counters_.reverse;
        for (std::size_t row = 0; row < point_size; ++row) {
            double value = 0.0;
            for (std::size_t input = 0; input < input_count; ++input) {
                double tangent = input == row ? 1.0 : 0.0;
                if (input >= point_size) {
                    tangent = site_sensitivity(
                        static_cast<Eigen::Index>(input - point_size),
                        static_cast<Eigen::Index>(row)
                    );
                }
                value += tangent * second[2 * input + 1];
            }
            result.hessian[row * point_size + column] = value;
        }
    }
    if (order == DerivativeOrder::PressureSlope) {
        result.volume_second_derivative = result.hessian[1 * point_size + 1];
        scope.probe("volume_second_derivative", result.volume_second_derivative);
        scope.probe("gradient_floor", result.gradient_floor);
        return result;
    }
    if (sites_count != 0) {
        double maximum = 0.0;
        for (std::size_t coordinate = 0; coordinate < point_size; ++coordinate) {
            double norm = 0.0;
            for (std::size_t site = 0; site < sites_count; ++site) {
                const auto component = static_cast<std::size_t>(parameters_.site_component[site]);
                const double diagonal = total * composition[component]
                    * parameters_.site_multiplicity[site] / sites[site];
                norm += std::abs(site_sensitivity(static_cast<Eigen::Index>(site),
                    static_cast<Eigen::Index>(coordinate)) * diagonal);
            }
            maximum = std::max(maximum, norm);
        }
        result.gradient_floor = maximum * policy::ASSOCIATION_RESIDUAL.value;
    }
    result.volume_second_derivative = result.hessian[1 * point_size + 1];
    scope.probe("hessian", result.hessian);
    scope.probe("volume_second_derivative", result.volume_second_derivative);
    scope.probe("site_fractions", result.site_fractions);
    scope.probe("gradient_floor", result.gradient_floor);
    return result;
}

PropertyActionBatch HelmholtzModel::directional_actions(
    double temperature, double molar_density, const std::vector<double>& mole_fractions,
    const std::vector<PropertyObservable>& observables,
    const std::vector<DensityActionDirection>& directions,
    const DensityRoot* root, const std::vector<double>& pressure_directions
) const {
    ActionPlan plan;
    for (const auto observable : observables) {
        const auto index = property_index(observable);
        if (is_ideal(observable) || is_total(observable)) plan.ideal[index] = true;
        if (!is_ideal(observable)) plan.residual[index] = true;
        if (!is_ideal(observable) && index == 3) plan.pressure = true;
    }
    std::vector<ActionDirection> all;
    all.reserve(directions.size());
    for (const auto& direction : directions)
        all.push_back({direction.temperature, direction.molar_density,
            direction.composition, project_parameters(*this, direction.parameter)});
    const bool pressure_chain = root != nullptr && (plan.ideal[2]
        || std::any_of(plan.residual.begin(), plan.residual.end(), [](bool value) { return value; }));
    const auto point = action_point(*this, temperature, molar_density, mole_fractions);
    const bool need_ideal = std::any_of(plan.ideal.begin(), plan.ideal.end(), [](bool value) { return value; });
    bool need_residual = pressure_chain || std::any_of(plan.residual.begin(), plan.residual.end(), [](bool value) { return value; });
    bool endpoint = false;
    const auto correlations = need_ideal
        ? selected_ideal_correlations(thermochemistry_, temperature, endpoint)
        : std::vector<IdealCorrelation>{};
    const bool ideal_available = need_ideal && correlations.size() == parameters_.component_count();
    const auto ideal_status = !need_ideal || ideal_available ? PropertyActionStatus::Available
        : thermochemistry_.has_value() ? PropertyActionStatus::OutsideInterval
        : PropertyActionStatus::MissingThermochemistry;
    PropertyActionStatus association_status = PropertyActionStatus::Available;
    std::vector<double> sites, sigma, diameter;
    if (need_residual) {
        sigma = segment_sigma(temperature, parameters_);
        require(positive_vector(sigma), "sigma correlation");
        diameter = segment_diameter(temperature, parameters_, sigma);
        require(positive_vector(diameter), "segment diameter");
        try {
            sites = solve_site_fractions(temperature, molar_density, mole_fractions,
                parameters_, sigma, diameter);
        } catch (const core::ConvergenceFailure&) {
            need_residual = false;
            association_status = PropertyActionStatus::CentralPropertyUnavailable;
        }
    }
    const auto signature = predicate_signature(temperature, molar_density, mole_fractions, need_residual, ideal_available);
    std::size_t recordings = 0;
    if (tape_ == nullptr || tape_->predicates != signature
        || tape_->residual_available != need_residual
        || tape_->ideal_available != ideal_available) {
        tape_ = detail::HelmholtzTape::record(
            parameters_, active_, point, sites, signature,
            ideal_available ? correlations : std::vector<IdealCorrelation>{},
            ideal_available ? thermochemistry_->standard_pressure : 0.0, need_residual);
        ++recordings;
        ++counters_.recordings;
    }
    tape_->site_values = sites;
    std::unique_ptr<Eigen::FullPivLU<Eigen::MatrixXd>> factorization;
    if (need_residual && tape_->site_count != 0) {
        factorization = std::make_unique<Eigen::FullPivLU<Eigen::MatrixXd>>(
            site_factorization(tape_->site_values, mole_fractions, parameters_,
                molar_density * NUMBER_DENSITY_PER_ANGSTROM3,
                association_delta(temperature, molar_density, mole_fractions,
                    parameters_, sigma, diameter), &association_status));
    }
    ActionContext context{*tape_, point, molar_density, all,
        association_status == PropertyActionStatus::Available ? factorization.get() : nullptr};
    context.maximum_order = directions.size()
        + (plan.residual[3] || plan.residual[4] ? 2 : 1);
    const auto result_support = support_status(*this, mole_fractions, all);
    const bool temperature_action = std::any_of(all.begin(), all.end(),
        [](const auto& direction) { return direction.temperature != 0.0; });
    PropertyActionBatch result;
    if (!pressure_chain || association_status != PropertyActionStatus::Available) {
        std::vector<std::size_t> ids(all.size());
        for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = i;
        context.planning = true;
        (void)evaluate_actions(context, ids, observables, plan, result_support,
            ideal_status, association_status, endpoint, temperature_action,
            root != nullptr, root != nullptr ? root->pressure_density_derivative : 0.0);
        context.run();
        result = evaluate_actions(context, ids, observables, plan, result_support,
            ideal_status, association_status, endpoint, temperature_action,
            root != nullptr, root != nullptr ? root->pressure_density_derivative : 0.0);
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
            const double density = density_solve(pressure_directions[i] - pressure_of({base_ids[i]}));
            lifted.push_back({directions[i].temperature, density,
                directions[i].composition, all[i].parameter});
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
        (void)property_jets(context, lifted_ids, plan.ideal, plan.residual, true);
        if (directions.size() >= 2)
            (void)property_jets(context, {radial_id}, plan.ideal, plan.residual, false);
        if (directions.size() == 3) {
            for (const auto id : lifted_ids)
                (void)property_jets(context, {id, radial_id}, plan.ideal, plan.residual, true);
        }
        context.run();
        result = evaluate_actions(context, lifted_ids, observables, plan, status,
            ideal_status, association_status, endpoint, temperature_action, true,
            root->pressure_density_derivative);
        std::array<double, 3> density_pairs{};
        if (directions.size() >= 2) {
            density_pairs[0] = density_solve(-pressure_of({lifted_ids[0], lifted_ids[1]}));
        }
        if (directions.size() == 3) {
            density_pairs[1] = density_solve(-pressure_of({lifted_ids[0], lifted_ids[2]}));
            density_pairs[2] = density_solve(-pressure_of({lifted_ids[1], lifted_ids[2]}));
        }
        PropertyActionBatch radial;
        std::vector<PropertyActionBatch> mixed;
        double density_123 = 0.0;
        if (directions.size() >= 2) {
            radial = evaluate_actions(context, {radial_id}, observables, plan,
                PropertyActionStatus::Available, ideal_status, association_status,
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
        }
        result.diagnostics.pressure_factorizations = 1;
        result.diagnostics.pressure_solves = pressure_solves;
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
    }
    result.diagnostics.recordings = recordings;
    result.diagnostics.thermal_evaluations = ideal_available ? correlations.size() : 0;
    result.diagnostics.directional_paths = context.paths.size();
    result.diagnostics.forward_sweeps = context.forward_sweeps;
    result.diagnostics.dense_hessian_calls = 0;
    result.diagnostics.site_factorizations = factorization ? 1 : 0;
    if (factorization) {
        result.diagnostics.site_rank = static_cast<std::size_t>(factorization->rank());
        const auto pivots = factorization->matrixLU().diagonal().cwiseAbs();
        if (pivots.size() != 0 && pivots.minCoeff() > 0.0)
            result.diagnostics.site_pivot_ratio = pivots.maxCoeff() / pivots.minCoeff();
    }
    result.diagnostics.site_solves = context.solves;
    result.diagnostics.association_residual = context.coefficient_residual;
    if (association_status == PropertyActionStatus::Available
        && context.coefficient_residual > policy::ASSOCIATION_ACTION_RESIDUAL.value) {
        for (auto& item : result.results) {
            if ((!is_ideal(item.observable) || (pressure_chain && property_index(item.observable) == 2))
                && item.status == PropertyActionStatus::Available) {
                item.action.reset();
                item.status = PropertyActionStatus::AssociationResidual;
            }
        }
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
    } catch (const std::runtime_error&) {
        return unavailable_batch(observables, PropertyActionStatus::DensityUnavailable, directions.size());
    }
    std::lock_guard<std::mutex> lock(mutex_);
    return directional_actions(temperature, root.molar_density, mole_fractions,
        observables, validation, &root, pressure_directions);
}

NeutralPropertyDerivatives HelmholtzModel::neutral_property(
    double temperature,
    double molar_density,
    NeutralProperty property,
    DerivativeOrder order
) const {
    core::Scope scope("neutral_property");
    std::lock_guard<std::mutex> lock(mutex_);
    require(order != DerivativeOrder::PressureSlope, "neutral property derivative order");
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
        tape_->function.Forward(0, point);
        ++counters_.forward;
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
    const std::vector<double>& mole_fractions
) const {
    core::Scope scope("state");
    validate_composition(mole_fractions);
    require(std::isfinite(molar_density) && molar_density > 0.0, "molar_density");
    const double volume = 1.0 / molar_density;
    std::vector<double> amounts = mole_fractions;
    const auto derivatives = evaluate(temperature, volume, amounts, DerivativeOrder::Hessian);
    const auto sigma = segment_sigma(temperature, parameters_);
    const auto diameter = segment_diameter(temperature, parameters_, sigma);
    const auto sites = solve_site_fractions(
        temperature, molar_density, mole_fractions, parameters_, sigma, diameter
    );
    const auto electro = electrolyte_contribution(
        temperature, molar_density, temperature, molar_density,
        mole_fractions, mole_fractions, parameters_
    );
    const auto born = born_contribution(
        temperature, mole_fractions, parameters_, electro.bulk_relative_permittivity
    );
    const auto terms = residual_helmholtz(
        temperature, molar_density, mole_fractions, parameters_, sigma, diameter,
        sites, electro.debye_huckel, born
    );
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
    result.hard_chain = terms.hard_chain;
    result.dispersion = terms.dispersion;
    result.association = terms.association;
    result.debye_huckel = terms.debye_huckel;
    result.born = terms.born;
    result.bulk_relative_permittivity = electro.bulk_relative_permittivity;
    result.residual_helmholtz = derivatives.value;
    result.compressibility_factor = z;
    result.pressure = pressure;
    result.residual_chemical_potential_over_rt.resize(mole_fractions.size());
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        result.residual_chemical_potential_over_rt[component] =
            derivatives.gradient[2 + component];
    }
    if (z > 0.0 && std::isfinite(z)) {
        result.log_fugacity_coefficient.emplace(mole_fractions.size());
        for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
            (*result.log_fugacity_coefficient)[component] =
                result.residual_chemical_potential_over_rt[component] - std::log(z);
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
            result.total_isochoric_heat_capacity = ideal.value->isochoric_heat_capacity
                + result.residual_isochoric_heat_capacity;
            if (result.residual_isobaric_heat_capacity.has_value()) {
                result.total_isobaric_heat_capacity = ideal.value->isobaric_heat_capacity
                    + *result.residual_isobaric_heat_capacity;
            }
        }
    }
    result.molar_density = molar_density;
    result.pressure_density_derivative = dpressure_drho;
    result.site_fractions = sites;
    return result;
}

State HelmholtzModel::state_at_pressure(
    double temperature,
    double pressure,
    const std::vector<double>& mole_fractions,
    DensityBranch branch,
    std::optional<double> anchor_density
) const {
    core::Scope scope("state_at_pressure");
    const auto root = solve_density(
        *this, temperature, pressure, mole_fractions, branch, anchor_density
    );
    State result = state(temperature, root.molar_density, mole_fractions);
    result.density_convergence = root.convergence;
    result.density_branch = root.branch;
    result.stable_root_count = root.stable_root_count;
    scope.finish(root.convergence);
    return result;
}

double HelmholtzModel::packing_fraction_coefficient(
    double temperature,
    const std::vector<double>& mole_fractions
) const {
    validate_composition(mole_fractions);
    require(std::isfinite(temperature) && temperature > 0.0, "temperature");
    const auto sigma = segment_sigma(temperature, parameters_);
    const auto diameter = segment_diameter(temperature, parameters_, sigma);
    require(positive_vector(sigma), "sigma correlation");
    require(positive_vector(diameter), "segment diameter");
    return packing_fraction_from_geometry(1.0, mole_fractions, parameters_, diameter);
}

double HelmholtzModel::phase_packing_fraction(
    double temperature,
    const std::vector<double>& amounts,
    double volume
) const {
    validate_point(temperature, volume, amounts);
    double total = 0.0;
    for (double amount : amounts) total += amount;
    std::vector<double> composition;
    composition.reserve(amounts.size());
    for (double amount : amounts) composition.push_back(amount / total);
    const auto sigma = segment_sigma(temperature, parameters_);
    const auto diameter = segment_diameter(temperature, parameters_, sigma);
    return packing_fraction_from_geometry(total / volume, composition, parameters_, diameter);
}

double HelmholtzModel::molar_density_upper_bound(
    double temperature,
    const std::vector<double>& mole_fractions
) const {
    const double coefficient = packing_fraction_coefficient(temperature, mole_fractions);
    require(std::isfinite(coefficient) && coefficient > 0.0, "packing coefficient");
    return policy::CLOSE_PACKING_FRACTION.value / coefficient;
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
